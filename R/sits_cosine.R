# sits_cosine.R
# Functions for fitting stable isotope time series with the cosine model

#' @title Prepare results for cosine model fit.
#' 
#' @description Performs the nonlinear least squares (NLS) regression method for the cosine 
#' model, with the proposed initialization for all the parameters. It fits the NLS method 
#' as required, and then computes different quantities for the birth seasonality estimates
#' corresponding to different individuals.
#' 
#' @param paths A list of data frames, where each frame contains the data for one individual. Every 
#' data frame should have two columns with names 'distance' and 'oxygen'.
#' 
#' @return A data frame with the estimated parameters in the model, birth seasonality estimate, 
#' predicted/observed minimum/maximum for the oxygen isotope variable, mean squared error
#' and Pearson's R^2 corresponding to the model fit for every individual.

makeFits = function(paths) {
  
  fits = c()
  for(i in 1:length(paths)) {
    data = paths[[i]]
    curve = sineFit(data)
    fit = convertParameters(curve)
    fit$predictedMin = fit$intercept - abs(fit$amplitude)
    fit$predictedMax = fit$intercept + abs(fit$amplitude)
    fit$observedMin = min(data$oxygen)
    fit$observedMax = max(data$oxygen)
    fit$MSE = mean((predict(curve) - data$oxygen)^2)
    fit$PearsonCorrelation = cor(predict(curve),paths[[i]]$oxygen,method = "pearson")
    fits = rbind(fits, as.numeric(fit))
  }
  fits = data.frame(fits)
  colnames(fits) = c("amplitude","intercept","x0","X","birth","predictedMin",
                     "predictedMax","observedMin","observedMax","MSE","Pearson")
  
  return(fits)
  
}

#' @title Fit Cosine Model using Nonlinear Least Squares
#' 
#' @description
#' Performs the updated nonlinear least squares (NLS) regression method for the cosine 
#' model proposed by Balasse et al. The method calculates with the proposed initial values at first, 
#' and then fits the NLS method as required.
#' 
#' @param data      A data frame containing the stable isotope values for a single specimen.
#' @param formula   Formula describing the independent and dependent variables, e.g. `oxygen~distance`. Will be converted to a cosine model.
#' 
#' @return 
#' A fitted model object. See [stats::nls()].
#' 
#' @examples
#' data(tsaghkahovit)
#' sits_cosine(tsaghkahovit[tsaghkahovit$id==6223,], oxygen~distance)
#' 
#' @export
sits_cosine <- function(data, formula) {
  # Rewrite formula as a sine/cosine model
  environment(formula) <- environment()
  formula %>% 
    formula.tools::rhs() %>% 
    rlang::eval_tidy(data) %>%
    {2 * pi / max(.)} ->
    frequency
  lm_model <- update.formula(formula, . ~ sin(frequency*.) + cos(frequency*.))
  
  # Fit linear model to provide starting estimates for NLS
  fit <- lm(lm_model, data = data) 
  coefs <- list(intercept = fit$coef[[1]],
                amplitude = sqrt(fit$coef[[2]]^2 + fit$coef[[3]]^2),
                phase = -atan(fit$coef[[2]] / fit$coef[[3]]),
                frequency = frequency)
  
  # Fit with Nonlinear Least Squares regression
  nls_model <- update.formula(formula, . ~ intercept + I(amplitude * cos(frequency*.+phase)))
  nls(nls_model, data = data, start = coefs, 
      control = nls.control(warnOnly = TRUE)) %>% 
    return()
}
