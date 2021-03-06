% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sineFitWrong.R
\name{sineFitWrong}
\alias{sineFitWrong}
\title{Cosine model fitting with given initialization for two parameters.}
\usage{
sineFitWrong(data, amplitude, intercept)
}
\arguments{
\item{data}{A data frame that contains the data for one individual. There should be two columns
with names 'distance' and 'oxygen'.}

\item{amplitude}{A number, corresponding to a given value for the amplitude parameter.}

\item{intercept}{A number, corresponding to a given value for the intercept parameter.}
}
\value{
A fitted model object from the nls function in R.
}
\description{
Performs the updated nonlinear least squares (NLS) regression method for the cosine
model proposed by Balasse et al. The method calculates with the proposed initial values at first,
and then fits the NLS method as required.
}
