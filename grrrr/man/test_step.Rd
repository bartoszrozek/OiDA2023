% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testing.R
\name{test_step}
\alias{test_step}
\title{Function that test stepwise regression method}
\usage{
test_step(X, y, true_betas, groups)
}
\arguments{
\item{X}{design matrix}

\item{y}{target variable}

\item{true_betas}{beta coefficients used in Y calculation}

\item{groups}{added only to keep function's shape}
}
\value{
instance of test_result class
}
\description{
Function that test stepwise regression method
}
