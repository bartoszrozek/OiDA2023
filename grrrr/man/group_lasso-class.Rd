% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\docType{class}
\name{group_lasso-class}
\alias{group_lasso-class}
\title{Class storing information about group lasso model}
\value{
instance of group_lasso class
}
\description{
Class storing information about group lasso model
}
\section{Slots}{

\describe{
\item{\code{X}}{matrix. Design matrix}

\item{\code{y}}{numeric. Target variable}

\item{\code{betas}}{numeric. Final beta coefficients}

\item{\code{betas_path}}{list. List of all beta coefficients obtain during calculations}

\item{\code{true_betas}}{numericOrNULL. Beta coefficients used in target variable calculations}

\item{\code{lambda_max}}{numeric. Maximum value of lambda}

\item{\code{lambda_best}}{numeric. Value of lambda used for final model}

\item{\code{Cp}}{numeric. Value of Cp}

\item{\code{Cp_path}}{list. List of values of Cp obtained during calculations}

\item{\code{model_error}}{numericOrNULL. Value of model_error for final model. Not null only if 
true_betas was supplied.}

\item{\code{me_path}}{listOrNULL. List of values of model_error obtained during calculations. Not null only if 
true_betas was supplied.}
}}

