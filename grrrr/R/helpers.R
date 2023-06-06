library(purrr)

# this is implementation of
# b_{-j}=(b^'_1, ..., b^'_{j-1}, 0', b^'_{j+1}, ..., b^'_{J})

#' A easier form of setting part of array to zero
#'
#' @param vector array which is meant to be used
#' @param index indexes were zeros will be inserted
#' 
#' @description
#' Implementation of b_{-j}=(b^'_1, ..., b^'_{j-1}, 0', b^'_{j+1}, ..., b^'_{J})
#' 
#'
#' @return array with zeros in selected indexes
`%-%` <- function(vector, index) {
    vector[index] <- 0
    return(vector)
}


#' Vector norm mentioned in the article
#'
#' @param vector array
#' @param p multiplier of identity matrix
#'
#' @return vector norm
norm_L <- function(vector, p) {
    K <- diag(p)
    norm_ <- sqrt(t(vector) %*% K %*% vector)[1]
    return(norm_)
}

#' Calculation of Cp value
#'
#' @param indexes array with factors chosen in the model
#' @param group_sizes array with sizes of consecutive groups
#' @param betas beta coefficients in the investigated model
#' @param betas_ls beta coefficients in the OLS model build on the same data 
#' @param X matrix of regressors
#' @param y target variable
#' @param df_function function that calculates degrees of freedom for specific model
#'
#' @return value of Cp statistic
calculate_cp <- function(indexes, group_sizes,
                         betas, betas_ls, X, y, df_function) {
    dg_f <- df_function(indexes, group_sizes, betas, betas_ls)

    mu <- X %*% betas

    Cp <- (sum((y - mu)**2) / var(y)[1]) - nrow(X) + 2 * dg_f
    return(Cp)
}

#' Calculation of model error value
#'
#' @param X matrix of regressors
#' @param beta_hat beta coefficients in the investigated model
#' @param beta original beta coefficients used to generate data set
#'
#' @return value of model error 
calculate_me <- function(X, beta_hat, beta) {
    me <- t((beta_hat - beta)) %*% t(X) %*% X %*% (beta_hat - beta)
    me <- t((beta_hat - beta)) %*% (beta_hat - beta)
    return(me[[1]])
}

#' Makes first letter of string uppercase
#'
#' @param x string
#'
#' @return transformed string
first_up <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}
