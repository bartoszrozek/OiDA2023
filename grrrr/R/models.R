library(MASS)
library(purrr)
library(pracma)

#' Trichotomization of values in the matrix
#'
#' @param Z matrix
#'
#' @return matrix with trichomizated values
categorize_matrix <- function(Z) {
    Z[Z > qnorm(2 / 3)] <- 2
    Z[Z > qnorm(1 / 3) & Z < qnorm(2 / 3)] <- 1
    Z[Z < qnorm(1 / 3)] <- 0
    return(Z)
}

#' Generating noise for target variable
#'
#' @param Y array
#' @param ratio signal-to-noise-ratio 
#'
#' @return array with noise
generate_noise <- function(Y, ratio) {
    Y_var <- var(Y)
    Y_noise <- rnorm(
        n = length(Y),
        mean = 0,
        sd = Y_var / ratio
    )
    return(Y_noise)
}

onehot_encoding <- function(X) {
    origin_colnames <- colnames(X)
    for (column_name in origin_colnames) {
        for (value in 0:1) {
            X[paste0(column_name, "_", value)] <-
                as.numeric(X[column_name] == value)
        }
    }
    X <- X[, !(names(X) %in% origin_colnames)]
    return(X)
}

columns_powers <- function(X, pow = 3) {
    origin_colnames <- colnames(X)
    for (column_name in origin_colnames) {
        for (power_ in 1:pow) {
            X[paste0(column_name, "_", power_)] <-
                X[column_name]**power_
        }
    }
    X <- X[, !(names(X) %in% origin_colnames)]
    return(X)
}

basic_continous_matrix <- function(n = 100, p = 16) {
    Z <- rnorm(n * (p + 1))
    dim(Z) <- c(n, (p + 1))
    X <- (Z + Z[, (p + 1)]) / sqrt(2)
    X <- X[, -(p + 1)] |> as.data.frame()
    colnames(X) <- paste0("X", 1:ncol(X))
    return(X)
}

#' Creation of type 1 data set 
#'
#' @param n number of observations
#' @param p number of variables
#'
#' @return list with three elements - X: design matrix, y: target variable,
#' betas: coefficients used to create y
#' @export
#'
create_model1 <- function(n = 50, p = 15) {
    cov_matrix <- matrix(
        nrow = p,
        ncol = p,
        data = 1
    )
    for (i in 1:p) {
        for (j in 1:p) {
            cov_matrix[i, j] <- 0.5**abs(i - j)
        }
    }

    Z <- mvrnorm(n, rep(0, p), cov_matrix)
    Z <- categorize_matrix(Z) |> as.data.frame()
    colnames(Z) <- paste0("Z", 1:ncol(Z))
    Z <- onehot_encoding(Z)

    colnames_ <- colnames(Z)
    gs <- gramSchmidt(as.matrix(Z))
    Z <- gs$Q
    colnames(Z) <- colnames_

    Y <- 1.8 * Z[, "Z1_1"] - 1.2 * Z[, "Z1_0"] +
        0.5 * Z[, "Z3_0"] + Z[, "Z5_0"] + Z[, "Z5_1"]
    Y_noise <- generate_noise(Y, 1.8)
    Y <- Y + Y_noise
    Z <- Z |> as.matrix()


    betas <- c(-1.2, 1.8, 0, 0, 0.5, 1, 0, 0, 1, 1)
    betas <- c(betas, rep(0, 2 * p - length(betas)))

    return(list(
        X = Z, Y = Y,
        groups = rep(1:p, each = 2),
        betas = betas
    ))
}

#' Creation of type 2 data set 
#'
#' @param n number of observations
#' @param p number of variables
#'
#' @return list with three elements - X: design matrix, y: target variable,
#' betas: coefficients used to create y
#' @export
#'
create_model2 <- function(n = 100, p = 4) {
    cov_matrix <- matrix(
        nrow = p,
        ncol = p,
        data = 1
    )
    for (i in 1:p) {
        for (j in 1:p) {
            cov_matrix[i, j] <- 0.5**abs(i - j)
        }
    }

    Z <- mvrnorm(n, rep(0, p), cov_matrix)
    Z <- categorize_matrix(Z) |> as.data.frame()
    colnames(Z) <- paste0("Z", 1:ncol(Z))
    Z <- onehot_encoding(Z)
    
    cols <- colnames(Z)
    for (idx1 in 1:length(cols)){
        for (idx2 in idx1:length(cols)){
            column1 <- cols[idx1]
            column2 <- cols[idx2]
            name1 <- substr(column1, 1, 2)
            name2 <- substr(column2, 1, 2)
            value1 <- substr(column1, 4, 4)
            value2 <- substr(column2, 4, 4)

            if (name1 != name2){
                Z[paste0(name1, name2, "_", value1, "_", value2)] = Z[column1] * Z[column2]
            }
        }    
    }
    Z_new <- Z[,(p*2+1):ncol(Z)]
    Z_new <- Z_new[,order(colnames(Z_new))]
    Z <- cbind(Z[,1:(p*2)], Z_new)

    colnames_ <- colnames(Z)
    gs <- gramSchmidt(as.matrix(Z))
    Z <- gs$Q
    colnames(Z) <- colnames_

    Y <- 3 * Z[, "Z1_1"] + 2 * Z[, "Z1_0"] +
        3 * Z[, "Z2_1"] + 2 * Z[, "Z2_0"] +
        Z[, "Z1Z2_1_1"] + 1.5 * Z[,"Z1Z2_1_0"] +
        2 * Z[, "Z1Z2_0_1"] + 2.5 * Z[, "Z1Z2_0_0"]
        

    Y_noise <- generate_noise(Y, 3)
    Y <- Y + Y_noise
    Z <- Z |> as.matrix()

    betas <- c(2, 3, 2, 3,
               0, 0, 0, 0,
               2.5, 2, 1.5, 1)
    betas <- c(betas, rep(0, ncol(Z) - length(betas)))

    return(list(
        X = Z, Y = Y,
        groups = c(rep(1:p, each = 2), rep(1:(p*(p - 1)/2), each = 4)),
        betas = betas
    ))
}

#' Creation of type 3 data set 
#'
#' @param n number of observations
#' @param p number of variables
#'
#' @return list with three elements - X: design matrix, y: target variable,
#' betas: coefficients used to create y
#' @export
#'
create_model3 <- function(n = 100, p = 16) {
    X <- basic_continous_matrix(n = n, p = p)
    X <- columns_powers(X)

    colnames_ <- colnames(X)
    gs <- gramSchmidt(as.matrix(X))
    X <- gs$Q
    colnames(X) <- colnames_

    Y <- X[, "X3_3"] + X[, "X3_2"] +
        (1 / 3) * X[, "X6_3"] - X[, "X6_2"] +
        (2 / 3) * X[, "X6_1"]

    Y <- Y + generate_noise(Y, 5)
    X <- X |> as.matrix()

    betas <- c(
        0, 0, 0,
        0, 0, 0,
        1, 1, 1,
        0, 0, 0,
        0, 0, 0,
        2 / 3, -1, 1 / 3
    )
    betas <- c(betas, rep(0, 3 * p - length(betas)))

    return(
        list(
            X = X, Y = Y,
            groups = rep(1:p, each = 3),
            betas = betas
        )
    )
}

#' Creation of type 4 data set 
#'
#' @param n number of observations
#' @param p1 number of continuous variables
#' @param p1 number of discrete variables
#'
#' @return list with three elements - X: design matrix, y: target variable,
#' betas: coefficients used to create y
#' @export
#'
create_model4 <- function(n = 100, p1 = 10, p2 = 10) {
    X1 <- basic_continous_matrix(n = n, p = p1)
    X1 <- columns_powers(X1)

    X2 <- basic_continous_matrix(n = n, p = p2)
    X2 <- categorize_matrix(X2) |> as.data.frame()
    colnames(X2) <- paste0("X", (1 + p2):(ncol(X2) + p2))
    X2 <- onehot_encoding(X2)

    X <- cbind(X1, X2)

    colnames_ <- colnames(X)
    gs <- gramSchmidt(as.matrix(X))
    X <- gs$Q
    colnames(X) <- colnames_

    Y <- X[, "X3_3"] + X[, "X3_2"] + X[, "X3_1"] +
        (1 / 3) * X[, "X6_3"] - X[, "X6_2"] +
        (2 / 3) * X[, "X6_1"] +
        2 * X[, "X11_0"] + X[, "X11_1"]

    Y <- Y + generate_noise(Y, 1)
    X <- X |> as.matrix()

    groups <- c(rep(1:p1, each = 3), rep(1:p2, each = 2))

    betas <- c(
        0, 0, 0,
        0, 0, 0,
        1, 1, 1,
        0, 0, 0,
        0, 0, 0,
        2 / 3, -1, 1 / 3,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        2, 1
    )
    betas <- c(betas, rep(0, 3 * p1 + 2 * p2 - length(betas)))

    return(list(X = X, Y = Y, groups = groups, betas = betas))
}
