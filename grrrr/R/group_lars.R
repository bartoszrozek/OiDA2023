# source("R/models.R")
source("R/helpers.R")
# source("R/classes.R")

#' Finds optimum for quadratic equation needed to find next factor included in the LARS algorithm. 
#'
#' @param X matrix with regressors
#' @param r current residuals
#' @param j candidate factor
#' @param mcs current "active set"
#' @param gamma_ current direction
#'
#' @return value of root which is in [0,1] interval
#'
find_alpha_lars <- function(X, r, j, mcs, gamma_) {
    sq_coef <- t(gamma_) %*% t(X) %*% X[, j] %*% t(X[, j]) %*% X %*% gamma_ -
        t(gamma_) %*% t(X) %*% X[, mcs] %*% t(X[, mcs]) %*% X %*% gamma_
    lin_coef <- t(r) %*% X[, mcs] %*% t(X[, mcs]) %*% X %*% gamma_ +
        t(gamma_) %*% t(X) %*% X[, mcs] %*% t(X[, mcs]) %*% r -
        t(r) %*% X[, j] %*% t(X[, j]) %*% X %*% gamma_ +
        t(gamma_) %*% t(X) %*% X[, j] %*% t(X[, j]) %*% r
    ww <- t(r) %*% X[, j] %*% t(X[, j]) %*% r -
        t(r) %*% X[, mcs] %*% t(X[, mcs]) %*% r
    roots <- quad_roots(sq_coef, lin_coef, ww)
    root <- roots[roots >= 0 & roots <= 1]
    if (length(root) != 1) {
        return(1)
    }
    return(root)
}

#' Calculates degrees of freedom for group lars model.
#'
#' @param indexes array with factors chosen in the model
#' @param group_sizes array with sizes of consecutive groups
#' @param betas beta coefficients in the investigated model
#' @param betas_ls beta coefficients in the OLS model build on the same data 
#'
#' @return number of degrees of freedom
#'
df_lars <- function(indexes, group_sizes, betas, betas_ls) {
    dg_f <- map2_dbl(
        indexes,
        group_sizes,
        \(j, p)
        as.integer(norm_L(betas[j], p) > 0) +
            sqrt(sum(betas[j]**2)) * (p - 1) / sqrt(sum(betas_ls[j]**
                2))
    ) |> sum()
    return(dg_f)
}


#' Creates a instance of group lars model 
#'
#' @param X matrix with regressors
#' @param y target variable
#' @param groups list of integers with a length equals to number of columns in X.
#' Indicates to which group given variable belongs to
#' @param result_indicator one of values ("cp", "me"). Indicates which of those two statistic 
#' should be used to select the final model. To use "me" also true_betas needs to be supplied.
#' @param true_betas array of true values of betas
#'
#' @return object of class group_lars
#' @export
#'
calc_group_lars <- function(X, y, groups,
                            result_indicator = "cp", true_betas = NULL) {
    
    if (!is.matrix(X)) {
        stop("X should be a matrix.")
    }
    
    if (!is.numeric(y) && is.null(dim(y))) {
        stop("y should be a numeric array.")
    }
    
    if (nrow(X) != length(y)) {
        stop("Dimensions of X and y does not match.")
    }
    
    if (!result_indicator %in% c("cp", "me")) {
        stop("result_indicator should be one of c('cp', 'me').")
    }
    if (is.null(true_betas) && result_indicator == "me") {
        stop("To use me as result_indicator use should supply true_betas.")
    }
    if (length(groups) != ncol(X)) {
        stop("groups should have the same length as number of X columns")
    }
    
    n <- nrow(X)
    n_var <- ncol(X)
    group_sizes <- table(groups) |> as.numeric()
    n_groups <- length(group_sizes)
    indexes <- list()
    ls <- lm(y ~ X+0)
    betas_ls <- as.numeric(ls$coefficients)

    for (value in unique(groups)) {
        indexes[[length(indexes) + 1]] <- which(groups == value)
    }

    r <- list()
    betas <- list()

    best_betas <- c()
    best_betas_me <- c()
    Cp_min <- Inf
    me_min <- Inf
    cp_path <- list()
    me_path <- list()

    betas[[1]] <- rep(0, n_var)
    r[[1]] <- y
    k <- 1

    direction <- purrr::map2_dbl(
        indexes,
        group_sizes,
        \(j, p) norm_L(t(X[, j]) %*% r[[k]], p) / p
    ) |> which.max()

    mcs <- indexes[direction] |> unlist()
    indexes_rest <- indexes[-direction]


    while (length(indexes_rest) > 0) {
        gamma_mcs <- solve(t(X[, mcs]) %*% X[, mcs]) %*%
            t(X[, mcs]) %*% r[[k]]
        gamma_ <- betas[[k]] * 0
        gamma_[mcs] <- gamma_mcs

        alphas <- map(
            indexes_rest,
            \(jp) find_alpha_lars(X, r[[k]], unlist(jp), mcs, gamma_)
        )

        min_idx <- which.min(alphas)
        alpha_ <- min(unlist(alphas))
        mcs <- c(mcs, indexes_rest[min_idx]) |> unlist()
        indexes_rest <- indexes_rest[-min_idx]

        k <- k + 1
        betas[[k]] <- betas[[k - 1]] + alpha_ * gamma_
        r[[k]] <- y - X %*% betas[[k]]

        cp_path[[k]] <- calculate_cp(
            indexes, group_sizes,
            betas[[k]], betas_ls, X, y,
            df_lars
        )

        if (cp_path[[k]] < Cp_min) {
            best_betas_cp <- betas[[k]]
            Cp_min <- cp_path[[k]]
        }
        
        if (!is.null(true_betas)) {
            me_path[[k]] <- calculate_me(
                X, betas[[k]], true_betas
            )
            
            if (me_path[[k]] < me_min) {
                best_betas_me <- betas[[k]]
                me_min <- me_path[[k]]
            }
        }
        
    }

    if (result_indicator == "cp") {
        model <- new("group_lars",
                     X = X,
                     y = y,
                     betas = best_betas_cp,
                     betas_path = betas,
                     true_betas = true_betas,
                     Cp = Cp_min,
                     Cp_path = cp_path,
                     model_error = me_min
        )
    } else if (result_indicator == "me"){
        model <- new("group_lars",
                     X = X,
                     y = y,
                     betas = best_betas_me,
                     betas_path = betas,
                     true_betas = true_betas,
                     Cp = Cp_min,
                     Cp_path = cp_path,
                     model_error = me_min
        )
    } else{
        stop("Wrong result_indicator value!")
    }
    return(model)
}


#' Very simple quadratic equation solver
#'
#' @param a quadratic coefficient
#' @param b linear coefficient
#' @param c constant coefficient
#'
#' @return array with two roots
quad_roots <- function(a, b, c) {
    return(c(((
        -b - sqrt(b ^ 2 - 4 * a * c)
    ) / (2 * a)), ((
        -b + sqrt(b ^ 2 -
                      4 * a * c)
    ) / (2 * a))))
}