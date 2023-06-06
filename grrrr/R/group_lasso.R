# source("R/models.R")
source("R/helpers.R")
# source("R/classes.R")

#' Creates a instance of group lasso model 
#'
#' @param X matrix with regressors
#' @param y target variable
#' @param groups list of integers with a length equals to number of columns in X.
#' Indicates to which group given variable belongs to
#' @param result_indicator one of values ("cp", "me"). Indicates which of those two statistic 
#' should be used to select the final model. To use "me" also true_betas needs to be supplied.
#' @param true_betas array of true values of betas
#'
#' @return object of class group_lasso
#' @export
#'
calc_group_lasso <- function(
    X, y, groups,
    result_indicator = "cp",
    true_betas = NULL) {
    
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
    
    n_var <- ncol(X)
    group_sizes <- table(groups) |> as.numeric()
    n_groups <- length(group_sizes)
    indexes <- list()
    for (value in unique(groups)) {
        indexes[[length(indexes) + 1]] <- which(groups == value)
    }


    max_lambda <- purrr::map_dbl(indexes, \(j) {
        (norm_L(
            t(X[, j]) %*% y,
            length(j)
        )/sqrt(length(j))) 
    }) |> max() 


    best_betas_cp <- c()
    best_betas_me <- c()
    Cp_min <- Inf
    me_min <- Inf

    ls <- lm(y ~ X+0)
    betas_ls <- as.numeric(ls$coefficients)
    betas <- list()
    cp_path <- list()
    me_path <- list()
    i <- 0

    for (lambda in seq(
        from = 0,
        to = max_lambda,
        length.out = 100
    )) {
        i <- i + 1
        betas[[i]] <- rep(0, n_var)
        betas_prev <- rep(-1, n_var)

        while (norm_L(betas[[i]] - betas_prev, n_var) > 0.00001) {
            betas_prev <- betas[[i]]
            for (q in 1:n_groups) {
                j <- which(groups == q)
                p <- length(j)

                S <- t(X[, j]) %*% (y - X %*% (betas[[i]] %-% j))
                betas[[i]][j] <-
                    max((1 - (lambda * sqrt(p)) / norm_L(S, p)), 0) * S
            }
        }

        cp_path[[i]] <- calculate_cp(
            indexes, group_sizes,
            betas[[i]], betas_ls, X, y,
            df_lasso
        )
        

        if (cp_path[[i]] < Cp_min) {
            best_betas_cp <- betas[[i]]
            best_lambda_cp <- lambda
            Cp_min <- cp_path[[i]]
        }

        if (!is.null(true_betas)) {
            me_path[[i]] <- calculate_me(
                X, betas[[i]], true_betas
            )

            if (me_path[[i]] < me_min) {
                best_betas_me <- betas[[i]]
                best_lambda_me <- lambda
                me_min <- me_path[[i]]
            }
        }
    }

    if (result_indicator == "cp") {
        model <- new("group_lasso",
            X = X,
            y = y,
            betas = best_betas_cp,
            betas_path = betas,
            true_betas = true_betas,
            lambda_max = max_lambda,
            lambda_best = best_lambda_cp,
            Cp = Cp_min,
            Cp_path = cp_path,
            model_error = me_min
        )
    } else if (result_indicator == "me"){
        model <- new("group_lasso",
            X = X,
            y = y,
            betas = best_betas_me,
            betas_path = betas,
            true_betas = true_betas,
            lambda_max = max_lambda,
            lambda_best = best_lambda_me,
            Cp = Cp_min,
            Cp_path = cp_path,
            model_error = me_min
        )
    }



    return(model)
}


#' Calculates degrees of freedom for group lasso model.
#'
#' @param indexes array with factors chosen in the model
#' @param group_sizes array with sizes of consecutive groups
#' @param betas beta coefficients in the investigated model
#' @param betas_ls beta coefficients in the OLS model build on the same data 
#'
#' @return number of degrees of freedom
#'
df_lasso <- function(indexes, group_sizes, betas, betas_ls) {
    dg_f <- map2_dbl(
        indexes,
        group_sizes,
        \(j, p)
        as.integer(norm_L(betas[j], p) > 0) +
            sqrt(sum(betas[j]**2)) * (p - 1) / sqrt(sum(betas_ls[j]**2))
    ) |> sum()
    return(dg_f)
}
