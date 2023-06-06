source("R/models.R")
source("R/helpers.R")
source("R/classes.R")

calc_group_lasso <- function(
    X, y, groups,
    result_indicator = "cp",
    true_betas = NULL) {
    
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
    best_betas_r2 <- c()
    Cp_min <- Inf
    me_min <- Inf
    r2_min <- -Inf

    ls <- lm(y ~ X+0)
    betas_ls <- as.numeric(ls$coefficients)
    betas <- list()
    cp_path <- list()
    me_path <- list()
    r2_path <- list()
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
        
        r2_path[[i]] <- r2(y, X %*% betas[[i]])[1]
        
        

        if (cp_path[[i]] < Cp_min) {
            best_betas_cp <- betas[[i]]
            best_lambda_cp <- lambda
            Cp_min <- cp_path[[i]]
        }
        
        if (r2_path[[i]] > r2_min) {
            best_betas_r2 <- betas[[i]]
            best_lambda_r2 <- lambda
            r2_min <- r2_path[[i]]
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
            model_error = me_min,
            r2 = r2_min,
            r2_path = r2_path
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
            model_error = me_min,
            r2 = r2_min,
            r2_path = r2_path
        )
    } else if (result_indicator == "r2"){
        model <- new("group_lasso",
                     X = X,
                     y = y,
                     betas = best_betas_r2,
                     betas_path = betas,
                     true_betas = true_betas,
                     lambda_max = max_lambda,
                     lambda_best = best_lambda_r2,
                     Cp = Cp_min,
                     Cp_path = cp_path,
                     model_error = me_min,
                     r2 = r2_min,
                     r2_path = r2_path
        )
    }



    return(model)
}


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
