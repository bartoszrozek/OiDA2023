source("R/models.R")
source("R/helpers.R")
source("R/classes.R")

calc_group_lasso <- function(
    X, y, groups,
    cp_indicator = TRUE, true_betas = NULL) {
    n_var <- ncol(X)
    group_sizes <- table(groups) |> as.numeric()
    n_groups <- length(group_sizes)
    indexes <- list()
    for (value in unique(groups)) {
        indexes[[length(indexes) + 1]] <- which(groups == value)
    }


    max_lambda <- purrr::map_dbl(indexes, \(j)
    (norm_L(
        t(X[, j]) %*% y,
        length(j)
    ) / sqrt(p))[1]) |> max()


    best_betas_cp <- c()
    best_betas_me <- c()
    Cp_min <- Inf
    me_min <- Inf

    ls <- lm(y ~ X)
    betas_ls <- as.numeric(ls$coefficients)
    betas <- list()
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

        Cp <- calculate_cp(
            indexes, group_sizes,
            betas[[i]], betas_ls, X, y,
            df_lasso
        )

        if (Cp < Cp_min) {
            best_betas_cp <- betas[[i]]
            best_lambda_cp <- lambda
            Cp_min <- Cp
        }

        if (!is.null(true_betas)) {
            me <- calculate_me(
                X, betas[[i]], true_betas
            )

            if (me < me_min) {
                best_betas_me <- betas[[i]]
                best_lambda_me <- lambda
                me_min <- me
            }
        }
    }

    if (cp_indicator) {
        model <- new("group_lasso",
            X = X,
            y = y,
            betas = best_betas_cp,
            true_betas = true_betas,
            lambda_max = max_lambda,
            lambda_best = best_lambda_cp,
            Cp = Cp_min,
            model_error = me_min
        )
    } else {
        model <- new("group_lasso",
            X = X,
            y = y,
            betas = best_betas_cp,
            true_betas = true_betas,
            lambda_max = max_lambda,
            lambda_best = best_lambda_cp,
            Cp = Cp_min,
            model_error = me_min
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
            sqrt(sum(betas[j]**2)) * (p - 1) / sqrt(sum(betas_ls[j]**
                2))
    ) |> sum()
    return(dg_f)
}
