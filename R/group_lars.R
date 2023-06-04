source("models.R")
source("helpers.R")
source("classes.R")

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



calc_group_lars <- function(X, y, groups,
                            result_indicator = "cp", true_betas = NULL) {
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


quad_roots <- function(a, b, c) {
    return(c(((
        -b - sqrt(b ^ 2 - 4 * a * c)
    ) / (2 * a)), ((
        -b + sqrt(b ^ 2 -
                      4 * a * c)
    ) / (2 * a))))
}