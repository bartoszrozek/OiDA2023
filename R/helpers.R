box::use(
    purrr[
        map, pmap, pmap_dfr,
        map2, map_dbl, map2_dbl
    ]
)


# this is implementation of
# b_{-j}=(b^'_1, ..., b^'_{j-1}, 0', b^'_{j+1}, ..., b^'_{J})
`%-%` <- function(vector, index) {
    vector[index] <- 0
    return(vector)
}

`%_%` <- function(object, index) {
    if (dim(object) |> is.null()) {
        out <- object * 0
        out[index] <- object[index]
    } else {
        out <- object * 0
        out[, index] <- object[, index]
    }

    return(out)
}

norm_L <- function(vector, p) {
    K <- diag(p)
    norm_ <- sqrt(t(vector) %*% K %*% vector)[1]
    return(norm_)
}

calculate_cp <- function(indexes, group_sizes,
                         betas, betas_ls, X, y, df_function) {
    dg_f <- df_function(indexes, group_sizes, betas, betas_ls)

    mu <- X %*% betas

    Cp <- (sum((y - mu)**2) / var(y - mu)[1]) - nrow(X) + 2 * dg_f
    return(Cp)
}
