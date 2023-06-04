source("helpers.R")
source("models.R")
source("group_lars.R")
source("group_lasso.R")
library(gglasso)

box::use(
    MASS[stepAIC]
)

calculate_test <- function(name, test_function, n, create_model, ...) {
    model_error <- c()
    n_factors <- c()
    cpu_time <- c()

    for (i in 1:n) {
        tryCatch(
            {
                model_data <- create_model()
                X <- model_data$X
                y <- model_data$Y
                true_betas <- model_data$betas
                groups <- model_data$groups

                result <- test_function(X, y, true_betas, groups, ...)
                model_error <- c(model_error, result@model_error)
                n_factors <- c(n_factors, result@n_factors)
                cpu_time <- c(cpu_time, result@cpu_time)
            },
            error = function(cond) {}
        )
    }

    agg_results <- new(
        "test_results",
        name = name,
        model_error = mean(model_error),
        model_error_list = model_error,
        model_error_std = sd(model_error),
        n_factors = mean(n_factors),
        n_factors_list = n_factors,
        n_factors_std = sd(n_factors),
        cpu_time = mean(cpu_time),
        cpu_time_list = cpu_time,
        cpu_time_std = sd(cpu_time)
    )

    return(agg_results)
}

test_ls <- function(X, y, true_betas, groups) {
    start_time <- Sys.time()
    ls_model <- lm(y ~ X + 0)
    betas <- ls_model$coefficients
    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    me <- calculate_me(X, betas, true_betas)

    model <- new("test_result",
        model_error = me,
        n_factors = length(betas),
        cpu_time = time
    )
    return(model)
}

test_step <- function(X, y, true_betas, groups) {
    start_time <- Sys.time()
    ls_model <- lm(y ~ X + 0)
    # Stepwise regression model
    step_model <- stepAIC(ls_model,
        direction = "both",
        trace = FALSE
    )

    betas <- step_model$coefficients

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    me <- calculate_me(X, betas, true_betas)

    model <- new("test_result",
        model_error = me,
        n_factors = (summary(step_model)$coefficients[, 4] < 0.05) |> sum(),
        cpu_time = time
    )
    return(model)
}

test_lars_group <- function(X, y, true_betas, groups, ...) {
    start_time <- Sys.time()
    lars_model <- calc_group_lars(X, y, groups,true_betas = true_betas, ...)

    betas <- lars_model@betas

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    me <- calculate_me(X, betas, true_betas)

    model <- new("test_result",
        model_error = me,
        n_factors = count_factors(betas, colnames(X)),
        cpu_time = time
    )
    return(model)
}

test_lars <- function(X, y, true_betas, groups, ...) {
    start_time <- Sys.time()
    lars_model <- calc_group_lars(X, y, 1:length(groups), true_betas = true_betas, ...)

    betas <- lars_model@betas

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    me <- calculate_me(X, betas, true_betas)

    model <- new("test_result",
        model_error = me,
        n_factors = count_factors(betas, colnames(X)),
        cpu_time = time
    )
    return(model)
}

test_lasso_group_library <- function(X, y, true_betas, groups) {
    start_time <- Sys.time()
    gr_cv <- cv.gglasso(
        x = X, y = y, group = groups,
        loss = "ls", pred.loss = "L2",
        intercept = F, nfolds = 5
    )

    gr <- gglasso(X, y,
        lambda = gr_cv$lambda.1se,
        group = groups, loss = "ls",
        intercept = F
    )

    betas <- gr$beta


    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    me <- calculate_me(X, betas, true_betas)

    model <- new("test_result",
        model_error = me,
        n_factors = count_factors(betas, colnames(X)),
        cpu_time = time
    )
    return(model)
}

test_lasso_group <- function(X, y, true_betas, groups, ...) {
    start_time <- Sys.time()
    lars_model <- calc_group_lasso(X, y, groups, true_betas = true_betas, ...)

    betas <- lars_model@betas

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    me <- calculate_me(X, betas, true_betas)

    model <- new("test_result",
        model_error = me,
        n_factors = count_factors(betas, colnames(X)),
        cpu_time = time
    )
    return(model)
}

count_factors <- function(betas, betas_names){
    non_zero_betas <- betas_names[betas != 0]
    factors_names <- map(non_zero_betas, 
                         \(name) substr(name, 1, nchar(name)-2)) |>
        unique()
    n_factors <- length(factors_names)
    return(n_factors)
}
