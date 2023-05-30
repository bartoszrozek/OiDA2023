source("R/helpers.R")

box::use(
    MASS[stepAIC]
)

tester <- function(test_function, n, ...) {
    model_error <- c()
    n_factors <- c()
    cpu_time <- c()

    for (i in 1:n) {
        result <- test_function(...)
        model_error <- c(model_error, result@model_error)
        n_factors <- c(n_factors, result@n_factors)
        cpu_time <- c(cpu_time, result@cpu_time)
    }

    agg_results <- new("test_results",
        model_error = mean(model_error),
        model_error_std = sd(model_error),
        n_factors = mean(n_factors),
        n_factors_std = sd(n_factors),
        cpu_time = mean(cpu_time),
        cpu_time_std = sd(cpu_time)
    )

    return(agg_results)
}

test_ls <- function(X, y, true_betas) {
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

test_step <- function(X, y, true_betas) {
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
        n_factors = length(betas),
        cpu_time = time
    )
    return(model)
}
