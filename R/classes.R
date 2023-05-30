# S4 classes to wrap up all parameters of models

setClassUnion("numericOrNULL", c("numeric", "NULL"))

setClass(
    "group_lasso",
    representation(
        X = "matrix",
        y = "numeric",
        betas = "numeric",
        true_betas = "numericOrNULL",
        lambda_max = "numeric",
        lambda_best = "numeric",
        Cp = "numeric",
        model_error = "numericOrNULL"
    )
)

setClass(
    "group_lars",
    representation(
        X = "matrix",
        y = "numeric",
        betas = "numeric",
        true_betas = "numericOrNULL",
        betas_path = "list",
        Cp = "numeric",
        model_error = "numericOrNULL"
    )
)


setClass(
    "test_results",
    representation(
        model_error = "numeric",
        model_error_std = "numeric",
        n_factors = "numeric",
        n_factors_std = "numeric",
        cpu_time = "numeric",
        cpu_time_std = "numeric"
    )
)

setClass(
    "test_result",
    representation(
        model_error = "numeric",
        n_factors = "numeric",
        cpu_time = "numeric"
    )
)
class(test_result)
setMethod("mean", signature("test_result"), function(x) 1)
mean(c(model, model))
