# S4 classes to wrap up all parameters of models

setClass(
    "group_lasso",
    representation(
        X = "matrix",
        y = "numeric",
        betas = "numeric",
        lambda_max = "numeric",
        lambda_best = "numeric",
        Cp = "numeric"
    )
)

setClass(
    "group_lars",
    representation(
        X = "matrix",
        y = "numeric",
        betas = "numeric",
        betas_path = "list",
        Cp = "numeric"
    )
)
