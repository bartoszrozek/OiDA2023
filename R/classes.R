# S4 classes to wrap up all parameters of models

library(ggplot2)

source("helpers.R")

setClassUnion("numericOrNULL", c("numeric", "NULL"))
setClassUnion("listOrNULL", c("numeric", "NULL"))

setClass(
    "group_lasso",
    representation(
        X = "matrix",
        y = "numeric",
        betas = "numeric",
        betas_path = "list",
        true_betas = "numericOrNULL",
        lambda_max = "numeric",
        lambda_best = "numeric",
        Cp = "numeric",
        Cp_path = "list",
        r2 = "numeric",
        r2_path = "list",
        model_error = "numericOrNULL",
        me_path = "listOrNULL"
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
        Cp_path = "list",
        model_error = "numericOrNULL",
        me_path = "listOrNULL"
    )
)


setClass(
    "test_results",
    representation(
        name = "character",
        model_error = "numeric",
        model_error_list = "numeric",
        model_error_std = "numeric",
        n_factors = "numeric",
        n_factors_list = "integer",
        n_factors_std = "numeric",
        cpu_time = "numeric",
        cpu_time_list = "numeric",
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

setClass(
    "test_container",
    representation(
        tests = "list"
    )
)

container <- new(
    "test_container"
)

setMethod("+", signature(e1 = "test_container", e2 = "test_results"), function(e1, e2) {
    new(
        "test_container",
        tests = append(e1@tests, e2)
    )
})

setGeneric("create_table", function(container) standardGeneric("create_table")) 

setMethod("create_table", signature("test_container"), function(container) {
    output_table <- data.frame()
    
    for (test in container@tests) {
        row <- data.frame(
            row.names = test@name,
            "Model error" = test@model_error,
            "Model error (std)" = test@model_error_std,
            "N factors" = test@n_factors,
            "N factors (std)" = test@n_factors_std,
            "CPU time" = test@cpu_time,
            "CPU time (std)" = test@cpu_time_std,
            check.names = FALSE
            
        )
        output_table <- rbind(output_table, row)
    }
    output_table <- output_table[order(row.names(output_table)), ]
    return(output_table)
})


setGeneric("create_boxplot", function(container, column) standardGeneric("create_boxplot")) 
setMethod("create_boxplot",
          signature(container = "test_container",
                    column = "character"),
          function(container, column) {
              results_table <- data.frame()
              
              column_name <- gsub("_", " ", column)
              column_name <- first_up(column_name)
              column_geter <- paste0(column, "_list")
              
              
              for (test in container@tests) {
                  row <- data.frame(
                      Method = test@name,
                      value = slot(test, column_geter),
                      check.names = FALSE
                      
                  )
                  results_table <- rbind(results_table, row)
              }
              results_table <- results_table[order(results_table[,1]), ]
             
              
              p <- ggplot(results_table, aes(x = Method, y = value)) +
                  geom_boxplot() +
                  xlab("") +
                  ylab(column_name) +
                  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
              
              return(p)
          })


setGeneric("get_test", function(container, name) standardGeneric("get_test")) 
setMethod("get_test",
          signature(container = "test_container",
                    name = "character"),
          function(container, name) {
              
              for (test in container@tests){
                  if (test@name == name) return(test)
              }
              stop("Test with a given name does not exist!")
          })

setGeneric("perform_ttest", function(container, tests_rows, tests_cols) standardGeneric("perform_ttest")) 
setMethod("perform_ttest",
          signature(container = "test_container",
                    tests_rows = "character",
                    tests_cols = "character"),
          function(container, tests_rows, tests_cols) {
              
              output_table <- data.frame()
              
              for (model_row in tests_rows){
                  test_row <- get_test(container, model_row)
                  results <- c()
                  for (model_col in tests_cols){
                      test_col <- get_test(container, model_col)
                      results <- c(results, t.test(test_row@model_error_list, 
                                                   test_col@model_error_list)$p.value)
                  }
                  output_table <- rbind(output_table, t(results))
              }
              colnames(output_table) <- tests_cols
              rownames(output_table) <- tests_rows
              return(output_table)
          })

