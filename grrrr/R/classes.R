# S4 classes to wrap up all parameters of models

library(ggplot2)

source("R/helpers.R")

setClassUnion("numericOrNULL", c("numeric", "NULL"))
setClassUnion("listOrNULL", c("numeric", "NULL"))

#' Class storing information about group lasso model 
#'
#' @slot X matrix. Design matrix
#' @slot y numeric. Target variable
#' @slot betas numeric. Final beta coefficients 
#' @slot betas_path list. List of all beta coefficients obtain during calculations
#' @slot true_betas numericOrNULL. Beta coefficients used in target variable calculations
#' @slot lambda_max numeric. Maximum value of lambda
#' @slot lambda_best numeric. Value of lambda used for final model
#' @slot Cp numeric. Value of Cp
#' @slot Cp_path list. List of values of Cp obtained during calculations 
#' @slot model_error numericOrNULL. Value of model_error for final model. Not null only if 
#' true_betas was supplied.
#' @slot me_path listOrNULL. List of values of model_error obtained during calculations. Not null only if 
#' true_betas was supplied.
#'
#' @return instance of group_lasso class
#' @export
#'
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
        model_error = "numericOrNULL",
        me_path = "listOrNULL"
    )
)

#' Class storing information about group lasso model 
#'
#' @slot X matrix. Design matrix
#' @slot y numeric. Target variable
#' @slot betas numeric. Final beta coefficients 
#' @slot betas_path list. List of all beta coefficients obtain during calculations
#' @slot true_betas numericOrNULL. Beta coefficients used in target variable calculations
#' @slot Cp numeric. Value of Cp
#' @slot Cp_path list. List of values of Cp obtained during calculations 
#' @slot model_error numericOrNULL. Value of model_error for final model. Not null only if 
#' true_betas was supplied.
#' @slot me_path listOrNULL. List of values of model_error obtained during calculations. Not null only if 
#' true_betas was supplied.
#'
#' @return instance of group_lars class
#' @export
#'
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


#' Class containing information from multiple tests runs
#'
#' @slot name character. Name of the model 
#' @slot model_error numeric. Mean model error
#' @slot model_error_list numeric. All model errors obtained during testing
#' @slot model_error_std numeric. Standard deviation of model error
#' @slot n_factors numeric. Mean number of factors
#' @slot n_factors_list integer. All numbers of factors obtained during testing
#' @slot n_factors_std numeric. Standard deviation of number of factors
#' @slot cpu_time numeric. Mean CPU time
#' @slot cpu_time_list numeric. All CPU times obtained during testing
#' @slot cpu_time_std numeric. Standard deviation of CPU time
#'
#' @return instance of test_results class
#' @export
#'
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

#' Title
#'
#' @slot model_error numeric. Model error obtained in the test
#' @slot n_factors numeric. Number of factors obtained in the test
#' @slot cpu_time numeric. CPU time obtained in the test
#'
#' @return instance of test_result class
#' @export
setClass(
    "test_result",
    representation(
        model_error = "numeric",
        n_factors = "numeric",
        cpu_time = "numeric"
    )
)

#' Object that stores instances of tests_results
#'
#' @slot tests list. List of tests_results instances
#'
#' @return object of class test_container
#' @export
#'
setClass(
    "test_container",
    representation(
        tests = "list"
    )
)

#' Adding test_results to container
#'
#' @param e1 test_container. Instance of class test_container
#' @param e2 test_results. Instance of class test_results
#'
#' @return instance of class test_container with added new test_results.
#' @export
#'
setMethod("+", signature(e1 = "test_container", e2 = "test_results"), function(e1, e2) {
    
    tests = append(e1@tests, e2)
    
})

setGeneric("create_table", function(container) standardGeneric("create_table")) 

#' Creates table with aggregated results of the tests
#'
#' @param container test_container. Instance of class test_container
#'
#' @return data frame with results of all tests. This table's shape is based on the results in the
#' article.
#' @export
#'
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

#' Creates boxplot for test_container.
#'
#' @param container test_container. Instance of test_container_class
#' @param column character. Column which values will be presented in the boxplot. 
#' One of the ("model_error", "n_factors", "cpu_time")
#'
#' @return ggplot2 object with boxplot
#' @export
#'
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

#' Test getter
#'
#' @param container test_container. Instance of test_container class 
#' @param name character. Name of the test to be returned
#'
#' @return instance of class test_results from the container. If there is no test with such a name
#' method will throw an error.
#' @export
#'
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

#' Performs check if results of the models are statistically different
#'
#' @param container test_container. Instance of test_container class 
#' @param tests_rows character. One group of tests (may be an array). Will be presented in the rows 
#' @param tests_cols character. Second group of tests (may be an array). 
#' Will be presented in the columns
#'
#' @return table with p-values of t-test.
#' @export
#'
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

