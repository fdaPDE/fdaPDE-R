## TODOs:
## TODO fdaPDE_Model: add support for function evaluation
##                    !! each statistical model must expose a getter for the Psi matrix !!
## TODO fdaPDE_Model: implement domain_data consistency check in sanity_check_domain
##                    - ncol(nodes) == ncol(locations)
##                    - locations \in domain
##                    - ??
## TODO fdaPDE_Model: add support for timing the execution times
## TODO fdaPDE_Model: add support for areal sampling
## TODO fdaPDE_Model: add support for space-time models and regularization_type detection


# fdaPDE_Model class ----

## interface for a generic fdaPDE model
fdaPDE_Base_Model <- R6::R6Class(
  classname = "fdaPDE_Base_Model",
  public = list(
    ## data
    domain = NULL,
    locations = NULL,
    ## results
    results = list(),
    ## utils
    model_traits = list(),
    ## options
    VERBOSE = FALSE
  ),
  private = list(
    ## model instance
    cpp_model = NULL,
    get_cpp_model = function() {
      return(cpp_model)
    },
    ## setters
    set_verbosity = function(VERBOSE = FALSE) {
      self$VERBOSE <- VERBOSE
    },
    ## utils
    display = function(message) {
      if (self$VERBOSE) {
        cat(paste(message, "\n", sep = ""))
      }
    },
    sanity_check_domain = function(data) {
      ## TODO: implement dimensions checks
      self$domain <- data$domain
      ## sampling type deduction
      if (is.null(data$locations)) {
        self$model_traits$sampling_type <- Sampling(0)
        self$locations <- as.matrix(data$domain$nodes)
      } else {
        self$model_traits$sampling_type <- Sampling(1)
        self$locations <- as.matrix(data$locations)
      } ## TODO: support for areal sampling
      ## regularization type deduction (forced to SpaceOnly)
      self$model_traits$regularization_type <- Regularization(0)
    }
  )
)

# helping functions

export_calibrator_results <- function(cpp_calibrator) {
  ## list initialization
  calibrator_results <- list()
  ## optimal lambda
  calibrator_results$lambda_opt <- cpp_calibrator$optimum()
  ## calibrator specific results
  if (cpp_calibrator$get_calibration_strategy() == Calibration("gcv")) {
    calibrator_results$edfs <- as.matrix(cpp_calibrator$edfs())
    calibrator_results$gcvs <- as.matrix(cpp_calibrator$gcvs())
  }
  if (cpp_calibrator$get_calibration_strategy() == Calibration("kcv")) {
    calibrator_results$avg_scores <- as.matrix(cpp_calibrator$avg_scores())
    calibrator_results$std_scores <- as.matrix(cpp_calibrator$std_scores())
    calibrator_results$scores <- as.matrix(cpp_calibrator$scores())
  }
  return(calibrator_results)
}
