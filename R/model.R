## TODOs:
## TODO fdaPDE_Model: add support for function evaluation
## TODO fdaPDE_Model: add support for timing of execution times
## TODO fdaPDE_Model: add support for areal sampling
## TODO fdaPDE_Regression_Model: implement regression data sanity checks:
##      - coherence between locations and observations dimensions
## TODO fdaPDE_Model: add support for space-time models
## TODO: fix error raised when multiple SRPDE model are initialized one after the other

# fdaPDE_Model class ----

## interface for a generic fdaPDE model
fdaPDE_Model <- R6::R6Class(
  classname = "fdaPDE_Model",
  public = list(
    ## data
    domain = NULL,
    locations = NULL,
    ## utils
    model_traits = list(),
    VERBOSE = FALSE,
    ## constructor
    initialize = function(VERBOSE = FALSE) {
      self$VERBOSE <- VERBOSE
    }
  ),
  private = list(
    ## model instance
    cpp_model = NULL,
    ## utils
    display = function(message) {
      if (self$VERBOSE) {
        cat(paste(message, "\n", sep = ""))
      }
    },
    sanity_check_domain = function(data) {
      self$domain <- data$domain
      if (is.null(data$locations)) {
        self$model_traits$sampling_type <- Sampling(0)
        self$locations <- data$domain$nodes
      } else {
        self$model_traits$sampling_type <- Sampling(1)
        self$locations <- data$locations
      }
      ## TODO: support for areal sampling
    }
  )
)

## fdaPDE_Regression_Model models ----

## interface for a generic regression model
fdaPDE_Regression_Model <- R6::R6Class(
  classname = "fdaPDE_Regression_Model",
  inherit = fdaPDE_Model,
  public = list(
    ## regression data
    observations = NULL,
    covariates = NULL,
    ## results
    results = list()
  ),
  private = list(
    ## calibrator instance
    cpp_calibrator = NULL,
    ## utilities
    init_calibrator = function(calibrator) {
      private$cpp_calibrator <- calibrators_factory(calibrator)
    },
    sanity_check_data = function(formula = NULL, data, penalty) {
      ## domain infotmations
      super$sanity_check_domain(data)
      ## regression data
      private$sanity_check_regression_data(formula, data)
      ## smoother parameters
      smoother_params <- smoother("SRPDE",
        penalty = penalty,
        sampling_type = self$model_traits$sampling_type
      )
      return(smoother_params)
    },
    init_model = function(smoother_params) {
      super$cpp_model <- regression_models_factory(
        self$domain,
        smoother_params
      )
      ## set locations
      super$cpp_model$set_spatial_locations(as.matrix(self$locations))
      super$display("  Locations have been set.")
      ## set data
      super$cpp_model$set_observations(as.matrix(self$observations))
      super$display("  Observations have been set.")
      if (!is.null(self$covariates)) {
        super$cpp_model$set_covariates(as.matrix(self$covariates))
        super$display("  Covariates have been set.")
      }
    },
    sanity_check_regression_data = function(formula, data) {
      self$observations <- data$observations
      self$covariates <- data$covariates
      ## TODO: check coherence between observations and locations dimensions
    },
    calibrate = function(lambda) {
      ## mondel initialization (necessary for GCV and KCV calibration strategies)
      super$cpp_model$init()
      ## calibration
      private$cpp_calibrator$set_lambda(lambda)
      if (private$cpp_calibrator$get_calibration_strategy() != 1) {
        lambda_opt <- private$cpp_calibrator$fit(super$cpp_model$get_view())
      } else {
        lambda_opt <- private$cpp_calibrator$fit(super$cpp_model$get_gcv())
      }
      ## statistical model preparation for fit with optimal lambda
      private$set_lambda(lambda_opt)
      ## save calibrator's results
      self$results$calibrator$lambda_opt <- lambda_opt
      if (private$cpp_calibrator$get_calibration_strategy() == 1) {
        self$results$calibrator$edfs <- as.matrix(private$cpp_calibrator$edfs())
        self$results$calibrator$gcvs <- as.matrix(private$cpp_calibrator$gcvs())
      }
      if (private$cpp_calibrator$get_calibration_strategy() == 2) {
        self$results$calibrator$avg_scores <- as.matrix(private$cpp_calibrator$avg_scores())
        self$results$calibrator$std_scores <- as.matrix(private$cpp_calibrator$std_scores())
        self$results$calibrator$scores <- as.matrix(private$cpp_calibrator$scores())
      }
    },
    set_lambda = function(lambda) {
      super$cpp_model$set_lambda(lambda[1])
    },
    fit = function() {
      ## fitting the statistical model using the optimal lambda set by the calibrator
      super$cpp_model$init()
      super$cpp_model$solve()
      ## save fit results
      self$results$f <- as.matrix(super$cpp_model$f())
      self$results$beta <- as.matrix(super$cpp_model$beta())
    }
  )
)


### SRPDE model implementation ----

#' @export
SRPDE <- R6::R6Class(
  classname = "SRPDE_model",
  inherit = fdaPDE_Regression_Model,
  public = list(
    ## constructor
    initialize = function(formula = NULL, data,
                          penalty = simple_laplacian_penalty(),
                          solver = NULL,
                          calibrator = off(),
                          VERBOSE = FALSE) {
      super$initialize(VERBOSE)
      super$display("\n\nSRPDE model\n")
      ## calibrator initialization
      super$display("- Calibrator initialization")
      super$init_calibrator(calibrator)
      ## inputs sanity check
      super$display("- Inputs sanity check")
      smoother_params <- super$sanity_check_data(formula, data, penalty)
      ## model initialization
      super$display("- Model initialization")
      super$init_model(smoother_params)
    },
    ## utilities
    fit = function(lambda = hyperparameters(1e-4),
                   calibrator = NULL) {
      ## clear previous results
      self$results <- list()
      ## calibrator update
      if (!is.null(calibrator)) {
        super$init_calibrator(calibrator)
      }
      ## calibration
      super$display("- Calibration")
      super$calibrate(lambda)
      ## fit
      super$display("- Fit")
      super$fit()
    }
  ),
  private = list()
)
