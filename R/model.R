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
      smoother <- smoothing("SRPDE",
        penalty = penalty,
        sampling_type = self$model_traits$sampling_type
      )
      return(smoother)
    },
    init_model = function(smoother) {
      super$cpp_model <- regression_models_factory(
        self$domain,
        smoother
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
      private$cpp_calibrator$configure_calibrator(lambda)
      if (private$cpp_calibrator$get_calibration_strategy() != 1) {
        lambda_opt <- private$cpp_calibrator$fit(super$cpp_model$get_view())
      } else {
        lambda_opt <- private$cpp_calibrator$fit(super$cpp_model$get_gcv())
      }
      ## statistical model preparation for fit with optimal lambda
      private$set_lambda(lambda_opt)
      ## save calibrator's results
      self$results$calibrator$lambda_opt <- private$cpp_calibrator$optimum()
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

## functional centering models ----
fCentering <- R6::R6Class(
  classname = "fCentering",
  inherit = fdaPDE_Model,
  public = list(
    ## data
    X = NULL, ## data to center
    w = NULL, ## weights
    ## results
    results = list(),
    ## constructor
    initialize = function(data,
                          smoother = smoothing(),
                          calibrator = off(),
                          VERBOSE = FALSE) {
      super$initialize(VERBOSE)
      super$display("\n\nCentering model\n")
      ## inputs sanity check
      super$display("- Inputs sanity check")
      private$sanity_check_data(data, penalty)
      ## Centering initialization
      super$display("- Centering initialization")
      private$init_centering()
      ## calibrator initialization
      super$display("- Calibrator initialization")
      private$init_calibrator(calibrator)
      ## smoother initialization
      super$display("- Smoother initialization")
      private$init_smoother(smoother)
    },
    ## utilities
    fit = function(lambda = hyperparameters(1e-4),
                   smoother = NULL, calibrator = NULL) {
      ## clear previous results
      self$results <- list()
      ## smoother update
      if (!is.null(smoother)) {
        super$display("- Smoother re-initialization")
        private$init_smoother(smoother)
      }
      ## calibrator update
      if (!is.null(calibrator)) {
        super$display("- Calibrator re-initialization")
        private$init_calibrator(calibrator)
      }
      ## calibrato configuration
      super$display("- Calibrator configuration")
      super$cpp_model$set_calibrator(
        private$cpp_calibrator$configure_calibrator(lambda)
      )
      # centering laveraging the selected calibrator
      super$display("- Centering")
      super$cpp_model$set_lambda(lambda)
      super$cpp_model$init()
      super$cpp_model$solve()
      ## save fit results
      self$results$mean <- as.matrix(super$cpp_model$mean())
      self$results$centered <- as.matrix(super$cpp_model$centered())
    }
  ),
  private = list(
    ## smoother instance
    cpp_smoother = NULL,
    ## calibrator instance
    cpp_calibrator = NULL,
    ## initialize utilities
    sanity_check_data = function(data, penalty) {
      ## domain informations
      super$sanity_check_domain(data)
      ## regression data
      private$sanity_check_centering_data(data)
    },
    sanity_check_centering_data = function(data) {
      self$X <- data$X
      self$w <- data$w
      ## TODO: check coherence between X and locations dimensions
      ## TODO: check coherence between X and w dimensions
    },
    init_centering = function() {
      super$cpp_model <- new(cpp_center)
      ## set data
      super$cpp_model$set_data(as.matrix(self$X))
      super$display("  Data have been set.")
      if (!is.null(self$w)) {
        super$cpp_model$set_weights(as.matrix(self$w))
        super$display("  Weights have been set.")
      }
    },
    init_calibrator = function(calibrator_params) {
      ## init calibrator
      private$cpp_calibrator <- calibrators_factory(calibrator_params)
    },
    init_smoother = function(smoother_params) {
      ## init calibrator
      smoother_params$sampling_type <- self$model_traits$sampling_type
      private$cpp_smoother <- regression_models_factory(
        self$domain,
        smoother_params
      )
      ## set locations
      private$cpp_smoother$set_spatial_locations(as.matrix(self$locations))
      super$display("  Locations have been set.")
      ## set calibrator into cpp_model
      super$cpp_model$set_smoother(private$cpp_smoother$get_view())
    },
    ## fit utilities
    save_calibration_results = function() {
      ## save calibrator's results
      self$results$calibrator$lambda_opt <- private$cpp_calibrator$optimum()
      if (private$cpp_calibrator$get_calibration_strategy() == 1) {
        self$results$calibrator$edfs <- as.matrix(private$cpp_calibrator$edfs())
        self$results$calibrator$gcvs <- as.matrix(private$cpp_calibrator$gcvs())
      }
      if (private$cpp_calibrator$get_calibration_strategy() == 2) {
        self$results$calibrator$avg_scores <- as.matrix(private$cpp_calibrator$avg_scores())
        self$results$calibrator$std_scores <- as.matrix(private$cpp_calibrator$std_scores())
        self$results$calibrator$scores <- as.matrix(private$cpp_calibrator$scores())
      }
    }
  )
)
