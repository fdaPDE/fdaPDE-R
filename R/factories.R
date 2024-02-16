calibrators_factory <- function(calibrator_params) {
  ## calibration startegy
  calibration_strategy <- calibrator_params$calibration_strategy
  ## calibrator initialization
  calibrator <- switch(calibration_strategy,
    "off" = {
      new(cpp_off, calibrator_params)
    },
    "gcv" = {
      new(cpp_gcv, calibrator_params)
    },
    "kcv" = {
      new(cpp_kcv, calibrator_params)
    }
  )
  return(calibrator)
}

regression_models_factory <- function(domain, smoother_params) {
  ## regression model parameters
  smoother_name <- smoother_params$smoother_name
  sampling_type <- Sampling(smoother_params$sampling_type)
  ## pde initialization
  pde_callable <- smoother_params$penalty$pde_callable
  pde_parameters <- smoother_params$penalty$pde_parameters
  pde <- pde_callable(domain, pde_parameters)
  ## model initialization
  model <- switch(smoother_name,
    "SRPDE" = {
      new(cpp_srpde, pde, sampling_type)
    }
  )
  return(model)
}
