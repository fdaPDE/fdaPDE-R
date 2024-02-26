## TODOs:
## TODO regression_model_factory: let this factory return also a calibrator


# calibrators ----

# calibrator_init_list
# - name
# - parameters

calibrators_factory <- function(calibrator_init_list) {
  calibrator <- switch(calibrator_init_list$name,
    "off" = {
      new(cpp_off, calibrator_init_list$parameters)
    },
    "gcv" = {
      new(cpp_gcv, calibrator_init_list$parameters)
    },
    "kcv" = {
      new(cpp_kcv, calibrator_init_list$parameters)
    }
  )
  return(calibrator)
}

# models ----

# smoother_init_list
# - name
# - penalty
#   - pde_callable
#   - parameters
# - regularization_type
# - sampling_type
# - parameters
# - calibrator
#   - name
#   - parameters
# - lambda

regression_models_factory <- function(domain, smoother_init_list) {
  ## pde initialization
  pde_callable <- smoother_init_list$penalty$pde_callable
  pde_parameters <- smoother_init_list$penalty$parameters
  pde <- pde_callable(domain, pde_parameters)
  ## model initialization
  cpp_pair <- switch(smoother_init_list$name,
    "SRPDE" = {
      list(
        cpp_model = new(cpp_srpde, pde, Sampling(smoother_init_list$sampling_type), smoother_init_list$parameters),
        cpp_calibrator = calibrators_factory(smoother_init_list$calibrator)
      )
    }
  )
  return(cpp_pair)
}

# center_init_list = smoother_init_list

centering_factory <- function(domain, center_init_list) {
  ## smoother & calibrator initialization
  cpp_pair <- regression_models_factory(domain, center_init_list)
  ## assembly
  cpp_touple <- list(
    cpp_model = new(cpp_center),
    cpp_smoother = cpp_pair$cpp_model,
    cpp_calibrator = cpp_pair$cpp_calibrator
  )
  return(cpp_touple)
}

# functional_model_init_list
# - name
# - penalty
# - regularization_type
# - sampling_type
# - center
# - .... (depending on the specific funcitonal model)
# - parameters
# - CENTER

functional_models_factory <- function(domain, fm_init_list) {
  ## pde initialization
  pde_callable <- fm_init_list$penalty$pde_callable
  parameters <- fm_init_list$penalty$parameters
  pde <- pde_callable(domain, parameters)
  ## model initialization
  functional_model <- switch(fm_init_list$name,
    "fPCA" = {
      ## calibrator recasting
      calibrator <- fm_init_list$solver$calibrator
      calibrator$calibration_strategy <- Calibration(calibrator$name)
      ## model initialization
      switch(fm_init_list$regularization_type,
        "SpaceOnly" = {
          ## statistical model initialization
          cpp_fPCA <- new(
            cpp_fpca_spaceonly,
            pde,
            Sampling(fm_init_list$sampling_type),
            fm_init_list$parameters
          )
          ## set solver module
          cpp_fPCA$set_solver(
            SolutionPolicy(fm_init_list$solver$policy),
            fm_init_list$solver$parameters,
            calibrator
          )
          return(cpp_fPCA)
        }
      )
    }
  )
  return(functional_model)
}
