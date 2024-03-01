## TODOs:
## TODO regression_model_factory: let this factory return also a calibrator


# calibrators ----

# calibrator_init_list
# - name
# - lambda
# - parameters

calibrators_factory <- function(calibrator_init_list) {
  ## calibrator initialization
  calibrator <- switch(calibrator_init_list$name,
    "off" = {
      new(
        ## cpp module
        cpp_off,
        ## contructor's arguments
        calibrator_init_list$parameters
      )
    },
    "gcv" = {
      new(
        ## cpp module
        cpp_gcv,
        ## contructor's arguments
        calibrator_init_list$parameters
      )
    },
    "kcv" = {
      new(
        ## cpp module
        cpp_kcv,
        ## contructor's arguments
        calibrator_init_list$parameters
      )
    }
  )
  ## calibrator configuration
  calibrator$configure_calibrator(calibrator_init_list$lambda)
  return(calibrator)
}

# models ----

# model_traits
# - regularization_type
# - sampling_type

# smoother_init_list
# - name
# - penalty
#   - pde_callable
#   - parameters
# - parameters
# - calibrator
#   - name
#   - parameters
# - lambda

regression_models_factory <- function(domain, model_traits, smoother_init_list) {
  ## pde initialization
  pde_callable <- smoother_init_list$penalty$pde_callable
  pde_parameters <- smoother_init_list$penalty$parameters
  pde <- pde_callable(domain, pde_parameters)
  ## model initialization
  cpp_pair <- switch(smoother_init_list$name,
    "SRPDE" = {
      list(
        cpp_model = new(
          ## cpp module
          cpp_srpde, pde,
          ## contructor's arguments
          Sampling(model_traits$sampling_type),
          smoother_init_list$parameters
        ),
        cpp_calibrator = calibrators_factory(
          smoother_init_list$calibrator
        )
      )
    }
  )
  return(cpp_pair)
}

# center_init_list = smoother_init_list

centering_factory <- function(domain, model_traits, center_init_list) {
  ## smoother & calibrator initialization
  cpp_pair <- regression_models_factory(domain, model_traits, center_init_list)
  ## assembly
  cpp_touple <- list(
    cpp_model = new(cpp_center),
    cpp_smoother = cpp_pair$cpp_model,
    cpp_calibrator = cpp_pair$cpp_calibrator
  )
  return(cpp_touple)
}

# model_traits
# - regularization_type
# - sampling_type

# functional_model_init_list
# - name
# - penalty
# - center
# - .... (depend on the specific funcitonal model)
# - parameters
# - CENTER

functional_models_factory <- function(domain, model_traits, fm_init_list) {
  ## pde initialization
  pde_callable <- fm_init_list$penalty$pde_callable
  parameters <- fm_init_list$penalty$parameters
  pde <- pde_callable(domain, parameters)
  ## model initialization
  functional_model <- switch(fm_init_list$name,
    "fPCA" = {
      ## model initialization
      switch(model_traits$regularization_type,
        "SpaceOnly" = {
          ## statistical model initialization
          cpp_fPCA <- new(
            cpp_fpca_spaceonly,
            ## contructor's arguments
            pde,
            Sampling(model_traits$sampling_type),
            fm_init_list$parameters
          )
          return(cpp_fPCA)
        }
      )
    },
    "fPLS" = {
      ## model initialization
      switch(model_traits$regularization_type,
        "SpaceOnly" = {
          ## statistical model initialization
          cpp_fPLS <- new(
            ## cpp module
            cpp_fpls_spaceonly,
            ## contructor's arguments
            pde,
            Sampling(model_traits$sampling_type),
            fm_init_list$parameters
          )
          return(cpp_fPLS)
        }
      )
    }
  )
  return(functional_model)
}
