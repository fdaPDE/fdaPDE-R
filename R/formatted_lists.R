## enumerators ----

match.list <- function(options_list, x) {
  if (is.numeric(x)) {
    return(names(options_list)[options_list == x])
  } else {
    return(options_list[[x]])
  }
}

Calibration <- function(x) {
  calibration_options <- list("off" = 0, "gcv" = 1, "kcv" = 2)
  match.list(calibration_options, x)
}

Sampling <- function(x) {
  sampling_options <- list("mesh_nodes" = 0, "pointwise" = 1, "areal" = 2)
  match.list(sampling_options, x)
}

# formatted lists ----

## lambda ----

#' @export
hyperparameters <- function(lambda_D, lambda_T = 0) {
  lambda <- list(
    space = lambda_D,
    time = lambda_T
  )
  return(lambda)
}

## data ----

#' @export
spatial_data <- function(domain, observations,
                         locations = NULL,
                         covariates = NULL) {
  data <- list(
    type = "spatial_data",
    domain = domain,
    locations = locations,
    observations = observations,
    covariates = covariates
  )
  return(data)
}

#' @export
functional_data <- function(domain, X,
                            locations = NULL,
                            w = NULL) {
  data <- list(
    type = "functional_data",
    domain = domain,
    locations = locations,
    X = X,
    w = w
  )
  return(data)
}

## calibrators ----

#' @export
off <- function() {
  off_params <- list(
    ## calibration strategy
    calibration_strategy = "off"
  )
  return(off_params)
}

#' @export
gcv <- function(optimizer = c("grid"),
                edf_computation = c("stochastic", "exact"),
                seed = -1, mc_samples = 100,
                max_iter = NULL, step = NULL,
                tolerance = NULL) {
  gcv_params <- list(
    ## calibration strategy
    calibration_strategy = "gcv",
    ## properties related to the computation of the expected degrees of freedom
    edf_computation = match.arg(edf_computation),
    seed = seed, ## meaningfull only for stochastic edf_computation
    mc_samples = mc_samples, ## meaningfull only for stochastic edf_computation
    ## properties related to the optimization algorithm
    optimizer = match.arg(optimizer),
    max_iter = max_iter,
    step = step,
    tolerance = tolerance
  )
  return(gcv_params)
}

#' @export
kcv <- function(n_folds = 10, shuffle = TRUE, seed = NULL) {
  kcv_params <- list(
    ## calibration strategy
    calibration_strategy = "kcv",
    ## kcv options
    K = n_folds,
    shuffle = shuffle,
    seed = seed
  )
  return(kcv_params)
}

#' @export
calibration <- function(strategy = c("off", "gcv", "kcv"), ...) {
  calibrator_params <- switch(match.arg(strategy),
    "off" = {
      off(...)
    },
    "gcv" = {
      gcv(...)
    },
    "kcv" = {
      kcv(...)
    },
    stop("The selected calibrator does not exist.")
  )
  return(calibrator_params)
}

## smoother ----

#' @export
smoothing <- function(name = c("SRPDE"), penalty = simple_laplacian_penalty(),
                      sampling_type = c("mesh_nodes", "pointwise", "areal")) {
  smoother_params <- list(
    smoother_name = match.arg(name),
    penalty = penalty,
    sampling_type = match.arg(sampling_type)
  )
  return(smoother_params)
}
