## TODOs:
## TODO: possible fault when setting a centering device in functional models
##       - NULL: default values are used
##       - logic: default values are used if true
##       - other: used as centering device
##       In this lasta case we have to check that it is a proper centering
##       device initialization list


# enumerators ----

match.list <- function(options_list, x) {
  if (is.numeric(x)) {
    return(names(options_list)[options_list == x])
  } else {
    return(options_list[[x]])
  }
}

Calibration <- function(x) {
  calibration_options <- list(
    "off" = 0,
    "gcv" = 1,
    "kcv" = 2
  )
  match.list(calibration_options, x)
}

Sampling <- function(x) {
  sampling_options <- list(
    "mesh_nodes" = 0,
    "pointwise" = 1,
    "areal" = 2
  )
  match.list(sampling_options, x)
}

Regularization <- function(x) {
  regularization_options <- list(
    "SpaceOnly" = 0,
    "SpaceTimeSeparable" = 1,
    "SpaceTimeParabolic" = 2
  )
  match.list(regularization_options, x)
}

SolutionPolicy <- function(x) {
  solution_policy_options <- list(
    "sequential" = 0,
    "monolithic" = 1
  )
  match.list(solution_policy_options, x)
}


# hyperparameters ----

#' @export
hyperparameters <- function(space, time = 0) {
  lambda <- list(
    space = space,
    time = time
  )
  return(lambda)
}

# data ----

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

# calibrators ----

#' @export
off <- function() {
  ## calibration strategy specific parameters
  off_params <- NULL
  ## init list assembly
  off_init_list <- list(
    ## calibration strategy name
    name = "off",
    ## parameters
    parameters = off_params
  )
  return(off_init_list)
}

#' @export
gcv <- function(optimizer = c("grid"),
                edf_computation = c("stochastic", "exact"),
                seed = -1, mc_samples = 100,
                max_iter = NULL, step = NULL,
                tolerance = NULL) {
  ## calibration strategy specific parameters
  gcv_params <- list(
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
  ## init list assembly
  gcv_init_list <- list(
    ## calibration strategy name
    name = "gcv",
    ## parameters
    parameters = gcv_params
  )
  return(gcv_init_list)
}

#' @export
kcv <- function(n_folds = 10, shuffle = TRUE, seed = NULL) {
  ## calibration strategy specific parameters
  kcv_params <- list(
    K = n_folds,
    shuffle = shuffle,
    seed = seed
  )
  ## init list assembly
  kcv_init_list <- list(
    ## calibration strategy name
    name = "kcv",
    ## parameters
    parameters = kcv_params
  )
  return(kcv_init_list)
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


# models ----


## Regularized SVD ----

### sequential ----

sequential_params <- function(tolerance = 1e-6, max_iter = 20, seed = NULL) {
  parameters <- list(
    tolerance = tolerance,
    max_iter = max_iter,
    seed = seed
  )
  return(parameters)
}

#' @export
sequential <- function(calibrator = off(), lambda = hyperparameters(space = 1e-4), ...) {
  sequential_init_list <- list(
    policy = "sequential",
    parameters = sequential_params(...),
    calibrator = calibrator,
    lambda = lambda
  )
  return(sequential_init_list)
}

### monolithic ----

monolithic_params <- function() {
  return(NULL)
}

#' @export
monolithic <- function(calibrator = off(), lambda = hyperparameters(space = 1e-4)) {
  ## check calibrator availability
  if (calibrator$name != "off") {
    stop("Only off calibrator is available for this solution policy")
  }
  ##
  monolithic_init_list <- list(
    policy = "monolithic",
    parameters = monolithic_params(),
    calibrator = calibrator,
    lambda = lambda
  )
  return(monolithic_init_list)
}

### generic ---

#' @export
RSVD <- function(policy = c("sequential", "monolithic"), ...) {
  RSVD_params <- switch(match.arg(policy),
    "sequential" = {
      sequential(...)
    },
    "monolithic" = {
      monolithic(...)
    },
    stop("The selected policy does not exist.")
  )
  return(RSVD_params)
}


## spatial ----


### regression ----

#### SRPDE ----

SRPDE_params <- function() {
  return(NULL)
}

#### generic ----

#' @export
smoothing <- function(name = c("SRPDE"),
                      penalty = simple_laplacian_penalty(),
                      regularization_type = c(
                        "SpaceOnly",
                        "SpaceTimeSeparable"
                      ),
                      sampling_type = c(
                        "mesh_nodes",
                        "pointwise",
                        "areal"
                      ),
                      calibrator = off(),
                      lambda = hyperparameters(space = 1e-4),
                      ...) {
  ## method specific parameters initialization
  parameters <- switch(match.arg(name),
    "SRPDE" = {
      SRPDE_params(...)
    },
    stop("The selected smoother does not exist.")
  )
  ## init list assembly
  smoother_init_list <- list(
    name = match.arg(name),
    ## regularization
    penalty = penalty,
    regularization_type = match.arg(regularization_type),
    sampling_type = match.arg(sampling_type),
    ## method specific parameters
    parameters = parameters,
    ## calibration
    calibrator = calibrator,
    lambda = lambda
  )
  return(smoother_init_list)
}


## functional ----


### fCentering ----

## this method has a slightly different behavior because in the library
## it is not impemented as a method but as a function (that wraps a smoothing method)
## the goal of this function is to guarantee the consistency of
## the interfaces on both the user (R) and library (R-cpp) sides
centering <- function(penalty = simple_laplacian_penalty(),
                      smoother = smoothing(),
                      calibrator = off(),
                      lambda = hyperparameters(space = 1e-4)) {
  ## overwriting smoother defaults with the required ones
  ## necessary to guarantee the consistency with the other statistical models
  ## (each statistical model is responsible of its own calibration and lambda)
  ## while providing a clean interface to the final user
  smoother$penalty <- penalty
  smoother$calibrator <- calibrator
  smoother$lambda <- lambda
  center_init_list <- smoother
  return(center_init_list)
}


### fPCA ----

fPCA_params <- function(n_pc = 3) {
  parameters <- list(
    n_pc = n_pc
  )
  return(parameters)
}

fPCA_init <- function(penalty = simple_laplacian_penalty(),
                      regularization_type = c(
                        "SpaceOnly",
                        "SpaceTimeSeparable"
                      ),
                      sampling_type = c(
                        "mesh_nodes",
                        "pointwise",
                        "areal"
                      ),
                      center = NULL, ## also defaults centering lambda
                      solver = sequential(), ## also defaults solver's lambda
                      ...) {
  fPCA_init_list <- list(
    name = "fPCA",
    ## regularization
    penalty = penalty,
    regularization_type = match.arg(regularization_type),
    sampling_type = match.arg(sampling_type),
    ## centering
    center = ifelse(
      is.null(center),
      centering(), ## default centering method
      ifelse(
        is.logical(center),
        centering(), ## if a logic value was provided it sets the default centering device
        center ## if a centering device was provided it uses that
      )
    ),
    ## analsysis
    solver = solver,
    ## method parameters
    parameters = fPCA_params(...),
    ## options
    CENTER = ifelse(
      is.null(center),
      TRUE, ## fPCA defaults center to TRUE
      ifelse(is.logical(center),
        center, ## if a logic was provided it is used as option
        TRUE ## if a calibrator was provided it sets the center option to true
      )
    )
  )
  return(fPCA_init_list)
}

### fPLS ----

fPLS_params <- function(n_comp = 3) {
  parameters <- list(
    n_comp = n_comp
  )
  return(parameters)
}

fPLS_init <- function(penalty = simple_laplacian_penalty(),
                      regularization_type = c(
                        "SpaceOnly",
                        "SpaceTimeSeparable"
                      ),
                      sampling_type = c(
                        "mesh_nodes",
                        "pointwise",
                        "areal"
                      ),
                      center = NULL, ## also defaults centering's lambda
                      rsvd = sequential(), ## also defaults covariance maximization step's lambda
                      smoother = smoothing(), ## also defaults regression step's lambda
                      CENTER = NULL,
                      ...) {
  fPLS_init_list <- list(
    name = "fPCA",
    ## regularization
    penalty = penalty,
    regularization_type = match.arg(regularization_type),
    sampling_type = match.arg(sampling_type),
    ## centering
    center = ifelse(
      is.null(center),
      centering(), ## default centering method
      ifelse(
        is.logical(center),
        centering(),
        center
      )
    ),
    ## analysis
    rsvd = rsvd,
    smoother = smoother,
    ## method parameters
    parameters = fPLS_params(...),
    ## options
    CENTER = ifelse(
      is.null(center),
      TRUE, ## fPCA defaults center to TRUE
      ifelse(is.logical(center),
        center,
        TRUE
      )
    )
  )
  return(fPCA_init_list)
}


### generic ----

functional_model <- function(name = c("fPCA", "fPLS"),
                             penalty = simple_laplacian_penalty(),
                             regularization_type = c(
                               "SpaceOnly",
                               "SpaceTimeSeparable"
                             ),
                             sampling_type = c(
                               "mesh_nodes",
                               "pointwise",
                               "areal"
                             ),
                             center = NULL,
                             solver = sequential(),
                             ...) {
  functional_model_init_list <- switch(match.arg(name),
    "fPCA" = {
      fPCA_init(
        penalty = penalty,
        regularization_type = match.arg(regularization_type),
        sampling_type = match.arg(sampling_type),
        solver = solver,
        center = center,
        ...
      )
    },
    "fPLS" = {
      fPLS_init(
        penalty = penalty,
        regularization_type = match.arg(regularization_type),
        sampling_type = match.arg(sampling_type),
        rsvd = solver,
        center = center,
        ...
      )
    },
    stop("The selected functional model does not exist.")
  )
  return(functional_model_init_list)
}
