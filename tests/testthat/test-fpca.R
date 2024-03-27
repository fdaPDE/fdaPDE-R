test_that("spaceonly_fixed_lambda", {
  unit_square <- MeshUnitSquare(n = 20)
  Vh <- FunctionalSpace(unit_square, type = "fe", order = 1)

  ## eigenfunctions of laplacian over square, neumann conditions
  square_eigenfunction <- function(locs, n, m) {
    return(as.matrix(cos(n * locs[, 1]) * cos(m * locs[, 2])))
  }
  locations <- unit_square$nodes
  n_locations <- nrow(locations)

  ## not normalized... ok...
  f_1 <- square_eigenfunction(locations, 1 * pi, 1 * pi)
  f_2 <- square_eigenfunction(locations, 1 * pi, 3 * pi)
  f_3 <- square_eigenfunction(locations, 4 * pi, 2 * pi)
  data_range <- max(c(f_1, f_2, f_3)) - min(c(f_1, f_2, f_3))

  ## set scores
  sd_score1 <- 0.4
  sd_score2 <- 0.3
  sd_score3 <- 0.2
  sd_error <- 0.1
  n_subjects <- 50 ## set number of statistical units
  set.seed(4513)

  score1 <- as.matrix(rnorm(n = n_subjects, sd = sd_score1 * data_range))
  score2 <- as.matrix(rnorm(n = n_subjects, sd = sd_score2 * data_range))
  score3 <- as.matrix(rnorm(n = n_subjects, sd = sd_score3 * data_range))

  exact_data <- score1 %*% t(f_1) + score2 %*% t(f_2) + score3 %*% t(f_3)
  error <- rnorm(n = n_subjects * n_locations, sd = sd_error * data_range)

  data <- exact_data + error

  model <- fPCA(data, Vh)
  model$fit(
    ncomp = 3,
    calibration = 1e-2
  )
})

test_that("spaceonly_gcv", {
  unit_square <- MeshUnitSquare(n = 20)
  Vh <- FunctionalSpace(unit_square, type = "fe", order = 1)

  ## eigenfunctions of laplacian over square, neumann conditions
  square_eigenfunction <- function(locs, n, m) {
    return(as.matrix(cos(n * locs[, 1]) * cos(m * locs[, 2])))
  }
  locations <- unit_square$nodes
  n_locations <- nrow(locations)

  ## not normalized... ok...
  f_1 <- square_eigenfunction(locations, 1 * pi, 1 * pi)
  f_2 <- square_eigenfunction(locations, 1 * pi, 3 * pi)
  f_3 <- square_eigenfunction(locations, 4 * pi, 2 * pi)
  data_range <- max(c(f_1, f_2, f_3)) - min(c(f_1, f_2, f_3))

  ## set scores
  sd_score1 <- 0.4
  sd_score2 <- 0.3
  sd_score3 <- 0.2
  sd_error <- 0.1
  n_subjects <- 50 ## set number of statistical units
  set.seed(4513)

  score1 <- as.matrix(rnorm(n = n_subjects, sd = sd_score1 * data_range))
  score2 <- as.matrix(rnorm(n = n_subjects, sd = sd_score2 * data_range))
  score3 <- as.matrix(rnorm(n = n_subjects, sd = sd_score3 * data_range))

  exact_data <- score1 %*% t(f_1) + score2 %*% t(f_2) + score3 %*% t(f_3)
  error <- rnorm(n = n_subjects * n_locations, sd = sd_error * data_range)

  data <- exact_data + error

  model <- fPCA(data, Vh)
  model$fit(
    ncomp = 3,
    calibration = gcv(lambda = 10^seq(-4, -2, by = 0.2))
  )
})

test_that("spacetime_fixed_lambda", {
  unit_square <- MeshUnitSquare(n = 20)
  unit_interval <- MeshUnitInterval(n = 10)
  ## space-time functional space
  Vs <- FunctionalSpace(unit_square, type = "fe", order = 1)
  Vt <- FunctionalSpace(unit_interval, type = "bs", order = 3)
  Vh <- Vs %X% Vt
  
  ## eigenfunctions of laplacian over cube, neumann conditions
  cube_eigenfunction <- function(locs, n, m, h) {
    return(as.matrix(cos(n*locs[,1])*cos(m*locs[,2])*cos(h*locs[,3])))
  }
  set.seed(4513)

  spatial_locations  <- unit_square$nodes
  temporal_locations <- as.matrix(unit_interval$nodes)
  n_space_locations  <- nrow(spatial_locations)
  n_time_locations   <- nrow(temporal_locations)
  n_locations        <- n_space_locations * n_time_locations
  n_stat_units       <- 50

  ## spatio-temporal locations where data are observed
  locations <- cbind(
    rep(spatial_locations[, 1], n_time_locations),
    rep(spatial_locations[, 2], n_time_locations),
    rep(temporal_locations, each = n_space_locations)
  )

  ## true PC function generation
  f_1 <- cube_eigenfunction(locations, 1 * pi, 1 * pi, 2 * pi)
  f_2 <- cube_eigenfunction(locations, 1 * pi, 3 * pi, 2 * pi)
  f_3 <- cube_eigenfunction(locations, 4 * pi, 2 * pi, 3 * pi)
  data_range <- max(c(f_1, f_2, f_3)) - min(c(f_1, f_2, f_3))
  ## true scores generation
  sd_score1 <- 0.4 sd_score2 <- 0.3 sd_score3 <- 0.2 sd_error <- 0.1
  score1 <- as.matrix(rnorm(n = n_stat_units, sd = sd_score1 * data_range))
  score2 <- as.matrix(rnorm(n = n_stat_units, sd = sd_score2 * data_range))
  score3 <- as.matrix(rnorm(n = n_stat_units, sd = sd_score3 * data_range))
  ## data generation
  exact_data <- score1 %*% t(f_1) + score2 %*% t(f_2) + score3 %*% t(f_3)
  error <- rnorm(n = n_stat_units * n_locations, sd = sd_error * data_range)
  data <- exact_data + error ## n_stat_units x n_locations matrix

  ## define model and fit
  model <- fPCA(data, Vh)
  model$fit(
    ncomp = 3,
    calibration = c(1e-2, 1e-2) ## fixed smoothing parameters
  )
  
  model$fit(
    ncomp = 3,
    calibration = gcv(lambda = expand.grid(10^seq(-4, -3, by = 0.2), 10^seq(-4, -3, by = 0.2)))
  )
})
