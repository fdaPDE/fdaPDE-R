test_that("gcv_grid_edf_exact", {
  unit_square <- MeshUnitSquare(n = 20)
  Vh <- FunctionalSpace(unit_square, type = "fe", order = 1)
  f <- Function(Vh)

  set.seed(42)
  ## data locations
  n_loc <- 800
  x_obs <- runif(min = 0, max = 1, n = n_loc)
  y_obs <- runif(min = 0, max = 1, n = n_loc)
  data_locs <- unit_square$nodes ## cbind(x_obs, y_obs)
  ## generating data
  a1 <- -1.5
  a2 <- 0.4
  z <- function(x, y) {
    a1 * sin(2 * pi * x) * cos(2 * pi * y) + a2 * sin(3 * pi * x) - 2
  }
  exact_data <- z(data_locs[, 1], data_locs[, 2])
  inv.link <- function(x) { ## gamma inv link function
    -1 / x
  }
  mu <- inv.link(exact_data)
  data <- data.frame(
    y = rgamma(length(data_locs[, 1]), shape = mu, scale = 1)
  )
  
  model <- SRPDE(y ~ f, data = data, family = "gamma")
  model$fit(
    calibration = 1e-2
  )
})
