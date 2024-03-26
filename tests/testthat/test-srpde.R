test_that("gcv_grid_edf_exact", {
  unit_square <- MeshUnitSquare(n = 20)
  Vh <- FunctionalSpace(unit_square, type = "fe", order = 1)
  f <- Function(Vh)

  ## test function
  test_function <- function(x, y, z = 1) {
    coe <- function(x, y) 1 / 2 * sin(5 * pi * x) * exp(-x^2) + 1
    return(sin(2 * pi * (coe(y, 1) * x * cos(z - 2) - y * sin(z - 2))) *
      cos(2 * pi * (coe(y, 1) * x * cos(z - 2 + pi / 2) + coe(x, 1) * y * sin((z - 2) * pi / 2))))
  }
  exact_data <- as.matrix(test_function(unit_square$nodes[, 1], unit_square$nodes[, 2]), ncol = 1)
  # simulate data by adding normal error
  set.seed(7893475)
  data <- data.frame(
    y = exact_data + rnorm(nrow(exact_data), mean = 0, sd = 0.05 * abs(diff(range(exact_data))))
  )
  
  model <- SRPDE(y ~ f, data = data)
  model$fit(
    calibration = gcv(lambda = 10^seq(-6, -3, length = 20), optimizer = "grid", edf_evaluation = "exact")
  )
})
