test_that("fixed_lambda", {
  unit_square <- MeshUnitSquare(n = 20)
  Vh <- FunctionalSpace(unit_square, type = "fe", order = 1)

  ## data locations
  data <- data.frame(
    y = rgamma(length(data_locs[, 1]), shape = mu, scale = 1)
  )
  
  model <- fPCA(data, Vh)
  model$fit(
    calibration = 1e-2
  )
})
