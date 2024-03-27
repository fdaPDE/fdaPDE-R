## This file is part of fdaPDE, an R library for physics-informed
## spatial and functional data analysis.

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

.fPCA <- R6::R6Class(
  "fPCAModel",
  private = list(
    model_ = NULL
  ),
  public = list(
    initialize = function(data, basis = NULL) {
      if (is.null(basis)) stop("missing functional basis.")
      if (!inherits(basis, "TensorProductSpace")) {
        ## instantiate laplacian, homogeneous neumann and dirichlet BCs, null forcing term
        f <- Function(basis)
        Lf <- -laplace(f)
        u <- function(nodes) matrix(rep(0, times = nrow(nodes)), ncol = 1)
        penalty <- PDE(Lf, u)
        ## TODO: need to derive type of sampling
        sampling <- 0
        ## instantiate model cpp backend
        private$model_ <- new(cpp_fpca_space, get_private(penalty)$pde_, sampling)
      } else {
        ## penalty in space
        Vs <- basis$lhs
        fs <- Function(Vs)
        spatial_penalty <- PDE(-laplace(fs))
        ## penalty in time
        Vt <- basis$rhs
        if(Vt$mesh$local_dim != 1 && Vt$mesh$embed_dim != 1) stop("wrong dimensions for time domain.")
        ft <- Function(Vt)
        temporal_penalty <- PDE(-bilaplace(ft)) ## TODO: questa è una presa in giro...

        ## TODO: need to derive type of sampling
        sampling <- 0
        ## instantiate model cpp backend
        private$model_ <- new(
          cpp_fpca_spacetime,
          get_private(spatial_penalty)$pde_, get_private(temporal_penalty)$pde_, sampling
        )
      }
      private$model_$set_observations(as.matrix(data))
    },
    fit = function(solver = c("sequential", "monolithic"), ncomp = 3, tolerance = 1e-6, max_iter = 20, ...) {
      args <- list(...)
      fit_data <- list()
      fit_data$tolerance <- tolerance
      fit_data$max_iter <- max_iter
      fit_data$ncomp <- ncomp
      fit_data[["solver"]] <- match.arg(solver)

      ## select calibration type
      if (is.numeric(args$calibration)) {
        fit_data$calibration <- "off"
        fit_data$lambda <- as.matrix(args$calibration)
      } else {
        if (is.list(args$calibration)) {
          fit_data <- c(fit_data, args$calibration)
        } else {
          stop("invalid type bounded to calibration argument.")
        }
      }
      private$model_$fit(fit_data)
    }
  ),
  ## active bindings
  active = list(
    scores   = function() private$model_$scores(),
    loadings = function() private$model_$loadings()
  )
)

## space-only regression model family
fPCA <- function(data, basis = NULL) {
  model <- .fPCA$new(data, basis)
  return(model)
}
