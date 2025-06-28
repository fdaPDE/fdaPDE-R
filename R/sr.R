## This file is part of fdaPDE, a C++ library for physics-informed
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

#' A triangulated spatial domain
#'
.sr <- R6::R6Class(
  "sr",
  private = list(
    model_ = NULL ## cpp backend
  ),
  public = list(
    initialize = function(formula, data, penalty) {
      f_symbol <- as.formula(formula)[[3]]
      f <- get(f_symbol, envir = globalenv())
      fe_type <- get_private(f)$type
      if (is.null(penalty)) {
        ## fallback to laplacian penalty
        private$model_ = new(cpp_sr_2_2, formula, get_private(data)$geoframe_, penalty)
      } else {
        domain <- get_private(get_private(data)$triangulation_)$mesh_
        params <- list()
        quad_nodes <- matrix()
        if (is.function(penalty$K) || is.function(penalty$b) || is.function(penalty$c) || is.function(penalty$u)) {
          quad_nodes <- fe_simplex_quad_nodes(domain, fe_type)
        }
        n_quad_nodes <- new(cpp_fe_space_2_2_p1, domain)$n_quad_nodes()
        embed_dim <- 2

        ## diffusion
        if (is.null(penalty$K)) {
          params$K = matrix(rep(0, times = (embed_dim * embed_dim) * n_quad_nodes), nrow = n_quad_nodes)
        } else {
          fdapde_assert(is.matrix(K) || is.function(K), "Invalid diffusion field type.")
          if (is.matrix(K)) {
            fdapde_assert(nrow(K) == embed_dim && ncol(K) == embed_dim, "Not a square matrix.")
            params$K = matrix(rep(c(K), times = n_quad_nodes), nrow = n_quad_nodes, byrow = TRUE)
          } else {
            q <- quad_nodes[1, ]
            fdapde_assert(
              nrow(K(q)) == embed_dim && ncol(K(q)) == embed_dim,
              "Not evaluates to a square matrix."
            )
            params$K = K(quad_nodes)
          }
        }
        ## transport
        if (is.null(penalty$b)) {
          params$b = matrix(rep(0, times = embed_dim * n_quad_nodes), nrow = n_quad_nodes)
        } else {
          fdapde_assert(is.matrix(b) || is.vector(b) || is.function(b), "Invalid transport field type.")
          if (is.matrix(b) || is.vector(b)) {
            fdapde_assert(nrow(b) == embed_dim && ncol(b) == 1, "Not a vector.")
            params$b = matrix(rep(c(b), times = n_quad_nodes), nrow = n_quad_nodes, byrow = TRUE)
          } else {
            q <- quad_nodes[1, ]
            fdapde_assert(nrow(b(q)) == embed_dim && ncol(b(q)) == 1, "Not evaluates to a vector.")
            params$b = b(quad_nodes)
          }
        }
        ## reaction
        if (is.null(penalty$c)) {
          params$c = matrix(rep(0, times = n_quad_nodes), nrow = n_quad_nodes)
        } else {
          fdapde_assert(is.numeric(c) || is.function(c), "Invalid reaction field type.")
          if (is.numeric(c)) {
            params$c = matrix(rep(c, times = n_quad_nodes), nrow = n_quad_nodes)
          } else {
            q <- quad_nodes[1, ]
            fdapde_assert(is.numeric(c(q)), "Not evaluates to a scalar.")
            params$c = c(quad_nodes)
          }
        }
        ## force
        if (is.null(penalty$u)) {
          params$u = matrix(rep(0, times = n_quad_nodes), nrow = n_quad_nodes)
        } else {
          fdapde_assert(is.numeric(u) || is.function(u), "Invalid forcing field type.")
          if (is.numeric(u)) {
            params$u = matrix(rep(u, times = n_quad_nodes), nrow = n_quad_nodes)
          } else {
            q <- quad_nodes[1, ]
            fdapde_assert(is.numeric(u(q)), "Not evaluates to a scalar.")
            params$c = c(quad_nodes)
          }
        }
        private$model_ = new(cpp_sr_2_2, formula, get_private(data)$geoframe_, params)
      }
    },
    fit = function(lambda = NULL, calibrator = NULL) {
      fdapde_assert(!is.null(calibrator) && !is.null(lambda), "Unable to select smoothing level.")
      if (is.null(calibrator)) {
        private$model_$fit(lambda)
      } else {
        private$model_$fit_gcv(calibrator)
      }
    }
  ),
  active = list(
    f = function() private$model_$f(),
    beta = function() private$model_$beta(),
    fitted = function() private$model_$fitted()
  )
)

#' @export
sr <- function(formula, data, penalty = NULL) {
  return(.sr$new(
    formula = deparse(formula),
    data = data,
    penalty = penalty
  ))
}

#' @export
fe_elliptic <- function(K = NULL, b = NULL, c = NULL, u = NULL) {
  fdapde_assert(!(is.null(K) && is.null(b) && is.null(c) && is.null(u)))
  return(list(K = K, b = b, c = c, u = u))
}
