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
.fe_function <- R6::R6Class(
  "fe_function",
  private = list(
    fe_function_ = NULL, ## cpp backend
    type = character()
  ),
  public = list(
    initialize = function(domain, type, coeff) {
      local_dim = domain$local_dim
      embed_dim = domain$embed_dim
      if (local_dim == 2 && embed_dim == 2) {
        private$fe_function_ = new(cpp_fe_function_2_2_p1, get_private(domain)$mesh_)
      }
      private$type = type
      if (!is.null(coeff)) {
        n_dofs <- private$fe_function_$n_dofs()
        fdapde_assert(nrow(coeff) == n_dofs, "invalid coefficient vector dimensions.")
        private$fe_function_$set_coeff(coeff)
      } else {
        private$fe_function_$set_coeff(matrix(rep(0, times = private$fe_function_$n_dofs()), ncol = 1))
      }
    },
    integral = function(marker = NULL) {
      if (is.null(marker)) {
        marker <- -1
      }
      return(private$fe_function_$cell_integrate_on(marker))
    },
    eval = function(locations) {
      return(private$fe_function_$grid_eval(as.matrix(locations)))
    }
  ),
  active = list(
    n_dofs = function() private$fe_function_$n_dofs(),
    l2_norm = function() private$fe_function_$l2_norm(),
    h1_norm = function() private$fe_function_$h1_norm(),
    l2_squared_norm = function() private$fe_function_$l2_squared_norm(),
    h1_squared_norm = function() private$fe_function_$h1_squared_norm(),
    coeff = function(c) {
      if (missing(c)) {
        return(private$fe_function_$coeff())
      } else {
        private$fe_function_$set_coeff(as.matrix(c))
      }
    }
  )
)

#' @export
fe_function <- function(domain, type, coeff = NULL) {
  ## check domain is of type triangulation

  ## check type is supported
  return(.fe_function$new(
    domain = domain,
    type = type,
    coeff = coeff
  ))
}
