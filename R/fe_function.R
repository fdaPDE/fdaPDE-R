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

#' R6 Class representing a function
#'
.fe_function <- R6::R6Class(
  "fe_function",
  private = list(
    fe_function_ = NULL, ## cpp backend
    mesh_ = NULL, ## geometry r6 handler
    type = character()
  ),
  public = list(
    #' @description
    #' Creates a new function object defined over a spatial domain.
    #'
    #' @param domain A triangulation of the spatial domain, created by [triangulation()].
    #' @param type A character string indicating the order of the finite element space.
    #' @param coeff A numeric vector of basis expansion coefficients.
    #'
    #' @return A function object defined over the given domain.
    initialize = function(domain, type, coeff) {
      private$mesh_ <- domain
      local_dim = domain$local_dim
      embed_dim = domain$embed_dim
      if (local_dim == 2 && embed_dim == 2) {
        private$fe_function_ = new(cpp_fe_function_2_2_p1, get_private(domain)$mesh_)
      }
      private$type = type
      if (!is.null(coeff)) {
        n_dofs <- private$fe_function_$n_dofs()
        fdapde_assert(length(coeff) == n_dofs, "invalid coefficient vector dimensions.")
        private$fe_function_$set_coeff(coeff)
      } else {
        private$fe_function_$set_coeff(matrix(rep(0, times = private$fe_function_$n_dofs()), ncol = 1))
      }
    },
    #' @description
    #' Computes the integral of the function.
    #'
    #' @param marker An integer denoting the subdomain over which to compute the integral.
    integral = function(marker = NULL) {
      if (is.null(marker)) {
        marker <- -1
      }
      return(private$fe_function_$cell_integrate_on(marker))
    },
    #' @description
    #' Evaluates the function over a set of points.
    #'
    #' @param locations A matrix containing the evaluatio points.
    eval = function(locations) {
      return(private$fe_function_$grid_eval(as.matrix(locations)))
    }
  ),
  active = list(
    #' @field n_dofs (`integer(1)`)\cr
    #' The number of degrees of freedom of the function.
    n_dofs = function() private$fe_function_$n_dofs(),
    #' @field l2_norm (`numeric(1)`)\cr
    #' The L2 norm of the function
    l2_norm = function() private$fe_function_$l2_norm(),
    #' @field h1_norm (`numeric(1)`)\cr
    #' The H1 norm of the function
    h1_norm = function() private$fe_function_$h1_norm(),
    #' @field l2_squared_norm (`numeric(1)`)\cr
    #' The L2 squared norm of the function
    l2_squared_norm = function() private$fe_function_$l2_squared_norm(),
    #' @field h1_squared_norm (`numeric(1)`)\cr
    #' The H1 squared norm of the function
    h1_squared_norm = function() private$fe_function_$h1_squared_norm(),
    #' @field coeff (`numeric(n_dofs)`)\cr
    #' Gets or sets the vector of basis expansion coefficients.
    #' If a value is provided, it assigns the given vector; otherwise, it returns the current coefficients.
    coeff = function(c) {
      if (missing(c)) {
        return(private$fe_function_$coeff())
      } else {
        private$fe_function_$set_coeff(as.matrix(c))
      }
    },
    #' @field geometry A \code{triangulation} object that defines the domain over which the finite element function is defined.
    geometry = function() private$mesh_
  )
)

#' Create a function object
#'
#' @param domain A \code{triangulation_2_2} R6 class object created by [triangulation()].
#' @param type A character string indicating the type of basis functions. Use "P1" or "P2" to select finite elements of order 1 or 2, respectively.
#' @param coeff A numeric vector of basis expansion coefficients. Defaults to \code{NULL}.
#'
#' @return An R6 object representing a function belonging to the specified function space.
#'
#' @export
#' @examples
#' \dontrun{
#' library(fdaPDE2)
#' unit_square <- MeshUnitSquare(n = 20)
#' f <- Function(unit_square, type="fe")
#' }
fe_function <- function(domain, type, coeff = NULL) {
  ## check domain is of type triangulation

  ## check type is supported
  return(.fe_function$new(
    domain = domain,
    type = type,
    coeff = coeff
  ))
}

#' Plot a finite element function
#'
#' Plots a \code{fe_function} object over its domain.
#'
#' @param x An object of class \code{fe_function}.
#' @param palette A function that takes a single integer (e.g., the number of levels or bins) and returns a character vector specifying a color palette.
#' @export
plot.fe_function <- function(x, palette = NULL, ...) {
  n_col <- 100
  if (is.null(palette)) {
    palette_ <- colorRampPalette(colors = c("lightyellow", "darkred"))(n_col)
  } else {
    palette_ <- palette(n_col)
  }

  nodes <- get_private(x)$mesh_$nodes
  x_grid <- seq(min(nodes[, 1]), max(nodes[, 1]), length.out = 250)
  y_grid <- seq(min(nodes[, 2]), max(nodes[, 2]), length.out = 250)
  xy_grid <- expand.grid(x_grid, y_grid)
  ## evaluate fe_function at fine grid
  vals <- x$eval(xy_grid)

  col <- palette_[as.numeric(cut(vals, breaks = n_col))]
  par(mar = c(1, 1, 1, 1))
  plot(
    xy_grid[, 1],
    xy_grid[, 2],
    xlab = "",
    ylab = "",
    pch = 15,
    col = col,
    asp = 1,
    cex = .6
  )
}
