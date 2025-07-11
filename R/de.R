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

#' An R6 class encapsulating the nonparametric density estimation method with partial differential equation regularization.
#'
#' @rdname de
#' @order 2
.ppe <- R6::R6Class(
  "ppe",
  private = list(
    model_ = NULL
  ),
  public = list(
    #' @description
    #' Creates a new \code{de} object.
    #' 
    #' @param data A \code{geoframe} containing the triangulation of the domain and the data locations (see also [geoframe()]).
    #' @param penalty A penalty object returned by [fe_elliptic()]. Default is \code{NULL}, which corresponds to using
    #'        the Laplacian operator in the penalty term (i.e., isotropic smoothing).
    initialize = function(data, penalty = NULL) {
      private$model_ <- new(cpp_de_2_2, get_private(data$gf__ptr__)$ptr_, penalty)
    },
    #' @description
    #' Fits the statistical model.
    #'
    #' @param lambda The smoothing parameter. Default is \code{NULL}.
    #' @param optimizer An Optimizer object. Available optimizers include [gradient_descent(), newton_fd(), bfgs()].
    fit = function(lambda, optimizer) {
      private$model_$fit(lambda, as.character(optimizer$opt_t), optimizer)
    }
  ),
  active = list(
    #' @field density The estimated density function.
    density = function() private$model_$density(),
    #' @field log_density The estimated log density function.
    log_density = function() private$model_$log_density(),
    #' @field fitted A numeric vector containing the estimated density values at the data locations.
    fitted = function() private$model_$fitted()
  )
)

#' Create a \code{de} object
#'
#' @param data A \code{geoframe} containing the triangulation of the domain and the data locations (see also [geoframe()]).
#' @param penalty A penalty object returned by [fe_elliptic()]. Default is \code{NULL}, which corresponds to using
#'        the Laplacian operator in the penalty term (i.e., isotropic smoothing).
#' @rdname de
#' @order 1
#' @references
#' \itemize{
#'    \item{Ferraccioli, F., Arnone, E., Finos, L., Ramsay, J. O., Sangalli, L. M. (2021). Nonparametric density estimation over complicated domains.
#'          Journal of the Royal Statistical Society: Series B (Statistical Methodology), 83(2), 346-368.}
#'    \item{Arnone, E., Ferraccioli, F., Pigolotti, C., Sangalli, L.M. (2021), A roughness penalty approach to estimate densities over two-dimensional manifolds, 
#'    Computational Statistics and Data Analysis, 174, 107527..}
#' }
#' 
#' @export
ppe <- function(data, penalty = NULL) {
  return(.ppe$new(
    data = data,
    penalty = penalty
  ))
}



