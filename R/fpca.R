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

#' An R6 class encapsulating the smooth functional principal component analysis.
#'
#' @rdname fpca
#' @order 2
.fpca <- R6::R6Class(
  "fpca",
  private = list(
    model_ = NULL ## cpp backend
  ),
  public = list(
    #' @description
    #' Creates a new \code{fpca} object.
    #' 
    #' @param column A string specifying the name of the data field in the geoframe layer.
    #' @param data A \code{geoframe} containing both the triangulation of the domain and the associated data (see also [geoframe()]).
    initialize = function(column, data) {
      private$model_ = new(cpp_fpca_laplace_2_2, column, get_private(data$gf__ptr__)$ptr_)
    },
    #' @description
    #' Fits the statistical model.
    
    #' @param npc integer denoting the number of principal components to be extracted.
    #' @param calibrator A calibrator object. Available calibrators include [gcv()].
    fit = function(npc = NULL, calibrator = NULL) {
      fdapde_assert(!is.null(calibrator), "Unable to select smoothing level.")
      private$model_$fit(npc, calibrator)
    }
  ),
  active = list(
    #' @field loadings A matrix containing the fPCs evaluated at the sampling locations.
    loadings = function() private$model_$loadings(),
    #' @field scores A matrix containing the scores for each fPC.
    scores = function() private$model_$scores(),
    #'@field pcs A matrix containing the coefficients of the basis expansions of the fPCs.
    pcs = function() private$model_$pcs()
  )
)

#' Create an \code{fpcs} object
#'
#' @param column A string specifying the name of the data field in the geoframe layer.
#' @param data A \code{geoframe} containing both the triangulation of the domain and the associated data (see also [geoframe()]).
#' @rdname fpca
#' @order 1
#' @references Lila, E., Aston, J.A.D.,  Sangalli, L.M., 2016a. Smooth Principal Component Analysis over two-dimensional
#' manifolds with an application to neuroimaging. Ann. Appl. Stat., 10(4), pp. 1854-1879.
#' @export
fpca <- function(column, data) {
  return(.fpca$new(
    column = column,
    data = data
  ))
}
