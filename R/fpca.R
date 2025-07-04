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

.fpca <- R6::R6Class(
  "fpca",
  private = list(
    model_ = NULL ## cpp backend
  ),
  public = list(
    initialize = function(column, data) {
      private$model_ = new(cpp_fpca_laplace_2_2, column, get_private(data$gf__ptr__)$ptr_)
    },
    fit = function(npc = NULL, calibrator = NULL) {
      fdapde_assert(!is.null(calibrator), "Unable to select smoothing level.")
      private$model_$fit(npc, calibrator)
    }
  ),
  active = list(
    loadings = function() private$model_$loadings(),
    scores = function() private$model_$scores(),
    pcs = function() private$model_$pcs()
  )
)

#' @export
fpca <- function(column, data) {
  return(.fpca$new(
    column = column,
    data = data
  ))
}
