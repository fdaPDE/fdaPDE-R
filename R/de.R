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

.de <- R6::R6Class(
  "de",
  private = list(
    model_ = NULL
  ),
  public = list(
    initialize = function(data, penalty = NULL) {
      private$model_ <- new(cpp_de_2_2, get_private(data$gf__ptr__)$ptr_, penalty)
    },
    fit = function(lambda, optimizer) {
      private$model_$fit(lambda, as.character(optimizer$opt_t), optimizer)
    }
  ),
  active = list(
    density = function() private$model_$density(),
    log_density = function() private$model_$log_density(),
    fitted = function() private$model_$fitted()
  )
)

#' @export
de <- function(data, penalty = NULL) {
  return(.de$new(
    data = data,
    penalty = penalty
  ))
}



