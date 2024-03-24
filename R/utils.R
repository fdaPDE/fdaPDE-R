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

cpp_aligned_index <- function(index) return(index - 1)
r_aligned_index   <- function(index) return(index + 1)

get_private <- function(x) {
    if(!inherits(x, "R6")) stop(deparse(substitute(x)), ": not an R6 class")
    return(x$.__enclos_env__$private)
}

set_private <- function(x, attribute, value) {
    if(!inherits(x, "R6")) stop(deparse(substitute(x)), ": not an R6 class")
    invisible(x$.__enclos_env__$private[[attribute]] <- value)
}
