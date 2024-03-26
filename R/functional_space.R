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

#' R6 Class representing a space of functions defined over a domain
#' @export
.FunctionalSpace <- R6::R6Class("FunctionalSpace",
  private = list(
    basis_ = "ANY",
    mesh_  = "Mesh",
    order_ = numeric(),
    integrate = function(f) {
      if (!is.function(f) && !is.vector(f)) {
          stop("object", deparse(substitute(f)), "is neither a function nor an expansion coefficient vector.")
          if(is.vector(f) && nrow(f) != private$basis_$size()) {
              stop("vector ", deparse(substitute(f)), " has wrong dimensions.")
          }
      }
      if (is.function(f)) { ## recover fem expansion and integrate
          return(private$basis_$integrate(matrix(f(private$basis_$dofs_coords()), ncol = 1)))
      } else {
          return(private$basis_$integrate(f))
      }
    }
  ),
  public = list(
    initialize = function(basis = NA, mesh = NA, order = NA) {
        private$basis_ <- basis
        private$mesh_  <- mesh
        private$order_ <- order
    },
    ## evaluates the basis system on a given set of locations
    eval = function(locations, type = c("pointwise", "areal")) {
      sampling_type <- match(match.arg(type), c("pointwise", "areal"))
      return(private$basis_$eval(cpp_aligned_index(sampling_type), as.matrix(locations)))
    }
  ),
  ## active binding
  active = list(
    size  = function() private$basis_$size(),
    mesh  = function() private$mesh_,
    order = function() private$order_
  )
)

#' @export
FunctionalSpace <- function(mesh, type, ...) {
  args <- list(...)
  if (type == "fe") { ## continous galerkin elements
    order <- if("order" %in% names(args)) args[["order"]] else 1 ## default to linear elements
    if(!(order == 1 || order == 2)) stop("not supported basis order.")
    cpp_backend <- new(
      eval(parse(text = paste("cpp_fe_space_lagrange",
        as.character(get_private(mesh)$local_dim_),
        as.character(get_private(mesh)$embed_dim_),
        as.character(order),
        sep = "_"
      ))),
      get_private(mesh)$mesh_, 0
    )
  }
  if (type == "bs") { ## BSpline basis
    order <- if("order" %in% names(args)) args[["order"]] else 3 ## default to cubic splines
    if(!(order == 3)) stop("not supported basis order.")
    if (get_private(mesh)$local_dim_ != 1 && get_private(mesh)$embed_dim_ != 1) {
      stop(deparse(substitute(mesh)), " is not a 1D interval.")
    }
    cpp_backend <- new(eval(parse(text = paste("cpp_bspline_space",sep = ""))), get_private(mesh)$mesh_, 0)
  }
  .FunctionalSpace$new(
    basis = cpp_backend,
    mesh  = mesh,
    order = order
  )
}

## tensor prodcut system of functional spaces
.TensorProductSpace <- R6::R6Class("TensorProductSpace",
  private = list(
      lhs_basis_ = "FunctionalSpace",
      rhs_basis_ = "FunctionalSpace"
  ),
  public = list(
      initialize = function(lhs_basis_ = NA, rhs_basis_ = NA) {
          private$lhs_basis_ <- lhs_basis_
          private$rhs_basis_ <- rhs_basis_
    },
    ## evaluates the basis system on a given set of locations
    eval = function(lhs_locations, rhs_locations, type = c("pointwise", "areal")) {
        sampling_type <- match(match.arg(type), c("pointwise", "areal"))
        lhs_eval <- private$lhs_basis_$eval(lhs_locations, type)
        rhs_eval <- private$rhs_basis_$eval(rhs_locations, type)
      return(kronecker_tensor_product(lhs_eval, rhs_eval))
    },
    size = function() {
      return(private$lhs_basis_$size() * private$rhs_basis_$size())
    }
  ),
  active = list(
    lhs_mesh = function() private$lhs_basis_$mesh,
    rhs_mesh = function() private$rhs_basis_$mesh      
  )
)

#' @export
`%X%` <- function(x, ...) { UseMethod("%X%", x) }

#' @export
`%X%.FunctionalSpace` <- function(rhs, lhs){
    if(!inherits(lhs, "FunctionalSpace")) stop(deparse(substitute(lhs)), ": not a FunctionalSpace object.")
    return(.TensorProductSpace$new(lhs, rhs))
}
