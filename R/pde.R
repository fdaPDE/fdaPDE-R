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

pde_type <- list("laplacian" = 1, "elliptic" = 2, "parabolic" = 3)
classify <- function(x, ...) UseMethod("classify", x)
classify.SymbolicExpression <- function(operator) {
  expr <- gsub("<(.*)>\\*laplace\\(<(.*)>\\)|laplace\\(<(.*)>\\)|[+|-]| ", "", operator$expr)
  if (expr == "") { ## detected as simple laplacian: \mu*laplace()
    return(pde_type$laplacian)
  }
  expr <- gsub(
    "div\\(<(.*)>\\*grad\\(<(.*)>\\)\\)|inner\\(<(.*)>,grad\\(<(.*)>\\)\\)|[<(.*)>\\*]?<(.*)>| [+|-] ",
    "", expr
  )
  if (expr == "") { ## detected as elliptic: [+|-]div(K*grad(f))[+|-]inner(b,grad(f))[+|-]c*f
    return(pde_type$elliptic)
  } else {
    if (regexpr("dt", expr)[1] != -1) { ## detected as parabolic
      return(pde_type$parabolic)
    }
    stop("Unrecognized differential operator.")
  }
}

diff_coeff_of <- function(operator) {
  if (!is.symbolic_expression(operator)) stop("Object ", deparse(substitute(operator)), " is not a valid operator.")
  if (classify(operator) == pde_type$laplacian) {
      r <- regexpr("<(.*)>\\*laplace\\(<(.*)>\\)", operator$expr)
      if (r[1] == -1) { ## special sintax for [+|-]laplacian(f)
          if(regexpr("-laplace\\(<(.*)>\\)", operator$expr)[1] != -1) return (-1.0)
          else return(1.0)
      }
      else { ## matched as [+|-] mu * laplace(f)
          match <- regmatches(operator$expr, r)
          return(operator$symbol_table[[substr(match, 1, 12)]])
      }
  }
  ## general diffusion operator [+|-]div(K * grad(f))
  r <- regexpr("div\\(<(.*)>\\*grad\\(<(.*)>\\)\\)", operator$expr)
  if (r[1] == -1) { ## no match found
    return(NULL)
  }
  match <- regmatches(operator$expr, r)
  return(operator$symbol_table[[substr(match, 5, 16)]])
}

tran_coeff_of <- function(operator) {
  if (!is.symbolic_expression(operator)) stop("Object ", deparse(substitute(operator)), " is not a valid operator.")
  r <- regexpr("inner\\(<(.*)>,grad\\(<(.*)>\\)\\)", operator$expr)
  if (r[1] == -1) { ## no match found
    return(NULL)
  }
  match <- regmatches(operator$expr, r)
  return(operator$symbol_table[[substr(match, 7, 18)]])
}

reac_coeff_of <- function(operator) {
  if (!is.symbolic_expression(operator)) stop("Object ", deparse(substitute(operator)), " is not a valid operator.")
  expr <- operator$expr
  ## clean operator from non-reactive terms
  expr <- gsub("div\\(<(.*)>\\*grad\\(<(.*)>\\)\\)|inner\\(<(.*)>,grad\\(<(.*)>\\)\\)|[+|-]| ", "", expr)
  if (regexpr("<(.*)>", expr)[1] == -1) { ## no match found
    return(NULL)
  }
  r <- regexpr("<(.*)>\\*<(.*)>", expr) ## matches reactive forms of type c*f
  if (r[1] == -1) {
    return(1.0) ## no c supplied, set c = 1
  } else {
    return(operator$symbol_table[[substr(expr, 1, 12)]])
  }
}

## Partial Differential Equation
.PartialDifferentialEquation <- R6::R6Class("PartialDifferentialEquation",
  private = list(
    pde_ = "ANY",
    operator_ = "SymbolicExpression",
    force_ = NULL
  ),
  public = list(
    initialize = function(operator = NULL, force = NULL) {
      private$operator_ <- operator
      private$force_ <- force

      ## recover the function on which the operator is applied, togheter with its functional space
      f <- NULL
      for (sym in operator$symbol_table) {
        if (is.symbolic_function(sym)) f <- sym
      }
      functional_space <- get_private(f)$functional_space_
      order <- functional_space$order
      quadrature_nodes <- as.matrix(get_private(functional_space)$basis_$quadrature_nodes())
      ## recover domain informations
      mesh <- functional_space$mesh
      local_dim <- mesh$local_dim
      embed_dim <- mesh$embed_dim

      ## infer PDE type
      pde_type_ <- classify(operator)
      ## parse operator expression and recover PDE parameters
      pde_parameters <- list()
      K <- diff_coeff_of(operator)
      b <- tran_coeff_of(operator)
      c <- reac_coeff_of(operator)
      space_varying <- FALSE
      if (is.function(K) || is.function(b) || is.function(c)) {
        space_varying <- TRUE
        zero <- function(nodes, ncol = 1) matrix(rep(0, times = nrow(nodes) * ncol), ncol = ncol)
        ## if any of the parameters is a function, the problem is considered space-varying
        pde_parameters$diff <- if (is.null(K)) zero(quadrature_nodes, ncol = local_dim^2) else K(quadrature_nodes)
        pde_parameters$tran <- if (is.null(b)) zero(quadrature_nodes, ncol = local_dim) else b(quadrature_nodes)
        pde_parameters$reac <- if (is.null(c)) zero(quadrature_nodes, ncol = 1) else c(quadrature_nodes)
      } else { ## constant coefficients case
        if (is.null(K)) {
          ## the laplacian case never returns NULL    
          pde_parameters$diff <- matrix(rep(0, times = local_dim^2), nrow = local_dim)
        } else {
          if (pde_type_ == pde_type$laplacian) {
            pde_parameters$diff <- as.numeric(K)
          } else {
            pde_parameters$diff <- if (is(K, "numeric")) { ## case -div(\mu*grad(f)), with \mu \in \mathbb{R}
              K * matrix(c(1, 0, 0, 1), nrow = local_dim, ncol = local_dim, byrow = T) 
            } else { ## general diffusion tensor
              as.matrix(K)
            }
          }
        }
        pde_parameters$tran <- if (is.null(b)) {
          matrix(rep(0, times = local_dim), nrow = local_dim, ncol = 1)
        } else {
          matrix(b, nrow = local_dim, ncol = 1)
        }
        pde_parameters$reac <- if (is.null(c)) 0.0 else as.numeric(c)
      }
      ## evaluate forcing function
      u <- as.matrix(force(quadrature_nodes))
      ## create cpp backend
      private$pde_ <- new(
        eval(parse(text = paste("cpp", "pde", local_dim, embed_dim, order, sep = "_"))),
        get_private(mesh)$mesh_,
        cpp_aligned_index(pde_type_),
        pde_parameters
      )
      private$pde_$set_forcing(u)
    },
    init = function() private$pde_$init()
  ),
  active = list(
    dofs_coordinates = function() private$pde_$dofs_coordinates(),
    mass  = function() private$pde_$mass (),
    stiff = function() private$pde_$stiff(),
    force = function() private$pde_$force(),
    operator = function() get_private(private$operator_)$symbolic_
  )
)

#' @export
PDE <- function(operator, forcing = NULL) {
  if (is.null(forcing)) {
    forcing <- function(nodes) matrix(rep(0, times = nrow(nodes)), ncol = 1)
  }
  pde <- .PartialDifferentialEquation$new(
    operator = operator,
    force = forcing
  )
  pde$init() ## initialize
  return(pde)
}
