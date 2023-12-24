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

pde_type_list <- list("laplacian" = 1, "elliptic" = 2, "parabolic" = 3)

## utility R wrapper around C++ PDE concept
.PdeCtr <- setRefClass(
    Class = "PdeObject",
    fields = c(
        domain = "MeshObject",  ## vector of time instants (parabolic problem only)
        pde_type = "integer",   ## pde type tag (one in pde_type_list)
        cpp_pde_handler = "ANY"
    ),
    methods = c(
        dofs_coordinates = function() { pde$get_dofs_coordinates() },
        set_dirichlet_bc = function(dirichlet_bc) {
            if (!is_space_time(domain)) {
                dirichlet_bc_ <- as.matrix(dirichlet_bc(pde$get_dofs_coordinates()))
            } else {
                dirichlet_bc_ <- as.matrix(dirichlet_bc(pde$get_dofs_coordinates(), domain$time_mesh))
            }
            pde$set_dirichlet_bc(dirichlet_bc_)
        },
        set_initial_condition = function(initial_condtion) {
            if (pde_type != pde_type_list$parabolic || !is_space_time(domain)) {
                stop("set_initial_condition is for space-time penalties only.")
            }
            pde$set_initial_condition(initial_condition(pde$get_dofs_coordinates()))
        },
        mass  = function() { pde$mass()  }, ## discretized mass matrix
        stiff = function() { pde$stiff() }, ## discretization of differential operator
        force = function() { pde$force() }, ## discretized forcing
        get_pde = function() { pde }
    )
)

## infers the type of a pde
extract_pde_type <- function(L) {
    if ("time" %in% names(L$params)) {
        return(pde_type_list$parabolic)
    }
    if ("diffusion" %in% names(L$params) && !is.matrix(L$params$diffusion)) {
        return(pde_type_list$laplacian)
    }
    return(pde_type_list$elliptic)
}

make_pde <- function(L, u, dirichlet_bc = NULL, initial_condition = NULL) {
    pde_type = extract_pde_type(L)

    time_domain <- vector()
    ## prepare pde parameters structure
    pde_parameters <- NULL
    pde_parameters$diffusion <- 1.0
    pde_parameters$transport <- matrix(0, nrow = 2, ncol = 1)
    pde_parameters$reaction  <- 0.0

    ## specific function bindings for pde types
    parse_parameters = list(
        laplacian = function() {
            pde_parameters$diffusion <- 1.0
        },
        elliptic  = function() {
            ## recover pde parameters from operator
            for (i in seq_len(L$tokens)) {
                pde_parameters[[paste(L$tokens[i])]] <- L$params[[paste(L$tokens[i])]]
            }
        },
        parabolic = function() {
            time_domain <- L$f$FunctionSpace$mesh$times ## time domain
            ## recover pde parameters from operator (if any)
            for (i in seq_len(L$tokens)) {
                pde_parameters[[paste(L$tokens[i])]] <- L$params[[paste(L$tokens[i])]]
            }
        }
        
    )
    parse_parameters[[pde_type]]()
    
    ## define Rcpp module
    D <- L$f$FunctionSpace$mesh$cpp_handler ## domain
    fe_order <- L$f$FunctionSpace$fe_order  ## finite element order
    pde_ <- NULL
    if (fe_order == 1) { ## linear finite elements
        pde_ <- new(cpp_pde_2d_fe1, D, pde_type - 1, pde_parameters)
    }
    if (fe_order == 2) { ## quadratic finite elements
        pde_ <- new(cpp_pde_2d_fe2, D, pde_type - 1, pde_parameters)
    }

    ## initialize and return
    quad_nodes <- as.matrix(pde_$get_quadrature_nodes())
    if(pde_type == pde_type_list$parabolic) {
        pde_$set_forcing(as.matrix(u(quad_nodes, times)))
        if (!is.null(initial_condition)) {
            pde_$set_initial_condition(initial_condition(pde_$get_dofs_coordinates()))
        }
    } else {
        pde_$set_forcing(as.matrix(u(quad_nodes)))
    }
    if (!is.null(dirichlet_bc)) {
        pde_$set_dirichlet_bc(dirichlet_bc(pde_$get_dofs_coordinates()))
    }
    pde_$init()
    return(pde_)
}

setGeneric("PDE", function(L, u, dirichlet_bc, initial_condition) standardGeneric("PDE"))
setMethod("PDE", ## space-only pde with boundary conditions
    signature = c(
        L = "DiffOpObject", u = "ANY", dirichlet_bc = "ANY",
        initial_condition = "missing"
    ),
    function(L, u, dirichlet_bc) {
        make_pde(L, u, dirichlet_bc, initial_condition = NULL)
    }
)
setMethod("PDE", ## parabolic pde
    signature = c(
        L = "DiffOpObject", u = "ANY", dirichlet_bc = "ANY",
        initial_condition = "ANY"
    ),
    function(L, u, dirichlet_bc, initial_condition) {
        make_pde(L, u, dirichlet_bc, initial_condition)
    }
)
setMethod("PDE", ## space-only pde with implicitly set homogeneous boundary conditions
    signature = c(
        L = "DiffOpObject", u = "ANY", dirichlet_bc = "missing",
        initial_condition = "missing"
    ),
    function(L, u, dirichlet_bc, initial_condition) {
        make_pde(L, u, dirichlet_bc = NULL, initial_condition = NULL)
    }
)
