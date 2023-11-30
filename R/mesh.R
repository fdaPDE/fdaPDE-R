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

.MeshCtr <- setRefClass(
    Class = "MeshObject",
    fields = c(
        cpp_handler = "ANY",       ## cpp backend
        m = "integer",             ## local dimension
        n = "integer",             ## embedding dimension
        time_mesh = "vector",      ## for spatio-temporal domains, the time dimension
        domain_type = "character",
        nodes = "matrix"           ## node of the triangulation
    )
)

## checks if a MeshObject represnts a spatio-temporal domain
is_space_time <- function(mesh) {
    return(length(mesh$time_mesh) != 0)
}

#' Create mesh object
#'
#' @param domain could be a \code{triangulation} returned by \code{\link[RTriangle]{triangulate}} or a named list containing:
#' \itemize{
#'    \item{\code{nodes}, a #nodes-by-2 matrix containing the x and y coordinates of the mesh nodes;}
#'    \item{\code{elements}, a #elements-by-3 matrix specifiying the triangles giving the row's indices in \code{nodes} of the triangles' vertices;}
#'    \item{\code{boundary}, a #nodes-by-1 matrix, with entries either '1' or '0'. An entry '1' indicates that the corresponding node is a boundary node; 
#'           an entry '0' indicates that the corresponding node is not a boundary node.}
#' }
#' 
#' @return An S4 object representing a Mesh.
#' @rdname MeshObject
#' @export
#' @examples
#' \dontrun{
#' library(RTriangle)
#' library(femR)
#' p <- pslg(P=rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)),
#' S=rbind(c(1, 2), c(2, 3), c(3, 4), c(4,1)))
#' unit_square <- triangulate(p, a = 0.00125, q=30)
#' mesh <- Mesh(unit_square)
#' }
setGeneric("Mesh", function(domain) standardGeneric("Mesh"))

#' @rdname MeshObject
setMethod("Mesh",
    signature = c(domain = "list"),
    function(domain) {
        domain$elements <- domain$elements - 1 ## perform index realignment for cpp handler
        storage.mode(domain$elements) <- "integer"
        ## extract local and embedding dimensions
        m <- ncol(domain$elements) - 1
        n <- ncol(domain$nodes)
        ## derive domain type (TODO: handle network domains)
        if (m == 2 && n == 2) {
            name_flag <- "2d"
        } else if (m == 2 && n == 3) {
            name_flag <- "surface"
        } else if (m == 3 && n == 3) {
            name_flag <- "3d"
        } else {
            stop("wrong input argument provided.")
        }
        ## construct mesh and return
        .MeshCtr(
            cpp_handler = new(eval(parse(text = paste("cpp_", name_flag, "_domain", sep = ""))), domain),
            m = as.integer(m),
            n = as.integer(n),
            time_mesh = vector(mode = "double"),
            domain_type = name_flag,
            nodes = domain$nodes
        )
    }
)

# @importClassesFrom RTriangle triangulation

#' @rdname MeshObject
setMethod("Mesh",
    signature = c(domain = "triangulation"),
    function(domain) {
        elements <- domain$T - 1 ## perform index realignment for cpp handler
        ## extract node and boundary informations
        nodes <- domain$P
        boundary <- matrix(0, nrow = nrow(nodes), ncol = 1)
        boundary[as.vector(domain$E[domain$EB == 1, ]), ] <- 1
        storage.mode(elements) <- "integer"
        storage.mode(nodes)    <- "numeric"
        storage.mode(boundary) <- "integer"
        ## extract local and embedding dimensions
        m <- ncol(elements) - 1
        n <- ncol(nodes)
        ## prepare list and construct mesh
        domain <- list(elements = elements, nodes = nodes, boundary = boundary)
        if (m == 2 & n == 2) {
            .MeshCtr(
                cpp_handler = new(cpp_2d_domain, domain),
                m = as.integer(m),
                n = as.integer(n),
                time_mesh = vector(mode = "double"),
                domain_type = "2d",
                nodes = domain$P
            )
        } else {
            stop("wrong input argument provided.")
        }
    }
)

#' create spatio-temporal domain
#'
#' @param op1 A mesh object created by \code{Mesh}.
#' @param op2 A numeric vector.
#' @return An S4 object representing a spatio-temporal domain.
#' @rdname MeshObject_times_vector
#' @export 
setGeneric("%X%", function(op1, op2) standardGeneric("%X%"))

#' @rdname MeshObject_times_vector
setMethod("%X%",
    signature = c(op1 = "MeshObject", op2 = "numeric"),
    function(op1, op2) {
        if (op2[1] > op2[length(op2)]) {
              stop("first time instant greater than last time instant.")
        }
        op1$time_mesh <- op2
        op1
    }
)

## Mesh - auxiliary methods
unroll_edges_aux <- function(mesh_) {
    mesh <- mesh_$cpp_handler
    edges <- matrix(nrow = 3 * nrow(mesh$elements()), ncol = 2)
    for (i in 1:nrow(mesh$elements())) {
        edges[(3 * (i - 1) + 1), ] <- mesh$elements()[i, c(1, 2)] + 1
        edges[(3 * (i - 1) + 2), ] <- mesh$elements()[i, c(2, 3)] + 1
        edges[(3 * (i - 1) + 3), ] <- mesh$elements()[i, c(3, 1)] + 1
    }
    edges
}

setGeneric("unroll_edges", function(Mesh) standardGeneric("unroll_edges"))
setMethod("unroll_edges", "MeshObject", function(Mesh) {
    unroll_edges_aux(Mesh)
})

plot_mesh_aux <- function(mesh_, ...) {
    mesh <- mesh_$cpp_handler
    edges <- unroll_edges(mesh)
    plot_ly(...) %>%
        add_markers(
            x = mesh$nodes()[, 1],
            y = mesh$nodes()[, 2],
            color = I("black"), size = I(1),
            hoverinfo = "text",
            text = paste(
                "</br><b> Coordinates:", round(mesh$nodes()[, 1], 2),
                round(mesh$nodes()[, 2], 2)
            ),
            showlegend = T,
            visible = T
        ) %>%
        add_segments(
            x = mesh$nodes()[edges[, 1], 1],
            y = mesh$nodes()[edges[, 1], 2],
            xend = mesh$nodes()[edges[, 2], 1],
            yend = mesh$nodes()[edges[, 2], 2],
            color = I("black"), size = I(1),
            showlegend = F
        ) %>%
        layout(
            xaxis = list(
                title = "",
                showgrid = F,
                zeroline = F,
                showticklabels = F
            ),
            yaxis = list(
                title = "",
                showgrid = F,
                zeroline = F,
                showticklabels = F
            )
        )
}


# setMethod("plot", signature=c(x="MeshObject"), function(x, ...){
#   plot_mesh_aux(x, ...)  
# })

#' Plot a Mesh object
#'
#' @param x A \code{MeshObject} object defining the triangular mesh, as generated by \code{Mesh}
#' @param ... Arguments representing graphical options to be passed to \code{\link[plotly]{plot_ly}}.
#' @return A plotly object
#' @export
#' @examples
#' \dontrun{
#' library(femR)
#' data("unit_square")
#' mesh <- Mesh(unit_square)
#' plot(mesh)
#' }
plot.MeshObject <-function(x, ...){
  plot_mesh_aux(x, ...)
}
