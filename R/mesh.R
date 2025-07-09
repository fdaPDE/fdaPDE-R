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

#' An R6 class representing a triangulated spatial domain.
#'
#' @rdname triangulation
.triangluation <- R6::R6Class(
  "triangulation_2_2",
  private = list(
    mesh_ = "ANY", ## cpp backend
    local_dim_ = 0L,
    embed_dim_ = 0L
  ),
  public = list(
    #' @description
    #' Creates a new triangulation object.
    #'
    #' @param mesh A C++ backend object created internally by [triangulation()].
    #' @param local_dim The local dimension of the space.
    #' @param embed_dim The dimension of the space in which the domain is embedded.
    initialize = function(mesh, local_dim, embed_dim) {
      private$mesh_ <- mesh
      private$local_dim_ <- local_dim
      private$embed_dim_ <- embed_dim
    },
    #' @description
    #' A utility to locate points over the domain.  
    #' For each point, the function returns the index of the cell that contains it, or \code{-1} if the point lies  
    #' outside the spatial domain.
    #'
    #' @param locations (`matrix`) A matrix of locations (each row is a point).
    #'
    #' @return An \code{nlocs}-by-1 matrix containing the index of the containing cell for each point, or \code{-1} if the point is outside the domain.
    locate = function(locations) {
      fdapde_assert(dim(locations)[1] > 0 && dim(locations)[2] == private$embed_dim_, "wrong matrix dimensions.")
      return(r_aligned_index(private$mesh_$locate(as.matrix(locations))))
    },
    #' @description
    #' Generates a uniform sample over the domain.
    #'
    #' @param n_samples (`integer(1)`) The number of points to sample.
    #' @param seed (`integer(1)`) A seed for the random number generator. Defaults to \code{NULL}.
    sample = function(n_samples, seed = NULL) {
      fdapde_assert(n_samples > 0, "number of samples must be positive.")
      if (is.null(seed)) {
        seed <- -1
      } ## triggers random seed at cpp layer
      return(private$mesh_$sample(n_samples, seed))
    },
    #' @description
    #' Marks a subdomain by marking all cells that satisfy a geometric condition.
    #'
    #' @param marker (`integer(1)`) An integer denoting the subdomain.
    #' @param predicate A function implementing the geometric condition that cells must satisfy to be marked as belonging to the subdomain \code{marker}.
    mark_cells = function(marker, predicate) {
      fdapde_assert(marker >= 0, "invalid marker value.")
      ## evaluate predicate at cells
      mask = matrix(rep(0, times = self$n_cells), ncol = 1)
      for (i in cpp_aligned_index(seq_len(self$n_cells))) {
        if (predicate(.cell$new(private$mesh_, i))) {
          mask[i + 1] = marker
        }
      }
      private$mesh_$mark_cells(marker, mask)
    },
    #' @description
    #' Marks a sub-boundary by marking all boundary nodes that satisfy a geometric condition.
    #'
    #' @param marker (`integer(1)`) An integer denoting the sub-boundary.
    #' @param predicate A function implementing the geometric condition that boundary nodes must satisfy to be marked as belonging to the sub-boundary \code{marker}.
    mark_boundary = function(marker, predicate) {
      fdapde_assert(marker >= 0, "invalid marker value.")
      ## evaluate predicate at cells
      mask = matrix(rep(0, times = self$n_edges), ncol = 1)
      for (i in cpp_aligned_index(seq_len(self$n_edges))) {
        if (predicate(.edge$new(private$mesh_, i))) {
          mask[i + 1] = marker
        }
      }
      private$mesh_$mark_boundary(marker, mask)
    },
    #' @description
    #' Clears the boundary markers previously set by \code{mark_boundary}.
    clear_boundary_markers = function() {
      private$mesh_$clear_boundary_markers()
    },
    #' @description
    #' Clears the cell markers previously set by \code{mark_cells}.
    clear_cells_markers = function() {
      private$mesh_$clear_cells_markers()
    },
    #' @description
    #' Filters the cells according to the provided marker.
    #'
    #' @param marker (`integer(1)`) An integer denoting the subdomain.
    #' @return ....
    filter_cells_by_marker = function(marker) {
      fdapde_assert(marker >= 0, "invalid marker value.")
      return(r_aligned_index(private$mesh_$filter_cells_by_marker(marker)))
    },
    #' @description
    #' Filters the boundary nodes according to the provided marker.
    #' 
    #' @param marker (`integer(1)`) An integer denoting the sub-boundary.
    #' @return ....
    filter_boundary_by_marker = function(marker) {
      fdapde_assert(marker >= 0, "invalid marker value.")
      return(r_aligned_index(private$mesh_$filter_boundary_by_marker(marker)))
    },
    #' @description
    #' Returns a cell.
    #'
    #' @param cell_id (`integer(1)`) An integer denoting the cell ID.
    #' @return An R6 object representing the \code{cell}.
    cell = function(cell_id) {
      fdapde_assert(cell_id > 0, "invalid cell index.")
      return(.cell$new(private$mesh_, cpp_aligned_index(cell_id)))
    }
  ),
  active = list(
    #' @field nodes (`matrix`)\cr
    #' Return the matrix of nodes
    nodes = function() private$mesh_$nodes(),
    #' @field cells (`matrix`)\cr
    #' Return the matrix of cells
    cells = function() r_aligned_index(private$mesh_$cells()),
    #' @field edges (`matrix`)\cr
    #' Return the matrix of edges
    edges = function() r_aligned_index(private$mesh_$edges()),
    #' @field neighbors (`matrix`)\cr
    #' Return the matrix of neighbors
    neighbors = function() {
      neigh_ <- r_aligned_index(private$mesh_$neighbors())
      neigh_[which(neigh_ == 0)] <- NA ## signal missing neighbor with NULL
      return(neigh_)
    },
    ## boundary
    #' @field boundary_nodes (`matrix`)\cr
    #' Return the matrix of boundary nodes
    boundary_nodes = function() private$mesh_$boundary_nodes(),
    #' @field boundary_edges (`matrix`)\cr
    #' Return the matrix of boundary_edges
    boundary_edges = function() private$mesh_$boundary_edges(),
    ## sizes
    #' @field n_nodes (`integer(1)`)\cr
    #' Return the number of nodes
    n_nodes = function() private$mesh_$n_nodes(),
    #' @field n_cells (`integer(1)`)\cr
    #' Return the number of cells
    n_cells = function() private$mesh_$n_cells(),
    #' @field n_edges (`integer(1)`)\cr
    #' Return the number of edges
    n_edges = function() private$mesh_$n_edges(),
    #' @field n_boundary_nodes (`integer(1)`)\cr
    #' Return the number of boundary nodes
    n_boundary_nodes = function() private$mesh_$n_boundary_nodes(),
    #' @field n_boundary_edges (`integer(1)`)\cr
    #' Return the number of boundary edges
    n_boundary_edges = function() private$mesh_$n_boundary_edges(),
    ## dimensions
    #' @field local_dim (`integer(1)`)\cr
    #' Return the local dimension of the domain 
    local_dim = function() private$local_dim_,
    #' @field embed_dim (`integer(1)`)\cr
    #' Return the dimension of the space where the domain is embedded 
    embed_dim = function() private$embed_dim_,
    ## utilities
    #' @field bbox (`integer(1)`)\cr
    #' Return the bounding box
    bbox = function() private$mesh_$bbox(),
    #' @field area (`integer(1)`)\cr
    #' Return the measure of the domain
    area = function() private$mesh_$measure(),
    #' @field cells_markers (`matrix`)\cr
    #' Return the cells markers
    cells_markers = function() private$mesh_$cells_markers(),
    #' @field edges_markers (`matrix`)\cr
    #' Return the edges markers
    edges_markers = function() private$mesh_$edges_markers()
  )
)

.cell <- R6::R6Class(
  "cell",
  private = list(
    mesh_ = "ANY", ## cpp backend
    id_ = integer()
  ),
  public = list(
    initialize = function(mesh, id) {
      private$mesh_ <- mesh
      private$id_ <- id
    }
  ),
  active = list(
    coords = function() private$mesh_$cell_coords(private$id_),
    measure = function() private$mesh_$cell_measure(private$id_),
    bbox = function() private$mesh_$cell_bbox(private$id_),
    barycenter = function() private$mesh_$cell_barycenter(private$id_),
    circumcenter = function() private$mesh_$cell_circumcenter(private$id_),
    diameter = function() private$mesh_$cell_diameter(private$id_)
  )
)

.edge <- R6::R6Class(
  "edge",
  private = list(
    mesh_ = "ANY", ## cpp backend
    id_ = integer()
  ),
  public = list(
    initialize = function(mesh, id) {
      private$mesh_ <- mesh
      private$id_ <- id
    }
  ),
  active = list(
    coords = function() private$mesh_$edge_coords(private$id_)
  )
)


#' Create a mesh object
#'
#' @param nodes A \code{#nodes}-by-2 matrix containing the coordinates of the mesh nodes.
#' @param cells A \code{#cells}-by-3 matrix specifying the cells by giving the row indices in \code{nodes} of the cells' vertices.
#' @param boundary A \code{#nodes}-by-1 matrix with entries either \code{1} or \code{0}. An entry of \code{1} indicates that the corresponding node is a boundary node; \code{0} indicates it is not.
#' 
#' @return An R6 object representing a triangulation.
#' 
#' @rdname triangulation
#' @export
#' @examples
#' \dontrun{
#' library(RTriangle)
#' library(fdaPDE2)
#' p <- pslg(P=rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)),
#' S=rbind(c(1, 2), c(2, 3), c(3, 4), c(4,1)))
#' mesh_data <- triangulate(p, a = 0.00125, q=30)
#' mesh <- triangulation(nodes = mesh_data$P, cells = mesh_data$T, boundary = mesh_data$PB)
#' }
triangulation <- function(nodes, cells, boundary) {
  ## check input dimensions
  fdapde_assert(
    nrow(nodes) > 0 &&
      nrow(cells) > 0 &&
      nrow(boundary) == nrow(nodes) &&
      ncol(nodes) %in% c(1, 2, 3) &&
      ncol(cells) %in% c(2, 3, 4) &&
      ncol(boundary) == 1,
    "invalid matrix dimensions"
  )

  local_dim <- ncol(cells) - 1
  embed_dim <- ncol(nodes)
  ## instantiate cpp backend
  data = list(
    nodes = as.matrix(nodes),
    cells = as.matrix(cells),
    boundary = as.matrix(boundary)
  )
  if (local_dim == 2 && embed_dim == 2) cpp_backend <- new(cpp_triangulation_2_2, data)

  ## construct mesh and return
  return(.triangluation$new(
    mesh = cpp_backend,
    local_dim = local_dim,
    embed_dim = embed_dim
  ))
}

#' @export
plot.triangulation_2_2 <- function(x, boundary_markers = NULL, boundary_palette = NULL, ...) {
  ## start plotting nodes
  plot(x$nodes, xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", pch = 19, cex = 0.5, ...)
  segments(
    x$nodes[x$edges[, 1], 1],
    x$nodes[x$edges[, 1], 2],
    x$nodes[x$edges[, 2], 1],
    x$nodes[x$edges[, 2], 2],
    ...
  )

  if (!is.null(boundary_markers)) {
    i <- 1
    sorted_boundary_markers <- sort(boundary_markers) ## sort, to give higher precedence to large labels
    for (marker in sorted_boundary_markers) {
      e = unit_square$filter_boundary_by_marker(marker)
      segments(
        x0 = x$nodes[x$edges[e, 1], 1],
        y0 = x$nodes[x$edges[e, 1], 2],
        x1 = x$nodes[x$edges[e, 2], 1],
        y1 = x$nodes[x$edges[e, 2], 2],
        col = boundary_palette[which(marker == boundary_markers)],
        lwd = 5
      )
      i <- i + 1
    }
  }
}
