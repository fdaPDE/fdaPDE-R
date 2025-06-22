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

#' A triangulated spatial domain
#'
.triangluation <- R6::R6Class(
  "triangulation_2_2",
  private = list(
    mesh_ = "ANY", ## cpp backend
    local_dim_ = 0L,
    embed_dim_ = 0L
  ),
  public = list(
    initialize = function(mesh, local_dim, embed_dim) {
      private$mesh_ <- mesh
      private$local_dim_ <- local_dim
      private$embed_dim_ <- embed_dim
    },
    locate = function(locations) {
      fdapde_assert(dim(locations)[1] > 0 && dim(locations)[2] == private$embed_dim_, "wrong matrix dimensions.")
      return(r_aligned_index(private$mesh_$locate(as.matrix(locations))))
    },
    sample = function(n_samples, seed = NULL) {
      fdapde_assert(n_samples > 0, "number of samples must be positive.")
      if (is.null(seed)) {
        seed <- -1
      } ## triggers random seed at cpp layer
      return(private$mesh_$sample(n_samples, seed))
    },
    ## mark
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
    clear_boundary_markers = function() {
      private$mesh_$clear_boundary_markers()
    },
    clear_cells_markers = function() {
      private$mesh_$clear_cells_markers()
    },
    ## filtering operations
    filter_cells_by_marker = function(marker) {
      fdapde_assert(marker >= 0, "invalid marker value.")
      return(r_aligned_index(private$mesh_$filter_cells_by_marker(marker)))
    },
    filter_boundary_by_marker = function(marker) {
      fdapde_assert(marker >= 0, "invalid marker value.")
      return(r_aligned_index(private$mesh_$filter_boundary_by_marker(marker)))
    },
    cell = function(cell_id) {
      fdapde_assert(cell_id > 0, "invalid cell index.")
      return(.cell$new(private$mesh_, cpp_aligned_index(cell_id)))
    }
  ),
  active = list(
    nodes = function() private$mesh_$nodes(),
    cells = function() r_aligned_index(private$mesh_$cells()),
    edges = function() r_aligned_index(private$mesh_$edges()),
    neighbors = function() {
      neigh_ <- r_aligned_index(private$mesh_$neighbors())
      neigh_[which(neigh_ == 0)] <- NA ## signal missing neighbor with NULL
      return(neigh_)
    },
    ## boundary
    boundary_nodes = function() private$mesh_$boundary_nodes(),
    boundary_edges = function() private$mesh_$boundary_edges(),
    ## sizes
    n_nodes = function() private$mesh_$n_nodes(),
    n_cells = function() private$mesh_$n_cells(),
    n_edges = function() private$mesh_$n_edges(),
    n_boundary_nodes = function() private$mesh_$n_boundary_nodes(),
    n_boundary_edges = function() private$mesh_$n_boundary_edges(),
    ## dimensions
    local_dim = function() private$local_dim_,
    embed_dim = function() private$embed_dim_,
    ## utilities
    bbox = function() private$mesh_$bbox(),
    area = function() private$mesh_$measure(),
    cells_markers = function() private$mesh_$cells_markers(),
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


#' @export
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

## the total number of intervals will be nx, and the overall number of nodes is nx + 1
#' @export
## triangluationInterval <- function(a, b, n = NULL, by = NULL) {
##   if (!is.null(n) && !is.null(by)) stop("too many arguments.")
##   mesh_data <- list()
##   by_x <- if (!is.null(n)) ((b - a) / (n - 1)) else by
##   mesh_data$nodes <- if (!is.null(by_x)) {
##     as.matrix(seq(from = a, to = b, by = by_x))
##   } else {
##     as.matrix(seq(from = a, to = b))
##   }
##   return(.triangluation$new(
##     mesh = new(cpp_mesh_1_1, mesh_data),
##     local_dim = 1,
##     embed_dim = 1
##   ))
## }

## #' @export
## triangluationUnitInterval <- function(n = NULL, by = NULL) {
##   return(triangluationInterval(0, 1, n))
## }

## #' @export
## triangluationRectangle <- function(a_x, b_x, a_y, b_y, nx = NULL, ny = NULL, by_x = NULL, by_y = NULL) {
##   if ((!is.null(nx) && !is.null(by_x)) || (!is.null(ny) && !is.null(by_y))) stop("too many arguments.")
##   mesh_data <- list()
##   by_x <- if (!is.null(nx)) ((b_x - a_x) / (nx - 1)) else by_x
##   by_y <- if (!is.null(ny)) ((b_y - a_y) / (ny - 1)) else by_y

##   grid_x <- if (!is.null(by_x)) as.matrix(seq(a_x, b_x, by = by_x)) else as.matrix(seq(a_x, b_X))
##   grid_y <- if (!is.null(by_y)) as.matrix(seq(a_y, b_y, by = by_y)) else as.matrix(seq(a_y, b_y))
##   mesh_data$nodes <- as.matrix(expand.grid(grid_x, grid_y))
##   ## build triangles (each subrectangle is split in 2 triangles)
##   triangles <- matrix(0, nrow = 2 * (nx - 1) * (ny - 1), ncol = 3)
##   j <- 1
##   for (y in seq_len(ny - 1)) {
##     for (x in seq_len(nx - 1)) {
##       ## build vector of vertices of j-th subrectangle
##       p <- x + (y - 1) * nx ## base point
##       v <- c(p, p + 1, p + nx, p + nx + 1)
##       ## compute vertices of each triangle in the subrectangle
##       triangles[j,     ] <- c(v[1], v[2], v[3])
##       triangles[j + 1, ] <- c(v[2], v[3], v[4])
##       j <- j + 2
##     }
##   }
##   mesh_data$elements <- cpp_aligned_index(triangles)
##   ## build boundary
##   boundary <- matrix(0, nrow(mesh_data$nodes))
##   for (i in 1:nrow(mesh_data$nodes)) {
##     if ((mesh_data$nodes[i, 1] == a_x || mesh_data$nodes[i, 1] == b_x) ||
##       (mesh_data$nodes[i, 2] == a_y || mesh_data$nodes[i, 2] == b_y)) {
##       boundary[i] <- 1
##     }
##   }
##   mesh_data$boundary <- boundary
##   return(.triangluation$new(
##     mesh = new(cpp_mesh_2_2, mesh_data),
##     local_dim = 2,
##     embed_dim = 2
##   ))
## }

## #' @export
## triangluationSquare <- function(a, b, n = NULL, by = NULL) {
##   return(triangluationRectangle(a, b, a, b, n, n, by, by))
## }

## #' @export
## triangluationUnitSquare <- function(n = NULL, by = NULL) {
##   return(triangluationSquare(0, 1, n, by))
## }

## #' @export
## triangluationCube <- function(a, b, n = NULL, by = NULL) {
##   if (!is.null(n) && !is.null(by)) stop("too many arguments.")
##   mesh_data <- list()
##   if (!is.null(n)) {
##     by_x <- ((b - a) / (n - 1))
##   } else {
##     by_x <- by
##     n <- ((b - a) / by_x) + 1
##   }
##   grid <- if (!is.null(by_x)) as.matrix(seq(a, b, by = by_x)) else as.matrix(seq(a, b))
##   mesh_data$nodes <- as.matrix(expand.grid(grid, grid, grid))
##   ## build tetrahedrons (each subcube can be split in 5 tetrahedrons)
##   tetrahedrons <- matrix(0, nrow = 6 * (n - 1)^3, ncol = 4)
##   j <- 1
##   for (z in seq_len(n - 1)) {
##     for (y in seq_len(n - 1)) {
##       for (x in seq_len(n - 1)) {
##         ## build vector of vertices of i-th subcube
##         p <- x + (y - 1) * n + (z - 1) * n^2 ## base point
##         v <- c(p, p + 1, p + n, p + n + 1, p + n^2, p + n^2 + 1, p + n^2 + n, p + n^2 + n + 1)
##         ## compute vertices of each thetraedron in the subcube
##         tetrahedrons[j,     ] <- c(v[7], v[3], v[2], v[1])
##         tetrahedrons[j + 1, ] <- c(v[7], v[5], v[2], v[1])
##         tetrahedrons[j + 2, ] <- c(v[7], v[6], v[8], v[2])
##         tetrahedrons[j + 3, ] <- c(v[7], v[4], v[8], v[2])
##         tetrahedrons[j + 4, ] <- c(v[7], v[4], v[3], v[2])
##         tetrahedrons[j + 5, ] <- c(v[7], v[6], v[5], v[2])
##         j <- j + 6
##       }
##     }
##   }
##   mesh_data$elements <- cpp_aligned_index(tetrahedrons)
##   ## build boundary
##   boundary <- matrix(0, nrow(mesh_data$nodes))
##   for (i in 1:nrow(mesh_data$nodes)) {
##     if ((mesh_data$nodes[i, 1] == a || mesh_data$nodes[i, 1] == b) ||
##       (mesh_data$nodes[i, 2] == a || mesh_data$nodes[i, 2] == b) ||
##       (mesh_data$nodes[i, 3] == a || mesh_data$nodes[i, 3] == b)) {
##       boundary[i] <- 1
##     }
##   }
##   mesh_data$boundary <- boundary
##   ##return(.triangluation$new(
##   ##  mesh = new(cpp_mesh_3_3, mesh_data),
##   ##  local_dim = 3,
##   ##  embed_dim = 3
##   ##))
##   return(mesh_data)
## }

## #' @export
## triangluationUnitCube <- function(n = NULL, by = NULL) {
##   return(triangluationCube(0, 1, n))
## }
