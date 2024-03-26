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
.Mesh <- R6::R6Class(
  "Mesh",
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
      return(r_aligned_index(private$mesh_$locate(as.matrix(locations))))
    }
  ),
  active = list(
    nodes = function() private$mesh_$nodes(),
    elements  = function() r_aligned_index(private$mesh_$elements()),
    boundary  = function() private$mesh_$boundary(),
    neighbors = function() {
      neigh_ <- r_aligned_index(private$mesh_$neighbors())
      neigh_[neigh_ == 0] <- -1 ## signal missing neighbor with -1
      return(neigh_)
    },
    local_dim = function() private$local_dim_,
    embed_dim = function() private$embed_dim_
  )
)

#' @export
Mesh <- function(data) {
  data$elements <- cpp_aligned_index(data$elements)
  storage.mode(data$elements) <- "integer"
  ## extract local and embedding dimensions
  local_dim <- ncol(data$elements) - 1
  embed_dim <- ncol(data$nodes)
  ## derive domain type
  cpp_backend <- new(
    eval(parse(text = paste("cpp", "mesh", as.character(local_dim), as.character(embed_dim), sep = "_"))), data
  )
  ## construct mesh and return
  return(.Mesh$new(
    mesh = cpp_backend,
    local_dim = local_dim,
    embed_dim = embed_dim
  ))
}

## the total number of intervals will be nx, and the overall number of nodes is nx + 1
#' @export
MeshInterval <- function(a, b, n = NULL, by = NULL) {
  if (!is.null(n) && !is.null(by)) stop("too many arguments.")
  mesh_data <- list()
  by_x <- if (!is.null(n)) ((b - a) / (n - 1)) else by
  mesh_data$nodes <- if (!is.null(by_x)) {
    as.matrix(seq(from = a, to = b, by = by_x))
  } else {
    as.matrix(seq(from = a, to = b))
  }
  return(.Mesh$new(
    mesh = new(cpp_mesh_1_1, mesh_data),
    local_dim = 1,
    embed_dim = 1
  ))
}

#' @export
MeshUnitInterval <- function(n = NULL, by = NULL) {
  return(IntervalMesh(0, 1, n))
}

#' @export
MeshRectangle <- function(a_x, b_x, a_y, b_y, nx = NULL, ny = NULL, by_x = NULL, by_y = NULL) {
  if ((!is.null(nx) && !is.null(by_x)) || (!is.null(ny) && !is.null(by_y))) stop("too many arguments.")
  mesh_data <- list()
  by_x <- if (!is.null(nx)) ((b_x - a_x) / (nx - 1)) else by_x
  by_y <- if (!is.null(ny)) ((b_y - a_y) / (ny - 1)) else by_y

  grid_x <- if (!is.null(by_x)) as.matrix(seq(a_x, b_x, by = by_x)) else as.matrix(seq(a_x, b_X))
  grid_y <- if (!is.null(by_y)) as.matrix(seq(a_y, b_y, by = by_y)) else as.matrix(seq(a_y, b_y))
  mesh_data$nodes <- as.matrix(expand.grid(grid_x, grid_y))
  ## build triangles (each subrectangle is split in 2 triangles)
  triangles <- matrix(0, nrow = 2 * (nx - 1) * (ny - 1), ncol = 3)
  j <- 1
  for (y in seq_len(ny - 1)) {
    for (x in seq_len(nx - 1)) {
      ## build vector of vertices of j-th subrectangle
      p <- x + (y - 1) * nx ## base point
      v <- c(p, p + 1, p + nx, p + nx + 1)
      ## compute vertices of each triangle in the subrectangle
      triangles[j,     ] <- c(v[1], v[2], v[3])
      triangles[j + 1, ] <- c(v[2], v[3], v[4])
      j <- j + 2
    }
  }
  mesh_data$elements <- cpp_aligned_index(triangles)
  ## build boundary
  boundary <- matrix(0, nrow(mesh_data$nodes))
  for (i in 1:nrow(mesh_data$nodes)) {
    if ((mesh_data$nodes[i, 1] == a_x || mesh_data$nodes[i, 1] == b_x) ||
      (mesh_data$nodes[i, 2] == a_y || mesh_data$nodes[i, 2] == b_y)) {
      boundary[i] <- 1
    }
  }
  mesh_data$boundary <- boundary
  return(.Mesh$new(
    mesh = new(cpp_mesh_2_2, mesh_data),
    local_dim = 2,
    embed_dim = 2
  ))
}

#' @export
MeshSquare <- function(a, b, n = NULL, by = NULL) {
  return(MeshRectangle(a, b, a, b, n, n, by, by))
}

#' @export
MeshUnitSquare <- function(n = NULL, by = NULL) {
  return(MeshSquare(0, 1, n, by))
}

#' @export
MeshCube <- function(a, b, n = NULL, by = NULL) {
  if (!is.null(n) && !is.null(by)) stop("too many arguments.")
  mesh_data <- list()
  if (!is.null(n)) {
    by_x <- ((b - a) / (n - 1))
  } else {
    by_x <- by
    n <- ((b - a) / by_x) + 1
  }
  grid <- if (!is.null(by_x)) as.matrix(seq(a, b, by = by_x)) else as.matrix(seq(a, b))
  mesh_data$nodes <- as.matrix(expand.grid(grid, grid, grid))
  ## build tetrahedrons (each subcube can be split in 5 tetrahedrons)
  tetrahedrons <- matrix(0, nrow = 5 * (n - 1)^3, ncol = 4)
  j <- 1
  for (z in seq_len(n - 1)) {
    for (y in seq_len(n - 1)) {
      for (x in seq_len(n - 1)) {
        ## build vector of vertices of i-th subcube
        p <- x + (y - 1) * n + (z - 1) * n^2 ## base point
        v <- c(p, p + 1, p + n, p + n + 1, p + n^2, p + n^2 + 1, p + n^2 + n, p + n^2 + n + 1)
        ## compute vertices of each thetraedron in the subcube
        tetrahedrons[j,     ] <- c(v[1], v[2], v[3], v[5])
        tetrahedrons[j + 1, ] <- c(v[2], v[3], v[4], v[8])
        tetrahedrons[j + 2, ] <- c(v[2], v[3], v[5], v[8])
        tetrahedrons[j + 3, ] <- c(v[2], v[5], v[6], v[8])
        tetrahedrons[j + 4, ] <- c(v[3], v[5], v[7], v[8])
        j <- j + 5
      }
    }
  }
  mesh_data$elements <- cpp_aligned_index(tetrahedrons)
  ## build boundary
  boundary <- matrix(0, nrow(mesh_data$nodes))
  for (i in 1:nrow(mesh_data$nodes)) {
    if ((mesh_data$nodes[i, 1] == a || mesh_data$nodes[i, 1] == b) ||
      (mesh_data$nodes[i, 2] == a || mesh_data$nodes[i, 2] == b) ||
      (mesh_data$nodes[i, 3] == a || mesh_data$nodes[i, 3] == b)) {
      boundary[i] <- 1
    }
  }
  mesh_data$boundary <- boundary
  return(.Mesh$new(
    mesh = new(cpp_mesh_3_3, mesh_data),
    local_dim = 3,
    embed_dim = 3
  ))
}

#' @export
MeshUnitCube <- function(n = NULL, by = NULL) {
  return(MeshCube(0, 1, n))
}

.TensorizedMesh <- R6::R6Class(
  "TensorizedMesh",
  inherit = .Mesh,
  private = list(
    rhs_ = vector(mode = "double", length = 0L)
  ),
  public = list(
    initialize = function(lhs, rhs) {
      super$initialize(lhs, lhs$local_dim, lhs$embed_dim)
      private$rhs_ <- rhs
    }
  ),
  active = list(
    rhs_nodes = function() rhs_
  )
)

#' @export
`%X%.Mesh` <- function(lhs, rhs) {
  if (!is(rhs, "vector") || !inheris(rsh, "Mesh") ||
    (inherits(rhs, "Mesh") && !(rhs$local_dim == 1 && rhs$embed_dim == 1))) {
    stop(deparse(substitute(rhs)), " must be a vector or a 1D Mesh object.")
  }
  return(.TensorizedMesh$new(lhs, rhs))
}

`%<<%` <- function(x, ...) UseMethod("%<<%", x)
`%<<%.file` <- function(x, what) {
  if (!isOpen(x)) stop(deparse(substitute(x)), ": file not open.")
  cat(what, file = x, append = TRUE)
}

#' @export
write_vtu <- function(mesh, file = NULL, data = NULL) {
  ## check that mesh is embedded in 3D
  if (!inherits(mesh, "Mesh")) stop(deparse(substitute(mesh)), ": not a valid Mesh object.")
  if (!(mesh$embed_dim == 3)) stop(deparse(substitute(mesh)), ": not embedded in 3D space.")
  n_nodes <- nrow(mesh$nodes)
  n_elems <- nrow(mesh$elements)
  dim <- mesh$embed_dim

  ## open file connector in write mode
  if (!is.character(file)) stop("invalid output file.")
  vtu <- file(file, open = "w")

  ## insert vtk header
  vtu %<<% '<?xml version="1.0"?>\n'
  vtu %<<% '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n'
  vtu %<<% "<UnstructuredGrid>\n"
  vtu %<<% paste0('<Piece NumberOfPoints="', n_nodes, '" NumberOfCells="', n_elems, '">\n')
  ## mesh nodes (vtk points)
  vtu %<<% "<Points>\n"
  vtu %<<% paste0('<DataArray type="Float64" Name="nodes" NumberOfComponents="', dim, '" format="ascii">\n')
  write.table(mesh$nodes, file = vtu, append = T, row.names = F, col.names = F)
  vtu %<<% "</DataArray>\n"
  vtu %<<% "</Points>\n"
  ## mesh elements (vtk cells)
  vtu %<<% "<Cells>\n"
  ## specifies the point connectivity. All the cells’ point lists are concatenated together.
  vtu %<<% '<DataArray type="Int32" Name="connectivity" format="ascii">\n'
  write.table(
    format(mesh$elements - 1, scientific = F),
    file = vtu, append = T, row.names = F, col.names = F, quote = F
  )
  vtu %<<% "</DataArray>\n"
  ## specifies the offset into the connectivity array for the end of each cell
  vtu %<<% '<DataArray type="Int32" Name="offsets" format="ascii">\n'
  offset <- if (mesh$local_dim == 2) 3 else 4 ## number of points per element for triangular and tetrahedral meshes
  write.table(
    format(matrix(seq(offset, offset * n_elems, by = offset), nrow = 1, byrow = T), scientific = F),
    file = vtu, append = T, row.names = F, col.names = F, quote = F
  )
  vtu %<<% "</DataArray>\n"
  ## specifies the type of each cell
  vtu %<<% '<DataArray type="Int32" Name="types" format="ascii">\n'
  vtk_cell_type <- if (mesh$local_dim == 2) 5 else 10 ## vtk code for (surface) triangular and tetrahedral elements
  write.table(
    format(matrix(rep(vtk_cell_type, n_elems), nrow = 1, byrow = T), scientific = F),
    file = vtu, append = T, row.names = F, col.names = F, quote = F
  )
  vtu %<<% "</DataArray>\n"
  vtu %<<% "</Cells>\n"

  if (!is.null(data)) {
    nvdata <- ncol(data) ## number of signals to plot
    vtu %<<% "<PointData>\n"
    for (i in seq_len(nvdata)) {
      vtu %<<% paste0('<DataArray type="Float64" Name="data', i, '" NumberOfComponents="1" format="ascii">\n')
      write.table(matrix(data[, i], nrow = 1), file = vtu, append = T, row.names = F, col.names = F)
      vtu %<<% "</DataArray>\n"
    }
    vtu %<<% "</PointData>\n"
  }
  ## footer
  vtu %<<% "</Piece>\n"
  vtu %<<% "</UnstructuredGrid>\n"
  vtu %<<% "</VTKFile>\n"

  ## close file
  close(vtu)
}

## missing data import and boundary import

#' @export
read_vtu <- function(file) {
  ## open file descriptor in read-only mode
  vtu <- file(file, open = "r")
  ## informations parsed from file
  n_points <- 0
  n_cells <- 0
  dim <- 0
  mesh_data <- list()

  read_xml_property <- function(line, name) {
    l <- regexpr(paste0(name, '=\"[0-9]+\"'), line)
    if (l[1] == -1) stop("XML property ", deparse(substitute(name)), " not found.")
    return(as.numeric(substr(line, l[1] + nchar(paste0(name, '=\"')), l[1] + attr(l, "match.length") - 2)))
  }

  linenum <- 1 ## current line number
  status <- 0 ## parsing status (1: <Points> detected, 2: <Cells> detected)
  while (TRUE) {
    line <- readLines(vtu, n = 1)
    if (length(line) == 0) { ## EOF
      break
    }
    if (grepl("VTKFile", line) && !grepl("/VTKFile", line)) {
      ## check if type is UnstructeredGrid, stop if not
      if (!grepl("UnstructuredGrid", line)) stop("only UnstructedGrid VTK format supported.")
    }
    if (grepl("Piece", line) && !grepl("/Piece", line)) {
      n_points <- read_xml_property(line, "NumberOfPoints")
      n_cells <- read_xml_property(line, "NumberOfCells")
    }
    ## <Points> ... </Points> block parsing logic
    if (status == 1) {
      if (grepl("DataArray", line) && !grepl("/DataArray", line)) {
        dim <- read_xml_property(line, "NumberOfComponents")
        ## read points matrix
        mesh_data$nodes <- as.matrix(read.table(vtu, nrows = n_points))
        colnames(mesh_data$nodes) <- NULL
        status <- 0
      }
    }
    ## <Cells> ... </Cells> block parsing logic
    if (status == 2) {
      if (grepl("DataArray", line) && grepl("connectivity", line)) {
        ## read elements matrix
        mesh_data$elements <- as.matrix(read.table(vtu, nrows = n_cells))
        colnames(mesh_data$elements) <- NULL
        status <- 0
      }
    }
    if (grepl("<Points>", line)) status <- 1
    if (grepl("<Cells>" , line)) status <- 2
    linenum <- linenum + 1
  }
  if (nrow(mesh_data$nodes) == 0 || nrow(mesh_data$elements) == 0) stop("VTK parsing error.")
  close(vtu) ## close file connection

  mesh_data$boundary <- matrix(0, nrow(mesh_data$nodes))
  local_dim <- if (ncol(mesh_data$elements) == 4) 3 else 2 
  return(.Mesh$new(
    mesh = new(eval(parse(text = paste("cpp_mesh", local_dim, dim, sep = "_"))), mesh_data),
    local_dim = local_dim,
    embed_dim = dim
  ))
}
