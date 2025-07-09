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

## supported layer and data types
layer_t <- list(point = 0, areal = 1)
data_t  <- list(flt64 = 0, flt32 = 1, int64 = 2, int32 = 3, bin = 4, str = 5)

.gf_inject_r6_to_s3 <- function(s3obj, r6obj, blacklist = c("initialize", "clone")) {
  methods <- names(r6obj)  
  for (m in methods) {
    ## blacklisted methods or marked as gf__ are not exposed at S3 level (except gf__ptr__)
    if (!m %in% blacklist && !m %in% names(s3obj) && (!grepl("gf__", m) || m == "gf__ptr__")) {
      local({
        method_name <- m
        if (is.function(r6obj[[method_name]])) {
          f <- function(...) {
            do.call(r6obj[[method_name]], list(...))
          }
          s3obj[[method_name]] <<- f
        }
      })
    }
  }
  return(s3obj)
}

#' An R6 class representing a geo referenced dataframe
#'
#'@rdname geoframe
.geoframe <- R6::R6Class(
  "cpp_gf",
  private = list(
    ptr_ = NULL,
    mesh_ = NULL,
    layer_map_ = list() ## maps layer names to layer types
  ),
  public = list(
    #' @description
    #' Creates a new geoframe object.
    #'
    #' @param ... Arguments created internally.
    initialize = function(...) {
      args <- list(...)
      if (length(args) == 1) {
        if ("triangulation_2_2" %in% class(args[[1]])) {
          private$mesh_ <- args[[1]]
          private$ptr_ <- new(cpp_geoframe_2_2, args[[1]]$.__enclos_env__$private$mesh_)
        } else {
          stop("Invalid argument to geoframe initializer.")
        }
      }
    },
    ## modifiers
    #' @description
    #' Adds a new data layer.
    #'
    #' @param layer A string indicating the name of the data layer.
    #' @param type A string indicating the type of the layer. It can be either \code{"point"} or \code{"areal"}.
    #' @param geo A character vector indicating the names of the columns in \code{data} to be used as coordinates.
    #' Default is \code{NULL}.
    #' @param data A data frame containing the data. Default is \code{NULL}.
    insert = function(layer, type, geo = NULL, data = NULL) {
      env <- new.env()
      env$data_ <- list()
      env$data_[["int_data"]] = list()
      env$data_[["dbl_data"]] = list()
      env$data_[["str_data"]] = list()

      load <- function(a, b, colname) {
        ## copies a inside b, depending on a's class
        if (is.numeric(a))   b$data_[["dbl_data"]][[colname]] = as.matrix(a)
        if (is.integer(a))   b$data_[["int_data"]][[colname]] = as.matrix(a)
        if (is.character(a)) b$data_[["str_data"]][[colname]] = as.matrix(a)
      }

      if (!is.null(data)) {
        if (is.matrix(data) || is.data.frame(data)) {
          n_col = dim(data)[2]
          for (i in seq(from = 1, to = n_col)) {
            colname = if (is.null(colnames(data))) paste("V", i, sep = "") else colnames(data)[i]
            if (is.character(geo)) {
              if (!(colname %in% geo)) load(data[, i], env, colname)
            } else {
              load(data[, i], env, colname)
            }
          }
        }
      }
      ## load geometrical informations
      if (is.null(geo)) {
        fdapde_assert(!is.null(geo) || type == "point", "Missing geometry.")
        private$ptr_$insert_scalar_point_layer_mesh_nodes(layer, 0, env$data_)
        private$layer_map_[[layer]] <- "point"
      } else {
        geo_ = matrix(0, nrow = 0, ncol = 0)
        if (is.character(geo)) {
          ## geometry is referenced as columns of data
          geo_ = as.matrix(data[, geo])
        } else {
          geo_ = as.matrix(geo)
        }
        if (type == "point") {
          private$ptr_$insert_scalar_point_layer(layer, geo_, env$data_)
          private$layer_map_[[layer]] <- "point"
        }
        if (type == "areal") {
          private$ptr_$insert_scalar_areal_layer(layer, geo_, env$data_)
          private$layer_map_[[layer]] <- "areal"
        }
      }
    },
    #' @description
    #' Adds a new data layer from a shapefile.
    #'
    #' @param layer A string indicating the name of the data layer.
    #' @param filename A string containing the path to the shapefile to be loaded.
    load_shp = function(layer, filename) {
      if (!file.exists(filename)) stop(paste("File", filename, "not found", sep = " "))
      private$ptr_$load_shp(layer, filename)
      private$layer_map_[[layer]] <- "areal"
    },
    #' @description
    #' Print method for the \code{geoframe} class.
    #'
    #' @param ... Ellipsis for compatibility with \code{plot.base}.
    gf__print__ = function(...) {
      n_layers <- length(private$layer_map_)
      layer_names <- names(private$layer_map_)
      bbox <- private$ptr_$bbox()
      n_nodes <- private$ptr_$n_nodes()
      n_cells <- private$ptr_$n_cells()
      cat("Geoframe with ", n_layers, " layers\n", sep = "")
      cat(
        "Bounding box:    xmin: ",
        bbox[1],
        " ymin: ",
        bbox[3],
        " xmax: ",
        bbox[2],
        " ymax: ",
        bbox[4],
        "\n",
        sep = ""
      )
      cat("Number of nodes: ", n_nodes, "\n", sep = "")
      cat("Number of cells: ", n_cells, "\n", sep = "")
      cat("\n")
      if (n_layers > 0) {
        for (i in seq(from = 1, to = n_layers)) {
          cat("Layer: ", layer_names[i], "\n", sep = "")
          cat("Type:  ", toupper(private$layer_map_[[i]]), "\n", sep = "")
          cat(
            "Dims: ",
            private$ptr_$rows(layer_names[i]),
            ", ",
            private$ptr_$cols(layer_names[i]),
            "\n",
            sep = ""
          )
          if (private$layer_map_[[i]] == "areal") {
            layer <- gf_areal(private$ptr_, layer_names[i])
            cat("First ", min(6, private$ptr_$rows(layer_names[i])), " data rows:\n", sep = "")
            print(layer)
          }
          if (private$layer_map_[[i]] == "point") {
            layer <- gf_point(private$ptr_, layer_names[i])
            cat("First ", min(6, private$ptr_$rows(layer_names[i])), " data rows:\n", sep = "")
            print(layer)
          }
        }
      }
      cat("")
    },
    ## subsetting
    #' @description
    #' Access method for the \code{geoframe} class.
    #'
    #' @param layer_name A string specifying the name of the data layer to be returned.
    #' @return The requested data layer.
    gf__layer__ = function(layer_name) {
      if (!layer_name %in% names(private$layer_map_)) stop(paste("Layer ", layer_name, " not found.", sep = ""))
      if (private$layer_map_[[layer_name]] == "areal") return(gf_areal(private$ptr_, layer_name))
      if (private$layer_map_[[layer_name]] == "point") return(gf_point(private$ptr_, layer_name))
    }
  ),
  active = list(
    #' @field colnames
    #' Returns the column names
    colnames = function() return(private$ptr_$colnames_all()),
    geometry = function() return(private$mesh_)
  )
)

#' Create a Geoframe
#'
#' @param domain An object of the triangulation class created by [triangulation()].
#' @return An R6 object representing a geo-referenced data frame.
#' @rdname geoframe
#' @export
#' @examples
#' \dontrun{
#' library(RTriangle)
#' library(fdaPDE2)
#' p <- pslg(P=rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)),
#' S=rbind(c(1, 2), c(2, 3), c(3, 4), c(4,1)))
#' mesh_data <- triangulate(p, a = 0.00125, q=30)
#' mesh <- triangulation(nodes = mesh_data$P, cells = mesh_data$T, boundary = mesh_data$PB)
#' points <- mesh.sample(100)
#' response <- rnorm(100)
#' data <- data.frame(p1=points[,1], p2=points[,2], y=response)
#' gf <- geoframe(domain = florida)
#' gf$insert(layer = "temp", type = "point", geo = c("p1", "p2"), data = data)
#' gf
#' }
geoframe <- function(domain) {
  ptr <- .geoframe$new(domain)
  obj <- list(gf__ptr__ = ptr)
  obj <- .gf_inject_r6_to_s3(obj, ptr)
  class(obj) <- "gf"
  return(obj)
}

#' @export
print.gf <- function(x) {
  x$gf__ptr__$gf__print__()
}

#' @export
`[[.gf` <- function(x, layer_name) {
  x$gf__ptr__$gf__layer__(layer_name)
}

#' @export
gf_geometry <- function(x) {
  fdapde_assert("gf" %in% class(g), "Not a geoframe.")
  return(x$gf__ptr__$geometry)
}

## low-level typed dispatch logic
.gf_cpp_access <- function(x, rows, col) {
  r_backend <- if ("gf__ptr__" %in% ls(x)) x[["gf__ptr__"]] else x
  cpp_backend <- get_private(r_backend)$ptr_
  layer_name <- r_backend$name
  dtype <- cpp_backend$dtype(layer_name, col)

  v <- NULL
  if (dtype == data_t$flt64) v <- cpp_backend$flt64_access(layer_name, rows, col)
  if (dtype == data_t$flt32) v <- cpp_backend$flt32_access(layer_name, rows, col)
  if (dtype == data_t$int64) v <- cpp_backend$int64_access(layer_name, rows, col)
  if (dtype == data_t$int32) v <- cpp_backend$int32_access(layer_name, rows, col)
  if (dtype == data_t$bin)   v <- cpp_backend$bin_access(layer_name, rows, col)
  if (dtype == data_t$str)   v <- cpp_backend$str_access(layer_name, rows, col)
  return(v)
}

.gf_cpp_assign <- function(x, rows, col, value) {
  r_backend <- if ("gf__ptr__" %in% ls(x)) x[["gf__ptr__"]] else x
  cpp_backend <- get_private(r_backend)$ptr_
  layer_name <- r_backend$name
  dtype <- cpp_backend$dtype(layer_name, col)

  if (dtype == data_t$flt64) cpp_backend$flt64_assign(layer_name, rows, col, as.numeric(value))
  if (dtype == data_t$flt32) cpp_backend$flt32_assign(layer_name, rows, col, as.numeric(value))
  if (dtype == data_t$int64) cpp_backend$int64_assign(layer_name, rows, col, as.integer(value))
  if (dtype == data_t$int32) cpp_backend$int32_assign(layer_name, rows, col, as.integer(value))
  if (dtype == data_t$bin)   cpp_backend$bin_assign  (layer_name, rows, col, as.logical(value))
  if (dtype == data_t$str)   cpp_backend$str_assign  (layer_name, rows, col, as.character(value))
}

.gf_cpp_insert <- function(x, col, value) {
  r_backend <- if ("gf__ptr__" %in% ls(x)) x[["gf__ptr__"]] else x
  cpp_backend <- get_private(r_backend)$ptr_
  layer_name <- r_backend$name
  nrows <- r_backend$rows
  if (length(value) == 1) value <- as.vector(rep(value, times = nrows))

  if (is.numeric(value))   cpp_backend$flt64_insert(layer_name, name, as.numeric(value))
  if (is.integer(value))   cpp_backend$int64_insert(layer_name, name, as.integer(value))
  if (is.character(value)) cpp_backend$str_insert  (layer_name, name, as.character(value))
}

## low-level data managment logic
.data_layer <- R6::R6Class(
  "cpp_gf_data",
  private = list(
    ptr_  = NULL, ## geoframe cpp backend
    name_ = NULL,
    type_ = NULL
  ),
  public = list(
    initialize = function(ptr, name, type) {
      private$ptr_ <- ptr
      private$name_ <- name
      fdapde_assert(type == "point" || type == "areal", "Invalid type.")
      private$type_ <- type
    },
    ## subsetting operation
    gf__get__ = function(rows, cols) {
      rows <- if (missing(rows)) {
        as.vector(seq(from = 0, to = (private$ptr_$rows(private$name_) - 1)))
      } else {
        if (is.logical(rows)) {
          as.vector(which(rows) - 1)
        } else {
          as.vector(rows - 1)
        }
      }
      if (missing(cols)) {
        cols <- private$ptr_$colnames(private$name_) ## take all columns
      } else {
        if (!is.character(cols)) {
          cols <- private$ptr_$colnames(private$name_)[cols]
        } else {
          cols <- as.vector(cols)
        }
      }
      ## create a new geoframe
      ptr <- .geoframe$new()
      ptr$.__enclos_env__$private$ptr_ <- new(
        cpp_geoframe_2_2,
        private$ptr_,
        private$name_,
        rows,
        cols
      )
      ptr$.__enclos_env__$private$layer_map_[[private$name_]] <- private$type_
      obj <- list(gf__ptr__ = ptr)
      obj <- .gf_inject_r6_to_s3(obj, ptr)
      class(obj) <- "gf"
      return(obj)
    },
    gf__set__ = function(rows, col, value) {
      dtype <- private$ptr_$dtype(private$name_, col)
      ## prepare row subsetting vector
      if (missing(rows)) {
        rows <- as.vector(seq(from = 0, to = (private$ptr_$rows(private$name_) - 1)))
      } else {
        rows <- as.vector(rows - 1) ## cpp aligned indexes
      }
      if (length(value) == 1) value <- as.vector(rep(value, times = length(rows)))
      .gf_cpp_assign(self, rows, col, value)
    },
    gf__print__ = function(...) {
      output <- list()
      colnames <- private$ptr_$colnames(private$name_)
      rows <- seq(from = 0, to = min(5, private$ptr_$rows(private$name_) - 1))
      for (i in seq(from = 1, to = length(colnames))) {
        v <- .gf_cpp_access(self, rows, colnames[i])

        output[[i]] <- rep("", times = 2 + length(rows))
        output[[i]][1] <- colnames[i]
        dtype <- private$ptr_$dtype(private$name_, colnames[i])
        if (dtype == data_t$flt64) output[[i]][2] <- "<flt64>"
        if (dtype == data_t$flt32) output[[i]][2] <- "<flt32>"
        if (dtype == data_t$int64) output[[i]][2] <- "<int64>"
        if (dtype == data_t$int32) output[[i]][2] <- "<int32>"
        if (dtype == data_t$bin)   output[[i]][2] <- "<bin>"
        if (dtype == data_t$str)   output[[i]][2] <- "<chr>"
        for (j in seq(from = 1, to = length(v))) {
          output[[i]][j + 2] <- format(v[j], digits = 6) ## rounds to 6 decimal points
        }
      }
      ## pretty format
      for (i in seq(from = 1, to = length(colnames))) {
        max_length <- max(nchar(output[[i]])) ## maximum string length per column
        output[[i]] <- format(output[[i]], width = (max_length + 1), justify = "left")
      }
      ## output to console
      for (i in seq(from = 1, to = length(colnames))) {
        cat(output[[i]][1])
      }
      cat("\n")
      for (i in seq(from = 1, to = length(colnames))) {
        cat(paste0("\033[0;", 31, "m", output[[i]][2], "\033[0m"))
      }
      cat("\n")
      for (i in seq(from = 3, to = (2 + length(rows)))) {
        for (j in seq(from = 1, to = length(colnames))) {
          cat(format(output[[j]][i], digits = 6))
        }
        cat("\n")
      }
    }
  ),
  active = list(
    name = function() return(private$name_),
    rows = function() return(private$ptr_$rows(private$name_)),
    cols = function() return(length(private$ptr_$colnames(private$name_))),
    colnames = function() return(private$ptr_$colnames(private$name_))
  )
)

#' @export
`$.gf_data` <- function(x, colname) {
  ## access column data from backend
  nrows <- x[["gf__ptr__"]]$rows
  rows <- as.vector(seq(from = 0, to = (nrows - 1)))
  .gf_cpp_access(x, rows, colname)
}

#' @export
`$<-.gf_data` <- function(x, colname, value) {
  cpp_backend <- get_private(x[["gf__ptr__"]])$ptr_
  layer_name <- x[["gf__ptr__"]]$name
  nrows <- x[["gf__ptr__"]]$rows

  if (length(value) == 1) value <- as.vector(rep(value, times = nrows))
  if (colname %in% x[["gf__ptr__"]]$colnames) {
    ## modify in place
    rows <- as.vector(seq(from = 0, to = (nrows - 1)))
    .gf_cpp_assign(x, rows, colname, value)
  } else {
    ## column insertion
    fdapde_assert(length(value) == 1 || length(value) == nrows, "Invalid assignment.")
    .gf_cpp_insert(x, colname, value)
  }
}


## areal layer
.areal_layer <- R6::R6Class(
  inherit = .data_layer,
  "cpp_gf_areal",
  active = list(
      polygons = function() {
          polygons <- private$ptr_$areal_polygons(private$name_)
          for (i in seq_len(length(polygons))) {
            polygons[[i]]$edges <- r_aligned_index(polygons[[i]]$edges)
          }
          return(polygons)
      }
  )
)

gf_areal <- function(geoframe, name) {
  ptr <- .areal_layer$new(geoframe, name, "areal")
  obj <- list(gf__ptr__ = ptr)
  obj <- .gf_inject_r6_to_s3(obj, ptr)
  class(obj) <- c("gf_areal", "gf_data")
  return(obj)
}

#' @export
`[.gf_areal` <- function(x, rows, cols) {
  return(x[["gf__ptr__"]]$gf__get__(rows, cols))
}

#' @export
`[<-.gf_areal` <- function(x, rows, cols, value) {
  x[["gf__ptr__"]]$gf__set__(rows, cols, value)
}

#' @export
print.gf_areal <- function(x) {
  x[["gf__ptr__"]]$gf__print__()
}

#' @export
names.gf_areal <- function(x, ...) {
  return(x[["gf__ptr__"]]$colnames)
}

#' @export
dim.gf_areal <- function(x, ...) {
  return(c(x[["gf__ptr__"]]$rows, x[["gf__ptr__"]]$cols))
}

#' @export
gf_polygons <- function(x) {
  fdapde_assert("gf_areal" %in% class(x), "Invalid layer type.")
  return(x[["gf__ptr__"]]$polygons)
}

## point layer
.point_layer <- R6::R6Class(
  inherit = .data_layer,
  "cpp_gf_point",
  public = list(
    gf__plot__ = function(covs = NULL, mesh = TRUE, ...) {
      coords = private$geoframe_handler()$point_coordinates(private$name_)
      x_range <- range(coords[, 1])
      y_range <- range(coords[, 2])

      par(mar = c(1, 1, 1, 1))
      plot(private$geoframe_$.__enclos_env__$private$triangulation_)
      col <- NULL
      if (is.null(covs)) {
        col <- "red"
      } else {
        n_col <- 50
        palette <- colorRampPalette(colors = c("lightyellow", "darkred"))(n_col)
        col <- palette[as.numeric(cut(vals, breaks = 50))]
      }
      points(coords, xlim = x_range, ylim = y_range, xlab = "", ylab = "", asp = 1, col = col, pch = 21, cex = 1)
    }
  ),
  active = list(
    coordinates = function() return(private$ptr_$point_coordinates(private$name_))
  )
)

gf_point <- function(geoframe, name) {
  ptr <- .point_layer$new(geoframe, name, "point")
  obj <- list(gf__ptr__ = ptr)
  obj <- .gf_inject_r6_to_s3(obj, ptr)
  class(obj) <- c("gf_point", "gf_data")
  return(obj)
}

#' @export
`[.gf_point` <- function(x, rows, cols) {
  return(x[["gf__ptr__"]]$gf__get__(rows, cols))
}

#' @export
`[<-.gf_point` <- function(x, rows, cols, value) {
  x[["gf__ptr__"]]$gf__set__(rows, cols, value)
}

#' @export
print.gf_point <- function(x) {
  x[["gf__ptr__"]]$gf__print__()
}

#' @export
dim.gf_point <- function(x) {
  return(c(x[["gf__ptr__"]]$rows, x[["gf__ptr__"]]$cols))
}

#' @export
names.gf_point <- function(x) {
  return(x[["gf__ptr__"]]$colnames)
}

#' @export
gf_locations <- function(x) {
  fdapde_assert("gf_point" %in% class(x), "Invalid layer type.")
  return(x[["gf__ptr__"]]$coordinates)
}
