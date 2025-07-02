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

## type of supported layers
layer_t <- list(point = 0, areal = 1)
data_t <- list(flt64 = 0, flt32 = 1, int64 = 2, int32 = 3, bin = 4, str = 5)

## .geoframe_ctr <- R6::R6Class(
##   "geoframe",
##   private = list(
##     geoframe_ = "geoframe_t",
##     triangulation_ = "triangulation_t",
##     layer_map_ = list() ## maps layer names to layer types
##   ),
##   public = list(
##     initialize = function(...) {
##       args <- list(...)
##       if (length(args) == 1) {
##         if ("triangulation_2_2" %in% class(args[[1]])) {
##           private$triangulation_ <- args[[1]]
##           private$geoframe_ <- new(cpp_geoframe_2_2, args[[1]]$.__enclos_env__$private$mesh_)
##         } else {
##           stop("Invalid argument to geoframe initializer.")
##         }
##       }
##       print("CREO GEOFRAME")
##     },
##     plot = function() {
##       n_layers = length(private$layer_map_)
##       n_col = ceiling(sqrt(n_layers))
##       n_row = ceiling(n_layers / n_col)
##       par(mfrow = c(n_row, n_col))
##       for (i in seq(from = 1, to = n_layers)) {
##         if (private$layer_map_[[i]] == "areal") {
##           plot(.areal_layer$new(self, names(private$layer_map_)[i]))
##           title(names(private$layer_map_)[i])
##         }
##         if (private$layer_map_[[i]] == "point") {
##           plot(.point_layer$new(self, names(private$layer_map_)[i]))
##           title(names(private$layer_map_)[i])
##         }
##       }
##     },
##     print = function(...) {
##       n_layers <- length(private$layer_map_)
##       layer_names <- names(private$layer_map_)
##       bbox <- private$geoframe_$bbox()
##       n_nodes <- private$geoframe_$n_nodes()
##       n_cells <- private$geoframe_$n_cells()
##       cat("Geoframe with ", n_layers, " layers\n", sep = "")
##       cat(
##         "Bounding box:    xmin: ",
##         bbox[1],
##         " ymin: ",
##         bbox[3],
##         " xmax: ",
##         bbox[2],
##         " ymax: ",
##         bbox[4],
##         "\n",
##         sep = ""
##       )
##       cat("Number of nodes: ", n_nodes, "\n", sep = "")
##       cat("Number of cells: ", n_cells, "\n", sep = "")
##       cat("\n")
##       if (n_layers > 0) {
##         for (i in seq(from = 1, to = n_layers)) {
##           cat("Layer: ", layer_names[i], "\n", sep = "")
##           cat("Type:  ", toupper(private$layer_map_[[i]]), "\n", sep = "")
##           cat(
##             "Dims: ",
##             private$geoframe_$rows(layer_names[i]),
##             ", ",
##             private$geoframe_$cols(layer_names[i]),
##             "\n",
##             sep = ""
##           )
##           if (private$layer_map_[[i]] == "areal") {
##             layer <- .areal_layer$new(self, layer_names[i])
##             cat("First ", min(6, private$geoframe_$rows(layer_names[i])), " data rows:\n", sep = "")
##             print(layer)
##           }
##           if (private$layer_map_[[i]] == "point") {
##             layer <- .point_layer$new(self, layer_names[i])
##             cat("First ", min(6, private$geoframe_$rows(layer_names[i])), " data rows:\n", sep = "")
##             print(layer)
##           }
##         }
##       }
##       cat("")
##     },
##     ## modifiers
##     insert = function(layer, type, geo, data) {
##         ## divide data by types
##         env <- new.env()
##         env$data_ <- list()
##         env$data_[["int_data"]] = list()
##         env$data_[["dbl_data"]] = list()
##         env$data_[["str_data"]] = list()

##         load <- function(a, b, colname) { ## copies a inside b, depending on a's class
##             if(is.numeric  (a)) b$data_[["dbl_data"]][[colname]] = as.matrix(a)
##             if(is.integer  (a)) b$data_[["int_data"]][[colname]] = as.matrix(a)
##             if(is.character(a)) b$data_[["str_data"]][[colname]] = as.matrix(a)
##         }

##         if (is.matrix(data) || is.data.frame(data)) {
##             n_col = dim(data)[2]
##             for (i in seq(from = 1, to = n_col)) {
##                 colname = if (is.null(colnames(data))) paste("V", i, sep = "") else colnames(data)[i]
##                 if (is.character(geo)) {
##                     if (!(colname %in% geo)) load(data[, i], env, colname)
##                 } else {
##                     load(data[, i], env, colname)
##                 }
##             }
##         }
##         ## load geometrical informations
##         geo_ = matrix(0, nrow = 0, ncol = 0)
##         if (is.character(geo)) {
##             ## geometry is referenced as columns of data
##             geo_ = as.matrix(data[, geo])
##         } else {
##             geo_ = geo
##         }
##         if(type == "point") {
##             private$geoframe_$insert_scalar_point_layer(layer, geo_, env$data_)
##             private$layer_map_[[layer]] <- "point"
##         }
##         if(type == "areal") {
##             private$geoframe_$insert_scalar_areal_layer(layer, geo_, env$data_)
##             private$layer_map_[[layer]] <- "areal"
##         }
##     },
##     load = function(layer, filename) {
##         private$geoframe_$load_shp(layer, filename)
##         private$layer_map_[[layer]] <- "areal"
##     },
##     ## subsetting
##     `[[` = function(layer_name) {
##         print("ACCESS LAYER")

##       if (!layer_name %in% names(private$layer_map_)) stop(paste("Layer ", layer_name, " not found.", sep = ""))
##       if (private$layer_map_[[layer_name]] == "areal") return(gf_areal_layer(private$geoframe_, layer_name)) ##return(invisible(.areal_layer$new(private$geoframe_, layer_name)))
##       if (private$layer_map_[[layer_name]] == "point") return(invisible(.point_layer$new(private$geoframe_, layer_name)))
##     },
##     rows = function(layer_name) {
##       private$geoframe_$rows(layer_name)
##     },
##     `==` = function(rhs) {
##       if (length(private$layer_map_) == 1) {
##         layer_name = names(private$layer_map_)[1]
##         return(.areal_layer$new(self, layer_name) == rhs)
##       }
##     },
##     `!=` = function(rhs) {
##       if (length(private$layer_map_) == 1) {
##         layer_name = names(private$layer_map_)[1]
##         return(.areal_layer$new(self, layer_name) != rhs)
##       }
##     },
##     `>=` = function(rhs) {
##       if (length(private$layer_map_) == 1) {
##         layer_name = names(private$layer_map_)[1]
##         return(.areal_layer$new(self, layer_name) >= rhs)
##       }
##     },
##     `<=` = function(rhs) {
##       if (length(private$layer_map_) == 1) {
##         layer_name = names(private$layer_map_)[1]
##         return(.areal_layer$new(self, layer_name) <= rhs)
##       }
##     },
##     `>` = function(rhs) {
##       if (length(private$layer_map_) == 1) {
##         layer_name = names(private$layer_map_)[1]
##         return(.areal_layer$new(self, layer_name) > rhs)
##       }
##     },
##     `<` = function(rhs) {
##       if (length(private$layer_map_) == 1) {
##         layer_name = names(private$layer_map_)[1]
##         return(.areal_layer$new(self, layer_name) < rhs)
##       }
##     }
##   ),
##   active = list(
##     colnames = function() private$geoframe_$colnames_all()
##   )
## )

## ## #' @export
## ## geoframe <- function(triangulation) {
## ##   structure(
## ##     list(ptr = .geoframe_ctr$new(triangulation)),
## ##     class = "gf"
## ##   )
## ## }

## .areal_layer <- R6::R6Class(
##   "geoframe_areal_layer",
##   private = list(
##     geoframe_ = "geoframe",
##     layer_name_ = "character",
##     geoframe_handler = function() {
##       private$geoframe_
##     }
##   ),
##   public = list(
##       initialize = function(geoframe, layer_name) {
##         private$geoframe_ <- geoframe
##       private$layer_name_ <- layer_name

##     },
##     ## plot = function() {
##     ##   ## recover all coordinates of polygons
##     ##   nodes <- private$geoframe_$.__enclos_env__$private$geoframe_$areal_poly_nodes(private$layer_name_)
##     ##   edges <- private$geoframe_$.__enclos_env__$private$geoframe_$areal_poly_edges(private$layer_name_)
##     ##   ## open plotting device
##     ##   x_range <- c(+Inf, -Inf)
##     ##   y_range <- c(+Inf, -Inf)
##     ##   for (i in seq(from = 1, to = length(nodes))) {
##     ##     x_r <- range(nodes[[i]][, 1])
##     ##     y_r <- range(nodes[[i]][, 2])
##     ##     if (x_range[1] > x_r[1]) x_range[1] <- x_r[1]
##     ##     if (x_range[2] < x_r[2]) x_range[2] <- x_r[2]
##     ##     if (y_range[1] > y_r[1]) y_range[1] <- y_r[1]
##     ##     if (y_range[2] < y_r[2]) y_range[2] <- y_r[2]
##     ##   }
##     ##   plot(NA, xlim = x_range, ylim = y_range, xlab = "", ylab = "", asp = 1)
##     ##   ## plot each polygon
##     ##   for (i in seq(from = 1, to = length(nodes))) {
##     ##     segments(
##     ##       x0 = nodes[[i]][edges[[i]][, 1], 1],
##     ##       y0 = nodes[[i]][edges[[i]][, 1], 2],
##     ##       x1 = nodes[[i]][edges[[i]][, 2], 1],
##     ##       y1 = nodes[[i]][edges[[i]][, 2], 2]
##     ##     )
##     ##   }
##     ## },
##     print = function() {
##         print("PRINTO")
##       output <- list()
##       colnames <- private$geoframe_handler()$colnames(private$layer_name_)
##       rows <- seq(from = 0, to = min(5, private$geoframe_handler()$rows(private$layer_name_) - 1))
##       for (i in seq(from = 1, to = length(colnames))) {
##         output[[i]] <- rep("", times = 2 + length(rows))
##         output[[i]][1] <- colnames[i]
##         dtype <- private$geoframe_handler()$ctype(private$layer_name_, colnames[i])
##         if (dtype == data_t$flt64) output[[i]][2] <- "<flt64>"
##         if (dtype == data_t$flt32) output[[i]][2] <- "<flt32>"
##         if (dtype == data_t$int64) output[[i]][2] <- "<int64>"
##         if (dtype == data_t$int32) output[[i]][2] <- "<int32>"
##         if (dtype == data_t$bin) output[[i]][2] <- "<bin>"
##         if (dtype == data_t$str) output[[i]][2] <- "<chr>"

##         if (dtype == data_t$flt64) v <- private$geoframe_handler()$flt64_access(private$layer_name_, rows, colnames[i])
##         if (dtype == data_t$flt32) v <- private$geoframe_handler()$flt32_access(private$layer_name_, rows, colnames[i])
##         if (dtype == data_t$int64) v <- private$geoframe_handler()$int64_access(private$layer_name_, rows, colnames[i])
##         if (dtype == data_t$int32) v <- private$geoframe_handler()$int32_access(private$layer_name_, rows, colnames[i])
##         if (dtype == data_t$bin) v <- private$geoframe_handler()$bin_access(private$layer_name_, rows, colnames[i])
##         if (dtype == data_t$str) v <- private$geoframe_handler()$str_access(private$layer_name_, rows, colnames[i])

##         for (j in seq(from = 1, to = length(v))) {
##           output[[i]][j + 2] <- as.character(v[j])
##         }
##       }
##                 print("PRINTO")

##       ## pretty format
##       for (i in seq(from = 1, to = length(colnames))) {
##         max_length <- max(nchar(output[[i]])) ## maximum string length per column
##         output[[i]] <- format(output[[i]], width = (max_length + 1), justify = "left")
##       }
##       ## output to console
##       for (i in seq(from = 1, to = length(colnames))) {
##         cat(output[[i]][1])
##       }
##       cat("\n")
##       for (i in seq(from = 1, to = length(colnames))) {
##         cat(paste0("\033[0;", 31, "m", output[[i]][2], "\033[0m"))
##       }
##       cat("\n")
##       for (i in seq(from = 3, to = (2 + length(rows)))) {
##         for (j in seq(from = 1, to = length(colnames))) {
##           cat(format(output[[j]][i], digits = 6))
##         }
##         cat("\n")
##       }
##                 print("PRINTO")

##     },
##     ## subsetting operation
##     get = function(rows, cols) {
##       rows <- if (missing(rows)) {
##         as.vector(seq(from = 0, to = (private$geoframe_$rows(private$layer_name_) - 1)))
##       } else {
##         if (is.logical(rows)) {
##           as.vector(which(rows) - 1)
##         } else {
##           as.vector(rows - 1)
##         }
##       }
##       if (missing(cols)) {
##         cols <- private$geoframe_handler()$colnames(private$layer_name_) ## take all columns
##       } else {
##         if (!is.character(cols)) {
##           cols <- private$geoframe_handler()$colnames(private$layer_name_)[cols]
##         } else {
##           cols <- as.vector(cols)
##         }
##       }
##       ## create a new geoframe
##       gf <- .geoframe_ctr$new()
##       gf$.__enclos_env__$private$geoframe_ <- new(
##         cpp_geoframe_2_2,
##         private$geoframe_handler(), ## there is a copy made here...
##         private$layer_name_,
##         rows,
##         cols
##       )
##       gf$.__enclos_env__$private$layer_map_[[private$layer_name_]] <- "areal"
##       return(gf)
##     },
##     set = function(rows, cols, value) {
##         print("SETTO")
##       ## ctype <- private$geoframe_handler()$ctype(private$layer_name_, cols)
##       ## ## prepare row subsetting vector
##       ## if (missing(rows)) {
##       ##   rows <- as.vector(seq(from = 0, to = (private$geoframe_handler()$rows(private$layer_name_) - 1)))
##       ## } else {
##       ##   rows <- as.vector(rows - 1) ## cpp aligned indexes
##       ## }
##       ## if (length(value) == 1) value <- as.vector(rep(value, times = length(rows)))
##       ## ## dispatch to typed assignment logic
##       ## if (ctype == data_t$flt64) private$geoframe_handler()$flt64_assign(layer_name, rows, cols, as.numeric(value))
##       ## if (ctype == data_t$flt32) private$geoframe_handler()$flt32_assign(layer_name, rows, cols, as.numeric(value))
##       ## if (ctype == data_t$int64) private$geoframe_handler()$int64_assign(layer_name, rows, cols, as.integer(value))
##       ## if (ctype == data_t$int32) private$geoframe_handler()$int32_assign(layer_name, rows, cols, as.integer(value))
##       ## if (ctype == data_t$bin)   private$geoframe_handler()$bin_assign  (layer_name, rows, cols, as.logical(value))
##       ## if (ctype == data_t$str)   private$geoframe_handler()$str_assign  (layer_name, rows, cols, as.character(value))
##     },
##     sample = function(n_samples, seed) {
##       seed_ <- if (missing(seed)) -1 else seed
##       return(private$geoframe_handler()$areal_sample(private$layer_name_, n_samples, seed_))
##     },
##     ## comparison operators
##     `==` = function(rhs) return(.gf_apply_cmp(private$geoframe_, private$layer_name_, rhs, function(x, y) x == y)),
##     `!=` = function(rhs) return(.gf_apply_cmp(private$geoframe_, private$layer_name_, rhs, function(x, y) x != y)),
##     `>=` = function(rhs) return(.gf_apply_cmp(private$geoframe_, private$layer_name_, rhs, function(x, y) x >= y)),
##     `<=` = function(rhs) return(.gf_apply_cmp(private$geoframe_, private$layer_name_, rhs, function(x, y) x <= y)),
##     `>`  = function(rhs) return(.gf_apply_cmp(private$geoframe_, private$layer_name_, rhs, function(x, y) x > y)),
##     `<`  = function(rhs) return(.gf_apply_cmp(private$geoframe_, private$layer_name_, rhs, function(x, y) x < y))
##   )
## )

## ## S3 wrapper
## ## required to have assignment-like semantic (R6 classes are locked environments)
## #' @export
## gf_areal_layer <- function(geoframe, layer_name) {
##   structure(
##     list(ptr = .areal_layer$new(geoframe, layer_name)),
##     class = "gf_areal_layer"
##   )
## }

## ## Assignment via [<-
## #' @export
## `[<-.gf_areal_layer` <- function(x, rows, cols, value) {
##     print("STO APPLICANDO QUESTO??")
##     x$ptr$set(rows, cols, value)
##     print("SETTATO")
## }

## #' @export
## `[.gf_areal_layer` <- function(x, rows, cols) {
##   x$ptr$get(rows, cols)
## }

## .gf_apply_cmp <- function(geoframe, layer_name, rhs, f) {
##   cpp_handler = geoframe$.__enclos_env__$private$geoframe_
##   n_rows = cpp_handler$rows(layer_name)
##   n_cols = cpp_handler$cols(layer_name)
##   output <- matrix(rep(FALSE, times = n_rows * n_cols), nrow = n_rows, ncol = n_cols)

##   rows = as.vector(seq(from = 0, to = (n_rows - 1)))
##   colnames <- cpp_handler$colnames(layer_name)

##   for (i in seq(from = 1, to = length(colnames))) {
##     ctype <- cpp_handler$ctype(layer_name, colnames[i])
##     if (ctype == data_t$flt64) v <- cpp_handler$flt64_access(layer_name, rows, colnames[i])
##     if (ctype == data_t$flt32) v <- cpp_handler$flt32_access(layer_name, rows, colnames[i])
##     if (ctype == data_t$int64) v <- cpp_handler$int64_access(layer_name, rows, colnames[i])
##     if (ctype == data_t$int32) v <- cpp_handler$int32_access(layer_name, rows, colnames[i])
##     if (ctype == data_t$bin) v <- cpp_handler$bin_access(layer_name, rows, colnames[i])
##     if (ctype == data_t$str) v <- cpp_handler$str_access(layer_name, rows, colnames[i])
##     output[, i] = f(v, rhs)
##   }
##   return(output)
## }

## .point_layer <- R6::R6Class(
##   "geoframe_point_layer",
##   private = list(
##     geoframe_ = "geoframe",
##     layer_name_ = "character",
##     geoframe_handler = function() private$geoframe_$.__enclos_env__$private$geoframe_
##   ),
##   public = list(
##     initialize = function(geoframe, layer_name) {
##       private$geoframe_ <- geoframe
##       private$layer_name_ <- layer_name
##     },
##     ## subsetting
##     `[` = function(rows, cols) {
##       rows <- if (missing(rows)) {
##         as.vector(seq(from = 0, to = (private$geoframe_$rows(private$layer_name_) - 1)))
##       } else {
##         if (is.logical(rows)) {
##           as.vector(which(rows) - 1)
##         } else {
##           as.vector(rows - 1)
##         }
##       }
##       if (missing(cols)) {
##         cols <- private$geoframe_handler()$colnames(private$layer_name_) ## take all columns
##       } else {
##         if (!is.character(cols)) {
##           cols <- private$geoframe_handler()$colnames(private$layer_name_)[cols]
##         } else {
##           cols <- as.vector(cols)
##         }
##       }
##       ## create a new geoframe
##       print("STO SUBSETTANDO QUI?")
##       gf <- .geoframe_ctr$new()
##       gf$.__enclos_env__$private$geoframe_ <- new(
##         cpp_geoframe_2_2,
##         private$geoframe_handler(),
##         private$layer_name_,
##         rows,
##         cols
##       )
##       gf$.__enclos_env__$private$layer_map_[[private$layer_name_]] <- "point"
##       print("SI, COGLIONE")

##       return(gf)
##     },
##     ## assignment
##     set = function(rows, cols, vals) {
##       layer_name_ <- private$layer_name_
##       ctype <- private$geoframe_handler()$ctype(layer_name_, cols)
##       ## prepare row subsetting vector
##       if (missing(rows)) {
##         rows <- as.vector(seq(from = 0, to = (private$geoframe_handler()$rows(layer_name_) - 1)))
##       } else {
##         rows <- as.vector(rows - 1) ## cpp aligned indexes
##       }
##       if (length(vals) == 1) vals <- as.vector(rep(vals, times = length(rows))) else {
##         fdapde_assert(length(vals) == length(rows))
##       }
##       ## dispatch to typed assignment logic
##       if (ctype == data_t$flt64) private$geoframe_handler()$flt64_assign(layer_name_, rows, cols, as.numeric(vals))
##       if (ctype == data_t$flt32) private$geoframe_handler()$flt32_assign(layer_name_, rows, cols, as.numeric(vals))
##       if (ctype == data_t$int64) private$geoframe_handler()$int64_assign(layer_name_, rows, cols, as.integer(vals))
##       if (ctype == data_t$int32) private$geoframe_handler()$int32_assign(layer_name_, rows, cols, as.integer(vals))
##       if (ctype == data_t$bin) private$geoframe_handler()$bin_assign(layer_name_, rows, cols, as.logical(vals))
##       if (ctype == data_t$str) private$geoframe_handler()$str_assign(layer_name_, rows, cols, as.character(vals))
##       invisible(self)
##     },
##     plot = function(covs = NULL, mesh = TRUE, ...) {
##       coords = private$geoframe_handler()$point_coordinates(private$layer_name_)
##       x_range <- range(coords[, 1])
##       y_range <- range(coords[, 2])

##       par(mar = c(1, 1, 1, 1))
##       plot(private$geoframe_$.__enclos_env__$private$triangulation_)
##       if (is.null(covs)) {
##         points(coords, xlim = x_range, ylim = y_range, xlab = "", ylab = "", asp = 1, col = "red", pch = 21, cex = 1)
##       } else {
##         n_col <- 50
##         palette <- colorRampPalette(colors = c("lightyellow", "darkred"))(n_col)
##         ## if you have strings, use as many colors as different values of strings
##         ## if you have binary, use 2 colors
##         ## else, if numeric
##         ## create value-palette mapping
##         points(coords, xlim = x_range, ylim = y_range, xlab = "", ylab = "", asp = 1, col = palette, pch = 21, cex = 1)
##       }
##     },
##     print = function() {
##       output <- list()
##       colnames <- private$geoframe_handler()$colnames(private$layer_name_)
##       rows <- seq(from = 0, to = min(5, private$geoframe_handler()$rows(private$layer_name_) - 1))
##       for (i in seq(from = 1, to = length(colnames))) {
##         output[[i]] <- rep("", times = 2 + length(rows))
##         output[[i]][1] <- colnames[i]
##         dtype <- private$geoframe_handler()$ctype(private$layer_name_, colnames[i])
##         if (dtype == data_t$flt64) output[[i]][2] <- "<flt64>"
##         if (dtype == data_t$flt32) output[[i]][2] <- "<flt32>"
##         if (dtype == data_t$int64) output[[i]][2] <- "<int64>"
##         if (dtype == data_t$int32) output[[i]][2] <- "<int32>"
##         if (dtype == data_t$bin) output[[i]][2] <- "<bin>"
##         if (dtype == data_t$str) output[[i]][2] <- "<chr>"

##         if (dtype == data_t$flt64) v <- private$geoframe_handler()$flt64_access(private$layer_name_, rows, colnames[i])
##         if (dtype == data_t$flt32) v <- private$geoframe_handler()$flt32_access(private$layer_name_, rows, colnames[i])
##         if (dtype == data_t$int64) v <- private$geoframe_handler()$int64_access(private$layer_name_, rows, colnames[i])
##         if (dtype == data_t$int32) v <- private$geoframe_handler()$int32_access(private$layer_name_, rows, colnames[i])
##         if (dtype == data_t$bin) v <- private$geoframe_handler()$bin_access(private$layer_name_, rows, colnames[i])
##         if (dtype == data_t$str) v <- private$geoframe_handler()$str_access(private$layer_name_, rows, colnames[i])

##         for (j in seq(from = 1, to = length(v))) {
##           output[[i]][j + 2] <- as.character(v[j])
##         }
##       }
##       ## pretty format
##       for (i in seq(from = 1, to = length(colnames))) {
##         max_length <- max(nchar(output[[i]])) ## maximum string length per column
##         output[[i]] <- format(output[[i]], width = (max_length + 1), justify = "left")
##       }
##       ## output to console
##       for (i in seq(from = 1, to = length(colnames))) {
##         cat(output[[i]][1])
##       }
##       cat("\n")
##       for (i in seq(from = 1, to = length(colnames))) {
##         cat(paste0("\033[0;", 31, "m", output[[i]][2], "\033[0m"))
##       }
##       cat("\n")
##       for (i in seq(from = 3, to = (2 + length(rows)))) {
##         for (j in seq(from = 1, to = length(colnames))) {
##           cat(output[[j]][i])
##         }
##         cat("\n")
##       }
##     },
##     ## comparison operators
##     `==` = function(rhs) return(.gf_apply_cmp(private$geoframe_, private$layer_name_, rhs, function(x, y) x == y)),
##     `!=` = function(rhs) return(.gf_apply_cmp(private$geoframe_, private$layer_name_, rhs, function(x, y) x != y)),
##     `>=` = function(rhs) return(.gf_apply_cmp(private$geoframe_, private$layer_name_, rhs, function(x, y) x >= y)),
##     `<=` = function(rhs) return(.gf_apply_cmp(private$geoframe_, private$layer_name_, rhs, function(x, y) x <= y)),
##     `>` = function(rhs) return(.gf_apply_cmp(private$geoframe_, private$layer_name_, rhs, function(x, y) x > y)),
##     `<` = function(rhs) return(.gf_apply_cmp(private$geoframe_, private$layer_name_, rhs, function(x, y) x < y))
##   )
## )

## #' @export
## geoframe <- function(triangulation) {
##   return(.geoframe_ctr$new(triangulation))
## }

## #' @export
## `[[.geoframe` <- function(x, ...) {
##   x$`[[`(...)
## }

## ## #' @export
## ## `[.geoframe_areal_layer` <- function(x, ...) {
## ##   x$`[`(...)
## ## }

## ## #' @export
## ## `[.geoframe_point_layer` <- function(x, ...) {
## ##   x$`[`(...)
## ## }

## #' @export
## `==.geoframe_areal_layer` <- function(x, ...) x$`==`(...)
## #' @export
## `!=.geoframe_areal_layer` <- function(x, ...) x$`!=`(...)
## #' @export
## `>=.geoframe_areal_layer` <- function(x, ...) x$`>=`(...)
## #' @export
## `<=.geoframe_areal_layer` <- function(x, ...) x$`<=`(...)
## #' @export
## `>.geoframe_areal_layer` <- function(x, ...) x$`>`(...)
## #' @export
## `<.geoframe_areal_layer` <- function(x, ...) x$`<`(...)

## #' @export
## `==.geoframe_point_layer` <- function(x, ...) x$`==`(...)
## #' @export
## `!=.geoframe_point_layer` <- function(x, ...) x$`!=`(...)
## #' @export
## `>=.geoframe_point_layer` <- function(x, ...) x$`>=`(...)
## #' @export
## `<=.geoframe_point_layer` <- function(x, ...) x$`<=`(...)
## #' @export
## `>.geoframe_point_layer` <- function(x, ...) x$`>`(...)
## #' @export
## `<.geoframe_point_layer` <- function(x, ...) x$`<`(...)

## ## comparison operators
## #' @export
## `==.geoframe` <- function(x, ...) x$`==`(...)
## #' @export
## `!=.geoframe` <- function(x, ...) x$`!=`(...)
## #' @export
## `>=.geoframe` <- function(x, ...) x$`>=`(...)
## #' @export
## `<=.geoframe` <- function(x, ...) x$`<=`(...)
## #' @export
## `>.geoframe` <- function(x, ...) x$`>`(...)
## #' @export
## `<.geoframe` <- function(x, ...) x$`<`(...)

## #' @export
## gf_load <- function(geoframe, layer, file) {
##   return(geoframe$load_shp(layer, file))
## }

## #' @export
## gf_sample <- function(geoframe, layer, size) {
##   return(.areal_layer$new(geoframe, layer)$sample(size))
## }

# Define R6 class
.geoframe <- R6::R6Class(
  "cpp_gf",
  private = list(
    ptr_ = NULL,
    mesh_ = NULL,
    layer_map_ = list() ## maps layer names to layer types
  ),
  public = list(
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
    gf__insert__ = function(layer, type, geo, data) {
      ## divide data by types
      env <- new.env()
      env$data_ <- list()
      env$data_[["int_data"]] = list()
      env$data_[["dbl_data"]] = list()
      env$data_[["str_data"]] = list()

      load <- function(a, b, colname) {
        ## copies a inside b, depending on a's class
        if (is.numeric(a)) b$data_[["dbl_data"]][[colname]] = as.matrix(a)
        if (is.integer(a)) b$data_[["int_data"]][[colname]] = as.matrix(a)
        if (is.character(a)) b$data_[["str_data"]][[colname]] = as.matrix(a)
      }

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
      ## load geometrical informations
      geo_ = matrix(0, nrow = 0, ncol = 0)
      if (is.character(geo)) {
        ## geometry is referenced as columns of data
        geo_ = as.matrix(data[, geo])
      } else {
        geo_ = geo
      }
      if (type == "point") {
        private$ptr_$insert_scalar_point_layer(layer, geo_, env$data_)
        private$layer_map_[[layer]] <- "point"
      }
      if (type == "areal") {
        private$ptr_$insert_scalar_areal_layer(layer, geo_, env$data_)
        private$layer_map_[[layer]] <- "areal"
      }
    },
    load_shp = function(layer, filename) {
      private$ptr_$load_shp(layer, filename)
      private$layer_map_[[layer]] <- "areal"
    },
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
    gf__layer__ = function(layer_name) {
      if (!layer_name %in% names(private$layer_map_)) stop(paste("Layer ", layer_name, " not found.", sep = ""))
      if (private$layer_map_[[layer_name]] == "areal") return(gf_areal(private$ptr_, layer_name))
      if (private$layer_map_[[layer_name]] == "point") return(gf_point(private$ptr_, layer_name))
    }
  )
)

inject_r6_to_s3 <- function(s3obj, r6obj, blacklist = c("initialize", "clone")) {
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
  s3obj
}

# S3 constructor
#' @export
geoframe <- function(domain) {
  ptr <- .geoframe$new(domain)
  obj <- list(gf__ptr__ = ptr)
  obj <- inject_r6_to_s3(obj, ptr)
  class(obj) <- "gf"
  return(obj)
}

#' @export
print.gf <- function(x) {
  x$gf__ptr__$gf__print__()
}

#' @export
`[[.gf` <- function(x, ...) {
    x$gf__ptr__$gf__layer__(...)
}

## areal layer
.areal_layer <- R6::R6Class(
  "cpp_gf_areal",
  private = list(
    ptr_ = NULL,
    name_ = NULL
  ),
  public = list(
    initialize = function(geoframe, layer_name) {
      private$ptr_ <- geoframe
      private$name_ <- layer_name
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
      ptr$.__enclos_env__$private$layer_map_[[private$name_]] <- "areal"
      obj <- list(gf__ptr__ = ptr)
      obj <- inject_r6_to_s3(obj, ptr)
      class(obj) <- "gf"
      return(obj)
    },
    gf__set__ = function(rows, cols, value) {
      ctype <- private$ptr_$ctype(private$name_, cols)
      ## prepare row subsetting vector
      if (missing(rows)) {
        rows <- as.vector(seq(from = 0, to = (private$ptr_$rows(private$name_) - 1)))
      } else {
        rows <- as.vector(rows - 1) ## cpp aligned indexes
      }
      if (length(value) == 1) value <- as.vector(rep(value, times = length(rows)))
      ## dispatch to typed assignment logic
      if (ctype == data_t$flt64) private$ptr_$flt64_assign(private$name_, rows, cols, as.numeric(value))
      if (ctype == data_t$flt32) private$ptr_$flt32_assign(private$name_, rows, cols, as.numeric(value))
      if (ctype == data_t$int64) private$ptr_$int64_assign(private$name_, rows, cols, as.integer(value))
      if (ctype == data_t$int32) private$ptr_$int32_assign(private$name_, rows, cols, as.integer(value))
      if (ctype == data_t$bin)   private$ptr_$bin_assign  (private$name_, rows, cols, as.logical(value))
      if (ctype == data_t$str)   private$ptr_$str_assign  (private$name_, rows, cols, as.character(value))
    },
    gf__print__ = function(...) {
      output <- list()
      colnames <- private$ptr_$colnames(private$name_)
      rows <- seq(from = 0, to = min(5, private$ptr_$rows(private$name_) - 1))
      for (i in seq(from = 1, to = length(colnames))) {
        output[[i]] <- rep("", times = 2 + length(rows))
        output[[i]][1] <- colnames[i]
        dtype <- private$ptr_$ctype(private$name_, colnames[i])
        if (dtype == data_t$flt64) output[[i]][2] <- "<flt64>"
        if (dtype == data_t$flt32) output[[i]][2] <- "<flt32>"
        if (dtype == data_t$int64) output[[i]][2] <- "<int64>"
        if (dtype == data_t$int32) output[[i]][2] <- "<int32>"
        if (dtype == data_t$bin)   output[[i]][2] <- "<bin>"
        if (dtype == data_t$str)   output[[i]][2] <- "<chr>"

        if (dtype == data_t$flt64) v <- private$ptr_$flt64_access(private$name_, rows, colnames[i])
        if (dtype == data_t$flt32) v <- private$ptr_$flt32_access(private$name_, rows, colnames[i])
        if (dtype == data_t$int64) v <- private$ptr_$int64_access(private$name_, rows, colnames[i])
        if (dtype == data_t$int32) v <- private$ptr_$int32_access(private$name_, rows, colnames[i])
        if (dtype == data_t$bin)   v <- private$ptr_$bin_access  (private$name_, rows, colnames[i])
        if (dtype == data_t$str)   v <- private$ptr_$str_access  (private$name_, rows, colnames[i])

        for (j in seq(from = 1, to = length(v))) {
          output[[i]][j + 2] <- as.character(v[j])
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
      name = function() { return(private$name_) },
      rows = function() { return(private$ptr_$rows(private$name_)) },
      colnames = function() { return(private$ptr_$colnames(private$name_)) }
  )
)

# Print method
gf_areal <- function(geoframe, name) {
  ptr <- .areal_layer$new(geoframe, name)
  obj <- list(gf__ptr__ = ptr)
  obj <- inject_r6_to_s3(obj, ptr)
  class(obj) <- "gf_areal"
  return(obj)
}

## S3 subsetting getter
#' @export
`[.gf_areal` <- function(x, rows, cols) {
  x$gf__ptr__$gf__get__(rows, cols)
}

## S3 subsetting setter
#' @export
`[<-.gf_areal` <- function(x, rows, cols, value) {
  x$gf__ptr__$gf__set__(rows, cols, value)
}

#' @export
print.gf_areal <- function(x) {
  x$gf__ptr__$gf__print__()
}

#' @export
`$.gf_areal` <- function(x, name) {
  if (name %in% names(x)) {
    return(x[[name]])
  } else {
      ## access column data from the backend
      cpp_backend <- get_private(x$gf__ptr__)$ptr_
      layer_name  <- x$gf__ptr__$name
      nrows       <- x$gf__ptr__$rows
      
      rows  <- as.vector(seq(from = 0, to = (nrows - 1)))
      dtype <- cpp_backend$ctype(layer_name, name)

      if (dtype == data_t$flt64) return(cpp_backend$flt64_access(layer_name, rows, name))
      if (dtype == data_t$flt32) return(cpp_backend$flt32_access(layer_name, rows, name))
      if (dtype == data_t$int64) return(cpp_backend$int64_access(layer_name, rows, name))
      if (dtype == data_t$int32) return(cpp_backend$int32_access(layer_name, rows, name))
      if (dtype == data_t$bin)   return(cpp_backend$bin_access  (layer_name, rows, name))
      if (dtype == data_t$str)   return(cpp_backend$str_access  (layer_name, rows, name))
  }
}

#' @export
`$<-.gf_areal` <- function(x, name, value) {
  if (name %in% names(x)) {
    x[[name]] <- value
  } else {
      cpp_backend <- get_private(x$gf__ptr__)$ptr_
      layer_name  <- x$gf__ptr__$name
      nrows       <- x$gf__ptr__$rows
      
      rows  <- as.vector(seq(from = 0, to = (nrows - 1)))
      dtype <- cpp_backend$ctype(layer_name, name)

      if (length(value) == 1) value <- as.vector(rep(value, times = length(rows)))
      ## dispatch to typed assignment logic
      if (dtype == data_t$flt64) cpp_backend$flt64_assign(layer_name, rows, name, as.numeric(value))
      if (dtype == data_t$flt32) cpp_backend$flt32_assign(layer_name, rows, name, as.numeric(value))
      if (dtype == data_t$int64) cpp_backend$int64_assign(layer_name, rows, name, as.integer(value))
      if (dtype == data_t$int32) cpp_backend$int32_assign(layer_name, rows, name, as.integer(value))
      if (dtype == data_t$bin)   cpp_backend$bin_assign  (layer_name, rows, name, as.logical(value))
      if (dtype == data_t$str)   cpp_backend$str_assign  (layer_name, rows, name, as.character(value))

  }
}

#' @export
gf_colnames <- function(x, ...) {
  return(x$gf__ptr__$colnames)
}


## bisogna dare nomi improponibili agli elementi R6
