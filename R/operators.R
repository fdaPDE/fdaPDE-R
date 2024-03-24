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

## base class for differential operators
## .DiffOpCtr <- setRefClass(
##     Class = "DiffOpObject",
##     fields = c(
##         tokens = "vector",
##         params = "list",
##         f = "FunctionObject"
##     ),
##     methods = c(
##         function_space = function() { return(f$function_space) }
##     )
## )

## ## sum of differential operators
## setMethod("+",
##     signature = c(e1 = "DiffOpObject", e2 = "DiffOpObject"),
##     function(e1, e2) {
##         if (tracemem(e1$f) != tracemem(e2$f)) {
##             stop("operator arguments must be the same")
##         }
##         ## check for duplicated operator tokens
##         if (any(duplicated(c(e1$tokens, e2$tokens)))) {
##             stop("duplicated operator")
##         }
##         .DiffOpCtr(
##             tokens = c(e1$tokens, e2$tokens),
##             params = c(e1$params, e2$params),
##             f = e1$f
##         )
##     }
## )

## # difference of differential operators
## setMethod("-",
##     signature = c(e1 = "DiffOpObject", e2 = "DiffOpObject"),
##     function(e1, e2) {
##         if (tracemem(e1$f) != tracemem(e2$f)) {
##             stop("operator arguments must be the same")
##         }
##         ## check for duplicated operator tokens
##         if (any(duplicated(c(e1$tokens, e2$tokens)))) {
##             stop("duplicated operator")
##         }
##         for (i in 1:length(e2$params)) {
##             e2$params[[i]] <- -e2$params[[i]]
##         }
##         .DiffOpCtr(
##             tokens = c(e1$tokens, e2$tokens),
##             params = c(e1$params, e2$params),
##             f = e1$f
##         )
##     }
## )
## ## minus (unary) operator for DiffOpObject
## setMethod(
##     "-", signature(e1 = "DiffOpObject", e2 = "missing"),
##     function(e1) {
##         e1$params[[1]] <- -e1$params[[1]]
##         e1
##     }
## )

## ## differential operator product by scalar
## setMethod("*",
##     signature = c(e1 = "numeric", e2 = "DiffOpObject"),
##     function(e1, e2) {
##         # if (!is.numeric(e1)) stop("bad product")
##         e2$params[[1]] <- e1 * e2$params[[1]]
##         e2
##     }
## )

## ## diffusion term
## .DiffusionCtr <- setRefClass(
##     Class = "DiffusionOperator",
##     contains = "DiffOpObject"
## )
## ## laplace() returns a special operator for the case of
## ## isotropic and stationary diffusion
## laplace <- function(f) {
##     if (!is(f, "FunctionObject")) {
##         stop("wrong argument type")
##     }
##     .DiffusionCtr(
##         tokens = "diffusion",
##         params = list(diffusion = 1),
##         f = f
##     )
## }
## ## the general non-isotrpic, non-stationary diffusion operator
## div <- function(f) {
##     if (is(f, "FunctionGradObject")) {
##         if (!is.null(f$K)) {
##             return(.DiffusionCtr(
##                 tokens = "diffusion",
##                 params = list(diffusion = f$K),
##                 f = f$f
##             ))
##         }
##     }
##     stop("wrong argument type")
## }

## ## transport term
## .TransportCtr <- setRefClass(
##     Class = "TransportOperator",
##     contains = "DiffOpObject"
## )
## setGeneric("dot", function(op1, op2) standardGeneric("dot"))
## setMethod(
##     "dot", signature(op1 = "vector", op2 = "FunctionGradObject"),
##     function(op1, op2) {
##         .TransportCtr(
##             tokens = "transport",
##             params = list(transport = as.matrix(op1)),
##             f = op2$f
##         )
##     }
## )

## ## reaction term
## .ReactionCtr <- setRefClass(
##     Class = "ReactionOperator",
##     contains = "DiffOpObject"
## )
## setMethod("*",
##     signature = c(e1 = "numeric", e2 = "FunctionObject"),
##     function(e1, e2) {
##         .ReactionCtr(
##             tokens = "reaction",
##             params = list(reaction = e1),
##             f = e2
##         )
##     }
## )

## .TimeDerivativeCtr <- setRefClass(
##     Class = "TimeDerivative",
##     contains = "DiffOpObject"
## )
## # overloading stats::dt function
## setMethod("dt",
##     signature = c(x = "FunctionObject", df = "missing", ncp = "missing"),
##     function(x, df, ncp) {
##         .TimeDerivativeCtr(
##             tokens = "time",
##             params = list(time = 1L),
##             f = f
##         )
##     }
## )
