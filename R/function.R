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

## a symbolic expression
.SymbolicExpression <- R6::R6Class("SymbolicExpression",
  private = list(
    expr_ = character(),
    symbolic_ = character(),
    symbol_table_ = list()
  ),
  public = list(
    initialize = function(expr = NA, symbolic = NA, symbol_table = list()) {
      private$expr_ <- expr
      private$symbol_table_ <- symbol_table
      private$symbolic_ <- symbolic
    },
    print = function(...) {
      cat("Symbolic expression: ", private$symbolic_, "\n", sep = "")
    }
  ),
  active = list(
    expr = function() private$expr_,
    symbol_table = function() private$symbol_table_,
    symbolic = function() private$symbolic_
  )
)

## local utilities
parenthesize <- function(expr, apply = TRUE) ifelse(apply, paste0("(", expr, ")"), expr)
is.symbolic_function <- function(obj) inherits(obj, "SymbolicFunction")
is.symbolic_expression <- function(obj) inherits(obj, "SymbolicExpression")
is.composite_expression <- function(obj) is.symbolic_expression(obj) && !is.symbolic_function(obj)
## returns false if the evaluation of expression e is invalid
try_eval <- function(e) {
  success <- TRUE
  tryCatch(e, error = function(e) success <<- FALSE)
  return(success)
}
random_name <- function(length = 10) {
  name <- paste0(sample(c(letters, "1", "2", "3", "4", "5", "6", "7", "8", "9"), length, TRUE),
    collapse = ""
  )
  return(paste0("<", name, ">"))
}

merge_symbolics <- function(obj1, obj2, expr1, expr2, OP) {
  ## check if operands require parenthesization
  p1 <- regexpr("^\\((.*)\\)", expr1)[1] != -1
  p2 <- regexpr("^\\((.*)\\)", expr2)[1] != -1
  ## obj1 and obj2 are both composite expressions or symbolic functions
  if ((is.composite_expression(obj1) && is.composite_expression(obj2)) ||
    (is.symbolic_function(obj1) && is.symbolic_function(obj2))) {
    return(
      .SymbolicExpression$new(
        expr = paste0(parenthesize(obj1$expr, p1), as.character(OP), parenthesize(obj2$expr, p2)),
        symbolic = if (is.symbolic_function(obj1)) {
          paste0(expr1, as.character(OP), expr2)
        } else {
          paste0(obj1$symbolic, as.character(OP), obj2$symbolic)
        },
        symbol_table = modifyList(obj1$symbol_table, obj2$symbol_table)
      )
    )
  } else {
    identifier <- list()
    value_name <- parenthesize(random_name(), p1)
    ## obj1 is either a value or a symbolic function, and obj2 is a symbolic expression
    if (!is.symbolic_expression(obj1) || is.symbolic_function(obj1)) {
      identifier[[value_name]] <- obj1
      e <- obj2 ## this must inherits from SymbolicExpression
      expr <- paste0(value_name, as.character(OP), parenthesize(e$expr, p2))
      symb <- if (is.symbolic_function(obj2)) {
        paste0(expr1, as.character(OP), expr2)
      } else {
        paste0(expr1, as.character(OP), e$symbolic)
      }
    } else { ## obj2 must be either a value or a symbolic function
      identifier[[value_name]] <- obj2
      e <- obj1 ## this must inherits from SymbolicExpression
      expr <- paste0(parenthesize(e$expr, p1), as.character(OP), value_name)
      symb <- if (is.symbolic_function(obj1)) {
        paste0(expr1, as.character(OP), expr2)
      } else {
        paste0(e$symbolic, as.character(OP), expr2)
      }
    }
    if (is.symbolic_function(e)) { ## create new symbolic expression
      return(
        .SymbolicExpression$new(
          expr = expr,
          symbolic = symb,
          symbol_table = modifyList(identifier, e$symbol_table)
        )
      )
    } else { ## modify symbolic in place
      set_private(e, "expr_", expr)
      set_private(e, "symbolic_", symb)
      set_private(e, "symbol_table_", modifyList(identifier, e$symbol_table))
      return(e)
    }
  }
}

#' addition between symbolic expressions
#' @export
`+.SymbolicExpression` <- function(op1, op2) {
  merge_symbolics(obj1 = op1, obj2 = op2, expr1 = deparse(substitute(op1)), expr2 = deparse(substitute(op2)), OP = "+")
}

#' unary negation and subtraction between symbolic expressions
#' @export
`-.SymbolicExpression` <- function(op1, op2 = NULL) {
  negate <- function(x) {
    if (is.function(x) || is.symbolic_function(x)) {
      return(x)
    } else {
      return(-x)
    }
  }
  expr1 <- deparse(substitute(op1))
  p1 <- regexpr("^\\((.*)\\)", expr1)[1] != -1
  if (is.null(op2)) { ## unary minus
    set_private(op1, "expr_", paste0("-", parenthesize(op1$expr, p1)))
    set_private(op1, "symbolic_", if (is.symbolic_function(op1)) paste0("-", expr1) else paste0("-", op1$symbolic))
    set_private(op1, "symbol_table_", lapply(op1$symbol_table, negate))
    return(op1)
  }
  expr2 <- deparse(substitute(op2))
  p2 <- regexpr("^\\((.*)\\)", expr2)[1] != -1
  if (is.symbolic_expression(op1) && is.symbolic_expression(op2)) {
    set_private(op2, "symbol_table_", lapply(op2$symbol_table, negate))
  } else { ## one between op1 and op2 is a value
    if (is.symbolic_expression(op2)) {
      set_private(op2, "symbol_table_", lapply(op2$symbol_table, negate))
    } else {
      op2 <- -op2
    }
  }
  merge_symbolics(obj1 = op1, obj2 = op2, expr1 = deparse(substitute(op1)), expr2 = deparse(substitute(op2)), OP = "-")
}

#' product between symbolic expressions
#' @export
`*.SymbolicExpression` <- function(op1, op2) {
  merge_symbolics(obj1 = op1, obj2 = op2, expr1 = deparse(substitute(op1)), expr2 = deparse(substitute(op2)), OP = "*")
}

unary_symbolic <- function(tag, obj, expr) {
  if (!inherits(obj, "SymbolicExpression")) {
    if (!try_eval(obj)) stop("object ", expr, " not found.")
  }
  if (is.symbolic_function(obj)) {
    ## create new symbolic expression
    return(.SymbolicExpression$new(
      expr = paste0(tag, parenthesize(obj$expr)),
      symbolic = paste0(tag, "(", expr, ")"),
      symbol_table = obj$symbol_table
    ))
  }
  set_private(obj, "expr_", paste0(tag, parenthesize(obj$expr)))
  set_private(obj, "symbolic_", expr)
  return(obj)
}

binary_symbolic <- function(tag, obj1, obj2, expr1, expr2) {
  if (is.symbolic_expression(obj1) && is.symbolic_expression(obj2)) {
    return(
      .SymbolicExpression$new(
        expr = paste0(tag, "(", obj1$expr, ",", obj2$expr, ")"),
        symbolic = paste0(tag, "(", expr1, ",", expr2, ")"),
        symbol_table = modifyList(obj1$symbol_table, obj2$symbol_table)
      )
    )
  } else { ## one between e1 and e2 is a value
    if (is.symbolic_expression(obj1)) {
      v <- obj2
      e <- obj1
    } else {
      v <- obj1
      e <- obj2
    }
    value_name <- random_name()
    expr <- if (!is.symbolic_expression(obj1)) {
      paste0(tag, "(", value_name, ",", e$expr, ")")
    } else {
      paste0(tag, "(", e$expr, ",", value_name, ")")
    }
    identifier <- list()
    identifier[[value_name]] <- v
    if (is.symbolic_function(e)) { ## create new symbolic expression
      return(.SymbolicExpression$new(
        expr = expr,
        symbolic = paste0(tag, "(", expr1, ",", expr2, ")"),
        symbol_table = modifyList(identifier, e$symbol_table)
      ))
    } else { ## modify symbolic in place
      set_private(e, "expr_", expr)
      set_private(e, "symbolic_", paste0(tag, "(", expr1, ",", expr2, ")"))
      set_private(e, "symbol_table_", modifyList(identifier, e$symbol_table))
      return(e)
    }
  }
}

grad <- function(op) unary_symbolic("grad", op, deparse(substitute(op)))
div <- function(op) unary_symbolic("div", op, deparse(substitute(op)))
laplace <- function(op, name = NULL) {
  if (is.null(name)) unary_symbolic("laplace", op, deparse(substitute(op))) else unary_symbolic("laplace", op, name)
}
dt <- function(op) unary_symbolic("dt", op, deparse(substitute(op)))
inner <- function(op1, op2) binary_symbolic("inner", op1, op2, deparse(substitute(op1)), deparse(substitute(op2)))

## symbolic mathematical function
.SymbolicFunction <- R6::R6Class("SymbolicFunction",
  inherit = .SymbolicExpression,
  private = list(
    functional_space_ = NULL,
    coeff_ = "matrix" ## expansion coefficient vector
  ),
  public = list(
    initialize = function(functional_space_ = NA) {
      private$functional_space_ <- functional_space_
      private$coeff_ <- matrix(ncol = 1, nrow = 0)
      private$expr_ <- random_name()
      private$symbol_table_[[private$expr_]] <- self
    },
    eval = function(locations) {
      return(ifelse(
        inherits(functional_space_, "TensorProductSpace"),
        private$functional_space_$eval(locations) %*% private$coeff_,
        private$functional_space_$eval_expansion(private$coeff_, locations)
      ))
    },
    print = function(...) {
      cat("SymbolicFunction\n", sep = "")
    }
  ),
  ## active binding
  active = list(
    coeff = function() private$coeff_
  )
)

#' @export
Function <- function(functional_space) {
  return(.SymbolicFunction$new(functional_space))
}
