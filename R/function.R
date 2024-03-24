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
    expr_ = "character",
    symbol_table_ = list()
  ),
  public = list(
    initialize = function(symbolic_expr = NA, symbol_table = list()) {
      private$expr_ <- symbolic_expr
      private$symbol_table_ <- symbol_table
    },
    print = function(...) {
      cat("<SymbolicExpression>: ", private$expr_, "\n", sep = "")
      cat("Symbols Table\n", sep = "")
      for (p in names(private$symbol_table_)) {
        cat(p, "\n", sep = "")
        print(private$symbol_table_[[p]])
      }
    }
  ),
  active = list(
    expr = function() private$expr_,
    symbol_table = function() private$symbol_table_
  )
)

## local utilities
parenthesize <- function(expr, apply = TRUE) ifelse(apply, paste("(", expr, ")", sep = ""), expr)
is.symbolic_function <- function(sym) inherits(sym, "SymbolicFunction")
is.symbolic_expression <- function(sym) inherits(sym, "SymbolicExpression")
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
  return(paste("<", name, ">", sep = ""))
}

merge_symbolics <- function(expr1, expr2, OP) {
  ## check if operands require parenthesization
  p1 <- regexpr("^\\((.*)\\)", expr1)[1] != -1
  p2 <- regexpr("^\\((.*)\\)", expr2)[1] != -1
  ## recursively evaluate operands
  e1 <- eval(parse(text = expr1))
  e2 <- eval(parse(text = expr2))
  if (is.symbolic_expression(e1) && is.symbolic_expression(e2)) {
    return(
      .SymbolicExpression$new(
        symbolic_expr = paste(parenthesize(e1$expr, p1), as.character(OP), parenthesize(e2$expr, p2)),
        symbol_table = modifyList(e1$symbol_table, e2$symbol_table)
      )
    )
  } else { ## one between e1 and e2 is a value
    if (is.symbolic_expression(e1)) {
      v <- e2
      e <- e1
    } else {
      v <- e1
      e <- e2
    }
    value_name <- parenthesize(random_name(), p1)
    symbolic_expr <- ifelse(!is.symbolic_expression(e1),
      paste(value_name, as.character(OP), parenthesize(e$expr, p2), sep = ""),
      paste(parenthesize(e$expr, p2), as.character(OP), value_name, sep = "")
    )
    identifier <- list()
    identifier[[value_name]] <- v

    if (is.symbolic_function(e)) { ## create new symbolic expression
      return(.SymbolicExpression$new(
        symbolic_expr = symbolic_expr,
        symbol_table = modifyList(identifier, e$symbol_table)
      ))
    } else { ## modify symbolic in place
      set_private(e, "expr_", symbolic_expr)
      set_private(e, "symbol_table_", modifyList(identifier, e$symbol_table))
      return(e)
    }
  }
}

#' addition between symbolic expressions
#' @export
`+.SymbolicExpression` <- function(op1, op2) {
  merge_symbolics(expr1 = deparse(substitute(op1)), expr2 = deparse(substitute(op2)), OP = "+")
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
  p1 <- regexpr("^\\((.*)\\)", deparse(substitute(op1)))[1] != -1
  if (is.null(op2)) { ## unary minus
    set_private(op1, "expr_", paste("-", parenthesize(op1$expr, p1)))
    set_private(op1, "symbol_table_", lapply(op1$symbol_table, negate))
    return(op1)
  }
  p2 <- regexpr("^\\((.*)\\)", deparse(substitute(op2)))[1] != -1
  if (is.symbolic_expression(op1) && is.symbolic_expression(op2)) {
    set_private(op2, "symbol_table_", lapply(op2$symbol_table, negate))
  } else { ## one between op1 and op2 is a value
    if (is.symbolic_expression(op2)) {
      set_private(op2, "symbol_table_", lapply(op2$symbol_table, negate))
    } else {
      op2 <- -op2
    }
  }
  merge_symbolics(expr1 = deparse(substitute(op1)), expr2 = deparse(substitute(op2)), OP = "-")
}

#' product between symbolic expressions
#' @export
`*.SymbolicExpression` <- function(op1, op2) {
  merge_symbolics(expr1 = deparse(substitute(op1)), expr2 = deparse(substitute(op2)), OP = "*")
}

unary_symbolic <- function(tag, expr) {
  e <- eval(parse(text = expr))
  if (!inherits(e, "SymbolicExpression")) {
    if (!try_eval(e)) stop("object ", e, " not found.")
  }
  if (is.symbolic_function(e)) {
    ## create new symbolic expression
    return(.SymbolicExpression$new(
      symbolic_expr = paste(tag, parenthesize(e$expr), sep = ""),
      symbol_table = e$symbol_table
    ))
  }
  set_private(e, "expr_", paste(tag, parenthesize(e$expr), sep = ""))
  return(e)
}

binary_symbolic <- function(tag, expr1, expr2) {
  e1 <- eval(parse(text = expr1))
  e2 <- eval(parse(text = expr2))
  if(is.symbolic_expression(e1) && is.symbolic_expression(e2)) {
      return(
      .SymbolicExpression$new(
        symbolic_expr = paste(tag, "(", e1$expr, ",", e2$expr, ")", sep = ""),
        symbol_table = modifyList(e1$symbol_table, e2$symbol_table)
      )
    )
  } else { ## one between e1 and e2 is a value
    if (is.symbolic_expression(e1)) {
      v <- e2
      e <- e1
    } else {
      v <- e1
      e <- e2
    }
    value_name <- random_name()
    symbolic_expr <- ifelse(!is.symbolic_expression(e1),
        paste(tag, "(", value_name, ",", e$expr, ")", sep = ""),
        paste(tag, "(", e$expr, ",", value_name, ")", sep = "")
    )
    identifier <- list()
    identifier[[value_name]] <- v
    if (is.symbolic_function(e)) { ## create new symbolic expression
      return(.SymbolicExpression$new(
        symbolic_expr = symbolic_expr,
        symbol_table = modifyList(identifier, e$symbol_table)
      ))
    } else { ## modify symbolic in place
      set_private(e, "expr_", symbolic_expr)
      set_private(e, "symbol_table_", modifyList(identifier, e$symbol_table))
      return(e)
    }
  }
}

grad    <- function(op) unary_symbolic("grad"   , deparse(substitute(op)))
div     <- function(op) unary_symbolic("div"    , deparse(substitute(op)))
laplace <- function(op) unary_symbolic("laplace", deparse(substitute(op)))
dt      <- function(op) unary_symbolic("dt"     , deparse(substitute(op)))
inner   <- function(op1, op2) binary_symbolic("inner", deparse(substitute(op1)), deparse(substitute(op2)))

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
      cat("<SymbolicFunction>\n", sep = "")
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
