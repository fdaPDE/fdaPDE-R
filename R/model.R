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

## colored output
colored_text <- function(color, text) cat(paste0("\033[0;", color, "m", text, "\033[0m"))
red <- function(text) colored_text(31, text)
green <- function(text) colored_text(32, text)
yellow <- function(text) colored_text(33, text)
blue <- function(text) colored_text(34, text)

lhs <- function(formula) setdiff(all.vars(formula), labels(terms(formula)))
rhs <- function(formula) labels(terms(formula))

#' @export
gcv <- function(lambda = NULL, optimizer = c("grid", "newton", "gd", "bfgs"),
                edf_evaluation = c("stochastic", "exact"), seed = NULL, n_mc_samples = 100, max_iter = 20,
                step = 1e-2, tolerance = 1e-4) {
  out <- list(
    calibration = "gcv",
    lambda = as.matrix(lambda),
    ## edf computation settings
    edf = match.arg(edf_evaluation),
    n_mc_samples = n_mc_samples,
    seed = if(is.null(seed)) -1 else seed,
    ## optimization settings
    optimizer = match.arg(optimizer),
    max_iter = max_iter,
    tolerance = tolerance,
    step = step
  )
  return(out)
}

model_parse_args <- function(x, ...) UseMethod("model_parse_args", x)
model_parse_args.SRPDE <- function(x, args) {
  return(list())
}
model_parse_args.GSRPDE <- function(x, args) {
  out <- list()
  ## defaults
  out$fpirls_tolerance <- 1e-4
  out$fprils_max_iter <- 200
  if ("fpirls_params" %in% names(args)) {
    out$fpirls_tolerance <- if ("tolerance" %in% args[["fpirls_params"]]) args[["fpirls_params"]]$tolerance
    out$fpirls_max_iter <- if ("max_iter" %in% args[["fpirls_params"]]) args[["fpirls_params"]]$max_iter
  }
  return(out)
}

.SpatialRegression <- R6::R6Class(
  "SpatialRegressionModel",
  private = list(
    model_   = NULL,
    formula_ = formula(),
    penalty_ = NULL,
    family_  = character()
  ),
  public = list(
    initialize = function(formula, data, penalty = NULL, family = NULL) {
      private$formula_ <- formula
      private$family_ <- family
      field_name <- NULL
      for (term in labels(terms(formula))) {
        if (!term %in% colnames(data)) {
          if (!is.null(field_name)) stop("ill-formed formula, missing field.")
          field_name <- term
        }
      }
      if (is.null(penalty)) {
        ## default to simple laplacian penalty, homogeneous neumann and dirichlet BCs, null forcing term
        f <- eval(parse(text = as.character(field_name)))
        if (!inherits(f, "SymbolicFunction")) {
          stop(as.charcter(field_name), ": not a SymbolicFunction object.")
        }
        Lf <- -laplace(f, field_name)
        u <- function(nodes) matrix(rep(0, times = nrow(nodes)), ncol = 1)
        penalty <- PDE(Lf, u)
      }
      private$penalty_ <- penalty
      ## TODO: need to derive type of sampling
      sampling <- 0

      ## instantiate model cpp backend
      if (family == "gaussian") {
        private$model_ <- new(cpp_srpde, get_private(penalty)$pde_, sampling)
      } else if (family %in% names(distributions)) {
        private$model_ <- new(
          cpp_gsrpde_space, get_private(penalty)$pde_, sampling, cpp_aligned_index(distributions[[family]])
        )
      }
      ## analyze formula
      if (length(lhs(formula)) == 0) stop("ill-formed formula, missing lhs")
      private$model_$set_observations(as.matrix(data[lhs(formula)]))
      ## check if model has covariates
      if (length(setdiff(all.vars(formula), c(lhs(formula), field_name))) != 0) { 
        covnames <- setdiff(labels(terms(formula)), field_name)
        private$model_$set_covariates(as.matrix(data[covnames]))
      }
    },
    fit = function(...) {
      args <- list(...)
      fit_data <- list()
      ## select calibration type
      if (is.numeric(args$calibration)) {
        fit_data$calibration <- "off"
        fit_data$lambda <- as.matrix(args$calibration)
      } else {
        if (is.list(args$calibration)) {
          fit_data <- c(fit_data, args$calibration)
        } else {
          stop("invalid type bounded to calibration argument.")
        }
      }
      fit_data <- modifyList(model_parse_args(self, args), fit_data) ## S3 dispatch to model specific fit logic
      private$model_$fit(fit_data)
    },
    print = function(...) {
      cat("Spatial regression model:", deparse(private$formula_), "\n")
      cat("PDE:", private$penalty_$operator, "\n")
      ## S3 dispatch for specific printing logic
    }
  ),
  ## active bindings
  active = list(
    f = function() as.matrix(private$model_$f()),
    fitted = function() as.matrix(private$model_$fitted()),
    gcvs = function() as.matrix(private$model_$gcvs()),
    edfs = function() as.matrix(private$model_$edfs()),
    optimal_lambda = function() as.matrix(private$model_$optimum())
  )
)

## these indexes, once cpp aligned, must match with the model instantiation at cpp side
distributions <- list(poisson = 1, bernulli = 2, exponential = 3, gamma = 4)

## space-only regression model family
SRPDE <- function(formula, data, penalty = NULL,
                  family = c("gaussian", "poisson", "bernulli", "exponential", "gamma", "quantile")) {
  family <- match.arg(family)
  model <- .SpatialRegression$new(formula, data, penalty, family)
  if (family == "gaussian") {
      class(model) <- append("SRPDE", class(model))
  }
  if (family %in% names(distributions)) {
    class(model) <- append("GSRPDE", class(model))
  }
  return(model)
}


## #' @export
## parse_eval_args.qsr <- function(x, ...) {
##     args <- list(...)
##     cpp_model_ <- args$cpp_model
##     ## parse and set FPIRLS parameters, if provided
##     parse_fpirls_args(cpp_model_, args)
##     ## set quantile level
##     if(is_in_list("alpha", args)) {
##         alpha = as.double(args[["alpha"]])
##         if ( alpha <= 0 || alpha > 1 ) stop("quantile level alpha should be between 0 and 1")
##         cpp_model_$set_alpha(as.double(args[["alpha"]]))
##     }
## }
