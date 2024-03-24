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

## .RegressionModel <- R6::R6Class(
##   "RegressionModel",
##   private = list(
##       model_ = NULL,
##       penalty_ = character(),
##       tag_ = list() ## can i do S3 dispatching on myself?? see below
##     ## calibrator
##   ),
##   public = list(
##       fit = function(lambda, ...) {
##           if()
##               args <- list(...)
##           ## check for calibration, if absent fallback to default
##           ## parse model-specific arguments (S3 dispatch to model processing)

##       }
##   ),
##   ## active bindings
##   active = list(
##     ## results
##     ## covariate/response
##   )
## )

## #' S3 dispatch to fit argument processing
## #' @export
## call_dispatch <- function(x, ...) UseMethod("call_dispatch")

## #' @export
## call_dispatch.mortis <- function(x, ...) {
##   print("mortis")
## }

## #' @export
## call_dispatch.fango <- function(x, ...) {
##   print("fango")
## }

## #' @export
## .Base <- R6::R6Class(
##   "Base",
##   private = list(
##     name = character()
##   ),
##   public = list(
##     initialize = function(name = NULL) {
##       private$name <- name
##     },
##     dispatch = function() {
##       self$derived_function()
##       # call_dispatch(self)
##     }
##   )
## )

## ## base model only provides a fit method, which performs basic initialization +
## ## discharge pack argument parse to derived instance + orchestrate calibration

## ## essentially fit recives the pack and consumes its arguments, customizing the model
## ## behaviour. finally it actually fits the model

## ## derived expose getters to its solution??
## ## we have to expose an R6 class for each model (is this necessary??) what about S3 dispatch??

## .Derived <- R6::R6Class(
##   "Derived",
##   inherit = .Base,
##   public = list(
##     derived_function = function() {
##       print("hasan")
##     }
##   )
## )

## use composition

is_nonparametric <- function(formula, data) {
  return(length(labels(terms(formula))) == 1 && !(labels(terms(formula))[1] %in% colnames(data)))
}

lhs <- function(formula) setdiff(all.vars(formula), labels(terms(formula)))
rhs <- function(formula) labels(terms(formula))

.SRPDE <- R6::R6Class(
  "LinearSpatialRegression",
  private = list(
    model_ = NULL,
    formula_ = formula()
  ),
  public = list(
    initialize = function(formula, data, penalty = NULL) {
      if (is.null(penalty)) {
        ## default to simple laplacian penalty, homogeneous neumann and dirichlet BCs, null forcing term
        field <- NULL
        for (term in labels(terms(formula))) {
          if (!term %in% colnames(data)) {
            if (!is.null(field)) stop("ill-formed formula, missing spatial field")
            field <- term
          }
        }
        f <- eval(parse(text = as.character(field)))
        Lf <- -laplace(f)
        penalty <- PDE(Lf) ## poisson problem over domain
      }
      ## TODO: need to derive type of sampling
      sampling <- 0
      ## ------------- do not depend on regression model type

      
      ## instantiate model cpp backend
      model_ <- new(cpp_srpde, penalty, sampling)

      ## analyze formula and prepare data
      if (length(lhs(formula)) == 0) stop("ill-formed formula, missing lhs")
      obs <- data[lhs(formula)]
      cov <- matrix()
      if (length(setdiff(all.vars(formula), c(lhs(formula), field))) != 0) { ## semiparametric model
        covnames <- setdiff(labels(terms(formula)), field)
        cov <- as.matrix(data[covnames])
      }
      ## set data to model
      model_$set_observations(obs)
      model_$set_covariates(cov)
      
    },
    fit = function(lambda = NULL, ...) {
        ## analyze pack arguments to customize model fitting behaviour

        ## analyze lambda type, if fixed impose otherwise
        
        ## if gcv, optimize via gcv calibration (add active bindings to get results? is this ok?)
        ## if kcv, optimize via kcv
        ## impose optimal lambda and refit
    }
  ),
  ## active bindings
  active = list(
      ## getters
  )
)

## the use of two parameters lambda/calibration in fit() is ambigous for various reasons
## lambda cannot be inserted outside calibration, since the meaning of lambda is different depending on the type of optimization
## suppose I say model$fit(lambda = 1e-2), this would mean a fixed lambda optimization strategy, or a defaulted gcv strategy
## (which uses newton optimization) where lambda is the starting point?? ambiguous

## if we remove lambda parameter and use only calibration, what about a fixed lambda?? doesn't seems good to say
## calibration = "a number", as by calibration we usually mean the strategy/act to select lambda.

## we take a single parameter, named lambda (rho??), which can be a number (it is semantically fine to say lambda = "a number"),
## or one between gcv and kcv, with the meaning: select the optimal lambda using a gcv/kcv strategy

## discard the possibility to say calibration = off(lambda = ), as would result in an heavy interface. If I have a lambda
## to impose, just give that number/vector, without need for extra code.

## we remove the lambda parameter from the model signature, as having lambda both in the model signature and in the fit method is
## again ambigous and redundant.

## space-only regression model family
SRPDE <- function(formula, data, penalty = NULL, family = NULL) {
  ## based on family, instantiate correct module providing penalty

  ## if penalty is null, default to simple laplacian

  ## set data
  ## set locations (how to handle areal data?? to ask)
}

## model is eventually inited at fit time (model inits only what it needs)


## need for something more than a simple data-frame, we can have covariates recorded at different locations


## simple_laplacian_penalty <- function(D, fe_order = 1L) {
##     ## generate simple laplacian regularization
##     L <- -laplace(Function(FunctionSpace(D, as.integer(fe_order))))
##     u <- function(points) { ## homogeneous forcing term
##         return(matrix(0, nrow = nrow(points), ncol = 1))
##     }
##     return(PDE(L, u))
## }

## ## seed and n_mc_samples meaningfull only for stochastic edf_computation
## gcv <- function(optimizer, lambda, edf_computation = c("stochastic", "exact"),
##                 seed = NULL, mc_samples = NULL, max_iter = NULL, step = NULL, tolerance = NULL) {
##     gcv_params <- list(
##         calibration_type = "gcv",
##         ## properties related to the computation of the expected degrees of freedom
##         edf_computation = match.arg(edf_computation),
##         seed = seed,
##         mc_samples = mc_samples,
##         ## properties related to the optimization algorithm
##         optimizer = optimizer,
##         lambda = lambda,
##         max_iter = max_iter,
##         step = step,
##         tolerance = tolerance
##     )
##     return(gcv_params)
## }

## is_nonparametric <- function(formula, data) {
##     return(length(labels(terms(formula))) == 1 &&
##         !(labels(terms(formula))[1] %in% colnames(data)))
## }

## is_in_list <- function(key, list) {
##     if(!is.character(key)) stop("invalid key to list argument")
##     return(key %in% names(list))
1 ## }

## parse_fpirls_args <- function(cpp_model, args) {
##     if (is_in_list("fpirls_params", args)) {
##         if (is_in_list("tolerance", args[["fpirls_params"]])) {
##             cpp_model$set_fpirls_tolerance(as.double(args[["fpirls_params"]]$tolerance))
##         }
##         if (is_in_list("max_iter", args[["fpirls_params"]])) {
##             cpp_model$set_fpirls_max_iter(as.integer(args[["fpirls_params"]]$max_iter))
##         }
##     }
## }

## #' S3 dispatch to fit argument processing
## #' @export
## parse_eval_args <- function(x, ...) UseMethod("parse_eval_args")

## #' @export
## parse_eval_args.sr <- function(x, ...) {
##     return(invisible(NULL))
## }

## #' @export
## parse_eval_args.gsr <- function(x, ...) {
##     args <- list(...)
##     cpp_model_ <- args$cpp_model
##     ## parse and set FPIRLS parameters, if provided
##     parse_fpirls_args(cpp_model_, args)
## }

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

## ## utility R wrapper around C++ PDE concept
## .RegressionModelCtr <- setRefClass(
##     Class = "FdaPDEModel",
##     fields = c(
##         cpp_model  = "ANY", ## C++ model backend
##         cpp_gcv    = "ANY", ## C++ gcv functor backend
##         lambda_fixed = "logical",
##         spatial_field = "FunctionObject",
##         model_tag  = "ANY"  ## model type tag
##     ),
##     methods = c(
##         fit = function(...) {
##             args <- list(...)
##             ## customize fitting behaviour, depending on model type tag
##             parse_eval_args(model_tag, cpp_model, args)
##             ## set up calibration approach
##             if (is_in_list("lambda", args)) {
##                 if (args[["lambda"]]$calibration_type == "gcv") {
##                     ## require GCV-based selection of smoothing parameter
##                     cpp_model$init()
##                     cpp_gcv$calibrate(args[["lambda"]]) ## forward arguments

##                     ## final fit with optimal lambda
##                     opt_lambda = args[["lambda"]]$lambda[which.min(cpp_gcv$gcvs())]
##                     cpp_model$set_lambda_D(opt_lambda)
##                     cpp_model$init()
##                     cpp_model$solve()
##                 }
##             } else if (lambda_fixed) { ## fit with provided lambda
##                 cpp_model$init()
##                 cpp_model$solve()
##             } else {
##                 ## default
##                 stop("using default")
##             }
##             spatial_field$coeff = as.matrix(cpp_model$f())
##         },
##         gcvs = function() { return(cpp_gcv$gcvs()) },
##         edfs = function() { return(cpp_gcv$edfs()) }
##     )
## )

## setGeneric("SRPDE", function(formula, data, domain, penalty, ...) { standardGeneric("SRPDE") })

## .SRPDE <- function(formula, data, domain, penalty, lambda,
##                    family = c("gaussian", "poisson", "exponential", "gamma", "bernulli", "quantile")) {
##     ## count number of terms in formula which are not covariates
##     count <- 0
##     field_name <- NULL
##     for (term in labels(terms(formula))) {
##         if (!term %in% colnames(data)) {
##             count <- count + 1
##             field_name <- term
##         }
##     }
##     if (count == 0) stop("no spatial field found in formula.")
##     if (count >= 2) stop("more than one spatial field supplied in formula.")
##     ## if not supplied, default to simple laplacian penalty
##     if (is.null(penalty)) {
##         penalty <- simple_laplacian_penalty(domain)
##         ## inject field in calling environment
##         field <- Function(FunctionSpace(domain))
##         eval(parse(text = paste(field_name, "<<- field", sep = "")))
##     } else {
##         ## TODO: check that name of the spatial field in formula coincide with name of field in penalty

##     }
##     ## extract data according to formula
##     obs <- data[all.vars(formula)[1]]
##     cov <- NULL
##     if (!is_nonparametric(formula, data)) {
##         covariates <- labels(terms(formula))
##         covariates <- covariates[!covariates == field_name]
##         cov <- as.matrix(data[covariates])
##     }
##     ## create model object
##     model_ <- NULL ## C++ model backend
##     model_tag <- list() ## model type tag
##     family_ <- match.arg(family)
##     model_ <- switch(family_,
##         "gaussian" = {
##             class(model_tag) <- "sr"
##             new(cpp_srpde, penalty, 0)
##         },
##         "quantile" = {
##             class(model_tag) <- "qsr"
##             new(cpp_qsrpde_s, penalty, 0)
##         },
##         { ## default case
##             class(model_tag) <- "gsr"
##             new(cpp_gsrpde_s, penalty, 0, family_)
##         }
##     )
##     model_$set_observations(as.matrix(obs))

##     if (!is.null(lambda)) model_$set_lambda_D(lambda)
##     if (!is.null(cov)) model_$set_covariates(as.matrix(cov))

##     return(.RegressionModelCtr(
##         cpp_model = model_,
##         model_tag = model_tag,
##         lambda_fixed = !is.null(lambda),
##         cpp_gcv = new(cpp_gcv, model_$get_gcv()),
##         spatial_field = field
##     ))
## }

## #' @export
## setMethod(
##     f = "SRPDE",
##     signature = c(
##         formula = "formula", data = "ANY", domain = "MeshObject", penalty = "missing"
##     ),
##     ## default to simple laplacian penalty
##     definition = function(formula, data, domain, lambda = NULL, family = "gaussian") {
##         .SRPDE(formula, data, domain, NULL, lambda, family)
##     }
## )

## #' @export
## setMethod(
##     f = "SRPDE",
##     signature = c(
##         formula = "formula", data = "ANY", domain = "missing", penalty = "ANY"
##     ),
##     ## general penalty (domain information indirectly supplied from penalty)
##     definition = function(formula, data, penalty, lambda = NULL, family = "gaussian") {
##         .SRPDE(formula, data, NULL, penalty, lambda, family)
##     }
## )
