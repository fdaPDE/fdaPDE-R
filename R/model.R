simple_laplacian_penalty <- function(D, fe_order = 1L) {
    ## generate simple laplacian regularization
    L <- -laplace(Function(FunctionSpace(D, as.integer(fe_order))))
    u <- function(points) { ## homogeneous forcing term
        return(matrix(0, nrow = nrow(points), ncol = 1))
    }
    return(Pde(L, u))
}

## seed and n_mc_samples meaningfull only for stochastic edf_computation
gcv <- function(optimizer, lambda, edf_computation = c("stochastic", "exact"),
                seed = NULL, mc_samples = NULL, max_iter = NULL, step = NULL, tolerance = NULL) {
    gcv_params <- list(
        calibration_type = "gcv",
        optimizer = optimizer,
        edf_computation = match.arg(edf_computation),
        seed = seed,
        mc_samples = mc_samples,
        lambda = lambda,
        max_iter = max_iter,
        step = step,
        tolerance = tolerance
    )
    return(gcv_params)
}

is.nonparametric <- function(formula, data) {
    return(length(labels(terms(formula))) == 1 &&
        !(labels(terms(formula))[1] %in% colnames(data)))
}

is_in_list <- function(key, list) {
    if(!is.character(key)) stop("invalid key to list argument")
    return(key %in% names(list))
}

#' S3 dispatch to fit argument processing
#' @export
parse_eval_args <- function(x, ...) UseMethod("parse_eval_args")

#' @export
parse_eval_args.sr <- function(x, ...) {
    return(invisible(NULL))
}

#' @export
parse_eval_args.gsr <- function(x, ...) {
    args <- list(...)
    cpp_model <- args$cpp_model
    ## parse and set FPIRLS parameters, if provided
    if (is_in_list("fpirls_params", args)) {
        if (is_in_list("tolerance", args[["fpirls_params"]])) {
            cpp_model_$set_fpirls_tolerance(as.double(args[["fpirls_params"]]$tolerance))
        }
        if (is_in_list("max_iter", args[["fpirls_params"]])) {
            cpp_model_$set_fpirls_max_iter(args[["fpirls_params"]]$max_iter)
        }
    }
}

## utility R wrapper around C++ PDE concept
.RegressionModelCtr <- setRefClass(
    Class = "FdaPDEModel",
    fields = c(
        cpp_model  = "ANY", ## C++ model backend
        cpp_gcv    = "ANY", ## C++ gcv functor backend
        lambda_fixed = "logical",
        model_tag  = "ANY" ## model type tag
    ),
    methods = c(
        fit = function(...) {
            args <- list(...)
            ## customize fitting behaviour, depending on model type tag
            parse_eval_args(model_tag, cpp_model, args)
            ## set up calibration approach
            if (is_in_list("lambda", args)) {
                if (args[["lambda"]]$calibration_type == "gcv") {
                    ## require GCV-based selection of smoothing parameter
                    cpp_model$init()
                    cpp_gcv$calibrate(args[["lambda"]]) ## forward arguments
                }
            } else if (lambda_fixed) { ## fit with provided lambda
                cpp_model$init()
                cpp_model$solve()
            } else {
                ## default
                stop("using default")
            }
        },
        gcvs = function() { return(cpp_gcv$gcvs()) },
        edfs = function() { return(cpp_gcv$edfs()) }
    )
)

setGeneric("SRPDE", function(formula, data, domain, penalty, ...) { standardGeneric("SRPDE") })

.SRPDE <- function(formula, data, domain, penalty, lambda,
                   family = c("gaussian", "poisson", "exponential", "gamma", "bernulli")) {
    ## count number of terms in formula which are not covariates
    count <- 0
    field_name <- NULL
    for (term in labels(terms(formula))) {
        if (!term %in% colnames(data)) {
            count <- count + 1
            field_name <- term
        }
    }
    if (count == 0) stop("no spatial field found in formula.")
    if (count >= 2) stop("more than one spatial field supplied in formula.")
    ## if not supplied, default to simple laplacian penalty
    if (is.null(penalty)) {
        penalty <- simple_laplacian_penalty(domain)
        ## inject field in calling environment
        eval(parse(text = paste(field_name, "<<- Function(FunctionSpace(domain))", sep = "")))
    } else {
        ## check that name of the spatial field in formula coincide with name of field in penalty
        
    }
    ## extract data according to formula
    obs <- data[all.vars(formula)[1]]
    cov <- NULL
    if (!is.nonparametric(formula, data)) {
        covariates <- labels(terms(formula))
        covariates <- covariates[!covariates == field_name]
        cov <- as.matrix(data[covariates])
    }
    ## create model object
    model_ <- NULL ## C++ model backend
    model_tag <- list() ## model type tag
    family_ <- match.arg(family)
    model_ <- switch(family_,
        "gaussian" = {
            class(model_tag) <- "sr"  ## linear regression
            new(cpp_srpde, penalty, 0)
        },
        { ## default case
            class(model_tag) <- "gsr" ## generalized regression
            new(cpp_gsrpde_s, penalty, 0, family_)
        }
    )
    model_$set_observations(as.matrix(obs))

    if (!is.null(lambda)) model_$set_lambda_D(lambda)
    if (!is.null(cov)) model_$set_covariates(as.matrix(cov))

    return(.RegressionModelCtr(
        cpp_model = model_,
        model_tag = model_tag,
        lambda_fixed = !is.null(lambda),
        cpp_gcv = new(cpp_gcv, model_$get_gcv())
    ))
}

#' @export
setMethod(
    f = "SRPDE",
    signature = c(
        formula = "formula", data = "ANY", domain = "MeshObject", penalty = "missing"
    ),
    ## default to simple laplacian penalty
    definition = function(formula, data, domain, lambda = NULL, family = "gaussian") {
        .SRPDE(formula, data, domain, NULL, lambda, family)
    }
)

#' @export
setMethod(
    f = "SRPDE",
    signature = c(
        formula = "formula", data = "ANY", domain = "missing", penalty = "ANY"
    ),
    ## general penalty case (domain information indirectly supplied from penalty)
    definition = function(formula, data, penalty, lambda = NULL, family = "gaussian") {
        .SRPDE(formula, data, NULL, penalty, lambda, family)
    }
)

is.formula <- function(x){
   inherits(x, "formula")
}

GSRPDE <- function(model, domain, data = NULL, lambda = NULL,
                   family = c("poisson", "exponential", "gamma", "bernulli"),
                   fpirls_params = NULL) {
    penalty_ <- simple_laplacian_penalty(domain)

    obs <- NULL
    cov <- NULL
    if (is.formula(model)) { ## parametric model (covariates are supplied)
        if (is.null(data)) stop("for parametric models, you must supply a valid data frame")
        obs <- data[all.vars(model)[1]]
        cov <- as.matrix(data[labels(terms(model))])
    } else { ## non-parametric model
        if (ncol(as.matrix(model)) != 1) stop("non parametric models require a vector of observations")
        obs <- as.matrix(model)[, 1]
    }
    ## create model object
    model_ <- new(R_GSRPDE, penalty_, 0, match.arg(family))
    model_$set_observations(as.matrix(obs))

    if (!is.null(fpirls_params)) {
        if ("tolerance" %in% fpirls_params) model_$set_fpirls_tolerance(as.double(fpirls_params$tolerance))
        if ("max_iter"  %in% fpirls_params) model_$set_fpirls_max_iter (fpirls_params$max_iter)
    }
    if (!is.null(lambda)) model_$set_lambda_D(lambda)
    if (!is.null(cov))    model_$set_covariates(as.matrix(cov))

    model_$init()
    return(model_)
}

STRPDE <- function(model, domain, data = NULL, lambda = NULL) {
    penalty_ <- simple_laplacian_penalty(domain)

    obs <- NULL
    cov <- NULL
    if (is.formula(model)) { ## parametric model (covariates are supplied)
        if (is.null(data)) stop("for parametric models, you must supply a valid data frame")
        obs <- data[all.vars(model)[1]]
        cov <- as.matrix(data[labels(terms(model))])
    } else { ## non-parametric model
        if (ncol(as.matrix(model)) != 1) stop("non parametric models require a vector of observations")
        obs <- as.matrix(model)[, 1]
    }
    ## create model object
    model_ <- new(R_STRPDE, penalty_, domain$times, 0)
    model_$set_observations(as.matrix(obs))

    if (!is.null(lambda)) {
        model_$set_lambda_D(lambda[1])
        model_$set_lambda_T(lambda[2])
    }
    if (!is.null(cov))      model_$set_covariates(as.matrix(cov))

    model_$init()
    return(model_)    
}
