simple_laplacian_penalty <- function(D, fe_order = 1L) {
    ## generate simple laplacian regularization
    L <- -laplace(Function(FunctionSpace(D, as.integer(fe_order))))
    u <- function(points) { ## homogeneous forcing term
        return(matrix(0, nrow = nrow(points), ncol = 1))
    }
    return(PDE(L, u))
}

## seed and n_mc_samples meaningfull only for stochastic edf_computation
gcv <- function(optimizer, lambda, edf_computation = c("stochastic", "exact"),
                seed = NULL, mc_samples = NULL, max_iter = NULL, step = NULL, tolerance = NULL) {
    gcv_params <- list(
        calibration_type = "gcv",
        ## properties related to the computation of the expected degrees of freedom
        edf_computation = match.arg(edf_computation),
        seed = seed,
        mc_samples = mc_samples,
        ## properties related to the optimization algorithm
        optimizer = optimizer,
        lambda = lambda,
        max_iter = max_iter,
        step = step,
        tolerance = tolerance
    )
    return(gcv_params)
}

is_nonparametric <- function(formula, data) {
    return(length(labels(terms(formula))) == 1 &&
        !(labels(terms(formula))[1] %in% colnames(data)))
}

is_in_list <- function(key, list) {
    if(!is.character(key)) stop("invalid key to list argument")
    return(key %in% names(list))
}

parse_fpirls_args <- function(cpp_model, args) {
    if (is_in_list("fpirls_params", args)) {
        if (is_in_list("tolerance", args[["fpirls_params"]])) {
            cpp_model$set_fpirls_tolerance(as.double(args[["fpirls_params"]]$tolerance))
        }
        if (is_in_list("max_iter", args[["fpirls_params"]])) {
            cpp_model$set_fpirls_max_iter(as.integer(args[["fpirls_params"]]$max_iter))
        }
    }
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
    cpp_model_ <- args$cpp_model
    ## parse and set FPIRLS parameters, if provided
    parse_fpirls_args(cpp_model_, args)
}

#' @export
parse_eval_args.qsr <- function(x, ...) {
    args <- list(...)
    cpp_model_ <- args$cpp_model
    ## parse and set FPIRLS parameters, if provided
    parse_fpirls_args(cpp_model_, args)
    ## set quantile level
    if(is_in_list("alpha", args)) {
        alpha = as.double(args[["alpha"]])
        if ( alpha <= 0 || alpha > 1 ) stop("quantile level alpha should be between 0 and 1")
        cpp_model_$set_alpha(as.double(args[["alpha"]]))
    }
}

## utility R wrapper around C++ PDE concept
.RegressionModelCtr <- setRefClass(
    Class = "FdaPDEModel",
    fields = c(
        cpp_model  = "ANY", ## C++ model backend
        cpp_gcv    = "ANY", ## C++ gcv functor backend
        lambda_fixed = "logical",
        spatial_field = "FunctionObject",
        model_tag  = "ANY"  ## model type tag
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

                    ## final fit with optimal lambda
                    opt_lambda = args[["lambda"]]$lambda[which.min(cpp_gcv$gcvs())]
                    cpp_model$set_lambda_D(opt_lambda)
                    cpp_model$init()
                    cpp_model$solve()
                }
            } else if (lambda_fixed) { ## fit with provided lambda
                cpp_model$init()
                cpp_model$solve()
            } else {
                ## default
                stop("using default")
            }
            spatial_field$coeff = as.matrix(cpp_model$f())
        },
        gcvs = function() { return(cpp_gcv$gcvs()) },
        edfs = function() { return(cpp_gcv$edfs()) }
    )
)

setGeneric("SRPDE", function(formula, data, domain, penalty, ...) { standardGeneric("SRPDE") })

.SRPDE <- function(formula, data, domain, penalty, lambda,
                   family = c("gaussian", "poisson", "exponential", "gamma", "bernulli", "quantile")) {
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
        field <- Function(FunctionSpace(domain))
        eval(parse(text = paste(field_name, "<<- field", sep = "")))
    } else {
        ## TODO: check that name of the spatial field in formula coincide with name of field in penalty
        
    }
    ## extract data according to formula
    obs <- data[all.vars(formula)[1]]
    cov <- NULL
    if (!is_nonparametric(formula, data)) {
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
            class(model_tag) <- "sr"
            new(cpp_srpde, penalty, 0)
        },
        "quantile" = {
            class(model_tag) <- "qsr"
            new(cpp_qsrpde_s, penalty, 0)
        },
        { ## default case
            class(model_tag) <- "gsr"
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
        cpp_gcv = new(cpp_gcv, model_$get_gcv()),
        spatial_field = field
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
    ## general penalty (domain information indirectly supplied from penalty)
    definition = function(formula, data, penalty, lambda = NULL, family = "gaussian") {
        .SRPDE(formula, data, NULL, penalty, lambda, family)
    }
)
