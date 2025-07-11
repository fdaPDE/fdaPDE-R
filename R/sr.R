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

# An R6 class encapsulating the spatial regression method with partial differential equation regularization.
#
#' @rdname sr
# @order 2
.sr <- R6::R6Class(
  "sr",
  private = list(
    model_ = NULL, ## cpp backend
    f_ = "fe_function"
  ),
  public = list(
    #' @description
    #' Creates a new \code{sr} object.
    #'
    #' @param formula R formula.
    #' @param data A \code{geoframe} containing the triangulation of the domain and the data (see also [geoframe()]).
    #' @param penalty A penalty term returned by [fe_elliptic()].
    initialize = function(formula, data, penalty) {
      ## recover name of non-parametric field
      vars <- all.vars(as.formula(formula)[[3]])
      f_symbol <- setdiff(vars, names(data))
      fdapde_assert(length(f_symbol) == 1, "Expected exactly one nonparametric term.")
      private$f_ <- get(f_symbol, envir = globalenv())
      fe_type <- get_private(private$f_)$type

      if (is.null(penalty)) {
        ## fallback to laplacian penalty
        private$model_ = new(cpp_sr_2_2, formula, get_private(data$gf__ptr__)$ptr_, penalty)
      } else {
        domain <- get_private(get_private(data$gf__ptr__)$mesh_)$mesh_
        params <- list()
        quad_nodes <- matrix()
        if (is.function(penalty$K) || is.function(penalty$b) || is.function(penalty$c) || is.function(penalty$u)) {
          quad_nodes <- fe_simplex_quad_nodes(domain, fe_type)
        }
        n_quad_nodes <- new(cpp_fe_space_2_2_p1, domain)$n_quad_nodes()
        embed_dim <- 2

        ## diffusion
        if (is.null(penalty$K)) {
          params$K = matrix(rep(0, times = (embed_dim * embed_dim) * n_quad_nodes), nrow = n_quad_nodes)
        } else {
          fdapde_assert(is.matrix(K) || is.function(K), "Invalid diffusion field type.")
          if (is.matrix(K)) {
            fdapde_assert(nrow(K) == embed_dim && ncol(K) == embed_dim, "Not a square matrix.")
            params$K = matrix(rep(c(K), times = n_quad_nodes), nrow = n_quad_nodes, byrow = TRUE)
          } else {
            q <- quad_nodes[1, ]
            fdapde_assert(
              nrow(K(q)) == embed_dim && ncol(K(q)) == embed_dim,
              "Not evaluates to a square matrix."
            )
            params$K = K(quad_nodes)
          }
        }
        ## transport
        if (is.null(penalty$b)) {
          params$b = matrix(rep(0, times = embed_dim * n_quad_nodes), nrow = n_quad_nodes)
        } else {
          fdapde_assert(is.matrix(b) || is.vector(b) || is.function(b), "Invalid transport field type.")
          if (is.matrix(b) || is.vector(b)) {
            fdapde_assert(nrow(b) == embed_dim && ncol(b) == 1, "Not a vector.")
            params$b = matrix(rep(c(b), times = n_quad_nodes), nrow = n_quad_nodes, byrow = TRUE)
          } else {
            q <- quad_nodes[1, ]
            fdapde_assert(nrow(b(q)) == embed_dim && ncol(b(q)) == 1, "Not evaluates to a vector.")
            params$b = b(quad_nodes)
          }
        }
        ## reaction
        if (is.null(penalty$c)) {
          params$c = matrix(rep(0, times = n_quad_nodes), nrow = n_quad_nodes)
        } else {
          fdapde_assert(is.numeric(c) || is.function(c), "Invalid reaction field type.")
          if (is.numeric(c)) {
            params$c = matrix(rep(c, times = n_quad_nodes), nrow = n_quad_nodes)
          } else {
            q <- quad_nodes[1, ]
            fdapde_assert(is.numeric(c(q)), "Not evaluates to a scalar.")
            params$c = c(quad_nodes)
          }
        }
        ## force
        if (is.null(penalty$u)) {
          params$u = matrix(rep(0, times = n_quad_nodes), nrow = n_quad_nodes)
        } else {
          fdapde_assert(is.numeric(u) || is.function(u), "Invalid forcing field type.")
          if (is.numeric(u)) {
            params$u = matrix(rep(u, times = n_quad_nodes), nrow = n_quad_nodes)
          } else {
            q <- quad_nodes[1, ]
            fdapde_assert(is.numeric(u(q)), "Not evaluates to a scalar.")
            params$c = c(quad_nodes)
          }
        }
        private$model_ = new(cpp_sr_2_2, formula, get_private(data$gf__ptr__)$ptr_, params)
      }
    },
    #' @description
    #' Fits the statistical model.
    #'
    #' @param lambda The smoothing parameter. Default is \code{NULL}.
    #' @param calibrator A calibrator object. Available calibrators include [gcv()].
    fit = function(lambda = NULL, calibrator = NULL) {
      fdapde_assert(!is.null(calibrator) || !is.null(lambda), "Unable to select smoothing level.")
      if (is.null(calibrator)) {
        private$model_$fit(lambda)
      } else {
        r <- private$model_$fit_gcv(as.character(calibrator$opt_t), calibrator)
      }
      ## update estimated field coefficient vector
      get_private(private$f_)$fe_function_$set_coeff(private$model_$f())
      if (!is.null(calibrator)) return(r)
    }
  ),
  active = list(
    #' @field f The estimated nonparametric component of the model.
    f = function() private$model_$f(),
    #' @field beta The estimated parametric component of the model.
    beta = function() private$model_$beta(),
    #' @field fitted The vector of values fitted by the model.
    fitted = function() private$model_$fitted()
  )
)

#' Create an \code{sr} object
#'
#' @param formula A formula object describing the relationship between the response variable and the model components (nonparametric or semiparametric).
#' @param data A \code{geoframe} containing both the triangulation of the domain and the associated data (see also [geoframe()]).
#' @param penalty A penalty term returned by [fe_elliptic()]. Default is \code{NULL}.
#' @rdname sr
#' @order 1
#' @references
#' \itemize{
#'    \item Sangalli, L. M., Ramsay, J. O., Ramsay, T. O. (2013). Spatial spline regression models.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 75(4), 681-703.
#'    \item Azzimonti, L., Sangalli, L. M., Secchi, P., Domanin, M., Nobile, F. (2015). Blood flow velocity field estimation
#' via spatial regression with PDE penalization. Journal of the American Statistical Association, 110(511), 1057-1071.
#' }
#' @export
sr <- function(formula, data, penalty = NULL) {
  return(.sr$new(
    formula = deparse(formula),
    data = data,
    penalty = penalty
  ))
}

# An R6 class encapsulating the generalized spatial regression method with partial differential equation regularization.
#
#' @rdname gsr
#' @order 2
.gsr <- R6::R6Class(
  "gsr",
  private = list(
    model_ = NULL, ## cpp backend
    f_ = "fe_function"
  ),
  public = list(
    #' @description
    #' Creates a new \code{gsr} object.
    #'
    #' @param formula R formula.
    #' @param data A \code{geoframe} containing the triangulation of the domain and the data (see also [geoframe()]).
    #' @param family A string specifying the exponential family to be used.
    #'        Avaialble option are \code{"bernoulli"}, \code{"poisson"}, \code{"exponential"} or \code{"gamma"}.
    #' @param penalty A penalty term returned by [fe_elliptic()].
    initialize = function(formula, data, family, penalty) {
      ## recover name of non-parametric field
      vars <- all.vars(as.formula(formula)[[3]])
      f_symbol <- setdiff(vars, names(data))
      fdapde_assert(length(f_symbol) == 1, "Expected exactly one nonparametric term.")
      private$f_ <- get(f_symbol, envir = globalenv())
      fe_type <- get_private(private$f_)$type

      if (is.null(penalty)) {
        ## fallback to laplacian penalty
        private$model_ = new(cpp_gsr_2_2, formula, get_private(data$gf__ptr__)$ptr_, family, penalty)
      } else {
        domain <- get_private(get_private(data$gf__ptr__)$mesh_)$mesh_
        params <- list()
        quad_nodes <- matrix()
        if (is.function(penalty$K) || is.function(penalty$b) || is.function(penalty$c) || is.function(penalty$u)) {
          quad_nodes <- fe_simplex_quad_nodes(domain, fe_type)
        }
        n_quad_nodes <- new(cpp_fe_space_2_2_p1, domain)$n_quad_nodes()
        embed_dim <- 2

        ## diffusion
        if (is.null(penalty$K)) {
          params$K = matrix(rep(0, times = (embed_dim * embed_dim) * n_quad_nodes), nrow = n_quad_nodes)
        } else {
          fdapde_assert(is.matrix(K) || is.function(K), "Invalid diffusion field type.")
          if (is.matrix(K)) {
            fdapde_assert(nrow(K) == embed_dim && ncol(K) == embed_dim, "Not a square matrix.")
            params$K = matrix(rep(c(K), times = n_quad_nodes), nrow = n_quad_nodes, byrow = TRUE)
          } else {
            q <- quad_nodes[1, ]
            fdapde_assert(
              nrow(K(q)) == embed_dim && ncol(K(q)) == embed_dim,
              "Not evaluates to a square matrix."
            )
            params$K = K(quad_nodes)
          }
        }
        ## transport
        if (is.null(penalty$b)) {
          params$b = matrix(rep(0, times = embed_dim * n_quad_nodes), nrow = n_quad_nodes)
        } else {
          fdapde_assert(is.matrix(b) || is.vector(b) || is.function(b), "Invalid transport field type.")
          if (is.matrix(b) || is.vector(b)) {
            fdapde_assert(nrow(b) == embed_dim && ncol(b) == 1, "Not a vector.")
            params$b = matrix(rep(c(b), times = n_quad_nodes), nrow = n_quad_nodes, byrow = TRUE)
          } else {
            q <- quad_nodes[1, ]
            fdapde_assert(nrow(b(q)) == embed_dim && ncol(b(q)) == 1, "Not evaluates to a vector.")
            params$b = b(quad_nodes)
          }
        }
        ## reaction
        if (is.null(penalty$c)) {
          params$c = matrix(rep(0, times = n_quad_nodes), nrow = n_quad_nodes)
        } else {
          fdapde_assert(is.numeric(c) || is.function(c), "Invalid reaction field type.")
          if (is.numeric(c)) {
            params$c = matrix(rep(c, times = n_quad_nodes), nrow = n_quad_nodes)
          } else {
            q <- quad_nodes[1, ]
            fdapde_assert(is.numeric(c(q)), "Not evaluates to a scalar.")
            params$c = c(quad_nodes)
          }
        }
        ## force
        if (is.null(penalty$u)) {
          params$u = matrix(rep(0, times = n_quad_nodes), nrow = n_quad_nodes)
        } else {
          fdapde_assert(is.numeric(u) || is.function(u), "Invalid forcing field type.")
          if (is.numeric(u)) {
            params$u = matrix(rep(u, times = n_quad_nodes), nrow = n_quad_nodes)
          } else {
            q <- quad_nodes[1, ]
            fdapde_assert(is.numeric(u(q)), "Not evaluates to a scalar.")
            params$c = c(quad_nodes)
          }
        }
        private$model_ = new(cpp_gsr_2_2, formula, get_private(data$gf__ptr__)$ptr_, family, params)
      }
    },
    #' @description
    #' Fits the statistical model.
    #'
    #' @param lambda The smoothing parameter. Default is \code{NULL}.
    #' @param calibrator A calibrator object. Available calibrators include [gcv()].
    fit = function(lambda = NULL, calibrator = NULL) {
      fdapde_assert(!is.null(calibrator) || !is.null(lambda), "Unable to select smoothing level.")
      if (is.null(calibrator)) {
        private$model_$fit(lambda)
      } else {
        r <- private$model_$fit_gcv(as.character(calibrator$opt_t), calibrator)
      }
      ## update estimated field coefficient vector
      get_private(private$f_)$fe_function_$set_coeff(private$model_$f())
      if (!is.null(calibrator)) return(r)
    }
  ),
  active = list(
    #' @field f The estimated nonparametric component of the model.
    f = function() private$model_$f(),
    #' @field beta The estimated parametric component of the model.
    beta = function() private$model_$beta(),
    #' @field fitted The vector of values fitted by the model.
    fitted = function() private$model_$fitted()
  )
)

#' Create a \code{gsr} object
#'
#' @inheritParams sr
#' @param family A string specifying the exponential family to be used.
#'        Avaialble option are \code{"bernoulli"}, \code{"poisson"}, \code{"exponential"} or \code{"gamma"}.
#' @references Wilhelm, M., and Sangalli, L.M. (2016), Generalized Spatial Regression with Differential Regularization,
#'             Journal of Statistical Computation and Simulation, 86 (13), 2497-2518.
#' @rdname gsr
#' @order 1
#' @export
gsr <- function(formula, data, family, penalty = NULL) {
  return(.gsr$new(
    formula = deparse(formula),
    data = data,
    family = family,
    penalty = penalty
  ))
}

# An R6 class encapsulating the quantile spatial regression method with partial differential equation regularization.
#
#' @rdname qsr
#' @order 2
.qsr <- R6::R6Class(
  "qsr",
  private = list(
    model_ = NULL, ## cpp backend
    f_ = "fe_function"
  ),
  public = list(
    #' @description
    #' Creates a new \code{sr} object.
    #'
    #' @param formula R formula.
    #' @param data A \code{geoframe} containing the triangulation of the domain and the data (see also [geoframe()]).
    #' @param level A numeric value in (0,1) denoting the quantile level to estimate.
    #' @param penalty A penalty term returned by [fe_elliptic()].
    initialize = function(formula, data, level, penalty) {
      ## recover name of non-parametric field
      vars <- all.vars(as.formula(formula)[[3]])
      f_symbol <- setdiff(vars, names(data))
      fdapde_assert(length(f_symbol) == 1, "Expected exactly one nonparametric term.")
      private$f_ <- get(f_symbol, envir = globalenv())
      fe_type <- get_private(private$f_)$type

      if (is.null(penalty)) {
        ## fallback to laplacian penalty
        private$model_ = new(cpp_qsr_2_2, formula, get_private(data$gf__ptr__)$ptr_, level, penalty)
      } else {
        domain <- get_private(get_private(data$gf__ptr__)$mesh_)$mesh_
        params <- list()
        quad_nodes <- matrix()
        if (is.function(penalty$K) || is.function(penalty$b) || is.function(penalty$c) || is.function(penalty$u)) {
          quad_nodes <- fe_simplex_quad_nodes(domain, fe_type)
        }
        n_quad_nodes <- new(cpp_fe_space_2_2_p1, domain)$n_quad_nodes()
        embed_dim <- 2

        ## diffusion
        if (is.null(penalty$K)) {
          params$K = matrix(rep(0, times = (embed_dim * embed_dim) * n_quad_nodes), nrow = n_quad_nodes)
        } else {
          fdapde_assert(is.matrix(K) || is.function(K), "Invalid diffusion field type.")
          if (is.matrix(K)) {
            fdapde_assert(nrow(K) == embed_dim && ncol(K) == embed_dim, "Not a square matrix.")
            params$K = matrix(rep(c(K), times = n_quad_nodes), nrow = n_quad_nodes, byrow = TRUE)
          } else {
            q <- quad_nodes[1, ]
            fdapde_assert(
              nrow(K(q)) == embed_dim && ncol(K(q)) == embed_dim,
              "Not evaluates to a square matrix."
            )
            params$K = K(quad_nodes)
          }
        }
        ## transport
        if (is.null(penalty$b)) {
          params$b = matrix(rep(0, times = embed_dim * n_quad_nodes), nrow = n_quad_nodes)
        } else {
          fdapde_assert(is.matrix(b) || is.vector(b) || is.function(b), "Invalid transport field type.")
          if (is.matrix(b) || is.vector(b)) {
            fdapde_assert(nrow(b) == embed_dim && ncol(b) == 1, "Not a vector.")
            params$b = matrix(rep(c(b), times = n_quad_nodes), nrow = n_quad_nodes, byrow = TRUE)
          } else {
            q <- quad_nodes[1, ]
            fdapde_assert(nrow(b(q)) == embed_dim && ncol(b(q)) == 1, "Not evaluates to a vector.")
            params$b = b(quad_nodes)
          }
        }
        ## reaction
        if (is.null(penalty$c)) {
          params$c = matrix(rep(0, times = n_quad_nodes), nrow = n_quad_nodes)
        } else {
          fdapde_assert(is.numeric(c) || is.function(c), "Invalid reaction field type.")
          if (is.numeric(c)) {
            params$c = matrix(rep(c, times = n_quad_nodes), nrow = n_quad_nodes)
          } else {
            q <- quad_nodes[1, ]
            fdapde_assert(is.numeric(c(q)), "Not evaluates to a scalar.")
            params$c = c(quad_nodes)
          }
        }
        ## force
        if (is.null(penalty$u)) {
          params$u = matrix(rep(0, times = n_quad_nodes), nrow = n_quad_nodes)
        } else {
          fdapde_assert(is.numeric(u) || is.function(u), "Invalid forcing field type.")
          if (is.numeric(u)) {
            params$u = matrix(rep(u, times = n_quad_nodes), nrow = n_quad_nodes)
          } else {
            q <- quad_nodes[1, ]
            fdapde_assert(is.numeric(u(q)), "Not evaluates to a scalar.")
            params$c = c(quad_nodes)
          }
        }
        private$model_ = new(cpp_qsr_2_2, formula, get_private(data$gf__ptr__)$ptr_, level, params)
      }
    },
    #' @description
    #' Fits the statistical model.
    #'
    #' @param lambda The smoothing parameter. Default is \code{NULL}.
    #' @param calibrator A calibrator object. Available calibrators include [gcv()].
    fit = function(lambda = NULL, calibrator = NULL) {
      fdapde_assert(!is.null(calibrator) || !is.null(lambda), "Unable to select smoothing level.")
      if (is.null(calibrator)) {
        private$model_$fit(lambda)
      } else {
        r <- private$model_$fit_gcv(as.character(calibrator$opt_t), calibrator)
      }
      ## update estimated field coefficient vector
      get_private(private$f_)$fe_function_$set_coeff(private$model_$f())
      if (!is.null(calibrator)) return(r)
    }
  ),
  active = list(
    #' @field f The estimated nonparametric component of the model.
    f = function() private$model_$f(),
    #' @field beta The estimated parametric component of the model.
    beta = function() private$model_$beta(),
    #' @field fitted The vector of values fitted by the model.
    fitted = function() private$model_$fitted()
  )
)

#' Create a \code{qsr} object
#'
#' @inheritParams sr
#' @param level A numeric value in (0,1) denoting the quantile level to estimate.
#' @references Castiglione, C., Arnone, E., Bernardi, M., Farcomeni, A., and Sangalli, L.M. (2025),
#'             PDE-regularised spatial quantile regression, Journal of Multivariate Analysis, DOI: 10.1016/j.jmva.2024.105381.
#' @rdname qsr
#' @order 1
#' @export
qsr <- function(formula, data, level, penalty = NULL) {
  return(.qsr$new(
    formula = deparse(formula),
    data = data,
    level = level,
    penalty = penalty
  ))
}


#' Space-only Second-Order Penalty Term
#'
#' A list encapsulating the coefficients of a second-order linear differential operator
#' typically used to define spatial smoothing penalties in regularization frameworks.
#'
#' @param K The diffusion tensor. A \code{local_dim} \eqn{\times} \code{local_dim} matrix encoding
#' spatially varying diffusion. Can also be \code{NULL} to omit this term. Default is \code{NULL}.
#' @param b The advection vector field. A \code{local_dim} \eqn{\times} 1 matrix representing
#' first-order derivatives. Can also be \code{NULL}. Default is \code{NULL}.
#' @param c The reaction term. A scalar \code{numeric} coefficient multiplying the zeroth-order
#' component. Can also be \code{NULL}. Default is \code{NULL}.
#' @param u The forcing term or known right-hand side function. Default is \code{NULL}. (... da completare ...)
#' @export
fe_elliptic <- function(K = NULL, b = NULL, c = NULL, u = NULL) {
  fdapde_assert(!(is.null(K) && is.null(b) && is.null(c) && is.null(u)))
  return(list(K = K, b = b, c = c, u = u))
}

#' Generalized Cross Validation (GCV) Calibrator
#'
#' A list encapsulating the Generalized Cross Validation (GCV) objective, which balances
#' goodness-of-fit and model complexity. This calibrator is typically used to select
#' smoothing parameters by minimizing the GCV score.
#'
#' @param optimizer An optimizer object. The library currently provides the brute-force
#'        grid search optimizer \code{grid_optimizer}, which can be created using [grid_optimizer()].
#' @param edf A string specifying the method to compute the equivalent degrees of freedom when minimizing the GCV.
#'        Available option: \code{"stochastic"}. Default is \code{"stochastic"}.
#' @param mc_samples An integer specifying the number of ... to use when
#'        estimating the stochastic equivalent degrees of freedom. Default is \code{100}.
#' @param seed An integer seed for the random number generator used in computing the
#'        \code{"stochastic"} EDF. Default is \code{NULL}, meaning ... (seed is set to 0???).
#' @export
gcv <- function(optimizer, edf = "stochastic", mc_samples = 100, seed = NULL) {
  args <- append(list(), optimizer)
  args$mc_samples = mc_samples
  args$seed = if (is.null(seed)) -1 else seed
  return(args)
}

#' A Grid Search Optimizer
#'
#' An optimizer that performs a brute-force search over a predefined set of candidate values
#' to minimize an objective function. This approach exhaustively evaluates the function at
#' each point in the grid and selects the value that yields the lowest objective.
#'
#' @param grid A numeric vector of values to be considered during the optimization.
#' @export
grid_search <- function(grid) {
  return(list(opt_t = "grid", grid = grid))
}

#' An Inexact Newton Optimizer
#' 
#' An optimizer that estimates the derivative of the objective function using finite differences. 
#' The minimization is performed through an iterative scheme.
#' 
#' @param max_iter An integer specifying the maximum number of iterations. Default is \code{100}.
#' @param tolerance A numeric value between 0 and 1 controlling the precision of the optimization: smaller values yield higher accuracy.
#' @param step A numeric value specifying the step size used during each update.
#' 
#' @export
newton_fd <- function(max_iter = 100, tolerance = 0.01, step = 0.01) {
  return(list(
    opt_t = "newton_fe",
    max_iter = max_iter,
    tolerance = tolerance,
    step = step
  ))
}

#' A Gradient Descent Optimizer
#' 
#' An optimizer that minimizes an objective function using its gradient. 
#' The minimization is performed through an iterative scheme.
#' 
#' @param max_iter An integer specifying the maximum number of iterations. Default is \code{100}.
#' @param tolerance A numeric value between 0 and 1 controlling the precision of the optimization: smaller values yield higher accuracy.
#' @param step A numeric value specifying the step size used during each update.
#' 
#' @export
gradient_descent <- function(max_iter = 100, tolerance = 0.01, step = 0.01) {
  return(list(
    opt_t = "gradient_descent",
    max_iter = max_iter,
    tolerance = tolerance,
    step = step
  ))
}

#' A BFGS Optimizer
#' 
#' An optimizer that minimizes an objective function using the Broyden–Fletcher–Goldfarb–Shanno (BFGS) quasi-Newton method. 
#' The optimization is performed through an iterative scheme that updates an approximation of the inverse Hessian matrix.
#' 
#' @param max_iter An integer specifying the maximum number of iterations. Default is \code{100}.
#' @param tolerance A numeric value between 0 and 1 controlling the precision of the optimization: smaller values yield higher accuracy.
#' @param step A numeric value specifying the initial step size used during the updates.
#' 
#' @export
bfgs <- function(max_iter = 100, tolerance = 0.01, step = 0.01) {
  return(list(
    opt_t = "bfgs",
    max_iter = max_iter,
    tolerance = tolerance,
    step = step
  ))
}
