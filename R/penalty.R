## TODOs:
## TODO: add the other funnctional spaces
## TODO: the forcing term u needs to be generalized to the case of fe_order > 1
## TODO: add the support for generic diffusion advection reaction regularization term
## TODO: add the support for time dependent regularization

# functional spaces ----

functions_space <- function(domain, order, basis) {
  ## function space
  Vh <- switch(basis,
    "FEM" = {
      ## the functional space of finite element
      ## functions of order fe_order, over the domain
      FunctionSpace(domain, fe_order = as.integer(order))
    },
    stop("This type of functional space is not implemented")
  )
  ## generic function of Vh
  f <- Function(Vh)
  ## homogeneous forcing term
  u <- function(points) {
    return(matrix(0, nrow = nrow(points), ncol = 1)) ## ?? correct number of rows if order > 1
  }
  return(list(
    Vh = Vh,
    f = f,
    u = u
  ))
}

# penalty ----

## callable simple lapalcian pde
simple_laplacian_pde_callable <- function(domain, parameters) {
  ## unpack parameters
  order <- parameters$order
  basis <- parameters$basis
  ## generate the functions space
  functions_space_result <- functions_space(domain, order, basis)
  f <- functions_space_result$f
  u <- functions_space_result$u
  ## differential operator
  L <- -laplace(f)
  ## pde
  pde <- PDE(L, u)
  return(pde)
}

#' @export
## simple laplacian penalty user interface
## it wraps together:
## - simple_laplacian_pde_callable
## - parameters
## this object is meant to contain all the information necessary to define
## a proper pde object once the actual domain is provided
simple_laplacian_penalty <- function(order = 1, basis = c("FEM")) {
  penalty <- list(
    pde_callable = simple_laplacian_pde_callable,
    parameters = list(
      order = order,
      basis = match.arg(basis)
    )
  )
  return(penalty)
}
