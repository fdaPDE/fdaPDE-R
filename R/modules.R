#' @useDynLib fdaPDE2
#' @import methods Rcpp
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL

## load required modules

## mesh and pde
Rcpp::loadModule("cpp_network_domain", TRUE)
Rcpp::loadModule("cpp_2d_domain", TRUE)
Rcpp::loadModule("cpp_surface_domain", TRUE)
Rcpp::loadModule("cpp_3d_domain", TRUE)
Rcpp::loadModule("cpp_pde_2d_fe1", TRUE)
Rcpp::loadModule("cpp_pde_2d_fe2", TRUE)
Rcpp::loadModule("cpp_lagrange_basis_2d_fe1", TRUE)
Rcpp::loadModule("cpp_lagrange_basis_2d_fe2", TRUE)

## calibrators
Rcpp::loadModule("cpp_off", TRUE)
# Rcpp::loadModule("cpp_gcv", TRUE)
Rcpp::loadModule("cpp_kcv", TRUE)

## regression models
Rcpp::loadModule("cpp_srpde", TRUE)

## functional models
Rcpp::loadModule("cpp_center", TRUE)
