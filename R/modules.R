#' @useDynLib fdaPDE2
#' @import methods Rcpp
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL

## load required modules
Rcpp::loadModule("cpp_network_domain", TRUE)
Rcpp::loadModule("cpp_2d_domain",      TRUE)
Rcpp::loadModule("cpp_surface_domain", TRUE)
Rcpp::loadModule("cpp_3d_domain",      TRUE)
Rcpp::loadModule("cpp_pde_2d_fe1",     TRUE)
Rcpp::loadModule("cpp_pde_2d_fe2",     TRUE)
Rcpp::loadModule("cpp_srpde",          TRUE)
Rcpp::loadModule("cpp_gsrpde_s",       TRUE)
Rcpp::loadModule("cpp_gcv",            TRUE)
Rcpp::loadModule("cpp_lagrange_basis_2d_fe1", TRUE)
Rcpp::loadModule("cpp_lagrange_basis_2d_fe2", TRUE)
#Rcpp::loadModule("R_GSRPDE", TRUE)
#Rcpp::loadModule("R_STRPDE", TRUE)
