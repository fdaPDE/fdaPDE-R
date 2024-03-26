#' @useDynLib fdaPDE2
#' @import methods Rcpp
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL

## load required modules
Rcpp::loadModule("cpp_mesh_1_1", TRUE)
Rcpp::loadModule("cpp_mesh_1_2", TRUE)
Rcpp::loadModule("cpp_mesh_2_2", TRUE)
Rcpp::loadModule("cpp_mesh_2_3", TRUE)
Rcpp::loadModule("cpp_mesh_3_3", TRUE)
Rcpp::loadModule("cpp_pde_2_2_1", TRUE)
Rcpp::loadModule("cpp_pde_3_3_1", TRUE)
Rcpp::loadModule("cpp_fe_space_lagrange_2_2_1", TRUE)
Rcpp::loadModule("cpp_fe_space_lagrange_3_3_1", TRUE)
## Rcpp::loadModule("cpp_fe_space_lagrange_2_2_2",       TRUE)
Rcpp::loadModule("cpp_bspline_space", TRUE)

Rcpp::loadModule("cpp_srpde", TRUE)
Rcpp::loadModule("cpp_gsrpde_space", TRUE)
