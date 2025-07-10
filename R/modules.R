#' @useDynLib fdaPDE2
#' @import methods Rcpp
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL

## load required modules
Rcpp::loadModule("cpp_triangulation_2_2", TRUE)


Rcpp::loadModule("cpp_fe_function_2_2_p1", TRUE)

Rcpp::loadModule("cpp_geoframe_2_2", TRUE)

Rcpp::loadModule("cpp_sr_2_2", TRUE)
Rcpp::loadModule("cpp_gsr_2_2", TRUE)
Rcpp::loadModule("cpp_qsr_2_2", TRUE)
Rcpp::loadModule("cpp_fpca_laplace_2_2", TRUE)
Rcpp::loadModule("cpp_de_2_2", TRUE)

Rcpp::loadModule("cpp_fe_space_2_2_p1", TRUE)
