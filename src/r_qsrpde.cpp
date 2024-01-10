// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "headers/r_qsrpde.h"
using fdapde::models::SpaceOnly;

#define CPP_QSRPDE_COMMON_METHODS(regularization_type)                                                                 \
    .method("f"               , &R_QSRPDE<regularization_type>::f)                                                     \
    .method("fitted"          , &R_QSRPDE<regularization_type>::fitted)                                                \
    .method("get_gcv"         , &R_QSRPDE<regularization_type>::get_gcv)                                               \
    .method("set_lambda_D"    , &R_QSRPDE<regularization_type>::set_lambda_D)                                          \
    .method("set_observations", &R_QSRPDE<regularization_type>::set_observations)                                      \
    .method("set_covariates"  , &R_QSRPDE<regularization_type>::set_covariates)                                        \
    .method("init"            , &R_QSRPDE<regularization_type>::init)                                                  \
    .method("solve"           , &R_QSRPDE<regularization_type>::solve)

// Rcpp modules definition
using cpp_qsrpde_s = R_QSRPDE<SpaceOnly>;
RCPP_MODULE(cpp_qsrpde_s) {
    Rcpp::class_<R_QSRPDE<SpaceOnly>>("cpp_qsrpde_s")
      .constructor<Rcpp::Environment, int>()
      CPP_QSRPDE_COMMON_METHODS(SpaceOnly)
      .method("set_alpha"     , &R_QSRPDE<SpaceOnly>::set_alpha);
}
