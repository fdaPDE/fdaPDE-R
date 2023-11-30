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
#include "headers/r_gsrpde.h"
using fdapde::models::SpaceOnly;

#define CPP_GSRPDE_COMMON_METHODS(regularization_type)                                                                 \
    .method("f"               , &R_GSRPDE<regularization_type>::f)                                                     \
    .method("fitted"          , &R_GSRPDE<regularization_type>::fitted)                                                \
    .method("get_gcv"         , &R_GSRPDE<regularization_type>::get_gcv)                                               \
    .method("set_lambda_D"    , &R_GSRPDE<regularization_type>::set_lambda_D)                                          \
    .method("set_observations", &R_GSRPDE<regularization_type>::set_observations)                                      \
    .method("set_covariates"  , &R_GSRPDE<regularization_type>::set_covariates)                                        \
    .method("init"            , &R_GSRPDE<regularization_type>::init)                                                  \
    .method("solve"           , &R_GSRPDE<regularization_type>::solve)

// Rcpp modules definition
using cpp_gsrpde_s = R_GSRPDE<SpaceOnly>;
RCPP_MODULE(cpp_gsrpde_s) {
    Rcpp::class_<R_GSRPDE<SpaceOnly>>("cpp_gsrpde_s")
      .constructor<Rcpp::Environment, int, std::string>()
      CPP_GSRPDE_COMMON_METHODS(SpaceOnly);
}
