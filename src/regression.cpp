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

#include "../inst/include/regression.h"

namespace fdapde {
namespace r {
  
#define regression_rcpp_interface(Model)                                                                               \
       method("set_lambda"      , &RegressionModel<Model>::set_lambda      )                                           \
      .method("set_covariates"  , &RegressionModel<Model>::set_covariates  )                                           \
      .method("set_observations", &RegressionModel<Model>::set_observations)                                           \
      .method("fit"             , &RegressionModel<Model>::fit             )                                           \
      .method("fitted"          , &RegressionModel<Model>::fitted          )                                           \
      .method("f"               , &RegressionModel<Model>::f               )                                           \
      .method("gcvs"            , &RegressionModel<Model>::gcvs            )                                           \
      .method("edfs"            , &RegressionModel<Model>::edfs            )                                           \
      .method("optimum"         , &RegressionModel<Model>::optimum         ) 

using cpp_srpde = SRPDE;
RCPP_MODULE(cpp_srpde) {
    Rcpp::class_<RegressionModel<models::SRPDE>>("cpp_regression").regression_rcpp_interface(models::SRPDE);
    Rcpp::class_<SRPDE>("cpp_srpde")
      .derives<RegressionModel<models::SRPDE>>("cpp_regression")
      .constructor<Rcpp::Environment, int>();
};

using cpp_gsrpde_space = GSRPDE<models::SpaceOnly>;
RCPP_MODULE(cpp_gsrpde_space) {
    Rcpp::class_<RegressionModel<models::GSRPDE<models::SpaceOnly>>>("cpp_regression")
      .regression_rcpp_interface(models::GSRPDE<models::SpaceOnly>);
    Rcpp::class_<GSRPDE<models::SpaceOnly>>("cpp_gsrpde_space")
      .derives<RegressionModel<models::GSRPDE<models::SpaceOnly>>>("cpp_regression")
      .constructor<Rcpp::Environment, int, int>()
      .method("set_fpirls_tolerance", &GSRPDE<models::SpaceOnly>::set_fpirls_tolerance)
      .method("set_fpirls_max_iter" , &GSRPDE<models::SpaceOnly>::set_fpirls_max_iter );
};

  
}   // namespace r
}   // namespace fdapde
