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
#include "headers/r_regression_model.h"

// Rcpp modules definition
using cpp_srpde = R_SRPDE;
RCPP_MODULE(cpp_srpde) {
    Rcpp::class_<R_SRPDE>("cpp_srpde") 
      .constructor<Rcpp::Environment, int, Rcpp::List>() 
      .method("get_view"               , &R_SRPDE::get_view               )
      .method("get_gcv"                , &R_SRPDE::get_gcv                )
      .method("f"                      , &R_SRPDE::f                      )
      .method("beta"                   , &R_SRPDE::beta                   ) 
      .method("set_lambda"             , &R_SRPDE::set_lambda             ) 
      .method("set_observations"       , &R_SRPDE::set_observations       )
      .method("set_covariates"         , &R_SRPDE::set_covariates         )
      .method("set_spatial_locations"  , &R_SRPDE::set_spatial_locations  ) 
      .method("init"                   , &R_SRPDE::init                   )
      .method("solve"                  , &R_SRPDE::solve                  );  
}
