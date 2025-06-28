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
#include "../inst/include/sr.h"

namespace fdapde {
namespace r {

// clang-format off
  
using cpp_sr_2_2 = SRPDE<2, 2>;
RCPP_MODULE(cpp_sr_2_2) {
  Rcpp::class_<SRPDE<2, 2>>("cpp_sr_2_2")
    .constructor<std::string, Rcpp::Environment, Rcpp::Nullable<Rcpp::List>>()
    .method("fit"   , &SRPDE<2, 2>::fit   )
    .method("f"     , &SRPDE<2, 2>::f     )
    .method("beta"  , &SRPDE<2, 2>::beta  )
    .method("fitted", &SRPDE<2, 2>::fitted);
}

  
// clang-format on

}   // namespace r
}   // namespace fdapde
