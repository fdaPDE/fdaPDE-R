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
#include "../inst/include/de.h"

namespace fdapde {
namespace r {

// clang-format on
  
#define fe_de_elliptic_rcpp_interface(LocalDim, EmbedDim)                                                              \
     method("fit"        , &de_elliptic<LocalDim, EmbedDim>::fit        )                                              \
    .method("density"    , &de_elliptic<LocalDim, EmbedDim>::density    )                                              \
    .method("log_density", &de_elliptic<LocalDim, EmbedDim>::log_density)                                              \
    .method("fitted"     , &de_elliptic<LocalDim, EmbedDim>::fitted     )

// spatial regression
using cpp_de_2_2 = de_elliptic<2, 2>;
RCPP_MODULE(cpp_de_2_2) {
    Rcpp::class_<de_elliptic<2, 2>>("cpp_de_2_2")
      .constructor<Rcpp::Environment, Rcpp::Nullable<Rcpp::List>>()
      .fe_de_elliptic_rcpp_interface(2, 2);
}
  
// clang-format on
  
}   // namespace r
}   // namespace fdapde
