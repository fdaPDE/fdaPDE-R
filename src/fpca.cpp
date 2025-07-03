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
#include "../inst/include/fpca.h"

namespace fdapde {
namespace r {

// clang-format on

#define fe_fpca_laplace_rcpp_interface(LocalDim, EmbedDim)                                                             \
     method("fit"      , &fpca_laplace<LocalDim, EmbedDim>::fit          )                                             \
    .method("loadings" , &fpca_laplace<LocalDim, EmbedDim>::Fn           )                                             \
    .method("scores"   , &fpca_laplace<LocalDim, EmbedDim>::S            )                                             \
    .method("pcs"      , &fpca_laplace<LocalDim, EmbedDim>::F            )                                             \
    .method("pcs_norm" , &fpca_laplace<LocalDim, EmbedDim>::loadings_norm)                                             \
    .method("lambda"   , &fpca_laplace<LocalDim, EmbedDim>::lambda       )
  
using cpp_fpca_laplace_2_2 = fpca_laplace<2, 2>;
RCPP_MODULE(cpp_fpca_laplace_2_2) {
    Rcpp::class_<fpca_laplace<2, 2>>("cpp_fpca_laplace_2_2")
      .constructor<std::string, Rcpp::Environment>()
      .fe_fpca_laplace_rcpp_interface(2, 2);
}

// clang-format on
  
}   // namespace r
}   // namespace fdapde
