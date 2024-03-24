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
#include "../inst/include/functional_space.h"

namespace fdapde {
namespace r {

#define fe_functional_space_rcpp_interface(LocalDim, EmbedDim, Order)                                                  \
       method("size",             &FEFunctionalSpace<core::Mesh<LocalDim, EmbedDim>, Order>::size            )         \
      .method("eval",             &FEFunctionalSpace<core::Mesh<LocalDim, EmbedDim>, Order>::Psi             )         \
      .method("integrate",        &FEFunctionalSpace<core::Mesh<LocalDim, EmbedDim>, Order>::integrate       )         \
      .method("dofs_coords",      &FEFunctionalSpace<core::Mesh<LocalDim, EmbedDim>, Order>::dofs_coords     )         \
      .method("eval_expansion",   &FEFunctionalSpace<core::Mesh<LocalDim, EmbedDim>, Order>::eval_expansion  )         \
      .method("quadrature_nodes", &FEFunctionalSpace<core::Mesh<LocalDim, EmbedDim>, Order>::quadrature_nodes)
  
using cpp_fe_space_lagrange_2_2_1 = FEFunctionalSpace<core::Mesh<2, 2>, 1>;
RCPP_MODULE(cpp_fe_space_lagrange_2_2_1) {
    Rcpp::class_<FEFunctionalSpace<core::Mesh<2, 2>, 1>>("cpp_fe_space_lagrange_2_2_1")
      .constructor<Rcpp::Environment, int>()
      .fe_functional_space_rcpp_interface(2, 2, 1);
}
// using cpp_fe_space_lagrange_2_2_2 = FEFunctionalSpace<core::Mesh<2, 2>, 2>;
// RCPP_MODULE(cpp_fe_space_lagrange_2_2_2) {
//     Rcpp::class_<FEFunctionalSpace<core::Mesh<2, 2>, 2>>("cpp_fe_space_lagrange_2_2_2")
//       .constructor<Rcpp::Environment, int>()
//       .fe_functional_space_rcpp_interface(2, 2, 2);
// }
using cpp_fe_space_lagrange_2_3_1 = FEFunctionalSpace<core::Mesh<2, 3>, 1>;
RCPP_MODULE(cpp_fe_space_lagrange_2_3_1) {
    Rcpp::class_<FEFunctionalSpace<core::Mesh<2, 3>, 1>>("cpp_fe_space_lagrange_2_3_1")
      .constructor<Rcpp::Environment, int>()
      .fe_functional_space_rcpp_interface(2, 3, 1);
}
// using cpp_fe_space_lagrange_2_3_2 = FEFunctionalSpace<core::Mesh<2, 3>, 2>;
// RCPP_MODULE(cpp_fe_space_lagrange_2_3_2) {
//     Rcpp::class_<FEFunctionalSpace<core::Mesh<2, 3>, 2>>("cpp_fe_space_lagrange_2_3_2")
//       .constructor<Rcpp::Environment, int>()
//       .fe_functional_space_rcpp_interface(2, 3, 2);
// }
using cpp_fe_space_lagrange_3_3_1 = FEFunctionalSpace<core::Mesh<3, 3>, 1>;
RCPP_MODULE(cpp_fe_space_lagrange_3_3_1) {
    Rcpp::class_<FEFunctionalSpace<core::Mesh<3, 3>, 1>>("cpp_fe_space_lagrange_3_3_1")
      .constructor<Rcpp::Environment, int>()
      .fe_functional_space_rcpp_interface(3, 3, 1);
}
// using cpp_fe_space_lagrange_3_3_2 = FEFunctionalSpace<core::Mesh<3, 3>, 2>;
// RCPP_MODULE(cpp_fe_space_lagrange_3_3_2) {
//     Rcpp::class_<FEFunctionalSpace<core::Mesh<3, 3>, 2>>("cpp_fe_space_lagrange_3_3_2")
//       .constructor<Rcpp::Environment, int>()
//       .fe_functional_space_rcpp_interface(3, 3, 2);
// }

using cpp_bspline_space = BSplineFunctionalSpace<3>;   // cubic BSpline basis function
RCPP_MODULE(cpp_bspline_space) {
    Rcpp::class_<BSplineFunctionalSpace<3>>("cpp_bspline_space")
      .constructor<Rcpp::Environment, int>()
      .method("eval_expansion",   &BSplineFunctionalSpace<3>::eval_expansion  )
      .method("quadrature_nodes", &BSplineFunctionalSpace<3>::quadrature_nodes)
      .method("size",             &BSplineFunctionalSpace<3>::size            )
      .method("eval",             &BSplineFunctionalSpace<3>::Psi             );
}
    
}   // namespace r
}   // namespace fdapde
