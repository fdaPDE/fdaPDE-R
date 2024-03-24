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

#include "../inst/include/pde.h"

namespace fdapde {
namespace r {
  
#define pde_rcpp_interface(LocalDim, EmbedDim, Order)                                                                  \
       constructor<Rcpp::Environment, int, Rcpp::Nullable<Rcpp::List>>()                                               \
      .method("get_quadrature_nodes" , &PDE_<LocalDim, EmbedDim, Order>::get_quadrature_nodes )                        \
      .method("get_dofs_coordinates" , &PDE_<LocalDim, EmbedDim, Order>::get_dofs_coordinates )                        \
      .method("mass"                 , &PDE_<LocalDim, EmbedDim, Order>::mass                 )                        \
      .method("stiff"                , &PDE_<LocalDim, EmbedDim, Order>::stiff                )                        \
      .method("force"                , &PDE_<LocalDim, EmbedDim, Order>::force                )                        \
      .method("set_dirichlet_bc"     , &PDE_<LocalDim, EmbedDim, Order>::set_dirichlet_bc     )                        \
      .method("set_forcing"          , &PDE_<LocalDim, EmbedDim, Order>::set_forcing          )                        \
      .method("set_initial_condition", &PDE_<LocalDim, EmbedDim, Order>::set_initial_condition)                        \
      .method("init"                 , &PDE_<LocalDim, EmbedDim, Order>::init                 )

using cpp_pde_2_2_1 = PDE_<2, 2, 1>;
RCPP_MODULE(cpp_pde_2_2_1) { Rcpp::class_<PDE_<2, 2, 1>>("cpp_pde_2_2_1").pde_rcpp_interface(2, 2, 1); }
// using cpp_pde_2_2_2 = PDE_<2, 2, 2>;
// RCPP_MODULE(cpp_pde_2_2_2) { Rcpp::class_<PDE_<2, 2, 2>>("cpp_pde_2_2_2").pde_rcpp_interface(2, 2, 2); }
// using cpp_pde_2_3_1 = PDE_<2, 3, 1>;
// RCPP_MODULE(cpp_pde_2_3_1) { Rcpp::class_<PDE_<2, 3, 1>>("cpp_pde_2_3_1").pde_rcpp_interface(2, 3, 1); }
// using cpp_pde_2_3_2 = PDE_<2, 3, 2>;
// RCPP_MODULE(cpp_pde_2_3_2) { Rcpp::class_<PDE_<2, 3, 2>>("cpp_pde_2_3_2").pde_rcpp_interface(2, 3, 2); }
using cpp_pde_3_3_1 = PDE_<3, 3, 1>;
RCPP_MODULE(cpp_pde_3_3_1) { Rcpp::class_<PDE_<3, 3, 1>>("cpp_pde_3_3_1").pde_rcpp_interface(3, 3, 1); }
// using cpp_pde_3_3_2 = PDE_<3, 3, 2>;
// RCPP_MODULE(cpp_pde_3_3_2) { Rcpp::class_<PDE_<3, 3, 2>>("cpp_pde_3_3_2").pde_rcpp_interface(3, 3, 2); }

}   // namespace r
}   // namespace fdapde
