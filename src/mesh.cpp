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
#include "../inst/include/mesh.h"

namespace fdapde {
namespace r {
  
#define mesh_rcpp_interface(LocalDim, EmbedDim)                                                                        \
       constructor<Rcpp::List>()                                                                                       \
      .method("nodes"    , &Mesh<LocalDim, EmbedDim>::nodes    )                                                       \
      .method("elements" , &Mesh<LocalDim, EmbedDim>::elements )                                                       \
      .method("neighbors", &Mesh<LocalDim, EmbedDim>::neighbors)                                                       \
      .method("boundary" , &Mesh<LocalDim, EmbedDim>::boundary )                                                       \
      .method("locate"   , &Mesh<LocalDim, EmbedDim>::locate   )

using cpp_mesh_1_1 = Mesh<1, 1>;
RCPP_MODULE(cpp_mesh_1_1) { Rcpp::class_<Mesh<1, 1>>("cpp_mesh_1_1").mesh_rcpp_interface(1, 1); }
using cpp_mesh_1_2 = Mesh<1, 2>;
RCPP_MODULE(cpp_mesh_1_2) { Rcpp::class_<Mesh<1, 2>>("cpp_mesh_1_2").mesh_rcpp_interface(1, 2); }
using cpp_mesh_2_2 = Mesh<2, 2>;
RCPP_MODULE(cpp_mesh_2_2) { Rcpp::class_<Mesh<2, 2>>("cpp_mesh_2_2").mesh_rcpp_interface(2, 2); }
using cpp_mesh_2_3 = Mesh<2, 3>;
RCPP_MODULE(cpp_mesh_2_3) { Rcpp::class_<Mesh<2, 3>>("cpp_mesh_2_3").mesh_rcpp_interface(2, 3); }
using cpp_mesh_3_3 = Mesh<3, 3>;
RCPP_MODULE(cpp_mesh_3_3) { Rcpp::class_<Mesh<3, 3>>("cpp_mesh_3_3").mesh_rcpp_interface(3, 3); }

}   // namespace r
}   // namespace fdapde
