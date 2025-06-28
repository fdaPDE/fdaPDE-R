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
#include "../inst/include/geoframe.h"

namespace fdapde {
namespace r {

// clang-format off
  
#define geoframe_rcpp_interface(Triangulation)                                                                         \
     constructor<Rcpp::Environment>()                                                                                  \
      .constructor<Rcpp::Environment, std::string, std::vector<int>, std::vector<std::string>>()                       \
      .method("flt64_assign"    , &GeoFrame<Triangulation>::assign<double>)                                            \
      .method("flt32_assign"    , &GeoFrame<Triangulation>::assign<float>)                                             \
      .method("int64_assign"    , &GeoFrame<Triangulation>::assign<std::int64_t>)                                      \
      .method("int32_assign"    , &GeoFrame<Triangulation>::assign<std::int32_t>)                                      \
      .method("bin_assign"      , &GeoFrame<Triangulation>::assign<bool>)                                              \
      .method("str_assign"      , &GeoFrame<Triangulation>::assign<std::string>)                                       \
      .method("flt64_access"    , &GeoFrame<Triangulation>::access<double>)                                            \
      .method("flt32_access"    , &GeoFrame<Triangulation>::access<float>)                                             \
      .method("int64_access"    , &GeoFrame<Triangulation>::access<std::int64_t>)                                      \
      .method("int32_access"    , &GeoFrame<Triangulation>::access<std::int32_t>)                                      \
      .method("bin_access"      , &GeoFrame<Triangulation>::access<bool>)                                              \
      .method("str_access"      , &GeoFrame<Triangulation>::access<std::string>)                                       \
      .method("ltype"           , &GeoFrame<Triangulation>::ltype)                                                     \
      .method("ctype"           , &GeoFrame<Triangulation>::ctype)                                                     \
      .method("rows"            , &GeoFrame<Triangulation>::rows)                                                      \
      .method("cols"            , &GeoFrame<Triangulation>::cols)                                                      \
      .method("colnames"        , &GeoFrame<Triangulation>::colnames)                                                  \
      .method("bbox"            , &GeoFrame<Triangulation>::bbox)                                                      \
      .method("n_nodes"         , &GeoFrame<Triangulation>::n_nodes)                                                   \
      .method("n_cells"         , &GeoFrame<Triangulation>::n_cells)                                                   \
      /* point layer */                                                                                                \
      .method("insert_scalar_point_layer" , &GeoFrame<Triangulation>::insert_scalar_point_layer)                       \
      .method("point_coordinates"         , &GeoFrame<Triangulation>::point_coordinates)
     
using cpp_geoframe_2_2 = GeoFrame<fdapde::Triangulation<2, 2>>;
RCPP_MODULE(cpp_geoframe_2_2) {
    Rcpp::class_<GeoFrame<fdapde::Triangulation<2, 2>>>("cpp_geoframe_2_2")
      .geoframe_rcpp_interface(fdapde::Triangulation<2 FDAPDE_COMMA 2>);
}

// clang-format on

}   // namespace r
}   // namespace fdapde
