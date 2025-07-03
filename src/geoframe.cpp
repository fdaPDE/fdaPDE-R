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

using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
  
#define geoframe_rcpp_interface(Triangulation)                                                                         \
     constructor<Rcpp::Environment>()                                                                                  \
      .constructor<Rcpp::Environment, std::string, std::vector<int>, std::vector<std::string>>()                       \
      .method("flt64_assign"    , &GeoFrame<Triangulation>::assign<double>)                                            \
      .method("flt32_assign"    , &GeoFrame<Triangulation>::assign<float>)                                             \
      .method("int64_assign"    , &GeoFrame<Triangulation>::assign<std::int64_t>)                                      \
      .method("int32_assign"    , &GeoFrame<Triangulation>::assign<std::int32_t>)                                      \
      .method("str_assign"      , &GeoFrame<Triangulation>::assign<std::string>)                                       \
      .method("flt64_access"    , &GeoFrame<Triangulation>::access<double>)                                            \
      .method("flt32_access"    , &GeoFrame<Triangulation>::access<float>)                                             \
      .method("int64_access"    , &GeoFrame<Triangulation>::access<std::int64_t>)                                      \
      .method("int32_access"    , &GeoFrame<Triangulation>::access<std::int32_t>)                                      \
      .method("str_access"      , &GeoFrame<Triangulation>::access<std::string>)                                       \
      .method("flt64_insert"    , &GeoFrame<Triangulation>::insert<double>)                                            \
      .method("flt32_insert"    , &GeoFrame<Triangulation>::insert<float>)                                             \
      .method("int64_insert"    , &GeoFrame<Triangulation>::insert<std::int64_t>)                                      \
      .method("int32_insert"    , &GeoFrame<Triangulation>::insert<std::int32_t>)                                      \
      .method("str_insert"      , &GeoFrame<Triangulation>::insert<std::string>)                                       \
      .method("ltype"           , &GeoFrame<Triangulation>::ltype)                                                     \
      .method("flt64_blk_insert", &GeoFrame<Triangulation>::blk_insert<double>)                                        \
      .method("int64_blk_insert", &GeoFrame<Triangulation>::blk_insert<std::int64_t>)                                  \
      .method("ltype"           , &GeoFrame<Triangulation>::ltype)                                                     \
      .method("ctype"           , &GeoFrame<Triangulation>::ctype)                                                     \
      .method("rows"            , &GeoFrame<Triangulation>::rows)                                                      \
      .method("cols"            , &GeoFrame<Triangulation>::cols)                                                      \
      .method("colnames"        , &GeoFrame<Triangulation>::colnames)                                                  \
      .method("colnames_all"    , &GeoFrame<Triangulation>::colnames_all)                                              \
      .method("bbox"            , &GeoFrame<Triangulation>::bbox)                                                      \
      .method("n_nodes"         , &GeoFrame<Triangulation>::n_nodes)                                                   \
      .method("n_cells"         , &GeoFrame<Triangulation>::n_cells)                                                   \
      /* point layer */                                                                                                \
      .method("insert_scalar_point_layer" , &GeoFrame<Triangulation>::insert_scalar_point_layer<matrix_t>)             \
      .method("insert_scalar_point_layer_mesh_nodes" , &GeoFrame<Triangulation>::insert_scalar_point_layer<int>)       \
      .method("point_coordinates"         , &GeoFrame<Triangulation>::point_coordinates)                               \
      /* areal layer */                                                                                                \
      .method("insert_scalar_areal_layer" , &GeoFrame<Triangulation>::insert_scalar_areal_layer)                       \
      .method("load_shp"                  , &GeoFrame<Triangulation>::load_shp)

// colnames_all can be done directly from R
	      
using cpp_geoframe_2_2 = GeoFrame<fdapde::Triangulation<2, 2>>;
RCPP_MODULE(cpp_geoframe_2_2) {
    Rcpp::class_<GeoFrame<fdapde::Triangulation<2, 2>>>("cpp_geoframe_2_2")
      .geoframe_rcpp_interface(fdapde::Triangulation<2 FDAPDE_COMMA 2>);
}

// clang-format on

}   // namespace r
}   // namespace fdapde
