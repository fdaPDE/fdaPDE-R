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

// clang-format off
  
#define triangulation_rcpp_interface(LocalDim, EmbedDim)                                                               \
       method("nodes"                 , &TriangulationBase<LocalDim, EmbedDim>::nodes                 )                \
      .method("cells"                 , &TriangulationBase<LocalDim, EmbedDim>::cells                 )                \
      .method("boundary_nodes"        , &TriangulationBase<LocalDim, EmbedDim>::boundary_nodes        )                \
      .method("n_nodes"               , &TriangulationBase<LocalDim, EmbedDim>::n_nodes               )                \
      .method("n_cells"               , &TriangulationBase<LocalDim, EmbedDim>::n_cells               )                \
      .method("n_boundary_nodes"      , &TriangulationBase<LocalDim, EmbedDim>::n_boundary_nodes      )                \
      .method("bbox"                  , &TriangulationBase<LocalDim, EmbedDim>::bbox                  )                \
      .method("measure"               , &TriangulationBase<LocalDim, EmbedDim>::measure               )                \
      .method("marked_measure"        , &TriangulationBase<LocalDim, EmbedDim>::marked_measure        )                \
      .method("sample"                , &TriangulationBase<LocalDim, EmbedDim>::sample                )                \
      .method("mark_cells"            , &TriangulationBase<LocalDim, EmbedDim>::mark_cells            )                \
      .method("cells_markers"         , &TriangulationBase<LocalDim, EmbedDim>::cells_markers         )                \
      .method("clear_cells_markers"   , &TriangulationBase<LocalDim, EmbedDim>::clear_cells_markers   )                \
      .method("filter_cells_by_marker", &TriangulationBase<LocalDim, EmbedDim>::filter_cells_by_marker)                \
      .method("cell_coords"           , &TriangulationBase<LocalDim, EmbedDim>::cell_coords           )                \
      .method("cell_measure"          , &TriangulationBase<LocalDim, EmbedDim>::cell_measure          )                \
      .method("cell_bbox"             , &TriangulationBase<LocalDim, EmbedDim>::cell_bbox             )                \
      .method("cell_barycenter"       , &TriangulationBase<LocalDim, EmbedDim>::cell_barycenter       )                \
      .method("cell_circumcenter"     , &TriangulationBase<LocalDim, EmbedDim>::cell_circumcenter     )                \
      .method("cell_diameter"         , &TriangulationBase<LocalDim, EmbedDim>::cell_diameter         )                \
      .method("quadrature_nodes"      , &TriangulationBase<LocalDim, EmbedDim>::quadrature_nodes      )
  
using cpp_triangulation_2_2 = Triangulation<2, 2>;
RCPP_MODULE(cpp_triangulation_2_2) {
    Rcpp::class_<TriangulationBase<2, 2>>("cpp_triangulation_base").triangulation_rcpp_interface(2, 2);
    Rcpp::class_<Triangulation<2, 2>>("cpp_triangulation_2_2")
      .derives<TriangulationBase<2, 2>>("cpp_triangulation_base")
      .constructor<Rcpp::List>()
      .method("neighbors"                , &Triangulation<2, 2>::neighbors                )
      .method("edges"                    , &Triangulation<2, 2>::edges                    )
      .method("n_edges"                  , &Triangulation<2, 2>::n_edges                  )
      .method("n_boundary_edges"         , &Triangulation<2, 2>::n_boundary_edges         )
      .method("boundary_edges"           , &Triangulation<2, 2>::boundary_edges           )
      .method("mark_boundary"            , &Triangulation<2, 2>::mark_boundary            )
      .method("edges_markers"            , &Triangulation<2, 2>::edges_markers            )
      .method("clear_boundary_markers"   , &Triangulation<2, 2>::clear_boundary_markers   )
      .method("filter_boundary_by_marker", &Triangulation<2, 2>::filter_boundary_by_marker)
      .method("edge_coords"              , &Triangulation<2, 2>::edge_coords              )
      .method("locate"                   , &Triangulation<2, 2>::locate                   );
}

// clang-format on

}   // namespace r
}   // namespace fdapde
