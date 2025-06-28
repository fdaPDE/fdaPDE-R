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
#include "../inst/include/fe_space.h"

namespace fdapde {
namespace r {

// clang-format off
    
using fe_space_2_2_p1 = FeSpace<2, 2, FeP<1, 1>>;
RCPP_MODULE(cpp_fe_space_2_2_p1) {
    Rcpp::class_<FeSpace<2, 2, FeP<1, 1>>>("cpp_fe_space_2_2_p1")
      .constructor<Rcpp::Environment>()
      .method("eval"        , &FeSpace<2, 2, FeP<1, 1>>::eval        )
      .method("n_dofs"      , &FeSpace<2, 2, FeP<1, 1>>::n_dofs      )
      .method("n_quad_nodes", &FeSpace<2, 2, FeP<1, 1>>::n_quad_nodes)
      .method("quad_nodes"  , &FeSpace<2, 2, FeP<1, 1>>::quad_nodes  );
}
    
// clang-format on

}   // namespace r
}   // namespace fdapde
