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
#include "../inst/include/fe_function.h"

namespace fdapde {
namespace r {

// clang-format off
    
using fe_function_2_2_p1 = FeFunction<2, 2, FeP<1, 1>>;
RCPP_MODULE(cpp_fe_function_2_2_p1) {
    Rcpp::class_<FeFunction<2, 2, FeP<1, 1>>>("cpp_fe_function_2_2_p1")
      .constructor<Rcpp::Environment>()
      .method("coeff"            , &FeFunction<2, 2, FeP<1, 1>>::coeff            )
      .method("eval"             , &FeFunction<2, 2, FeP<1, 1>>::eval             )
      .method("grid_eval"        , &FeFunction<2, 2, FeP<1, 1>>::grid_eval        )
      .method("n_dofs"           , &FeFunction<2, 2, FeP<1, 1>>::n_dofs           )
      .method("l2_squared_norm"  , &FeFunction<2, 2, FeP<1, 1>>::l2_squared_norm  )
      .method("h1_squared_norm"  , &FeFunction<2, 2, FeP<1, 1>>::h1_squared_norm  )
      .method("l2_norm"          , &FeFunction<2, 2, FeP<1, 1>>::l2_norm          )
      .method("h1_norm"          , &FeFunction<2, 2, FeP<1, 1>>::h1_norm          )
      .method("cell_integrate_on", &FeFunction<2, 2, FeP<1, 1>>::cell_integrate_on)
      .method("set_coeff"        , &FeFunction<2, 2, FeP<1, 1>>::set_coeff        );
}
    
// clang-format on

}   // namespace r
}   // namespace fdapde
