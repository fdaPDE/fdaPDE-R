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

#include "../inst/include/fpca.h"

namespace fdapde {
namespace r {
  
using cpp_fpca_space = FPCA<models::SpaceOnly>;
RCPP_MODULE(cpp_fpca_space) {
    Rcpp::class_<FPCA<models::SpaceOnly>>("cpp_fpca_space")
      .constructor<Rcpp::Environment, int>()
      .method("fit"                , &FPCA<models::SpaceOnly>::fit              )
      .method("loadings"           , &FPCA<models::SpaceOnly>::loadings         )
      .method("scores"             , &FPCA<models::SpaceOnly>::scores           )
      .method("set_observations"   , &FPCA<models::SpaceOnly>::set_observations );
};

using cpp_fpca_spacetime = FPCA<models::SpaceTimeSeparable>;
RCPP_MODULE(cpp_fpca_spacetime) {
    Rcpp::class_<FPCA<models::SpaceTimeSeparable>>("cpp_fpca_spacetime")
      .constructor<Rcpp::Environment, Rcpp::Environment, int>()
      .method("fit"                , &FPCA<models::SpaceTimeSeparable>::fit              )
      .method("loadings"           , &FPCA<models::SpaceTimeSeparable>::loadings         )
      .method("scores"             , &FPCA<models::SpaceTimeSeparable>::scores           )
      .method("set_observations"   , &FPCA<models::SpaceTimeSeparable>::set_observations );
};

  
}   // namespace r
}   // namespace fdapde
