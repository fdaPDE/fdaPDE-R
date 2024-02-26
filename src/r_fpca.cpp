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
#include "headers/r_fpca.h"


using cpp_fpca_spaceonly_base = R_FPCA<SpaceOnly>;
using cpp_fpca_spaceonly = R_FPCA_SpaceOnly;
RCPP_MODULE(cpp_fpca_spaceonly) {
    // inherited methods
    Rcpp::class_<cpp_fpca_spaceonly_base>("cpp_fpca_spaceonly_base")
      // getters
      .method("scores",                   &cpp_fpca_spaceonly_base::scores                 )
      .method("loadings",                 &cpp_fpca_spaceonly_base::loadings               )
      .method("Psi",                      &cpp_fpca_spaceonly_base::Psi                    )
      // setters
      .method("set_data",                 &cpp_fpca_spaceonly_base::set_data               )
      .method("set_spatial_locations",    &cpp_fpca_spaceonly_base::set_spatial_locations  )
      .method("set_npc",                  &cpp_fpca_spaceonly_base::set_npc                )
      .method("set_lambda",               &cpp_fpca_spaceonly_base::set_lambda             )
      .method("set_solver",               &cpp_fpca_spaceonly_base::set_solver             )
      // utilities
      .method("init",                     &cpp_fpca_spaceonly_base::init                   )
      .method("solve",                    &cpp_fpca_spaceonly_base::solve                  );
    // specific methods
    Rcpp::class_<R_FPCA_SpaceOnly>("cpp_fpca_spaceonly")
      .derives<cpp_fpca_spaceonly_base>("cpp_fpca_spaceonly_base")
      .constructor<Rcpp::Environment,
                   int,
                   Rcpp::List>();
}

using cpp_fpca_spacetimeseparable_base = R_FPCA<SpaceTimeSeparable>;
using cpp_fpca_spacetimeseparable = R_FPCA_SpaceTimeSeparable;
RCPP_MODULE(cpp_fpca_spacetimeseparable) {
    // inherited methods
    Rcpp::class_<cpp_fpca_spacetimeseparable_base>("cpp_fpca_spacetimeseparable_base") 
      // getters
      .method("scores",                   &cpp_fpca_spacetimeseparable_base::scores                 )
      .method("loadings",                 &cpp_fpca_spacetimeseparable_base::loadings               )
      .method("Psi",                      &cpp_fpca_spacetimeseparable_base::Psi                    )
      // setters
      .method("set_data",                 &cpp_fpca_spacetimeseparable_base::set_data               )
      .method("set_spatial_locations",    &cpp_fpca_spacetimeseparable_base::set_spatial_locations  )
      .method("set_npc",                  &cpp_fpca_spacetimeseparable_base::set_npc                )
      .method("set_solver",               &cpp_fpca_spacetimeseparable_base::set_solver             )
      // utilities
      .method("init",                     &cpp_fpca_spacetimeseparable_base::init                   )
      .method("solve",                    &cpp_fpca_spacetimeseparable_base::solve                  );
    // specific methods
    Rcpp::class_<R_FPCA_SpaceTimeSeparable>("cpp_fpca_spacetimeseparable")
      .derives<cpp_fpca_spacetimeseparable_base>("cpp_fpca_spacetimeseparable_base")
      .constructor<Rcpp::Environment,
                   Rcpp::Environment,
                   int,
                   Rcpp::List>();
}