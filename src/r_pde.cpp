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
#include "headers/r_pde.h"

using cpp_pde_2d_fe1 = R_PDE<2,2,1>;
RCPP_MODULE(cpp_pde_2d_fe1) {
    Rcpp::class_<R_PDE<2,2,1>>("cpp_pde_2d_fe1")
      .constructor<Rcpp::Environment, int, Rcpp::Nullable<Rcpp::List>>()
      .method("get_quadrature_nodes" , &R_PDE<2,2,1>::get_quadrature_nodes )
      .method("get_dofs_coordinates" , &R_PDE<2,2,1>::get_dofs_coordinates )
      .method("mass"                 , &R_PDE<2,2,1>::R0                   )
      .method("stiff"                , &R_PDE<2,2,1>::R1                   )
      .method("force"                , &R_PDE<2,2,1>::u                    )
      .method("set_dirichlet_bc"     , &R_PDE<2,2,1>::set_dirichlet_bc     )
      .method("set_forcing"          , &R_PDE<2,2,1>::set_forcing          )
      .method("set_initial_condition", &R_PDE<2,2,1>::set_initial_condition)
      .method("init"                 , &R_PDE<2,2,1>::init                 );
}
using cpp_pde_2d_fe2 = R_PDE<2,2,2>;
RCPP_MODULE(cpp_pde_2d_fe2) {
    Rcpp::class_<R_PDE<2,2,2>>("cpp_pde_2d_fe2")
      .constructor<Rcpp::Environment, int, Rcpp::Nullable<Rcpp::List>>()
      .method("get_quadrature_nodes" , &R_PDE<2,2,2>::get_quadrature_nodes )
      .method("get_dofs_coordinates" , &R_PDE<2,2,2>::get_dofs_coordinates )
      .method("mass"                 , &R_PDE<2,2,2>::R0                   )
      .method("stiff"                , &R_PDE<2,2,2>::R1                   )
      .method("force"                , &R_PDE<2,2,2>::u                    )
      .method("set_dirichlet_bc"     , &R_PDE<2,2,2>::set_dirichlet_bc     )
      .method("set_forcing"          , &R_PDE<2,2,2>::set_forcing          )
      .method("set_initial_condition", &R_PDE<2,2,2>::set_initial_condition)
      .method("init"                 , &R_PDE<2,2,2>::init                 );
}
