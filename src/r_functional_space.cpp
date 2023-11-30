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
#include "headers/r_functional_space.h"

using cpp_lagrange_basis_2d_fe1 = R_FunctionalSpace<2,2,1>;
RCPP_MODULE(cpp_lagrange_basis_2d_fe1) {
    Rcpp::class_<R_FunctionalSpace<2,2,1>>("cpp_lagrange_basis_2d_fe1")
      .constructor<Rcpp::Environment, int>()
      .method("size"      , &R_FunctionalSpace<2,2,1>::size     )
      .method("eval"      , &R_FunctionalSpace<2,2,1>::eval     )
      .method("integrate" , &R_FunctionalSpace<2,2,1>::integrate);
}

using cpp_lagrange_basis_2d_fe2 = R_FunctionalSpace<2,2,2>;
RCPP_MODULE(cpp_lagrange_basis_2d_fe2) {
    Rcpp::class_<R_FunctionalSpace<2,2,2>>("cpp_lagrange_basis_2d_fe2")
      .constructor<Rcpp::Environment, int>()
      .method("size"      , &R_FunctionalSpace<2,2,2>::size     )
      .method("eval"      , &R_FunctionalSpace<2,2,2>::eval     )
      .method("integrate" , &R_FunctionalSpace<2,2,2>::integrate);
}
