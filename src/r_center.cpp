
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
#include "headers/r_center.h"

// Rcpp modules definition

using cpp_center = R_CENTER;
RCPP_MODULE(cpp_center) { 
    Rcpp::class_<R_CENTER>("cpp_center")
      .constructor() 
      .method("centered",        &R_CENTER::centered        ) 
      .method("mean",            &R_CENTER::mean            ) 
      .method("set_data",        &R_CENTER::set_data        )
      .method("set_weights",     &R_CENTER::set_weights     )
      .method("set_smoother",    &R_CENTER::set_smoother    )
      .method("set_calibrator",  &R_CENTER::set_calibrator  )
      .method("init",            &R_CENTER::init            )
      .method("solve",           &R_CENTER::solve           );
}