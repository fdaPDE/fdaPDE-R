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
#include "headers/r_calibrator.h"

// Rcpp modules definition
 
using cpp_off = R_OFF; 
RCPP_MODULE(cpp_off) { 
    Rcpp::class_<R_OFF>("cpp_off")
      .constructor<Rcpp::List>() 
      .method("get_calibration_strategy",   &R_OFF::get_calibration_strategy      )
      .method("set_lambda",                 &R_OFF::set_lambda                    )
      .method("fit",                        &R_OFF::fit                           );    
}
 
using cpp_gcv = R_GCV; 
RCPP_MODULE(cpp_gcv) { 
    Rcpp::class_<R_GCV>("cpp_gcv")
      .constructor<Rcpp::List>()
      .method("get_calibration_strategy",   &R_GCV::get_calibration_strategy      ) 
      .method("set_lambda",                 &R_GCV::set_lambda                    )
      .method("fit",                        &R_GCV::fit                           )
      .method("gcvs",                       &R_GCV::gcvs                          ) 
      .method("edfs",                       &R_GCV::edfs                          );
}

using cpp_kcv = R_KCV;
RCPP_MODULE(cpp_kcv) {
    Rcpp::class_<R_KCV>("cpp_kcv")
      .constructor<Rcpp::List>()
      .method("get_calibration_strategy",   &R_KCV::get_calibration_strategy      )
      .method("set_lambda",                 &R_KCV::set_lambda                    )
      .method("fit",                        &R_KCV::fit                           )
      .method("avg_scores",                 &R_KCV::avg_scores                    )
      .method("std_scores",                 &R_KCV::std_scores                    )
      .method("scores",                     &R_KCV::scores                        )
      .method("optimum",                    &R_KCV::optimum                       )
      .method("set_n_folds",                &R_KCV::set_n_folds                   ); 
}
