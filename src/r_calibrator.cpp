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

using cpp_off_base = R_CALIBRATOR<fdapde::calibration::Off>;
using cpp_off = R_OFF;
RCPP_MODULE(cpp_off) {
    // inherited methods
    Rcpp::class_<cpp_off_base>("cpp_off_base")
      .method("get_calibration_strategy",   &cpp_off_base::get_calibration_strategy      )
      .method("optimum",                    &cpp_off_base::optimum                       );
    // specific methods
    Rcpp::class_<R_OFF>("cpp_off")
      .derives<cpp_off_base>("cpp_off_base")
      .constructor<Rcpp::List>() 
      .method("configure_calibrator",       &R_OFF::configure_calibrator          )
      .method("fit",                        &R_OFF::fit                           );    
}

using cpp_gcv_base = R_CALIBRATOR<fdapde::calibration::GCV<void>>;
using cpp_gcv = R_GCV; 
RCPP_MODULE(cpp_gcv) { 
    // inherited methods
    Rcpp::class_<cpp_gcv_base>("cpp_gcv_base")
      .method("get_calibration_strategy",   &cpp_gcv_base::get_calibration_strategy      )
      .method("optimum",                    &cpp_gcv_base::optimum                       );
    // specific methods
    Rcpp::class_<R_GCV>("cpp_gcv")
      .derives<cpp_gcv_base>("cpp_gcv_base")
      .constructor<Rcpp::List>()
      .method("configure_calibrator",       &R_GCV::configure_calibrator          )
      .method("fit",                        &R_GCV::fit                           )
      .method("gcvs",                       &R_GCV::gcvs                          ) 
      .method("edfs",                       &R_GCV::edfs                          );
}

using cpp_kcv_base = R_CALIBRATOR<fdapde::calibration::KCV>;
using cpp_kcv = R_KCV;
RCPP_MODULE(cpp_kcv) {
    // inherited methods
    Rcpp::class_<cpp_kcv_base>("cpp_kcv_base")
      .method("get_calibration_strategy",   &cpp_kcv_base::get_calibration_strategy      )
      .method("optimum",                    &cpp_kcv_base::optimum                       );
    // specific methods
    Rcpp::class_<R_KCV>("cpp_kcv")
      .derives<cpp_kcv_base>("cpp_kcv_base")
      .constructor<Rcpp::List>()
      .method("configure_calibrator",       &R_KCV::configure_calibrator          )
      .method("fit",                        &R_KCV::fit                           )
      .method("avg_scores",                 &R_KCV::avg_scores                    )
      .method("std_scores",                 &R_KCV::std_scores                    )
      .method("scores",                     &R_KCV::scores                        )
      .method("set_n_folds",                &R_KCV::set_n_folds                   ); 
}
