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
#include "../inst/include/sr.h"
#include "../inst/include/gsr.h"
#include "../inst/include/qsr.h"

namespace fdapde {
namespace r {

// clang-format on
  
#define fe_ls_elliptic_rcpp_interface(LocalDim, EmbedDim, Model)                                                       \
     method("fit"    , &fe_ls_elliptic<LocalDim, EmbedDim, Model>::fit    )                                            \
    .method("fit_gcv", &fe_ls_elliptic<LocalDim, EmbedDim, Model>::fit_gcv)                                            \
    .method("f"      , &fe_ls_elliptic<LocalDim, EmbedDim, Model>::f      )                                            \
    .method("beta"   , &fe_ls_elliptic<LocalDim, EmbedDim, Model>::beta   )                                            \
    .method("fitted" , &fe_ls_elliptic<LocalDim, EmbedDim, Model>::fitted )

// spatial regression
using cpp_sr_2_2 = sr_elliptic<2, 2>;
RCPP_MODULE(cpp_sr_2_2) {
    Rcpp::class_<fe_ls_elliptic<2, 2, fdapde::SRPDE<internals::fe_ls_elliptic>>>("cpp_fe_ls_elliptic")
      .fe_ls_elliptic_rcpp_interface(2, 2, fdapde::SRPDE<internals::fe_ls_elliptic>);
    Rcpp::class_<sr_elliptic<2, 2>>("cpp_sr_2_2")
      .derives<fe_ls_elliptic<2, 2, fdapde::SRPDE<internals::fe_ls_elliptic>>>("cpp_fe_ls_elliptic")
      .constructor<std::string, Rcpp::Environment, Rcpp::Nullable<Rcpp::List>>();
}

// generalized regression
using cpp_gsr_2_2 = gsr_elliptic<2, 2>;
RCPP_MODULE(cpp_gsr_2_2) {
    Rcpp::class_<fe_ls_elliptic<2, 2, fdapde::GSRPDE<internals::fe_ls_elliptic>>>("cpp_fe_ls_elliptic")
      .fe_ls_elliptic_rcpp_interface(2, 2, fdapde::GSRPDE<internals::fe_ls_elliptic>);
    Rcpp::class_<gsr_elliptic<2, 2>>("cpp_gsr_2_2")
      .derives<fe_ls_elliptic<2, 2, fdapde::GSRPDE<internals::fe_ls_elliptic>>>("cpp_fe_ls_elliptic")
      .constructor<std::string, Rcpp::Environment, std::string, Rcpp::Nullable<Rcpp::List>>();
}

// quantile regression
using cpp_qsr_2_2 = qsr_elliptic<2, 2>;
RCPP_MODULE(cpp_qsr_2_2) {
    Rcpp::class_<fe_ls_elliptic<2, 2, fdapde::QSRPDE<internals::fe_ls_elliptic>>>("cpp_fe_ls_elliptic")
      .fe_ls_elliptic_rcpp_interface(2, 2, fdapde::QSRPDE<internals::fe_ls_elliptic>);
    Rcpp::class_<qsr_elliptic<2, 2>>("cpp_qsr_2_2")
      .derives<fe_ls_elliptic<2, 2, fdapde::QSRPDE<internals::fe_ls_elliptic>>>("cpp_fe_ls_elliptic")
      .constructor<std::string, Rcpp::Environment, double, Rcpp::Nullable<Rcpp::List>>();
}
  
// clang-format on
  
}   // namespace r
}   // namespace fdapde
