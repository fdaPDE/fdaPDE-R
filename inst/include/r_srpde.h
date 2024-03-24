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

#ifndef __R_SRPDE_H__
#define __R_SRPDE_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/utils/symbols.h>
#include <fdaPDE/models.h>
#include "r_pde.h"
using fdapde::models::SRPDE;
using fdapde::models::Sampling;
using fdapde::models::RegressionView;
using fdapde::models::GCV;

namespace fdapde {
namespace r {
  
class SRPDE {
   private:
    SRPDE model_;              // statistical model to wrap
    GCV* gcv_ptr_ = nullptr;   // model's view pointer to be wrapped as external pointer
    BlockFrame<double, int> data_;
   public:
    R_SRPDE(Rcpp::Environment pde, int sampling_type) {
        // recover pointer to penalty
        SEXP pdeptr = pde[".pointer"];
        PDEWrapper* ptr = reinterpret_cast<PDEWrapper*>(R_ExternalPtrAddr(pdeptr));
        // set model instance
        model_ = SRPDE(ptr->get_pde(), Sampling(sampling_type));
        gcv_ptr_ = new GCV(model_);   // set up pointer to gcv functor
    }
    // setters
    void set_lambda_D(double lambda_D) { model_.set_lambda_D(lambda_D); }
    void set_observations(const DMatrix<double>& y) { data_.template insert<double>(OBSERVATIONS_BLK, y); }
    void set_covariates(const DMatrix<double>& X) { data_.template insert<double>(DESIGN_MATRIX_BLK, X); }
    // getters
    const DVector<double>& f() const { return model_.f(); }
    DMatrix<double> fitted() const { return model_.fitted(); }
    // utilities
    void init() {
        model_.set_data(data_);
        model_.init();
    }
    void fit() { model_.solve(); }
    Rcpp::XPtr<GCV> get_gcv() { return Rcpp::XPtr<GCV>(gcv_ptr_); }
    // destructor
    ~R_SRPDE() { delete gcv_ptr_; };
};

}   // namespace r
}   // namespace fdapde

#endif // __R_SRPDE_H__
