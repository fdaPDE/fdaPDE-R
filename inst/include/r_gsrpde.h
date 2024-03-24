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

#ifndef __R_GSRPDE_H__
#define __R_GSRPDE_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/utils/symbols.h>
#include <fdaPDE/models.h>
#include "r_pde.h"
using fdapde::models::GSRPDE;
using fdapde::models::Distribution;
using fdapde::models::Sampling;
using fdapde::models::GCV;

template <typename RegularizationType>
class R_GSRPDE {
private:
  using ModelType = GSRPDE<RegularizationType>;
  ModelType model_; // statistical model to wrap
  GCV* gcv_ptr_ = nullptr; // model's view pointer to be wrapped as external pointer
  BlockFrame<double, int> data_;
public:
  R_GSRPDE(Rcpp::Environment pde, int sampling_type, std::string distribution) {
    // recover pointer to penalty
    SEXP pdeptr = pde[".pointer"];
    PDEWrapper* ptr = reinterpret_cast<PDEWrapper*>(R_ExternalPtrAddr(pdeptr));
    // instantiate probability distribution object
    Distribution distr;
    if(distribution == "poisson")     distr = fdapde::models::Poisson();
    if(distribution == "gamma"  )     distr = fdapde::models::Gamma();
    if(distribution == "exponential") distr = fdapde::models::Exponential();
    if(distribution == "bernulli")    distr = fdapde::models::Bernulli();
    // set model instance
    model_ = ModelType(ptr->get_pde(), Sampling(sampling_type), distr);
    gcv_ptr_ = new GCV(model_);   // set up pointer to gcv functor
  }
  // setters
  void set_lambda_D(double lambda_D) { model_.set_lambda_D(lambda_D); }
  void set_observations(const DMatrix<double>& y) { data_.template insert<double>(OBSERVATIONS_BLK, y); }
  void set_covariates(const DMatrix<double>& X) { data_.template insert<double>(DESIGN_MATRIX_BLK, X); }
  void set_fpirls_tolerance(double tol) { model_.set_fpirls_tolerance(tol); }
  void set_fpirls_max_iter(std::size_t max_iter) { model_.set_fpirls_max_iter(max_iter); }
  // getters
  const DVector<double>& f() const { return model_.f(); }
  DMatrix<double> fitted() const { return model_.fitted(); }
  // utilities
  void init() {
    model_.set_data(data_);
    model_.init();
  }
  void solve() { model_.solve(); }
  Rcpp::XPtr<GCV> get_gcv() { return Rcpp::XPtr<GCV>(gcv_ptr_); }
  // destructor
  ~R_GSRPDE() {
    delete gcv_ptr_;
  };
};

#endif // __R_GSRPDE_H__
