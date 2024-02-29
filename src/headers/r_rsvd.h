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

#ifndef __R_RSVD_H__
#define __R_RSVD_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// utils
#include <fdaPDE/utils/symbols.h>
#include <fdaPDE/models.h>

// calibration
#include <fdapde/calibration/symbols.h>
using fdapde::calibration::Calibration;

// RSVD
#include <fdaPDE/models/functional/regularized_svd.h>
using fdapde::models::RegularizedSVD;
using fdapde::models::RSVDType;

enum RSVDSolutionPolicy{ sequential, monolithic };

template<typename ModelType> 
class R_RSVD {
  private:
    RegularizedSVD<fdapde::sequential> rsvd_sequential_;
    RegularizedSVD<fdapde::monolithic> rsvd_monolithic_;
    RSVDType<ModelType> type_erased_rsvd_;
    RSVDSolutionPolicy policy_;
  public:
    R_RSVD() : policy_(RSVDSolutionPolicy::sequential) {}
    R_RSVD(RSVDSolutionPolicy policy, Rcpp::List sequential_params, Calibration calibration_strategy, Rcpp::List calibrator_params) : policy_(policy) {
      if (policy_ == RSVDSolutionPolicy::sequential) {
        // rsvd and calibrator initialization
        rsvd_sequential_ = RegularizedSVD<fdapde::sequential>{calibration_strategy};
        if(calibration_strategy == Calibration::kcv){
          rsvd_sequential_.set_nfolds(calibrator_params["n_folds"]);
        }
        if(calibration_strategy != Calibration::off ){
          if(calibrator_params["seed"] != R_NilValue) rsvd_sequential_.set_seed(calibrator_params["seed"]);
        }
        // rsvd configuration
        rsvd_sequential_.set_tolerance(Rcpp::as<double>(sequential_params["tolerance"]));
        rsvd_sequential_.set_max_iter(sequential_params["max_iter"]);
        type_erased_rsvd_ = RSVDType<ModelType>(rsvd_sequential_);
      } else {
        type_erased_rsvd_ = RSVDType<ModelType>(rsvd_monolithic_);
      }
    }
    Calibration calibration() { return type_erased_rsvd_.calibration(); }
    RSVDType<ModelType> operator()(std::vector<DVector<double>> lambda_grid_) {
      switch (policy_) {
        case RSVDSolutionPolicy::sequential : {
          type_erased_rsvd_ = RSVDType<ModelType>(rsvd_sequential_.set_lambda(lambda_grid_));
          return type_erased_rsvd_;
        }
        case RSVDSolutionPolicy::monolithic : {
          fdapde_assert("the monolithic solver does not implement a grid based calibration approach.");
        }
      }
    }
    RSVDType<ModelType> operator()() {
      switch (policy_) {
        case RSVDSolutionPolicy::sequential : {
          type_erased_rsvd_ = RSVDType<ModelType>(rsvd_sequential_);
          return type_erased_rsvd_;
        }
        case RSVDSolutionPolicy::monolithic : {
          type_erased_rsvd_ = RSVDType<ModelType>(rsvd_monolithic_);
          return type_erased_rsvd_;
        }
      }
    }
};

#endif // __R_RSVD_H__