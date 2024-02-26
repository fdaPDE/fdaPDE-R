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

#ifndef __R_CENTER_H__
#define __R_CENTER_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// utils
#include <fdaPDE/utils/symbols.h>
#include "r_pde.h"

// models
#include <fdaPDE/models.h>
using fdapde::models::Sampling;
using fdapde::models::RegressionView;

// calibrators
#include <fdaPDE/calibration/calibration_base.h>
using fdapde::calibration::Calibrator;

// center method
#include <fdaPDE/models/functional/center.h>
using fdapde::models::CenterReturnType;


// implementation of CENTER wrapper
class R_CENTER {

private:
    RegressionView<void> smoother_view_;
    Calibrator<RegressionView<void>> calibrator_;
    DMatrix<double> X_;
    DMatrix<double> w_;
    bool weighted_ = FALSE;
    CenterReturnType results_;

public:

    // constructor
    R_CENTER() = default;

    // setters
    void set_data(const DMatrix<double>& X) { X_ = X; }
    void set_weights(const DVector<double>& w) {
      w_ = w;
      weighted_ = TRUE;
    }
    void set_smoother(Rcpp::XPtr<RegressionView<void>> smoother_view_ptr){
      smoother_view_ = *smoother_view_ptr.get();
    }
    void set_calibrator(Rcpp::XPtr<Calibrator<RegressionView<void>>> calibrator_ptr){
      calibrator_ = *calibrator_ptr.get();
    }

    // utilities
    void init() { return; }
    void solve() {
      if(weighted_){
        results_ = center(X_, w_, smoother_view_, calibrator_);
      }
      else {
        results_ = center(X_, smoother_view_, calibrator_);
      }
    }

    // getters
    DMatrix<double> centered() {return results_.fitted; }
    DMatrix<double> mean() {return results_.mean; }

};

#endif // __R_CENTER_H__
