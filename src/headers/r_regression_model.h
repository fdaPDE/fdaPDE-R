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

#ifndef __R_REGRESSION_MODEL_H__
#define __R_REGRESSION_MODEL_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// pde utilities
#include "r_pde.h"

// RegressionModels
#include <fdaPDE/models/regression/regression_base.h>
#include <fdaPDE/models/regression/srpde.h>
#include <fdaPDE/models/regression/strpde.h>
#include <fdaPDE/models/regression/gsrpde.h>

using fdapde::models::RegressionView;
using fdapde::models::SRPDE;
using fdapde::models::STRPDE;
using fdapde::models::GSRPDE;

// implementation of the wrapper for SRPDE regression model
class R_SRPDE {
    private:
        // statistical model to wrap 
        SRPDE model_;
        // data
        BlockFrame<double, int> data_;
        // calibration utilities
        RegressionView<void> model_view_;
    public:
        // constructor
        R_SRPDE(Rcpp::Environment pde, int sampling_type, const Rcpp::List & smoother_params) {
            // recover pointer to penalty
            SEXP pdeptr = pde[".pointer"];
            PDEWrapper* ptr = reinterpret_cast<PDEWrapper*>(R_ExternalPtrAddr(pdeptr));
            // set model instance
            model_ = SRPDE(ptr->get_pde(), Sampling(sampling_type));
            // smoother customization
            // ... SRPDE does not have any smoother parameter yet
            // create a model view
            model_view_ = model_;
        }
        // setters
        void set_lambda(Rcpp::List lambda) {
            model_.set_lambda_D(lambda["space"]);
            // model_.set_lambda_T(lambda["time"]);
        }
        void set_spatial_locations(const DMatrix<double>& locs) { model_.set_spatial_locations(locs); }
        void set_observations(const DMatrix<double>& y) { data_.template insert<double>(OBSERVATIONS_BLK, y); }
        void set_covariates(const DMatrix<double>& X) { data_.template insert<double>(DESIGN_MATRIX_BLK, X); }
        // getters
        Rcpp::XPtr<RegressionView<void>> get_view() { return Rcpp::XPtr<RegressionView<void>>(&model_view_); }
        const DVector<double>& f() const { return model_.f(); }
        const DVector<double>& beta() const { return model_.beta(); }
        // utilities
        void init() {
            model_.set_data(data_);
            model_.init();
        }
        void solve() { model_.solve(); }
};

#endif // __R_REGRESSION_MODEL_H__