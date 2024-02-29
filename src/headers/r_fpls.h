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

#ifndef __R_FPLS_H__
#define __R_FPLS_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// utils
#include <fdaPDE/utils/symbols.h>
#include <fdaPDE/models.h>
using fdapde::models::SpaceOnly;
using fdapde::models::SpaceTimeSeparable;

// pde
#include "r_pde.h"

// fpca
#include <fdaPDE/models/functional/fpls.h>
using fdapde::models::FPLS;

// rsvd
#include <fdaPDE/models/functional/regularized_svd.h>
#include "r_rsvd.h"
using fdapde::models::RSVDType;

// generic fPCA model wrapper signature
template<typename RegularizationType> class R_FPLS {
  public:
    using ModelType = std::decay_t<FPLS<RegularizationType>>;
    using SolverType = RSVDType<ModelType>;
    // constructor
    R_FPLS() : calibration_strategy_(Calibration::off) {}
    // getters
    const SpMatrix<double>&  Psi() const { return model_.Psi(); }
    DMatrix<double> B(){ return model_.B(); }
    DMatrix<double> fitted(){ return model_.fitted(); }
    DMatrix<double> reconstructed(){ return model_.reconstructed(); }
    // setters
    void set_data(const Rcpp::List & data) {
      BlockFrame<double, int> df;
      df.template insert<double>(OBSERVATIONS_BLK, Rcpp::as<DMatrix<double>>(data["Y"]));
      df.template insert<double>(DESIGN_MATRIX_BLK, Rcpp::as<DMatrix<double>>(data["X"]));
      model_.set_data(df);
    }
    void set_spatial_locations(const DMatrix<double>& locs) { model_.set_spatial_locations(locs); }
    void set_ncomp(std::size_t n_comp) { model_.set_ncomp(n_comp); }
    void set_lambda(Rcpp::List R_lambda){
      // parse the R input
      Rcpp::NumericVector lambda_D_grid = R_lambda["space"];
      // Rcpp::NumericVector lambda_T_grid = R_lambda["time"];
      // move Rcpp::NumericVector into something understandable by cpp layer
      lambda_grid_.resize(lambda_D_grid.size());
      for(std::size_t i = 0; i < lambda_grid_.size(); ++i){
        lambda_grid_[i].resize(1); // 2
        lambda_grid_[i][0] = lambda_D_grid[i];
        // lambda_grid_[i][1] = lambda_T_grid[i];
      }
    }
    void set_rsvd(int policy, Rcpp::List rsvd_params, int calibration_strategy, Rcpp::List calibrator_params) {
      rsvd_ = R_RSVD<ModelType>(RSVDSolutionPolicy(policy), rsvd_params, Calibration(calibration_strategy), calibrator_params);
      calibration_strategy_ = rsvd_.calibration();
    }
    // utilities
    void init() { 
      load_solver();
      model_.init();
    }
    void solve() { model_.solve(); }
  protected:
    ModelType model_;
    R_RSVD<ModelType> rsvd_;
    Calibration calibration_strategy_;
    std::vector<DVector<double>> lambda_grid_;
    // utilities
    void load_solver() {
      if(calibration_strategy_ == Calibration::off){
        model_.set_rsvd(rsvd_());
        model_.set_lambda_D(lambda_grid_.front()[0]);
      } else {
        model_.set_rsvd(rsvd_(lambda_grid_));
      }
    }
};

// implementation of the fPCA model wrapper for SpaceOnly regulatization
class R_FPLS_SpaceOnly : public R_FPLS<SpaceOnly> {
  public:
    R_FPLS_SpaceOnly(Rcpp::Environment pde,
                     int sampling_type,
                     Rcpp::List fPCA_params) {
        // recover pointer to penalty
        SEXP pdeptr = pde[".pointer"];
        PDEWrapper* ptr = reinterpret_cast<PDEWrapper*>(R_ExternalPtrAddr(pdeptr));
        // set model instance
        model_ = ModelType(ptr->get_pde(), Sampling(sampling_type));
        // model configuration
        model_.set_ncomp(fPCA_params["n_comp"]);
    }
};

#endif // __R_FPLS_H__