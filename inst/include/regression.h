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

#ifndef __R_REGRESSION_H__
#define __R_REGRESSION_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <variant>

#include <fdaPDE/regression.h>
#include <fdaPDE/optimization.h>
#include "fdaPDE/fdaPDE/calibration/kfold_cv.h"
#include "fdaPDE/fdaPDE/calibration/gcv.h"
#include "fdaPDE/fdaPDE/calibration/rmse.h"
#include "utils.h"
#include "pde.h"

namespace fdapde {
namespace r {

template <typename Objective> core::Optimizer<Objective> parse_optimizer(const Rcpp::List& r_input) {
    core::Optimizer<Objective> opt;
    if (r_input["optimizer"] == R_NilValue) { return opt; }   // return empty optimizer
    std::string opt_flag = r_input["optimizer"];   // optimizer name
    // newton/quasi-netwon optimization parameters
    std::size_t max_iter = Rcpp::as<int>   (r_input["max_iter" ]);
    double tolerance     = Rcpp::as<double>(r_input["tolerance"]);
    double step          = Rcpp::as<double>(r_input["step"     ]);
    // instantiate optimization algorithms
    if (opt_flag == "grid")   { opt = core::Grid<fdapde::Dynamic> {}; }
    if (opt_flag == "newton") { opt = core::Newton<fdapde::Dynamic> {max_iter, tolerance, step}; }
    if (opt_flag == "gd")     { opt = core::GradientDescent<fdapde::Dynamic> {max_iter, tolerance, step}; }
    if (opt_flag == "bfgs")   { opt = core::BFGS<fdapde::Dynamic> {max_iter, tolerance, step}; }

    // handle step adaptivity in optimizers
    
    return opt;
}
  
template <typename ModelType> class RegressionModel {
   protected:
    using RegularizationType = typename ModelType::RegularizationType;
    using DataContainer = BlockFrame<double, int>;
    ModelType model_;      // wrapped model
    DataContainer data_;   // data storage (can be removed if we initialize the data stack in init())
    // calibration
    DVector<double> optim_lambda_;   // optimal selected level of smoothing
    calibration::GCV<
      std::conditional_t<std::is_same_v<RegularizationType, models::SpaceOnly>, models::SpaceOnly, models::SpaceTime>>
      gcv_;
    calibration::KCV kcv_;

    // calibrate regression model according to provided input
    void calibrate(const Rcpp::List& r_input) {
        // recover calibration type
        std::string calibration = r_input["calibration"];
        if (calibration == "off") { optim_lambda_ = Rcpp::as<DVector<double>>(r_input["lambda"]); }
        if (calibration == "gcv") {
            std::string edf = r_input["edf"];
            if (edf == "exact") {
                gcv_.set_edf_strategy(models::ExactEDF());
            } else {   // stochastic approximation of edf
                int n_mc_samples = Rcpp::as<int>(r_input["n_mc_samples"]);
                int seed = r_input["seed"] == R_NilValue ? fdapde::random_seed : Rcpp::as<int>(r_input["seed"]);
                gcv_.set_edf_strategy(models::StochasticEDF(n_mc_samples, seed));
            }
            gcv_.set_optimizer(parse_optimizer<models::GCV>(r_input));
            optim_lambda_ = gcv_(Rcpp::as<DMatrix<double>>(r_input["lambda"])).fit(model_);
        }
        if (calibration == "kcv") {
            int n_folds = Rcpp::as<int>(r_input["n_folds"]);
            int seed = r_input["seed"] == R_NilValue ? fdapde::random_seed : Rcpp::as<int>(r_input["seed"]);
            kcv_ = calibration::KCV(n_folds, seed);
            optim_lambda_ = kcv_(Rcpp::as<DMatrix<double>>(r_input["lambda"]), calibration::RMSE()).fit(model_);
        }
        model_.set_lambda(optim_lambda_);
        return;
    }
  // TODO:
  // need to supply the spatial locations and the temporal locations
  
   public:
    // constructor
    RegressionModel() = default;
    // space-only, space-time parabolic constructor
    RegressionModel(Rcpp::Environment r_env, int sampling_type)
        requires(models::is_space_only<ModelType>::value || models::is_space_time_parabolic<ModelType>::value) {
        // set model instance
        model_ = ModelType(get_env_as<r::PDE>(r_env)->pde, Sampling(sampling_type));
    }
    // space-time separable constructor
    RegressionModel(Rcpp::Environment r_env1, Rcpp::Environment r_env2, int sampling_type)
        requires(models::is_space_time_separable<ModelType>::value) {
        // set model instance
        model_ =
          ModelType(get_env_as<r::PDE>(r_env1)->pde, get_env_as<r::PDE>(r_env2)->pde, Sampling(sampling_type));
    }

    void fit(const Rcpp::List& r_input) {
        model_.set_data(data_);
        calibrate(r_input);   // find optimal tuning parameter according to desired strategy
        model_.init();        // initialize model with optimal parameters
        model_.solve();       // final fitting
    }
    // setters
    void set_observations(const DMatrix<double>& y) { data_.template insert<double>(OBSERVATIONS_BLK, y); }
    void set_covariates(const DMatrix<double>& X) { data_.template insert<double>(DESIGN_MATRIX_BLK, X); }
    void set_lambda(const DVector<double>& lambda) { model_.set_lambda(lambda); }
    // getters
    const DVector<double>& f() const { return model_.f(); }
    DMatrix<double> fitted() const { return model_.fitted(); }
    // destructor
    ~RegressionModel() = default;
};

struct SRPDE : public RegressionModel<models::SRPDE> {
    SRPDE(Rcpp::Environment pde, int sampling_type) : RegressionModel<models::SRPDE>(pde, sampling_type) {}
};

}   // namespace r
}   // namespace fdapde

#endif   // __R_REGRESSION_H__
