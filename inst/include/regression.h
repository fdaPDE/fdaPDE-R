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

#include <fdaPDE/optimization.h>
#include <fdaPDE/regression.h>

#include <variant>

#include "fdaPDE/fdaPDE/calibration/gcv.h"
#include "fdaPDE/fdaPDE/calibration/kfold_cv.h"
#include "fdaPDE/fdaPDE/calibration/rmse.h"
#include "pde.h"
#include "utils.h"

namespace fdapde {
namespace r {

template <typename Objective> core::Optimizer<Objective> parse_optimizer(const Rcpp::List& r_input) {
    core::Optimizer<Objective> opt;
    if (r_input["optimizer"] == R_NilValue) { return opt; }   // return empty optimizer
    std::string opt_flag = r_input["optimizer"];              // optimizer name

    // newton/quasi-netwon optimization parameters
    std::size_t max_iter = Rcpp::as<std::size_t>(r_input["max_iter"]);
    double tolerance = Rcpp::as<double>(r_input["tolerance"]);
    double step = Rcpp::as<double>(r_input["step"]);
    // instantiate optimization algorithms
    if (opt_flag == "grid") { opt = core::Grid<fdapde::Dynamic> {}; }
    if (opt_flag == "newton") { opt = core::Newton<fdapde::Dynamic> {max_iter, tolerance, step}; }
    if (opt_flag == "gd") { opt = core::GradientDescent<fdapde::Dynamic> {max_iter, tolerance, step}; }
    if (opt_flag == "bfgs") { opt = core::BFGS<fdapde::Dynamic> {max_iter, tolerance, step}; }

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
                int seed = Rcpp::as<int>(r_input["seed"]);
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
    const std::vector<double>& gcvs() const { return gcv_.gcvs(); }
    const std::vector<double>& edfs() const { return gcv_.edfs(); }
    const DVector<double>& optimum() const { return optim_lambda_; }

    // destructor
    ~RegressionModel() = default;
};

struct SRPDE : public RegressionModel<models::SRPDE> {
    using ModelType = models::SRPDE;
    using Base = RegressionModel<ModelType>;
    SRPDE(Rcpp::Environment r_env, int sampling_type) {
        Base::model_ = ModelType(get_env_as<r::PDE>(r_env)->pde, Sampling(sampling_type));
    }
};

template <typename RegularizationType> struct GSRPDE : public RegressionModel<models::GSRPDE<RegularizationType>> {
    using ModelType = models::GSRPDE<RegularizationType>;
    using Base = RegressionModel<ModelType>;
    // space-only, space-time parabolic
    GSRPDE(Rcpp::Environment r_env, int sampling_type, int distribution)
    requires(models::is_space_only<ModelType>::value || models::is_space_time_parabolic<ModelType>::value) {
        init_(distribution, get_env_as<r::PDE>(r_env)->pde, Sampling(sampling_type));
    }
    // space-time separable
    GSRPDE(Rcpp::Environment r_env1, Rcpp::Environment r_env2, int sampling_type, int distribution)
    requires(models::is_space_time_separable<ModelType>::value) {
        init_(distribution, get_env_as<r::PDE>(r_env1)->pde, get_env_as<r::PDE>(r_env2)->pde, Sampling(sampling_type));
    }
    void set_fpirls_tolerance(double tol) { Base::model_.set_fpirls_tolerance(tol); }
    void set_fpirls_max_iter(std::size_t max_iter) { Base::model_.set_fpirls_max_iter(max_iter); }
   private:
    // instantiate model depending on distribution type
    template <typename... Args> void init_(int distribution, Args&&... args) {
        if (distribution == 0) Base::model_ = ModelType(std::forward<Args>(args)..., models::Poisson());
        if (distribution == 1) Base::model_ = ModelType(std::forward<Args>(args)..., models::Bernulli());
        if (distribution == 2) Base::model_ = ModelType(std::forward<Args>(args)..., models::Exponential());
        if (distribution == 3) Base::model_ = ModelType(std::forward<Args>(args)..., models::Gamma());
    }
};

template <typename RegularizationType, typename SolutionPolicy>
struct STRPDE : public RegressionModel<models::STRPDE<RegularizationType, SolutionPolicy>> {
    using ModelType = models::STRPDE<RegularizationType, SolutionPolicy>;
    using Base = RegressionModel<ModelType>;
    // space-time parabolic
    STRPDE(Rcpp::Environment r_env, int sampling_type, int distribution)
    requires(models::is_space_time_parabolic<ModelType>::value) {
        Base::model_ = ModelType(get_env_as<r::PDE>(r_env)->pde, Sampling(sampling_type));
    }
    // space-time separable
    STRPDE(Rcpp::Environment r_env1, Rcpp::Environment r_env2, int sampling_type)
    requires(models::is_space_time_separable<ModelType>::value) {
        Base::model_ =
          ModelType(get_env_as<r::PDE>(r_env1)->pde, get_env_as<r::PDE>(r_env2)->pde, Sampling(sampling_type));
    }
};

}   // namespace r
}   // namespace fdapde

#endif   // __R_REGRESSION_H__
