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


#ifndef __R_CALIBRATOR_H__
#define __R_CALIBRATOR_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// CalibrationStrategy enumerator
#include <fdaPDE/calibration/symbols.h>
using CalibrationStrategy = fdapde::calibration::Calibration;

// calibrators classes
#include <fdaPDE/calibration/off.h> 
#include <fdaPDE/calibration/gcv.h> 
#include <fdaPDE/models/regression/gcv.h> 
#include <fdaPDE/calibration/kfold_cv.h>
#include <fdaPDE/calibration/rmse.h>
using fdapde::calibration::RMSE;
using fdapde::models::GCV;

// optimization
#include <fdaPDE/optimization.h>
using fdapde::models::EDFStrategy;
using fdapde::models::StochasticEDF;
using fdapde::models::ExactEDF;

// RegressionModel
#include <fdaPDE/models/regression/regression_type_erasure.h>
using fdapde::models::RegressionModel;
using ModelType = RegressionModel<void>; // generic RegressionModel

/* TODOs:
// TODO gcv: generalize the GCV to work with other optimizers, by now the optimizer is forced to be opt_grid.
// TODO gcv: generalize the GCV to work with other EDF strategies, by now the EDF strategies is forced to be stochastic.
// TODO kcv: the fit method needs to be generalized, by now the ScoreType is set ot RMSE
// TODO gcv/kcv: the set_lambda method needs to be generalized to space-time
*/

// implementation of calibrator wrapper for Off calibration strategy
class R_OFF{
    typedef fdapde::calibration::Off CalibratorType;
    protected:
        CalibratorType calibrator_;
    public:
        // constructor
        R_OFF() = default;
        R_OFF(Rcpp::List off_params) { return; };
        // getters
        int get_calibration_strategy() { return CalibrationStrategy::off; }
        // setters
        void set_lambda(Rcpp::List R_lambda){
            DVector<double> lambda(2);
            lambda[0] = Rcpp::as<double>(R_lambda["space"]);
            lambda[1] = Rcpp::as<double>(R_lambda["time"]);
            calibrator_(lambda);
        }
        // fit
        DVector<double> fit(Rcpp::XPtr<fdapde::models::RegressionView<void>> model_view_ptr){
            return calibrator_.fit(*model_view_ptr.get());
        }

};

// implementation of calibrator wrapper for GCV calibration strategy
class R_GCV{
    typedef fdapde::calibration::GCV CalibratorType;
    protected:
        CalibratorType calibrator_;
        std::vector<DVector<double>> lambda_grid_;
    public:
        // constructor
        R_GCV(Rcpp::List gcv_params) {
            // grid optimizer (fixed to grid)
            fdapde::core::Grid<fdapde::Dynamic> opt;
            // EDF strategy (set to stochastic)
            StochasticEDF edf{};
            edf.set_seed(Rcpp::as<int>(gcv_params["seed"]));
            edf.set_n_mc_samples(Rcpp::as<int>(gcv_params["mc_samples"]));
            // calibrator
            calibrator_ = fdapde::calibration::GCV {opt, edf};
        }
        // getters
        int get_calibration_strategy() { return CalibrationStrategy::gcv; }
        std::vector<double> gcvs() const { return calibrator_.gcvs(); }
        std::vector<double> edfs() const { return calibrator_.edfs(); }
        // setters
        void set_lambda(Rcpp::List R_lambda){
            Rcpp::NumericVector lambda_D_grid = R_lambda["space"];
            // move Rcpp::NumericVector into something understandable by cpp layer
            lambda_grid_.resize(lambda_D_grid.size());
            for(std::size_t i = 0; i < lambda_grid_.size(); ++i){
              lambda_grid_[i].resize(1);
              lambda_grid_[i][0] = lambda_D_grid[i];
            }
        }
        void set_step(double step) { calibrator_.set_step(step); }
        // fit
        DVector<double> fit(Rcpp::XPtr<GCV> model_gcv_ptr){
            return calibrator_.fit(*model_gcv_ptr.get(), lambda_grid_);
        }
};

// implementation of calibrator wrapper for KCV calibration strategy
class R_KCV {
    typedef fdapde::calibration::KCV CalibratorType;
    protected:
        CalibratorType calibrator_;
        std::vector<DVector<double>> lambda_grid_;
    public:
        // constructor
        R_KCV(Rcpp::List kcv_params) {
            // number of folds
            if (kcv_params["K"] == R_NilValue){
                calibrator_ = fdapde::calibration::KCV {};
                return;
            }
            std::size_t K = kcv_params["K"];
            // folds randomization
            bool shuffle = kcv_params["shuffle"];
            if (kcv_params["seed"] != R_NilValue){
                calibrator_ = fdapde::calibration::KCV{K, kcv_params["seed"], shuffle};
                return;
            } else{
                calibrator_ = fdapde::calibration::KCV {K, shuffle};
                return;
            }
        }
        // getters
        int get_calibration_strategy() { return CalibrationStrategy::kcv; }
        const DVector<double>& avg_scores() const { return calibrator_.avg_scores(); }
        const DVector<double>& std_scores() const { return calibrator_.std_scores(); }
        const DMatrix<double>& scores() const { return calibrator_.scores(); }
        const DVector<double>& optimum() const { return calibrator_.optimum(); }
        // setters
        void set_lambda(Rcpp::List R_lambda){
            Rcpp::NumericVector lambda_D_grid = R_lambda["space"];
            // move Rcpp::NumericVector into something understandable by cpp layer
            lambda_grid_.resize(lambda_D_grid.size());
            for(std::size_t i = 0; i < lambda_grid_.size(); ++i){
              lambda_grid_[i].resize(1);
              lambda_grid_[i][0] = lambda_D_grid[i];
            }
        }
        void set_n_folds(std::size_t K) { calibrator_.set_n_folds(K); }
        // fit
        DVector<double> fit(Rcpp::XPtr<fdapde::models::RegressionView<void>> model_view_ptr){
            return calibrator_.fit(*model_view_ptr.get(), lambda_grid_, RMSE(*model_view_ptr.get()));
        }
};

#endif   // __R_CALIBRATOR_H__