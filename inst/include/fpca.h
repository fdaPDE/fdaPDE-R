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

#ifndef __R_FPCA_H__
#define __R_FPCA_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/functional.h>

#include "fdaPDE/fdaPDE/calibration/calibration_base.h"
#include "fdaPDE/fdaPDE/calibration/gcv.h"
#include "fdaPDE/fdaPDE/calibration/kfold_cv.h"
#include "fdaPDE/fdaPDE/calibration/off.h"
#include "fdaPDE/fdaPDE/models/functional/fpca.h"

#include "pde.h"
#include "utils.h"

namespace fdapde {
namespace r {

template <typename RegularizationType> class FPCA {
   private:
    using DataContainer = BlockFrame<double, int>;
    using ModelType = models::FPCA<RegularizationType>;
    ModelType model_;
    DataContainer data_;
   public:
    // space-only
    FPCA(Rcpp::Environment r_env, int sampling_type)
    requires(models::is_space_only<ModelType>::value) {
        model_ = ModelType(
          get_env_as<r::PDE>(r_env)->pde, Sampling(sampling_type),
          models::RegularizedSVD<fdapde::sequential> {Calibration::off});
    }
    // space-time separable
    FPCA(Rcpp::Environment r_env1, Rcpp::Environment r_env2, int sampling_type)
    requires(models::is_space_time_separable<ModelType>::value) {
        model_ = ModelType(
          get_env_as<r::PDE>(r_env1)->pde, get_env_as<r::PDE>(r_env2)->pde, Sampling(sampling_type),
          models::RegularizedSVD<fdapde::sequential> {Calibration::off});
    }

    void fit(const Rcpp::List& r_input) {
        // set desired resolution strategy
        if (r_input["solver"] == "monolithic") {
            model_.set_solver(models::RegularizedSVD<fdapde::monolithic>());
            model_.set_lambda(Rcpp::as<DVector<double>>(r_input["lambda"]));
        } else {
            double tolerance = Rcpp::as<double>(r_input["tolerance"]);
            int max_iter     = Rcpp::as<int>(r_input["max_iter"]);
            // configure calibration strategy
            std::string calibration = r_input["calibration"];
            if (calibration == "off") {
                models::RegularizedSVD<fdapde::sequential> rsvd(Calibration::off);
                rsvd.set_tolerance(tolerance);
                rsvd.set_max_iter(max_iter);
                model_.set_solver(rsvd);
                model_.set_lambda(Rcpp::as<DVector<double>>(r_input["lambda"]));
            }
            if (calibration == "gcv") {
                models::RegularizedSVD<fdapde::sequential> rsvd(Calibration::gcv);
                rsvd.set_tolerance(tolerance);
                rsvd.set_max_iter(max_iter);
                rsvd.set_seed(Rcpp::as<int>(r_input["seed"]));
                rsvd.set_lambda(Rcpp::as<DMatrix<double>>(r_input["lambda"]));
                model_.set_solver(rsvd);
            }
            if (calibration == "kcv") {
                models::RegularizedSVD<fdapde::sequential> rsvd(Calibration::kcv);
                rsvd.set_tolerance(tolerance);
                rsvd.set_max_iter(max_iter);
                rsvd.set_seed(Rcpp::as<int>(r_input["seed"]));
                rsvd.set_nfolds(Rcpp::as<int>(r_input["n_folds"]));
                rsvd.set_lambda(Rcpp::as<DMatrix<double>>(r_input["lambda"]));
                model_.set_solver(rsvd);
            }
        }
        // set number of components
        model_.set_npc(Rcpp::as<int>(r_input["ncomp"]));
        model_.set_data(data_);
        model_.init();
        model_.solve();
    }
    // getters
    const DMatrix<double>& loadings() const { return model_.loadings(); }
    const DMatrix<double>& scores() const { return model_.scores(); }
    // setters
    void set_observations(const DMatrix<double>& y) { data_.template insert<double>(OBSERVATIONS_BLK, y); }

    // destructor
    ~FPCA() = default;
};

}   // namespace r
}   // namespace fdapde

#endif   // __R_FPCA_H__
