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

#ifndef __R_GCV_H__
#define __R_GCV_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/utils/symbols.h>
#include <fdaPDE/optimization.h>
#include <fdaPDE/models.h>
#include "r_pde.h"
using fdapde::models::EDFStrategy;
using fdapde::models::StochasticEDF;
using fdapde::models::ExactEDF;
using fdapde::models::GCV;

class R_GCV {
private:
  GCV* gcv_; // gcv functor
public:
  R_GCV(Rcpp::XPtr<GCV> gcv) {
    // recover pointer to model's GCV functor
    gcv_ = gcv.get();
  }
  // calibrate model_ using a gcv approach, customized with the given parameters
  void calibrate(Rcpp::List gcv_params) {
    // parse gcv_params list
    EDFStrategy edf;
    std::string edf_type = gcv_params["edf_computation"];
    if(edf_type == "exact") {
      edf = ExactEDF();
    } else {
      StochasticEDF edf_{};
      if (gcv_params["seed"] != R_NilValue) edf_.set_seed(Rcpp::as<int>(gcv_params["seed"]));
      if (gcv_params["mc_samples"] != R_NilValue) edf_.set_n_mc_samples(Rcpp::as<int>(gcv_params["mc_samples"]));
      edf = edf_;
    }
    gcv_->set_edf_strategy(edf);
    // customize optimization
    if(gcv_params.containsElementNamed("optimizer")) {
      std::string opt_name = gcv_params["optimizer"];
      if(opt_name == "grid") {
	fdapde::core::Grid<fdapde::Dynamic> opt;
	// grid of lambda values assumed to be present (must be guaranteed by R wrapper)
	Rcpp::NumericVector lambda_grid = gcv_params["lambda"];
	// move Rcpp::NumericVector into something understandable by cpp layer
	std::vector<DVector<double>> cpp_grid;
	cpp_grid.resize(lambda_grid.size());
	for(std::size_t i = 0; i < cpp_grid.size(); ++i){
	  cpp_grid[i].resize(1);
	  cpp_grid[i][0] = lambda_grid[i];
	}
	// optimize
	opt.optimize(*gcv_, cpp_grid);
	return;
      }
      // iterative optimization methods
      int max_iter = 10, tolerance = 1e-2, step = 1e-1; // default values
      if (gcv_params["max_iter"]  != R_NilValue) max_iter = Rcpp::as<int>(gcv_params["max_iter"]);
      if (gcv_params["tolerance"] != R_NilValue) tolerance = Rcpp::as<double>(gcv_params["tolerance"]);
      if (gcv_params["step"]      != R_NilValue) step = Rcpp::as<double>(gcv_params["step"]);
      // optimization starting point
      DVector<double> lambda;
      lambda.resize(1);
      lambda[0] = Rcpp::as<double>(gcv_params["lambda"]);
      // iterative optimizers
      if (opt_name == "newton") {
	fdapde::core::Newton<fdapde::Dynamic> opt(max_iter, tolerance, step);
	opt.optimize(*gcv_, lambda);
	return;
      }
      if (opt_name == "gradient_descent") {
	fdapde::core::GradientDescent<fdapde::Dynamic> opt(max_iter, tolerance, step);
	opt.optimize(*gcv_, lambda);
	return;
      }
      if (opt_name == "bfgs") {
	fdapde::core::BFGS<fdapde::Dynamic> opt(max_iter, tolerance, step);
	opt.optimize(*gcv_, lambda);
	return;
      }      
    } else {
      // default optimization
    }
  }
  // getters
  std::vector<double> gcvs() const { return gcv_->gcvs(); }
  std::vector<double> edfs() const { return gcv_->edfs(); }
  
  // destructor
  ~R_GCV() = default;
};

#endif // __R_SRPDE_H__
