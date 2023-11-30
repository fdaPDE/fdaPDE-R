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

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/utils/symbols.h>
#include <fdaPDE/optimization.h>
#include <fdaPDE/models.h>
using fdapde::models::SRPDE;
using fdapde::models::GSRPDE;
using fdapde::models::STRPDE;
using fdapde::models::SpaceOnly;
using fdapde::models::SpaceTimeSeparable;
using fdapde::models::Sampling;
using fdapde::models::Distribution;
using fdapde::models::GCV;
#include "headers/r_pde.h"

// class R_SRPDE {
// private:
//   SRPDE model_;
//   BlockFrame<double, int> data_;

//   std::vector<double> edfs_;
//   std::vector<double> gcvs_;
//   DVector<double> optimal_lambda_;
  
// public:
//   R_SRPDE(Rcpp::Environment pde, int sampling_type) {
//     // recover pointer to penalty
//     SEXP pdeptr = pde[".pointer"];
//     PDEWrapper* ptr = reinterpret_cast<PDEWrapper*>(R_ExternalPtrAddr(pdeptr));
//     // set model instance
//     model_ = SRPDE(ptr->get_pde(), Sampling(sampling_type));
//   }

//   // setters
//   void set_lambda_D(double lambda_D) { model_.set_lambda_D(lambda_D); }
//   void set_observations(const DMatrix<double>& y) { data_.template insert<double>(OBSERVATIONS_BLK, y); }
//   void set_covariates(const DMatrix<double>& X) { data_.template insert<double>(DESIGN_MATRIX_BLK, X); }
//   // getters
//   const DVector<double>& f() const { return model_.f(); }

//   void init() {
//     model_.set_data(data_);
//     model_.init();
//   }
//   void solve() { model_.solve(); }

//   RegressionView<SpaceOnly> get_regression_view() { return *this; }
  
//   void gcv(const std::vector<double>& lambdas_, int edf_type) {
//     std::vector<DVector<double>> lambdas;
//     for (double l : lambdas_) lambdas.push_back(SVector<1>(l));

//     // select optimizer
//     fdapde::core::Grid<fdapde::Dynamic> opt;
//     if (edf_type == 0) { // exact optimization
//       auto GCV = model_.gcv<fdapde::models::ExactEDF>();
//       opt.optimize(GCV, lambdas);
//       edfs_ = GCV.edfs_;
//       gcvs_ = GCV.gcvs_;
//       optimal_lambda_ = opt.optimum();
//     } else {
//       auto GCV = model_.gcv<fdapde::models::StochasticEDF>(100, 123); // custom seed and iteration
//       opt.optimize(GCV, lambdas);
//       edfs_ = GCV.edfs_;
//       gcvs_ = GCV.gcvs_;
//       optimal_lambda_ = opt.optimum();
//     }
//     return;
//   }

//   const std::vector<double>& edfs() const { return edfs_; }
//   const std::vector<double>& gcvs() const { return gcvs_; }
//   const DVector<double>& optim_lambda() const { return optimal_lambda_; }
//   // destructor
//   ~R_SRPDE() = default;
// };

// // Rcpp modules definition
// RCPP_MODULE(R_SRPDE) {
//     Rcpp::class_<R_SRPDE>("R_SRPDE")
//       .constructor<Rcpp::Environment, int>()
//       // getters
//       .method("f", &R_SRPDE::f)
//       // setters
//       .method("set_lambda_D", &R_SRPDE::set_lambda_D)
//       .method("set_observations", &R_SRPDE::set_observations)
//       .method("set_covariates", &R_SRPDE::set_covariates)
//       .method("gcv", &R_SRPDE::gcv)
//       .method("edfs", &R_SRPDE::edfs)
//       .method("gcvs", &R_SRPDE::gcvs)
//       .method("optimal_lambda", &R_SRPDE::optim_lambda)
//       // init
//       .method("init",  &R_SRPDE::init)
//       .method("solve", &R_SRPDE::solve);
// }

// class R_GSRPDE {
// private:
//   GSRPDE<SpaceOnly> model_;
//   BlockFrame<double, int> data_;
// public:
//   R_GSRPDE(Rcpp::Environment pde, int sampling_type, std::string distribution) {
//     // recover pointer to penalty
//     SEXP pdeptr = pde[".pointer"];
//     PDEWrapper* ptr = reinterpret_cast<PDEWrapper*>(R_ExternalPtrAddr(pdeptr));
//     // instantiate probability distribution object
//     Distribution distr;
//     if(distribution == "poisson")     distr = fdapde::models::Poisson();
//     if(distribution == "gamma"  )     distr = fdapde::models::Gamma();
//     if(distribution == "exponential") distr = fdapde::models::Exponential();
//     if(distribution == "bernulli")    distr = fdapde::models::Bernulli();
//     // set model instance
//     model_ = GSRPDE<SpaceOnly>(ptr->get_pde(), Sampling(sampling_type), distr);
//   }

//   // setters
//   void set_lambda_D(double lambda_D) { model_.set_lambda_D(lambda_D); }
//   void set_observations(const DMatrix<double>& y) { data_.template insert<double>(OBSERVATIONS_BLK, y); }
//   void set_covariates(const DMatrix<double>& X) { data_.template insert<double>(DESIGN_MATRIX_BLK, X); }
//   void set_fpirls_tolerance(double tol) { model_.set_fpirls_tolerance(tol); }
//   void set_fpirls_max_iter(std::size_t max_iter) { model_.set_fpirls_max_iter(max_iter); }
//   // getters
//   const DVector<double>& f() const { return model_.f(); }
  
//   void init() {
//     model_.set_data(data_);
//     model_.init();
//   }
//   void solve() { model_.solve(); }

//   // destructor
//   ~R_GSRPDE() = default;
// };

// // Rcpp modules definition
// RCPP_MODULE(R_GSRPDE) {
//     Rcpp::class_<R_GSRPDE>("R_GSRPDE")
//       .constructor<Rcpp::Environment, int, std::string>()
//       // getters
//       .method("f", &R_GSRPDE::f)
//       // setters
//       .method("set_lambda_D", &R_GSRPDE::set_lambda_D)
//       .method("set_observations", &R_GSRPDE::set_observations)
//       .method("set_covariates", &R_GSRPDE::set_covariates)
//       .method("set_fpirls_tolerance", &R_GSRPDE::set_fpirls_tolerance)
//       .method("set_fpirls_max_iter", &R_GSRPDE::set_fpirls_max_iter)
//       // init
//       .method("init",  &R_GSRPDE::init)
//       .method("solve", &R_GSRPDE::solve);
// }

class R_STRPDE {
private:
  STRPDE<SpaceTimeSeparable, fdapde::monolithic> model_;
  BlockFrame<double, int> data_;
public:
  R_STRPDE(Rcpp::Environment pde, const DVector<double>& time_mesh, int sampling_type) {
    // recover pointer to penalty
    SEXP pdeptr = pde[".pointer"];
    PDEWrapper* ptr = reinterpret_cast<PDEWrapper*>(R_ExternalPtrAddr(pdeptr));
    // set model instance
    model_ = STRPDE<SpaceTimeSeparable, fdapde::monolithic>(ptr->get_pde(), Sampling(sampling_type), time_mesh);
  }

  // setters
  void set_lambda_D(double lambda_D) { model_.set_lambda_D(lambda_D); }
  void set_lambda_T(double lambda_T) { model_.set_lambda_T(lambda_T); }
  void set_observations(const DMatrix<double>& y) { data_.template insert<double>(OBSERVATIONS_BLK, y); }
  void set_covariates(const DMatrix<double>& X) { data_.template insert<double>(DESIGN_MATRIX_BLK, X); }
  // getters
  const DVector<double>& f() const { return model_.f(); }
  
  void init() {
    model_.set_data(data_);
    model_.init();
  }
  void solve() { model_.solve(); }

  // destructor
  ~R_STRPDE() = default;
};

// Rcpp modules definition
RCPP_MODULE(R_STRPDE) {
    Rcpp::class_<R_STRPDE>("R_STRPDE")
      .constructor<Rcpp::Environment, DVector<double>, int>()
      // getters
      .method("f", &R_STRPDE::f)
      // setters
      .method("set_lambda_D", &R_STRPDE::set_lambda_D)
      .method("set_lambda_T", &R_STRPDE::set_lambda_T)
      .method("set_observations", &R_STRPDE::set_observations)
      .method("set_covariates", &R_STRPDE::set_covariates)
      // init
      .method("init",  &R_STRPDE::init)
      .method("solve", &R_STRPDE::solve);
}
