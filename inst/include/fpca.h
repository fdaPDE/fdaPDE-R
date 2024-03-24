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

#include <fdaPDE/utils/symbols.h>
#include <fdaPDE/models.h>
#include <fdaPDE/splines.h>
#include <fdaPDE/finite_elements.h>
#include "fdaPDE/fdaPDE/models/functional/fpca.h"
#include "fdaPDE/fdaPDE/calibration/calibration_base.h"
#include "fdaPDE/fdaPDE/calibration/off.h"
#include "fdaPDE/fdaPDE/calibration/gcv.h"
#include "fdaPDE/fdaPDE/calibration/kfold_cv.h"
#include "fdaPDE/fdaPDE/models/functional/center.h"
using fdapde::models::FPCA;
using fdapde::models::Sampling;
using fdapde::calibration::Calibration;
using fdapde::core::advection;
using fdapde::core::diffusion;
using fdapde::core::dt;
using fdapde::core::FEM;
using fdapde::core::fem_order;
using fdapde::core::Integrator;
using fdapde::core::laplacian;
using fdapde::core::Mesh;
using fdapde::core::PDE;
using fdapde::core::pde_ptr;
using fdapde::core::reaction;


#include "fdaPDE/fdaPDE/models/regression/strpde.h"
#include "fdaPDE/fdaPDE/models/sampling_design.h"
using fdapde::models::STRPDE;
using fdapde::models::SpaceTimeSeparable;
using fdapde::models::SpaceTimeParabolic;
using fdapde::models::Sampling;


#include "mesh.h"

namespace fdapde{
  namespace r{

template <int M, int N>
class R_FPCA {
private:
  pde_ptr space_pen_;
  pde_ptr time_pen_;
  core::BlockFrame<double, int> data_;
  models::Sampling sampling_;
  double lambda_D_;
  double lambda_T_;
  DMatrix<double> loadings_;
  DMatrix<double> scores_;
  std::vector<DVector<double>> lambda_grid_;
public:
  R_FPCA(Rcpp::Environment mesh, double a, double b, int n, int sampling_type) {
    // recover pointer to penalty
    using RDomainType = r::Mesh<M, N>;
    SEXP meshptr = mesh[".pointer"];
    RDomainType* ptr = reinterpret_cast<RDomainType*>(R_ExternalPtrAddr(meshptr));
    core::Mesh<M, N> domain = ptr->domain();
    // define time domain
    core::Mesh<1,1> time_mesh(a, b, n);

    sampling_ = models::Sampling(sampling_type);
    
    // define regularizing PDE in space
    auto Ld = -core::laplacian<fdapde::core::FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.n_elements() * 3 * time_mesh.n_nodes(), 1);
    space_pen_ = core::PDE<core::Mesh<M, N>, decltype(Ld), DMatrix<double>, core::FEM, core::fem_order<1>>(domain, Ld, u);

    // define regularizing PDE in time
    auto Lt = -core::bilaplacian<core::SPLINE>();
    time_pen_ = core::PDE<core::Mesh<1, 1>, decltype(Lt), DMatrix<double>, core::SPLINE, core::spline_order<3>>(time_mesh, Lt);
  }
  // setters
  void set_lambda_D(double lambda_D) { lambda_D_ = lambda_D; }
  void set_lambda_T(double lambda_T) { lambda_T_ = lambda_T; }
  void set_observations(const DMatrix<double>& y) { data_.template insert<double>(OBSERVATIONS_BLK, y); }
  // getters
  DMatrix<double> loadings() const { return loadings_; }
  DMatrix<double> scores() const { return scores_; }

  void set_lambda_grid(const std::vector<double>& lambda_s, const std::vector<double>& lambda_t) {
    for(int i = 0; i < lambda_s.size(); ++i) {
      lambda_grid_.push_back(SVector<2>(lambda_s[i], lambda_t[i]));
    }
  }

  // DMatrix<double> center(const DMatrix<double>& X, double lambda_D, double lambda_T) {
  // auto centered_data =
  //   models::center(X, fdapde::models::STRPDE<models::SpaceTimeSeparable, fdapde::monolithic> (space_pen_, time_pen_ sampling_),
  // 		   fdapde::calibration::Off {}(SVector<2>(lambda_D, lambda_T)));
  // return centered_data.fitted;
  // }

  
  // utilities
  void solve() {

    models::FPCA<models::SpaceTimeSeparable> model(space_pen_, time_pen_, sampling_, models::RegularizedSVD<fdapde::sequential> {Calibration::off});

    // RegularizedSVD<fdapde::sequential> rsvd(Calibration::gcv);
    // rsvd.set_lambda(lambda_grid_);
    // rsvd.set_seed(78965);   // for reproducibility purposes in testing
    // models::FPCA<models::SpaceTimeSeparable> model(space_pen, time_pen_, sampling_, rsvd);
    
  model.set_lambda_D(lambda_D_);
  model.set_lambda_T(lambda_T_);
  model.set_data(data_);
  model.init();
  model.solve();
  loadings_ = model.loadings();
  scores_ = model.scores();
  }
};

  }}
    
#endif // __R_FPCA_H__
