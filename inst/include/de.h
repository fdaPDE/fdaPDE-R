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

#ifndef __R_FE_DE_H__
#define __R_FE_DE_H__

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/models.h>
#include "utility.h"

namespace fdapde {
namespace r {

template <int LocalDim, int EmbedDim> class de_elliptic {
    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;

    using Triangulation = fdapde::Triangulation<local_dim, embed_dim>;
    using GeoFrame = fdapde::GeoFrame<Triangulation>;
    using Model = fdapde::DEPDE<internals::fe_de_elliptic>;
   public:
    de_elliptic() noexcept = default;
    de_elliptic(
      const Rcpp::Environment& geoframe, const Rcpp::Nullable<Rcpp::List>& penalty) {
        const GeoFrame& gf = get_env_as<GeoFrame>(geoframe);
	const Triangulation& D = get_env_as<GeoFrame>(geoframe).template triangulation<0>();
	FeSpace Vh(D, P1<1>);
	n_dofs_ = Vh.n_dofs();
	measure_ = D.measure();
	
	// discretize
	TrialFunction f(Vh);
        TestFunction  v(Vh);
        if (penalty.isNotNull()) {   // general elliptic operator
            Rcpp::List ls(penalty);
            // bilinear form
            FeCoeff<local_dim, local_dim, local_dim, matrix_t> K(Rcpp::as<matrix_t>(ls["K"]));
            FeCoeff<local_dim, local_dim, 1, matrix_t> b(Rcpp::as<matrix_t>(ls["b"]));
            FeCoeff<local_dim, 1, 1, vector_t> c(Rcpp::as<matrix_t>(ls["c"]));
            auto a = integral(D)(dot(K * grad(f), grad(v)) + dot(b, grad(f)) * v + c * f * v);
            // linear form
            FeCoeff<local_dim, 1, 1, vector_t> u(Rcpp::as<matrix_t>(ls["u"]));
            auto F = integral(D)(u * v);

            model_.discretize(gf, std::pair {a, F});
        } else {   // fallback to isotropic laplacian penalty
            auto a = integral(D)(dot(grad(f), grad(v)));
            ScalarField<local_dim, decltype([](const vector_t&) { return 0; })> u;
            auto F = integral(D)(u * v);
	    
            model_.discretize(gf, std::pair {a, F});
        }
	model_.analyze_data(gf);
    }

    // fitting
    void fit(double lambda, const std::string& opt_t, const Rcpp::List& params) {
        // unpack optimization parameters
        int max_iter = params["max_iter"];
        double tol = params["tolerance"], step = params["step"];

	vector_t g_init(n_dofs_);
	for(int i = 0; i < n_dofs_; ++i) { g_init[i] = 1.0/measure_; }
	
        if (opt_t == "gradient_descent" || opt_t == "bfgs") {
            if (opt_t == "gradient_descent") {
                model_.fit(lambda, g_init, GradientDescent<Dynamic> {max_iter, tol, step});
            }
            if (opt_t == "bfgs") { model_.fit(lambda, g_init, BFGS<Dynamic> {max_iter, tol, step}); }
        }
        return;
    }
    // observers
    vector_t density() const { return model_.density(); }
    vector_t log_density() const { return model_.log_density(); }
    vector_t fitted() const { return model_.fn(); }
   protected:
    Model model_;
    using edf_cache_t = std::unordered_map<std::array<double, 1>, double, internals::std_array_hash<double, 1>>;
    edf_cache_t edf_cache_;

    double measure_ = 0;   // domain measure, just to set initial density to 1/|D|
    int n_dofs_;
};

}   // namespace r
}   // namespace fdapde

#endif   // __R_FE_DE_H__
