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

#ifndef __R_FE_LS_ELLIPTIC_H__
#define __R_FE_LS_ELLIPTIC_H__

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/models.h>
#include "utility.h"

namespace fdapde {
namespace r {

template <int LocalDim, int EmbedDim, typename Model> class fe_ls_elliptic {
    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;

    using Triangulation = fdapde::Triangulation<local_dim, embed_dim>;
    using GeoFrame = fdapde::GeoFrame<Triangulation>;
   public:
    fe_ls_elliptic() noexcept = default;
    fe_ls_elliptic(
      const std::string& formula, const Rcpp::Environment& geoframe, const Rcpp::Nullable<Rcpp::List>& penalty) {
        const GeoFrame& gf = get_env_as<GeoFrame>(geoframe);
	const Triangulation& D = get_env_as<GeoFrame>(geoframe).template triangulation<0>();
	FeSpace Vh(D, P1<1>);
	
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

            model_.discretize(std::pair {a, F});
        } else {   // fallback to isotropic laplacian penalty
            auto a = integral(D)(dot(grad(f), grad(v)));
            ScalarField<local_dim, decltype([](const vector_t&) { return 0; })> u;
            auto F = integral(D)(u * v);
	    
            model_.discretize(std::pair {a, F});
        }
        model_.analyze_data(formula, gf);
    }

    // fitting
    void fit(double lambda) { model_.fit(lambda); }
    Rcpp::List fit_gcv(const std::string& opt_t, const Rcpp::List& params) {
        int mc_samples = params["mc_samples"];
        int seed = params["seed"];
	auto gcv = model_.gcv(edf_cache_, mc_samples, seed);
	double optimum;
	std::vector<double> values;
	std::vector<double> points;
        if (opt_t == "grid") {
            // unpack optimization parameters
            std::vector<double> lambda_grid = params["grid"];
	    points = lambda_grid;
	    
            GridOptimizer<1> optimizer;
            optimizer.optimize(gcv, lambda_grid);
	    optimum = optimizer.optimum()[0];
	    values  = optimizer.values();
        }
        edf_cache_.insert(gcv.edf_cache().begin(), gcv.edf_cache().end());
        model_.fit(optimum);
	
        Rcpp::List result = Rcpp::List::create(
          Rcpp::Named("optimum") = optimum, Rcpp::Named("values") = values, Rcpp::Named("points") = points);
        return result;
    }
    // observers
    const vector_t& f() const { return model_.f(); }
    const vector_t& beta() const { return model_.beta(); }
    vector_t fitted() const { return model_.fitted(); }
  
   protected:
    Model model_;
    using edf_cache_t = std::unordered_map<std::array<double, 1>, double, internals::std_array_hash<double, 1>>;
    edf_cache_t edf_cache_;
};

}   // namespace r
}   // namespace fdapde

#endif   // __R_FE_LS_ELLIPTIC_H__
