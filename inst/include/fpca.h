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

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include "fe_ls_elliptic.h"

namespace fdapde {
namespace r {

template <int LocalDim, int EmbedDim> class fpca_laplace {
    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;

    using Triangulation = fdapde::Triangulation<local_dim, embed_dim>;
    using GeoFrame = fdapde::GeoFrame<Triangulation>;
    using Model = fdapde::fPCA<internals::fe_ls_elliptic>;
   public:
    fpca_laplace() noexcept = default;
    fpca_laplace(const std::string& colname, const Rcpp::Environment& geoframe) {
        const GeoFrame& gf = get_env_as<GeoFrame>(geoframe);
        const Triangulation& D = get_env_as<GeoFrame>(geoframe).template triangulation<0>();
        FeSpace Vh(D, P1<1>);

        // discretize
        TrialFunction f(Vh);
        TestFunction  v(Vh);
        auto a = integral(D)(dot(grad(f), grad(v)));
        ScalarField<local_dim, decltype([](const vector_t&) { return 0; })> u;
        auto F = integral(D)(u * v);

        model_.discretize(std::pair {a, F});
	model_.analyze_data(colname, gf);
    }
    void fit(int rank, const Rcpp::List& params) {
        std::vector<double> lambda_grid = params["grid"];
        model_.fit(rank, lambda_grid, ComputeRandSVD);
    }
    // observers
    const matrix_t& S() const { return model_.S(); }   // scoring matrix
    const matrix_t& F() const { return model_.F(); }   // loading matrix
    matrix_t Fn() const { return model_.Fn(); }
    const std::vector<double>& loadings_norm() const { return model_.loadings_norm(); }
    const matrix_t& lambda() const { return model_.lambda(); }
   private:
    Model model_;
};

}   // namespace r
}   // namespace fdapde

#endif   // __R_FPCA_H__
