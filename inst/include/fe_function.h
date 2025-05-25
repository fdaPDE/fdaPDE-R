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

#ifndef __R_FE_FUNCTION_H__
#define __R_FE_FUNCTION_H__

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/finite_elements.h>
#include "utility.h"

namespace fdapde {
namespace r {

template <int LocalDim, int EmbedDim, typename FeType> class FeFunction {
    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;
    using TriangulationType = fdapde::Triangulation<local_dim, embed_dim>;
    using FeSpaceType = fdapde::FeSpace<TriangulationType, FeType>;
    using DofHandlerType = typename FeSpaceType::DofHandlerType;

    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
   public:
    FeFunction() noexcept = default;
    FeFunction(const Rcpp::Environment& triangulation) :
        fe_space_(get_env_as<TriangulationType>(triangulation), FeType {}), fe_function_(fe_space_) { }

    // observers
    const vector_t& coeff() const { return fe_function_.coeff(); }
    double eval(const vector_t& p) const { return fe_function_.operator()(p); }
    vector_t grid_eval(const matrix_t& ps) const {
        vector_t res(ps.rows());
        for (int i = 0, n = ps.rows(); i < n; ++i) { res[i] = fe_function_(ps.row(i)); }
        return res;
    }
    double l2_squared_norm() { return fe_function_.l2_squared_norm(); }
    double l2_norm() { return fe_function_.l2_norm(); }
    double h1_squared_norm() const { return fe_function_.h1_squared_norm(); }
    double h1_norm() const { return fe_function_.h1_norm(); }
    // integration
    double cell_integrate_on(int marker) const {
        if (marker == BoundaryAll) {
            return fe_function_.integrate_on(
              fe_space_.triangulation().cells_begin(), fe_space_.triangulation().cells_end());
        } else {
            return fe_function_.integrate_on(
              fe_space_.triangulation().cells_begin(marker), fe_space_.triangulation().cells_end(marker));
        }
    }
    int n_dofs() const { return fe_space_.n_dofs(); }
    // modifiers
    void set_coeff(const vector_t& coeff) { fe_function_.set_coeff(coeff); }
   private:
    FeSpaceType fe_space_;
    fdapde::FeFunction<FeSpaceType> fe_function_;
};

}   // namespace r
}   // namespace fdapde

#endif   // __R_FE_FUNCTION_H__
