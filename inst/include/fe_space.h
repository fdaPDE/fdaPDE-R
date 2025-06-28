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

template <int LocalDim, int EmbedDim, typename FeType> class FeSpace {
    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;
    using Triangulation = fdapde::Triangulation<local_dim, embed_dim>;
    using FeSpaceType = fdapde::FeSpace<Triangulation, FeType>;
    using Quadrature = typename FeType::template cell_quadrature_t<local_dim>;
    static constexpr int n_quadrature_nodes = Quadrature::order;

    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
   public:
    FeSpace() noexcept = default;
    FeSpace(const Rcpp::Environment& triangulation) :
        fe_space_(get_env_as<Triangulation>(triangulation), FeType {}) { }

    // observers
    sparse_matrix_t eval(const matrix_t& ps) const {
        TestFunction v(fe_space_);
        auto a = eval_at(ps)(v);
        return a.assemble();
    }
    int n_dofs() const { return fe_space_.n_dofs(); }
    int n_quad_nodes() const { return fe_space_.triangulation().n_cells() * n_quadrature_nodes; }
    matrix_t quad_nodes() { return simplex_quadrature_nodes(fe_space_.triangulation(), Quadrature {}); }
   private:
    FeSpaceType fe_space_;
};

}   // namespace r
}   // namespace fdapde

#endif   // __R_FE_FUNCTION_H__
