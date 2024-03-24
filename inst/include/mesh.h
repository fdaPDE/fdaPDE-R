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

#ifndef __R_MESH_H__
#define __R_MESH_H__

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/finite_elements.h>
#include <fdaPDE/geometry.h>
#include <fdaPDE/utils/symbols.h>

namespace fdapde {
namespace r {

// Rcpp wrapper for mesh. M : local dimension, N : embedding dimension
template <int M, int N> class Mesh {
   private:
    using NeighborsContainerType = typename core::Mesh<M, N>::NeighborsContainerType;
    using NodesType = typename std::conditional_t<M == 1 && N == 1, DVector<double>, DMatrix<double>>;
    core::Mesh<M, N> domain_ {};
   public:
    // constructor
    Mesh(const Rcpp::List& mesh_data) {
        if constexpr (M == 1 && N == 1) {
            domain_ = core::Mesh<1, 1>(Rcpp::as<DMatrix<double>>(mesh_data["nodes"]));
        } else {
            domain_ = core::Mesh<M, N>(
              Rcpp::as<DMatrix<double>>(mesh_data["nodes"]), Rcpp::as<DMatrix<int>>(mesh_data["elements"]),
              Rcpp::as<DMatrix<int>>(mesh_data["boundary"]));
        }
    };
    // getters
    const core::Mesh<M, N>& domain() const { return domain_; }
    const NodesType& nodes() const { return domain_.nodes(); }
    const DMatrix<int, Eigen::RowMajor>& elements() const { return domain_.elements(); }
    const NeighborsContainerType& neighbors() const { return domain_.neighbors(); }
    DMatrix<int> boundary() const { return domain_.boundary(); }
    // solves point location problem over triangulation
    DVector<int> locate(const DMatrix<double>& points) { return domain_.locate(points); }
    // destructor
    ~Mesh() = default;
};

}   // namespace r
}   // namespace fdapde

#endif   // __R_MESH_H__
