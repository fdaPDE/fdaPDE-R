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

#include <fdaPDE/finite_elements.h>
#include <fdaPDE/geometry.h>
#include <fdaPDE/utils/symbols.h>
using fdapde::core::TreeSearch;
using fdapde::core::Mesh;
using fdapde::core::neighboring_structure;

// Rcpp wrapper for mesh. M : local dimension, N : embedding dimension
template <int M, int N> class R_Mesh {
   private:
    using DomainType = Mesh<M, N>;
    DomainType domain_ {};
   public:
    // constructor
    R_Mesh(const Rcpp::List& mesh_data) :
        domain_(
          Rcpp::as<DMatrix<double>>(mesh_data["nodes"]), Rcpp::as<DMatrix<int>>(mesh_data["elements"]),
          Rcpp::as<DMatrix<int>>(mesh_data["boundary"])) {};
    // getters
    const DomainType& domain() const { return domain_; }
    const DMatrix<double>& nodes() const { return domain_.nodes(); }
    const DMatrix<int, Eigen::RowMajor>& elements() const { return domain_.elements(); }
    const typename neighboring_structure<M, N>::type& neighbors() const { return domain_.neighbors(); }
    const DMatrix<int>& boundary() const { return domain_.boundary(); }
    // destructor
    ~R_Mesh() = default;
};

#endif   // __R_MESH_H__
