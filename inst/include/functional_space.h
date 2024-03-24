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

#ifndef __R_FUNCTIONAL_SPACE_H__
#define __R_FUNCTIONAL_SPACE_H__

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/finite_elements.h>
#include <fdaPDE/geometry.h>
#include <fdaPDE/utils.h>
#include <fdaPDE/splines.h>

#include "mesh.h"

namespace fdapde {
namespace r {

// functional basis concept
struct FunctionalSpace__ {
    template <typename T>
    using fn_ptrs = fdapde::mem_fn_ptrs<
      &T::size, &T::operator(), &T::template eval<core::pointwise_evaluation>,
      &T::template eval<core::areal_evaluation>, &T::dofs_coords>;
    // forwardings
    std::size_t size() const { return fdapde::invoke<std::size_t, 0>(*this); }
    DVector<double> operator()(const DVector<double>& c, const DMatrix<double>& locs) const {
        return fdapde::invoke<DVector<double>, 1>(*this, c, locs);
    }
    SpMatrix<double> eval(int evaluation_type, const DMatrix<double>& locs) const {
        using RetType = std::pair<SpMatrix<double>, DVector<double>>;
        switch (evaluation_type) {
        case 0:   // pointwise evaluation
            return fdapde::invoke<RetType, 2>(*this, locs).first;
        case 1:   // areal evaluation
            return fdapde::invoke<RetType, 3>(*this, locs).first;
        }
        return SpMatrix<double> {};
    }
    DMatrix<double> dofs_coords() { return fdapde::invoke<DMatrix<double>, 4>(*this); }
};

// TODO: we could merge the two templates above in one, required: finite elements should work also on 1D domains,

// space of finite element basis functions
template <typename DomainType, int Order> class FEFunctionalSpace {
   private:
    using FunctionalSpaceType = fdapde::erase<fdapde::heap_storage, FunctionalSpace__>;
    using QuadratureRule = core::Integrator<core::FEM, DomainType::local_dimension, Order>;
    using VectorType = DVector<double>;
    using MatrixType = DMatrix<double>;

    FunctionalSpaceType fun_space_;
    DomainType domain_;   // pointer to physical domain
    QuadratureRule integrator_ {};   // quadrature rule exact for the requested Order_
    MatrixType quadrature_nodes_;    // quadrature nodes on phyisical domain
    MatrixType dofs_coords_;         // dofs phyisical coordinates
   public:
    // constructor
    FEFunctionalSpace(Rcpp::Environment mesh, int type) {
        // set domain
        using RDomainType = r::Mesh<DomainType::local_dimension, DomainType::embedding_dimension>;
        SEXP meshptr = mesh[".pointer"];
        RDomainType* ptr = reinterpret_cast<RDomainType*>(R_ExternalPtrAddr(meshptr));
        domain_ = ptr->domain();
        // define functional space
        switch (type) {
        case 0: {   // continous Galerkin finite elements (CG, Lagrange)
            fun_space_ = core::LagrangianBasis<DomainType, Order>(domain_);
        } break;
        }
    }
    std::size_t size() const { return fun_space_.size(); }
    const MatrixType& dofs_coords() {
        if (is_empty(dofs_coords_)) dofs_coords_ = fun_space_.dofs_coords();
        return dofs_coords_;
    }
    const MatrixType& quadrature_nodes() {
        if (is_empty(quadrature_nodes_)) quadrature_nodes_ = integrator_.quadrature_nodes(domain_);
        return quadrature_nodes_;
    }
    // pointwise evaluates the basis expansion of c \in \mathbb{R}^size_ at locs
    VectorType eval_expansion(const VectorType& c, const MatrixType& locs) const { return fun_space_(c, locs); }
    SpMatrix<double> Psi(int eval_type, const MatrixType& locs) const {
        fdapde_assert(eval_type == 0 || eval_type == 1);
        return fun_space_.eval(eval_type, locs);
    }
    // evaluates the integral of the basis expasion of c \in \mathbb{R}^size_ over the whole domain
    double integrate(const VectorType& c) {
        if (is_empty(quadrature_nodes_)) quadrature_nodes_ = integrator_.quadrature_nodes(domain_);
        return integrator_.integrate(domain_, fun_space_(c, quadrature_nodes_));
    }
    // destructor
    ~FEFunctionalSpace() = default;
};

// space of BSpline basis functions
template <int Order> class BSplineFunctionalSpace {
   private:
    using DomainType = core::Mesh<1, 1>;
    using FunctionalSpaceType = fdapde::erase<fdapde::heap_storage, FunctionalSpace__>;
    using QuadratureRule = core::Integrator<core::SPLINE, 1, Order>;
    using VectorType = DVector<double>;
    using MatrixType = DMatrix<double>;

    FunctionalSpaceType fun_space_;
    DomainType domain_;   // pointer to physical domain
    QuadratureRule integrator_ {};   // quadrature rule exact for the requested Order_
    MatrixType quadrature_nodes_;    // quadrature nodes on phyisical domain
    MatrixType dofs_coords_;         // dofs phyisical coordinates
   public:
    // constructor
    BSplineFunctionalSpace(Rcpp::Environment mesh, int type) {
        // set domain
        using RDomainType = r::Mesh<DomainType::local_dimension, DomainType::embedding_dimension>;
        SEXP meshptr = mesh[".pointer"];
        RDomainType* ptr = reinterpret_cast<RDomainType*>(R_ExternalPtrAddr(meshptr));
        domain_ = ptr->domain();
        // define functional space
        switch (type) {
        case 0: { // cubic BSpline basis function
             fun_space_ = core::SplineBasis<Order>(domain_.nodes());
        } break;
        }
    }
    std::size_t size() const { return fun_space_.size(); }
    const MatrixType& dofs_coords() {
        if (is_empty(dofs_coords_)) dofs_coords_ = fun_space_.dofs_coords();
        return dofs_coords_;
    }
    const MatrixType& quadrature_nodes() {
        if (is_empty(quadrature_nodes_)) quadrature_nodes_ = integrator_.quadrature_nodes(domain_);
        return quadrature_nodes_;
    }
    // pointwise evaluates the basis expansion of c \in \mathbb{R}^size_ at locs
    VectorType eval_expansion(const VectorType& c, const MatrixType& locs) const { return fun_space_(c, locs); }
    SpMatrix<double> Psi(int eval_type, const MatrixType& locs) const {
        fdapde_assert(eval_type == 0 || eval_type == 1);
        return fun_space_.eval(eval_type, locs);
    }
    // destructor
    ~BSplineFunctionalSpace() = default;
};
  

}   // namespace r
}   // namespace fdapde

#endif   // __R_FUNCTIONAL_SPACE_H__
