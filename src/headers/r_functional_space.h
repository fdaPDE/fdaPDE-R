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

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/utils.h>
#include <fdaPDE/finite_elements.h>
#include <fdaPDE/mesh.h>
#include "r_mesh.h"
using fdapde::core::Mesh;
using fdapde::core::Integrator;
using fdapde::core::LagrangianBasis;
using fdapde::core::pointwise_evaluation;
using fdapde::core::areal_evaluation;
using fdapde::core::eval;

// type-erasure wrapper for a function space object
struct I_FunctionalSpace {
    template <typename T>
    using fn_ptrs = fdapde::mem_fn_ptrs<
      &T::size, &T::operator(), &T::template eval<pointwise_evaluation>, &T::template eval<areal_evaluation>>;
    // forwardings
    decltype(auto) size() const { return fdapde::invoke<std::size_t, 0>(*this); }
    decltype(auto) operator()(const DVector<double>& c, const DMatrix<double>& locs) const {
      return fdapde::invoke<DVector<double>, 1>(*this, c, locs);
    }
    decltype(auto) eval(int evaluation_type, const DMatrix<double>& locs) const {
        using RetType = std::pair<SpMatrix<double>, DVector<double>>;
        switch (evaluation_type) {   // run-time switch based on evaluation strategy
        case eval::pointwise: // 0
            return fdapde::invoke<RetType, 2>(*this, locs).first;
        case eval::areal: // 1
            return fdapde::invoke<RetType, 3>(*this, locs).first;
        }
        return SpMatrix<double>{};
    }
};
using FunctionalSpace = fdapde::erase<fdapde::heap_storage, I_FunctionalSpace>;

// supported basis function types
enum space_type { fem_lagrange };

template <int M, int N, int R> class R_FunctionalSpace {
   private:
    using DomainType = Mesh<M, N>;
    using FunctionalSpaceType = FunctionalSpace;
    using QuadratureRule = Integrator<DomainType::local_dimension, R>; // exact for FEM spaces (TODO: generalize)

    // internal data
    DomainType domain_ {};
    FunctionalSpaceType fun_space_ {}; // wrapped functional space
    QuadratureRule integrator_ {};     // quadrature rule (exact for the provided fem order)
    int space_type_;
   public:
    // constructor
    R_FunctionalSpace(Rcpp::Environment mesh, int space_type) : space_type_(space_type) {
        // set domain
        SEXP meshptr = mesh[".pointer"];
        R_Mesh<M, N>* ptr = reinterpret_cast<R_Mesh<M, N>*>(R_ExternalPtrAddr(meshptr));
	domain_ = ptr->domain();
	// define functional space
	std::size_t size = (R == 1 ? domain_.n_nodes() : domain_.n_nodes() + domain_.n_edges());
	switch(space_type_) {
	case space_type::fem_lagrange: {
	  fun_space_ = LagrangianBasis<DomainType, R>(domain_, size);
	} break;
        }
    }
    std::size_t size() const { return fun_space_.size(); }   // number of basis functions on phyiscal domain
    // given a coefficient vector c \in \mathbb{R}^size_, evaluates the corresponding basis expansion at locs
    DVector<double> operator()(const DVector<double>& c, const DMatrix<double>& locs) const {
        return fun_space_(c, locs);
    }
    // returns the \Psi is the matrix of basis functions evaluations according
    // to the given policy (see the specific policy for details)
    SpMatrix<double> eval(int evaluation_type, const DMatrix<double>& locs) const {
        return fun_space_.eval(evaluation_type, locs);
    }
    // given a coefficient vector c \in \mathbb{R}^size_, evaluates the integral of the corresponding basis expasion
    // over the whole domain
    double integrate(const DVector<double>& c) const {
      return integrator_.integrate(domain_, fun_space_(c, integrator_.quadrature_nodes(domain_)));
    }
    // destructor
    ~R_FunctionalSpace() = default;
};

#endif   // __R_FUNCTIONAL_SPACE_H__
