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

#ifndef __R_PDE_H__
#define __R_PDE_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/pde.h>
#include <fdaPDE/utils/symbols.h>

#include "r_mesh.h"
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

// supported pde instantiations
enum pde_type {
    simple_laplacian,        // \mu * \Delta f
    second_order_elliptic,   // div(K * \nabla f) + dot(b, \nabla f) + c * f
    second_order_parabolic   // df/dt + div(K * \nabla f) + dot(b, \nabla f) + c * f
};

// base non-templated class holding a type-erasure wrapper for a general pde object
struct PDEWrapper {
   protected:
    pde_ptr pde_;   // type-erased pde wrapper
   public:
    pde_ptr& get_pde() { return pde_; }
    const pde_ptr& get_pde() const { return pde_; }
};

template <int M, int N, int R> class R_PDE : public PDEWrapper {
   private:
    using DomainType = Mesh<M, N>;
    using QuadratureRule = Integrator<FEM, DomainType::local_dimension, R>;
    template <typename L> using PDEType = PDE<DomainType, L, DMatrix<double>, FEM, fem_order<R>>;

    // internal data
    DomainType domain_ {};           // triangulation
    QuadratureRule integrator_ {};   // quadrature rule (exact for the provided fem order)
    int pde_type_;                   // one of the pde_type enum values
   public:
    // constructor
    R_PDE(Rcpp::Environment mesh, int pde_type, const Rcpp::Nullable<Rcpp::List>& pde_parameters) :
        pde_type_(pde_type) {
        // set domain
        SEXP meshptr = mesh[".pointer"];
        R_Mesh<M, N>* ptr = reinterpret_cast<R_Mesh<M, N>*>(R_ExternalPtrAddr(meshptr));
        domain_ = ptr->domain();

        // unpack pde parameters
        Rcpp::List pde_parameters_(pde_parameters);
        SVector<M> b = Rcpp::as<DMatrix<double>>(pde_parameters_["transport"]);
        double c = Rcpp::as<double>(pde_parameters_["reaction"]);
        // instantiate pde template based on pde type
        switch (pde_type_) {
        case pde_type::simple_laplacian: {   // special acceleration for simple laplacian penalty
            double mu = Rcpp::as<double>(pde_parameters_["diffusion"]);
            auto L = mu * laplacian<FEM>();
            pde_ = PDEType<decltype(L)>(domain_, L);
        } break;
        case pde_type::second_order_elliptic: {
            SMatrix<M> K = Rcpp::as<DMatrix<double>>(pde_parameters_["diffusion"]);
            auto L = diffusion<FEM>(K) + advection<FEM>(b) + reaction<FEM>(c);
            pde_ = PDEType<decltype(L)>(domain_, L);
        } break;
        case pde_type::second_order_parabolic: {
            SMatrix<M> K = Rcpp::as<DMatrix<double>>(pde_parameters_["diffusion"]);
            auto L = dt<FEM>() + diffusion<FEM>(K) + advection<FEM>(b) + reaction<FEM>(c);
            pde_ = PDEType<decltype(L)>(domain_, DVector<double>(), L);
        } break;
        }
    }
    // setters
    void set_dirichlet_bc(const DMatrix<double>& data) { pde_.set_dirichlet_bc(data); }
    void set_forcing(const DMatrix<double>& data) { pde_.set_forcing(data); }
    void set_initial_condition(const DVector<double>& data) {
        if (pde_type_ != pde_type::second_order_parabolic) throw "set_initial_condition is for space-time pdes only.";
        pde_.set_initial_condition(data);
    }
    // getters
    DMatrix<double> get_quadrature_nodes() const { return integrator_.quadrature_nodes(domain_); };
    DMatrix<double> get_dofs_coordinates() const { return pde_.dof_coords(); };
    // avaiable only after initialization
    const SpMatrix<double>& R0() const { return pde_.mass(); }
    const SpMatrix<double>& R1() const { return pde_.stiff(); }
    const DMatrix<double>& u() const { return pde_.force(); }
    // initialize internal pde status
    void init() { pde_.init(); }
    // destructor
    ~R_PDE() = default;
};

#endif   // __R_PDE_H__
