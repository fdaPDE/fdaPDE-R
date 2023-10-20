// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/utils/symbols.h>
#include "R_Mesh.h"
#include <fdaPDE/pde.h>
using fdapde::core::advection;
using fdapde::core::diffusion;
using fdapde::core::FEM;
using fdapde::core::FiniteElementBasis;
using fdapde::core::Integrator;
using fdapde::core::LagrangianElement;
using fdapde::core::laplacian;
using fdapde::core::dt;
using fdapde::core::Mesh;
using fdapde::core::PDE;
using fdapde::core::PDEBase;
using fdapde::core::reaction;
using fdapde::core::fem_order;

// M : local dimension, N : embedding dimension, R : FEM order
template <int M, int N, int R> class R_PDE {
   private:
    typedef Mesh<M, N> DomainType;
    typedef Integrator<DomainType::local_dimension, R> QuadratureRule;
    typedef LagrangianElement<DomainType::local_dimension, R> FunctionSpace;
    typedef FiniteElementBasis<FunctionSpace> FunctionBasis;

    // internal data
    DomainType domain_ {};
    QuadratureRule integrator_ {};
    std::unique_ptr<PDEBase> pde_;                      
    int pde_type_;
    Rcpp::Nullable<Rcpp::List> pde_parameters_;   // eventual PDE coefficients
    DMatrix<double> u_;                           // discretized forcing term
    FunctionBasis basis_ {};                      // functional basis
    Rcpp::Environment solution_;
    bool is_parabolic_ = false;

   public:
    // constructor
    R_PDE(
      Rcpp::Environment mesh, int pde_type, const Rcpp::Nullable<Rcpp::List>& pde_parameters,
      Rcpp::Environment solution) :
        pde_type_(pde_type), pde_parameters_(pde_parameters), solution_(solution) {
        SEXP meshptr = mesh[".pointer"];
        R_Mesh<M, N>* ptr = (R_Mesh<M, N>*)R_ExternalPtrAddr(meshptr);
        domain_ = ptr->domain();   // assign domain (re-enumerate if R>1)
    }

    // setters
    void set_dirichlet_bc(const DMatrix<double>& data) { pde_->set_dirichlet_bc(data); }
    void set_forcing(const DMatrix<double>& data) { u_ = data; }
    void set_initial_condition(const DVector<double>& data) {
        if(!is_parabolic_) 
            throw std::runtime_error("set initial condition for elliptic problem.");
        pde_->set_initial_condition(data);}

    // getters
    DMatrix<double> get_quadrature_nodes() const { return integrator_.quadrature_nodes(domain_); };
    DMatrix<double> get_dofs_coordinates() const { return pde_->dof_coords(); };
    SpMatrix<double> R0() const { return pde_->R0(); }
    SpMatrix<double> R1() const { return pde_->R1(); }
    DMatrix<double> solution() const { return pde_->solution(); }

    // init pde_ pointer
    void init() {
        switch (pde_type_) {
        case 1: {   // L = K*Laplacian + b*grad(f) + c*f 
            Rcpp::List pde_parameters(pde_parameters_);
            double K = Rcpp::as<double>(pde_parameters["diffusion"]);
            SVector<M> b = Rcpp::as<DMatrix<double>>(pde_parameters["transport"]);
            double c = Rcpp::as<double>(pde_parameters["reaction"]);
            // define bilinear form
            auto L = K*laplacian<FEM>() + advection<FEM>(b) + reaction<FEM>(c);
            pde_ = std::make_unique<PDE<DomainType, decltype(L), DMatrix<double>, FEM, fem_order<R>>>(domain_, L, u_); 
        } break;
        case 2: {   // L = div(K*grad(f)) + b*grad(f) + c*f
            if (pde_parameters_.isNotNull()) {
                // extract parameters from pack
                Rcpp::List pde_parameters(pde_parameters_);
                SMatrix<M> K = Rcpp::as<DMatrix<double>>(pde_parameters["diffusion"]);
                SVector<M> b = Rcpp::as<DMatrix<double>>(pde_parameters["transport"]);
                double c = Rcpp::as<double>(pde_parameters["reaction"]);
                // define bilinear form
                auto L = diffusion<FEM>(K) + advection<FEM>(b) + reaction<FEM>(c);
                pde_ = std::make_unique<PDE<DomainType, decltype(L), DMatrix<double>, FEM, fem_order<R>>>(domain_, L, u_);
            } 
            
            else {
                throw std::runtime_error("pde parameters not supplied.");
            }
        } break;
        case 3: {
            Rcpp::List pde_parameters(pde_parameters_);
            double K = Rcpp::as<double>(pde_parameters["diffusion"]);
            SVector<M> b = Rcpp::as<DMatrix<double>>(pde_parameters["transport"]);
            double c = Rcpp::as<double>(pde_parameters["reaction"]);
            this-> is_parabolic_ = true;
            // define bilinear form
            auto L = dt<FEM>() + K*laplacian<FEM>() + advection<FEM>(b) + reaction<FEM>(c);
            pde_ = std::make_unique<PDE<DomainType, decltype(L), DMatrix<double>, FEM, fem_order<R>>>(domain_, L, u_); 
        } break;
        case 4: {
            if (pde_parameters_.isNotNull()) {
                // extract parameters from pack
                Rcpp::List pde_parameters(pde_parameters_);
                SMatrix<M> K = Rcpp::as<DMatrix<double>>(pde_parameters["diffusion"]);
                SVector<M> b = Rcpp::as<DMatrix<double>>(pde_parameters["transport"]);
                double c = Rcpp::as<double>(pde_parameters["reaction"]);
                this-> is_parabolic_ = true;
                // define bilinear form
                auto L = dt<FEM>() + diffusion<FEM>(K) + advection<FEM>(b) + reaction<FEM>(c);
                pde_ = std::make_unique<PDE<DomainType, decltype(L), DMatrix<double>, FEM, fem_order<R>>>(domain_, L, u_);
            } 
            
            else {
                throw std::runtime_error("pde parameters not supplied.");
            }
        } break;
        }
        // initialize pointer
        pde_->init();
        return;
    }

    void solve() {
        pde_->solve();
        solution_["coeff"] = pde_->solution();
    }

    // eval solution at point
    DMatrix<double> eval(Rcpp::Environment mesh, const DMatrix<double>& coeff, const DMatrix<double>& point) {
        SEXP meshptr = mesh[".pointer"];
        R_Mesh<M, N>* ptr = (R_Mesh<M, N>*)R_ExternalPtrAddr(meshptr);
        if (point.cols() != M) {
            throw std::runtime_error("evaluation point must be a " + std::to_string(M) + "-dimensional vector");
        }
        DMatrix<double> evals;
        evals.resize(point.rows(), 1);
        for (std::size_t j = 0; j < point.rows(); ++j) {
            SVector<N> p;
            for (int i = 0; i < N; ++i) p[i] = point(j, i);
            // locate element containing point
            auto e = ptr->point_locator().locate(p);
            // compute value of function at point
            double v = std::numeric_limits<double>::quiet_NaN();
            if (e != nullptr) {
                v = 0;
                for (std::size_t k = 0; k < FunctionSpace::n_basis; ++k) {
                    v += coeff(ptr->domain().elements()(e->ID(), k), 0) * basis_(*e, k)(p);
                }
            }
            evals(j, 0) = v;
        }
        return evals;
    }

    // destructor
    ~R_PDE() = default;
};

// Rcpp modules definition

typedef R_PDE<2, 2, 1> PDE_2D_ORDER_1;
RCPP_MODULE(PDE_2D_ORDER_1) {
    Rcpp::class_<PDE_2D_ORDER_1>("PDE_2D_ORDER_1")
      .constructor<Rcpp::Environment, int, Rcpp::Nullable<Rcpp::List>, Rcpp::Environment>()
      // getters
      .method("get_quadrature_nodes", &PDE_2D_ORDER_1::get_quadrature_nodes)
      .method("get_dofs_coordinates", &PDE_2D_ORDER_1::get_dofs_coordinates)
      .method("get_mass", &PDE_2D_ORDER_1::R0)
      .method("get_stiff", &PDE_2D_ORDER_1::R1)
      .method("solution", &PDE_2D_ORDER_1::solution)
      .method("eval", &PDE_2D_ORDER_1::eval)
      // setters
      .method("set_dirichlet_bc", &PDE_2D_ORDER_1::set_dirichlet_bc)
      .method("set_forcing", &PDE_2D_ORDER_1::set_forcing)
      // init and solve
      .method("init", &PDE_2D_ORDER_1::init)
      .method("solve", &PDE_2D_ORDER_1::solve)
      .method("set_initial_condition", &PDE_2D_ORDER_1::set_initial_condition);
}

typedef R_PDE<2, 2, 2> PDE_2D_ORDER_2;
RCPP_MODULE(PDE_2D_ORDER_2) {
    Rcpp::class_<PDE_2D_ORDER_2>("PDE_2D_ORDER_2")
      .constructor<Rcpp::Environment, int, Rcpp::Nullable<Rcpp::List>, Rcpp::Environment>()
      // getters
      .method("get_quadrature_nodes", &PDE_2D_ORDER_2::get_quadrature_nodes)
      .method("get_dofs_coordinates", &PDE_2D_ORDER_2::get_dofs_coordinates)
      .method("get_mass", &PDE_2D_ORDER_2::R0)
      .method("get_stiff", &PDE_2D_ORDER_2::R1)
      .method("solution", &PDE_2D_ORDER_2::solution)
      .method("eval", &PDE_2D_ORDER_2::eval)
      // setters
      .method("set_dirichlet_bc", &PDE_2D_ORDER_2::set_dirichlet_bc)
      .method("set_forcing", &PDE_2D_ORDER_2::set_forcing)
      // init and solve
      .method("init", &PDE_2D_ORDER_2::init)
      .method("solve", &PDE_2D_ORDER_2::solve)
      .method("set_initial_condition", &PDE_2D_ORDER_2::set_initial_condition);
}
