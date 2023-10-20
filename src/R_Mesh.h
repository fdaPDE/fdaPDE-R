#ifndef __R_MESH_H__
#define __R_MESH_H__

#include <fdaPDE/utils/symbols.h>
#include <fdaPDE/mesh.h>
#include<fdaPDE/finite_elements.h>
#include <memory>
using fdapde::core::ADT;
using fdapde::core::Mesh;

// M : local dimension, N : embedding dimension
template <int M, int N> class R_Mesh {
   private:
    typedef Mesh<M, N> DomainType;
    typedef ADT<M, N> PointLocatorType;
    
    // internal data (use shared_ptr since Rcpp::wrap will make a copy when passed as argument)
    DomainType domain_ {};
    PointLocatorType point_locator_;
    
    // FunctionBasis basis_ {};  
    
   public:
    // constructor
    R_Mesh(const Rcpp::List& mesh_data) :
        domain_(
          Rcpp::as<DMatrix<double>>(mesh_data["nodes"]), Rcpp::as<DMatrix<int>>(mesh_data["elements"]),
          Rcpp::as<DMatrix<int>>(mesh_data["boundary"])),  point_locator_(domain_) {};
    
    //static constexpr int fe_order_ = DomainType::order;
    
    // getters
    const DomainType& domain() const { return domain_; }
    const PointLocatorType& point_locator() const { return point_locator_; }
    const DMatrix<double>& nodes() const { return domain_.nodes(); }
    const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& elements() const {
        return domain_.elements();
    }
    const DMatrix<int>& neighbors() const { return domain_.neighbors(); }
    const DMatrix<int>& boundary() const { return domain_.boundary(); }

    // destructor
    ~R_Mesh() = default;
};

typedef R_Mesh<2, 2> Mesh_2D;
  
#endif   // __R_MESH_H__