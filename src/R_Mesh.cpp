#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "R_Mesh.h"

// Rcpp module definition
RCPP_MODULE(Mesh_2D) {
    Rcpp::class_<Mesh_2D>("Mesh_2D")
      .constructor<Rcpp::List>()
      .method("nodes", &Mesh_2D::nodes)
      .method("elements", &Mesh_2D::elements)
      .method("neighbors", &Mesh_2D::neighbors)
      .method("boundary", &Mesh_2D::boundary);
}
