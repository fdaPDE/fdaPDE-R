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

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <fdaPDE/linear_algebra.h>
  
// [[Rcpp::export]]
Eigen::SparseMatrix<double>
kronecker_tensor_product(const Eigen::SparseMatrix<double>& lhs, const Eigen::SparseMatrix<double>& rhs) {
    return Eigen::SparseMatrix<double>(fdapde::core::Kronecker(lhs, rhs));
}
