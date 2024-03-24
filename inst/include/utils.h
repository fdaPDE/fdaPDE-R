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

#ifndef __R_UTILS_H__
#define __R_UTILS_H__

namespace fdapde {
namespace r {

template <typename T> T* get_env_as(const Rcpp::Environment& r_env) {
    SEXP pde_ptr = r_env[".pointer"];
    return reinterpret_cast<T*>(R_ExternalPtrAddr(pde_ptr));
}

}   // namespace r
}   // namespace fdapde

#endif   // __R_UTILS_H__
