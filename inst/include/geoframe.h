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

#ifndef __R_GEOFRAME_H__
#define __R_GEOFRAME_H__

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/geoframe.h>
#include "mesh.h"
#include "utility.h"

namespace fdapde {
namespace r {

template <typename Triangulation> class GeoFrame {
    using int_t = int;
    using str_t = std::string;
    using dbl_t = double;
   public:
    static constexpr int local_dim = Triangulation::local_dim;
    static constexpr int embed_dim = Triangulation::embed_dim;
    static constexpr int Order = 1;
    using TriangulationType = fdapde::Triangulation<local_dim, embed_dim>;
    using geoframe_t = fdapde::GeoFrame<TriangulationType>;

    GeoFrame(const Rcpp::Environment& triangulation) : data_(get_env_as<TriangulationType>(triangulation)) { }
    // constrcut from an existing geoframe using layer subsetting
    GeoFrame(
      const Rcpp::Environment& geoframe, const std::string& layer_name, const std::vector<int>& rows,
      const std::vector<std::string>& cols) {
        geoframe_t& gf = get_env_as<geoframe_t>(geoframe);
        data_ = geoframe_t(gf.template triangulation<0>());
        auto make_ = [&]<typename GeoInfo>(GeoInfo) {
            auto row_filter = geo_cast<GeoInfo>(gf[layer_name]).select(rows.begin(), rows.end());
            data_.template insert_scalar_layer<GeoInfo>(layer_name, row_filter, cols);
        };
        fdapde::ltype ltype = gf[layer_name].category()[0];
        if (ltype == ltype::point) { make_(POINT {}); }
        if (ltype == ltype::areal) { make_(POLYGON {}); }
    }
    // layer insertion
    template <typename GeoDescriptor>
    void insert_scalar_point_layer(const std::string& layer_name, const GeoDescriptor& locs, const Rcpp::List& data) {
        auto& l = data_.template insert_scalar_layer<POINT>(layer_name, locs);
        // copy data from R list
        auto copy_ = [&]<typename T>(const std::string& field) {
            if (data.containsElementNamed(field.data())) {
                Rcpp::List lst = Rcpp::as<Rcpp::List>(data[field]);
                if (lst.size() != 0) {
                    std::vector<std::string> names = lst.names();
                    for (const std::string& name : names) { l.load_vec(name, Rcpp::as<std::vector<T>>(lst[name])); }
                }
            }
        };
        copy_.template operator()<int_t>("int_data");
        copy_.template operator()<dbl_t>("dbl_data");
        copy_.template operator()<str_t>("str_data");
    }
    void
    insert_scalar_areal_layer(const std::string& layer_name, const std::vector<int>& regions, const Rcpp::List& data) {
        auto& l = data_.template insert_scalar_layer<POLYGON>(layer_name, regions);
        // copy data from R list
        auto copy_ = [&]<typename T>(const std::string& field) {
            if (data.containsElementNamed(field.data())) {
                Rcpp::List lst = Rcpp::as<Rcpp::List>(data[field]);
                if (lst.size() != 0) {
                    std::vector<std::string> names = lst.names();
                    for (const std::string& name : names) { l.load_vec(name, Rcpp::as<std::vector<T>>(lst[name])); }
                }
            }
        };
        copy_.template operator()<int_t>("int_data");
        copy_.template operator()<dbl_t>("dbl_data");
        copy_.template operator()<str_t>("str_data");
    }
    void load_shp(const std::string& layer_name, const std::string& filename) { data_.load_shp(layer_name, filename); }
    template <typename T>
    void insert(const std::string& layer_name, const std::string& colname, const std::vector<T>& data) {
        data_[layer_name].add_column(colname, data);
    }
    template <typename T>
    void blk_insert(
      const std::string& layer_name, const std::string& colname, const Eigen::Matrix<T, Dynamic, Dynamic>& data) {
        data_[layer_name].add_block(colname, data);
    }
    // observers
    Eigen::Matrix<double, Dynamic, Dynamic> bbox() { return data_.template triangulation<0>().bbox(); }
    int n_nodes() { return data_.template triangulation<0>().n_nodes(); }
    int n_cells() { return data_.template triangulation<0>().n_cells(); }

    template <typename T>
    void assign(
      const std::string& layer_name, const std::vector<int>& rows, const std::string& column,
      const std::vector<T>& values) {
        fdapde::ltype ltype = data_[layer_name].category()[0];
        auto assign_ = [&]<typename GeoInfo>(GeoInfo, geoframe_t& gf) {
            auto row_filter =
              geo_cast<GeoInfo>(gf[layer_name]).select(rows.begin(), rows.end()).template col<T>(column);
            for (int i = 0, n = rows.size(); i < n; ++i) { row_filter(i, 0) = values[i]; }
        };
        if (ltype == ltype::point) { assign_(POINT   {}, data_); }
	if (ltype == ltype::areal) { assign_(POLYGON {}, data_); }
    }
    template <typename T>
    std::vector<T> access(const std::string& layer_name, const std::vector<int>& rows, const std::string& column) {
        fdapde::ltype ltype = data_[layer_name].category()[0];
        auto access_ = [&]<typename GeoInfo>(GeoInfo, std::vector<T>& buff) {
            auto row_filter = geo_cast<GeoInfo>(data_[layer_name]).select(rows.begin(), rows.end());
            buff.reserve(rows.size());
            for (int i = 0, n = rows.size(); i < n; ++i) {
                buff.emplace_back(row_filter.template col<T>(column)(i, 0));
            }
            return buff;
        };
        std::vector<T> buff;
        if (ltype == ltype::point) { buff = access_(POINT   {}, buff); }
	if (ltype == ltype::areal) { buff = access_(POLYGON {}, buff); }
	return buff;
    }

    int ltype(const std::string& layer_name) const { return static_cast<int>(data_[layer_name].category()[0]); }
    int ctype(const std::string& layer_name, const std::string& col_name) {
        fdapde::ltype ltype = data_[layer_name].category()[0];
        int ctype_ = 0;
        if (ltype == ltype::areal) {
            ctype_ = int(geo_cast<POLYGON>(data_[layer_name]).data().field_descriptor(col_name).type_id());
        }
        if (ltype == ltype::point) {
            ctype_ = int(geo_cast<POINT  >(data_[layer_name]).data().field_descriptor(col_name).type_id());
        }
        return ctype_;
    }
    std::vector<std::string> colnames_all() const { return data_.colnames(); }
    std::vector<std::string> colnames(const std::string& layer_name) {
        fdapde::ltype ltype = data_[layer_name].category()[0];
        std::vector<std::string> cols_;
        if (ltype == ltype::areal) { cols_ = geo_cast<POLYGON>(data_[layer_name]).data().colnames(); }
        if (ltype == ltype::point) { cols_ = geo_cast<POINT  >(data_[layer_name]).data().colnames(); }
	return cols_;
    }
    int rows(const std::string& layer_name) {
        fdapde::ltype ltype = data_[layer_name].category()[0];
        int rows = 0;
        if (ltype == ltype::areal) { rows = geo_cast<POLYGON>(data_[layer_name]).rows(); }
	if (ltype == ltype::point) { rows = geo_cast<POINT  >(data_[layer_name]).rows(); }
        return rows;
    }
    int cols(const std::string& layer_name) {
        fdapde::ltype ltype = data_[layer_name].category()[0];
        int cols = 0;
        if (ltype == ltype::areal) { cols = geo_cast<POLYGON>(data_[layer_name]).cols(); }
	if (ltype == ltype::point) { cols = geo_cast<POINT  >(data_[layer_name]).cols(); }
        return cols;
    }  
    // areal layer
    // Eigen::Matrix<double, Dynamic, Dynamic>
    // areal_sample(const std::string& layer_name, int n_samples, int seed = fdapde::random_seed) {
    //     const internals::areal_layer<GeoFrame_>& layer = data_.get_as(layer_t::areal, layer_name);
    //     return layer.sample(n_samples, seed);
    // }
    // std::vector<Eigen::Matrix<double, Dynamic, Dynamic>> areal_poly_nodes(const std::string& layer_name) {
    //     const internals::areal_layer<GeoFrame_>& layer = data_.get_as(layer_t::areal, layer_name);
    //     int n_regions = layer.n_regions();
    //     std::vector<Eigen::Matrix<double, Dynamic, Dynamic>> nodes;
    //     nodes.resize(n_regions);
    //     for (int i = 0; i < n_regions; ++i) { nodes[i] = layer.geometry(i).nodes(); }
    //     return nodes;
    // }
    // std::vector<Eigen::Matrix<int, Dynamic, Dynamic>> areal_poly_edges(const std::string& layer_name) {
    //     const internals::areal_layer<GeoFrame_>& layer = data_.get_as(layer_t::areal, layer_name);
    //     int n_regions = layer.n_regions();
    //     std::vector<Eigen::Matrix<int, Dynamic, Dynamic>> edges;
    //     edges.resize(n_regions);
    //     for (int i = 0; i < n_regions; ++i) { edges[i] = layer.geometry(i).boundary_edges().array() + 1; }
    //     return edges;
    // }

    // point layer
    Eigen::Matrix<double, Dynamic, Dynamic> point_coordinates(const std::string& layer_name) {
        // assert layer is pointwise
        return geo_index_cast<0, POINT>(data_[layer_name]).coordinates();
    }

    // Eigen::Matrix<double, Dynamic, Dynamic> sample(const std::string& layer_name, int n_sample) {
    //     const internals::areal_layer<GeoFrame_>& layer = data_.get_as(layer_t::areal, layer_name);
    //     return layer.sample(n_sample);
    // }

    geoframe_t& data() { return data_; }
   private:
    geoframe_t data_;
};



  
}   // namespace r
}   // namespace fdapde

#endif // __R_GEOFRAME_H__
