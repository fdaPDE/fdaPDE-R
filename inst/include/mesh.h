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

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/geometry.h>

namespace fdapde {
namespace r {

template <int LocalDim, int EmbedDim> class Triangulation;
template <int LocalDim, int EmbedDim> class TriangulationBase {
   public:
    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;
    using double_mtx = Eigen::Matrix<double, Dynamic, Dynamic>;
    using int_mtx = Eigen::Matrix<int, Dynamic, Dynamic>;
    using triangulation_t = fdapde::Triangulation<local_dim, embed_dim>;
    static constexpr int n_nodes_per_cell = triangulation_t::n_nodes_per_cell;
    using node_t = Eigen::Matrix<double, embed_dim, 1>;

    TriangulationBase() noexcept = default;
    TriangulationBase(const double_mtx& nodes, const int_mtx& cells, const int_mtx& boundary) noexcept :
        triangulation_(nodes, cells, boundary, cache_cells) { }   // always enable cell caching
    // observers
    const double_mtx& nodes() const { return triangulation_.nodes(); }
    const Eigen::Matrix<int, Dynamic, Dynamic, Eigen::RowMajor>& cells() const { return triangulation_.cells(); }
    int_mtx boundary_nodes() const {
        int_mtx boundary_(triangulation_.n_nodes(), 1);
        for (int i = 0, n = triangulation_.n_nodes(); i < n; ++i) {
            boundary_(i, 0) = triangulation_.is_node_on_boundary(i) ? 1 : 0;
        }
        return boundary_;
    }
    int n_cells() const { return triangulation_.n_cells(); }
    int n_nodes() const { return triangulation_.n_nodes(); }
    int n_boundary_nodes() const { return triangulation_.n_boundary_nodes(); }
    double_mtx bbox() const { return triangulation_.bbox(); }
    double measure() const { return triangulation_.measure(); }
    double marked_measure(int marker) const {
        double measure_ = 0;
        for (auto it = triangulation_.cells_begin(marker); it != triangulation_.cells_end(marker); ++it) {
            measure_ += it->measure();
        }
        return measure_;
    }
    // random sample in triangulation
    double_mtx sample(int n_samples, int seed = random_seed) const { return triangulation_.sample(n_samples, seed); }
    // cell marker handling
    void clear_cells_markers() { triangulation_.clear_cell_markers(); }
    const std::vector<int>& cells_markers() const { return triangulation_.cells_markers(); }
    void mark_cells(int marker, const std::vector<bool>& mask) {
        triangulation_.mark_cells(marker, [&](const auto& cell) { return mask[cell.id()]; });
    }
    // all indices of cells having marker = m
    std::vector<int> filter_cells_by_marker(int m) const {
        std::vector<int> cell_ids;
        if (!triangulation_.cells_markers().empty()) {
            for (int i = 0, n = triangulation_.n_cells(); i < n; ++i) {
                if (triangulation_.cells_markers()[i] == m) { cell_ids.push_back(i); }
            }
        }
        return cell_ids;
    }
    // cell properties
    const auto& cell_coords(int cell_id) { return triangulation_.cell(cell_id).nodes(); }
    double cell_measure(int cell_id) const { return triangulation_.cell(cell_id).measure(); }
    double_mtx cell_bbox(int cell_id) {
        auto [ll, ur] = triangulation_.cell(cell_id).bounding_box();
	double_mtx bbox(2, embed_dim);
	bbox.row(0) = ll;
	bbox.row(1) = ur;
        return bbox;
    }
    node_t cell_barycenter(int cell_id)   const { return triangulation_.cell(cell_id).barycenter(); }
    node_t cell_circumcenter(int cell_id) const { return triangulation_.cell(cell_id).circumcenter(); }
    double cell_diameter(int cell_id) const { return triangulation_.cell(cell_id).diameter(); }

    ~TriangulationBase() = default;
   protected:
    triangulation_t triangulation_;
};
  
// planar and surface triangulations
template <int EmbedDim> class Triangulation<2, EmbedDim> : public TriangulationBase<2, EmbedDim> {
   public:
    using Base = TriangulationBase<2, EmbedDim>;
    static constexpr int local_dim = Base::local_dim;
    static constexpr int embed_dim = Base::embed_dim;
    using Base::triangulation_;
    using double_mtx = typename Base::double_mtx;
    using int_mtx = typename Base::int_mtx;
    using triangulation_t = fdapde::Triangulation<local_dim, embed_dim>;
    static constexpr int n_nodes_per_edge = triangulation_t::n_nodes_per_edge;

    Triangulation() noexcept = default;
    Triangulation(const Rcpp::List& data) noexcept :
        Base(
          Rcpp::as<double_mtx>(data["nodes"]), Rcpp::as<int_mtx>(data["cells"]), Rcpp::as<int_mtx>(data["boundary"])) {
    }
    // observers
    const Eigen::Matrix<int, Dynamic, Dynamic, Eigen::RowMajor>& neighbors() const {
        return triangulation_.neighbors(); }
    int_mtx edges() const { return triangulation_.edges(); }
    int n_edges() const { return triangulation_.n_edges(); }
    int n_boundary_edges() const { return triangulation_.n_boundary_edges(); }
    int_mtx boundary_edges() const {
        int_mtx boundary_(triangulation_.n_edges(), 1);
        for (int i = 0, n = triangulation_.n_edges(); i < n; ++i) {
            boundary_(i, 0) = triangulation_.is_edge_on_boundary(i) ? 1 : 0;
        }
        return boundary_;
    }
    // edge marker handling
    const std::vector<int>& edges_markers() const { return triangulation_.edges_markers(); }
    void mark_boundary(int marker, const std::vector<bool>& mask) {
      triangulation_.mark_boundary(marker, [&](const auto& edge) { return mask[edge.id()]; });
    }
    void clear_boundary_markers() { triangulation_.clear_boundary_markers(); }
    // all indices of boundary edges having marker = m
    std::vector<int> filter_boundary_by_marker(int m) const {
        std::vector<int> edge_ids;
        if (!triangulation_.edges_markers().empty()) {
            for (int i = 0, n = triangulation_.n_edges(); i < n; ++i) {
                if (triangulation_.is_edge_on_boundary(i) && triangulation_.edges_markers()[i] == m) {
                    edge_ids.push_back(i);
                }
            }
        }
        return edge_ids;
    }
    // edge properties
    double_mtx edge_coords(int edge_id) {
        double_mtx coords(n_nodes_per_edge, embed_dim);
        for (int i = 0; i < n_nodes_per_edge; ++i) {
            coords.row(i) = triangulation_.nodes().row(triangulation_.edges()(edge_id, i));
        }
        return coords;
    }
    // point location
    Eigen::Matrix<int, Dynamic, 1> locate(const double_mtx& points) { return triangulation_.locate(points); }
    ~Triangulation() = default;
};

}   // namespace r
}   // namespace fdapde

#endif   // __R_MESH_H__
