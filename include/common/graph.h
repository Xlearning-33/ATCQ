#pragma once

#include "common/types.h"

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace tacc {

struct BipartiteSubgraph {
    std::vector<std::vector<VertexID>> u_adj;
    std::vector<std::vector<VertexID>> v_adj;
    std::vector<std::uint8_t> u_selected;
    std::vector<std::uint8_t> v_selected;

    BipartiteSubgraph(std::size_t num_u = 0, std::size_t num_v = 0)
        : u_adj(num_u),
          v_adj(num_v),
          u_selected(num_u, 0),
          v_selected(num_v, 0) {}
};

struct WindowAdjacencyStats {
    std::size_t incident_work_u{};
    std::size_t incident_work_v{};
    std::size_t unique_edges_examined{};
    std::size_t window_edges{};
    bool traversed_from_u{};
    bool adaptive_seed_recovery{};
    bool reversed_seed_scan{};
    std::size_t seed_scan_work{};
    std::size_t opposite_scan_work{};
};

struct GroupedVTimeAdjacency {
    std::vector<std::size_t> v_groups;
    std::vector<std::size_t> starts;
    std::vector<Timestamp> times;
    std::vector<VertexID> neighbors;
    double build_seconds{};
    std::size_t scratch_peak_bytes{};

    std::size_t bytes() const;
    std::pair<std::size_t, std::size_t> range(VertexID v, TimeWindow window) const;
};

class TemporalBipartiteGraph {
public:
    bool load_from_files(const std::string& ulabel_file,
                         const std::string& vlabel_file,
                         const std::string& edge_file);

    std::size_t num_u_vertices() const;
    std::size_t num_v_vertices() const;
    std::size_t num_edges() const;
    std::size_t num_attributes() const;

    Timestamp min_timestamp() const;
    Timestamp max_timestamp() const;

    AttributeID u_attribute(VertexID u) const;
    AttributeID v_attribute(VertexID v) const;
    std::optional<AttributeID> attribute_id(const std::string& name) const;
    const std::string& attribute_name(AttributeID id) const;

    VertexID original_u_id(VertexID u) const;
    VertexID original_v_id(VertexID v) const;

    const std::vector<TemporalEdge>& edges() const;
    const std::vector<Timestamp>& time_points() const;

    WindowAdjacencyStats populate_window_adjacency(
        BipartiteSubgraph& subgraph,
        TimeWindow window,
        const std::vector<VertexID>& active_u,
        const std::vector<VertexID>& active_v) const;
    WindowAdjacencyStats populate_seed_window_adjacency(
        BipartiteSubgraph& subgraph,
        TimeWindow window,
        const std::vector<VertexID>& seeds,
        bool seed_is_u,
        AttributeMask other_attributes,
        std::vector<VertexID>& other_vertices,
        bool adaptive = true) const;
    std::size_t attribute_incident_work(bool side_is_u,
                                        AttributeMask attributes) const;
    std::size_t num_unique_edges() const;
    std::size_t num_unique_temporal_edges() const;
    VertexID unique_edge_u(std::size_t edge_id) const;
    VertexID unique_edge_v(std::size_t edge_id) const;
    std::size_t unique_edge_timestamp_offset(std::size_t edge_id) const;
    std::uint32_t unique_edge_timestamp_count(std::size_t edge_id) const;
    Timestamp unique_edge_timestamp(std::size_t edge_id,
                                    std::size_t offset) const;
    std::size_t static_u_degree(VertexID u) const;
    std::size_t static_v_degree(VertexID v) const;
    VertexID static_u_neighbor(VertexID u, std::size_t offset) const;
    VertexID static_v_neighbor(VertexID v, std::size_t offset) const;
    std::size_t temporal_query_aux_bytes() const;
    const GroupedVTimeAdjacency& grouped_v_time_adjacency() const;
    std::size_t grouped_v_time_aux_bytes() const;

private:
    friend std::uint64_t graph_fingerprint(const TemporalBipartiteGraph& graph);

    using UniqueEdgeID = std::uint32_t;
    struct UniqueTemporalEdge {
        VertexID u{};
        VertexID v{};
        std::size_t timestamp_offset{};
        std::uint32_t timestamp_count{};
    };

    bool rebuild_temporal_adjacency();
    bool active_in_window(const UniqueTemporalEdge& edge,
                          TimeWindow window) const;

    std::size_t num_u_vertices_{};
    std::size_t num_v_vertices_{};
    Timestamp min_timestamp_{};
    Timestamp max_timestamp_{};

    std::vector<TemporalEdge> edges_;
    std::vector<Timestamp> time_points_;
    std::vector<UniqueTemporalEdge> unique_edges_;
    std::vector<Timestamp> unique_edge_timestamps_;
    std::vector<std::size_t> u_edge_offsets_;
    std::vector<std::size_t> v_edge_offsets_;
    std::vector<UniqueEdgeID> u_incident_edges_;
    std::vector<UniqueEdgeID> v_incident_edges_;
    std::vector<std::size_t> u_attribute_incident_work_;
    std::vector<std::size_t> v_attribute_incident_work_;
    std::vector<AttributeID> u_attributes_;
    std::vector<AttributeID> v_attributes_;
    std::vector<VertexID> original_u_ids_;
    std::vector<VertexID> original_v_ids_;

    std::unordered_map<VertexID, VertexID> u_id_to_internal_;
    std::unordered_map<VertexID, VertexID> v_id_to_internal_;
    std::unordered_map<std::string, AttributeID> attribute_name_to_id_;
    std::vector<std::string> attribute_names_;
    mutable std::shared_ptr<const GroupedVTimeAdjacency> grouped_v_time_;
    mutable std::shared_ptr<const std::uint64_t> fingerprint_;
};

}
