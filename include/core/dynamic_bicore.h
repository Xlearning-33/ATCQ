#pragma once

#include "common/graph.h"
#include "common/types.h"
#include "core/bicore_decomposition.h"
#include "index/index_types.h"

#include <cstddef>
#include <unordered_set>
#include <vector>

namespace tacc {

struct UFrontierChange {
    VertexID vertex{};
    Count alpha{};
    Count old_beta_max{};
    Count new_beta_max{};
};

struct VFrontierChange {
    VertexID vertex{};
    Count beta{};
    Count old_alpha_max{};
    Count new_alpha_max{};
};

struct DynamicBicoreUpdate {
    std::vector<VertexID> changed_u;
    std::vector<VertexID> changed_v;
    std::vector<UFrontierChange> u_frontier_changes;
    std::vector<VFrontierChange> v_frontier_changes;
    std::size_t affected_u{};
    std::size_t affected_v{};
    std::size_t changed_memberships{};
};

struct DynamicBicoreProfile {
    double axis_seconds{};
    std::size_t axis_adjacency_entries{};
    std::size_t axis_frontier_changes{};
    double layered_seconds{};
    std::size_t layered_sweeps{};
    std::size_t layered_region_vertices{};
    std::size_t layered_removed_vertices{};
    std::size_t layered_adjacency_entries{};
    std::size_t layered_anchor_edges{};
    double candidate_seconds{};
    double support_seconds{};
    double cascade_seconds{};
    double merge_seconds{};
    std::size_t candidate_planes{};
    std::size_t seeded_planes{};
    std::size_t changed_planes{};
    std::size_t support_adjacency_entries{};
    std::size_t cascade_queue_pops{};
    std::size_t cascade_removed_vertices{};
};

class DynamicBicoreMaintainer {
public:
    explicit DynamicBicoreMaintainer(BipartiteSubgraph graph);

    bool contains_edge(VertexID u, VertexID v) const;
    DynamicBicoreUpdate delete_edge_bucket(
        const std::vector<std::pair<VertexID, VertexID>>& edges,
        DynamicBicoreProfile* profile = nullptr);

    Count max_alpha() const { return skyline_.max_alpha; }
    Count max_beta() const { return skyline_.max_beta; }
    const BicoreSkyline& skyline() const { return skyline_; }
    std::size_t state_count() const { return frontier_entries_; }
    std::size_t membership_count() const { return membership_entries_; }
    std::size_t approximate_bytes() const;
    bool layered_enabled() const { return layered_enabled_; }
    std::size_t layered_scratch_peak_bytes() const { return layered_scratch_peak_bytes_; }

private:
    DynamicBicoreUpdate delete_edge(VertexID u, VertexID v);
    bool u_in_core(const BicoreSkyline& skyline, VertexID u,
                   Count alpha, Count beta) const;
    bool v_in_core(const BicoreSkyline& skyline, VertexID v,
                   Count alpha, Count beta) const;
    Count u_support(const BicoreSkyline& skyline, VertexID u,
                    Count alpha, Count beta,
                    const std::unordered_set<VertexID>& removed_v) const;
    Count v_support(const BicoreSkyline& skyline, VertexID v,
                    Count alpha, Count beta,
                    const std::unordered_set<VertexID>& removed_u) const;
    bool deactivate_edge(VertexID u, VertexID v);
    void record_frontier(bool is_u, VertexID id, Count parameter, Count value,
                         std::vector<std::pair<VertexID, Count>>& touched_u,
                         std::vector<std::pair<VertexID, Count>>& touched_v);
    Count next_unwritten(bool is_u, VertexID id, Count parameter);
    std::uint32_t next_plane_epoch();
    void maintain_layered(
        const std::vector<std::pair<VertexID, VertexID>>& edges,
        std::vector<std::pair<VertexID, Count>>& touched_u,
        std::vector<std::pair<VertexID, Count>>& touched_v,
        DynamicBicoreProfile* profile);
    BicoreSkyline skyline_;
    std::vector<std::vector<VertexID>> u_neighbors_;
    std::vector<std::vector<VertexID>> v_neighbors_;
    std::vector<std::vector<std::uint8_t>> u_edge_active_;
    std::vector<std::vector<std::uint8_t>> v_edge_active_;
    std::vector<Count> live_u_degree_;
    std::vector<Count> live_v_degree_;
    std::vector<std::uint32_t> axis_u_epoch_;
    std::vector<std::uint32_t> axis_v_epoch_;
    std::vector<std::uint32_t> removed_u_epoch_;
    std::vector<std::uint32_t> removed_v_epoch_;
    std::vector<std::uint32_t> queued_u_epoch_;
    std::vector<std::uint32_t> queued_v_epoch_;
    std::vector<std::uint32_t> support_u_epoch_;
    std::vector<std::uint32_t> support_v_epoch_;
    std::vector<Count> support_u_value_;
    std::vector<Count> support_v_value_;
    std::vector<Count> endpoint_u_beta_max_;
    std::vector<Count> endpoint_v_beta_max_;
    std::vector<std::size_t> u_frontier_offsets_;
    std::vector<std::size_t> v_frontier_offsets_;
    std::vector<std::uint32_t> pending_u_epoch_;
    std::vector<std::uint32_t> pending_v_epoch_;
    std::vector<Count> pending_u_new_beta_;
    std::vector<Count> pending_v_new_alpha_;
    std::vector<Count> successor_u_;
    std::vector<Count> successor_v_;
    bool layered_enabled_{};
    Count delta_bound_{};
    std::size_t layered_scratch_peak_bytes_{};
    std::uint32_t plane_epoch_{};
    std::uint32_t bucket_epoch_{};
    std::size_t frontier_entries_{};
    std::size_t membership_entries_{};
};

}
