#pragma once

#include "common/types.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace tacc {

enum class Side : std::uint8_t { U = 0, V = 1 };

enum class IndexScope : std::uint8_t { UV = 1 };

inline const char* index_scope_name(IndexScope scope) {
    (void)scope;
    return "uv";
}

struct IndexMembership {
    VertexID vertex{};
    AttributeMask attributes{};
    Timestamp start{};
    Timestamp end{};
    Count alpha{};
    Count beta{};
};

struct IndexBlock {
    AttributeMask attributes{};
    Timestamp start{};
    Timestamp end{};
    Count alpha{};
    Count beta{};
    std::size_t vertex_offset{};
    std::size_t vertex_count{};
};

struct IndexPostingRange {
    std::size_t offset{};
    std::size_t count{};
};

enum class BuildMethod : std::uint8_t { Basic = 0, Opt = 1 };

struct MATStreamingProfile {
    bool enabled{};
    std::size_t sampled_calls{};
    std::size_t rejected_records{};
    std::size_t erased_records{};
    std::size_t reject_buckets{};
    std::size_t reject_bounds_passed{};
    std::size_t reject_comparisons{};
    std::size_t reject_time_groups{};
    std::size_t ordered_reject_buckets{};
    std::size_t start_group_build_records{};
    std::size_t start_group_capacity_bytes{};
    std::size_t remove_buckets{};
    std::size_t remove_bounds_passed{};
    std::size_t remove_comparisons{};
    std::size_t remove_time_groups{};
    std::size_t refreshed_records{};
    std::size_t max_bucket_records{};
    double sampled_add_seconds{};
    double sampled_locator_seconds{};
    double sampled_reject_seconds{};
    double sampled_remove_seconds{};
    double sampled_append_seconds{};
};

struct IndexBuildStats {
    IndexScope scope{IndexScope::UV};
    std::size_t uat_peak_flat_bytes{};
    std::size_t uat_peak_core_directory_bytes{};
    std::size_t uat_core_groups{};
    std::size_t uat_completed_directions{};
    BuildMethod method{BuildMethod::Basic};
    bool attribute_pruning_enabled{};
    bool dynamic_profile_enabled{};
    bool dynamic_layered_enabled{};
    bool dynamic_kernel_known{};
    MATStreamingProfile mat_streaming_profile;
    bool mat_ordered_streaming_enabled{};
    bool mat_streaming_mode_known{};
    double build_seconds{};
    double attribute_planning_seconds{};
    double decomposition_seconds{};
    double dynamic_update_seconds{};
    double dynamic_candidate_seconds{};
    double dynamic_axis_seconds{};
    std::size_t dynamic_axis_adjacency_entries{};
    std::size_t dynamic_axis_frontier_changes{};
    double dynamic_layered_seconds{};
    std::size_t dynamic_layered_sweeps{};
    std::size_t dynamic_layered_region_vertices{};
    std::size_t dynamic_layered_removed_vertices{};
    std::size_t dynamic_layered_adjacency_entries{};
    std::size_t dynamic_layered_anchor_edges{};
    double dynamic_support_seconds{};
    double dynamic_cascade_seconds{};
    double dynamic_merge_seconds{};
    double streaming_seconds{};
    double finalize_seconds{};
    std::size_t attribute_sets_total{};
    std::size_t attribute_sets_no_edge_pruned{};
    std::size_t attribute_sets_equivalent_pruned{};
    std::size_t attribute_sets_nonempty{};
    std::size_t attribute_sets_empty{};
    std::size_t windows{};
    std::size_t decompositions{};
    std::size_t fallback_decompositions{};
    std::size_t full_window_memberships{};
    std::size_t streaming_records_avoided{};
    std::size_t generated_records{};
    std::size_t kept_records{};
    std::size_t attribute_time_removed{};
    std::size_t cross_core_removed{};
    std::size_t grouped_blocks{};
    std::size_t vertex_references{};
    std::size_t approx_bytes{};
    std::size_t inserted_unique_edges{};
    std::size_t deleted_unique_edges{};
    std::size_t dynamic_edge_updates{};
    std::size_t dynamic_candidate_planes{};
    std::size_t dynamic_seeded_planes{};
    std::size_t dynamic_changed_planes{};
    std::size_t dynamic_support_adjacency_entries{};
    std::size_t dynamic_cascade_queue_pops{};
    std::size_t dynamic_cascade_removed_vertices{};
    std::size_t dynamic_peak_bytes{};
    std::size_t affected_u_vertices_total{};
    std::size_t affected_v_vertices_total{};
    std::size_t max_affected_region_vertices{};
    std::size_t changed_u_vertices{};
    std::size_t changed_v_vertices{};
    std::size_t changed_frontier_entries{};
    std::size_t max_changed_vertices_per_bucket{};
    std::size_t nonzero_changed_buckets{};
};

inline bool attribute_subset(AttributeMask lhs, AttributeMask rhs) {
    return (lhs & rhs) == lhs;
}

}
