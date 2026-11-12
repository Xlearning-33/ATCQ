#pragma once

#include <algorithm>
#include <cstddef>

namespace tacc {

struct QueryProfile {
    double attribute_filter_us{};
    double directory_lookup_us{};
    double candidate_collect_us{};
    double candidate_intersection_us{};
    double window_recovery_us{};
    double peeling_us{};
    double grouped_aux_build_us{};
    double workspace_clear_us{};

    std::size_t directory_keys{};
    std::size_t directory_hits{};
    std::size_t posting_entries{};
    std::size_t vertex_refs{};
    std::size_t directional_candidates{};
    std::size_t candidate_u{};
    std::size_t candidate_v{};
    std::size_t incident_work_u{};
    std::size_t incident_work_v{};
    std::size_t unique_edges_examined{};
    std::size_t window_edges{};
    std::size_t result_u{};
    std::size_t result_v{};
    std::size_t online_fallbacks{};
    std::size_t adaptive_seed_recoveries{};
    std::size_t reversed_seed_scans{};
    std::size_t seed_scan_work{};
    std::size_t opposite_scan_work{};
    std::size_t grouped_recoveries{};
    std::size_t grouped_window_events{};
    std::size_t grouped_duplicate_events{};
    std::size_t grouped_budget_fallbacks{};
    std::size_t grouped_aux_bytes{};
    std::size_t workspace_fixed_bytes{};
    std::size_t workspace_retained_bytes{};
};

inline QueryProfile& operator+=(QueryProfile& lhs, const QueryProfile& rhs) {
    lhs.attribute_filter_us += rhs.attribute_filter_us;
    lhs.directory_lookup_us += rhs.directory_lookup_us;
    lhs.candidate_collect_us += rhs.candidate_collect_us;
    lhs.candidate_intersection_us += rhs.candidate_intersection_us;
    lhs.window_recovery_us += rhs.window_recovery_us;
    lhs.peeling_us += rhs.peeling_us;
    lhs.grouped_aux_build_us += rhs.grouped_aux_build_us;
    lhs.workspace_clear_us += rhs.workspace_clear_us;
    lhs.directory_keys += rhs.directory_keys;
    lhs.directory_hits += rhs.directory_hits;
    lhs.posting_entries += rhs.posting_entries;
    lhs.vertex_refs += rhs.vertex_refs;
    lhs.directional_candidates += rhs.directional_candidates;
    lhs.candidate_u += rhs.candidate_u;
    lhs.candidate_v += rhs.candidate_v;
    lhs.incident_work_u += rhs.incident_work_u;
    lhs.incident_work_v += rhs.incident_work_v;
    lhs.unique_edges_examined += rhs.unique_edges_examined;
    lhs.window_edges += rhs.window_edges;
    lhs.result_u += rhs.result_u;
    lhs.result_v += rhs.result_v;
    lhs.online_fallbacks += rhs.online_fallbacks;
    lhs.adaptive_seed_recoveries += rhs.adaptive_seed_recoveries;
    lhs.reversed_seed_scans += rhs.reversed_seed_scans;
    lhs.seed_scan_work += rhs.seed_scan_work;
    lhs.opposite_scan_work += rhs.opposite_scan_work;
    lhs.grouped_recoveries += rhs.grouped_recoveries;
    lhs.grouped_window_events += rhs.grouped_window_events;
    lhs.grouped_duplicate_events += rhs.grouped_duplicate_events;
    lhs.grouped_budget_fallbacks += rhs.grouped_budget_fallbacks;
    lhs.grouped_aux_bytes = std::max(lhs.grouped_aux_bytes, rhs.grouped_aux_bytes);
    lhs.workspace_fixed_bytes = std::max(lhs.workspace_fixed_bytes, rhs.workspace_fixed_bytes);
    lhs.workspace_retained_bytes = std::max(lhs.workspace_retained_bytes, rhs.workspace_retained_bytes);
    return lhs;
}

}
