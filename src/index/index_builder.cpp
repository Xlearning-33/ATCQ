#include "index/index_builder.h"

#include "core/bicore_decomposition.h"
#include "core/dynamic_bicore.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <type_traits>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tacc {
namespace {

using Clock = std::chrono::steady_clock;


AttributeMask u_universe(const TemporalBipartiteGraph& graph) {
    AttributeMask universe = 0;
    for (VertexID u = 0; u < graph.num_u_vertices(); ++u)
        universe |= AttributeMask{1} << graph.u_attribute(u);
    return universe;
}
struct AttributeMaskPlan {
    std::vector<AttributeMask> masks;
    std::size_t total{};
    std::size_t no_edge_pruned{};
    std::size_t equivalent_pruned{};
};

bool attribute_pruning_enabled() {
    const char* value = std::getenv("TAC_ATTRIBUTE_PRUNING");
    if (value == nullptr) return false;
    const std::string text(value);
    return text == "1" || text == "true" || text == "on";
}
bool dynamic_profile_enabled() {
    const char* value = std::getenv("TAC_PROFILE_DYNAMIC");
    if (value == nullptr) return false;
    const std::string text(value);
    return text == "1" || text == "true" || text == "on";
}
bool dynamic_bucket_check_enabled() {
    const char* value = std::getenv("TAC_VERIFY_DYNAMIC_BUCKET");
    if (value == nullptr) return false;
    const std::string text(value);
    return text == "1" || text == "true" || text == "on";
}

bool mat_boundary_only_enabled() {
    static const bool enabled = [] {
        const char* value = std::getenv("TAC_MAT_BOUNDARY_ONLY");
        if (value == nullptr) return true;
        const std::string text(value);
        return text != "0" && text != "false" && text != "off";
    }();
    return enabled;
}

std::vector<AttributeMask> all_nonempty_subsets(AttributeMask universe) {
    std::vector<AttributeMask> result;
    for (AttributeMask mask = universe; mask != 0;
         mask = (mask - 1) & universe) {
        result.push_back(mask);
    }
    std::sort(result.begin(), result.end(),
              [](AttributeMask lhs, AttributeMask rhs) {
                  const int lhs_size = __builtin_popcountll(lhs);
                  const int rhs_size = __builtin_popcountll(rhs);
                  if (lhs_size != rhs_size) return lhs_size < rhs_size;
                  return lhs < rhs;
              });
    return result;
}

using TemporalEdgeSignature =
    std::tuple<VertexID, VertexID, Timestamp>;

AttributeMaskPlan make_attribute_plan(const TemporalBipartiteGraph& graph,
                                      bool enable_pruning) {
    AttributeMaskPlan plan;
    plan.masks = all_nonempty_subsets(u_universe(graph));
    plan.total = plan.masks.size();
    if (!enable_pruning) return plan;

    std::sort(plan.masks.begin(), plan.masks.end(),
              [](AttributeMask lhs, AttributeMask rhs) {
                  const int lhs_size = __builtin_popcountll(lhs);
                  const int rhs_size = __builtin_popcountll(rhs);
                  if (lhs_size != rhs_size) return lhs_size < rhs_size;
                  return lhs < rhs;
              });

    std::map<std::vector<TemporalEdgeSignature>,
             std::vector<AttributeMask>> representatives;
    std::vector<AttributeMask> kept;
    for (AttributeMask mask : plan.masks) {
        std::vector<TemporalEdgeSignature> signature;
        for (const TemporalEdge& edge : graph.edges()) {
            const AttributeID attribute = graph.u_attribute(edge.u);
            if ((mask & (AttributeMask{1} << attribute)) != 0) {
                signature.emplace_back(edge.u, edge.v, edge.timestamp);
            }
        }
        if (signature.empty()) {
            ++plan.no_edge_pruned;
            continue;
        }

        bool covered = false;
        auto& same_signature = representatives[signature];
        for (AttributeMask representative : same_signature) {
            if (attribute_subset(representative, mask)) {
                covered = true;
                break;
            }
        }
        if (covered) {
            ++plan.equivalent_pruned;
            continue;
        }
        same_signature.push_back(mask);
        kept.push_back(mask);
    }
    plan.masks = std::move(kept);
    return plan;
}

struct TemporalProjection {
    std::vector<std::pair<VertexID, VertexID>> pairs;
    std::vector<std::vector<std::uint32_t>> pair_ids_by_time;
    std::vector<std::uint32_t> full_counts;
};

BipartiteSubgraph snapshot_from_counts(
    const TemporalBipartiteGraph& graph, AttributeMask mask, const TemporalProjection& projection,
    const std::vector<std::uint32_t>& counts) {
    BipartiteSubgraph snapshot(graph.num_u_vertices(), graph.num_v_vertices());
    std::fill(snapshot.u_selected.begin(), snapshot.u_selected.end(), 1);
    std::fill(snapshot.v_selected.begin(), snapshot.v_selected.end(), 1);

        for (VertexID u = 0; u < graph.num_u_vertices(); ++u) {
            snapshot.u_selected[u] =
                (mask & (AttributeMask{1} << graph.u_attribute(u))) != 0;
        }
    for (std::size_t pair_id = 0; pair_id < counts.size(); ++pair_id) {
        if (counts[pair_id] == 0) continue;
        const auto [u, v] = projection.pairs[pair_id];
        snapshot.u_adj[u].push_back(v);
        snapshot.v_adj[v].push_back(u);
    }
    return snapshot;
}

struct AttributeProjectionCatalog {
    std::vector<std::uint32_t> time_indices;
    std::vector<std::vector<std::uint32_t>> u_edge_ids_by_attribute;
};

AttributeProjectionCatalog make_projection_catalog(
    const TemporalBipartiteGraph& graph) {
    const auto& times = graph.time_points();
    if (graph.num_unique_edges() >
            std::numeric_limits<std::uint32_t>::max() ||
        times.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw std::runtime_error(
            "projection catalog requires 32-bit edge and time IDs");
    }

    AttributeProjectionCatalog catalog;
    catalog.time_indices.resize(graph.num_unique_temporal_edges());
    catalog.u_edge_ids_by_attribute.resize(graph.num_attributes());
    for (std::size_t edge_id = 0;
         edge_id < graph.num_unique_edges(); ++edge_id) {
        const VertexID u = graph.unique_edge_u(edge_id);
        const auto compact_edge_id = static_cast<std::uint32_t>(edge_id);
        catalog.u_edge_ids_by_attribute[graph.u_attribute(u)].push_back(
            compact_edge_id);

        const std::size_t timestamp_offset =
            graph.unique_edge_timestamp_offset(edge_id);
        const std::uint32_t timestamp_count =
            graph.unique_edge_timestamp_count(edge_id);
        for (std::uint32_t i = 0; i < timestamp_count; ++i) {
            const Timestamp timestamp =
                graph.unique_edge_timestamp(edge_id, i);
            const auto found = std::lower_bound(
                times.begin(), times.end(), timestamp);
            if (found == times.end() || *found != timestamp) {
                throw std::runtime_error(
                    "unique-edge timestamp missing from time points");
            }
            catalog.time_indices[timestamp_offset + i] =
                static_cast<std::uint32_t>(found - times.begin());
        }
    }
    return catalog;
}

TemporalProjection make_temporal_projection(
    const TemporalBipartiteGraph& graph,
    const AttributeProjectionCatalog& catalog,
    AttributeMask mask) {
    TemporalProjection projection;
    projection.pair_ids_by_time.resize(graph.time_points().size());
    const auto& buckets = catalog.u_edge_ids_by_attribute;

    std::size_t pair_count = 0;
    for (std::size_t attribute = 0; attribute < buckets.size(); ++attribute) {
        if ((mask & (AttributeMask{1} << attribute)) != 0) {
            pair_count += buckets[attribute].size();
        }
    }
    projection.pairs.reserve(pair_count);
    projection.full_counts.reserve(pair_count);

    for (std::size_t attribute = 0; attribute < buckets.size(); ++attribute) {
        if ((mask & (AttributeMask{1} << attribute)) == 0) continue;
        for (std::uint32_t edge_id : buckets[attribute]) {
            const std::uint32_t pair_id = static_cast<std::uint32_t>(
                projection.pairs.size());
            projection.pairs.emplace_back(
                graph.unique_edge_u(edge_id), graph.unique_edge_v(edge_id));
            const std::uint32_t timestamp_count =
                graph.unique_edge_timestamp_count(edge_id);
            projection.full_counts.push_back(timestamp_count);
            const std::size_t timestamp_offset =
                graph.unique_edge_timestamp_offset(edge_id);
            for (std::uint32_t i = 0; i < timestamp_count; ++i) {
                projection.pair_ids_by_time[
                    catalog.time_indices[timestamp_offset + i]].push_back(
                        pair_id);
            }
        }
    }
    return projection;
}

Count u_frontier(const BicoreSkyline& skyline, VertexID u, Count alpha) {
    if (u >= skyline.u_beta_max.size() ||
        alpha >= skyline.u_beta_max[u].size()) {
        return 0;
    }
    return skyline.u_beta_max[u][alpha];
}

Count v_frontier(const BicoreSkyline& skyline, VertexID v, Count beta) {
    if (v >= skyline.v_alpha_max.size() ||
        beta >= skyline.v_alpha_max[v].size()) {
        return 0;
    }
    return skyline.v_alpha_max[v][beta];
}

std::size_t skyline_membership_count(const BicoreSkyline& skyline) {
    std::size_t count = 0;
    for (const auto& row : skyline.u_beta_max) {
        for (std::size_t alpha = 1; alpha < row.size(); ++alpha) {
            count += row[alpha];
            if (alpha == 1 && row[alpha] != 0) --count;
        }
    }
    for (const auto& row : skyline.v_alpha_max) {
        for (std::size_t beta = 1; beta < row.size(); ++beta) {
            count += row[beta];
            if (beta == 1 && row[beta] != 0) --count;
        }
    }
    return count;
}

void verify_dynamic_bucket(const BicoreSkyline& actual,
                           const BicoreSkyline& expected,
                           AttributeMask mask, Timestamp start,
                           Timestamp removed_time) {
    const std::string context =
        std::string("dynamic bucket mismatch side=u") +
        " mask=" + std::to_string(mask) +
        " start=" + std::to_string(start) +
        " removed_time=" + std::to_string(removed_time);
    if (actual.u_beta_max.size() != expected.u_beta_max.size() ||
        actual.v_alpha_max.size() != expected.v_alpha_max.size()) {
        throw std::runtime_error(context + " vertex counts");
    }
    for (VertexID u = 0; u < actual.u_beta_max.size(); ++u) {
        const std::size_t length = std::max(
            actual.u_beta_max[u].size(), expected.u_beta_max[u].size());
        for (std::size_t alpha = 1; alpha < length; ++alpha) {
            if (u_frontier(actual, u, static_cast<Count>(alpha)) !=
                u_frontier(expected, u, static_cast<Count>(alpha))) {
                throw std::runtime_error(context + " u=" +
                    std::to_string(u) + " alpha=" +
                    std::to_string(alpha));
            }
        }
    }
    for (VertexID v = 0; v < actual.v_alpha_max.size(); ++v) {
        const std::size_t length = std::max(
            actual.v_alpha_max[v].size(), expected.v_alpha_max[v].size());
        for (std::size_t beta = 1; beta < length; ++beta) {
            if (v_frontier(actual, v, static_cast<Count>(beta)) !=
                v_frontier(expected, v, static_cast<Count>(beta))) {
                throw std::runtime_error(context + " v=" +
                    std::to_string(v) + " beta=" +
                    std::to_string(beta));
            }
        }
    }
}

template <typename Index>
void emit_frontier_growth(Index& index,
                          AttributeMask mask, Timestamp start, Timestamp end,
                          const BicoreSkyline& previous,
                          const BicoreSkyline& current) {
    const auto streaming_start = Clock::now();
    const bool boundary_only =
        std::is_same_v<Index, MATIndex> && mat_boundary_only_enabled();
    for (VertexID v = 0; v < current.v_alpha_max.size(); ++v) {
        for (Count beta = 1; beta < current.v_alpha_max[v].size(); ++beta) {
            const Count old_alpha = v_frontier(previous, v, beta);
            const Count new_alpha = current.v_alpha_max[v][beta];
            if (boundary_only) {
                if (new_alpha > old_alpha && !(new_alpha == 1 && beta == 1)) {
                    index.add({v, mask, start, end,
                               new_alpha, beta});
                }
                continue;
            }
            for (Count alpha = old_alpha + 1; alpha <= new_alpha; ++alpha) {
                if (alpha == 1 && beta == 1) continue;
                index.add({v, mask, start, end,
                           alpha, beta});
            }
        }
    }
    index.mutable_stats().streaming_seconds +=
        std::chrono::duration<double>(Clock::now() - streaming_start).count();
}

template <typename Index>
void emit_frontier_drop(Index& index,
                        AttributeMask mask, Timestamp start, Timestamp end,
                        const DynamicBicoreUpdate& update) {
    const auto streaming_start = Clock::now();
    const bool boundary_only =
        std::is_same_v<Index, MATIndex> && mat_boundary_only_enabled();
    for (const VFrontierChange& change : update.v_frontier_changes) {
        if (boundary_only) {
            if (change.old_alpha_max > change.new_alpha_max &&
                !(change.old_alpha_max == 1 && change.beta == 1)) {
                index.add({change.vertex, mask,
                           start, end, change.old_alpha_max, change.beta});
            }
            continue;
        }
        for (Count alpha = change.new_alpha_max + 1;
             alpha <= change.old_alpha_max; ++alpha) {
            if (alpha == 1 && change.beta == 1) continue;
            index.add({change.vertex, mask,
                       start, end, alpha, change.beta});
        }
    }
    index.mutable_stats().streaming_seconds +=
        std::chrono::duration<double>(Clock::now() - streaming_start).count();
}

template <typename Index>
void build_basic(
    Index& index, const TemporalBipartiteGraph& graph,
    const AttributeProjectionCatalog& catalog,
    const std::vector<AttributeMask>& masks,
    const std::function<void(AttributeMask)>& complete_attribute_set) {
    const auto& times = graph.time_points();
    IndexBuildStats& stats = index.mutable_stats();
    for (AttributeMask mask : masks) {
        bool has_record = false;
        const TemporalProjection projection =
            make_temporal_projection(graph, catalog, mask);
        for (std::size_t s = 0; s < times.size(); ++s) {
            std::vector<std::uint32_t> counts(projection.pairs.size(), 0);
            BicoreSkyline previous;
            for (std::size_t e = s; e < times.size(); ++e) {
                ++stats.windows;
                for (std::uint32_t pair_id :
                     projection.pair_ids_by_time[e]) {
                    ++counts[pair_id];
                }
                const BipartiteSubgraph snapshot =
                    snapshot_from_counts(graph, mask,
                                         projection, counts);
                const auto decomposition_start = Clock::now();
                BicoreSkyline current = bicore_skyline_decomposition(snapshot);
                stats.decomposition_seconds += std::chrono::duration<double>(
                    Clock::now() - decomposition_start).count();
                ++stats.decompositions;
                stats.full_window_memberships +=
                    skyline_membership_count(current);
                const std::size_t before = stats.generated_records;
                emit_frontier_growth(index, mask, times[s], times[e],
                                     previous, current);
                has_record = has_record || stats.generated_records > before;
                previous = std::move(current);
            }
        }
        if (has_record) ++stats.attribute_sets_nonempty;
        else ++stats.attribute_sets_empty;
        complete_attribute_set(mask);
    }
}

template <typename Index>
void build_opt(
    Index& index, const TemporalBipartiteGraph& graph,
    const AttributeProjectionCatalog& catalog,
    const std::vector<AttributeMask>& masks,
    const std::function<void(AttributeMask)>& complete_attribute_set) {
    const auto& times = graph.time_points();
    IndexBuildStats& stats = index.mutable_stats();
    const bool measure_dynamic = dynamic_profile_enabled();
    const bool verify_dynamic = dynamic_bucket_check_enabled();
    stats.dynamic_profile_enabled = measure_dynamic;
    DynamicBicoreProfile dynamic_profile;
    for (AttributeMask mask : masks) {
        bool has_record = false;
        const TemporalProjection projection =
            make_temporal_projection(graph, catalog, mask);
        std::vector<std::uint32_t> suffix_counts = projection.full_counts;
        for (std::size_t s = 0; s < times.size(); ++s) {
            auto counts = suffix_counts;
            stats.inserted_unique_edges += static_cast<std::size_t>(
                std::count_if(counts.begin(), counts.end(),
                              [](std::uint32_t count) { return count != 0; }));

            BipartiteSubgraph broad_snapshot =
                snapshot_from_counts(graph, mask,
                                     projection, counts);
            const auto initialization_start = Clock::now();
            DynamicBicoreMaintainer maintainer(std::move(broad_snapshot));
            stats.dynamic_layered_enabled = maintainer.layered_enabled();
            stats.dynamic_kernel_known = true;
            const auto maintainer_bytes = maintainer.approximate_bytes();
            stats.dynamic_peak_bytes = std::max(
                stats.dynamic_peak_bytes, maintainer_bytes);
            stats.decomposition_seconds += std::chrono::duration<double>(
                Clock::now() - initialization_start).count();
            ++stats.decompositions;

            for (std::size_t e = times.size(); e-- > s;) {
                ++stats.windows;
                stats.full_window_memberships += maintainer.membership_count();
                std::vector<std::pair<VertexID, VertexID>> disappeared;
                disappeared.reserve(
                    projection.pair_ids_by_time[e].size());
                for (std::uint32_t pair_id :
                     projection.pair_ids_by_time[e]) {
                    if (counts[pair_id] == 0 || --counts[pair_id] != 0) {
                        continue;
                    }
                    disappeared.push_back(projection.pairs[pair_id]);
                }
                stats.deleted_unique_edges += disappeared.size();
                stats.dynamic_edge_updates += disappeared.size();
                const auto update_start = Clock::now();
                const DynamicBicoreUpdate update =
                    maintainer.delete_edge_bucket(
                        disappeared, measure_dynamic ? &dynamic_profile : nullptr);
                stats.dynamic_update_seconds += std::chrono::duration<double>(
                    Clock::now() - update_start).count();
                stats.dynamic_peak_bytes = std::max(stats.dynamic_peak_bytes,
                    maintainer_bytes + maintainer.layered_scratch_peak_bytes());
                if (verify_dynamic) {
                    const BicoreSkyline expected =
                        bicore_skyline_decomposition(
                            snapshot_from_counts(
                                graph, mask, projection, counts));
                    verify_dynamic_bucket(maintainer.skyline(), expected,
                                          mask, times[s], times[e]);
                }

                const std::size_t before = stats.generated_records;
                emit_frontier_drop(index, mask, times[s], times[e],
                                   update);
                has_record = has_record || stats.generated_records > before;

                const std::size_t changed_vertices =
                    update.changed_u.size() + update.changed_v.size();
                stats.affected_u_vertices_total += update.affected_u;
                stats.affected_v_vertices_total += update.affected_v;
                stats.max_affected_region_vertices = std::max(
                    stats.max_affected_region_vertices,
                    update.affected_u + update.affected_v);
                stats.changed_u_vertices += update.changed_u.size();
                stats.changed_v_vertices += update.changed_v.size();
                stats.changed_frontier_entries +=
                    update.u_frontier_changes.size() +
                    update.v_frontier_changes.size();
                stats.max_changed_vertices_per_bucket = std::max(
                    stats.max_changed_vertices_per_bucket, changed_vertices);
                if (changed_vertices != 0) ++stats.nonzero_changed_buckets;
            }

            for (std::uint32_t pair_id :
                 projection.pair_ids_by_time[s]) {
                if (suffix_counts[pair_id] != 0) --suffix_counts[pair_id];
            }
        }
        if (has_record) ++stats.attribute_sets_nonempty;
        else ++stats.attribute_sets_empty;
        complete_attribute_set(mask);
    }
    stats.dynamic_candidate_seconds += dynamic_profile.candidate_seconds;
    stats.dynamic_axis_seconds += dynamic_profile.axis_seconds;
    stats.dynamic_axis_adjacency_entries += dynamic_profile.axis_adjacency_entries;
    stats.dynamic_axis_frontier_changes += dynamic_profile.axis_frontier_changes;
    stats.dynamic_layered_seconds += dynamic_profile.layered_seconds;
    stats.dynamic_layered_sweeps += dynamic_profile.layered_sweeps;
    stats.dynamic_layered_region_vertices += dynamic_profile.layered_region_vertices;
    stats.dynamic_layered_removed_vertices += dynamic_profile.layered_removed_vertices;
    stats.dynamic_layered_adjacency_entries += dynamic_profile.layered_adjacency_entries;
    stats.dynamic_layered_anchor_edges += dynamic_profile.layered_anchor_edges;
    stats.dynamic_support_seconds += dynamic_profile.support_seconds;
    stats.dynamic_cascade_seconds += dynamic_profile.cascade_seconds;
    stats.dynamic_merge_seconds += dynamic_profile.merge_seconds;
    stats.dynamic_candidate_planes += dynamic_profile.candidate_planes;
    stats.dynamic_seeded_planes += dynamic_profile.seeded_planes;
    stats.dynamic_changed_planes += dynamic_profile.changed_planes;
    stats.dynamic_support_adjacency_entries +=
        dynamic_profile.support_adjacency_entries;
    stats.dynamic_cascade_queue_pops += dynamic_profile.cascade_queue_pops;
    stats.dynamic_cascade_removed_vertices +=
        dynamic_profile.cascade_removed_vertices;
}
template <typename Index>
Index build_index(const TemporalBipartiteGraph& graph, BuildMethod method,
                  const IndexBuildProgressCallback& progress) {
    if (graph.num_attributes() > 64) {
        throw std::runtime_error("formal indexes support at most 64 attributes");
    }
    Index index;
    IndexBuildStats& stats = index.mutable_stats();
    stats.scope = configured_index_scope();
    if constexpr (std::is_same_v<Index, VATIndex>)
        index.initialize_vertices(graph.num_v_vertices());
    stats.method = method;
    stats.attribute_pruning_enabled = attribute_pruning_enabled();
    const auto build_start = Clock::now();
    const auto planning_start = Clock::now();
    const AttributeMaskPlan u_plan = make_attribute_plan(
        graph, stats.attribute_pruning_enabled);
    stats.attribute_planning_seconds = std::chrono::duration<double>(
        Clock::now() - planning_start).count();
    for (const AttributeMaskPlan* plan : {&u_plan}) {
        stats.attribute_sets_total += plan->total;
        stats.attribute_sets_no_edge_pruned += plan->no_edge_pruned;
        stats.attribute_sets_equivalent_pruned += plan->equivalent_pruned;
    }

    const AttributeProjectionCatalog projection_catalog =
        make_projection_catalog(graph);

    const std::size_t total_attribute_sets =
        u_plan.masks.size();
    std::size_t completed_attribute_sets = 0;
    auto build_masks = [&](const AttributeMaskPlan& plan) {
        const auto complete = [&](AttributeMask mask) {
            ++completed_attribute_sets;
            if (progress) {
                progress({completed_attribute_sets, total_attribute_sets,
                          Side::U, mask,
                          std::chrono::duration<double>(
                              Clock::now() - build_start).count()});
            }
        };
        if (method == BuildMethod::Basic) {
            build_basic(index, graph, projection_catalog,
                                  plan.masks, complete);
        } else {
            build_opt(index, graph, projection_catalog,
                                plan.masks, complete);
        }
    };
    build_masks(u_plan);
    if (stats.full_window_memberships >= stats.generated_records) {
        stats.streaming_records_avoided =
            stats.full_window_memberships - stats.generated_records;
    }
    const auto finalize_start = Clock::now();
    index.finalize();
    stats.finalize_seconds += std::chrono::duration<double>(
        Clock::now() - finalize_start).count();
    stats.build_seconds = std::chrono::duration<double>(
        Clock::now() - build_start).count();
    return index;
}
}

IndexScope configured_index_scope() {
    const char* value = std::getenv("TAC_INDEX_SCOPE");
    if (value == nullptr || std::string(value) == "uv") return IndexScope::UV;
    throw std::invalid_argument("This build supports UV only; unset TAC_INDEX_SCOPE or use uv");
}

BuildMethod configured_build_method() {
    const char* value = std::getenv("TAC_INDEX_BUILD_METHOD");
    if (value == nullptr || *value == '\0' || std::string(value) == "opt") {
        return BuildMethod::Opt;
    }
    if (std::string(value) == "basic") return BuildMethod::Basic;
    throw std::invalid_argument(
        "TAC_INDEX_BUILD_METHOD must be 'basic' or 'opt'");
}

const char* build_method_name(BuildMethod method) {
    return method == BuildMethod::Basic ? "basic" : "opt";
}

VATIndex build_vat_index(
    const TemporalBipartiteGraph& graph, BuildMethod method,
    const IndexBuildProgressCallback& progress) {
    return build_index<VATIndex>(graph, method, progress);
}

UATIndex build_uat_index(
    const TemporalBipartiteGraph& graph, BuildMethod method,
    const IndexBuildProgressCallback& progress) {
    return build_index<UATIndex>(graph, method, progress);
}

MATIndex build_mat_index(
    const TemporalBipartiteGraph& graph, BuildMethod method,
    const IndexBuildProgressCallback& progress) {
    return build_index<MATIndex>(graph, method, progress);
}

}
