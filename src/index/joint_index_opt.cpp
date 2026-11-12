#include "index/joint_index.h"
#include "core/dynamic_bicore.h"

#include <algorithm>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace tacc::joint {
namespace {
using CoreKey = std::pair<Count, Count>;

std::vector<AttributeMask> subsets(AttributeMask universe) {
    std::vector<AttributeMask> result;
    for (auto mask = universe; mask; mask = (mask - 1) & universe) result.push_back(mask);
    std::sort(result.begin(), result.end());
    return result;
}

struct Catalog {
    AttributePair all;
    std::map<std::pair<AttributeID, AttributeID>, std::vector<std::size_t>> by_attributes;
    std::vector<std::size_t> time_indices;
    explicit Catalog(const TemporalBipartiteGraph& graph) {
        const auto add = [](AttributeMask& mask, AttributeID id) {
            if (id >= 64) throw std::runtime_error("joint OPT supports global attribute IDs 0--63");
            mask |= AttributeMask{1} << id;
        };
        for (VertexID u = 0; u < graph.num_u_vertices(); ++u) add(all.u, graph.u_attribute(u));
        for (VertexID v = 0; v < graph.num_v_vertices(); ++v) add(all.v, graph.v_attribute(v));
        time_indices.resize(graph.num_unique_temporal_edges());
        const auto& times = graph.time_points();
        for (std::size_t id = 0; id < graph.num_unique_edges(); ++id) {
            by_attributes[{graph.u_attribute(graph.unique_edge_u(id)),
                           graph.v_attribute(graph.unique_edge_v(id))}].push_back(id);
            const auto offset = graph.unique_edge_timestamp_offset(id);
            for (std::uint32_t i = 0; i < graph.unique_edge_timestamp_count(id); ++i) {
                const auto time = graph.unique_edge_timestamp(id, i);
                const auto found = std::lower_bound(times.begin(), times.end(), time);
                if (found == times.end() || *found != time || offset + i >= time_indices.size())
                    throw std::runtime_error("invalid joint OPT time catalog");
                time_indices[offset + i] = static_cast<std::size_t>(found - times.begin());
            }
        }
    }
};

struct Projection {
    std::vector<std::pair<VertexID, VertexID>> pairs;
    std::vector<std::vector<std::size_t>> by_time;
    std::vector<std::uint32_t> full_counts;
};

Projection project(const TemporalBipartiteGraph& graph, const Catalog& catalog, AttributePair attributes) {
    Projection result;
    result.by_time.resize(graph.time_points().size());
    for (const auto& bucket : catalog.by_attributes) {
        if (!(attributes.u & (AttributeMask{1} << bucket.first.first)) ||
            !(attributes.v & (AttributeMask{1} << bucket.first.second))) continue;
        for (auto edge : bucket.second) {
            const auto local = result.pairs.size();
            result.pairs.emplace_back(graph.unique_edge_u(edge), graph.unique_edge_v(edge));
            const auto count = graph.unique_edge_timestamp_count(edge);
            result.full_counts.push_back(count);
            const auto offset = graph.unique_edge_timestamp_offset(edge);
            for (std::uint32_t i = 0; i < count; ++i)
                result.by_time[catalog.time_indices[offset + i]].push_back(local);
        }
    }
    return result;
}

BipartiteSubgraph snapshot(const TemporalBipartiteGraph& graph, AttributePair attributes,
                           const Projection& projection, const std::vector<std::uint32_t>& counts) {
    BipartiteSubgraph result(graph.num_u_vertices(), graph.num_v_vertices());
    for (VertexID u = 0; u < graph.num_u_vertices(); ++u)
        result.u_selected[u] = (attributes.u & (AttributeMask{1} << graph.u_attribute(u))) != 0;
    for (VertexID v = 0; v < graph.num_v_vertices(); ++v)
        result.v_selected[v] = (attributes.v & (AttributeMask{1} << graph.v_attribute(v))) != 0;
    for (std::size_t i = 0; i < counts.size(); ++i) if (counts[i]) {
        const auto [u, v] = projection.pairs[i];
        result.u_adj[u].push_back(v);
        result.v_adj[v].push_back(u);
    }
    return result;
}

std::size_t core_count(const BicoreSkyline& skyline) {
    std::vector<Count> maximum;
    for (const auto& row : skyline.u_beta_max) {
        if (maximum.size() < row.size()) maximum.resize(row.size());
        for (std::size_t a = 1; a < row.size(); ++a) maximum[a] = std::max(maximum[a], row[a]);
    }
    std::size_t result = 0;
    for (std::size_t a = 1; a < maximum.size(); ++a)
        result += maximum[a] - (a == 1 && maximum[a] ? 1 : 0);
    return result;
}

template<class Init, class Emit>
void sweep(const TemporalBipartiteGraph& graph, BuildStats& stats, const BuildProgress& progress,
           Init initialize, Emit emit) {
    const Catalog catalog(graph);
    const auto u_masks = subsets(catalog.all.u), v_masks = subsets(catalog.all.v);
    if (!v_masks.empty() && u_masks.size() > std::numeric_limits<std::size_t>::max() / v_masks.size())
        throw std::runtime_error("joint OPT attribute pair count overflow");
    stats.attribute_pairs = u_masks.size() * v_masks.size();
    const auto& times = graph.time_points();
    if (times.empty()) { stats.attribute_pairs = 0; return; }
    for (auto u : u_masks) for (auto v : v_masks) {
        const AttributePair attributes{u, v};
        const auto projection = project(graph, catalog, attributes);
        auto suffix_counts = projection.full_counts;
        std::size_t suffix_edges = projection.pairs.size();
        if (!suffix_edges) ++stats.empty_attribute_pairs;
        for (std::size_t s = 0; s < times.size() && suffix_edges; ++s) {
            auto counts = suffix_counts;
            DynamicBicoreMaintainer maintainer(snapshot(graph, attributes, projection, counts));
            ++stats.decompositions;
            if (s == 0) {
                const auto valid = core_count(maintainer.skyline());
                stats.joint_valid_cores += valid;
                if (u == catalog.all.u && v == catalog.all.v) stats.global_valid_cores = valid;
                if (!valid) { ++stats.empty_attribute_pairs; break; }
                initialize(attributes, maintainer.skyline());
            }
            if (!maintainer.membership_count()) break;
            for (std::size_t e = times.size(); e-- > s;) {
                if (!maintainer.membership_count()) break;
                ++stats.windows;
                std::vector<std::pair<VertexID, VertexID>> disappeared;
                for (auto id : projection.by_time[e])
                    if (counts[id] && --counts[id] == 0) disappeared.push_back(projection.pairs[id]);
                if (disappeared.empty()) continue;
                ++stats.dynamic_bucket_updates;
                stats.deleted_unique_edges += disappeared.size();
                const auto update = maintainer.delete_edge_bucket(disappeared);
                stats.frontier_changes += update.u_frontier_changes.size() + update.v_frontier_changes.size();
                for (const auto& change : update.u_frontier_changes)
                    for (std::uint64_t b = std::uint64_t(change.new_beta_max) + 1; b <= change.old_beta_max; ++b) {
                        if (change.alpha == 1 && b == 1) continue;
                        emit(attributes, Side::U, change.vertex, change.alpha, static_cast<Count>(b), times[s], times[e]);
                        ++stats.emitted_memberships;
                    }
                for (const auto& change : update.v_frontier_changes)
                    for (std::uint64_t a = std::uint64_t(change.new_alpha_max) + 1; a <= change.old_alpha_max; ++a) {
                        if (a == 1 && change.beta == 1) continue;
                        emit(attributes, Side::V, change.vertex, static_cast<Count>(a), change.beta, times[s], times[e]);
                        ++stats.emitted_memberships;
                    }
            }
            for (auto id : projection.by_time[s])
                if (suffix_counts[id] && --suffix_counts[id] == 0) --suffix_edges;
        }
        ++stats.completed_attribute_pairs;
        if (progress) progress(stats.completed_attribute_pairs, stats.attribute_pairs, attributes);
    }
}

void prepare_vertex(VATVertex& vertex, AttributePair attributes, const std::vector<Count>& beta_max) {
    if (beta_max.size() <= 1) return;
    std::size_t first = 1, last = beta_max.size();
    const auto valid = [&](std::size_t a) { return beta_max[a] >= (a == 1 ? 2u : 1u); };
    while (first < last && !valid(first)) ++first;
    while (last > first && !valid(last - 1)) --last;
    if (first == last) return;
    VATAttributes node{attributes, static_cast<Count>(first), {}};
    node.alpha_rows.resize(last - first);
    for (std::size_t a = first; a < last; ++a) if (valid(a)) {
        auto& row = node.alpha_rows[a - first];
        row.beta_begin = a == 1 ? 2 : 1;
        row.beta_lists.resize(std::size_t(beta_max[a]) - row.beta_begin + 1);
    }
    vertex.attributes.push_back(std::move(node));
}
}

VAT build_vat_opt(const TemporalBipartiteGraph& graph, BuildStats& stats, const BuildProgress& progress) {
    stats = {};
    VAT index;
    for (VertexID u = 0; u < graph.num_u_vertices(); ++u) index.u_vertices.push_back({graph.original_u_id(u), {}});
    for (VertexID v = 0; v < graph.num_v_vertices(); ++v) index.v_vertices.push_back({graph.original_v_id(v), {}});
    sweep(graph, stats, progress, [&](AttributePair attributes, const BicoreSkyline& skyline) {
        for (VertexID u = 0; u < skyline.u_beta_max.size(); ++u)
            prepare_vertex(index.u_vertices[u], attributes, skyline.u_beta_max[u]);
        for (VertexID v = 0; v < skyline.v_alpha_max.size(); ++v) {
            const auto& row = skyline.v_alpha_max[v];
            const auto max_alpha = row.empty() ? Count{0} : *std::max_element(row.begin(), row.end());
            std::vector<Count> beta_max(std::size_t(max_alpha) + 1);
            for (std::size_t b = 1; b < row.size(); ++b)
                for (std::size_t a = 1; a <= row[b]; ++a) beta_max[a] = static_cast<Count>(b);
            prepare_vertex(index.v_vertices[v], attributes, beta_max);
        }
    }, [&](AttributePair, Side side, VertexID id, Count a, Count b, Timestamp start, Timestamp end) {
        auto& vertex = side == Side::U ? index.u_vertices[id] : index.v_vertices[id];
        auto& node = vertex.attributes.back();
        auto& row = node.alpha_rows.at(a - node.alpha_begin);
        auto& intervals = row.beta_lists.at(b - row.beta_begin);
        if (intervals.empty()) ++stats.vat_core_records;
        if (!intervals.empty() && intervals.back().end == end) intervals.back().start = start;
        else { intervals.push_back({start, end}); ++stats.vat_intervals; }
    });
    build_vat_directory(index);
    return index;
}

UAT build_uat_opt(const TemporalBipartiteGraph& graph, BuildStats& stats, const BuildProgress& progress) {
    stats = {};
    UAT index;
    index.time_points = graph.time_points();
    std::map<CoreKey, std::vector<UATAttributes>> cores;
    sweep(graph, stats, progress, [](AttributePair, const BicoreSkyline&) {},
          [&](AttributePair attributes, Side side, VertexID id, Count a, Count b, Timestamp start, Timestamp end) {
        auto& entries = cores[{a, b}];
        if (entries.empty() || entries.back().attributes.u != attributes.u || entries.back().attributes.v != attributes.v)
            entries.push_back({attributes, {}});
        auto& starts = entries.back().starts;
        if (starts.empty() || starts.back().start != start) starts.push_back({start, {}});
        auto& ends = starts.back().ends;
        if (ends.empty() || ends.back().end != end) { ends.push_back({end, {}, {}}); ++stats.uat_blocks; }
        auto& vertices = side == Side::U ? ends.back().u_vertices : ends.back().v_vertices;
        vertices.push_back(side == Side::U ? graph.original_u_id(id) : graph.original_v_id(id));
        ++stats.vertex_refs;
    });
    for (auto& core : cores) {
        for (auto& attributes : core.second) for (auto& start : attributes.starts) {
            std::reverse(start.ends.begin(), start.ends.end());
            for (auto& end : start.ends) {
                std::sort(end.u_vertices.begin(), end.u_vertices.end());
                std::sort(end.v_vertices.begin(), end.v_vertices.end());
            }
        }
        index.cores.push_back({core.first.first, core.first.second, std::move(core.second)});
    }
    return index;
}
}
