#include "index/joint_index.h"
#include "query/online_query.h"

#include <algorithm>
#include <limits>
#include <iterator>
#include <stdexcept>
#include <unordered_map>

namespace tacc::joint {
namespace {
using CoreKey = std::pair<Count, Count>;

std::vector<AttributeID> attributes(AttributeMask mask) {
    std::vector<AttributeID> result;
    for (AttributeID a = 0; a < 64; ++a)
        if (mask & (AttributeMask{1} << a)) result.push_back(a);
    return result;
}

std::vector<AttributeMask> subsets(AttributeMask mask) {
    std::vector<AttributeMask> result;
    for (auto sub = mask; sub; sub = (sub - 1) & mask) result.push_back(sub);
    std::sort(result.begin(), result.end());
    return result;
}

QueryResult compute(const TemporalBipartiteGraph& graph, AttributePair pair,
                    CoreKey core, TimeWindow window, BuildStats& stats, bool precheck = false) {
    auto answer = online_query(graph, {window, attributes(pair.u), attributes(pair.v), core.first, core.second}, nullptr);
    ++stats.online_calls;
    if (precheck) ++stats.full_window_prechecks;
    std::sort(answer.u_vertices.begin(), answer.u_vertices.end());
    std::sort(answer.v_vertices.begin(), answer.v_vertices.end());
    return answer;
}

struct Plan {
    AttributePair all;
    std::vector<AttributeMask> u_masks, v_masks;
    std::vector<CoreKey> cores;
};

Plan plan(const TemporalBipartiteGraph& graph, BuildStats& stats) {
    Plan result;
    std::size_t max_a = 0, max_b = 0;
    const auto add = [](AttributeMask& mask, AttributeID id) {
        if (id >= 64) throw std::runtime_error("joint BASIC supports global attribute IDs 0--63");
        mask |= AttributeMask{1} << id;
    };
    for (VertexID u = 0; u < graph.num_u_vertices(); ++u) {
        add(result.all.u, graph.u_attribute(u));
        max_a = std::max(max_a, graph.static_u_degree(u));
    }
    for (VertexID v = 0; v < graph.num_v_vertices(); ++v) {
        add(result.all.v, graph.v_attribute(v));
        max_b = std::max(max_b, graph.static_v_degree(v));
    }
    if (graph.time_points().empty()) return result;
    if (max_a > std::numeric_limits<Count>::max() || max_b > std::numeric_limits<Count>::max())
        throw std::runtime_error("joint BASIC degree bound overflow");
    const TimeWindow full{graph.time_points().front(), graph.time_points().back()};
    for (std::size_t a = 1; a <= max_a; ++a) {
        for (std::size_t b = 1; b <= max_b; ++b) {
            if (a == 1 && b == 1) continue;
            const CoreKey core{static_cast<Count>(a), static_cast<Count>(b)};
            if (!compute(graph, result.all, core, full, stats, true).u_vertices.empty())
                result.cores.push_back(core);
        }
    }
    stats.global_valid_cores = result.cores.size();
    result.u_masks = subsets(result.all.u);
    result.v_masks = subsets(result.all.v);
    if (!result.v_masks.empty() && result.u_masks.size() > std::numeric_limits<std::size_t>::max() / result.v_masks.size())
        throw std::runtime_error("joint attribute pair count overflow");
    stats.attribute_pairs = result.u_masks.size() * result.v_masks.size();
    return result;
}

template<class Consume>
void windows(const TemporalBipartiteGraph& graph, AttributePair pair, CoreKey core,
             QueryResult full, BuildStats& stats, Consume consume) {
    const auto& times = graph.time_points();
    for (std::size_t s = 0; s < times.size(); ++s) {
        for (std::size_t e = s; e < times.size(); ++e) {
            auto answer = s == 0 && e + 1 == times.size() ? std::move(full)
                : compute(graph, pair, core, {times[s], times[e]}, stats);
            consume(s, e, std::move(answer));
        }
    }
}

void completed_pair(BuildStats& stats, const BuildProgress& progress, AttributePair pair) {
    ++stats.completed_attribute_pairs;
    if (progress) progress(stats.completed_attribute_pairs, stats.attribute_pairs, pair);
}

void append(VATVertex& vertex, AttributePair pair, CoreKey core, std::vector<Interval> intervals) {
    if (vertex.attributes.empty() || vertex.attributes.back().attributes.u != pair.u ||
        vertex.attributes.back().attributes.v != pair.v)
        vertex.attributes.push_back({pair, core.first, {}});
    auto& attr = vertex.attributes.back();
    attr.alpha_rows.resize(std::size_t(core.first) - attr.alpha_begin + 1);
    auto& row = attr.alpha_rows[core.first - attr.alpha_begin];
    if (row.beta_lists.empty()) row.beta_begin = core.second;
    row.beta_lists.resize(std::size_t(core.second) - row.beta_begin + 1);
    row.beta_lists[core.second - row.beta_begin] = std::move(intervals);
}

std::vector<VertexID> difference(const std::vector<VertexID>& current, const std::vector<VertexID>& previous) {
    std::vector<VertexID> result;
    std::set_difference(current.begin(), current.end(), previous.begin(), previous.end(), std::back_inserter(result));
    return result;
}
}

VAT build_vat_basic(const TemporalBipartiteGraph& graph, BuildStats& stats, const BuildProgress& progress) {
    stats = {};
    VAT index;
    std::unordered_map<VertexID, std::size_t> u_slot, v_slot;
    for (VertexID u = 0; u < graph.num_u_vertices(); ++u) {
        u_slot.emplace(graph.original_u_id(u), index.u_vertices.size());
        index.u_vertices.push_back({graph.original_u_id(u), {}});
    }
    const auto nu = index.u_vertices.size();
    for (VertexID v = 0; v < graph.num_v_vertices(); ++v) {
        v_slot.emplace(graph.original_v_id(v), nu + index.v_vertices.size());
        index.v_vertices.push_back({graph.original_v_id(v), {}});
    }
    const auto build_plan = plan(graph, stats);
    const auto& times = graph.time_points();
    if (times.empty()) { build_vat_directory(index); return index; }
    for (const auto u : build_plan.u_masks) for (const auto v : build_plan.v_masks) {
        const AttributePair pair{u, v};
        for (const auto& core : build_plan.cores) {
            auto full = compute(graph, pair, core, {times.front(), times.back()}, stats, true);
            if (full.u_vertices.empty()) continue;
            ++stats.joint_valid_cores;
            std::vector<std::vector<Interval>> lists(nu + index.v_vertices.size());
            std::vector<std::size_t> seen(lists.size(), std::numeric_limits<std::size_t>::max());
            windows(graph, pair, core, std::move(full), stats, [&](std::size_t s, std::size_t e, QueryResult answer) {
                const auto record = [&](std::size_t slot) {
                    if (seen[slot] == s) return;
                    seen[slot] = s;
                    auto& intervals = lists[slot];
                    if (!intervals.empty() && intervals.back().end == times[e]) intervals.back().start = times[s];
                    else intervals.push_back({times[s], times[e]});
                };
                for (auto id : answer.u_vertices) record(u_slot.at(id));
                for (auto id : answer.v_vertices) record(v_slot.at(id));
            });
            for (std::size_t slot = 0; slot < lists.size(); ++slot) {
                if (lists[slot].empty()) continue;
                ++stats.vat_core_records;
                stats.vat_intervals += lists[slot].size();
                append(slot < nu ? index.u_vertices[slot] : index.v_vertices[slot - nu], pair, core, std::move(lists[slot]));
            }
        }
        completed_pair(stats, progress, pair);
    }
    build_vat_directory(index);
    return index;
}

UAT build_uat_basic(const TemporalBipartiteGraph& graph, BuildStats& stats, const BuildProgress& progress) {
    stats = {};
    UAT index;
    index.time_points = graph.time_points();
    const auto build_plan = plan(graph, stats);
    const auto& times = index.time_points;
    if (times.empty()) return index;
    for (const auto& core : build_plan.cores) index.cores.push_back({core.first, core.second, {}});
    for (const auto u : build_plan.u_masks) for (const auto v : build_plan.v_masks) {
        const AttributePair pair{u, v};
        for (auto& core : index.cores) {
            const CoreKey key{core.alpha, core.beta};
            auto full = compute(graph, pair, key, {times.front(), times.back()}, stats, true);
            if (full.u_vertices.empty()) continue;
            ++stats.joint_valid_cores;
            UATAttributes attr{pair, {}};
            QueryResult previous;
            windows(graph, pair, key, std::move(full), stats, [&](std::size_t s, std::size_t e, QueryResult answer) {
                if (s == e) previous = {};
                UATEnd end{times[e], difference(answer.u_vertices, previous.u_vertices),
                                      difference(answer.v_vertices, previous.v_vertices)};
                if (!end.u_vertices.empty() || !end.v_vertices.empty()) {
                    if (attr.starts.empty() || attr.starts.back().start != times[s]) attr.starts.push_back({times[s], {}});
                    ++stats.uat_blocks;
                    stats.vertex_refs += end.u_vertices.size() + end.v_vertices.size();
                    attr.starts.back().ends.push_back(std::move(end));
                }
                previous = std::move(answer);
            });
            if (!attr.starts.empty()) core.attributes.push_back(std::move(attr));
        }
        completed_pair(stats, progress, pair);
    }
    index.cores.erase(std::remove_if(index.cores.begin(), index.cores.end(),
        [](const auto& core) { return core.attributes.empty(); }), index.cores.end());
    return index;
}
}
