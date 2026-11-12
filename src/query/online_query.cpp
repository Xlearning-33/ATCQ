#include "query/online_query.h"

#include "core/bicore_decomposition.h"

#include <algorithm>
#include <chrono>

namespace tacc {

QueryResult online_query(const TemporalBipartiteGraph& graph,
                         const Query& query,
                         QueryProfile* profile,
                         QueryTimer* timer) {
    if (profile != nullptr) *profile = {};
    if (query.window.start > query.window.end ||
        query.alpha == 0 || query.beta == 0 ||
        query.u_attributes.empty() || query.v_attributes.empty()) {
        stop_query_timer(timer);
        return {};
    }

    std::chrono::steady_clock::time_point recovery_start;
    if (profile != nullptr) recovery_start = std::chrono::steady_clock::now();
    const auto& events = graph.edges();
    const auto first = std::lower_bound(
        events.begin(), events.end(), query.window.start,
        [](const TemporalEdge& edge, Timestamp t) { return edge.timestamp < t; });
    const auto last = std::upper_bound(
        first, events.end(), query.window.end,
        [](Timestamp t, const TemporalEdge& edge) { return t < edge.timestamp; });

    BipartiteSubgraph subgraph(graph.num_u_vertices(), graph.num_v_vertices());
    for (auto it = first; it != last; ++it) {
        if (std::find(query.u_attributes.begin(), query.u_attributes.end(),
                      graph.u_attribute(it->u)) == query.u_attributes.end() ||
            std::find(query.v_attributes.begin(), query.v_attributes.end(),
                      graph.v_attribute(it->v)) == query.v_attributes.end()) {
            continue;
        }
        subgraph.u_selected[it->u] = 1;
        subgraph.v_selected[it->v] = 1;
        subgraph.u_adj[it->u].push_back(it->v);
    }
    std::size_t window_edges = 0;
    for (VertexID u = 0; u < subgraph.u_adj.size(); ++u) {
        auto& neighbors = subgraph.u_adj[u];
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
        window_edges += neighbors.size();
        for (VertexID v : neighbors) subgraph.v_adj[v].push_back(u);
    }
    if (profile != nullptr) {
        profile->window_recovery_us =
            std::chrono::duration<double, std::micro>(
                std::chrono::steady_clock::now() - recovery_start).count();
        profile->window_edges = window_edges;
        profile->candidate_u = std::count(subgraph.u_selected.begin(), subgraph.u_selected.end(), 1);
        profile->candidate_v = std::count(subgraph.v_selected.begin(), subgraph.v_selected.end(), 1);
    }

    std::chrono::steady_clock::time_point peeling_start;
    if (profile != nullptr) peeling_start = std::chrono::steady_clock::now();
    QueryResult internal = bicore_peeling(subgraph, query.alpha, query.beta);
    if (profile != nullptr) {
        profile->peeling_us =
            std::chrono::duration<double, std::micro>(
                std::chrono::steady_clock::now() - peeling_start).count();
    }

    QueryResult result;
    result.u_vertices.reserve(internal.u_vertices.size());
    result.v_vertices.reserve(internal.v_vertices.size());
    for (VertexID u : internal.u_vertices) {
        result.u_vertices.push_back(graph.original_u_id(u));
    }
    for (VertexID v : internal.v_vertices) {
        result.v_vertices.push_back(graph.original_v_id(v));
    }
    if (profile != nullptr) {
        profile->result_u = result.u_vertices.size();
        profile->result_v = result.v_vertices.size();
    }
    stop_query_timer(timer);
    return result;
}

}
