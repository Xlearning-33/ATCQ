#include "core/bicore_decomposition.h"

#include <algorithm>
#include <functional>
#include <queue>
#include <utility>

namespace tacc {
namespace {

using DegreeEntry = std::pair<Count, VertexID>;

Count symmetric_degeneracy(const BipartiteSubgraph& graph) {
    const VertexID v_offset = static_cast<VertexID>(graph.u_adj.size());
    const std::size_t total = graph.u_adj.size() + graph.v_adj.size();
    std::vector<Count> degree(total, 0);
    std::vector<std::uint8_t> active(total, 0);
    std::priority_queue<DegreeEntry, std::vector<DegreeEntry>,
                        std::greater<DegreeEntry>> heap;
    for (VertexID u = 0; u < graph.u_adj.size(); ++u) {
        if (!graph.u_selected[u]) continue;
        degree[u] = static_cast<Count>(graph.u_adj[u].size());
        active[u] = 1;
        heap.emplace(degree[u], u);
    }
    for (VertexID v = 0; v < graph.v_adj.size(); ++v) {
        if (!graph.v_selected[v]) continue;
        const VertexID id = v_offset + v;
        degree[id] = static_cast<Count>(graph.v_adj[v].size());
        active[id] = 1;
        heap.emplace(degree[id], id);
    }

    Count delta = 0;
    while (!heap.empty()) {
        const auto [d, id] = heap.top();
        heap.pop();
        if (!active[id] || degree[id] != d) continue;
        active[id] = 0;
        delta = std::max(delta, d);
        if (id < v_offset) {
            for (VertexID v : graph.u_adj[id]) {
                const VertexID neighbor = v_offset + v;
                if (!active[neighbor]) continue;
                --degree[neighbor];
                heap.emplace(degree[neighbor], neighbor);
            }
        } else {
            const VertexID v = id - v_offset;
            for (VertexID u : graph.v_adj[v]) {
                if (!active[u]) continue;
                --degree[u];
                heap.emplace(degree[u], u);
            }
        }
    }
    return delta;
}

void compute_beta_max(const BipartiteSubgraph& graph, Count alpha,
                      BicoreSkyline& skyline) {
    std::vector<Count> u_degree(graph.u_adj.size(), 0);
    std::vector<Count> v_degree(graph.v_adj.size(), 0);
    std::vector<std::uint8_t> u_active = graph.u_selected;
    std::vector<std::uint8_t> v_active = graph.v_selected;
    std::queue<VertexID> remove_u;
    std::queue<VertexID> remove_v;
    std::priority_queue<DegreeEntry, std::vector<DegreeEntry>,
                        std::greater<DegreeEntry>> v_heap;

    for (VertexID u = 0; u < graph.u_adj.size(); ++u) {
        if (!u_active[u]) continue;
        u_degree[u] = static_cast<Count>(graph.u_adj[u].size());
        if (u_degree[u] < alpha) remove_u.push(u);
    }
    for (VertexID v = 0; v < graph.v_adj.size(); ++v) {
        if (!v_active[v]) continue;
        v_degree[v] = static_cast<Count>(graph.v_adj[v].size());
    }
    while (!remove_u.empty()) {
        const VertexID u = remove_u.front();
        remove_u.pop();
        if (!u_active[u] || u_degree[u] >= alpha) continue;
        u_active[u] = 0;
        for (VertexID v : graph.u_adj[u]) {
            if (v_active[v] && v_degree[v] != 0) --v_degree[v];
        }
    }
    for (VertexID v = 0; v < graph.v_adj.size(); ++v) {
        if (v_active[v]) v_heap.emplace(v_degree[v], v);
    }

    while (!v_heap.empty()) {
        while (!v_heap.empty() &&
               (!v_active[v_heap.top().second] ||
                v_degree[v_heap.top().second] != v_heap.top().first)) {
            v_heap.pop();
        }
        if (v_heap.empty()) break;
        const Count beta_level = v_heap.top().first;
        while (!v_heap.empty()) {
            const auto [degree, v] = v_heap.top();
            if (degree > beta_level) break;
            v_heap.pop();
            if (v_active[v] && v_degree[v] <= beta_level) remove_v.push(v);
        }

        while (!remove_v.empty() || !remove_u.empty()) {
            while (!remove_v.empty()) {
                const VertexID v = remove_v.front();
                remove_v.pop();
                if (!v_active[v] || v_degree[v] > beta_level) continue;
                v_active[v] = 0;
                const Count limit = std::min<Count>(
                    beta_level,
                    skyline.v_alpha_max[v].empty()
                        ? 0
                        : static_cast<Count>(
                              skyline.v_alpha_max[v].size() - 1));
                for (Count beta = 1; beta <= limit; ++beta) {
                    skyline.v_alpha_max[v][beta] = std::max(
                        skyline.v_alpha_max[v][beta], alpha);
                }
                for (VertexID u : graph.v_adj[v]) {
                    if (!u_active[u]) continue;
                    if (u_degree[u] != 0) --u_degree[u];
                    if (u_degree[u] < alpha) remove_u.push(u);
                }
            }
            while (!remove_u.empty()) {
                const VertexID u = remove_u.front();
                remove_u.pop();
                if (!u_active[u] || u_degree[u] >= alpha) continue;
                u_active[u] = 0;
                if (alpha < skyline.u_beta_max[u].size()) {
                    skyline.u_beta_max[u][alpha] = beta_level;
                }
                for (VertexID v : graph.u_adj[u]) {
                    if (!v_active[v]) continue;
                    if (v_degree[v] != 0) --v_degree[v];
                    v_heap.emplace(v_degree[v], v);
                    if (v_degree[v] <= beta_level) remove_v.push(v);
                }
            }
        }
    }
}

void compute_alpha_max(const BipartiteSubgraph& graph, Count beta,
                       BicoreSkyline& skyline) {
    std::vector<Count> u_degree(graph.u_adj.size(), 0);
    std::vector<Count> v_degree(graph.v_adj.size(), 0);
    std::vector<std::uint8_t> u_active = graph.u_selected;
    std::vector<std::uint8_t> v_active = graph.v_selected;
    std::queue<VertexID> remove_u;
    std::queue<VertexID> remove_v;
    std::priority_queue<DegreeEntry, std::vector<DegreeEntry>,
                        std::greater<DegreeEntry>> u_heap;

    for (VertexID u = 0; u < graph.u_adj.size(); ++u) {
        if (!u_active[u]) continue;
        u_degree[u] = static_cast<Count>(graph.u_adj[u].size());
    }
    for (VertexID v = 0; v < graph.v_adj.size(); ++v) {
        if (!v_active[v]) continue;
        v_degree[v] = static_cast<Count>(graph.v_adj[v].size());
        if (v_degree[v] < beta) remove_v.push(v);
    }
    while (!remove_v.empty()) {
        const VertexID v = remove_v.front();
        remove_v.pop();
        if (!v_active[v] || v_degree[v] >= beta) continue;
        v_active[v] = 0;
        for (VertexID u : graph.v_adj[v]) {
            if (u_active[u] && u_degree[u] != 0) --u_degree[u];
        }
    }
    for (VertexID u = 0; u < graph.u_adj.size(); ++u) {
        if (u_active[u]) u_heap.emplace(u_degree[u], u);
    }

    while (!u_heap.empty()) {
        while (!u_heap.empty() &&
               (!u_active[u_heap.top().second] ||
                u_degree[u_heap.top().second] != u_heap.top().first)) {
            u_heap.pop();
        }
        if (u_heap.empty()) break;
        const Count alpha_level = u_heap.top().first;
        while (!u_heap.empty()) {
            const auto [degree, u] = u_heap.top();
            if (degree > alpha_level) break;
            u_heap.pop();
            if (u_active[u] && u_degree[u] <= alpha_level) remove_u.push(u);
        }

        while (!remove_u.empty() || !remove_v.empty()) {
            while (!remove_u.empty()) {
                const VertexID u = remove_u.front();
                remove_u.pop();
                if (!u_active[u] || u_degree[u] > alpha_level) continue;
                u_active[u] = 0;
                const Count limit = std::min<Count>(
                    alpha_level,
                    skyline.u_beta_max[u].empty()
                        ? 0
                        : static_cast<Count>(
                              skyline.u_beta_max[u].size() - 1));
                for (Count alpha = 1; alpha <= limit; ++alpha) {
                    skyline.u_beta_max[u][alpha] = std::max(
                        skyline.u_beta_max[u][alpha], beta);
                }
                for (VertexID v : graph.u_adj[u]) {
                    if (!v_active[v]) continue;
                    if (v_degree[v] != 0) --v_degree[v];
                    if (v_degree[v] < beta) remove_v.push(v);
                }
            }
            while (!remove_v.empty()) {
                const VertexID v = remove_v.front();
                remove_v.pop();
                if (!v_active[v] || v_degree[v] >= beta) continue;
                v_active[v] = 0;
                if (beta < skyline.v_alpha_max[v].size()) {
                    skyline.v_alpha_max[v][beta] = alpha_level;
                }
                for (VertexID u : graph.v_adj[v]) {
                    if (!u_active[u]) continue;
                    if (u_degree[u] != 0) --u_degree[u];
                    u_heap.emplace(u_degree[u], u);
                    if (u_degree[u] <= alpha_level) remove_u.push(u);
                }
            }
        }
    }
}

}

QueryResult bicore_peeling(const BipartiteSubgraph& graph,
                           Count alpha,
                           Count beta) {
    const std::size_t num_u = graph.u_adj.size();
    const std::size_t num_v = graph.v_adj.size();
    std::vector<Count> u_degree(num_u, 0);
    std::vector<Count> v_degree(num_v, 0);
    std::vector<std::uint8_t> u_active = graph.u_selected;
    std::vector<std::uint8_t> v_active = graph.v_selected;
    std::queue<std::pair<bool, VertexID>> peel_queue;

    for (VertexID u = 0; u < num_u; ++u) {
        if (!u_active[u]) continue;
        u_degree[u] = static_cast<Count>(graph.u_adj[u].size());
        if (u_degree[u] < alpha) peel_queue.emplace(true, u);
    }
    for (VertexID v = 0; v < num_v; ++v) {
        if (!v_active[v]) continue;
        v_degree[v] = static_cast<Count>(graph.v_adj[v].size());
        if (v_degree[v] < beta) peel_queue.emplace(false, v);
    }

    while (!peel_queue.empty()) {
        const auto [is_u_side, id] = peel_queue.front();
        peel_queue.pop();
        if (is_u_side) {
            if (!u_active[id]) continue;
            u_active[id] = 0;
            for (VertexID v : graph.u_adj[id]) {
                if (!v_active[v]) continue;
                if (v_degree[v] > 0) --v_degree[v];
                if (v_degree[v] < beta) peel_queue.emplace(false, v);
            }
        } else {
            if (!v_active[id]) continue;
            v_active[id] = 0;
            for (VertexID u : graph.v_adj[id]) {
                if (!u_active[u]) continue;
                if (u_degree[u] > 0) --u_degree[u];
                if (u_degree[u] < alpha) peel_queue.emplace(true, u);
            }
        }
    }

    QueryResult result;
    for (VertexID u = 0; u < num_u; ++u) {
        if (u_active[u]) result.u_vertices.push_back(u);
    }
    for (VertexID v = 0; v < num_v; ++v) {
        if (v_active[v]) result.v_vertices.push_back(v);
    }
    return result;
}

QueryResult bicore_peeling_candidates(
    BipartiteSubgraph& graph,
    Count alpha,
    Count beta,
    const std::vector<VertexID>& candidate_u,
    const std::vector<VertexID>& candidate_v) {
    thread_local std::vector<Count> u_degree;
    thread_local std::vector<Count> v_degree;
    thread_local std::vector<std::pair<bool, VertexID>> peel_queue;
    if (u_degree.size() < graph.u_adj.size()) {
        u_degree.resize(graph.u_adj.size());
    }
    if (v_degree.size() < graph.v_adj.size()) {
        v_degree.resize(graph.v_adj.size());
    }
    peel_queue.clear();
    if (peel_queue.capacity() <
        candidate_u.size() + candidate_v.size()) {
        peel_queue.reserve(candidate_u.size() + candidate_v.size());
    }

    for (VertexID u : candidate_u) {
        if (!graph.u_selected[u]) continue;
        u_degree[u] = static_cast<Count>(graph.u_adj[u].size());
        if (u_degree[u] < alpha) peel_queue.emplace_back(true, u);
    }
    for (VertexID v : candidate_v) {
        if (!graph.v_selected[v]) continue;
        v_degree[v] = static_cast<Count>(graph.v_adj[v].size());
        if (v_degree[v] < beta) peel_queue.emplace_back(false, v);
    }

    std::size_t head = 0;
    while (head < peel_queue.size()) {
        const auto [is_u_side, id] = peel_queue[head++];
        if (is_u_side) {
            if (!graph.u_selected[id]) continue;
            graph.u_selected[id] = 0;
            for (VertexID v : graph.u_adj[id]) {
                if (!graph.v_selected[v]) continue;
                const Count previous = v_degree[v];
                if (previous != 0) --v_degree[v];
                if (previous >= beta && v_degree[v] < beta) {
                    peel_queue.emplace_back(false, v);
                }
            }
        } else {
            if (!graph.v_selected[id]) continue;
            graph.v_selected[id] = 0;
            for (VertexID u : graph.v_adj[id]) {
                if (!graph.u_selected[u]) continue;
                const Count previous = u_degree[u];
                if (previous != 0) --u_degree[u];
                if (previous >= alpha && u_degree[u] < alpha) {
                    peel_queue.emplace_back(true, u);
                }
            }
        }
    }

    QueryResult result;
    result.u_vertices.reserve(candidate_u.size());
    result.v_vertices.reserve(candidate_v.size());
    for (VertexID u : candidate_u) {
        if (graph.u_selected[u]) result.u_vertices.push_back(u);
    }
    for (VertexID v : candidate_v) {
        if (graph.v_selected[v]) result.v_vertices.push_back(v);
    }
    return result;
}

BicoreSkyline bicore_skyline_decomposition(
    const BipartiteSubgraph& graph) {
    BicoreSkyline skyline;
    for (VertexID u = 0; u < graph.u_adj.size(); ++u) {
        if (graph.u_selected[u]) {
            skyline.max_alpha = std::max(
                skyline.max_alpha,
                static_cast<Count>(graph.u_adj[u].size()));
        }
    }
    for (VertexID v = 0; v < graph.v_adj.size(); ++v) {
        if (graph.v_selected[v]) {
            skyline.max_beta = std::max(
                skyline.max_beta,
                static_cast<Count>(graph.v_adj[v].size()));
        }
    }
    skyline.u_beta_max.resize(graph.u_adj.size());
    skyline.v_alpha_max.resize(graph.v_adj.size());
    for (VertexID u = 0; u < graph.u_adj.size(); ++u) {
        if (graph.u_selected[u]) {
            skyline.u_beta_max[u].assign(graph.u_adj[u].size() + 1, 0);
        }
    }
    for (VertexID v = 0; v < graph.v_adj.size(); ++v) {
        if (graph.v_selected[v]) {
            skyline.v_alpha_max[v].assign(graph.v_adj[v].size() + 1, 0);
        }
    }

    const Count delta = symmetric_degeneracy(graph);
    for (Count alpha = 1; alpha <= delta; ++alpha) {
        compute_beta_max(graph, alpha, skyline);
    }
    for (Count beta = 1; beta <= delta; ++beta) {
        compute_alpha_max(graph, beta, skyline);
    }
    return skyline;
}

}
