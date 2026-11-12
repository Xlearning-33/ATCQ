#include "query/mat_recovery.h"

#include "core/bicore_decomposition.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace tacc {

MATRecovery configured_mat_recovery() {
    static const MATRecovery mode = [] {
        const char* value = std::getenv("TAC_MAT_RECOVERY");
        if (value == nullptr || std::string(value) == "grouped") return MATRecovery::Grouped;
        if (std::string(value) == "legacy") return MATRecovery::Legacy;
        throw std::invalid_argument("Invalid TAC_MAT_RECOVERY; expected grouped|legacy");
    }();
    return mode;
}

const char* mat_recovery_name(MATRecovery mode) {
    return mode == MATRecovery::Grouped ? "grouped" : "legacy";
}

namespace {
using Clock = std::chrono::steady_clock;
constexpr std::size_t retained_adjacency_cap = 64 * 1024 * 1024;
constexpr std::size_t repeated_event_limit = 8 * 1024 * 1024;

template <class T> std::size_t bytes(const std::vector<T>& v) {
    return v.capacity() * sizeof(T);
}

struct Workspace {
    BipartiteSubgraph snapshot;
    std::vector<VertexID> u, v;
    std::vector<std::size_t> last_u_bytes, last_v_bytes;
    std::size_t retained{};

    Workspace(std::size_t nu = 0, std::size_t nv = 0)
        : snapshot(nu, nv), last_u_bytes(nu), last_v_bytes(nv) {}

    void clear() noexcept {
        const auto reset = [&](auto& adj, auto& marks, auto& previous, const auto& ids) {
            for (VertexID id : ids) {
                adj[id].clear();
                marks[id] = 0;
                retained -= previous[id];
                const auto now = bytes(adj[id]);
                if (now > retained_adjacency_cap - retained) {
                    std::vector<VertexID>().swap(adj[id]);
                    previous[id] = 0;
                } else {
                    retained += now;
                    previous[id] = now;
                }
            }
        };
        reset(snapshot.u_adj, snapshot.u_selected, last_u_bytes, u);
        reset(snapshot.v_adj, snapshot.v_selected, last_v_bytes, v);
        u.clear();
        v.clear();
    }

    std::size_t fixed_bytes() const {
        return bytes(snapshot.u_adj) + bytes(snapshot.v_adj) +
            bytes(snapshot.u_selected) + bytes(snapshot.v_selected) +
            bytes(u) + bytes(v) +
            bytes(last_u_bytes) + bytes(last_v_bytes);
    }
};

Workspace& workspace(const TemporalBipartiteGraph& graph) {
    thread_local Workspace work;
    if (work.snapshot.u_adj.size() != graph.num_u_vertices() ||
        work.snapshot.v_adj.size() != graph.num_v_vertices())
        work = Workspace(graph.num_u_vertices(), graph.num_v_vertices());
    return work;
}

struct ClearOnExit {
    Workspace& work;
    QueryProfile* profile;
    ~ClearOnExit() {
        const auto start = profile ? Clock::now() : Clock::time_point{};
        work.clear();
        if (profile) {
            profile->workspace_clear_us += std::chrono::duration<double, std::micro>(Clock::now() - start).count();
            profile->workspace_fixed_bytes = work.fixed_bytes();
            profile->workspace_retained_bytes = work.retained;
        }
    }
};
}

QueryResult query_mat_grouped_seeds(const TemporalBipartiteGraph& graph,
                                   const Query& query,
                                   std::vector<VertexID> seeds,
                                   AttributeMask u_mask,
                                   QueryProfile* profile,
                                   QueryTimer* timer) {
    if (seeds.empty()) { stop_query_timer(timer); return {}; }
    const auto start = profile ? Clock::now() : Clock::time_point{};
    const bool cached = profile && graph.grouped_v_time_aux_bytes() != 0;
    const auto& aux = graph.grouped_v_time_adjacency();
    if (profile) {
        profile->grouped_aux_bytes = aux.bytes();
        if (!cached) profile->grouped_aux_build_us =
            std::chrono::duration<double, std::micro>(Clock::now() - start).count();
    }
    auto& work = workspace(graph);
    ClearOnExit clear{work, profile};
    work.v = std::move(seeds);
    auto& snapshot = work.snapshot;
    for (VertexID v : work.v) snapshot.v_selected[v] = 1;

    const bool unique = graph.num_edges() == graph.num_unique_edges();
    bool fallback = false;
    if (!unique && graph.num_edges() > repeated_event_limit) {
        std::size_t events = 0;
        for (VertexID v : work.v) {
            const auto range = aux.range(v, query.window);
            const auto count = range.second - range.first;
            if (count > repeated_event_limit - events) { fallback = true; break; }
            events += count;
        }
    }
    std::size_t edges = 0, accepted = 0;
    if (fallback) {
        if (profile) ++profile->grouped_budget_fallbacks;
        WindowAdjacencyStats stats;
        try {
            stats = graph.populate_seed_window_adjacency(snapshot, query.window,
                work.v, false, u_mask, work.u, true);
        } catch (...) {
            work = Workspace{};
            throw;
        }
        edges = stats.window_edges;
        if (profile) {
            profile->unique_edges_examined = stats.unique_edges_examined;
            profile->adaptive_seed_recoveries = 1;
            profile->reversed_seed_scans = stats.reversed_seed_scan ? 1 : 0;
            profile->seed_scan_work = stats.seed_scan_work;
            profile->opposite_scan_work = stats.opposite_scan_work;
            profile->incident_work_u = stats.incident_work_u;
            profile->incident_work_v = stats.incident_work_v;
        }
    } else {
        if (profile) ++profile->grouped_recoveries;
        for (VertexID v : work.v) {
            const auto range = aux.range(v, query.window);
            if (profile) profile->grouped_window_events += range.second - range.first;
            for (auto i = range.first; i < range.second; ++i) {
                const VertexID u = aux.neighbors[i];
                const auto attr = graph.u_attribute(u);
                if (attr >= 64 || (u_mask & (AttributeMask{1} << attr)) == 0) continue;
                ++accepted;
                if (!snapshot.u_selected[u]) {
                    work.u.push_back(u);
                    snapshot.u_selected[u] = 1;
                }
                snapshot.u_adj[u].push_back(v);
                if (unique) { snapshot.v_adj[v].push_back(u); ++edges; }
            }
        }
        if (!unique) {
            for (VertexID u : work.u) {
                auto& row = snapshot.u_adj[u];
                std::sort(row.begin(), row.end());
                row.erase(std::unique(row.begin(), row.end()), row.end());
                for (VertexID v : row) { snapshot.v_adj[v].push_back(u); ++edges; }
            }
            if (profile) profile->grouped_duplicate_events = accepted - edges;
        }
        std::sort(work.u.begin(), work.u.end());
    }
    if (profile) {
        profile->window_recovery_us += std::chrono::duration<double, std::micro>(Clock::now() - start).count();
        profile->candidate_u = work.u.size();
        profile->candidate_v = work.v.size();
        profile->window_edges = edges;
    }
    if (work.u.empty()) { stop_query_timer(timer); return {}; }
    const auto peel_start = profile ? Clock::now() : Clock::time_point{};
    const auto internal = bicore_peeling_candidates(snapshot, query.alpha, query.beta, work.u, work.v);
    if (profile) profile->peeling_us +=
        std::chrono::duration<double, std::micro>(Clock::now() - peel_start).count();
    QueryResult result;
    result.u_vertices.reserve(internal.u_vertices.size());
    result.v_vertices.reserve(internal.v_vertices.size());
    for (VertexID u : internal.u_vertices) result.u_vertices.push_back(graph.original_u_id(u));
    for (VertexID v : internal.v_vertices) result.v_vertices.push_back(graph.original_v_id(v));
    if (profile) {
        profile->result_u = result.u_vertices.size();
        profile->result_v = result.v_vertices.size();
    }
    stop_query_timer(timer);
    return result;
}

}
