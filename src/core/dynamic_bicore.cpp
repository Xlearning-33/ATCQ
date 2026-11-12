#include "core/dynamic_bicore.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <functional>
#include <limits>
#include <queue>
#include <string>
#include <tuple>
#include <utility>

namespace tacc {

namespace {
using Clock = std::chrono::steady_clock;
using Seconds = std::chrono::duration<double>;
}

DynamicBicoreMaintainer::DynamicBicoreMaintainer(
    BipartiteSubgraph graph)
    : skyline_(bicore_skyline_decomposition(graph)),
      u_neighbors_(std::move(graph.u_adj)),
      v_neighbors_(std::move(graph.v_adj)),
      u_edge_active_(u_neighbors_.size()),
      v_edge_active_(v_neighbors_.size()) {
    const char* layered = std::getenv("TAC_DYNAMIC_LAYERED");
    layered_enabled_ = layered == nullptr || std::string(layered) != "0";
    removed_u_epoch_.assign(u_neighbors_.size(), 0);
    removed_v_epoch_.assign(v_neighbors_.size(), 0);
    queued_u_epoch_.assign(u_neighbors_.size(), 0);
    queued_v_epoch_.assign(v_neighbors_.size(), 0);
    support_u_epoch_.assign(u_neighbors_.size(), 0);
    support_v_epoch_.assign(v_neighbors_.size(), 0);
    support_u_value_.assign(u_neighbors_.size(), 0);
    support_v_value_.assign(v_neighbors_.size(), 0);
    endpoint_u_beta_max_.assign(u_neighbors_.size(), 0);
    endpoint_v_beta_max_.assign(v_neighbors_.size(), 0);
    live_u_degree_.resize(u_neighbors_.size());
    live_v_degree_.resize(v_neighbors_.size());
    axis_u_epoch_.assign(u_neighbors_.size(), 0);
    axis_v_epoch_.assign(v_neighbors_.size(), 0);
    u_frontier_offsets_.reserve(skyline_.u_beta_max.size() + 1);
    u_frontier_offsets_.push_back(0);
    for (const auto& row : skyline_.u_beta_max) {
        u_frontier_offsets_.push_back(u_frontier_offsets_.back() + row.size());
    }
    v_frontier_offsets_.reserve(skyline_.v_alpha_max.size() + 1);
    v_frontier_offsets_.push_back(0);
    for (const auto& row : skyline_.v_alpha_max) {
        v_frontier_offsets_.push_back(v_frontier_offsets_.back() + row.size());
    }
    pending_u_epoch_.assign(u_frontier_offsets_.back(), 0);
    pending_v_epoch_.assign(v_frontier_offsets_.back(), 0);
    pending_u_new_beta_.resize(u_frontier_offsets_.back());
    pending_v_new_alpha_.resize(v_frontier_offsets_.back());
    if (layered_enabled_) {
        successor_u_.resize(u_frontier_offsets_.back());
        successor_v_.resize(v_frontier_offsets_.back());
    }
    for (VertexID u = 0; u < u_neighbors_.size(); ++u) {
        auto& row = u_neighbors_[u];
        std::sort(row.begin(), row.end());
        row.erase(std::unique(row.begin(), row.end()), row.end());
        u_edge_active_[u].assign(row.size(), 1);
        live_u_degree_[u] = static_cast<Count>(row.size());
    }
    for (VertexID v = 0; v < v_neighbors_.size(); ++v) {
        auto& row = v_neighbors_[v];
        std::sort(row.begin(), row.end());
        row.erase(std::unique(row.begin(), row.end()), row.end());
        v_edge_active_[v].assign(row.size(), 1);
        live_v_degree_[v] = static_cast<Count>(row.size());
    }
    for (const auto& row : skyline_.u_beta_max) {
        for (std::size_t alpha = 1; alpha < row.size(); ++alpha) {
            if (row[alpha] != 0) ++frontier_entries_;
            membership_entries_ += row[alpha];
            if (alpha == 1 && row[alpha] != 0) --membership_entries_;
            if (row[alpha] >= alpha)
                delta_bound_ = std::max(delta_bound_, static_cast<Count>(alpha));
        }
    }
    for (const auto& row : skyline_.v_alpha_max) {
        for (std::size_t beta = 1; beta < row.size(); ++beta) {
            if (row[beta] != 0) ++frontier_entries_;
            membership_entries_ += row[beta];
            if (beta == 1 && row[beta] != 0) --membership_entries_;
        }
    }
}

bool DynamicBicoreMaintainer::contains_edge(VertexID u, VertexID v) const {
    if (u >= u_neighbors_.size()) return false;
    const auto& row = u_neighbors_[u];
    const auto it = std::lower_bound(row.begin(), row.end(), v);
    return it != row.end() && *it == v &&
           u_edge_active_[u][static_cast<std::size_t>(it - row.begin())];
}

bool DynamicBicoreMaintainer::deactivate_edge(VertexID u, VertexID v) {
    if (u >= u_neighbors_.size() || v >= v_neighbors_.size()) return false;
    const auto u_it = std::lower_bound(
        u_neighbors_[u].begin(), u_neighbors_[u].end(), v);
    if (u_it == u_neighbors_[u].end() || *u_it != v) return false;
    const std::size_t u_index =
        static_cast<std::size_t>(u_it - u_neighbors_[u].begin());
    if (!u_edge_active_[u][u_index]) return false;
    const auto v_it = std::lower_bound(
        v_neighbors_[v].begin(), v_neighbors_[v].end(), u);
    if (v_it == v_neighbors_[v].end() || *v_it != u) return false;
    u_edge_active_[u][u_index] = 0;
    v_edge_active_[v][static_cast<std::size_t>(
        v_it - v_neighbors_[v].begin())] = 0;
    --live_u_degree_[u];
    --live_v_degree_[v];
    return true;
}

void DynamicBicoreMaintainer::record_frontier(
    bool is_u, VertexID id, Count parameter, Count value,
    std::vector<std::pair<VertexID, Count>>& touched_u,
    std::vector<std::pair<VertexID, Count>>& touched_v) {
    const auto& row = is_u ? skyline_.u_beta_max[id] : skyline_.v_alpha_max[id];
    if (parameter >= row.size() || value >= row[parameter]) return;
    const auto slot = (is_u ? u_frontier_offsets_[id] : v_frontier_offsets_[id]) +
        parameter;
    auto& epochs = is_u ? pending_u_epoch_ : pending_v_epoch_;
    auto& values = is_u ? pending_u_new_beta_ : pending_v_new_alpha_;
    auto& touched = is_u ? touched_u : touched_v;
    if (epochs[slot] != bucket_epoch_) {
        epochs[slot] = bucket_epoch_;
        values[slot] = value;
        if (layered_enabled_)
            (is_u ? successor_u_ : successor_v_)[slot] = parameter + 1;
        touched.emplace_back(id, parameter);
    } else {
        values[slot] = std::min(values[slot], value);
    }
}

std::uint32_t DynamicBicoreMaintainer::next_plane_epoch() {
    if (plane_epoch_ == std::numeric_limits<std::uint32_t>::max()) {
        for (auto* marks : {&removed_u_epoch_, &removed_v_epoch_,
                            &queued_u_epoch_, &queued_v_epoch_,
                            &support_u_epoch_, &support_v_epoch_}) {
            std::fill(marks->begin(), marks->end(), 0);
        }
        plane_epoch_ = 1;
    } else {
        ++plane_epoch_;
    }
    return plane_epoch_;
}

Count DynamicBicoreMaintainer::next_unwritten(bool is_u, VertexID id,
                                            Count parameter) {
    const auto& row = is_u ? skyline_.u_beta_max[id] : skyline_.v_alpha_max[id];
    const auto base = is_u ? u_frontier_offsets_[id] : v_frontier_offsets_[id];
    auto& next = is_u ? successor_u_ : successor_v_;
    const auto& epochs = is_u ? pending_u_epoch_ : pending_v_epoch_;
    Count root = parameter;
    while (root < row.size() && epochs[base + root] == bucket_epoch_)
        root = next[base + root];
    while (parameter != root) {
        const Count old_next = next[base + parameter];
        next[base + parameter] = root;
        parameter = old_next;
    }
    return root;
}

void DynamicBicoreMaintainer::maintain_layered(
    const std::vector<std::pair<VertexID, VertexID>>& edges,
    std::vector<std::pair<VertexID, Count>>& touched_u,
    std::vector<std::pair<VertexID, Count>>& touched_v,
    DynamicBicoreProfile* profile) {
    const auto start = profile != nullptr ? Clock::now() : Clock::time_point{};
    std::vector<std::pair<bool, VertexID>> region, queue;
    using HeapEntry = std::pair<Count, VertexID>;
    std::vector<HeapEntry> heap;
    const auto heap_compare = std::greater<HeapEntry>{};
    const auto sweep = [&](bool transposed, Count fixed, Count minimum) {
        const auto epoch = next_plane_epoch();
        const auto& fa = transposed ? v_neighbors_ : u_neighbors_;
        const auto& va = transposed ? u_neighbors_ : v_neighbors_;
        const auto& fe = transposed ? v_edge_active_ : u_edge_active_;
        const auto& ve = transposed ? u_edge_active_ : v_edge_active_;
        const auto& fr = transposed ? skyline_.v_alpha_max : skyline_.u_beta_max;
        const auto& vr = transposed ? skyline_.u_beta_max : skyline_.v_alpha_max;
        auto& fc_epoch = transposed ? queued_v_epoch_ : queued_u_epoch_;
        auto& vc_epoch = transposed ? queued_u_epoch_ : queued_v_epoch_;
        auto& fc = transposed ? endpoint_v_beta_max_ : endpoint_u_beta_max_;
        auto& vc = transposed ? endpoint_u_beta_max_ : endpoint_v_beta_max_;
        auto& region_f = transposed ? support_v_epoch_ : support_u_epoch_;
        auto& region_v = transposed ? support_u_epoch_ : support_v_epoch_;
        auto& fd = transposed ? support_v_value_ : support_u_value_;
        auto& vd = transposed ? support_u_value_ : support_v_value_;
        auto& removed_f = transposed ? removed_v_epoch_ : removed_u_epoch_;
        auto& removed_v = transposed ? removed_u_epoch_ : removed_v_epoch_;
        const auto cap_f = [&](VertexID id) {
            if (fc_epoch[id] != epoch) {
                fc_epoch[id] = epoch;
                fc[id] = fixed < fr[id].size() ? fr[id][fixed] : 0;
            }
            return fc[id];
        };
        const auto cap_v = [&](VertexID id) {
            if (vc_epoch[id] != epoch) {
                vc_epoch[id] = epoch;
                vc[id] = 0;
                const auto& row = vr[id];
                if (row.size() > 1 && row[1] >= fixed) {
                    Count lo = 1, hi = static_cast<Count>(row.size() - 1);
                    while (lo < hi) {
                        const Count mid = lo + (hi - lo + 1) / 2;
                        if (row[mid] >= fixed) lo = mid;
                        else hi = mid - 1;
                    }
                    vc[id] = lo;
                }
            }
            return vc[id];
        };
        Count boundary = 0;
        for (const auto& [u, v] : edges) {
            const auto f = transposed ? v : u, var = transposed ? u : v;
            boundary = std::max(boundary, std::min(cap_f(f), cap_v(var)));
        }
        if (boundary < minimum) return;
        if (profile != nullptr) ++profile->layered_sweeps;
        region.clear();
        queue.clear();
        heap.clear();
        const auto visit_f = [&](VertexID id) {
            const Count cap = cap_f(id);
            if (cap < minimum || cap > boundary || region_f[id] == epoch) return;
            region_f[id] = epoch;
            fd[id] = 0;
            region.emplace_back(true, id);
        };
        const auto visit_v = [&](VertexID id) {
            const Count cap = cap_v(id);
            if (cap < minimum || cap > boundary || region_v[id] == epoch) return;
            region_v[id] = epoch;
            vd[id] = 0;
            region.emplace_back(false, id);
        };
        for (const auto& [u, v] : edges) {
            const auto f = transposed ? v : u, var = transposed ? u : v;
            if (cap_f(f) < minimum || cap_v(var) < minimum) continue;
            visit_f(f);
            visit_v(var);
        }
        for (std::size_t head = 0; head < region.size(); ++head) {
            const auto [is_fixed, id] = region[head];
            if (is_fixed) {
                for (std::size_t i = 0; i < fa[id].size(); ++i) {
                    if (profile != nullptr) ++profile->layered_adjacency_entries;
                    if (!fe[id][i]) continue;
                    const VertexID v = fa[id][i];
                    const Count cap = cap_v(v);
                    if (cap < minimum) continue;
                    ++fd[id];
                    if (cap <= boundary) visit_v(v);
                    else if (profile != nullptr) ++profile->layered_anchor_edges;
                }
            } else {
                for (std::size_t i = 0; i < va[id].size(); ++i) {
                    if (profile != nullptr) ++profile->layered_adjacency_entries;
                    if (!ve[id][i]) continue;
                    const VertexID f = va[id][i];
                    const Count cap = cap_f(f);
                    if (cap < minimum) continue;
                    ++vd[id];
                    if (cap <= boundary) visit_f(f);
                    else if (profile != nullptr) ++profile->layered_anchor_edges;
                }
            }
        }
        if (profile != nullptr) profile->layered_region_vertices += region.size();
        const auto push_heap = [&](VertexID id) {
            heap.emplace_back(vd[id], id);
            std::push_heap(heap.begin(), heap.end(), heap_compare);
        };
        for (const auto& [is_fixed, id] : region) {
            if (is_fixed) {
                if (fd[id] < fixed) queue.emplace_back(true, id);
            } else {
                if (vd[id] < minimum) queue.emplace_back(false, id);
                push_heap(id);
            }
        }
        const auto drain = [&](Count level) {
            for (std::size_t head = 0; head < queue.size(); ++head) {
                const auto [is_fixed, id] = queue[head];
                if (is_fixed) {
                    if (removed_f[id] == epoch || fd[id] >= fixed) continue;
                    removed_f[id] = epoch;
                    if (profile != nullptr) ++profile->layered_removed_vertices;
                    record_frontier(!transposed, id, fixed, level, touched_u, touched_v);
                    for (std::size_t i = 0; i < fa[id].size(); ++i) {
                        if (profile != nullptr) ++profile->layered_adjacency_entries;
                        if (!fe[id][i]) continue;
                        const VertexID v = fa[id][i];
                        if (region_v[v] != epoch || removed_v[v] == epoch) continue;
                        --vd[v];
                        if (vd[v] <= level) queue.emplace_back(false, v);
                        push_heap(v);
                    }
                } else {
                    if (removed_v[id] == epoch || vd[id] > level) continue;
                    removed_v[id] = epoch;
                    if (profile != nullptr) ++profile->layered_removed_vertices;
                    Count parameter = next_unwritten(transposed, id, level + 1);
                    while (parameter <= vc[id]) {
                        record_frontier(transposed, id, parameter, fixed - 1,
                                        touched_u, touched_v);
                        parameter = next_unwritten(transposed, id, parameter);
                    }
                    for (std::size_t i = 0; i < va[id].size(); ++i) {
                        if (profile != nullptr) ++profile->layered_adjacency_entries;
                        if (!ve[id][i]) continue;
                        const VertexID f = va[id][i];
                        if (region_f[f] != epoch || removed_f[f] == epoch) continue;
                        --fd[f];
                        if (fd[f] < fixed) queue.emplace_back(true, f);
                    }
                }
            }
            queue.clear();
        };
        drain(minimum - 1);
        while (!heap.empty()) {
            std::pop_heap(heap.begin(), heap.end(), heap_compare);
            const auto [degree, id] = heap.back();
            heap.pop_back();
            if (removed_v[id] == epoch || degree != vd[id]) continue;
            if (degree >= boundary) break;
            queue.emplace_back(false, id);
            drain(degree);
        }
        layered_scratch_peak_bytes_ = std::max(layered_scratch_peak_bytes_,
            (region.capacity() + queue.capacity()) * sizeof(region[0]) +
            heap.capacity() * sizeof(HeapEntry));
    };
    for (Count k = 2; k <= delta_bound_; ++k) {
        sweep(false, k, k);
        sweep(true, k, k + 1);
    }
    if (profile != nullptr)
        profile->layered_seconds += Seconds(Clock::now() - start).count();
}

bool DynamicBicoreMaintainer::u_in_core(const BicoreSkyline& skyline,
                                        VertexID u, Count alpha,
                                        Count beta) const {
    return u < skyline.u_beta_max.size() &&
           alpha < skyline.u_beta_max[u].size() &&
           skyline.u_beta_max[u][alpha] >= beta;
}

bool DynamicBicoreMaintainer::v_in_core(const BicoreSkyline& skyline,
                                        VertexID v, Count alpha,
                                        Count beta) const {
    return v < skyline.v_alpha_max.size() &&
           beta < skyline.v_alpha_max[v].size() &&
           skyline.v_alpha_max[v][beta] >= alpha;
}

Count DynamicBicoreMaintainer::u_support(
    const BicoreSkyline& skyline, VertexID u, Count alpha, Count beta,
    const std::unordered_set<VertexID>& removed_v) const {
    Count support = 0;
    for (std::size_t i = 0; i < u_neighbors_[u].size(); ++i) {
        if (!u_edge_active_[u][i]) continue;
        const VertexID v = u_neighbors_[u][i];
        if (removed_v.find(v) == removed_v.end() &&
            v_in_core(skyline, v, alpha, beta)) {
            ++support;
        }
    }
    return support;
}

Count DynamicBicoreMaintainer::v_support(
    const BicoreSkyline& skyline, VertexID v, Count alpha, Count beta,
    const std::unordered_set<VertexID>& removed_u) const {
    Count support = 0;
    for (std::size_t i = 0; i < v_neighbors_[v].size(); ++i) {
        if (!v_edge_active_[v][i]) continue;
        const VertexID u = v_neighbors_[v][i];
        if (removed_u.find(u) == removed_u.end() &&
            u_in_core(skyline, u, alpha, beta)) {
            ++support;
        }
    }
    return support;
}

DynamicBicoreUpdate DynamicBicoreMaintainer::delete_edge(VertexID u,
                                                          VertexID v) {
    DynamicBicoreUpdate update;
    if (!contains_edge(u, v)) return update;

    const BicoreSkyline& old = skyline_;
    BicoreSkyline next = skyline_;
    deactivate_edge(u, v);

    std::unordered_set<VertexID> touched_u;
    std::unordered_set<VertexID> touched_v;
    for (Count alpha = 1; alpha < old.u_beta_max[u].size(); ++alpha) {
        const Count endpoint_beta_max = old.u_beta_max[u][alpha];
        for (Count beta = 1;
             beta <= endpoint_beta_max && beta <= old.max_beta; ++beta) {
            if (alpha == 1 && beta == 1) continue;
            if (!v_in_core(old, v, alpha, beta)) continue;

            std::unordered_set<VertexID> removed_u;
            std::unordered_set<VertexID> removed_v;
            std::queue<std::pair<Side, VertexID>> queue;
            if (u_support(old, u, alpha, beta, removed_v) < alpha) {
                queue.emplace(Side::U, u);
            }
            if (v_support(old, v, alpha, beta, removed_u) < beta) {
                queue.emplace(Side::V, v);
            }

            while (!queue.empty()) {
                const auto [side, id] = queue.front();
                queue.pop();
                if (side == Side::U) {
                    if (removed_u.find(id) != removed_u.end() ||
                        !u_in_core(old, id, alpha, beta) ||
                        u_support(old, id, alpha, beta, removed_v) >= alpha) {
                        continue;
                    }
                    removed_u.insert(id);
                    touched_u.insert(id);
                    for (std::size_t i = 0; i < u_neighbors_[id].size(); ++i) {
                        if (!u_edge_active_[id][i]) continue;
                        const VertexID neighbor = u_neighbors_[id][i];
                        if (removed_v.find(neighbor) == removed_v.end() &&
                            v_in_core(old, neighbor, alpha, beta)) {
                            queue.emplace(Side::V, neighbor);
                        }
                    }
                } else {
                    if (removed_v.find(id) != removed_v.end() ||
                        !v_in_core(old, id, alpha, beta) ||
                        v_support(old, id, alpha, beta, removed_u) >= beta) {
                        continue;
                    }
                    removed_v.insert(id);
                    touched_v.insert(id);
                    for (std::size_t i = 0; i < v_neighbors_[id].size(); ++i) {
                        if (!v_edge_active_[id][i]) continue;
                        const VertexID neighbor = v_neighbors_[id][i];
                        if (removed_u.find(neighbor) == removed_u.end() &&
                            u_in_core(old, neighbor, alpha, beta)) {
                            queue.emplace(Side::U, neighbor);
                        }
                    }
                }
            }

            for (VertexID removed : removed_u) {
                next.u_beta_max[removed][alpha] = std::min(
                    next.u_beta_max[removed][alpha], beta - 1);
            }
            for (VertexID removed : removed_v) {
                next.v_alpha_max[removed][beta] = std::min(
                    next.v_alpha_max[removed][beta], alpha - 1);
            }
        }
    }

    update.affected_u = touched_u.size();
    update.affected_v = touched_v.size();
    for (VertexID x : touched_u) {
        bool changed = false;
        for (Count alpha = 1; alpha < old.u_beta_max[x].size(); ++alpha) {
            if (old.u_beta_max[x][alpha] == next.u_beta_max[x][alpha]) continue;
            if (old.u_beta_max[x][alpha] != 0) --frontier_entries_;
            if (next.u_beta_max[x][alpha] != 0) ++frontier_entries_;
            std::size_t lost = old.u_beta_max[x][alpha] -
                next.u_beta_max[x][alpha];
            if (alpha == 1 && next.u_beta_max[x][alpha] == 0 &&
                old.u_beta_max[x][alpha] != 0) {
                --lost;
            }
            membership_entries_ -= lost;
            ++update.changed_memberships;
            changed = true;
        }
        if (changed) update.changed_u.push_back(x);
    }
    for (VertexID x : touched_v) {
        bool changed = false;
        for (Count beta = 1; beta < old.v_alpha_max[x].size(); ++beta) {
            if (old.v_alpha_max[x][beta] == next.v_alpha_max[x][beta]) continue;
            if (old.v_alpha_max[x][beta] != 0) --frontier_entries_;
            if (next.v_alpha_max[x][beta] != 0) ++frontier_entries_;
            std::size_t lost = old.v_alpha_max[x][beta] -
                next.v_alpha_max[x][beta];
            if (beta == 1 && next.v_alpha_max[x][beta] == 0 &&
                old.v_alpha_max[x][beta] != 0) {
                --lost;
            }
            membership_entries_ -= lost;
            ++update.changed_memberships;
            changed = true;
        }
        if (changed) update.changed_v.push_back(x);
    }
    std::sort(update.changed_u.begin(), update.changed_u.end());
    std::sort(update.changed_v.begin(), update.changed_v.end());
    skyline_ = std::move(next);
    return update;
}

DynamicBicoreUpdate DynamicBicoreMaintainer::delete_edge_bucket(
    const std::vector<std::pair<VertexID, VertexID>>& edges,
    DynamicBicoreProfile* profile) {
    DynamicBicoreUpdate bucket_update;
    if (edges.empty()) return bucket_update;
    const auto candidate_start = profile != nullptr ? Clock::now()
                                                    : Clock::time_point{};

    const BicoreSkyline& before = skyline_;
    struct RemovedEdge {
        VertexID u;
        VertexID v;
        Count alpha_limit;
    };
    std::vector<RemovedEdge> removed_edges;
    std::vector<VertexID> endpoints_u;
    std::vector<VertexID> endpoints_v;
    std::vector<std::pair<VertexID, VertexID>> layered_edges;

    if (bucket_epoch_ == std::numeric_limits<std::uint32_t>::max()) {
        std::fill(pending_u_epoch_.begin(), pending_u_epoch_.end(), 0);
        std::fill(pending_v_epoch_.begin(), pending_v_epoch_.end(), 0);
        std::fill(axis_u_epoch_.begin(), axis_u_epoch_.end(), 0);
        std::fill(axis_v_epoch_.begin(), axis_v_epoch_.end(), 0);
        bucket_epoch_ = 1;
    } else {
        ++bucket_epoch_;
    }

    for (const auto& [u, v] : edges) {
        if (!deactivate_edge(u, v)) continue;
        if (layered_enabled_) layered_edges.emplace_back(u, v);
        if (axis_u_epoch_[u] != bucket_epoch_) {
            axis_u_epoch_[u] = bucket_epoch_;
            endpoints_u.push_back(u);
            endpoint_u_beta_max_[u] = 0;
        }
        if (axis_v_epoch_[v] != bucket_epoch_) {
            axis_v_epoch_[v] = bucket_epoch_;
            endpoints_v.push_back(v);
            endpoint_v_beta_max_[v] = 0;
        }
        ++endpoint_u_beta_max_[u];
        ++endpoint_v_beta_max_[v];
        const auto& v_frontier = before.v_alpha_max[v];
        if (before.u_beta_max[u].size() > 1 && v_frontier.size() > 1) {
            removed_edges.push_back({u, v, std::min(
                static_cast<Count>(before.u_beta_max[u].size() - 1),
                v_frontier[1])});
        }
    }
    if (removed_edges.empty()) return bucket_update;
    auto sort_unique = [](auto& values) {
        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end()), values.end());
    };
    sort_unique(endpoints_u);
    sort_unique(endpoints_v);
    std::sort(removed_edges.begin(), removed_edges.end(),
              [](const RemovedEdge& lhs, const RemovedEdge& rhs) {
                  return lhs.alpha_limit > rhs.alpha_limit;
              });
    if (profile != nullptr) {
        profile->candidate_seconds +=
            Seconds(Clock::now() - candidate_start).count();
    }

    std::vector<std::pair<VertexID, Count>> touched_u;
    std::vector<std::pair<VertexID, Count>> touched_v;

    const auto axis_start = profile != nullptr ? Clock::now() : Clock::time_point{};
    std::vector<VertexID> axis_u = endpoints_u;
    std::vector<VertexID> axis_v = endpoints_v;
    for (VertexID u : endpoints_u) {
        const Count old_degree = live_u_degree_[u] + endpoint_u_beta_max_[u];
        for (Count a = live_u_degree_[u] + 1; a <= old_degree; ++a)
            record_frontier(true, u, a, 0, touched_u, touched_v);
        for (std::size_t i = 0; i < u_neighbors_[u].size(); ++i) {
            if (profile != nullptr) ++profile->axis_adjacency_entries;
            if (!u_edge_active_[u][i]) continue;
            const VertexID v = u_neighbors_[u][i];
            if (axis_v_epoch_[v] != bucket_epoch_ &&
                before.v_alpha_max[v][1] == old_degree) {
                axis_v_epoch_[v] = bucket_epoch_;
                axis_v.push_back(v);
            }
        }
    }
    for (VertexID v : endpoints_v) {
        const Count old_degree = live_v_degree_[v] + endpoint_v_beta_max_[v];
        for (Count b = live_v_degree_[v] + 1; b <= old_degree; ++b)
            record_frontier(false, v, b, 0, touched_u, touched_v);
        for (std::size_t i = 0; i < v_neighbors_[v].size(); ++i) {
            if (profile != nullptr) ++profile->axis_adjacency_entries;
            if (!v_edge_active_[v][i]) continue;
            const VertexID u = v_neighbors_[v][i];
            if (axis_u_epoch_[u] != bucket_epoch_ &&
                before.u_beta_max[u][1] == old_degree) {
                axis_u_epoch_[u] = bucket_epoch_;
                axis_u.push_back(u);
            }
        }
    }
    for (VertexID u : axis_u) {
        Count bound = 0;
        for (std::size_t i = 0; i < u_neighbors_[u].size(); ++i) {
            if (profile != nullptr) ++profile->axis_adjacency_entries;
            if (!u_edge_active_[u][i]) continue;
            bound = std::max(bound, live_v_degree_[u_neighbors_[u][i]]);
            if (bound == before.u_beta_max[u][1]) break;
        }
        record_frontier(true, u, 1, bound, touched_u, touched_v);
    }
    for (VertexID v : axis_v) {
        Count bound = 0;
        for (std::size_t i = 0; i < v_neighbors_[v].size(); ++i) {
            if (profile != nullptr) ++profile->axis_adjacency_entries;
            if (!v_edge_active_[v][i]) continue;
            bound = std::max(bound, live_u_degree_[v_neighbors_[v][i]]);
            if (bound == before.v_alpha_max[v][1]) break;
        }
        record_frontier(false, v, 1, bound, touched_u, touched_v);
    }
    if (profile != nullptr) {
        profile->axis_seconds += Seconds(Clock::now() - axis_start).count();
        profile->axis_frontier_changes += touched_u.size() + touched_v.size();
    }

    std::vector<VertexID> removed_u;
    std::vector<VertexID> removed_v;
    std::vector<std::pair<Side, VertexID>> queue;
    if (layered_enabled_) {
        maintain_layered(layered_edges, touched_u, touched_v, profile);
    }
    for (Count alpha = 2;
         !layered_enabled_ && alpha <= removed_edges.front().alpha_limit; ++alpha) {
        const auto endpoint_start = profile != nullptr ? Clock::now()
                                                       : Clock::time_point{};
        for (VertexID u : endpoints_u) endpoint_u_beta_max_[u] = 0;
        for (VertexID v : endpoints_v) endpoint_v_beta_max_[v] = 0;
        Count affected_beta_max = 0;
        for (const RemovedEdge& edge : removed_edges) {
            if (edge.alpha_limit < alpha) break;
            const auto& v_frontier = before.v_alpha_max[edge.v];
            const Count beta_limit = std::min({
                before.u_beta_max[edge.u][alpha], before.max_beta,
                static_cast<Count>(v_frontier.size() - 1)});
            if (beta_limit == 0) continue;
            Count low = 1;
            Count high = beta_limit;
            while (low < high) {
                const Count middle = low + (high - low + 1) / 2;
                if (v_frontier[middle] >= alpha) low = middle;
                else high = middle - 1;
            }
            endpoint_u_beta_max_[edge.u] =
                std::max(endpoint_u_beta_max_[edge.u], low);
            endpoint_v_beta_max_[edge.v] =
                std::max(endpoint_v_beta_max_[edge.v], low);
            affected_beta_max = std::max(affected_beta_max, low);
        }
        if (profile != nullptr) {
            if (affected_beta_max > 1)
                profile->candidate_planes += affected_beta_max - 1;
            profile->candidate_seconds +=
                Seconds(Clock::now() - endpoint_start).count();
        }
        for (Count beta = 2; beta <= affected_beta_max; ++beta) {
            const auto support_start = profile != nullptr ? Clock::now()
                                                          : Clock::time_point{};
            const std::uint32_t epoch = next_plane_epoch();
            removed_u.clear();
            removed_v.clear();
            queue.clear();
    
            const auto is_removed_u = [&](VertexID u) {
                return removed_u_epoch_[u] == epoch;
            };
            const auto is_removed_v = [&](VertexID v) {
                return removed_v_epoch_[v] == epoch;
            };
            const auto current_u_support = [&](VertexID u) {
                if (support_u_epoch_[u] == epoch) return support_u_value_[u];
                Count support = 0;
                for (std::size_t i = 0; i < u_neighbors_[u].size(); ++i) {
                    if (profile != nullptr) ++profile->support_adjacency_entries;
                    if (!u_edge_active_[u][i]) continue;
                    const VertexID v = u_neighbors_[u][i];
                    if (!is_removed_v(v) &&
                        v_in_core(before, v, alpha, beta)) {
                        ++support;
                    }
                }
                support_u_epoch_[u] = epoch;
                support_u_value_[u] = support;
                return support;
            };
            const auto current_v_support = [&](VertexID v) {
                if (support_v_epoch_[v] == epoch) return support_v_value_[v];
                Count support = 0;
                for (std::size_t i = 0; i < v_neighbors_[v].size(); ++i) {
                    if (profile != nullptr) ++profile->support_adjacency_entries;
                    if (!v_edge_active_[v][i]) continue;
                    const VertexID u = v_neighbors_[v][i];
                    if (!is_removed_u(u) &&
                        u_in_core(before, u, alpha, beta)) {
                        ++support;
                    }
                }
                support_v_epoch_[v] = epoch;
                support_v_value_[v] = support;
                return support;
            };
            const auto u_below_threshold = [&](VertexID u) {
                Count support = 0;
                for (std::size_t i = 0; i < u_neighbors_[u].size(); ++i) {
                    if (profile != nullptr) ++profile->support_adjacency_entries;
                    if (!u_edge_active_[u][i]) continue;
                    const VertexID v = u_neighbors_[u][i];
                    if (v_in_core(before, v, alpha, beta) && ++support >= alpha) {
                        return false;
                    }
                }
                support_u_epoch_[u] = epoch;
                support_u_value_[u] = support;
                return true;
            };
            const auto v_below_threshold = [&](VertexID v) {
                Count support = 0;
                for (std::size_t i = 0; i < v_neighbors_[v].size(); ++i) {
                    if (profile != nullptr) ++profile->support_adjacency_entries;
                    if (!v_edge_active_[v][i]) continue;
                    const VertexID u = v_neighbors_[v][i];
                    if (u_in_core(before, u, alpha, beta) && ++support >= beta) {
                        return false;
                    }
                }
                support_v_epoch_[v] = epoch;
                support_v_value_[v] = support;
                return true;
            };
            const auto enqueue_u = [&](VertexID u) {
                if (is_removed_u(u) || queued_u_epoch_[u] == epoch) return;
                queued_u_epoch_[u] = epoch;
                queue.emplace_back(Side::U, u);
            };
            const auto enqueue_v = [&](VertexID v) {
                if (is_removed_v(v) || queued_v_epoch_[v] == epoch) return;
                queued_v_epoch_[v] = epoch;
                queue.emplace_back(Side::V, v);
            };
    
            for (VertexID u : endpoints_u) {
                if (endpoint_u_beta_max_[u] >= beta &&
                    u_below_threshold(u)) {
                    enqueue_u(u);
                }
            }
            for (VertexID v : endpoints_v) {
                if (endpoint_v_beta_max_[v] >= beta &&
                    v_below_threshold(v)) {
                    enqueue_v(v);
                }
            }
    
            if (profile != nullptr) {
                profile->support_seconds +=
                    Seconds(Clock::now() - support_start).count();
                if (!queue.empty()) ++profile->seeded_planes;
            }
            const auto cascade_start = profile != nullptr ? Clock::now()
                                                          : Clock::time_point{};
            std::size_t queue_head = 0;
            while (queue_head < queue.size()) {
                if (profile != nullptr) ++profile->cascade_queue_pops;
                const auto [side, id] = queue[queue_head++];
                if (side == Side::U) {
                    queued_u_epoch_[id] = 0;
                    if (is_removed_u(id) ||
                        !u_in_core(before, id, alpha, beta) ||
                        current_u_support(id) >= alpha) {
                        continue;
                    }
                    removed_u_epoch_[id] = epoch;
                    removed_u.push_back(id);
                    for (std::size_t i = 0; i < u_neighbors_[id].size(); ++i) {
                        if (!u_edge_active_[id][i]) continue;
                        const VertexID neighbor = u_neighbors_[id][i];
                        if (is_removed_v(neighbor) ||
                            !v_in_core(before, neighbor, alpha, beta)) {
                            continue;
                        }
                        if (support_v_epoch_[neighbor] == epoch &&
                            support_v_value_[neighbor] != 0) {
                            --support_v_value_[neighbor];
                        }
                        enqueue_v(neighbor);
                    }
                } else {
                    queued_v_epoch_[id] = 0;
                    if (is_removed_v(id) ||
                        !v_in_core(before, id, alpha, beta) ||
                        current_v_support(id) >= beta) {
                        continue;
                    }
                    removed_v_epoch_[id] = epoch;
                    removed_v.push_back(id);
                    for (std::size_t i = 0; i < v_neighbors_[id].size(); ++i) {
                        if (!v_edge_active_[id][i]) continue;
                        const VertexID neighbor = v_neighbors_[id][i];
                        if (is_removed_u(neighbor) ||
                            !u_in_core(before, neighbor, alpha, beta)) {
                            continue;
                        }
                        if (support_u_epoch_[neighbor] == epoch &&
                            support_u_value_[neighbor] != 0) {
                            --support_u_value_[neighbor];
                        }
                        enqueue_u(neighbor);
                    }
                }
            }
    
            if (profile != nullptr) {
                profile->cascade_seconds +=
                    Seconds(Clock::now() - cascade_start).count();
                if (!removed_u.empty() || !removed_v.empty()) {
                    ++profile->changed_planes;
                }
                profile->cascade_removed_vertices +=
                    removed_u.size() + removed_v.size();
            }
            for (VertexID u : removed_u) {
                record_frontier(true, u, alpha, beta - 1, touched_u, touched_v);
            }
            for (VertexID v : removed_v) {
                record_frontier(false, v, beta, alpha - 1, touched_u, touched_v);
            }
        }
    }

    const auto merge_start = profile != nullptr ? Clock::now()
                                                : Clock::time_point{};
    std::sort(touched_u.begin(), touched_u.end());
    for (const auto& [u, alpha] : touched_u) {
        const Count old_beta = skyline_.u_beta_max[u][alpha];
        const Count new_beta =
            pending_u_new_beta_[u_frontier_offsets_[u] + alpha];
        if (new_beta >= old_beta) continue;
        if (old_beta != 0) --frontier_entries_;
        if (new_beta != 0) ++frontier_entries_;
        std::size_t lost = old_beta - new_beta;
        if (alpha == 1 && new_beta == 0 && old_beta != 0) --lost;
        membership_entries_ -= lost;
        skyline_.u_beta_max[u][alpha] = new_beta;
        if (bucket_update.changed_u.empty() ||
            bucket_update.changed_u.back() != u) {
            bucket_update.changed_u.push_back(u);
        }
        bucket_update.u_frontier_changes.push_back(
            {u, alpha, old_beta, new_beta});
        bucket_update.changed_memberships += old_beta - new_beta;
    }

    std::sort(touched_v.begin(), touched_v.end());
    for (const auto& [v, beta] : touched_v) {
        const Count old_alpha = skyline_.v_alpha_max[v][beta];
        const Count new_alpha =
            pending_v_new_alpha_[v_frontier_offsets_[v] + beta];
        if (new_alpha >= old_alpha) continue;
        if (old_alpha != 0) --frontier_entries_;
        if (new_alpha != 0) ++frontier_entries_;
        std::size_t lost = old_alpha - new_alpha;
        if (beta == 1 && new_alpha == 0 && old_alpha != 0) --lost;
        membership_entries_ -= lost;
        skyline_.v_alpha_max[v][beta] = new_alpha;
        if (bucket_update.changed_v.empty() ||
            bucket_update.changed_v.back() != v) {
            bucket_update.changed_v.push_back(v);
        }
        bucket_update.v_frontier_changes.push_back(
            {v, beta, old_alpha, new_alpha});
        bucket_update.changed_memberships += old_alpha - new_alpha;
    }
    bucket_update.affected_u = bucket_update.changed_u.size();
    bucket_update.affected_v = bucket_update.changed_v.size();
    if (profile != nullptr) {
        profile->merge_seconds += Seconds(Clock::now() - merge_start).count();
    }
    return bucket_update;
}

std::size_t DynamicBicoreMaintainer::approximate_bytes() const {
    std::size_t bytes = sizeof(*this);
    for (const auto& row : u_neighbors_)
        bytes += sizeof(row) + row.size() * sizeof(VertexID);
    for (const auto& row : v_neighbors_)
        bytes += sizeof(row) + row.size() * sizeof(VertexID);
    for (const auto& row : u_edge_active_)
        bytes += sizeof(row) + row.capacity() * sizeof(std::uint8_t);
    for (const auto& row : v_edge_active_)
        bytes += sizeof(row) + row.capacity() * sizeof(std::uint8_t);
    for (const auto& row : skyline_.u_beta_max)
        bytes += sizeof(row) + row.capacity() * sizeof(Count);
    for (const auto& row : skyline_.v_alpha_max)
        bytes += sizeof(row) + row.capacity() * sizeof(Count);
    bytes += (removed_u_epoch_.capacity() + removed_v_epoch_.capacity() +
              queued_u_epoch_.capacity() + queued_v_epoch_.capacity() +
              support_u_epoch_.capacity() + support_v_epoch_.capacity()) *
        sizeof(std::uint32_t);
    bytes += (axis_u_epoch_.capacity() + axis_v_epoch_.capacity()) *
        sizeof(std::uint32_t);
    bytes += (live_u_degree_.capacity() + live_v_degree_.capacity()) * sizeof(Count);
    bytes += (support_u_value_.capacity() + support_v_value_.capacity() +
              endpoint_u_beta_max_.capacity() + endpoint_v_beta_max_.capacity() +
              pending_u_new_beta_.capacity() +
              pending_v_new_alpha_.capacity()) * sizeof(Count);
    bytes += (pending_u_epoch_.capacity() + pending_v_epoch_.capacity()) *
        sizeof(std::uint32_t);
    bytes += (u_frontier_offsets_.capacity() + v_frontier_offsets_.capacity()) *
        sizeof(std::size_t);
    bytes += (successor_u_.capacity() + successor_v_.capacity()) * sizeof(Count);
    return bytes;
}

}
