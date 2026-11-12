#include "query/index_query_utils.h"

#include "core/bicore_decomposition.h"
#include "query/online_query.h"
#include "query/mat_recovery.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace tacc {

SeedRecovery configured_seed_recovery() {
    static const SeedRecovery recovery = [] {
        const char* env = std::getenv("TAC_SEED_RECOVERY");
        if (env == nullptr || std::string(env) == "adaptive")
            return SeedRecovery::Adaptive;
        if (std::string(env) == "seed") return SeedRecovery::Seed;
        throw std::invalid_argument("Invalid TAC_SEED_RECOVERY; expected adaptive|seed");
    }();
    return recovery;
}

const char* seed_recovery_name(SeedRecovery recovery) {
    return recovery == SeedRecovery::Adaptive ? "adaptive" : "seed";
}

UATRecovery configured_uat_recovery() {
    static const UATRecovery recovery=[] {
        const char* env=std::getenv("TAC_UAT_RECOVERY");
        if (env==nullptr || std::string(env)=="adaptive") return UATRecovery::Adaptive;
        if (std::string(env)=="grouped") return UATRecovery::Grouped;
        throw std::invalid_argument("Invalid TAC_UAT_RECOVERY; expected adaptive|grouped");
    }();
    return recovery;
}

const char* uat_recovery_name(UATRecovery recovery) {
    return recovery==UATRecovery::Grouped ? "grouped" : "adaptive";
}

bool index_uses_grouped_recovery(const std::string& kind) {
    if (kind=="uat") return configured_uat_recovery()==UATRecovery::Grouped;
    return (kind=="vat" || kind=="mat") && configured_mat_recovery()==MATRecovery::Grouped;
}

QueryRoute configured_query_route() {
    static const QueryRoute route = [] {
        const char* env = std::getenv("TAC_QUERY_ROUTE");
        if (env == nullptr) return QueryRoute::UV;
        const std::string value(env);
        if (value == "auto") return QueryRoute::Auto;
        if (value == "uv") return QueryRoute::UV;
        throw std::invalid_argument("Invalid TAC_QUERY_ROUTE='" + value +
            "'; this build supports auto|uv only");
    }();
    return route;
}

const char* query_route_name(QueryRoute route) {
    switch (route) {
    case QueryRoute::UV: return "uv";
    case QueryRoute::Auto: return "auto";
    }
    throw std::invalid_argument("Unknown query route");
}

QueryRoute effective_query_route(IndexScope) {
    (void)configured_query_route();
    return QueryRoute::UV;
}

namespace {

std::chrono::steady_clock::time_point profile_start(
    const QueryProfile* profile) {
    return profile == nullptr
        ? std::chrono::steady_clock::time_point{}
        : std::chrono::steady_clock::now();
}

AttributeMask attribute_mask(const std::vector<AttributeID>& attributes) {
    AttributeMask mask = 0;
    for (AttributeID attribute : attributes) {
        if (attribute < 64) mask |= AttributeMask{1} << attribute;
    }
    return mask;
}

struct MarkWorkspace {
    std::vector<std::uint32_t> v_marks;
    std::uint32_t generation = 0;

    std::pair<std::vector<std::uint32_t>*, std::uint32_t> fresh(
        std::size_t size) {
        if (generation == std::numeric_limits<std::uint32_t>::max()) {
            std::fill(v_marks.begin(), v_marks.end(), 0);
            generation = 0;
        }
        ++generation;
        auto& marks = v_marks;
        if (marks.size() < size) marks.resize(size, 0);
        return {&marks, generation};
    }
};

MarkWorkspace& mark_workspace() {
    thread_local MarkWorkspace workspace;
    return workspace;
}

std::pair<std::vector<std::uint32_t>*, std::uint32_t> fresh_marks(
    std::size_t size) {
    return mark_workspace().fresh(size);
}

void mark_vertex(std::vector<std::uint32_t>& marks,
                 std::uint32_t generation,
                 VertexID vertex,
                 bool track_sparse_hits,
                 std::vector<VertexID>& selected) {
    if (marks[vertex] == generation) return;
    marks[vertex] = generation;
    if (track_sparse_hits) selected.push_back(vertex);
}

std::vector<VertexID> finish_vertices(
    std::vector<VertexID> selected,
    const std::vector<std::uint32_t>& marks,
    std::uint32_t generation,
    std::size_t vertex_count) {
    std::size_t log_factor = 0;
    for (std::size_t value = selected.size(); value > 1; value >>= 1) {
        ++log_factor;
    }
    if (log_factor != 0 &&
        selected.size() > vertex_count / log_factor) {
        selected.clear();
        for (VertexID vertex = 0; vertex < vertex_count; ++vertex) {
            if (marks[vertex] == generation) selected.push_back(vertex);
        }
    } else {
        std::sort(selected.begin(), selected.end());
    }
    return selected;
}

class VertexCollector {
public:
    VertexCollector(std::size_t vertex_count)
        : vertex_count_(vertex_count),
          sparse_(vertex_count > 4096) {
        if (sparse_) {
            auto fresh = fresh_marks(vertex_count);
            marks_ = fresh.first;
            generation_ = fresh.second;
        } else {
            dense_.assign(vertex_count, 0);
        }
    }

    void add(VertexID vertex) {
        if (sparse_) {
            mark_vertex(*marks_, generation_, vertex, true, selected_);
        } else {
            dense_[vertex] = 1;
        }
    }

    std::vector<VertexID> finish() {
        if (sparse_) {
            return finish_vertices(std::move(selected_), *marks_,
                                   generation_, vertex_count_);
        }
        std::vector<VertexID> result;
        for (VertexID vertex = 0; vertex < vertex_count_; ++vertex) {
            if (dense_[vertex]) result.push_back(vertex);
        }
        return result;
    }

private:
    std::size_t vertex_count_{};
    bool sparse_{};
    std::vector<std::uint8_t> dense_;
    std::vector<std::uint32_t>* marks_{};
    std::uint32_t generation_{};
    std::vector<VertexID> selected_;
};

template <typename Directory, typename Key>
const IndexPostingRange* find_posting_range(const Directory& directory,
                                            const Key& key) {
    const auto found = std::lower_bound(
        directory.begin(), directory.end(), key,
        [](const typename Directory::value_type& entry, const Key& target) {
            return entry.first < target;
        });
    if (found == directory.end() || found->first != key) return nullptr;
    return &found->second;
}

QueryResult peel_snapshot(const TemporalBipartiteGraph& graph,
                          const Query& query,
                          BipartiteSubgraph& snapshot,
                          const std::vector<VertexID>& candidate_u,
                          const std::vector<VertexID>& candidate_v,
                          QueryProfile* profile,
                          QueryTimer* timer = nullptr) {
    const auto peeling_start = profile_start(profile);
    const QueryResult internal = bicore_peeling_candidates(
        snapshot, query.alpha, query.beta, candidate_u, candidate_v);
    if (profile != nullptr) {
        profile->peeling_us +=
            std::chrono::duration<double, std::micro>(
                std::chrono::steady_clock::now() - peeling_start).count();
    }
    QueryResult result;
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

QueryResult peel_seed_candidates(const TemporalBipartiteGraph& graph,
                                 const Query& query,
                                 std::vector<VertexID> seeds,
                                 AttributeMask u_mask, AttributeMask v_mask,
                                 QueryProfile* profile, bool mat_uv, bool force_grouped,
                                 QueryTimer* timer = nullptr) {
    if (seeds.empty()) { stop_query_timer(timer); return {}; }
    {
        const auto attribute_start = profile_start(profile);
        const AttributeMask allowed = v_mask;
        seeds.erase(std::remove_if(seeds.begin(), seeds.end(), [&](VertexID id) {
            const AttributeID attr = graph.v_attribute(id);
            return attr >= 64 || (allowed & (AttributeMask{1} << attr)) == 0;
        }), seeds.end());
        if (profile != nullptr) profile->attribute_filter_us +=
            std::chrono::duration<double, std::micro>(
                std::chrono::steady_clock::now() - attribute_start).count();
        if (seeds.empty()) { stop_query_timer(timer); return {}; }
    }
    if (force_grouped || (mat_uv && configured_mat_recovery() == MATRecovery::Grouped))
        return query_mat_grouped_seeds(graph, query, std::move(seeds), u_mask, profile, timer);
    const auto recovery_start = profile_start(profile);
    BipartiteSubgraph snapshot(graph.num_u_vertices(), graph.num_v_vertices());
    auto& marks = snapshot.v_selected;
    for (VertexID seed : seeds) marks[seed] = 1;
    std::vector<VertexID> neighbors;
    const auto stats = graph.populate_seed_window_adjacency(
        snapshot, query.window, seeds, false,
        u_mask, neighbors,
        configured_seed_recovery() == SeedRecovery::Adaptive);
    const auto& candidate_u = neighbors;
    const auto& candidate_v = seeds;
    if (profile != nullptr) {
        profile->candidate_u = candidate_u.size();
        profile->candidate_v = candidate_v.size();
        profile->window_recovery_us += std::chrono::duration<double, std::micro>(
            std::chrono::steady_clock::now() - recovery_start).count();
        profile->incident_work_u = stats.incident_work_u;
        profile->incident_work_v = stats.incident_work_v;
        profile->unique_edges_examined = stats.unique_edges_examined;
        profile->window_edges = stats.window_edges;
        profile->adaptive_seed_recoveries = stats.adaptive_seed_recovery ? 1 : 0;
        profile->reversed_seed_scans = stats.reversed_seed_scan ? 1 : 0;
        profile->seed_scan_work = stats.seed_scan_work;
        profile->opposite_scan_work = stats.opposite_scan_work;
    }
    if (candidate_u.empty() || candidate_v.empty()) { stop_query_timer(timer); return {}; }
    return peel_snapshot(graph, query, snapshot, candidate_u, candidate_v, profile, timer);
}

template <typename Collect>
QueryResult query_uv_index(const TemporalBipartiteGraph& graph,
                               const Query& query, IndexScope scope, Collect collect,
                               QueryProfile* profile, bool mat_uv = false, bool force_grouped = false,
                               const EncodedAttributes* attributes = nullptr, QueryTimer* timer = nullptr) {
    (void)effective_query_route(scope);
    if (query.window.start > query.window.end || query.alpha == 0 || query.beta == 0) {
        stop_query_timer(timer);
        return {};
    }
    if (query.alpha == 1 && query.beta == 1) {
        QueryResult result = online_query(graph, query, profile, timer);
        if (profile != nullptr) profile->online_fallbacks = 1;
        return result;
    }
    const AttributeMask u_mask = attributes ? attributes->u : attribute_mask(query.u_attributes);
    const AttributeMask v_mask = attributes ? attributes->v : attribute_mask(query.v_attributes);
    if (u_mask == 0 || v_mask == 0) { stop_query_timer(timer); return {}; }
    auto seeds = collect(u_mask);
    if (profile != nullptr) profile->directional_candidates = seeds.size();
    return peel_seed_candidates(graph, query, std::move(seeds), u_mask, v_mask, profile, mat_uv, force_grouped, timer);
}

}

QueryResult query_membership_index(
    const VATIndex& index, const TemporalBipartiteGraph& graph,
    const Query& query, QueryProfile* profile) {
    if (profile != nullptr) *profile = {};
    if (!index.ready()) throw std::logic_error("VAT query before finalize");
    const auto& intervals = index.intervals();
    return query_uv_index(
        graph, query, index.stats().scope, [&](AttributeMask attrs) {
            const auto& tree = index.vertices();
            std::vector<VertexID> result;
            const auto allowed = attribute_mask(query.v_attributes);
            const auto key = std::make_pair(query.alpha, query.beta);
            const auto collect_start = profile_start(profile);
            std::vector<std::uint8_t> selected(graph.num_v_vertices(), 0);
            const auto& cores = index.query_cores();
            const auto& refs = index.query_refs();
            for (const auto& m : index.query_masks()) {
                if (profile != nullptr) ++profile->directory_keys;
                if (!attribute_subset(m.mask, attrs)) continue;
                const auto first = cores.begin() + m.offset, last = first + m.count;
                const auto c = std::lower_bound(first, last, key,
                    [](const auto& c, const auto& k) { return std::make_pair(c.alpha, c.beta) < k; });
                if (c == last || c->alpha != query.alpha || c->beta != query.beta) continue;
                if (profile != nullptr) ++profile->directory_hits;
                for (std::size_t i = c->offset; i < c->offset + c->count; ++i) {
                    const auto& ref = refs[i];
                    if (profile != nullptr) ++profile->posting_entries;
                    if (selected[ref.vertex]) continue;
                    const auto attr = graph.v_attribute(ref.vertex);
                    if (attr >= 64 || !(allowed & (AttributeMask{1} << attr))) continue;
                    const auto& leaf = tree[ref.vertex][ref.attribute_slot].cores[ref.core_slot];
                    if (profile != nullptr) ++profile->vertex_refs;
                    const auto begin = intervals.begin() +
                        static_cast<std::ptrdiff_t>(leaf.interval_offset);
                    const auto end = begin +
                        static_cast<std::ptrdiff_t>(leaf.interval_count);
                    const auto first = std::lower_bound(
                        begin, end, query.window.start,
                        [](const VATInterval& interval, Timestamp start) {
                            return interval.start < start;
                        });
                    if (first != end &&
                        first->end <= query.window.end) {
                        selected[ref.vertex] = 1;
                        result.push_back(ref.vertex);
                    }
                }
            }
            std::sort(result.begin(), result.end());
            if (profile != nullptr) {
                profile->candidate_collect_us +=
                    std::chrono::duration<double, std::micro>(
                        std::chrono::steady_clock::now() - collect_start).count();
            }
            return result;
        }, profile, true);
}

QueryResult query_uat_block_index(
    const UATIndex& index, const TemporalBipartiteGraph& graph,
    const Query& query, QueryProfile* profile) {
    if (profile != nullptr) *profile = {};
    (void)effective_query_route(index.stats().scope);
    const auto recovery=configured_uat_recovery();
    if (!index.ready()) throw std::logic_error("UAT query before finalize");
    const auto& cores = index.cores();
    const auto& vertex_pool = index.vertex_pool();
    const auto& times = graph.time_points();
    const auto effective =
        std::lower_bound(times.begin(), times.end(), query.window.start);
    if (effective == times.end() || *effective > query.window.end) return {};
    const Timestamp effective_start = *effective;
    return query_uv_index(
        graph, query, index.stats().scope, [&](AttributeMask attrs) {
            const std::size_t vertex_count = graph.num_v_vertices();
            VertexCollector selected(vertex_count);
            const auto lookup_start = profile_start(profile);
            const auto key=std::make_pair(query.alpha,query.beta);
            const auto core=std::lower_bound(cores.begin(),cores.end(),key,
                [](const UATCoreEntry& c,const auto& k) {
                    return std::make_pair(c.alpha,c.beta)<k;
                });
            if (profile != nullptr) {
                profile->directory_lookup_us +=
                    std::chrono::duration<double, std::micro>(
                        std::chrono::steady_clock::now() - lookup_start).count();
            }
            const auto collect_start = profile_start(profile);
            if (core==cores.end() || core->alpha!=query.alpha || core->beta!=query.beta)
                return selected.finish();
            for (const auto& attribute : core->attributes) {
                if (profile) ++profile->directory_keys;
                if (!attribute_subset(attribute.mask,attrs)) continue;
                if (profile) ++profile->directory_hits;
                const auto start=std::lower_bound(attribute.starts.begin(),attribute.starts.end(),effective_start,
                    [](const UATStartEntry& s,Timestamp t) {return s.start<t;});
                if (start==attribute.starts.end() || start->start!=effective_start) continue;
                const auto eligible_end = std::upper_bound(
                    start->ends.begin(),start->ends.end(),query.window.end,
                    [](Timestamp t,const UATEndEntry& e) {return t<e.end;});
                if (profile != nullptr) {
                    profile->posting_entries += static_cast<std::size_t>(
                        std::distance(start->ends.begin(), eligible_end));
                }
                for (auto it = start->ends.begin(); it != eligible_end; ++it) {
                    const UATEndEntry& entry = *it;
                    if (profile != nullptr) {
                        profile->vertex_refs += entry.vertex_count;
                    }
                    const std::size_t end =
                        entry.vertex_offset + entry.vertex_count;
                    for (std::size_t i = entry.vertex_offset; i < end; ++i) {
                        selected.add(vertex_pool[i]);
                    }
                }
            }
            std::vector<VertexID> result = selected.finish();
            if (profile != nullptr) {
                profile->candidate_collect_us +=
                    std::chrono::duration<double, std::micro>(
                        std::chrono::steady_clock::now() - collect_start).count();
            }
            return result;
        }, profile, false, recovery==UATRecovery::Grouped);
}

QueryResult query_mat_block_index(
    const MATIndex& index, const TemporalBipartiteGraph& graph,
    const Query& query, QueryProfile* profile,
    const EncodedAttributes* attributes, QueryTimer* timer) {
    if (profile != nullptr) *profile = {};
    const auto& blocks = index.blocks();
    const auto& vertex_pool = index.vertex_pool();
    const auto& directory = index.query_directory();
    const auto& time_groups = index.query_time_groups();
    return query_uv_index(
        graph, query, index.stats().scope, [&](AttributeMask attrs) {
            const std::size_t vertex_count = graph.num_v_vertices();
            VertexCollector selected(vertex_count);
            std::vector<IndexPostingRange> ranges;
            const auto lookup_start = profile_start(profile);
            for (AttributeMask subset = attrs; subset != 0;
                 subset = (subset - 1) & attrs) {
                if (profile != nullptr) ++profile->directory_keys;
                const Count fixed_core = query.beta;
                const MATIndex::QueryKey key{
                    subset, fixed_core};
                const IndexPostingRange* found =
                    find_posting_range(directory, key);
                if (found == nullptr) continue;
                ranges.push_back(*found);
                if (profile != nullptr) ++profile->directory_hits;
            }
            if (profile != nullptr) {
                profile->directory_lookup_us +=
                    std::chrono::duration<double, std::micro>(
                        std::chrono::steady_clock::now() - lookup_start).count();
            }
            const auto collect_start = profile_start(profile);
            for (const IndexPostingRange range : ranges) {
                const auto groups_begin = time_groups.begin() +
                    static_cast<std::ptrdiff_t>(range.offset);
                const auto groups_end = groups_begin +
                    static_cast<std::ptrdiff_t>(range.count);
                const auto first_group = std::lower_bound(
                    groups_begin, groups_end, query.window.start,
                    [](const MATIndex::TimeGroup& group, Timestamp start) {
                        return group.start < start;
                    });
                const auto last_group = std::upper_bound(
                    first_group, groups_end, query.window.end,
                    [](Timestamp end, const MATIndex::TimeGroup& group) {
                        return end < group.start;
                    });
                const Count skyline_threshold = query.alpha;
                for (auto group = first_group;
                     group != last_group; ++group) {
                    if (group->end > query.window.end) continue;
                    const auto entries_begin = blocks.begin() +
                        static_cast<std::ptrdiff_t>(group->postings.offset);
                    const auto entries_end = entries_begin +
                        static_cast<std::ptrdiff_t>(group->postings.count);
                    const auto first_entry = std::lower_bound(
                        entries_begin, entries_end, skyline_threshold,
                        [&](const IndexBlock& block, Count threshold) {
                            const Count skyline = block.alpha;
                            return skyline < threshold;
                        });
                    if (profile != nullptr) {
                        profile->posting_entries += static_cast<std::size_t>(
                            std::distance(first_entry, entries_end));
                    }
                    for (auto entry = first_entry;
                         entry != entries_end; ++entry) {
                        const IndexBlock& block = *entry;
                        if (profile != nullptr) {
                            profile->vertex_refs += block.vertex_count;
                        }
                        const std::size_t vertex_end =
                            block.vertex_offset + block.vertex_count;
                        for (std::size_t i = block.vertex_offset;
                             i < vertex_end; ++i) {
                            selected.add(vertex_pool[i]);
                        }
                    }
                }
            }
            std::vector<VertexID> result = selected.finish();
            if (profile != nullptr) {
                profile->candidate_collect_us +=
                    std::chrono::duration<double, std::micro>(
                        std::chrono::steady_clock::now() - collect_start).count();
            }
            return result;
        }, profile, true, false, attributes, timer);
}

}
