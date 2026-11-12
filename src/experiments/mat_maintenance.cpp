#include "experiments/mat_maintenance.h"
#include "core/bicore_decomposition.h"
#include "core/dynamic_bicore.h"
#include "index/index_builder.h"
#include "index/index_io.h"
#include "query/mat_query.h"
#include "query/online_query.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <sys/resource.h>

namespace tacc {
namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;
namespace {
void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}
double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}
bool profile_enabled{}, component_pruning{}, window_pruning{}, compact_graph{};
Clock::time_point profile_start() { return profile_enabled ? Clock::now() : Clock::time_point{}; }
void profile_stop(double& total, Clock::time_point start) {
    if (profile_enabled) total += seconds(start);
}
std::uint64_t number(const char* value) {
    const std::string text(value);
    require(!text.empty() && text.find_first_not_of("0123456789") == std::string::npos,
            "expected unsigned integer");
    return std::stoull(text);
}
std::uint64_t bounded(std::mt19937_64& random, std::uint64_t limit) {
    require(limit != 0, "empty random range");
    const auto threshold = (std::uint64_t{0} - limit) % limit;
    std::uint64_t value;
    do { value = random(); } while (value < threshold);
    return value % limit;
}
AttributeMask universe(const TemporalBipartiteGraph& graph, bool u_side) {
    AttributeMask mask = 0;
    const auto size = u_side ? graph.num_u_vertices() : graph.num_v_vertices();
    for (VertexID id = 0; id < size; ++id) {
        const auto a = u_side ? graph.u_attribute(id) : graph.v_attribute(id);
        require(a < 64, "attribute mask overflow");
        mask |= AttributeMask{1} << a;
    }
    return mask;
}
std::vector<AttributeMask> subsets(AttributeMask all) {
    std::vector<AttributeMask> masks;
    for (auto mask = all; mask; mask = (mask - 1) & all) masks.push_back(mask);
    std::sort(masks.begin(), masks.end(), [](auto a, auto b) {
        return std::make_pair(__builtin_popcountll(a), a) <
               std::make_pair(__builtin_popcountll(b), b);
    });
    return masks;
}
std::vector<AttributeID> attributes(AttributeMask mask) {
    std::vector<AttributeID> values;
    for (AttributeID a = 0; a < 64; ++a)
        if (mask & (AttributeMask{1} << a)) values.push_back(a);
    return values;
}
void normalize(QueryResult& result) {
    std::sort(result.u_vertices.begin(), result.u_vertices.end());
    std::sort(result.v_vertices.begin(), result.v_vertices.end());
}
TemporalBipartiteGraph load(const fs::path& data, const fs::path& edges) {
    TemporalBipartiteGraph graph;
    require(graph.load_from_files((data / "graph_ulabel.txt").string(),
            (data / "graph_vlabel.txt").string(), edges.string()), "graph load failed");
    return graph;
}

struct PairEvents {
    VertexID u{}, v{};
    std::vector<Timestamp> times;
    std::vector<std::uint8_t> present;
};
struct Event {
    std::size_t pair{}, occurrence{};
    VertexID original_u{}, original_v{};
    Timestamp time{};
};
struct MutableEvents {
    const TemporalBipartiteGraph& labels;
    std::vector<PairEvents> pairs;
    std::vector<Event> events;
    std::size_t active{};
    explicit MutableEvents(const TemporalBipartiteGraph& graph) : labels(graph) {
        for (std::size_t id = 0; id < graph.num_unique_edges(); ++id) {
            PairEvents pair{graph.unique_edge_u(id), graph.unique_edge_v(id), {}, {}};
            for (std::size_t i = 0; i < graph.unique_edge_timestamp_count(id); ++i) {
                const auto t = graph.unique_edge_timestamp(id, i);
                pair.times.push_back(t);
                pair.present.push_back(1);
                events.push_back({pairs.size(), i, graph.original_u_id(pair.u), graph.original_v_id(pair.v), t});
            }
            pairs.push_back(std::move(pair));
        }
        std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
            return std::tie(a.original_u, a.original_v, a.time) < std::tie(b.original_u, b.original_v, b.time);
        });
        active = events.size();
    }
    void apply(const Event& event, bool insertion) {
        auto& flag = pairs.at(event.pair).present.at(event.occurrence);
        require(bool(flag) != insertion, "duplicate insertion or missing deletion");
        flag = insertion;
        if (insertion) ++active; else --active;
    }
    BipartiteSubgraph snapshot(AttributeMask mask, const std::vector<std::uint32_t>& counts) const {
        BipartiteSubgraph graph(labels.num_u_vertices(), labels.num_v_vertices());
        for (VertexID u = 0; u < graph.u_selected.size(); ++u)
            graph.u_selected[u] = (mask & (AttributeMask{1} << labels.u_attribute(u))) != 0;
        std::fill(graph.v_selected.begin(), graph.v_selected.end(), 1);
        for (std::size_t i = 0; i < pairs.size(); ++i) if (counts[i]) {
            const auto& pair = pairs[i];
            require(graph.u_selected[pair.u], "projection contains an unselected U vertex");
            graph.u_adj[pair.u].push_back(pair.v);
            graph.v_adj[pair.v].push_back(pair.u);
        }
        return graph;
    }
    void write(const fs::path& file) const {
        std::ofstream out(file);
        require(bool(out), "cannot write graph checkpoint");
        for (const auto& event : events)
            if (pairs[event.pair].present[event.occurrence])
                out << event.original_u << ' ' << event.original_v << ' ' << event.time << '\n';
        out.close();
        require(bool(out), "graph checkpoint write failed");
    }
};

struct Work {
    std::size_t masks{}, windows{}, decompositions{}, bucket_updates{};
    std::size_t surviving_memberships{}, removed_memberships{}, emitted{};
    std::size_t region_vertices{}, projected_pairs{}, component_mark_bytes{};
    std::size_t skipped_empty_starts{}, core_early_stops{}, skipped_shorter_windows{}, compact_vertices{};
    double component_prepare_s{}, retain_s{}, recompute_s{}, finalize_s{};
    double projection_s{}, state_copy_s{}, snapshot_s{}, initialize_s{}, bucket_prepare_s{}, dynamic_s{}, emit_s{}, release_s{};
};
struct Projection {
    const MutableEvents& source;
    AttributeMask mask;
    std::vector<std::pair<VertexID, VertexID>> local_pairs;
    std::vector<VertexID> global_v;
    std::vector<std::vector<std::size_t>> by_time;
    std::vector<std::uint32_t> suffix_counts;
    std::size_t num_u{};

    Projection(const MutableEvents& graph, AttributeMask attributes,
               const std::vector<std::uint8_t>* region, Work& work)
        : source(graph), mask(attributes), by_time(graph.labels.time_points().size()) {
        const auto invalid = std::numeric_limits<VertexID>::max();
        std::vector<VertexID> u_map, v_map;
        if (compact_graph) {
            u_map.assign(graph.labels.num_u_vertices(), invalid);
            v_map.assign(graph.labels.num_v_vertices(), invalid);
        } else suffix_counts.resize(graph.pairs.size());
        const auto& times = graph.labels.time_points();
        for (std::size_t id = 0; id < graph.pairs.size(); ++id) {
            const auto& pair = graph.pairs[id];
            if (!(mask & (AttributeMask{1} << graph.labels.u_attribute(pair.u)))) continue;
            if (region && !region->at(pair.v)) continue;
            ++work.projected_pairs;
            const auto active = std::count(pair.present.begin(), pair.present.end(), std::uint8_t{1});
            if (!active) continue;
            auto slot = id;
            if (compact_graph) {
                if (u_map[pair.u] == invalid) u_map[pair.u] = num_u++;
                if (v_map[pair.v] == invalid) {
                    v_map[pair.v] = global_v.size();
                    global_v.push_back(pair.v);
                }
                slot = local_pairs.size();
                local_pairs.emplace_back(u_map[pair.u], v_map[pair.v]);
                suffix_counts.push_back(static_cast<std::uint32_t>(active));
            } else suffix_counts[slot] = static_cast<std::uint32_t>(active);
            for (std::size_t i = 0; i < pair.times.size(); ++i) if (pair.present[i]) {
                const auto ti = std::lower_bound(times.begin(), times.end(), pair.times[i]) - times.begin();
                by_time[ti].push_back(slot);
            }
        }
        if (compact_graph) work.compact_vertices += num_u + global_v.size();
    }
    std::pair<VertexID, VertexID> endpoints(std::size_t slot) const {
        if (compact_graph) return local_pairs.at(slot);
        const auto& pair = source.pairs.at(slot);
        return {pair.u, pair.v};
    }
    VertexID original_v(VertexID local) const {
        return compact_graph ? global_v.at(local) : local;
    }
    BipartiteSubgraph snapshot(const std::vector<std::uint32_t>& counts) const {
        if (!compact_graph) return source.snapshot(mask, counts);
        BipartiteSubgraph result(num_u, global_v.size());
        std::fill(result.u_selected.begin(), result.u_selected.end(), 1);
        std::fill(result.v_selected.begin(), result.v_selected.end(), 1);
        for (std::size_t id = 0; id < local_pairs.size(); ++id) if (counts[id]) {
            const auto [u, v] = local_pairs[id];
            result.u_adj[u].push_back(v);
            result.v_adj[v].push_back(u);
        }
        return result;
    }
};
std::map<AttributeMask, std::vector<std::uint8_t>> regions(
    const MutableEvents& graph, const std::vector<Event>& changes, const std::vector<AttributeMask>& masks,
    AttributeMask bit, Work& work) {
    std::map<AttributeMask, std::vector<std::uint8_t>> result;
    const auto nu = graph.labels.num_u_vertices(), nv = graph.labels.num_v_vertices();
    std::vector<std::uint8_t> active(graph.pairs.size());
    for (const auto& event : changes) active[event.pair] = 1;
    for (std::size_t id = 0; id < graph.pairs.size(); ++id)
        active[id] = active[id] || std::any_of(graph.pairs[id].present.begin(), graph.pairs[id].present.end(),
                                                   [](auto flag) { return flag != 0; });
    for (auto mask : masks) if (mask & bit) {
        std::vector<std::size_t> parent(nu + nv), sizes(nu + nv, 1);
        std::iota(parent.begin(), parent.end(), 0);
        const auto root = [&](std::size_t node) {
            while (parent[node] != node) { parent[node] = parent[parent[node]]; node = parent[node]; }
            return node;
        };
        for (std::size_t id = 0; id < graph.pairs.size(); ++id) {
            const auto& pair = graph.pairs[id];
            if (!active[id] || !(mask & (AttributeMask{1} << graph.labels.u_attribute(pair.u)))) continue;
            auto a = root(pair.u), b = root(nu + pair.v);
            if (a == b) continue;
            if (sizes[a] < sizes[b]) std::swap(a, b);
            parent[b] = a; sizes[a] += sizes[b];
        }
        std::vector<std::uint8_t> changed_roots(nu + nv);
        for (const auto& event : changes) {
            const auto u = graph.pairs[event.pair].u;
            if (!(mask & (AttributeMask{1} << graph.labels.u_attribute(u)))) continue;
            const auto changed_root = root(u);
            if (!changed_roots[changed_root]) work.region_vertices += sizes[changed_root];
            changed_roots[changed_root] = 1;
        }
        auto& marked = result[mask];
        marked.resize(nv);
        for (VertexID v = 0; v < nv; ++v) marked[v] = changed_roots[root(nu + v)];
        work.component_mark_bytes += marked.size();
    }
    return result;
}
MATIndex refresh(const MATIndex& old, const MutableEvents& graph, const std::vector<Event>& changes,
                 const std::vector<AttributeMask>& masks, Work& work) {
    MATIndex next;
    next.mutable_stats().scope = IndexScope::UV;
    next.mutable_stats().method = BuildMethod::Opt;
    auto start = Clock::now();
    AttributeMask bit = 0;
    for (const auto& event : changes)
        bit |= AttributeMask{1} << graph.labels.u_attribute(graph.pairs[event.pair].u);
    std::map<AttributeMask, std::vector<Timestamp>> affected_times;
    for (auto mask : masks) if (mask & bit) {
        auto& changed_times = affected_times[mask];
        for (const auto& event : changes)
            if (mask & (AttributeMask{1} << graph.labels.u_attribute(graph.pairs[event.pair].u)))
                changed_times.push_back(event.time);
        std::sort(changed_times.begin(), changed_times.end());
        changed_times.erase(std::unique(changed_times.begin(), changed_times.end()), changed_times.end());
    }
    const auto affected_regions = component_pruning ? regions(graph, changes, masks, bit, work)
        : std::map<AttributeMask, std::vector<std::uint8_t>>{};
    work.component_prepare_s = seconds(start);
    start = Clock::now();
    for (const auto& block : old.blocks()) {
        bool affected = false;
        const auto domain = affected_times.find(block.attributes);
        if (domain != affected_times.end()) {
            const auto first = std::lower_bound(domain->second.begin(), domain->second.end(), block.start);
            affected = first != domain->second.end() && *first <= block.end;
        }
        for (std::size_t i = 0; i < block.vertex_count; ++i) {
            const auto vertex = old.vertex_pool().at(block.vertex_offset + i);
            const bool in_region = !component_pruning ||
                (affected && affected_regions.at(block.attributes).at(vertex));
            if (affected && in_region) { ++work.removed_memberships; continue; }
            next.add({vertex, block.attributes,
                      block.start, block.end, block.alpha, block.beta});
            ++work.surviving_memberships;
        }
    }
    work.retain_s = seconds(start);
    start = Clock::now();
    const auto& times = graph.labels.time_points();
    for (auto mask : masks) if (mask & bit) {
        ++work.masks;
        auto phase = profile_start();
        std::vector<std::size_t> changed_indices;
        for (auto t : affected_times.at(mask))
            changed_indices.push_back(std::lower_bound(times.begin(), times.end(), t) - times.begin());
        Projection projection(graph, mask, component_pruning ? &affected_regions.at(mask) : nullptr, work);
        const auto& by_time = projection.by_time;
        auto& suffix_counts = projection.suffix_counts;
        profile_stop(work.projection_s, phase);
        for (std::size_t s = 0; s <= changed_indices.back(); ++s) {
            const auto t_index = *std::lower_bound(changed_indices.begin(), changed_indices.end(), s);
            if (window_pruning && s < t_index && by_time[s].empty()) {
                ++work.skipped_empty_starts;
                continue;
            }
            phase = profile_start();
            auto counts = suffix_counts;
            profile_stop(work.state_copy_s, phase);
            phase = profile_start();
            auto snapshot = projection.snapshot(counts);
            profile_stop(work.snapshot_s, phase);
            phase = profile_start();
            bool stop_starts = false;
            {
            DynamicBicoreMaintainer maintainer(std::move(snapshot));
            profile_stop(work.initialize_s, phase);
            ++work.decompositions;
            stop_starts = window_pruning && maintainer.membership_count() == 0;
            for (std::size_t e = times.size(); e-- > t_index;) {
                if (window_pruning && maintainer.membership_count() == 0) {
                    ++work.core_early_stops;
                    work.skipped_shorter_windows += e - t_index + 1;
                    break;
                }
                ++work.windows;
                phase = profile_start();
                std::vector<std::pair<VertexID, VertexID>> removed;
                for (auto id : by_time[e]) if (--counts[id] == 0)
                    removed.push_back(projection.endpoints(id));
                profile_stop(work.bucket_prepare_s, phase);
                if (removed.empty()) continue;
                phase = profile_start();
                const auto update = maintainer.delete_edge_bucket(removed);
                profile_stop(work.dynamic_s, phase);
                ++work.bucket_updates;
                phase = profile_start();
                for (const auto& change : update.v_frontier_changes) {
                    if (!(change.old_alpha_max == 1 && change.beta == 1)) {
                        next.add({projection.original_v(change.vertex), mask, times[s], times[e], change.old_alpha_max, change.beta});
                        ++work.emitted;
                    }
                }
                profile_stop(work.emit_s, phase);
            }
            phase = profile_start();
            }
            profile_stop(work.release_s, phase);
            if (stop_starts) break;
            for (auto id : by_time[s]) --suffix_counts[id];
        }
    }
    work.recompute_s = seconds(start);
    start = Clock::now();
    next.finalize();
    work.finalize_s = seconds(start);
    return next;
}
void equal_index(const MATIndex& a, const MATIndex& b) {
    require(a.blocks().size() == b.blocks().size(), "block counts differ from full rebuild");
    for (std::size_t i = 0; i < a.blocks().size(); ++i) {
        const auto& x = a.blocks()[i]; const auto& y = b.blocks()[i];
        require(std::tie(x.attributes, x.start, x.end, x.alpha, x.beta, x.vertex_count) ==
                std::tie(y.attributes, y.start, y.end, y.alpha, y.beta, y.vertex_count), "block conditions differ");
        for (std::size_t j = 0; j < x.vertex_count; ++j)
            require(a.vertex_pool().at(x.vertex_offset + j) == b.vertex_pool().at(y.vertex_offset + j),
                    "block members differ");
    }
}
std::vector<Query> queries(const TemporalBipartiteGraph& graph, std::uint64_t seed) {
    std::mt19937_64 random(seed ^ 0x9e3779b97f4a7c15ULL);
    const auto us = subsets(universe(graph, true)), vs = subsets(universe(graph, false));
    const auto& times = graph.time_points();
    std::size_t max_a = 1, max_b = 1;
    for (VertexID u = 0; u < graph.num_u_vertices(); ++u) max_a = std::max(max_a, graph.static_u_degree(u));
    for (VertexID v = 0; v < graph.num_v_vertices(); ++v) max_b = std::max(max_b, graph.static_v_degree(v));
    std::vector<Query> result;
    for (std::size_t i = 0; i < 500; ++i) {
        auto s = bounded(random, times.size()), e = bounded(random, times.size());
        if (s > e) std::swap(s, e);
        const auto um = i % 10 ? us[bounded(random, us.size())] : universe(graph, true);
        const auto vm = i % 10 ? vs[bounded(random, vs.size())] : universe(graph, false);
        Count a = 1 + bounded(random, i % 4 ? std::min<std::size_t>(6, max_a + 1) : max_a + 1);
        Count b = 1 + bounded(random, i % 4 ? std::min<std::size_t>(6, max_b + 1) : max_b + 1);
        if (a == 1 && b == 1) a = 2;
        if (i % 10 == 0) { s = 0; e = times.size() - 1; a = 2; b = 1; }
        result.push_back({{times[s], times[e]}, attributes(um), attributes(vm), a, b});
    }
    return result;
}
void check_answers(const MATIndex& index, const TemporalBipartiteGraph& graph,
                   const std::vector<Query>& conditions) {
    for (const auto& q : conditions) {
        auto actual = mat_query(index, graph, q, nullptr), expected = online_query(graph, q, nullptr);
        normalize(actual); normalize(expected);
        require(actual.u_vertices == expected.u_vertices && actual.v_vertices == expected.v_vertices,
                "MAT/Online full U/V mismatch");
    }
}
struct Result { std::string operation; std::size_t step{}, events{}; double maintenance{}, rebuild{}; };
}

MATMaintenanceOptions configure_mat_maintenance(std::size_t updates, std::uint64_t seed) {
        MATMaintenanceOptions options;
        options.updates = updates;
        options.seed = seed;
        const auto batch_size = std::getenv("TACC_MAINT_BATCH_SIZE") ? number(std::getenv("TACC_MAINT_BATCH_SIZE")) : 1;
        require(updates > 0, "update count must be positive");
        require(batch_size > 0 && batch_size <= updates, "batch size must be within [1, COUNT]");
        options.batch_size = static_cast<std::size_t>(batch_size);
        const auto enabled = [](const char* name) {
            const auto* value = std::getenv(name);
            return value && std::string(value) == "1";
        };
        options.profile = enabled("TACC_MAINT_PROFILE");
        options.component_pruning = enabled("TACC_MAINT_COMPONENT_PRUNING");
        options.window_pruning = enabled("TACC_MAINT_WINDOW_PRUNING");
        options.compact_graph = enabled("TACC_MAINT_COMPACT_GRAPH");
        for (const auto* key : {"TAC_QUERY_PROFILE", "TAC_PROFILE_DYNAMIC", "TAC_PROFILE_MAT_STREAMING",
                                "TAC_VERIFY_DYNAMIC_BUCKET", "TAC_ATTRIBUTE_PRUNING"})
            require(setenv(key, "0", 1) == 0, "setenv failed");
        for (const auto* key : {"TAC_DYNAMIC_LAYERED", "TAC_MAT_BOUNDARY_ONLY", "TAC_MAT_ORDERED_STREAMING"})
            require(setenv(key, "1", 1) == 0, "setenv failed");
        require(setenv("TAC_INDEX_SCOPE", "uv", 1) == 0 && setenv("TAC_QUERY_ROUTE", "uv", 1) == 0 &&
                setenv("TAC_MAT_RECOVERY", "grouped", 1) == 0, "setenv failed");
        return options;
}

void run_mat_maintenance(const TemporalBipartiteGraph& base, MATIndex index,
                         const MATMaintenanceOptions& options, std::ostream& log) {
        const auto& data = options.data;
        const auto& output = options.output;
        const auto count = options.updates;
        const auto seed = options.seed;
        const auto batch_size = options.batch_size;
        profile_enabled = options.profile;
        component_pruning = options.component_pruning;
        window_pruning = options.window_pruning;
        compact_graph = options.compact_graph;
        require(count > 0 && batch_size > 0 && batch_size <= count, "invalid maintenance update/batch count");
        require(!base.time_points().empty() && base.num_attributes() <= 64, "unsupported graph");
        require(base.num_edges() == base.num_unique_temporal_edges(),
                "maintenance requires deduplicated (u,v,t) input; do not silently discard duplicate events");
        require(index.stats().method == BuildMethod::Opt && index.stats().scope == IndexScope::UV,
                "initial index must be opt/UV");
        MutableEvents events(base);
        require(count < events.events.size(), "update count must leave at least one temporal event");
        require(!options.prefix.empty(), "maintenance output prefix is empty");
        for (const auto* directory : {"updates", "graphs", "indexes"}) fs::create_directories(output / directory);
        const auto trace_file = output / "updates" / (options.prefix + "_updates.tsv");
        const auto working_edges = output / "graphs" / (options.prefix + "_working_graph_timestamp.txt");
        const auto deleted_edges = output / "graphs" / (options.prefix + "_deleted_graph_timestamp.txt");
        for (const auto& path : {trace_file, working_edges, deleted_edges,
                output / "indexes" / (options.prefix + "_deleted_mat_opt_uv_index.bin"),
                output / "indexes" / (options.prefix + "_restored_mat_opt_uv_index.bin")})
            require(!fs::exists(path) && !fs::is_symlink(path), "refusing to overwrite maintenance output: " + path.string());
        log << "MAINTENANCE_CONFIG method=affected_window_dynamic_recompute fine_grained_incremental=0 count=" << count
            << " detailed_profile=" << profile_enabled << " component_pruning=" << component_pruning
            << " window_pruning=" << window_pruning << " compact_graph=" << compact_graph
            << " batch_size=" << batch_size << " queryable_state=batch_end_only"
            << " seed=" << seed << " sampling=uniform_distinct_existing_temporal_events data=" << std::quoted(data.string()) << '\n';
        log << "TIMING maintenance=event_apply+working_graph_write_load+affected_window_recompute+global_recompression_directory"
            << " rebuild=same_event_apply_and_graph_refresh+full_opt_build validation=excluded checkpoint_save=excluded"
            << " warmup=0 order=maintenance_then_rebuild process_peak_includes_validation=1\n";
        log << "BATCH_TIMING step_times=whole_batch summary_avg_ms=amortized_per_event"
            << " summary_avg_batch_ms=whole_batch not_individual_event_latency=1"
            << " intermediate_states_within_batch_validated=0 step_event_fields=first_event\n";
        auto start = Clock::now();
        { auto oracle = build_mat_index(base, BuildMethod::Opt); equal_index(index, oracle); }
        const auto initial_hash = formal_index_hash(index);
        const auto conditions = queries(base, seed);
        const auto masks = subsets(universe(base, true));
        std::vector<std::size_t> order(events.events.size());
        std::iota(order.begin(), order.end(), 0);
        std::mt19937_64 random(seed);
        for (std::size_t i = 0; i < count; ++i)
            std::swap(order[i], order[i + bounded(random, order.size() - i)]);
        order.resize(count);
        std::ofstream trace(trace_file);
        require(bool(trace), "cannot open update trace");
        trace << "# seed=" << seed << " deletion_order=ascending_id insertion_order=descending_id\n"
              << "id\tu\tv\tt\n";
        for (std::size_t i = 0; i < order.size(); ++i) {
            const auto& event = events.events[order[i]];
            trace << i + 1 << '\t' << event.original_u << '\t' << event.original_v << '\t' << event.time << '\n';
        }
        trace.close();
        require(bool(trace), "cannot save update trace");
        std::vector<Result> results;
        log << "ARTIFACT updates=" << std::quoted(trace_file.string())
            << " working_graph=" << std::quoted(working_edges.string()) << '\n';
        const auto full_windows = masks.size() * base.time_points().size() * (base.time_points().size() + 1) / 2;
        for (std::size_t operation = 0; operation < 2; ++operation) {
            const bool insertion = operation == 1;
            const std::string name = insertion ? "insert" : "delete";
            for (std::size_t step = 0; step < count; step += batch_size) {
                const auto end = std::min<std::size_t>(count, step + batch_size);
                const auto event_id = insertion ? count - 1 - step : step;
                const auto& event = events.events[order[event_id]];
                start = Clock::now();
                std::vector<Event> changes;
                changes.reserve(end - step);
                for (auto offset = step; offset < end; ++offset) {
                    const auto id = insertion ? count - 1 - offset : offset;
                    changes.push_back(events.events[order[id]]);
                    events.apply(changes.back(), insertion);
                }
                const double apply_s = seconds(start);
                start = Clock::now();
                events.write(working_edges);
                auto current_graph = load(data, working_edges);
                const double graph_refresh_s = seconds(start);
                require(current_graph.num_unique_temporal_edges() == events.active, "event count mismatch");
                require(current_graph.num_u_vertices() == base.num_u_vertices() &&
                        current_graph.num_v_vertices() == base.num_v_vertices(), "vertex mapping changed");
                Work work;
                start = Clock::now();
                index = refresh(index, events, changes, masks, work);
                const double patch_s = seconds(start);
                const double maintenance_s = apply_s + graph_refresh_s + patch_s;
                start = Clock::now();
                auto rebuilt = build_mat_index(current_graph, BuildMethod::Opt);
                const double full_build_s = seconds(start);
                const double rebuild_s = apply_s + graph_refresh_s + full_build_s;
                equal_index(index, rebuilt);
                check_answers(index, current_graph, conditions);
                const auto batch_id = step / batch_size + 1;
                results.push_back({name, batch_id, changes.size(), maintenance_s, rebuild_s});
                log << "STEP operation=" << name << " step=" << batch_id << " trace_id=" << event_id + 1
                    << " batch_events=" << changes.size() << " events_completed=" << end
                    << " u=" << event.original_u << " v=" << event.original_v << " t=" << event.time
                    << " active_events=" << events.active << " affected_masks=" << work.masks
                    << " total_masks=" << masks.size() << " recomputed_windows=" << work.windows
                    << " decompositions=" << work.decompositions << " dynamic_bucket_updates=" << work.bucket_updates
                    << " full_windows=" << full_windows << " retained_memberships=" << work.surviving_memberships
                    << " removed_memberships=" << work.removed_memberships << " candidates=" << work.emitted
                    << " event_apply_s=" << apply_s << " graph_refresh_s=" << graph_refresh_s
                    << " component_prepare_s=" << work.component_prepare_s
                    << " retain_s=" << work.retain_s << " recompute_s=" << work.recompute_s
                    << " projection_s=" << work.projection_s << " state_copy_s=" << work.state_copy_s
                    << " snapshot_s=" << work.snapshot_s << " initialize_s=" << work.initialize_s
                    << " bucket_prepare_s=" << work.bucket_prepare_s << " dynamic_s=" << work.dynamic_s
                    << " emit_s=" << work.emit_s << " release_s=" << work.release_s
                    << " projected_pairs=" << work.projected_pairs << " region_vertices=" << work.region_vertices
                    << " component_mark_bytes=" << work.component_mark_bytes
                    << " skipped_empty_starts=" << work.skipped_empty_starts
                    << " core_early_stops=" << work.core_early_stops
                    << " skipped_shorter_windows=" << work.skipped_shorter_windows
                    << " compact_vertices=" << work.compact_vertices
                    << " finalize_s=" << work.finalize_s << " index_update_s=" << patch_s
                    << " maintenance_s=" << maintenance_s << " full_opt_build_s=" << full_build_s
                    << " rebuild_s=" << rebuild_s << " speedup=" << rebuild_s / maintenance_s
                    << " block_count=" << index.blocks().size() << " full_index_equal=1 query_pass=500 status=ok\n" << std::flush;
                std::cout << "PROGRESS operation=" << name << " completed=" << end << '/' << count
                          << " batch_events=" << changes.size()
                          << " maintenance_s=" << maintenance_s << " rebuild_s=" << rebuild_s << " status=ok\n" << std::flush;
                if (end == count) {
                    const auto prefix = insertion ? "restored" : "deleted";
                    const auto checkpoint = output / "indexes" / (options.prefix + '_' + prefix + "_mat_opt_uv_index.bin");
                    std::string error;
                    require(save_index_binary(checkpoint.string(), current_graph, index, &error),
                            "checkpoint save failed: " + error);
                    if (!insertion) {
                        fs::copy_file(working_edges, deleted_edges);
                        MATIndex reloaded;
                        require(load_index_binary(checkpoint.string(), current_graph, reloaded, &error),
                                "deleted checkpoint load failed: " + error);
                        equal_index(index, reloaded);
                        index = std::move(reloaded);
                        log << "INSERT_INITIAL source=deleted_checkpoint round_trip_equal=1\n";
                    } else {
                        require(events.active == events.events.size() && formal_index_hash(index) == initial_hash,
                                "restore did not recover initial index");
                        require(graph_fingerprint(current_graph) == graph_fingerprint(base),
                                "restored graph fingerprint differs");
                    }
                }
            }
        }
        log << "\n# SUMMARY time_unit=ms speedup=sum_rebuild/sum_maintenance\n"
            << "Operation Batch Events MaintenanceBatch(ms) RebuildBatch(ms) Speedup\n";
        for (const auto& result : results)
            log << result.operation << ' ' << result.step << ' ' << result.events << ' ' << result.maintenance * 1000 << ' '
                << result.rebuild * 1000 << ' ' << result.rebuild / result.maintenance << '\n';
        std::ostringstream table;
        table << std::fixed << std::setprecision(3) << std::left
            << "\n# MAINTENANCE_TABLE time_unit=ms average_scope=per_event batch_size=" << batch_size << '\n'
            << "# Batch>1: per-event amortized cost, not individual-event latency; validation excluded from timings.\n"
            << std::setw(12) << "Operation" << std::setw(10) << "Updates" << std::setw(10) << "Batches"
            << std::setw(20) << "Maintenance(ms)" << std::setw(20) << "Rebuild(ms)" << "Speedup\n";
        for (const std::string name : {"delete", "insert"}) {
            double maintenance_s = 0, rebuild_s = 0;
            std::size_t batches = 0;
            for (const auto& result : results) if (result.operation == name) {
                maintenance_s += result.maintenance; rebuild_s += result.rebuild;
                ++batches;
            }
            log << "SUMMARY operation=" << name << " updates=" << count << " batches=" << batches
                << " maintenance_total_s=" << maintenance_s << " rebuild_total_s=" << rebuild_s
                << " maintenance_avg_ms=" << maintenance_s * 1000 / count
                << " rebuild_avg_ms=" << rebuild_s * 1000 / count
                << " maintenance_avg_batch_ms=" << maintenance_s * 1000 / batches
                << " rebuild_avg_batch_ms=" << rebuild_s * 1000 / batches
                << " speedup=" << rebuild_s / maintenance_s << '\n';
            table << std::setw(12) << name << std::setw(10) << count << std::setw(10) << batches
                << std::setw(20) << maintenance_s * 1000 / count
                << std::setw(20) << rebuild_s * 1000 / count << rebuild_s / maintenance_s << "x\n";
        }
        log << table.str();
        rusage usage{};
        require(getrusage(RUSAGE_SELF, &usage) == 0, "getrusage failed");
        log << "RESULT exp=exp4 method=mat status=ok updates=" << count * 2 << " full_index_checks=" << results.size()
            << " query_checks=" << results.size() * conditions.size() << " restored_initial_index=1"
            << " process_peak_mb=" << usage.ru_maxrss / 1024.0 << " peak_includes_oracles=1\n";
        require(bool(log), "cannot write maintenance result log");
}

}
