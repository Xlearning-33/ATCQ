#include "common/graph.h"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <utility>

namespace tacc {

namespace {

using AttributeRows = std::map<VertexID, std::string>;

bool read_attribute_file(const std::string& path,
                         const char* side,
                         AttributeRows& rows) {
    std::ifstream input(path);
    if (!input) {
        std::cerr << "Cannot open " << side << " attribute file: " << path << '\n';
        return false;
    }

    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#' || line[0] == '%') {
            continue;
        }

        std::istringstream iss(line);
        VertexID original_id{};
        std::string attribute;
        if (!(iss >> original_id >> attribute)) {
            std::cerr << "Invalid " << side << " attribute row at " << path
                      << ':' << line_number << '\n';
            return false;
        }
        if (original_id == 0) {
            std::cerr << side << " vertex IDs must start from 1 at " << path
                      << ':' << line_number << '\n';
            return false;
        }
        if (!rows.emplace(original_id, std::move(attribute)).second) {
            std::cerr << "Duplicate " << side << " vertex ID at " << path
                      << ':' << line_number << '\n';
            return false;
        }
    }

    if (rows.empty()) {
        std::cerr << "No attributes found in: " << path << '\n';
        return false;
    }
    return true;
}

}

bool TemporalBipartiteGraph::load_from_files(
    const std::string& ulabel_file,
    const std::string& vlabel_file,
    const std::string& edge_file) {
    std::atomic_store(&fingerprint_, std::shared_ptr<const std::uint64_t>{});
    num_u_vertices_ = 0;
    num_v_vertices_ = 0;
    min_timestamp_ = 0;
    max_timestamp_ = 0;
    edges_.clear();
    std::atomic_store(&grouped_v_time_, std::shared_ptr<const GroupedVTimeAdjacency>{});
    time_points_.clear();
    unique_edges_.clear();
    unique_edge_timestamps_.clear();
    u_edge_offsets_.clear();
    v_edge_offsets_.clear();
    u_incident_edges_.clear();
    v_incident_edges_.clear();
    u_attribute_incident_work_.clear();
    v_attribute_incident_work_.clear();
    u_attributes_.clear();
    v_attributes_.clear();
    original_u_ids_.clear();
    original_v_ids_.clear();
    u_id_to_internal_.clear();
    v_id_to_internal_.clear();
    attribute_name_to_id_.clear();
    attribute_names_.clear();

    AttributeRows u_rows;
    AttributeRows v_rows;
    if (!read_attribute_file(ulabel_file, "U", u_rows) ||
        !read_attribute_file(vlabel_file, "V", v_rows)) {
        return false;
    }

    auto get_attribute_id = [this](const std::string& name) {
        auto found = attribute_name_to_id_.find(name);
        if (found != attribute_name_to_id_.end()) {
            return found->second;
        }
        const AttributeID id = static_cast<AttributeID>(attribute_names_.size());
        attribute_names_.push_back(name);
        attribute_name_to_id_.emplace(name, id);
        return id;
    };

    for (const auto& [original_id, attribute] : u_rows) {
        const VertexID internal_id = static_cast<VertexID>(original_u_ids_.size());
        u_id_to_internal_.emplace(original_id, internal_id);
        original_u_ids_.push_back(original_id);
        u_attributes_.push_back(get_attribute_id(attribute));
    }
    for (const auto& [original_id, attribute] : v_rows) {
        const VertexID internal_id = static_cast<VertexID>(original_v_ids_.size());
        v_id_to_internal_.emplace(original_id, internal_id);
        original_v_ids_.push_back(original_id);
        v_attributes_.push_back(get_attribute_id(attribute));
    }
    num_u_vertices_ = original_u_ids_.size();
    num_v_vertices_ = original_v_ids_.size();

    std::ifstream input(edge_file);
    if (!input) {
        std::cerr << "Cannot open temporal edge file: " << edge_file << '\n';
        return false;
    }

    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#' || line[0] == '%') {
            continue;
        }

        std::istringstream iss(line);
        VertexID original_u{};
        VertexID original_v{};
        Timestamp timestamp{};
        if (!(iss >> original_u >> original_v >> timestamp)) {
            std::cerr << "Invalid temporal edge row at " << edge_file
                      << ':' << line_number << '\n';
            return false;
        }

        const auto u = u_id_to_internal_.find(original_u);
        const auto v = v_id_to_internal_.find(original_v);
        if (u == u_id_to_internal_.end() || v == v_id_to_internal_.end()) {
            std::cerr << "Unknown endpoint at " << edge_file << ':' << line_number
                      << " (U=" << original_u << ", V=" << original_v << ")\n";
            return false;
        }
        edges_.push_back({u->second, v->second, timestamp});
    }

    if (edges_.empty()) {
        std::cerr << "No temporal edges found in: " << edge_file << '\n';
        return false;
    }

    std::sort(edges_.begin(), edges_.end(),
              [](const TemporalEdge& lhs, const TemporalEdge& rhs) {
                  if (lhs.timestamp != rhs.timestamp) {
                      return lhs.timestamp < rhs.timestamp;
                  }
                  if (lhs.u != rhs.u) return lhs.u < rhs.u;
                  return lhs.v < rhs.v;
              });
    min_timestamp_ = edges_.front().timestamp;
    max_timestamp_ = edges_.back().timestamp;
    time_points_.reserve(edges_.size());
    for (const TemporalEdge& edge : edges_) {
        if (time_points_.empty() || time_points_.back() != edge.timestamp) {
            time_points_.push_back(edge.timestamp);
        }
    }
    return rebuild_temporal_adjacency();
}

bool TemporalBipartiteGraph::rebuild_temporal_adjacency() {
    if (edges_.size() >
        std::numeric_limits<UniqueEdgeID>::max()) {
        std::cerr << "Too many temporal events for 32-bit sort IDs\n";
        return false;
    }
    std::vector<UniqueEdgeID> by_pair(edges_.size());
    std::iota(by_pair.begin(), by_pair.end(), UniqueEdgeID{0});
    std::sort(by_pair.begin(), by_pair.end(),
              [&](UniqueEdgeID lhs_id, UniqueEdgeID rhs_id) {
                  const TemporalEdge& lhs = edges_[lhs_id];
                  const TemporalEdge& rhs = edges_[rhs_id];
                  if (lhs.u != rhs.u) return lhs.u < rhs.u;
                  if (lhs.v != rhs.v) return lhs.v < rhs.v;
                  return lhs.timestamp < rhs.timestamp;
              });

    unique_edges_.clear();
    unique_edge_timestamps_.clear();
    unique_edges_.reserve(by_pair.size());
    unique_edge_timestamps_.reserve(by_pair.size());
    for (std::size_t begin = 0; begin < by_pair.size();) {
        const TemporalEdge& first = edges_[by_pair[begin]];
        std::size_t end = begin + 1;
        while (end < by_pair.size() &&
               edges_[by_pair[end]].u == first.u &&
               edges_[by_pair[end]].v == first.v) {
            ++end;
        }
        if (unique_edges_.size() >=
            std::numeric_limits<UniqueEdgeID>::max()) {
            std::cerr << "Too many unique edges for 32-bit edge IDs\n";
            return false;
        }
        const std::size_t timestamp_offset =
            unique_edge_timestamps_.size();
        Timestamp previous = 0;
        bool has_previous = false;
        for (std::size_t i = begin; i < end; ++i) {
            const Timestamp timestamp = edges_[by_pair[i]].timestamp;
            if (has_previous &&
                timestamp == previous) {
                continue;
            }
            previous = timestamp;
            has_previous = true;
            unique_edge_timestamps_.push_back(previous);
        }
        const std::size_t timestamp_count =
            unique_edge_timestamps_.size() - timestamp_offset;
        if (timestamp_count >
            std::numeric_limits<std::uint32_t>::max()) {
            std::cerr << "Too many timestamps on one unique edge\n";
            return false;
        }
        unique_edges_.push_back(
            {first.u, first.v, timestamp_offset,
             static_cast<std::uint32_t>(timestamp_count)});
        begin = end;
    }
    std::vector<UniqueEdgeID>().swap(by_pair);

    u_edge_offsets_.assign(num_u_vertices_ + 1, 0);
    v_edge_offsets_.assign(num_v_vertices_ + 1, 0);
    for (const UniqueTemporalEdge& edge : unique_edges_) {
        ++u_edge_offsets_[static_cast<std::size_t>(edge.u) + 1];
        ++v_edge_offsets_[static_cast<std::size_t>(edge.v) + 1];
    }
    std::partial_sum(u_edge_offsets_.begin(), u_edge_offsets_.end(),
                     u_edge_offsets_.begin());
    std::partial_sum(v_edge_offsets_.begin(), v_edge_offsets_.end(),
                     v_edge_offsets_.begin());
    u_incident_edges_.resize(unique_edges_.size());
    v_incident_edges_.resize(unique_edges_.size());
    std::vector<std::size_t> u_cursor = u_edge_offsets_;
    std::vector<std::size_t> v_cursor = v_edge_offsets_;
    for (std::size_t i = 0; i < unique_edges_.size(); ++i) {
        const UniqueEdgeID edge_id = static_cast<UniqueEdgeID>(i);
        const UniqueTemporalEdge& edge = unique_edges_[i];
        u_incident_edges_[u_cursor[edge.u]++] = edge_id;
        v_incident_edges_[v_cursor[edge.v]++] = edge_id;
    }
    u_attribute_incident_work_.assign(num_attributes(), 0);
    v_attribute_incident_work_.assign(num_attributes(), 0);
    for (VertexID u = 0; u < num_u_vertices_; ++u) {
        u_attribute_incident_work_[u_attributes_[u]] += static_u_degree(u);
    }
    for (VertexID v = 0; v < num_v_vertices_; ++v) {
        v_attribute_incident_work_[v_attributes_[v]] += static_v_degree(v);
    }
    return true;
}

std::size_t TemporalBipartiteGraph::attribute_incident_work(
    bool side_is_u, AttributeMask attributes) const {
    const auto& totals = side_is_u
        ? u_attribute_incident_work_ : v_attribute_incident_work_;
    std::size_t work = 0;
    for (std::size_t a = 0; a < totals.size() && a < 64; ++a) {
        if ((attributes & (AttributeMask{1} << a)) != 0) work += totals[a];
    }
    return work;
}

bool TemporalBipartiteGraph::active_in_window(
    const UniqueTemporalEdge& edge, TimeWindow window) const {
    const auto begin = unique_edge_timestamps_.begin() +
        static_cast<std::ptrdiff_t>(edge.timestamp_offset);
    const auto end = begin +
        static_cast<std::ptrdiff_t>(edge.timestamp_count);
    const auto found = std::lower_bound(begin, end, window.start);
    return found != end && *found <= window.end;
}

WindowAdjacencyStats TemporalBipartiteGraph::populate_window_adjacency(
    BipartiteSubgraph& subgraph, TimeWindow window,
    const std::vector<VertexID>& active_u,
    const std::vector<VertexID>& active_v) const {
    WindowAdjacencyStats stats;
    if (window.start > window.end ||
        subgraph.u_selected.size() != num_u_vertices_ ||
        subgraph.v_selected.size() != num_v_vertices_) {
        return stats;
    }

    std::size_t u_work = 0;
    for (VertexID u : active_u) {
        if (u < num_u_vertices_ && subgraph.u_selected[u]) {
            u_work += u_edge_offsets_[u + 1] - u_edge_offsets_[u];
        }
    }
    std::size_t v_work = 0;
    for (VertexID v : active_v) {
        if (v < num_v_vertices_ && subgraph.v_selected[v]) {
            v_work += v_edge_offsets_[v + 1] - v_edge_offsets_[v];
        }
    }
    stats.incident_work_u = u_work;
    stats.incident_work_v = v_work;

    auto add_edge = [&](UniqueEdgeID edge_id) {
        ++stats.unique_edges_examined;
        const UniqueTemporalEdge& edge = unique_edges_[edge_id];
        if (!subgraph.u_selected[edge.u] ||
            !subgraph.v_selected[edge.v] ||
            !active_in_window(edge, window)) {
            return;
        }
        ++stats.window_edges;
        subgraph.u_adj[edge.u].push_back(edge.v);
        subgraph.v_adj[edge.v].push_back(edge.u);
    };

    if (u_work <= v_work) {
        stats.traversed_from_u = true;
        for (VertexID u : active_u) {
            if (u >= num_u_vertices_ || !subgraph.u_selected[u]) continue;
            for (std::size_t i = u_edge_offsets_[u];
                 i < u_edge_offsets_[u + 1]; ++i) {
                add_edge(u_incident_edges_[i]);
            }
        }
    } else {
        stats.traversed_from_u = false;
        for (VertexID v : active_v) {
            if (v >= num_v_vertices_ || !subgraph.v_selected[v]) continue;
            for (std::size_t i = v_edge_offsets_[v];
                 i < v_edge_offsets_[v + 1]; ++i) {
                add_edge(v_incident_edges_[i]);
            }
        }
    }
    return stats;
}

WindowAdjacencyStats TemporalBipartiteGraph::populate_seed_window_adjacency(
    BipartiteSubgraph& subgraph, TimeWindow window,
    const std::vector<VertexID>& seeds, bool seed_is_u,
    AttributeMask other_attributes, std::vector<VertexID>& other_vertices,
    bool adaptive) const {
    WindowAdjacencyStats stats;
    stats.traversed_from_u = seed_is_u;
    other_vertices.clear();
    if (window.start > window.end ||
        subgraph.u_selected.size() != num_u_vertices_ ||
        subgraph.v_selected.size() != num_v_vertices_ ||
        subgraph.u_adj.size() != num_u_vertices_ ||
        subgraph.v_adj.size() != num_v_vertices_) return stats;

    if (adaptive) {
        stats.adaptive_seed_recovery = true;
        stats.opposite_scan_work = attribute_incident_work(
            !seed_is_u, other_attributes);
        for (VertexID seed : seeds) {
            if (seed_is_u) {
                if (seed < num_u_vertices_ && subgraph.u_selected[seed])
                    stats.seed_scan_work += static_u_degree(seed);
            } else {
                if (seed < num_v_vertices_ && subgraph.v_selected[seed])
                    stats.seed_scan_work += static_v_degree(seed);
            }
        }
        if (stats.opposite_scan_work < stats.seed_scan_work) {
            stats.reversed_seed_scan = true;
            stats.traversed_from_u = !seed_is_u;
            const auto add_from_opposite = [&](UniqueEdgeID edge_id) {
                ++stats.unique_edges_examined;
                const auto& edge = unique_edges_[edge_id];
                if (!(seed_is_u ? subgraph.u_selected[edge.u]
                               : subgraph.v_selected[edge.v]) ||
                    !active_in_window(edge, window)) return;
                auto& mark = seed_is_u ? subgraph.v_selected[edge.v]
                                       : subgraph.u_selected[edge.u];
                if (!mark) {
                    mark = 1;
                    other_vertices.push_back(seed_is_u ? edge.v : edge.u);
                }
                ++stats.window_edges;
                subgraph.u_adj[edge.u].push_back(edge.v);
                subgraph.v_adj[edge.v].push_back(edge.u);
            };
            const auto& attributes = seed_is_u ? v_attributes_ : u_attributes_;
            const auto& offsets = seed_is_u ? v_edge_offsets_ : u_edge_offsets_;
            const auto& incident = seed_is_u ? v_incident_edges_ : u_incident_edges_;
            for (VertexID vertex = 0; vertex < attributes.size(); ++vertex) {
                const AttributeID attr = attributes[vertex];
                if (attr >= 64 ||
                    (other_attributes & (AttributeMask{1} << attr)) == 0) continue;
                for (std::size_t p = offsets[vertex]; p < offsets[vertex + 1]; ++p)
                    add_from_opposite(incident[p]);
            }
            if (seed_is_u) {
                stats.incident_work_u = stats.seed_scan_work;
                for (VertexID v : other_vertices) stats.incident_work_v += static_v_degree(v);
            } else {
                stats.incident_work_v = stats.seed_scan_work;
                for (VertexID u : other_vertices) stats.incident_work_u += static_u_degree(u);
            }
            std::sort(other_vertices.begin(), other_vertices.end());
            return stats;
        }
    }

    auto add_edge = [&](UniqueEdgeID edge_id) {
        ++stats.unique_edges_examined;
        const auto& edge = unique_edges_[edge_id];
        const AttributeID attr = seed_is_u ? v_attributes_[edge.v] : u_attributes_[edge.u];
        if (attr >= 64 || (other_attributes & (AttributeMask{1} << attr)) == 0 ||
            !active_in_window(edge, window)) return;
        auto& mark = seed_is_u ? subgraph.v_selected[edge.v] : subgraph.u_selected[edge.u];
        if (!mark) {
            mark = 1;
            other_vertices.push_back(seed_is_u ? edge.v : edge.u);
        }
        ++stats.window_edges;
        subgraph.u_adj[edge.u].push_back(edge.v);
        subgraph.v_adj[edge.v].push_back(edge.u);
    };
    if (seed_is_u) {
        for (VertexID u : seeds) {
            if (u >= num_u_vertices_ || !subgraph.u_selected[u]) continue;
            stats.incident_work_u += u_edge_offsets_[u+1] - u_edge_offsets_[u];
            for (std::size_t p = u_edge_offsets_[u]; p < u_edge_offsets_[u+1]; ++p)
                add_edge(u_incident_edges_[p]);
        }
        for (VertexID v : other_vertices)
            stats.incident_work_v += v_edge_offsets_[v+1] - v_edge_offsets_[v];
    } else {
        for (VertexID v : seeds) {
            if (v >= num_v_vertices_ || !subgraph.v_selected[v]) continue;
            stats.incident_work_v += v_edge_offsets_[v+1] - v_edge_offsets_[v];
            for (std::size_t p = v_edge_offsets_[v]; p < v_edge_offsets_[v+1]; ++p)
                add_edge(v_incident_edges_[p]);
        }
        for (VertexID u : other_vertices)
            stats.incident_work_u += u_edge_offsets_[u+1] - u_edge_offsets_[u];
    }
    std::sort(other_vertices.begin(), other_vertices.end());
    return stats;
}

std::size_t TemporalBipartiteGraph::num_u_vertices() const {
    return num_u_vertices_;
}

std::size_t TemporalBipartiteGraph::num_v_vertices() const {
    return num_v_vertices_;
}

std::size_t TemporalBipartiteGraph::num_edges() const {
    return edges_.size();
}

std::size_t TemporalBipartiteGraph::num_attributes() const {
    return attribute_names_.size();
}

Timestamp TemporalBipartiteGraph::min_timestamp() const {
    return min_timestamp_;
}

Timestamp TemporalBipartiteGraph::max_timestamp() const {
    return max_timestamp_;
}

AttributeID TemporalBipartiteGraph::u_attribute(VertexID u) const {
    return u_attributes_.at(u);
}

AttributeID TemporalBipartiteGraph::v_attribute(VertexID v) const {
    return v_attributes_.at(v);
}

std::optional<AttributeID> TemporalBipartiteGraph::attribute_id(
    const std::string& name) const {
    const auto found = attribute_name_to_id_.find(name);
    if (found == attribute_name_to_id_.end()) return std::nullopt;
    return found->second;
}

const std::string& TemporalBipartiteGraph::attribute_name(AttributeID id) const {
    return attribute_names_.at(id);
}

VertexID TemporalBipartiteGraph::original_u_id(VertexID u) const {
    return original_u_ids_.at(u);
}

VertexID TemporalBipartiteGraph::original_v_id(VertexID v) const {
    return original_v_ids_.at(v);
}

const std::vector<TemporalEdge>& TemporalBipartiteGraph::edges() const {
    return edges_;
}

const std::vector<Timestamp>& TemporalBipartiteGraph::time_points() const {
    return time_points_;
}

std::size_t TemporalBipartiteGraph::num_unique_edges() const {
    return unique_edges_.size();
}

std::size_t TemporalBipartiteGraph::num_unique_temporal_edges() const {
    return unique_edge_timestamps_.size();
}

VertexID TemporalBipartiteGraph::unique_edge_u(std::size_t edge_id) const {
    return unique_edges_[edge_id].u;
}

VertexID TemporalBipartiteGraph::unique_edge_v(std::size_t edge_id) const {
    return unique_edges_[edge_id].v;
}

std::size_t TemporalBipartiteGraph::unique_edge_timestamp_offset(
    std::size_t edge_id) const {
    return unique_edges_[edge_id].timestamp_offset;
}

std::uint32_t TemporalBipartiteGraph::unique_edge_timestamp_count(
    std::size_t edge_id) const {
    return unique_edges_[edge_id].timestamp_count;
}

Timestamp TemporalBipartiteGraph::unique_edge_timestamp(
    std::size_t edge_id, std::size_t offset) const {
    const UniqueTemporalEdge& edge = unique_edges_[edge_id];
    return unique_edge_timestamps_[edge.timestamp_offset + offset];
}

std::size_t TemporalBipartiteGraph::static_u_degree(VertexID u) const {
    return u_edge_offsets_[u + 1] - u_edge_offsets_[u];
}

std::size_t TemporalBipartiteGraph::static_v_degree(VertexID v) const {
    return v_edge_offsets_[v + 1] - v_edge_offsets_[v];
}

VertexID TemporalBipartiteGraph::static_u_neighbor(
    VertexID u, std::size_t offset) const {
    return unique_edges_[u_incident_edges_[u_edge_offsets_[u] + offset]].v;
}

VertexID TemporalBipartiteGraph::static_v_neighbor(
    VertexID v, std::size_t offset) const {
    return unique_edges_[v_incident_edges_[v_edge_offsets_[v] + offset]].u;
}

std::size_t TemporalBipartiteGraph::temporal_query_aux_bytes() const {
    return unique_edges_.size() * sizeof(UniqueTemporalEdge) +
        unique_edge_timestamps_.size() * sizeof(Timestamp) +
        (u_edge_offsets_.size() + v_edge_offsets_.size()) *
            sizeof(std::size_t) +
        (u_incident_edges_.size() + v_incident_edges_.size()) *
            sizeof(UniqueEdgeID) +
        (u_attribute_incident_work_.size() + v_attribute_incident_work_.size()) *
            sizeof(std::size_t) + grouped_v_time_aux_bytes();
}

std::size_t GroupedVTimeAdjacency::bytes() const {
    return (v_groups.capacity() + starts.capacity()) * sizeof(std::size_t) +
        times.capacity() * sizeof(Timestamp) + neighbors.capacity() * sizeof(VertexID);
}

std::pair<std::size_t, std::size_t> GroupedVTimeAdjacency::range(
    VertexID v, TimeWindow window) const {
    const auto begin = times.begin() + v_groups[v];
    const auto end = times.begin() + v_groups[v + 1];
    const auto first = std::lower_bound(begin, end, window.start);
    const auto last = std::upper_bound(first, end, window.end);
    return {starts[static_cast<std::size_t>(first - times.begin())],
            starts[static_cast<std::size_t>(last - times.begin())]};
}

std::size_t TemporalBipartiteGraph::grouped_v_time_aux_bytes() const {
    const auto cached = std::atomic_load(&grouped_v_time_);
    return cached ? cached->bytes() : 0;
}

const GroupedVTimeAdjacency& TemporalBipartiteGraph::grouped_v_time_adjacency() const {
    auto cached = std::atomic_load(&grouped_v_time_);
    if (cached) return *cached;
    const auto start = std::chrono::steady_clock::now();
    auto data = std::make_shared<GroupedVTimeAdjacency>();
    data->v_groups.assign(num_v_vertices_ + 1, 0);
    for (const auto& edge : edges_) ++data->v_groups[edge.v + 1];
    std::partial_sum(data->v_groups.begin(), data->v_groups.end(), data->v_groups.begin());
    data->neighbors.resize(edges_.size());
    {
        std::vector<std::size_t> cursor(data->v_groups.begin(), data->v_groups.end() - 1);
        data->scratch_peak_bytes = cursor.capacity() * sizeof(std::size_t);
        for (std::size_t i = 0; i < edges_.size(); ++i)
            data->neighbors[cursor[edges_[i].v]++] = static_cast<VertexID>(i);
    }
    std::size_t groups = 0;
    for (VertexID v = 0; v < num_v_vertices_; ++v) {
        const auto begin = data->v_groups[v], end = data->v_groups[v + 1];
        for (auto i = begin; i < end; ++i)
            groups += i == begin || edges_[data->neighbors[i]].timestamp !=
                edges_[data->neighbors[i - 1]].timestamp;
    }
    data->times.resize(groups);
    data->starts.resize(groups + 1);
    std::size_t group = 0;
    for (VertexID v = 0; v < num_v_vertices_; ++v) {
        const auto begin = data->v_groups[v], end = data->v_groups[v + 1];
        data->v_groups[v] = group;
        Timestamp previous{};
        for (auto i = begin; i < end; ++i) {
            const auto& edge = edges_[data->neighbors[i]];
            if (i == begin || edge.timestamp != previous) {
                data->times[group] = edge.timestamp;
                data->starts[group++] = i;
            }
            previous = edge.timestamp;
            data->neighbors[i] = edge.u;
        }
    }
    data->v_groups[num_v_vertices_] = group;
    data->starts[group] = data->neighbors.size();
    data->build_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    std::shared_ptr<const GroupedVTimeAdjacency> prepared = data;
    if (std::atomic_compare_exchange_strong(&grouped_v_time_, &cached, prepared))
        cached = std::move(prepared);
    return *cached;
}

}
