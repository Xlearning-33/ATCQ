#include "experiments/dataset_stats.h"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace tacc {
namespace {

std::filesystem::path raw_event_file(const std::filesystem::path& data_dir) {
    const auto normalized = data_dir.lexically_normal();
    if (normalized.parent_path().filename() != "real_span") {
        return {};
    }
    const auto path = normalized.parent_path().parent_path() /
        "real_raw" / normalized.filename();
    std::error_code error;
    const auto status = std::filesystem::status(path, error);
    if (error == std::errc::no_such_file_or_directory) {
        return {};
    }
    if (error) {
        throw std::runtime_error("cannot_stat_raw_file");
    }
    return std::filesystem::is_regular_file(status) ? path : std::filesystem::path{};
}

std::string time_granularity(const std::string& dataset) {
    if (dataset == "04_mu") return "month";
    if (dataset == "06_na" || dataset == "09_lb" ||
        dataset == "10_ls") return "halfyear";
    if (dataset == "01_br" || dataset == "02_sg" ||
        dataset == "03_sv" || dataset == "05_ta" ||
        dataset == "07_so" || dataset == "08_wu") return "quarter";
    return "unspecified";
}

std::uint64_t count_raw_events(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot_open_raw_file");
    }
    std::uint64_t count = 0;
    std::string line;
    while (std::getline(input, line)) {
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos ||
            line[first] == '%' || line[first] == '#') {
            continue;
        }
        ++count;
    }
    if (input.bad()) {
        throw std::runtime_error("cannot_read_raw_file");
    }
    return count;
}

std::size_t symmetric_degeneracy(
    const TemporalBipartiteGraph& graph,
    std::size_t max_degree) {
    const std::size_t num_u = graph.num_u_vertices();
    const std::size_t num_v = graph.num_v_vertices();
    const std::size_t total = num_u + num_v;
    if (total == 0) return 0;

    std::vector<std::size_t> degree(total);
    std::vector<std::size_t> bin(max_degree + 1, 0);
    for (std::size_t u = 0; u < num_u; ++u) {
        degree[u] = graph.static_u_degree(static_cast<VertexID>(u));
        ++bin[degree[u]];
    }
    for (std::size_t v = 0; v < num_v; ++v) {
        degree[num_u + v] =
            graph.static_v_degree(static_cast<VertexID>(v));
        ++bin[degree[num_u + v]];
    }

    std::size_t next = 0;
    for (std::size_t d = 0; d < bin.size(); ++d) {
        const auto count = bin[d];
        bin[d] = next;
        next += count;
    }
    std::vector<std::size_t> position(total);
    std::vector<std::size_t> order(total);
    for (std::size_t vertex = 0; vertex < total; ++vertex) {
        const auto slot = bin[degree[vertex]]++;
        order[slot] = vertex;
        position[vertex] = slot;
    }
    for (std::size_t d = max_degree; d > 0; --d) {
        bin[d] = bin[d - 1];
    }
    bin[0] = 0;

    std::size_t maximum_core = 0;
    for (std::size_t i = 0; i < total; ++i) {
        const std::size_t vertex = order[i];
        const std::size_t current_degree = degree[vertex];
        maximum_core = std::max(maximum_core, current_degree);
        auto decrement_higher_neighbor = [&](std::size_t neighbor) {
            if (degree[neighbor] <= current_degree) return;
            const auto old_degree = degree[neighbor];
            const auto neighbor_slot = position[neighbor];
            const auto first_slot = bin[old_degree];
            const auto first_vertex = order[first_slot];
            if (neighbor != first_vertex) {
                std::swap(order[neighbor_slot], order[first_slot]);
                position[neighbor] = first_slot;
                position[first_vertex] = neighbor_slot;
            }
            ++bin[old_degree];
            --degree[neighbor];
        };
        if (vertex < num_u) {
            const auto u = static_cast<VertexID>(vertex);
            for (std::size_t j = 0; j < graph.static_u_degree(u); ++j) {
                decrement_higher_neighbor(
                    num_u + graph.static_u_neighbor(u, j));
            }
        } else {
            const auto v = static_cast<VertexID>(vertex - num_u);
            for (std::size_t j = 0; j < graph.static_v_degree(v); ++j) {
                decrement_higher_neighbor(graph.static_v_neighbor(v, j));
            }
        }
    }
    return maximum_core;
}

DatasetStatistics calculate_statistics(
    const TemporalBipartiteGraph& graph,
    const std::filesystem::path& raw_file,
    const std::string& granularity) {
    DatasetStatistics stats;
    stats.raw_available = !raw_file.empty();
    if (stats.raw_available) {
        stats.raw_events = count_raw_events(raw_file);
    }
    stats.temporal_edges = graph.num_unique_temporal_edges();
    stats.static_edges = graph.num_unique_edges();
    stats.time_buckets = graph.time_points().size();
    stats.first_time = graph.min_timestamp();
    stats.last_time = graph.max_timestamp();
    stats.time_granularity = granularity;

    std::vector<std::uint8_t> u_seen(graph.num_attributes(), 0);
    std::vector<std::uint8_t> v_seen(graph.num_attributes(), 0);
    for (std::size_t i = 0; i < graph.num_u_vertices(); ++i) {
        const auto u = static_cast<VertexID>(i);
        const auto degree = graph.static_u_degree(u);
        if (degree == 0) continue;
        ++stats.active_u_vertices;
        stats.max_u_degree = std::max(stats.max_u_degree, degree);
        const auto attribute = graph.u_attribute(u);
        if (!u_seen[attribute]) {
            u_seen[attribute] = 1;
            ++stats.u_attributes;
        }
    }
    for (std::size_t i = 0; i < graph.num_v_vertices(); ++i) {
        const auto v = static_cast<VertexID>(i);
        const auto degree = graph.static_v_degree(v);
        if (degree == 0) continue;
        ++stats.active_v_vertices;
        stats.max_v_degree = std::max(stats.max_v_degree, degree);
        const auto attribute = graph.v_attribute(v);
        if (!v_seen[attribute]) {
            v_seen[attribute] = 1;
            ++stats.v_attributes;
        }
    }
    stats.symmetric_d = symmetric_degeneracy(
        graph, std::max(stats.max_u_degree, stats.max_v_degree));
    return stats;
}

}

std::size_t balanced_core_threshold(const TemporalBipartiteGraph& graph) {
    std::size_t maximum_degree = 0;
    for (VertexID u = 0; u < graph.num_u_vertices(); ++u)
        maximum_degree = std::max(maximum_degree, graph.static_u_degree(u));
    for (VertexID v = 0; v < graph.num_v_vertices(); ++v)
        maximum_degree = std::max(maximum_degree, graph.static_v_degree(v));
    return symmetric_degeneracy(graph, maximum_degree);
}

DatasetStatsReport collect_dataset_stats(
    const TemporalBipartiteGraph& graph,
    const std::filesystem::path& data_dir) {
    const auto start = std::chrono::steady_clock::now();
    const auto raw_file = raw_event_file(data_dir);
    const std::string granularity =
        time_granularity(data_dir.lexically_normal().filename().string());
    DatasetStatsReport report;
    report.stats = calculate_statistics(graph, raw_file, granularity);
    report.stats_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    return report;
}

void write_dataset_stats_fields(std::ostream& output,
                                const DatasetStatsReport& report) {
    const auto& stats = report.stats;
    output << "TAC_DATASET_STATS\nraw_file_available " << static_cast<int>(stats.raw_available)
           << "\nraw_events ";
    if (stats.raw_available) output << stats.raw_events;
    else output << "NA";
    output << '\n'
           << "temporal_edges " << stats.temporal_edges << '\n'
           << "static_edges " << stats.static_edges << '\n'
           << "active_u_vertices " << stats.active_u_vertices << '\n'
           << "active_v_vertices " << stats.active_v_vertices << '\n'
           << "u_attributes " << stats.u_attributes << '\n'
           << "v_attributes " << stats.v_attributes << '\n'
           << "max_u_degree " << stats.max_u_degree << '\n'
           << "max_v_degree " << stats.max_v_degree << '\n'
           << "symmetric_d " << stats.symmetric_d << '\n'
           << "time_granularity " << stats.time_granularity << '\n'
           << "time_buckets " << stats.time_buckets << '\n'
           << "first_time " << stats.first_time << '\n'
           << "last_time " << stats.last_time << '\n';
}

void write_dataset_stats(std::ostream& output,
                         const DatasetStatsReport& report) {
    const auto& stats = report.stats;
    output << "RESULT exp=dataset_stats status=ok"
           << " stats_s=" << report.stats_seconds
           << " raw_available=" << (stats.raw_available ? 1 : 0)
           << " raw_events=";
    if (stats.raw_available) output << stats.raw_events;
    else output << "NA";
    output << " temporal_edges=" << stats.temporal_edges
           << " static_edges=" << stats.static_edges
           << " active_u_vertices=" << stats.active_u_vertices
           << " active_v_vertices=" << stats.active_v_vertices
           << " u_attributes=" << stats.u_attributes
           << " v_attributes=" << stats.v_attributes
           << " max_u_degree=" << stats.max_u_degree
           << " max_v_degree=" << stats.max_v_degree
           << " symmetric_d=" << stats.symmetric_d
           << " time_granularity=" << stats.time_granularity
           << " time_buckets=" << stats.time_buckets
           << " first_time=" << stats.first_time
           << " last_time=" << stats.last_time
           << '\n';
}

}
