#pragma once

#include "common/graph.h"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iosfwd>
#include <string>

namespace tacc {

struct DatasetStatistics {
    bool raw_available{};
    std::uint64_t raw_events{};
    std::size_t temporal_edges{};
    std::size_t static_edges{};
    std::size_t active_u_vertices{};
    std::size_t active_v_vertices{};
    std::size_t u_attributes{};
    std::size_t v_attributes{};
    std::size_t max_u_degree{};
    std::size_t max_v_degree{};
    std::size_t symmetric_d{};
    std::size_t time_buckets{};
    Timestamp first_time{};
    Timestamp last_time{};
    std::string time_granularity;
};

struct DatasetStatsReport {
    DatasetStatistics stats;
    double stats_seconds{};
};

std::size_t balanced_core_threshold(const TemporalBipartiteGraph& graph);

DatasetStatsReport collect_dataset_stats(
    const TemporalBipartiteGraph& graph,
    const std::filesystem::path& data_dir);

void write_dataset_stats_fields(std::ostream& output,
                                const DatasetStatsReport& report);

void write_dataset_stats(std::ostream& output,
                         const DatasetStatsReport& report);

}
