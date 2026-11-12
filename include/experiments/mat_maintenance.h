#pragma once

#include "common/graph.h"
#include "index/mat_index.h"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <ostream>
#include <string>

namespace tacc {

struct MATMaintenanceOptions {
    std::size_t updates{10}, batch_size{1};
    std::uint64_t seed{20261031};
    bool profile{}, component_pruning{}, window_pruning{}, compact_graph{};
    std::filesystem::path data, output;
    std::string prefix;
};

MATMaintenanceOptions configure_mat_maintenance(std::size_t updates, std::uint64_t seed);

void run_mat_maintenance(const TemporalBipartiteGraph& graph, MATIndex index,
                         const MATMaintenanceOptions& options, std::ostream& log);

}
