#pragma once

#include "common/graph.h"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <tuple>
#include <vector>

namespace tacc::experiments {

constexpr std::size_t queries_per_setting = 500;
constexpr std::uint64_t workload_seed = 20261023;

struct Parameters {
    Count alpha{}, beta{};
    std::size_t window_points{}, u_attributes{}, v_attributes{};
    auto key() const {
        return std::make_tuple(alpha, beta, window_points, u_attributes, v_attributes);
    }
    bool operator<(const Parameters& other) const { return key() < other.key(); }
    bool excluded() const { return alpha == 1 && beta == 1; }
};

struct ParameterSetting {
    std::string family;
    unsigned percent{};
    Parameters parameters;
};

Parameters default_parameters(const TemporalBipartiteGraph& graph, Count delta);
std::vector<Query> generate_queries(const TemporalBipartiteGraph& graph,
                                    const Parameters& parameters,
                                    std::uint64_t seed = workload_seed);
std::vector<Query> generate_parameter_queries(const TemporalBipartiteGraph& graph,
                                              const Parameters& parameters);
std::vector<ParameterSetting> parameter_settings(const TemporalBipartiteGraph& graph,
                                                Count delta);
std::uint64_t workload_hash(const std::vector<Query>& queries);
std::size_t distinct_queries(const std::vector<Query>& queries);
std::vector<Query> load_queries(const std::filesystem::path& path,
                                const TemporalBipartiteGraph& graph,
                                const Parameters& parameters);

}
