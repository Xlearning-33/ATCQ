#pragma once

#include "common/graph.h"
#include "common/types.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace tacc {

struct EncodedAttributes;

struct TABQueryKey {
    AttributeMask u_attributes{};
    AttributeMask v_attributes{};
    Timestamp start{};
    Timestamp end{};
    Count alpha{};
    Count beta{};

    bool operator==(const TABQueryKey& other) const;
};

struct TABQueryKeyHash {
    std::size_t operator()(const TABQueryKey& key) const;
};

struct TABRecord {
    TABQueryKey key;
    QueryResult result;
};

struct TABBuildStats {
    double build_seconds{};
    std::size_t enumerated_queries{};
    std::size_t nonempty_records{};
    std::size_t vertex_references{};
    std::size_t attribute_pairs{};
    std::size_t windows{};
    std::size_t approx_bytes{};
};

struct TABIndex {
    std::vector<TABRecord> records;
    std::unordered_map<TABQueryKey, std::size_t, TABQueryKeyHash> lookup;
    std::vector<Timestamp> time_points;
    AttributeMask u_attribute_universe{};
    AttributeMask v_attribute_universe{};
    Count max_alpha{};
    Count max_beta{};
    std::optional<std::uint64_t> source_graph_fingerprint;
    TABBuildStats stats;
};

using TABBuildProgress = std::function<void(std::size_t completed, std::size_t total)>;
TABIndex build_tab_index_exhaustive(const TemporalBipartiteGraph& graph,
                                   const TABBuildProgress& progress = {});

std::size_t tab_index_structure_bytes(const TABIndex& index);

QueryResult query_tab_index(const TABIndex& index,
                            const TemporalBipartiteGraph& graph,
                            const Query& query,
                            const EncodedAttributes* attributes = nullptr);

std::uint64_t tab_index_hash(const TABIndex& index);

bool save_tab_index_binary(const std::string& path,
                           const TABIndex& index,
                           std::string* error = nullptr);

bool load_tab_index_binary(const std::string& path,
                           TABIndex& index,
                           std::string* error = nullptr);

bool save_tab_index_text(const std::string& path,
                         const TemporalBipartiteGraph& graph,
                         const TABIndex& index,
                         std::string* error = nullptr,
                         const std::string& generated_at = {},
                         const std::string& data_directory = {});

}
