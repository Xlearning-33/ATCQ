#pragma once

#include "common/graph.h"
#include "index/mat_index.h"
#include "index/uat_index.h"
#include "index/vat_index.h"

#include <cstdint>
#include <string>

namespace tacc {

std::uint64_t graph_fingerprint(const TemporalBipartiteGraph& graph);

std::uint64_t formal_index_hash(const VATIndex& index);
std::uint64_t formal_index_hash(const UATIndex& index);
std::uint64_t formal_index_hash(const MATIndex& index);

bool read_index_scope(const std::string& path, IndexScope& scope,
                      std::string* error = nullptr);

bool save_index_binary(const std::string& path,
                       const TemporalBipartiteGraph& graph,
                       const VATIndex& index,
                       std::string* error = nullptr);
bool save_index_binary(const std::string& path,
                       const TemporalBipartiteGraph& graph,
                       const UATIndex& index,
                       std::string* error = nullptr);
bool save_index_binary(const std::string& path,
                       const TemporalBipartiteGraph& graph,
                       const MATIndex& index,
                       std::string* error = nullptr);

bool load_index_binary(const std::string& path,
                       const TemporalBipartiteGraph& graph,
                       VATIndex& index,
                       std::string* error = nullptr);
bool load_index_binary(const std::string& path,
                       const TemporalBipartiteGraph& graph,
                       UATIndex& index,
                       std::string* error = nullptr);
bool load_index_binary(const std::string& path,
                       const TemporalBipartiteGraph& graph,
                       MATIndex& index,
                       std::string* error = nullptr,
                       bool verify_content = true);

bool save_index_text(const std::string& path,
                     const TemporalBipartiteGraph& graph,
                     const VATIndex& index,
                     std::string* error = nullptr,
                     const std::string& data_directory = "",
                     std::uint64_t build_peak_rss_bytes = 0);
bool save_index_intervals_text(const std::string& path,
                              const TemporalBipartiteGraph& graph,
                              const VATIndex& index,
                              std::string* error = nullptr,
                              const std::string& data_directory = "",
                              std::uint64_t build_peak_rss_bytes = 0);
bool save_index_text(const std::string& path,
                     const TemporalBipartiteGraph& graph,
                     const UATIndex& index,
                     std::string* error = nullptr);
bool save_index_text(const std::string& path,
                     const TemporalBipartiteGraph& graph,
                     const MATIndex& index,
                     std::string* error = nullptr,
                     const std::string& data_directory = "");

}
