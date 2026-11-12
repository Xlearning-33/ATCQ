#pragma once

#include "common/graph.h"

#include <filesystem>
#include <functional>
#include <iosfwd>
#include <string>
#include <utility>

namespace tacc::joint {

struct AttributePair { AttributeMask u{}, v{}; };
struct Interval { Timestamp start{}, end{}; };

struct AlphaRow { Count beta_begin{}; std::vector<std::vector<Interval>> beta_lists; };
struct VATAttributes { AttributePair attributes; Count alpha_begin{}; std::vector<AlphaRow> alpha_rows; };
struct VATVertex { VertexID original_id{}; std::vector<VATAttributes> attributes; };
struct VATDirectoryRef {
    std::uint32_t vertex_side{};
    std::uint32_t attribute_slot{};
};
static_assert(sizeof(VATDirectoryRef) == 8, "VAT directory references must occupy 8 bytes");
inline constexpr std::uint32_t vat_upper_flag = std::uint32_t{1} << 31;
using VATDirectoryGrid = std::vector<std::vector<std::vector<VATDirectoryRef>>>;
struct VATDirectory {
    std::vector<AttributePair> pairs;
    std::vector<VATDirectoryGrid> groups;
};
struct VAT {
    std::vector<VATVertex> u_vertices, v_vertices;
    VATDirectory directory;
    bool directory_ready{};
};

void build_vat_directory(VAT&);
std::size_t vat_directory_bytes(const VAT&);
std::size_t vat_directory_references(const VAT&);

struct UATEnd { Timestamp end{}; std::vector<VertexID> u_vertices, v_vertices; };
struct UATStart { Timestamp start{}; std::vector<UATEnd> ends; };
struct UATAttributes { AttributePair attributes; std::vector<UATStart> starts; };
struct UATCore { Count alpha{}, beta{}; std::vector<UATAttributes> attributes; };
struct UAT { std::vector<UATCore> cores; std::vector<Timestamp> time_points; };

struct BuildStats {
    std::size_t online_calls{}, full_window_prechecks{}, global_valid_cores{}, joint_valid_cores{};
    std::size_t attribute_pairs{}, completed_attribute_pairs{};
    std::size_t vat_core_records{}, vat_intervals{}, uat_blocks{}, vertex_refs{};
    std::size_t decompositions{}, windows{}, dynamic_bucket_updates{}, deleted_unique_edges{};
    std::size_t frontier_changes{}, emitted_memberships{}, empty_attribute_pairs{};
};

using BuildProgress = std::function<void(std::size_t completed, std::size_t total, AttributePair current)>;

VAT build_vat_basic(const TemporalBipartiteGraph&, BuildStats&, const BuildProgress& = {});
UAT build_uat_basic(const TemporalBipartiteGraph&, BuildStats&, const BuildProgress& = {});
VAT build_vat_opt(const TemporalBipartiteGraph&, BuildStats&, const BuildProgress& = {});
UAT build_uat_opt(const TemporalBipartiteGraph&, BuildStats&, const BuildProgress& = {});

QueryResult query(const VAT&, AttributePair, Count alpha, Count beta, TimeWindow, bool sort_result = true);
QueryResult query(const UAT&, AttributePair, Count alpha, Count beta, TimeWindow, bool sort_result = true);

std::size_t structure_bytes(const VAT&);
std::size_t structure_bytes(const UAT&);

void save_binary(const std::filesystem::path&, const TemporalBipartiteGraph&, const VAT&);
void save_binary(const std::filesystem::path&, const TemporalBipartiteGraph&, const UAT&);
void load_binary(const std::filesystem::path&, const TemporalBipartiteGraph&, VAT&);
void load_binary(const std::filesystem::path&, const TemporalBipartiteGraph&, UAT&);
void save_text(const std::filesystem::path&, const TemporalBipartiteGraph&, const VAT&,
               const std::filesystem::path& data, const std::string& generated_at,
               const std::string& build_method = "basic");
void save_text(const std::filesystem::path&, const TemporalBipartiteGraph&, const UAT&,
               const std::filesystem::path& data, const std::string& generated_at,
               const std::string& build_method = "basic");

}
