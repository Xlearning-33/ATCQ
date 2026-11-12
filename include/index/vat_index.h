#pragma once

#include "index/index_types.h"
#include <memory>
#include <vector>

namespace tacc {
struct FormalIndexIOAccess;
struct VATInterval { Timestamp start{}, end{}; };

struct VATCoreEntry {
    Count alpha{}, beta{};
    std::size_t interval_offset{}, interval_count{};
};
struct VATAttributeEntry {
    AttributeMask mask{};
    std::vector<VATCoreEntry> cores;
};

struct VATRecordRef { VertexID vertex{}; std::uint32_t attribute_slot{}, core_slot{}; };
struct VATCoreRange { Count alpha{}, beta{}; std::size_t offset{}, count{}; };
struct VATMaskRange { AttributeMask mask{}; std::size_t offset{}, count{}; };

class VATIndex {
public:
    using VertexTree = std::vector<std::vector<VATAttributeEntry>>;
    bool add(IndexMembership record);
    void finalize();
    void clear();
    void initialize_vertices(std::size_t num_v);
    bool ready() const { return finalized_; }
    const VertexTree& vertices() const { return uv_vertices_; }
    const std::vector<VATInterval>& intervals() const { return intervals_; }
    const std::vector<VATMaskRange>& query_masks() const { return query_masks_; }
    const std::vector<VATCoreRange>& query_cores() const { return query_cores_; }
    const std::vector<VATRecordRef>& query_refs() const { return query_refs_; }
    std::size_t tree_bytes() const;
    std::size_t query_directory_bytes() const;
    void validate() const;
    const IndexBuildStats& stats() const { return stats_; }
    IndexBuildStats& mutable_stats() { return stats_; }

private:
    friend struct FormalIndexIOAccess;
    struct BuildCore { Count alpha{}, beta{}; std::vector<VATInterval> intervals; };
    struct BuildAttribute { AttributeMask mask{}; std::vector<BuildCore> cores; };
    using BuildTree = std::vector<std::vector<BuildAttribute>>;
    void recount();
    void build_query_directory();

    VertexTree uv_vertices_;
    std::unique_ptr<BuildTree> build_tree_;
    std::vector<VATInterval> intervals_;
    std::vector<VATMaskRange> query_masks_;
    std::vector<VATCoreRange> query_cores_;
    std::vector<VATRecordRef> query_refs_;
    IndexBuildStats stats_;
    bool finalized_{};
};
}
