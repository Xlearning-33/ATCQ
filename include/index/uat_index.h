#pragma once

#include "index/index_types.h"
#include "index/pool_interner.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

namespace tacc {

struct FormalIndexIOAccess;

struct UATEndEntry {
    Timestamp end{};
    std::size_t vertex_offset{}, vertex_count{};
};
struct UATStartEntry {
    Timestamp start{};
    std::vector<UATEndEntry> ends;
};
struct UATAttributeEntry {
    AttributeMask mask{};
    std::vector<UATStartEntry> starts;
};
struct UATCoreEntry {
    Count alpha{}, beta{};
    std::vector<UATAttributeEntry> attributes;
};

class UATIndex {
public:
    bool add(IndexMembership record);
    void finalize();
    void clear();
    const std::vector<UATCoreEntry>& cores() const { return cores_; }
    const std::vector<VertexID>& vertex_pool() const { return vertex_pool_; }
    bool ready() const { return finalized_; }
    void validate() const;
    std::size_t tree_bytes() const;
    template<class Visit> void for_each_block(Visit visit) const {
        for (const auto& c : cores_) for (const auto& a : c.attributes)
            for (const auto& s : a.starts) for (const auto& e : s.ends)
                visit(IndexBlock{a.mask,s.start,e.end,c.alpha,c.beta,
                                 e.vertex_offset,e.vertex_count});
    }
    const IndexBuildStats& stats() const { return stats_; }
    IndexBuildStats& mutable_stats() { return stats_; }

private:
    friend struct FormalIndexIOAccess;

    using BuildAttributeID = std::uint32_t;

    struct BuildMembership {
        BuildAttributeID attribute_id{};
        Timestamp start{};
        Timestamp end{};
    };
    struct CoreBucket {
        Count alpha{};
        Count beta{};
        VertexID vertex{};
        std::vector<BuildMembership> memberships;
    };
    struct VertexBucket {
        VertexID vertex{};
        std::vector<CoreBucket> cores;
    };
    struct BuildBuckets {
        std::vector<std::uint32_t> vertex_positions;
        std::vector<VertexBucket> active_vertices;
    };
    BuildAttributeID build_attribute_id(AttributeMask attributes);
    void clear_build_buckets();
    void finalize_buckets();
    void recount();

    BuildBuckets build_buckets_;
    std::vector<AttributeMask> build_attribute_masks_;
    AttributeMask last_build_attribute_mask_{};
    BuildAttributeID last_build_attribute_id_{};
    bool has_last_build_attribute_{};
    bool finalized_{};
    std::size_t logical_references_{};
    std::optional<PoolInterner> build_interner_;
    std::vector<UATCoreEntry> cores_;
    std::vector<VertexID> vertex_pool_;
    IndexBuildStats stats_;
};

}
