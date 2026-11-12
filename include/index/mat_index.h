#pragma once

#include "index/index_types.h"

#include <cstddef>
#include <cstdint>
#include <tuple>
#include <utility>
#include <vector>

namespace tacc {

struct FormalIndexIOAccess;

class MATIndex {
public:
    MATIndex();
    using QueryKey = std::tuple<AttributeMask, Count>;
    using QueryDirectory =
        std::vector<std::pair<QueryKey, IndexPostingRange>>;
    struct TimeGroup {
        Timestamp start{};
        Timestamp end{};
        IndexPostingRange postings;
    };
    using QueryTimeGroups = std::vector<TimeGroup>;

    bool add(IndexMembership record);
    void finalize();
    void clear();
    const std::vector<IndexBlock>& blocks() const { return blocks_; }
    const std::vector<VertexID>& vertex_pool() const { return vertex_pool_; }
    const QueryDirectory& query_directory() const { return query_directory_; }
    const QueryTimeGroups& query_time_groups() const {
        return query_time_groups_;
    }
    const IndexBuildStats& stats() const { return stats_; }
    IndexBuildStats& mutable_stats() { return stats_; }

private:
    friend struct FormalIndexIOAccess;

    struct BucketStats {
        Timestamp min_start{};
        Timestamp max_start{};
        Timestamp min_end{};
        Timestamp max_end{};
        Count min_variable{};
        Count max_variable{};
        bool valid{};
    };
    struct MaskBucket {
        struct StartGroup {
            Timestamp start{};
            std::uint32_t begin{};
            std::uint32_t end{};
        };
        AttributeMask attributes{};
        std::vector<IndexMembership> memberships;
        BucketStats stats;
        bool starts_sorted{true};
        bool staircase_sorted{true};
        std::vector<StartGroup> start_groups;
        bool start_groups_valid{};
    };
    struct CoreBucket {
        Count fixed{};
        std::vector<MaskBucket> masks;
    };
    struct VertexBucket {
        VertexID vertex{};
        std::vector<std::uint32_t> core_positions;
        std::vector<CoreBucket> active_cores;
    };
    struct BuildBuckets {
        std::vector<std::uint32_t> vertex_positions;
        std::vector<VertexBucket> active_vertices;
    };
    template <bool Profile>
    bool add_impl(IndexMembership record);
    static void update_bucket_stats(BucketStats& stats,
                                    const IndexMembership& record);
    static void recompute_bucket_stats(MaskBucket& bucket);
    static void update_bucket_order(MaskBucket& bucket,
                                    const IndexMembership& record);
    void clear_build_buckets();
    void rebuild_query_directory();
    std::size_t query_directory_bytes() const;

    BuildBuckets build_buckets_;
    std::vector<IndexBlock> blocks_;
    std::vector<VertexID> vertex_pool_;
    QueryDirectory query_directory_;
    QueryTimeGroups query_time_groups_;
    IndexBuildStats stats_;
    bool streaming_profile_enabled_{};
    bool ordered_streaming_enabled_{true};
};

}
