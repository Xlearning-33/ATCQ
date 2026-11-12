#include "index/mat_index.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

namespace tacc {
namespace {
using Clock = std::chrono::steady_clock;

template <bool Enabled>
struct SampleTimer {
    double* seconds{};
    Clock::time_point start;
    SampleTimer(double& accumulator, bool sample) {
        if constexpr (Enabled) {
            if (sample) { seconds = &accumulator; start = Clock::now(); }
        }
    }
    void stop() {
        if constexpr (Enabled) {
            if (seconds != nullptr) {
                *seconds += std::chrono::duration<double>(Clock::now() - start).count();
                seconds = nullptr;
            }
        }
    }
    ~SampleTimer() { stop(); }
};
bool covers(const IndexMembership& lhs, const IndexMembership& rhs) {
    if (!attribute_subset(lhs.attributes, rhs.attributes) ||
        lhs.start < rhs.start || lhs.end > rhs.end) {
        return false;
    }
    return lhs.beta == rhs.beta && lhs.alpha >= rhs.alpha;
}
using BlockKey = std::tuple<AttributeMask, Timestamp,
                            Timestamp, Count, Count>;

BlockKey block_key(const IndexMembership& record) {
    return {record.attributes,
            record.start, record.end, record.alpha, record.beta};
}

Count fixed_core(const IndexBlock& block) {
    return block.beta;
}

Count skyline_core(const IndexBlock& block) {
    return block.alpha;
}

Count fixed_core(const IndexMembership& record) {
    return record.beta;
}

Count skyline_core(const IndexMembership& record) {
    return record.alpha;
}

auto query_sort_key(const IndexBlock& block) {
    return std::make_tuple(block.attributes, fixed_core(block),
                           block.start, block.end, skyline_core(block));
}

std::uint64_t vertex_set_hash(const std::vector<VertexID>& vertices) {
    std::uint64_t hash = 1469598103934665603ULL;
    hash ^= static_cast<std::uint64_t>(vertices.size()) +
        0x9e3779b97f4a7c15ULL + (hash << 6U) + (hash >> 2U);
    for (VertexID vertex : vertices) {
        hash ^= static_cast<std::uint64_t>(vertex) +
            0x9e3779b97f4a7c15ULL + (hash << 6U) + (hash >> 2U);
    }
    return hash;
}

bool same_vertex_set(const std::vector<VertexID>& pool,
                     IndexPostingRange range,
                     const std::vector<VertexID>& vertices) {
    return range.count == vertices.size() &&
        range.offset <= pool.size() &&
        range.count <= pool.size() - range.offset &&
        std::equal(vertices.begin(), vertices.end(),
                   pool.begin() + static_cast<std::ptrdiff_t>(range.offset));
}
}

MATIndex::MATIndex() {
    const char* value = std::getenv("TAC_PROFILE_MAT_STREAMING");
    streaming_profile_enabled_ = value != nullptr && value[0] == '1' && value[1] == '\0';
    const char* ordered = std::getenv("TAC_MAT_ORDERED_STREAMING");
    ordered_streaming_enabled_ = !(ordered != nullptr && ordered[0] == '0' && ordered[1] == '\0');
    stats_.mat_ordered_streaming_enabled = ordered_streaming_enabled_;
    stats_.mat_streaming_mode_known = true;
}

void MATIndex::update_bucket_stats(
    BucketStats& stats, const IndexMembership& record) {
    const Count variable = skyline_core(record);
    if (!stats.valid) {
        stats.min_start = stats.max_start = record.start;
        stats.min_end = stats.max_end = record.end;
        stats.min_variable = stats.max_variable = variable;
        stats.valid = true;
        return;
    }
    stats.min_start = std::min(stats.min_start, record.start);
    stats.max_start = std::max(stats.max_start, record.start);
    stats.min_end = std::min(stats.min_end, record.end);
    stats.max_end = std::max(stats.max_end, record.end);
    stats.min_variable = std::min(stats.min_variable, variable);
    stats.max_variable = std::max(stats.max_variable, variable);
}

void MATIndex::recompute_bucket_stats(MaskBucket& bucket) {
    bucket.stats = {};
    for (const IndexMembership& membership : bucket.memberships) {
        update_bucket_stats(bucket.stats, membership);
    }
}

void MATIndex::update_bucket_order(MaskBucket& bucket,
                                  const IndexMembership& record) {
    bucket.start_groups_valid = false;
    if (bucket.memberships.empty()) {
        bucket.starts_sorted = bucket.staircase_sorted = true;
        return;
    }
    const auto& last = bucket.memberships.back();
    if (last.start > record.start) {
        bucket.starts_sorted = bucket.staircase_sorted = false;
    } else if (last.start == record.start &&
               (last.end < record.end || skyline_core(last) < skyline_core(record))) {
        bucket.staircase_sorted = false;
    }
}

void MATIndex::clear_build_buckets() {
    auto& buckets = build_buckets_;

        std::vector<std::uint32_t>().swap(buckets.vertex_positions);
        std::vector<VertexBucket>().swap(buckets.active_vertices);

}

bool MATIndex::add(IndexMembership record) {
    if (streaming_profile_enabled_) return add_impl<true>(record);
    return add_impl<false>(record);
}

template <bool Profile>
bool MATIndex::add_impl(IndexMembership record) {
    ++stats_.generated_records;
    auto& diagnostics = stats_.mat_streaming_profile;
    const bool sample = Profile &&
        ((stats_.generated_records * std::uint64_t{11400714819323198485ULL}) >> 54U) == 0;
    if constexpr (Profile) {
        diagnostics.enabled = true;
        if (sample) ++diagnostics.sampled_calls;
    }
    SampleTimer<Profile> total_timer(diagnostics.sampled_add_seconds, sample);
    SampleTimer<Profile> locator_timer(diagnostics.sampled_locator_seconds, sample);
    BuildBuckets& buckets =
        build_buckets_;
    constexpr std::uint32_t kMissing =
        std::numeric_limits<std::uint32_t>::max();
    if (record.vertex >= buckets.vertex_positions.size()) {
        buckets.vertex_positions.resize(
            static_cast<std::size_t>(record.vertex) + 1, kMissing);
    }
    std::uint32_t& position = buckets.vertex_positions[record.vertex];
    if (position == kMissing) {
        if (buckets.active_vertices.size() >= kMissing) {
            throw std::overflow_error("MAT active vertex locator exceeds uint32_t");
        }
        position = static_cast<std::uint32_t>(
            buckets.active_vertices.size());
        buckets.active_vertices.push_back({record.vertex, {}, {}});
    }
    VertexBucket& vertex = buckets.active_vertices[position];
    const Count fixed = fixed_core(record);
    if (fixed >= vertex.core_positions.size()) {
        vertex.core_positions.resize(
            static_cast<std::size_t>(fixed) + 1, kMissing);
    }
    std::uint32_t& core_position = vertex.core_positions[fixed];
    if (core_position == kMissing) {
        if (vertex.active_cores.size() >= kMissing) {
            throw std::overflow_error("MAT active core locator exceeds uint32_t");
        }
        core_position = static_cast<std::uint32_t>(
            vertex.active_cores.size());
        vertex.active_cores.push_back({fixed, {}});
    }
    CoreBucket& core = vertex.active_cores[core_position];
    const Count variable = skyline_core(record);
    locator_timer.stop();

    SampleTimer<Profile> reject_timer(diagnostics.sampled_reject_seconds, sample);
    for (MaskBucket& bucket : core.masks) {
        if constexpr (Profile) ++diagnostics.reject_buckets;
        if (!attribute_subset(bucket.attributes, record.attributes)) continue;
        const BucketStats& bounds = bucket.stats;
        if (!bounds.valid || bounds.max_start < record.start ||
            bounds.min_end > record.end ||
            bounds.max_variable < variable) {
            continue;
        }
        if constexpr (Profile) ++diagnostics.reject_bounds_passed;
        const auto covered_by = [&](const IndexMembership& kept) {
            if constexpr (Profile) ++diagnostics.reject_comparisons;
            if (!covers(kept, record)) return false;
            if (kept.alpha != record.alpha || kept.beta != record.beta) {
                ++stats_.cross_core_removed;
            } else {
                ++stats_.attribute_time_removed;
            }
            if constexpr (Profile) ++diagnostics.rejected_records;
            return true;
        };
        const auto& memberships = bucket.memberships;
        if (ordered_streaming_enabled_ && bucket.staircase_sorted && !memberships.empty() &&
            record.start == memberships.back().start &&
            record.end < memberships.back().end) {
            continue;
        }
        if (ordered_streaming_enabled_ && bucket.staircase_sorted &&
            bucket.attributes != record.attributes &&
            memberships.size() <= std::numeric_limits<std::uint32_t>::max()) {
            if constexpr (Profile) ++diagnostics.ordered_reject_buckets;
            if (!bucket.start_groups_valid) {
                const auto old_capacity = bucket.start_groups.capacity();
                bucket.start_groups.clear();
                for (std::size_t begin = 0; begin < memberships.size();) {
                    std::size_t end = begin + 1;
                    while (end < memberships.size() &&
                           memberships[end].start == memberships[begin].start) ++end;
                    bucket.start_groups.push_back({memberships[begin].start,
                        static_cast<std::uint32_t>(begin), static_cast<std::uint32_t>(end)});
                    begin = end;
                }
                bucket.start_groups_valid = true;
                if constexpr (Profile) {
                    diagnostics.start_group_build_records += memberships.size();
                    diagnostics.start_group_capacity_bytes +=
                        (bucket.start_groups.capacity() - old_capacity) *
                        sizeof(MaskBucket::StartGroup);
                }
            }
            auto group = std::lower_bound(bucket.start_groups.begin(),
                bucket.start_groups.end(), record.start,
                [](const MaskBucket::StartGroup& group, Timestamp start) {
                    return group.start < start;
                });
            for (; group != bucket.start_groups.end(); ++group) {
                if constexpr (Profile) ++diagnostics.reject_time_groups;
                if (memberships[group->end - 1].end > record.end ||
                    skyline_core(memberships[group->begin]) < variable) continue;
                const auto last = memberships.begin() + group->end;
                const auto begin = memberships.begin() + group->begin;
                const auto candidate = begin->end <= record.end ? begin :
                    std::lower_bound(begin, last, record.end,
                    [](const IndexMembership& kept, Timestamp end) {
                        return kept.end > end;
                    });
                if (candidate != last && covered_by(*candidate)) return false;
            }
            continue;
        }
        auto first = memberships.begin();
        if (ordered_streaming_enabled_ && bucket.starts_sorted) {
            first = std::lower_bound(first, memberships.end(), record.start,
                [](const IndexMembership& kept, Timestamp start) {
                    return kept.start < start;
                });
        }
        if (ordered_streaming_enabled_ && bucket.staircase_sorted) {
            if constexpr (Profile) ++diagnostics.ordered_reject_buckets;
            while (first != memberships.end()) {
                if constexpr (Profile) ++diagnostics.reject_time_groups;
                const auto last = std::upper_bound(first, memberships.end(), first->start,
                    [](Timestamp start, const IndexMembership& kept) {
                        return start < kept.start;
                    });
                if ((last - 1)->end > record.end || skyline_core(*first) < variable) {
                    first = last;
                    continue;
                }
                const auto candidate = std::lower_bound(first, last, record.end,
                    [](const IndexMembership& kept, Timestamp end) {
                        return kept.end > end;
                    });
                if (candidate != last && covered_by(*candidate)) return false;
                first = last;
            }
        } else {
            for (; first != memberships.end(); ++first)
                if (covered_by(*first)) return false;
        }
    }
    reject_timer.stop();

    SampleTimer<Profile> remove_timer(diagnostics.sampled_remove_seconds, sample);
    for (MaskBucket& bucket : core.masks) {
        if constexpr (Profile) ++diagnostics.remove_buckets;
        if (!attribute_subset(record.attributes, bucket.attributes)) continue;
        const BucketStats& bounds = bucket.stats;
        if (!bounds.valid || bounds.min_start > record.start ||
            bounds.max_end < record.end ||
            bounds.min_variable > variable) {
            continue;
        }
        if constexpr (Profile) ++diagnostics.remove_bounds_passed;
        auto& memberships = bucket.memberships;
        const auto remove_record = [&](const IndexMembership& kept) {
            if constexpr (Profile) ++diagnostics.remove_comparisons;
            if (!covers(record, kept)) return false;
            if (record.alpha != kept.alpha || record.beta != kept.beta) {
                ++stats_.cross_core_removed;
            } else {
                ++stats_.attribute_time_removed;
            }
            if constexpr (Profile) ++diagnostics.erased_records;
            return true;
        };
        const auto new_end = std::remove_if(
            memberships.begin(), memberships.end(), remove_record);
        if (new_end != memberships.end()) {
            bucket.start_groups_valid = false;
            memberships.erase(new_end, memberships.end());
            if constexpr (Profile) diagnostics.refreshed_records += memberships.size();
            recompute_bucket_stats(bucket);
        }
    }
    remove_timer.stop();

    SampleTimer<Profile> timer(diagnostics.sampled_append_seconds, sample);
    MaskBucket* target = nullptr;
    for (MaskBucket& bucket : core.masks) {
        if (bucket.attributes == record.attributes) {
            target = &bucket;
            break;
        }
    }
    if (target == nullptr) {
        core.masks.emplace_back();
        core.masks.back().attributes = record.attributes;
        target = &core.masks.back();
    }
    if (ordered_streaming_enabled_) update_bucket_order(*target, record);
    target->memberships.push_back(std::move(record));
    update_bucket_stats(target->stats, target->memberships.back());
    if constexpr (Profile) diagnostics.max_bucket_records = std::max(
        diagnostics.max_bucket_records, target->memberships.size());
    return true;
}

void MATIndex::finalize() {
    std::size_t memberships = 0;
    { const auto& buckets = build_buckets_;
        for (const VertexBucket& vertex : buckets.active_vertices) {
            for (const CoreBucket& core : vertex.active_cores) {
                for (const MaskBucket& bucket : core.masks) {
                    memberships += bucket.memberships.size();
                }
            }
        }
    }

    std::vector<const IndexMembership*> flat;
    flat.reserve(memberships);
    { const auto& buckets = build_buckets_;
        for (const VertexBucket& vertex : buckets.active_vertices) {
            for (const CoreBucket& core : vertex.active_cores) {
                for (const MaskBucket& bucket : core.masks) {
                    for (const IndexMembership& membership :
                         bucket.memberships) {
                        flat.push_back(&membership);
                    }
                }
            }
        }
    }
    std::sort(flat.begin(), flat.end(),
              [](const IndexMembership* lhs,
                 const IndexMembership* rhs) {
                  return std::tie(lhs->attributes, lhs->start, lhs->end,
                                  lhs->alpha, lhs->beta, lhs->vertex) <
                         std::tie(rhs->attributes, rhs->start, rhs->end,
                                  rhs->alpha, rhs->beta, rhs->vertex);
              });

    std::size_t block_count = 0;
    for (std::size_t i = 0; i < flat.size();) {
        ++block_count;
        const BlockKey key = block_key(*flat[i]);
        do {
            ++i;
        } while (i < flat.size() && block_key(*flat[i]) == key);
    }

    blocks_.clear();
    vertex_pool_.clear();
    blocks_.reserve(block_count);
    vertex_pool_.reserve(flat.size());
    constexpr std::size_t kNoRange =
        std::numeric_limits<std::size_t>::max();
    std::unordered_map<std::uint64_t, std::size_t> range_heads;
    range_heads.reserve(block_count);
    std::vector<IndexPostingRange> interned_ranges;
    std::vector<std::size_t> next_same_hash;
    interned_ranges.reserve(block_count);
    next_same_hash.reserve(block_count);
    std::vector<VertexID> group_vertices;
    std::size_t logical_vertex_references = 0;
    for (std::size_t begin = 0; begin < flat.size();) {
        const BlockKey key = block_key(*flat[begin]);
        std::size_t group_end = begin + 1;
        while (group_end < flat.size() &&
               block_key(*flat[group_end]) == key) {
            ++group_end;
        }
        const auto [attrs, start, end, alpha, beta] = key;
        group_vertices.clear();
        group_vertices.reserve(group_end - begin);
        VertexID previous = 0;
        bool has_previous = false;
        for (std::size_t i = begin; i < group_end; ++i) {
            if (has_previous && flat[i]->vertex == previous) continue;
            previous = flat[i]->vertex;
            has_previous = true;
            group_vertices.push_back(previous);
        }
        logical_vertex_references += group_vertices.size();
        const std::uint64_t hash = vertex_set_hash(group_vertices);
        IndexPostingRange canonical{};
        bool reused = false;
        auto head = range_heads.find(hash);
        for (std::size_t candidate = head == range_heads.end()
                 ? kNoRange : head->second;
             candidate != kNoRange;
             candidate = next_same_hash[candidate]) {
            if (!same_vertex_set(vertex_pool_, interned_ranges[candidate],
                                 group_vertices)) {
                continue;
            }
            canonical = interned_ranges[candidate];
            reused = true;
            break;
        }
        if (!reused) {
            canonical = {vertex_pool_.size(), group_vertices.size()};
            vertex_pool_.insert(vertex_pool_.end(), group_vertices.begin(),
                                group_vertices.end());
            const std::size_t position = interned_ranges.size();
            interned_ranges.push_back(canonical);
            next_same_hash.push_back(
                head == range_heads.end() ? kNoRange : head->second);
            if (head == range_heads.end()) {
                range_heads.emplace(hash, position);
            } else {
                head->second = position;
            }
        }
        blocks_.push_back({attrs, start, end,
                           alpha, beta, canonical.offset, canonical.count});
        begin = group_end;
    }
    clear_build_buckets();
    std::vector<const IndexMembership*>().swap(flat);
    rebuild_query_directory();
    stats_.kept_records = logical_vertex_references;
    stats_.grouped_blocks = blocks_.size();
    stats_.vertex_references = vertex_pool_.size();
    stats_.approx_bytes = blocks_.size() * sizeof(IndexBlock) +
        vertex_pool_.size() * sizeof(VertexID) +
        query_directory_bytes();
}

void MATIndex::rebuild_query_directory() {
    query_directory_.clear();
    query_time_groups_.clear();
    std::sort(blocks_.begin(), blocks_.end(),
              [](const IndexBlock& left, const IndexBlock& right) {
                  return query_sort_key(left) < query_sort_key(right);
              });
    for (std::size_t key_begin = 0; key_begin < blocks_.size();) {
        const IndexBlock& first = blocks_[key_begin];
        const QueryKey key{first.attributes, fixed_core(first)};
        std::size_t key_end = key_begin;
        while (key_end < blocks_.size()) {
            const IndexBlock& next = blocks_[key_end];
            if (QueryKey{next.attributes, fixed_core(next)} != key) {
                break;
            }
            ++key_end;
        }
        const std::size_t time_group_offset = query_time_groups_.size();
        for (std::size_t group_begin = key_begin;
             group_begin < key_end;) {
            const IndexBlock& group_first =
                blocks_[group_begin];
            std::size_t group_end = group_begin + 1;
            while (group_end < key_end) {
                const IndexBlock& next = blocks_[group_end];
                if (next.start != group_first.start ||
                    next.end != group_first.end) {
                    break;
                }
                ++group_end;
            }
            query_time_groups_.push_back(
                {group_first.start, group_first.end,
                 {group_begin, group_end - group_begin}});
            group_begin = group_end;
        }
        query_directory_.push_back(
            {key, IndexPostingRange{
                time_group_offset,
                query_time_groups_.size() - time_group_offset}});
        key_begin = key_end;
    }
    query_directory_.shrink_to_fit();
    query_time_groups_.shrink_to_fit();
}

std::size_t MATIndex::query_directory_bytes() const {
    return query_directory_.size() *
        sizeof(QueryDirectory::value_type) +
        query_time_groups_.size() * sizeof(TimeGroup);
}

void MATIndex::clear() {
    clear_build_buckets();
    blocks_.clear();
    vertex_pool_.clear();
    query_directory_.clear();
    query_time_groups_.clear();
    stats_ = {};
    stats_.mat_ordered_streaming_enabled = ordered_streaming_enabled_;
    stats_.mat_streaming_mode_known = true;
}
}
