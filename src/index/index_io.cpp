#include "index/index_io.h"
#include "index/index_space.h"

#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <set>
#include <stdexcept>
#include <limits>
#include <ostream>
#include <optional>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace tacc {

struct VATRecord {
    VertexID vertex{};
    AttributeMask attributes{};
    Count alpha{}, beta{};
    std::size_t interval_offset{}, interval_count{};
};

struct FormalIndexIOAccess {
    static void restore(VATIndex& index,
                        std::vector<VATRecord> records,
                        std::vector<VATInterval> intervals,
                        IndexBuildStats stats,
                        const TemporalBipartiteGraph& graph) {
        index.clear();
        index.stats_ = stats;
        index.initialize_vertices(graph.num_v_vertices());
        std::sort(records.begin(), records.end(), [](const auto& x, const auto& y) {
            return std::tie(x.vertex, x.attributes, x.alpha, x.beta) <
                   std::tie(y.vertex, y.attributes, y.alpha, y.beta);
        });
        for (const auto& r : records) {
            auto& attrs = index.uv_vertices_[r.vertex];
            if (attrs.empty() || attrs.back().mask != r.attributes)
                attrs.push_back({r.attributes, {}});
            attrs.back().cores.push_back({r.alpha, r.beta, r.interval_offset, r.interval_count});
        }
        index.intervals_ = std::move(intervals);
        finish(index);
    }
    static VATIndex::VertexTree& tree(VATIndex& index) {
        return index.uv_vertices_;
    }
    static void pool(VATIndex& index, std::vector<VATInterval> intervals) {
        index.intervals_ = std::move(intervals);
    }
    static std::vector<VATMaskRange>& masks(VATIndex& index) { return index.query_masks_; }
    static std::vector<VATCoreRange>& cores(VATIndex& index) { return index.query_cores_; }
    static std::vector<VATRecordRef>& refs(VATIndex& index) { return index.query_refs_; }
    static void finish(VATIndex& index, bool legacy_without_directory = true) {
        if (legacy_without_directory) index.build_query_directory();
        index.finalized_ = true;
        index.recount();
        index.validate();
    }

    static void restore(UATIndex& index,
                        std::vector<UATCoreEntry> cores,
                        std::vector<VertexID> vertex_pool,
                        IndexBuildStats stats) {
        index.clear_build_buckets();
        index.finalized_ = true;
        index.cores_ = std::move(cores);
        index.vertex_pool_ = std::move(vertex_pool);
        index.stats_ = stats;
        index.recount();
        index.validate();
    }

    static void restore(MATIndex& index,
                        std::vector<IndexBlock> blocks,
                        std::vector<VertexID> vertex_pool,
                        MATIndex::QueryDirectory query_directory,
                        MATIndex::QueryTimeGroups query_time_groups,
                        IndexBuildStats stats) {
        index.clear_build_buckets();
        index.blocks_ = std::move(blocks);
        index.vertex_pool_ = std::move(vertex_pool);
        index.query_directory_ = std::move(query_directory);
        index.query_time_groups_ = std::move(query_time_groups);
        index.stats_ = stats;
        const std::size_t references = index.vertex_pool_.size();
        std::size_t bytes = index.blocks_.size() * sizeof(IndexBlock) +
            index.vertex_pool_.size() * sizeof(VertexID);
        bytes += index.query_directory_bytes();
        index.stats_.grouped_blocks = index.blocks_.size();
        index.stats_.vertex_references = references;
        index.stats_.approx_bytes = bytes;
    }
};

namespace {

constexpr std::uint64_t kUVBinaryMagic = 0x54414355564f4e31ULL;
constexpr std::uint32_t kVATBinaryVersion = 11;
constexpr std::uint32_t kUATBinaryVersion = 9;
constexpr std::uint32_t kMATBinaryVersion = 13;
constexpr std::size_t kStatsIntegerCount = 29;

enum class IndexKind : std::uint8_t { VAT = 1, UAT = 2, MAT = 3 };

std::uint32_t binary_version(IndexKind kind) {
    switch (kind) {
        case IndexKind::VAT: return kVATBinaryVersion;
        case IndexKind::UAT: return kUATBinaryVersion;
        case IndexKind::MAT: return kMATBinaryVersion;
    }
    return 0;
}

void set_error(std::string* error, const std::string& message) {
    if (error != nullptr) *error = message;
}

template <typename T>
bool write_pod(std::ostream& output, const T& value) {
    static_assert(std::is_trivially_copyable<T>::value,
                  "write_pod requires a trivially copyable type");
    output.write(reinterpret_cast<const char*>(&value), sizeof(T));
    return output.good();
}

template <typename T>
bool read_pod(std::istream& input, T& value) {
    static_assert(std::is_trivially_copyable<T>::value,
                  "read_pod requires a trivially copyable type");
    input.read(reinterpret_cast<char*>(&value), sizeof(T));
    return input.good();
}

bool remaining_bytes_at_least(std::istream& input, std::uint64_t required) {
    const std::streampos current = input.tellg();
    if (current == std::streampos(-1)) return false;
    input.seekg(0, std::ios::end);
    const std::streampos end = input.tellg();
    input.seekg(current);
    if (end == std::streampos(-1) || end < current) return false;
    return required <= static_cast<std::uint64_t>(end - current);
}

bool read_count(std::istream& input,
                std::size_t& count,
                std::size_t minimum_item_bytes) {
    std::uint64_t encoded = 0;
    if (!read_pod(input, encoded) ||
        encoded > static_cast<std::uint64_t>(
            std::numeric_limits<std::size_t>::max())) {
        return false;
    }
    if (minimum_item_bytes != 0 &&
        (encoded > std::numeric_limits<std::uint64_t>::max() /
                       minimum_item_bytes ||
         !remaining_bytes_at_least(
             input, encoded * static_cast<std::uint64_t>(minimum_item_bytes)))) {
        return false;
    }
    count = static_cast<std::size_t>(encoded);
    return true;
}

std::uint64_t mix_hash(std::uint64_t hash, std::uint64_t value) {
    hash ^= value + 0x9e3779b97f4a7c15ULL + (hash << 6U) + (hash >> 2U);
    return hash;
}

std::uint64_t mix_string(std::uint64_t hash, const std::string& value) {
    hash = mix_hash(hash, value.size());
    for (unsigned char ch : value) hash = mix_hash(hash, ch);
    return hash;
}

std::array<std::uint64_t, kStatsIntegerCount> stats_integers(
    const IndexBuildStats& stats) {
    return {{
        stats.attribute_sets_total,
        stats.attribute_sets_no_edge_pruned,
        stats.attribute_sets_equivalent_pruned,
        stats.attribute_sets_nonempty,
        stats.attribute_sets_empty,
        stats.windows,
        stats.decompositions,
        stats.fallback_decompositions,
        stats.full_window_memberships,
        stats.streaming_records_avoided,
        stats.generated_records,
        stats.kept_records,
        stats.attribute_time_removed,
        stats.cross_core_removed,
        stats.grouped_blocks,
        stats.vertex_references,
        stats.approx_bytes,
        stats.inserted_unique_edges,
        stats.deleted_unique_edges,
        stats.dynamic_edge_updates,
        stats.dynamic_peak_bytes,
        stats.affected_u_vertices_total,
        stats.affected_v_vertices_total,
        stats.max_affected_region_vertices,
        stats.changed_u_vertices,
        stats.changed_v_vertices,
        stats.changed_frontier_entries,
        stats.max_changed_vertices_per_bucket,
        stats.nonzero_changed_buckets
    }};
}

bool write_stats(std::ostream& output, const IndexBuildStats& stats) {
    const std::uint8_t method = static_cast<std::uint8_t>(stats.method);
    const std::uint8_t pruning = stats.attribute_pruning_enabled ? 1 : 0;
    if (!write_pod(output, method) || !write_pod(output, pruning) ||
        !write_pod(output, stats.build_seconds) ||
        !write_pod(output, stats.attribute_planning_seconds) ||
        !write_pod(output, stats.decomposition_seconds) ||
        !write_pod(output, stats.dynamic_update_seconds) ||
        !write_pod(output, stats.streaming_seconds) ||
        !write_pod(output, stats.finalize_seconds)) {
        return false;
    }
    for (std::uint64_t value : stats_integers(stats)) {
        if (!write_pod(output, value)) return false;
    }
    return true;
}

bool read_stats(std::istream& input, IndexBuildStats& stats) {
    std::uint8_t method = 0;
    std::uint8_t pruning = 0;
    if (!read_pod(input, method) || !read_pod(input, pruning) ||
        method > static_cast<std::uint8_t>(BuildMethod::Opt) || pruning > 1 ||
        !read_pod(input, stats.build_seconds) ||
        !read_pod(input, stats.attribute_planning_seconds) ||
        !read_pod(input, stats.decomposition_seconds) ||
        !read_pod(input, stats.dynamic_update_seconds) ||
        !read_pod(input, stats.streaming_seconds) ||
        !read_pod(input, stats.finalize_seconds)) {
        return false;
    }
    stats.method = static_cast<BuildMethod>(method);
    stats.attribute_pruning_enabled = pruning != 0;

    std::array<std::uint64_t, kStatsIntegerCount> values{};
    for (std::uint64_t& value : values) {
        if (!read_pod(input, value) ||
            value > static_cast<std::uint64_t>(
                std::numeric_limits<std::size_t>::max())) {
            return false;
        }
    }
    std::size_t i = 0;
    stats.attribute_sets_total = values[i++];
    stats.attribute_sets_no_edge_pruned = values[i++];
    stats.attribute_sets_equivalent_pruned = values[i++];
    stats.attribute_sets_nonempty = values[i++];
    stats.attribute_sets_empty = values[i++];
    stats.windows = values[i++];
    stats.decompositions = values[i++];
    stats.fallback_decompositions = values[i++];
    stats.full_window_memberships = values[i++];
    stats.streaming_records_avoided = values[i++];
    stats.generated_records = values[i++];
    stats.kept_records = values[i++];
    stats.attribute_time_removed = values[i++];
    stats.cross_core_removed = values[i++];
    stats.grouped_blocks = values[i++];
    stats.vertex_references = values[i++];
    stats.approx_bytes = values[i++];
    stats.inserted_unique_edges = values[i++];
    stats.deleted_unique_edges = values[i++];
    stats.dynamic_edge_updates = values[i++];
    stats.dynamic_peak_bytes = values[i++];
    stats.affected_u_vertices_total = values[i++];
    stats.affected_v_vertices_total = values[i++];
    stats.max_affected_region_vertices = values[i++];
    stats.changed_u_vertices = values[i++];
    stats.changed_v_vertices = values[i++];
    stats.changed_frontier_entries = values[i++];
    stats.max_changed_vertices_per_bucket = values[i++];
    stats.nonzero_changed_buckets = values[i++];
    return true;
}

bool read_legacy_uv_tag(std::istream& input) {
    std::uint8_t u=0, v=0;
    return read_pod(input,u) && read_pod(input,v) && u==0 && v==1;
}

bool read_vat_record(std::istream& input, VATRecord& record) {
    std::uint64_t count = 0, offset = 0;
    const bool ok = read_legacy_uv_tag(input) &&
           read_pod(input, record.vertex) &&
           read_pod(input, record.attributes) &&
           read_pod(input, record.alpha) &&
           read_pod(input, record.beta) &&
           read_pod(input, count) && read_pod(input, offset);
    if (!ok || count > std::numeric_limits<std::size_t>::max() ||
        offset > std::numeric_limits<std::size_t>::max()) return false;
    record.interval_count = static_cast<std::size_t>(count);
    record.interval_offset = static_cast<std::size_t>(offset);
    return true;
}

bool write_vat_interval(std::ostream& output, const VATInterval& interval) {
    return write_pod(output, interval.start) &&
           write_pod(output, interval.end);
}

bool read_vat_interval(std::istream& input, VATInterval& interval) {
    return read_pod(input, interval.start) &&
           read_pod(input, interval.end);
}

bool write_size(std::ostream& output, std::size_t value) {
    return write_pod(output, static_cast<std::uint64_t>(value));
}

bool read_size(std::istream& input, std::size_t& value) {
    std::uint64_t encoded = 0;
    if (!read_pod(input, encoded) ||
        encoded > static_cast<std::uint64_t>(
            std::numeric_limits<std::size_t>::max())) {
        return false;
    }
    value = static_cast<std::size_t>(encoded);
    return true;
}

bool write_vat_query_aux(std::ostream& output, const VATIndex& index) {
    if (!write_size(output, index.query_masks().size())) return false;
    for (const auto& m : index.query_masks())
        if (!write_pod(output, m.mask) || !write_size(output, m.offset) || !write_size(output, m.count)) return false;
    if (!write_size(output, index.query_cores().size())) return false;
    for (const auto& c : index.query_cores())
        if (!write_pod(output, c.alpha) || !write_pod(output, c.beta) ||
            !write_size(output, c.offset) || !write_size(output, c.count)) return false;
    if (!write_size(output, index.query_refs().size())) return false;
    for (const auto& r : index.query_refs())
        if (!write_pod(output, r.vertex) || !write_pod(output, r.attribute_slot) || !write_pod(output, r.core_slot)) return false;
    return true;
}

bool read_vat_query_aux(std::istream& input, VATIndex& index, std::size_t leaves) {
    std::size_t n = 0;
    if (!read_count(input, n, 24) || n > leaves) return false;
    auto& masks = FormalIndexIOAccess::masks(index); masks.resize(n);
    for (auto& m : masks)
        if (!read_pod(input, m.mask) || !read_size(input, m.offset) || !read_size(input, m.count)) return false;
    if (!read_count(input, n, 24) || n > leaves) return false;
    auto& cores = FormalIndexIOAccess::cores(index); cores.resize(n);
    for (auto& c : cores)
        if (!read_pod(input, c.alpha) || !read_pod(input, c.beta) ||
            !read_size(input, c.offset) || !read_size(input, c.count)) return false;
    if (!read_count(input, n, 12) || n != leaves) return false;
    auto& refs = FormalIndexIOAccess::refs(index); refs.resize(n);
    for (auto& r : refs)
        if (!read_pod(input, r.vertex) || !read_pod(input, r.attribute_slot) || !read_pod(input, r.core_slot)) return false;
    return true;
}

bool write_shared_block(std::ostream& output, const IndexBlock& block) {
    return write_pod(output, block.attributes) &&
           write_pod(output, block.start) &&
           write_pod(output, block.end) &&
           write_pod(output, block.alpha) &&
           write_pod(output, block.beta) &&
           write_size(output, block.vertex_offset) &&
           write_size(output, block.vertex_count);
}

bool read_shared_block(std::istream& input, IndexBlock& block, bool legacy) {
    return (!legacy || read_legacy_uv_tag(input)) &&
           read_pod(input, block.attributes) &&
           read_pod(input, block.start) &&
           read_pod(input, block.end) &&
           read_pod(input, block.alpha) &&
           read_pod(input, block.beta) &&
           read_size(input, block.vertex_offset) &&
           read_size(input, block.vertex_count);
}

bool write_vertex_pool(std::ostream& output,
                       const std::vector<VertexID>& pool) {
    if (!write_size(output, pool.size())) return false;
    if (!pool.empty()) {
        output.write(reinterpret_cast<const char*>(pool.data()),
                     static_cast<std::streamsize>(
                         pool.size() * sizeof(VertexID)));
    }
    return output.good();
}

bool read_vertex_pool(std::istream& input,
                      std::vector<VertexID>& pool) {
    std::size_t count = 0;
    if (!read_count(input, count, sizeof(VertexID)) ||
        count > pool.max_size()) {
        return false;
    }
    pool.resize(count);
    if (count != 0) {
        input.read(reinterpret_cast<char*>(pool.data()),
                   static_cast<std::streamsize>(
                       count * sizeof(VertexID)));
    }
    return input.good();
}

using CoreQueryKey = std::tuple<Count, Count, AttributeMask>;
using CoreQueryDirectory =
    std::vector<std::pair<CoreQueryKey, IndexPostingRange>>;
using MATQueryKey = std::tuple<AttributeMask, Count>;
using MATQueryDirectory =
    std::vector<std::pair<MATQueryKey, IndexPostingRange>>;

bool read_core_query_aux(std::istream& input,
                         CoreQueryDirectory& directory,
                         bool legacy = false) {
    const std::size_t kMinimumEntryBytes =
        (legacy ? 2 : 0) + 2 * sizeof(Count) + sizeof(AttributeMask) +
        2 * sizeof(std::uint64_t);
    std::size_t directory_count = 0;
    if (!read_count(input, directory_count, kMinimumEntryBytes)) return false;
    directory.resize(directory_count);
    for (auto& entry : directory) {
        Count alpha = 0;
        Count beta = 0;
        AttributeMask attributes = 0;
        if ((legacy && !read_legacy_uv_tag(input)) ||
            !read_pod(input, alpha) || !read_pod(input, beta) ||
            !read_pod(input, attributes) ||
            !read_size(input, entry.second.offset) ||
            !read_size(input, entry.second.count)) {
            return false;
        }
        entry.first = {alpha, beta, attributes};
    }
    return true;
}

bool write_mat_query_aux(
    std::ostream& output,
    const MATQueryDirectory& directory,
    const MATIndex::QueryTimeGroups& time_groups) {
    if (!write_size(output, directory.size())) return false;
    for (const auto& entry : directory) {
        const auto& [attributes, fixed_core] =
            entry.first;
        if (!write_pod(output, attributes) ||
            !write_pod(output, fixed_core) ||
            !write_size(output, entry.second.offset) ||
            !write_size(output, entry.second.count)) {
            return false;
        }
    }
    if (!write_size(output, time_groups.size())) return false;
    for (const MATIndex::TimeGroup& group : time_groups) {
        if (!write_pod(output, group.start) ||
            !write_pod(output, group.end) ||
            !write_size(output, group.postings.offset) ||
            !write_size(output, group.postings.count)) {
            return false;
        }
    }
    return true;
}

bool read_mat_query_aux(std::istream& input,
                        MATQueryDirectory& directory,
                        MATIndex::QueryTimeGroups& time_groups,
                        bool legacy = false) {
    const std::size_t kMinimumEntryBytes =
        (legacy ? 2 : 0) + sizeof(AttributeMask) + sizeof(Count) +
        2 * sizeof(std::uint64_t);
    std::size_t directory_count = 0;
    if (!read_count(input, directory_count, kMinimumEntryBytes)) return false;
    directory.resize(directory_count);
    for (auto& entry : directory) {
        AttributeMask attributes = 0;
        Count fixed_core = 0;
        if ((legacy && !read_legacy_uv_tag(input)) ||
            !read_pod(input, attributes) ||
            !read_pod(input, fixed_core) ||
            !read_size(input, entry.second.offset) ||
            !read_size(input, entry.second.count)) {
            return false;
        }
        entry.first = {attributes, fixed_core};
    }
    constexpr std::size_t kMinimumTimeGroupBytes =
        2 * sizeof(Timestamp) + 2 * sizeof(std::uint64_t);
    std::size_t time_group_count = 0;
    if (!read_count(input, time_group_count, kMinimumTimeGroupBytes)) {
        return false;
    }
    time_groups.resize(time_group_count);
    for (MATIndex::TimeGroup& group : time_groups) {
        if (!read_pod(input, group.start) ||
            !read_pod(input, group.end) ||
            !read_size(input, group.postings.offset) ||
            !read_size(input, group.postings.count)) {
            return false;
        }
    }
    return true;
}

template <typename Directory, typename KeyOf, typename Less>
bool valid_query_aux(const Directory& directory,
                     std::size_t item_count,
                     KeyOf key_of,
                     Less less) {
    const auto count = item_count;
    if (count != item_count || directory.empty() != (count == 0)) {
        return false;
    }
    std::size_t expected_offset = 0;
    for (std::size_t d = 0; d < directory.size(); ++d) {
        const auto& [key, range] = directory[d];
        if (range.count == 0 || range.offset != expected_offset ||
            range.offset > count ||
            range.count > count - range.offset ||
            (d != 0 && !(directory[d - 1].first < key))) {
            return false;
        }
        const std::size_t end = range.offset + range.count;
        for (std::size_t i = range.offset; i < end; ++i) {
            const std::size_t posting = i;
            if (posting >= item_count ||
                key_of(posting) != key ||
                (i != range.offset &&
                 !less(i - 1, posting))) {
                return false;
            }
        }
        expected_offset = end;
    }
    return expected_offset == count;
}

bool valid_vat_query_aux(const std::vector<VATRecord>& records,
                         const CoreQueryDirectory& directory) {
    return valid_query_aux(
        directory, records.size(),
        [&](std::size_t i) {
            const VATRecord& record = records[i];
            return CoreQueryKey{record.alpha, record.beta,
                                record.attributes};
        },
        [&](std::size_t lhs, std::size_t rhs) {
            return records[lhs].vertex < records[rhs].vertex;
        });
}

bool valid_uat_query_aux(const std::vector<IndexBlock>& blocks,
                         const CoreQueryDirectory& directory) {
    return valid_query_aux(
        directory, blocks.size(),
        [&](std::size_t i) {
            const IndexBlock& block = blocks[i];
            return CoreQueryKey{block.alpha, block.beta,
                                block.attributes};
        },
        [&](std::size_t lhs, std::size_t rhs) {
            return std::tie(blocks[lhs].start, blocks[lhs].end) <
                   std::tie(blocks[rhs].start, blocks[rhs].end);
        });
}

bool valid_mat_query_aux(const std::vector<IndexBlock>& blocks,
                         const MATQueryDirectory& directory,
                         const MATIndex::QueryTimeGroups& time_groups) {
    const auto count = blocks.size();
    if (count != blocks.size() ||
        directory.empty() != time_groups.empty() ||
        time_groups.empty() != (count == 0)) {
        return false;
    }
    std::size_t expected_group = 0;
    std::size_t expected_posting = 0;
    for (std::size_t d = 0; d < directory.size(); ++d) {
        const auto& [key, range] = directory[d];
        if (range.count == 0 || range.offset != expected_group ||
            range.offset > time_groups.size() ||
            range.count > time_groups.size() - range.offset ||
            (d != 0 && !(directory[d - 1].first < key))) {
            return false;
        }
        const std::size_t group_end = range.offset + range.count;
        for (std::size_t g = range.offset; g < group_end; ++g) {
            const MATIndex::TimeGroup& group = time_groups[g];
            if (group.start > group.end || group.postings.count == 0 ||
                group.postings.offset != expected_posting ||
                group.postings.offset > count ||
                group.postings.count >
                    count - group.postings.offset ||
                (g != range.offset &&
                 std::tie(time_groups[g - 1].start,
                          time_groups[g - 1].end) >=
                     std::tie(group.start, group.end))) {
                return false;
            }
            Count previous_skyline = 0;
            bool has_previous = false;
            const std::size_t posting_end =
                group.postings.offset + group.postings.count;
            for (std::size_t p = group.postings.offset;
                 p < posting_end; ++p) {
                const std::size_t posting = p;
                if (posting >= blocks.size()) {
                    return false;
                }
                const IndexBlock& block = blocks[posting];
                const Count fixed_core = block.beta;
                const Count skyline_core = block.alpha;
                const MATQueryKey block_key{
                    block.attributes, fixed_core};
                if (block_key != key || block.start != group.start ||
                    block.end != group.end ||
                    (has_previous && skyline_core <= previous_skyline)) {
                    return false;
                }
                    previous_skyline = skyline_core;
                has_previous = true;
            }
            expected_posting = posting_end;
        }
        expected_group = group_end;
    }
    return expected_group == time_groups.size() &&
           expected_posting == count;
}

AttributeMask graph_attribute_universe(const TemporalBipartiteGraph& graph) {
    if (graph.num_attributes() >= 64) {
        return std::numeric_limits<AttributeMask>::max();
    }
    return graph.num_attributes() == 0
        ? 0
        : (AttributeMask{1} << graph.num_attributes()) - 1;
}

bool valid_vat_record(const TemporalBipartiteGraph& graph,
                      const VATRecord& record,
                      const std::vector<VATInterval>& intervals) {
    const AttributeMask universe = graph_attribute_universe(graph);
    if (record.vertex >= graph.num_v_vertices() ||
        record.attributes == 0 ||
        (record.attributes & ~universe) != 0 ||
        record.alpha == 0 || record.beta == 0 ||
        record.interval_count == 0 ||
        record.interval_offset > intervals.size() ||
        record.interval_count > intervals.size() - record.interval_offset) {
        return false;
    }
    const std::size_t end = record.interval_offset + record.interval_count;
    for (std::size_t i = record.interval_offset; i < end; ++i) {
        if (intervals[i].start > intervals[i].end) return false;
        if (i != record.interval_offset &&
            (intervals[i].start <= intervals[i - 1].start ||
             intervals[i].end <= intervals[i - 1].end)) {
            return false;
        }
    }
    return true;
}

bool valid_block(const TemporalBipartiteGraph& graph,
                 const IndexBlock& block,
                 const std::vector<VertexID>& vertex_pool,
                 bool verify_members = true) {
    const AttributeMask universe = graph_attribute_universe(graph);
    if (block.attributes == 0 || (block.attributes & ~universe) != 0 ||
        block.start > block.end || block.alpha == 0 || block.beta == 0 ||
        block.vertex_count == 0 ||
        block.vertex_offset > vertex_pool.size() ||
        block.vertex_count > vertex_pool.size() - block.vertex_offset) {
        return false;
    }
    if (!verify_members) return true;
    const std::size_t end = block.vertex_offset + block.vertex_count;
    for (std::size_t i = block.vertex_offset; i < end; ++i) {
        const VertexID vertex = vertex_pool[i];
        if (vertex >= graph.num_v_vertices()) {
            return false;
        }
    }
    return true;
}

bool write_header(std::ostream& output,
                  IndexKind kind,
                  const TemporalBipartiteGraph& graph,
                  const IndexBuildStats& stats,
                  std::size_t item_count) {
    const std::uint8_t encoded_kind = static_cast<std::uint8_t>(kind);
    const std::uint32_t version = binary_version(kind);
    const std::uint64_t count = item_count;
    const std::uint64_t fingerprint = graph_fingerprint(graph);
    const std::uint64_t magic = kUVBinaryMagic;
    return write_pod(output, magic) &&
           write_pod(output, version) &&
           write_pod(output, encoded_kind) &&
           write_pod(output, fingerprint) &&
           write_stats(output, stats) &&
           write_pod(output, count);
}

bool read_header(std::istream& input,
                 IndexKind expected_kind,
                 const TemporalBipartiteGraph& graph,
                 IndexBuildStats& stats,
                 std::size_t& item_count,
                 std::size_t minimum_item_bytes,
                 std::string* error,
                 std::uint32_t* loaded_version = nullptr) {
    std::uint64_t magic = 0;
    std::uint32_t version = 0;
    std::uint8_t kind = 0;
    std::uint64_t fingerprint = 0;
    if (!read_pod(input, magic) || !read_pod(input, version) ||
        !read_pod(input, kind) || !read_pod(input, fingerprint)) {
        set_error(error, "read_failed");
        return false;
    }
    if (magic != kUVBinaryMagic) {
        set_error(error, "uv_only_index_required"); return false;
    }
    const bool legacy = (expected_kind==IndexKind::VAT && (version==9 || version==10)) ||
                        (expected_kind==IndexKind::UAT && (version==7 || version==8)) ||
                        (expected_kind==IndexKind::MAT && version==12);
    if (version != binary_version(expected_kind) && !legacy) {
        set_error(error, "invalid_format_or_version"); return false;
    }
    if (kind != static_cast<std::uint8_t>(expected_kind)) {
        set_error(error, "index_type_mismatch");
        return false;
    }
    if (fingerprint != graph_fingerprint(graph)) {
        set_error(error, "graph_fingerprint_mismatch");
        return false;
    }
    if (!read_stats(input, stats) ||
        !read_count(input, item_count, minimum_item_bytes)) {
        set_error(error, "read_failed");
        return false;
    }
    stats.scope = IndexScope::UV;
    if (loaded_version != nullptr) *loaded_version = version;
    return true;
}

template <typename Index>
bool finish_scoped_load(std::istream& input, const Index& index,
                        std::string* error, std::optional<std::uint64_t> computed_hash = std::nullopt,
                        bool verify_content = true) {
    std::uint64_t stored_hash = 0;
    if (!read_pod(input, stored_hash) ||
        (verify_content && stored_hash != (computed_hash ? *computed_hash : formal_index_hash(index))) ||
        input.peek() != std::char_traits<char>::eof()) {
        set_error(error, "uv_content_checksum_or_trailing_data_mismatch");
        return false;
    }
    return true;
}

std::uint64_t hash_vat_record(std::uint64_t hash,
                              const VATRecord& record) {
    hash = mix_hash(mix_hash(hash, 0), 1);
    hash = mix_hash(hash, record.vertex);
    hash = mix_hash(hash, record.attributes);
    hash = mix_hash(hash, record.alpha);
    hash = mix_hash(hash, record.beta);
    return mix_hash(hash, record.interval_count);
}

std::uint64_t hash_block(std::uint64_t hash, const IndexBlock& block,
                         const std::vector<VertexID>& vertex_pool, bool legacy = false) {
    if (legacy) hash = mix_hash(mix_hash(hash, 0), 1);
    hash = mix_hash(hash, block.attributes);
    hash = mix_hash(hash, block.start);
    hash = mix_hash(hash, block.end);
    hash = mix_hash(hash, block.alpha);
    hash = mix_hash(hash, block.beta);
    hash = mix_hash(hash, block.vertex_count);
    const std::size_t end = block.vertex_offset + block.vertex_count;
    for (std::size_t i = block.vertex_offset; i < end; ++i) {
        hash = mix_hash(hash, vertex_pool[i]);
    }
    return hash;
}

template <typename Index>
std::uint64_t legacy_uv_blocks_hash(const Index& index) {
    const auto kind = std::is_same_v<Index, UATIndex> ? IndexKind::UAT : IndexKind::MAT;
    auto hash = mix_hash(1469598103934665603ULL, static_cast<std::uint8_t>(kind));
    hash = mix_hash(hash, index.blocks().size());
    for (const auto& block : index.blocks())
        hash = hash_block(hash, block, index.vertex_pool(), true);
    return hash;
}

void write_attributes(std::ostream& output,
                      const TemporalBipartiteGraph& graph,
                      AttributeMask mask) {
    output << '{';
    bool first = true;
    for (AttributeID attribute = 0;
         attribute < graph.num_attributes() && attribute < 64; ++attribute) {
        if ((mask & (AttributeMask{1} << attribute)) == 0) continue;
        if (!first) output << ',';
        output << graph.attribute_name(attribute);
        first = false;
    }
    output << '}';
}

void write_generated_at(std::ostream& output) {
    const auto now = std::time(nullptr);
    std::tm local {};
    if (!localtime_r(&now, &local)) throw std::runtime_error("cannot read local time");
    output << "Generated at: " << std::put_time(&local, "%Y-%m-%d %H:%M:%S") << '\n';
}

void write_text_header(std::ostream& output,
                       const TemporalBipartiteGraph& graph,
                       const MATIndex& index,
                       const std::string& data_directory) {
    const auto& stats = index.stats();
    const auto space = formal_index_space(index);
    auto data = std::filesystem::path(data_directory).lexically_normal();
    if (!data.empty() && data.filename().empty()) data = data.parent_path();
    output << std::setprecision(12) << "TACC MAT single-sided attribute index log\n";
    write_generated_at(output);
    output << "Dataset: " << (data.empty() ? "unknown" : data.filename().string())
           << "\nConstruction algorithm: "
           << (stats.method == BuildMethod::Basic ? "basic (fresh snapshot decomposition)"
                                                  : "opt (dynamic bicore maintenance)") << '\n';
    if (!data.empty()) {
        output << "Label files: " << (data / "graph_ulabel.txt").string() << ", "
               << (data / "graph_vlabel.txt").string()
               << "\nEdge file: " << (data / "graph_timestamp.txt").string() << '\n';
    }
    output << "Graph summary: |U|=" << graph.num_u_vertices() << ", |V|=" << graph.num_v_vertices()
           << ", temporal_events=" << graph.num_edges() << ", unique_pairs=" << graph.num_unique_edges()
           << ", |T|=" << graph.time_points().size()
           << "\nIndex summary: index_mb=" << space_mb(space.index_bytes())
           << ", blocks=" << index.blocks().size() << ", physical_vertex_refs=" << index.vertex_pool().size()
           << "\nSpace breakdown (MB): root=" << space_mb(space.root_bytes)
           << ", blocks_and_pool=" << space_mb(space.body_bytes())
           << ", query_aux=" << space_mb(space.query_aux_bytes())
           << "\nContent hash: " << formal_index_hash(index)
           << "\nSemantics: constrain U attributes and project onto V; stored vertices are candidates for the joint-attribute query, not the final answer.\n"
           << "Storage: vertex blocks with shared vertex lists; persisted attribute/beta directory and time groups. Lists below are expanded for readability.\n"
           << "Each block stores an alpha threshold at a fixed beta; it also supports queries with a smaller alpha.\n"
           << "Core (1,1) is excluded from stored records. Original V vertex IDs shown.\n\n[ATTRIBUTES]\n";
    for (AttributeID a = 0; a < graph.num_attributes() && a < 64; ++a)
        output << "attribute_id=" << a << " name=" << graph.attribute_name(a)
               << " bit=" << (AttributeMask{1} << a) << '\n';
    output << "\n[TIME_POINTS]";
    for (const auto time : graph.time_points()) output << ' ' << time;
    output << "\n\n[MAT] (U attribute set, time interval, alpha threshold, beta) -> V vertex set\n\n";
}

bool finish_text(std::ostream& output, std::string* error) {
    if (!output.good()) {
        set_error(error, "write_failed");
        return false;
    }
    return true;
}

}

std::uint64_t graph_fingerprint(const TemporalBipartiteGraph& graph) {
    if (const auto cached = std::atomic_load(&graph.fingerprint_)) return *cached;
    std::uint64_t hash = 1469598103934665603ULL;
    hash = mix_hash(hash, graph.num_u_vertices());
    hash = mix_hash(hash, graph.num_v_vertices());
    hash = mix_hash(hash, graph.num_edges());
    hash = mix_hash(hash, graph.num_attributes());
    for (AttributeID id = 0; id < graph.num_attributes(); ++id) {
        hash = mix_string(hash, graph.attribute_name(id));
    }
    for (VertexID u = 0; u < graph.num_u_vertices(); ++u) {
        hash = mix_hash(hash, graph.original_u_id(u));
        hash = mix_hash(hash, graph.u_attribute(u));
    }
    for (VertexID v = 0; v < graph.num_v_vertices(); ++v) {
        hash = mix_hash(hash, graph.original_v_id(v));
        hash = mix_hash(hash, graph.v_attribute(v));
    }
    for (const TemporalEdge& edge : graph.edges()) {
        hash = mix_hash(hash, edge.u);
        hash = mix_hash(hash, edge.v);
        hash = mix_hash(hash, edge.timestamp);
    }
    std::atomic_store(&graph.fingerprint_, std::make_shared<const std::uint64_t>(hash));
    return hash;
}

static std::uint64_t legacy_vat_v10_hash(const VATIndex& index) {
    std::uint64_t hash = mix_hash(1469598103934665603ULL,
                                  static_cast<std::uint8_t>(IndexKind::VAT));
    hash = mix_hash(hash, static_cast<std::uint8_t>(index.stats().scope));
    {
        hash = mix_hash(hash, 1);
        const auto& tree = index.vertices();
        hash = mix_hash(hash, tree.size());
        for (const auto& attrs : tree) {
            hash = mix_hash(hash, attrs.size());
            for (const auto& a : attrs) {
                hash = mix_hash(hash, a.mask);
                hash = mix_hash(hash, a.cores.size());
                for (const auto& c : a.cores) {
                    hash = mix_hash(hash, c.alpha); hash = mix_hash(hash, c.beta);
                    hash = mix_hash(hash, c.interval_offset); hash = mix_hash(hash, c.interval_count);
                }
            }
        }
    }
    hash = mix_hash(hash, index.intervals().size());
    for (const auto& t : index.intervals()) {
        hash = mix_hash(hash, t.start); hash = mix_hash(hash, t.end);
    }
    return hash;
}

std::uint64_t formal_index_hash(const VATIndex& index) {
    auto hash = legacy_vat_v10_hash(index);
    hash = mix_hash(hash, index.query_masks().size());
    for (const auto& m : index.query_masks()) {
        hash = mix_hash(hash, m.mask); hash = mix_hash(hash, m.offset); hash = mix_hash(hash, m.count);
    }
    hash = mix_hash(hash, index.query_cores().size());
    for (const auto& c : index.query_cores()) {
        hash = mix_hash(hash, c.alpha); hash = mix_hash(hash, c.beta);
        hash = mix_hash(hash, c.offset); hash = mix_hash(hash, c.count);
    }
    hash = mix_hash(hash, index.query_refs().size());
    for (const auto& r : index.query_refs()) {
        hash = mix_hash(hash, r.vertex); hash = mix_hash(hash, r.attribute_slot); hash = mix_hash(hash, r.core_slot);
    }
    return hash;
}

std::uint64_t formal_index_hash(const UATIndex& index) {
    std::uint64_t hash = mix_hash(1469598103934665603ULL,
                                  static_cast<std::uint8_t>(IndexKind::UAT));
    hash = mix_hash(hash,index.cores().size());
    for (const auto& c : index.cores()) {
        hash=mix_hash(hash,c.alpha); hash=mix_hash(hash,c.beta);
        hash=mix_hash(hash,c.attributes.size());
        for (const auto& a : c.attributes) {
            hash=mix_hash(hash,a.mask); hash=mix_hash(hash,a.starts.size());
            for (const auto& s : a.starts) {
                hash=mix_hash(hash,s.start); hash=mix_hash(hash,s.ends.size());
                for (const auto& e : s.ends) {
                    hash=mix_hash(hash,e.end); hash=mix_hash(hash,e.vertex_offset);
                    hash=mix_hash(hash,e.vertex_count);
                }
            }
        }
    }
    hash=mix_hash(hash,index.vertex_pool().size());
    for (auto v : index.vertex_pool()) hash=mix_hash(hash,v);
    return hash;
}

std::uint64_t formal_index_hash(const MATIndex& index) {
    std::uint64_t hash = mix_hash(1469598103934665603ULL,
                                  static_cast<std::uint8_t>(IndexKind::MAT));
    hash = mix_hash(hash, index.blocks().size());
    for (const IndexBlock& block : index.blocks()) {
        hash = hash_block(hash, block, index.vertex_pool());
    }
    return hash;
}

bool read_index_scope(const std::string& path, IndexScope& scope,
                      std::string* error) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        set_error(error, "open_failed");
        return false;
    }
    std::uint64_t magic = 0;
    if (!read_pod(input, magic) ||
        magic != kUVBinaryMagic) {
        set_error(error, "invalid_format_or_version");
        return false;
    }
    scope = IndexScope::UV;
    return true;
}

bool save_index_binary(const std::string& path,
                       const TemporalBipartiteGraph& graph,
                       const VATIndex& index,
                       std::string* error) {
    try { index.validate(); } catch (const std::exception& e) {
        set_error(error, e.what()); return false;
    }
    std::ofstream output(path, std::ios::binary);
    if (!output) {
        set_error(error, "open_failed");
        return false;
    }
    bool ok = write_header(output, IndexKind::VAT, graph, index.stats(),
                           index.stats().grouped_blocks);
    {
        const auto& tree = index.vertices();
        if (ok) ok = write_pod(output, std::uint64_t(tree.size()));
        for (const auto& attrs : tree) {
            if (ok) ok = write_pod(output, std::uint64_t(attrs.size()));
            for (const auto& a : attrs) {
                if (ok) ok = write_pod(output, a.mask) && write_pod(output, std::uint64_t(a.cores.size()));
                for (const auto& c : a.cores) if (ok)
                    ok = write_pod(output, c.alpha) && write_pod(output, c.beta) &&
                         write_pod(output, std::uint64_t(c.interval_offset)) && write_pod(output, std::uint64_t(c.interval_count));
            }
        }
    }
    if (ok) ok = write_pod(output, std::uint64_t(index.intervals().size()));
    for (const VATInterval& interval : index.intervals()) {
        if (ok) ok = write_vat_interval(output, interval);
    }
    if (ok) ok = write_vat_query_aux(output, index);
    if (ok) ok = write_pod(output, formal_index_hash(index));
    output.flush();
    ok = ok && output.good();
    if (!ok) set_error(error, "write_failed");
    return ok;
}

bool save_index_binary(const std::string& path,
                       const TemporalBipartiteGraph& graph,
                       const UATIndex& index,
                       std::string* error) {
    if (!index.ready()) {set_error(error,"uat_not_finalized");return false;}
    std::ofstream output(path, std::ios::binary);
    if (!output) {
        set_error(error, "open_failed");
        return false;
    }
    bool ok = write_header(output, IndexKind::UAT, graph, index.stats(), index.cores().size()) &&
              write_vertex_pool(output, index.vertex_pool());
    for (const auto& c : index.cores()) {
        if (ok) ok=write_pod(output,c.alpha) && write_pod(output,c.beta) && write_size(output,c.attributes.size());
        for (const auto& a : c.attributes) {
            if (ok) ok=write_pod(output,a.mask) && write_size(output,a.starts.size());
            for (const auto& s : a.starts) {
                if (ok) ok=write_pod(output,s.start) && write_size(output,s.ends.size());
                for (const auto& e : s.ends)
                    if (ok) ok=write_pod(output,e.end) && write_size(output,e.vertex_offset) && write_size(output,e.vertex_count);
            }
        }
    }
    if (ok)
        ok = write_pod(output, formal_index_hash(index));
    output.flush();
    ok = ok && output.good();
    if (!ok) set_error(error, "write_failed");
    return ok;
}

bool save_index_binary(const std::string& path,
                       const TemporalBipartiteGraph& graph,
                       const MATIndex& index,
                       std::string* error) {
    std::ofstream output(path, std::ios::binary);
    if (!output) {
        set_error(error, "open_failed");
        return false;
    }
    bool ok = write_header(output, IndexKind::MAT, graph, index.stats(),
                           index.blocks().size()) &&
        write_vertex_pool(output, index.vertex_pool());
    for (const IndexBlock& block : index.blocks()) {
        if (ok) ok = write_shared_block(output, block);
    }
    if (ok) {
        ok = write_mat_query_aux(output, index.query_directory(),
                                 index.query_time_groups());
    }
    if (ok)
        ok = write_pod(output, formal_index_hash(index));
    output.flush();
    ok = ok && output.good();
    if (!ok) set_error(error, "write_failed");
    return ok;
}

bool load_index_binary(const std::string& path,
                       const TemporalBipartiteGraph& graph,
                       VATIndex& index,
                       std::string* error) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        set_error(error, "open_failed");
        return false;
    }
    IndexBuildStats stats;
    std::size_t count = 0;
    std::uint32_t version = 0;
    if (!read_header(input, IndexKind::VAT, graph, stats, count,
                     24, error, &version)) {
        return false;
    }
    if (version == 10 || version == 11) {
        VATIndex loaded;
        loaded.mutable_stats() = stats;
        loaded.initialize_vertices(graph.num_v_vertices());
        std::size_t cores = 0, logical = 0;
        {
            auto& tree = FormalIndexIOAccess::tree(loaded);
            std::size_t n = 0;
            if (!read_count(input, n, 8) || n != tree.size()) {
                set_error(error, "invalid_vertex_dimension"); return false;
            }
            for (auto& attrs : tree) {
                std::size_t na = 0;
                if (!read_count(input, na, 16)) { set_error(error, "invalid_attribute_count"); return false; }
                attrs.resize(na);
                for (auto& a : attrs) {
                    std::size_t nc = 0;
                    if (!read_pod(input, a.mask) || !a.mask || (a.mask & ~graph_attribute_universe(graph)) ||
                        !read_count(input, nc, 24) || nc > count - cores) {
                        set_error(error, "invalid_attribute_or_core_count"); return false;
                    }
                    a.cores.resize(nc); cores += nc;
                    for (auto& c : a.cores) {
                        if (!read_pod(input, c.alpha) || !read_pod(input, c.beta) ||
                            !read_count(input, c.interval_offset, 0) || !read_count(input, c.interval_count, 0) ||
                            c.interval_count > std::numeric_limits<std::size_t>::max() - logical) {
                            set_error(error, "invalid_core"); return false;
                        }
                        logical += c.interval_count;
                    }
                }
            }
        }
        std::size_t np = 0;
        if (!read_count(input, np, 8) || cores != count || cores != stats.grouped_blocks ||
            logical != stats.kept_records || np != stats.vertex_references) {
            set_error(error, "invalid_nested_totals"); return false;
        }
        std::vector<VATInterval> pool(np);
        for (auto& t : pool) if (!read_vat_interval(input, t) ||
            !std::binary_search(graph.time_points().begin(), graph.time_points().end(), t.start) ||
            !std::binary_search(graph.time_points().begin(), graph.time_points().end(), t.end)) {
            set_error(error, "invalid_interval"); return false;
        }
        FormalIndexIOAccess::pool(loaded, std::move(pool));
        if (version == 11 && !read_vat_query_aux(input, loaded, count)) {
            set_error(error, "invalid_vat_query_directory_encoding"); return false;
        }
        try { FormalIndexIOAccess::finish(loaded, version == 10); }
        catch (const std::exception& e) { set_error(error, e.what()); return false; }
        if (version == 11 && loaded.stats().approx_bytes != stats.approx_bytes) {
            set_error(error, "vat_structural_size_mismatch"); return false;
        }
        std::uint64_t stored_hash = 0;
        if (!read_pod(input, stored_hash) ||
            stored_hash != (version == 10 ? legacy_vat_v10_hash(loaded) : formal_index_hash(loaded)) ||
            input.peek() != std::char_traits<char>::eof()) {
            set_error(error, "vat_content_checksum_or_trailing_data_mismatch"); return false;
        }
        index = std::move(loaded);
        return true;
    }
    std::size_t interval_count = 0;
    if (!read_count(input, interval_count, 2 * sizeof(Timestamp))) {
        set_error(error, "read_failed");
        return false;
    }
    std::vector<VATRecord> records(count);
    for (VATRecord& record : records) {
        if (!read_vat_record(input, record)) {
            set_error(error, "invalid_record");
            return false;
        }
    }
    std::vector<VATInterval> intervals(interval_count);
    for (VATInterval& interval : intervals) {
        if (!read_vat_interval(input, interval)) {
            set_error(error, "invalid_interval");
            return false;
        }
    }
    std::size_t expected_offset = 0, logical_intervals = 0;
    for (VATRecord& record : records) {
        if (!valid_vat_record(graph, record, intervals) ||
            record.interval_count > std::numeric_limits<std::size_t>::max() - logical_intervals) {
            set_error(error, "invalid_record");
            return false;
        }
        expected_offset = std::max(expected_offset, record.interval_offset + record.interval_count);
        logical_intervals += record.interval_count;
    }
    if (expected_offset != intervals.size() || stats.kept_records != logical_intervals ||
        stats.vertex_references != intervals.size()) {
        set_error(error, "invalid_interval_pool");
        return false;
    }
    CoreQueryDirectory query_directory;
    if (!read_core_query_aux(input, query_directory, true) ||
        !valid_vat_query_aux(records, query_directory)) {
        set_error(error, "invalid_query_aux");
        return false;
    }
    VATIndex loaded;
    {
        std::uint64_t hash = mix_hash(1469598103934665603ULL, static_cast<std::uint8_t>(IndexKind::VAT));
        hash = mix_hash(hash, records.size()); hash = mix_hash(hash, logical_intervals);
        for (const auto& r : records) {
            hash = hash_vat_record(hash, r);
            for (std::size_t i = r.interval_offset; i < r.interval_offset + r.interval_count; ++i) {
                hash = mix_hash(hash, intervals[i].start); hash = mix_hash(hash, intervals[i].end);
            }
        }
        std::uint64_t stored_hash = 0;
        if (!read_pod(input, stored_hash) || stored_hash != hash || input.peek() != std::char_traits<char>::eof()) {
            set_error(error, "uv_content_checksum_or_trailing_data_mismatch"); return false;
        }
    }
    CoreQueryDirectory().swap(query_directory);
    try {
    FormalIndexIOAccess::restore(loaded, std::move(records),
                                 std::move(intervals),
                                 stats, graph);
    } catch (const std::exception& e) { set_error(error, e.what()); return false; }
    index = std::move(loaded);
    return true;
}

bool load_index_binary(const std::string& path,
                       const TemporalBipartiteGraph& graph,
                       UATIndex& index,
                       std::string* error) try {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        set_error(error, "open_failed");
        return false;
    }
    IndexBuildStats stats;
    std::size_t count = 0;
    std::uint32_t version = 0;
    if (!read_header(input, IndexKind::UAT, graph, stats, count,
                     2*sizeof(Count)+sizeof(std::uint64_t), error, &version)) {
        return false;
    }
    std::vector<VertexID> vertex_pool;
    if (!read_vertex_pool(input, vertex_pool) || vertex_pool.size() != stats.vertex_references) {
        set_error(error, "invalid_vertex_pool");
        return false;
    }
    std::vector<UATCoreEntry> cores;
    if (version==9) {
        cores.resize(count);
        for (auto& c : cores) {
            std::size_t n=0;
            if (!read_pod(input,c.alpha) || !read_pod(input,c.beta) || !read_count(input,n,16)) {
                set_error(error,"invalid_uat_core");return false;
            }
            c.attributes.resize(n);
            for (auto& a : c.attributes) {
                if (!read_pod(input,a.mask) || !read_count(input,n,12)) {
                    set_error(error,"invalid_uat_attribute");return false;
                }
                a.starts.resize(n);
                for (auto& s : a.starts) {
                    if (!read_pod(input,s.start) || !read_count(input,n,20)) {
                        set_error(error,"invalid_uat_start");return false;
                    }
                    s.ends.resize(n);
                    for (auto& e : s.ends)
                        if (!read_pod(input,e.end) || !read_size(input,e.vertex_offset) || !read_size(input,e.vertex_count)) {
                            set_error(error,"invalid_uat_end");return false;
                        }
                }
            }
        }
    } else {
        if (count>std::numeric_limits<std::uint64_t>::max()/40 ||
            !remaining_bytes_at_least(input,std::uint64_t(count)*40)) {
            set_error(error,"invalid_legacy_uat_count");return false;
        }
        std::vector<IndexBlock> blocks(count);
        auto hash=mix_hash(1469598103934665603ULL,static_cast<std::uint8_t>(IndexKind::UAT));
        hash=mix_hash(hash,count);
        for (auto& b : blocks) {
            if (!read_shared_block(input,b,version==7) || !valid_block(graph,b,vertex_pool)) {
                set_error(error,"invalid_legacy_uat_block");return false;
            }
            hash=hash_block(hash,b,vertex_pool,version==7);
        }
        CoreQueryDirectory directory;
        if (!read_core_query_aux(input,directory,version==7) || !valid_uat_query_aux(blocks,directory)) {
            set_error(error,"invalid_legacy_uat_directory");return false;
        }
        std::uint64_t stored_hash=0;
        if (!read_pod(input,stored_hash) || stored_hash!=hash || input.peek()!=std::char_traits<char>::eof()) {
            set_error(error,"uv_content_checksum_or_trailing_data_mismatch");return false;
        }
        for (const auto& b : blocks) {
            if (cores.empty() || cores.back().alpha!=b.alpha || cores.back().beta!=b.beta)
                cores.push_back({b.alpha,b.beta,{}});
            auto& c=cores.back();
            if (c.attributes.empty() || c.attributes.back().mask!=b.attributes)
                c.attributes.push_back({b.attributes,{}});
            auto& a=c.attributes.back();
            if (a.starts.empty() || a.starts.back().start!=b.start) a.starts.push_back({b.start,{}});
            a.starts.back().ends.push_back({b.end,b.vertex_offset,b.vertex_count});
        }
    }
    UATIndex loaded;
    FormalIndexIOAccess::restore(loaded,std::move(cores),std::move(vertex_pool),stats);
    if (loaded.stats().kept_records!=stats.kept_records || loaded.stats().grouped_blocks!=stats.grouped_blocks ||
        (version==9 && loaded.stats().approx_bytes!=stats.approx_bytes)) {
        set_error(error,"invalid_uat_statistics");return false;
    }
    bool valid=true;
    loaded.for_each_block([&](const IndexBlock& b) {
        valid=valid && valid_block(graph,b,loaded.vertex_pool()) &&
              b.start>=graph.min_timestamp() && b.end<=graph.max_timestamp();
    });
    if (!valid) {set_error(error,"invalid_uat_graph_membership");return false;}
    if (version==9 && !finish_scoped_load(input,loaded,error)) return false;
    index = std::move(loaded);
    return true;
} catch (const std::exception& e) {
    set_error(error,std::string("invalid_uat_index: ")+e.what());return false;
}

bool load_index_binary(const std::string& path,
                       const TemporalBipartiteGraph& graph,
                       MATIndex& index,
                       std::string* error,
                       bool verify_content) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        set_error(error, "open_failed");
        return false;
    }
    IndexBuildStats stats;
    std::size_t count = 0;
    std::uint32_t version = 0;
    constexpr std::size_t kMinimumBlockBytes =
        2 + sizeof(AttributeMask) + 4 * sizeof(std::uint32_t) +
        2 * sizeof(std::uint64_t);
    if (!read_header(input, IndexKind::MAT, graph, stats, count,
                     kMinimumBlockBytes, error, &version)) {
        return false;
    }
    std::vector<VertexID> vertex_pool;
    if (!read_vertex_pool(input, vertex_pool) ||
        vertex_pool.size() != stats.vertex_references) {
        set_error(error, "invalid_vertex_pool");
        return false;
    }
    if (std::any_of(vertex_pool.begin(), vertex_pool.end(), [&](VertexID vertex) {
            return vertex >= graph.num_v_vertices();
        })) {
        set_error(error, "invalid_vertex_pool_member");
        return false;
    }
    std::vector<IndexBlock> blocks(count);
    for (IndexBlock& block : blocks) {
        if (!read_shared_block(input, block, version == 12) ||
            !valid_block(graph, block, vertex_pool, false)) {
            set_error(error, "invalid_block");
            return false;
        }
    }
    MATQueryDirectory query_directory;
    MATIndex::QueryTimeGroups query_time_groups;
    if (!read_mat_query_aux(input, query_directory, query_time_groups, version == 12) ||
        !valid_mat_query_aux(blocks, query_directory, query_time_groups)) {
        set_error(error, "invalid_query_aux");
        return false;
    }
    MATIndex loaded;
    FormalIndexIOAccess::restore(loaded, std::move(blocks),
                                 std::move(vertex_pool),
                                 std::move(query_directory),
                                 std::move(query_time_groups),
                                 stats);
    if (!finish_scoped_load(input, loaded, error,
            verify_content && version == 12 ? std::make_optional(legacy_uv_blocks_hash(loaded)) : std::nullopt,
            verify_content)) return false;
    index = std::move(loaded);
    return true;
}

static bool save_vat_text(const std::string& path,
                     const TemporalBipartiteGraph& graph,
                     const VATIndex& index,
                     std::string* error,
                     const std::string& data_directory,
                     std::uint64_t build_peak_rss_bytes,
                     bool expanded_intervals) {
    std::ofstream output(path);
    if (!output) {
        set_error(error, "open_failed");
        return false;
    }
    try { index.validate(); } catch (const std::exception& e) { set_error(error, e.what()); return false; }
    const auto& s = index.stats();
    const auto now = std::time(nullptr); std::tm local{};
    localtime_r(&now, &local);
    std::set<AttributeID> ua, va;
    for (VertexID v = 0; v < graph.num_u_vertices(); ++v) ua.insert(graph.u_attribute(v));
    for (VertexID v = 0; v < graph.num_v_vertices(); ++v) va.insert(graph.v_attribute(v));
    const auto dataset = data_directory.empty() ? std::string("unknown") :
        data_directory.substr(data_directory.find_last_of("/\\") + 1);
    output << std::setprecision(12)
           << (expanded_intervals ? "TACC VAT expanded interval view\nGenerated at: " : "TACC VAT index log\nGenerated at: ")
           << std::put_time(&local, "%Y-%m-%d %H:%M:%S") << "\nDataset: " << dataset
           << "\nConstruction algorithm: " << (s.method == BuildMethod::Opt ? "opt" : "basic");
    for (const auto& file : {std::make_pair("U attribute file", "graph_ulabel.txt"),
                            std::make_pair("V attribute file", "graph_vlabel.txt"),
                            std::make_pair("Edge file", "graph_timestamp.txt")})
        output << '\n' << file.first << ": " << (data_directory.empty() ? "unknown" : data_directory + '/' + file.second);
    output << "\nSpan mode: active\nDirection: " << "U -> V"
           << "\nGraph summary: |U|=" << graph.num_u_vertices() << ", |V|=" << graph.num_v_vertices()
           << ", |E|=" << graph.num_edges() << ", |T|=" << graph.time_points().size()
           << ", |L_U|=" << ua.size() << ", |L_V|=" << va.size()
           << "\nIndex summary: status=ok, build_s=" << s.build_seconds
           << ", index_mb=" << space_mb(formal_index_space(index).index_bytes()) << ", build_peak_rss_mb=";
    if (build_peak_rss_bytes) output << space_mb(build_peak_rss_bytes); else output << "NA";
    output << ", core_entries=" << s.grouped_blocks << ", logical_intervals=" << s.kept_records
           << ", physical_intervals=" << index.intervals().size()
           << ", directory_refs=" << index.query_refs().size() << "\n\n[ATTRIBUTES]\n";
    for (AttributeID a = 0; a < graph.num_attributes() && a < 64; ++a)
        output << "attribute_id=" << a << " attribute_name=" << graph.attribute_name(a)
               << " bit=" << (AttributeMask{1} << a) << '\n';
    if (expanded_intervals)
        output << "\nView: expanded intervals for reading only; actual storage remains the shared pool.\n"
               << "record_id is display-only vertex-major numbering, not an extra stored index.\n"
               << "Directory record_offset/count refer to positions in that attribute's record_ids list.\n";
    output << '\n' << std::string(80, '=')
           << (expanded_intervals
               ? "\nVAT expanded view: vertex -> attribute set (MAS) -> (alpha,beta) -> MAI intervals\n"
               : "\nVAT index: vertex -> attribute set (MAS) -> (alpha,beta) -> MAI interval pool range\n")
           << std::string(80, '=') << "\ncore_entries = " << s.grouped_blocks
           << "\nlogical_intervals = " << s.kept_records << "\npooled_intervals = " << index.intervals().size()
           << "\nStorage: shared complete interval lists; "
           << "no direction dimension; persisted sparse lookup directory.\n";
    if (expanded_intervals) {
        const auto& tree = index.vertices();
        std::vector<std::vector<std::size_t>> first_record(tree.size());
        std::size_t record_id = 0;
        output << "N = " << tree.size() << '\n';
        for (std::size_t v = 0; v < tree.size(); ++v) {
            output << "\nvertex " << graph.original_v_id(VertexID(v)) << ":\n";
            first_record[v].resize(tree[v].size());
            for (std::size_t a = 0; a < tree[v].size(); ++a) {
                const auto& attr = tree[v][a]; first_record[v][a] = record_id;
                output << "  MAS=" << attr.mask << ' '; write_attributes(output, graph, attr.mask); output << '\n';
                for (const auto& c : attr.cores) {
                    output << "    core=(" << c.alpha << ',' << c.beta << ") record_id=" << record_id++ << ": ";
                    for (std::size_t i = 0; i < c.interval_count; ++i) {
                        if (i) output << ", ";
                        const auto& t = index.intervals()[c.interval_offset + i];
                        output << '[' << t.start << ',' << t.end << ']';
                    }
                    output << '\n';
                }
            }
        }
        output << "\nlookup_by_attribute\n";
        for (const auto& m : index.query_masks()) {
            output << "  "; write_attributes(output, graph, m.mask);
            output << " mask=" << m.mask << " record_ids=[";
            const auto first = index.query_cores()[m.offset].offset;
            const auto& last = index.query_cores()[m.offset + m.count - 1];
            for (std::size_t i = first; i < last.offset + last.count; ++i) {
                if (i != first) output << ',';
                const auto& r = index.query_refs()[i];
                output << first_record[r.vertex][r.attribute_slot] + r.core_slot;
            }
            output << "]\n";
        }
        output << "\ncore_start_by_attribute\n";
        for (const auto& m : index.query_masks()) {
            output << "  "; write_attributes(output, graph, m.mask);
            output << " mask=" << m.mask << '\n';
            const auto first = index.query_cores()[m.offset].offset;
            for (std::size_t i = m.offset; i < m.offset + m.count; ++i) {
                const auto& c = index.query_cores()[i];
                output << "    core=(" << c.alpha << ',' << c.beta << ") record_offset="
                       << c.offset - first << " record_count=" << c.count << '\n';
            }
        }
        return finish_text(output, error);
    }
    {
        const auto& tree = index.vertices();
        output << "N = " << tree.size() << '\n';
        for (std::size_t v = 0; v < tree.size(); ++v) {
            output << '\n' << 'V' << '='
                   << graph.original_v_id(VertexID(v)) << " internal=" << v << '\n';
            for (const auto& a : tree[v]) {
                output << "  MAS=" << a.mask << ' '; write_attributes(output, graph, a.mask); output << '\n';
                for (const auto& c : a.cores)
                    output << "    core=(" << c.alpha << ',' << c.beta << ") interval_offset=" << c.interval_offset
                           << " interval_count=" << c.interval_count << '\n';
            }
        }
    }
    output << "\n[interval_pool]\n";
    for (std::size_t i = 0; i < index.intervals().size(); ++i)
        output << "interval_pool[" << i << "] = (" << index.intervals()[i].start << ',' << index.intervals()[i].end << ")\n";
    output << "\nlookup_by_attribute\n";
    for (const auto& m : index.query_masks()) {
        output << "  "; write_attributes(output, graph, m.mask);
        output << " mask=" << m.mask << " core_offset=" << m.offset << " core_count=" << m.count << '\n';
    }
    output << "\ncore_start_by_attribute\n";
    for (std::size_t i = 0; i < index.query_cores().size(); ++i) {
        const auto& c = index.query_cores()[i];
        output << "  cores[" << i << "] core=(" << c.alpha << ',' << c.beta
               << ") ref_offset=" << c.offset << " ref_count=" << c.count << '\n';
    }
    output << "\nrecord_refs\n";
    for (std::size_t i = 0; i < index.query_refs().size(); ++i) {
        const auto& r = index.query_refs()[i];
        output << "  record_refs[" << i << "] V=" << graph.original_v_id(r.vertex)
               << " internal=" << r.vertex << " attribute_slot=" << r.attribute_slot
               << " core_slot=" << r.core_slot << '\n';
    }
    return finish_text(output, error);
}

bool save_index_text(const std::string& path, const TemporalBipartiteGraph& graph,
                     const VATIndex& index, std::string* error,
                     const std::string& data_directory, std::uint64_t build_peak_rss_bytes) {
    return save_vat_text(path, graph, index, error, data_directory, build_peak_rss_bytes, false);
}

bool save_index_intervals_text(const std::string& path, const TemporalBipartiteGraph& graph,
                              const VATIndex& index, std::string* error,
                              const std::string& data_directory, std::uint64_t build_peak_rss_bytes) {
    return save_vat_text(path, graph, index, error, data_directory, build_peak_rss_bytes, true);
}

bool save_block_index_text(const std::string& path,
                           const TemporalBipartiteGraph& graph,
                           const MATIndex& index,
                           bool core_first,
                           std::string* error,
                           const std::string& data_directory) {
    std::ofstream output(path);
    if (!output) {
        set_error(error, "open_failed");
        return false;
    }
    write_text_header(output, graph, index, data_directory);
    for (std::size_t i = 0; i < index.blocks().size(); ++i) {
        const IndexBlock& block = index.blocks()[i];
        output << "block=" << i + 1;
        if (core_first) {
            output << " core=(" << block.alpha << ',' << block.beta << ')'
                   << " attributes=";
            write_attributes(output, graph, block.attributes);
            output << " window=[" << block.start << ',' << block.end << ']';
        } else {
            output << " attributes=";
            write_attributes(output, graph, block.attributes);
            output << " window=[" << block.start << ',' << block.end << ']'
                   << " core=(" << block.alpha << ',' << block.beta << ')';
        }
        output << " vertices={";
        const std::size_t vertex_end =
            block.vertex_offset + block.vertex_count;
        for (std::size_t j = block.vertex_offset; j < vertex_end; ++j) {
            if (j != block.vertex_offset) output << ',';
            output << graph.original_v_id(index.vertex_pool()[j]);
        }
        output << "}\n";
    }
    output << "\n[QUERY_DIRECTORY] (U attribute set, beta) -> time group -> block range\n"
           << "Block ranges are inclusive and use the block numbers above.\n";
    for (const auto& entry : index.query_directory()) {
        output << "\nattributes=";
        write_attributes(output, graph, std::get<0>(entry.first));
        output << " beta=" << std::get<1>(entry.first)
               << " time_groups=" << entry.second.count << '\n';
        const auto group_end = entry.second.offset + entry.second.count;
        for (auto g = entry.second.offset; g < group_end; ++g) {
            const auto& group = index.query_time_groups()[g];
            output << "  window=[" << group.start << ',' << group.end << ']'
                   << " block_range=" << group.postings.offset + 1 << ".."
                   << group.postings.offset + group.postings.count << '\n';
        }
    }
    return finish_text(output, error);
}

bool save_index_text(const std::string& path,
                     const TemporalBipartiteGraph& graph,
                     const UATIndex& index,
                     std::string* error) {
    if (!index.ready()) {set_error(error,"uat_not_finalized");return false;}
    std::ofstream output(path);
    if (!output) {set_error(error,"open_failed");return false;}
    output << "TACC UAT index log\n";
    write_generated_at(output);
    output << "Scope: UV only (U-side attribute constraints, V-side stored vertices).\n"
           << "Actual structure: (alpha,beta) -> MAS -> ts -> MSE te -> shared vertex-pool range.\n"
           << "Construction algorithm: " << (index.stats().method==BuildMethod::Opt?"opt":"basic") << '\n'
           << "Index summary: cores=" << index.cores().size() << " records=" << index.stats().grouped_blocks
           << " logical_refs=" << index.stats().kept_records << " physical_refs=" << index.vertex_pool().size()
           << " build_s=" << index.stats().build_seconds
           << " index_mb=" << space_mb(formal_index_space(index).index_bytes()) << '\n'
           << "\n[LABELS]\n";
    for (AttributeID a=0;a<graph.num_attributes() && a<64;++a)
        output << "label_id=" << a << " label_name=" << graph.attribute_name(a)
               << " bit=" << (AttributeMask{1}<<a) << '\n';
    output << "\n[INDEX]\n";
    for (const auto& c : index.cores()) {
        output << "\ncore=(" << c.alpha << ',' << c.beta << ")\n";
        for (const auto& a : c.attributes) {
            output << "  MAS=" << a.mask << ' ';write_attributes(output,graph,a.mask);output << '\n';
            for (const auto& s : a.starts) {
                output << "    ts=" << s.start << '\n';
                for (const auto& e : s.ends)
                    output << "      te=" << e.end << " vertex_offset=" << e.vertex_offset
                           << " vertex_count=" << e.vertex_count << '\n';
            }
        }
    }
    output << "\n[SHARED_VERTEX_POOL]\n";
    for (std::size_t i=0;i<index.vertex_pool().size();++i)
        output << "vertex_pool[" << i << "] V=" << graph.original_v_id(index.vertex_pool()[i])
               << " internal=" << index.vertex_pool()[i] << '\n';
    return finish_text(output,error);
}

bool save_index_text(const std::string& path,
                     const TemporalBipartiteGraph& graph,
                     const MATIndex& index,
                     std::string* error,
                     const std::string& data_directory) {
    return save_block_index_text(path, graph, index, false, error, data_directory);
}

}
