#include "index/uat_index.h"
#include "index/pool_interner.h"

#include <algorithm>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace tacc {
namespace {
using CompactBlockKey =
    std::tuple<Count, Count, std::uint32_t, Timestamp, Timestamp>;
}

UATIndex::BuildAttributeID UATIndex::build_attribute_id(
    AttributeMask attributes) {
    if (has_last_build_attribute_ &&
        attributes == last_build_attribute_mask_) {
        return last_build_attribute_id_;
    }
    const auto found = std::find(build_attribute_masks_.begin(),
                                 build_attribute_masks_.end(), attributes);
    std::size_t position = static_cast<std::size_t>(
        std::distance(build_attribute_masks_.begin(), found));
    if (found == build_attribute_masks_.end()) {
        if (build_attribute_masks_.size() >=
            std::numeric_limits<BuildAttributeID>::max()) {
            throw std::overflow_error(
                "UAT build attribute dictionary exceeds uint32_t");
        }
        position = build_attribute_masks_.size();
        build_attribute_masks_.push_back(attributes);
    }
    last_build_attribute_mask_ = attributes;
    last_build_attribute_id_ = static_cast<BuildAttributeID>(position);
    has_last_build_attribute_ = true;
    return last_build_attribute_id_;
}

void UATIndex::clear_build_buckets() {
    auto& buckets = build_buckets_;

    std::vector<std::uint32_t>().swap(buckets.vertex_positions);
    std::vector<VertexBucket>().swap(buckets.active_vertices);

    std::vector<AttributeMask>().swap(build_attribute_masks_);
    last_build_attribute_mask_ = 0;
    last_build_attribute_id_ = 0;
    has_last_build_attribute_ = false;
    finalized_ = false;
    logical_references_ = 0;
    build_interner_.reset();
}

bool UATIndex::add(IndexMembership record) {
    if (finalized_)
        throw std::logic_error("UAT add after constraint side was finalized");
    if (!record.attributes || !record.alpha || !record.beta || record.start > record.end)
        throw std::invalid_argument("invalid UAT membership");
    ++stats_.generated_records;
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
            throw std::overflow_error("UAT active vertex locator exceeds uint32_t");
        }
        position = static_cast<std::uint32_t>(
            buckets.active_vertices.size());
        buckets.active_vertices.push_back({record.vertex, {}});
    }
    VertexBucket& vertex = buckets.active_vertices[position];
    auto core_it = std::find_if(
        vertex.cores.begin(), vertex.cores.end(),
        [&](const CoreBucket& core) {
            return core.alpha == record.alpha && core.beta == record.beta;
        });
    if (core_it == vertex.cores.end()) {
        vertex.cores.push_back({record.alpha, record.beta, record.vertex, {}});
        core_it = std::prev(vertex.cores.end());
    }
    auto& bucket = core_it->memberships;
    const BuildAttributeID attribute_id =
        build_attribute_id(record.attributes);
    const BuildMembership candidate{
        attribute_id, record.start, record.end};
    for (const BuildMembership& kept : bucket) {
        if (attribute_subset(build_attribute_masks_[kept.attribute_id],
                             record.attributes) &&
            kept.start == record.start && kept.end <= record.end) {
            ++stats_.attribute_time_removed;
            return false;
        }
    }
    bucket.erase(std::remove_if(bucket.begin(), bucket.end(),
        [&](const BuildMembership& kept) {
            if (!attribute_subset(
                    record.attributes,
                    build_attribute_masks_[kept.attribute_id]) ||
                record.start != kept.start || record.end > kept.end) {
                return false;
            }
            ++stats_.attribute_time_removed;
            return true;
        }), bucket.end());
    bucket.push_back(candidate);
    return true;
}

void UATIndex::finalize_buckets() {
    using FlatRef = std::uint64_t;
    static_assert(sizeof(BuildMembership) == 12,
                  "UAT BuildMembership must remain compact");
    if (!build_interner_) build_interner_.emplace();
    std::vector<VertexID> vertices;

    std::vector<BuildAttributeID> sorted_attribute_ids(
        build_attribute_masks_.size());
    std::iota(sorted_attribute_ids.begin(), sorted_attribute_ids.end(), 0);
    std::sort(sorted_attribute_ids.begin(), sorted_attribute_ids.end(),
              [&](BuildAttributeID lhs, BuildAttributeID rhs) {
                  return build_attribute_masks_[lhs] <
                         build_attribute_masks_[rhs];
              });
    std::vector<BuildAttributeID> attribute_ranks(
        build_attribute_masks_.size());
    for (BuildAttributeID rank = 0;
         rank < sorted_attribute_ids.size(); ++rank) {
        attribute_ranks[sorted_attribute_ids[rank]] = rank;
    }

    auto& buckets = build_buckets_;
    std::size_t core_count = 0;
    for (const auto& vertex : buckets.active_vertices)
        core_count += vertex.cores.size();
    std::vector<CoreBucket*> core_table;
    core_table.reserve(core_count);
    for (auto& vertex : buckets.active_vertices)
        for (auto& core : vertex.cores) core_table.push_back(&core);
    stats_.uat_peak_core_directory_bytes = std::max(
        stats_.uat_peak_core_directory_bytes, core_table.capacity() * sizeof(CoreBucket*));
    std::sort(core_table.begin(), core_table.end(), [](const CoreBucket* a, const CoreBucket* b) {
        return std::tie(a->alpha, a->beta, a->vertex) < std::tie(b->alpha, b->beta, b->vertex);
    });
    std::vector<FlatRef> flat;
    for (std::size_t core_begin = 0; core_begin < core_table.size();) {
        auto core_end = core_begin + 1;
        while (core_end < core_table.size() &&
               core_table[core_end]->alpha == core_table[core_begin]->alpha &&
               core_table[core_end]->beta == core_table[core_begin]->beta) ++core_end;
        if (core_end - core_begin > std::numeric_limits<std::uint32_t>::max())
            throw std::overflow_error("UAT core GROUP locator exceeds uint32_t");
        std::size_t group_memberships = 0;
        for (auto c = core_begin; c < core_end; ++c) {
            const auto count = core_table[c]->memberships.size();
            if (count > std::numeric_limits<std::uint32_t>::max() ||
                count > flat.max_size() - group_memberships)
                throw std::overflow_error("UAT group membership locator exceeds capacity");
            group_memberships += count;
        }
        flat.clear();
        flat.reserve(group_memberships);
        for (auto c = core_begin; c < core_end; ++c) {
            for (std::size_t m = 0; m < core_table[c]->memberships.size(); ++m)
                flat.push_back((static_cast<std::uint64_t>(c - core_begin) << 32U) | m);
        }
        stats_.uat_peak_flat_bytes = std::max(stats_.uat_peak_flat_bytes,
                                              flat.capacity() * sizeof(FlatRef));
        ++stats_.uat_core_groups;

        const auto core_of = [&](FlatRef ref) -> const CoreBucket& {
            return *core_table[core_begin + static_cast<std::uint32_t>(ref >> 32U)];
        };
        const auto membership_of = [&](FlatRef ref,
                                       const CoreBucket& core)
            -> const BuildMembership& {
            return core.memberships[static_cast<std::uint32_t>(ref)];
        };
        const auto key_of = [&](FlatRef ref) {
            const CoreBucket& core = core_of(ref);
            const BuildMembership& membership = membership_of(ref, core);
            return CompactBlockKey{
                core.alpha, core.beta,
                attribute_ranks[membership.attribute_id],
                membership.start, membership.end};
        };
        std::sort(flat.begin(), flat.end(),
                  [&](FlatRef lhs, FlatRef rhs) {
                      const CoreBucket& lhs_core = core_of(lhs);
                      const CoreBucket& rhs_core = core_of(rhs);
                      const BuildMembership& lhs_membership =
                          membership_of(lhs, lhs_core);
                      const BuildMembership& rhs_membership =
                          membership_of(rhs, rhs_core);
                      return std::tie(
                                 lhs_core.alpha, lhs_core.beta,
                                 attribute_ranks[
                                     lhs_membership.attribute_id],
                                 lhs_membership.start, lhs_membership.end,
                                 lhs_core.vertex) <
                             std::tie(
                                 rhs_core.alpha, rhs_core.beta,
                                 attribute_ranks[
                                     rhs_membership.attribute_id],
                                 rhs_membership.start, rhs_membership.end,
                                 rhs_core.vertex);
                  });

        for (std::size_t begin = 0; begin < flat.size();) {
            const CompactBlockKey key = key_of(flat[begin]);
            std::size_t group_end = begin + 1;
            while (group_end < flat.size() &&
                   key_of(flat[group_end]) == key) {
                ++group_end;
            }
            const auto [alpha, beta, attribute_rank, start, end] = key;
            const BuildAttributeID attribute_id =
                sorted_attribute_ids[attribute_rank];
            const AttributeMask attrs =
                build_attribute_masks_[attribute_id];
            vertices.clear();
            VertexID previous = 0;
            bool has_previous = false;
            for (std::size_t i = begin; i < group_end; ++i) {
                const VertexID vertex = core_of(flat[i]).vertex;
                if (has_previous && vertex == previous) continue;
                previous = vertex;
                has_previous = true;
                vertices.push_back(vertex);
            }
            const auto range = build_interner_->intern(vertex_pool_, vertices,
                [](VertexID vertex) { return std::uint64_t{vertex}; });
            logical_references_ += range.count;
            if (cores_.empty() || cores_.back().alpha != alpha || cores_.back().beta != beta)
                cores_.push_back({alpha,beta,{}});
            auto& core = cores_.back();
            if (core.attributes.empty() || core.attributes.back().mask != attrs)
                core.attributes.push_back({attrs,{}});
            auto& attribute = core.attributes.back();
            if (attribute.starts.empty() || attribute.starts.back().start != start)
                attribute.starts.push_back({start,{}});
            attribute.starts.back().ends.push_back({end,range.offset,range.count});
            begin = group_end;
        }
        flat.clear();
        for (auto c = core_begin; c < core_end; ++c)
            std::vector<BuildMembership>().swap(core_table[c]->memberships);
        core_begin = core_end;
    }
    std::vector<std::uint32_t>().swap(buckets.vertex_positions);
    std::vector<VertexBucket>().swap(buckets.active_vertices);
    ++stats_.uat_completed_directions;
}

void UATIndex::finalize() {
    if (finalized_) return;
    finalize_buckets();
    std::vector<AttributeMask>().swap(build_attribute_masks_);
    has_last_build_attribute_ = false;
    build_interner_.reset();
    for (auto& c : cores_) {
        for (auto& a : c.attributes) {
            for (auto& s : a.starts) s.ends.shrink_to_fit();
            a.starts.shrink_to_fit();
        }
        c.attributes.shrink_to_fit();
    }
    cores_.shrink_to_fit(); vertex_pool_.shrink_to_fit();
    finalized_ = true;
    recount();
    validate();
}

std::size_t UATIndex::tree_bytes() const {
    std::size_t bytes = cores_.size()*sizeof(UATCoreEntry);
    for (const auto& c : cores_) {
        bytes += c.attributes.size()*sizeof(UATAttributeEntry);
        for (const auto& a : c.attributes) {
            bytes += a.starts.size()*sizeof(UATStartEntry);
            for (const auto& s : a.starts) bytes += s.ends.size()*sizeof(UATEndEntry);
        }
    }
    return bytes;
}

void UATIndex::recount() {
    stats_.grouped_blocks = stats_.kept_records = 0;
    for_each_block([&](const IndexBlock& b) {
        ++stats_.grouped_blocks;
        if (b.vertex_count > std::numeric_limits<std::size_t>::max()-stats_.kept_records)
            throw std::overflow_error("UAT logical references overflow");
        stats_.kept_records += b.vertex_count;
    });
    logical_references_ = stats_.kept_records;
    stats_.vertex_references = vertex_pool_.size();
    stats_.approx_bytes = tree_bytes()+vertex_pool_.size()*sizeof(VertexID);
}

void UATIndex::validate() const {
    if (!finalized_) throw std::logic_error("UAT not finalized");
    std::pair<Count,Count> last_core{};
    std::size_t covered=0;
    for (const auto& c : cores_) {
        const auto key=std::make_pair(c.alpha,c.beta);
        if (!c.alpha || !c.beta || key<=last_core || c.attributes.empty())
            throw std::runtime_error("invalid UAT core ordering");
        last_core=key; AttributeMask last_mask=0;
        for (const auto& a : c.attributes) {
            if (!a.mask || a.mask<=last_mask || a.starts.empty())
                throw std::runtime_error("invalid UAT attribute ordering");
            last_mask=a.mask;
            for (std::size_t i=0;i<a.starts.size();++i) {
                const auto& s=a.starts[i];
                if ((i && a.starts[i-1].start>=s.start) || s.ends.empty())
                    throw std::runtime_error("invalid UAT start ordering");
                for (std::size_t j=0;j<s.ends.size();++j) {
                    const auto& e=s.ends[j];
                    if (e.end<s.start || (j && s.ends[j-1].end>=e.end) || !e.vertex_count ||
                        e.vertex_offset>vertex_pool_.size() || e.vertex_count>vertex_pool_.size()-e.vertex_offset)
                        throw std::runtime_error("invalid UAT end/pool range");
                    covered=std::max(covered,e.vertex_offset+e.vertex_count);
                    for (std::size_t k=1;k<e.vertex_count;++k)
                        if (vertex_pool_[e.vertex_offset+k-1]>=vertex_pool_[e.vertex_offset+k])
                            throw std::runtime_error("invalid UAT vertex ordering");
                }
            }
        }
    }
    if (covered!=vertex_pool_.size()) throw std::runtime_error("unreferenced UAT pool suffix");
}

void UATIndex::clear() {
    clear_build_buckets();
    cores_.clear();
    vertex_pool_.clear();
    stats_ = {};
}
}
