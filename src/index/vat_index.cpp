#include "index/vat_index.h"
#include "index/pool_interner.h"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

namespace tacc {
namespace {
std::uint64_t word(const VATInterval& t) { return (std::uint64_t(t.start) << 32U) | t.end; }
bool equal_lists(const std::vector<VATInterval>& a, const std::vector<VATInterval>& b) {
    return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin(),
        [](const auto& x, const auto& y) { return word(x) == word(y); });
}
void require(bool ok, const char* why) { if (!ok) throw std::runtime_error(why); }
struct DirectoryKey {
    AttributeMask mask; Count alpha, beta;
    bool operator==(const DirectoryKey& b) const {
        return mask == b.mask && alpha == b.alpha && beta == b.beta;
    }
    bool operator<(const DirectoryKey& b) const {
        return std::tie(mask, alpha, beta) < std::tie(b.mask, b.alpha, b.beta);
    }
};
struct DirectoryKeyHash {
    static std::uint64_t mix(std::uint64_t n) {
        n ^= n >> 30; n *= 0xbf58476d1ce4e5b9ULL;
        n ^= n >> 27; n *= 0x94d049bb133111ebULL; return n ^ (n >> 31);
    }
    std::size_t operator()(const DirectoryKey& k) const {
        return mix(k.mask) ^ mix((std::uint64_t(k.alpha) << 32) | k.beta);
    }
};
}
void VATIndex::initialize_vertices(std::size_t v) {
    if (finalized_ || build_tree_) throw std::logic_error("VAT initialize after construction");
    uv_vertices_.resize(v);
}
bool VATIndex::add(IndexMembership r) {
    if (finalized_) throw std::logic_error("VAT add after finalize");
    if (!r.attributes || !r.alpha || !r.beta || r.start > r.end)
        throw std::invalid_argument("invalid VAT membership");
    if (!build_tree_) build_tree_ = std::make_unique<BuildTree>();
    auto& tree = *build_tree_;
    if (tree.size() <= r.vertex) tree.resize(std::size_t(r.vertex) + 1);
    auto& attrs = tree[r.vertex];
    const auto find_core = [&](auto& cores) -> BuildCore* {
        const auto it = std::find_if(cores.begin(), cores.end(), [&](const auto& c) {
            return c.alpha == r.alpha && c.beta == r.beta;
        });
        return it == cores.end() ? nullptr : &*it;
    };
    ++stats_.generated_records;
    for (auto& a : attrs) if (attribute_subset(a.mask, r.attributes)) {
        if (auto* c = find_core(a.cores)) for (const auto& t : c->intervals)
            if (t.start >= r.start && t.end <= r.end) {
                ++stats_.attribute_time_removed; return false;
            }
    }
    for (auto& a : attrs) if (attribute_subset(r.attributes, a.mask)) {
        if (auto* c = find_core(a.cores)) {
            const auto old = c->intervals.size();
            c->intervals.erase(std::remove_if(c->intervals.begin(), c->intervals.end(),
                [&](const auto& t) { return r.start >= t.start && r.end <= t.end; }), c->intervals.end());
            stats_.attribute_time_removed += old - c->intervals.size();
        }
    }
    auto a = std::find_if(attrs.begin(), attrs.end(), [&](const auto& x) { return x.mask == r.attributes; });
    if (a == attrs.end()) { attrs.push_back({r.attributes, {}}); a = attrs.end() - 1; }
    auto* c = find_core(a->cores);
    if (!c) { a->cores.push_back({r.alpha, r.beta, {}}); c = &a->cores.back(); }
    c->intervals.push_back({r.start, r.end});
    return true;
}
void VATIndex::finalize() {
    if (finalized_) return;
    if (build_tree_) {
        for (auto& attrs : *build_tree_) {
            for (auto& a : attrs) {
                a.cores.erase(std::remove_if(a.cores.begin(), a.cores.end(),
                    [](const auto& c) { return c.intervals.empty(); }), a.cores.end());
                std::sort(a.cores.begin(), a.cores.end(), [](const auto& x, const auto& y) {
                    return std::tie(x.alpha, x.beta) < std::tie(y.alpha, y.beta);
                });
                for (auto& c : a.cores) {
                    std::sort(c.intervals.begin(), c.intervals.end(), [](const auto& x, const auto& y) {
                        return std::tie(x.start, x.end) < std::tie(y.start, y.end);
                    });
                    c.intervals.erase(std::unique(c.intervals.begin(), c.intervals.end(),
                        [](const auto& x, const auto& y) { return word(x) == word(y); }), c.intervals.end());
                }
            }
            attrs.erase(std::remove_if(attrs.begin(), attrs.end(),
                [](const auto& a) { return a.cores.empty(); }), attrs.end());
            std::sort(attrs.begin(), attrs.end(), [](const auto& a, const auto& b) { return a.mask < b.mask; });
        }
        const auto limit = PoolInterner::configured_limit();
        std::size_t planned = 0;
        {
            std::unordered_multimap<std::uint64_t, const std::vector<VATInterval>*> cache;
            const auto visit = [&](const auto& list) {
                std::uint64_t h = 1469598103934665603ULL;
                for (const auto& t : list) { h ^= word(t); h *= 1099511628211ULL; }
                const auto range = cache.equal_range(h);
                for (auto p = range.first; p != range.second; ++p)
                    if (equal_lists(*p->second, list)) return;
                if (list.size() > std::numeric_limits<std::size_t>::max() - planned)
                    throw std::overflow_error("VAT pool size overflow");
                planned += list.size();
                if (cache.size() < limit) cache.emplace(h, &list);
            };
            for (const auto& attrs : *build_tree_)
                for (const auto& a : attrs) for (const auto& c : a.cores) visit(c.intervals);
        }
        intervals_.reserve(planned);
        {
            PoolInterner interner(limit);
            {
                auto& source = *build_tree_;
                auto& tree = uv_vertices_;
                if (tree.size() < source.size()) tree.resize(source.size());
                for (std::size_t v = 0; v < source.size(); ++v) {
                    auto& target = tree[v]; target.reserve(source[v].size());
                    for (auto& a : source[v]) {
                        target.push_back({a.mask, {}});
                        auto& cores = target.back().cores; cores.reserve(a.cores.size());
                        for (auto& c : a.cores) {
                            const auto range = interner.intern(intervals_, c.intervals, word);
                            cores.push_back({c.alpha, c.beta, range.offset, range.count});
                            std::vector<VATInterval>().swap(c.intervals);
                        }
                        std::vector<BuildCore>().swap(a.cores);
                    }
                    std::vector<BuildAttribute>().swap(source[v]);
                }
            }
        }
        build_tree_.reset();
        require(intervals_.size() == planned, "VAT planning/packing disagreement");
    }
    build_query_directory();
    finalized_ = true;
    recount();
    validate();
}
void VATIndex::build_query_directory() {
    std::unordered_map<DirectoryKey, std::size_t, DirectoryKeyHash> positions;
    std::size_t records = 0;
    for (const auto& attrs : uv_vertices_) for (const auto& a : attrs) for (const auto& c : a.cores) {
        require(records < std::numeric_limits<std::size_t>::max(), "VAT directory size overflow");
        ++positions[{a.mask, c.alpha, c.beta}]; ++records;
    }
    std::vector<DirectoryKey> keys;
    keys.reserve(positions.size());
    for (const auto& entry : positions) keys.push_back(entry.first);
    std::sort(keys.begin(), keys.end());
    std::size_t masks = 0;
    for (std::size_t i = 0; i < keys.size(); ++i)
        masks += i == 0 || keys[i].mask != keys[i - 1].mask;
    query_masks_.resize(masks); query_cores_.resize(keys.size()); query_refs_.resize(records);
    std::size_t cursor = 0, mi = 0;
    for (std::size_t i = 0; i < keys.size(); ++i) {
        const auto k = keys[i]; auto& count = positions.at(k);
        if (!i || k.mask != keys[i - 1].mask) {
            if (i) ++mi;
            query_masks_[mi] = {k.mask, i, 0};
        }
        ++query_masks_[mi].count;
        query_cores_[i] = {k.alpha, k.beta, cursor, count};
        const auto n = count; count = cursor; cursor += n;
    }
    for (std::size_t v = 0; v < uv_vertices_.size(); ++v) {
        const auto& attrs = uv_vertices_[v];
        require(v <= std::numeric_limits<VertexID>::max() &&
                attrs.size() <= std::numeric_limits<std::uint32_t>::max(), "VAT directory vertex/attribute overflow");
        for (std::size_t a = 0; a < attrs.size(); ++a) {
            require(attrs[a].cores.size() <= std::numeric_limits<std::uint32_t>::max(), "VAT directory core slot overflow");
            for (std::size_t c = 0; c < attrs[a].cores.size(); ++c) {
                const auto& leaf = attrs[a].cores[c];
                auto& pos = positions.at({attrs[a].mask, leaf.alpha, leaf.beta});
                require(pos < query_refs_.size(), "VAT directory reference overflow");
                query_refs_[pos++] = {VertexID(v), std::uint32_t(a), std::uint32_t(c)};
            }
        }
    }
    for (std::size_t i = 0; i < keys.size(); ++i)
        require(positions.at(keys[i]) == query_cores_[i].offset + query_cores_[i].count,
                "incomplete VAT directory");
}
std::size_t VATIndex::query_directory_bytes() const {
    return query_masks_.size() * sizeof(VATMaskRange) +
        query_cores_.size() * sizeof(VATCoreRange) + query_refs_.size() * sizeof(VATRecordRef);
}
std::size_t VATIndex::tree_bytes() const {
    std::size_t bytes = 0;
    {
        const auto& tree = uv_vertices_;
        bytes += tree.size() * sizeof(VertexTree::value_type);
        for (const auto& attrs : tree) {
            bytes += attrs.size() * sizeof(VATAttributeEntry);
            for (const auto& a : attrs) bytes += a.cores.size() * sizeof(VATCoreEntry);
        }
    }
    return bytes;
}
void VATIndex::recount() {
    stats_.grouped_blocks = stats_.kept_records = 0;
    {
        for (const auto& attrs : uv_vertices_)
            for (const auto& a : attrs) for (const auto& c : a.cores) {
                ++stats_.grouped_blocks;
                stats_.kept_records += c.interval_count;
            }
    }
    stats_.vertex_references = intervals_.size();
    stats_.approx_bytes = tree_bytes() + intervals_.size() * sizeof(VATInterval) + query_directory_bytes();
}
void VATIndex::validate() const {
    require(finalized_ && !build_tree_, "VAT is not finalized");
    std::size_t leaves = 0, logical = 0;
    {
        for (const auto& attrs : uv_vertices_) {
            AttributeMask previous = 0;
            for (const auto& a : attrs) {
                require(a.mask > previous && !a.cores.empty(), "invalid VAT MAS order"); previous = a.mask;
                std::pair<Count, Count> last{};
                for (const auto& c : a.cores) {
                    require(leaves < std::numeric_limits<std::size_t>::max() &&
                            c.interval_count <= std::numeric_limits<std::size_t>::max() - logical, "VAT count overflow");
                    ++leaves; logical += c.interval_count;
                    const auto key = std::make_pair(c.alpha, c.beta);
                    require(c.alpha && c.beta && key > last && c.interval_count, "invalid VAT core order"); last = key;
                    require(c.interval_offset <= intervals_.size() &&
                        c.interval_count <= intervals_.size() - c.interval_offset, "invalid VAT pool range");
                    for (std::size_t i = 0; i < c.interval_count; ++i) {
                        const auto& t = intervals_[c.interval_offset + i];
                        require(t.start <= t.end, "invalid VAT interval");
                        if (i) {
                            const auto& p = intervals_[c.interval_offset + i - 1];
                            require(p.start < t.start && p.end < t.end, "invalid VAT MAI staircase");
                        }
                    }
                }
            }
        }
    }
    std::size_t ci = 0, ri = 0; AttributeMask previous = 0;
    for (const auto& m : query_masks_) {
        require(m.mask > previous && m.offset == ci && m.count &&
                ci <= query_cores_.size() && m.count <= query_cores_.size() - ci, "invalid VAT directory mask range");
        previous = m.mask; std::pair<Count, Count> last{};
        for (std::size_t i = m.offset; i < m.offset + m.count; ++i) {
            const auto& c = query_cores_[i]; const auto key = std::make_pair(c.alpha, c.beta);
            require(c.alpha && c.beta && key > last && c.offset == ri && c.count &&
                    ri <= query_refs_.size() && c.count <= query_refs_.size() - ri, "invalid VAT directory core range");
            last = key; VertexID last_v = 0; bool first = true;
            for (std::size_t j = c.offset; j < c.offset + c.count; ++j) {
                const auto& r = query_refs_[j];
                require(r.vertex < uv_vertices_.size() && (first || r.vertex > last_v), "invalid/duplicate VAT directory vertex");
                first = false; last_v = r.vertex;
                const auto& attrs = uv_vertices_[r.vertex];
                require(r.attribute_slot < attrs.size(), "invalid VAT directory attribute slot");
                const auto& a = attrs[r.attribute_slot];
                require(a.mask == m.mask && r.core_slot < a.cores.size(), "invalid VAT directory core slot");
                const auto& leaf = a.cores[r.core_slot];
                require(leaf.alpha == c.alpha && leaf.beta == c.beta, "VAT directory references a different leaf");
            }
            ri += c.count;
        }
        ci += m.count;
    }
    require(ci == query_cores_.size() && ri == query_refs_.size() && ri == leaves,
            "VAT directory does not cover every leaf exactly once");
    require(leaves == stats_.grouped_blocks && logical == stats_.kept_records &&
            intervals_.size() == stats_.vertex_references &&
            stats_.approx_bytes == tree_bytes() + intervals_.size() * sizeof(VATInterval) + query_directory_bytes(),
            "VAT index statistics disagree with storage");
}
void VATIndex::clear() {
    VertexTree().swap(uv_vertices_); build_tree_.reset();
    std::vector<VATInterval>().swap(intervals_);
    std::vector<VATMaskRange>().swap(query_masks_);
    std::vector<VATCoreRange>().swap(query_cores_);
    std::vector<VATRecordRef>().swap(query_refs_);
    stats_ = {}; finalized_ = false;
}
}
