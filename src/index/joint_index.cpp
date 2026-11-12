#include "index/joint_index.h"
#include "index/index_io.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace tacc::joint {
namespace {
bool less(AttributePair a, AttributePair b) { return a.u < b.u || (a.u == b.u && a.v < b.v); }
bool same(AttributePair a, AttributePair b) { return a.u == b.u && a.v == b.v; }
void check_core(Count a, Count b) {
    if (a == 1 && b == 1) throw std::invalid_argument("joint index does not store core (1,1)");
}
QueryResult canonical(QueryResult result) {
    for (auto* side : {&result.u_vertices, &result.v_vertices}) {
        std::sort(side->begin(), side->end());
        side->erase(std::unique(side->begin(), side->end()), side->end());
    }
    return result;
}
template<class T> std::size_t bytes(const std::vector<T>& values) { return sizeof(T) * values.size(); }

void integer(std::ostream& out, std::uint64_t value) {
    for (unsigned i = 0; i < 8; ++i) { out.put(static_cast<char>(value & 255)); value >>= 8; }
}
void pair(std::ostream& out, AttributePair value) { integer(out, value.u); integer(out, value.v); }
void ids(std::ostream& out, const std::vector<VertexID>& values) {
    integer(out, values.size());
    for (auto id : values) integer(out, id);
}
void close(std::ofstream& out) {
    out.close();
    if (!out) throw std::runtime_error("cannot finish joint index output");
}
std::ofstream binary(const std::filesystem::path& path, const char* magic, const TemporalBipartiteGraph& graph) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot create joint index binary: " + path.string());
    out.write(magic, 8);
    integer(out, graph_fingerprint(graph));
    return out;
}

class BinaryReader {
public:
    BinaryReader(const std::filesystem::path& path, const char* magic, const TemporalBipartiteGraph& graph,
                 bool allow_legacy_vat = false)
        : input_(path, std::ios::binary) {
        check(bool(input_), "cannot open joint index: " + path.string());
        remaining_ = std::filesystem::file_size(path);
        check(remaining_ >= 16, "truncated joint index header");
        char header[8];
        input_.read(header, 8);
        remaining_ -= 8;
        const std::string format(header, 8);
        legacy_vat_ = allow_legacy_vat && format == "TACJVAT1";
        check(bool(input_) && (format == std::string(magic, 8) || legacy_vat_), "joint index type/version mismatch");
        check(integer() == graph_fingerprint(graph), "joint index graph fingerprint mismatch");
    }
    static void check(bool valid, const std::string& message) {
        if (!valid) throw std::runtime_error(message);
    }
    static void check(bool valid, const char* message) {
        if (!valid) throw std::runtime_error(message);
    }
    std::uint64_t integer() {
        check(remaining_ >= 8, "truncated joint index integer");
        if (buffer_pos_ == buffer_size_) {
            buffer_size_ = static_cast<std::size_t>(std::min<std::uintmax_t>(remaining_, buffer_.size()));
            input_.read(reinterpret_cast<char*>(buffer_.data()), static_cast<std::streamsize>(buffer_size_));
            check(bool(input_), "cannot read joint index integer");
            buffer_pos_ = 0;
        }
        const auto* raw = buffer_.data() + buffer_pos_;
        buffer_pos_ += 8;
        remaining_ -= 8;
        std::uint64_t value = 0;
        for (unsigned i = 0; i < 8; ++i) value |= std::uint64_t(raw[i]) << (8 * i);
        return value;
    }
    template<class T> T value() {
        const auto number = integer();
        check(number <= std::numeric_limits<T>::max(), "joint index integer overflow");
        return static_cast<T>(number);
    }
    std::size_t count(std::size_t minimum_bytes) {
        const auto number = integer();
        check(number <= std::numeric_limits<std::size_t>::max() && number <= remaining_ / minimum_bytes,
              "invalid joint index container length");
        return static_cast<std::size_t>(number);
    }
    void finish() { check(remaining_ == 0, "unexpected trailing joint index bytes"); }
    bool legacy_vat() const { return legacy_vat_; }
private:
    std::ifstream input_;
    std::array<unsigned char, 64 * 1024> buffer_;
    std::size_t buffer_pos_{}, buffer_size_{};
    std::uintmax_t remaining_{};
    bool legacy_vat_{};
};

struct ReadContext {
    AttributePair all;
    std::size_t max_alpha{}, max_beta{};
    std::size_t num_u{}, num_v{};
    const std::vector<Timestamp>& times;
    explicit ReadContext(const TemporalBipartiteGraph& graph)
        : num_u(graph.num_u_vertices()), num_v(graph.num_v_vertices()), times(graph.time_points()) {
        for (VertexID u = 0; u < graph.num_u_vertices(); ++u) {
            BinaryReader::check(graph.u_attribute(u) < 64, "joint index attribute overflow");
            all.u |= AttributeMask{1} << graph.u_attribute(u);
            max_alpha = std::max(max_alpha, graph.static_u_degree(u));
        }
        for (VertexID v = 0; v < graph.num_v_vertices(); ++v) {
            BinaryReader::check(graph.v_attribute(v) < 64, "joint index attribute overflow");
            all.v |= AttributeMask{1} << graph.v_attribute(v);
            max_beta = std::max(max_beta, graph.static_v_degree(v));
        }
    }
    AttributePair attributes(BinaryReader& in) const {
        const AttributePair result{in.integer(), in.integer()};
        BinaryReader::check(result.u && result.v && !(result.u & ~all.u) && !(result.v & ~all.v),
                            "invalid joint index attribute pair");
        return result;
    }
};

std::vector<VertexID> read_ids(BinaryReader& in, std::size_t max_count) {
    const auto count = in.count(8);
    BinaryReader::check(count <= max_count, "invalid joint index vertex count");
    std::vector<VertexID> values;
    values.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        values.push_back(in.value<VertexID>());
    }
    return values;
}

std::string names(const TemporalBipartiteGraph& graph, AttributeMask mask) {
    std::ostringstream out;
    out << '{';
    bool first = true;
    for (AttributeID id = 0; id < 64; ++id) if (mask & (AttributeMask{1} << id)) {
        if (!first) out << ',';
        first = false;
        out << graph.attribute_name(id);
    }
    out << '}';
    return out.str();
}
void text_header(std::ostream& out, const char* method, const TemporalBipartiteGraph& graph,
                 const std::filesystem::path& data, const std::string& generated_at, std::size_t size,
                 const std::string& build_method) {
    out << "TACC exact joint " << method << " index\nGenerated at: " << generated_at
        << "\nDataset: " << data.filename().string() << "\nData directory: " << data.string()
        << "\nBuild method: " << build_method << "\nScope: exact_joint_attributes\nAttribute compression: disabled"
        << "\nCore (1,1): excluded\nGraph fingerprint: " << graph_fingerprint(graph)
        << "\nVertices: U=" << graph.num_u_vertices() << " V=" << graph.num_v_vertices()
        << "\nTemporal events: " << graph.num_edges() << "\nTime points: " << graph.time_points().size()
        << "\nindex_mb=" << std::fixed << std::setprecision(9) << size / 1048576.0 << '\n';
}

QueryResult directory_query(const VAT& index, AttributePair key, Count a, Count b, TimeWindow window) {
    QueryResult result;
    const auto& directory = index.directory;
    const auto pair_it = std::lower_bound(directory.pairs.begin(), directory.pairs.end(), key, less);
    if (pair_it == directory.pairs.end() || !same(*pair_it, key)) return result;
    const auto& grid = directory.groups[static_cast<std::size_t>(pair_it - directory.pairs.begin())];
    if (a >= grid.size() || b >= grid[a].size()) return result;
    for (const auto ref : grid[a][b]) {
        const bool upper = ref.vertex_side & vat_upper_flag;
        const auto& vertex = (upper ? index.u_vertices : index.v_vertices)[ref.vertex_side & ~vat_upper_flag];
        const auto& attr = vertex.attributes[ref.attribute_slot];
        const auto& row = attr.alpha_rows[a - attr.alpha_begin];
        const auto& intervals = row.beta_lists[b - row.beta_begin];
        if (intervals.empty() || window.end < intervals.front().end || window.start > intervals.back().start) continue;
        const auto hit = std::lower_bound(intervals.begin(), intervals.end(), window.start,
            [](Interval interval, Timestamp start) { return interval.start < start; });
        if (hit != intervals.end() && hit->end <= window.end)
            (upper ? result.u_vertices : result.v_vertices).push_back(vertex.original_id);
    }
    return result;
}
}

void build_vat_directory(VAT& index) {
    VATDirectory directory;
    std::set<std::pair<AttributeMask, AttributeMask>> keys;
    for (const auto* vertices : {&index.u_vertices, &index.v_vertices})
        for (const auto& vertex : *vertices) for (const auto& attr : vertex.attributes)
            keys.emplace(attr.attributes.u, attr.attributes.v);
    directory.pairs.reserve(keys.size());
    for (const auto& key : keys) directory.pairs.push_back({key.first, key.second});
    directory.groups.resize(directory.pairs.size());
    for (bool upper : {true, false}) {
        const auto& vertices = upper ? index.u_vertices : index.v_vertices;
        if (vertices.size() > vat_upper_flag) throw std::runtime_error("VAT directory vertex slot overflow");
        for (std::size_t vi = 0; vi < vertices.size(); ++vi) {
            const auto& vertex = vertices[vi];
            if (vertex.attributes.size() > std::numeric_limits<std::uint32_t>::max())
                throw std::runtime_error("VAT directory attribute slot overflow");
            for (std::size_t slot = 0; slot < vertex.attributes.size(); ++slot) {
                const auto& attr = vertex.attributes[slot];
                const auto key = std::lower_bound(directory.pairs.begin(), directory.pairs.end(), attr.attributes, less);
                auto& grid = directory.groups[static_cast<std::size_t>(key - directory.pairs.begin())];
                for (std::size_t ai = 0; ai < attr.alpha_rows.size(); ++ai) {
                    const auto a = std::size_t(attr.alpha_begin) + ai;
                    const auto& row = attr.alpha_rows[ai];
                    for (std::size_t bi = 0; bi < row.beta_lists.size(); ++bi) {
                        if (row.beta_lists[bi].empty()) continue;
                        const auto b = std::size_t(row.beta_begin) + bi;
                        if (grid.size() <= a) grid.resize(a + 1);
                        if (grid[a].size() <= b) grid[a].resize(b + 1);
                        grid[a][b].push_back({static_cast<std::uint32_t>(vi) | (upper ? vat_upper_flag : 0),
                                             static_cast<std::uint32_t>(slot)});
                    }
                }
            }
        }
    }
    index.directory = std::move(directory);
    index.directory_ready = true;
}

std::size_t vat_directory_bytes(const VAT& index) {
    if (!index.directory_ready) return 0;
    const auto& directory = index.directory;
    std::size_t total = sizeof(directory) + bytes(directory.pairs) + bytes(directory.groups);
    for (const auto& grid : directory.groups) {
        total += bytes(grid);
        for (const auto& row : grid) {
            total += bytes(row);
            for (const auto& refs : row) total += bytes(refs);
        }
    }
    return total;
}

std::size_t vat_directory_references(const VAT& index) {
    std::size_t count = 0;
    for (const auto& grid : index.directory.groups) for (const auto& row : grid)
        for (const auto& refs : row) count += refs.size();
    return count;
}

QueryResult query(const VAT& index, AttributePair pair, Count a, Count b, TimeWindow window, bool sort_result) {
    check_core(a, b);
    QueryResult result;
    if (!a || !b || !pair.u || !pair.v || window.start > window.end) return result;

    if (index.directory_ready) {
        auto answer = directory_query(index, pair, a, b, window);
        return sort_result ? canonical(std::move(answer)) : std::move(answer);
    }

    for (bool u_side : {true, false}) {
        const auto& vertices = u_side ? index.u_vertices : index.v_vertices;
        auto& output = u_side ? result.u_vertices : result.v_vertices;

        for (const auto& vertex : vertices) {
            if (vertex.attributes.empty() || less(pair, vertex.attributes.front().attributes) ||
                less(vertex.attributes.back().attributes, pair)) continue;

            const auto attr = std::lower_bound(vertex.attributes.begin(), vertex.attributes.end(), pair,
                [](const VATAttributes& node, AttributePair key) { return less(node.attributes, key); });

            if (attr == vertex.attributes.end() || !same(attr->attributes, pair) || a < attr->alpha_begin ||
                std::size_t(a - attr->alpha_begin) >= attr->alpha_rows.size()) continue;

            const auto& row = attr->alpha_rows[a - attr->alpha_begin];
            if (b < row.beta_begin || std::size_t(b - row.beta_begin) >= row.beta_lists.size()) continue;
            const auto& intervals = row.beta_lists[b - row.beta_begin];

            if (intervals.empty() || window.end < intervals.front().end || window.start > intervals.back().start) continue;

            const auto hit = std::lower_bound(intervals.begin(), intervals.end(), window.start,
                [](Interval interval, Timestamp start) { return interval.start < start; });

            if (hit != intervals.end() && hit->end <= window.end) output.push_back(vertex.original_id);
        }
    }

    if (sort_result) return canonical(std::move(result));
    return result;
}

QueryResult query(const UAT& index, AttributePair pair, Count a, Count b, TimeWindow window, bool sort_result) {
    check_core(a, b);
    QueryResult result;
    if (!a || !b || !pair.u || !pair.v || window.start > window.end) return result;
    const auto core = std::lower_bound(index.cores.begin(), index.cores.end(), std::make_pair(a, b),
        [](const UATCore& node, const auto& key) { return std::make_pair(node.alpha, node.beta) < key; });
    if (core == index.cores.end() || core->alpha != a || core->beta != b) return result;
    const auto attr = std::lower_bound(core->attributes.begin(), core->attributes.end(), pair,
        [](const UATAttributes& node, AttributePair key) { return less(node.attributes, key); });
    if (attr == core->attributes.end() || !same(attr->attributes, pair)) return result;
    const auto effective = std::lower_bound(index.time_points.begin(), index.time_points.end(), window.start);
    if (effective == index.time_points.end() || *effective > window.end) return result;
    const auto start = std::lower_bound(attr->starts.begin(), attr->starts.end(), *effective,
        [](const UATStart& node, Timestamp time) { return node.start < time; });
    if (start == attr->starts.end() || start->start != *effective) return result;
    const auto last = std::upper_bound(start->ends.begin(), start->ends.end(), window.end,
        [](Timestamp time, const UATEnd& node) { return time < node.end; });
    for (auto end = start->ends.begin(); end != last; ++end) {
        result.u_vertices.insert(result.u_vertices.end(), end->u_vertices.begin(), end->u_vertices.end());
        result.v_vertices.insert(result.v_vertices.end(), end->v_vertices.begin(), end->v_vertices.end());
    }
    if (sort_result) return canonical(std::move(result));
    return result;
}

std::size_t structure_bytes(const VAT& index) {
    std::size_t total = sizeof(index);
    for (const auto* side : {&index.u_vertices, &index.v_vertices}) {
        total += bytes(*side);
        for (const auto& vertex : *side) {
            total += bytes(vertex.attributes);
            for (const auto& attr : vertex.attributes) {
                total += bytes(attr.alpha_rows);
                for (const auto& row : attr.alpha_rows) {
                    total += bytes(row.beta_lists);
                    for (const auto& intervals : row.beta_lists) total += bytes(intervals);
                }
            }
        }
    }
    if (index.directory_ready) total += vat_directory_bytes(index) - sizeof(index.directory);
    return total;
}
std::size_t structure_bytes(const UAT& index) {
    std::size_t total = sizeof(index) + bytes(index.cores) + bytes(index.time_points);
    for (const auto& core : index.cores) {
        total += bytes(core.attributes);
        for (const auto& attr : core.attributes) {
            total += bytes(attr.starts);
            for (const auto& start : attr.starts) {
                total += bytes(start.ends);
                for (const auto& end : start.ends) total += bytes(end.u_vertices) + bytes(end.v_vertices);
            }
        }
    }
    return total;
}

void save_binary(const std::filesystem::path& path, const TemporalBipartiteGraph& graph, const VAT& index) {
    auto out = binary(path, index.directory_ready ? "TACJVAT2" : "TACJVAT1", graph);
    for (const auto* side : {&index.u_vertices, &index.v_vertices}) {
        integer(out, side->size());
        for (const auto& vertex : *side) {
            integer(out, vertex.original_id); integer(out, vertex.attributes.size());
            for (const auto& attr : vertex.attributes) {
                pair(out, attr.attributes); integer(out, attr.alpha_begin); integer(out, attr.alpha_rows.size());
                for (const auto& row : attr.alpha_rows) {
                    integer(out, row.beta_begin); integer(out, row.beta_lists.size());
                    for (const auto& intervals : row.beta_lists) {
                        integer(out, intervals.size());
                        for (auto interval : intervals) { integer(out, interval.start); integer(out, interval.end); }
                    }
                }
            }
        }
    }
    if (index.directory_ready) {
        integer(out, index.directory.pairs.size());
        for (const auto key : index.directory.pairs) pair(out, key);
        for (const auto& grid : index.directory.groups) {
            integer(out, grid.size());
            for (const auto& row : grid) {
                integer(out, row.size());
                for (const auto& refs : row) {
                    integer(out, refs.size());
                    for (const auto ref : refs)
                        integer(out, std::uint64_t(ref.vertex_side) | (std::uint64_t(ref.attribute_slot) << 32));
                }
            }
        }
    }
    close(out);
}
void save_binary(const std::filesystem::path& path, const TemporalBipartiteGraph& graph, const UAT& index) {
    auto out = binary(path, "TACJUAT1", graph);
    integer(out, index.time_points.size());
    for (auto time : index.time_points) integer(out, time);
    integer(out, index.cores.size());
    for (const auto& core : index.cores) {
        integer(out, core.alpha); integer(out, core.beta); integer(out, core.attributes.size());
        for (const auto& attr : core.attributes) {
            pair(out, attr.attributes); integer(out, attr.starts.size());
            for (const auto& start : attr.starts) {
                integer(out, start.start); integer(out, start.ends.size());
                for (const auto& end : start.ends) {
                    integer(out, end.end); ids(out, end.u_vertices); ids(out, end.v_vertices);
                }
            }
        }
    }
    close(out);
}

void load_binary(const std::filesystem::path& path, const TemporalBipartiteGraph& graph, VAT& index) {
    BinaryReader in(path, "TACJVAT2", graph, true);
    const ReadContext context(graph);
    VAT loaded;
    for (bool u_side : {true, false}) {
        auto& vertices = u_side ? loaded.u_vertices : loaded.v_vertices;
        const auto count = in.count(16);
        BinaryReader::check(count == (u_side ? graph.num_u_vertices() : graph.num_v_vertices()),
                            "joint VAT vertex count mismatch");
        vertices.resize(count);
        for (std::size_t i = 0; i < count; ++i) {
            auto& vertex = vertices[i];
            vertex.original_id = in.value<VertexID>();
            BinaryReader::check(vertex.original_id == (u_side ? graph.original_u_id(i) : graph.original_v_id(i)),
                                "joint VAT vertex ID mismatch");
            vertex.attributes.resize(in.count(32));
            AttributePair previous{};
            for (auto& attr : vertex.attributes) {
                attr.attributes = context.attributes(in);
                BinaryReader::check(less(previous, attr.attributes), "unordered joint VAT attributes");
                previous = attr.attributes;
                attr.alpha_begin = in.value<Count>();
                const auto rows = in.count(16);
                BinaryReader::check(attr.alpha_begin && rows &&
                                    std::uint64_t(attr.alpha_begin) + rows - 1 <= context.max_alpha,
                                    "invalid joint VAT alpha range");
                attr.alpha_rows.resize(rows);
                for (std::size_t ai = 0; ai < rows; ++ai) {
                    auto& row = attr.alpha_rows[ai];
                    row.beta_begin = in.value<Count>();
                    const auto lists = in.count(8);
                    BinaryReader::check(lists ? row.beta_begin &&
                                        std::uint64_t(row.beta_begin) + lists - 1 <= context.max_beta : row.beta_begin == 0,
                                        "invalid joint VAT beta range");
                    row.beta_lists.resize(lists);
                    for (std::size_t bi = 0; bi < lists; ++bi) {
                        auto& intervals = row.beta_lists[bi];
                        intervals.resize(in.count(16));
                        BinaryReader::check(intervals.empty() || attr.alpha_begin + ai != 1 || row.beta_begin + bi != 1,
                                            "joint VAT stores excluded core (1,1)");
                        for (auto& window : intervals)
                            window = {in.value<Timestamp>(), in.value<Timestamp>()};
                    }
                }
            }
        }
    }
    if (!in.legacy_vat()) {
        auto& directory = loaded.directory;
        const auto pairs = in.count(24);
        directory.pairs.reserve(pairs);
        directory.groups.resize(pairs);
        AttributePair previous{};
        for (std::size_t i = 0; i < pairs; ++i) {
            const auto key = context.attributes(in);
            BinaryReader::check(less(previous, key), "unordered VAT directory attributes");
            directory.pairs.push_back(key);
            previous = key;
        }
        for (std::size_t k = 0; k < pairs; ++k) {
            auto& grid = directory.groups[k];
            const auto rows = in.count(8);
            BinaryReader::check(rows <= context.max_alpha + 1, "invalid VAT directory alpha range");
            grid.resize(rows);
            for (std::size_t a = 0; a < rows; ++a) {
                auto& row = grid[a];
                const auto lists = in.count(8);
                BinaryReader::check(lists <= context.max_beta + 1, "invalid VAT directory beta range");
                row.resize(lists);
                for (std::size_t b = 0; b < lists; ++b) {
                    auto& refs = row[b];
                    const auto n = in.count(8);
                    refs.reserve(n);
                    for (std::size_t r = 0; r < n; ++r) {
                        const auto packed = in.integer();
                        const VATDirectoryRef ref{static_cast<std::uint32_t>(packed), static_cast<std::uint32_t>(packed >> 32)};
                        const auto& vertices = (ref.vertex_side & vat_upper_flag) ? loaded.u_vertices : loaded.v_vertices;
                        const auto slot = ref.vertex_side & ~vat_upper_flag;
                        BinaryReader::check(slot < vertices.size() && ref.attribute_slot < vertices[slot].attributes.size(),
                                            "VAT directory position out of range");
                        const auto& attr = vertices[slot].attributes[ref.attribute_slot];
                        BinaryReader::check(a >= attr.alpha_begin && a - attr.alpha_begin < attr.alpha_rows.size(),
                                            "VAT directory alpha position out of range");
                        const auto& alpha = attr.alpha_rows[a - attr.alpha_begin];
                        BinaryReader::check(b >= alpha.beta_begin && b - alpha.beta_begin < alpha.beta_lists.size(),
                                            "VAT directory beta position out of range");
                        refs.push_back(ref);
                    }
                }
            }
        }
        loaded.directory_ready = true;
    }
    in.finish();
    index = std::move(loaded);
}

void load_binary(const std::filesystem::path& path, const TemporalBipartiteGraph& graph, UAT& index) {
    BinaryReader in(path, "TACJUAT1", graph);
    const ReadContext context(graph);
    UAT loaded;
    const auto times = in.count(8);
    BinaryReader::check(times == context.times.size(), "joint UAT time-point count mismatch");
    loaded.time_points.resize(times);
    for (std::size_t i = 0; i < times; ++i) {
        loaded.time_points[i] = in.value<Timestamp>();
        BinaryReader::check(loaded.time_points[i] == context.times[i], "joint UAT time axis mismatch");
    }
    loaded.cores.resize(in.count(24));
    std::pair<Count, Count> previous_core{};
    for (auto& core : loaded.cores) {
        core.alpha = in.value<Count>(); core.beta = in.value<Count>();
        const auto key = std::make_pair(core.alpha, core.beta);
        BinaryReader::check(core.alpha && core.beta && core.alpha <= context.max_alpha && core.beta <= context.max_beta &&
                            key != std::make_pair(Count{1}, Count{1}) && previous_core < key,
                            "invalid joint UAT core order");
        previous_core = key;
        core.attributes.resize(in.count(24));
        BinaryReader::check(!core.attributes.empty(), "empty joint UAT core");
        AttributePair previous_attr{};
        for (auto& attr : core.attributes) {
            attr.attributes = context.attributes(in);
            BinaryReader::check(less(previous_attr, attr.attributes), "unordered joint UAT attributes");
            previous_attr = attr.attributes;
            attr.starts.resize(in.count(16));
            BinaryReader::check(!attr.starts.empty(), "empty joint UAT attributes");
            for (std::size_t si = 0; si < attr.starts.size(); ++si) {
                auto& start = attr.starts[si];
                start.start = in.value<Timestamp>();
                BinaryReader::check(!si || attr.starts[si-1].start < start.start, "unordered joint UAT starts");
                start.ends.resize(in.count(24));
                BinaryReader::check(!start.ends.empty(), "empty joint UAT start");
                for (std::size_t ei = 0; ei < start.ends.size(); ++ei) {
                    auto& end = start.ends[ei];
                    end.end = in.value<Timestamp>();
                    BinaryReader::check(end.end >= start.start && (!ei || start.ends[ei-1].end < end.end),
                                        "unordered joint UAT ends");
                    end.u_vertices = read_ids(in, context.num_u);
                    end.v_vertices = read_ids(in, context.num_v);
                    BinaryReader::check(!end.u_vertices.empty() || !end.v_vertices.empty(), "empty joint UAT end block");
                }
            }
        }
    }
    in.finish();
    index = std::move(loaded);
}

void save_text(const std::filesystem::path& path, const TemporalBipartiteGraph& graph, const VAT& index,
               const std::filesystem::path& data, const std::string& generated_at, const std::string& build_method) {
    std::ofstream out(path);
    text_header(out, "VAT", graph, data, generated_at, structure_bytes(index), build_method);
    out << "\n[VAT] vertex -> attribute pair -> alpha array -> beta array -> minimal intervals\n";
    for (bool u_side : {true, false}) {
        out << '\n' << (u_side ? "[U_VERTICES]" : "[V_VERTICES]") << '\n';
        for (const auto& vertex : u_side ? index.u_vertices : index.v_vertices) {
            out << (u_side ? "U=" : "V=") << vertex.original_id << '\n';
            for (const auto& attr : vertex.attributes) {
                out << "  AU=" << names(graph, attr.attributes.u) << " AV=" << names(graph, attr.attributes.v)
                    << " alpha_begin=" << attr.alpha_begin << " alpha_slots=" << attr.alpha_rows.size() << '\n';
                for (std::size_t a = 0; a < attr.alpha_rows.size(); ++a) {
                    const auto& row = attr.alpha_rows[a];
                    out << "    alpha=" << attr.alpha_begin + a << " beta_begin=" << row.beta_begin
                        << " beta_slots=" << row.beta_lists.size() << '\n';
                    for (std::size_t b = 0; b < row.beta_lists.size(); ++b) {
                        out << "      core=(" << attr.alpha_begin + a << ',' << row.beta_begin + b << ") intervals={";
                        bool first = true;
                        for (auto interval : row.beta_lists[b]) {
                            if (!first) out << ',';
                            first = false;
                            out << '[' << interval.start << ',' << interval.end << ']';
                        }
                        out << "}\n";
                    }
                }
            }
        }
    }
    out << "\n[VAT_QUERY_DIRECTORY] ready=" << index.directory_ready
        << " reference_bytes=" << sizeof(VATDirectoryRef) << " references=" << vat_directory_references(index)
        << " directory_bytes=" << vat_directory_bytes(index) << '\n'
        << "# U[i]/V[i] decodes vertex_side; slots are zero-based.\n"
        << "# Text after # is display-only; no original IDs or ATF intervals are copied into references.\n";
    for (std::size_t k = 0; k < index.directory.pairs.size(); ++k) {
        const auto key = index.directory.pairs[k];
        out << "AU=" << names(graph, key.u) << " AL=" << names(graph, key.v) << '\n';
        const auto& grid = index.directory.groups[k];
        for (std::size_t a = 0; a < grid.size(); ++a)
            for (std::size_t b = 0; b < grid[a].size(); ++b) {
                if (grid[a][b].empty()) continue;
                out << "  core=(" << a << ',' << b << ")\n";
                for (const auto ref : grid[a][b]) {
                    const bool upper = ref.vertex_side & vat_upper_flag;
                    const auto slot = ref.vertex_side & ~vat_upper_flag;
                    const auto& vertex = (upper ? index.u_vertices : index.v_vertices)[slot];
                    out << "    (" << (upper ? 'U' : 'V') << '[' << slot << "], attribute_slot=" << ref.attribute_slot
                        << ")  # -> " << (upper ? 'u' : 'v') << vertex.original_id << "'s existing ATF\n";
                }
            }
    }
    close(out);
}
void save_text(const std::filesystem::path& path, const TemporalBipartiteGraph& graph, const UAT& index,
               const std::filesystem::path& data, const std::string& generated_at, const std::string& build_method) {
    std::ofstream out(path);
    text_header(out, "UAT", graph, data, generated_at, structure_bytes(index), build_method);
    out << "\n[UAT] exact core -> attribute pair -> start -> critical end -> U/V vertex blocks\n";
    for (const auto& core : index.cores) {
        out << "core=(" << core.alpha << ',' << core.beta << ")\n";
        for (const auto& attr : core.attributes) {
            out << "  AU=" << names(graph, attr.attributes.u) << " AV=" << names(graph, attr.attributes.v) << '\n';
            for (const auto& start : attr.starts) {
                out << "    start=" << start.start << '\n';
                for (const auto& end : start.ends) {
                    out << "      end=" << end.end;
                    for (bool u_side : {true, false}) {
                        out << (u_side ? " U={" : " V={");
                        const auto& vertices = u_side ? end.u_vertices : end.v_vertices;
                        for (std::size_t i = 0; i < vertices.size(); ++i) { if (i) out << ','; out << vertices[i]; }
                        out << '}';
                    }
                    out << '\n';
                }
            }
        }
    }
    close(out);
}
}
