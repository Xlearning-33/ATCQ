#include "index/tab_index.h"
#include "query/query_execution.h"
#include "index/index_io.h"

#include "query/online_query.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <utility>

namespace tacc {

namespace {

constexpr std::uint64_t kTABBinaryMagic = 0x544142494e444558ULL;
constexpr std::uint32_t kTABBinaryVersion = 2;

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

template <typename T>
bool write_vector(std::ostream& output, const std::vector<T>& values) {
    static_assert(std::is_trivially_copyable<T>::value,
                  "write_vector requires trivially copyable elements");
    const std::uint64_t size = values.size();
    if (!write_pod(output, size)) return false;
    if (!values.empty()) {
        output.write(reinterpret_cast<const char*>(values.data()),
                     static_cast<std::streamsize>(values.size() * sizeof(T)));
    }
    return output.good();
}

template <typename T>
bool read_vector(std::istream& input, std::vector<T>& values) {
    static_assert(std::is_trivially_copyable<T>::value,
                  "read_vector requires trivially copyable elements");
    std::uint64_t size = 0;
    if (!read_pod(input, size) ||
        size > static_cast<std::uint64_t>(
                   std::numeric_limits<std::size_t>::max() / sizeof(T))) {
        return false;
    }
    values.resize(static_cast<std::size_t>(size));
    if (!values.empty()) {
        input.read(reinterpret_cast<char*>(values.data()),
                   static_cast<std::streamsize>(values.size() * sizeof(T)));
    }
    return input.good();
}

bool write_key(std::ostream& output, const TABQueryKey& key) {
    return write_pod(output, key.u_attributes) &&
           write_pod(output, key.v_attributes) &&
           write_pod(output, key.start) && write_pod(output, key.end) &&
           write_pod(output, key.alpha) && write_pod(output, key.beta);
}

bool read_key(std::istream& input, TABQueryKey& key) {
    return read_pod(input, key.u_attributes) &&
           read_pod(input, key.v_attributes) &&
           read_pod(input, key.start) && read_pod(input, key.end) &&
           read_pod(input, key.alpha) && read_pod(input, key.beta);
}

void build_lookup(TABIndex& index) {
    index.lookup.clear();
    index.lookup.reserve(index.records.size() * 2 + 1);
    for (std::size_t i = 0; i < index.records.size(); ++i) {
        index.lookup.emplace(index.records[i].key, i);
    }
}

std::uint64_t mix_hash(std::uint64_t hash, std::uint64_t value) {
    hash ^= value + 0x9e3779b97f4a7c15ULL + (hash << 6U) + (hash >> 2U);
    return hash;
}

AttributeMask make_attribute_mask(const std::vector<AttributeID>& attributes) {
    AttributeMask mask = 0;
    for (AttributeID attribute : attributes) {
        if (attribute >= 64) {
            throw std::runtime_error(
                "TAB-Index supports at most 64 distinct attributes");
        }
        mask |= AttributeMask{1} << attribute;
    }
    return mask;
}

std::vector<AttributeID> attributes_from_mask(AttributeMask mask) {
    std::vector<AttributeID> attributes;
    for (AttributeID attribute = 0; attribute < 64; ++attribute) {
        if ((mask & (AttributeMask{1} << attribute)) != 0) {
            attributes.push_back(attribute);
        }
    }
    return attributes;
}

std::vector<std::pair<AttributeMask, std::vector<AttributeID>>>
enumerate_nonempty_subsets(AttributeMask universe) {
    std::vector<std::pair<AttributeMask, std::vector<AttributeID>>> subsets;
    for (AttributeMask mask = universe; mask != 0;
         mask = (mask - 1) & universe) {
        subsets.emplace_back(mask, attributes_from_mask(mask));
    }
    return subsets;
}

std::pair<Count, Count> maximum_static_degrees(
    const TemporalBipartiteGraph& graph) {
    std::vector<std::unordered_set<VertexID>> u_neighbors(
        graph.num_u_vertices());
    std::vector<std::unordered_set<VertexID>> v_neighbors(
        graph.num_v_vertices());
    for (const TemporalEdge& edge : graph.edges()) {
        u_neighbors[edge.u].insert(edge.v);
        v_neighbors[edge.v].insert(edge.u);
    }

    Count max_alpha = 0;
    Count max_beta = 0;
    for (const auto& neighbors : u_neighbors) {
        max_alpha = std::max(max_alpha,
                             static_cast<Count>(neighbors.size()));
    }
    for (const auto& neighbors : v_neighbors) {
        max_beta = std::max(max_beta,
                            static_cast<Count>(neighbors.size()));
    }
    return {max_alpha, max_beta};
}

}

bool TABQueryKey::operator==(const TABQueryKey& other) const {
    return u_attributes == other.u_attributes &&
           v_attributes == other.v_attributes && start == other.start &&
           end == other.end && alpha == other.alpha && beta == other.beta;
}

std::size_t TABQueryKeyHash::operator()(const TABQueryKey& key) const {
    std::uint64_t hash = 1469598103934665603ULL;
    hash = mix_hash(hash, key.u_attributes);
    hash = mix_hash(hash, key.v_attributes);
    hash = mix_hash(hash, key.start);
    hash = mix_hash(hash, key.end);
    hash = mix_hash(hash, key.alpha);
    hash = mix_hash(hash, key.beta);
    return static_cast<std::size_t>(hash);
}

std::size_t tab_index_structure_bytes(const TABIndex& index) {
    std::size_t bytes = sizeof(TABIndex) + index.time_points.size() * sizeof(Timestamp) +
        index.records.size() * sizeof(TABRecord) +
        index.lookup.size() * sizeof(decltype(index.lookup)::value_type) +
        index.lookup.bucket_count() * sizeof(void*);
    for (const auto& record : index.records)
        bytes += (record.result.u_vertices.size() + record.result.v_vertices.size()) * sizeof(VertexID);
    return bytes;
}

TABIndex build_tab_index_exhaustive(const TemporalBipartiteGraph& graph,
                                   const TABBuildProgress& progress) {
    if (graph.num_attributes() > 64) {
        throw std::runtime_error(
            "TAB-Index supports at most 64 distinct attributes");
    }

    TABIndex index;
    index.source_graph_fingerprint = graph_fingerprint(graph);
    index.time_points = graph.time_points();

    for (VertexID u = 0; u < graph.num_u_vertices(); ++u) {
        index.u_attribute_universe |=
            AttributeMask{1} << graph.u_attribute(u);
    }
    for (VertexID v = 0; v < graph.num_v_vertices(); ++v) {
        index.v_attribute_universe |=
            AttributeMask{1} << graph.v_attribute(v);
    }

    const auto [max_alpha, max_beta] = maximum_static_degrees(graph);
    index.max_alpha = max_alpha;
    index.max_beta = max_beta;

    const auto u_subsets =
        enumerate_nonempty_subsets(index.u_attribute_universe);
    const auto v_subsets =
        enumerate_nonempty_subsets(index.v_attribute_universe);
    index.stats.attribute_pairs = u_subsets.size() * v_subsets.size();
    index.stats.windows =
        index.time_points.size() * (index.time_points.size() + 1) / 2;

    const auto build_start = std::chrono::steady_clock::now();
    std::size_t completed_pairs = 0;
    for (const auto& [u_mask, u_attributes] : u_subsets) {
        for (const auto& [v_mask, v_attributes] : v_subsets) {
            for (std::size_t start_index = 0;
                 start_index < index.time_points.size(); ++start_index) {
                for (std::size_t end_index = start_index;
                     end_index < index.time_points.size(); ++end_index) {
                    const TimeWindow window{index.time_points[start_index],
                                            index.time_points[end_index]};
                    for (Count alpha = 1; alpha <= index.max_alpha; ++alpha) {
                        for (Count beta = 1; beta <= index.max_beta; ++beta) {
                            if (alpha == 1 && beta == 1) continue;
                            ++index.stats.enumerated_queries;
                            const Query query{window, u_attributes, v_attributes,
                                              alpha, beta};
                            QueryResult result = online_query(graph, query);
                            if (result.u_vertices.empty() ||
                                result.v_vertices.empty()) {
                                continue;
                            }

                            TABQueryKey key{u_mask, v_mask, window.start,
                                            window.end, alpha, beta};
                            index.stats.vertex_references +=
                                result.u_vertices.size() +
                                result.v_vertices.size();
                            index.records.push_back(
                                TABRecord{key, std::move(result)});
                        }
                    }
                }
            }
            ++completed_pairs;
            if (progress) progress(completed_pairs, index.stats.attribute_pairs);
        }
    }

    build_lookup(index);
    index.stats.nonempty_records = index.records.size();
    index.stats.build_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - build_start).count();
    index.stats.approx_bytes = tab_index_structure_bytes(index);
    return index;
}

QueryResult query_tab_index(const TABIndex& index,
                            const TemporalBipartiteGraph& graph,
                            const Query& query,
                            const EncodedAttributes* attributes) {
    if (query.alpha == 1 && query.beta == 1) {
        return online_query(graph, query);
    }
    if (query.alpha == 0 || query.beta == 0 ||
        query.alpha > index.max_alpha || query.beta > index.max_beta ||
        query.window.start > query.window.end || index.time_points.empty()) {
        return {};
    }

    AttributeMask u_mask = (attributes ? attributes->u : make_attribute_mask(query.u_attributes)) &
                           index.u_attribute_universe;
    AttributeMask v_mask = (attributes ? attributes->v : make_attribute_mask(query.v_attributes)) &
                           index.v_attribute_universe;
    if (u_mask == 0 || v_mask == 0) return {};

    const auto start = std::lower_bound(index.time_points.begin(),
                                        index.time_points.end(),
                                        query.window.start);
    const auto after_end = std::upper_bound(index.time_points.begin(),
                                            index.time_points.end(),
                                            query.window.end);
    if (start == index.time_points.end() || after_end == index.time_points.begin()) {
        return {};
    }
    const auto end = std::prev(after_end);
    if (*start > *end) return {};

    const TABQueryKey key{u_mask, v_mask, *start, *end,
                          query.alpha, query.beta};
    const auto found = index.lookup.find(key);
    if (found == index.lookup.end()) return {};
    return index.records[found->second].result;
}

std::uint64_t tab_index_hash(const TABIndex& index) {
    std::uint64_t hash = 1469598103934665603ULL;
    hash = mix_hash(hash, index.u_attribute_universe);
    hash = mix_hash(hash, index.v_attribute_universe);
    hash = mix_hash(hash, index.max_alpha);
    hash = mix_hash(hash, index.max_beta);
    hash = mix_hash(hash, index.time_points.size());
    for (Timestamp timestamp : index.time_points) {
        hash = mix_hash(hash, timestamp);
    }
    hash = mix_hash(hash, index.records.size());
    for (const TABRecord& record : index.records) {
        hash = mix_hash(hash, record.key.u_attributes);
        hash = mix_hash(hash, record.key.v_attributes);
        hash = mix_hash(hash, record.key.start);
        hash = mix_hash(hash, record.key.end);
        hash = mix_hash(hash, record.key.alpha);
        hash = mix_hash(hash, record.key.beta);
        hash = mix_hash(hash, record.result.u_vertices.size());
        for (VertexID u : record.result.u_vertices) hash = mix_hash(hash, u);
        hash = mix_hash(hash, record.result.v_vertices.size());
        for (VertexID v : record.result.v_vertices) hash = mix_hash(hash, v);
    }
    return hash;
}

bool save_tab_index_binary(const std::string& path,
                           const TABIndex& index,
                           std::string* error) {
    std::ofstream output(path, std::ios::binary);
    if (!output) {
        set_error(error, "open_failed");
        return false;
    }

    const std::uint64_t enumerated_queries = index.stats.enumerated_queries;
    const std::uint64_t nonempty_records = index.stats.nonempty_records;
    const std::uint64_t vertex_references = index.stats.vertex_references;
    const std::uint64_t attribute_pairs = index.stats.attribute_pairs;
    const std::uint64_t windows = index.stats.windows;
    const std::uint64_t approx_bytes = index.stats.approx_bytes;

    bool ok = write_pod(output, kTABBinaryMagic) &&
              write_pod(output, index.source_graph_fingerprint ? kTABBinaryVersion : std::uint32_t{1});
    if (ok && index.source_graph_fingerprint)
        ok = write_pod(output, *index.source_graph_fingerprint);
    ok = ok &&
              write_pod(output, index.u_attribute_universe) &&
              write_pod(output, index.v_attribute_universe) &&
              write_pod(output, index.max_alpha) &&
              write_pod(output, index.max_beta) &&
              write_pod(output, index.stats.build_seconds) &&
              write_pod(output, enumerated_queries) &&
              write_pod(output, nonempty_records) &&
              write_pod(output, vertex_references) &&
              write_pod(output, attribute_pairs) &&
              write_pod(output, windows) &&
              write_pod(output, approx_bytes) &&
              write_vector(output, index.time_points);

    const std::uint64_t record_count = index.records.size();
    ok = ok && write_pod(output, record_count);
    if (ok) {
        for (const TABRecord& record : index.records) {
            if (!write_key(output, record.key) ||
                !write_vector(output, record.result.u_vertices) ||
                !write_vector(output, record.result.v_vertices)) {
                ok = false;
                break;
            }
        }
    }

    output.close();
    ok = ok && !output.fail();
    if (!ok) set_error(error, "write_failed");
    return ok;
}

bool load_tab_index_binary(const std::string& path,
                           TABIndex& index,
                           std::string* error) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        set_error(error, "open_failed");
        return false;
    }

    TABIndex loaded;
    std::uint64_t magic = 0;
    std::uint32_t version = 0;
    std::uint64_t record_count = 0;
    std::uint64_t enumerated_queries = 0;
    std::uint64_t nonempty_records = 0;
    std::uint64_t vertex_references = 0;
    std::uint64_t attribute_pairs = 0;
    std::uint64_t windows = 0;
    std::uint64_t approx_bytes = 0;
    bool ok = read_pod(input, magic) && read_pod(input, version);
    if (!ok || magic != kTABBinaryMagic || (version != 1 && version != kTABBinaryVersion)) {
        set_error(error, "invalid_format_or_version");
        return false;
    }

    if (version == 2) {
        std::uint64_t fingerprint{};
        ok = read_pod(input, fingerprint);
        if (ok) loaded.source_graph_fingerprint = fingerprint;
    }
    ok = ok && read_pod(input, loaded.u_attribute_universe) &&
         read_pod(input, loaded.v_attribute_universe) &&
         read_pod(input, loaded.max_alpha) &&
         read_pod(input, loaded.max_beta) &&
         read_pod(input, loaded.stats.build_seconds) &&
         read_pod(input, enumerated_queries) &&
         read_pod(input, nonempty_records) &&
         read_pod(input, vertex_references) &&
         read_pod(input, attribute_pairs) &&
         read_pod(input, windows) && read_pod(input, approx_bytes) &&
         read_vector(input, loaded.time_points) &&
         read_pod(input, record_count);
    if (!ok || record_count >
                   static_cast<std::uint64_t>(
                       std::numeric_limits<std::size_t>::max())) {
        set_error(error, "read_failed");
        return false;
    }

    loaded.stats.enumerated_queries =
        static_cast<std::size_t>(enumerated_queries);
    loaded.stats.nonempty_records =
        static_cast<std::size_t>(nonempty_records);
    loaded.stats.vertex_references =
        static_cast<std::size_t>(vertex_references);
    loaded.stats.attribute_pairs = static_cast<std::size_t>(attribute_pairs);
    loaded.stats.windows = static_cast<std::size_t>(windows);
    loaded.stats.approx_bytes = static_cast<std::size_t>(approx_bytes);

    loaded.records.resize(static_cast<std::size_t>(record_count));
    for (TABRecord& record : loaded.records) {
        if (!read_key(input, record.key) ||
            !read_vector(input, record.result.u_vertices) ||
            !read_vector(input, record.result.v_vertices)) {
            set_error(error, "read_failed");
            return false;
        }
    }

    std::size_t actual_refs = 0;
    for (const auto& record : loaded.records)
        actual_refs += record.result.u_vertices.size() + record.result.v_vertices.size();
    if (input.peek() != std::char_traits<char>::eof() || input.bad() ||
        loaded.records.size() != loaded.stats.nonempty_records ||
        actual_refs != loaded.stats.vertex_references) {
        set_error(error, "invalid_TAB_content_or_trailing_bytes");
        return false;
    }

    build_lookup(loaded);
    if (loaded.lookup.size() != loaded.records.size()) {
        set_error(error, "duplicate_TAB_query_keys");
        return false;
    }
    loaded.stats.nonempty_records = loaded.records.size();
    loaded.stats.approx_bytes = tab_index_structure_bytes(loaded);
    index = std::move(loaded);
    return true;
}

bool save_tab_index_text(const std::string& path,
                         const TemporalBipartiteGraph& graph,
                         const TABIndex& index,
                         std::string* error,
                         const std::string& generated_at,
                         const std::string& data_directory) {
    std::ofstream output(path);
    if (!output) {
        set_error(error, "open_failed");
        return false;
    }

    auto write_attributes = [&](AttributeMask mask) {
        output << '{';
        bool first = true;
        for (AttributeID attribute = 0;
             attribute < graph.num_attributes(); ++attribute) {
            if ((mask & (AttributeMask{1} << attribute)) == 0) continue;
            if (!first) output << ',';
            output << graph.attribute_name(attribute);
            first = false;
        }
        output << '}';
    };
    auto write_result_vertices = [&](const std::vector<VertexID>& vertices) {
        output << '{';
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            if (i != 0) output << ',';
            output << vertices[i];
        }
        output << '}';
    };

    output << "# TACC TAB-Index (human-readable)\n";
    if (!generated_at.empty()) output << "Generated at: " << generated_at << '\n';
    if (!data_directory.empty()) output << "Data directory: " << data_directory << '\n';
    output << "Build method: exhaustive\n"
           << "Structure bytes: " << tab_index_structure_bytes(index) << '\n'
           << "# Non-empty Online answers are stored; (1,1) is handled directly.\n"
           << "records=" << index.records.size()
           << " enumerated_queries=" << index.stats.enumerated_queries
           << " windows=" << index.stats.windows
           << " max_alpha=" << index.max_alpha
           << " max_beta=" << index.max_beta << "\n\n";

    for (std::size_t i = 0; i < index.records.size(); ++i) {
        const TABRecord& record = index.records[i];
        output << "record=" << i + 1 << " u_attributes=";
        write_attributes(record.key.u_attributes);
        output << " v_attributes=";
        write_attributes(record.key.v_attributes);
        output << " window=[" << record.key.start << ',' << record.key.end
               << "] alpha=" << record.key.alpha
               << " beta=" << record.key.beta << '\n';
        output << "  U=";
        write_result_vertices(record.result.u_vertices);
        output << '\n' << "  V=";
        write_result_vertices(record.result.v_vertices);
        output << "\n\n";
    }

    output.close();
    if (!output.good()) {
        set_error(error, "write_failed");
        return false;
    }
    return true;
}

}
