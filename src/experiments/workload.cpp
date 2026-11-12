#include "experiments/workload.h"
#include "index/index_io.h"

#include <algorithm>
#include <fstream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>

namespace tacc::experiments {
namespace {


void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}

std::size_t scaled(std::size_t total, unsigned percent) {
    require(total && total <= (std::numeric_limits<std::size_t>::max() - 50) / 100,
            "invalid parameter scale");
    return std::max<std::size_t>(1, std::min(total, (total * percent + 50) / 100));
}

std::vector<AttributeID> side_attributes(const TemporalBipartiteGraph& graph, bool u_side) {
    std::set<AttributeID> attributes;
    const auto count = u_side ? graph.num_u_vertices() : graph.num_v_vertices();
    for (VertexID v = 0; v < count; ++v) {
        const auto a = u_side ? graph.u_attribute(v) : graph.v_attribute(v);
        require(a < 64, "query file masks support attribute IDs 0--63");
        attributes.insert(a);
    }
    return {attributes.begin(), attributes.end()};
}

AttributeMask attribute_mask(const std::vector<AttributeID>& attributes) {
    AttributeMask mask = 0;
    for (const auto a : attributes) {
        require(a < 64, "attribute mask overflow");
        mask |= AttributeMask{1} << a;
    }
    return mask;
}

auto key(const Query& q) {
    return std::make_tuple(q.window.start, q.window.end, q.alpha, q.beta,
                          attribute_mask(q.u_attributes), attribute_mask(q.v_attributes));
}

std::uint64_t unsigned_number(const std::string& text) {
    require(!text.empty() && std::all_of(text.begin(), text.end(),
                [](char c) { return c >= '0' && c <= '9'; }), "invalid query integer");
    return std::stoull(text);
}

std::vector<AttributeID> selected_attributes(AttributeMask mask,
                                             const std::vector<AttributeID>& universe) {
    std::vector<AttributeID> selected;
    for (const auto a : universe) {
        if (mask & (AttributeMask{1} << a)) {
            selected.push_back(a);
            mask &= ~(AttributeMask{1} << a);
        }
    }
    require(!mask, "query contains an attribute from the wrong side or graph");
    return selected;
}

}

Parameters default_parameters(const TemporalBipartiteGraph& graph, Count delta) {
    require(delta > 0, "graph has no nonempty balanced core");
    return {static_cast<Count>(scaled(delta, 60)), static_cast<Count>(scaled(delta, 40)),
            scaled(graph.time_points().size(), 60),
            scaled(side_attributes(graph, true).size(), 40),
            scaled(side_attributes(graph, false).size(), 40)};
}

std::vector<Query> generate_queries(const TemporalBipartiteGraph& graph, const Parameters& p,
                                    std::uint64_t seed) {
    auto u = side_attributes(graph, true), v = side_attributes(graph, false);
    const auto& times = graph.time_points();
    require(p.alpha > 0 && p.beta > 0 && !p.excluded(), "invalid batch core thresholds");
    require(p.u_attributes > 0 && p.u_attributes <= u.size() &&
            p.v_attributes > 0 && p.v_attributes <= v.size() &&
            p.window_points > 0 && p.window_points <= times.size(), "invalid batch parameter sizes");
    std::mt19937_64 random(seed);
    std::uniform_int_distribution<std::size_t> start(0, times.size() - p.window_points);
    std::vector<Query> queries;
    queries.reserve(queries_per_setting);
    for (std::size_t i = 0; i < queries_per_setting; ++i) {
        const auto first = start(random);
        std::shuffle(u.begin(), u.end(), random);
        std::shuffle(v.begin(), v.end(), random);
        std::vector<AttributeID> au(u.begin(), u.begin() + p.u_attributes);
        std::vector<AttributeID> av(v.begin(), v.begin() + p.v_attributes);
        std::sort(au.begin(), au.end());
        std::sort(av.begin(), av.end());
        queries.push_back({{times[first], times[first + p.window_points - 1]},
                           std::move(au), std::move(av), p.alpha, p.beta});
    }
    return queries;
}

std::vector<Query> generate_parameter_queries(const TemporalBipartiteGraph& graph, const Parameters& p) {
    auto u = side_attributes(graph, true), v = side_attributes(graph, false);
    const auto& times = graph.time_points();
    require(p.alpha > 0 && p.beta > 0 && !p.excluded(), "invalid parameter core thresholds");
    require(p.u_attributes > 0 && p.u_attributes <= u.size() &&
            p.v_attributes > 0 && p.v_attributes <= v.size() &&
            p.window_points > 0 && p.window_points <= times.size(), "invalid parameter sizes");
    const auto starts = times.size() - p.window_points + 1;
    std::mt19937_64 random(workload_seed);
    std::vector<Query> queries;
    queries.reserve(queries_per_setting);
    for (std::size_t i = 0; i < queries_per_setting; ++i) {
        const double quantile = std::generate_canonical<double, 53>(random);
        const auto first = std::min(starts - 1, static_cast<std::size_t>(quantile * starts));
        std::shuffle(u.begin(), u.end(), random);
        std::shuffle(v.begin(), v.end(), random);
        std::vector<AttributeID> au(u.begin(), u.begin() + p.u_attributes);
        std::vector<AttributeID> av(v.begin(), v.begin() + p.v_attributes);
        std::sort(au.begin(), au.end());
        std::sort(av.begin(), av.end());
        queries.push_back({{times[first], times[first + p.window_points - 1]},
                           std::move(au), std::move(av), p.alpha, p.beta});
    }
    return queries;
}

std::vector<ParameterSetting> parameter_settings(const TemporalBipartiteGraph& graph,
                                                Count delta) {
    const auto base = default_parameters(graph, delta);
    std::vector<ParameterSetting> settings;
    for (const std::string family : {"alpha_beta", "alpha", "beta"}) {
        for (unsigned percent : {10, 30, 50, 70, 90}) {
            auto p = base;
            if (family != "beta") p.alpha = static_cast<Count>(scaled(delta, percent));
            if (family != "alpha") p.beta = static_cast<Count>(scaled(delta, percent));
            settings.push_back({family, percent, p});
        }
    }
    for (const std::string family : {"u_attributes", "v_attributes", "uv_attributes"}) {
        for (unsigned percent : {20, 40, 60, 80}) {
            auto p = base;
            if (family != "v_attributes")
                p.u_attributes = scaled(side_attributes(graph, true).size(), percent);
            if (family != "u_attributes")
                p.v_attributes = scaled(side_attributes(graph, false).size(), percent);
            settings.push_back({family, percent, p});
        }
    }
    for (unsigned percent : {20, 40, 60, 80}) {
        auto p = base;
        p.window_points = scaled(graph.time_points().size(), percent);
        settings.push_back({"window", percent, p});
    }
    return settings;
}

std::uint64_t workload_hash(const std::vector<Query>& queries) {
    std::uint64_t hash = 14695981039346656037ULL;
    for (const auto& q : queries) {
        for (std::uint64_t value : {std::uint64_t(q.window.start), std::uint64_t(q.window.end),
                std::uint64_t(q.alpha), std::uint64_t(q.beta),
                attribute_mask(q.u_attributes), attribute_mask(q.v_attributes)}) {
            for (unsigned i = 0; i < 8; ++i) {
                hash ^= value & 255;
                hash *= 1099511628211ULL;
                value >>= 8;
            }
        }
    }
    return hash;
}

std::size_t distinct_queries(const std::vector<Query>& queries) {
    std::set<decltype(key(Query{}))> distinct;
    for (const auto& q : queries) distinct.insert(key(q));
    return distinct.size();
}

std::vector<Query> load_queries(const std::filesystem::path& path,
                                const TemporalBipartiteGraph& graph,
                                const Parameters& p) {
    std::ifstream input(path);
    std::string line, token;
    require(bool(std::getline(input, line)), "cannot read queries: " + path.string());
    std::istringstream metadata(line);
    metadata >> token;
    require(token == "TACC_QUERY1", "expected a new TACC_QUERY1 query file");
    std::map<std::string, std::string> fields;
    while (metadata >> token) {
        const auto separator = token.find('=');
        require(separator != std::string::npos &&
                fields.emplace(token.substr(0, separator), token.substr(separator + 1)).second,
                "invalid query metadata");
    }
    require(unsigned_number(fields.at("graph_fingerprint")) == graph_fingerprint(graph) &&
            unsigned_number(fields.at("queries")) == queries_per_setting,
            "query graph fingerprint or count mismatch");
    require(bool(std::getline(input, line)) &&
            line == "id\tts\tte\talpha\tbeta\tu_mask\tv_mask", "invalid query header");
    const auto u = side_attributes(graph, true), v = side_attributes(graph, false);
    const auto& times = graph.time_points();
    std::vector<Query> queries;
    while (std::getline(input, line)) {
        std::istringstream row(line);
        std::vector<std::uint64_t> values;
        while (row >> token) values.push_back(unsigned_number(token));
        require(values.size() == 7 && values[0] == queries.size() + 1 &&
                values[0] <= queries_per_setting, "invalid query row");
        require(values[1] <= values[2] && values[2] <= std::numeric_limits<Timestamp>::max() &&
                values[3] == p.alpha && values[4] == p.beta && !p.excluded(),
                "query thresholds do not match the requested setting");
        const auto first = std::lower_bound(times.begin(), times.end(), values[1]);
        const auto last = std::lower_bound(times.begin(), times.end(), values[2]);
        require(first != times.end() && last != times.end() && *first == values[1] &&
                *last == values[2] && std::size_t(last - first + 1) == p.window_points,
                "query window does not match the requested setting");
        auto au = selected_attributes(values[5], u), av = selected_attributes(values[6], v);
        require(au.size() == p.u_attributes && av.size() == p.v_attributes,
                "query attribute counts do not match the requested setting");
        queries.push_back({{static_cast<Timestamp>(values[1]), static_cast<Timestamp>(values[2])},
                           std::move(au), std::move(av), p.alpha, p.beta});
    }
    require(input.eof() && queries.size() == queries_per_setting, "incomplete query file");
    require(unsigned_number(fields.at("query_hash")) == workload_hash(queries), "query hash mismatch");
    return queries;
}

}
