#include "experiments/online_baseline.h"
#include "index/index_io.h"
#include "query/online_query.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unistd.h>

namespace tacc::experiments {
namespace {
namespace fs = std::filesystem;
constexpr const char* version = "TACC_ONLINE_BASELINE2";
constexpr const char* online_version = "plain_time_scan_batch_500_v3";

void require(bool ok, const std::string& reason) {
    if (!ok) throw std::runtime_error(reason);
}
void require(bool ok, const char* reason) {
    if (!ok) throw std::runtime_error(reason);
}
void hash_bytes(std::uint64_t& hash, std::string_view bytes) {
    for (unsigned char ch : bytes) {
        hash ^= ch;
        hash *= 1099511628211ULL;
    }
}
std::string host() {
    std::array<char, 256> name{};
    require(gethostname(name.data(), name.size() - 1) == 0, "cannot identify baseline host");
    return name.data();
}
std::uint64_t number(const std::string& text) {
    require(!text.empty() && std::all_of(text.begin(), text.end(),
        [](char c) { return c >= '0' && c <= '9'; }), "invalid baseline integer");
    return std::stoull(text);
}
double duration(const std::string& text) {
    std::size_t used = 0;
    const double value = std::stod(text, &used);
    require(used == text.size() && std::isfinite(value) && value >= 0,
            "invalid baseline duration");
    return value;
}
std::uint64_t file_hash(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    require(bool(input), "missing baseline file: " + path.string());
    std::array<char, 65536> buffer{};
    std::uint64_t hash = 14695981039346656037ULL;
    while (input) {
        input.read(buffer.data(), buffer.size());
        hash_bytes(hash, std::string_view(buffer.data(), static_cast<std::size_t>(input.gcount())));
    }
    require(input.eof(), "cannot hash baseline file");
    return hash;
}
OnlineBaseline paths_for(const fs::path& queries) {
    const auto text = queries.string();
    const std::string suffix = "queries.tsv";
    require(text.size() >= suffix.size() &&
            text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0,
            "baseline query filename must end with queries.tsv");
    const auto prefix = text.substr(0, text.size() - suffix.size());
    OnlineBaseline baseline;
    baseline.queries_path = queries;
    baseline.answers_path = prefix + "answers.tsv";
    baseline.samples_path = prefix + "samples.csv";
    baseline.complete_path = prefix + "online.complete";
    return baseline;
}
std::vector<VertexID> read_ids(std::string_view text) {
    if (text == "-") return {};
    require(!text.empty() && text.back() != ',', "invalid baseline vertex list");
    std::vector<VertexID> vertices;
    const char* current = text.data();
    const char* end = current + text.size();
    while (current != end) {
        VertexID value = 0;
        const auto parsed = std::from_chars(current, end, value);
        require(parsed.ec == std::errc{} && (parsed.ptr == end || *parsed.ptr == ','),
                "invalid baseline vertex ID or overflow");
        require(vertices.empty() || vertices.back() < value, "baseline vertex list must be sorted and unique");
        vertices.push_back(value);
        current = parsed.ptr == end ? end : parsed.ptr + 1;
    }
    return vertices;
}

void write_ids(std::ostream& output, const std::vector<VertexID>& vertices) {
    if (vertices.empty()) output << '-';
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        if (i) output << ',';
        output << vertices[i];
    }
}

void close_checked(std::ofstream& output) {
    output.close();
    require(bool(output), "cannot finish writing Online baseline");
}

}

OnlineBaseline create_online_baseline(const fs::path& queries_path,
    const TemporalBipartiteGraph& graph, const Parameters& parameters,
    const std::vector<Query>& prepared_queries) {
    auto paths = paths_for(queries_path);
    for (const auto& path : {paths.queries_path, paths.answers_path, paths.samples_path, paths.complete_path})
        require(!fs::exists(path) && !fs::is_symlink(path), "refusing to overwrite Online baseline: " + path.string());
    const auto queries = prepared_queries.empty() ? generate_queries(graph, parameters) : prepared_queries;
    require(queries.size() == queries_per_setting, "need exactly 500 prepared queries");
    const auto fingerprint = graph_fingerprint(graph), hash = workload_hash(queries);
    fs::create_directories(queries_path.parent_path());
    std::ofstream conditions(paths.queries_path);
    require(bool(conditions), "cannot create query conditions");
    conditions << "TACC_QUERY1 graph_fingerprint=" << fingerprint << " query_hash=" << hash
               << " queries=" << queries.size() << " seed=" << workload_seed << '\n'
               << "id\tts\tte\talpha\tbeta\tu_mask\tv_mask\n";
    for (std::size_t i = 0; i < queries.size(); ++i) {
        const auto& query = queries[i];
        AttributeMask u_mask = 0, v_mask = 0;
        for (const auto a : query.u_attributes) u_mask |= AttributeMask{1} << a;
        for (const auto a : query.v_attributes) v_mask |= AttributeMask{1} << a;
        conditions << i + 1 << '\t' << query.window.start << '\t' << query.window.end
                   << '\t' << query.alpha << '\t' << query.beta << '\t' << u_mask << '\t' << v_mask << '\n';
    }
    close_checked(conditions);
    std::ofstream answers(paths.answers_path), samples(paths.samples_path);
    require(bool(answers) && bool(samples), "cannot create Online answers/timings");
    paths.queries = queries;
    paths.answers_hash = 14695981039346656037ULL;
    const auto write_answer = [&](const std::string& text) {
        answers.write(text.data(), static_cast<std::streamsize>(text.size()));
        hash_bytes(paths.answers_hash, text);
    };
    std::ostringstream header;
    header << "TACC_ANSWERS2 " << fingerprint << ' ' << hash << ' ' << queries.size()
           << "\nid\tu_ids\tv_ids\n";
    write_answer(header.str());
    samples << "id,empty,u_vertices,v_vertices\n";
    std::cout << "PROGRESS phase=online_baseline status=start queries=" << queries.size() << '\n' << std::flush;
    auto batch = run_timed_query_batch(queries.size(), [&](std::size_t i) {
        return online_query(graph, queries[i], nullptr, nullptr);
    });
    paths.groups[0].total_us = batch.total_us;
    for (std::size_t i = 0; i < queries.size(); ++i) {
        auto core = std::move(batch.answers[i]);
        std::sort(core.u_vertices.begin(), core.u_vertices.end());
        std::sort(core.v_vertices.begin(), core.v_vertices.end());
        require(core.u_vertices.empty() == core.v_vertices.empty(), "invalid Online core");
        std::ostringstream row;
        row << i + 1 << '\t';
        write_ids(row, core.u_vertices);
        row << '\t';
        write_ids(row, core.v_vertices);
        row << '\n';
        write_answer(row.str());
        samples << i + 1 << ',' << core.u_vertices.empty() << ','
                << core.u_vertices.size() << ',' << core.v_vertices.size() << '\n';
        paths.answer_counts.push_back({core.u_vertices.size(), core.v_vertices.size()});
        for (const auto group : {0u, core.u_vertices.empty() ? 1u : 2u}) {
            ++paths.groups[group].count;
        }
        if ((i + 1) % 100 == 0)
            std::cout << "PROGRESS phase=online_baseline completed=" << i + 1
                      << " total=" << queries.size() << '\n' << std::flush;
    }
    close_checked(answers);
    close_checked(samples);
    std::ofstream complete(paths.complete_path);
    require(bool(complete), "cannot create Online completion marker");
    complete << std::setprecision(std::numeric_limits<double>::max_digits10)
             << version << " status=ok graph=" << fingerprint << " query_hash=" << hash
             << " queries=" << queries.size() << " seed=" << workload_seed << " host=" << host()
             << " online_version=" << online_version << " online_batch_total_us=" << batch.total_us
             << " clock_reads=2 queries_hash=" << file_hash(paths.queries_path)
             << " answers_hash=" << paths.answers_hash
             << " samples_hash=" << file_hash(paths.samples_path) << '\n';
    close_checked(complete);
    paths.freshly_measured = true;
    return paths;
}

OnlineBaseline load_online_baseline(const fs::path& queries_path,
    const TemporalBipartiteGraph& graph, const Parameters& parameters) {
    auto baseline = paths_for(queries_path);
    std::ifstream complete(baseline.complete_path);
    std::string tag, token;
    require(bool(complete >> tag), "Online baseline is missing or incomplete");
    require(tag == version, "Online baseline timing model mismatch; regenerate Online timings for batch_500_v3");
    std::map<std::string, std::string> fields;
    while (complete >> token) {
        const auto pos = token.find('=');
        require(pos != std::string::npos &&
                fields.emplace(token.substr(0, pos), token.substr(pos + 1)).second,
                "invalid baseline completion marker");
    }
    require(complete.eof() && fields.at("status") == "ok" &&
            number(fields.at("graph")) == graph_fingerprint(graph) &&
            number(fields.at("queries")) == queries_per_setting &&
            number(fields.at("seed")) == workload_seed, "Online baseline metadata mismatch");
    require(fields.at("online_version") == online_version,
            "Online baseline timing model mismatch; regenerate Online timings for batch_500_v3");
    require(fields.at("clock_reads") == "2", "invalid Online batch clock count");
    baseline.groups[0].total_us = duration(fields.at("online_batch_total_us"));
    require(fields.at("host") == host(), "baseline host differs; prepare a new exp1 baseline on this host");
    baseline.queries = load_queries(queries_path, graph, parameters);
    require(number(fields.at("query_hash")) == workload_hash(baseline.queries) &&
            number(fields.at("queries_hash")) == file_hash(queries_path) &&
            number(fields.at("samples_hash")) == file_hash(baseline.samples_path),
            "Online baseline file integrity mismatch");
    baseline.answers_hash = number(fields.at("answers_hash"));
    require(fs::is_regular_file(baseline.answers_path), "missing baseline file: " + baseline.answers_path.string());
    std::ifstream samples(baseline.samples_path);
    std::string line;
    require(bool(std::getline(samples, line)) && line == "id,empty,u_vertices,v_vertices",
            "invalid Online samples header");
    while (std::getline(samples, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        std::string id, empty, u, v, extra;
        require(bool(row >> id >> empty >> u >> v) && !(row >> extra) &&
                number(id) == baseline.answer_counts.size() + 1 && number(id) <= queries_per_setting,
                "invalid Online sample row");
        BaselineAnswerCounts timing{static_cast<std::size_t>(number(u)),
                              static_cast<std::size_t>(number(v))};
        require(timing.u_count <= graph.num_u_vertices() && timing.v_count <= graph.num_v_vertices(),
                "baseline answer count exceeds the graph");
        require((empty == "0" || empty == "1") && (timing.u_count == 0) == (timing.v_count == 0) &&
                (empty == "1") == (timing.u_count == 0), "baseline sample class mismatch");
        baseline.answer_counts.push_back(timing);
        for (const auto group : {0u, empty == "1" ? 1u : 2u}) {
            ++baseline.groups[group].count;
        }
    }
    require(samples.eof() && baseline.answer_counts.size() == queries_per_setting, "truncated Online samples");
    for (const auto& group : baseline.groups)
        require(std::isfinite(group.total_us), "baseline duration sum overflow");
    return baseline;
}

BaselineAnswerReader::BaselineAnswerReader(const OnlineBaseline& baseline,
                                          const TemporalBipartiteGraph& graph)
    : baseline_(baseline), input_(baseline.answers_path, std::ios::binary) {
    std::string tag, line, extra;
    std::uint64_t fingerprint = 0, hash = 0;
    std::size_t count = 0;
    require(read_line(line), "missing Online answers header");
    std::istringstream metadata(line);
    require(bool(metadata >> tag >> fingerprint >> hash >> count) && !(metadata >> extra) &&
            tag == "TACC_ANSWERS2" && fingerprint == graph_fingerprint(graph) &&
            hash == workload_hash(baseline.queries) && count == queries_per_setting,
            "Online answer metadata mismatch");
    require(read_line(line) && line == "id\tu_ids\tv_ids",
            "invalid Online answer columns");
}

bool BaselineAnswerReader::read_line(std::string& line) {
    if (!std::getline(input_, line)) return false;
    hash_bytes(answers_hash_, line);
    if (!input_.eof()) hash_bytes(answers_hash_, "\n");
    return true;
}

QueryResult BaselineAnswerReader::next() {
    std::string line, id, u, v, extra;
    require(position_ < baseline_.answer_counts.size() && read_line(line),
            "truncated Online answers");
    std::istringstream row(line);
    require(bool(row >> id >> u >> v) && !(row >> extra) && number(id) == position_ + 1,
            "invalid Online answer row");
    const auto& timing = baseline_.answer_counts[position_++];
    QueryResult answer{read_ids(u), read_ids(v)};
    require(answer.u_vertices.size() == timing.u_count &&
            answer.v_vertices.size() == timing.v_count, "Online answer/sample mismatch");
    return answer;
}

void BaselineAnswerReader::finish() {
    std::string extra;
    require(position_ == queries_per_setting && !read_line(extra) && input_.eof() &&
            answers_hash_ == baseline_.answers_hash,
            "Online baseline file integrity mismatch: checksum or extra answers");
}

}
