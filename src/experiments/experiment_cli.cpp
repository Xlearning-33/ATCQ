#include "experiments/experiment_cli.h"
#include "experiments/experiment.h"
#include "experiments/dataset_catalog.h"
#include "experiments/dataset_stats.h"
#include "experiments/workload.h"
#include "experiments/online_baseline.h"
#include "experiments/mat_maintenance.h"
#include "common/memory_guard.h"
#include "index/index_builder.h"
#include "index/index_io.h"
#include "index/index_space.h"
#include "index/joint_index.h"
#include "index/tab_index.h"
#include "query/index_query_utils.h"
#include "query/mat_query.h"
#include "query/online_query.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <charconv>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <fcntl.h>
#include <sys/resource.h>
#include <unistd.h>

namespace tacc {
namespace {
namespace fs = std::filesystem;
namespace workload = experiments;
using Clock = std::chrono::steady_clock;


void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}
std::string wall_time(std::time_t value = std::time(nullptr)) {
    std::tm local {};
    require(localtime_r(&value, &local) != nullptr, "cannot read local time");
    std::ostringstream output;
    output << std::put_time(&local, "%Y-%m-%d %H:%M:%S");
    return output.str();
}
std::uint64_t peak_bytes() {
    struct rusage usage {};
    require(getrusage(RUSAGE_SELF, &usage) == 0, "cannot read process peak memory");
    return static_cast<std::uint64_t>(usage.ru_maxrss) * 1024;
}
struct Options {
    std::string experiment, task, selector, method{"mat"}, build_method{"opt"}, family{"all"};
    fs::path data, output, index, queries;
    std::string memory_percent;
    Query single_query;
    std::string u_names, v_names;
    std::string u_label_file{"graph_ulabel.txt"};
    std::string time_file{"graph_timestamp.txt"}, time_granularity;
    std::size_t u_label_count{}, v_label_count{};
    bool dump_text{};
    bool query_only{};
    std::string balanced_stage;
    std::size_t updates{10};
    std::uint64_t update_seed{20261031};
};

Options parse(int argc, char** argv) {
    Options o;
    require(argc >= 2, "missing experiment");
    o.experiment = argv[1];
    require(o.experiment.size() == 4 && o.experiment.substr(0, 3) == "exp" &&
            o.experiment[3] >= '0' && o.experiment[3] <= '9', "unknown experiment: " + o.experiment);
    if (o.experiment == "exp6") {
        o.method = "tab";
        o.build_method = "exhaustive";
    }

    int next = 3;
    if (o.experiment == "exp0") {
        require(argc >= 4, "expected exp0 <stats|online> <dataset>");
        o.task = argv[2]; o.selector = argv[3]; next = 4;
        require(o.task == "stats" || o.task == "online", "exp0 task must be stats or online");
        if (o.task == "online") {
            require(argc >= 10, "expected exp0 online <dataset> <alpha> <beta> <U_attrs> <V_attrs> <ts> <te>");
            const auto number = [](const std::string& text) {
                require(!text.empty() && std::all_of(text.begin(), text.end(),
                    [](char c) { return c >= '0' && c <= '9'; }), "query parameters must be unsigned integers");
                const auto value = std::stoull(text);
                require(value <= std::numeric_limits<std::uint32_t>::max(), "query integer overflow");
                return static_cast<std::uint32_t>(value);
            };
            o.single_query.alpha = number(argv[4]);
            o.single_query.beta = number(argv[5]);
            o.u_names = argv[6]; o.v_names = argv[7];
            o.single_query.window = {number(argv[8]), number(argv[9])};
            require(o.single_query.alpha > 0 && o.single_query.beta > 0, "alpha and beta must be positive");
            require(o.single_query.window.start <= o.single_query.window.end, "time window requires ts <= te");
            next = 10;
        }
    } else if (o.experiment == "exp5") {
        require(argc >= 5, "expected exp5 L <dataset> <U label count> or exp5 T <dataset> <time granularity>");
        o.task = argv[2]; o.selector = argv[3]; next = 4;
        require(o.task == "L" || o.task == "T", "exp5 task must be L or T");
        if (o.task == "L") {
            const std::string count = argv[4];
            require(!count.empty() && std::all_of(count.begin(), count.end(),
                [](char c) { return c >= '0' && c <= '9'; }), "U label count must be an integer from 1 to 64");
            o.u_label_count = std::stoull(count);
            require(o.u_label_count >= 1 && o.u_label_count <= 64, "U label count must be from 1 to 64");
            o.u_label_file = "graph_ulabel_" + std::to_string(o.u_label_count) + ".txt";
        } else {
            o.time_granularity = argv[4];
            const std::array<std::string, 9> grains{{"year", "halfyear", "quarter", "halfquarter",
                "month", "halfmonth", "week", "halfweek", "day"}};
            require(std::find(grains.begin(), grains.end(), o.time_granularity) != grains.end(),
                    "time granularity must be year, halfyear, quarter, halfquarter, month, halfmonth, week, halfweek or day");
            o.time_file = "graph_timestamp_" + o.time_granularity + ".txt";
        }
        next = 5;
    } else if (o.experiment == "exp6") {
        require(argc >= 3, "expected exp6 <build|query> <dataset>");
        const std::string task = argv[2];
        if (task == "build" || task == "query") {
            require(argc >= 4, "expected exp6 " + task + " <dataset>");
            o.task = task; o.selector = argv[3]; next = 4;
        } else {
            o.task = "build"; o.selector = task;
        }
    } else {
        require(argc >= 3, "expected " + o.experiment + " <dataset>");
        o.selector = argv[2];
        o.task = o.experiment == "exp3" ? "construction" : o.experiment == "exp4" ? "maintenance" : "query";
    }
    const bool tab_query = o.experiment == "exp6" && o.task == "query";
    const bool query_index = o.experiment == "exp1" || o.experiment == "exp2" || tab_query;

    std::set<std::string> supplied;
    for (int i = next; i < argc; ++i) {
        const std::string flag = argv[i];
        require(supplied.insert(flag).second, "duplicate option: " + flag);
        if (flag == "--query-only") {
            require(o.experiment == "exp5", "--query-only is only valid for exp5");
            o.query_only = true;
            continue;
        }
        if (flag == "--dump-text") {
            require(o.experiment == "exp3" || (o.experiment == "exp6" && o.task == "build"),
                    "--dump-text is only valid for exp3 or exp6 build");
            o.dump_text = true;
            continue;
        }
        require(flag == "--data-dir" || flag == "--output-root" || flag == "--index" ||
                flag == "--queries" || flag == "--method" || flag == "--build-method" ||
                flag == "--family" || flag == "--memory-percent" || flag == "--updates" || flag == "--seed" ||
                flag == "--balanced",
                "unknown option: " + flag);
        require(i + 1 < argc, "missing value: " + flag);
        const std::string value = argv[++i];
        require(!value.empty() && value.rfind("--", 0) != 0, "missing value: " + flag);
        if (flag == "--balanced") {
            require(o.experiment == "exp1" && (value == "prepare" || value == "query"),
                    "--balanced requires exp1 and prepare/query");
            o.balanced_stage = value;
        } else if (flag == "--data-dir") o.data = value;
        else if (flag == "--output-root") o.output = value;
        else if (flag == "--index") {
            require(query_index || o.experiment == "exp4" || o.experiment == "exp5",
                    "--index is only valid for query experiments or exp4");
            o.index = value;
        } else if (flag == "--queries") {
            require(o.experiment == "exp1" || tab_query || o.experiment == "exp5",
                    "--queries is only valid for exp1, exp5 --query-only or exp6 query");
            o.queries = value;
        } else if (flag == "--method") {
            require(o.experiment != "exp0", "exp0 does not take --method");
            if (o.experiment == "exp6") require(value == "tab", "exp6 method must be tab");
            else require(value == "mat" || value == "vat" || value == "uat", "method must be mat, vat or uat");
            require(o.experiment != "exp4" || value == "mat", "exp4 currently supports MAT maintenance only");
            o.method = value;
        } else if (flag == "--build-method") {
            require(o.experiment == "exp3" || o.experiment == "exp5",
                    "--build-method is only valid for exp3 or exp5");
            require(value == "opt" || value == "basic", "build method must be opt or basic");
            o.build_method = value;
        } else if (flag == "--updates" || flag == "--seed") {
            require(o.experiment == "exp4", flag + " is only valid for exp4");
            require(std::all_of(value.begin(), value.end(), [](char c) { return c >= '0' && c <= '9'; }),
                    flag + " must be an unsigned integer");
            const auto number = std::stoull(value);
            if (flag == "--updates") {
                require(number > 0 && number <= std::numeric_limits<std::size_t>::max(), "updates must be positive and fit size_t");
                o.updates = static_cast<std::size_t>(number);
            } else o.update_seed = number;
        } else if (flag == "--family") {
            require(o.experiment == "exp2", "--family is only valid for exp2");
            require(value == "all" || value == "alpha_beta" || value == "alpha" ||
                    value == "beta" || value == "u_attributes" || value == "v_attributes" ||
                    value == "uv_attributes" || value == "window", "unknown parameter family: " + value);
            o.family = value;
        } else {
            require(std::all_of(value.begin(), value.end(), [](char c) { return c >= '0' && c <= '9'; }) &&
                    std::stoul(value) <= 95, "memory percent must be 0--95");
            o.memory_percent = value;
        }
    }
    if (o.experiment == "exp5") {
        require(o.query_only || (o.index.empty() && o.queries.empty()),
                "exp5 --index/--queries require --query-only");
        require(!o.query_only || !supplied.count("--build-method"),
                "exp5 --query-only does not take --build-method");
    }
    require(fs::is_regular_file("CMakeLists.txt") && fs::is_directory("include/experiments"),
            "run from the project root directory");
    const auto& dataset = find_dataset(o.selector);
    o.selector = dataset.id;
    if (o.balanced_stage == "prepare") {
        if (!supplied.count("--method"))
            o.method = o.selector == "06" || o.selector == "08" || o.selector == "09" || o.selector == "10"
                ? "mat" : "uat";
        require(o.method == "mat" || o.method == "uat", "balanced sampling uses MAT or UAT");
        require(o.queries.empty(), "balanced prepare does not take --queries");
        o.task = "prepare";
    }
    o.data = fs::absolute(o.data.empty() ? fs::path(dataset.relative_directory) : o.data);
    o.output = fs::absolute(o.output.empty() ? fs::path("results/official") : o.output);
    return o;
}

std::string scalability_variant(const Options& o) {
    return o.task == "L" ? "u" + std::to_string(o.u_label_count) : o.time_granularity;
}

std::string experiment_directory(const Options& o) {
    if (o.experiment == "exp0") return "exp0_" + o.task;
    if (o.experiment == "exp1") return o.balanced_stage.empty() ? "exp1_query_default" : "exp1_query_5050";
    if (o.experiment == "exp2") return "exp2_query_vary";
    if (o.experiment == "exp3") return "exp3_construction";
    if (o.experiment == "exp4") return "exp4_maintenance";
    if (o.experiment == "exp5") return "exp5_scalability_" + o.task +
        (o.query_only ? "/query/" + scalability_variant(o) : "");
    if (o.experiment == "exp6") return "exp6_tab_comparison";
    throw std::invalid_argument("unsupported experiment output directory");
}

struct OutputPaths {
    fs::path directory;
    std::string prefix;
    std::string generated_at;
    std::string log_tag;

    fs::path log_file() const {
        return directory / (prefix + log_tag + '_' + directory.filename().string() + ".log");
    }
};

OutputPaths new_run(const Options& o, const DatasetEntry& dataset) {
    const auto base = o.output / dataset.name / experiment_directory(o);
    fs::create_directories(base);
    std::string dataset_tag = o.experiment == "exp0" && o.task == "stats" ? "_" + std::string(dataset.name) : "";
    if (o.experiment == "exp5")
        dataset_tag = "_" + std::string(dataset.name) + "_" +
            (o.task == "L" ? "u" + std::to_string(o.u_label_count) : o.time_granularity);
    for (unsigned attempt = 0; attempt < 100000; ++attempt) {
        const auto now = std::time(nullptr);
        std::tm local {};
        require(localtime_r(&now, &local) != nullptr, "cannot read local time");
        std::ostringstream name;
        name << std::put_time(&local, "%Y%m%d_%H%M%S") << dataset_tag;
        auto tag = o.experiment != "exp0"
            ? '_' + o.method : std::string{};
        if ((o.experiment == "exp6" && o.task == "query") || o.query_only) tag += "_query";
        if (o.balanced_stage == "prepare") tag += "_prepare";
        OutputPaths paths{base, name.str(), wall_time(now), tag};
        const auto log_path = paths.log_file();
        const int descriptor = ::open(log_path.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0644);
        if (descriptor < 0) {
            if (errno == EEXIST) {
                std::this_thread::sleep_until(std::chrono::system_clock::from_time_t(now + 1));
                continue;
            }
            throw std::runtime_error("cannot reserve output prefix: " + log_path.string());
        }
        require(::close(descriptor) == 0, "cannot close output reservation");
        if (o.experiment == "exp3" || (o.experiment == "exp5" && !o.query_only) ||
            (o.experiment == "exp6" && o.task == "build"))
            fs::create_directories(base / "indexes");
        return paths;
    }
    throw std::runtime_error("too many output name collisions");
}

TemporalBipartiteGraph load_graph(const fs::path& data, std::ostream& log,
                                const std::string& u_label_file = "graph_ulabel.txt",
                                const std::string& time_file = "graph_timestamp.txt") {
    const auto start = Clock::now();
    TemporalBipartiteGraph graph;
    require(graph.load_from_files((data / u_label_file).string(),
            (data / "graph_vlabel.txt").string(), (data / time_file).string()),
            "graph load failed");
    log << "GRAPH load_s=" << seconds(start) << " fingerprint=" << graph_fingerprint(graph) << '\n';
    return graph;
}

std::vector<AttributeID> parse_attributes(const TemporalBipartiteGraph& graph,
                                        const std::string& text, bool u_side) {
    require(!text.empty() && text.back() != ',', "attribute set must not be empty or end with a comma");
    std::set<AttributeID> allowed, selected;
    const auto count = u_side ? graph.num_u_vertices() : graph.num_v_vertices();
    for (VertexID vertex = 0; vertex < count; ++vertex)
        allowed.insert(u_side ? graph.u_attribute(vertex) : graph.v_attribute(vertex));
    std::istringstream input(text);
    std::string name;
    while (std::getline(input, name, ',')) {
        const auto first = name.find_first_not_of(" \t");
        require(first != std::string::npos, "empty attribute name");
        name = name.substr(first, name.find_last_not_of(" \t") - first + 1);
        const auto id = graph.attribute_id(name);
        require(id && allowed.count(*id), std::string(u_side ? "invalid U attribute: " : "invalid V attribute: ") + name);
        selected.insert(*id);
    }
    return {selected.begin(), selected.end()};
}

struct GroupTotal {
    std::size_t count{}, mismatches{};
    double online{}, index_us{};
};
struct Measurement {
    std::array<GroupTotal, 3> groups{};
    std::size_t distinct{}, mismatches{};
    std::size_t index_calls{};
    std::size_t online_calls{};
};

std::string average(double total, std::size_t count) {
    if (!count) return "NA";
    std::ostringstream value;
    value << std::fixed << std::setprecision(9) << total / count;
    return value.str();
}
std::string speedup(const GroupTotal& total) {
    if (!total.count || total.index_us <= 0) return "NA";
    std::ostringstream value;
    value << std::fixed << std::setprecision(9) << total.online / total.index_us;
    return value.str();
}

QueryResult index_answer(const MATIndex& index, const TemporalBipartiteGraph& graph, const Query& query,
                         const EncodedAttributes& attributes, QueryTimer* timer) {
    return mat_query(index, graph, query, nullptr, &attributes, timer);
}
QueryResult index_answer(const TABIndex& index, const TemporalBipartiteGraph& graph, const Query& query,
                         const EncodedAttributes& attributes, QueryTimer*) {
    return query_tab_index(index, graph, query, &attributes);
}
template<class JointIndex>
QueryResult index_answer(const JointIndex& index, const TemporalBipartiteGraph&, const Query& query,
                         const EncodedAttributes& attributes, QueryTimer*) {
    return joint::query(index, {attributes.u, attributes.v},
                        query.alpha, query.beta, query.window, false);
}

template<class QueryIndex>
Measurement measure(const TemporalBipartiteGraph& graph,
                    const workload::OnlineBaseline& baseline, const QueryIndex* index,
                    std::ostream& log, const std::string& case_id) {
    const auto& queries = baseline.queries;
    require(index != nullptr, "index stage requires an index");
    require(queries.size() == workload::queries_per_setting, "need exactly 500 queries");
    Measurement result;
    result.online_calls = baseline.freshly_measured ? queries.size() : 0;
    result.distinct = workload::distinct_queries(queries);
    result.groups[0].online = baseline.groups[0].total_us;

    std::vector<EncodedAttributes> attributes;
    attributes.reserve(queries.size());
    for (const auto& query : queries) attributes.push_back(encode_query_attributes(query));
    workload::BaselineAnswerReader answers(baseline, graph);
    auto batch = run_timed_query_batch(queries.size(), [&](std::size_t i) {
        return index_answer(*index, graph, queries[i], attributes[i], nullptr);
    });
    result.index_calls = queries.size();
    result.groups[0].index_us = batch.total_us;

    for (std::size_t i = 0; i < queries.size(); ++i) {
        auto reference = answers.next();
        auto actual = std::move(batch.answers[i]);
        std::sort(actual.u_vertices.begin(), actual.u_vertices.end());
        std::sort(actual.v_vertices.begin(), actual.v_vertices.end());
        require(reference.u_vertices.empty() == reference.v_vertices.empty(), "invalid cached Online core");
        const bool empty = reference.u_vertices.empty();
        const bool same = reference.u_vertices == actual.u_vertices &&
                          reference.v_vertices == actual.v_vertices;
        result.mismatches += !same;
        for (auto group : {0u, empty ? 1u : 2u}) {
            auto& total = result.groups[group];
            ++total.count;
            total.mismatches += !same;
        }
        if (!same) log << "MISMATCH case=" << case_id << " query_id=" << i + 1
                       << " expected_u=" << reference.u_vertices.size()
                       << " expected_v=" << reference.v_vertices.size()
                       << " actual_u=" << actual.u_vertices.size()
                       << " actual_v=" << actual.v_vertices.size() << '\n';
    }
    answers.finish();
    log << "BATCH_TIMING case=" << case_id << " queries=" << queries.size()
        << " clock_reads=2 per_query_timers=disabled total_us=" << batch.total_us << '\n';
    require(!result.mismatches, "index differs from saved Online answers; see MISMATCH lines in the query log");
    return result;
}

std::vector<fs::path> newest_files(const fs::path& directory, const std::string& suffix) {
    std::vector<fs::path> files;
    if (!fs::is_directory(directory)) return files;
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (!entry.is_regular_file()) continue;
        const auto name = entry.path().filename().string();
        if (name.size() >= suffix.size() &&
            name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0)
            files.push_back(entry.path());
    }
    std::sort(files.begin(), files.end(), [](const auto& a, const auto& b) {
        const auto ta = fs::last_write_time(a), tb = fs::last_write_time(b);
        return ta != tb ? ta > tb : a.filename() > b.filename();
    });
    return files;
}

std::string query_index_suffix(const std::string& method) {
    if (method == "tab") return "_tab_index.bin";
    return method == "mat" ? "_mat_opt_uv_index.bin" : '_' + method + "_opt_index.bin";
}

fs::path query_index_directory(const Options& o, const DatasetEntry& dataset) {
    return o.output / dataset.name / (o.experiment == "exp5"
        ? "exp5_scalability_" + o.task + "/indexes"
        : o.method == "tab" ? "exp6_tab_comparison/indexes" : "exp3_construction/indexes");
}

std::vector<fs::path> query_index_files(const Options& o, const DatasetEntry& dataset) {
    const auto suffix = (o.experiment == "exp5"
        ? "_" + std::string(dataset.name) + "_" + scalability_variant(o) : "") + query_index_suffix(o.method);
    return newest_files(query_index_directory(o, dataset), suffix);
}

workload::OnlineBaseline saved_baseline(const Options& o, const DatasetEntry& dataset,
    const TemporalBipartiteGraph& graph, const workload::Parameters& parameters, std::ostream& log,
    const fs::path& new_queries, const std::vector<Query>& prepared_queries = {}) {
    if (!o.queries.empty()) {
        auto baseline = workload::load_online_baseline(o.queries, graph, parameters);
        if (!prepared_queries.empty())
            require(workload::workload_hash(baseline.queries) == workload::workload_hash(prepared_queries),
                    "explicit baseline does not match selected queries");
        log << "BASELINE selection=explicit source=cached path=" << std::quoted(o.queries.string()) << '\n';
        return baseline;
    }
    std::vector<fs::path> directories;
    if (!new_queries.empty()) directories.push_back(new_queries.parent_path());
    if (o.experiment == "exp2" || o.experiment == "exp6")
        directories.push_back(o.output / dataset.name / "exp1_query_default/query_data");
    if (o.experiment == "exp1" && o.balanced_stage.empty())
        directories.push_back(o.output / dataset.name / "exp6_tab_comparison/query_data");
    if (o.experiment != "exp5" && o.balanced_stage.empty())
        directories.push_back(o.output / dataset.name / "exp0_online/query_data");
    const auto expected_hash = prepared_queries.empty() ? 0 : workload::workload_hash(prepared_queries);
    for (const auto& directory : directories) {
        for (const auto& candidate : newest_files(directory, "_queries.tsv")) {
            try {
                auto baseline = workload::load_online_baseline(candidate, graph, parameters);
                if (!prepared_queries.empty() && workload::workload_hash(baseline.queries) != expected_hash)
                    continue;
                log << "BASELINE selection=automatic source=cached path=" << std::quoted(candidate.string()) << '\n';
                return baseline;
            } catch (const std::exception&) {
            }
        }
    }
    if (!new_queries.empty()) {
        auto baseline = workload::create_online_baseline(new_queries, graph, parameters, prepared_queries);
        log << "BASELINE selection=generated source=fresh path=" << std::quoted(new_queries.string())
            << " queries=" << baseline.queries.size() << " online_executed_queries=" << baseline.queries.size() << '\n';
        return baseline;
    }
    throw std::runtime_error("missing output path for a new Online baseline");
}

void load_query_index(Options& o, const DatasetEntry& dataset,
                      const TemporalBipartiteGraph& graph, MATIndex& index, std::ostream& log) {
    const bool verify_content = o.experiment != "exp1" && o.experiment != "exp2" && !o.query_only;
    if (!o.index.empty()) {
        std::string error;
        require(load_index_binary(o.index.string(), graph, index, &error, verify_content), "cannot load MAT index: " + error);
        require(index.stats().method == BuildMethod::Opt && index.stats().scope == IndexScope::UV,
                "exp1/exp2 require an opt MAT index with UV scope");
        log << "INDEX selection=explicit path=" << std::quoted(o.index.string())
            << " content_validation=" << (verify_content ? "checksum" : "skipped") << '\n';
        return;
    }
    for (const auto& candidate : query_index_files(o, dataset)) {
        std::string error;
        if (load_index_binary(candidate.string(), graph, index, &error, verify_content) &&
            index.stats().method == BuildMethod::Opt && index.stats().scope == IndexScope::UV) {
            o.index = candidate;
            log << "INDEX selection=automatic path=" << std::quoted(candidate.string())
                << " content_validation=" << (verify_content ? "checksum" : "skipped") << '\n';
            return;
        }
        log << "INDEX_REJECT path=" << std::quoted(candidate.string()) << " reason=" << std::quoted(error) << '\n';
    }
    if (o.query_only) throw std::runtime_error("no matching opt MAT index for this scalability variant");
    throw std::runtime_error("no matching opt MAT index; run tacc exp3 " + std::string(dataset.id) + " first");
}

void load_query_index(Options& o, const DatasetEntry& dataset,
                      const TemporalBipartiteGraph& graph, TABIndex& index, std::ostream& log) {
    const auto expected = graph_fingerprint(graph);
    const auto load_one = [&](const fs::path& path, const char* selection) {
        TABIndex loaded;
        std::string error;
        const bool ok = load_tab_index_binary(path.string(), loaded, &error);
        require(ok, "cannot load TAB index: " + error);
        const bool legacy = !loaded.source_graph_fingerprint.has_value();
        if (legacy) {
            const auto name = path.filename().string();
            const auto prefix = name.substr(0, name.size() - std::string("_tab_index.bin").size());
            const auto build_log = path.parent_path().parent_path() / (prefix + "_tab_exp6_tab_comparison.log");
            std::ifstream input(build_log);
            std::optional<std::uint64_t> fingerprint;
            bool completed = false;
            std::string line;
            while (std::getline(input, line)) {
                if (line.rfind("GRAPH ", 0) == 0) {
                    std::istringstream fields(line);
                    std::string field;
                    while (fields >> field) {
                        if (field.rfind("fingerprint=", 0) != 0) continue;
                        const auto value = field.substr(12);
                        require(!value.empty() && std::all_of(value.begin(), value.end(),
                            [](char c) { return c >= '0' && c <= '9'; }), "invalid TAB build log fingerprint");
                        fingerprint = std::stoull(value);
                    }
                }
                if (line.rfind("RESULT exp=exp6 method=tab ", 0) == 0 &&
                    line.find(" status=ok") != std::string::npos) completed = true;
            }
            require(fingerprint.has_value() && completed,
                    "legacy TAB v1 requires its completed build log: " + build_log.string() +
                    "; otherwise rebuild with exp6 build " + dataset.id);
            require(*fingerprint == expected, "TAB graph_fingerprint_mismatch (legacy build log)");
        } else {
            require(*loaded.source_graph_fingerprint == expected, "TAB graph_fingerprint_mismatch");
        }
        require(loaded.time_points == graph.time_points(), "TAB time points do not match graph");
        index = std::move(loaded);
        o.index = path;
        log << "INDEX selection=" << selection << " path=" << std::quoted(path.string())
            << " graph_validation=" << (legacy ? "legacy_build_log" : "binary_fingerprint") << '\n';
    };
    if (!o.index.empty()) {
        load_one(o.index, "explicit");
        return;
    }
    const auto directory = o.output / dataset.name / "exp6_tab_comparison/indexes";
    for (const auto& candidate : newest_files(directory, "_tab_index.bin")) {
        try {
            load_one(candidate, "automatic");
            return;
        } catch (const std::exception& error) {
            log << "INDEX_REJECT path=" << std::quoted(candidate.string())
                << " reason=" << std::quoted(error.what()) << '\n';
        }
    }
    throw std::runtime_error("no matching TAB index; run tacc exp6 build " + std::string(dataset.id) + " first");
}

template<class JointIndex>
void load_query_index(Options& o, const DatasetEntry& dataset,
                      const TemporalBipartiteGraph& graph, JointIndex& index, std::ostream& log) {
    if (!o.index.empty()) {
        joint::load_binary(o.index, graph, index);
        log << "INDEX selection=explicit path=" << std::quoted(o.index.string()) << '\n';
        return;
    }
    for (const auto& candidate : query_index_files(o, dataset)) {
        try {
            joint::load_binary(candidate, graph, index);
            o.index = candidate;
            log << "INDEX selection=automatic path=" << std::quoted(candidate.string()) << '\n';
            return;
        } catch (const std::exception& error) {
            log << "INDEX_REJECT path=" << std::quoted(candidate.string())
                << " reason=" << std::quoted(error.what()) << '\n';
        }
    }
    if (o.query_only) throw std::runtime_error("no matching opt " + o.method + " index for this scalability variant");
    throw std::runtime_error("no matching opt " + o.method + " index; run tacc exp3 " + dataset.id +
                             " --method " + o.method + " first");
}

void write_result(std::ostream& log, const Options& o, const workload::Parameters& p,
                  const Measurement& result, const std::string& family,
                  unsigned percent, const std::string& case_id, bool reused) {
    const auto& all = result.groups[0];
    const auto& empty = result.groups[1];
    const auto& nonempty = result.groups[2];
    const unsigned alpha_ratio = family == "alpha" || family == "alpha_beta" ? percent : 60;
    const unsigned beta_ratio = family == "beta" || family == "alpha_beta" ? percent : 40;
    const unsigned window_ratio = family == "window" ? percent : 60;
    const unsigned u_ratio = family == "u_attributes" || family == "uv_attributes" ? percent : 40;
    const unsigned v_ratio = family == "v_attributes" || family == "uv_attributes" ? percent : 40;
    std::ostringstream heading;
    heading << "RESULT exp=" << o.experiment << " dataset=" << find_dataset(o.selector).name
        << " family=" << family << " percent=" << percent << " case=" << case_id
        << " alpha_ratio=" << alpha_ratio << " beta_ratio=" << beta_ratio
        << " window_ratio=" << window_ratio << " u_attribute_ratio=" << u_ratio
        << " v_attribute_ratio=" << v_ratio
        << " alpha=" << p.alpha << " beta=" << p.beta << " window_points=" << p.window_points
        << " u_attributes=" << p.u_attributes << " v_attributes=" << p.v_attributes
        << " seed=" << workload::workload_seed
        << " queries=" << all.count << " distinct_conditions=" << result.distinct
        << " empty=" << empty.count << " nonempty=" << nonempty.count;
    if (o.query_only) {
        heading << " task=query axis=" << o.task << " variant=" << scalability_variant(o)
                << " u_labels=" << o.u_label_count << " v_labels=" << o.v_label_count;
    }
    if (o.balanced_stage == "query") heading << " workload=balanced_5050_v1";
    for (const bool online : {true, false}) {
        log << heading.str() << " method=" << (online ? "Online" : o.method)
            << " source=" << (online ? (result.online_calls && !reused ? "fresh" : "cached") : "loaded")
            << " avg_total_us=" << average(online ? all.online : all.index_us, all.count)
            << " speedup_vs_online=" << (online ? "1" : speedup(all))
            << " empty_us=NA speedup_empty=NA nonempty_us=NA speedup_nonempty=NA"
            << " timing_model=" << batch_query_timing_model
            << " batch_total_us=" << (online ? all.online : all.index_us)
            << " clock_reads=" << (reused ? 0 : online ? (result.online_calls ? 2 : 0) : 2)
            << " group_timing=not_measured";
        if (!online)
            log << " pass=" << all.count - result.mismatches << " mismatches=" << result.mismatches;
        log << " reused=" << reused
            << " online_executed_queries=" << (online && !reused ? result.online_calls : 0)
            << " index_executed_queries=" << (online || reused ? 0 : result.index_calls)
            << " status=ok\n";
    }
}

std::array<std::string, 6> performance_table_row(
    const workload::ParameterSetting& setting, const Measurement* result) {
    const auto& p = setting.parameters;
    std::string actual;
    if (setting.family == "alpha") actual = std::to_string(p.alpha);
    else if (setting.family == "beta") actual = std::to_string(p.beta);
    else if (setting.family == "u_attributes") actual = std::to_string(p.u_attributes);
    else if (setting.family == "v_attributes") actual = std::to_string(p.v_attributes);
    else if (setting.family == "uv_attributes")
        actual = '(' + std::to_string(p.u_attributes) + ',' + std::to_string(p.v_attributes) + ')';
    else if (setting.family == "window") actual = std::to_string(p.window_points);
    else actual = '(' + std::to_string(p.alpha) + ',' + std::to_string(p.beta) + ')';
    std::array<std::string, 6> row{{setting.family == "default" ? "default" : std::to_string(setting.percent) + '%',
                                   actual, "-", "NA", "NA", "SKIP"}};
    if (!result) return row;
    const auto number = [](double value, int digits) {
        std::ostringstream output;
        output << std::fixed << std::setprecision(digits) << value;
        return output.str();
    };
    const auto& all = result->groups[0];
    row[2] = std::to_string(result->groups[1].count) + '/' + std::to_string(result->groups[2].count);
    if (all.count) {
        row[3] = number(all.online / all.count, 3);
        row[4] = number(all.index_us / all.count, 3);
    }
    row[5] = all.count && all.index_us > 0 ? number(all.online / all.index_us, 2) + 'x' : "NA";
    return row;
}

void write_performance_table(std::ostream& log, const std::string& family,
                             const std::vector<std::array<std::string, 6>>& rows,
                             const std::string& method = "mat") {
    const std::string label = method == "tab" ? "TAB" : method == "vat" ? "VAT" : method == "uat" ? "UAT" : "MAT";
    std::string parameter = "Core(a,b)";
    if (family == "alpha") parameter = "Alpha";
    else if (family == "beta") parameter = "Beta";
    else if (family == "u_attributes") parameter = "U_count";
    else if (family == "v_attributes") parameter = "V_count";
    else if (family == "uv_attributes") parameter = "U/V_count";
    else if (family == "window") parameter = "Window_points";
    const auto line = [&](const std::array<std::string, 6>& cells) {
        std::ostringstream output;
        const std::array<int, 6> widths{{9, 16, 18, 15, 15, 12}};
        for (std::size_t i = 0; i < cells.size(); ++i)
            output << std::left << std::setw(widths[i]) << cells[i];
        log << output.str() << '\n';
    };
    log << "\n# PERFORMANCE_TABLE family=" << family << " time_unit=us\n"
        << "# Online/" << label << ": mean query time; Speedup: Online total / " << label << " total.\n"
        << (family == "default" ? "# Default parameters are listed in PARAMETERS above.\n"
                               : "# SKIP: core(1,1); repeated configurations reuse the same measured result.\n");
    line({{family == "default" ? "Setting" : "Ratio", parameter, "Empty/Nonempty", "Online(us)", label + "(us)", "Speedup"}});
    log << std::string(85, '-') << '\n';
    for (const auto& row : rows) line(row);
    log << '\n';
}

template<class QueryIndex>
void default_queries(const Options& o, const TemporalBipartiteGraph& graph,
                     const workload::OnlineBaseline& baseline, const QueryIndex& index,
                     const workload::Parameters& parameters, std::ostream& log) {
    const auto result = measure(graph, baseline, &index, log, "default");
    write_result(log, o, parameters, result, "default", 0, "default", false);
    log << "QUERIES count=" << workload::queries_per_setting << " distinct=" << result.distinct
        << " empty=" << result.groups[1].count << " nonempty=" << result.groups[2].count
        << " query_hash=" << workload::workload_hash(baseline.queries)
        << " online_executed_queries=" << result.online_calls << '\n';
    write_performance_table(log, "default", {performance_table_row({"default", 0, parameters}, &result)}, o.method);
}

template<class QueryIndex>
void parameter_experiment(const Options& o, const TemporalBipartiteGraph& graph,
                          const QueryIndex& index,
                          const std::vector<workload::ParameterSetting>& settings,
                          const std::map<workload::Parameters, workload::OnlineBaseline>& baselines,
                          const workload::Parameters& defaults,
                          const OutputPaths& paths,
                          std::ostream& log) {
    std::map<workload::Parameters, std::pair<std::string, Measurement>> measured;
    std::string family;
    std::ofstream family_log;
    std::vector<std::array<std::string, 6>> table_rows;
    const auto close_family = [&] {
        if (!family_log.is_open()) return;
        write_performance_table(family_log, family, table_rows, o.method);
        table_rows.clear();
        family_log << "PROGRESS phase=parameter_family family=" << family
                   << " status=finished wall_time=" << wall_time() << " exit_code=0\n";
        family_log.close();
        require(bool(family_log), "cannot finish parameter family log: " + family);
    };
    for (const auto& setting : settings) {
        if (family != setting.family) {
            close_family();
            family = setting.family;
            const auto path = paths.directory / (paths.prefix + paths.log_tag + "_exp2_query_vary_" + family + ".log");
            require(!fs::exists(path) && !fs::is_symlink(path), "refusing to overwrite family log: " + path.string());
            family_log.open(path);
            require(bool(family_log), "cannot open parameter family log: " + path.string());
            const auto& base = defaults;
            family_log << std::fixed << std::setprecision(9) << std::unitbuf
                    << "TACC exp2 parameter variation\nGenerated at: " << paths.generated_at
                    << "\nCONFIG exp=exp2 dataset=" << find_dataset(o.selector).name << " family=" << family
                    << " method=" << o.method
                    << " queries_per_configuration=" << workload::queries_per_setting
                    << " seed=" << workload::workload_seed << " profile=0 rounds=1 warmup=0\n"
                    << "RUN_LOG path=" << std::quoted(paths.log_file().string()) << '\n'
                    << "INDEX source=shared built=0 path=" << std::quoted(o.index.string()) << '\n'
                    << "DEFAULT_PARAMETERS alpha=" << base.alpha << " beta=" << base.beta
                    << " window_points=" << base.window_points << " u_attributes=" << base.u_attributes
                    << " v_attributes=" << base.v_attributes << '\n'
                    << "# avg_total_us is in microseconds; repeated configurations reuse the same case.\n"
                    << "PROGRESS phase=parameter_family family=" << family
                    << " status=start wall_time=" << wall_time() << '\n';
            log << "QUERY_LOG family=" << family << " path=" << std::quoted(path.string()) << '\n';
            std::cout << "QUERY_LOG family=" << family << " path=" << std::quoted(path.string()) << '\n';
            family_log << "\n# ===== family=" << family << " =====\n";
        }
        const auto& p = setting.parameters;
        if (p.excluded()) {
            family_log << "RESULT exp=" << o.experiment << " family=" << setting.family
                << " percent=" << setting.percent << " alpha=1 beta=1 queries=0 status=skipped_11\n";
            table_rows.push_back(performance_table_row(setting, nullptr));
            continue;
        }
        auto found = measured.find(p);
        const bool reused = found != measured.end();
        if (!reused) {
            std::ostringstream name;
            name << "case_" << std::setw(3) << std::setfill('0') << measured.size() + 1;
            const auto result = measure(graph, baselines.at(p), &index, family_log, name.str());
            found = measured.emplace(p, std::make_pair(name.str(), result)).first;
            std::cout << "PROGRESS exp=" << o.experiment << " case=" << name.str() << " family=" << setting.family
                      << " status=ok\n" << std::flush;
        }
        const auto& result = found->second.second;
        write_result(family_log, o, p, result, setting.family, setting.percent,
                     found->second.first, reused);
        table_rows.push_back(performance_table_row(setting, &result));
        if (!reused)
            family_log << "QUERIES case=" << found->second.first
                << " count=" << workload::queries_per_setting << " distinct=" << result.distinct
                << " empty=" << result.groups[1].count << " nonempty=" << result.groups[2].count
                << " query_hash=" << workload::workload_hash(baselines.at(p).queries)
                << " online_executed_queries=" << result.online_calls << '\n';
    }
    close_family();
}

void write_construction_table(std::ostream& log, const Options& o, double build_s,
                              std::size_t index_bytes, std::uint64_t build_peak_bytes) {
    std::ostringstream table;
    const int build_width = o.method == "tab" ? 12 : 10;
    table << std::fixed << std::setprecision(3) << std::left
          << "\n# CONSTRUCTION_TABLE time_unit=s space_unit=MiB\n"
          << "# Index: structural space, not binary file size; Peak: process peak RSS through build, including graph.\n"
          << std::setw(16) << "Dataset" << std::setw(8) << "Index" << std::setw(build_width) << "Build"
          << std::setw(16) << "Build(s)" << std::setw(16) << "Index(MiB)" << std::setw(16) << "Peak(MiB)" << '\n'
          << std::string(72 + build_width, '-') << '\n'
          << std::setw(16) << find_dataset(o.selector).name
          << std::setw(8) << (o.method == "tab" ? "TAB" : o.method == "vat" ? "VAT" : o.method == "uat" ? "UAT" : "MAT")
          << std::setw(build_width) << o.build_method << std::setw(16) << build_s
          << std::setw(16) << space_mb(index_bytes) << std::setw(16) << space_mb(build_peak_bytes) << "\n\n";
    log << table.str();
}

void tab_construction(const Options& o, const TemporalBipartiteGraph& graph,
                      const OutputPaths& paths, std::ostream& log) {
    log << "PROGRESS phase=construction status=build_start wall_time=" << wall_time() << '\n';
    const char* flag = std::getenv("TAC_BUILD_PROGRESS");
    const bool show_progress = !flag || (std::string(flag) != "0" && std::string(flag) != "false" &&
                                         std::string(flag) != "FALSE");
    unsigned next_percent = 10;
    const auto start = Clock::now();
    auto index = build_tab_index_exhaustive(graph, [&](std::size_t completed, std::size_t total) {
        if (!show_progress || total == 0) return;
        while (next_percent <= 100 && 100.0L * completed / total >= next_percent) {
            std::ostringstream line;
            line << "PROGRESS method=tab build_method=exhaustive phase=attribute_pairs"
                 << " completed=" << completed << " total=" << total << " percent=" << next_percent
                 << " elapsed_s=" << seconds(start) << " wall_time=" << wall_time() << '\n';
            log << line.str();
            std::cout << line.str() << std::flush;
            next_percent += 10;
        }
    });
    const double build_s = seconds(start);
    const auto peak = peak_bytes();
    index.stats.build_seconds = build_s;
    const auto size = tab_index_structure_bytes(index);
    index.stats.approx_bytes = size;
    log << "PROGRESS phase=construction status=build_done wall_time=" << wall_time() << '\n';

    const auto stem = paths.prefix + "_tab_index";
    const auto binary = paths.directory / "indexes" / (stem + ".bin");
    std::string error;
    const bool saved = save_tab_index_binary(binary.string(), index, &error);
    require(saved, "TAB index save failed: " + error);
    if (o.dump_text || o.selector == "00") {
        const auto text = paths.directory / "indexes" / (stem + "_text.log");
        const bool text_saved = save_tab_index_text(text.string(), graph, index, &error,
                                                    paths.generated_at, o.data.string());
        require(text_saved, "TAB index text save failed: " + error);
    }
    log << "INDEX selection=constructed path=" << std::quoted(binary.string()) << '\n'
        << "RESULT exp=exp6 method=tab build_method=exhaustive build_s=" << build_s
        << " index_mb=" << space_mb(size) << " index_G=" << space_G(size)
        << " build_peak_mb=" << space_mb(peak) << " build_peak_G=" << space_G(peak)
        << " attribute_pairs=" << index.stats.attribute_pairs << " windows=" << index.stats.windows
        << " max_alpha=" << index.max_alpha << " max_beta=" << index.max_beta
        << " enumerated_queries=" << index.stats.enumerated_queries
        << " records=" << index.records.size() << " vertex_refs=" << index.stats.vertex_references
        << " nonempty_rate=" << (index.stats.enumerated_queries == 0 ? 0.0 :
            static_cast<double>(index.records.size()) / index.stats.enumerated_queries)
        << " index_bytes=" << size << " size_model=logical_elements_and_hash_buckets"
        << " peak_scope=process_through_build dump_excluded=1 status=ok\n";
    write_construction_table(log, o, build_s, size, peak);
    std::cout << "INDEX_FILE=" << binary.string() << '\n';
}

void construction(const Options& o, const TemporalBipartiteGraph& graph,
                  const OutputPaths& paths, std::ostream& log) {
    const auto method = o.build_method == "opt" ? BuildMethod::Opt : BuildMethod::Basic;
    const auto progress = make_index_build_progress_reporter("mat", method, log);
    log << "PROGRESS phase=construction status=build_start wall_time=" << wall_time() << '\n';
    const auto start = Clock::now();
    auto index = build_mat_index(graph, method, progress);
    const auto build_seconds = seconds(start);
    const auto peak = peak_bytes();
    log << "PROGRESS phase=construction status=build_done wall_time=" << wall_time() << '\n';
    const auto space = formal_index_space(index);
    std::string error;
    const auto index_stem = paths.prefix + "_mat_" + o.build_method + "_uv_index";
    const auto binary = paths.directory / "indexes" / (index_stem + ".bin");
    require(save_index_binary(binary.string(), graph, index, &error), "index save failed: " + error);
    if (o.dump_text || (o.experiment == "exp3" && o.selector == "00")) {
        require(save_index_text((paths.directory / "indexes" / (index_stem + "_text.log")).string(), graph, index, &error, o.data.string()),
                "index text save failed: " + error);
    }
    log << "RESULT exp=" << o.experiment << " method=mat";
    if (o.experiment == "exp5") {
        log << " axis=" << o.task << " u_labels=" << o.u_label_count << " v_labels=" << o.v_label_count;
        if (o.task == "T")
            log << " time_granularity=" << o.time_granularity
                << " time_points=" << graph.time_points().size() << " temporal_events=" << graph.num_edges();
    }
    log << " build_method=" << o.build_method << " build_s=" << build_seconds
        << " index_mb=" << space_mb(space.index_bytes()) << " index_G=" << space_G(space.index_bytes())
        << " build_peak_mb=" << space_mb(peak) << " build_peak_G=" << space_G(peak)
        << " peak_scope=process_through_build dump_excluded=1 status=ok\n";
    write_construction_table(log, o, build_seconds, space.index_bytes(), peak);
    std::cout << "INDEX_FILE=" << binary.string() << '\n';
}

void joint_construction(const Options& o, const TemporalBipartiteGraph& graph,
                        const OutputPaths& paths, std::ostream& log) {
    joint::BuildStats stats;
    const auto reporter = make_index_build_progress_reporter(o.method.c_str(),
        o.build_method == "opt" ? BuildMethod::Opt : BuildMethod::Basic, log);
    log << "PROGRESS phase=construction status=build_start wall_time=" << wall_time() << '\n';
    const auto construct_one = [&](auto builder) {
        const auto start = Clock::now();
        joint::BuildProgress progress;
        if (reporter) {
            progress = [&, start](std::size_t completed, std::size_t total, joint::AttributePair pair) {
                reporter({completed, total, Side::U, pair.u, seconds(start), true, pair.v});
            };
        }
        auto index = builder(graph, stats, progress);
        const double build_s = seconds(start);
        const auto peak = peak_bytes();
        log << "PROGRESS phase=construction status=build_done wall_time=" << wall_time() << '\n';
        const auto size = joint::structure_bytes(index);
        const auto stem = paths.prefix + '_' + o.method + '_' + o.build_method + "_index";
        const auto binary = paths.directory / "indexes" / (stem + ".bin");
        joint::save_binary(binary, graph, index);
        if (o.dump_text || (o.experiment == "exp3" && o.selector == "00"))
            joint::save_text(paths.directory / "indexes" / (stem + "_text.log"), graph, index, o.data, paths.generated_at, o.build_method);
        log << "RESULT exp=" << o.experiment << " method=" << o.method;
        if (o.experiment == "exp5") {
            log << " axis=" << o.task << " u_labels=" << o.u_label_count << " v_labels=" << o.v_label_count;
            if (o.task == "T")
                log << " time_granularity=" << o.time_granularity
                    << " time_points=" << graph.time_points().size() << " temporal_events=" << graph.num_edges();
        }
        log << " build_method=" << o.build_method << " scope=exact_joint_attributes"
            << " built_indexes=1 attribute_compression=0 core_1_1_stored=0"
            << " build_s=" << build_s << " index_mb=" << space_mb(size) << " index_G=" << space_G(size)
            << " build_peak_mb=" << space_mb(peak) << " build_peak_G=" << space_G(peak)
            << " peak_scope=process_through_build dump_excluded=1 online_calls=" << stats.online_calls
            << " full_window_prechecks=" << stats.full_window_prechecks
            << " global_valid_cores=" << stats.global_valid_cores << " joint_valid_cores=" << stats.joint_valid_cores
            << " attribute_pairs=" << stats.attribute_pairs << " vat_core_records=" << stats.vat_core_records
            << " vat_intervals=" << stats.vat_intervals << " uat_blocks=" << stats.uat_blocks
            << " vertex_refs=" << stats.vertex_refs << " decompositions=" << stats.decompositions
            << " windows=" << stats.windows << " dynamic_bucket_updates=" << stats.dynamic_bucket_updates
            << " deleted_unique_edges=" << stats.deleted_unique_edges << " frontier_changes=" << stats.frontier_changes
            << " emitted_memberships=" << stats.emitted_memberships << " empty_attribute_pairs=" << stats.empty_attribute_pairs;
        if constexpr (std::is_same_v<std::decay_t<decltype(index)>, joint::VAT>)
            log << " query_directory=" << index.directory_ready << " directory_ref_bytes=" << sizeof(joint::VATDirectoryRef)
                << " directory_refs=" << joint::vat_directory_references(index)
                << " directory_mb=" << space_mb(joint::vat_directory_bytes(index)) << " directory_in_build=1";
        log << " status=ok\n";
        write_construction_table(log, o, build_s, size, peak);
        std::cout << "INDEX_FILE=" << binary.string() << '\n';
    };
    if (o.method == "vat") {
        if (o.build_method == "opt") construct_one(joint::build_vat_opt);
        else construct_one(joint::build_vat_basic);
    } else {
        if (o.build_method == "opt") construct_one(joint::build_uat_opt);
        else construct_one(joint::build_uat_basic);
    }
}

struct BalancedSample {
    std::vector<Query> queries;
    std::string classes;
};

fs::path balanced_sample_path(const OutputPaths& paths) {
    return paths.directory / "samples/selected_queries.tsv";
}

BalancedSample read_balanced_sample(const fs::path& path, const TemporalBipartiteGraph& graph,
                                    const workload::Parameters& parameters) {
    BalancedSample sample;
    sample.queries = workload::load_queries(path, graph, parameters);
    std::ifstream input(path);
    std::string line, token;
    std::getline(input, line);
    std::istringstream header(line);
    std::map<std::string, std::string> fields;
    header >> token;
    while (header >> token) {
        const auto separator = token.find('=');
        fields[token.substr(0, separator)] = token.substr(separator + 1);
    }
    require(fields["workload"] == "balanced_5050_v1", "not a balanced 5050 sample");
    sample.classes = fields["classes"];
    require(sample.classes.size() == workload::queries_per_setting &&
            std::all_of(sample.classes.begin(), sample.classes.end(), [](char c) { return c == '0' || c == '1'; }) &&
            std::count(sample.classes.begin(), sample.classes.end(), '1') == 250,
            "balanced sample must contain 250 empty and 250 nonempty labels");
    return sample;
}

void prepare_balanced_sample(Options& o, const DatasetEntry& dataset,
                            const TemporalBipartiteGraph& graph, const OutputPaths& paths, std::ostream& log) {
    const auto parameters = workload::default_parameters(graph, static_cast<Count>(balanced_core_threshold(graph)));
    require(!parameters.excluded(), "balanced sampling excludes core (1,1)");
    const auto sample_path = balanced_sample_path(paths);
    if (fs::exists(sample_path)) {
        const auto sample = read_balanced_sample(sample_path, graph, parameters);
        log << "BALANCED_SAMPLE source=cached path=" << std::quoted(sample_path.string())
            << " empty=250 nonempty=250 queries=500 distinct=" << workload::distinct_queries(sample.queries)
            << " screening_calls=0 online_executed_queries=0 status=ok\n";
        return;
    }
    constexpr std::size_t candidate_limit = 100000;
    std::vector<Query> empty, nonempty;
    std::map<std::tuple<Timestamp, Timestamp, AttributeMask, AttributeMask>, bool> classifications;
    std::size_t candidates = 0, calls = 0;
    const auto collect = [&](const auto& index) {
        for (std::uint64_t batch = 0; candidates < candidate_limit &&
             (empty.size() < 250 || nonempty.size() < 250); ++batch) {
            const auto queries = workload::generate_queries(graph, parameters, workload::workload_seed + batch);
            for (const auto& query : queries) {
                if (candidates == candidate_limit || (empty.size() == 250 && nonempty.size() == 250)) break;
                ++candidates;
                const auto attributes = encode_query_attributes(query);
                const auto key = std::make_tuple(query.window.start, query.window.end, attributes.u, attributes.v);
                auto found = classifications.find(key);
                if (found == classifications.end()) {
                    auto answer = index_answer(index, graph, query, attributes, nullptr);
                    found = classifications.emplace(key, !answer.u_vertices.empty()).first;
                    ++calls;
                }
                auto& selected = found->second ? nonempty : empty;
                if (selected.size() < 250) selected.push_back(query);
            }
            if (candidates % 10000 == 0)
                log << "PROGRESS phase=balanced_sampling candidates=" << candidates << " screening_calls=" << calls
                    << " empty=" << empty.size() << " nonempty=" << nonempty.size() << '\n';
        }
    };
    const auto load_start = Clock::now();
    if (o.method == "mat") {
        MATIndex index;
        load_query_index(o, dataset, graph, index, log);
        log << "SAMPLING_INDEX method=mat load_s=" << seconds(load_start) << '\n';
        if (index_uses_grouped_recovery("mat")) (void)graph.grouped_v_time_adjacency();
        collect(index);
    } else {
        joint::UAT index;
        load_query_index(o, dataset, graph, index, log);
        log << "SAMPLING_INDEX method=uat load_s=" << seconds(load_start) << '\n';
        collect(index);
    }
    log << "BALANCED_SAMPLE source=generated screening_method=" << o.method << " candidate_limit=" << candidate_limit
        << " candidates=" << candidates << " screening_calls=" << calls << " distinct_candidates=" << classifications.size()
        << " empty=" << empty.size() << " nonempty=" << nonempty.size()
        << " online_executed_queries=0 status=" << (empty.size() == 250 && nonempty.size() == 250 ? "ready" : "insufficient") << '\n';
    require(empty.size() == 250 && nonempty.size() == 250,
            "insufficient balanced samples under unchanged default parameters; no formal queries were run");
    std::vector<std::pair<Query, char>> selected;
    for (auto& q : empty) selected.emplace_back(std::move(q), '0');
    for (auto& q : nonempty) selected.emplace_back(std::move(q), '1');
    std::mt19937_64 random(workload::workload_seed);
    std::shuffle(selected.begin(), selected.end(), random);
    BalancedSample sample;
    for (auto& entry : selected) {
        sample.queries.push_back(std::move(entry.first));
        sample.classes.push_back(entry.second);
    }
    fs::create_directories(sample_path.parent_path());
    const auto temporary = sample_path.parent_path() / (paths.prefix + "_selected.tmp");
    std::ofstream output(temporary);
    require(bool(output), "cannot save balanced sample");
    output << "TACC_QUERY1 graph_fingerprint=" << graph_fingerprint(graph)
           << " query_hash=" << workload::workload_hash(sample.queries) << " queries=500 seed=" << workload::workload_seed
           << " workload=balanced_5050_v1 screening_method=" << o.method << " classes=" << sample.classes
           << "\nid\tts\tte\talpha\tbeta\tu_mask\tv_mask\n";
    for (std::size_t i = 0; i < sample.queries.size(); ++i) {
        const auto& q = sample.queries[i];
        const auto attributes = encode_query_attributes(q);
        output << i + 1 << '\t' << q.window.start << '\t' << q.window.end << '\t' << q.alpha << '\t' << q.beta
               << '\t' << attributes.u << '\t' << attributes.v << '\n';
    }
    output.close();
    require(bool(output), "cannot finish balanced sample");
    fs::rename(temporary, sample_path);
    log << "BALANCED_SAMPLE source=saved path=" << std::quoted(sample_path.string())
        << " queries=500 empty=250 nonempty=250 distinct=" << workload::distinct_queries(sample.queries) << " status=ok\n";
    std::cout << "SAMPLE_FILE=" << sample_path.string() << '\n';
}

bool query_experiment(Options& o, const DatasetEntry& dataset,
                      const TemporalBipartiteGraph& graph,
                      const OutputPaths& paths, std::ostream& log, bool vary_parameters) {
    const auto delta_start = Clock::now();
    const auto delta = balanced_core_threshold(graph);
    log << "PARAMETER_SETUP delta_compute_s=" << seconds(delta_start) << '\n';
    require(delta <= std::numeric_limits<Count>::max(), "balanced core threshold overflow");
    const auto parameters = workload::default_parameters(graph, static_cast<Count>(delta));
    log << "PARAMETERS delta=" << delta << " alpha=" << parameters.alpha
        << " beta=" << parameters.beta << " window_points=" << parameters.window_points
        << " u_attributes=" << parameters.u_attributes << " v_attributes=" << parameters.v_attributes << '\n';
    workload::OnlineBaseline baseline;
    std::vector<workload::ParameterSetting> settings;
    std::map<workload::Parameters, workload::OnlineBaseline> baselines;
    if (vary_parameters) {
        settings = workload::parameter_settings(graph, static_cast<Count>(delta));
        settings.erase(std::remove_if(settings.begin(), settings.end(), [&](const auto& setting) {
            return o.family != "all" && setting.family != o.family;
        }), settings.end());
        log << "PARAMETER_PLAN family=" << o.family << " settings=" << settings.size()
            << " queries_per_configuration=" << workload::queries_per_setting
            << " sampling=paired_parameter_v1 duplicate_configurations=reused\n";
    }
    const bool need_queries = vary_parameters
        ? std::any_of(settings.begin(), settings.end(), [](const auto& setting) { return !setting.parameters.excluded(); })
        : !parameters.excluded();
    if (need_queries && o.index.empty()) {
        const auto files = query_index_files(o, dataset);
        if (o.method == "tab")
            require(!files.empty(), "no matching TAB index; run tacc exp6 build " + std::string(dataset.id) + " first");
        if (o.query_only)
            require(!files.empty(), "no matching opt " + o.method + " index for exp5 " + o.task +
                    " " + dataset.id + " " + scalability_variant(o) + "; construct this variant first");
        require(!files.empty(),
                "no matching opt " + (o.method == "mat" ? std::string("MAT") : o.method) +
                " index; run tacc exp3 " + dataset.id + " --method " + o.method + " first");
    }
    if (!vary_parameters && !parameters.excluded()) {
        BalancedSample sample;
        if (o.balanced_stage == "query") {
            sample = read_balanced_sample(balanced_sample_path(paths), graph, parameters);
            log << "BALANCED_SAMPLE source=selected queries=500 empty=250 nonempty=250 query_hash="
                << workload::workload_hash(sample.queries) << '\n';
        }
        baseline = saved_baseline(o, dataset, graph, parameters, log,
            paths.directory / "query_data" / (paths.prefix + "_queries.tsv"), sample.queries);
        if (o.balanced_stage == "query") {
            require(baseline.groups[1].count == 250 && baseline.groups[2].count == 250,
                    "Online does not confirm the selected 250/250 sample");
            for (std::size_t i = 0; i < sample.classes.size(); ++i)
                require((baseline.answer_counts[i].u_count > 0) == (sample.classes[i] == '1'),
                        "Online classification differs from the screening index");
            log << "BALANCED_VERIFY online_empty=250 online_nonempty=250 labels_match=500 status=ok\n";
        }
    }
    if (vary_parameters) {
        for (const auto& setting : settings) {
            const auto& p = setting.parameters;
            if (p.excluded() || baselines.count(p)) continue;
            std::ostringstream name;
            name << "case_" << std::setw(3) << std::setfill('0') << baselines.size() + 1;
            log << "BASELINE_CONFIG case=" << name.str() << " family=" << setting.family
                << " percent=" << setting.percent << " alpha=" << p.alpha << " beta=" << p.beta
                << " window_points=" << p.window_points << " u_attributes=" << p.u_attributes
                << " v_attributes=" << p.v_attributes << '\n';
            const auto queries = workload::generate_parameter_queries(graph, p);
            std::cout << "PROGRESS exp=exp2 case=" << name.str() << " family=" << setting.family
                      << " percent=" << setting.percent << " method=Online status=prepare\n" << std::flush;
            baselines.emplace(p, saved_baseline(o, dataset, graph, p, log,
                paths.directory / "query_data" / (paths.prefix + '_' + name.str() + "_queries.tsv"), queries));
        }
        log << "PARAMETER_BASELINES unique_configurations=" << baselines.size()
            << " status=ready\n";
    }
    if (o.method == "tab") {
        require(!vary_parameters, "TAB parameter experiments are not connected");
        if (parameters.excluded()) {
            log << "RESULT exp=exp6 status=skipped_11\n";
            return false;
        }
        TABIndex index;
        const auto start = Clock::now();
        load_query_index(o, dataset, graph, index, log);
        const auto load_s = seconds(start);
        log << "INDEX source=loaded built=0 load_s=" << load_s
            << " index_mb=" << space_mb(tab_index_structure_bytes(index)) << '\n';
        default_queries(o, graph, baseline, index, parameters, log);
        return true;
    }
    if (o.method != "mat") {
        if (!vary_parameters && parameters.excluded()) {
            log << "RESULT exp=" << o.experiment << " status=skipped_11\n";
            return false;
        }
        const auto query_joint = [&](auto& index) {
            const auto start = Clock::now();
            load_query_index(o, dataset, graph, index, log);
            const double load_s = seconds(start);
            log << "INDEX source=loaded built=0 load_s=" << load_s
                << " index_mb=" << space_mb(joint::structure_bytes(index))
                << " required_build_method=opt build_method=not_encoded scope=exact_joint_attributes attribute_compression=0";
            if constexpr (std::is_same_v<std::decay_t<decltype(index)>, joint::VAT>)
                log << " query_directory=" << index.directory_ready
                    << " directory_ref_bytes=" << sizeof(joint::VATDirectoryRef)
                    << " directory_mb=" << space_mb(joint::vat_directory_bytes(index));
            log << '\n';
            if (vary_parameters)
                parameter_experiment(o, graph, index, settings, baselines, parameters, paths, log);
            else
                default_queries(o, graph, baseline, index, parameters, log);
        };
        if (o.method == "vat") { joint::VAT index; query_joint(index); }
        else { joint::UAT index; query_joint(index); }
        return true;
    }
    MATIndex index;
    {
        const auto start = Clock::now();
        load_query_index(o, dataset, graph, index, log);
        log << "INDEX source=loaded built=0 load_s=" << seconds(start)
            << " index_mb=" << space_mb(formal_index_space(index).index_bytes()) << '\n';
        if (index_uses_grouped_recovery("mat")) {
            const auto prep_start = Clock::now();
            const auto& aux = graph.grouped_v_time_adjacency();
            log << "QUERY_AUX prepare_s=" << seconds(prep_start)
                << " aux_mb=" << space_mb(grouped_aux_structural_bytes(aux))
                << " index_and_aux_mb=" << space_mb(formal_index_space(index).index_bytes() +
                                                   grouped_aux_structural_bytes(aux))
                << " scratch_peak_mb=" << space_mb(aux.scratch_peak_bytes) << '\n';
        }
    }
    if (vary_parameters) {
        parameter_experiment(o, graph, index, settings, baselines, parameters, paths, log);
    } else if (parameters.excluded()) {
        log << "RESULT exp=" << o.experiment << " status=skipped_11\n";
        return false;
    } else {
        default_queries(o, graph, baseline, index, parameters, log);
    }
    return true;
}

void check_inputs(const Options& o) {
    for (const auto& name : {o.u_label_file, std::string("graph_vlabel.txt"), o.time_file})
        require(fs::is_regular_file(o.data / name), "missing input: " + (o.data / name).string());
    require(o.index.empty() || fs::is_regular_file(o.index), "index file does not exist");
    if (!o.index.empty()) {
        const auto name = o.index.filename().string();
        const auto suffix = query_index_suffix(o.method);
        require(name.size() >= suffix.size() &&
                name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0,
                (o.method == "tab" ? "TAB query requires an index named *" :
                    (o.experiment == "exp4" ? "exp4 requires an opt " : "exp1/exp2 require an opt ") + o.method + " index named *") + suffix);
    }
    require(o.queries.empty() || fs::is_regular_file(o.queries), "query file does not exist");
}

void open_run_log(const Options& o, const DatasetEntry& dataset, const MemoryGuardInfo& guard,
                  OutputPaths& paths, std::ofstream& log) {
    for (const char* name : {"TAC_QUERY_PROFILE", "TAC_PROFILE_DYNAMIC", "TAC_PROFILE_MAT_STREAMING"})
        require(setenv(name, "0", 1) == 0, "cannot disable profiling");
    if (paths.directory.empty()) paths = new_run(o, dataset);
    std::cout << "RESULT_DIRECTORY=" << paths.directory.string()
              << "\nRESULT_PREFIX=" << paths.prefix
              << "\nRESULT_LOG=" << paths.log_file().string() << '\n' << std::flush;
    if (!log.is_open()) log.open(paths.log_file(), std::ios::app);
    require(bool(log), "cannot open reserved run log");
    if (o.experiment == "exp0" && o.task == "stats") log << "\n\n";
    const bool query_index = (o.experiment == "exp1" && o.balanced_stage != "prepare") || o.experiment == "exp2" || o.query_only ||
                             (o.experiment == "exp6" && o.task == "query");
    log << std::fixed << std::setprecision(9);
    log << std::unitbuf << "TACC experiment\nGenerated at: " << paths.generated_at
        << "\nCONFIG exp=" << o.experiment << " task=" << o.task
        << " dataset=" << dataset.name << " data=" << std::quoted(o.data.string())
        << " index=" << std::quoted(o.index.string()) << " queries=" << std::quoted(o.queries.string())
        << " profile=0 rounds=1 warmup=0 probe_queries=0"
        << " index_implementation=" << (o.experiment == "exp0" ? "none" :
            o.method == "tab" ? "tab_exhaustive" :
            o.method == "mat" ? "production_mat_uv" : "production_joint_" + o.method) << '\n';
    if (!o.balanced_stage.empty()) log << "WORKLOAD model=balanced_5050_v1 stage=" << o.balanced_stage << '\n';
    if (query_index)
        log << "TIMING model=" << batch_query_timing_model
            << " clock=steady_wall batch_queries=500 clock_reads=2 per_query_timers=disabled"
            << " attribute_encoding=excluded internal_cleanup=included returned_answer_cleanup=excluded"
            << " answer_generation=included verification=excluded group_timing=not_measured\n";
    else if (o.experiment == "exp0" && o.task == "online")
        log << "TIMING model=" << query_timing_model
            << " clock=steady_wall attribute_encoding=excluded post_answer_cleanup=excluded"
            << " answer_generation=included verification=excluded\n";
    if (o.experiment == "exp5") {
        log << "INPUT axis=" << o.task << " operation=" << (o.query_only ? "query" : "build");
        if (o.task == "L") log << " requested_u_labels=" << o.u_label_count;
        else log << " time_granularity=" << o.time_granularity;
        log << " u_label_file=" << std::quoted((o.data / o.u_label_file).string())
            << " v_label_file=" << std::quoted((o.data / "graph_vlabel.txt").string())
            << " edge_file=" << std::quoted((o.data / o.time_file).string()) << '\n';
    }
    if (o.experiment == "exp4")
        log << "# Maintenance/rebuild means are in SUMMARY and MAINTENANCE_TABLE, in milliseconds per event; validation excluded.\n";
    else if (o.experiment == "exp2")
        log << "# Parameter results are in QUERY_LOG files; avg_total_us is in microseconds.\n";
    else if (o.experiment == "exp3" || (o.experiment == "exp5" && !o.query_only) ||
             (o.experiment == "exp6" && o.task == "build"))
        log << "# Main results are the RESULT lines below; build_s is in seconds; index_mb is structural MiB.\n";
    else
        log << "# Main results are the RESULT lines below; avg_total_us is query time in microseconds.\n";
    log << "PROGRESS phase=run status=start wall_time=" << paths.generated_at << '\n';
    write_memory_guard(log, guard);
}

int finish_run(const Options& o, const OutputPaths& paths, std::ofstream& log, bool completed = true) {
    if (!completed) {
        log << "PROGRESS phase=run status=finished wall_time=" << wall_time() << " exit_code=0\n";
        log.close();
        require(bool(log), "output write failed");
        std::cout << "RESULT exp=" << o.experiment << " status=skipped_11\n";
        return 0;
    }
    log << "RESULT exp=" << o.experiment << " status=ok process_peak_mb=" << space_mb(peak_bytes()) << '\n';
    log << "PROGRESS phase=run status=finished wall_time=" << wall_time() << " exit_code=0\n";
    log.close();
    require(bool(log), "output write failed");
    std::cout << "RESULT exp=" << o.experiment << " status=ok log="
              << std::quoted(paths.log_file().string()) << '\n';
    return 0;
}

int report_failure(const std::exception& error, const OutputPaths& paths, std::ofstream& log) {
    if (log.is_open()) {
        log << "RESULT status=failed reason=" << std::quoted(error.what()) << '\n';
        log << "PROGRESS phase=run status=failed wall_time=" << wall_time() << " exit_code=2\n";
    }
    std::cerr << "RESULT status=failed reason=" << std::quoted(error.what());
    if (!paths.directory.empty())
        std::cerr << " result_directory=" << std::quoted(paths.directory.string())
                  << " result_prefix=" << paths.prefix;
    std::cerr << '\n';
    return 2;
}

int not_implemented(int argc, char** argv) {
    try {
        const auto o = parse(argc, argv);
        (void)find_dataset(o.selector);
        std::cerr << "RESULT exp=" << o.experiment << " status=not_implemented\n";
        return 3;
    } catch (const std::exception& error) {
        std::cerr << "RESULT status=failed reason=" << std::quoted(error.what()) << '\n';
        return 2;
    }
}

}

int run_exp0(int argc, char** argv) {
    OutputPaths paths;
    std::ofstream log;
    try {
        auto o = parse(argc, argv);
        const auto& dataset = find_dataset(o.selector);
        check_inputs(o);

        if (!o.memory_percent.empty())
            require(setenv("TAC_MEMORY_GUARD_PERCENT", o.memory_percent.c_str(), 1) == 0, "cannot set memory guard");
        const auto guard = configure_memory_guard();
        require(memory_guard_ready(guard), "memory guard could not be activated");

        if (o.task == "stats") {
            paths = new_run(o, dataset);
            std::ostringstream graph_log;
            auto graph = load_graph(o.data, graph_log);
            const auto report = collect_dataset_stats(graph, o.data);
            log.open(paths.log_file(), std::ios::app);
            require(bool(log), "cannot open reserved run log");
            write_dataset_stats_fields(log, report);
            open_run_log(o, dataset, guard, paths, log);
            log << graph_log.str();
            write_dataset_stats(log, report);
        } else {
            open_run_log(o, dataset, guard, paths, log);
            auto graph = load_graph(o.data, log);
            auto query = o.single_query;
            query.u_attributes = parse_attributes(graph, o.u_names, true);
            query.v_attributes = parse_attributes(graph, o.v_names, false);
            QueryTimer timer;
            auto core = online_query(graph, query, nullptr, &timer);
            timer.stop();
            const double online_us = timer.microseconds();
            std::sort(core.u_vertices.begin(), core.u_vertices.end());
            std::sort(core.v_vertices.begin(), core.v_vertices.end());
            for (auto* output : {static_cast<std::ostream*>(&std::cout), static_cast<std::ostream*>(&log)}) {
                *output << "RESULT exp=exp0 task=online alpha=" << query.alpha << " beta=" << query.beta
                        << " U_attrs=" << o.u_names << " V_attrs=" << o.v_names
                        << " ts=" << query.window.start << " te=" << query.window.end
                        << " queries=1 online_us=" << online_us
                        << " u_vertices=" << core.u_vertices.size() << " v_vertices=" << core.v_vertices.size()
                        << " status=ok\n";
                for (const auto& side : {std::make_pair("U", &core.u_vertices), std::make_pair("V", &core.v_vertices)}) {
                    *output << side.first << " = {";
                    for (std::size_t i = 0; i < side.second->size(); ++i) {
                        if (i) *output << ',';
                        *output << (*side.second)[i];
                    }
                    *output << "}\n";
                }
            }
        }
        return finish_run(o, paths, log);
    } catch (const std::exception& error) {
        return report_failure(error, paths, log);
    }
}

int run_exp1(int argc, char** argv) {
    OutputPaths paths;
    std::ofstream log;
    try {
        auto o = parse(argc, argv);
        const auto& dataset = find_dataset(o.selector);
        check_inputs(o);
        if (!o.memory_percent.empty())
            require(setenv("TAC_MEMORY_GUARD_PERCENT", o.memory_percent.c_str(), 1) == 0, "cannot set memory guard");
        const auto guard = configure_memory_guard();
        require(memory_guard_ready(guard), "memory guard could not be activated");
        if (o.method == "mat") {
            (void)effective_query_route(configured_index_scope());
            (void)index_uses_grouped_recovery("mat");
        }
        open_run_log(o, dataset, guard, paths, log);
        auto graph = load_graph(o.data, log);
        if (o.balanced_stage == "prepare") {
            prepare_balanced_sample(o, dataset, graph, paths, log);
            return finish_run(o, paths, log);
        }
        const bool completed = query_experiment(o, dataset, graph, paths, log,  false);
        return finish_run(o, paths, log, completed);
    } catch (const std::exception& error) {
        return report_failure(error, paths, log);
    }
}

int run_exp2(int argc, char** argv) {
    OutputPaths paths;
    std::ofstream log;
    try {
        auto o = parse(argc, argv);
        const auto& dataset = find_dataset(o.selector);
        check_inputs(o);
        if (!o.memory_percent.empty())
            require(setenv("TAC_MEMORY_GUARD_PERCENT", o.memory_percent.c_str(), 1) == 0, "cannot set memory guard");
        const auto guard = configure_memory_guard();
        require(memory_guard_ready(guard), "memory guard could not be activated");
        if (o.method == "mat") {
            (void)effective_query_route(configured_index_scope());
            (void)index_uses_grouped_recovery("mat");
        }
        open_run_log(o, dataset, guard, paths, log);
        auto graph = load_graph(o.data, log);
        const bool completed = query_experiment(o, dataset, graph, paths, log,  true);
        return finish_run(o, paths, log, completed);
    } catch (const std::exception& error) {
        return report_failure(error, paths, log);
    }
}

int run_exp3(int argc, char** argv) {
    OutputPaths paths;
    std::ofstream log;
    try {
        auto o = parse(argc, argv);
        const auto& dataset = find_dataset(o.selector);
        check_inputs(o);
        if (!o.memory_percent.empty())
            require(setenv("TAC_MEMORY_GUARD_PERCENT", o.memory_percent.c_str(), 1) == 0, "cannot set memory guard");
        const auto guard = configure_memory_guard();
        require(memory_guard_ready(guard), "memory guard could not be activated");
        if (o.method == "mat") {
            (void)effective_query_route(configured_index_scope());
            (void)index_uses_grouped_recovery("mat");
        }
        open_run_log(o, dataset, guard, paths, log);
        auto graph = load_graph(o.data, log);
        if (o.method == "mat") construction(o, graph, paths, log);
        else joint_construction(o, graph, paths, log);
        return finish_run(o, paths, log);
    } catch (const std::exception& error) {
        return report_failure(error, paths, log);
    }
}

int run_exp4(int argc, char** argv) {
    OutputPaths paths;
    std::ofstream log;
    try {
        auto o = parse(argc, argv);
        const auto& dataset = find_dataset(o.selector);
        check_inputs(o);
        auto settings = configure_mat_maintenance(o.updates, o.update_seed);
        if (!o.memory_percent.empty())
            require(setenv("TAC_MEMORY_GUARD_PERCENT", o.memory_percent.c_str(), 1) == 0, "cannot set memory guard");
        const auto guard = configure_memory_guard();
        require(memory_guard_ready(guard), "memory guard could not be activated");
        open_run_log(o, dataset, guard, paths, log);
        const auto graph = load_graph(o.data, log);
        require(!graph.time_points().empty() && graph.num_attributes() <= 64, "unsupported maintenance graph");
        require(graph.num_edges() == graph.num_unique_temporal_edges(), "maintenance requires deduplicated (u,v,t) input");
        require(o.updates < graph.num_unique_temporal_edges(), "update count must leave at least one temporal event");

        const auto start = Clock::now();
        MATIndex initial;
        const auto directory = o.output / dataset.name / "exp3_construction/indexes";
        const bool reuse = !o.index.empty() || !newest_files(directory, "_mat_opt_uv_index.bin").empty();
        if (reuse) load_query_index(o, dataset, graph, initial, log);
        else initial = build_mat_index(graph, BuildMethod::Opt);
        log << "INITIAL source=" << (reuse ? "loaded" : "built") << " prepare_s=" << seconds(start)
            << " fingerprint=" << graph_fingerprint(graph) << " unique_events=" << graph.num_unique_temporal_edges() << '\n';
        settings.data = o.data;
        settings.output = paths.directory;
        settings.prefix = paths.prefix;
        run_mat_maintenance(graph, std::move(initial), settings, log);
        return finish_run(o, paths, log);
    } catch (const std::exception& error) {
        return report_failure(error, paths, log);
    }
}

int run_exp5(int argc, char** argv) {
    OutputPaths paths;
    std::ofstream log;
    try {
        auto o = parse(argc, argv);
        const auto& dataset = find_dataset(o.selector);
        check_inputs(o);
        if (!o.memory_percent.empty())
            require(setenv("TAC_MEMORY_GUARD_PERCENT", o.memory_percent.c_str(), 1) == 0, "cannot set memory guard");
        const auto guard = configure_memory_guard();
        require(memory_guard_ready(guard), "memory guard could not be activated");
        (void)effective_query_route(configured_index_scope());
        if (o.method == "mat") (void)index_uses_grouped_recovery("mat");
        open_run_log(o, dataset, guard, paths, log);
        auto graph = load_graph(o.data, log, o.u_label_file, o.time_file);
        std::set<AttributeID> u_labels, v_labels;
        for (VertexID u = 0; u < graph.num_u_vertices(); ++u) u_labels.insert(graph.u_attribute(u));
        for (VertexID v = 0; v < graph.num_v_vertices(); ++v) v_labels.insert(graph.v_attribute(v));
        if (o.task == "L")
            require(u_labels.size() == o.u_label_count,
                    "U label count mismatch: requested=" + std::to_string(o.u_label_count) +
                    " actual=" + std::to_string(u_labels.size()));
        o.u_label_count = u_labels.size();
        o.v_label_count = v_labels.size();
        log << "SCALABILITY axis=" << o.task;
        if (o.task == "T") log << " time_granularity=" << o.time_granularity;
        log << " u_labels=" << o.u_label_count << " v_labels=" << o.v_label_count
            << " u_vertices=" << graph.num_u_vertices() << " v_vertices=" << graph.num_v_vertices()
            << " temporal_events=" << graph.num_edges() << " time_points=" << graph.time_points().size() << '\n';
        if (o.query_only) {
            const bool completed = query_experiment(o, dataset, graph, paths, log, false);
            return finish_run(o, paths, log, completed);
        }
        if (o.method == "mat") construction(o, graph, paths, log);
        else joint_construction(o, graph, paths, log);
        return finish_run(o, paths, log);
    } catch (const std::exception& error) {
        return report_failure(error, paths, log);
    }
}

int run_exp6(int argc, char** argv) {
    OutputPaths paths;
    std::ofstream log;
    try {
        auto o = parse(argc, argv);
        const auto& dataset = find_dataset(o.selector);
        check_inputs(o);
        if (!o.memory_percent.empty())
            require(setenv("TAC_MEMORY_GUARD_PERCENT", o.memory_percent.c_str(), 1) == 0, "cannot set memory guard");
        const auto guard = configure_memory_guard();
        require(memory_guard_ready(guard), "memory guard could not be activated");
        open_run_log(o, dataset, guard, paths, log);
        const auto graph = load_graph(o.data, log);
        if (o.task == "build") {
            tab_construction(o, graph, paths, log);
            return finish_run(o, paths, log);
        }
        const bool completed = query_experiment(o, dataset, graph, paths, log, false);
        return finish_run(o, paths, log, completed);
    } catch (const std::exception& error) {
        return report_failure(error, paths, log);
    }
}

int run_exp7(int argc, char** argv) {
    return not_implemented(argc, argv);
}

namespace {
namespace dblp_case {
namespace fs = std::filesystem;
using Number = std::uint32_t;
using Names = std::map<Number, std::string>;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

Number number(const std::string& text, const char* field) {
    Number value{};
    const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
    require(!text.empty() && result.ec == std::errc{} && result.ptr == text.data() + text.size(),
            std::string(field) + " must be an unsigned 32-bit integer: " + text);
    return value;
}

bool is_graph_directory(const fs::path& path) {
    return fs::is_regular_file(path / "graph_ulabel.txt") &&
           fs::is_regular_file(path / "graph_vlabel.txt") &&
           fs::is_regular_file(path / "graph_timestamp.txt");
}

fs::path default_data_directory() {
    const fs::path relative = "data/case_dblp/case_study_data";
    for (auto root = fs::current_path(); !root.empty();) {
        if (is_graph_directory(root / relative)) return root / relative;
        const auto parent = root.parent_path();
        if (parent == root) break;
        root = parent;
    }
    const fs::path source = __FILE__;
    if (source.is_absolute()) {
        const auto root = source.parent_path().parent_path().parent_path();
        if (is_graph_directory(root / relative)) return root / relative;
    }
    throw std::runtime_error("cannot locate DBLP case data; specify --data-dir <directory>");
}

void usage(const char* program) {
    std::cout << "exp8 DBLP Online (same algorithm as exp0 online)\n"
              << "Usage: " << program << " <alpha> <beta> <U_attribute_IDs> <V_attribute_IDs> <ts> <te>\n"
              << "       [--data-dir DIR] [--note-dir DIR]\n"
              << "Use comma-separated attribute IDs, or '*' / all for every attribute of that side.\n"
              << "Times are inclusive publication years. alpha/beta must be positive.\n"
              << "Example: " << program << " 2 2 15 57 2010 2016\n"
              << "Example: " << program << " 2 2 4,15 42,57 2004 2017\n"
              << "Default data: data/case_dblp/case_study_data; names: DATA/note/*.txt\n"
              << "Results are printed to stdout; no project files are written.\n";
}

struct Options {
    Number alpha{}, beta{}, start{}, end{};
    std::string u_text, v_text;
    fs::path data, notes;
};

Options parse(int argc, char** argv) {
    require(argc >= 7, "expected alpha beta U_attribute_IDs V_attribute_IDs ts te; use --help");
    Options options;
    options.alpha = number(argv[1], "alpha");
    options.beta = number(argv[2], "beta");
    options.u_text = argv[3];
    options.v_text = argv[4];
    options.start = number(argv[5], "ts");
    options.end = number(argv[6], "te");
    require(options.alpha > 0 && options.beta > 0, "alpha and beta must be positive");
    require(options.start <= options.end, "time window requires ts <= te");
    bool data_seen = false, notes_seen = false;
    for (int i = 7; i < argc; ++i) {
        const std::string flag = argv[i];
        require(i + 1 < argc, "missing value for option: " + flag);
        if (flag == "--data-dir") {
            require(!data_seen, "--data-dir specified twice");
            options.data = argv[++i];
            data_seen = true;
        } else if (flag == "--note-dir") {
            require(!notes_seen, "--note-dir specified twice");
            options.notes = argv[++i];
            notes_seen = true;
        } else {
            throw std::runtime_error("unknown option: " + flag);
        }
    }
    if (options.data.empty()) options.data = default_data_directory();
    require(is_graph_directory(options.data), "missing one or more graph TXT files in: " + options.data.string());
    if (options.notes.empty()) options.notes = options.data / "note";
    return options;
}

Names read_names(const fs::path& path) {
    Names names;
    if (!fs::exists(path)) return names;
    std::ifstream input(path);
    require(input.good(), "cannot open name list: " + path.string());
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line_number == 1 && line.compare(0, 3, "\xef\xbb\xbf") == 0) line.erase(0, 3);
        line = trim(line);
        if (line.empty() || line[0] == '#' || line[0] == '%') continue;
        const auto separator = line.find_first_of(" \t");
        require(separator != std::string::npos, "invalid name list row: " + path.string() + ":" + std::to_string(line_number));
        const auto id = number(line.substr(0, separator), "name-list ID");
        const auto name = trim(line.substr(separator + 1));
        require(!name.empty(), "empty name in: " + path.string());
        require(names.emplace(id, name).second, "duplicate name-list ID in: " + path.string());
    }
    require(!input.bad(), "error reading: " + path.string());
    return names;
}

std::string name_or_id(const Names& names, Number id) {
    const auto found = names.find(id);
    return found == names.end() ? std::to_string(id) : found->second;
}

struct Vertex {
    Number external_attribute{};
    tacc::AttributeID internal_attribute{};
};
using Vertices = std::map<tacc::VertexID, Vertex>;
using AttributeDomain = std::map<Number, tacc::AttributeID>;

struct Side {
    Vertices vertices;
    AttributeDomain attributes;
};

Side read_side(const tacc::TemporalBipartiteGraph& graph, bool upper) {
    Side side;
    const auto count = upper ? graph.num_u_vertices() : graph.num_v_vertices();
    require(count <= std::numeric_limits<tacc::VertexID>::max(), "too many vertices");
    for (std::size_t i = 0; i < count; ++i) {
        const auto internal_vertex = static_cast<tacc::VertexID>(i);
        const auto original = upper ? graph.original_u_id(internal_vertex) : graph.original_v_id(internal_vertex);
        const auto internal_attr = upper ? graph.u_attribute(internal_vertex) : graph.v_attribute(internal_vertex);
        const auto external_attr = number(graph.attribute_name(internal_attr), "attribute ID in graph TXT");
        require(side.vertices.emplace(original, Vertex{external_attr, internal_attr}).second, "duplicate vertex ID");
        const auto existing = side.attributes.emplace(external_attr, internal_attr);
        require(existing.second || existing.first->second == internal_attr, "inconsistent attribute mapping");
    }
    return side;
}

struct Selection {
    std::vector<Number> external;
    std::vector<tacc::AttributeID> internal;
};

Selection select_attributes(std::string text, const AttributeDomain& domain, const char* side) {
    std::set<Number> wanted;
    text = trim(text);
    if (text == "*" || text == "all") {
        for (const auto& entry : domain) wanted.insert(entry.first);
    } else {
        for (auto& c : text) if (c == ';' || c == '|') c = ',';
        require(!text.empty() && text.back() != ',', std::string(side) + " attributes must not be empty");
        std::istringstream input(text);
        std::string token;
        while (std::getline(input, token, ',')) wanted.insert(number(trim(token), "query attribute ID"));
    }
    require(!wanted.empty(), std::string(side) + " attributes must not be empty");
    Selection selection;
    for (const auto external : wanted) {
        const auto found = domain.find(external);
        require(found != domain.end(), std::string("unknown ") + side + "-side attribute ID: " + std::to_string(external));
        selection.external.push_back(external);
        selection.internal.push_back(found->second);
    }
    return selection;
}

template<class T>
void print_ids(const char* label, const std::vector<T>& values) {
    std::cout << label << " = {";
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << values[i];
    }
    std::cout << "}\n";
}

int run(const Options& options) {
    tacc::TemporalBipartiteGraph graph;
    require(graph.load_from_files((options.data / "graph_ulabel.txt").string(),
                                  (options.data / "graph_vlabel.txt").string(),
                                  (options.data / "graph_timestamp.txt").string()), "failed to load graph");
    const auto upper = read_side(graph, true);
    const auto lower = read_side(graph, false);
    const auto u_selection = select_attributes(options.u_text, upper.attributes, "U");
    const auto v_selection = select_attributes(options.v_text, lower.attributes, "V");
    tacc::Query query{{options.start, options.end}, u_selection.internal, v_selection.internal, options.alpha, options.beta};

    tacc::QueryProfile profile;
    tacc::QueryTimer timer;
    auto answer = tacc::online_query(graph, query, &profile, &timer);
    const auto query_us = timer.microseconds();
    std::sort(answer.u_vertices.begin(), answer.u_vertices.end());
    std::sort(answer.v_vertices.begin(), answer.v_vertices.end());
    const std::set<tacc::VertexID> active_u(answer.u_vertices.begin(), answer.u_vertices.end());
    const std::set<tacc::VertexID> active_v(answer.v_vertices.begin(), answer.v_vertices.end());
    std::set<std::tuple<Number, Number, Number>> retained_events;
    std::set<std::pair<Number, Number>> retained_pairs;
    std::map<Number, std::size_t> u_degree, v_degree;
    for (const auto& event : graph.edges()) {
        if (event.timestamp < options.start || event.timestamp > options.end) continue;
        const auto u = graph.original_u_id(event.u), v = graph.original_v_id(event.v);
        if (active_u.count(u) == 0 || active_v.count(v) == 0) continue;
        retained_events.emplace(u, v, event.timestamp);
        if (retained_pairs.emplace(u, v).second) { ++u_degree[u]; ++v_degree[v]; }
    }
    for (const auto u : answer.u_vertices) require(u_degree[u] >= options.alpha, "result violates alpha");
    for (const auto v : answer.v_vertices) require(v_degree[v] >= options.beta, "result violates beta");

    const auto author_names = read_names(options.notes / "author_list.txt");
    const auto paper_names = read_names(options.notes / "paper_list.txt");
    const auto direction_names = read_names(options.notes / "research_direction_list.txt");
    const auto venue_names = read_names(options.notes / "venue_list.txt");
    std::cout << "RESULT exp=exp8 status=ok method=Online data_dir=" << std::quoted(options.data.string())
              << " alpha=" << options.alpha << " beta=" << options.beta
              << " ts=" << options.start << " te=" << options.end
              << " graph_u=" << graph.num_u_vertices() << " graph_v=" << graph.num_v_vertices()
              << " graph_events=" << graph.num_edges() << " graph_attributes=" << graph.num_attributes()
              << " candidate_u=" << profile.candidate_u << " candidate_v=" << profile.candidate_v
              << " window_edges=" << profile.window_edges
              << " result_u=" << answer.u_vertices.size() << " result_v=" << answer.v_vertices.size()
              << " result_edges=" << retained_pairs.size() << " result_events=" << retained_events.size()
              << " query_us=" << std::fixed << std::setprecision(3) << query_us << '\n';
    print_ids("U_ATTRIBUTES", u_selection.external);
    print_ids("V_ATTRIBUTES", v_selection.external);
    print_ids("U", answer.u_vertices);
    print_ids("V", answer.v_vertices);
    std::cout << "\n[U_VERTICES]\nauthor_id\tresearch_direction_id\tdegree\tauthor_name\tresearch_direction\n";
    for (const auto u : answer.u_vertices) {
        const auto attr = upper.vertices.at(u).external_attribute;
        std::cout << u << '\t' << attr << '\t' << u_degree[u] << '\t'
                  << name_or_id(author_names, u) << '\t' << name_or_id(direction_names, attr) << '\n';
    }
    std::cout << "\n[V_VERTICES]\npaper_id\tvenue_id\tdegree\tpaper_title\tvenue_name\n";
    for (const auto v : answer.v_vertices) {
        const auto attr = lower.vertices.at(v).external_attribute;
        std::cout << v << '\t' << attr << '\t' << v_degree[v] << '\t'
                  << name_or_id(paper_names, v) << '\t' << name_or_id(venue_names, attr) << '\n';
    }
    std::cout << "\n[TEMPORAL_EDGES]\nauthor_id\tpaper_id\tyear\n";
    for (const auto& event : retained_events) {
        std::cout << std::get<0>(event) << '\t' << std::get<1>(event) << '\t' << std::get<2>(event) << '\n';
    }
    require(std::cout.good(), "failed to write results");
    return 0;
}
}
}

int run_exp8(int argc, char** argv) {
    try {
        if (argc == 3 && (std::string(argv[2]) == "--help" || std::string(argv[2]) == "-h")) {
            const auto command = std::string(argv[0]) + " exp8";
            dblp_case::usage(command.c_str());
            return 0;
        }
        return dblp_case::run(dblp_case::parse(argc - 1, argv + 1));
    } catch (const std::exception& error) {
        std::cerr << "RESULT exp=exp8 status=failed reason=" << std::quoted(error.what()) << '\n';
        return 2;
    }
}

int run_exp9(int argc, char** argv) {
    return not_implemented(argc, argv);
}

}
