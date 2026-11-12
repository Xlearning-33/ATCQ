#include "experiments/experiment.h"

#include "common/graph.h"
#include "common/memory_guard.h"
#include "index/tab_index.h"
#include "index/index_space.h"
#include "query/index_query_utils.h"
#include "query/mat_recovery.h"
#include "query/online_query.h"

#include <chrono>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <new>
#include <optional>
#include <random>
#include <set>
#include <type_traits>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>
#include <sys/resource.h>

namespace tacc {

namespace {

MemoryGuardInfo& active_memory_guard() {
    static MemoryGuardInfo info;
    return info;
}

std::optional<std::vector<AttributeID>> resolve_attributes(
    const TemporalBipartiteGraph& graph,
    const std::vector<std::string>& names) {
    std::vector<AttributeID> attributes;
    attributes.reserve(names.size());
    for (const std::string& name : names) {
        const auto id = graph.attribute_id(name);
        if (!id) {
            std::cerr << "Unknown query attribute: " << name << '\n';
            return std::nullopt;
        }
        attributes.push_back(*id);
    }
    return attributes;
}

void write_vertices(std::ostream& output,
                    const char* side,
                    const std::vector<VertexID>& vertices) {
    output << side << " = {";
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        if (i != 0) output << ',';
        output << vertices[i];
    }
    output << "}\n";
}

void print_vertices(const char* side, const std::vector<VertexID>& vertices) {
    write_vertices(std::cout, side, vertices);
}

std::vector<std::string> split_attributes(std::string text) {
    for (char& ch : text) {
        if (ch == ';' || ch == '|') ch = ',';
    }

    std::vector<std::string> names;
    std::istringstream input(text);
    std::string name;
    while (std::getline(input, name, ',')) {
        if (!name.empty()) names.push_back(name);
    }
    return names;
}

std::optional<std::uint32_t> parse_uint32(const char* text) {
    try {
        std::size_t parsed = 0;
        const unsigned long long value = std::stoull(text, &parsed);
        if (text[parsed] != '\0' ||
            value > std::numeric_limits<std::uint32_t>::max()) {
            return std::nullopt;
        }
        return static_cast<std::uint32_t>(value);
    } catch (...) {
        return std::nullopt;
    }
}

void print_usage(const char* program) {
    std::cerr
        << "Usage:\n"
        << "  " << program << "\n"
        << "  " << program << " <data_dir>\n"
        << "  " << program
        << " <data_dir> <ts> <te> <alpha> <beta> <u_attrs> <v_attrs>\n\n"
        << "  " << program << " tab <data_dir>\n"
        << "  " << program
        << " tab <data_dir> <ts> <te> <alpha> <beta> <u_attrs> <v_attrs>\n\n"
        << "Attribute lists use commas, for example: a,b  x\n";
}

bool run_query_case(const TemporalBipartiteGraph& graph,
                    const char* name,
                    const std::vector<std::string>& u_names,
                    const std::vector<std::string>& v_names,
                    TimeWindow window,
                    Count alpha,
                    Count beta,
                    const QueryResult& expected) {
    const auto u_attributes = resolve_attributes(graph, u_names);
    const auto v_attributes = resolve_attributes(graph, v_names);
    if (!u_attributes || !v_attributes) return false;

    Query query{window, *u_attributes, *v_attributes, alpha, beta};
    const QueryResult result = online_query(graph, query);

    std::cout << "\n[" << name << "] window=[" << window.start << ','
              << window.end << "], alpha=" << alpha << ", beta=" << beta
              << '\n';
    print_vertices("U", result.u_vertices);
    print_vertices("V", result.v_vertices);

    const bool passed = result.u_vertices == expected.u_vertices &&
                        result.v_vertices == expected.v_vertices;
    std::cout << "check=" << (passed ? "PASS" : "FAIL") << '\n';
    return passed;
}

std::string current_filename_timestamp() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t time = std::chrono::system_clock::to_time_t(now);
    std::tm local_time{};
    localtime_r(&time, &local_time);
    std::ostringstream output;
    output << std::put_time(&local_time, "%Y%m%d_%H%M%S");
    return output.str();
}

std::string current_wall_time() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t time = std::chrono::system_clock::to_time_t(now);
    std::tm local_time{};
    localtime_r(&time, &local_time);
    std::ostringstream output;
    output << std::put_time(&local_time, "%Y-%m-%d %H:%M:%S");
    return output.str();
}

std::string result_dataset_name(const std::filesystem::path& data_dir) {
    std::string name = data_dir.lexically_normal().filename().string();
    if (name.empty() || name == "." || name == "..") name = "dataset";
    for (char& ch : name) {
        const bool safe = (ch >= 'a' && ch <= 'z') ||
                          (ch >= 'A' && ch <= 'Z') ||
                          (ch >= '0' && ch <= '9') ||
                          ch == '_' || ch == '-';
        if (!safe) ch = '_';
    }
    return name;
}

struct DatasetResultPaths {
    std::string dataset;
    std::filesystem::path root;
    std::filesystem::path check;
    std::filesystem::path indexes;
};

DatasetResultPaths dataset_result_paths(
    const std::filesystem::path& data_dir) {
    DatasetResultPaths paths;
    paths.dataset = result_dataset_name(data_dir);
    paths.root = std::filesystem::path("results") / paths.dataset;
    paths.check = paths.root / "check";
    paths.indexes = paths.root / "indexes";
    return paths;
}

bool create_result_directories(const DatasetResultPaths& paths,
                               std::error_code& error) {
    std::filesystem::create_directories(paths.root, error);
    if (!error) std::filesystem::create_directories(paths.check, error);
    if (!error) std::filesystem::create_directories(paths.indexes, error);
    return !error;
}

double bytes_to_mb(std::uintmax_t bytes) {
    return static_cast<double>(bytes) / (1024.0 * 1024.0);
}

struct MemorySnapshot {
    std::size_t rss_kb{};
    std::size_t hwm_kb{};
    std::size_t vmsize_kb{};
    std::size_t vmpeak_kb{};
    bool complete{};
};

MemorySnapshot read_memory_snapshot() {
    MemorySnapshot memory;
    unsigned found = 0;
    std::ifstream input("/proc/self/status");
    for (std::string line; std::getline(input, line);) {
        std::istringstream fields(line);
        std::string key;
        std::size_t value = 0;
        std::string unit;
        if (!(fields >> key >> value >> unit)) continue;
        if (key == "VmRSS:") { memory.rss_kb = value; found |= 1; }
        if (key == "VmHWM:") { memory.hwm_kb = value; found |= 2; }
        if (key == "VmSize:") { memory.vmsize_kb = value; found |= 4; }
        if (key == "VmPeak:") { memory.vmpeak_kb = value; found |= 8; }
    }
    memory.complete = found == 15;
    struct rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) == 0 && usage.ru_maxrss > 0) {
        memory.hwm_kb = std::max(memory.hwm_kb,
                              static_cast<std::size_t>(usage.ru_maxrss));
    }
    return memory;
}

double kilobytes_to_mb(std::size_t kilobytes) {
    return static_cast<double>(kilobytes) / 1024.0;
}

void write_bad_alloc_memory(std::ostream& output, const char* stage) {
    MemorySnapshot memory;
    try {
        memory = read_memory_snapshot();
    } catch (const std::bad_alloc&) {
    }
    output << "RESULT exp=memory_at_failure"
           << " stage=" << stage
           << " snapshot_status=" << (memory.complete ? "ok" : "unavailable")
           << " rss_mb=" << kilobytes_to_mb(memory.rss_kb)
           << " hwm_mb=" << kilobytes_to_mb(memory.hwm_kb)
           << " vmsize_mb=" << kilobytes_to_mb(memory.vmsize_kb)
           << " vmpeak_mb=" << kilobytes_to_mb(memory.vmpeak_kb)
           << " guard_limit_bytes=" << active_memory_guard().effective_bytes
           << '\n';
}

int finish_log(const std::filesystem::path& log_path, int exit_code) {
    {
        std::ofstream output(log_path, std::ios::app);
        output << "PROGRESS phase=run status=finished"
               << " wall_time=" << current_wall_time()
               << " exit_code=" << exit_code << '\n';
    }
    std::ifstream input(log_path);
    std::cout << input.rdbuf();
    std::cout << "\nLog saved to " << log_path.string() << '\n';
    return exit_code;
}

int run_tab_experiment(int argc,
                       char** argv,
                       const std::filesystem::path& data_dir,
                       const TemporalBipartiteGraph& graph) {
    const std::string stamp = current_filename_timestamp();
    const DatasetResultPaths result_paths = dataset_result_paths(data_dir);
    std::error_code directory_error;
    if (!create_result_directories(result_paths, directory_error)) {
        std::cerr << "Cannot create result directory: "
                  << result_paths.root.string() << ": "
                  << directory_error.message() << '\n';
        return 1;
    }

    const std::filesystem::path log_path =
        result_paths.root /
        (stamp + "_" + result_paths.dataset + "_tab_index.log");
    const std::filesystem::path index_path =
        result_paths.indexes /
        (stamp + "_" + result_paths.dataset + "_tab_index.bin");
    const std::filesystem::path text_index_path =
        result_paths.indexes /
        (stamp + "_" + result_paths.dataset + "_tab_index_text.log");
    const bool write_text_index =
        data_dir.lexically_normal().filename() == "small";
    std::ofstream log(log_path);
    if (!log) {
        std::cerr << "Cannot create log file: " << log_path.string() << '\n';
        return 1;
    }
    log << std::unitbuf << std::fixed << std::setprecision(6);
    log << "# TACC TAB-Index experiment\n"
        << "Generated at: " << current_wall_time() << '\n'
        << "RESULT exp=tab_dataset"
        << " dataset=" << data_dir.string()
        << " u_vertices=" << graph.num_u_vertices()
        << " v_vertices=" << graph.num_v_vertices()
        << " temporal_edges=" << graph.num_edges()
        << " unique_edges=" << graph.num_unique_edges()
        << " graph_query_aux_bytes=" << graph.temporal_query_aux_bytes()
        << " timestamps=" << graph.time_points().size()
        << " min_timestamp=" << graph.min_timestamp()
        << " max_timestamp=" << graph.max_timestamp() << '\n'
        << "PROGRESS phase=tab_index status=build_start"
        << " wall_time=" << current_wall_time() << '\n';

    TABIndex index;
    try {
        index = build_tab_index_exhaustive(graph);
    } catch (const std::exception& error) {
        log << "RESULT exp=tab_index status=failed reason="
            << error.what() << '\n';
        log.close();
        return finish_log(log_path, 1);
    }
    log << "PROGRESS phase=tab_index status=build_done"
        << " wall_time=" << current_wall_time() << '\n';

    bool text_dump_ok = true;
    double text_dump_seconds = 0.0;
    std::uintmax_t text_bytes = 0;
    std::error_code text_size_error;
    if (write_text_index) {
        std::string text_error;
        const auto text_start = std::chrono::steady_clock::now();
        text_dump_ok = save_tab_index_text(
            text_index_path.string(), graph, index, &text_error);
        text_dump_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - text_start).count();
        if (text_dump_ok) {
            text_bytes = std::filesystem::file_size(text_index_path,
                                                    text_size_error);
        }
        log << "RESULT exp=tab_index_text_dump"
            << " dataset_type=small"
            << " status=" << (text_dump_ok ? "ok" : "failed")
            << " dump_s=" << text_dump_seconds
            << " text_bytes=" << (text_size_error ? 0 : text_bytes)
            << " text_mb="
            << (text_size_error ? 0.0 : bytes_to_mb(text_bytes))
            << " path=" << text_index_path.string();
        if (!text_dump_ok) log << " reason=" << text_error;
        log << '\n';
    } else {
        log << "RESULT exp=tab_index_text_dump"
            << " dataset_type=real"
            << " status=skipped"
            << " reason=binary_only\n";
    }

    const std::uint64_t hash_before = tab_index_hash(index);
    std::string dump_error;
    const auto dump_start = std::chrono::steady_clock::now();
    const bool dump_ok =
        save_tab_index_binary(index_path.string(), index, &dump_error);
    const double dump_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - dump_start).count();

    log << "RESULT exp=tab_index_dump"
        << " format=binary"
        << " status=" << (dump_ok ? "ok" : "failed")
        << " dump_s=" << dump_seconds
        << " path=" << index_path.string();
    if (!dump_ok) log << " reason=" << dump_error;
    log << '\n';

    TABIndex loaded;
    std::string load_error;
    const auto load_start = std::chrono::steady_clock::now();
    const bool load_ok = dump_ok &&
        load_tab_index_binary(index_path.string(), loaded, &load_error);
    const double load_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - load_start).count();
    const std::uint64_t hash_after = load_ok ? tab_index_hash(loaded) : 0;
    const bool load_same = load_ok && hash_before == hash_after &&
        index.records.size() == loaded.records.size() &&
        index.stats.vertex_references == loaded.stats.vertex_references;

    log << "RESULT exp=tab_index_load"
        << " status=" << (load_ok ? "ok" : "failed")
        << " load_s=" << load_seconds
        << " same=" << (load_same ? 1 : 0)
        << " records=" << (load_ok ? loaded.records.size() : 0)
        << " hash_before=" << hash_before
        << " hash_after=" << hash_after;
    if (!load_ok) log << " reason=" << load_error;
    log << '\n';

    bool query_pass = true;
    if (argc == 9 && load_same) {
        const auto start = parse_uint32(argv[3]);
        const auto end = parse_uint32(argv[4]);
        const auto alpha = parse_uint32(argv[5]);
        const auto beta = parse_uint32(argv[6]);
        const auto u_attributes =
            resolve_attributes(graph, split_attributes(argv[7]));
        const auto v_attributes =
            resolve_attributes(graph, split_attributes(argv[8]));
        if (!start || !end || !alpha || !beta || *start > *end ||
            *alpha == 0 || *beta == 0 || !u_attributes || !v_attributes ||
            u_attributes->empty() || v_attributes->empty()) {
            log << "RESULT exp=tab_index_query status=failed"
                << " reason=invalid_query\n";
            query_pass = false;
        } else {
            const Query query{{*start, *end}, *u_attributes, *v_attributes,
                              *alpha, *beta};
            const auto query_start = std::chrono::steady_clock::now();
            const QueryResult tab_result = query_tab_index(loaded, graph, query);
            const double query_us = std::chrono::duration<double, std::micro>(
                std::chrono::steady_clock::now() - query_start).count();
            const QueryResult online_result = online_query(graph, query);
            query_pass =
                tab_result.u_vertices == online_result.u_vertices &&
                tab_result.v_vertices == online_result.v_vertices;

            log << "RESULT exp=tab_index_query"
                << " start=" << *start
                << " end=" << *end
                << " alpha=" << *alpha
                << " beta=" << *beta
                << " source="
                << ((*alpha == 1 && *beta == 1) ? "direct_11" : "tab")
                << " query_us=" << query_us
                << " u_result=" << tab_result.u_vertices.size()
                << " v_result=" << tab_result.v_vertices.size()
                << " online_same=" << (query_pass ? 1 : 0)
                << " status=" << (query_pass ? "ok" : "failed") << '\n';
            write_vertices(log, "U", tab_result.u_vertices);
            write_vertices(log, "V", tab_result.v_vertices);
        }
    } else if (argc == 9) {
        query_pass = false;
        log << "RESULT exp=tab_index_query status=failed"
            << " reason=index_load_validation_failed\n";
    }

    const double nonempty_rate = index.stats.enumerated_queries == 0
        ? 0.0
        : static_cast<double>(index.stats.nonempty_records) /
              static_cast<double>(index.stats.enumerated_queries);
    const bool success = text_dump_ok && dump_ok && load_same && query_pass;
    log << "RESULT exp=tab_index"
        << " method=tab-exhaustive"
        << " attribute_pairs=" << index.stats.attribute_pairs
        << " windows=" << index.stats.windows
        << " max_alpha=" << index.max_alpha
        << " max_beta=" << index.max_beta
        << " omitted_11_queries="
        << index.stats.attribute_pairs * index.stats.windows
        << " enumerated_queries=" << index.stats.enumerated_queries
        << " records=" << index.stats.nonempty_records
        << " vertex_refs=" << index.stats.vertex_references
        << " nonempty_rate=" << nonempty_rate
        << " build_s=" << index.stats.build_seconds
        << " approx_memory_bytes=" << index.stats.approx_bytes
        << " approx_memory_mb=" << bytes_to_mb(index.stats.approx_bytes)
        << " text_dump_status="
        << (write_text_index ? (text_dump_ok ? "ok" : "failed") : "skipped")
        << " text_bytes="
        << (write_text_index && !text_size_error ? text_bytes : 0)
        << " text_mb="
        << (write_text_index && !text_size_error ? bytes_to_mb(text_bytes) : 0.0)
        << " dump_s=" << dump_seconds
        << " load_s=" << load_seconds
        << " dump_status=" << (dump_ok ? "ok" : "failed")
        << " load_status=" << (load_ok ? "ok" : "failed")
        << " load_same=" << (load_same ? 1 : 0)
        << " status=" << (success ? "ok" : "failed")
        << " index_path=" << index_path.string() << '\n';
    log.close();
    return finish_log(log_path, success ? 0 : 2);
}

}

int run_experiment(int argc, char** argv) {
    if (argc >= 2 && (std::string(argv[1]) == "quota-online" ||
                      std::string(argv[1]) == "quota-query"))
        return run_quota_experiment(argc, argv);
    if (argc >= 2 &&
        (std::string(argv[1]) == "index" ||
         std::string(argv[1]) == "index-query")) {
        try {
            return run_index_experiment(argc, argv);
        } catch (const std::bad_alloc&) {
            std::cerr << "RESULT exp=memory_guard status=failed"
                      << " reason=bad_alloc stage=outside_builder\n";
            write_bad_alloc_memory(std::cerr, "outside_builder");
            return 2;
        }
    }
    const bool tab_mode = argc >= 2 && std::string(argv[1]) == "tab";
    const bool valid_online_args = argc == 1 || argc == 2 || argc == 8;
    const bool valid_tab_args = tab_mode && (argc == 3 || argc == 9);
    if ((!tab_mode && !valid_online_args) || (tab_mode && !valid_tab_args)) {
        print_usage(argv[0]);
        return 1;
    }

    const std::filesystem::path data_dir =
        tab_mode ? std::filesystem::path(argv[2])
                 : (argc >= 2 ? std::filesystem::path(argv[1])
                              : std::filesystem::path("data/small"));

    TemporalBipartiteGraph graph;
    if (!graph.load_from_files((data_dir / "graph_ulabel.txt").string(),
                               (data_dir / "graph_vlabel.txt").string(),
                               (data_dir / "graph_timestamp.txt").string())) {
        return 1;
    }

    if (tab_mode) {
        return run_tab_experiment(argc, argv, data_dir, graph);
    }

    std::cout << "dataset=" << data_dir.string()
              << " U=" << graph.num_u_vertices()
              << " V=" << graph.num_v_vertices()
              << " temporal_edges=" << graph.num_edges()
              << " unique_edges=" << graph.num_unique_edges()
              << " graph_query_aux_bytes="
              << graph.temporal_query_aux_bytes()
              << " time=[" << graph.min_timestamp() << ','
              << graph.max_timestamp() << "]\n";

    if (argc == 8) {
        const auto start = parse_uint32(argv[2]);
        const auto end = parse_uint32(argv[3]);
        const auto alpha = parse_uint32(argv[4]);
        const auto beta = parse_uint32(argv[5]);
        if (!start || !end || !alpha || !beta || *start > *end ||
            *alpha == 0 || *beta == 0) {
            std::cerr << "Invalid time window or (alpha, beta) parameters.\n";
            print_usage(argv[0]);
            return 1;
        }

        const std::vector<std::string> u_names = split_attributes(argv[6]);
        const std::vector<std::string> v_names = split_attributes(argv[7]);
        const auto u_attributes = resolve_attributes(graph, u_names);
        const auto v_attributes = resolve_attributes(graph, v_names);
        if (!u_attributes || !v_attributes || u_attributes->empty() ||
            v_attributes->empty()) {
            return 1;
        }

        const Query query{{*start, *end}, *u_attributes, *v_attributes,
                          *alpha, *beta};
        const QueryResult result = online_query(graph, query);
        std::cout << "\n[online query] window=[" << *start << ',' << *end
                  << "], alpha=" << *alpha << ", beta=" << *beta << '\n';
        print_vertices("U", result.u_vertices);
        print_vertices("V", result.v_vertices);
        return 0;
    }

    const bool empty_case = run_query_case(
        graph, "attribute constrained empty case", {"a", "b"}, {"x"},
        {1, 3}, 3, 3, {});

    const bool nonempty_case = run_query_case(
        graph, "non-empty correctness case", {"a", "b", "c"},
        {"x", "y", "z"}, {1, 3}, 2, 2, {{1, 2}, {1, 2}});

    return empty_case && nonempty_case ? 0 : 2;
}

}

#include "common/graph.h"
#include "index/index_builder.h"
#include "index/index_io.h"
#include "query/mat_query.h"
#include "query/online_query.h"
#include "query/uat_query.h"
#include "query/vat_query.h"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace tacc {
namespace {

void formal_usage(const char* program) {
    std::cerr
        << "Usage:\n  " << program
        << " index <vat|uat|mat> <data_dir>"
           " [ts te alpha beta u_attrs v_attrs]\n  "
        << program
        << " index-query <vat|uat|mat> <data_dir> <index_file>"
           " <ts> <te> <alpha> <beta> <u_attrs> <v_attrs>\n"
        << "Build method: TAC_INDEX_BUILD_METHOD=basic|opt "
           "(default: opt)\n"
        << "Build scope: UV only (TAC_INDEX_SCOPE=uv or unset)\n"
        << "Query route: UV only (TAC_QUERY_ROUTE=auto|uv or unset)\n"
        << "500-query verification: TAC_INDEX_500CHECK=1\n"
        << "Quota batch query diagnostics: TAC_QUERY_PROFILE=0|1 (default: 0; 1 is diagnostic timing)\n"
        << "Seed recovery: TAC_SEED_RECOVERY=adaptive|seed (default: adaptive)\n"
        << "UAT recovery: TAC_UAT_RECOVERY=adaptive|grouped (default: adaptive)\n"
        << "MAT UV recovery: TAC_MAT_RECOVERY=grouped|legacy (default: grouped)\n"
        << "Build progress: TAC_BUILD_PROGRESS=0|1 (default: 1)\n";
}

std::optional<Count> formal_number(const char* text) {
    try {
        std::size_t used = 0;
        const auto value = std::stoull(text, &used);
        if (text[used] != '\0' || value > UINT32_MAX) return std::nullopt;
        return static_cast<Count>(value);
    } catch (...) {
        return std::nullopt;
    }
}

std::vector<std::string> formal_split(std::string text) {
    std::vector<std::string> output;
    std::istringstream input(text);
    for (std::string item; std::getline(input, item, ',');) {
        if (!item.empty()) output.push_back(item);
    }
    return output;
}

std::optional<std::vector<AttributeID>> formal_attributes(
    const TemporalBipartiteGraph& graph, const char* text) {
    std::vector<AttributeID> output;
    for (const std::string& name : formal_split(text)) {
        const auto id = graph.attribute_id(name);
        if (!id) return std::nullopt;
        output.push_back(*id);
    }
    return output;
}

std::optional<Query> parse_formal_query(
    const TemporalBipartiteGraph& graph, char** argv, int offset) {
    const auto ts = formal_number(argv[offset]);
    const auto te = formal_number(argv[offset + 1]);
    const auto alpha = formal_number(argv[offset + 2]);
    const auto beta = formal_number(argv[offset + 3]);
    const auto u_attrs = formal_attributes(graph, argv[offset + 4]);
    const auto v_attrs = formal_attributes(graph, argv[offset + 5]);
    if (!ts || !te || !alpha || !beta || !u_attrs || !v_attrs ||
        *ts > *te || *alpha == 0 || *beta == 0 ||
        u_attrs->empty() || v_attrs->empty()) {
        return std::nullopt;
    }
    return Query{{*ts, *te}, *u_attrs, *v_attrs, *alpha, *beta};
}

bool same_result(QueryResult lhs, QueryResult rhs) {
    std::sort(lhs.u_vertices.begin(), lhs.u_vertices.end());
    std::sort(lhs.v_vertices.begin(), lhs.v_vertices.end());
    std::sort(rhs.u_vertices.begin(), rhs.u_vertices.end());
    std::sort(rhs.v_vertices.begin(), rhs.v_vertices.end());
    return lhs.u_vertices == rhs.u_vertices &&
           lhs.v_vertices == rhs.v_vertices;
}

bool is_small_dataset(const std::filesystem::path& data_dir) {
    for (const auto& part : data_dir.lexically_normal()) {
        const std::string name = part.string();
        if (name.rfind("small", 0) == 0) return true;
    }
    return false;
}

void write_formal_stats(std::ostream& output,
                        const char* kind,
                        const IndexBuildStats& stats) {
    output << "RESULT exp=formal_index"
           << " index=" << kind
           << " index_scope=" << index_scope_name(stats.scope)
           << " method=" << build_method_name(stats.method)
           << " attribute_pruning="
           << (stats.attribute_pruning_enabled ? 1 : 0)
           << " attribute_sets=" << stats.attribute_sets_total
           << " no_edge_attribute_sets_pruned="
           << stats.attribute_sets_no_edge_pruned
           << " equivalent_attribute_sets_pruned="
           << stats.attribute_sets_equivalent_pruned
           << " nonempty_attribute_sets=" << stats.attribute_sets_nonempty
           << " empty_attribute_sets=" << stats.attribute_sets_empty
           << " windows=" << stats.windows
           << " decompositions=" << stats.decompositions
           << " fallback_decompositions=" << stats.fallback_decompositions
           << " full_window_memberships="
           << stats.full_window_memberships
           << " streaming_records_avoided="
           << stats.streaming_records_avoided
           << " generated_records=" << stats.generated_records
           << " kept_records=" << stats.kept_records
           << " uat_peak_flat_bytes=" << stats.uat_peak_flat_bytes
           << " uat_peak_core_directory_bytes=" << stats.uat_peak_core_directory_bytes
           << " uat_core_groups=" << stats.uat_core_groups
           << " uat_completed_directions=" << stats.uat_completed_directions
           << " attribute_time_removed=" << stats.attribute_time_removed
           << " cross_core_removed=" << stats.cross_core_removed
           << " blocks=" << stats.grouped_blocks
           << " vertex_refs=" << stats.vertex_references
           << " inserted_unique_edges=" << stats.inserted_unique_edges
           << " deleted_unique_edges=" << stats.deleted_unique_edges
           << " dynamic_edge_updates=" << stats.dynamic_edge_updates
           << " dynamic_kernel=" << (stats.method == BuildMethod::Basic ? "none" :
                (!stats.dynamic_kernel_known ? "not_recorded" :
                 (stats.dynamic_layered_enabled ? "axis_layered" : "axis_planes")))
           << " dynamic_peak_bytes=" << stats.dynamic_peak_bytes
           << " affected_u_vertices_total="
           << stats.affected_u_vertices_total
           << " affected_v_vertices_total="
           << stats.affected_v_vertices_total
           << " max_affected_region_vertices="
           << stats.max_affected_region_vertices
           << " changed_u_vertices=" << stats.changed_u_vertices
           << " changed_v_vertices=" << stats.changed_v_vertices
           << " changed_frontier_entries="
           << stats.changed_frontier_entries
           << " max_changed_vertices_per_bucket="
           << stats.max_changed_vertices_per_bucket
           << " nonzero_changed_buckets="
           << stats.nonzero_changed_buckets
           << " decomposition_s=" << stats.decomposition_seconds
           << " attribute_planning_s=" << stats.attribute_planning_seconds
           << " dynamic_update_s=" << stats.dynamic_update_seconds
           << " streaming_s=" << stats.streaming_seconds
           << " mat_streaming=" << (std::string(kind) != "mat" ? "none" :
                (!stats.mat_streaming_mode_known ? "not_recorded" :
                 (stats.mat_ordered_streaming_enabled ? "ordered" : "scan")))
           << " finalize_s=" << stats.finalize_seconds
           << " build_s=" << stats.build_seconds
           << " approx_bytes=" << stats.approx_bytes
           << " approx_bytes_scope=heap_without_root\n";
    if (stats.mat_streaming_profile.enabled) {
        const auto& p = stats.mat_streaming_profile;
        output << "RESULT exp=mat_streaming_profile"
               << " generated_records=" << stats.generated_records
               << " sampled_calls=" << p.sampled_calls
               << " rejected_records=" << p.rejected_records
               << " erased_records=" << p.erased_records
               << " reject_buckets=" << p.reject_buckets
               << " reject_bounds_passed=" << p.reject_bounds_passed
               << " reject_comparisons=" << p.reject_comparisons
               << " reject_time_groups=" << p.reject_time_groups
               << " ordered_reject_buckets=" << p.ordered_reject_buckets
               << " start_group_build_records=" << p.start_group_build_records
               << " start_group_capacity_bytes=" << p.start_group_capacity_bytes
               << " remove_buckets=" << p.remove_buckets
               << " remove_bounds_passed=" << p.remove_bounds_passed
               << " remove_comparisons=" << p.remove_comparisons
               << " remove_time_groups=" << p.remove_time_groups
               << " refreshed_records=" << p.refreshed_records
               << " max_bucket_records=" << p.max_bucket_records
               << " sampled_add_s=" << p.sampled_add_seconds
               << " sampled_locator_s=" << p.sampled_locator_seconds
               << " sampled_reject_s=" << p.sampled_reject_seconds
               << " sampled_remove_s=" << p.sampled_remove_seconds
               << " sampled_append_s=" << p.sampled_append_seconds << '\n';
    }
    if (stats.dynamic_profile_enabled) {
        output << "RESULT exp=formal_index_dynamic_profile"
               << " index=" << kind
               << " method=" << build_method_name(stats.method)
               << " candidate_planes=" << stats.dynamic_candidate_planes
               << " seeded_planes=" << stats.dynamic_seeded_planes
               << " changed_planes=" << stats.dynamic_changed_planes
               << " support_adjacency_entries="
               << stats.dynamic_support_adjacency_entries
               << " cascade_queue_pops="
               << stats.dynamic_cascade_queue_pops
               << " cascade_removed_vertices="
               << stats.dynamic_cascade_removed_vertices
               << " candidate_s=" << stats.dynamic_candidate_seconds
               << " axis_s=" << stats.dynamic_axis_seconds
               << " axis_adjacency_entries=" << stats.dynamic_axis_adjacency_entries
               << " axis_frontier_changes=" << stats.dynamic_axis_frontier_changes
               << " layered_s=" << stats.dynamic_layered_seconds
               << " layered_sweeps=" << stats.dynamic_layered_sweeps
               << " layered_region_vertices=" << stats.dynamic_layered_region_vertices
               << " layered_removed_vertices=" << stats.dynamic_layered_removed_vertices
               << " layered_adjacency_entries=" << stats.dynamic_layered_adjacency_entries
               << " layered_anchor_edges=" << stats.dynamic_layered_anchor_edges
               << " support_s=" << stats.dynamic_support_seconds
               << " cascade_s=" << stats.dynamic_cascade_seconds
               << " merge_s=" << stats.dynamic_merge_seconds << '\n';
    }
}



bool formal_build_progress_enabled() {
    const char* value = std::getenv("TAC_BUILD_PROGRESS");
    if (value == nullptr || *value == '\0') return true;
    const std::string text(value);
    return text != "0" && text != "false" && text != "FALSE";
}

IndexBuildProgressCallback make_formal_progress_reporter(
    const char* kind, BuildMethod method, std::ostream& log) {
    if (!formal_build_progress_enabled()) return {};
    return [kind = std::string(kind),
            method_name = std::string(build_method_name(method)),
            &log, next_percent = 10](const IndexBuildProgress& progress)
               mutable {
        if (progress.total_attribute_sets == 0) return;
        while (next_percent <= 100 &&
               progress.completed_attribute_sets * 100 >=
                   progress.total_attribute_sets *
                       static_cast<std::size_t>(next_percent)) {
            std::ostringstream line;
            line << std::fixed << std::setprecision(6)
                 << "PROGRESS phase=index_enumeration"
                 << " index=" << kind
                 << " method=" << method_name
                 << " percent=" << next_percent;
            if (progress.joint_attributes) {
                line << " completed_attribute_pairs=" << progress.completed_attribute_sets
                     << " total_attribute_pairs=" << progress.total_attribute_sets
                     << " side=UV current_u_attribute_mask=" << progress.current_attribute_mask
                     << " current_v_attribute_mask=" << progress.current_v_attribute_mask;
            } else {
                line << " completed_attribute_sets=" << progress.completed_attribute_sets
                     << " total_attribute_sets=" << progress.total_attribute_sets
                     << " side=" << (progress.side == Side::U ? "U" : "V")
                     << " current_attribute_mask=" << progress.current_attribute_mask;
            }
            line << " elapsed_s=" << progress.elapsed_seconds
                 << " wall_time=" << current_wall_time();
            const MemorySnapshot memory = read_memory_snapshot();
            line << " rss_mb=" << kilobytes_to_mb(memory.rss_kb)
                 << " hwm_mb=" << kilobytes_to_mb(memory.hwm_kb)
                 << '\n';
            log << line.str() << std::flush;
            std::cout << line.str() << std::flush;
            next_percent += 10;
        }
    };
}

bool formal_500check_enabled() {
    const char* value = std::getenv("TAC_INDEX_500CHECK");
    if (value == nullptr) return false;
    const std::string text(value);
    return text == "1" || text == "true" || text == "on";
}

std::vector<AttributeID> used_attributes(
    const TemporalBipartiteGraph& graph, Side side) {
    std::set<AttributeID> unique;
    if (side == Side::U) {
        for (VertexID u = 0; u < graph.num_u_vertices(); ++u) {
            unique.insert(graph.u_attribute(u));
        }
    } else {
        for (VertexID v = 0; v < graph.num_v_vertices(); ++v) {
            unique.insert(graph.v_attribute(v));
        }
    }
    return {unique.begin(), unique.end()};
}

AttributeMask check_attribute_mask(
    const std::vector<AttributeID>& attributes) {
    AttributeMask mask = 0;
    for (AttributeID attribute : attributes) {
        if (attribute < 64) mask |= AttributeMask{1} << attribute;
    }
    return mask;
}

std::vector<AttributeID> check_attributes(
    const std::vector<AttributeID>& universe, std::uint64_t local_mask) {
    std::vector<AttributeID> attributes;
    for (std::size_t i = 0; i < universe.size(); ++i) {
        if ((local_mask & (std::uint64_t{1} << i)) != 0) {
            attributes.push_back(universe[i]);
        }
    }
    return attributes;
}

std::string check_query_key(const Query& query) {
    std::ostringstream output;
    output << query.window.start << ':' << query.window.end << ':'
           << query.alpha << ':' << query.beta << ':'
           << check_attribute_mask(query.u_attributes) << ':'
           << check_attribute_mask(query.v_attributes);
    return output.str();
}

struct CheckQuery {
    Query query;
    QueryResult expected;
    bool nonempty_group{};
    double online_us{};
    QueryProfile online_profile;
};

void write_query_profile(std::ostream& output, const QueryProfile& profile) {
    output << ' ' << profile.attribute_filter_us
           << ' ' << profile.directory_lookup_us
           << ' ' << profile.candidate_collect_us
           << ' ' << profile.candidate_intersection_us
           << ' ' << profile.window_recovery_us
           << ' ' << profile.peeling_us
           << ' ' << profile.directory_keys
           << ' ' << profile.directory_hits
           << ' ' << profile.posting_entries
           << ' ' << profile.vertex_refs
           << ' ' << profile.directional_candidates
           << ' ' << profile.candidate_u
           << ' ' << profile.candidate_v
           << ' ' << profile.incident_work_u
           << ' ' << profile.incident_work_v
           << ' ' << profile.unique_edges_examined
           << ' ' << profile.window_edges
           << ' ' << profile.result_u
           << ' ' << profile.result_v
           << ' ' << profile.online_fallbacks;
}

bool read_query_profile(std::istream& input, QueryProfile& profile) {
    return static_cast<bool>(
        input >> profile.attribute_filter_us
              >> profile.directory_lookup_us
              >> profile.candidate_collect_us
              >> profile.candidate_intersection_us
              >> profile.window_recovery_us
              >> profile.peeling_us
              >> profile.directory_keys
              >> profile.directory_hits
              >> profile.posting_entries
              >> profile.vertex_refs
              >> profile.directional_candidates
              >> profile.candidate_u
              >> profile.candidate_v
              >> profile.incident_work_u
              >> profile.incident_work_v
              >> profile.unique_edges_examined
              >> profile.window_edges
              >> profile.result_u
              >> profile.result_v
              >> profile.online_fallbacks);
}

double profiled_phase_total_us(const QueryProfile& profile) {
    return profile.attribute_filter_us + profile.directory_lookup_us +
        profile.candidate_collect_us + profile.candidate_intersection_us +
        profile.window_recovery_us + profile.peeling_us + profile.workspace_clear_us;
}

void write_query_profile_summary(std::ostream& log,
                                 const char* engine,
                                 const QueryProfile& profile,
                                 double measured_total_us,
                                 std::size_t query_count,
                                 QueryRoute route = QueryRoute::UV) {
    const double divisor = query_count == 0
        ? 1.0 : static_cast<double>(query_count);
    const double phase_total = profiled_phase_total_us(profile);
    log << "RESULT exp=500check_profile"
        << " engine=" << engine
        << " query_route=" << (std::string(engine) == "online" ? "online" :
                               query_route_name(route))
        << " queries=" << query_count
        << " measured_avg_us=" << measured_total_us / divisor
        << " attribute_filter_avg_us="
        << profile.attribute_filter_us / divisor
        << " directory_lookup_avg_us="
        << profile.directory_lookup_us / divisor
        << " candidate_collect_avg_us="
        << profile.candidate_collect_us / divisor
        << " candidate_intersection_avg_us="
        << profile.candidate_intersection_us / divisor
        << " window_recovery_avg_us="
        << profile.window_recovery_us / divisor
        << " peeling_avg_us=" << profile.peeling_us / divisor
        << " workspace_clear_avg_us=" << profile.workspace_clear_us / divisor
        << " grouped_aux_build_us=" << profile.grouped_aux_build_us
        << " unprofiled_avg_us="
        << std::max(0.0, measured_total_us - phase_total) / divisor
        << " directory_keys_avg=" << profile.directory_keys / divisor
        << " directory_hits_avg=" << profile.directory_hits / divisor
        << " posting_entries_avg=" << profile.posting_entries / divisor
        << " vertex_refs_avg=" << profile.vertex_refs / divisor
        << " directional_candidates_avg="
        << profile.directional_candidates / divisor
        << " candidate_u_avg=" << profile.candidate_u / divisor
        << " candidate_v_avg=" << profile.candidate_v / divisor
        << " incident_work_u_avg=" << profile.incident_work_u / divisor
        << " incident_work_v_avg=" << profile.incident_work_v / divisor
        << " unique_edges_examined_avg="
        << profile.unique_edges_examined / divisor
        << " window_edges_avg=" << profile.window_edges / divisor
        << " result_u_avg=" << profile.result_u / divisor
        << " result_v_avg=" << profile.result_v / divisor
        << " online_fallbacks=" << profile.online_fallbacks
        << " adaptive_seed_recoveries=" << profile.adaptive_seed_recoveries
        << " reversed_seed_scans=" << profile.reversed_seed_scans
        << " seed_scan_work_avg=" << (profile.seed_scan_work / divisor)
        << " opposite_scan_work_avg=" << (profile.opposite_scan_work / divisor)
        << " grouped_recoveries=" << profile.grouped_recoveries
        << " grouped_window_events_avg=" << profile.grouped_window_events / divisor
        << " grouped_duplicate_events_avg=" << profile.grouped_duplicate_events / divisor
        << " grouped_budget_fallbacks=" << profile.grouped_budget_fallbacks
        << " grouped_aux_bytes=" << profile.grouped_aux_bytes
        << " workspace_fixed_bytes=" << profile.workspace_fixed_bytes
        << " workspace_retained_bytes=" << profile.workspace_retained_bytes
        << " status=ok\n";
}

void prepare_uv_query_aux(const char* kind,
                           const TemporalBipartiteGraph& graph, std::ostream& log) {
    if (!index_uses_grouped_recovery(kind)) return;
    const bool cached = graph.grouped_v_time_aux_bytes() != 0;
    const auto start = std::chrono::steady_clock::now();
    const auto& aux = graph.grouped_v_time_adjacency();
    const double prepare_s = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    const auto memory = read_memory_snapshot();
    log << "RESULT exp=" << kind << "_query_aux recovery=grouped source="
        << (cached ? "cached" : "built") << " prepare_s=" << prepare_s
        << " aux_build_s=" << aux.build_seconds << " aux_bytes=" << aux.bytes()
        << " scratch_peak_bytes=" << aux.scratch_peak_bytes
        << " groups=" << aux.times.size() << " events=" << aux.neighbors.size()
        << " rss_mb=" << kilobytes_to_mb(memory.rss_kb)
        << " hwm_mb=" << kilobytes_to_mb(memory.hwm_kb)
        << " status=ok\n";
}

std::optional<Query> make_random_check_query(
    const TemporalBipartiteGraph& graph,
    const std::vector<AttributeID>& u_universe,
    const std::vector<AttributeID>& v_universe,
    Count alpha_cap, Count beta_cap, bool nonempty_group,
    std::mt19937_64& random) {
    const auto& times = graph.time_points();
    if (times.empty() || u_universe.empty() || v_universe.empty() ||
        u_universe.size() >= 64 || v_universe.size() >= 64) {
        return std::nullopt;
    }

    std::uniform_int_distribution<std::size_t> start_pick(
        0, times.size() - 1);
    const std::size_t start_index = start_pick(random);
    std::uniform_int_distribution<std::size_t> end_pick(
        start_index, times.size() - 1);
    const std::size_t end_index = end_pick(random);
    std::uniform_int_distribution<std::uint64_t> u_mask_pick(
        1, (std::uint64_t{1} << u_universe.size()) - 1);
    std::uniform_int_distribution<std::uint64_t> v_mask_pick(
        1, (std::uint64_t{1} << v_universe.size()) - 1);

    Count alpha = 1;
    Count beta = 1;
    if (nonempty_group) {
        static const std::vector<std::pair<Count, Count>> preferred{
            {2, 2}, {2, 1}, {1, 2}, {3, 2}, {2, 3}, {3, 3}};
        std::vector<std::pair<Count, Count>> available;
        for (const auto& pair : preferred) {
            if (pair.first <= alpha_cap && pair.second <= beta_cap) {
                available.push_back(pair);
            }
        }
        if (available.empty()) return std::nullopt;
        std::uniform_int_distribution<std::size_t> pair_pick(
            0, available.size() - 1);
        std::tie(alpha, beta) = available[pair_pick(random)];
    } else {
        std::uniform_int_distribution<Count> alpha_pick(1, alpha_cap);
        std::uniform_int_distribution<Count> beta_pick(1, beta_cap);
        alpha = alpha_pick(random);
        beta = beta_pick(random);
    }

    return Query{{times[start_index], times[end_index]},
                 check_attributes(u_universe, u_mask_pick(random)),
                 check_attributes(v_universe, v_mask_pick(random)),
                 alpha, beta};
}

std::optional<std::vector<CheckQuery>> make_500check_queries(
    const TemporalBipartiteGraph& graph, std::ostream& log) {
    constexpr std::size_t nonempty_target = 250;
    constexpr std::size_t total_target = 500;
    constexpr std::size_t max_attempts = 50000;
    constexpr std::uint64_t seed = 20260918ULL;

    const auto u_universe = used_attributes(graph, Side::U);
    const auto v_universe = used_attributes(graph, Side::V);
    const Count alpha_cap = std::min<Count>(
        10, static_cast<Count>(graph.num_v_vertices()));
    const Count beta_cap = std::min<Count>(
        10, static_cast<Count>(graph.num_u_vertices()));
    if (alpha_cap == 0 || beta_cap == 0) {
        log << "RESULT exp=500check_query_generation status=failed"
            << " reason=empty_graph\n";
        return std::nullopt;
    }

    std::mt19937_64 random(seed);
    std::unordered_set<std::string> seen;
    std::vector<CheckQuery> queries;
    queries.reserve(total_target);
    std::size_t attempts = 0;

    while (queries.size() < nonempty_target && attempts < max_attempts) {
        ++attempts;
        const auto candidate = make_random_check_query(
            graph, u_universe, v_universe, alpha_cap, beta_cap, true, random);
        if (!candidate ||
            !seen.insert(check_query_key(*candidate)).second) {
            continue;
        }
        QueryProfile online_profile;
        const auto start = std::chrono::steady_clock::now();
        QueryResult expected = online_query(graph, *candidate, &online_profile);
        const double online_us = std::chrono::duration<double, std::micro>(
            std::chrono::steady_clock::now() - start).count();
        if (expected.u_vertices.empty() || expected.v_vertices.empty()) {
            continue;
        }
        queries.push_back(
            {*candidate, std::move(expected), true, online_us,
             online_profile});
    }
    if (queries.size() != nonempty_target) {
        log << "RESULT exp=500check_query_generation status=failed"
            << " reason=insufficient_nonempty_queries"
            << " generated=" << queries.size()
            << " attempts=" << attempts << '\n';
        return std::nullopt;
    }

    while (queries.size() < total_target) {
        const auto candidate = make_random_check_query(
            graph, u_universe, v_universe, alpha_cap, beta_cap, false, random);
        if (!candidate ||
            !seen.insert(check_query_key(*candidate)).second) {
            continue;
        }
        QueryProfile online_profile;
        const auto start = std::chrono::steady_clock::now();
        QueryResult expected = online_query(graph, *candidate, &online_profile);
        const double online_us = std::chrono::duration<double, std::micro>(
            std::chrono::steady_clock::now() - start).count();
        queries.push_back(
            {*candidate, std::move(expected), false, online_us,
             online_profile});
    }

    log << "RESULT exp=500check_query_generation"
        << " seed=" << seed
        << " nonempty_queries=250 random_queries=250"
        << " nonempty_attempts=" << attempts
        << " alpha_cap=" << alpha_cap
        << " beta_cap=" << beta_cap
        << " status=ok\n";
    return queries;
}


std::vector<AttributeID> check_attributes_from_global_mask(
    const TemporalBipartiteGraph& graph, Side side, AttributeMask mask) {
    std::vector<AttributeID> attributes;
    for (AttributeID attribute : used_attributes(graph, side)) {
        if (attribute < 64 &&
            (mask & (AttributeMask{1} << attribute)) != 0) {
            attributes.push_back(attribute);
        }
    }
    return attributes;
}

bool save_500check_workload(
    const std::filesystem::path& path,
    const TemporalBipartiteGraph& graph,
    const std::vector<CheckQuery>& queries,
    std::string& error) {
    const std::filesystem::path temporary =
        std::filesystem::path(path.string() + ".tmp");
    std::ofstream output(temporary, std::ios::trunc);
    if (!output) {
        error = "open_failed";
        return false;
    }
    output << std::setprecision(17)
           << "TACC_500CHECK 5\n"
           << "GRAPH_FINGERPRINT " << graph_fingerprint(graph) << '\n'
           << "QUERY_COUNT " << queries.size() << '\n';
    for (std::size_t i = 0; i < queries.size(); ++i) {
        const CheckQuery& item = queries[i];
        output << "QUERY " << (i + 1)
               << ' ' << (item.nonempty_group ? 'N' : 'R')
               << ' ' << item.query.window.start
               << ' ' << item.query.window.end
               << ' ' << item.query.alpha
               << ' ' << item.query.beta
               << ' ' << check_attribute_mask(item.query.u_attributes)
               << ' ' << check_attribute_mask(item.query.v_attributes)
               << ' ' << item.online_us;
        write_query_profile(output, item.online_profile);
        output << ' ' << item.expected.u_vertices.size();
        for (VertexID vertex : item.expected.u_vertices) {
            output << ' ' << vertex;
        }
        output << ' ' << item.expected.v_vertices.size();
        for (VertexID vertex : item.expected.v_vertices) {
            output << ' ' << vertex;
        }
        output << '\n';
    }
    output.close();
    if (!output) {
        error = "write_failed";
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        return false;
    }

    std::error_code filesystem_error;
    std::filesystem::remove(path, filesystem_error);
    filesystem_error.clear();
    std::filesystem::rename(temporary, path, filesystem_error);
    if (filesystem_error) {
        error = "rename_failed:" + filesystem_error.message();
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        return false;
    }
    return true;
}

bool load_500check_workload(
    const std::filesystem::path& path,
    const TemporalBipartiteGraph& graph,
    std::vector<CheckQuery>& queries,
    std::string& error) {
    std::ifstream input(path);
    if (!input) {
        error = "open_failed";
        return false;
    }

    std::string magic;
    std::uint32_t version = 0;
    std::string field;
    std::uint64_t fingerprint = 0;
    std::size_t query_count = 0;
    if (!(input >> magic >> version) ||
        magic != "TACC_500CHECK" || version != 5) {
        error = "format_or_version_mismatch";
        return false;
    }
    if (!(input >> field >> fingerprint) ||
        field != "GRAPH_FINGERPRINT") {
        error = "missing_graph_fingerprint";
        return false;
    }
    if (fingerprint != graph_fingerprint(graph)) {
        error = "graph_fingerprint_mismatch";
        return false;
    }
    if (!(input >> field >> query_count) || field != "QUERY_COUNT" ||
        query_count != 500) {
        error = "query_count_mismatch";
        return false;
    }

    std::vector<CheckQuery> loaded;
    loaded.reserve(query_count);
    for (std::size_t i = 0; i < query_count; ++i) {
        std::string tag;
        std::size_t query_id = 0;
        char group = '\0';
        Query query;
        AttributeMask u_mask = 0;
        AttributeMask v_mask = 0;
        double online_us = 0.0;
        QueryProfile online_profile;
        if (!(input >> tag >> query_id >> group
                    >> query.window.start >> query.window.end
                    >> query.alpha >> query.beta >> u_mask >> v_mask
                    >> online_us) ||
            tag != "QUERY" || query_id != i + 1 ||
            (group != 'N' && group != 'R')) {
            error = "invalid_query_record";
            return false;
        }
        if (!read_query_profile(input, online_profile)) {
            error = "invalid_online_profile";
            return false;
        }
        std::size_t u_count = 0;
        if (!(input >> u_count)) {
            error = "missing_u_answer_size";
            return false;
        }
        query.u_attributes =
            check_attributes_from_global_mask(graph, Side::U, u_mask);
        query.v_attributes =
            check_attributes_from_global_mask(graph, Side::V, v_mask);
        if (query.window.start > query.window.end ||
            query.alpha == 0 || query.beta == 0 ||
            query.u_attributes.empty() || query.v_attributes.empty() ||
            check_attribute_mask(query.u_attributes) != u_mask ||
            check_attribute_mask(query.v_attributes) != v_mask) {
            error = "invalid_query_condition";
            return false;
        }

        QueryResult expected;
        expected.u_vertices.resize(u_count);
        for (VertexID& vertex : expected.u_vertices) {
            if (!(input >> vertex)) {
                error = "truncated_u_answer";
                return false;
            }
        }
        std::size_t v_count = 0;
        if (!(input >> v_count)) {
            error = "missing_v_answer_size";
            return false;
        }
        expected.v_vertices.resize(v_count);
        for (VertexID& vertex : expected.v_vertices) {
            if (!(input >> vertex)) {
                error = "truncated_v_answer";
                return false;
            }
        }
        if (group == 'N' &&
            (expected.u_vertices.empty() || expected.v_vertices.empty())) {
            error = "nonempty_group_contains_empty_answer";
            return false;
        }
        loaded.push_back(
            {std::move(query), std::move(expected), group == 'N', online_us,
             online_profile});
    }
    queries = std::move(loaded);
    return true;
}

std::optional<std::vector<CheckQuery>> load_or_create_500check_workload(
    const std::filesystem::path& path,
    const TemporalBipartiteGraph& graph,
    std::ostream& log,
    bool& loaded_existing) {
    loaded_existing = false;
    std::vector<CheckQuery> queries;
    std::string load_error;
    if (std::filesystem::exists(path)) {
        if (load_500check_workload(path, graph, queries, load_error)) {
            loaded_existing = true;
            log << "RESULT exp=500check_workload"
                << " source=loaded"
                << " queries=" << queries.size()
                << " graph_fingerprint=" << graph_fingerprint(graph)
                << " path=" << path.string()
                << " status=ok\n";
            return queries;
        }
        log << "RESULT exp=500check_workload_load"
            << " path=" << path.string()
            << " status=failed reason=" << load_error
            << " action=regenerate\n";
    } else {
        log << "RESULT exp=500check_workload_load"
            << " path=" << path.string()
            << " status=missing action=generate\n";
    }

    auto generated = make_500check_queries(graph, log);
    if (!generated) return std::nullopt;
    std::string save_error;
    if (!save_500check_workload(path, graph, *generated, save_error)) {
        log << "RESULT exp=500check_workload"
            << " source=generated"
            << " path=" << path.string()
            << " status=failed reason=" << save_error << '\n';
        return std::nullopt;
    }
    log << "RESULT exp=500check_workload"
        << " source=generated"
        << " queries=" << generated->size()
        << " graph_fingerprint=" << graph_fingerprint(graph)
        << " path=" << path.string()
        << " status=ok\n";
    return generated;
}

bool has_filename_timestamp(const std::string& name) {
    if (name.size() < 15 || name[8] != '_') return false;
    for (std::size_t i = 0; i < 15; ++i) {
        if (i == 8) continue;
        if (name[i] < '0' || name[i] > '9') return false;
    }
    return true;
}

std::string available_formal_stamp(const DatasetResultPaths& paths,
                                   const std::string& suffix, bool check) {
    const std::string timestamp = current_filename_timestamp();
    for (std::size_t run = 0;; ++run) {
        const std::string stamp = timestamp +
            (run == 0 ? "" : "_run" + std::to_string(run));
        const auto log = (check ? paths.check : paths.root) / (stamp + suffix + ".log");
        const auto binary = paths.indexes / (stamp + suffix + ".bin");
        const auto text = paths.indexes / (stamp + suffix + "_text.log");
        const auto intervals_text = paths.indexes / (stamp + suffix + "_intervals_text.log");
        if (!std::filesystem::exists(log) && !std::filesystem::exists(binary) &&
            !std::filesystem::exists(text) && !std::filesystem::exists(intervals_text)) return stamp;
    }
}

std::vector<std::filesystem::path> timestamped_index_candidates(
    const DatasetResultPaths& paths,
    const char* kind,
    const std::string& method_name) {
    const std::string suffix = "_" + paths.dataset + "_" + kind + "_" +
        method_name + "_uv" + "_index.bin";
    std::vector<std::filesystem::path> candidates;
    std::error_code error;
    for (std::filesystem::directory_iterator iterator(paths.indexes, error),
         end;
         !error && iterator != end;
         iterator.increment(error)) {
        std::error_code type_error;
        if (!iterator->is_regular_file(type_error) || type_error) continue;
        const std::string name = iterator->path().filename().string();
        if (!has_filename_timestamp(name) || name.size() < 15 + suffix.size()) {
            continue;
        }
        const std::string extra = name.substr(15, name.size() - 15 - suffix.size());
        if (!extra.empty() && (extra.size() <= 4 || extra.compare(0, 4, "_run") != 0 ||
            !std::all_of(extra.begin() + 4, extra.end(), [](char ch) {
                return ch >= '0' && ch <= '9';
            }))) continue;
        if (name.compare(name.size() - suffix.size(), suffix.size(), suffix) ==
            0) {
            candidates.push_back(iterator->path());
        }
    }
    std::sort(candidates.begin(), candidates.end(),
              [&suffix](const auto& lhs, const auto& rhs) {
                  const std::string a = lhs.filename().string();
                  const std::string b = rhs.filename().string();
                  const int timestamp_order = a.compare(0, 15, b, 0, 15);
                  if (timestamp_order != 0) return timestamp_order > 0;
                  const auto run = [&](const std::string& name) {
                      const auto length = name.size() - 15 - suffix.size();
                      std::string value = length == 0 ? "0" : name.substr(19, length - 4);
                      const auto first = value.find_first_not_of('0');
                      return first == std::string::npos ? std::string("0") : value.substr(first);
                  };
                  const auto ar = run(a), br = run(b);
                  if (ar.size() != br.size()) return ar.size() > br.size();
                  return ar > br;
              });
    return candidates;
}

template <typename Index, typename Builder, typename QueryFunction>
int run_formal_500check(const char* kind,
                        const std::filesystem::path& data_dir,
                        const TemporalBipartiteGraph& graph,
                        BuildMethod method,
                        Builder builder,
                        QueryFunction query_function) {
    const DatasetResultPaths result_paths = dataset_result_paths(data_dir);
    const std::string& dataset = result_paths.dataset;
    const std::string method_name = build_method_name(method);
    const IndexScope scope = configured_index_scope();
    const QueryRoute resolved_route = effective_query_route(scope);
    const std::string route = query_route_name(resolved_route);
    std::error_code directory_error;
    if (!create_result_directories(result_paths, directory_error)) {
        std::cerr << "Cannot create 500check output directory: "
                  << directory_error.message() << '\n';
        return 1;
    }

    const std::string check_suffix = ("_route_" + route) +
        "_scope_uv" + "_" + dataset +
        "_" + kind + "_" + method_name + "_500check";
    const std::string stamp = available_formal_stamp(result_paths, check_suffix, true);
    const std::string base = stamp + check_suffix;
    const std::filesystem::path log_path =
        result_paths.check / (base + ".log");
    const std::string index_suffix = "_" + dataset + "_" + kind + "_" +
        method_name + "_uv" + "_index.bin";
    const std::filesystem::path workload_path =
        result_paths.check / (dataset + "_500check_online.log");
    std::ofstream log(log_path);
    if (!log) {
        std::cerr << "Cannot create log file: " << log_path.string() << '\n';
        return 1;
    }
    log << std::unitbuf << std::fixed << std::setprecision(6);
    log << "# TACC 500-query index verification\n"
        << "Generated at: " << current_wall_time() << '\n'
        << "RESULT exp=500check_dataset"
        << " index=" << kind
        << " index_scope=" << index_scope_name(scope)
        << " query_route=" << route
        << " method=" << method_name
        << " dataset=" << data_dir.string()
        << " graph_fingerprint=" << graph_fingerprint(graph)
        << " u_vertices=" << graph.num_u_vertices()
        << " v_vertices=" << graph.num_v_vertices()
        << " temporal_edges=" << graph.num_edges()
        << " unique_edges=" << graph.num_unique_edges()
        << " graph_query_aux_bytes=" << graph.temporal_query_aux_bytes()
        << " timestamps=" << graph.time_points().size()
        << " index_dir=" << result_paths.indexes.string()
        << " index_suffix=" << index_suffix
        << " workload_path=" << workload_path.string() << '\n';

    log << "PROGRESS phase=500check_workload status=load_or_generate_start"
        << " wall_time=" << current_wall_time() << '\n';
    write_memory_guard(log, active_memory_guard());
    bool workload_loaded = false;
    auto queries = load_or_create_500check_workload(
        workload_path, graph, log, workload_loaded);
    if (!queries) {
        log.close();
        return finish_log(log_path, 3);
    }
    log << "PROGRESS phase=500check_workload status=ready"
        << " source=" << (workload_loaded ? "loaded" : "generated")
        << " queries=" << queries->size()
        << " wall_time=" << current_wall_time() << '\n';

    Index index;
    bool cache_hit = false;
    bool built = false;
    double load_seconds = 0.0;
    double dump_seconds = 0.0;
    std::filesystem::path index_path;
    const auto candidates = timestamped_index_candidates(
        result_paths, kind, method_name);
    log << "PROGRESS phase=500check_index status=discover_done"
        << " candidates=" << candidates.size()
        << " wall_time=" << current_wall_time() << '\n';
    for (const auto& candidate : candidates) {
        Index loaded;
        std::string candidate_error;
        const auto start = std::chrono::steady_clock::now();
        IndexScope marker;
        bool compatible = read_index_scope(candidate.string(), marker, &candidate_error);
        if (compatible) compatible = load_index_binary(
            candidate.string(), graph, loaded, &candidate_error);
        const double candidate_load_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        load_seconds += candidate_load_seconds;
        if (compatible && loaded.stats().method != method) {
            compatible = false;
            candidate_error = "build_method_mismatch";
        }
        log << "RESULT exp=500check_index_candidate"
            << " status=" << (compatible ? "ok" : "rejected")
            << " load_s=" << candidate_load_seconds
            << " path=" << candidate.string();
        if (!compatible) log << " reason=" << candidate_error;
        log << '\n';
        if (compatible) {
            index = std::move(loaded);
            index_path = candidate;
            cache_hit = true;
            break;
        }
    }
    if (cache_hit) {
        log << "RESULT exp=500check_index_load status=ok"
            << " candidates=" << candidates.size()
            << " load_s=" << load_seconds
            << " path=" << index_path.string() << '\n';
    } else {
        const std::string index_stamp = available_formal_stamp(
            result_paths, index_suffix.substr(0, index_suffix.size() - 4), false);
        index_path = result_paths.indexes /
            (index_stamp + index_suffix);
        log << "RESULT exp=500check_index_load status="
            << (candidates.empty() ? "missing" : "no_compatible")
            << " candidates=" << candidates.size()
            << " build_path=" << index_path.string() << '\n';
    }

    if (!cache_hit) {
        MemorySnapshot build_memory;
        log << "PROGRESS phase=500check_index status=build_start"
            << " wall_time=" << current_wall_time() << '\n';
        try {
            const auto progress =
                make_formal_progress_reporter(kind, method, log);
            index = builder(graph, method, progress);
            build_memory = read_memory_snapshot();
            built = true;
        } catch (const std::bad_alloc&) {
            log << "RESULT exp=500check_index_build status=failed reason=bad_alloc"
                << " memory_guard_percent=" << active_memory_guard().percent << '\n';
            write_bad_alloc_memory(log, "500check_index_build");
            log.close();
            return finish_log(log_path, 2);
        } catch (const std::exception& error) {
            log << "RESULT exp=500check_index_build status=failed reason="
                << error.what() << '\n';
            log.close();
            return finish_log(log_path, 2);
        }
        log << "PROGRESS phase=500check_index status=build_done"
            << " wall_time=" << current_wall_time() << '\n';
        write_formal_stats(log, kind, index.stats());

        log << "RESULT exp=formal_index_build_memory index=" << kind
            << " build_peak_rss_bytes=" << std::uint64_t(build_memory.hwm_kb) * 1024
            << " build_peak_rss_mb=" << kilobytes_to_mb(build_memory.hwm_kb)
            << " build_peak_rss_G=" << space_G(std::uint64_t(build_memory.hwm_kb) * 1024)
            << " peak_scope=process_through_build peak_sample=before_dump_reload"
            << " prior_workload_generation_included=1 snapshot_status="
            << (build_memory.complete ? "ok" : "partial") << '\n';

        std::string dump_error;
        const auto dump_start = std::chrono::steady_clock::now();
        const bool dump_ok = save_index_binary(
            index_path.string(), graph, index, &dump_error);
        dump_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - dump_start).count();
        log << "RESULT exp=500check_index_dump"
            << " status=" << (dump_ok ? "ok" : "failed")
            << " dump_s=" << dump_seconds
            << " path=" << index_path.string();
        if (!dump_ok) log << " reason=" << dump_error;
        log << '\n';
        if (!dump_ok) {
            log.close();
            return finish_log(log_path, 2);
        }
    } else {
        write_formal_stats(log, kind, index.stats());
        log << "RESULT exp=formal_index_build_memory index=" << kind
            << " build_peak_rss_status=not_measured source=loaded\n";
    }
    const IndexSpaceBreakdown index_space = formal_index_space(index);
    write_formal_index_space(log, kind, index_space,
                             index.stats().approx_bytes);

    log << "RESULT exp=500check_index_ready"
        << " source=" << (cache_hit ? "loaded" : "built")
        << " method=" << build_method_name(index.stats().method)
        << " hash=" << formal_index_hash(index)
        << " load_s=" << load_seconds
        << " dump_s=" << dump_seconds
        << " path=" << index_path.string()
        << " status=ok\n"
        << "PROGRESS phase=500check_queries status=query_start"
        << " workload_source="
        << (workload_loaded ? "loaded" : "generated")
        << " workload_path=" << workload_path.string()
        << " wall_time=" << current_wall_time() << '\n';

    std::size_t mismatches = 0;
    std::size_t empty_random = 0;
    double total_online_us = 0.0;
    double total_index_us = 0.0;
    QueryProfile total_online_profile;
    QueryProfile total_index_profile;
    prepare_uv_query_aux(kind, graph, log);
    write_query_index_space(log, kind, index_space, graph,
        index_uses_grouped_recovery(kind));
    for (std::size_t i = 0; i < queries->size(); ++i) {
        CheckQuery& item = (*queries)[i];
        QueryProfile fresh_online_profile;
        const auto online_start = std::chrono::steady_clock::now();
        const QueryResult fresh_online = online_query(graph, item.query, &fresh_online_profile);
        item.online_us = std::chrono::duration<double, std::micro>(
            std::chrono::steady_clock::now() - online_start).count();
        item.online_profile = fresh_online_profile;
        if (!same_result(fresh_online, item.expected)) {
            log << "RESULT exp=500check status=failed reason=online_reference_mismatch"
                << " id=" << (i + 1) << '\n';
            log.close();
            return finish_log(log_path, 3);
        }
        total_online_us += item.online_us;
        total_online_profile += item.online_profile;
        if (!item.nonempty_group &&
            (item.expected.u_vertices.empty() ||
             item.expected.v_vertices.empty())) {
            ++empty_random;
        }
        QueryProfile index_profile;
        const auto start = std::chrono::steady_clock::now();
        QueryResult actual = query_function(
            index, graph, item.query, &index_profile);
        const double index_us = std::chrono::duration<double, std::micro>(
            std::chrono::steady_clock::now() - start).count();
        total_index_us += index_us;
        total_index_profile += index_profile;
        const bool same = same_result(actual, item.expected);
        if (!same) ++mismatches;
        log << "CHECK id=" << (i + 1)
            << " query_route=" << route
            << " group="
            << (item.nonempty_group ? "nonempty" : "random")
            << " ts=" << item.query.window.start
            << " te=" << item.query.window.end
            << " alpha=" << item.query.alpha
            << " beta=" << item.query.beta
            << " u_mask=" << check_attribute_mask(item.query.u_attributes)
            << " v_mask=" << check_attribute_mask(item.query.v_attributes)
            << " expected_u=" << item.expected.u_vertices.size()
            << " expected_v=" << item.expected.v_vertices.size()
            << " actual_u=" << actual.u_vertices.size()
            << " actual_v=" << actual.v_vertices.size()
            << " online_us=" << item.online_us
            << " online_implementation=plain_time_scan"
            << " online_timing=fresh"
            << " index_us=" << index_us
            << " online_attribute_us="
            << item.online_profile.attribute_filter_us
            << " online_window_us="
            << item.online_profile.window_recovery_us
            << " online_peeling_us=" << item.online_profile.peeling_us
            << " index_attribute_us=" << index_profile.attribute_filter_us
            << " index_directory_us="
            << index_profile.directory_lookup_us
            << " index_collect_us="
            << index_profile.candidate_collect_us
            << " index_intersection_us="
            << index_profile.candidate_intersection_us
            << " index_window_us=" << index_profile.window_recovery_us
            << " index_peeling_us=" << index_profile.peeling_us
            << " index_candidate_u=" << index_profile.candidate_u
            << " index_candidate_v=" << index_profile.candidate_v
            << " index_edges_examined="
            << index_profile.unique_edges_examined
            << " index_window_edges=" << index_profile.window_edges
            << " index_adaptive_recovery=" << index_profile.adaptive_seed_recoveries
            << " index_reversed_scan=" << index_profile.reversed_seed_scans
            << " index_seed_scan_work=" << index_profile.seed_scan_work
            << " index_opposite_scan_work=" << index_profile.opposite_scan_work
            << " index_grouped_recovery=" << index_profile.grouped_recoveries
            << " index_grouped_window_events=" << index_profile.grouped_window_events
            << " index_grouped_duplicate_events=" << index_profile.grouped_duplicate_events
            << " index_grouped_budget_fallback=" << index_profile.grouped_budget_fallbacks
            << " index_workspace_clear_us=" << index_profile.workspace_clear_us
            << " index_workspace_retained_bytes=" << index_profile.workspace_retained_bytes
            << " same=" << (same ? 1 : 0) << '\n';
    }

    const bool success = mismatches == 0;
    const MemorySnapshot check_memory = read_memory_snapshot();
    write_query_profile_summary(
        log, "online", total_online_profile, total_online_us,
        queries->size());
    write_query_profile_summary(
        log, kind, total_index_profile, total_index_us, queries->size(), resolved_route);
    log << "RESULT exp=500check"
        << " online_implementation=plain_time_scan"
        << " online_timing=fresh"
        << " index=" << kind
        << " index_scope=" << index_scope_name(index.stats().scope)
        << " query_route=" << route
        << " method=" << method_name
        << " source=" << (cache_hit ? "loaded" : "built")
        << " seed_recovery=" << seed_recovery_name(configured_seed_recovery())
        << " mat_recovery=" << (std::string(kind) == "mat"
            ? mat_recovery_name(configured_mat_recovery()) : "not_applicable")
        << " uat_recovery=" << (std::string(kind) == "uat"
            ? uat_recovery_name(configured_uat_recovery()) : "not_applicable")
        << " built=" << (built ? 1 : 0)
        << " workload_source="
        << (workload_loaded ? "loaded" : "generated")
        << " online_executed=1"
        << " online_executed_queries=" << queries->size()
        << " workload_path=" << workload_path.string()
        << " queries=" << queries->size()
        << " nonempty_queries=250 random_queries=250"
        << " empty_random_queries=" << empty_random
        << " mismatches=" << mismatches
        << " baseline_online_total_us=" << total_online_us
        << " index_total_us=" << total_index_us
        << " baseline_online_avg_us="
        << (total_online_us / queries->size())
        << " index_avg_us=" << (total_index_us / queries->size())
        << " rss_mb=" << kilobytes_to_mb(check_memory.rss_kb)
        << " hwm_mb=" << kilobytes_to_mb(check_memory.hwm_kb)
        << " peak_rss_bytes="
        << static_cast<std::uint64_t>(check_memory.hwm_kb) * 1024ULL
        << " status=" << (success ? "ok" : "failed")
        << " index_path=" << index_path.string() << '\n';
    log.close();
    return finish_log(log_path, success ? 0 : 4);
}

template <typename Index, typename Builder, typename QueryFunction>
int run_formal_build(const char* kind,
                     const std::filesystem::path& data_dir,
                     const TemporalBipartiteGraph& graph,
                     BuildMethod method,
                     const std::optional<Query>& query,
                     Builder builder,
                     QueryFunction query_function) {
    const DatasetResultPaths result_paths = dataset_result_paths(data_dir);
    const std::string& dataset_name = result_paths.dataset;
    const std::string method_name = build_method_name(method);
    const IndexScope scope = configured_index_scope();
    const QueryRoute route = effective_query_route(scope);
    std::error_code directory_error;
    if (!create_result_directories(result_paths, directory_error)) {
        std::cerr << "Cannot create result directory: "
                  << result_paths.root.string() << ": "
                  << directory_error.message() << '\n';
        return 1;
    }

    const std::string suffix = "_" + dataset_name + "_" + kind +
        "_" + method_name + "_uv" + "_index";
    const std::string base = available_formal_stamp(result_paths, suffix, false) + suffix;
    const std::filesystem::path log_path =
        result_paths.root / (base + ".log");
    const std::filesystem::path binary_path = result_paths.indexes /
        (base + ".bin");
    const std::filesystem::path text_path =
        result_paths.indexes / (base + "_text.log");
    const std::filesystem::path intervals_text_path =
        result_paths.indexes / (base + "_intervals_text.log");
    std::ofstream log(log_path);
    if (!log) {
        std::cerr << "Cannot create log file: " << log_path.string() << '\n';
        return 1;
    }
    log << std::unitbuf << std::fixed << std::setprecision(6);
    const std::uint64_t dataset_fingerprint = graph_fingerprint(graph);
    log << "# TACC formal index experiment\n"
        << "Generated at: " << current_wall_time() << '\n'
        << "RESULT exp=formal_index_dataset"
        << " index=" << kind
        << " index_scope=" << index_scope_name(scope)
        << " query_route=" << query_route_name(route)
        << " method=" << build_method_name(method)
        << " dataset=" << dataset_name
        << " dataset_path=" << data_dir.string()
        << " graph_fingerprint=" << dataset_fingerprint
        << " u_vertices=" << graph.num_u_vertices()
        << " v_vertices=" << graph.num_v_vertices()
        << " temporal_edges=" << graph.num_edges()
        << " unique_edges=" << graph.num_unique_edges()
        << " graph_query_aux_bytes=" << graph.temporal_query_aux_bytes()
        << " timestamps=" << graph.time_points().size() << '\n';
    write_memory_guard(log, active_memory_guard());
    log << "PROGRESS phase=formal_index status=build_start"
        << " wall_time=" << current_wall_time() << '\n';
    Index index;
    MemorySnapshot build_memory;
    try {
        const auto progress =
            make_formal_progress_reporter(kind, method, log);
        index = builder(graph, method, progress);
        build_memory = read_memory_snapshot();
    } catch (const std::bad_alloc&) {
        log << "RESULT exp=formal_index status=failed reason=bad_alloc"
            << " memory_guard_percent=" << active_memory_guard().percent << '\n';
        write_bad_alloc_memory(log, "formal_index_build");
        log.close();
        return finish_log(log_path, 2);
    } catch (const std::exception& error) {
        log << "RESULT exp=formal_index status=failed reason="
            << error.what() << '\n';
        log.close();
        return finish_log(log_path, 2);
    }
    log << "PROGRESS phase=formal_index status=build_done"
        << " wall_time=" << current_wall_time()
        << " rss_mb=" << kilobytes_to_mb(build_memory.rss_kb)
        << " hwm_mb=" << kilobytes_to_mb(build_memory.hwm_kb) << '\n';
    write_formal_stats(log, kind, index.stats());
    const IndexBuildStats build_stats = index.stats();
    const IndexSpaceBreakdown index_space = formal_index_space(index);
    write_formal_index_space(log, kind, index_space,
                             index.stats().approx_bytes);
    log << "RESULT exp=formal_index_build_memory"
        << " index=" << kind
        << " dataset=" << dataset_name
        << " build_end_rss_bytes="
        << static_cast<std::uint64_t>(build_memory.rss_kb) * 1024ULL
        << " build_end_rss_mb=" << kilobytes_to_mb(build_memory.rss_kb)
        << " build_peak_rss_bytes="
        << static_cast<std::uint64_t>(build_memory.hwm_kb) * 1024ULL
        << " build_peak_rss_mb=" << kilobytes_to_mb(build_memory.hwm_kb)
        << " build_peak_rss_G=" << space_G(std::uint64_t(build_memory.hwm_kb) * 1024)
        << " peak_scope=process_through_build peak_sample=before_dump_reload"
        << " snapshot_status=" << (build_memory.complete ? "ok" : "partial")
        << '\n';

    const std::uint64_t hash_before = formal_index_hash(index);
    const bool small_dataset = is_small_dataset(data_dir);
    const bool write_text = small_dataset || dataset_name == "01_br";
    const bool write_intervals_text = small_dataset && std::is_same_v<Index, VATIndex>;
    bool intervals_text_ok = true;
    bool text_ok = true;
    double text_seconds = 0.0;
    std::uintmax_t text_bytes = 0;
    std::error_code text_size_error;
    std::string text_error;
    if (write_text) {
        const auto start = std::chrono::steady_clock::now();
        if constexpr (std::is_same_v<Index, VATIndex>)
            text_ok = save_index_text(text_path.string(), graph, index, &text_error,
                                     data_dir.string(), std::uint64_t(build_memory.hwm_kb) * 1024);
        else
            text_ok = save_index_text(text_path.string(), graph, index, &text_error);
        text_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        if (text_ok) {
            text_bytes = std::filesystem::file_size(
                text_path, text_size_error);
        }
        log << "RESULT exp=formal_index_text_dump"
            << " index=" << kind
            << " dataset_type=" << (small_dataset ? "small" : "real")
            << " dataset=" << dataset_name
            << " status=" << (text_ok ? "ok" : "failed")
            << " dump_s=" << text_seconds
            << " text_bytes=" << (text_size_error ? 0 : text_bytes)
            << " text_mb="
            << (text_size_error ? 0.0 : bytes_to_mb(text_bytes))
            << " path=" << text_path.string();
        if (!text_ok) log << " reason=" << text_error;
        log << '\n';
    } else {
        log << "RESULT exp=formal_index_text_dump"
            << " index=" << kind
            << " dataset_type=real status=skipped"
            << " reason=binary_only\n";
    }

    if constexpr (std::is_same_v<Index, VATIndex>) {
        if (write_intervals_text) {
            std::string error;
            const auto start = std::chrono::steady_clock::now();
            intervals_text_ok = save_index_intervals_text(intervals_text_path.string(), graph, index,
                &error, data_dir.string(), std::uint64_t(build_memory.hwm_kb) * 1024);
            const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            log << "RESULT exp=formal_index_intervals_text_dump index=" << kind
                << " dataset=" << dataset_name << " status=" << (intervals_text_ok ? "ok" : "failed")
                << " view=expanded_intervals physical_storage_unchanged=1 dump_s=" << seconds
                << " path=" << intervals_text_path.string();
            if (!intervals_text_ok) log << " reason=" << error;
            log << '\n';
        }
    }

    std::string dump_error;
    const auto dump_start = std::chrono::steady_clock::now();
    const bool dump_ok = save_index_binary(
        binary_path.string(), graph, index, &dump_error);
    const double dump_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - dump_start).count();
    log << "RESULT exp=formal_index_dump"
        << " index=" << kind
        << " format=binary"
        << " status=" << (dump_ok ? "ok" : "failed")
        << " dump_s=" << dump_seconds
        << " path=" << binary_path.string();
    if (!dump_ok) log << " reason=" << dump_error;
    log << '\n';

    index = Index{};
    Index loaded;
    std::string load_error;
    const auto load_start = std::chrono::steady_clock::now();
    const bool load_ok = dump_ok && load_index_binary(
        binary_path.string(), graph, loaded, &load_error);
    const double load_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - load_start).count();
    const std::uint64_t hash_after =
        load_ok ? formal_index_hash(loaded) : 0;
    const bool load_same = load_ok && hash_before == hash_after &&
        build_stats.vertex_references ==
            loaded.stats().vertex_references &&
        build_stats.kept_records == loaded.stats().kept_records;

    log << "RESULT exp=formal_index_load"
        << " index=" << kind
        << " status=" << (load_ok ? "ok" : "failed")
        << " load_s=" << load_seconds
        << " same=" << (load_same ? 1 : 0)
        << " hash_before=" << hash_before
        << " hash_after=" << hash_after
        << " vertex_refs="
        << (load_ok ? loaded.stats().vertex_references : 0);
    if (!load_ok) {
        log << " reason="
            << (dump_ok ? load_error : "binary_dump_failed");
    }
    log << '\n';

    bool query_ok = true;
    if (query && load_same) {
        prepare_uv_query_aux(kind, graph, log);
        write_query_index_space(log, kind, index_space, graph,
            index_uses_grouped_recovery(kind));
        const auto query_start = std::chrono::steady_clock::now();
        const QueryResult result = query_function(loaded, graph, *query);
        const double query_us =
            std::chrono::duration<double, std::micro>(
                std::chrono::steady_clock::now() - query_start).count();
        const QueryResult expected = online_query(graph, *query);
        query_ok = same_result(result, expected);
        log << "RESULT exp=formal_index_query"
            << " index=" << kind
            << " mat_recovery=" << (std::string(kind) == "mat"
                ? mat_recovery_name(configured_mat_recovery()) : "not_applicable")
            << " uat_recovery=" << (std::string(kind) == "uat"
                ? uat_recovery_name(configured_uat_recovery()) : "not_applicable")
            << " index_scope=" << index_scope_name(loaded.stats().scope)
            << " query_route=" << query_route_name(effective_query_route(loaded.stats().scope))
            << " source=loaded_binary"
            << " query_us=" << query_us
            << " u_vertices=" << result.u_vertices.size()
            << " v_vertices=" << result.v_vertices.size()
            << " online_same=" << (query_ok ? 1 : 0)
            << " status=" << (query_ok ? "ok" : "failed") << '\n';
        write_vertices(log, "U", result.u_vertices);
        write_vertices(log, "V", result.v_vertices);
    } else if (query) {
        query_ok = false;
        log << "RESULT exp=formal_index_query"
            << " index=" << kind
            << " status=failed reason=index_load_validation_failed\n";
    }

    const bool success = text_ok && intervals_text_ok && dump_ok && load_same && query_ok;
    log << "RESULT exp=formal_index_pipeline"
        << " index=" << kind
        << " index_scope=" << index_scope_name(build_stats.scope)
        << " dataset=" << dataset_name
        << " method=" << build_method_name(method)
        << " build_peak_rss_bytes="
        << static_cast<std::uint64_t>(build_memory.hwm_kb) * 1024ULL
        << " build_peak_rss_mb=" << kilobytes_to_mb(build_memory.hwm_kb)
        << " build_peak_rss_G=" << space_G(std::uint64_t(build_memory.hwm_kb) * 1024)
        << " build_s=" << build_stats.build_seconds
        << " index_bytes=" << index_space.index_bytes()
        << " index_mb=" << bytes_to_mb(index_space.index_bytes())
        << " index_G=" << space_G(index_space.index_bytes())
        << " size_model=structural_v2 unit_mb_bytes=1048576 unit_G_bytes=1073741824"
        << " dump_s=" << dump_seconds
        << " load_s=" << load_seconds
        << " text_dump_status="
        << (write_text ? (text_ok ? "ok" : "failed") : "skipped")
        << " intervals_text_dump_status="
        << (write_intervals_text ? (intervals_text_ok ? "ok" : "failed") : "skipped")
        << " dump_status=" << (dump_ok ? "ok" : "failed")
        << " load_status=" << (load_ok ? "ok" : "failed")
        << " load_same=" << (load_same ? 1 : 0)
        << " status=" << (success ? "ok" : "failed")
        << " index_path=" << binary_path.string() << '\n';
    log.close();
    return finish_log(log_path, success ? 0 : 3);
}

template <typename Index, typename QueryFunction>
int run_formal_loaded_query(const char* kind,
                            const std::filesystem::path& index_path,
                            const TemporalBipartiteGraph& graph,
                            const Query& query,
                            QueryFunction query_function) {
    Index index;
    std::string error;
    const auto load_start = std::chrono::steady_clock::now();
    const bool load_ok = load_index_binary(
        index_path.string(), graph, index, &error);
    const double load_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - load_start).count();
    if (!load_ok) {
        std::cerr << "RESULT exp=formal_index_query_load"
                  << " index=" << kind
                  << " status=failed reason=" << error
                  << " path=" << index_path.string() << '\n';
        return 2;
    }

    QueryRoute route;
    try {
        route = effective_query_route(index.stats().scope);
    } catch (const std::invalid_argument& route_error) {
        std::cerr << "RESULT exp=formal_index_query_load index=" << kind
                  << " index_scope=" << index_scope_name(index.stats().scope)
                  << " status=failed reason=" << route_error.what() << '\n';
        return 2;
    }
    prepare_uv_query_aux(kind, graph, std::cout);
    write_formal_index_space(std::cout, kind, formal_index_space(index), index.stats().approx_bytes);
    write_query_index_space(std::cout, kind, formal_index_space(index), graph,
        index_uses_grouped_recovery(kind));
    const auto query_start = std::chrono::steady_clock::now();
    const QueryResult result = query_function(index, graph, query);
    const double query_us = std::chrono::duration<double, std::micro>(
        std::chrono::steady_clock::now() - query_start).count();
    const QueryResult expected = online_query(graph, query);
    const bool ok = same_result(result, expected);
    std::cout << std::fixed << std::setprecision(6)
              << "RESULT exp=formal_index_query_load"
              << " index=" << kind
              << " mat_recovery=" << (std::string(kind) == "mat"
                  ? mat_recovery_name(configured_mat_recovery()) : "not_applicable")
              << " uat_recovery=" << (std::string(kind) == "uat"
                  ? uat_recovery_name(configured_uat_recovery()) : "not_applicable")
              << " index_scope=" << index_scope_name(index.stats().scope)
              << " query_route=" << query_route_name(route)
              << " method=" << build_method_name(index.stats().method)
              << " load_s=" << load_seconds
              << " query_us=" << query_us
              << " hash=" << formal_index_hash(index)
              << " u_vertices=" << result.u_vertices.size()
              << " v_vertices=" << result.v_vertices.size()
              << " online_same=" << (ok ? 1 : 0)
              << " status=" << (ok ? "ok" : "failed") << '\n';
    write_vertices(std::cout, "U", result.u_vertices);
    write_vertices(std::cout, "V", result.v_vertices);
    return ok ? 0 : 3;
}

}

IndexBuildProgressCallback make_index_build_progress_reporter(
    const char* kind, BuildMethod method, std::ostream& log) {
    return make_formal_progress_reporter(kind, method, log);
}

int run_index_experiment(int argc, char** argv) {
    const bool query_only =
        argc >= 2 && std::string(argv[1]) == "index-query";
    if ((!query_only && argc != 4 && argc != 10) ||
        (query_only && argc != 11)) {
        formal_usage(argv[0]);
        return 1;
    }

    const std::string kind = argv[2];
    if (kind != "vat" && kind != "uat" && kind != "mat") {
        formal_usage(argv[0]);
        return 1;
    }
    try {
        (void)configured_query_route();
        (void)configured_seed_recovery();
        if (kind == "mat" || kind == "vat") (void)configured_mat_recovery();
        if (kind == "uat") (void)configured_uat_recovery();
        if (!query_only) (void)effective_query_route(configured_index_scope());
    } catch (const std::invalid_argument& error) {
        std::cerr << "Query route configuration error: " << error.what() << '\n';
        return 1;
    }
    try {
        active_memory_guard() = configure_memory_guard();
    } catch (const std::invalid_argument& error) {
        std::cerr << "Memory guard configuration error: " << error.what() << '\n';
        return 1;
    }
    write_memory_guard(std::cerr, active_memory_guard());
    if (!memory_guard_ready(active_memory_guard())) {
        std::cerr << "Memory guard could not be activated; refusing to run"
                  << " (set TAC_MEMORY_GUARD_PERCENT=0 to disable it).\n";
        return 1;
    }
    const std::filesystem::path data_dir = argv[3];
    TemporalBipartiteGraph graph;
    if (!graph.load_from_files((data_dir / "graph_ulabel.txt").string(),
                               (data_dir / "graph_vlabel.txt").string(),
                               (data_dir / "graph_timestamp.txt").string())) {
        return 2;
    }

    if (query_only) {
        const std::filesystem::path index_path = argv[4];
        const auto query = parse_formal_query(graph, argv, 5);
        if (!query) {
            formal_usage(argv[0]);
            return 1;
        }
        if (kind == "vat") {
            return run_formal_loaded_query<VATIndex>(
                "vat", index_path, graph, *query,
                [](const VATIndex& index,
                   const TemporalBipartiteGraph& current_graph,
                   const Query& current_query) {
                    return vat_query(index, current_graph, current_query);
                });
        }
        if (kind == "uat") {
            return run_formal_loaded_query<UATIndex>(
                "uat", index_path, graph, *query,
                [](const UATIndex& index,
                   const TemporalBipartiteGraph& current_graph,
                   const Query& current_query) {
                    return uat_query(index, current_graph, current_query);
                });
        }
        return run_formal_loaded_query<MATIndex>(
            "mat", index_path, graph, *query,
            [](const MATIndex& index,
               const TemporalBipartiteGraph& current_graph,
               const Query& current_query) {
                return mat_query(index, current_graph, current_query);
            });
    }

    BuildMethod method;
    try {
        method = configured_build_method();
    } catch (const std::invalid_argument& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    const bool check500 = formal_500check_enabled();
    if (check500 && argc != 4) {
        std::cerr << "TAC_INDEX_500CHECK=1 cannot be combined with a single query.\n";
        return 1;
    }

    std::optional<Query> query;
    if (argc == 10) {
        query = parse_formal_query(graph, argv, 4);
        if (!query) {
            formal_usage(argv[0]);
            return 1;
        }
    }

    if (kind == "vat") {
        if (check500) {
            return run_formal_500check<VATIndex>(
                "vat", data_dir, graph, method, build_vat_index,
                [](const VATIndex& index,
                   const TemporalBipartiteGraph& current_graph,
                   const Query& current_query,
                   QueryProfile* profile) {
                    return vat_query(
                        index, current_graph, current_query, profile);
                });
        }
        return run_formal_build<VATIndex>(
            "vat", data_dir, graph, method, query, build_vat_index,
            [](const VATIndex& index,
               const TemporalBipartiteGraph& current_graph,
               const Query& current_query) {
                return vat_query(index, current_graph, current_query);
            });
    }
    if (kind == "uat") {
        if (check500) {
            return run_formal_500check<UATIndex>(
                "uat", data_dir, graph, method, build_uat_index,
                [](const UATIndex& index,
                   const TemporalBipartiteGraph& current_graph,
                   const Query& current_query,
                   QueryProfile* profile) {
                    return uat_query(
                        index, current_graph, current_query, profile);
                });
        }
        return run_formal_build<UATIndex>(
            "uat", data_dir, graph, method, query, build_uat_index,
            [](const UATIndex& index,
               const TemporalBipartiteGraph& current_graph,
               const Query& current_query) {
                return uat_query(index, current_graph, current_query);
            });
    }
    if (check500) {
        return run_formal_500check<MATIndex>(
            "mat", data_dir, graph, method, build_mat_index,
            [](const MATIndex& index,
               const TemporalBipartiteGraph& current_graph,
               const Query& current_query,
               QueryProfile* profile) {
                return mat_query(
                    index, current_graph, current_query, profile);
            });
    }
    return run_formal_build<MATIndex>(
        "mat", data_dir, graph, method, query, build_mat_index,
        [](const MATIndex& index,
           const TemporalBipartiteGraph& current_graph,
           const Query& current_query) {
            return mat_query(index, current_graph, current_query);
        });
}

}
