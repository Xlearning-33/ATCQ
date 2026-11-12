#pragma once

#include "experiments/workload.h"

#include <array>
#include <fstream>

namespace tacc::experiments {

struct BaselineGroup {
    std::size_t count{};
    double total_us{};
};

struct BaselineAnswerCounts {
    std::size_t u_count{}, v_count{};
};

struct OnlineBaseline {
    std::filesystem::path queries_path, answers_path, samples_path, complete_path;
    std::vector<Query> queries;
    std::vector<BaselineAnswerCounts> answer_counts;
    std::array<BaselineGroup, 3> groups{};
    std::uint64_t answers_hash{};
    bool freshly_measured{};
};

OnlineBaseline create_online_baseline(const std::filesystem::path& queries_path,
    const TemporalBipartiteGraph& graph, const Parameters& parameters,
    const std::vector<Query>& prepared_queries = {});

OnlineBaseline load_online_baseline(const std::filesystem::path& queries_path,
    const TemporalBipartiteGraph& graph, const Parameters& parameters);

class BaselineAnswerReader {
public:
    BaselineAnswerReader(const OnlineBaseline& baseline,
                         const TemporalBipartiteGraph& graph);
    QueryResult next();
    void finish();
private:
    bool read_line(std::string& line);
    const OnlineBaseline& baseline_;
    std::ifstream input_;
    std::size_t position_{};
    std::uint64_t answers_hash_{14695981039346656037ULL};
};

}
