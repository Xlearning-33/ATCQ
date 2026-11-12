#pragma once

#include "common/types.h"

#include <chrono>
#include <stdexcept>

namespace tacc {

struct EncodedAttributes { AttributeMask u{}, v{}; };

inline EncodedAttributes encode_query_attributes(const Query& query) {
    const auto encode = [](const std::vector<AttributeID>& attributes) {
        AttributeMask mask = 0;
        for (auto id : attributes) {
            if (id >= 64) throw std::invalid_argument("query attribute ID exceeds 63");
            mask |= AttributeMask{1} << id;
        }
        return mask;
    };
    return {encode(query.u_attributes), encode(query.v_attributes)};
}

inline constexpr const char* query_timing_model = "answer_ready_v2";
inline constexpr const char* batch_query_timing_model = "batch_500_v3";

struct TimedQueryBatch {
    std::vector<QueryResult> answers;
    double total_us{};
};

template<class Clock = std::chrono::steady_clock, class Execute>
TimedQueryBatch run_timed_query_batch(std::size_t count, Execute&& execute) {
    TimedQueryBatch batch;
    batch.answers.resize(count);
    const auto start = Clock::now();
    for (std::size_t i = 0; i < count; ++i) batch.answers[i] = execute(i);
    const auto end = Clock::now();
    batch.total_us = std::chrono::duration<double, std::micro>(end - start).count();
    return batch;
}

class QueryTimer {
public:
    using Clock = std::chrono::steady_clock;
    QueryTimer() : start_(Clock::now()) {}
    void stop() {
        if (!stopped_) { end_ = Clock::now(); stopped_ = true; }
    }
    bool stopped() const { return stopped_; }
    double microseconds() const {
        return std::chrono::duration<double, std::micro>(end_ - start_).count();
    }
private:
    Clock::time_point start_, end_{};
    bool stopped_{};
};

inline void stop_query_timer(QueryTimer* timer) { if (timer) timer->stop(); }

}
