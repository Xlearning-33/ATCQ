#pragma once

#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>

namespace tacc {

struct MemoryGuardInfo {
    int percent{95};
    std::uint64_t memtotal_bytes{};
    std::uint64_t cgroup_bytes{};
    std::uint64_t usable_bytes{};
    std::uint64_t requested_bytes{};
    std::uint64_t effective_bytes{};
    int system_error{};
    std::string status{"unconfigured"};
};

inline std::uint64_t read_memory_guard_limit(const char* path) {
    std::ifstream input(path);
    std::string text;
    if (!(input >> text) || text == "max") return 0;
    try {
        std::size_t end = 0;
        const auto value = std::stoull(text, &end);
        return end == text.size() && value > 0 &&
                       value < (std::uint64_t{1} << 60)
            ? static_cast<std::uint64_t>(value) : 0;
    } catch (const std::exception&) {
        return 0;
    }
}

inline std::uint64_t read_memory_guard_memtotal() {
    std::ifstream input("/proc/meminfo");
    std::string key;
    std::uint64_t value_kb = 0;
    std::string unit;
    while (input >> key >> value_kb >> unit) {
        if (key == "MemTotal:") return value_kb * 1024ULL;
    }
    return 0;
}

inline MemoryGuardInfo configure_memory_guard() {
    MemoryGuardInfo info;
    if (const char* value = std::getenv("TAC_MEMORY_GUARD_PERCENT")) {
        const std::string text(value);
        std::size_t end = 0;
        int parsed = -1;
        try {
            parsed = std::stoi(text, &end);
        } catch (const std::exception&) {
            throw std::invalid_argument(
                "TAC_MEMORY_GUARD_PERCENT must be 0 or an integer from 1 to 95");
        }
        if (end != text.size() || parsed < 0 || parsed > 95) {
            throw std::invalid_argument(
                "TAC_MEMORY_GUARD_PERCENT must be 0 or an integer from 1 to 95");
        }
        info.percent = parsed;
    }
    if (info.percent == 0) {
        info.status = "disabled";
        return info;
    }

    info.memtotal_bytes = read_memory_guard_memtotal();
    info.cgroup_bytes = read_memory_guard_limit("/sys/fs/cgroup/memory.max");
    if (!info.cgroup_bytes) {
        info.cgroup_bytes = read_memory_guard_limit(
            "/sys/fs/cgroup/memory/memory.limit_in_bytes");
    }
    if (info.memtotal_bytes && info.cgroup_bytes) {
        info.usable_bytes = info.memtotal_bytes < info.cgroup_bytes
            ? info.memtotal_bytes : info.cgroup_bytes;
    } else {
        info.usable_bytes = info.memtotal_bytes
            ? info.memtotal_bytes : info.cgroup_bytes;
    }
    if (!info.usable_bytes) {
        info.status = "no_memory_limit_detected";
        return info;
    }
    info.requested_bytes =
        (info.usable_bytes / 100) * info.percent +
        ((info.usable_bytes % 100) * info.percent) / 100;

    struct rlimit limit {};
    if (getrlimit(RLIMIT_AS, &limit) != 0) {
        info.status = "getrlimit_failed";
        info.system_error = errno;
        return info;
    }
    if (limit.rlim_cur != RLIM_INFINITY &&
        static_cast<std::uint64_t>(limit.rlim_cur) <= info.requested_bytes) {
        info.status = "kept_existing_lower_limit";
        info.effective_bytes = static_cast<std::uint64_t>(limit.rlim_cur);
        return info;
    }
    if (info.requested_bytes >
        static_cast<std::uint64_t>(std::numeric_limits<rlim_t>::max())) {
        info.status = "limit_not_representable";
        return info;
    }
    limit.rlim_cur = static_cast<rlim_t>(info.requested_bytes);
    if (setrlimit(RLIMIT_AS, &limit) != 0) {
        info.status = "setrlimit_failed";
        info.system_error = errno;
        return info;
    }
    info.status = "set";
    info.effective_bytes = info.requested_bytes;
    return info;
}

inline bool memory_guard_ready(const MemoryGuardInfo& info) {
    return info.status == "set" ||
           info.status == "kept_existing_lower_limit" ||
           info.status == "disabled";
}

inline void write_memory_guard(std::ostream& output,
                               const MemoryGuardInfo& info) {
    output << "PROGRESS phase=memory_guard"
           << " status=" << info.status
           << " guard_percent=" << info.percent
           << " memtotal_bytes=" << info.memtotal_bytes
           << " cgroup_limit_bytes=" << info.cgroup_bytes
           << " usable_bytes=" << info.usable_bytes
           << " requested_limit_bytes=" << info.requested_bytes
           << " effective_limit_bytes=" << info.effective_bytes
           << " errno=" << info.system_error << '\n';
}

}
