#pragma once

#include "index/index_types.h"

#include <algorithm>
#include <charconv>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace tacc {

class PoolInterner {
public:
    static std::size_t configured_limit() {
        const char* value = std::getenv("TAC_POOL_INTERN_MAX_ENTRIES");
        if (value == nullptr) return 262144;
        const std::string text(value);
        std::size_t limit{};
        const auto parsed = std::from_chars(text.data(), text.data() + text.size(), limit);
        if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size() ||
            limit > std::numeric_limits<std::size_t>::max() / sizeof(Entry) / 4)
            throw std::invalid_argument("invalid TAC_POOL_INTERN_MAX_ENTRIES");
        return limit;
    }

    explicit PoolInterner(std::size_t limit = configured_limit()) : limit_(limit) {}

    template <typename T, typename Word>
    IndexPostingRange intern(std::vector<T>& pool, const std::vector<T>& list,
                             Word word) {
        std::uint64_t hash = 1469598103934665603ULL;
        for (const auto& value : list) {
            hash ^= word(value);
            hash *= 1099511628211ULL;
        }
        return intern_hashed(pool, list, word, hash);
    }

    template <typename T, typename Word>
    IndexPostingRange intern_hashed(std::vector<T>& pool, const std::vector<T>& list,
                                    Word word, std::uint64_t hash) {
        if (!heads_.empty()) {
            for (auto id = heads_[hash % heads_.size()]; id != missing;
                 id = entries_[id].next) {
                const auto& entry = entries_[id];
                if (entry.hash != hash || entry.count != list.size()) continue;
                bool equal = true;
                for (std::size_t i = 0; i < list.size(); ++i) {
                    if (word(pool[entry.offset + i]) != word(list[i])) {
                        equal = false;
                        break;
                    }
                }
                if (equal) return {entry.offset, entry.count};
            }
        }
        if (list.size() > pool.max_size() - pool.size())
            throw std::overflow_error("shared index pool exceeds max_size");
        const IndexPostingRange result{pool.size(), list.size()};
        pool.insert(pool.end(), list.begin(), list.end());
        if (entries_.size() < limit_) {
            if (entries_.size() == entries_.capacity()) {
                const auto next = entries_.capacity() == 0
                    ? std::min<std::size_t>(1024, limit_)
                    : std::min(limit_, entries_.capacity() * 2);
                entries_.reserve(next);
            }
            if (heads_.empty() || entries_.size() >= heads_.size()) {
                const auto count = heads_.empty() ? std::min<std::size_t>(1024, limit_)
                    : std::min(limit_, heads_.size() * 2);
                std::vector<std::size_t> heads(count, missing);
                for (std::size_t i = 0; i < entries_.size(); ++i) {
                    auto& entry = entries_[i];
                    const auto bucket = entry.hash % count;
                    entry.next = heads[bucket];
                    heads[bucket] = i;
                }
                heads_.swap(heads);
            }
            const auto bucket = hash % heads_.size();
            entries_.push_back({hash, result.offset, result.count, heads_[bucket]});
            heads_[bucket] = entries_.size() - 1;
        }
        return result;
    }

private:
    static constexpr std::size_t missing = std::numeric_limits<std::size_t>::max();
    struct Entry { std::uint64_t hash; std::size_t offset, count, next; };
    std::size_t limit_;
    std::vector<Entry> entries_;
    std::vector<std::size_t> heads_;
};
}
