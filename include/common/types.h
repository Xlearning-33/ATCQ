#pragma once

#include <cstdint>
#include <vector>

namespace tacc {

using VertexID = std::uint32_t;
using Timestamp = std::uint32_t;
using AttributeID = std::uint32_t;
using AttributeMask = std::uint64_t;
using Count = std::uint32_t;

struct TemporalEdge {
    VertexID u{};
    VertexID v{};
    Timestamp timestamp{};
};

struct TimeWindow {
    Timestamp start{};
    Timestamp end{};
};

struct Query {
    TimeWindow window;
    std::vector<AttributeID> u_attributes;
    std::vector<AttributeID> v_attributes;
    std::uint32_t alpha{};
    std::uint32_t beta{};
};

struct QueryResult {
    std::vector<VertexID> u_vertices;
    std::vector<VertexID> v_vertices;
};

}
