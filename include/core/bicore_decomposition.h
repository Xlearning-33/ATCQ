#pragma once

#include "common/graph.h"
#include "common/types.h"

#include <vector>

namespace tacc {

QueryResult bicore_peeling(const BipartiteSubgraph& graph,
                           Count alpha,
                           Count beta);

QueryResult bicore_peeling_candidates(
    BipartiteSubgraph& graph,
    Count alpha,
    Count beta,
    const std::vector<VertexID>& candidate_u,
    const std::vector<VertexID>& candidate_v);

struct BicoreSkyline {
    Count max_alpha{};
    Count max_beta{};
    std::vector<std::vector<Count>> u_beta_max;
    std::vector<std::vector<Count>> v_alpha_max;
};

BicoreSkyline bicore_skyline_decomposition(const BipartiteSubgraph& graph);

}
