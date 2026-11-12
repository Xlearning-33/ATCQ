#pragma once

#include "common/graph.h"
#include "query/query_profile.h"
#include "query/query_execution.h"

namespace tacc {

enum class MATRecovery { Grouped, Legacy };
MATRecovery configured_mat_recovery();
const char* mat_recovery_name(MATRecovery recovery);

QueryResult query_mat_grouped_seeds(const TemporalBipartiteGraph& graph,
                                   const Query& query,
                                   std::vector<VertexID> seeds,
                                   AttributeMask u_mask,
                                   QueryProfile* profile,
                                   QueryTimer* timer = nullptr);

}
