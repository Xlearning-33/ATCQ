#pragma once

#include "common/graph.h"
#include "common/types.h"
#include "index/uat_index.h"
#include "query/query_profile.h"

namespace tacc {
QueryResult uat_query(const UATIndex& index,
                      const TemporalBipartiteGraph& graph,
                      const Query& query,
                      QueryProfile* profile = nullptr);
}
