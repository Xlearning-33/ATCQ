#pragma once

#include "common/graph.h"
#include "query/query_profile.h"
#include "query/query_execution.h"
#include "common/types.h"

namespace tacc {

QueryResult online_query(const TemporalBipartiteGraph& graph,
                         const Query& query,
                         QueryProfile* profile = nullptr,
                         QueryTimer* timer = nullptr);

}
