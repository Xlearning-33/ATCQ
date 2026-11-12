#pragma once

#include "common/graph.h"
#include "common/types.h"
#include "index/mat_index.h"
#include "query/query_profile.h"
#include "query/query_execution.h"

namespace tacc {
QueryResult mat_query(const MATIndex& index,
                      const TemporalBipartiteGraph& graph,
                      const Query& query,
                      QueryProfile* profile = nullptr,
                      const EncodedAttributes* attributes = nullptr,
                      QueryTimer* timer = nullptr);
}
