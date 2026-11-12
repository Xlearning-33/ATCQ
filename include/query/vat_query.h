#pragma once

#include "common/graph.h"
#include "common/types.h"
#include "index/vat_index.h"
#include "query/query_profile.h"

namespace tacc {
QueryResult vat_query(const VATIndex& index,
                      const TemporalBipartiteGraph& graph,
                      const Query& query,
                      QueryProfile* profile = nullptr);
}
