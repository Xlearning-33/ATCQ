#include "query/uat_query.h"
#include "query/index_query_utils.h"

namespace tacc {
QueryResult uat_query(const UATIndex& index,
                      const TemporalBipartiteGraph& graph,
                      const Query& query,
                      QueryProfile* profile) {
    return query_uat_block_index(index, graph, query, profile);
}
}
