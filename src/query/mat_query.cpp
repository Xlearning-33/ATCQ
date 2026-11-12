#include "query/mat_query.h"
#include "query/index_query_utils.h"

namespace tacc {
QueryResult mat_query(const MATIndex& index,
                      const TemporalBipartiteGraph& graph,
                      const Query& query,
                      QueryProfile* profile,
                      const EncodedAttributes* attributes,
                      QueryTimer* timer) {
    return query_mat_block_index(index, graph, query, profile, attributes, timer);
}
}
