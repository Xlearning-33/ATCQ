#include "query/vat_query.h"
#include "query/index_query_utils.h"

namespace tacc {
QueryResult vat_query(const VATIndex& index,
                      const TemporalBipartiteGraph& graph,
                      const Query& query,
                      QueryProfile* profile) {
    return query_membership_index(index, graph, query, profile);
}
}
