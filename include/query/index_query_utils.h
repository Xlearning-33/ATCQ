#pragma once

#include "common/graph.h"
#include "index/mat_index.h"
#include "index/uat_index.h"
#include "index/vat_index.h"
#include "query/query_profile.h"
#include "query/query_execution.h"
#include <string>

namespace tacc {

enum class QueryRoute { UV, Auto };
enum class SeedRecovery { Adaptive, Seed };
enum class UATRecovery { Adaptive, Grouped };

SeedRecovery configured_seed_recovery();
const char* seed_recovery_name(SeedRecovery recovery);

UATRecovery configured_uat_recovery();
const char* uat_recovery_name(UATRecovery recovery);
bool index_uses_grouped_recovery(const std::string& kind);

QueryRoute configured_query_route();
const char* query_route_name(QueryRoute route);
QueryRoute effective_query_route(IndexScope scope);

QueryResult query_membership_index(
    const VATIndex& index,
    const TemporalBipartiteGraph& graph, const Query& query,
    QueryProfile* profile = nullptr);

QueryResult query_uat_block_index(const UATIndex& index,
                                  const TemporalBipartiteGraph& graph,
                                  const Query& query,
                                  QueryProfile* profile = nullptr);

QueryResult query_mat_block_index(const MATIndex& index,
                                  const TemporalBipartiteGraph& graph,
                                  const Query& query,
                                  QueryProfile* profile = nullptr,
                                  const EncodedAttributes* attributes = nullptr,
                                  QueryTimer* timer = nullptr);

}
