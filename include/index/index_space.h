#pragma once

#include "common/graph.h"
#include "index/mat_index.h"
#include "index/uat_index.h"
#include "index/vat_index.h"
#include <cstdint>
#include <ostream>

namespace tacc {

inline double space_mb(std::uint64_t bytes) { return double(bytes) / 1048576.0; }
inline double space_G(std::uint64_t bytes) { return double(bytes) / 1073741824.0; }

struct IndexSpaceBreakdown {
    const char* primary_component{};
    const char* secondary_component{};
    std::uint64_t root_bytes{};
    std::uint64_t primary_bytes{};
    std::uint64_t secondary_bytes{};
    std::uint64_t directory_bytes{};
    std::uint64_t time_group_bytes{};
    std::uint64_t postings_bytes{};

    std::uint64_t body_bytes() const { return primary_bytes + secondary_bytes; }
    std::uint64_t query_aux_bytes() const {
        return directory_bytes + time_group_bytes + postings_bytes;
    }
    std::uint64_t heap_bytes() const { return body_bytes() + query_aux_bytes(); }
    std::uint64_t index_bytes() const { return root_bytes + heap_bytes(); }
};

inline IndexSpaceBreakdown formal_index_space(const VATIndex& x) {
    IndexSpaceBreakdown s;
    s.primary_component = "vat_vertex_tree"; s.secondary_component = "mai_pool";
    s.root_bytes = sizeof(x);
    s.primary_bytes = x.tree_bytes();
    s.secondary_bytes = x.intervals().size() * sizeof(VATInterval);
    s.directory_bytes = x.query_masks().size() * sizeof(VATMaskRange) +
        x.query_cores().size() * sizeof(VATCoreRange);
    s.postings_bytes = x.query_refs().size() * sizeof(VATRecordRef);
    return s;
}
inline IndexSpaceBreakdown formal_index_space(const UATIndex& x) {
    IndexSpaceBreakdown s;
    s.primary_component = "uat_core_tree"; s.secondary_component = "vertex_pool";
    s.root_bytes = sizeof(x);
    s.primary_bytes = x.tree_bytes();
    s.secondary_bytes = x.vertex_pool().size() * sizeof(VertexID);
    return s;
}
inline IndexSpaceBreakdown formal_index_space(const MATIndex& x) {
    IndexSpaceBreakdown s;
    s.primary_component = "mat_blocks"; s.secondary_component = "vertex_pool";
    s.root_bytes = sizeof(x);
    s.primary_bytes = x.blocks().size() * sizeof(IndexBlock);
    s.secondary_bytes = x.vertex_pool().size() * sizeof(VertexID);
    s.directory_bytes = x.query_directory().size() * sizeof(MATIndex::QueryDirectory::value_type);
    s.time_group_bytes = x.query_time_groups().size() * sizeof(MATIndex::TimeGroup);
    return s;
}

inline void write_formal_index_space(std::ostream& out, const char* kind,
                                     const IndexSpaceBreakdown& s,
                                     std::size_t recorded_heap_bytes) {
    out << "RESULT exp=formal_index_space index=" << kind
        << " size_model=structural_v2 unit_mb_bytes=1048576 unit_G_bytes=1073741824"
        << " primary_component=" << s.primary_component
        << " secondary_component=" << s.secondary_component
        << " index_bytes=" << s.index_bytes()
        << " index_mb=" << space_mb(s.index_bytes())
        << " index_G=" << space_G(s.index_bytes())
        << " root_bytes=" << s.root_bytes << " heap_bytes=" << s.heap_bytes()
        << " body_bytes=" << s.body_bytes() << " body_mb=" << space_mb(s.body_bytes())
        << " query_aux_bytes=" << s.query_aux_bytes()
        << " query_aux_mb=" << space_mb(s.query_aux_bytes())
        << " primary_bytes=" << s.primary_bytes << " secondary_bytes=" << s.secondary_bytes
        << " query_directory_bytes=" << s.directory_bytes
        << " query_time_group_bytes=" << s.time_group_bytes
        << " query_postings_bytes=" << s.postings_bytes
        << " recorded_approx_bytes=" << recorded_heap_bytes
        << " recorded_approx_scope=heap_without_root size_matches_recorded="
        << (s.heap_bytes() == recorded_heap_bytes ? 1 : 0) << '\n';
}

inline std::uint64_t grouped_aux_structural_bytes(const GroupedVTimeAdjacency& x) {
    return sizeof(x) + (x.v_groups.size() + x.starts.size()) * sizeof(std::size_t)
        + x.times.size() * sizeof(Timestamp) + x.neighbors.size() * sizeof(VertexID);
}

inline void write_query_index_space(std::ostream& out, const char* kind,
                                    const IndexSpaceBreakdown& s,
                                    const TemporalBipartiteGraph& graph,
                                    bool uses_grouped) {
    const std::uint64_t support = uses_grouped
        ? grouped_aux_structural_bytes(graph.grouped_v_time_adjacency()) : 0;
    out << "RESULT exp=query_index_space index=" << kind
        << " size_model=structural_v2 unit_mb_bytes=1048576 unit_G_bytes=1073741824"
        << " index_bytes=" << s.index_bytes()
        << " query_support_bytes=" << support
        << " query_support_mb=" << space_mb(support)
        << " query_support_G=" << space_G(support)
        << " index_with_query_support_bytes=" << s.index_bytes() + support
        << " index_with_query_support_mb=" << space_mb(s.index_bytes() + support)
        << " index_with_query_support_G=" << space_G(s.index_bytes() + support)
        << " input_graph_included=0 transient_workspace_included=0\n";
}
}
