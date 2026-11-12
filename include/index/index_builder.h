#pragma once

#include "common/graph.h"
#include "index/index_types.h"
#include "index/mat_index.h"
#include "index/uat_index.h"
#include "index/vat_index.h"

#include <functional>

namespace tacc {

struct IndexBuildProgress {
    std::size_t completed_attribute_sets{};
    std::size_t total_attribute_sets{};
    Side side{Side::U};
    AttributeMask current_attribute_mask{};
    double elapsed_seconds{};
    bool joint_attributes{};
    AttributeMask current_v_attribute_mask{};
};

using IndexBuildProgressCallback =
    std::function<void(const IndexBuildProgress&)>;

VATIndex build_vat_index(
    const TemporalBipartiteGraph& graph,
    BuildMethod method = BuildMethod::Opt,
    const IndexBuildProgressCallback& progress = {});
UATIndex build_uat_index(
    const TemporalBipartiteGraph& graph,
    BuildMethod method = BuildMethod::Opt,
    const IndexBuildProgressCallback& progress = {});
MATIndex build_mat_index(
    const TemporalBipartiteGraph& graph,
    BuildMethod method = BuildMethod::Opt,
    const IndexBuildProgressCallback& progress = {});

BuildMethod configured_build_method();
const char* build_method_name(BuildMethod method);
IndexScope configured_index_scope();

}

