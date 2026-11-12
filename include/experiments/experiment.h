#pragma once

#include "index/index_builder.h"
#include <iosfwd>

namespace tacc {

IndexBuildProgressCallback make_index_build_progress_reporter(
    const char* kind, BuildMethod method, std::ostream& log);

int run_experiment(int argc, char** argv);
int run_index_experiment(int argc, char** argv);
int run_quota_experiment(int argc, char** argv);

}
