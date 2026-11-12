#pragma once

#include <string>
#include <vector>

namespace tacc {

struct DatasetEntry {
    const char* id;
    const char* abbreviation;
    const char* name;
    const char* relative_directory;
};

const std::vector<DatasetEntry>& dataset_catalog();
const DatasetEntry& find_dataset(const std::string& selector);

}
