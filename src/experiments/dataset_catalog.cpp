#include "experiments/dataset_catalog.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace tacc {

const std::vector<DatasetEntry>& dataset_catalog() {
    static const std::vector<DatasetEntry> entries{
        {"00", "SMALL", "small", "data/small"},
        {"01", "BR", "01_br", "data/real_span/01_br"},
        {"02", "SG", "02_sg", "data/real_span/02_sg"},
        {"03", "SV", "03_sv", "data/real_span/03_sv"},
        {"04", "MU", "04_mu", "data/real_span/04_mu"},
        {"05", "TA", "05_ta", "data/real_span/05_ta"},
        {"06", "NA", "06_na", "data/real_span/06_na"},
        {"07", "SO", "07_so", "data/real_span/07_so"},
        {"08", "WU", "08_wu", "data/real_span/08_wu"},
        {"09", "LB", "09_lb", "data/real_span/09_lb"},
        {"10", "LS", "10_ls", "data/real_span/10_ls"},
    };
    return entries;
}

const DatasetEntry& find_dataset(const std::string& selector) {
    std::string normalized = selector;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (normalized.size() == 1 && std::isdigit(static_cast<unsigned char>(normalized[0])))
        normalized.insert(normalized.begin(), '0');

    for (const auto& entry : dataset_catalog()) {
        std::string abbreviation = entry.abbreviation;
        std::transform(abbreviation.begin(), abbreviation.end(), abbreviation.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

        if (normalized == entry.id || normalized == entry.name || normalized == abbreviation)
            return entry;
    }

    throw std::invalid_argument("unknown dataset: " + selector + " (see exp.md)");
}

}
