#ifndef ECOEVOLITY_UTILS_FOR_TESTING
#define ECOEVOLITY_UTILS_FOR_TESTING

#include <iostream>
#include <algorithm>

#include "ecoevolity/path.hpp"
#include "ecoevolity/string_util.hpp"

inline void write_r_script(const std::vector<unsigned int> & counts,
        const std::string & path,
        const std::vector<unsigned int> num_topologies_by_ndivs = {}) {
    std::pair<std::string, std::string> prefix_ext = path::splitext(path::basename(path));
    std::ofstream os;
    os.open(path);
    os << "#! /usr/bin/env Rscript\n\n"
       << "plot_binom_on_hist <- function(counts, p) {\n"
       << "    k = seq(from = min(counts), to = max(counts), by = 1)\n"
       << "    n = sum(counts)\n"
       << "    binom_probs = dbinom(k, n, p)\n"
       << "    hist(counts, freq = F)\n"
       << "    lines(k, binom_probs, type = 'l')\n"
       << "}\n\n";

    if (num_topologies_by_ndivs.size() > 0) {
        os << "num_topologies_by_ndivs = c(";
        for (unsigned int i = 0; i < num_topologies_by_ndivs.size(); ++i) {
            if (i == 0) {
                os << num_topologies_by_ndivs.at(i);
            }
            else {
                os << ", " << num_topologies_by_ndivs.at(i);
            }
        }
        os << ")\n";
    }
    os << "counts = c(";
    for (unsigned int i = 0; i < counts.size(); ++i) {
        if (i == 0) {
            os << counts.at(i);
        }
        else {
            if ((i + 1) % 1000 == 0) {
                os << ",\n        " << counts.at(i);
            }
            else {
                os << ", " << counts.at(i);
            }
        }
    }
    os << ")\n"
       << "number_of_topologies = length(counts)\n"
       << "binomial_prob = 1.0 / number_of_topologies\n"
       << "topology = seq(from = 1, to = number_of_topologies, by = 1)\n"
       << "d = data.frame(topology = topology, count = counts)\n\n";
    std::string plot_path = prefix_ext.first + "-topo-vs-count.pdf";
    os << "pdf(\"" << plot_path << "\")\n"
       << "plot(x = d$topology, y = d$count, ylim = c(0, max(counts)))\n"
       << "dev.off()\n\n";
    plot_path = prefix_ext.first + "-count-distribution.pdf";
    os << "pdf(\"" << plot_path << "\")\n"
       << "plot_binom_on_hist(counts = counts, p = binomial_prob)\n"
       << "dev.off()\n\n"
       << "chisq.test(counts)\n";
    os.close();
}

inline void write_r_script(
        const std::map< std::set< std::set<Split> >, unsigned int> & split_counts,
        const unsigned int number_of_leaves,
        const std::string & path) {
    std::vector<unsigned int> counts;
    std::vector<unsigned int> num_topologies_by_ndivs(number_of_leaves - 1, 0);
    counts.reserve(split_counts.size());
    for (auto splitset_count : split_counts) {
        counts.push_back(splitset_count.second);
        ++num_topologies_by_ndivs[splitset_count.first.size()];
    }
    write_r_script(counts, path, num_topologies_by_ndivs);
}

inline int get_number_of_subsets_from_model_string(const std::string & model) {
    std::vector<std::string> indices_str = string_util::split(model, ',');
    std::vector<int> indices;
    indices.reserve(indices_str.size());
    for (const auto & index_str : indices_str) {
        indices.push_back(std::stoi(index_str));
    }
    auto max_iter = std::max_element(indices.begin(), indices.end());
    int max_index = *max_iter;
    return max_index + 1;
}

inline void write_sampled_models_tsv(
        const std::string & path,
        const std::map<std::string, int> & model_counts,
        const std::map<std::string, double> & expected_model_probs) {
    std::ofstream os;
    os.open(path);
    os.precision(12);
    os << "model\tnum_samples\texpected_prob" << std::endl;
    for (const auto & pair : model_counts) {
        os << pair.first
           << "\t"
           << pair.second
           << "\t"
           << expected_model_probs.at(pair.first)
           << std::endl;
    }
}

#endif
