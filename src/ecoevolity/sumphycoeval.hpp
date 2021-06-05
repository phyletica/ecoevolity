/******************************************************************************
 * Copyright (C) 2015-2016 Jamie R. Oaks.
 *
 * This file is part of Ecoevolity.
 *
 * Ecoevolity is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Ecoevolity is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Ecoevolity.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

#ifndef SUMPHYCOEVAL_HPP
#define SUMPHYCOEVAL_HPP

#include <limits>
#include <time.h>

#include "cpp-optparse/OptionParser.h"

#include "version.hpp"
#include "error.hpp"
#include "path.hpp"
#include "node.hpp"
#include "treesum.hpp"


void write_sum_phy_splash(std::ostream& out);

void check_sumphy_output_path(const std::string& path);

template <class NodeType>
int sumphycoeval_main(int argc, char * argv[]) {

    write_sum_phy_splash(std::cerr);
    std::cerr << "\n";

    const std::string usage = 
        "usage: %prog [OPTIONS] PHYCOEVAL-TREE-LOG-FILE-1 [PHYCOEVAL-TREE-LOG-FILE-2 ...]";

    std::ostringstream version_ss;
    version_ss << PROJECT_NAME << " version " << PROJECT_DETAILED_VERSION;
    const std::string version = version_ss.str();

    const std::string description =
        "Sumphycoeval: Summarizing phylogenetic coevality";
    // const std::string epilog =
    //     "Epilog goes here...";

    optparse::OptionParser parser = optparse::OptionParser()
            .usage(usage)
            .version(version)
            .description(description);
            // .epilog(epilog);

    parser.add_option("-b", "--burnin")
            .action("store")
            .type("unsigned int")
            .dest("burnin")
            .set_default("0")
            .help("Number of samples from the beginning of each tree log file "
                  "to ignore as burn in. "
                  "Default: 0.");
    parser.add_option("-t", "--target-tree")
            .action("store")
            .dest("target_tree_path")
            .set_default("")
            .help("Path to a file containing a tree on to which the MCMC "
                  "samples of trees will be summarized.");
    parser.add_option("--to", "--target-tree-out")
            .action("store")
            .dest("target_tree_out_path")
            .set_default("")
            .help("Path to a file where the annotated target tree "
                  "will be written in nexus format. Default: Do not write a "
                  "nexus-formatted file of the annotated target tree.");
    parser.add_option("--mo", "--map-tree-out")
            .action("store")
            .dest("map_tree_out_path")
            .set_default("")
            .help("Path to a file where the maximum a posteriori (MAP) tree(s) "
                  "will be written in nexus format. Default: Do not write a "
                  "nexus-formatted file of the MAP tree(s).");
    parser.add_option("--median-heights")
            .action("store_true")
            .dest("use_median_heights")
            .help("Use median (rather than mean) of MCMC samples for the node "
                  "heights of output MAP trees. Default: Use mean.");
    parser.add_option("--min-split-freq")
            .action("store")
            .type("double")
            .dest("min_split_freq")
            .set_default("0.1")
            .help("Minimum frequency of splits for calculating the average "
                  "standard deviation of split frequencies among MCMC "
                  "samples. Default: 0.1. This is only used if more than "
                  "one tree log file is provided.");
    parser.add_option("-m" "--multiplier")
            .action("store")
            .type("double")
            .dest("multiplier")
            .set_default("-1.0")
            .help("Scale all trees by multiplying their branch lengths by this "
                  "number. Default: Do not rescale trees.");
    parser.add_option("-n", "--newick-target")
            .action("store_true")
            .dest("newick_target_tree")
            .help("The target tree is in newick format. By default, the target "
                  "tree must be in nexus format.");
    parser.add_option("--include-merged-target-heights")
            .action("store_true")
            .dest("include_merged_target_heights")
            .help("Include a summary of merged heights from the target tree. "
                  "If a target tree is not provided, this option is ignored.");
    parser.add_option("-f", "--force")
            .action("store_true")
            .dest("force")
            .help("Force overwriting of existing output files.");

    optparse::Values& options = parser.parse_args(argc, argv);
    std::vector<std::string> log_paths = parser.args();

    if (log_paths.size() < 1) {
        throw EcoevolityError("Too few positional arguments; "
                "a path to at least one tree log file is required.");
    }

    for (unsigned int i = 0; i < log_paths.size(); ++i) {
        std::string log_path = log_paths.at(i);
        if (! path::exists(log_path)) {
            throw EcoevolityError("Log path \'" + log_path +
                    "\' does not exist");
        }
        if (! path::isfile(log_path)) {
            throw EcoevolityError("Log path \'" + log_path +
                    "\' is not a regular file");
        }
    }

    // Parse options
    const bool prevent_overwrite = (! options.get("force"));
    unsigned int burnin = options.get("burnin");

    std::string target_tree_path;
    bool target_tree_provided = false;
    if (options.is_set_by_user("target_tree_path")) {
        target_tree_path = options.get("target_tree_path").get_str();
        target_tree_provided = true;
        if (! path::exists(target_tree_path)) {
            throw EcoevolityError("Target tree path \'" + target_tree_path +
                    "\' does not exist");
        }
        if (! path::isfile(target_tree_path)) {
            throw EcoevolityError("Target tree path \'" + target_tree_path +
                    "\' is not a regular file");
        }
    }

    std::string target_tree_out_path;
    bool writing_target_to_nexus = false;
    if (options.is_set_by_user("target_tree_out_path")) {
        target_tree_out_path = options.get("target_tree_out_path").get_str();
        writing_target_to_nexus = true;
        if (prevent_overwrite && path::exists(target_tree_out_path)) {
            throw EcoevolityError("Target tree output path \'" +
                    target_tree_out_path +
                    "\' already exists. Please specify a different path or use "
                    "the \'--force\' option to overwrite the file.");
        }
    }

    if (writing_target_to_nexus && (! target_tree_provided)) {
        throw EcoevolityError(
                "Target tree output path was specified, but a target tree was "
                "not provided.");
    }

    std::string map_tree_out_path;
    bool writing_map_to_nexus = false;
    if (options.is_set_by_user("map_tree_out_path")) {
        map_tree_out_path = options.get("map_tree_out_path").get_str();
        writing_map_to_nexus = true;
        if (prevent_overwrite && path::exists(map_tree_out_path)) {
            throw EcoevolityError("MAP tree output path \'" +
                    map_tree_out_path +
                    "\' already exists. Please specify a different path or use "
                    "the \'--force\' option to overwrite the file.");
        }
    }

    std::string target_tree_format = "nexus";
    if (options.get("newick_target_tree")) {
        target_tree_format = "relaxedphyliptree";
    }

    const bool use_median_heights = options.get("use_median_heights");
    const double multiplier = options.get("multiplier");
    const double min_split_freq = options.get("min_split_freq");
    if ((min_split_freq < 0.0) || (min_split_freq >=1.0)) {
        throw EcoevolityError("\'--min-split-freq\' must be between 0 and 1\n");
    }
    const bool include_merged_target_heights = (
            options.get("include_merged_target_heights") && target_tree_provided);
    const double precision = 18;

    std::ofstream target_tree_out_stream;
    if (writing_target_to_nexus) {
        target_tree_out_stream.open(target_tree_out_path);
        if (! target_tree_out_stream.is_open()) {
            std::ostringstream message;
            message << "ERROR: Could not open target-tree output file \'"
                    << target_tree_out_path
                    << "\'\n";
            throw EcoevolityError(message.str());
        }
    }
    std::ofstream map_tree_out_stream;
    if (writing_map_to_nexus) {
        map_tree_out_stream.open(map_tree_out_path);
        if (! map_tree_out_stream.is_open()) {
            std::ostringstream message;
            message << "ERROR: Could not open MAP tree output file \'"
                    << map_tree_out_path
                    << "\'\n";
            throw EcoevolityError(message.str());
        }
    }

    const double ultrametricity_tolerance = 1e-6;

    time_t start;
    time_t finish;
    time(&start);

    std::cerr << "Parsing trees from files..." << std::endl;
    treesum::TreeSample<PopulationNode> tree_sample;
    if (target_tree_provided) {
        tree_sample = treesum::TreeSample<PopulationNode>(
                target_tree_path,
                log_paths,
                target_tree_format,
                "nexus",
                burnin,
                ultrametricity_tolerance,
                multiplier);
    }
    else {
        tree_sample = treesum::TreeSample<PopulationNode>(
                log_paths,
                "nexus",
                burnin,
                ultrametricity_tolerance,
                multiplier);
    }

    if (writing_target_to_nexus) {
        std::cerr << "Writing annotated target tree to:\n"
                  << "  " << target_tree_out_path << std::endl;
        tree_sample.write_target_tree_to_nexus(
                target_tree_out_stream,
                precision);
        target_tree_out_stream.close();
    }
    if (writing_map_to_nexus) {
        std::cerr << "Writing annotated MAP trees to:\n"
                  << "  " << map_tree_out_path << std::endl;
        tree_sample.write_map_trees_to_nexus(
                map_tree_out_stream,
                use_median_heights,
                precision);
        map_tree_out_stream.close();
    }

    std::cerr << "Writing YAML-formatted summary to standard output..." << std::endl;
    tree_sample.write_summary(std::cout,
            use_median_heights,
            min_split_freq,
            "",
            precision);
    if (include_merged_target_heights) {
        tree_sample.write_summary_of_merged_target_heights(std::cout, "", precision);
    }

    time(&finish);
    double duration = difftime(finish, start);
    std::cerr << "Runtime: " << duration << " seconds." << std::endl;

    return 0;
}

#endif
