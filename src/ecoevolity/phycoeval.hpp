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

#ifndef PHYCOEVAL_HPP
#define PHYCOEVAL_HPP

#include <limits>
#include <time.h>

#include "cpp-optparse/OptionParser.h"

#include "version.hpp"
#include "error.hpp"
#include "rng.hpp"
#include "path.hpp"
#include "string_util.hpp"
#include "general_tree_settings.hpp"
#include "settings_io.hpp"
#include "mcmc.hpp"


void write_phy_splash(std::ostream& out);

void update_log_paths(
        std::string & tree_log_path,
        std::string & state_log_path,
        std::string & operator_log_path,
        unsigned int max_number_of_attempts = 10000);

void increment_log_path(std::string & log_path);

template <class TreeType>
int phycoeval_main(int argc, char * argv[]) {

    write_phy_splash(std::cout);
    std::cout << "\n";

    const std::string usage = 
        "usage: %prog [OPTIONS] YAML-CONFIG-FILE";

    std::ostringstream version_ss;
    version_ss << PROJECT_NAME << " version " << PROJECT_DETAILED_VERSION;
    const std::string version = version_ss.str();

    const std::string description =
        "Phycoeval: Estimating phylogenetic coevality";
    // const std::string epilog =
    //     "Epilog goes here...";

    optparse::OptionParser parser = optparse::OptionParser()
            .usage(usage)
            .version(version)
            .description(description);
            // .epilog(epilog);

    parser.set_defaults("ignore_data", "false");

    parser.add_option("--seed")
            .action("store")
            .type("long")
            .dest("seed")
            .help("Seed for random number generator. Default: Set from clock.");
    parser.add_option("--ignore-data")
            .action("store_true")
            .dest("ignore_data")
            .help("Ignore data to sample from the prior distribution. Default: "
                  "Use data to sample from the posterior distribution");
#ifdef BUILD_WITH_THREADS
    parser.add_option("--nthreads")
            .action("store")
            .type("unsigned int")
            .dest("nthreads")
            .set_default("1")
            .help("Number of threads to use for likelihood calculations. "
                  "Default: 1 (no multithreading). If you are using "
                  "the \'--ignore-data\' option, no likelihood calculations "
                  "will be performed, and so no multithreading is used.");
#endif
    parser.add_option("--prefix")
            .action("store")
            .dest("prefix")
            .set_default("")
            .help("Optional string to prefix all output files.");
    parser.add_option("--relax-constant-sites")
            .action("store_true")
            .dest("relax_constant_sites")
            .help("By default, if you specify \'constant_sites_removed = true\' "
                  "and constant sites are found, phycoeval throws an error. "
                  "With this option, phycoeval will automatically ignore the "
                  "constant sites and only issue a warning (and correct for "
                  "constant sites in the likelihood calculation). Please make sure "
                  "you understand what you are doing when you use this option."
                );
    parser.add_option("--relax-missing-sites")
            .action("store_true")
            .dest("relax_missing_sites")
            .help("By default, if a column is found for which there is no data "
                  "for at least one population, phycoeval throws an error. "
                  "With this option, phycoeval will automatically ignore such "
                  "sites and only issue a warning."
                );
    parser.add_option("--relax-triallelic-sites")
            .action("store_true")
            .dest("relax_triallelic_sites")
            .help("By default, if a DNA site is found for which there is more "
                  "than two nucleotide states, phycoeval throws an error. "
                  "With this option, phycoeval will automatically recode such "
                  "sites as biallelic and only issue a warning. These sites "
                  "are recoded by assigning state 0 to the first nucleotide "
                  "found and state 1 to all others. If you do not wish to "
                  "recode such sites and prefer to ignore them, please remove "
                  "all sites with more than two nucleotide states from your "
                  "DNA alignments. NOTE: only alignments of nucleotides are "
                  "affected by this option, not alignments of standard "
                  "characters (i.e., 0, 1, 2)."
                );
    parser.add_option("--dry-run")
            .action("store_true")
            .dest("dry_run")
            .help("Do not run analysis; only process and report settings.");

    optparse::Values& options = parser.parse_args(argc, argv);
    std::vector<std::string> args = parser.args();

    RandomNumberGenerator rng;
    long seed_opt;
    if (options.is_set_by_user("seed")) {
        seed_opt = options.get("seed");
    }
    else {
        seed_opt = rng.uniform_int(1, std::numeric_limits<int>::max() - 1);
    }
    const long seed = seed_opt;
    rng.set_seed(seed);
    std::cout << "Seed: " << seed << std::endl;

    const bool dry_run = options.get("dry_run");

    const bool ignore_data = options.get("ignore_data");
    if (ignore_data) {
        std::cout << "Ignoring data in order to sample from the prior distribution..." << std::endl;
    }
    else {
        std::cout << "Using data in order to sample from the posterior distribution..." << std::endl;
    }

    const bool strict_on_constant_sites = (! options.get("relax_constant_sites"));
    const bool strict_on_missing_sites = (! options.get("relax_missing_sites"));
    const bool strict_on_triallelic_sites = (! options.get("relax_triallelic_sites"));

#ifdef BUILD_WITH_THREADS 
    unsigned int nthreads = options.get("nthreads");
#else
    unsigned int nthreads = 1;
#endif

    if (args.size() < 1) {
        throw EcoevolityError("Path to YAML-formatted config file is required");
    }
    if (args.size() > 1) {
        throw EcoevolityError("Too many arguments; only one config file is allowed");
    }
    const std::string config_path = args.at(0);
    if (! path::exists(config_path)) {
        throw EcoevolityError("Config file \'" + config_path +
                "\' does not exist");
    }
    if (! path::isfile(config_path)) {
        throw EcoevolityError("Config path \'" + config_path +
                "\' is not a regular file");
    }
    std::cout << "Config path: " << config_path << std::endl;

    std::cout << "Parsing config file..." << std::endl;
    PopulationTreeSettings settings(config_path);

    std::cout << "Configuring model..." << std::endl;

    TreeType tree(
            settings,
            rng,
            strict_on_constant_sites,
            strict_on_missing_sites,
            strict_on_triallelic_sites,
            false // store_seq_loci_info
            );

    GeneralTreeOperatorSchedule<BasePopulationTree> operator_schedule(
            settings.operator_settings, tree.get_leaf_node_count());

    std::cout << "\n" << string_util::banner('-') << "\n";
    write_settings(std::cout, settings, operator_schedule);
    std::cout << string_util::banner('-') << "\n\n";

    unsigned int n_moves_per_generation = tree.get_leaf_node_count();

    if (ignore_data) {
        tree.ignore_data();
    }
    else {
        tree.use_data();
    }

    unsigned int logging_precision = 18;

    std::string tree_log_path = settings.get_tree_log_path();
    std::string state_log_path = settings.get_state_log_path();
    std::string operator_log_path = settings.get_operator_log_path();
    if (options.is_set_by_user("prefix")) {
        std::string output_prefix = options.get("prefix").get_str();
        tree_log_path = output_prefix + path::basename(tree_log_path);
        state_log_path = output_prefix + path::basename(state_log_path);
        operator_log_path = output_prefix + path::basename(operator_log_path);
    }
    update_log_paths(tree_log_path, state_log_path, operator_log_path);

    std::cout << "\n" << string_util::banner('-') << "\n";
    tree.write_data_summary(std::cout);
    std::cout << string_util::banner('-') << "\n\n";

    std::cout << "Number of threads: " << nthreads << std::endl;

    if (dry_run) {
        return 0;
    }

    if (path::exists(tree_log_path)) {
        std::ostringstream message;
        message << "ERROR: The tree log file \'"
                << tree_log_path
                << "\' already exists!\n";
        throw EcoevolityError(message.str());
    }
    if (path::exists(state_log_path)) {
        std::ostringstream message;
        message << "ERROR: The state log file \'"
                << state_log_path
                << "\' already exists!\n";
        throw EcoevolityError(message.str());
    }
    if (path::exists(operator_log_path)) {
        std::ostringstream message;
        message << "ERROR: The operator log file \'"
                << operator_log_path
                << "\' already exists!\n";
        throw EcoevolityError(message.str());
    }

    std::ofstream tree_log_stream;
    std::ofstream state_log_stream;
    std::ofstream operator_log_stream;

    tree_log_stream.open(tree_log_path);
    state_log_stream.open(state_log_path);
    operator_log_stream.open(operator_log_path);
    
    if (! tree_log_stream.is_open()) {
        std::ostringstream message;
        message << "ERROR: Could not open tree log file \'"
                << tree_log_path
                << "\'\n";
        throw EcoevolityError(message.str());
    }
    if (! state_log_stream.is_open()) {
        std::ostringstream message;
        message << "ERROR: Could not open state log file \'"
                << state_log_path
                << "\'\n";
        throw EcoevolityError(message.str());
    }
    if (! operator_log_stream.is_open()) {
        std::ostringstream message;
        message << "ERROR: Could not open operator log file \'"
                << operator_log_path
                << "\'\n";
        throw EcoevolityError(message.str());
    }

    std::cout << "Tree log path: " << tree_log_path << std::endl;
    std::cout << "State log path: " << state_log_path << std::endl;
    std::cout << "Operator log path: " << operator_log_path << std::endl;

    time_t start;
    time_t finish;
    time(&start);

    std::cout << "Firing up MCMC..." << std::endl;
    mcmc<TreeType>(
            rng,
            tree,
            operator_schedule,
            settings.get_chain_length(),
            settings.get_sample_frequency(),
            n_moves_per_generation,
            tree_log_stream,
            state_log_stream,
            operator_log_stream,
            std::cout,
            "\t",
            logging_precision,
            nthreads);

    tree_log_stream.close();
    state_log_stream.close();
    operator_log_stream.close();

    time(&finish);
    double duration = difftime(finish, start);
    std::cout << "Runtime: " << duration << " seconds." << std::endl;

    return 0;
}

#endif
