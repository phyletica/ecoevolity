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

#ifndef SIMCOEVOLITY_HPP
#define SIMCOEVOLITY_HPP

#include <limits>
#include <time.h>

#include "cpp-optparse/OptionParser.h"

#include "version.hpp"
#include "error.hpp"
#include "rng.hpp"
#include "path.hpp"
#include "settings.hpp"
#include "collection.hpp"


void write_sim_splash(std::ostream& out);

void check_output_path(const std::string& path);


template <class SettingsType, class CollectionType>
int simcoevolity_main(int argc, char * argv[]) {

    write_sim_splash(std::cout);
    std::cout << "\n";

    const std::string usage = 
        "usage: %prog [--seed SEED] [--ignore-data] YAML-CONFIG-FILE";

    std::ostringstream version_ss;
    version_ss << PROJECT_NAME << " version " << PROJECT_DETAILED_VERSION;
    const std::string version = version_ss.str();

    const std::string description =
        "Simcoevolity: Simulating evolutionary coevality";
    // const std::string epilog =
    //     "Epilog goes here...";

    optparse::OptionParser parser = optparse::OptionParser()
            .usage(usage)
            .version(version)
            .description(description);
            // .epilog(epilog);

    parser.add_option("--seed")
            .action("store")
            .type("long")
            .dest("seed")
            .help("Random number seed. Default: Set from clock.");
    parser.add_option("-n", "--number-of-replicates")
            .action("store")
            .type("unsigned int")
            .dest("number_of_replicates")
            .set_default("100")
            .help("Number of simulation replicates. Default: 100.");
    parser.add_option("-o", "--output-directory")
            .action("store")
            .dest("output_directory")
            .set_default("")
            .help("The directory into which simulated alignments and "
                  "associated config files will be written. Default: "
                  "Use directory of YAML config file.");
    parser.add_option("-p", "--prior")
            .action("store")
            .dest("prior")
            .set_default("")
            .help("The path to the configuration file that contains the "
                  "priors you would like to use when you analyse the "
                  "simulated datasets. By default, the same priors will "
                  "specified in your subsequent analyses as were used to "
                  "simulate the datasets.");
    parser.add_option("--prefix")
            .action("store")
            .dest("prefix")
            .set_default("")
            .help("Optional string to prefix all output files.");
    parser.add_option("--relax-constant-sites")
            .action("store_true")
            .dest("relax_constant_sites")
            .help("By default, if you specify \'constant_sites_removed = true\' "
                  "and constant sites are found, Simcoevolity throws an error. "
                  "With this option, Simcoevolity will automatically ignore the "
                  "constant sites and only issue a warning (and correct for "
                  "constant sites in the likelihood calculation). Please make sure "
                  "you understand what you are doing when you use this option."
                );
    parser.add_option("--relax-missing-sites")
            .action("store_true")
            .dest("relax_missing_sites")
            .help("By default, if a column is found for which there is no data "
                  "for at least one population, Simcoevolity throws an error. "
                  "With this option, Simcoevolity will automatically ignore such "
                  "sites and only issue a warning."
                );

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

    unsigned int nreps = options.get("number_of_replicates");
    if (nreps < 1) {
        throw EcoevolityError(
                "Number of simulation replicates must be 1 or greater");
    }
    std::cout << "Number of simulation replicates: " << nreps << std::endl;

    const bool strict_on_constant_sites = (! options.get("relax_constant_sites"));
    const bool strict_on_missing_sites = (! options.get("relax_missing_sites"));

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

    std::string prior_config_path = config_path;
    bool using_prior_config = false;
    if (options.is_set_by_user("prior")) {
        using_prior_config = true;
        prior_config_path = options.get("prior").get_str();
        if (! path::exists(prior_config_path)) {
            throw EcoevolityError("Config file \'" + prior_config_path +
                    "\' does not exist");
        }
        if (! path::isfile(prior_config_path)) {
            throw EcoevolityError("Config path \'" + prior_config_path +
                    "\' is not a regular file");
        }
    }

    std::string output_dir = path::dirname(config_path);
    if (options.is_set_by_user("output_directory")) {
        output_dir = options.get("output_directory").get_str();
        if (! path::exists(output_dir)) {
            throw EcoevolityError("Output directory \'" + output_dir +
                    "\' does not exist");
        }
        if (! path::isdir(output_dir)) {
            throw EcoevolityError("Output directory \'" + output_dir +
                    "\' is not a regular directory");
        }
    }

    std::string output_prefix = "";
    if (options.is_set_by_user("prefix")) {
        output_prefix = options.get("prefix").get_str() + "-";
    }
    output_prefix += "simcoevolity-";

    std::cout << "Prior config path: " << prior_config_path << std::endl;

    std::cout << "Parsing config file..." << std::endl;
    SettingsType settings = SettingsType(config_path);

    SettingsType prior_settings = SettingsType(prior_config_path);

    if (using_prior_config) {
        if (! settings.same_comparison_paths(prior_settings)) {
            throw EcoevolityError(
                    "The comparison files specified in \'" + config_path +
                    "\' and \'" + prior_config_path +
                    "\' do not match");
        }
    }

    std::string sim_settings_path = path::join(
            output_dir,
            output_prefix + "model-used-for-sims.yml");
    check_output_path(sim_settings_path);
    std::ofstream sim_settings_stream;
    sim_settings_stream.open(sim_settings_path);
    settings.write_settings(sim_settings_stream);
    sim_settings_stream.close();

    std::cout << "Configuring model for simulations..." << std::endl;
    CollectionType comparisons = CollectionType(
            settings,
            rng,
            strict_on_constant_sites,
            strict_on_missing_sites);

    if (using_prior_config) {
        // Not used but creating instance to vet settings
        std::cout << "Vetting model for analyses of simulated data sets..." << std::endl;
        CollectionType prior_comparisons = CollectionType(
                prior_settings,
                rng,
                strict_on_constant_sites,
                strict_on_missing_sites);
    }

    std::cout << "\n" << string_util::banner('-') << "\n";
    comparisons.write_summary(std::cout);
    std::cout << string_util::banner('-') << "\n\n";

    time_t start;
    time_t finish;
    time(&start);

    unsigned int pad_width = std::to_string(nreps).size();
    std::string sim_prefix = path::join(output_dir,
            output_prefix + "sim-");
    std::map<std::string, BiallelicData> sim_alignments;
    for (unsigned int i = 0; i < nreps; ++i) {
        std::string rep_str = string_util::pad_int(i, pad_width);
        std::string analysis_config_path = sim_prefix + rep_str + "-config.yml";
        check_output_path(analysis_config_path);
        std::string true_state_path = sim_prefix + rep_str + "-true-values.txt";
        check_output_path(true_state_path);

        comparisons.draw_from_prior(rng);
        sim_alignments = comparisons.simulate_biallelic_data_sets(rng, true);

        std::ofstream true_state_stream;
        true_state_stream.open(true_state_path);
        true_state_stream.precision(comparisons.get_logging_precision());
        comparisons.write_state_log_header(true_state_stream);
        comparisons.log_state(true_state_stream, 0);
        true_state_stream.close();

        std::ofstream sim_alignment_stream;
        for (auto const & k_v: sim_alignments) {
            std::string sim_alignment_path = sim_prefix + rep_str + "-" + path::basename(k_v.first);
            check_output_path(sim_alignment_path);

            char delim = prior_settings.get_population_name_delimiter(k_v.first);
            prior_settings.replace_comparison_path(k_v.first, path::basename(sim_alignment_path));

            sim_alignment_stream.open(sim_alignment_path);
            k_v.second.write_nexus(sim_alignment_stream, delim);
            sim_alignment_stream.close();
        }
        prior_settings.blanket_set_population_name_is_prefix(true);
        std::ofstream analysis_settings_stream;
        analysis_settings_stream.open(analysis_config_path);
        prior_settings.write_settings(analysis_settings_stream);
        analysis_settings_stream.close();
        for (auto const & k_v: sim_alignments) {
            std::string sim_alignment_path = sim_prefix + rep_str + "-" + path::basename(k_v.first);
            prior_settings.replace_comparison_path(path::basename(sim_alignment_path), k_v.first);
        }
    }

    time(&finish);
    double duration = difftime(finish, start);
    std::cout << "Runtime: " << duration << " seconds." << std::endl;

    return 0;
}

#endif
