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

    write_sim_splash(std::cerr);
    std::cerr << "\n";

    const std::string usage = 
        "usage: %prog [OPTIONS] YAML-CONFIG-FILE";

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
            .help("Seed for random number generator. Default: Set from clock.");
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
    parser.add_option("--singleton-sample-probability")
            .action("store")
            .type("double")
            .dest("singleton_sample_probability")
            .set_default("1.0")
            .help("The probability of sampling singleton site patterns. This "
                  "is used to simulate data acquisition biases against "
                  "character patterns where only one gene copy has an allele "
                  "that is different from all the others. "
                  "Default: 1.0 (no acquisition bias).");
    parser.add_option("-l", "--locus-size")
            .action("store")
            .type("unsigned int")
            .dest("locus_size")
            .set_default("1")
            .help("Number of sites simulated on each gene tree. Default: 1 "
                  "(every site is unlinked). NOTE: When locus size is greater "
                  "than 1, 'complete' datasets are simulated (i.e., no missing "
                  "data). When the locus size is 1 (default), the pattern of "
                  "missing data in every simulated aligment will match the "
                  "corresponding empirical alignment");
    parser.add_option("--max-one-variable-site-per-locus")
            .action("store_true")
            .dest("max_one_variable_site_per_locus")
            .help("When locus size is greater than one, this option forces "
                  "only one variable site per locus to be retained. "
                  "This is for simulating datasets that mimic retaining only "
                  "one SNP per RADseq or sequence-capture locus, for example. "
                  "NOTE: this only works when the locus size is greater than "
                  "one, and each comparison specified in the config file has "
                  "not had constant sites removed (i.e., "
                  "\'constant_sites_removed: false\').");
    parser.add_option("-c", "--charsets")
            .action("store_true")
            .dest("charsets")
            .help("Use charsets defined in the nexus-formatted alignment files. "
                  "When this option is used, and there is a \'sets\' block in "
                  "every comparison\'s nexus-formatted alignment file that has "
                  "\'charsets\' delineating all the locus boundaries, then "
                  "multi-locus datasets are simulated that match the empirical "
                  "dataset (i.e., the same number and length of loci, and the "
                  "same patterns of missing data for each locus). NOTE: This "
                  "option will override the \'--locus-size\' option, but can "
                  "be used in combination with the "
                  "\'--max-one-variable-site-per-locus\' option.");
    parser.add_option("--parameters-only")
            .action("store_true")
            .dest("parameters_only")
            .help("By default, sequence alignments and associated config files "
                  "are generated for each replicate. When this option is "
                  "specified, only the parameter values drawn for each "
                  "replicate are written to standard output (no sequence data "
                  "are simulated). Because no data or configs are created, the "
                  "settings for the '-o/--output-directory', '-l/--locus-size' "
                  "and '-p/--prior' options will be ignored."
                );
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
                  "constant sites and only issue a warning. Simulated data sets "
                  "will have fewer sites than the original alignment."
                );
    parser.add_option("--relax-missing-sites")
            .action("store_true")
            .dest("relax_missing_sites")
            .help("By default, if a column is found for which there is no data "
                  "across all populations, Simcoevolity throws an error. "
                  "With this option, Simcoevolity will automatically ignore such "
                  "sites and only issue a warning. Simulated data sets will have "
                  "fewer sites than the original alignment."
                );
    parser.add_option("--relax-triallelic-sites")
            .action("store_true")
            .dest("relax_triallelic_sites")
            .help("By default, if a DNA site is found for which there is more "
                  "than two nucleotide states, Simcoevolity throws an error. "
                  "With this option, Simcoevolity will only issue a warning. "
                  "The number of sites with more than two nucleotides will be "
                  "retained, but all simulated sites will have at most two "
                  "character states."
                );
    parser.add_option("--nexus")
            .action("store_true")
            .dest("output_nexus")
            .help("Output simulated data in nexus format, rather than the "
                  "default YAML format."
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
    std::cerr << "Seed: " << seed << std::endl;

    const unsigned int nreps = options.get("number_of_replicates");
    if (nreps < 1) {
        throw EcoevolityError(
                "Number of simulation replicates must be 1 or greater");
    }
    std::cerr << "Number of simulation replicates: " << nreps << std::endl;

    const double singleton_sample_probability = options.get(
            "singleton_sample_probability");

    const unsigned int locus_size = options.get("locus_size");
    if (locus_size < 1) {
        throw EcoevolityError(
                "Number of sites simulated per locus must be 1 or greater");
    }
    std::cerr << "Number of sites simulated per locus: " << locus_size << std::endl;

    const bool use_charsets = options.get("charsets");

    const bool max_one_variable_site_per_locus = options.get(
            "max_one_variable_site_per_locus");
    if ((! use_charsets) && (max_one_variable_site_per_locus && (locus_size < 2))) {
        throw EcoevolityError(
                "Locus length must be 1 or greater to use "
                "\'--max-one-variable-site-per-locus\'");
    }

    const bool strict_on_constant_sites = (! options.get("relax_constant_sites"));
    const bool strict_on_missing_sites = (! options.get("relax_missing_sites"));
    const bool strict_on_triallelic_sites = (! options.get("relax_triallelic_sites"));
    const bool simulate_sequences = (! options.get("parameters_only"));
    const bool output_nexus = options.get("output_nexus");

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
    std::cerr << "Config path: " << config_path << std::endl;

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
        output_prefix = options.get("prefix").get_str();
    }
    output_prefix += "simcoevolity-";

    std::cerr << "Prior config path: " << prior_config_path << std::endl;

    std::cerr << "Parsing config file..." << std::endl;
    SettingsType settings = SettingsType(config_path);

    if (max_one_variable_site_per_locus && (! settings.all_comparions_have_constant_sites())) {
        throw EcoevolityError(
                "All comparisons must not have constant sites removed when "
                "using the \'--max-one-variable-site-per-locus\' option.");
    }

    SettingsType prior_settings = SettingsType(prior_config_path);

    if (using_prior_config) {
        if (! settings.same_comparison_paths(prior_settings)) {
            throw EcoevolityError(
                    "The comparison files specified in \'" + config_path +
                    "\' and \'" + prior_config_path +
                    "\' do not match");
        }
    }

    if (simulate_sequences) {
        std::string sim_settings_path = path::join(
                output_dir,
                output_prefix + "model-used-for-sims.yml");
        check_output_path(sim_settings_path);
        std::ofstream sim_settings_stream;
        sim_settings_stream.open(sim_settings_path);
        settings.write_settings(sim_settings_stream);
        sim_settings_stream.close();
    }
    else {
        settings.write_settings(std::cerr);
    }

    std::cerr << "Configuring model for simulations..." << std::endl;
    CollectionType comparisons = CollectionType(
            settings,
            rng,
            strict_on_constant_sites,
            strict_on_missing_sites,
            strict_on_triallelic_sites,
            use_charsets);

    if (using_prior_config) {
        // Not used but creating instance to vet settings
        std::cerr << "Vetting model for analyses of simulated data sets..." << std::endl;
        CollectionType prior_comparisons = CollectionType(
                prior_settings,
                rng,
                strict_on_constant_sites,
                strict_on_missing_sites,
                strict_on_triallelic_sites,
                use_charsets);
    }

    if (use_charsets) {
        if (! comparisons.has_seq_loci_info()) {
            throw EcoevolityError(
                    "All comparisons must have charsets defined when "
                    "using the \'--charsets\' option.");
        }
        std::cerr << "Using charset info to simulate multi-locus data..." << std::endl;
    }

    std::cerr << "\n" << string_util::banner('-') << "\n";
    comparisons.write_summary(std::cerr);
    std::cerr << string_util::banner('-') << "\n\n";

    time_t start;
    time_t finish;
    time(&start);

    std::cerr << "Starting simulations..." << std::endl;
    if (simulate_sequences) {
        unsigned int pad_width = std::to_string(nreps).size();
        std::string sim_prefix = path::join(output_dir,
                output_prefix + "sim-");
        std::map<std::string, BiallelicData> sim_alignments;
        for (unsigned int i = 0; i < nreps; ++i) {
            std::cerr << "Simulating data set " << (i + 1) << " of " << nreps << "\n";
            std::string rep_str = string_util::pad_int(i, pad_width);
            std::string analysis_config_path = sim_prefix + rep_str + "-config.yml";
            check_output_path(analysis_config_path);
            std::string true_state_path = sim_prefix + rep_str + "-true-values.txt";
            check_output_path(true_state_path);

            comparisons.draw_from_prior(rng);
            if (use_charsets) {
                sim_alignments = comparisons.simulate_linked_biallelic_data_sets(rng,
                        singleton_sample_probability,
                        max_one_variable_site_per_locus,
                        true);
            }
            else if (locus_size < 2) {
                sim_alignments = comparisons.simulate_biallelic_data_sets(rng,
                        singleton_sample_probability,
                        true);
            }
            else {
                sim_alignments = comparisons.simulate_complete_biallelic_data_sets(rng,
                        locus_size,
                        singleton_sample_probability,
                        max_one_variable_site_per_locus,
                        true);
            }

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
                if (output_nexus) {
                    k_v.second.write_nexus(sim_alignment_stream, delim);
                }
                else {
                    k_v.second.write_yaml(sim_alignment_stream);
                }
                sim_alignment_stream.close();
            }
            prior_settings.blanket_set_population_name_is_prefix(true);
            prior_settings.blanket_set_genotypes_are_diploid(false);
            if (max_one_variable_site_per_locus) {
                prior_settings.blanket_set_constant_sites_removed(true);
            }
            std::ofstream analysis_settings_stream;
            analysis_settings_stream.open(analysis_config_path);
            prior_settings.write_settings(analysis_settings_stream);
            analysis_settings_stream.close();
            for (auto const & k_v: sim_alignments) {
                std::string sim_alignment_path = sim_prefix + rep_str + "-" + path::basename(k_v.first);
                prior_settings.replace_comparison_path(path::basename(sim_alignment_path), k_v.first);
            }
        }
    }
    else {
        std::ostream & state_stream = std::cout;
        state_stream.precision(comparisons.get_logging_precision());
        comparisons.write_state_log_header(state_stream);
        std::cerr << "Only drawing samples of parameters and writing to stdout." << std::endl;
        for (unsigned int i = 0; i < nreps; ++i) {
            comparisons.draw_from_prior(rng);
            comparisons.log_state(state_stream, i + 1);
        }
    }

    time(&finish);
    double duration = difftime(finish, start);
    std::cerr << "Runtime: " << duration << " seconds." << std::endl;

    return 0;
}

#endif
