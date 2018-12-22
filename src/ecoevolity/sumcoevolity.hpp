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

#ifndef SUMCOEVOLITY_HPP
#define SUMCOEVOLITY_HPP

#include <limits>
#include <time.h>

#include "cpp-optparse/OptionParser.h"

#include "version.hpp"
#include "error.hpp"
#include "rng.hpp"
#include "util.hpp"
#include "math_util.hpp"
#include "probability.hpp"
#include "string_util.hpp"
#include "settings.hpp"
#include "spreadsheet.hpp"


void write_sumcoevolity_splash(std::ostream& out);

void check_sumcoevolity_output_path(const std::string& path);

template <class SettingsType>
int sumcoevolity_main(int argc, char * argv[]) {

    write_sumcoevolity_splash(std::cerr);
    std::cerr << "\n";

    const std::string usage = 
        "usage: %prog [OPTIONS] ECOEVOLITY-STATE-LOG-FILE-1 [ECOEVOLITY-STATE-LOG-FILE-2 ...]";

    std::ostringstream version_ss;
    version_ss << PROJECT_NAME << " version " << PROJECT_DETAILED_VERSION;
    const std::string version = version_ss.str();

    const std::string description =
        "Sumcoevolity: Summarizing evolutionary coevality";
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
            .help("Number of samples from the beginning of each log file to "
                  "ignore as burn in. "
                  "Default: 0.");
    parser.add_option("-p", "--prefix")
            .action("store")
            .dest("prefix")
            .set_default("")
            .help("Optional string to prefix all output files.");
    parser.add_option("-c", "--config")
            .action("store")
            .dest("config")
            .set_default("")
            .help("Path to the YAML config file used to generate the provided "
                  "log files. If provided, prior probabilities and Bayes "
                  "factors will be approximated via simulation.");
    parser.add_option("-n", "--number-of-samples")
            .action("store")
            .type("unsigned int")
            .dest("number_of_samples")
            .set_default("100000")
            .help("Number of simulation samples. "
                  "Only used if YAML config file is provided. "
                  "Default: 100000.");
    parser.add_option("--seed")
            .action("store")
            .type("long")
            .dest("seed")
            .help("Seed for random number generator. "
                  "Only used if YAML config file is provided. "
                  "Default: Set from clock.");
    parser.add_option("-f", "--force")
            .action("store_true")
            .dest("force")
            .help("Force overwriting of existing output files, if they exist.");

    optparse::Values& options = parser.parse_args(argc, argv);
    std::vector<std::string> log_paths = parser.args();

    if (log_paths.size() < 1) {
        throw EcoevolityError("Too few positional arguments; "
                "a path to at least one log file is required.");
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

    unsigned int nreps = options.get("number_of_samples");
    if (nreps < 1) {
        throw EcoevolityError(
                "Number of samples must be 1 or greater");
    }

    unsigned int burnin = options.get("burnin");
    // Using unsigned int for burnin, so no need to check if negative
    // if (burnin < 0) {
    //     throw EcoevolityError(
    //             "Burn in must be 0 or greater");
    // }

    std::string config_path = options.get("config").get_str();
    bool running_sims = false;
    if (options.is_set_by_user("config")) {
        running_sims = true;
        if (! path::exists(config_path)) {
            throw EcoevolityError("Config path \'" + config_path +
                    "\' does not exist");
        }
        if (! path::isfile(config_path)) {
            throw EcoevolityError("Config path \'" + config_path +
                    "\' is not a regular file");
        }
    }

    std::string prefix = "";
    if (options.is_set_by_user("prefix")) {
        prefix = options.get("prefix").get_str();
    }

    bool prevent_overwrite = (! options.get("force"));


    // Start real work
    time_t start;
    time_t finish;
    time(&start);

    std::cerr << "Parsing log files...\n";
    spreadsheet::Spreadsheet posterior_sample;
    posterior_sample.update(log_paths, burnin);
    std::vector<int> nevents = posterior_sample.get<int>("number_of_events");

    unsigned int number_of_posterior_samples = nevents.size();
    if (number_of_posterior_samples < 1) {
        throw EcoevolityError(
                "No samples were parsed from the log files. "
                "Perhaps you specified a burn in value larger than the number "
                "of samples in each log file?"
                );
    }

    std::cerr << "Parsed " << number_of_posterior_samples
              << " total samples from "
              << log_paths.size() << " log files.\n";

    std::vector<std::string> keys = posterior_sample.get_keys();
    std::vector<std::vector<unsigned int> > indices;
    for (auto const & k: keys) { 
        if (string_util::startswith(k, "root_height_index_")) {
            indices.push_back(posterior_sample.get<unsigned int>(k));
        }
    }
    unsigned int number_of_comparisons = indices.size();

    std::map<unsigned int, unsigned int> nevents_counts;
    std::map<unsigned int, unsigned int> prior_nevents_counts;
    for (unsigned int i = 1; i <= number_of_comparisons; ++i) {
        nevents_counts[i] = 0;
        prior_nevents_counts[i] = 0;
    }

    std::map<std::vector<unsigned int>, unsigned int> model_counts;
    std::map<std::vector<unsigned int>, unsigned int> prior_model_counts;

    for (unsigned int i = 0; i < nevents.size(); ++i) {
        std::vector<unsigned int> model;
        model.reserve(number_of_comparisons);
        for (unsigned int j = 0; j < number_of_comparisons; ++j) {
            model.push_back(indices.at(j).at(i));
        }
        if (model_counts.count(model) < 1) {
            model_counts[model] = 1;
            prior_model_counts[model] = 0;
        }
        else {
            ++model_counts[model];
        }
        ++nevents_counts[nevents.at(i)];
    }
    unsigned int tally = 0;
    for (auto const & kv: nevents_counts) {
        tally += kv.second;
    }
    ECOEVOLITY_ASSERT(tally == number_of_posterior_samples);
    tally = 0;
    for (auto const & kv: model_counts) {
        tally += kv.second;
    }
    ECOEVOLITY_ASSERT(tally == number_of_posterior_samples);

    // Sort descending by counts
    std::vector< std::pair<unsigned int, unsigned int> > nevents_count_pairs;
    for (auto const & kv: nevents_counts) {
        nevents_count_pairs.push_back(kv);
    }
    sort_pairs(nevents_count_pairs, false, true);
    std::vector< std::pair<std::vector<unsigned int>, unsigned int> > model_count_pairs;
    for (auto const & kv: model_counts) {
        model_count_pairs.push_back(kv);
    }
    sort_pairs(model_count_pairs, false, true);

    // Run simulations
    if (running_sims) {
        SettingsType settings = SettingsType(config_path);
        if (settings.get_number_of_comparisons() != number_of_comparisons) {
            throw EcoevolityError("Number of comparisons found in log files "
                    "and config do not match");
        }
        if (settings.event_model_is_fixed()) {
            throw EcoevolityError("The model is fixed in the config file, so "
                    "prior probabilities under the specified Dirichlet process "
                    "settings will not be comparable to the posterior "
                    "probabilities.");
        }
        if (settings.get_model_prior() == EcoevolityOptions::ModelPrior::uniform) {
            throw EcoevolityError("Approximating prior probabilites via "
                    "simulation is not supported for the uniform model "
                    "prior");
        }
        const ModelOperatorSettings & model_settings = settings.get_operator_schedule_settings().get_model_operator_settings();
        if (model_settings.get_weight() <= 0.0) {
            throw EcoevolityError("The model operator in the config file "
                    "was turned off (given no weight), so prior probabilities "
                    "under the specified Dirichlet process settings will not "
                    "be comparable to the posterior probabilities (because the "
                    "model was fixed in the analysis).");
        }

        std::cerr << "Approximating prior probabilities via simulations...\n";
        const PositiveRealParameterSettings & concentration_settings = settings.get_concentration_settings();
        std::shared_ptr<ContinuousProbabilityDistribution> concentration_prior = concentration_settings.get_prior_settings().get_instance();

        std::cerr << "Simulations settings:\n";
        std::cerr << "\tseed = " << seed << "\n";
        double concentration = 1.0;
        bool concentration_is_fixed = false;
        if (concentration_settings.is_fixed()) {
            concentration = concentration_settings.get_value();
            concentration_is_fixed = true;
            std::cerr << "\tconcentration = " << concentration << "\n";
        }
        else {
            std::cerr << "\tconcentration ~ " << concentration_prior->to_string() << "\n";
        }

        for (unsigned int i = 0; i < nreps; ++i) {
            if (! concentration_is_fixed) {
                concentration = concentration_prior->draw(rng);
            }
            std::vector<unsigned int> model(number_of_comparisons, 0);
            unsigned int number_of_categories = rng.dirichlet_process(model, concentration);
            ++prior_nevents_counts[number_of_categories];
            if (prior_model_counts.count(model) < 1) {
                prior_model_counts[model] = 1;
            }
            else {
                ++prior_model_counts[model];
            }
        }
        tally = 0;
        for (auto const & kv: prior_nevents_counts) {
            tally += kv.second;
        }
        ECOEVOLITY_ASSERT(tally == nreps);
        tally = 0;
        for (auto const & kv: prior_model_counts) {
            tally += kv.second;
        }
        ECOEVOLITY_ASSERT(tally == nreps);
    }

    double min_post_prob = 1.0 / (double)number_of_posterior_samples;

    double min_prior_prob = 1.0 / (double)nreps;

    double max_post_prob = (number_of_posterior_samples - 1) / (double)number_of_posterior_samples;

    double max_prior_prob = (nreps - 1) / (double)nreps;

    // output nevents results
    std::string nevents_path = prefix + "sumcoevolity-results-nevents.txt";
    if (prevent_overwrite) {
        check_sumcoevolity_output_path(nevents_path);
    }

    std::cerr << "Writing summary for number of events...\n";
    std::ofstream nevents_stream;
    nevents_stream.open(nevents_path);
    if (! nevents_stream.is_open()) {
        std::ostringstream message;
        message << "ERROR: Could not output file \'"
                << nevents_path
                << "\'\n";
        throw EcoevolityError(message.str());
    }
    try {
        nevents_stream << "number_of_events\tpost_prob\tcumulative_post_prob\tprior_prob\tbf\n";
        double cumulative_post_prob = 0.0;
        for (auto const & nc: nevents_count_pairs) {
            bool min_post = false;
            bool max_post = false;
            bool min_prior = false;
            bool max_prior = false;

            double post_prob = nc.second / (double)number_of_posterior_samples;
            cumulative_post_prob += post_prob;

            nevents_stream << nc.first << "\t";
            if (nc.second < 1) {
                post_prob = min_post_prob;
                nevents_stream << "<" << post_prob << "\t";
                min_post = true;
            }
            else if (nc.second == number_of_posterior_samples) {
                post_prob = max_post_prob;
                nevents_stream << ">" << post_prob << "\t";
                max_post = true;
            }
            else {
                nevents_stream << post_prob << "\t";
            }
            nevents_stream << cumulative_post_prob << "\t";
            if (running_sims) {
                double prior_prob = prior_nevents_counts[nc.first] / (double)nreps;
                if (prior_nevents_counts[nc.first] < 1) {
                    prior_prob = min_prior_prob;
                    nevents_stream << "<" << prior_prob << "\t";
                    min_prior = true;
                }
                else if (prior_nevents_counts[nc.first] == nreps) {
                    prior_prob = max_prior_prob;
                    nevents_stream << ">" << prior_prob << "\t";
                    max_prior = true;
                }
                else {
                    nevents_stream << prior_prob << "\t";
                }
                double bf = (
                        (post_prob / (1.0 - post_prob)) /
                        (prior_prob / (1.0 - prior_prob))
                        );
                if ((max_post && max_prior) || (min_post && min_prior)) {
                    // Cannot know the BF under these conditions
                    nevents_stream << "NA\n";
                }
                else if (max_post || min_prior) {
                    nevents_stream << ">" << bf << "\n";
                }
                else if (min_post || max_prior) {
                    nevents_stream << "<" << bf << "\n";
                }
                else {
                    nevents_stream << bf << "\n";
                }
            }
            else {
                nevents_stream << "NA\tNA\n";
            }
        }
        nevents_stream.close();
    }
    catch (...) {
        nevents_stream.close();
        throw;
    }
    std::cerr << "Summary written to \'" << nevents_path << "\'\n";


    // output model results
    std::string model_path = prefix + "sumcoevolity-results-model.txt";
    if (prevent_overwrite) {
        check_sumcoevolity_output_path(model_path);
    }

    std::cerr << "Writing summary for event models...\n";
    std::ofstream model_stream;
    model_stream.open(model_path);
    if (! model_stream.is_open()) {
        std::ostringstream message;
        message << "ERROR: Could not output file \'"
                << model_path
                << "\'\n";
        throw EcoevolityError(message.str());
    }
    try {
        model_stream << "model\tpost_prob\tcumulative_post_prob\tprior_prob\tbf\n";
        double cumulative_post_prob = 0.0;
        for (auto const & mc: model_count_pairs) {
            bool min_post = false;
            bool max_post = false;
            bool min_prior = false;
            bool max_prior = false;

            double post_prob = mc.second / (double)number_of_posterior_samples;
            cumulative_post_prob += post_prob;

            model_stream << mc.first.at(0);
            for (unsigned int idx = 1; idx < mc.first.size(); ++idx) {
                model_stream << "," << mc.first.at(idx);
            }
            model_stream << "\t";
            if (mc.second < 1) {
                post_prob = min_post_prob;
                model_stream << "<" << post_prob << "\t";
                min_post = true;
            }
            else if (mc.second == number_of_posterior_samples) {
                post_prob = max_post_prob;
                model_stream << ">" << post_prob << "\t";
                max_post = true;
            }
            else {
                model_stream << post_prob << "\t";
            }
            model_stream << cumulative_post_prob << "\t";
            if (running_sims) {
                double prior_prob = prior_model_counts[mc.first] / (double)nreps;
                if (prior_model_counts[mc.first] < 1) {
                    prior_prob = min_prior_prob;
                    model_stream << "<" << prior_prob << "\t";
                    min_prior = true;
                }
                else if (prior_model_counts[mc.first] == nreps) {
                    prior_prob = max_prior_prob;
                    model_stream << ">" << prior_prob << "\t";
                    max_prior = true;
                }
                else {
                    model_stream << prior_prob << "\t";
                }
                double bf = (
                        (post_prob / (1.0 - post_prob)) /
                        (prior_prob / (1.0 - prior_prob))
                        );
                if ((max_post && max_prior) || (min_post && min_prior)) {
                    // Cannot know the BF under these conditions
                    model_stream << "NA\n";
                }
                else if (max_post || min_prior) {
                    model_stream << ">" << bf << "\n";
                }
                else if (min_post || max_prior) {
                    model_stream << "<" << bf << "\n";
                }
                else {
                    model_stream << bf << "\n";
                }
            }
            else {
                model_stream << "NA\tNA\n";
            }
        }
        model_stream.close();
    }
    catch (...) {
        model_stream.close();
        throw;
    }
    std::cerr << "Summary written to \'" << model_path << "\'\n";

    time(&finish);
    double duration = difftime(finish, start);
    std::cerr << "Runtime: " << duration << " seconds." << std::endl;

    return 0;
}

#endif
