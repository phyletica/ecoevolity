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

#include "simcoevolity.hpp"

void write_sim_splash(std::ostream& out) {
    std::string v = "Version ";
    v += PROJECT_DETAILED_VERSION;
    out << string_util::banner('=') << "\n" 
        << string_util::center("Simcoevolity") << "\n"
        << string_util::center("Simulating evolutionary coevality") << "\n"
        << string_util::center("Part of:") << "\n"
        << string_util::center(PROJECT_NAME) << "\n"
        << string_util::center(v) << "\n"
        << string_util::banner('=') << "\n";
}

int simcoevolity_main(int argc, char * argv[]) {

    write_splash(std::cout);
    std::cout << "\n";

    const std::string usage = 
        "usage: %prog [--seed SEED] [--ignore-data] YAML-CONFIG-FILE";

    std::ostringstream version_ss;
    version_ss << PROJECT_NAME << " version " << PROJECT_DETAILED_VERSION;
    const std::string version = version_ss.str();

    const std::string description =
        "Ecoevolity: Estimating evolutionary coevality";
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
            .help("Random number seed. Default: Set from clock.");
    parser.add_option("-p", "--prior")
            .set_default("")
            .dest("prior")
            .help("The path to the configuration file that contains the "
                  "priors you would like to use when you analyse the "
                  "simulated datasets. By default, the same priors will "
                  "specified in your subsequent analyses as were used to "
                  "simulate the datasets.");
    parser.add_option("--relax-constant-sites")
            .action("store_true")
            .dest("relax_constant_sites")
            .help("By default, if you specify \'constant_sites_removed = true\' "
                  "and constant sites are found, Ecoevolity throws an error. "
                  "With this option, Ecoevolity will automatically ignore the "
                  "constant sites and only issue a warning (and correct for "
                  "constant sites in the likelihood calculation). Please make sure "
                  "you understand what you are doing when you use this option."
                );
    parser.add_option("--relax-missing-sites")
            .action("store_true")
            .dest("relax_missing_sites")
            .help("By default, if a column is found for which there is no data "
                  "for at least one population, Ecoevolity throws an error. "
                  "With this option, Ecoevolity will automatically ignore such "
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
    if (options.is_set_by_user("prior")) {
        prior_config_path = options.get("prior");
        if (! path::exists(prior_config_path)) {
            throw EcoevolityError("Config file \'" + config_path +
                    "\' does not exist");
        }
        if (! path::isfile(prior_config_path)) {
            throw EcoevolityError("Config path \'" + config_path +
                    "\' is not a regular file");
        }
    }
    std::cout << "Prior config path: " << prior_config_path << std::endl;


    std::cout << "Parsing config file..." << std::endl;
    CollectionSettings settings = CollectionSettings(config_path);

    std::cout << "\n" << string_util::banner('-') << "\n";
    settings.write_settings(std::cout);
    std::cout << string_util::banner('-') << "\n\n";

    std::cout << "Configuring model..." << std::endl;
    ComparisonPopulationTreeCollection comparisons =
            ComparisonPopulationTreeCollection(
                    settings,
                    rng,
                    strict_on_constant_sites,
                    strict_on_missing_sites);

    std::cout << "\n" << string_util::banner('-') << "\n";
    comparisons.write_summary(std::cout);
    std::cout << string_util::banner('-') << "\n\n";

    time_t start;
    time_t finish;
    time(&start);

    time(&finish);
    double duration = difftime(finish, start);
    std::cout << "Runtime: " << duration << " seconds." << std::endl;

    return 0;
}
