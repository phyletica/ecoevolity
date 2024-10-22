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

#ifndef NEX2YML_HPP
#define NEX2YML_HPP

#include <limits>
#include <time.h>

#include "cpp-optparse/OptionParser.h"

#include "version.hpp"
#include "error.hpp"
#include "settings.hpp"
#include "collection.hpp"
#include "general_tree_settings.hpp"
#include "settings_io.hpp"


void write_nex2yml_splash(std::ostream& out);


template <class SettingsType, class CollectionType, class TreeType>
int nex2yml_main(int argc, char * argv[]) {

    write_nex2yml_splash(std::cout);
    std::cout << "\n";

    const std::string usage = 
        "usage: %prog [OPTIONS] YAML-CONFIG-FILE";

    std::ostringstream version_ss;
    version_ss << PROJECT_NAME << " version " << PROJECT_DETAILED_VERSION;
    const std::string version = version_ss.str();

    const std::string description =
        "nex2yml: convert nexus alignments to YAML files of allele count patterns";
    // const std::string epilog =
    //     "Epilog goes here...";

    optparse::OptionParser parser = optparse::OptionParser()
            .usage(usage)
            .version(version)
            .description(description);
            // .epilog(epilog);

    parser.add_option("--relax-constant-sites")
            .action("store_true")
            .dest("relax_constant_sites")
            .help("By default, if you specify \'constant_sites_removed = true\' "
                  "and constant sites are found, nex2yml throws an error. "
                  "With this option, nex2yml will automatically ignore the "
                  "constant sites and only issue a warning."
                );
    parser.add_option("--relax-missing-sites")
            .action("store_true")
            .dest("relax_missing_sites")
            .help("By default, if a column is found for which there is no data "
                  "across all populations, nex2yml throws an error. "
                  "With this option, nex2yml will automatically ignore such "
                  "sites and only issue a warning."
                );
    parser.add_option("--relax-triallelic-sites")
            .action("store_true")
            .dest("relax_triallelic_sites")
            .help("By default, if a DNA site is found for which there is more "
                  "than two nucleotide states, nex2yml throws an error. "
                  "With this option, nex2yml will automatically recode such "
                  "sites as biallelic and only issue a warning. These sites "
                  "are recoded by assigning state 0 to the first nucleotide "
                  "found and state 1 to all others. If you do not wish to "
                  "recode such sites and prefer to ignore them, please remove "
                  "all sites with more than two nucleotide states from your "
                  "DNA alignments. NOTE: only alignments of nucleotides are "
                  "affected by this option, not alignments of standard "
                  "characters (i.e., 0, 1, 2)."
                );

    optparse::Values& options = parser.parse_args(argc, argv);
    std::vector<std::string> args = parser.args();

    RandomNumberGenerator rng;

    const bool strict_on_constant_sites = (! options.get("relax_constant_sites"));
    const bool strict_on_missing_sites = (! options.get("relax_missing_sites"));
    const bool strict_on_triallelic_sites = (! options.get("relax_triallelic_sites"));

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
    SettingsType settings;
    PopulationTreeSettings phyco_settings;
    bool using_phyco_settings = false;
    try {
        settings = SettingsType(config_path);
    }
    catch (...) {
        using_phyco_settings = true;
        try {
            phyco_settings = PopulationTreeSettings(config_path);
        }
        catch (...) {
            std::cerr << "ERROR: Problem parsing config file: "
                    << config_path << "\n";
            throw;
        }
    }

    if (! using_phyco_settings) {
        std::cout << "Configuring comparisons..." << std::endl;
        CollectionType comparisons = CollectionType(
                settings,
                rng,
                strict_on_constant_sites,
                strict_on_missing_sites,
                strict_on_triallelic_sites);

        std::cout << "\n" << string_util::banner('-') << "\n";
        comparisons.write_summary(std::cout);
        std::cout << string_util::banner('-') << "\n\n";

        std::cout << "Writing YAML files of allele count patterns..." << std::endl;
        std::ofstream yml_out_stream;
        for (unsigned int i = 0;
                i < comparisons.get_number_of_trees();
                ++i) {
            const BiallelicData & d = comparisons.get_tree(i)->get_data();
            std::string yml_out_path = d.get_path() + ".yml";

            std::cout << "    Writing \'" << yml_out_path << "\'" << std::endl;
            try {
                yml_out_stream.open(yml_out_path);
            }
            catch (...) {
                std::cerr << "An error occurred when trying to open \'"
                          << yml_out_path
                          << "\'\n";
                throw;
            }
            try {
                d.write_yaml(yml_out_stream);
            }
            catch (...) {
                yml_out_stream.close();
                std::cerr << "An error occurred when trying to write to \'"
                          << yml_out_path
                          << "\'\n";
                throw;
            }
            yml_out_stream.close();
        }
    }
    else {
        std::cout << "Configuring model..." << std::endl;

        TreeType tree(
                phyco_settings,
                rng,
                strict_on_constant_sites,
                strict_on_missing_sites,
                strict_on_triallelic_sites,
                false // store_seq_loci_info
                );

        GeneralTreeOperatorSchedule<BasePopulationTree> operator_schedule(
                phyco_settings.operator_settings, tree.get_leaf_node_count());

        std::cout << "\n" << string_util::banner('-') << "\n";
        write_settings(std::cout, phyco_settings, operator_schedule);
        std::cout << string_util::banner('-') << "\n\n";

        std::cout << "Writing YAML files of allele count patterns..." << std::endl;
        std::ofstream yml_out_stream;
        const BiallelicData & d = tree.get_data();
        std::string yml_out_path = d.get_path() + ".yml";

        std::cout << "    Writing \'" << yml_out_path << "\'" << std::endl;
        try {
            yml_out_stream.open(yml_out_path);
        }
        catch (...) {
            std::cerr << "An error occurred when trying to open \'"
                      << yml_out_path
                      << "\'\n";
            throw;
        }
        try {
            d.write_yaml(yml_out_stream);
        }
        catch (...) {
            yml_out_stream.close();
            std::cerr << "An error occurred when trying to write to \'"
                      << yml_out_path
                      << "\'\n";
            throw;
        }
        yml_out_stream.close();
    }
    return 0;
}

#endif
