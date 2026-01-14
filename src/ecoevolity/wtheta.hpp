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

#ifndef WTHETA_HPP
#define WTHETA_HPP

#include <limits>
#include <time.h>

#include "cpp-optparse/OptionParser.h"

#include "version.hpp"
#include "error.hpp"
#include "settings.hpp"
#include "collection.hpp"
#include "general_tree_settings.hpp"
#include "settings_io.hpp"


void write_wtheta_splash(std::ostream& out);


template <class SettingsType, class CollectionType, class TreeType>
int wtheta_main(int argc, char * argv[]) {

    write_wtheta_splash(std::cerr);
    std::cerr << "\n";

    const std::string usage = 
        "usage: %prog [OPTIONS] YAML-CONFIG-FILE";

    std::ostringstream version_ss;
    version_ss << PROJECT_NAME << " version " << PROJECT_DETAILED_VERSION;
    const std::string version = version_ss.str();

    const std::string description =
        "wtheta: Estimate Watterson's theta for each population";
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
                  "and constant sites are found, wtheta throws an error. "
                  "With this option, wtheta will automatically ignore the "
                  "constant sites and only issue a warning."
                );
    parser.add_option("--relax-missing-sites")
            .action("store_true")
            .dest("relax_missing_sites")
            .help("By default, if a column is found for which there is no data "
                  "across all populations, wtheta throws an error. "
                  "With this option, wtheta will automatically ignore such "
                  "sites and only issue a warning."
                );
    parser.add_option("--relax-triallelic-sites")
            .action("store_true")
            .dest("relax_triallelic_sites")
            .help("By default, if a DNA site is found for which there is more "
                  "than two nucleotide states, wtheta throws an error. "
                  "With this option, wtheta will automatically recode such "
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
    std::cerr << "Config path: " << config_path << std::endl;

    std::cerr << "Parsing config file..." << std::endl;
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

    std::cout.precision(12);

    if (! using_phyco_settings) {
        std::cerr << "Configuring comparisons..." << std::endl;
        CollectionType comparisons = CollectionType(
                settings,
                rng,
                strict_on_constant_sites,
                strict_on_missing_sites,
                strict_on_triallelic_sites);

        std::cerr << "\n" << string_util::banner('-') << "\n";
        comparisons.write_summary(std::cerr);
        std::cerr << string_util::banner('-') << "\n\n";

        std::cerr << "Estimating Watterson's theta per population..." << std::endl;
        std::cout << "population\twattersons_theta\tmax_number_of_alleles" << std::endl; 
        for (unsigned int i = 0;
                i < comparisons.get_number_of_trees();
                ++i) {
            const BiallelicData & d = comparisons.get_tree(i)->get_data();

            for (unsigned int pop_idx = 0;
                    pop_idx < d.get_number_of_populations();
                    ++pop_idx) {
                const std::string & pop_label = d.get_population_label(pop_idx);
                const unsigned int max_allele_count = d.get_max_allele_count(pop_idx);
                const double w = d.get_wattersons_theta(pop_idx);
                std::cout << pop_label << "\t" << w << "\t" << max_allele_count << std::endl;
            }
        }
    }
    else {
        std::cerr << "Configuring model..." << std::endl;

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

        std::cerr << "\n" << string_util::banner('-') << "\n";
        write_settings(std::cerr, phyco_settings, operator_schedule);
        std::cerr << string_util::banner('-') << "\n\n";

        std::cerr << "Estimating Watterson's theta per population..." << std::endl;
        const BiallelicData & d = tree.get_data();

        std::cout << "population\twattersons_theta\tmax_number_of_alleles" << std::endl; 
        for (unsigned int pop_idx = 0;
                pop_idx < d.get_number_of_populations();
                ++pop_idx) {
            const std::string & pop_label = d.get_population_label(pop_idx);
            const unsigned int max_allele_count = d.get_max_allele_count(pop_idx);
            const double w = d.get_wattersons_theta(pop_idx);
            std::cout << pop_label << "\t" << w << "\t" << max_allele_count << std::endl;
        }
    }
    return 0;
}

#endif
