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
#include "rng.hpp"
#include "path.hpp"
#include "string_util.hpp"
#include "general_tree_settings.hpp"
#include "settings_io.hpp"
#include "mcmc.hpp"


void write_sum_phy_splash(std::ostream& out);

void check_sumphy_output_path(const std::string& path);

template <class TreeType>
int sumphycoeval_main(int argc, char * argv[]) {

    write_sum_phy_splash(std::cout);
    std::cout << "\n";

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
    parser.add_option("-p", "--prefix")
            .action("store")
            .dest("prefix")
            .set_default("")
            .help("Optional string to prefix all output files.");
    parser.add_option("-f", "--force")
            .action("store_true")
            .dest("force")
            .help("Force overwriting of existing output files, if they exist.");

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
    unsigned int burnin = options.get("burnin");

    std::string prefix = "";
    if (options.is_set_by_user("prefix")) {
        prefix = options.get("prefix").get_str();
    }

    bool prevent_overwrite = (! options.get("force"));


    time_t start;
    time_t finish;
    time(&start);

    time(&finish);
    double duration = difftime(finish, start);
    std::cerr << "Runtime: " << duration << " seconds." << std::endl;

    return 0;
}

#endif
