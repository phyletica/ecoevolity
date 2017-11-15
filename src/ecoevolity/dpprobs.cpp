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

#include "dpprobs.hpp"


void write_dpprobs_splash(std::ostream& out) {
    std::string v = "Version ";
    v += PROJECT_DETAILED_VERSION;
    out << string_util::banner('=') << "\n" 
        << string_util::center("DPprobs") << "\n"
        << string_util::center("Simulating Dirichlet process probabilities") << "\n\n"
        << string_util::center("Part of:") << "\n"
        << string_util::center(PROJECT_NAME) << "\n"
        << string_util::center(v) << "\n"
        << string_util::banner('=') << "\n";
}

int dpprobs_main(int argc, char * argv[]) {

    write_dpprobs_splash(std::cerr);
    std::cerr << "\n";

    const std::string usage = 
        "usage: %prog [OPTIONS] PARAMETER-VALUE NUMBER-OF-ELEMENTS";

    std::ostringstream version_ss;
    version_ss << PROJECT_NAME << " version " << PROJECT_DETAILED_VERSION;
    const std::string version = version_ss.str();

    const std::string description =
        "DPprobs: Simulating Dirichlet Process probabilites";
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
    parser.add_option("-n", "--number-of-samples")
            .action("store")
            .type("unsigned int")
            .dest("number_of_samples")
            .set_default("100000")
            .help("Number of simulation samples. Default: 100000.");
    parser.add_option("-p", "--parameter").choices({"concentration", "mean"})
            .dest("parameter")
            .set_default("mean")
            .help("The parameter specified in the first argument. E.g., the "
                  "command "
                  "\'dpprobs -p concentration 2.0 6\' "
                  "would specify that the argument \'2.0\' is the value of the "
                  "concentration parameter, whereas "
                  "\'dpprobs -p mean 2.0 6\' "
                  "specifies that \'2.0\' is the desired mean number of "
                  "categories. "
                  "Default: mean.");
    parser.add_option("--shape")
            .action("store")
            .type("double")
            .dest("shape")
            .set_default("2.0")
            .help("Shape parameter of the gamma distributed prior on the "
                  "concentration parameter of the Dirichlet process. "
                  "If provided, the program will calculate the corresponding "
                  "scale parameter for the gamma distribution such that the "
                  "mean of the gamma prior is equal to the specified value of "
                  "the concentration parameter (or the value the concentration "
                  "parameter that corresponds to the specified mean number of "
                  "categories). If not provided, the concetration parameter "
                  "is simply fixed to the specifed value.");
    parser.add_option("--scale")
            .action("store")
            .type("double")
            .dest("scale")
            .set_default("1.0")
            .help("Scale parameter of the gamma-distributed prior on the "
                  "concentration parameter of the Dirichlet process. "
                  "If provided, the parameter value (the first argument) is "
                  "ignored, and the concentration parameter is treated as a "
                  "gamma-distriuted random variable with shape and scale "
                  "specified by the \'--shape\' and \'--scale\' options.");

    optparse::Values& options = parser.parse_args(argc, argv);
    std::vector<std::string> args = parser.args();

    // Parse number of elements argument
    if (args.size() < 2) {
        throw EcoevolityError("Too few positional arguments; "
                "a parameter value and the number of elements is required.");
    }
    if (args.size() > 2) {
        throw EcoevolityError("Too many positional arguments; "
                "a parameter value and the number of elements is required.");
    }

    double parameter_value = 1.0;
    const std::string parameter_value_arg = args.at(0);
    std::istringstream p_reader(parameter_value_arg);
    p_reader >> parameter_value;

    if (parameter_value <= 0.0) {
        throw EcoevolityError("Parameter value must be positive.");
    }

    unsigned int number_of_elements = 2;
    const std::string number_of_elements_arg = args.at(1);
    std::istringstream n_reader(number_of_elements_arg);
    n_reader >> number_of_elements;
    
    if (number_of_elements < 2) {
        throw EcoevolityError("Number of elements should be greater than 1.");
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

    double concentration;
    double mean_ncats;
    const char* p = options.get("parameter");
    std::string parameter(p);
    if (parameter == "concentration") {
        concentration = parameter_value;
        mean_ncats = get_dpp_expected_number_of_categories(concentration,
                number_of_elements);
    }
    else if (parameter == "mean") {
        mean_ncats = parameter_value;
        concentration = get_dpp_concentration(mean_ncats, number_of_elements);
    }
    else {
        throw EcoevolityError("Unrecognized parameter type: " + parameter);
    }

    bool concentration_is_fixed = true; 
    if (options.is_set_by_user("shape")) {
        concentration_is_fixed = false; 
    }
    if (options.is_set_by_user("scale")) {
        concentration_is_fixed = false; 
    }

    double shape = options.get("shape");
    if (shape <= 0.0) {
        throw EcoevolityError("Shape must be positive.");
    }
    
    double scale = concentration / shape;
    if (options.is_set_by_user("scale")) {
        scale = options.get("scale");
        // Need to override specified parameter value argument
        concentration = shape * scale;
        mean_ncats = get_dpp_expected_number_of_categories(concentration,
                number_of_elements);
    }
    if (scale <= 0.0) {
        throw EcoevolityError("Scale must be positive.");
    }

    GammaDistribution gamma_concentration(shape, scale);

    std::cerr << "Seed = " << seed << std::endl;
    std::cerr << "Number of samples = " << nreps << std::endl;
    std::cerr << "Number of elements = " << number_of_elements << std::endl;
    if (concentration_is_fixed) {
        std::cerr << "Concentration = " << concentration << std::endl;
        std::cerr << "Mean number of categories = " << mean_ncats << std::endl;
    }
    else {
        std::cerr << "Concentration ~ "
                  << gamma_concentration.to_string() << std::endl;
        std::cerr << "Mean number of categories given mean concentration = "
                  << mean_ncats << std::endl;
    }

    std::string number_of_elements_str = std::to_string(number_of_elements);
    unsigned int ncats_padding = number_of_elements_str.size();

    // Start real work
    time_t start;
    time_t finish;
    time(&start);

    std::vector<unsigned int> elements (number_of_elements, 0);
    std::map<unsigned int, unsigned int> number_of_categories_counts;
    for (unsigned int i = 1; i <= number_of_elements; ++i) {
        number_of_categories_counts[i] = 0;
    }
    unsigned int number_of_categories;
    for (unsigned int i = 0; i < nreps; ++i) {
        if (! concentration_is_fixed) {
            concentration = rng.gamma(shape, scale);
        }
        number_of_categories = rng.dirichlet_process(elements, concentration);
        ++number_of_categories_counts[number_of_categories];
    }

    unsigned int tally = 0;
    for (auto const & kv: number_of_categories_counts) {
        tally += kv.second;
    }
    ECOEVOLITY_ASSERT(tally == nreps);

    std::cerr << "\nEstimated probabilities of the number of categories:\n";
    std::cerr << string_util::banner('-') << "\n";
    for (auto const & kv: number_of_categories_counts) {
        std::cout << "p(ncats = "
                  << std::setw(ncats_padding) << std::right << kv.first
                  << ") = "
                  << std::setw(12) << std::left << kv.second / (double)nreps
                  << " (n = " << stirling2_float(number_of_elements, kv.first)
                  << ")\n";
    }
    std::cerr << string_util::banner('-') << "\n\n";

    time(&finish);
    double duration = difftime(finish, start);
    std::cerr << "Runtime: " << duration << " seconds." << std::endl;

    return 0;
}
