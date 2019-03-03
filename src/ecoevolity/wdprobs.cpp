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

#include "wdprobs.hpp"


void write_wdprobs_splash(std::ostream& out) {
    std::string v = "Version ";
    v += PROJECT_DETAILED_VERSION;
    out << string_util::banner('=') << "\n" 
        << string_util::center("WDprobs") << "\n"
        << string_util::center("Simulating weighted-discount process probabilities") << "\n\n"
        << string_util::center("Part of:") << "\n"
        << string_util::center(PROJECT_NAME) << "\n"
        << string_util::center(v) << "\n"
        << string_util::banner('=') << "\n";
}

int wdprobs_main(int argc, char * argv[]) {

    write_wdprobs_splash(std::cerr);
    std::cerr << "\n";

    const std::string usage = 
        "usage: %prog [OPTIONS] PARAMETER-VALUE NUMBER-OF-ELEMENTS";

    std::ostringstream version_ss;
    version_ss << PROJECT_NAME << " version " << PROJECT_DETAILED_VERSION;
    const std::string version = version_ss.str();

    const std::string description =
        "WDprobs: Simulating weighted-discount process probabilites";
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
    parser.add_option("-c", "--concentration")
            .action("store")
            .type("double")
            .dest("concentration")
            .set_default("1.0")
            .help("The value of the concentration parameter of the weighted-discount "
                  "process. ");
    parser.add_option("-d", "--discount")
            .action("store")
            .type("double")
            .dest("discount")
            .set_default("0.0")
            .help("The value of the discount parameter of the weighted-discount "
                  "process. ");
    parser.add_option("--concentration-shape")
            .action("store")
            .type("double")
            .dest("concentration_shape")
            .set_default("2.0")
            .help("Shape parameter of the gamma-distributed prior on the "
                  "concentration parameter of the weighted-discount process.");
    parser.add_option("--concentration-scale")
            .action("store")
            .type("double")
            .dest("concentration_scale")
            .set_default("1.0")
            .help("Scale parameter of the gamma-distributed prior on the "
                  "concentration parameter of the weighted-discount process.");
    parser.add_option("--discount-alpha")
            .action("store")
            .type("double")
            .dest("discount_alpha")
            .set_default("1.0")
            .help("Alpha parameter of the beta-distributed prior on the "
                  "discount parameter of the weighted-discount process.");
    parser.add_option("--discount-beta")
            .action("store")
            .type("double")
            .dest("discount_beta")
            .set_default("1.0")
            .help("Beta parameter of the beta-distributed prior on the "
                  "discount parameter of the weighted-discount process.");

    optparse::Values& options = parser.parse_args(argc, argv);
    std::vector<std::string> args = parser.args();

    // Parse number of elements argument
    if (args.size() < 1) {
        throw EcoevolityError("Too few positional arguments; "
                "a parameter value and the number of elements is required.");
    }
    if (args.size() > 1) {
        throw EcoevolityError("Too many positional arguments; "
                "a parameter value and the number of elements is required.");
    }

    // Parse discount options
    double discount = options.get("discount");
    double discount_mean = discount;
    double discount_alpha = options.get("discount_alpha");
    double discount_beta = options.get("discount_beta");
    bool discount_is_fixed = true; 
    if ((options.is_set_by_user("discount_alpha")) &&
            (options.is_set_by_user("discount_beta"))) {
        discount_is_fixed = false; 
        discount_mean = discount_alpha / (discount_alpha + discount_beta);
    }
    if ((discount < 0.0) || (discount >= 1.0)) {
        std::ostringstream message;
        message << "Invalid discount value ("
                << discount
                << "); it must be 0.0 <= discount < 1.0\n";
        throw EcoevolityError(message.str());
    }
    BetaDistribution beta_dist_discount(discount_alpha, discount_beta);

    double concentration = options.get("concentration");
    double concentration_mean = concentration;
    double concentration_shape = options.get("concentration_shape");
    double concentration_scale = options.get("concentration_scale");
    bool concentration_is_fixed = true; 
    if ((options.is_set_by_user("concentration_shape")) &&
            (options.is_set_by_user("concentration_scale"))) {
        concentration_is_fixed = false; 
        concentration_mean = concentration_shape * concentration_scale;
    }
    if ((concentration <= 0.0)) {
        std::ostringstream message;
        message << "Invalid concentration value ("
                << concentration
                << "); it must be 0.0 <= concentration < 1.0\n";
        throw EcoevolityError(message.str());
    }
    GammaDistribution gamma_concentration(concentration_shape, concentration_scale);

    unsigned int number_of_elements = 2;
    const std::string number_of_elements_arg = args.at(0);
    std::istringstream n_reader(number_of_elements_arg);
    n_reader >> number_of_elements;
    
    if (number_of_elements < 2) {
        throw EcoevolityError("Number of elements must be greater than 1.");
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

    std::cerr << "Seed = " << seed << std::endl;
    std::cerr << "Number of samples = " << nreps << std::endl;
    std::cerr << "Number of elements = " << number_of_elements << std::endl;
    if (concentration_is_fixed) {
        std::cerr << "Concentration = " << concentration << std::endl;
    }
    else {
        std::cerr << "Concentration ~ "
                  << gamma_concentration.to_string() << std::endl;
    }
    if (discount_is_fixed) {
        std::cerr << "Discount = " << discount << std::endl;
    }
    else {
        std::cerr << "Discount ~ "
                  << beta_dist_discount.to_string() << std::endl;
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
            concentration = rng.gamma(concentration_shape, concentration_scale);
        }
        if (! discount_is_fixed) {
            discount = rng.beta(discount_alpha, discount_beta);
        }
        number_of_categories = rng.weighted_discount_process(elements, concentration, discount);
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
    std::vector<unsigned int> part1(number_of_elements, 0);
    std::vector<unsigned int> part2 = part1;
    for (unsigned int i = 0; i < number_of_elements; ++i) {
        part2.at(i) = i;
    }

    double part1_ln_prob = get_wdp_log_prior_probability<unsigned int>(part1, concentration, discount);
    double part2_ln_prob = get_wdp_log_prior_probability<unsigned int>(part2, concentration, discount);
    double part1_prob = std::exp(part1_ln_prob);
    double part2_prob = std::exp(part2_ln_prob);
    std::cout << part1_prob << "\n";
    std::cout << part2_prob << "\n";

    std::cerr << string_util::banner('-') << "\n\n";

    time(&finish);
    double duration = difftime(finish, start);
    std::cerr << "Runtime: " << duration << " seconds." << std::endl;

    return 0;
}
