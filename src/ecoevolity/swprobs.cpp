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

#include "swprobs.hpp"


void write_swprobs_splash(std::ostream& out) {
    std::string v = "Version ";
    v += PROJECT_DETAILED_VERSION;
    out << string_util::banner('=') << "\n" 
        << string_util::center("swprobs") << "\n"
        << string_util::center("Calculating probabilities under the split-weight model prior") << "\n\n"
        << string_util::center("Part of:") << "\n"
        << string_util::center(PROJECT_NAME) << "\n"
        << string_util::center(v) << "\n"
        << string_util::banner('=') << "\n";
}

int swprobs_main(int argc, char * argv[]) {

    write_swprobs_splash(std::cerr);
    std::cerr << "\n";

    const std::string usage = 
        "usage: %prog [OPTIONS] SPLIT-WEIGHT-VALUE NUMBER-OF-ELEMENTS";

    std::ostringstream version_ss;
    version_ss << PROJECT_NAME << " version " << PROJECT_DETAILED_VERSION;
    const std::string version = version_ss.str();

    const std::string description =
        "swprobs: Simulating split-weight probabilites";
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
            .help("Number of simulation samples. Default: 100000. "
                  "Simulations are only used if a prior distribution on the "
                  "split weight is specified (i.e., using "
                  "\'--shape\'/\'--scale\' options.");
    parser.add_option("--shape")
            .action("store")
            .type("double")
            .dest("shape")
            .set_default("2.0")
            .help("Shape parameter of the gamma-distributed prior on the "
                  "split-weight parameter of the uniform distribution over "
                  "event models. If provided, the program will calculate the "
                  "corresponding scale parameter for the gamma distribution "
                  "such that the mean of the gamma prior is equal to the "
                  "specified value of the split-weight parameter. Also, "
                  "simulations will be performed to approximate the prior "
                  "probabilities for numbers of events. If not provided, the "
                  "split-weight parameter is simply fixed to the specifed "
                  "value, and the prior probabilities are calculated "
                  "exactly.");
    parser.add_option("--scale")
            .action("store")
            .type("double")
            .dest("scale")
            .set_default("1.0")
            .help("Scale parameter of the gamma distributed prior on the "
                  "split-weight parameter of the uniform distribution over "
                  "event models. If provided, the split-weight value provided "
                  "is ignored, and the split-weight parameter is treated as a "
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

    EcoevolityOptions::ModelPrior model_prior = EcoevolityOptions::ModelPrior::uniform;

    // Parse split-weight option
    double split_weight = 1.0;
    const std::string split_weight_arg = args.at(0);
    std::istringstream sw_reader(split_weight_arg);
    sw_reader >> split_weight;

    if (split_weight <= 0.0) {
        throw EcoevolityError("Split weight value must be positive.");
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

    bool split_weight_is_fixed = true; 
    if (options.is_set_by_user("shape")) {
        split_weight_is_fixed = false; 
    }
    if (options.is_set_by_user("scale")) {
        split_weight_is_fixed = false; 
    }

    double shape = options.get("shape");
    if (shape <= 0.0) {
        throw EcoevolityError("Shape must be positive.");
    }
    
    double scale = split_weight / shape;
    if (options.is_set_by_user("scale")) {
        scale = options.get("scale");
        // Need to override specified parameter value argument
    }
    if (scale <= 0.0) {
        throw EcoevolityError("Scale must be positive.");
    }

    GammaDistribution gamma_split_weight(shape, scale);

    std::cerr << "Prior = split-weighted uniform\n";
    std::cerr << "Number of elements = " << number_of_elements << std::endl;
    if (split_weight_is_fixed) {
        std::cerr << "Split weight = " << split_weight << std::endl;
    }
    else {
        std::cerr << "Split weight ~ "
                  << gamma_split_weight.to_string() << std::endl;
        std::cerr << "Seed = " << seed << std::endl;
        std::cerr << "Number of samples = " << nreps << std::endl;
    }
    std::string number_of_elements_str = std::to_string(number_of_elements);
    unsigned int ncats_padding = number_of_elements_str.size();

    // Start real work
    time_t start;
    time_t finish;
    time(&start);

    std::vector<long double> number_of_cat_probs(number_of_elements, 0.0);
    if (split_weight_is_fixed) {
        get_number_of_subset_probs(number_of_cat_probs, split_weight);
        std::vector<double> number_of_cats(number_of_elements, 0.0);
        for (unsigned int i = 0; i < number_of_elements; ++i) {
            number_of_cats.at(i) = i + 1.0;
        }
        double sample_mean_ncats = weighted_mean<double>(number_of_cats, number_of_cat_probs);
        std::cerr << "Mean number of categories = " << sample_mean_ncats << std::endl;

    }
    // else {
    //     unsigned int total = 0;
    //     std::vector<unsigned int> number_of_categories_counts(number_of_elements, 0);
    //     unsigned int number_of_categories;
    //     for (unsigned int i = 0; i < nreps; ++i) {
    //         split_weight = rng.gamma(shape, scale);
    //         number_of_categories = rng.random_number_of_subsets(number_of_elements, split_weight);
    //         ++number_of_categories_counts.at(number_of_categories - 1);
    //         total += number_of_categories;
    //     }
    //     double sample_mean_ncats = (double)total / (double)nreps;
    //     std::cerr << "Sample mean number of categories = " << sample_mean_ncats << std::endl;
    //     unsigned int tally = 0;
    //     double double_nreps = (double)nreps;
    //     for (unsigned int i = 0; i < number_of_elements; ++i) {
    //         tally += number_of_categories_counts.at(i);
    //         number_of_cat_probs.at(i) = number_of_categories_counts.at(i) / double_nreps; 
    //     }
    //     ECOEVOLITY_ASSERT(tally == nreps);
    // }
    //
    // It is much more efficient to only simulate split_weight and then average
    // the probabilities of the number of subsets (below) than it is to
    // simulate both split_weight and the number of subsets and then summarize
    // counts of the latter (commented out above).
    else {
        std::vector<long double> number_of_cat_probs_sum(number_of_elements, 0.0);
        for (unsigned int i = 0; i < nreps; ++i) {
            split_weight = rng.gamma(shape, scale);
            get_number_of_subset_probs(number_of_cat_probs, split_weight);
            for (unsigned int j = 0; j < number_of_elements; ++j) {
                number_of_cat_probs_sum.at(j) += number_of_cat_probs.at(j);
            }
        }
        for (unsigned int j = 0; j < number_of_elements; ++j) {
            number_of_cat_probs.at(j) = number_of_cat_probs_sum.at(j) / nreps;
        }
        std::vector<double> number_of_cats(number_of_elements, 0.0);
        for (unsigned int i = 0; i < number_of_elements; ++i) {
            number_of_cats.at(i) = i + 1.0;
        }
        double sample_mean_ncats = weighted_mean<double>(number_of_cats, number_of_cat_probs);
        std::cerr << "Mean number of categories = " << sample_mean_ncats << std::endl;
    }


    if (split_weight_is_fixed) {
        std::cerr << "\nProbabilities of the number of categories:\n";
    }
    else {
        std::cerr << "\nApproximated probabilities of the number of categories:\n";
    }
    std::cerr << string_util::banner('-') << "\n";
    for (unsigned int i = 0; i < number_of_elements; ++i) {
        std::cout << "p(ncats = "
                  << std::setw(ncats_padding) << std::right << i + 1
                  << ") = "
                  << std::setw(12) << std::left << number_of_cat_probs.at(i)
                  << " (n = " << stirling2_float(number_of_elements, i + 1)
                  << ")\n";
    }
    std::cerr << string_util::banner('-') << "\n\n";

    time(&finish);
    double duration = difftime(finish, start);
    std::cerr << "Runtime: " << duration << " seconds." << std::endl;

    return 0;
}
