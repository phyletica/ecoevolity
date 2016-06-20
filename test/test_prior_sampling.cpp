#include "catch.hpp"
#include "ecoevolity/ecoevolity.hpp"

#include "ecoevolity/rng.hpp"
#include "ecoevolity/path.hpp"
#include "ecoevolity/spreadsheet.hpp"

RandomNumberGenerator _PRIOR_SAMPLING_RNG = RandomNumberGenerator();

TEST_CASE("Testing sampling from prior with ComparisonHeightScaler", "[SamplingPrior]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        double shape = 10.0;
        double scale = 0.1;
        std::string auto_optimize = "false";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << shape << "\n";
        os << "        scale: " << scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 0.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: true\n";
        os << "    constrain_mutation_rates: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        u_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        time_multiplier:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "1234";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);
        std::vector<double> samples = prior_sample.get<double>("root_height_kya");
        REQUIRE(samples.size() == 10001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("root_height_kya");
        REQUIRE(summary.mean() == Approx(shape * scale).epsilon(0.01));
        REQUIRE(summary.variance() == Approx(shape * scale * scale).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(summary.mean() == 0.0);
        REQUIRE(summary.variance() == 0.0);
        summary = prior_sample.summarize<double>("ln_likelihood_kya");
        REQUIRE(summary.mean() == 0.0);
        REQUIRE(summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "root_height_kya",
                "ln_prior",
                "ln_prior_kya"};
        for (auto const &kv: prior_sample.data) {
            bool test = true;
            for (auto const &p: columns_to_ignore) {
                if (kv.first == p) {
                    test = false;
                }
            }
            if (test) {
                summary = prior_sample.summarize<double>(kv.first);
                REQUIRE(summary.variance() == 0.0);
            }
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing sampling from prior with ComparisonHeightScaler with optimizing",
        "[SamplingPrior]") {

    SECTION("Testing gamma(10.0, 0.1) prior and optimizing") {
        double shape = 10.0;
        double scale = 0.1;
        std::string auto_optimize = "false";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << shape << "\n";
        os << "        scale: " << scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 0.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: true\n";
        os << "    constrain_mutation_rates: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        u_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        time_multiplier:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "4321";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);
        std::vector<double> samples = prior_sample.get<double>("root_height_kya");
        REQUIRE(samples.size() == 10001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("root_height_kya");
        REQUIRE(summary.mean() == Approx(shape * scale).epsilon(0.01));
        REQUIRE(summary.variance() == Approx(shape * scale * scale).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(summary.mean() == 0.0);
        REQUIRE(summary.variance() == 0.0);
        summary = prior_sample.summarize<double>("ln_likelihood_kya");
        REQUIRE(summary.mean() == 0.0);
        REQUIRE(summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "root_height_kya",
                "ln_prior",
                "ln_prior_kya"};
        for (auto const &kv: prior_sample.data) {
            bool test = true;
            for (auto const &p: columns_to_ignore) {
                if (kv.first == p) {
                    test = false;
                }
            }
            if (test) {
                summary = prior_sample.summarize<double>(kv.first);
                REQUIRE(summary.variance() == 0.0);
            }
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing sampling from prior with RootPopulationSizeScaler", "[SamplingPrior]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        double shape = 10.0;
        double scale = 0.1;
        std::string auto_optimize = "false";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: 10.0\n";
        os << "        scale: 0.001\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 0.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: false\n";
        os << "    constrain_mutation_rates: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << shape << "\n";
        os << "                    scale: " << scale << "\n";
        os << "        u_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        time_multiplier:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "1111";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);
        std::vector<double> samples = prior_sample.get<double>("pop_size_root_kya");
        REQUIRE(samples.size() == 10001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("pop_size_root_kya");
        REQUIRE(summary.mean() == Approx(shape * scale).epsilon(0.01));
        REQUIRE(summary.variance() == Approx(shape * scale * scale).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(summary.mean() == 0.0);
        REQUIRE(summary.variance() == 0.0);
        summary = prior_sample.summarize<double>("ln_likelihood_kya");
        REQUIRE(summary.mean() == 0.0);
        REQUIRE(summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {"generation",
            "ln_prior",
            "ln_prior_kya",
            "pop_size_root_kya"};
        for (auto const &kv: prior_sample.data) {
            bool test = true;
            for (auto const &p: columns_to_ignore) {
                if (kv.first == p) {
                    test = false;
                }
            }
            if (test) {
                summary = prior_sample.summarize<double>(kv.first);
                REQUIRE(summary.variance() == 0.0);
            }
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing sampling from prior with RootPopulationSizeScaler with optimizing", "[SamplingPrior]") {

    SECTION("Testing gamma(10.0, 0.1) prior and optimizing") {
        double shape = 10.0;
        double scale = 0.1;
        std::string auto_optimize = "true";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: 10.0\n";
        os << "        scale: 0.001\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 0.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: false\n";
        os << "    constrain_mutation_rates: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << shape << "\n";
        os << "                    scale: " << scale << "\n";
        os << "        u_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        time_multiplier:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "2222";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);
        std::vector<double> samples = prior_sample.get<double>("pop_size_root_kya");
        REQUIRE(samples.size() == 10001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("pop_size_root_kya");
        REQUIRE(summary.mean() == Approx(shape * scale).epsilon(0.01));
        REQUIRE(summary.variance() == Approx(shape * scale * scale).epsilon(0.01));

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "ln_prior",
                "ln_prior_kya",
                "pop_size_root_kya"};
        for (auto const &kv: prior_sample.data) {
            bool test = true;
            for (auto const &p: columns_to_ignore) {
                if (kv.first == p) {
                    test = false;
                }
            }
            if (test) {
                summary = prior_sample.summarize<double>(kv.first);
                REQUIRE(summary.variance() == 0.0);
            }
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing sampling from prior with ChildPopulationSizeScaler", "[SamplingPrior]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        double shape = 10.0;
        double scale = 0.1;
        std::string auto_optimize = "false";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: 10.0\n";
        os << "        scale: 0.001\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 0.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: false\n";
        os << "    constrain_mutation_rates: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << shape << "\n";
        os << "                    scale: " << scale << "\n";
        os << "        u_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        time_multiplier:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "1111";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);
        std::vector<double> samples1 = prior_sample.get<double>("pop_size_kya");
        REQUIRE(samples1.size() == 10001);
        std::vector<double> samples2 = prior_sample.get<double>("pop_size_fas");
        REQUIRE(samples2.size() == 10001);

        for (size_t i = 0; i < samples1.size(); ++i) {
            if (i > 0) {
                REQUIRE(samples1.at(i) != samples2.at(i));
            }
        }

        SampleSummarizer<double> summary1 = prior_sample.summarize<double>("pop_size_kya");
        REQUIRE(summary1.mean() == Approx(shape * scale).epsilon(0.01));
        REQUIRE(summary1.variance() == Approx(shape * scale * scale).epsilon(0.01));

        SampleSummarizer<double> summary2 = prior_sample.summarize<double>("pop_size_fas");
        REQUIRE(summary2.mean() == Approx(shape * scale).epsilon(0.01));
        REQUIRE(summary2.variance() == Approx(shape * scale * scale).epsilon(0.01));

        SampleSummarizer<double> summary;
        std::vector<std::string> columns_to_ignore = {
                "generation",
                "ln_prior",
                "ln_prior_kya",
                "pop_size_kya",
                "pop_size_fas"};
        for (auto const &kv: prior_sample.data) {
            bool test = true;
            for (auto const &p: columns_to_ignore) {
                if (kv.first == p) {
                    test = false;
                }
            }
            if (test) {
                summary = prior_sample.summarize<double>(kv.first);
                REQUIRE(summary.variance() == 0.0);
            }
        }

        delete[] cfg_path;
    }
}


TEST_CASE("Testing sampling from prior with ChildPopulationSizeScaler with optimizing", "[SamplingPrior]") {

    SECTION("Testing gamma(10.0, 0.1) prior and optimizing") {
        double shape = 10.0;
        double scale = 0.1;
        std::string auto_optimize = "true";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: 10.0\n";
        os << "        scale: 0.001\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 0.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: false\n";
        os << "    constrain_mutation_rates: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << shape << "\n";
        os << "                    scale: " << scale << "\n";
        os << "        u_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        time_multiplier:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "2222";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);
        std::vector<double> samples1 = prior_sample.get<double>("pop_size_kya");
        REQUIRE(samples1.size() == 10001);
        std::vector<double> samples2 = prior_sample.get<double>("pop_size_fas");
        REQUIRE(samples2.size() == 10001);

        for (size_t i = 0; i < samples1.size(); ++i) {
            if (i > 0) {
                REQUIRE(samples1.at(i) != samples2.at(i));
            }
        }

        SampleSummarizer<double> summary1 = prior_sample.summarize<double>("pop_size_kya");
        REQUIRE(summary1.mean() == Approx(shape * scale).epsilon(0.01));
        REQUIRE(summary1.variance() == Approx(shape * scale * scale).epsilon(0.01));

        SampleSummarizer<double> summary2 = prior_sample.summarize<double>("pop_size_fas");
        REQUIRE(summary2.mean() == Approx(shape * scale).epsilon(0.01));
        REQUIRE(summary2.variance() == Approx(shape * scale * scale).epsilon(0.01));

        SampleSummarizer<double> summary;
        std::vector<std::string> columns_to_ignore = {
                "generation",
                "ln_prior",
                "ln_prior_kya",
                "pop_size_kya",
                "pop_size_fas"};
        for (auto const &kv: prior_sample.data) {
            bool test = true;
            for (auto const &p: columns_to_ignore) {
                if (kv.first == p) {
                    test = false;
                }
            }
            if (test) {
                summary = prior_sample.summarize<double>(kv.first);
                REQUIRE(summary.variance() == 0.0);
            }
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing sampling from prior with RootPopulationSizeScaler and ChildPopulationSizeScaler on constrained sizes",
        "[SamplingPrior]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        double shape = 10.0;
        double scale = 0.1;
        std::string auto_optimize = "false";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: 10.0\n";
        os << "        scale: 0.001\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 0.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: true\n";
        os << "    constrain_mutation_rates: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << shape << "\n";
        os << "                    scale: " << scale << "\n";
        os << "        u_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        time_multiplier:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "1111";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);
        std::vector<double> samples1 = prior_sample.get<double>("pop_size_kya");
        REQUIRE(samples1.size() == 10001);
        std::vector<double> samples2 = prior_sample.get<double>("pop_size_fas");
        REQUIRE(samples2.size() == 10001);
        std::vector<double> samples3 = prior_sample.get<double>("pop_size_root_kya");
        REQUIRE(samples2.size() == 10001);

        for (size_t i = 0; i < samples1.size(); ++i) {
            REQUIRE(samples1.at(i) == samples2.at(i));
            REQUIRE(samples1.at(i) == samples3.at(i));
        }

        SampleSummarizer<double> summary = prior_sample.summarize<double>("pop_size_kya");
        REQUIRE(summary.mean() == Approx(shape * scale).epsilon(0.01));
        REQUIRE(summary.variance() == Approx(shape * scale * scale).epsilon(0.01));

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "ln_prior",
                "ln_prior_kya",
                "pop_size_kya",
                "pop_size_fas",
                "pop_size_root_kya"};
        for (auto const &kv: prior_sample.data) {
            bool test = true;
            for (auto const &p: columns_to_ignore) {
                if (kv.first == p) {
                    test = false;
                }
            }
            if (test) {
                summary = prior_sample.summarize<double>(kv.first);
                REQUIRE(summary.variance() == 0.0);
            }
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing sampling from prior with RootPopulationSizeScaler and ChildPopulationSizeScaler on constrained sizes with optimizing",
        "[SamplingPrior]") {

    SECTION("Testing gamma(10.0, 0.1) prior and optimizing") {
        double shape = 10.0;
        double scale = 0.1;
        std::string auto_optimize = "true";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: 10.0\n";
        os << "        scale: 0.001\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 0.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: true\n";
        os << "    constrain_mutation_rates: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << shape << "\n";
        os << "                    scale: " << scale << "\n";
        os << "        u_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        time_multiplier:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "1111";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);
        std::vector<double> samples1 = prior_sample.get<double>("pop_size_kya");
        REQUIRE(samples1.size() == 10001);
        std::vector<double> samples2 = prior_sample.get<double>("pop_size_fas");
        REQUIRE(samples2.size() == 10001);
        std::vector<double> samples3 = prior_sample.get<double>("pop_size_root_kya");
        REQUIRE(samples2.size() == 10001);

        for (size_t i = 0; i < samples1.size(); ++i) {
            REQUIRE(samples1.at(i) == samples2.at(i));
            REQUIRE(samples1.at(i) == samples3.at(i));
        }

        SampleSummarizer<double> summary = prior_sample.summarize<double>("pop_size_kya");
        REQUIRE(summary.mean() == Approx(shape * scale).epsilon(0.01));
        REQUIRE(summary.variance() == Approx(shape * scale * scale).epsilon(0.01));

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "ln_prior",
                "ln_prior_kya",
                "pop_size_kya",
                "pop_size_fas",
                "pop_size_root_kya"};
        for (auto const &kv: prior_sample.data) {
            bool test = true;
            for (auto const &p: columns_to_ignore) {
                if (kv.first == p) {
                    test = false;
                }
            }
            if (test) {
                summary = prior_sample.summarize<double>(kv.first);
                REQUIRE(summary.variance() == 0.0);
            }
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing sampling from prior with MutationRateScaler without offset",
        "[SamplingPrior]") {

    SECTION("Testing gamma(1.0, 1.0) prior and no optimizing") {
        double shape = 1.0;
        double scale = 1.0;
        std::string auto_optimize = "false";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: 10.0\n";
        os << "        scale: 0.001\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 0.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: true\n";
        os << "    constrain_mutation_rates: false\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        u_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << shape << "\n";
        os << "                    scale: " << scale << "\n";
        os << "        time_multiplier:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "1234";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);
        std::vector<double> samples = prior_sample.get<double>("u_kya");
        REQUIRE(samples.size() == 10001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("u_kya");
        // This should *not* match prior, because u_rate is bounded to be
        // greater than 0.5
        REQUIRE(summary.mean() != Approx(shape * scale).epsilon(0.01));
        REQUIRE(summary.variance() != Approx(shape * scale * scale).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(summary.mean() == 0.0);
        REQUIRE(summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "ln_prior",
                "ln_prior_kya",
                "u_kya",
                "v_kya"};
        for (auto const &kv: prior_sample.data) {
            bool test = true;
            for (auto const &p: columns_to_ignore) {
                if (kv.first == p) {
                    test = false;
                }
            }
            if (test) {
                summary = prior_sample.summarize<double>(kv.first);
                REQUIRE(summary.variance() == 0.0);
            }
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing sampling from prior with MutationRateScaler with offset",
        "[SamplingPrior]") {

    SECTION("Testing gamma(1.0, 0.5, offset=0.5) prior and no optimizing") {
        double shape = 1.0;
        double scale = 0.5;
        double offset = 0.5;
        std::string auto_optimize = "false";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: 10.0\n";
        os << "        scale: 0.001\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 0.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: true\n";
        os << "    constrain_mutation_rates: false\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        u_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    offset: " << offset << "\n";
        os << "                    shape: " << shape << "\n";
        os << "                    scale: " << scale << "\n";
        os << "        time_multiplier:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "1234";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);
        std::vector<double> samples = prior_sample.get<double>("u_kya");
        REQUIRE(samples.size() == 10001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("u_kya");
        // This should match prior, because offset matches the bound on u_rate to be
        // greater than 0.5.
        REQUIRE(summary.mean() == Approx((shape * scale) + offset).epsilon(0.01));
        REQUIRE(summary.variance() == Approx(shape * scale * scale).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(summary.mean() == 0.0);
        REQUIRE(summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "ln_prior",
                "ln_prior_kya",
                "u_kya",
                "v_kya"};
        for (auto const &kv: prior_sample.data) {
            bool test = true;
            for (auto const &p: columns_to_ignore) {
                if (kv.first == p) {
                    test = false;
                }
            }
            if (test) {
                summary = prior_sample.summarize<double>(kv.first);
                REQUIRE(summary.variance() == 0.0);
            }
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing sampling from prior with MutationRateScaler with offset and optimizing",
        "[SamplingPrior]") {

    SECTION("Testing gamma(1.0, 0.5, offset=0.5) prior and optimizing") {
        double shape = 1.0;
        double scale = 0.5;
        double offset = 0.5;
        std::string auto_optimize = "true";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: 10.0\n";
        os << "        scale: 0.001\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 0.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: true\n";
        os << "    constrain_mutation_rates: false\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        u_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    offset: " << offset << "\n";
        os << "                    shape: " << shape << "\n";
        os << "                    scale: " << scale << "\n";
        os << "        time_multiplier:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "1234";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("u_kya");
        REQUIRE(summary.sample_size() == 10001);
        // This should match prior, because offset matches the bound on u_rate to be
        // greater than 0.5.
        REQUIRE(summary.mean() == Approx((shape * scale) + offset).epsilon(0.01));
        REQUIRE(summary.variance() == Approx(shape * scale * scale).epsilon(0.01));

        std::vector<double> u_sample = prior_sample.get<double>("u_kya");
        std::vector<double> v_sample = prior_sample.get<double>("v_kya");
        for (size_t i = 0; i < u_sample.size(); ++i) {
            double u = u_sample.at(i);
            double v = v_sample.at(i);
            REQUIRE(v == Approx(u / ((2.0 * u) - 1.0)));
        }

        // Make sure the rest of the prior sample is as expected
        summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(summary.mean() == 0.0);
        REQUIRE(summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "ln_prior",
                "ln_prior_kya",
                "u_kya",
                "v_kya"};
        for (auto const &kv: prior_sample.data) {
            bool test = true;
            for (auto const &p: columns_to_ignore) {
                if (kv.first == p) {
                    test = false;
                }
            }
            if (test) {
                summary = prior_sample.summarize<double>(kv.first);
                REQUIRE(summary.variance() == 0.0);
            }
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing sampling from prior with ComparisonHeightMultiplierScaler",
        "[SamplingPrior]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        double shape = 10.0;
        double scale = 0.1;
        std::string auto_optimize = "false";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: 10.0\n";
        os << "        scale: 0.001\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 0.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: false\n";
        os << "    constrain_mutation_rates: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.001\n";
        os << "            estimate: false\n";
        os << "        u_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        time_multiplier:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << shape << "\n";
        os << "                    scale: " << scale << "\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "7777";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);
        std::vector<double> samples = prior_sample.get<double>("time_multiplier_kya");
        REQUIRE(samples.size() == 10001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("time_multiplier_kya");
        REQUIRE(summary.mean() == Approx(shape * scale).epsilon(0.01));
        REQUIRE(summary.variance() == Approx(shape * scale * scale).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(summary.mean() == 0.0);
        REQUIRE(summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "ln_prior",
                "ln_prior_kya",
                "time_multiplier_kya"};
        for (auto const &kv: prior_sample.data) {
            bool test = true;
            for (auto const &p: columns_to_ignore) {
                if (kv.first == p) {
                    test = false;
                }
            }
            if (test) {
                summary = prior_sample.summarize<double>(kv.first);
                REQUIRE(summary.variance() == 0.0);
            }
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing sampling from prior with ComparisonHeightMultiplierScaler with optimizing",
        "[SamplingPrior]") {

    SECTION("Testing gamma(10.0, 0.1) prior and optimizing") {
        double shape = 10.0;
        double scale = 0.1;
        std::string auto_optimize = "true";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: 10.0\n";
        os << "        scale: 0.001\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 0.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: false\n";
        os << "    constrain_mutation_rates: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.001\n";
        os << "            estimate: false\n";
        os << "        u_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        time_multiplier:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << shape << "\n";
        os << "                    scale: " << scale << "\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "7777";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);
        std::vector<double> samples = prior_sample.get<double>("time_multiplier_kya");
        REQUIRE(samples.size() == 10001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("time_multiplier_kya");
        REQUIRE(summary.mean() == Approx(shape * scale).epsilon(0.01));
        REQUIRE(summary.variance() == Approx(shape * scale * scale).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(summary.mean() == 0.0);
        REQUIRE(summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "ln_prior",
                "ln_prior_kya",
                "time_multiplier_kya"};
        for (auto const &kv: prior_sample.data) {
            bool test = true;
            for (auto const &p: columns_to_ignore) {
                if (kv.first == p) {
                    test = false;
                }
            }
            if (test) {
                summary = prior_sample.summarize<double>(kv.first);
                REQUIRE(summary.variance() == 0.0);
            }
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing fully parameterized model for one pair",
        "[SamplingPrior]") {

    SECTION("Testing with no optimizing") {
        double height_shape = 1.0;
        double height_scale = 0.1;
        double size_shape = 10.0;
        double size_scale = 0.0001;
        double u_shape = 1.0;
        double u_scale = 0.5;
        double u_offset = 0.5;
        double mult_shape = 10.0;
        double mult_scale = 0.1;
        std::string auto_optimize = "false";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 500000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 0.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: false\n";
        os << "    constrain_mutation_rates: false\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shape << "\n";
        os << "                    scale: " << size_scale << "\n";
        os << "        u_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    offset: " << u_offset << "\n";
        os << "                    shape: " << u_shape << "\n";
        os << "                    scale: " << u_scale << "\n";
        os << "        time_multiplier:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult_shape << "\n";
        os << "                    scale: " << mult_scale << "\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "9876";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);

        SampleSummarizer<double> height_summary = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> size_summary1 = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size_summary2 = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size_summary3 = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> u_summary = prior_sample.summarize<double>("u_kya");
        SampleSummarizer<double> mult_summary = prior_sample.summarize<double>("time_multiplier_kya");

        REQUIRE(height_summary.sample_size() == 50001);
        REQUIRE(size_summary1.sample_size() == 50001);
        REQUIRE(size_summary2.sample_size() == 50001);
        REQUIRE(size_summary3.sample_size() == 50001);
        REQUIRE(u_summary.sample_size() == 50001);
        REQUIRE(mult_summary.sample_size() == 50001);

        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(size_summary1.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary1.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(u_summary.mean() == Approx((u_shape * u_scale) + u_offset).epsilon(0.01));
        REQUIRE(u_summary.variance() == Approx(u_shape * u_scale * u_scale).epsilon(0.03));
        REQUIRE(mult_summary.mean() == Approx(mult_shape * mult_scale).epsilon(0.01));
        REQUIRE(mult_summary.variance() == Approx(mult_shape * mult_scale * mult_scale).epsilon(0.01));

        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.sample_size() == 50001);
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        std::vector<double> u_sample = prior_sample.get<double>("u_kya");
        std::vector<double> v_sample = prior_sample.get<double>("v_kya");
        for (size_t i = 0; i < u_sample.size(); ++i) {
            double u = u_sample.at(i);
            double v = v_sample.at(i);
            REQUIRE(v == Approx(u / ((2.0 * u) - 1.0)));
        }

        std::vector<double> size_sample1 = prior_sample.get<double>("pop_size_kya");
        std::vector<double> size_sample2 = prior_sample.get<double>("pop_size_fas");
        std::vector<double> size_sample3 = prior_sample.get<double>("pop_size_root_kya");
        for (size_t i = 0; i < size_sample1.size(); ++i) {
            if (i > 99) {
                REQUIRE(size_sample1.at(i) != size_sample2.at(i));
                REQUIRE(size_sample1.at(i) != size_sample3.at(i));
                REQUIRE(size_sample2.at(i) != size_sample3.at(i));
            }
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing fully parameterized model for one pair with optimization",
        "[SamplingPrior]") {

    SECTION("Testing with optimizing") {
        double height_shape = 1.0;
        double height_scale = 0.1;
        double size_shape = 10.0;
        double size_scale = 0.0001;
        double u_shape = 1.0;
        double u_scale = 0.5;
        double u_offset = 0.5;
        double mult_shape = 10.0;
        double mult_scale = 0.1;
        std::string auto_optimize = "true";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 500000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 0.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: false\n";
        os << "    constrain_mutation_rates: false\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shape << "\n";
        os << "                    scale: " << size_scale << "\n";
        os << "        u_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    offset: " << u_offset << "\n";
        os << "                    shape: " << u_shape << "\n";
        os << "                    scale: " << u_scale << "\n";
        os << "        time_multiplier:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult_shape << "\n";
        os << "                    scale: " << mult_scale << "\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "9876";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);

        SampleSummarizer<double> height_summary = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> size_summary1 = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size_summary2 = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size_summary3 = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> u_summary = prior_sample.summarize<double>("u_kya");
        SampleSummarizer<double> mult_summary = prior_sample.summarize<double>("time_multiplier_kya");

        REQUIRE(height_summary.sample_size() == 50001);
        REQUIRE(size_summary1.sample_size() == 50001);
        REQUIRE(size_summary2.sample_size() == 50001);
        REQUIRE(size_summary3.sample_size() == 50001);
        REQUIRE(u_summary.sample_size() == 50001);
        REQUIRE(mult_summary.sample_size() == 50001);

        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(size_summary1.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary1.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(u_summary.mean() == Approx((u_shape * u_scale) + u_offset).epsilon(0.01));
        REQUIRE(u_summary.variance() == Approx(u_shape * u_scale * u_scale).epsilon(0.01));
        REQUIRE(mult_summary.mean() == Approx(mult_shape * mult_scale).epsilon(0.01));
        REQUIRE(mult_summary.variance() == Approx(mult_shape * mult_scale * mult_scale).epsilon(0.01));

        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.sample_size() == 50001);
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        std::vector<double> u_sample = prior_sample.get<double>("u_kya");
        std::vector<double> v_sample = prior_sample.get<double>("v_kya");
        for (size_t i = 0; i < u_sample.size(); ++i) {
            double u = u_sample.at(i);
            double v = v_sample.at(i);
            REQUIRE(v == Approx(u / ((2.0 * u) - 1.0)));
        }

        std::vector<double> size_sample1 = prior_sample.get<double>("pop_size_kya");
        std::vector<double> size_sample2 = prior_sample.get<double>("pop_size_fas");
        std::vector<double> size_sample3 = prior_sample.get<double>("pop_size_root_kya");
        for (size_t i = 0; i < size_sample1.size(); ++i) {
            if (i > 99) {
                REQUIRE(size_sample1.at(i) != size_sample2.at(i));
                REQUIRE(size_sample1.at(i) != size_sample3.at(i));
                REQUIRE(size_sample2.at(i) != size_sample3.at(i));
            }
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing DPP with 2 pairs and alpha 1.0", "[SamplingPrior]") {

    SECTION("Testing alpha fixed to 1.0") {
        double height_shape = 10.0;
        double height_scale = 0.1;
        std::string auto_optimize = "true";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                value: 1.0\n";
        os << "                estimate: false\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: true\n";
        os << "    constrain_mutation_rates: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        u_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        time_multiplier:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "1234";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        REQUIRE(height_summary1.sample_size() == 10001);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == 10001);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        std::vector<int> nevents = prior_sample.get<int>("number_of_events");
        std::vector<int> event_indices1 = prior_sample.get<int>("root_height_index_kya");
        std::vector<int> event_indices2 = prior_sample.get<int>("root_height_index_pop1");
        std::vector<double> heights1 = prior_sample.get<double>("root_height_kya");
        std::vector<double> heights2 = prior_sample.get<double>("root_height_pop1");

        int nshared = 0;
        int total = 0;
        for (size_t i = 0; i < nevents.size(); ++i) {
            REQUIRE((nevents.at(i) == 1 || nevents.at(i) == 2));
            REQUIRE((event_indices1.at(i) == 0 || event_indices1.at(i) == 1));
            REQUIRE((event_indices2.at(i) == 0 || event_indices2.at(i) == 1));
            if (nevents.at(i) == 1) {
                ++nshared;
                REQUIRE(event_indices1.at(i) == event_indices2.at(i));
                REQUIRE(heights1.at(i) == heights2.at(i));
            }
            else {
                REQUIRE(event_indices1.at(i) != event_indices2.at(i));
                REQUIRE(heights1.at(i) != heights2.at(i));
            }
            ++total;
        }
        REQUIRE(total == 10001);
        REQUIRE((nshared / 10001.0) == Approx(0.5).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "ln_prior",
                "ln_prior_kya",
                "ln_prior_pop1",
                "number_of_events",
                "root_height_kya",
                "root_height_pop1",
                "root_height_index_kya",
                "root_height_index_pop1"
        };
        SampleSummarizer<double> summary;
        for (auto const &kv: prior_sample.data) {
            bool test = true;
            for (auto const &p: columns_to_ignore) {
                if (kv.first == p) {
                    test = false;
                }
            }
            if (test) {
                summary = prior_sample.summarize<double>(kv.first);
                REQUIRE(summary.variance() == 0.0);
            }
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing DPP with 2 pairs and alpha 2.0", "[SamplingPrior]") {

    SECTION("Testing alpha fixed to 2.0") {
        double height_shape = 10.0;
        double height_scale = 0.1;
        std::string auto_optimize = "true";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                value: 2.0\n";
        os << "                estimate: false\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: true\n";
        os << "    constrain_mutation_rates: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        u_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        time_multiplier:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "1234";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        REQUIRE(height_summary1.sample_size() == 10001);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == 10001);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        std::vector<int> nevents = prior_sample.get<int>("number_of_events");
        std::vector<int> event_indices1 = prior_sample.get<int>("root_height_index_kya");
        std::vector<int> event_indices2 = prior_sample.get<int>("root_height_index_pop1");
        std::vector<double> heights1 = prior_sample.get<double>("root_height_kya");
        std::vector<double> heights2 = prior_sample.get<double>("root_height_pop1");

        int nshared = 0;
        int total = 0;
        for (size_t i = 0; i < nevents.size(); ++i) {
            REQUIRE((nevents.at(i) == 1 || nevents.at(i) == 2));
            REQUIRE((event_indices1.at(i) == 0 || event_indices1.at(i) == 1));
            REQUIRE((event_indices2.at(i) == 0 || event_indices2.at(i) == 1));
            if (nevents.at(i) == 1) {
                ++nshared;
                REQUIRE(event_indices1.at(i) == event_indices2.at(i));
                REQUIRE(heights1.at(i) == heights2.at(i));
            }
            else {
                REQUIRE(event_indices1.at(i) != event_indices2.at(i));
                REQUIRE(heights1.at(i) != heights2.at(i));
            }
            ++total;
        }
        REQUIRE(total == 10001);
        REQUIRE((nshared / 10001.0) == Approx(1.0/3.0).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "ln_prior",
                "ln_prior_kya",
                "ln_prior_pop1",
                "number_of_events",
                "root_height_kya",
                "root_height_pop1",
                "root_height_index_kya",
                "root_height_index_pop1"
        };
        SampleSummarizer<double> summary;
        for (auto const &kv: prior_sample.data) {
            bool test = true;
            for (auto const &p: columns_to_ignore) {
                if (kv.first == p) {
                    test = false;
                }
            }
            if (test) {
                summary = prior_sample.summarize<double>(kv.first);
                REQUIRE(summary.variance() == 0.0);
            }
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing DPP with 2 pairs and alpha 0.5", "[SamplingPrior]") {

    SECTION("Testing alpha fixed to 0.5") {
        double height_shape = 10.0;
        double height_scale = 0.1;
        std::string auto_optimize = "true";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                value: 0.5\n";
        os << "                estimate: false\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: true\n";
        os << "    constrain_mutation_rates: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        u_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        time_multiplier:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "1234";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        REQUIRE(height_summary1.sample_size() == 10001);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == 10001);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        std::vector<int> nevents = prior_sample.get<int>("number_of_events");
        std::vector<int> event_indices1 = prior_sample.get<int>("root_height_index_kya");
        std::vector<int> event_indices2 = prior_sample.get<int>("root_height_index_pop1");
        std::vector<double> heights1 = prior_sample.get<double>("root_height_kya");
        std::vector<double> heights2 = prior_sample.get<double>("root_height_pop1");

        int nshared = 0;
        int total = 0;
        for (size_t i = 0; i < nevents.size(); ++i) {
            REQUIRE((nevents.at(i) == 1 || nevents.at(i) == 2));
            REQUIRE((event_indices1.at(i) == 0 || event_indices1.at(i) == 1));
            REQUIRE((event_indices2.at(i) == 0 || event_indices2.at(i) == 1));
            if (nevents.at(i) == 1) {
                ++nshared;
                REQUIRE(event_indices1.at(i) == event_indices2.at(i));
                REQUIRE(heights1.at(i) == heights2.at(i));
            }
            else {
                REQUIRE(event_indices1.at(i) != event_indices2.at(i));
                REQUIRE(heights1.at(i) != heights2.at(i));
            }
            ++total;
        }
        REQUIRE(total == 10001);
        REQUIRE((nshared / 10001.0) == Approx(2.0/3.0).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "ln_prior",
                "ln_prior_kya",
                "ln_prior_pop1",
                "number_of_events",
                "root_height_kya",
                "root_height_pop1",
                "root_height_index_kya",
                "root_height_index_pop1"
        };
        SampleSummarizer<double> summary;
        for (auto const &kv: prior_sample.data) {
            bool test = true;
            for (auto const &p: columns_to_ignore) {
                if (kv.first == p) {
                    test = false;
                }
            }
            if (test) {
                summary = prior_sample.summarize<double>(kv.first);
                REQUIRE(summary.variance() == 0.0);
            }
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing DPP with 3 pairs and alpha 1.0", "[SamplingPrior]") {

    SECTION("Testing alpha fixed to 1.0") {
        double height_shape = 10.0;
        double height_scale = 0.1;
        std::string auto_optimize = "true";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                value: 1.0\n";
        os << "                estimate: false\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightMultiplierScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \"_\"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_mutation_rate_starting_values: false\n";
        os << "    constrain_population_sizes: true\n";
        os << "    constrain_mutation_rates: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        u_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        time_multiplier:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "1234";
        char arg3[] = "--ignore-data";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        ret = ecoevolity_main(argc, argv);
        REQUIRE(ret == 0);
        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");
        REQUIRE(height_summary1.sample_size() == 10001);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == 10001);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == 10001);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        std::vector<int> nevents = prior_sample.get<int>("number_of_events");
        std::vector<int> event_indices1 = prior_sample.get<int>("root_height_index_kya");
        std::vector<int> event_indices2 = prior_sample.get<int>("root_height_index_pop1");
        std::vector<int> event_indices3 = prior_sample.get<int>("root_height_index_pop1b");
        std::vector<double> heights1 = prior_sample.get<double>("root_height_kya");
        std::vector<double> heights2 = prior_sample.get<double>("root_height_pop1");
        std::vector<double> heights3 = prior_sample.get<double>("root_height_pop1b");

        std::map<std::string, int> model_counts = {
                {"000", 0},
                {"001", 0},
                {"010", 0},
                {"011", 0},
                {"012", 0}
        };
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0}
        };
        for (size_t i = 0; i < nevents.size(); ++i) {
            std::ostringstream stream;
            stream << event_indices1.at(i);
            stream << event_indices2.at(i);
            stream << event_indices3.at(i);
            std::string model_str = stream.str();
            std::cout << model_str << "\n";
            REQUIRE(model_counts.count(model_str) == 1);
            REQUIRE(nevent_counts.count(nevents.at(i)) == 1);
            ++model_counts[model_str];
            ++nevent_counts[nevents.at(i)];
            if (nevents.at(i) == 1) {
                REQUIRE(event_indices1.at(i) == event_indices2.at(i));
                REQUIRE(event_indices1.at(i) == event_indices3.at(i));
                REQUIRE(heights1.at(i) == heights2.at(i));
                REQUIRE(heights1.at(i) == heights3.at(i));
            }
            else if (nevents.at(i) == 3) {
                REQUIRE(event_indices1.at(i) != event_indices2.at(i));
                REQUIRE(event_indices1.at(i) != event_indices3.at(i));
                REQUIRE(event_indices2.at(i) != event_indices3.at(i));
                REQUIRE(heights1.at(i) != heights2.at(i));
                REQUIRE(heights1.at(i) != heights3.at(i));
                REQUIRE(heights2.at(i) != heights3.at(i));
            }
        }
        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == 10001);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == 10001);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));
        REQUIRE((model_counts.at("000") / 10001.0) == Approx(1.0/3.0).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "ln_prior",
                "ln_prior_kya",
                "ln_prior_pop1",
                "number_of_events",
                "root_height_kya",
                "root_height_pop1",
                "root_height_pop1b",
                "root_height_index_kya",
                "root_height_index_pop1",
                "root_height_index_pop1b"
        };
        SampleSummarizer<double> summary;
        for (auto const &kv: prior_sample.data) {
            bool test = true;
            for (auto const &p: columns_to_ignore) {
                if (kv.first == p) {
                    test = false;
                }
            }
            if (test) {
                summary = prior_sample.summarize<double>(kv.first);
                REQUIRE(summary.variance() == 0.0);
            }
        }

        delete[] cfg_path;
    }
}
