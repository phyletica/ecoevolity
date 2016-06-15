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
        char arg2[] = "111";
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
        REQUIRE(summary.mean() == Approx(shape * scale).epsilon(0.001));
        REQUIRE(summary.variance() == Approx(shape * scale * scale).epsilon(0.001));

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
        char arg2[] = "111";
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
        REQUIRE(summary.mean() == Approx(shape * scale).epsilon(0.001));
        REQUIRE(summary.variance() == Approx(shape * scale * scale).epsilon(0.001));

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
        char arg2[] = "111";
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

        for (auto v: samples) {
            std::cout << v << "\n";
        }

        SampleSummarizer<double> summary = prior_sample.summarize<double>("pop_size_root_kya");
        REQUIRE(summary.mean() == Approx(shape * scale).epsilon(0.001));
        REQUIRE(summary.variance() == Approx(shape * scale * scale).epsilon(0.001));

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
        char arg2[] = "111";
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

        for (auto v: samples) {
            std::cout << v << "\n";
        }

        SampleSummarizer<double> summary = prior_sample.summarize<double>("pop_size_root_kya");
        REQUIRE(summary.mean() == Approx(shape * scale).epsilon(0.001));
        REQUIRE(summary.variance() == Approx(shape * scale * scale).epsilon(0.001));

        delete[] cfg_path;
    }
}
