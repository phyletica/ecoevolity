#include "catch.hpp"
#include "ecoevolity/ecoevolity.hpp"

#include "ecoevolity/rng.hpp"
#include "ecoevolity/path.hpp"
#include "ecoevolity/spreadsheet.hpp"

RandomNumberGenerator _PRIOR_SAMPLING_RNG = RandomNumberGenerator();

TEST_CASE("Testing sampling from prior with CompositeHeightSizeRateMixer with 6 pairs", "[SamplingPrior]") {

    SECTION("Testing 6 pairs with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        double size1_shape = 10.0;
        double size1_scale = 0.1;
        double size2_shape = 2.0;
        double size2_scale = 0.2;
        double size3_shape = 5.0;
        double size3_scale = 0.2;
        double size4_shape = 4.0;
        double size4_scale = 0.5;
        double size5_shape = 8.0;
        double size5_scale = 0.2;
        double size6_shape = 6.0;
        double size6_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-collection-scaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-collection-scaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 2000000\n";
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size4_shape << "\n";
        os << "                    scale: " << size4_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname4.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size5_shape << "\n";
        os << "                    scale: " << size5_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname5.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size6_shape << "\n";
        os << "                    scale: " << size6_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "123456";
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

        unsigned int expected_sample_size = 200001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");
        SampleSummarizer<double> height_summary4 = prior_sample.summarize<double>("root_height_pop1c");
        SampleSummarizer<double> height_summary5 = prior_sample.summarize<double>("root_height_pop1d");
        SampleSummarizer<double> height_summary6 = prior_sample.summarize<double>("root_height_pop1e");

        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary4.sample_size() == expected_sample_size);
        REQUIRE(height_summary4.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary4.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary5.sample_size() == expected_sample_size);
        REQUIRE(height_summary5.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary5.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary6.sample_size() == expected_sample_size);
        REQUIRE(height_summary6.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary6.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        SampleSummarizer<double> size1_summary_a = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size1_summary_b = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size1_summary_c = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> size2_summary_a = prior_sample.summarize<double>("pop_size_pop1");
        SampleSummarizer<double> size2_summary_b = prior_sample.summarize<double>("pop_size_pop2");
        SampleSummarizer<double> size2_summary_c = prior_sample.summarize<double>("pop_size_root_pop1");
        SampleSummarizer<double> size3_summary_a = prior_sample.summarize<double>("pop_size_pop1b");
        SampleSummarizer<double> size3_summary_b = prior_sample.summarize<double>("pop_size_pop2b");
        SampleSummarizer<double> size3_summary_c = prior_sample.summarize<double>("pop_size_root_pop1b");
        SampleSummarizer<double> size4_summary_a = prior_sample.summarize<double>("pop_size_pop1c");
        SampleSummarizer<double> size4_summary_b = prior_sample.summarize<double>("pop_size_pop2c");
        SampleSummarizer<double> size4_summary_c = prior_sample.summarize<double>("pop_size_root_pop1c");
        SampleSummarizer<double> size5_summary_a = prior_sample.summarize<double>("pop_size_pop1d");
        SampleSummarizer<double> size5_summary_b = prior_sample.summarize<double>("pop_size_pop2d");
        SampleSummarizer<double> size5_summary_c = prior_sample.summarize<double>("pop_size_root_pop1d");
        SampleSummarizer<double> size6_summary_a = prior_sample.summarize<double>("pop_size_pop1e");
        SampleSummarizer<double> size6_summary_b = prior_sample.summarize<double>("pop_size_pop2e");
        SampleSummarizer<double> size6_summary_c = prior_sample.summarize<double>("pop_size_root_pop1e");

        REQUIRE(size1_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));

        REQUIRE(size2_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));

        REQUIRE(size3_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));

        REQUIRE(size4_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size4_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size4_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size4_summary_a.mean() == Approx(size4_shape * size4_scale).epsilon(0.01));
        REQUIRE(size4_summary_a.variance() == Approx(size4_shape * size4_scale * size4_scale).epsilon(0.01));
        REQUIRE(size4_summary_b.mean() == Approx(size4_shape * size4_scale).epsilon(0.01));
        REQUIRE(size4_summary_b.variance() == Approx(size4_shape * size4_scale * size4_scale).epsilon(0.01));
        REQUIRE(size4_summary_c.mean() == Approx(size4_shape * size4_scale).epsilon(0.01));
        REQUIRE(size4_summary_c.variance() == Approx(size4_shape * size4_scale * size4_scale).epsilon(0.01));

        REQUIRE(size5_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size5_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size5_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size5_summary_a.mean() == Approx(size5_shape * size5_scale).epsilon(0.01));
        REQUIRE(size5_summary_a.variance() == Approx(size5_shape * size5_scale * size5_scale).epsilon(0.01));
        REQUIRE(size5_summary_b.mean() == Approx(size5_shape * size5_scale).epsilon(0.01));
        REQUIRE(size5_summary_b.variance() == Approx(size5_shape * size5_scale * size5_scale).epsilon(0.01));
        REQUIRE(size5_summary_c.mean() == Approx(size5_shape * size5_scale).epsilon(0.01));
        REQUIRE(size5_summary_c.variance() == Approx(size5_shape * size5_scale * size5_scale).epsilon(0.01));

        REQUIRE(size6_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size6_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size6_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size6_summary_a.mean() == Approx(size6_shape * size6_scale).epsilon(0.01));
        REQUIRE(size6_summary_a.variance() == Approx(size6_shape * size6_scale * size6_scale).epsilon(0.01));
        REQUIRE(size6_summary_b.mean() == Approx(size6_shape * size6_scale).epsilon(0.01));
        REQUIRE(size6_summary_b.variance() == Approx(size6_shape * size6_scale * size6_scale).epsilon(0.01));
        REQUIRE(size6_summary_c.mean() == Approx(size6_shape * size6_scale).epsilon(0.01));
        REQUIRE(size6_summary_c.variance() == Approx(size6_shape * size6_scale * size6_scale).epsilon(0.01));

        delete[] cfg_path;
    }
}

TEST_CASE("Testing sampling from prior with CompositeHeightSizeRateMixer", "[SamplingPrior]") {

    SECTION("Testing gamma(10.0, 0.1) and gamma(5.0, 0.5) prior and no optimizing") {
        double time_shape = 10.0;
        double time_scale = 0.1;
        double size_shape = 5.0;
        double size_scale = 0.5;
        std::string auto_optimize = "false";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-collection-scaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-collection-scaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << time_shape << "\n";
        os << "        scale: " << time_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 1000000\n";
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shape << "\n";
        os << "                    scale: " << size_scale << "\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "654154";
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
        REQUIRE(samples.size() == 100001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("root_height_kya");
        REQUIRE(summary.mean() == Approx(time_shape * time_scale).epsilon(0.01));
        REQUIRE(summary.variance() == Approx(time_shape * time_scale * time_scale).epsilon(0.01));

        SampleSummarizer<double> size_summary1 = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size_summary2 = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size_summary3 = prior_sample.summarize<double>("pop_size_root_kya");

        REQUIRE(size_summary1.sample_size() == 100001);
        REQUIRE(size_summary2.sample_size() == 100001);
        REQUIRE(size_summary3.sample_size() == 100001);

        REQUIRE(size_summary1.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary1.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));

        SampleSummarizer<double> f_summary = prior_sample.summarize<double>("freq_1_kya");
        REQUIRE(f_summary.mean() == Approx(0.5));
        REQUIRE(f_summary.variance() == Approx(0.0));


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
                "pop_size_kya",
                "pop_size_fas",
                "pop_size_root_kya",
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

    SECTION("Testing gamma(10.0, 0.1) and gamma(5.0, 0.5) prior with optimizing") {
        double time_shape = 10.0;
        double time_scale = 0.1;
        double size_shape = 5.0;
        double size_scale = 0.5;
        std::string auto_optimize = "true";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-collection-scaler-test2-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-collection-scaler-test2-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << time_shape << "\n";
        os << "        scale: " << time_scale << "\n";
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shape << "\n";
        os << "                    scale: " << size_scale << "\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "123456";
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
        REQUIRE(samples.size() == 50001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("root_height_kya");
        REQUIRE(summary.mean() == Approx(time_shape * time_scale).epsilon(0.01));
        REQUIRE(summary.variance() == Approx(time_shape * time_scale * time_scale).epsilon(0.01));

        SampleSummarizer<double> size_summary1 = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size_summary2 = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size_summary3 = prior_sample.summarize<double>("pop_size_root_kya");

        REQUIRE(size_summary1.sample_size() == 50001);
        REQUIRE(size_summary2.sample_size() == 50001);
        REQUIRE(size_summary3.sample_size() == 50001);

        REQUIRE(size_summary1.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary1.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));

        SampleSummarizer<double> f_summary = prior_sample.summarize<double>("freq_1_kya");
        REQUIRE(f_summary.mean() == Approx(0.5));
        REQUIRE(f_summary.variance() == Approx(0.0));


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
                "pop_size_kya",
                "pop_size_fas",
                "pop_size_root_kya",
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

    SECTION("Testing gamma(10.0, 0.1) and gamma(5.0, 0.5) prior for singleton with optimizing") {
        double time_shape = 10.0;
        double time_scale = 0.1;
        double size_shape = 5.0;
        double size_scale = 0.5;
        std::string auto_optimize = "true";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-collection-scaler-test2-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-collection-scaler-test2-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << time_shape << "\n";
        os << "        scale: " << time_scale << "\n";
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shape << "\n";
        os << "                    scale: " << size_scale << "\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129-singleton.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "934028";
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
        REQUIRE(samples.size() == 50001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("root_height_kya");
        REQUIRE(summary.mean() == Approx(time_shape * time_scale).epsilon(0.01));
        REQUIRE(summary.variance() == Approx(time_shape * time_scale * time_scale).epsilon(0.01));

        SampleSummarizer<double> size_summary1 = prior_sample.summarize<double>("pop_size_kya");
        REQUIRE_THROWS(prior_sample.summarize<double>("pop_size_fas"));
        SampleSummarizer<double> size_summary3 = prior_sample.summarize<double>("pop_size_root_kya");

        REQUIRE(size_summary1.sample_size() == 50001);
        REQUIRE(size_summary3.sample_size() == 50001);

        REQUIRE(size_summary1.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary1.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));

        SampleSummarizer<double> f_summary = prior_sample.summarize<double>("freq_1_kya");
        REQUIRE(f_summary.mean() == Approx(0.5));
        REQUIRE(f_summary.variance() == Approx(0.0));


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
                "pop_size_kya",
                "pop_size_root_kya",
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

    SECTION("Testing gamma(10.0, 0.1) and gamma(5.0, 0.5) prior with constrained pop sizes and optimizing") {
        double time_shape = 10.0;
        double time_scale = 0.1;
        double size_shape = 5.0;
        double size_scale = 0.5;
        std::string auto_optimize = "true";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-collection-scaler-test3-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-collection-scaler-test3-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << time_shape << "\n";
        os << "        scale: " << time_scale << "\n";
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shape << "\n";
        os << "                    scale: " << size_scale << "\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "123456";
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
        REQUIRE(samples.size() == 50001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("root_height_kya");
        REQUIRE(summary.mean() == Approx(time_shape * time_scale).epsilon(0.01));
        REQUIRE(summary.variance() == Approx(time_shape * time_scale * time_scale).epsilon(0.01));

        SampleSummarizer<double> size_summary1 = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size_summary2 = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size_summary3 = prior_sample.summarize<double>("pop_size_root_kya");

        std::vector<double> samples1 = prior_sample.get<double>("pop_size_kya");
        std::vector<double> samples2 = prior_sample.get<double>("pop_size_fas");
        std::vector<double> samples3 = prior_sample.get<double>("pop_size_root_kya");

        for (size_t i = 0; i < samples1.size(); ++i) {
            REQUIRE(samples1.at(i) == samples2.at(i));
            REQUIRE(samples1.at(i) == samples3.at(i));
        }

        REQUIRE(size_summary1.sample_size() == 50001);
        REQUIRE(size_summary2.sample_size() == 50001);
        REQUIRE(size_summary3.sample_size() == 50001);

        REQUIRE(size_summary1.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary1.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));

        SampleSummarizer<double> f_summary = prior_sample.summarize<double>("freq_1_kya");
        REQUIRE(f_summary.mean() == Approx(0.5));
        REQUIRE(f_summary.variance() == Approx(0.0));


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
                "pop_size_kya",
                "pop_size_fas",
                "pop_size_root_kya",
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

    SECTION("Testing gamma(10.0, 0.1) and gamma(5.0, 0.5) prior with fixed pop sizes and optimizing") {
        double time_shape = 10.0;
        double time_scale = 0.1;
        double size_shape = 5.0;
        double size_scale = 0.5;
        std::string auto_optimize = "true";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-collection-scaler-test3-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-collection-scaler-test3-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << time_shape << "\n";
        os << "        scale: " << time_scale << "\n";
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "123456";
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
        // Should be no operators with any weight; should throw error
        REQUIRE_THROWS(ecoevolity_main(argc, argv));

        delete[] cfg_path;
    }

    SECTION("Testing gamma(10.0, 0.1) and gamma(5.0, 0.5) prior with fixed pop sizes, optimizing, 2 pairs") {
        double time_shape = 10.0;
        double time_scale = 0.1;
        double size_shape = 5.0;
        double size_scale = 0.5;
        std::string auto_optimize = "true";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-collection-scaler-test3-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-collection-scaler-test3-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 1]\n";
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << time_shape << "\n";
        os << "        scale: " << time_scale << "\n";
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shape << "\n";
        os << "                    scale: " << size_scale << "\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "123456";
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
        REQUIRE(samples.size() == 50001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("root_height_kya");
        REQUIRE(summary.mean() == Approx(time_shape * time_scale).epsilon(0.01));
        REQUIRE(summary.variance() == Approx(time_shape * time_scale * time_scale).epsilon(0.01));

        SampleSummarizer<double> summary2 = prior_sample.summarize<double>("root_height_pop1");
        REQUIRE(summary2.mean() == Approx(time_shape * time_scale).epsilon(0.01));
        REQUIRE(summary2.variance() == Approx(time_shape * time_scale * time_scale).epsilon(0.01));

        SampleSummarizer<double> size_summary1 = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size_summary2 = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size_summary3 = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> size_summary21 = prior_sample.summarize<double>("pop_size_pop1");
        SampleSummarizer<double> size_summary22 = prior_sample.summarize<double>("pop_size_pop2");
        SampleSummarizer<double> size_summary23 = prior_sample.summarize<double>("pop_size_root_pop1");

        std::vector<double> samples1 = prior_sample.get<double>("pop_size_kya");
        std::vector<double> samples2 = prior_sample.get<double>("pop_size_fas");
        std::vector<double> samples3 = prior_sample.get<double>("pop_size_root_kya");

        std::vector<double> samples21 = prior_sample.get<double>("pop_size_pop1");
        std::vector<double> samples22 = prior_sample.get<double>("pop_size_pop2");
        std::vector<double> samples23 = prior_sample.get<double>("pop_size_root_pop1");

        for (size_t i = 0; i < samples1.size(); ++i) {
            REQUIRE(samples1.at(i) == samples2.at(i));
            REQUIRE(samples1.at(i) == samples3.at(i));
            if (i > 10) {
                REQUIRE(samples21.at(i) != samples22.at(i));
                REQUIRE(samples21.at(i) != samples23.at(i));
                REQUIRE(samples22.at(i) != samples23.at(i));
            }
        }

        REQUIRE(size_summary1.sample_size() == 50001);
        REQUIRE(size_summary2.sample_size() == 50001);
        REQUIRE(size_summary3.sample_size() == 50001);
        REQUIRE(size_summary21.sample_size() == 50001);
        REQUIRE(size_summary22.sample_size() == 50001);
        REQUIRE(size_summary23.sample_size() == 50001);

        REQUIRE(size_summary1.mean() == Approx(0.005));
        REQUIRE(size_summary1.variance() == Approx(0.0));
        REQUIRE(size_summary2.mean() == Approx(0.005));
        REQUIRE(size_summary2.variance() == Approx(0.0));
        REQUIRE(size_summary3.mean() == Approx(0.005));
        REQUIRE(size_summary3.variance() == Approx(0.0));

        REQUIRE(size_summary21.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary21.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary22.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary22.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary23.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary23.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));

        SampleSummarizer<double> f_summary = prior_sample.summarize<double>("freq_1_kya");
        REQUIRE(f_summary.mean() == Approx(0.5));
        REQUIRE(f_summary.variance() == Approx(0.0));


        // Make sure the rest of the prior sample is as expected
        summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(summary.mean() == 0.0);
        REQUIRE(summary.variance() == 0.0);
        summary = prior_sample.summarize<double>("ln_likelihood_kya");
        REQUIRE(summary.mean() == 0.0);
        REQUIRE(summary.variance() == 0.0);
        summary = prior_sample.summarize<double>("ln_likelihood_pop1");
        REQUIRE(summary.mean() == 0.0);
        REQUIRE(summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "root_height_kya",
                "root_height_pop1",
                "pop_size_kya",
                "pop_size_pop1",
                "pop_size_fas",
                "pop_size_pop2",
                "pop_size_root_kya",
                "pop_size_root_pop1",
                "ln_prior",
                "ln_prior_pop1",
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

    SECTION("Testing gamma(10.0, 10.0) and gamma(2.0, 0.5) prior with optimizing") {
        double time_shape = 10.0;
        double time_scale = 10.0;
        double size_shape = 2.0;
        double size_scale = 0.5;
        std::string auto_optimize = "true";
        std::string tag = _PRIOR_SAMPLING_RNG.random_string(10);
        std::string test_path = "data/tmp-config-collection-scaler-test4-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-collection-scaler-test4-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << time_shape << "\n";
        os << "        scale: " << time_scale << "\n";
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shape << "\n";
        os << "                    scale: " << size_scale << "\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "123456";
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
        REQUIRE(samples.size() == 50001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("root_height_kya");
        REQUIRE(summary.mean() == Approx(time_shape * time_scale).epsilon(0.01));
        REQUIRE(summary.variance() == Approx(time_shape * time_scale * time_scale).epsilon(0.01));

        SampleSummarizer<double> size_summary1 = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size_summary2 = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size_summary3 = prior_sample.summarize<double>("pop_size_root_kya");

        REQUIRE(size_summary1.sample_size() == 50001);
        REQUIRE(size_summary2.sample_size() == 50001);
        REQUIRE(size_summary3.sample_size() == 50001);

        REQUIRE(size_summary1.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary1.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));

        SampleSummarizer<double> f_summary = prior_sample.summarize<double>("freq_1_kya");
        REQUIRE(f_summary.mean() == Approx(0.5));
        REQUIRE(f_summary.variance() == Approx(0.0));


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
                "pop_size_kya",
                "pop_size_fas",
                "pop_size_root_kya",
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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

        SampleSummarizer<double> f_summary = prior_sample.summarize<double>("freq_1_kya");
        REQUIRE(f_summary.mean() == Approx(0.5));
        REQUIRE(f_summary.variance() == Approx(0.0));


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
        std::string test_path = "data/tmp-config-comp-height-scaler-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-comp-height-scaler-" + tag + "-state-run-1.log";
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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

        SampleSummarizer<double> f_summary = prior_sample.summarize<double>("freq_1_kya");
        REQUIRE(f_summary.mean() == Approx(0.5));
        REQUIRE(f_summary.variance() == Approx(0.0));

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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << shape << "\n";
        os << "                    scale: " << scale << "\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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

        SampleSummarizer<double> f_summary = prior_sample.summarize<double>("freq_1_kya");
        REQUIRE(f_summary.mean() == Approx(0.5));
        REQUIRE(f_summary.variance() == Approx(0.0));

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
        std::string test_path = "data/tmp-config-root-pop-size-scaler-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-root-pop-size-scaler-" + tag + "-state-run-1.log";
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << shape << "\n";
        os << "                    scale: " << scale << "\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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

        SampleSummarizer<double> f_summary = prior_sample.summarize<double>("freq_1_kya");
        REQUIRE(f_summary.mean() == Approx(0.5));
        REQUIRE(f_summary.variance() == Approx(0.0));

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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << shape << "\n";
        os << "                    scale: " << scale << "\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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

        SampleSummarizer<double> f_summary = prior_sample.summarize<double>("freq_1_kya");
        REQUIRE(f_summary.mean() == Approx(0.5));
        REQUIRE(f_summary.variance() == Approx(0.0));

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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << shape << "\n";
        os << "                    scale: " << scale << "\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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

        SampleSummarizer<double> f_summary = prior_sample.summarize<double>("freq_1_kya");
        REQUIRE(f_summary.mean() == Approx(0.5));
        REQUIRE(f_summary.variance() == Approx(0.0));

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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << shape << "\n";
        os << "                    scale: " << scale << "\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << shape << "\n";
        os << "                    scale: " << scale << "\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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

TEST_CASE("Testing sampling from beta(1.5, 2.5) prior with FreqMover",
        "[SamplingPrior]") {

    SECTION("Testing beta(1.5, 2.5) prior and no optimizing") {
        double a = 1.5;
        double b = 2.5;
        double expected_mean = a / (a + b);
        double expected_variance = (a * b) / ((a + b) * (a + b) * (a + b + 1.0));

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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: false\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << a << "\n";
        os << "                    beta: " << b << "\n";
        os << "        mutation_rate:\n";
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
        std::vector<double> samples = prior_sample.get<double>("freq_1_kya");
        REQUIRE(samples.size() == 10001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("freq_1_kya");
        REQUIRE(summary.mean() == Approx(expected_mean).epsilon(0.01));
        REQUIRE(summary.variance() == Approx(expected_variance).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(summary.mean() == 0.0);
        REQUIRE(summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "ln_prior",
                "ln_prior_kya",
                "freq_1_kya"
                };
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

TEST_CASE("Testing sampling from beta(1.5, 2.5) prior with FreqMover and optimizing",
        "[SamplingPrior]") {

    SECTION("Testing beta(1.5, 2.5) prior and optimizing") {
        double a = 1.5;
        double b = 2.5;
        double expected_mean = a / (a + b);
        double expected_variance = (a * b) / ((a + b) * (a + b) * (a + b + 1.0));

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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: false\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << a << "\n";
        os << "                    beta: " << b << "\n";
        os << "        mutation_rate:\n";
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
        std::vector<double> samples = prior_sample.get<double>("freq_1_kya");
        REQUIRE(samples.size() == 10001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("freq_1_kya");
        REQUIRE(summary.mean() == Approx(expected_mean).epsilon(0.001));
        REQUIRE(summary.variance() == Approx(expected_variance).epsilon(0.001));

        // Make sure the rest of the prior sample is as expected
        summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(summary.mean() == 0.0);
        REQUIRE(summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "ln_prior",
                "ln_prior_kya",
                "freq_1_kya"
                };
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

TEST_CASE("Testing sampling from beta(2.5, 1.5) prior with FreqMover and optimizing",
        "[SamplingPrior]") {

    SECTION("Testing beta(2.5, 1.5) prior and optimizing") {
        double a = 2.5;
        double b = 1.5;
        double expected_mean = a / (a + b);
        double expected_variance = (a * b) / ((a + b) * (a + b) * (a + b + 1.0));

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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: false\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << a << "\n";
        os << "                    beta: " << b << "\n";
        os << "        mutation_rate:\n";
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
        std::vector<double> samples = prior_sample.get<double>("freq_1_kya");
        REQUIRE(samples.size() == 10001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("freq_1_kya");
        REQUIRE(summary.mean() == Approx(expected_mean).epsilon(0.001));
        REQUIRE(summary.variance() == Approx(expected_variance).epsilon(0.001));

        // Make sure the rest of the prior sample is as expected
        summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(summary.mean() == 0.0);
        REQUIRE(summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "ln_prior",
                "ln_prior_kya",
                "freq_1_kya"
                };
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

TEST_CASE("Testing sampling from prior with ComparisonMutationRateScaler",
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.001\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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
        std::vector<double> samples = prior_sample.get<double>("mutation_rate_kya");
        REQUIRE(samples.size() == 10001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("mutation_rate_kya");
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
                "mutation_rate_kya"};
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

TEST_CASE("Testing sampling from prior with ComparisonMutationRateScaler with optimizing",
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.001\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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
        std::vector<double> samples = prior_sample.get<double>("mutation_rate_kya");
        REQUIRE(samples.size() == 10001);

        SampleSummarizer<double> summary = prior_sample.summarize<double>("mutation_rate_kya");
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
                "mutation_rate_kya"};
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
        double f_a = 1.0;
        double f_b = 0.5;
        double expected_f_mean = f_a / (f_a + f_b);
        double expected_f_variance = (f_a * f_b) / ((f_a + f_b) * (f_a + f_b) * (f_a + f_b + 1.0));
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shape << "\n";
        os << "                    scale: " << size_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f_a << "\n";
        os << "                    beta: " << f_b << "\n";
        os << "        mutation_rate:\n";
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
        SampleSummarizer<double> f_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> mult_summary = prior_sample.summarize<double>("mutation_rate_kya");

        REQUIRE(height_summary.sample_size() == 50001);
        REQUIRE(size_summary1.sample_size() == 50001);
        REQUIRE(size_summary2.sample_size() == 50001);
        REQUIRE(size_summary3.sample_size() == 50001);
        REQUIRE(f_summary.sample_size() == 50001);
        REQUIRE(mult_summary.sample_size() == 50001);

        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(size_summary1.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary1.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(f_summary.mean() == Approx(expected_f_mean).epsilon(0.005));
        REQUIRE(f_summary.variance() == Approx(expected_f_variance).epsilon(0.005));
        REQUIRE(mult_summary.mean() == Approx(mult_shape * mult_scale).epsilon(0.01));
        REQUIRE(mult_summary.variance() == Approx(mult_shape * mult_scale * mult_scale).epsilon(0.01));

        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.sample_size() == 50001);
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

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
        double f_a = 1.0;
        double f_b = 0.5;
        double expected_f_mean = f_a / (f_a + f_b);
        double expected_f_variance = (f_a * f_b) / ((f_a + f_b) * (f_a + f_b) * (f_a + f_b + 1.0));
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shape << "\n";
        os << "                    scale: " << size_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f_a << "\n";
        os << "                    beta: " << f_b << "\n";
        os << "        mutation_rate:\n";
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
        SampleSummarizer<double> f_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> mult_summary = prior_sample.summarize<double>("mutation_rate_kya");

        REQUIRE(height_summary.sample_size() == 50001);
        REQUIRE(size_summary1.sample_size() == 50001);
        REQUIRE(size_summary2.sample_size() == 50001);
        REQUIRE(size_summary3.sample_size() == 50001);
        REQUIRE(f_summary.sample_size() == 50001);
        REQUIRE(mult_summary.sample_size() == 50001);

        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(size_summary1.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary1.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(f_summary.mean() == Approx(expected_f_mean).epsilon(0.005));
        REQUIRE(f_summary.variance() == Approx(expected_f_variance).epsilon(0.005));
        REQUIRE(mult_summary.mean() == Approx(mult_shape * mult_scale).epsilon(0.01));
        REQUIRE(mult_summary.variance() == Approx(mult_shape * mult_scale * mult_scale).epsilon(0.01));

        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.sample_size() == 50001);
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

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

TEST_CASE("Testing fully parameterized model for one pair with optimization and CompositeHeightSizeRateMixer",
        "[SamplingPrior]") {

    SECTION("Testing with optimizing") {
        double height_shape = 1.0;
        double height_scale = 0.1;
        double size_shape = 10.0;
        double size_scale = 0.0001;
        double f_a = 1.0;
        double f_b = 0.5;
        double expected_f_mean = f_a / (f_a + f_b);
        double expected_f_variance = (f_a * f_b) / ((f_a + f_b) * (f_a + f_b) * (f_a + f_b + 1.0));
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shape << "\n";
        os << "                    scale: " << size_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f_a << "\n";
        os << "                    beta: " << f_b << "\n";
        os << "        mutation_rate:\n";
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
        char arg2[] = "9876987345";
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
        SampleSummarizer<double> f_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> mult_summary = prior_sample.summarize<double>("mutation_rate_kya");

        REQUIRE(height_summary.sample_size() == 50001);
        REQUIRE(size_summary1.sample_size() == 50001);
        REQUIRE(size_summary2.sample_size() == 50001);
        REQUIRE(size_summary3.sample_size() == 50001);
        REQUIRE(f_summary.sample_size() == 50001);
        REQUIRE(mult_summary.sample_size() == 50001);

        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(size_summary1.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary1.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary2.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.mean() == Approx(size_shape * size_scale).epsilon(0.01));
        REQUIRE(size_summary3.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.01));
        REQUIRE(f_summary.mean() == Approx(expected_f_mean).epsilon(0.005));
        REQUIRE(f_summary.variance() == Approx(expected_f_variance).epsilon(0.005));
        REQUIRE(mult_summary.mean() == Approx(mult_shape * mult_scale).epsilon(0.01));
        REQUIRE(mult_summary.variance() == Approx(mult_shape * mult_scale * mult_scale).epsilon(0.01));

        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.sample_size() == 50001);
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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
        double concentration = 1.0;
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
        os << "                value: " << concentration << "\n";
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / 10001.0) == Approx(std::exp(
                    get_dpp_log_prior_probability(m, concentration))).epsilon(0.01));
        }

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

TEST_CASE("Testing DPP with 3 pairs and alpha 4.0", "[SamplingPrior]") {

    SECTION("Testing alpha fixed to 4.0") {
        double height_shape = 10.0;
        double height_scale = 0.1;
        double concentration = 4.0;
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
        os << "                value: " << concentration << "\n";
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / 10001.0) == Approx(std::exp(
                    get_dpp_log_prior_probability(m, concentration))).epsilon(0.01));
        }

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

TEST_CASE("Testing DPP with 6 pairs and alpha 1.7", "[SamplingPrior]") {

    SECTION("Testing alpha fixed to 1.7") {
        double height_shape = 10.0;
        double height_scale = 0.1;
        double concentration = 1.7;
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
        os << "                value: " << concentration << "\n";
        os << "                estimate: false\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 1000000\n";
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
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname4.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname5.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "64551529";
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
        SampleSummarizer<double> height_summary4 = prior_sample.summarize<double>("root_height_pop1c");
        SampleSummarizer<double> height_summary5 = prior_sample.summarize<double>("root_height_pop1d");
        SampleSummarizer<double> height_summary6 = prior_sample.summarize<double>("root_height_pop1e");
        REQUIRE(height_summary1.sample_size() == 100001);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary2.sample_size() == 100001);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary3.sample_size() == 100001);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.005));
        REQUIRE(height_summary4.sample_size() == 100001);
        REQUIRE(height_summary4.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary4.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary5.sample_size() == 100001);
        REQUIRE(height_summary5.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary5.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary6.sample_size() == 100001);
        REQUIRE(height_summary6.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary6.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));

        std::vector<int> nevents = prior_sample.get<int>("number_of_events");
        std::vector<int> event_indices1 = prior_sample.get<int>("root_height_index_kya");
        std::vector<int> event_indices2 = prior_sample.get<int>("root_height_index_pop1");
        std::vector<int> event_indices3 = prior_sample.get<int>("root_height_index_pop1b");
        std::vector<int> event_indices4 = prior_sample.get<int>("root_height_index_pop1c");
        std::vector<int> event_indices5 = prior_sample.get<int>("root_height_index_pop1d");
        std::vector<int> event_indices6 = prior_sample.get<int>("root_height_index_pop1e");
        std::vector<double> heights1 = prior_sample.get<double>("root_height_kya");
        std::vector<double> heights2 = prior_sample.get<double>("root_height_pop1");
        std::vector<double> heights3 = prior_sample.get<double>("root_height_pop1b");
        std::vector<double> heights4 = prior_sample.get<double>("root_height_pop1c");
        std::vector<double> heights5 = prior_sample.get<double>("root_height_pop1d");
        std::vector<double> heights6 = prior_sample.get<double>("root_height_pop1e");

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
                {4, 0},
                {5, 0},
                {6, 0}
        };
        for (size_t i = 0; i < nevents.size(); ++i) {
            std::ostringstream stream;
            stream << event_indices1.at(i);
            stream << event_indices2.at(i);
            stream << event_indices3.at(i);
            stream << event_indices4.at(i);
            stream << event_indices5.at(i);
            stream << event_indices6.at(i);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
            REQUIRE(nevent_counts.count(nevents.at(i)) == 1);
            ++nevent_counts[nevents.at(i)];
            if (nevents.at(i) == 1) {
                REQUIRE(event_indices1.at(i) == event_indices2.at(i));
                REQUIRE(event_indices1.at(i) == event_indices3.at(i));
                REQUIRE(event_indices1.at(i) == event_indices4.at(i));
                REQUIRE(event_indices1.at(i) == event_indices5.at(i));
                REQUIRE(event_indices1.at(i) == event_indices6.at(i));
                REQUIRE(heights1.at(i) == heights2.at(i));
                REQUIRE(heights1.at(i) == heights3.at(i));
                REQUIRE(heights1.at(i) == heights4.at(i));
                REQUIRE(heights1.at(i) == heights5.at(i));
                REQUIRE(heights1.at(i) == heights6.at(i));
            }
            else if (nevents.at(i) == 6) {
                REQUIRE(event_indices1.at(i) != event_indices2.at(i));
                REQUIRE(event_indices1.at(i) != event_indices3.at(i));
                REQUIRE(event_indices1.at(i) != event_indices4.at(i));
                REQUIRE(event_indices1.at(i) != event_indices5.at(i));
                REQUIRE(event_indices1.at(i) != event_indices6.at(i));
                REQUIRE(event_indices2.at(i) != event_indices3.at(i));
                REQUIRE(event_indices2.at(i) != event_indices4.at(i));
                REQUIRE(event_indices2.at(i) != event_indices5.at(i));
                REQUIRE(event_indices2.at(i) != event_indices6.at(i));
                REQUIRE(event_indices3.at(i) != event_indices4.at(i));
                REQUIRE(event_indices3.at(i) != event_indices5.at(i));
                REQUIRE(event_indices3.at(i) != event_indices6.at(i));
                REQUIRE(event_indices4.at(i) != event_indices5.at(i));
                REQUIRE(event_indices4.at(i) != event_indices6.at(i));
                REQUIRE(event_indices5.at(i) != event_indices6.at(i));
                REQUIRE(heights1.at(i) != heights2.at(i));
                REQUIRE(heights1.at(i) != heights3.at(i));
                REQUIRE(heights1.at(i) != heights4.at(i));
                REQUIRE(heights1.at(i) != heights5.at(i));
                REQUIRE(heights1.at(i) != heights6.at(i));
                REQUIRE(heights2.at(i) != heights3.at(i));
                REQUIRE(heights2.at(i) != heights4.at(i));
                REQUIRE(heights2.at(i) != heights5.at(i));
                REQUIRE(heights2.at(i) != heights6.at(i));
                REQUIRE(heights3.at(i) != heights4.at(i));
                REQUIRE(heights3.at(i) != heights5.at(i));
                REQUIRE(heights3.at(i) != heights6.at(i));
                REQUIRE(heights4.at(i) != heights5.at(i));
                REQUIRE(heights4.at(i) != heights6.at(i));
                REQUIRE(heights5.at(i) != heights6.at(i));
            }
        }
        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == 100001);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == 100001);

        REQUIRE(model_counts.at("000000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012345") == nevent_counts.at(6));

        for (std::string m: {"000000", "012345", "012344", "012340", "001234", "012314"}) {
            REQUIRE((model_counts.at(m) / 100001.0) == Approx(std::exp(
                    get_dpp_log_prior_probability(m, concentration))).epsilon(0.001));
        }

        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing DPP with 3 pairs and fully parameterized", "[SamplingPrior]") {

    SECTION("Testing alpha integrated") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_a = 2.0;
        double f1_b = 1.1;
        double f2_a = 1.0;
        double f2_b = 0.5;
        double f3_a = 1.5;
        double f3_b = 1.8;
        double expected_f1_mean = f1_a / (f1_a + f1_b);
        double expected_f1_variance = (f1_a * f1_b) / ((f1_a + f1_b) * (f1_a + f1_b) * (f1_a + f1_b + 1.0));
        double expected_f2_mean = f2_a / (f2_a + f2_b);
        double expected_f2_variance = (f2_a * f2_b) / ((f2_a + f2_b) * (f2_a + f2_b) * (f2_a + f2_b + 1.0));
        double expected_f3_mean = f3_a / (f3_a + f3_b);
        double expected_f3_variance = (f3_a * f3_b) / ((f3_a + f3_b) * (f3_a + f3_b) * (f3_a + f3_b + 1.0));
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
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
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 2000000\n";
        os << "    sample_frequency: 100\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_a << "\n";
        os << "                    beta: " << f1_b << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_a << "\n";
        os << "                    beta: " << f2_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_a << "\n";
        os << "                    beta: " << f3_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
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

        unsigned int expected_sample_size = 20001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");
        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        SampleSummarizer<double> size1_summary_a = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size1_summary_b = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size1_summary_c = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> size2_summary_a = prior_sample.summarize<double>("pop_size_pop1");
        SampleSummarizer<double> size2_summary_b = prior_sample.summarize<double>("pop_size_pop2");
        SampleSummarizer<double> size2_summary_c = prior_sample.summarize<double>("pop_size_root_pop1");
        SampleSummarizer<double> size3_summary_a = prior_sample.summarize<double>("pop_size_pop1b");
        SampleSummarizer<double> size3_summary_b = prior_sample.summarize<double>("pop_size_pop2b");
        SampleSummarizer<double> size3_summary_c = prior_sample.summarize<double>("pop_size_root_pop1b");

        REQUIRE(size1_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));

        REQUIRE(size2_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));

        REQUIRE(size3_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));

        SampleSummarizer<double> f1_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> f2_summary = prior_sample.summarize<double>("freq_1_pop1");
        SampleSummarizer<double> f3_summary = prior_sample.summarize<double>("freq_1_pop1b");
        REQUIRE(f1_summary.sample_size() == expected_sample_size);
        REQUIRE(f2_summary.sample_size() == expected_sample_size);
        REQUIRE(f3_summary.sample_size() == expected_sample_size);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.01));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.01));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.01));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.01));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.01));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.01));

        SampleSummarizer<double> mult1_summary = prior_sample.summarize<double>("mutation_rate_kya");
        SampleSummarizer<double> mult2_summary = prior_sample.summarize<double>("mutation_rate_pop1");
        SampleSummarizer<double> mult3_summary = prior_sample.summarize<double>("mutation_rate_pop1b");
        REQUIRE(mult1_summary.sample_size() == expected_sample_size);
        REQUIRE(mult2_summary.sample_size() == expected_sample_size);
        REQUIRE(mult3_summary.sample_size() == expected_sample_size);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.01));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.01));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.01));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.01));

        SampleSummarizer<double> conc_summary = prior_sample.summarize<double>("concentration");
        REQUIRE(conc_summary.sample_size() == expected_sample_size);
        REQUIRE(conc_summary.mean() == Approx(concentration_shape * concentration_scale).epsilon(0.01));
        REQUIRE(conc_summary.variance() == Approx(concentration_shape * concentration_scale * concentration_scale).epsilon(0.01));

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        REQUIRE((model_counts.at("000") / (double)expected_sample_size) == Approx(0.367).epsilon(0.01));
        REQUIRE((model_counts.at("012") / (double)expected_sample_size) == Approx(0.163).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing DPP with 3 pairs and fully parameterized and CompositeHeightSizeRateMixer",
        "[SamplingPrior]") {

    SECTION("Testing alpha integrated") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_a = 2.0;
        double f1_b = 1.1;
        double f2_a = 1.0;
        double f2_b = 0.5;
        double f3_a = 1.5;
        double f3_b = 1.8;
        double expected_f1_mean = f1_a / (f1_a + f1_b);
        double expected_f1_variance = (f1_a * f1_b) / ((f1_a + f1_b) * (f1_a + f1_b) * (f1_a + f1_b + 1.0));
        double expected_f2_mean = f2_a / (f2_a + f2_b);
        double expected_f2_variance = (f2_a * f2_b) / ((f2_a + f2_b) * (f2_a + f2_b) * (f2_a + f2_b + 1.0));
        double expected_f3_mean = f3_a / (f3_a + f3_b);
        double expected_f3_variance = (f3_a * f3_b) / ((f3_a + f3_b) * (f3_a + f3_b) * (f3_a + f3_b + 1.0));
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
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
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 2000000\n";
        os << "    sample_frequency: 100\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_a << "\n";
        os << "                    beta: " << f1_b << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_a << "\n";
        os << "                    beta: " << f2_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_a << "\n";
        os << "                    beta: " << f3_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
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

        unsigned int expected_sample_size = 20001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");
        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        SampleSummarizer<double> size1_summary_a = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size1_summary_b = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size1_summary_c = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> size2_summary_a = prior_sample.summarize<double>("pop_size_pop1");
        SampleSummarizer<double> size2_summary_b = prior_sample.summarize<double>("pop_size_pop2");
        SampleSummarizer<double> size2_summary_c = prior_sample.summarize<double>("pop_size_root_pop1");
        SampleSummarizer<double> size3_summary_a = prior_sample.summarize<double>("pop_size_pop1b");
        SampleSummarizer<double> size3_summary_b = prior_sample.summarize<double>("pop_size_pop2b");
        SampleSummarizer<double> size3_summary_c = prior_sample.summarize<double>("pop_size_root_pop1b");

        REQUIRE(size1_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));

        REQUIRE(size2_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));

        REQUIRE(size3_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));

        SampleSummarizer<double> f1_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> f2_summary = prior_sample.summarize<double>("freq_1_pop1");
        SampleSummarizer<double> f3_summary = prior_sample.summarize<double>("freq_1_pop1b");
        REQUIRE(f1_summary.sample_size() == expected_sample_size);
        REQUIRE(f2_summary.sample_size() == expected_sample_size);
        REQUIRE(f3_summary.sample_size() == expected_sample_size);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.01));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.01));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.01));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.01));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.01));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.01));

        SampleSummarizer<double> mult1_summary = prior_sample.summarize<double>("mutation_rate_kya");
        SampleSummarizer<double> mult2_summary = prior_sample.summarize<double>("mutation_rate_pop1");
        SampleSummarizer<double> mult3_summary = prior_sample.summarize<double>("mutation_rate_pop1b");
        REQUIRE(mult1_summary.sample_size() == expected_sample_size);
        REQUIRE(mult2_summary.sample_size() == expected_sample_size);
        REQUIRE(mult3_summary.sample_size() == expected_sample_size);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.01));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.01));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.01));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.01));

        SampleSummarizer<double> conc_summary = prior_sample.summarize<double>("concentration");
        REQUIRE(conc_summary.sample_size() == expected_sample_size);
        REQUIRE(conc_summary.mean() == Approx(concentration_shape * concentration_scale).epsilon(0.01));
        REQUIRE(conc_summary.variance() == Approx(concentration_shape * concentration_scale * concentration_scale).epsilon(0.01));

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        REQUIRE((model_counts.at("000") / (double)expected_sample_size) == Approx(0.367).epsilon(0.01));
        REQUIRE((model_counts.at("012") / (double)expected_sample_size) == Approx(0.163).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing sampling of small concentration", "[SamplingPrior]") {

    SECTION("Testing concentration sampling") {
        double height_shape = 10.0;
        double height_scale = 0.1;
        double concentration_shape = 10.0;
        double concentration_scale = 0.01;
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
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 500000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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

        unsigned int expected_sample_size = 50001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");
        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        SampleSummarizer<double> conc_summary = prior_sample.summarize<double>("concentration");
        REQUIRE(conc_summary.sample_size() == expected_sample_size);
        REQUIRE(conc_summary.mean() == Approx(concentration_shape * concentration_scale).epsilon(0.01));
        REQUIRE(conc_summary.variance() == Approx(concentration_shape * concentration_scale * concentration_scale).epsilon(0.01));

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        REQUIRE((model_counts.at("000") / (double)expected_sample_size) == Approx(0.866).epsilon(0.01));
        REQUIRE((model_counts.at("012") / (double)expected_sample_size) == Approx(0.0046).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "concentration",
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

TEST_CASE("Testing sampling of large concentration", "[SamplingPrior]") {

    SECTION("Testing concentration sampling") {
        double height_shape = 10.0;
        double height_scale = 0.1;
        double concentration_shape = 10.0;
        double concentration_scale = 2.0;
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
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 500000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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
        char arg2[] = "239472";
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

        unsigned int expected_sample_size = 50001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");
        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        SampleSummarizer<double> conc_summary = prior_sample.summarize<double>("concentration");
        REQUIRE(conc_summary.sample_size() == expected_sample_size);
        REQUIRE(conc_summary.mean() == Approx(concentration_shape * concentration_scale).epsilon(0.01));
        REQUIRE(conc_summary.variance() == Approx(concentration_shape * concentration_scale * concentration_scale).epsilon(0.1));

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        REQUIRE((model_counts.at("000") / (double)expected_sample_size) == Approx(0.0057).epsilon(0.01));
        REQUIRE((model_counts.at("012") / (double)expected_sample_size) == Approx(0.854).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "concentration",
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

TEST_CASE("Testing sampling of diffuse concentration", "[SamplingPrior]") {

    SECTION("Testing concentration sampling") {
        double height_shape = 10.0;
        double height_scale = 0.1;
        double concentration_shape = 0.5;
        double concentration_scale = 2.0;
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
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 1000000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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
        char arg2[] = "54866549";
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

        unsigned int expected_sample_size = 100001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");
        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.1));

        SampleSummarizer<double> conc_summary = prior_sample.summarize<double>("concentration");
        REQUIRE(conc_summary.sample_size() == expected_sample_size);
        REQUIRE(conc_summary.mean() == Approx(concentration_shape * concentration_scale).epsilon(0.01));
        REQUIRE(conc_summary.variance() == Approx(concentration_shape * concentration_scale * concentration_scale).epsilon(0.01));


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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        // Expected prior probs estimated by 10 million simulations using the
        // 'dmc_estimate_prior_probs.py' tool from the PyMsBayes package
        REQUIRE((model_counts.at("000") / (double)expected_sample_size) == Approx(0.554).epsilon(0.01));
        REQUIRE((model_counts.at("012") / (double)expected_sample_size) == Approx(0.140).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        std::vector<std::string> columns_to_ignore = {
                "generation",
                "concentration",
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

TEST_CASE("Testing DPP with 3 pairs, fully parameterized, and multithreading",
        "[SamplingPrior]") {

    SECTION("Testing multithreading") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_a = 2.0;
        double f1_b = 1.1;
        double f2_a = 1.0;
        double f2_b = 0.5;
        double f3_a = 1.5;
        double f3_b = 1.8;
        double expected_f1_mean = f1_a / (f1_a + f1_b);
        double expected_f1_variance = (f1_a * f1_b) / ((f1_a + f1_b) * (f1_a + f1_b) * (f1_a + f1_b + 1.0));
        double expected_f2_mean = f2_a / (f2_a + f2_b);
        double expected_f2_variance = (f2_a * f2_b) / ((f2_a + f2_b) * (f2_a + f2_b) * (f2_a + f2_b + 1.0));
        double expected_f3_mean = f3_a / (f3_a + f3_b);
        double expected_f3_variance = (f3_a * f3_b) / ((f3_a + f3_b) * (f3_a + f3_b) * (f3_a + f3_b + 1.0));
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
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
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 2000000\n";
        os << "    sample_frequency: 100\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_a << "\n";
        os << "                    beta: " << f1_b << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_a << "\n";
        os << "                    beta: " << f2_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_a << "\n";
        os << "                    beta: " << f3_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "--ignore-data";
        char arg4[] = "--nthreads";
        char arg5[] = "3";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
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

        unsigned int expected_sample_size = 20001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");
        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        SampleSummarizer<double> size1_summary_a = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size1_summary_b = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size1_summary_c = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> size2_summary_a = prior_sample.summarize<double>("pop_size_pop1");
        SampleSummarizer<double> size2_summary_b = prior_sample.summarize<double>("pop_size_pop2");
        SampleSummarizer<double> size2_summary_c = prior_sample.summarize<double>("pop_size_root_pop1");
        SampleSummarizer<double> size3_summary_a = prior_sample.summarize<double>("pop_size_pop1b");
        SampleSummarizer<double> size3_summary_b = prior_sample.summarize<double>("pop_size_pop2b");
        SampleSummarizer<double> size3_summary_c = prior_sample.summarize<double>("pop_size_root_pop1b");

        REQUIRE(size1_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));

        REQUIRE(size2_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));

        REQUIRE(size3_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));

        SampleSummarizer<double> f1_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> f2_summary = prior_sample.summarize<double>("freq_1_pop1");
        SampleSummarizer<double> f3_summary = prior_sample.summarize<double>("freq_1_pop1b");
        REQUIRE(f1_summary.sample_size() == expected_sample_size);
        REQUIRE(f2_summary.sample_size() == expected_sample_size);
        REQUIRE(f3_summary.sample_size() == expected_sample_size);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.01));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.01));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.01));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.01));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.01));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.01));

        SampleSummarizer<double> mult1_summary = prior_sample.summarize<double>("mutation_rate_kya");
        SampleSummarizer<double> mult2_summary = prior_sample.summarize<double>("mutation_rate_pop1");
        SampleSummarizer<double> mult3_summary = prior_sample.summarize<double>("mutation_rate_pop1b");
        REQUIRE(mult1_summary.sample_size() == expected_sample_size);
        REQUIRE(mult2_summary.sample_size() == expected_sample_size);
        REQUIRE(mult3_summary.sample_size() == expected_sample_size);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.01));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.01));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.01));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.01));

        SampleSummarizer<double> conc_summary = prior_sample.summarize<double>("concentration");
        REQUIRE(conc_summary.sample_size() == expected_sample_size);
        REQUIRE(conc_summary.mean() == Approx(concentration_shape * concentration_scale).epsilon(0.01));
        REQUIRE(conc_summary.variance() == Approx(concentration_shape * concentration_scale * concentration_scale).epsilon(0.01));

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        REQUIRE((model_counts.at("000") / (double)expected_sample_size) == Approx(0.367).epsilon(0.01));
        REQUIRE((model_counts.at("012") / (double)expected_sample_size) == Approx(0.163).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing DPP with 3 pairs, fully parameterized, and 2 threads",
        "[SamplingPrior]") {

    SECTION("Testing nthreads 2") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_a = 2.0;
        double f1_b = 1.1;
        double f2_a = 1.0;
        double f2_b = 0.5;
        double f3_a = 1.5;
        double f3_b = 1.8;
        double expected_f1_mean = f1_a / (f1_a + f1_b);
        double expected_f1_variance = (f1_a * f1_b) / ((f1_a + f1_b) * (f1_a + f1_b) * (f1_a + f1_b + 1.0));
        double expected_f2_mean = f2_a / (f2_a + f2_b);
        double expected_f2_variance = (f2_a * f2_b) / ((f2_a + f2_b) * (f2_a + f2_b) * (f2_a + f2_b + 1.0));
        double expected_f3_mean = f3_a / (f3_a + f3_b);
        double expected_f3_variance = (f3_a * f3_b) / ((f3_a + f3_b) * (f3_a + f3_b) * (f3_a + f3_b + 1.0));
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
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
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 2000000\n";
        os << "    sample_frequency: 100\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_a << "\n";
        os << "                    beta: " << f1_b << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_a << "\n";
        os << "                    beta: " << f2_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_a << "\n";
        os << "                    beta: " << f3_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "--ignore-data";
        char arg4[] = "--nthreads";
        char arg5[] = "2";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
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

        unsigned int expected_sample_size = 20001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");
        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        SampleSummarizer<double> size1_summary_a = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size1_summary_b = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size1_summary_c = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> size2_summary_a = prior_sample.summarize<double>("pop_size_pop1");
        SampleSummarizer<double> size2_summary_b = prior_sample.summarize<double>("pop_size_pop2");
        SampleSummarizer<double> size2_summary_c = prior_sample.summarize<double>("pop_size_root_pop1");
        SampleSummarizer<double> size3_summary_a = prior_sample.summarize<double>("pop_size_pop1b");
        SampleSummarizer<double> size3_summary_b = prior_sample.summarize<double>("pop_size_pop2b");
        SampleSummarizer<double> size3_summary_c = prior_sample.summarize<double>("pop_size_root_pop1b");

        REQUIRE(size1_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));

        REQUIRE(size2_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));

        REQUIRE(size3_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));

        SampleSummarizer<double> f1_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> f2_summary = prior_sample.summarize<double>("freq_1_pop1");
        SampleSummarizer<double> f3_summary = prior_sample.summarize<double>("freq_1_pop1b");
        REQUIRE(f1_summary.sample_size() == expected_sample_size);
        REQUIRE(f2_summary.sample_size() == expected_sample_size);
        REQUIRE(f3_summary.sample_size() == expected_sample_size);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.01));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.01));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.01));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.01));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.01));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.01));

        SampleSummarizer<double> mult1_summary = prior_sample.summarize<double>("mutation_rate_kya");
        SampleSummarizer<double> mult2_summary = prior_sample.summarize<double>("mutation_rate_pop1");
        SampleSummarizer<double> mult3_summary = prior_sample.summarize<double>("mutation_rate_pop1b");
        REQUIRE(mult1_summary.sample_size() == expected_sample_size);
        REQUIRE(mult2_summary.sample_size() == expected_sample_size);
        REQUIRE(mult3_summary.sample_size() == expected_sample_size);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.01));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.01));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.01));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.01));

        SampleSummarizer<double> conc_summary = prior_sample.summarize<double>("concentration");
        REQUIRE(conc_summary.sample_size() == expected_sample_size);
        REQUIRE(conc_summary.mean() == Approx(concentration_shape * concentration_scale).epsilon(0.01));
        REQUIRE(conc_summary.variance() == Approx(concentration_shape * concentration_scale * concentration_scale).epsilon(0.01));

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        REQUIRE((model_counts.at("000") / (double)expected_sample_size) == Approx(0.367).epsilon(0.01));
        REQUIRE((model_counts.at("012") / (double)expected_sample_size) == Approx(0.163).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing ReversibleJumpSampler with 2 pairs", "[SamplingPrior]") {

    SECTION("Testing rjMCMC with 2 pairs") {
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
        os << "    uniform:\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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
        char arg2[] = "6623";
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

        unsigned int expected_sample_size = 10001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");

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
        REQUIRE(total == expected_sample_size);
        REQUIRE((nshared / (double)expected_sample_size) == Approx(0.5).epsilon(0.01));

        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

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

TEST_CASE("Testing ReversibleJumpSampler with 3 pairs", "[SamplingPrior]") {

    SECTION("Testing rjMCMC with 3 pairs") {
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
        os << "    uniform:\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 100000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
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

        unsigned int expected_sample_size = 10001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)expected_sample_size << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)expected_sample_size << "\n";
        }
        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)expected_sample_size) == Approx(1.0/model_counts.size()).epsilon(0.01));
        }

        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

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

TEST_CASE("Testing ReversibleJumpSampler with 6 pairs", "[SamplingPrior]") {

    SECTION("Testing rjMCMC with 6 pairs") {
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
        os << "    uniform:\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 2000000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname4.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname5.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "68165454";
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

        unsigned int expected_sample_size = 200001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");
        SampleSummarizer<double> height_summary4 = prior_sample.summarize<double>("root_height_pop1c");
        SampleSummarizer<double> height_summary5 = prior_sample.summarize<double>("root_height_pop1d");
        SampleSummarizer<double> height_summary6 = prior_sample.summarize<double>("root_height_pop1e");

        std::vector<int> nevents = prior_sample.get<int>("number_of_events");
        std::vector<int> event_indices1 = prior_sample.get<int>("root_height_index_kya");
        std::vector<int> event_indices2 = prior_sample.get<int>("root_height_index_pop1");
        std::vector<int> event_indices3 = prior_sample.get<int>("root_height_index_pop1b");
        std::vector<int> event_indices4 = prior_sample.get<int>("root_height_index_pop1c");
        std::vector<int> event_indices5 = prior_sample.get<int>("root_height_index_pop1d");
        std::vector<int> event_indices6 = prior_sample.get<int>("root_height_index_pop1e");
        std::vector<double> heights1 = prior_sample.get<double>("root_height_kya");
        std::vector<double> heights2 = prior_sample.get<double>("root_height_pop1");
        std::vector<double> heights3 = prior_sample.get<double>("root_height_pop1b");
        std::vector<double> heights4 = prior_sample.get<double>("root_height_pop1c");
        std::vector<double> heights5 = prior_sample.get<double>("root_height_pop1d");
        std::vector<double> heights6 = prior_sample.get<double>("root_height_pop1e");

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
                {4, 0},
                {5, 0},
                {6, 0}
        };
        for (size_t i = 0; i < nevents.size(); ++i) {
            std::ostringstream stream;
            stream << event_indices1.at(i);
            stream << event_indices2.at(i);
            stream << event_indices3.at(i);
            stream << event_indices4.at(i);
            stream << event_indices5.at(i);
            stream << event_indices6.at(i);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
            REQUIRE(nevent_counts.count(nevents.at(i)) == 1);
            ++nevent_counts[nevents.at(i)];
            if (nevents.at(i) == 1) {
                REQUIRE(event_indices1.at(i) == event_indices2.at(i));
                REQUIRE(event_indices1.at(i) == event_indices3.at(i));
                REQUIRE(event_indices1.at(i) == event_indices4.at(i));
                REQUIRE(event_indices1.at(i) == event_indices5.at(i));
                REQUIRE(event_indices1.at(i) == event_indices6.at(i));
                REQUIRE(heights1.at(i) == heights2.at(i));
                REQUIRE(heights1.at(i) == heights3.at(i));
                REQUIRE(heights1.at(i) == heights4.at(i));
                REQUIRE(heights1.at(i) == heights5.at(i));
                REQUIRE(heights1.at(i) == heights6.at(i));
            }
            else if (nevents.at(i) == 6) {
                REQUIRE(event_indices1.at(i) != event_indices2.at(i));
                REQUIRE(event_indices1.at(i) != event_indices3.at(i));
                REQUIRE(event_indices1.at(i) != event_indices4.at(i));
                REQUIRE(event_indices1.at(i) != event_indices5.at(i));
                REQUIRE(event_indices1.at(i) != event_indices6.at(i));
                REQUIRE(event_indices2.at(i) != event_indices3.at(i));
                REQUIRE(event_indices2.at(i) != event_indices4.at(i));
                REQUIRE(event_indices2.at(i) != event_indices5.at(i));
                REQUIRE(event_indices2.at(i) != event_indices6.at(i));
                REQUIRE(event_indices3.at(i) != event_indices4.at(i));
                REQUIRE(event_indices3.at(i) != event_indices5.at(i));
                REQUIRE(event_indices3.at(i) != event_indices6.at(i));
                REQUIRE(event_indices4.at(i) != event_indices5.at(i));
                REQUIRE(event_indices4.at(i) != event_indices6.at(i));
                REQUIRE(event_indices5.at(i) != event_indices6.at(i));
                REQUIRE(heights1.at(i) != heights2.at(i));
                REQUIRE(heights1.at(i) != heights3.at(i));
                REQUIRE(heights1.at(i) != heights4.at(i));
                REQUIRE(heights1.at(i) != heights5.at(i));
                REQUIRE(heights1.at(i) != heights6.at(i));
                REQUIRE(heights2.at(i) != heights3.at(i));
                REQUIRE(heights2.at(i) != heights4.at(i));
                REQUIRE(heights2.at(i) != heights5.at(i));
                REQUIRE(heights2.at(i) != heights6.at(i));
                REQUIRE(heights3.at(i) != heights4.at(i));
                REQUIRE(heights3.at(i) != heights5.at(i));
                REQUIRE(heights3.at(i) != heights6.at(i));
                REQUIRE(heights4.at(i) != heights5.at(i));
                REQUIRE(heights4.at(i) != heights6.at(i));
                REQUIRE(heights5.at(i) != heights6.at(i));
            }
        }
        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        REQUIRE(model_counts.at("000000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012345") == nevent_counts.at(6));

        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)expected_sample_size << "\n";
        }
        for (std::string m: {"000000", "012345", "012344", "012340", "001234", "012314"}) {
            std::cout << m << ": " << model_counts.at(m) / (double)expected_sample_size << "\n";
        }

        for (std::string m: {"000000", "012345", "012344", "012340", "001234", "012314"}) {
            REQUIRE((model_counts.at(m) / (double)expected_sample_size) == Approx(1.0/bell_float(6)).epsilon(0.002));
        }
        for (auto const & kv: nevent_counts) {
            REQUIRE((kv.second / (double)expected_sample_size) == Approx(stirling2_float(6, kv.first)/bell_float(6)).epsilon(0.002));
        }

        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary4.sample_size() == expected_sample_size);
        REQUIRE(height_summary4.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary4.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary5.sample_size() == expected_sample_size);
        REQUIRE(height_summary5.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary5.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary6.sample_size() == expected_sample_size);
        REQUIRE(height_summary6.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary6.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));

        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing ReversibleJumpSampler with 6 pairs and diffuse gamma", "[SamplingPrior]") {

    SECTION("Testing rjMCMC with 6 pairs and diffuse gamma") {
        double height_shape = 1.0;
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
        os << "    uniform:\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 1000000\n";
        os << "    sample_frequency: 10\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.3\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.2\n";
        os << "            weight: 0.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 0.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    equal_state_frequencies: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.005\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname4.nex\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname5.nex\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "4826624";
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

        unsigned int expected_sample_size = 100001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");
        SampleSummarizer<double> height_summary4 = prior_sample.summarize<double>("root_height_pop1c");
        SampleSummarizer<double> height_summary5 = prior_sample.summarize<double>("root_height_pop1d");
        SampleSummarizer<double> height_summary6 = prior_sample.summarize<double>("root_height_pop1e");

        std::vector<int> nevents = prior_sample.get<int>("number_of_events");
        std::vector<int> event_indices1 = prior_sample.get<int>("root_height_index_kya");
        std::vector<int> event_indices2 = prior_sample.get<int>("root_height_index_pop1");
        std::vector<int> event_indices3 = prior_sample.get<int>("root_height_index_pop1b");
        std::vector<int> event_indices4 = prior_sample.get<int>("root_height_index_pop1c");
        std::vector<int> event_indices5 = prior_sample.get<int>("root_height_index_pop1d");
        std::vector<int> event_indices6 = prior_sample.get<int>("root_height_index_pop1e");
        std::vector<double> heights1 = prior_sample.get<double>("root_height_kya");
        std::vector<double> heights2 = prior_sample.get<double>("root_height_pop1");
        std::vector<double> heights3 = prior_sample.get<double>("root_height_pop1b");
        std::vector<double> heights4 = prior_sample.get<double>("root_height_pop1c");
        std::vector<double> heights5 = prior_sample.get<double>("root_height_pop1d");
        std::vector<double> heights6 = prior_sample.get<double>("root_height_pop1e");

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
                {4, 0},
                {5, 0},
                {6, 0}
        };
        for (size_t i = 0; i < nevents.size(); ++i) {
            std::ostringstream stream;
            stream << event_indices1.at(i);
            stream << event_indices2.at(i);
            stream << event_indices3.at(i);
            stream << event_indices4.at(i);
            stream << event_indices5.at(i);
            stream << event_indices6.at(i);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
            REQUIRE(nevent_counts.count(nevents.at(i)) == 1);
            ++nevent_counts[nevents.at(i)];
            if (nevents.at(i) == 1) {
                REQUIRE(event_indices1.at(i) == event_indices2.at(i));
                REQUIRE(event_indices1.at(i) == event_indices3.at(i));
                REQUIRE(event_indices1.at(i) == event_indices4.at(i));
                REQUIRE(event_indices1.at(i) == event_indices5.at(i));
                REQUIRE(event_indices1.at(i) == event_indices6.at(i));
                REQUIRE(heights1.at(i) == heights2.at(i));
                REQUIRE(heights1.at(i) == heights3.at(i));
                REQUIRE(heights1.at(i) == heights4.at(i));
                REQUIRE(heights1.at(i) == heights5.at(i));
                REQUIRE(heights1.at(i) == heights6.at(i));
            }
            else if (nevents.at(i) == 6) {
                REQUIRE(event_indices1.at(i) != event_indices2.at(i));
                REQUIRE(event_indices1.at(i) != event_indices3.at(i));
                REQUIRE(event_indices1.at(i) != event_indices4.at(i));
                REQUIRE(event_indices1.at(i) != event_indices5.at(i));
                REQUIRE(event_indices1.at(i) != event_indices6.at(i));
                REQUIRE(event_indices2.at(i) != event_indices3.at(i));
                REQUIRE(event_indices2.at(i) != event_indices4.at(i));
                REQUIRE(event_indices2.at(i) != event_indices5.at(i));
                REQUIRE(event_indices2.at(i) != event_indices6.at(i));
                REQUIRE(event_indices3.at(i) != event_indices4.at(i));
                REQUIRE(event_indices3.at(i) != event_indices5.at(i));
                REQUIRE(event_indices3.at(i) != event_indices6.at(i));
                REQUIRE(event_indices4.at(i) != event_indices5.at(i));
                REQUIRE(event_indices4.at(i) != event_indices6.at(i));
                REQUIRE(event_indices5.at(i) != event_indices6.at(i));
                REQUIRE(heights1.at(i) != heights2.at(i));
                REQUIRE(heights1.at(i) != heights3.at(i));
                REQUIRE(heights1.at(i) != heights4.at(i));
                REQUIRE(heights1.at(i) != heights5.at(i));
                REQUIRE(heights1.at(i) != heights6.at(i));
                REQUIRE(heights2.at(i) != heights3.at(i));
                REQUIRE(heights2.at(i) != heights4.at(i));
                REQUIRE(heights2.at(i) != heights5.at(i));
                REQUIRE(heights2.at(i) != heights6.at(i));
                REQUIRE(heights3.at(i) != heights4.at(i));
                REQUIRE(heights3.at(i) != heights5.at(i));
                REQUIRE(heights3.at(i) != heights6.at(i));
                REQUIRE(heights4.at(i) != heights5.at(i));
                REQUIRE(heights4.at(i) != heights6.at(i));
                REQUIRE(heights5.at(i) != heights6.at(i));
            }
        }
        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        REQUIRE(model_counts.at("000000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012345") == nevent_counts.at(6));

        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)expected_sample_size << "\n";
        }
        for (std::string m: {"000000", "012345", "012344", "012340", "001234", "012314"}) {
            std::cout << m << ": " << model_counts.at(m) / (double)expected_sample_size << "\n";
        }

        for (std::string m: {"000000", "012345", "012344", "012340", "001234", "012314"}) {
            REQUIRE((model_counts.at(m) / (double)expected_sample_size) == Approx(1.0/bell_float(6)).epsilon(0.002));
        }
        for (auto const & kv: nevent_counts) {
            REQUIRE((kv.second / (double)expected_sample_size) == Approx(stirling2_float(6, kv.first)/bell_float(6)).epsilon(0.002));
        }

        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary4.sample_size() == expected_sample_size);
        REQUIRE(height_summary4.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary4.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary5.sample_size() == expected_sample_size);
        REQUIRE(height_summary5.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary5.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary6.sample_size() == expected_sample_size);
        REQUIRE(height_summary6.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary6.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));

        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing ReversibleJumpSampler with 3 pairs and fully parameterized", "[SamplingPrior]") {

    SECTION("Testing rjMCMC with 3 pairs") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_a = 2.0;
        double f1_b = 1.1;
        double f2_a = 1.0;
        double f2_b = 0.5;
        double f3_a = 1.5;
        double f3_b = 1.8;
        double expected_f1_mean = f1_a / (f1_a + f1_b);
        double expected_f1_variance = (f1_a * f1_b) / ((f1_a + f1_b) * (f1_a + f1_b) * (f1_a + f1_b + 1.0));
        double expected_f2_mean = f2_a / (f2_a + f2_b);
        double expected_f2_variance = (f2_a * f2_b) / ((f2_a + f2_b) * (f2_a + f2_b) * (f2_a + f2_b + 1.0));
        double expected_f3_mean = f3_a / (f3_a + f3_b);
        double expected_f3_variance = (f3_a * f3_b) / ((f3_a + f3_b) * (f3_a + f3_b) * (f3_a + f3_b + 1.0));
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
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
        os << "    uniform:\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 2000000\n";
        os << "    sample_frequency: 100\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_a << "\n";
        os << "                    beta: " << f1_b << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_a << "\n";
        os << "                    beta: " << f2_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_a << "\n";
        os << "                    beta: " << f3_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "58961543";
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

        unsigned int expected_sample_size = 20001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");

        SampleSummarizer<double> size1_summary_a = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size1_summary_b = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size1_summary_c = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> size2_summary_a = prior_sample.summarize<double>("pop_size_pop1");
        SampleSummarizer<double> size2_summary_b = prior_sample.summarize<double>("pop_size_pop2");
        SampleSummarizer<double> size2_summary_c = prior_sample.summarize<double>("pop_size_root_pop1");
        SampleSummarizer<double> size3_summary_a = prior_sample.summarize<double>("pop_size_pop1b");
        SampleSummarizer<double> size3_summary_b = prior_sample.summarize<double>("pop_size_pop2b");
        SampleSummarizer<double> size3_summary_c = prior_sample.summarize<double>("pop_size_root_pop1b");

        SampleSummarizer<double> f1_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> f2_summary = prior_sample.summarize<double>("freq_1_pop1");
        SampleSummarizer<double> f3_summary = prior_sample.summarize<double>("freq_1_pop1b");

        SampleSummarizer<double> mult1_summary = prior_sample.summarize<double>("mutation_rate_kya");
        SampleSummarizer<double> mult2_summary = prior_sample.summarize<double>("mutation_rate_pop1");
        SampleSummarizer<double> mult3_summary = prior_sample.summarize<double>("mutation_rate_pop1b");

        REQUIRE_THROWS(prior_sample.summarize<double>("concentration"));

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)expected_sample_size << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)expected_sample_size << "\n";
        }

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (auto const & kv: model_counts) {
            REQUIRE((kv.second / (double)expected_sample_size) == Approx(1.0/bell_float(3)).epsilon(0.01));
        }
        for (auto const & kv: nevent_counts) {
            REQUIRE((kv.second / (double)expected_sample_size) == Approx(stirling2_float(3, kv.first)/bell_float(3)).epsilon(0.005));
        }

        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(size1_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));

        REQUIRE(size2_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));

        REQUIRE(size3_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));

        REQUIRE(f1_summary.sample_size() == expected_sample_size);
        REQUIRE(f2_summary.sample_size() == expected_sample_size);
        REQUIRE(f3_summary.sample_size() == expected_sample_size);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.01));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.01));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.01));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.01));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.01));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.01));

        REQUIRE(mult1_summary.sample_size() == expected_sample_size);
        REQUIRE(mult2_summary.sample_size() == expected_sample_size);
        REQUIRE(mult3_summary.sample_size() == expected_sample_size);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.01));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.01));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.01));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing ReversibleJumpSampler with 3 pairs and fully parameterized and CompositeHeightSizeRateMixer",
        "[SamplingPrior]") {

    SECTION("Testing rjMCMC with 3 pairs") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_a = 2.0;
        double f1_b = 1.1;
        double f2_a = 1.0;
        double f2_b = 0.5;
        double f3_a = 1.5;
        double f3_b = 1.8;
        double expected_f1_mean = f1_a / (f1_a + f1_b);
        double expected_f1_variance = (f1_a * f1_b) / ((f1_a + f1_b) * (f1_a + f1_b) * (f1_a + f1_b + 1.0));
        double expected_f2_mean = f2_a / (f2_a + f2_b);
        double expected_f2_variance = (f2_a * f2_b) / ((f2_a + f2_b) * (f2_a + f2_b) * (f2_a + f2_b + 1.0));
        double expected_f3_mean = f3_a / (f3_a + f3_b);
        double expected_f3_variance = (f3_a * f3_b) / ((f3_a + f3_b) * (f3_a + f3_b) * (f3_a + f3_b + 1.0));
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
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
        os << "    uniform:\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 2000000\n";
        os << "    sample_frequency: 20\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_a << "\n";
        os << "                    beta: " << f1_b << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_a << "\n";
        os << "                    beta: " << f2_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_a << "\n";
        os << "                    beta: " << f3_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "589615439";
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

        unsigned int expected_sample_size = 100001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");

        SampleSummarizer<double> size1_summary_a = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size1_summary_b = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size1_summary_c = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> size2_summary_a = prior_sample.summarize<double>("pop_size_pop1");
        SampleSummarizer<double> size2_summary_b = prior_sample.summarize<double>("pop_size_pop2");
        SampleSummarizer<double> size2_summary_c = prior_sample.summarize<double>("pop_size_root_pop1");
        SampleSummarizer<double> size3_summary_a = prior_sample.summarize<double>("pop_size_pop1b");
        SampleSummarizer<double> size3_summary_b = prior_sample.summarize<double>("pop_size_pop2b");
        SampleSummarizer<double> size3_summary_c = prior_sample.summarize<double>("pop_size_root_pop1b");

        SampleSummarizer<double> f1_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> f2_summary = prior_sample.summarize<double>("freq_1_pop1");
        SampleSummarizer<double> f3_summary = prior_sample.summarize<double>("freq_1_pop1b");

        SampleSummarizer<double> mult1_summary = prior_sample.summarize<double>("mutation_rate_kya");
        SampleSummarizer<double> mult2_summary = prior_sample.summarize<double>("mutation_rate_pop1");
        SampleSummarizer<double> mult3_summary = prior_sample.summarize<double>("mutation_rate_pop1b");

        REQUIRE_THROWS(prior_sample.summarize<double>("concentration"));

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)expected_sample_size << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)expected_sample_size << "\n";
        }

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (auto const & kv: model_counts) {
            REQUIRE((kv.second / (double)expected_sample_size) == Approx(1.0/bell_float(3)).epsilon(0.01));
        }
        for (auto const & kv: nevent_counts) {
            REQUIRE((kv.second / (double)expected_sample_size) == Approx(stirling2_float(3, kv.first)/bell_float(3)).epsilon(0.005));
        }

        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(size1_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));

        REQUIRE(size2_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));

        REQUIRE(size3_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));

        REQUIRE(f1_summary.sample_size() == expected_sample_size);
        REQUIRE(f2_summary.sample_size() == expected_sample_size);
        REQUIRE(f3_summary.sample_size() == expected_sample_size);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.01));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.01));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.01));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.01));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.01));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.01));

        REQUIRE(mult1_summary.sample_size() == expected_sample_size);
        REQUIRE(mult2_summary.sample_size() == expected_sample_size);
        REQUIRE(mult3_summary.sample_size() == expected_sample_size);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.01));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.01));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.01));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing ReversibleJumpSampler with 3 pairs, fully parameterized, and 2 threads",
        "[SamplingPrior]") {

    SECTION("Testing rjMCMC with 3 pairs and 2 threads") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_a = 2.0;
        double f1_b = 1.1;
        double f2_a = 1.0;
        double f2_b = 0.5;
        double f3_a = 1.5;
        double f3_b = 1.8;
        double expected_f1_mean = f1_a / (f1_a + f1_b);
        double expected_f1_variance = (f1_a * f1_b) / ((f1_a + f1_b) * (f1_a + f1_b) * (f1_a + f1_b + 1.0));
        double expected_f2_mean = f2_a / (f2_a + f2_b);
        double expected_f2_variance = (f2_a * f2_b) / ((f2_a + f2_b) * (f2_a + f2_b) * (f2_a + f2_b + 1.0));
        double expected_f3_mean = f3_a / (f3_a + f3_b);
        double expected_f3_variance = (f3_a * f3_b) / ((f3_a + f3_b) * (f3_a + f3_b) * (f3_a + f3_b + 1.0));
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
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
        os << "    uniform:\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 2000000\n";
        os << "    sample_frequency: 100\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_a << "\n";
        os << "                    beta: " << f1_b << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_a << "\n";
        os << "                    beta: " << f2_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_a << "\n";
        os << "                    beta: " << f3_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "571518456";
        char arg3[] = "--ignore-data";
        char arg4[] = "--nthreads";
        char arg5[] = "2";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
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

        unsigned int expected_sample_size = 20001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");

        SampleSummarizer<double> size1_summary_a = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size1_summary_b = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size1_summary_c = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> size2_summary_a = prior_sample.summarize<double>("pop_size_pop1");
        SampleSummarizer<double> size2_summary_b = prior_sample.summarize<double>("pop_size_pop2");
        SampleSummarizer<double> size2_summary_c = prior_sample.summarize<double>("pop_size_root_pop1");
        SampleSummarizer<double> size3_summary_a = prior_sample.summarize<double>("pop_size_pop1b");
        SampleSummarizer<double> size3_summary_b = prior_sample.summarize<double>("pop_size_pop2b");
        SampleSummarizer<double> size3_summary_c = prior_sample.summarize<double>("pop_size_root_pop1b");

        SampleSummarizer<double> f1_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> f2_summary = prior_sample.summarize<double>("freq_1_pop1");
        SampleSummarizer<double> f3_summary = prior_sample.summarize<double>("freq_1_pop1b");

        SampleSummarizer<double> mult1_summary = prior_sample.summarize<double>("mutation_rate_kya");
        SampleSummarizer<double> mult2_summary = prior_sample.summarize<double>("mutation_rate_pop1");
        SampleSummarizer<double> mult3_summary = prior_sample.summarize<double>("mutation_rate_pop1b");

        REQUIRE_THROWS(prior_sample.summarize<double>("concentration"));

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)expected_sample_size << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)expected_sample_size << "\n";
        }

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (auto const & kv: model_counts) {
            REQUIRE((kv.second / (double)expected_sample_size) == Approx(1.0/bell_float(3)).epsilon(0.005));
        }
        for (auto const & kv: nevent_counts) {
            REQUIRE((kv.second / (double)expected_sample_size) == Approx(stirling2_float(3, kv.first)/bell_float(3)).epsilon(0.005));
        }

        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(size1_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));

        REQUIRE(size2_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));

        REQUIRE(size3_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));

        REQUIRE(f1_summary.sample_size() == expected_sample_size);
        REQUIRE(f2_summary.sample_size() == expected_sample_size);
        REQUIRE(f3_summary.sample_size() == expected_sample_size);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.01));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.01));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.01));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.01));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.01));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.01));

        REQUIRE(mult1_summary.sample_size() == expected_sample_size);
        REQUIRE(mult2_summary.sample_size() == expected_sample_size);
        REQUIRE(mult3_summary.sample_size() == expected_sample_size);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.01));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.01));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.01));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing DPP with 2 singletons and 1 pair, fully parameterized",
        "[SamplingPrior]") {

    SECTION("Testing alpha integrated") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_a = 2.0;
        double f1_b = 1.1;
        double f2_a = 1.0;
        double f2_b = 0.5;
        double f3_a = 1.5;
        double f3_b = 1.8;
        double expected_f1_mean = f1_a / (f1_a + f1_b);
        double expected_f1_variance = (f1_a * f1_b) / ((f1_a + f1_b) * (f1_a + f1_b) * (f1_a + f1_b + 1.0));
        double expected_f2_mean = f2_a / (f2_a + f2_b);
        double expected_f2_variance = (f2_a * f2_b) / ((f2_a + f2_b) * (f2_a + f2_b) * (f2_a + f2_b + 1.0));
        double expected_f3_mean = f3_a / (f3_a + f3_b);
        double expected_f3_variance = (f3_a * f3_b) / ((f3_a + f3_b) * (f3_a + f3_b) * (f3_a + f3_b + 1.0));
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
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
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 2000000\n";
        os << "    sample_frequency: 100\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_a << "\n";
        os << "                    beta: " << f1_b << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_a << "\n";
        os << "                    beta: " << f2_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_a << "\n";
        os << "                    beta: " << f3_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
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

        unsigned int expected_sample_size = 20001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");
        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        SampleSummarizer<double> size1_summary_a = prior_sample.summarize<double>("pop_size_kya");
        REQUIRE_THROWS(prior_sample.summarize<double>("pop_size_fas"));
        SampleSummarizer<double> size1_summary_c = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> size2_summary_a = prior_sample.summarize<double>("pop_size_pop1");
        REQUIRE_THROWS(prior_sample.summarize<double>("pop_size_pop2"));
        SampleSummarizer<double> size2_summary_c = prior_sample.summarize<double>("pop_size_root_pop1");
        SampleSummarizer<double> size3_summary_a = prior_sample.summarize<double>("pop_size_pop1b");
        SampleSummarizer<double> size3_summary_b = prior_sample.summarize<double>("pop_size_pop2b");
        SampleSummarizer<double> size3_summary_c = prior_sample.summarize<double>("pop_size_root_pop1b");

        REQUIRE(size1_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));

        REQUIRE(size2_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));

        REQUIRE(size3_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));

        SampleSummarizer<double> f1_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> f2_summary = prior_sample.summarize<double>("freq_1_pop1");
        SampleSummarizer<double> f3_summary = prior_sample.summarize<double>("freq_1_pop1b");
        REQUIRE(f1_summary.sample_size() == expected_sample_size);
        REQUIRE(f2_summary.sample_size() == expected_sample_size);
        REQUIRE(f3_summary.sample_size() == expected_sample_size);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.01));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.01));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.01));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.01));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.01));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.01));

        SampleSummarizer<double> mult1_summary = prior_sample.summarize<double>("mutation_rate_kya");
        SampleSummarizer<double> mult2_summary = prior_sample.summarize<double>("mutation_rate_pop1");
        SampleSummarizer<double> mult3_summary = prior_sample.summarize<double>("mutation_rate_pop1b");
        REQUIRE(mult1_summary.sample_size() == expected_sample_size);
        REQUIRE(mult2_summary.sample_size() == expected_sample_size);
        REQUIRE(mult3_summary.sample_size() == expected_sample_size);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.01));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.01));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.01));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.01));

        SampleSummarizer<double> conc_summary = prior_sample.summarize<double>("concentration");
        REQUIRE(conc_summary.sample_size() == expected_sample_size);
        REQUIRE(conc_summary.mean() == Approx(concentration_shape * concentration_scale).epsilon(0.01));
        REQUIRE(conc_summary.variance() == Approx(concentration_shape * concentration_scale * concentration_scale).epsilon(0.01));

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        REQUIRE((model_counts.at("000") / (double)expected_sample_size) == Approx(0.367).epsilon(0.01));
        REQUIRE((model_counts.at("012") / (double)expected_sample_size) == Approx(0.163).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing DPP with 2 singletons and 1 pair, fully parameterized, and CompositeHeightSizeRateMixer",
        "[SamplingPrior]") {

    SECTION("Testing alpha integrated") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_a = 2.0;
        double f1_b = 1.1;
        double f2_a = 1.0;
        double f2_b = 0.5;
        double f3_a = 1.5;
        double f3_b = 1.8;
        double expected_f1_mean = f1_a / (f1_a + f1_b);
        double expected_f1_variance = (f1_a * f1_b) / ((f1_a + f1_b) * (f1_a + f1_b) * (f1_a + f1_b + 1.0));
        double expected_f2_mean = f2_a / (f2_a + f2_b);
        double expected_f2_variance = (f2_a * f2_b) / ((f2_a + f2_b) * (f2_a + f2_b) * (f2_a + f2_b + 1.0));
        double expected_f3_mean = f3_a / (f3_a + f3_b);
        double expected_f3_variance = (f3_a * f3_b) / ((f3_a + f3_b) * (f3_a + f3_b) * (f3_a + f3_b + 1.0));
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
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
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 2000000\n";
        os << "    sample_frequency: 100\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_a << "\n";
        os << "                    beta: " << f1_b << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_a << "\n";
        os << "                    beta: " << f2_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_a << "\n";
        os << "                    beta: " << f3_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
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

        unsigned int expected_sample_size = 20001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");
        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        SampleSummarizer<double> size1_summary_a = prior_sample.summarize<double>("pop_size_kya");
        REQUIRE_THROWS(prior_sample.summarize<double>("pop_size_fas"));
        SampleSummarizer<double> size1_summary_c = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> size2_summary_a = prior_sample.summarize<double>("pop_size_pop1");
        REQUIRE_THROWS(prior_sample.summarize<double>("pop_size_pop2"));
        SampleSummarizer<double> size2_summary_c = prior_sample.summarize<double>("pop_size_root_pop1");
        SampleSummarizer<double> size3_summary_a = prior_sample.summarize<double>("pop_size_pop1b");
        SampleSummarizer<double> size3_summary_b = prior_sample.summarize<double>("pop_size_pop2b");
        SampleSummarizer<double> size3_summary_c = prior_sample.summarize<double>("pop_size_root_pop1b");

        REQUIRE(size1_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));

        REQUIRE(size2_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));

        REQUIRE(size3_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));

        SampleSummarizer<double> f1_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> f2_summary = prior_sample.summarize<double>("freq_1_pop1");
        SampleSummarizer<double> f3_summary = prior_sample.summarize<double>("freq_1_pop1b");
        REQUIRE(f1_summary.sample_size() == expected_sample_size);
        REQUIRE(f2_summary.sample_size() == expected_sample_size);
        REQUIRE(f3_summary.sample_size() == expected_sample_size);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.01));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.01));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.01));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.01));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.01));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.01));

        SampleSummarizer<double> mult1_summary = prior_sample.summarize<double>("mutation_rate_kya");
        SampleSummarizer<double> mult2_summary = prior_sample.summarize<double>("mutation_rate_pop1");
        SampleSummarizer<double> mult3_summary = prior_sample.summarize<double>("mutation_rate_pop1b");
        REQUIRE(mult1_summary.sample_size() == expected_sample_size);
        REQUIRE(mult2_summary.sample_size() == expected_sample_size);
        REQUIRE(mult3_summary.sample_size() == expected_sample_size);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.01));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.01));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.01));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.01));

        SampleSummarizer<double> conc_summary = prior_sample.summarize<double>("concentration");
        REQUIRE(conc_summary.sample_size() == expected_sample_size);
        REQUIRE(conc_summary.mean() == Approx(concentration_shape * concentration_scale).epsilon(0.01));
        REQUIRE(conc_summary.variance() == Approx(concentration_shape * concentration_scale * concentration_scale).epsilon(0.01));

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        REQUIRE((model_counts.at("000") / (double)expected_sample_size) == Approx(0.367).epsilon(0.01));
        REQUIRE((model_counts.at("012") / (double)expected_sample_size) == Approx(0.163).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing ReversibleJumpSampler with 2 singletons, 1 pair, and fully parameterized",
        "[SamplingPrior]") {

    SECTION("Testing rjMCMC with 2 singletons and a pair") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_a = 2.0;
        double f1_b = 1.1;
        double f2_a = 1.0;
        double f2_b = 0.5;
        double f3_a = 1.5;
        double f3_b = 1.8;
        double expected_f1_mean = f1_a / (f1_a + f1_b);
        double expected_f1_variance = (f1_a * f1_b) / ((f1_a + f1_b) * (f1_a + f1_b) * (f1_a + f1_b + 1.0));
        double expected_f2_mean = f2_a / (f2_a + f2_b);
        double expected_f2_variance = (f2_a * f2_b) / ((f2_a + f2_b) * (f2_a + f2_b) * (f2_a + f2_b + 1.0));
        double expected_f3_mean = f3_a / (f3_a + f3_b);
        double expected_f3_variance = (f3_a * f3_b) / ((f3_a + f3_b) * (f3_a + f3_b) * (f3_a + f3_b + 1.0));
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
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
        os << "    uniform:\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 2000000\n";
        os << "    sample_frequency: 100\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_a << "\n";
        os << "                    beta: " << f1_b << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_a << "\n";
        os << "                    beta: " << f2_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_a << "\n";
        os << "                    beta: " << f3_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "68412635";
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

        unsigned int expected_sample_size = 20001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");

        SampleSummarizer<double> size1_summary_a = prior_sample.summarize<double>("pop_size_kya");
        REQUIRE_THROWS(prior_sample.summarize<double>("pop_size_fas"));
        SampleSummarizer<double> size1_summary_c = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> size2_summary_a = prior_sample.summarize<double>("pop_size_pop1");
        REQUIRE_THROWS(prior_sample.summarize<double>("pop_size_pop2"));
        SampleSummarizer<double> size2_summary_c = prior_sample.summarize<double>("pop_size_root_pop1");
        SampleSummarizer<double> size3_summary_a = prior_sample.summarize<double>("pop_size_pop1b");
        SampleSummarizer<double> size3_summary_b = prior_sample.summarize<double>("pop_size_pop2b");
        SampleSummarizer<double> size3_summary_c = prior_sample.summarize<double>("pop_size_root_pop1b");

        SampleSummarizer<double> f1_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> f2_summary = prior_sample.summarize<double>("freq_1_pop1");
        SampleSummarizer<double> f3_summary = prior_sample.summarize<double>("freq_1_pop1b");

        SampleSummarizer<double> mult1_summary = prior_sample.summarize<double>("mutation_rate_kya");
        SampleSummarizer<double> mult2_summary = prior_sample.summarize<double>("mutation_rate_pop1");
        SampleSummarizer<double> mult3_summary = prior_sample.summarize<double>("mutation_rate_pop1b");

        REQUIRE_THROWS(prior_sample.summarize<double>("concentration"));

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)expected_sample_size << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)expected_sample_size << "\n";
        }

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (auto const & kv: model_counts) {
            REQUIRE((kv.second / (double)expected_sample_size) == Approx(1.0/bell_float(3)).epsilon(0.01));
        }
        for (auto const & kv: nevent_counts) {
            REQUIRE((kv.second / (double)expected_sample_size) == Approx(stirling2_float(3, kv.first)/bell_float(3)).epsilon(0.005));
        }

        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(size1_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));

        REQUIRE(size2_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));

        REQUIRE(size3_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));

        REQUIRE(f1_summary.sample_size() == expected_sample_size);
        REQUIRE(f2_summary.sample_size() == expected_sample_size);
        REQUIRE(f3_summary.sample_size() == expected_sample_size);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.01));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.01));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.01));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.01));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.01));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.01));

        REQUIRE(mult1_summary.sample_size() == expected_sample_size);
        REQUIRE(mult2_summary.sample_size() == expected_sample_size);
        REQUIRE(mult3_summary.sample_size() == expected_sample_size);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.01));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.01));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.01));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing ReversibleJumpSampler with 2 singletons, 1 pair, CompositeHeightSizeRateMixer, and fully parameterized",
        "[SamplingPrior]") {

    SECTION("Testing rjMCMC with 2 singletons and a pair") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_a = 2.0;
        double f1_b = 1.1;
        double f2_a = 1.0;
        double f2_b = 0.5;
        double f3_a = 1.5;
        double f3_b = 1.8;
        double expected_f1_mean = f1_a / (f1_a + f1_b);
        double expected_f1_variance = (f1_a * f1_b) / ((f1_a + f1_b) * (f1_a + f1_b) * (f1_a + f1_b + 1.0));
        double expected_f2_mean = f2_a / (f2_a + f2_b);
        double expected_f2_variance = (f2_a * f2_b) / ((f2_a + f2_b) * (f2_a + f2_b) * (f2_a + f2_b + 1.0));
        double expected_f3_mean = f3_a / (f3_a + f3_b);
        double expected_f3_variance = (f3_a * f3_b) / ((f3_a + f3_b) * (f3_a + f3_b) * (f3_a + f3_b + 1.0));
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
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
        os << "    uniform:\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 2000000\n";
        os << "    sample_frequency: 100\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_a << "\n";
        os << "                    beta: " << f1_b << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_a << "\n";
        os << "                    beta: " << f2_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_a << "\n";
        os << "                    beta: " << f3_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "68412635";
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

        unsigned int expected_sample_size = 20001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");

        SampleSummarizer<double> size1_summary_a = prior_sample.summarize<double>("pop_size_kya");
        REQUIRE_THROWS(prior_sample.summarize<double>("pop_size_fas"));
        SampleSummarizer<double> size1_summary_c = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> size2_summary_a = prior_sample.summarize<double>("pop_size_pop1");
        REQUIRE_THROWS(prior_sample.summarize<double>("pop_size_pop2"));
        SampleSummarizer<double> size2_summary_c = prior_sample.summarize<double>("pop_size_root_pop1");
        SampleSummarizer<double> size3_summary_a = prior_sample.summarize<double>("pop_size_pop1b");
        SampleSummarizer<double> size3_summary_b = prior_sample.summarize<double>("pop_size_pop2b");
        SampleSummarizer<double> size3_summary_c = prior_sample.summarize<double>("pop_size_root_pop1b");

        SampleSummarizer<double> f1_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> f2_summary = prior_sample.summarize<double>("freq_1_pop1");
        SampleSummarizer<double> f3_summary = prior_sample.summarize<double>("freq_1_pop1b");

        SampleSummarizer<double> mult1_summary = prior_sample.summarize<double>("mutation_rate_kya");
        SampleSummarizer<double> mult2_summary = prior_sample.summarize<double>("mutation_rate_pop1");
        SampleSummarizer<double> mult3_summary = prior_sample.summarize<double>("mutation_rate_pop1b");

        REQUIRE_THROWS(prior_sample.summarize<double>("concentration"));

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)expected_sample_size << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)expected_sample_size << "\n";
        }

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (auto const & kv: model_counts) {
            REQUIRE((kv.second / (double)expected_sample_size) == Approx(1.0/bell_float(3)).epsilon(0.01));
        }
        for (auto const & kv: nevent_counts) {
            REQUIRE((kv.second / (double)expected_sample_size) == Approx(stirling2_float(3, kv.first)/bell_float(3)).epsilon(0.005));
        }

        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(size1_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));

        REQUIRE(size2_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));

        REQUIRE(size3_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));

        REQUIRE(f1_summary.sample_size() == expected_sample_size);
        REQUIRE(f2_summary.sample_size() == expected_sample_size);
        REQUIRE(f3_summary.sample_size() == expected_sample_size);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.01));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.01));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.01));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.01));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.01));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.01));

        REQUIRE(mult1_summary.sample_size() == expected_sample_size);
        REQUIRE(mult2_summary.sample_size() == expected_sample_size);
        REQUIRE(mult3_summary.sample_size() == expected_sample_size);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.01));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.01));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.01));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.01));

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing fixed 012 and fully parameterized", "[SamplingPrior]") {

    SECTION("Testing fixed 012 pairs") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_a = 2.0;
        double f1_b = 1.1;
        double f2_a = 1.0;
        double f2_b = 0.5;
        double f3_a = 1.5;
        double f3_b = 1.8;
        double expected_f1_mean = f1_a / (f1_a + f1_b);
        double expected_f1_variance = (f1_a * f1_b) / ((f1_a + f1_b) * (f1_a + f1_b) * (f1_a + f1_b + 1.0));
        double expected_f2_mean = f2_a / (f2_a + f2_b);
        double expected_f2_variance = (f2_a * f2_b) / ((f2_a + f2_b) * (f2_a + f2_b) * (f2_a + f2_b + 1.0));
        double expected_f3_mean = f3_a / (f3_a + f3_b);
        double expected_f3_variance = (f3_a * f3_b) / ((f3_a + f3_b) * (f3_a + f3_b) * (f3_a + f3_b + 1.0));
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
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
        os << "    fixed: [0, 1, 2]\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 500000\n";
        os << "    sample_frequency: 100\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_a << "\n";
        os << "                    beta: " << f1_b << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_a << "\n";
        os << "                    beta: " << f2_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_a << "\n";
        os << "                    beta: " << f3_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "58961543";
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

        unsigned int expected_sample_size = 5001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");
        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.05));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.05));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.05));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.05));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.05));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.05));

        SampleSummarizer<double> size1_summary_a = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size1_summary_b = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size1_summary_c = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> size2_summary_a = prior_sample.summarize<double>("pop_size_pop1");
        SampleSummarizer<double> size2_summary_b = prior_sample.summarize<double>("pop_size_pop2");
        SampleSummarizer<double> size2_summary_c = prior_sample.summarize<double>("pop_size_root_pop1");
        SampleSummarizer<double> size3_summary_a = prior_sample.summarize<double>("pop_size_pop1b");
        SampleSummarizer<double> size3_summary_b = prior_sample.summarize<double>("pop_size_pop2b");
        SampleSummarizer<double> size3_summary_c = prior_sample.summarize<double>("pop_size_root_pop1b");

        REQUIRE(size1_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_b.mean() == Approx(size1_shape * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_b.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.05));

        REQUIRE(size2_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_b.mean() == Approx(size2_shape * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_b.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.05));

        REQUIRE(size3_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.05));

        SampleSummarizer<double> f1_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> f2_summary = prior_sample.summarize<double>("freq_1_pop1");
        SampleSummarizer<double> f3_summary = prior_sample.summarize<double>("freq_1_pop1b");
        REQUIRE(f1_summary.sample_size() == expected_sample_size);
        REQUIRE(f2_summary.sample_size() == expected_sample_size);
        REQUIRE(f3_summary.sample_size() == expected_sample_size);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.05));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.05));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.05));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.05));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.05));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.05));

        SampleSummarizer<double> mult1_summary = prior_sample.summarize<double>("mutation_rate_kya");
        SampleSummarizer<double> mult2_summary = prior_sample.summarize<double>("mutation_rate_pop1");
        SampleSummarizer<double> mult3_summary = prior_sample.summarize<double>("mutation_rate_pop1b");
        REQUIRE(mult1_summary.sample_size() == expected_sample_size);
        REQUIRE(mult2_summary.sample_size() == expected_sample_size);
        REQUIRE(mult3_summary.sample_size() == expected_sample_size);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.05));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.05));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.05));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.05));

        REQUIRE_THROWS(prior_sample.summarize<double>("concentration"));

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        REQUIRE(model_counts.at("000") == 0);
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == 0);

        REQUIRE(model_counts.at("012") == expected_sample_size);
        REQUIRE(nevent_counts.at(3) == expected_sample_size);

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing fixed 012 and fully parameterized and CompositeHeightSizeRateMixer",
        "[SamplingPrior]") {

    SECTION("Testing fixed 012 pairs") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_a = 2.0;
        double f1_b = 1.1;
        double f2_a = 1.0;
        double f2_b = 0.5;
        double f3_a = 1.5;
        double f3_b = 1.8;
        double expected_f1_mean = f1_a / (f1_a + f1_b);
        double expected_f1_variance = (f1_a * f1_b) / ((f1_a + f1_b) * (f1_a + f1_b) * (f1_a + f1_b + 1.0));
        double expected_f2_mean = f2_a / (f2_a + f2_b);
        double expected_f2_variance = (f2_a * f2_b) / ((f2_a + f2_b) * (f2_a + f2_b) * (f2_a + f2_b + 1.0));
        double expected_f3_mean = f3_a / (f3_a + f3_b);
        double expected_f3_variance = (f3_a * f3_b) / ((f3_a + f3_b) * (f3_a + f3_b) * (f3_a + f3_b + 1.0));
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
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
        os << "    fixed: [0, 1, 2]\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 500000\n";
        os << "    sample_frequency: 100\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_a << "\n";
        os << "                    beta: " << f1_b << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_a << "\n";
        os << "                    beta: " << f2_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_a << "\n";
        os << "                    beta: " << f3_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "58961543";
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

        unsigned int expected_sample_size = 5001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");
        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.05));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.05));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.05));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.05));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.05));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.05));

        SampleSummarizer<double> size1_summary_a = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size1_summary_b = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size1_summary_c = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> size2_summary_a = prior_sample.summarize<double>("pop_size_pop1");
        SampleSummarizer<double> size2_summary_b = prior_sample.summarize<double>("pop_size_pop2");
        SampleSummarizer<double> size2_summary_c = prior_sample.summarize<double>("pop_size_root_pop1");
        SampleSummarizer<double> size3_summary_a = prior_sample.summarize<double>("pop_size_pop1b");
        SampleSummarizer<double> size3_summary_b = prior_sample.summarize<double>("pop_size_pop2b");
        SampleSummarizer<double> size3_summary_c = prior_sample.summarize<double>("pop_size_root_pop1b");

        REQUIRE(size1_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_b.mean() == Approx(size1_shape * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_b.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.05));

        REQUIRE(size2_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_b.mean() == Approx(size2_shape * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_b.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.05));

        REQUIRE(size3_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.05));

        SampleSummarizer<double> f1_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> f2_summary = prior_sample.summarize<double>("freq_1_pop1");
        SampleSummarizer<double> f3_summary = prior_sample.summarize<double>("freq_1_pop1b");
        REQUIRE(f1_summary.sample_size() == expected_sample_size);
        REQUIRE(f2_summary.sample_size() == expected_sample_size);
        REQUIRE(f3_summary.sample_size() == expected_sample_size);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.05));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.05));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.05));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.05));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.05));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.05));

        SampleSummarizer<double> mult1_summary = prior_sample.summarize<double>("mutation_rate_kya");
        SampleSummarizer<double> mult2_summary = prior_sample.summarize<double>("mutation_rate_pop1");
        SampleSummarizer<double> mult3_summary = prior_sample.summarize<double>("mutation_rate_pop1b");
        REQUIRE(mult1_summary.sample_size() == expected_sample_size);
        REQUIRE(mult2_summary.sample_size() == expected_sample_size);
        REQUIRE(mult3_summary.sample_size() == expected_sample_size);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.05));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.05));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.05));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.05));

        REQUIRE_THROWS(prior_sample.summarize<double>("concentration"));

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        REQUIRE(model_counts.at("000") == 0);
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == 0);

        REQUIRE(model_counts.at("012") == expected_sample_size);
        REQUIRE(nevent_counts.at(3) == expected_sample_size);

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing fixed 000 and fully parameterized", "[SamplingPrior]") {

    SECTION("Testing fixed 000 pairs") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_a = 2.0;
        double f1_b = 1.1;
        double f2_a = 1.0;
        double f2_b = 0.5;
        double f3_a = 1.5;
        double f3_b = 1.8;
        double expected_f1_mean = f1_a / (f1_a + f1_b);
        double expected_f1_variance = (f1_a * f1_b) / ((f1_a + f1_b) * (f1_a + f1_b) * (f1_a + f1_b + 1.0));
        double expected_f2_mean = f2_a / (f2_a + f2_b);
        double expected_f2_variance = (f2_a * f2_b) / ((f2_a + f2_b) * (f2_a + f2_b) * (f2_a + f2_b + 1.0));
        double expected_f3_mean = f3_a / (f3_a + f3_b);
        double expected_f3_variance = (f3_a * f3_b) / ((f3_a + f3_b) * (f3_a + f3_b) * (f3_a + f3_b + 1.0));
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
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
        os << "    fixed: [0, 0, 0]\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 500000\n";
        os << "    sample_frequency: 100\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 0.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_a << "\n";
        os << "                    beta: " << f1_b << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_a << "\n";
        os << "                    beta: " << f2_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_a << "\n";
        os << "                    beta: " << f3_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "58961543";
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

        unsigned int expected_sample_size = 5001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");
        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.05));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.05));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.05));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.05));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.05));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.05));

        SampleSummarizer<double> size1_summary_a = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size1_summary_b = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size1_summary_c = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> size2_summary_a = prior_sample.summarize<double>("pop_size_pop1");
        SampleSummarizer<double> size2_summary_b = prior_sample.summarize<double>("pop_size_pop2");
        SampleSummarizer<double> size2_summary_c = prior_sample.summarize<double>("pop_size_root_pop1");
        SampleSummarizer<double> size3_summary_a = prior_sample.summarize<double>("pop_size_pop1b");
        SampleSummarizer<double> size3_summary_b = prior_sample.summarize<double>("pop_size_pop2b");
        SampleSummarizer<double> size3_summary_c = prior_sample.summarize<double>("pop_size_root_pop1b");

        REQUIRE(size1_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_b.mean() == Approx(size1_shape * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_b.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.05));

        REQUIRE(size2_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_b.mean() == Approx(size2_shape * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_b.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.05));

        REQUIRE(size3_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.05));

        SampleSummarizer<double> f1_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> f2_summary = prior_sample.summarize<double>("freq_1_pop1");
        SampleSummarizer<double> f3_summary = prior_sample.summarize<double>("freq_1_pop1b");
        REQUIRE(f1_summary.sample_size() == expected_sample_size);
        REQUIRE(f2_summary.sample_size() == expected_sample_size);
        REQUIRE(f3_summary.sample_size() == expected_sample_size);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.05));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.05));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.05));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.05));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.05));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.05));

        SampleSummarizer<double> mult1_summary = prior_sample.summarize<double>("mutation_rate_kya");
        SampleSummarizer<double> mult2_summary = prior_sample.summarize<double>("mutation_rate_pop1");
        SampleSummarizer<double> mult3_summary = prior_sample.summarize<double>("mutation_rate_pop1b");
        REQUIRE(mult1_summary.sample_size() == expected_sample_size);
        REQUIRE(mult2_summary.sample_size() == expected_sample_size);
        REQUIRE(mult3_summary.sample_size() == expected_sample_size);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.05));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.05));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.05));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.05));

        REQUIRE_THROWS(prior_sample.summarize<double>("concentration"));

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == 0);
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == 0);

        REQUIRE(model_counts.at("000") == expected_sample_size);
        REQUIRE(nevent_counts.at(1) == expected_sample_size);

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing fixed 000 and fully parameterized and CompositeHeightSizeRateMixer",
        "[SamplingPrior]") {

    SECTION("Testing fixed 000 pairs") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_a = 2.0;
        double f1_b = 1.1;
        double f2_a = 1.0;
        double f2_b = 0.5;
        double f3_a = 1.5;
        double f3_b = 1.8;
        double expected_f1_mean = f1_a / (f1_a + f1_b);
        double expected_f1_variance = (f1_a * f1_b) / ((f1_a + f1_b) * (f1_a + f1_b) * (f1_a + f1_b + 1.0));
        double expected_f2_mean = f2_a / (f2_a + f2_b);
        double expected_f2_variance = (f2_a * f2_b) / ((f2_a + f2_b) * (f2_a + f2_b) * (f2_a + f2_b + 1.0));
        double expected_f3_mean = f3_a / (f3_a + f3_b);
        double expected_f3_variance = (f3_a * f3_b) / ((f3_a + f3_b) * (f3_a + f3_b) * (f3_a + f3_b + 1.0));
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
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
        os << "    fixed: [0, 0, 0]\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 500000\n";
        os << "    sample_frequency: 100\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        CompositeHeightSizeRateMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonHeightScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ComparisonMutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        ChildPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    use_empirical_starting_value_for_freq_1: false\n";
        os << "    equal_population_sizes: false\n";
        os << "    equal_state_frequencies: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_a << "\n";
        os << "                    beta: " << f1_b << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_a << "\n";
        os << "                    beta: " << f2_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_a << "\n";
        os << "                    beta: " << f3_b << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "58961543";
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

        unsigned int expected_sample_size = 5001;

        SampleSummarizer<double> height_summary1 = prior_sample.summarize<double>("root_height_kya");
        SampleSummarizer<double> height_summary2 = prior_sample.summarize<double>("root_height_pop1");
        SampleSummarizer<double> height_summary3 = prior_sample.summarize<double>("root_height_pop1b");
        REQUIRE(height_summary1.sample_size() == expected_sample_size);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.05));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.05));
        REQUIRE(height_summary2.sample_size() == expected_sample_size);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.05));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.05));
        REQUIRE(height_summary3.sample_size() == expected_sample_size);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.05));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.05));

        SampleSummarizer<double> size1_summary_a = prior_sample.summarize<double>("pop_size_kya");
        SampleSummarizer<double> size1_summary_b = prior_sample.summarize<double>("pop_size_fas");
        SampleSummarizer<double> size1_summary_c = prior_sample.summarize<double>("pop_size_root_kya");
        SampleSummarizer<double> size2_summary_a = prior_sample.summarize<double>("pop_size_pop1");
        SampleSummarizer<double> size2_summary_b = prior_sample.summarize<double>("pop_size_pop2");
        SampleSummarizer<double> size2_summary_c = prior_sample.summarize<double>("pop_size_root_pop1");
        SampleSummarizer<double> size3_summary_a = prior_sample.summarize<double>("pop_size_pop1b");
        SampleSummarizer<double> size3_summary_b = prior_sample.summarize<double>("pop_size_pop2b");
        SampleSummarizer<double> size3_summary_c = prior_sample.summarize<double>("pop_size_root_pop1b");

        REQUIRE(size1_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_b.mean() == Approx(size1_shape * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_b.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.05));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.05));

        REQUIRE(size2_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_b.mean() == Approx(size2_shape * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_b.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.05));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.05));

        REQUIRE(size3_summary_a.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_b.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_c.sample_size() == expected_sample_size);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.05));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.05));

        SampleSummarizer<double> f1_summary = prior_sample.summarize<double>("freq_1_kya");
        SampleSummarizer<double> f2_summary = prior_sample.summarize<double>("freq_1_pop1");
        SampleSummarizer<double> f3_summary = prior_sample.summarize<double>("freq_1_pop1b");
        REQUIRE(f1_summary.sample_size() == expected_sample_size);
        REQUIRE(f2_summary.sample_size() == expected_sample_size);
        REQUIRE(f3_summary.sample_size() == expected_sample_size);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.05));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.05));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.05));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.05));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.05));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.05));

        SampleSummarizer<double> mult1_summary = prior_sample.summarize<double>("mutation_rate_kya");
        SampleSummarizer<double> mult2_summary = prior_sample.summarize<double>("mutation_rate_pop1");
        SampleSummarizer<double> mult3_summary = prior_sample.summarize<double>("mutation_rate_pop1b");
        REQUIRE(mult1_summary.sample_size() == expected_sample_size);
        REQUIRE(mult2_summary.sample_size() == expected_sample_size);
        REQUIRE(mult3_summary.sample_size() == expected_sample_size);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.05));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.05));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.05));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.05));

        REQUIRE_THROWS(prior_sample.summarize<double>("concentration"));

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
        REQUIRE(total == expected_sample_size);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == expected_sample_size);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == 0);
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == 0);

        REQUIRE(model_counts.at("000") == expected_sample_size);
        REQUIRE(nevent_counts.at(1) == expected_sample_size);

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);

        delete[] cfg_path;
    }
}
