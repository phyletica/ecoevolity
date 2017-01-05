#include "catch.hpp"
#include "ecoevolity/collection.hpp"
#include "ecoevolity/stats_util.hpp"

#include "ecoevolity/rng.hpp"
#include "ecoevolity/spreadsheet.hpp"

RandomNumberGenerator _COLLECTION_RNG = RandomNumberGenerator();

TEST_CASE("Testing draw_heights_from_prior for DPP with 3 pairs and alpha 1.0",
        "[ComparisonPopulationTreeCollection]") {

    SECTION("Testing draw_heights_from_prior for DPP with 3 pairs and alpha 1.0") {
        std::string cfg_path = "data/dummy.yml";

        double height_shape = 10.0;
        double height_scale = 0.1;
        double concentration = 1.0;
        std::string auto_optimize = "true";

        std::stringstream cfg;
        cfg << "event_time_prior:\n";
        cfg << "    gamma_distribution:\n";
        cfg << "        shape: " << height_shape << "\n";
        cfg << "        scale: " << height_scale << "\n";
        cfg << "event_model_prior:\n";
        cfg << "    dirichlet_process:\n";
        cfg << "        parameters:\n";
        cfg << "            concentration:\n";
        cfg << "                value: " << concentration << "\n";
        cfg << "                estimate: false\n";
        cfg << "mcmc_settings:\n";
        cfg << "    chain_length: 100000\n";
        cfg << "    sample_frequency: 10\n";
        cfg << "operator_settings:\n";
        cfg << "    auto_optimize: " << auto_optimize << "\n";
        cfg << "    auto_optimize_delay: 10000\n";
        cfg << "    operators:\n";
        cfg << "        ModelOperator:\n";
        cfg << "            number_of_auxiliary_categories: 5\n";
        cfg << "            weight: 1.0\n";
        cfg << "        ConcentrationScaler:\n";
        cfg << "            scale: 0.2\n";
        cfg << "            weight: 0.0\n";
        cfg << "        ComparisonHeightScaler:\n";
        cfg << "            scale: 0.3\n";
        cfg << "            weight: 1.0\n";
        cfg << "        ComparisonMutationRateScaler:\n";
        cfg << "            scale: 0.5\n";
        cfg << "            weight: 0.0\n";
        cfg << "        RootPopulationSizeScaler:\n";
        cfg << "            scale: 0.2\n";
        cfg << "            weight: 0.0\n";
        cfg << "        ChildPopulationSizeScaler:\n";
        cfg << "            scale: 0.2\n";
        cfg << "            weight: 0.0\n";
        cfg << "        FreqMover:\n";
        cfg << "            window: 0.1\n";
        cfg << "            weight: 0.0\n";
        cfg << "global_comparison_settings:\n";
        cfg << "    genotypes_are_diploid: true\n";
        cfg << "    markers_are_dominant: false\n";
        cfg << "    population_name_delimiter: \" \"\n";
        cfg << "    population_name_is_prefix: true\n";
        cfg << "    constant_sites_removed: true\n";
        cfg << "    use_empirical_starting_value_for_freq_1: false\n";
        cfg << "    equal_population_sizes: true\n";
        cfg << "    equal_state_frequencies: true\n";
        cfg << "    parameters:\n";
        cfg << "        population_size:\n";
        cfg << "            value: 0.005\n";
        cfg << "            estimate: false\n";
        cfg << "        freq_1:\n";
        cfg << "            value: 1.0\n";
        cfg << "            estimate: false\n";
        cfg << "        mutation_rate:\n";
        cfg << "            value: 1.0\n";
        cfg << "            estimate: false\n";
        cfg << "comparisons:\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129.nex\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129-altname1.nex\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129-altname2.nex\n";

        CollectionSettings settings = CollectionSettings(cfg, cfg_path);
        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        ComparisonPopulationTreeCollection collection = ComparisonPopulationTreeCollection(
                settings,
                rng);

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        SampleSummarizer<double> height_summary1;
        SampleSummarizer<double> height_summary2;
        SampleSummarizer<double> height_summary3;

        unsigned int nsamples = 100000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            collection.draw_heights_from_prior(rng);
            height_summary1.add_sample(collection.get_height_of_tree(0));
            height_summary2.add_sample(collection.get_height_of_tree(1));
            height_summary3.add_sample(collection.get_height_of_tree(2));
            unsigned int nevents = collection.get_number_of_events();
            std::vector<unsigned int> h_indices = collection.get_standardized_height_indices();
            std::ostringstream stream;
            for (auto idx: h_indices) {
                stream << idx;
            }
            if (nevent_counts.count(nevents) < 1) {
                nevent_counts[nevents] = 1;
            }
            else {
                ++nevent_counts.at(nevents);
            }
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }

            if (nevents == 1) {
                REQUIRE(h_indices.at(0) == h_indices.at(1));
                REQUIRE(h_indices.at(0) == h_indices.at(2));
                REQUIRE(collection.get_height_of_tree(0) == collection.get_height_of_tree(1));
                REQUIRE(collection.get_height_of_tree(0) == collection.get_height_of_tree(2));
            }
            if (nevents == 3) {
                REQUIRE(h_indices.at(0) != h_indices.at(1));
                REQUIRE(h_indices.at(0) != h_indices.at(2));
                REQUIRE(h_indices.at(1) != h_indices.at(2));
                REQUIRE(collection.get_height_of_tree(0) != collection.get_height_of_tree(1));
                REQUIRE(collection.get_height_of_tree(0) != collection.get_height_of_tree(2));
                REQUIRE(collection.get_height_of_tree(1) != collection.get_height_of_tree(2));
            }
        }
        REQUIRE(height_summary1.sample_size() == nsamples);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary2.sample_size() == nsamples);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary3.sample_size() == nsamples);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_dpp_log_prior_probability(m, concentration))).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing draw_from_prior for DPP with 3 pairs and alpha 1.0",
        "[ComparisonPopulationTreeCollection]") {

    SECTION("Testing draw_from_prior for DPP with 3 pairs and alpha 1.0") {
        std::string cfg_path = "data/dummy.yml";

        double height_shape = 10.0;
        double height_scale = 0.1;
        double concentration = 1.0;

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

        std::stringstream cfg;
        cfg << "event_time_prior:\n";
        cfg << "    gamma_distribution:\n";
        cfg << "        shape: " << height_shape << "\n";
        cfg << "        scale: " << height_scale << "\n";
        cfg << "event_model_prior:\n";
        cfg << "    dirichlet_process:\n";
        cfg << "        parameters:\n";
        cfg << "            concentration:\n";
        cfg << "                value: " << concentration << "\n";
        cfg << "                estimate: false\n";
        cfg << "mcmc_settings:\n";
        cfg << "    chain_length: 100000\n";
        cfg << "    sample_frequency: 10\n";
        cfg << "operator_settings:\n";
        cfg << "    auto_optimize: " << auto_optimize << "\n";
        cfg << "    auto_optimize_delay: 10000\n";
        cfg << "    operators:\n";
        cfg << "        ModelOperator:\n";
        cfg << "            number_of_auxiliary_categories: 5\n";
        cfg << "            weight: 1.0\n";
        cfg << "        ConcentrationScaler:\n";
        cfg << "            scale: 0.2\n";
        cfg << "            weight: 0.0\n";
        cfg << "        ComparisonHeightScaler:\n";
        cfg << "            scale: 0.3\n";
        cfg << "            weight: 1.0\n";
        cfg << "        ComparisonMutationRateScaler:\n";
        cfg << "            scale: 0.5\n";
        cfg << "            weight: 0.0\n";
        cfg << "        RootPopulationSizeScaler:\n";
        cfg << "            scale: 0.2\n";
        cfg << "            weight: 0.0\n";
        cfg << "        ChildPopulationSizeScaler:\n";
        cfg << "            scale: 0.2\n";
        cfg << "            weight: 0.0\n";
        cfg << "        FreqMover:\n";
        cfg << "            window: 0.1\n";
        cfg << "            weight: 0.0\n";
        cfg << "global_comparison_settings:\n";
        cfg << "    genotypes_are_diploid: true\n";
        cfg << "    markers_are_dominant: false\n";
        cfg << "    population_name_delimiter: \" \"\n";
        cfg << "    population_name_is_prefix: true\n";
        cfg << "    constant_sites_removed: true\n";
        cfg << "    use_empirical_starting_value_for_freq_1: false\n";
        cfg << "    equal_population_sizes: false\n";
        cfg << "    equal_state_frequencies: false\n";
        cfg << "comparisons:\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129.nex\n";
        cfg << "    parameters:\n";
        cfg << "        population_size:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << size1_shape << "\n";
        cfg << "                    scale: " << size1_scale << "\n";
        cfg << "        freq_1:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                beta_distribution:\n";
        cfg << "                    alpha: " << f1_a << "\n";
        cfg << "                    beta: " << f1_b << "\n";
        cfg << "        mutation_rate:\n";
        cfg << "            value: 1.0\n";
        cfg << "            estimate: false\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129-altname1.nex\n";
        cfg << "    parameters:\n";
        cfg << "        population_size:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << size2_shape << "\n";
        cfg << "                    scale: " << size2_scale << "\n";
        cfg << "        freq_1:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                beta_distribution:\n";
        cfg << "                    alpha: " << f2_a << "\n";
        cfg << "                    beta: " << f2_b << "\n";
        cfg << "        mutation_rate:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << mult2_shape << "\n";
        cfg << "                    scale: " << mult2_scale << "\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129-altname2.nex\n";
        cfg << "    parameters:\n";
        cfg << "        population_size:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << size3_shape << "\n";
        cfg << "                    scale: " << size3_scale << "\n";
        cfg << "        freq_1:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                beta_distribution:\n";
        cfg << "                    alpha: " << f3_a << "\n";
        cfg << "                    beta: " << f3_b << "\n";
        cfg << "        mutation_rate:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << mult3_shape << "\n";
        cfg << "                    scale: " << mult3_scale << "\n";

        CollectionSettings settings = CollectionSettings(cfg, cfg_path);
        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        ComparisonPopulationTreeCollection collection = ComparisonPopulationTreeCollection(
                settings,
                rng);

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        SampleSummarizer<double> height_summary1;
        SampleSummarizer<double> height_summary2;
        SampleSummarizer<double> height_summary3;
        SampleSummarizer<double> size1_summary_a;
        SampleSummarizer<double> size1_summary_b;
        SampleSummarizer<double> size1_summary_c;
        SampleSummarizer<double> size2_summary_a;
        SampleSummarizer<double> size2_summary_b;
        SampleSummarizer<double> size2_summary_c;
        SampleSummarizer<double> size3_summary_a;
        SampleSummarizer<double> size3_summary_b;
        SampleSummarizer<double> size3_summary_c;
        SampleSummarizer<double> f1_summary;
        SampleSummarizer<double> f2_summary;
        SampleSummarizer<double> f3_summary;
        SampleSummarizer<double> mult1_summary;
        SampleSummarizer<double> mult2_summary;
        SampleSummarizer<double> mult3_summary;
        SampleSummarizer<double> conc_summary;

        unsigned int nsamples = 100000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            collection.draw_from_prior(rng);
            height_summary1.add_sample(collection.get_height_of_tree(0));
            height_summary2.add_sample(collection.get_height_of_tree(1));
            height_summary3.add_sample(collection.get_height_of_tree(2));
            unsigned int nevents = collection.get_number_of_events();
            std::vector<unsigned int> h_indices = collection.get_standardized_height_indices();
            std::ostringstream stream;
            for (auto idx: h_indices) {
                stream << idx;
            }
            if (nevent_counts.count(nevents) < 1) {
                nevent_counts[nevents] = 1;
            }
            else {
                ++nevent_counts.at(nevents);
            }
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }

            if (nevents == 1) {
                REQUIRE(h_indices.at(0) == h_indices.at(1));
                REQUIRE(h_indices.at(0) == h_indices.at(2));
                REQUIRE(collection.get_height_of_tree(0) == collection.get_height_of_tree(1));
                REQUIRE(collection.get_height_of_tree(0) == collection.get_height_of_tree(2));
            }
            if (nevents == 3) {
                REQUIRE(h_indices.at(0) != h_indices.at(1));
                REQUIRE(h_indices.at(0) != h_indices.at(2));
                REQUIRE(h_indices.at(1) != h_indices.at(2));
                REQUIRE(collection.get_height_of_tree(0) != collection.get_height_of_tree(1));
                REQUIRE(collection.get_height_of_tree(0) != collection.get_height_of_tree(2));
                REQUIRE(collection.get_height_of_tree(1) != collection.get_height_of_tree(2));
            }

            size1_summary_a.add_sample(collection.get_tree(0).get_root_population_size());
            size1_summary_b.add_sample(collection.get_tree(0).get_child_population_size(0));
            size1_summary_c.add_sample(collection.get_tree(0).get_child_population_size(1));
            size2_summary_a.add_sample(collection.get_tree(1).get_root_population_size());
            size2_summary_b.add_sample(collection.get_tree(1).get_child_population_size(0));
            size2_summary_c.add_sample(collection.get_tree(1).get_child_population_size(1));
            size3_summary_a.add_sample(collection.get_tree(2).get_root_population_size());
            size3_summary_b.add_sample(collection.get_tree(2).get_child_population_size(0));
            size3_summary_c.add_sample(collection.get_tree(2).get_child_population_size(1));
            f1_summary.add_sample(collection.get_tree(0).get_freq_1());
            f2_summary.add_sample(collection.get_tree(1).get_freq_1());
            f3_summary.add_sample(collection.get_tree(2).get_freq_1());
            mult1_summary.add_sample(collection.get_tree(0).get_mutation_rate());
            mult2_summary.add_sample(collection.get_tree(1).get_mutation_rate());
            mult3_summary.add_sample(collection.get_tree(2).get_mutation_rate());
            conc_summary.add_sample(collection.get_concentration());
        }

        REQUIRE(height_summary1.sample_size() == nsamples);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary2.sample_size() == nsamples);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary3.sample_size() == nsamples);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_dpp_log_prior_probability(m, concentration))).epsilon(0.01));
        }

        REQUIRE(size1_summary_a.sample_size() == nsamples);
        REQUIRE(size1_summary_b.sample_size() == nsamples);
        REQUIRE(size1_summary_c.sample_size() == nsamples);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));

        REQUIRE(size2_summary_a.sample_size() == nsamples);
        REQUIRE(size2_summary_b.sample_size() == nsamples);
        REQUIRE(size2_summary_c.sample_size() == nsamples);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));

        REQUIRE(size3_summary_a.sample_size() == nsamples);
        REQUIRE(size3_summary_b.sample_size() == nsamples);
        REQUIRE(size3_summary_c.sample_size() == nsamples);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));

        REQUIRE(f1_summary.sample_size() == nsamples);
        REQUIRE(f2_summary.sample_size() == nsamples);
        REQUIRE(f3_summary.sample_size() == nsamples);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.01));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.01));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.01));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.01));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.01));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.01));

        REQUIRE(mult1_summary.sample_size() == nsamples);
        REQUIRE(mult2_summary.sample_size() == nsamples);
        REQUIRE(mult3_summary.sample_size() == nsamples);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.01));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.01));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.01));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.01));
        REQUIRE(conc_summary.sample_size() == nsamples);
        REQUIRE(conc_summary.mean() == 1.0);
        REQUIRE(conc_summary.variance() == 0.0);
    }
}

TEST_CASE("Testing draw_from_prior for DPP with variable alpha",
        "[ComparisonPopulationTreeCollection]") {

    SECTION("Testing draw_from_prior for DPP with variable alpha") {
        std::string cfg_path = "data/dummy.yml";

        double height_shape = 10.0;
        double height_scale = 0.1;
        double concentration = 1.0;

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

        std::stringstream cfg;
        cfg << "event_time_prior:\n";
        cfg << "    gamma_distribution:\n";
        cfg << "        shape: " << height_shape << "\n";
        cfg << "        scale: " << height_scale << "\n";
        cfg << "event_model_prior:\n";
        cfg << "    dirichlet_process:\n";
        cfg << "        parameters:\n";
        cfg << "            concentration:\n";
        cfg << "                estimate: true\n";
        cfg << "                prior:\n";
        cfg << "                    gamma_distribution:\n";
        cfg << "                        shape: " << concentration_shape << "\n";
        cfg << "                        scale: " << concentration_scale << "\n";
        cfg << "mcmc_settings:\n";
        cfg << "    chain_length: 100000\n";
        cfg << "    sample_frequency: 10\n";
        cfg << "operator_settings:\n";
        cfg << "    auto_optimize: " << auto_optimize << "\n";
        cfg << "    auto_optimize_delay: 10000\n";
        cfg << "    operators:\n";
        cfg << "        ModelOperator:\n";
        cfg << "            number_of_auxiliary_categories: 5\n";
        cfg << "            weight: 1.0\n";
        cfg << "        ConcentrationScaler:\n";
        cfg << "            scale: 0.2\n";
        cfg << "            weight: 0.0\n";
        cfg << "        ComparisonHeightScaler:\n";
        cfg << "            scale: 0.3\n";
        cfg << "            weight: 1.0\n";
        cfg << "        ComparisonMutationRateScaler:\n";
        cfg << "            scale: 0.5\n";
        cfg << "            weight: 0.0\n";
        cfg << "        RootPopulationSizeScaler:\n";
        cfg << "            scale: 0.2\n";
        cfg << "            weight: 0.0\n";
        cfg << "        ChildPopulationSizeScaler:\n";
        cfg << "            scale: 0.2\n";
        cfg << "            weight: 0.0\n";
        cfg << "        FreqMover:\n";
        cfg << "            window: 0.1\n";
        cfg << "            weight: 0.0\n";
        cfg << "global_comparison_settings:\n";
        cfg << "    genotypes_are_diploid: true\n";
        cfg << "    markers_are_dominant: false\n";
        cfg << "    population_name_delimiter: \" \"\n";
        cfg << "    population_name_is_prefix: true\n";
        cfg << "    constant_sites_removed: true\n";
        cfg << "    use_empirical_starting_value_for_freq_1: false\n";
        cfg << "    equal_population_sizes: false\n";
        cfg << "    equal_state_frequencies: false\n";
        cfg << "comparisons:\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129.nex\n";
        cfg << "    parameters:\n";
        cfg << "        population_size:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << size1_shape << "\n";
        cfg << "                    scale: " << size1_scale << "\n";
        cfg << "        freq_1:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                beta_distribution:\n";
        cfg << "                    alpha: " << f1_a << "\n";
        cfg << "                    beta: " << f1_b << "\n";
        cfg << "        mutation_rate:\n";
        cfg << "            value: 1.0\n";
        cfg << "            estimate: false\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129-altname1.nex\n";
        cfg << "    parameters:\n";
        cfg << "        population_size:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << size2_shape << "\n";
        cfg << "                    scale: " << size2_scale << "\n";
        cfg << "        freq_1:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                beta_distribution:\n";
        cfg << "                    alpha: " << f2_a << "\n";
        cfg << "                    beta: " << f2_b << "\n";
        cfg << "        mutation_rate:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << mult2_shape << "\n";
        cfg << "                    scale: " << mult2_scale << "\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129-altname2.nex\n";
        cfg << "    parameters:\n";
        cfg << "        population_size:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << size3_shape << "\n";
        cfg << "                    scale: " << size3_scale << "\n";
        cfg << "        freq_1:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                beta_distribution:\n";
        cfg << "                    alpha: " << f3_a << "\n";
        cfg << "                    beta: " << f3_b << "\n";
        cfg << "        mutation_rate:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << mult3_shape << "\n";
        cfg << "                    scale: " << mult3_scale << "\n";

        CollectionSettings settings = CollectionSettings(cfg, cfg_path);
        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        ComparisonPopulationTreeCollection collection = ComparisonPopulationTreeCollection(
                settings,
                rng);

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        SampleSummarizer<double> height_summary1;
        SampleSummarizer<double> height_summary2;
        SampleSummarizer<double> height_summary3;
        SampleSummarizer<double> size1_summary_a;
        SampleSummarizer<double> size1_summary_b;
        SampleSummarizer<double> size1_summary_c;
        SampleSummarizer<double> size2_summary_a;
        SampleSummarizer<double> size2_summary_b;
        SampleSummarizer<double> size2_summary_c;
        SampleSummarizer<double> size3_summary_a;
        SampleSummarizer<double> size3_summary_b;
        SampleSummarizer<double> size3_summary_c;
        SampleSummarizer<double> f1_summary;
        SampleSummarizer<double> f2_summary;
        SampleSummarizer<double> f3_summary;
        SampleSummarizer<double> mult1_summary;
        SampleSummarizer<double> mult2_summary;
        SampleSummarizer<double> mult3_summary;
        SampleSummarizer<double> conc_summary;

        unsigned int nsamples = 100000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            collection.draw_from_prior(rng);
            height_summary1.add_sample(collection.get_height_of_tree(0));
            height_summary2.add_sample(collection.get_height_of_tree(1));
            height_summary3.add_sample(collection.get_height_of_tree(2));
            unsigned int nevents = collection.get_number_of_events();
            std::vector<unsigned int> h_indices = collection.get_standardized_height_indices();
            std::ostringstream stream;
            for (auto idx: h_indices) {
                stream << idx;
            }
            if (nevent_counts.count(nevents) < 1) {
                nevent_counts[nevents] = 1;
            }
            else {
                ++nevent_counts.at(nevents);
            }
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }

            if (nevents == 1) {
                REQUIRE(h_indices.at(0) == h_indices.at(1));
                REQUIRE(h_indices.at(0) == h_indices.at(2));
                REQUIRE(collection.get_height_of_tree(0) == collection.get_height_of_tree(1));
                REQUIRE(collection.get_height_of_tree(0) == collection.get_height_of_tree(2));
            }
            if (nevents == 3) {
                REQUIRE(h_indices.at(0) != h_indices.at(1));
                REQUIRE(h_indices.at(0) != h_indices.at(2));
                REQUIRE(h_indices.at(1) != h_indices.at(2));
                REQUIRE(collection.get_height_of_tree(0) != collection.get_height_of_tree(1));
                REQUIRE(collection.get_height_of_tree(0) != collection.get_height_of_tree(2));
                REQUIRE(collection.get_height_of_tree(1) != collection.get_height_of_tree(2));
            }

            size1_summary_a.add_sample(collection.get_tree(0).get_root_population_size());
            size1_summary_b.add_sample(collection.get_tree(0).get_child_population_size(0));
            size1_summary_c.add_sample(collection.get_tree(0).get_child_population_size(1));
            size2_summary_a.add_sample(collection.get_tree(1).get_root_population_size());
            size2_summary_b.add_sample(collection.get_tree(1).get_child_population_size(0));
            size2_summary_c.add_sample(collection.get_tree(1).get_child_population_size(1));
            size3_summary_a.add_sample(collection.get_tree(2).get_root_population_size());
            size3_summary_b.add_sample(collection.get_tree(2).get_child_population_size(0));
            size3_summary_c.add_sample(collection.get_tree(2).get_child_population_size(1));
            f1_summary.add_sample(collection.get_tree(0).get_freq_1());
            f2_summary.add_sample(collection.get_tree(1).get_freq_1());
            f3_summary.add_sample(collection.get_tree(2).get_freq_1());
            mult1_summary.add_sample(collection.get_tree(0).get_mutation_rate());
            mult2_summary.add_sample(collection.get_tree(1).get_mutation_rate());
            mult3_summary.add_sample(collection.get_tree(2).get_mutation_rate());
            conc_summary.add_sample(collection.get_concentration());
        }

        REQUIRE(height_summary1.sample_size() == nsamples);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary2.sample_size() == nsamples);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));
        REQUIRE(height_summary3.sample_size() == nsamples);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.001));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.001));

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        REQUIRE(size1_summary_a.sample_size() == nsamples);
        REQUIRE(size1_summary_b.sample_size() == nsamples);
        REQUIRE(size1_summary_c.sample_size() == nsamples);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));

        REQUIRE(size2_summary_a.sample_size() == nsamples);
        REQUIRE(size2_summary_b.sample_size() == nsamples);
        REQUIRE(size2_summary_c.sample_size() == nsamples);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));

        REQUIRE(size3_summary_a.sample_size() == nsamples);
        REQUIRE(size3_summary_b.sample_size() == nsamples);
        REQUIRE(size3_summary_c.sample_size() == nsamples);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));

        REQUIRE(f1_summary.sample_size() == nsamples);
        REQUIRE(f2_summary.sample_size() == nsamples);
        REQUIRE(f3_summary.sample_size() == nsamples);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.01));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.01));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.01));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.01));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.01));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.01));

        REQUIRE(mult1_summary.sample_size() == nsamples);
        REQUIRE(mult2_summary.sample_size() == nsamples);
        REQUIRE(mult3_summary.sample_size() == nsamples);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.01));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.01));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.01));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.01));
        REQUIRE(conc_summary.sample_size() == nsamples);
        REQUIRE(conc_summary.mean() == Approx(concentration_shape * concentration_scale).epsilon(0.01));
        REQUIRE(conc_summary.variance() == Approx(concentration_shape * concentration_scale * concentration_scale).epsilon(0.01));
    }
}

TEST_CASE("Testing draw_from_prior for reversible jump with 3 pairs",
        "[ComparisonPopulationTreeCollection]") {

    SECTION("Testing draw_from_prior for reversible jump with 3 pairs") {
        std::string cfg_path = "data/dummy.yml";

        double height_shape = 10.0;
        double height_scale = 0.1;

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

        std::stringstream cfg;
        cfg << "event_time_prior:\n";
        cfg << "    gamma_distribution:\n";
        cfg << "        shape: " << height_shape << "\n";
        cfg << "        scale: " << height_scale << "\n";
        cfg << "event_model_prior:\n";
        cfg << "    uniform:\n";
        cfg << "mcmc_settings:\n";
        cfg << "    chain_length: 100000\n";
        cfg << "    sample_frequency: 10\n";
        cfg << "operator_settings:\n";
        cfg << "    auto_optimize: " << auto_optimize << "\n";
        cfg << "    auto_optimize_delay: 10000\n";
        cfg << "    operators:\n";
        cfg << "        ModelOperator:\n";
        cfg << "            weight: 1.0\n";
        cfg << "        ComparisonHeightScaler:\n";
        cfg << "            scale: 0.3\n";
        cfg << "            weight: 1.0\n";
        cfg << "        ComparisonMutationRateScaler:\n";
        cfg << "            scale: 0.5\n";
        cfg << "            weight: 0.0\n";
        cfg << "        RootPopulationSizeScaler:\n";
        cfg << "            scale: 0.2\n";
        cfg << "            weight: 0.0\n";
        cfg << "        ChildPopulationSizeScaler:\n";
        cfg << "            scale: 0.2\n";
        cfg << "            weight: 0.0\n";
        cfg << "        FreqMover:\n";
        cfg << "            window: 0.1\n";
        cfg << "            weight: 0.0\n";
        cfg << "global_comparison_settings:\n";
        cfg << "    genotypes_are_diploid: true\n";
        cfg << "    markers_are_dominant: false\n";
        cfg << "    population_name_delimiter: \" \"\n";
        cfg << "    population_name_is_prefix: true\n";
        cfg << "    constant_sites_removed: true\n";
        cfg << "    use_empirical_starting_value_for_freq_1: false\n";
        cfg << "    equal_population_sizes: false\n";
        cfg << "    equal_state_frequencies: false\n";
        cfg << "comparisons:\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129.nex\n";
        cfg << "    parameters:\n";
        cfg << "        population_size:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << size1_shape << "\n";
        cfg << "                    scale: " << size1_scale << "\n";
        cfg << "        freq_1:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                beta_distribution:\n";
        cfg << "                    alpha: " << f1_a << "\n";
        cfg << "                    beta: " << f1_b << "\n";
        cfg << "        mutation_rate:\n";
        cfg << "            value: 1.0\n";
        cfg << "            estimate: false\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129-altname1.nex\n";
        cfg << "    parameters:\n";
        cfg << "        population_size:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << size2_shape << "\n";
        cfg << "                    scale: " << size2_scale << "\n";
        cfg << "        freq_1:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                beta_distribution:\n";
        cfg << "                    alpha: " << f2_a << "\n";
        cfg << "                    beta: " << f2_b << "\n";
        cfg << "        mutation_rate:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << mult2_shape << "\n";
        cfg << "                    scale: " << mult2_scale << "\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129-altname2.nex\n";
        cfg << "    parameters:\n";
        cfg << "        population_size:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << size3_shape << "\n";
        cfg << "                    scale: " << size3_scale << "\n";
        cfg << "        freq_1:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                beta_distribution:\n";
        cfg << "                    alpha: " << f3_a << "\n";
        cfg << "                    beta: " << f3_b << "\n";
        cfg << "        mutation_rate:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << mult3_shape << "\n";
        cfg << "                    scale: " << mult3_scale << "\n";

        CollectionSettings settings = CollectionSettings(cfg, cfg_path);
        RandomNumberGenerator rng = RandomNumberGenerator(98714654);
        ComparisonPopulationTreeCollection collection = ComparisonPopulationTreeCollection(
                settings,
                rng);

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        SampleSummarizer<double> height_summary1;
        SampleSummarizer<double> height_summary2;
        SampleSummarizer<double> height_summary3;
        SampleSummarizer<double> size1_summary_a;
        SampleSummarizer<double> size1_summary_b;
        SampleSummarizer<double> size1_summary_c;
        SampleSummarizer<double> size2_summary_a;
        SampleSummarizer<double> size2_summary_b;
        SampleSummarizer<double> size2_summary_c;
        SampleSummarizer<double> size3_summary_a;
        SampleSummarizer<double> size3_summary_b;
        SampleSummarizer<double> size3_summary_c;
        SampleSummarizer<double> f1_summary;
        SampleSummarizer<double> f2_summary;
        SampleSummarizer<double> f3_summary;
        SampleSummarizer<double> mult1_summary;
        SampleSummarizer<double> mult2_summary;
        SampleSummarizer<double> mult3_summary;

        unsigned int nsamples = 10000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            collection.draw_from_prior(rng);
            height_summary1.add_sample(collection.get_height_of_tree(0));
            height_summary2.add_sample(collection.get_height_of_tree(1));
            height_summary3.add_sample(collection.get_height_of_tree(2));
            unsigned int nevents = collection.get_number_of_events();
            std::vector<unsigned int> h_indices = collection.get_standardized_height_indices();
            std::ostringstream stream;
            for (auto idx: h_indices) {
                stream << idx;
            }
            if (nevent_counts.count(nevents) < 1) {
                nevent_counts[nevents] = 1;
            }
            else {
                ++nevent_counts.at(nevents);
            }
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }

            if (nevents == 1) {
                REQUIRE(h_indices.at(0) == h_indices.at(1));
                REQUIRE(h_indices.at(0) == h_indices.at(2));
                REQUIRE(collection.get_height_of_tree(0) == collection.get_height_of_tree(1));
                REQUIRE(collection.get_height_of_tree(0) == collection.get_height_of_tree(2));
            }
            if (nevents == 3) {
                REQUIRE(h_indices.at(0) != h_indices.at(1));
                REQUIRE(h_indices.at(0) != h_indices.at(2));
                REQUIRE(h_indices.at(1) != h_indices.at(2));
                REQUIRE(collection.get_height_of_tree(0) != collection.get_height_of_tree(1));
                REQUIRE(collection.get_height_of_tree(0) != collection.get_height_of_tree(2));
                REQUIRE(collection.get_height_of_tree(1) != collection.get_height_of_tree(2));
            }

            size1_summary_a.add_sample(collection.get_tree(0).get_root_population_size());
            size1_summary_b.add_sample(collection.get_tree(0).get_child_population_size(0));
            size1_summary_c.add_sample(collection.get_tree(0).get_child_population_size(1));
            size2_summary_a.add_sample(collection.get_tree(1).get_root_population_size());
            size2_summary_b.add_sample(collection.get_tree(1).get_child_population_size(0));
            size2_summary_c.add_sample(collection.get_tree(1).get_child_population_size(1));
            size3_summary_a.add_sample(collection.get_tree(2).get_root_population_size());
            size3_summary_b.add_sample(collection.get_tree(2).get_child_population_size(0));
            size3_summary_c.add_sample(collection.get_tree(2).get_child_population_size(1));
            f1_summary.add_sample(collection.get_tree(0).get_freq_1());
            f2_summary.add_sample(collection.get_tree(1).get_freq_1());
            f3_summary.add_sample(collection.get_tree(2).get_freq_1());
            mult1_summary.add_sample(collection.get_tree(0).get_mutation_rate());
            mult2_summary.add_sample(collection.get_tree(1).get_mutation_rate());
            mult3_summary.add_sample(collection.get_tree(2).get_mutation_rate());
        }

        REQUIRE(height_summary1.sample_size() == nsamples);
        REQUIRE(height_summary1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.sample_size() == nsamples);
        REQUIRE(height_summary2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.sample_size() == nsamples);
        REQUIRE(height_summary3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            std::cout << m << ": " << model_counts.at(m) / (double)nsamples << "\n";
        }

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(
                    1.0/5.0).epsilon(0.01));
        }

        REQUIRE(size1_summary_a.sample_size() == nsamples);
        REQUIRE(size1_summary_b.sample_size() == nsamples);
        REQUIRE(size1_summary_c.sample_size() == nsamples);
        REQUIRE(size1_summary_a.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_a.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_b.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.mean() == Approx(size1_shape * size1_scale).epsilon(0.01));
        REQUIRE(size1_summary_c.variance() == Approx(size1_shape * size1_scale * size1_scale).epsilon(0.01));

        REQUIRE(size2_summary_a.sample_size() == nsamples);
        REQUIRE(size2_summary_b.sample_size() == nsamples);
        REQUIRE(size2_summary_c.sample_size() == nsamples);
        REQUIRE(size2_summary_a.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_a.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_b.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.mean() == Approx(size2_shape * size2_scale).epsilon(0.01));
        REQUIRE(size2_summary_c.variance() == Approx(size2_shape * size2_scale * size2_scale).epsilon(0.01));

        REQUIRE(size3_summary_a.sample_size() == nsamples);
        REQUIRE(size3_summary_b.sample_size() == nsamples);
        REQUIRE(size3_summary_c.sample_size() == nsamples);
        REQUIRE(size3_summary_a.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_a.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_b.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.mean() == Approx(size3_shape * size3_scale).epsilon(0.01));
        REQUIRE(size3_summary_c.variance() == Approx(size3_shape * size3_scale * size3_scale).epsilon(0.01));

        REQUIRE(f1_summary.sample_size() == nsamples);
        REQUIRE(f2_summary.sample_size() == nsamples);
        REQUIRE(f3_summary.sample_size() == nsamples);
        REQUIRE(f1_summary.mean() ==     Approx(expected_f1_mean).epsilon(0.01));
        REQUIRE(f2_summary.mean() ==     Approx(expected_f2_mean).epsilon(0.01));
        REQUIRE(f3_summary.mean() ==     Approx(expected_f3_mean).epsilon(0.01));
        REQUIRE(f1_summary.variance() == Approx(expected_f1_variance).epsilon(0.01));
        REQUIRE(f2_summary.variance() == Approx(expected_f2_variance).epsilon(0.01));
        REQUIRE(f3_summary.variance() == Approx(expected_f3_variance).epsilon(0.01));

        REQUIRE(mult1_summary.sample_size() == nsamples);
        REQUIRE(mult2_summary.sample_size() == nsamples);
        REQUIRE(mult3_summary.sample_size() == nsamples);
        REQUIRE(mult1_summary.mean() == 1.0);
        REQUIRE(mult1_summary.variance() == 0.0);
        REQUIRE(mult2_summary.mean() == Approx(mult2_shape * mult2_scale).epsilon(0.01));
        REQUIRE(mult2_summary.variance() == Approx(mult2_shape * mult2_scale * mult2_scale).epsilon(0.01));
        REQUIRE(mult3_summary.mean() == Approx(mult3_shape * mult3_scale).epsilon(0.01));
        REQUIRE(mult3_summary.variance() == Approx(mult3_shape * mult3_scale * mult3_scale).epsilon(0.01));
    }
}

TEST_CASE("Testing draw_from_prior logging for reversible jump with 3 pairs",
        "[ComparisonPopulationTreeCollection]") {

    SECTION("Testing draw_from_prior logging for reversible jump with 3 pairs") {
        std::string cfg_path = "data/dummy.yml";

        double height_shape = 10.0;
        double height_scale = 0.1;

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

        std::string tag = _COLLECTION_RNG.random_string(10);
        std::string log_path = "data/tmp-prior-sampling-" + tag + ".log";

        std::stringstream cfg;
        cfg << "event_time_prior:\n";
        cfg << "    gamma_distribution:\n";
        cfg << "        shape: " << height_shape << "\n";
        cfg << "        scale: " << height_scale << "\n";
        cfg << "event_model_prior:\n";
        cfg << "    uniform:\n";
        cfg << "mcmc_settings:\n";
        cfg << "    chain_length: 100000\n";
        cfg << "    sample_frequency: 10\n";
        cfg << "operator_settings:\n";
        cfg << "    auto_optimize: " << auto_optimize << "\n";
        cfg << "    auto_optimize_delay: 10000\n";
        cfg << "    operators:\n";
        cfg << "        ModelOperator:\n";
        cfg << "            weight: 1.0\n";
        cfg << "        ComparisonHeightScaler:\n";
        cfg << "            scale: 0.3\n";
        cfg << "            weight: 1.0\n";
        cfg << "        ComparisonMutationRateScaler:\n";
        cfg << "            scale: 0.5\n";
        cfg << "            weight: 0.0\n";
        cfg << "        RootPopulationSizeScaler:\n";
        cfg << "            scale: 0.2\n";
        cfg << "            weight: 0.0\n";
        cfg << "        ChildPopulationSizeScaler:\n";
        cfg << "            scale: 0.2\n";
        cfg << "            weight: 0.0\n";
        cfg << "        FreqMover:\n";
        cfg << "            window: 0.1\n";
        cfg << "            weight: 0.0\n";
        cfg << "global_comparison_settings:\n";
        cfg << "    genotypes_are_diploid: true\n";
        cfg << "    markers_are_dominant: false\n";
        cfg << "    population_name_delimiter: \" \"\n";
        cfg << "    population_name_is_prefix: true\n";
        cfg << "    constant_sites_removed: true\n";
        cfg << "    use_empirical_starting_value_for_freq_1: false\n";
        cfg << "    equal_population_sizes: false\n";
        cfg << "    equal_state_frequencies: false\n";
        cfg << "comparisons:\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129.nex\n";
        cfg << "    parameters:\n";
        cfg << "        population_size:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << size1_shape << "\n";
        cfg << "                    scale: " << size1_scale << "\n";
        cfg << "        freq_1:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                beta_distribution:\n";
        cfg << "                    alpha: " << f1_a << "\n";
        cfg << "                    beta: " << f1_b << "\n";
        cfg << "        mutation_rate:\n";
        cfg << "            value: 1.0\n";
        cfg << "            estimate: false\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129-altname1.nex\n";
        cfg << "    parameters:\n";
        cfg << "        population_size:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << size2_shape << "\n";
        cfg << "                    scale: " << size2_scale << "\n";
        cfg << "        freq_1:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                beta_distribution:\n";
        cfg << "                    alpha: " << f2_a << "\n";
        cfg << "                    beta: " << f2_b << "\n";
        cfg << "        mutation_rate:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << mult2_shape << "\n";
        cfg << "                    scale: " << mult2_scale << "\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129-altname2.nex\n";
        cfg << "    parameters:\n";
        cfg << "        population_size:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << size3_shape << "\n";
        cfg << "                    scale: " << size3_scale << "\n";
        cfg << "        freq_1:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                beta_distribution:\n";
        cfg << "                    alpha: " << f3_a << "\n";
        cfg << "                    beta: " << f3_b << "\n";
        cfg << "        mutation_rate:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << mult3_shape << "\n";
        cfg << "                    scale: " << mult3_scale << "\n";

        CollectionSettings settings = CollectionSettings(cfg, cfg_path);
        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        ComparisonPopulationTreeCollection collection = ComparisonPopulationTreeCollection(
                settings,
                rng);

        std::ofstream out_stream;
        out_stream.open(log_path);
        collection.write_state_log_header(out_stream);

        unsigned int nsamples = 10000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            collection.draw_from_prior(rng);
            collection.log_state(out_stream, i);
        }
        out_stream.close();
        

        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);

        unsigned int expected_sample_size = nsamples;

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
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)expected_sample_size << "\n";
        }

        for (auto const & kv: model_counts) {
            REQUIRE((kv.second / (double)expected_sample_size) == Approx(1.0/bell_float(3)).epsilon(0.01));
        }
        for (auto const & kv: nevent_counts) {
            REQUIRE((kv.second / (double)expected_sample_size) == Approx(stirling2_float(3, kv.first)/bell_float(3)).epsilon(0.005));
        }

        // Make sure the rest of the prior sample is as expected
        SampleSummarizer<double> lnl_summary = prior_sample.summarize<double>("ln_likelihood");
        REQUIRE(lnl_summary.mean() == 0.0);
        REQUIRE(lnl_summary.variance() == 0.0);
    }
}

TEST_CASE("Testing draw_from_prior logging for DPP with variable alpha",
        "[ComparisonPopulationTreeCollection]") {

    SECTION("Testing draw_from_prior logging for DPP with variable alpha") {
        std::string cfg_path = "data/dummy.yml";

        double height_shape = 10.0;
        double height_scale = 0.1;
        double concentration = 1.0;

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

        std::string tag = _COLLECTION_RNG.random_string(10);
        std::string log_path = "data/tmp-prior-sampling-" + tag + ".log";

        std::stringstream cfg;
        cfg << "event_time_prior:\n";
        cfg << "    gamma_distribution:\n";
        cfg << "        shape: " << height_shape << "\n";
        cfg << "        scale: " << height_scale << "\n";
        cfg << "event_model_prior:\n";
        cfg << "    dirichlet_process:\n";
        cfg << "        parameters:\n";
        cfg << "            concentration:\n";
        cfg << "                estimate: true\n";
        cfg << "                prior:\n";
        cfg << "                    gamma_distribution:\n";
        cfg << "                        shape: " << concentration_shape << "\n";
        cfg << "                        scale: " << concentration_scale << "\n";
        cfg << "mcmc_settings:\n";
        cfg << "    chain_length: 100000\n";
        cfg << "    sample_frequency: 10\n";
        cfg << "operator_settings:\n";
        cfg << "    auto_optimize: " << auto_optimize << "\n";
        cfg << "    auto_optimize_delay: 10000\n";
        cfg << "    operators:\n";
        cfg << "        ModelOperator:\n";
        cfg << "            number_of_auxiliary_categories: 5\n";
        cfg << "            weight: 1.0\n";
        cfg << "        ConcentrationScaler:\n";
        cfg << "            scale: 0.2\n";
        cfg << "            weight: 0.0\n";
        cfg << "        ComparisonHeightScaler:\n";
        cfg << "            scale: 0.3\n";
        cfg << "            weight: 1.0\n";
        cfg << "        ComparisonMutationRateScaler:\n";
        cfg << "            scale: 0.5\n";
        cfg << "            weight: 0.0\n";
        cfg << "        RootPopulationSizeScaler:\n";
        cfg << "            scale: 0.2\n";
        cfg << "            weight: 0.0\n";
        cfg << "        ChildPopulationSizeScaler:\n";
        cfg << "            scale: 0.2\n";
        cfg << "            weight: 0.0\n";
        cfg << "        FreqMover:\n";
        cfg << "            window: 0.1\n";
        cfg << "            weight: 0.0\n";
        cfg << "global_comparison_settings:\n";
        cfg << "    genotypes_are_diploid: true\n";
        cfg << "    markers_are_dominant: false\n";
        cfg << "    population_name_delimiter: \" \"\n";
        cfg << "    population_name_is_prefix: true\n";
        cfg << "    constant_sites_removed: true\n";
        cfg << "    use_empirical_starting_value_for_freq_1: false\n";
        cfg << "    equal_population_sizes: false\n";
        cfg << "    equal_state_frequencies: false\n";
        cfg << "comparisons:\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129.nex\n";
        cfg << "    parameters:\n";
        cfg << "        population_size:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << size1_shape << "\n";
        cfg << "                    scale: " << size1_scale << "\n";
        cfg << "        freq_1:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                beta_distribution:\n";
        cfg << "                    alpha: " << f1_a << "\n";
        cfg << "                    beta: " << f1_b << "\n";
        cfg << "        mutation_rate:\n";
        cfg << "            value: 1.0\n";
        cfg << "            estimate: false\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129-altname1.nex\n";
        cfg << "    parameters:\n";
        cfg << "        population_size:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << size2_shape << "\n";
        cfg << "                    scale: " << size2_scale << "\n";
        cfg << "        freq_1:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                beta_distribution:\n";
        cfg << "                    alpha: " << f2_a << "\n";
        cfg << "                    beta: " << f2_b << "\n";
        cfg << "        mutation_rate:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << mult2_shape << "\n";
        cfg << "                    scale: " << mult2_scale << "\n";
        cfg << "- comparison:\n";
        cfg << "    path: hemi129-altname2.nex\n";
        cfg << "    parameters:\n";
        cfg << "        population_size:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << size3_shape << "\n";
        cfg << "                    scale: " << size3_scale << "\n";
        cfg << "        freq_1:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                beta_distribution:\n";
        cfg << "                    alpha: " << f3_a << "\n";
        cfg << "                    beta: " << f3_b << "\n";
        cfg << "        mutation_rate:\n";
        cfg << "            estimate: true\n";
        cfg << "            prior:\n";
        cfg << "                gamma_distribution:\n";
        cfg << "                    shape: " << mult3_shape << "\n";
        cfg << "                    scale: " << mult3_scale << "\n";

        CollectionSettings settings = CollectionSettings(cfg, cfg_path);
        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        ComparisonPopulationTreeCollection collection = ComparisonPopulationTreeCollection(
                settings,
                rng);

        std::ofstream out_stream;
        out_stream.open(log_path);
        collection.write_state_log_header(out_stream);

        unsigned int nsamples = 10000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            collection.draw_from_prior(rng);
            collection.log_state(out_stream, i);
        }
        out_stream.close();
        

        REQUIRE(path::exists(log_path));

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);

        unsigned int expected_sample_size = nsamples;

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
    }
}
