#include "catch.hpp"
#include "ecoevolity/operator.hpp"

#include <limits>
#include <memory>

#include "ecoevolity/probability.hpp"
#include "ecoevolity/parameter.hpp"
#include "ecoevolity/stats_util.hpp"
#include "ecoevolity/rng.hpp"

RandomNumberGenerator _TEST_OPERATOR_RNG = RandomNumberGenerator();

TEST_CASE("Testing HeightSizeScaler with 4 pairs",
        "[HeightSizeScaler]") {

    SECTION("Testing 4 pairs with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "                    shape: " << size_shapes.at(0) << "\n";
        os << "                    scale: " << size_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf2_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_root_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == ntrees);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            REQUIRE(size_leaf2_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf2_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf2_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            REQUIRE(size_root_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_root_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_root_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing HeightSizeScaler with 4 pairs with constrained sizes",
        "[HeightSizeScaler]") {

    SECTION("Testing 4 pairs with constrained sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "                    shape: " << size_shapes.at(0) << "\n";
        os << "                    scale: " << size_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        unsigned int niterations = 400000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == ntrees);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing HeightSizeScaler with 4 pairs with fixed sizes",
        "[HeightSizeScaler]") {

    SECTION("Testing 4 pairs with fixed sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        unsigned int niterations = 400000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == ntrees);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                        REQUIRE(size_leaf1 == 0.01);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(0.01));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing HeightSizeScaler with 4 singletons",
        "[HeightSizeScaler]") {

    SECTION("Testing 4 singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "    path: hemi129-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(0) << "\n";
        os << "                    scale: " << size_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_root_summaries(ntrees);

        double size_leaf1;
        double size_root;
        unsigned int niterations = 800000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == ntrees);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                    size_leaf1 = tree.get_child_population_size(0);
                    size_root = tree.get_root_population_size();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(tree.get_leaf_node_count() == 1);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            REQUIRE(size_root_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_root_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_root_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing HeightSizeScaler with mix of pairs and singletons",
        "[HeightSizeScaler]") {

    SECTION("Testing mix of pairs and singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "                    shape: " << size_shapes.at(0) << "\n";
        os << "                    scale: " << size_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname4-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 5);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        SampleSummarizer<double> height_summary_pair1;
        SampleSummarizer<double> height_summary_pair2;
        SampleSummarizer<double> height_summary_pair3;
        SampleSummarizer<double> height_summary_single1;
        SampleSummarizer<double> height_summary_single2;
        SampleSummarizer<double> size_root_summary_pair1;
        SampleSummarizer<double> size_root_summary_pair2;
        SampleSummarizer<double> size_root_summary_pair3;
        SampleSummarizer<double> size_root_summary_single1;
        SampleSummarizer<double> size_root_summary_single2;
        SampleSummarizer<double> size_leaf_summary_pair1;
        SampleSummarizer<double> size_leaf_summary_single1;
        SampleSummarizer<double> size_leaf_summary_single2;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                REQUIRE(comparisons.get_number_of_events() == ntrees);

                height_summary_pair1.add_sample(comparisons.get_tree(0).get_height());
                height_summary_pair2.add_sample(comparisons.get_tree(1).get_height());
                height_summary_pair3.add_sample(comparisons.get_tree(2).get_height());
                height_summary_single1.add_sample(comparisons.get_tree(3).get_height());
                height_summary_single2.add_sample(comparisons.get_tree(4).get_height());

                double  size_root_pair1 = comparisons.get_tree(0).get_root_population_size();
                double size_leaf1_pair1 = comparisons.get_tree(0).get_child_population_size(0);
                double size_leaf2_pair1 = comparisons.get_tree(0).get_child_population_size(1);

                double  size_root_pair2 = comparisons.get_tree(1).get_root_population_size();
                double size_leaf1_pair2 = comparisons.get_tree(1).get_child_population_size(0);
                double size_leaf2_pair2 = comparisons.get_tree(1).get_child_population_size(1);

                double  size_root_pair3 = comparisons.get_tree(2).get_root_population_size();
                double size_leaf1_pair3 = comparisons.get_tree(2).get_child_population_size(0);
                double size_leaf2_pair3 = comparisons.get_tree(2).get_child_population_size(1);

                double size_root_single1 = comparisons.get_tree(3).get_root_population_size();
                double size_leaf_single1 = comparisons.get_tree(3).get_child_population_size(0);

                double size_root_single2 = comparisons.get_tree(4).get_root_population_size();
                double size_leaf_single2 = comparisons.get_tree(4).get_child_population_size(0);

                if (i > (niterations - (sample_freq * 5))) {
                    REQUIRE(size_leaf1_pair1 != size_leaf2_pair1);
                    REQUIRE(size_leaf1_pair1 != size_root_pair1);
                    REQUIRE(size_leaf2_pair1 != size_root_pair1);

                    REQUIRE(size_leaf1_pair2 == size_leaf2_pair2);
                    REQUIRE(size_leaf1_pair2 == size_root_pair2);

                    REQUIRE(size_leaf1_pair3 == 0.01);
                    REQUIRE(size_leaf1_pair3 == size_leaf2_pair3);
                    REQUIRE(size_leaf1_pair3 == size_root_pair3);

                    REQUIRE(comparisons.get_tree(3).get_leaf_node_count() == 1);
                    REQUIRE(size_leaf_single1 != size_root_single1);

                    REQUIRE(comparisons.get_tree(4).get_leaf_node_count() == 1);
                    REQUIRE(size_leaf_single2 != size_root_single2);
                }

                size_leaf_summary_pair1.add_sample(size_leaf1_pair1);
                size_leaf_summary_pair1.add_sample(size_leaf2_pair1);

                size_leaf_summary_single1.add_sample(size_leaf_single1);
                size_leaf_summary_single2.add_sample(size_leaf_single2);

                size_root_summary_pair1.add_sample(size_root_pair1);
                size_root_summary_pair2.add_sample(size_root_pair2);
                size_root_summary_pair3.add_sample(size_root_pair3);
                size_root_summary_single1.add_sample(size_root_single1);
                size_root_summary_single2.add_sample(size_root_single2);
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        REQUIRE(height_summary_pair1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_pair1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(height_summary_pair2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_pair2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(height_summary_pair3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_pair3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(height_summary_single1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_single1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(height_summary_single2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_single2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        double size_sh = size_shapes.at(0);
        double size_sc = size_scales.at(0);
        REQUIRE(size_leaf_summary_pair1.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_leaf_summary_pair1.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_pair1.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_pair1.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));

        size_sh = size_shapes.at(1);
        size_sc = size_scales.at(1);
        REQUIRE(size_root_summary_pair2.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_pair2.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));

        REQUIRE(size_root_summary_pair3.mean() == Approx(0.01));
        REQUIRE(size_root_summary_pair3.variance() == Approx(0.0));

        size_sh = size_shapes.at(2);
        size_sc = size_scales.at(2);
        REQUIRE(size_leaf_summary_single1.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_leaf_summary_single1.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_single1.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_single1.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));

        size_sh = size_shapes.at(3);
        size_sc = size_scales.at(3);
        REQUIRE(size_leaf_summary_single2.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_leaf_summary_single2.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_single2.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_single2.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
    }
}

TEST_CASE("Testing CompositeHeightSizeScaler with 4 pairs",
        "[CompositeHeightSizeScaler]") {

    SECTION("Testing 4 pairs with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "                    shape: " << size_shapes.at(0) << "\n";
        os << "                    scale: " << size_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf2_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_root_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == ntrees);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            REQUIRE(size_leaf2_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf2_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf2_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            REQUIRE(size_root_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_root_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_root_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeScaler with 4 pairs with constrained sizes",
        "[CompositeHeightSizeScaler]") {

    SECTION("Testing 4 pairs with constrained sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "                    shape: " << size_shapes.at(0) << "\n";
        os << "                    scale: " << size_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        unsigned int niterations = 400000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == ntrees);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeScaler with 4 pairs with fixed sizes",
        "[CompositeHeightSizeScaler]") {

    SECTION("Testing 4 pairs with fixed sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        unsigned int niterations = 400000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == ntrees);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                        REQUIRE(size_leaf1 == 0.01);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(0.01));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeScaler with 4 singletons",
        "[CompositeHeightSizeScaler]") {

    SECTION("Testing 4 singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "    path: hemi129-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(0) << "\n";
        os << "                    scale: " << size_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_root_summaries(ntrees);

        double size_leaf1;
        double size_root;
        unsigned int niterations = 800000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == ntrees);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                    size_leaf1 = tree.get_child_population_size(0);
                    size_root = tree.get_root_population_size();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(tree.get_leaf_node_count() == 1);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            REQUIRE(size_root_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_root_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_root_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeScaler with mix of pairs and singletons",
        "[CompositeHeightSizeScaler]") {

    SECTION("Testing mix of pairs and singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "                    shape: " << size_shapes.at(0) << "\n";
        os << "                    scale: " << size_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname4-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 5);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        SampleSummarizer<double> height_summary_pair1;
        SampleSummarizer<double> height_summary_pair2;
        SampleSummarizer<double> height_summary_pair3;
        SampleSummarizer<double> height_summary_single1;
        SampleSummarizer<double> height_summary_single2;
        SampleSummarizer<double> size_root_summary_pair1;
        SampleSummarizer<double> size_root_summary_pair2;
        SampleSummarizer<double> size_root_summary_pair3;
        SampleSummarizer<double> size_root_summary_single1;
        SampleSummarizer<double> size_root_summary_single2;
        SampleSummarizer<double> size_leaf_summary_pair1;
        SampleSummarizer<double> size_leaf_summary_single1;
        SampleSummarizer<double> size_leaf_summary_single2;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                REQUIRE(comparisons.get_number_of_events() == ntrees);

                height_summary_pair1.add_sample(comparisons.get_tree(0).get_height());
                height_summary_pair2.add_sample(comparisons.get_tree(1).get_height());
                height_summary_pair3.add_sample(comparisons.get_tree(2).get_height());
                height_summary_single1.add_sample(comparisons.get_tree(3).get_height());
                height_summary_single2.add_sample(comparisons.get_tree(4).get_height());

                double  size_root_pair1 = comparisons.get_tree(0).get_root_population_size();
                double size_leaf1_pair1 = comparisons.get_tree(0).get_child_population_size(0);
                double size_leaf2_pair1 = comparisons.get_tree(0).get_child_population_size(1);

                double  size_root_pair2 = comparisons.get_tree(1).get_root_population_size();
                double size_leaf1_pair2 = comparisons.get_tree(1).get_child_population_size(0);
                double size_leaf2_pair2 = comparisons.get_tree(1).get_child_population_size(1);

                double  size_root_pair3 = comparisons.get_tree(2).get_root_population_size();
                double size_leaf1_pair3 = comparisons.get_tree(2).get_child_population_size(0);
                double size_leaf2_pair3 = comparisons.get_tree(2).get_child_population_size(1);

                double size_root_single1 = comparisons.get_tree(3).get_root_population_size();
                double size_leaf_single1 = comparisons.get_tree(3).get_child_population_size(0);

                double size_root_single2 = comparisons.get_tree(4).get_root_population_size();
                double size_leaf_single2 = comparisons.get_tree(4).get_child_population_size(0);

                if (i > (niterations - (sample_freq * 5))) {
                    REQUIRE(size_leaf1_pair1 != size_leaf2_pair1);
                    REQUIRE(size_leaf1_pair1 != size_root_pair1);
                    REQUIRE(size_leaf2_pair1 != size_root_pair1);

                    REQUIRE(size_leaf1_pair2 == size_leaf2_pair2);
                    REQUIRE(size_leaf1_pair2 == size_root_pair2);

                    REQUIRE(size_leaf1_pair3 == 0.01);
                    REQUIRE(size_leaf1_pair3 == size_leaf2_pair3);
                    REQUIRE(size_leaf1_pair3 == size_root_pair3);

                    REQUIRE(comparisons.get_tree(3).get_leaf_node_count() == 1);
                    REQUIRE(size_leaf_single1 != size_root_single1);

                    REQUIRE(comparisons.get_tree(4).get_leaf_node_count() == 1);
                    REQUIRE(size_leaf_single2 != size_root_single2);
                }

                size_leaf_summary_pair1.add_sample(size_leaf1_pair1);
                size_leaf_summary_pair1.add_sample(size_leaf2_pair1);

                size_leaf_summary_single1.add_sample(size_leaf_single1);
                size_leaf_summary_single2.add_sample(size_leaf_single2);

                size_root_summary_pair1.add_sample(size_root_pair1);
                size_root_summary_pair2.add_sample(size_root_pair2);
                size_root_summary_pair3.add_sample(size_root_pair3);
                size_root_summary_single1.add_sample(size_root_single1);
                size_root_summary_single2.add_sample(size_root_single2);
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        REQUIRE(height_summary_pair1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_pair1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(height_summary_pair2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_pair2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(height_summary_pair3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_pair3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(height_summary_single1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_single1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(height_summary_single2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_single2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        double size_sh = size_shapes.at(0);
        double size_sc = size_scales.at(0);
        REQUIRE(size_leaf_summary_pair1.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_leaf_summary_pair1.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_pair1.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_pair1.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));

        size_sh = size_shapes.at(1);
        size_sc = size_scales.at(1);
        REQUIRE(size_root_summary_pair2.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_pair2.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));

        REQUIRE(size_root_summary_pair3.mean() == Approx(0.01));
        REQUIRE(size_root_summary_pair3.variance() == Approx(0.0));

        size_sh = size_shapes.at(2);
        size_sc = size_scales.at(2);
        REQUIRE(size_leaf_summary_single1.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_leaf_summary_single1.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_single1.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_single1.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));

        size_sh = size_shapes.at(3);
        size_sc = size_scales.at(3);
        REQUIRE(size_leaf_summary_single2.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_leaf_summary_single2.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_single2.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_single2.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
    }
}

TEST_CASE("Testing SmartHeightSizeMixer with 4 pairs",
        "[SmartHeightSizeMixer]") {

    SECTION("Testing 4 pairs with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "                    shape: " << size_shapes.at(0) << "\n";
        os << "                    scale: " << size_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(12345678);
        std::shared_ptr<OperatorInterface> op = std::make_shared<SmartHeightSizeMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf2_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_root_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        unsigned int niterations = 4000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == ntrees);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            REQUIRE(size_leaf2_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf2_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf2_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            REQUIRE(size_root_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_root_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_root_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing SmartHeightSizeMixer with 4 pairs with constrained sizes",
        "[SmartHeightSizeMixer]") {

    SECTION("Testing 4 pairs with constrained sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "                    shape: " << size_shapes.at(0) << "\n";
        os << "                    scale: " << size_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<SmartHeightSizeMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        unsigned int niterations = 400000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == ntrees);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing SmartHeightSizeMixer with 4 pairs with fixed sizes",
        "[SmartHeightSizeMixer]") {

    SECTION("Testing 4 pairs with fixed sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<SmartHeightSizeMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        unsigned int niterations = 400000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == ntrees);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                        REQUIRE(size_leaf1 == 0.01);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(0.01));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing SmartHeightSizeMixer with 4 singletons",
        "[SmartHeightSizeMixer]") {

    SECTION("Testing 4 singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "    path: hemi129-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(0) << "\n";
        os << "                    scale: " << size_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456789);
        std::shared_ptr<OperatorInterface> op = std::make_shared<SmartHeightSizeMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_root_summaries(ntrees);

        double size_leaf1;
        double size_root;
        unsigned int niterations = 1000000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == ntrees);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                    size_leaf1 = tree.get_child_population_size(0);
                    size_root = tree.get_root_population_size();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(tree.get_leaf_node_count() == 1);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            REQUIRE(size_root_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_root_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_root_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing SmartHeightSizeMixer with mix of pairs and singletons",
        "[SmartHeightSizeMixer]") {

    SECTION("Testing mix of pairs and singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "                    shape: " << size_shapes.at(0) << "\n";
        os << "                    scale: " << size_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname4-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<SmartHeightSizeMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 5);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        SampleSummarizer<double> height_summary_pair1;
        SampleSummarizer<double> height_summary_pair2;
        SampleSummarizer<double> height_summary_pair3;
        SampleSummarizer<double> height_summary_single1;
        SampleSummarizer<double> height_summary_single2;
        SampleSummarizer<double> size_root_summary_pair1;
        SampleSummarizer<double> size_root_summary_pair2;
        SampleSummarizer<double> size_root_summary_pair3;
        SampleSummarizer<double> size_root_summary_single1;
        SampleSummarizer<double> size_root_summary_single2;
        SampleSummarizer<double> size_leaf_summary_pair1;
        SampleSummarizer<double> size_leaf_summary_single1;
        SampleSummarizer<double> size_leaf_summary_single2;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                REQUIRE(comparisons.get_number_of_events() == ntrees);

                height_summary_pair1.add_sample(comparisons.get_tree(0).get_height());
                height_summary_pair2.add_sample(comparisons.get_tree(1).get_height());
                height_summary_pair3.add_sample(comparisons.get_tree(2).get_height());
                height_summary_single1.add_sample(comparisons.get_tree(3).get_height());
                height_summary_single2.add_sample(comparisons.get_tree(4).get_height());

                double  size_root_pair1 = comparisons.get_tree(0).get_root_population_size();
                double size_leaf1_pair1 = comparisons.get_tree(0).get_child_population_size(0);
                double size_leaf2_pair1 = comparisons.get_tree(0).get_child_population_size(1);

                double  size_root_pair2 = comparisons.get_tree(1).get_root_population_size();
                double size_leaf1_pair2 = comparisons.get_tree(1).get_child_population_size(0);
                double size_leaf2_pair2 = comparisons.get_tree(1).get_child_population_size(1);

                double  size_root_pair3 = comparisons.get_tree(2).get_root_population_size();
                double size_leaf1_pair3 = comparisons.get_tree(2).get_child_population_size(0);
                double size_leaf2_pair3 = comparisons.get_tree(2).get_child_population_size(1);

                double size_root_single1 = comparisons.get_tree(3).get_root_population_size();
                double size_leaf_single1 = comparisons.get_tree(3).get_child_population_size(0);

                double size_root_single2 = comparisons.get_tree(4).get_root_population_size();
                double size_leaf_single2 = comparisons.get_tree(4).get_child_population_size(0);

                if (i > (niterations - (sample_freq * 5))) {
                    REQUIRE(size_leaf1_pair1 != size_leaf2_pair1);
                    REQUIRE(size_leaf1_pair1 != size_root_pair1);
                    REQUIRE(size_leaf2_pair1 != size_root_pair1);

                    REQUIRE(size_leaf1_pair2 == size_leaf2_pair2);
                    REQUIRE(size_leaf1_pair2 == size_root_pair2);

                    REQUIRE(size_leaf1_pair3 == 0.01);
                    REQUIRE(size_leaf1_pair3 == size_leaf2_pair3);
                    REQUIRE(size_leaf1_pair3 == size_root_pair3);

                    REQUIRE(comparisons.get_tree(3).get_leaf_node_count() == 1);
                    REQUIRE(size_leaf_single1 != size_root_single1);

                    REQUIRE(comparisons.get_tree(4).get_leaf_node_count() == 1);
                    REQUIRE(size_leaf_single2 != size_root_single2);
                }

                size_leaf_summary_pair1.add_sample(size_leaf1_pair1);
                size_leaf_summary_pair1.add_sample(size_leaf2_pair1);

                size_leaf_summary_single1.add_sample(size_leaf_single1);
                size_leaf_summary_single2.add_sample(size_leaf_single2);

                size_root_summary_pair1.add_sample(size_root_pair1);
                size_root_summary_pair2.add_sample(size_root_pair2);
                size_root_summary_pair3.add_sample(size_root_pair3);
                size_root_summary_single1.add_sample(size_root_single1);
                size_root_summary_single2.add_sample(size_root_single2);
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        REQUIRE(height_summary_pair1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_pair1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(height_summary_pair2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_pair2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(height_summary_pair3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_pair3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(height_summary_single1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_single1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(height_summary_single2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_single2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        double size_sh = size_shapes.at(0);
        double size_sc = size_scales.at(0);
        REQUIRE(size_leaf_summary_pair1.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_leaf_summary_pair1.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_pair1.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_pair1.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));

        size_sh = size_shapes.at(1);
        size_sc = size_scales.at(1);
        REQUIRE(size_root_summary_pair2.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_pair2.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));

        REQUIRE(size_root_summary_pair3.mean() == Approx(0.01));
        REQUIRE(size_root_summary_pair3.variance() == Approx(0.0));

        size_sh = size_shapes.at(2);
        size_sc = size_scales.at(2);
        REQUIRE(size_leaf_summary_single1.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_leaf_summary_single1.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_single1.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_single1.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));

        size_sh = size_shapes.at(3);
        size_sc = size_scales.at(3);
        REQUIRE(size_leaf_summary_single2.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_leaf_summary_single2.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_single2.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_single2.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
    }
}

TEST_CASE("Testing CompositeSmartHeightSizeMixer with 4 pairs",
        "[CompositeSmartHeightSizeMixer]") {

    SECTION("Testing 4 pairs with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "                    shape: " << size_shapes.at(0) << "\n";
        os << "                    scale: " << size_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeSmartHeightSizeMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf2_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_root_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == ntrees);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            REQUIRE(size_leaf2_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf2_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf2_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            REQUIRE(size_root_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_root_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_root_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeSmartHeightSizeMixer with 4 pairs with constrained sizes",
        "[CompositeSmartHeightSizeMixer]") {

    SECTION("Testing 4 pairs with constrained sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "                    shape: " << size_shapes.at(0) << "\n";
        os << "                    scale: " << size_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeSmartHeightSizeMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        unsigned int niterations = 400000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == ntrees);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeSmartHeightSizeMixer with 4 pairs with fixed sizes",
        "[CompositeSmartHeightSizeMixer]") {

    SECTION("Testing 4 pairs with fixed sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeSmartHeightSizeMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        unsigned int niterations = 400000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == ntrees);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                        REQUIRE(size_leaf1 == 0.01);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(0.01));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing CompositeSmartHeightSizeMixer with 4 singletons",
        "[CompositeSmartHeightSizeMixer]") {

    SECTION("Testing 4 singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "    path: hemi129-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(0) << "\n";
        os << "                    scale: " << size_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeSmartHeightSizeMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_root_summaries(ntrees);

        double size_leaf1;
        double size_root;
        unsigned int niterations = 800000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == ntrees);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                    size_leaf1 = tree.get_child_population_size(0);
                    size_root = tree.get_root_population_size();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(tree.get_leaf_node_count() == 1);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            REQUIRE(size_root_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_root_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_root_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeSmartHeightSizeMixer with mix of pairs and singletons",
        "[CompositeSmartHeightSizeMixer]") {

    SECTION("Testing mix of pairs and singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
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
        os << "                    shape: " << size_shapes.at(0) << "\n";
        os << "                    scale: " << size_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname4-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeSmartHeightSizeMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 5);
        REQUIRE(comparisons.get_number_of_events() == ntrees);
        SampleSummarizer<double> height_summary_pair1;
        SampleSummarizer<double> height_summary_pair2;
        SampleSummarizer<double> height_summary_pair3;
        SampleSummarizer<double> height_summary_single1;
        SampleSummarizer<double> height_summary_single2;
        SampleSummarizer<double> size_root_summary_pair1;
        SampleSummarizer<double> size_root_summary_pair2;
        SampleSummarizer<double> size_root_summary_pair3;
        SampleSummarizer<double> size_root_summary_single1;
        SampleSummarizer<double> size_root_summary_single2;
        SampleSummarizer<double> size_leaf_summary_pair1;
        SampleSummarizer<double> size_leaf_summary_single1;
        SampleSummarizer<double> size_leaf_summary_single2;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                REQUIRE(comparisons.get_number_of_events() == ntrees);

                height_summary_pair1.add_sample(comparisons.get_tree(0).get_height());
                height_summary_pair2.add_sample(comparisons.get_tree(1).get_height());
                height_summary_pair3.add_sample(comparisons.get_tree(2).get_height());
                height_summary_single1.add_sample(comparisons.get_tree(3).get_height());
                height_summary_single2.add_sample(comparisons.get_tree(4).get_height());

                double  size_root_pair1 = comparisons.get_tree(0).get_root_population_size();
                double size_leaf1_pair1 = comparisons.get_tree(0).get_child_population_size(0);
                double size_leaf2_pair1 = comparisons.get_tree(0).get_child_population_size(1);

                double  size_root_pair2 = comparisons.get_tree(1).get_root_population_size();
                double size_leaf1_pair2 = comparisons.get_tree(1).get_child_population_size(0);
                double size_leaf2_pair2 = comparisons.get_tree(1).get_child_population_size(1);

                double  size_root_pair3 = comparisons.get_tree(2).get_root_population_size();
                double size_leaf1_pair3 = comparisons.get_tree(2).get_child_population_size(0);
                double size_leaf2_pair3 = comparisons.get_tree(2).get_child_population_size(1);

                double size_root_single1 = comparisons.get_tree(3).get_root_population_size();
                double size_leaf_single1 = comparisons.get_tree(3).get_child_population_size(0);

                double size_root_single2 = comparisons.get_tree(4).get_root_population_size();
                double size_leaf_single2 = comparisons.get_tree(4).get_child_population_size(0);

                if (i > (niterations - (sample_freq * 5))) {
                    REQUIRE(size_leaf1_pair1 != size_leaf2_pair1);
                    REQUIRE(size_leaf1_pair1 != size_root_pair1);
                    REQUIRE(size_leaf2_pair1 != size_root_pair1);

                    REQUIRE(size_leaf1_pair2 == size_leaf2_pair2);
                    REQUIRE(size_leaf1_pair2 == size_root_pair2);

                    REQUIRE(size_leaf1_pair3 == 0.01);
                    REQUIRE(size_leaf1_pair3 == size_leaf2_pair3);
                    REQUIRE(size_leaf1_pair3 == size_root_pair3);

                    REQUIRE(comparisons.get_tree(3).get_leaf_node_count() == 1);
                    REQUIRE(size_leaf_single1 != size_root_single1);

                    REQUIRE(comparisons.get_tree(4).get_leaf_node_count() == 1);
                    REQUIRE(size_leaf_single2 != size_root_single2);
                }

                size_leaf_summary_pair1.add_sample(size_leaf1_pair1);
                size_leaf_summary_pair1.add_sample(size_leaf2_pair1);

                size_leaf_summary_single1.add_sample(size_leaf_single1);
                size_leaf_summary_single2.add_sample(size_leaf_single2);

                size_root_summary_pair1.add_sample(size_root_pair1);
                size_root_summary_pair2.add_sample(size_root_pair2);
                size_root_summary_pair3.add_sample(size_root_pair3);
                size_root_summary_single1.add_sample(size_root_single1);
                size_root_summary_single2.add_sample(size_root_single2);
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        REQUIRE(height_summary_pair1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_pair1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(height_summary_pair2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_pair2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(height_summary_pair3.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_pair3.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(height_summary_single1.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_single1.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        REQUIRE(height_summary_single2.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary_single2.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        double size_sh = size_shapes.at(0);
        double size_sc = size_scales.at(0);
        REQUIRE(size_leaf_summary_pair1.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_leaf_summary_pair1.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_pair1.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_pair1.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));

        size_sh = size_shapes.at(1);
        size_sc = size_scales.at(1);
        REQUIRE(size_root_summary_pair2.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_pair2.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));

        REQUIRE(size_root_summary_pair3.mean() == Approx(0.01));
        REQUIRE(size_root_summary_pair3.variance() == Approx(0.0));

        size_sh = size_shapes.at(2);
        size_sc = size_scales.at(2);
        REQUIRE(size_leaf_summary_single1.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_leaf_summary_single1.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_single1.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_single1.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));

        size_sh = size_shapes.at(3);
        size_sc = size_scales.at(3);
        REQUIRE(size_leaf_summary_single2.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_leaf_summary_single2.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_single2.mean() == Approx(size_sh * size_sc).epsilon(0.01));
        REQUIRE(size_root_summary_single2.variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
    }
}

TEST_CASE("Testing ComparisonHeightScaler", "[ComparisonHeightScaler]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        double shape = 10.0;
        double scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream out;
        out.open(test_path);
        out << "event_time_prior:\n";
        out << "    gamma_distribution:\n";
        out << "        shape: " << shape << "\n";
        out << "        scale: " << scale << "\n";
        out << "global_comparison_settings:\n";
        out << "    genotypes_are_diploid: true\n";
        out << "    markers_are_dominant: false\n";
        out << "    population_name_delimiter: \" \"\n";
        out << "    population_name_is_prefix: true\n";
        out << "    constant_sites_removed: true\n";
        out << "    use_empirical_starting_value_for_freq_1: false\n";
        out << "    equal_population_sizes: true\n";
        out << "    equal_state_frequencies: true\n";
        out << "    parameters:\n";
        out << "        population_size:\n";
        out << "            value: 0.005\n";
        out << "            estimate: false\n";
        out << "        freq_1:\n";
        out << "            value: 0.5\n";
        out << "            estimate: false\n";
        out << "        mutation_rate:\n";
        out << "            value: 1.0\n";
        out << "            estimate: false\n";
        out << "comparisons:\n";
        out << "- comparison:\n";
        out << "    path: hemi129.nex\n";
        out.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ComparisonHeightScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        os.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(os);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double h;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            o.perform_collection_move(rng, comparisons, 1);
            h = comparisons.get_height(0);
            mn = std::min(mn, h);
            mx = std::max(mx, h);
            ++n;
            d = h - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        
        REQUIRE(mean == Approx(shape * scale).epsilon(0.001));
        REQUIRE(variance == Approx(shape * scale * scale).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }

    SECTION("Testing gamma(10.0, 0.1) prior and with optimizing") {
        double shape = 10.0;
        double scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream out;
        out.open(test_path);
        out << "event_time_prior:\n";
        out << "    gamma_distribution:\n";
        out << "        shape: " << shape << "\n";
        out << "        scale: " << scale << "\n";
        out << "global_comparison_settings:\n";
        out << "    genotypes_are_diploid: true\n";
        out << "    markers_are_dominant: false\n";
        out << "    population_name_delimiter: \" \"\n";
        out << "    population_name_is_prefix: true\n";
        out << "    constant_sites_removed: true\n";
        out << "    use_empirical_starting_value_for_freq_1: false\n";
        out << "    equal_population_sizes: true\n";
        out << "    equal_state_frequencies: true\n";
        out << "    parameters:\n";
        out << "        population_size:\n";
        out << "            value: 0.005\n";
        out << "            estimate: false\n";
        out << "        freq_1:\n";
        out << "            value: 0.5\n";
        out << "            estimate: false\n";
        out << "        mutation_rate:\n";
        out << "            value: 1.0\n";
        out << "            estimate: false\n";
        out << "comparisons:\n";
        out << "- comparison:\n";
        out << "    path: hemi129.nex\n";
        out.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ComparisonHeightScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(os);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double h;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            o.perform_collection_move(rng, comparisons, 1);
            h = comparisons.get_height(0);
            mn = std::min(mn, h);
            mx = std::max(mx, h);
            ++n;
            d = h - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        
        REQUIRE(mean == Approx(shape * scale).epsilon(0.001));
        REQUIRE(variance == Approx(shape * scale * scale).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonHeightMover", "[ComparisonHeightMover]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        double shape = 10.0;
        double scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream out;
        out.open(test_path);
        out << "event_time_prior:\n";
        out << "    gamma_distribution:\n";
        out << "        shape: " << shape << "\n";
        out << "        scale: " << scale << "\n";
        out << "global_comparison_settings:\n";
        out << "    genotypes_are_diploid: true\n";
        out << "    markers_are_dominant: false\n";
        out << "    population_name_delimiter: \" \"\n";
        out << "    population_name_is_prefix: true\n";
        out << "    constant_sites_removed: true\n";
        out << "    use_empirical_starting_value_for_freq_1: false\n";
        out << "    equal_population_sizes: true\n";
        out << "    equal_state_frequencies: true\n";
        out << "    parameters:\n";
        out << "        population_size:\n";
        out << "            value: 0.005\n";
        out << "            estimate: false\n";
        out << "        freq_1:\n";
        out << "            value: 0.5\n";
        out << "            estimate: false\n";
        out << "        mutation_rate:\n";
        out << "            value: 1.0\n";
        out << "            estimate: false\n";
        out << "comparisons:\n";
        out << "- comparison:\n";
        out << "    path: hemi129.nex\n";
        out.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ComparisonHeightMover>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        os.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(os);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double h;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            o.perform_collection_move(rng, comparisons, 1);
            h = comparisons.get_height(0);
            mn = std::min(mn, h);
            mx = std::max(mx, h);
            ++n;
            d = h - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        
        REQUIRE(mean == Approx(shape * scale).epsilon(0.001));
        REQUIRE(variance == Approx(shape * scale * scale).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }

    SECTION("Testing gamma(10.0, 0.1) prior and with optimizing") {
        double shape = 10.0;
        double scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream out;
        out.open(test_path);
        out << "event_time_prior:\n";
        out << "    gamma_distribution:\n";
        out << "        shape: " << shape << "\n";
        out << "        scale: " << scale << "\n";
        out << "global_comparison_settings:\n";
        out << "    genotypes_are_diploid: true\n";
        out << "    markers_are_dominant: false\n";
        out << "    population_name_delimiter: \" \"\n";
        out << "    population_name_is_prefix: true\n";
        out << "    constant_sites_removed: true\n";
        out << "    use_empirical_starting_value_for_freq_1: false\n";
        out << "    equal_population_sizes: true\n";
        out << "    equal_state_frequencies: true\n";
        out << "    parameters:\n";
        out << "        population_size:\n";
        out << "            value: 0.005\n";
        out << "            estimate: false\n";
        out << "        freq_1:\n";
        out << "            value: 0.5\n";
        out << "            estimate: false\n";
        out << "        mutation_rate:\n";
        out << "            value: 1.0\n";
        out << "            estimate: false\n";
        out << "comparisons:\n";
        out << "- comparison:\n";
        out << "    path: hemi129.nex\n";
        out.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ComparisonHeightMover>(1.0, 0.2);
        OperatorSchedule os = OperatorSchedule();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(os);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double h;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            o.perform_collection_move(rng, comparisons, 1);
            h = comparisons.get_height(0);
            mn = std::min(mn, h);
            mx = std::max(mx, h);
            ++n;
            d = h - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        
        REQUIRE(mean == Approx(shape * scale).epsilon(0.001));
        REQUIRE(variance == Approx(shape * scale * scale).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}


TEST_CASE("Testing RootPopulationSizeScaler", "[RootPopulationSizeScaler]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(12933);
        std::shared_ptr<OperatorInterface> op = std::make_shared<RootPopulationSizeScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_root_population_size(1.0);
        tree.set_population_size_prior(prior);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            // REQUIRE(tree.is_dirty());
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_root_population_size();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
        
    }

    SECTION("Testing gamma(10.0, 0.1) prior with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(12933);
        std::shared_ptr<OperatorInterface> op = std::make_shared<RootPopulationSizeScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        //os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_root_population_size(1.0);
        tree.set_population_size_prior(prior);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            // REQUIRE(tree.is_dirty());
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_root_population_size();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
        
    }
}

TEST_CASE("Testing ChildPopulationSizeScaler", "[ChildPopulationSizeScaler]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(12933);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ChildPopulationSizeScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.constrain_population_sizes();
        tree.set_child_population_size(0, 1.0);
        tree.set_population_size_prior(prior);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            // REQUIRE(tree.is_dirty());
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_child_population_size(0);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
    }

    SECTION("Testing gamma(10.0, 0.1) prior with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(12933);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ChildPopulationSizeScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        // os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.constrain_population_sizes();
        tree.set_child_population_size(0, 1.0);
        tree.set_population_size_prior(prior);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            // REQUIRE(tree.is_dirty());
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_child_population_size(0);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
    }
}

TEST_CASE("Testing ComparisonMutationRateScaler", "[ComparisonMutationRateScaler]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(928374);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ComparisonMutationRateScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_state_frequencies();
        tree.fix_population_sizes();
        tree.set_mutation_rate(1.0);
        tree.set_mutation_rate_prior(prior);
        tree.estimate_mutation_rate();
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double old_v = tree.get_mutation_rate();
            double hastings = o.propose(rng, tree);
            double new_v = tree.get_mutation_rate();
            REQUIRE(tree.is_dirty());
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            /* double prior_ratio = */
            /*     prior->relative_ln_pdf(new_v) - */
            /*     prior->relative_ln_pdf(old_v); */
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_mutation_rate();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
        
    }

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(928374);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ComparisonMutationRateScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        // os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_state_frequencies();
        tree.fix_population_sizes();
        tree.set_mutation_rate(1.0);
        tree.set_mutation_rate_prior(prior);
        tree.estimate_mutation_rate();
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double old_v = tree.get_mutation_rate();
            double hastings = o.propose(rng, tree);
            double new_v = tree.get_mutation_rate();
            REQUIRE(tree.is_dirty());
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            /* double prior_ratio = */
            /*     prior->relative_ln_pdf(new_v) - */
            /*     prior->relative_ln_pdf(old_v); */
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_mutation_rate();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
        
    }
}

TEST_CASE("Testing FreqMover", "[FreqMover]") {

    SECTION("testing beta(1.0, 1.0) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(3648);
        std::shared_ptr<OperatorInterface> op = std::make_shared<FreqMover>(1.0, 0.1);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        double a = 1.0;
        double b = 1.0;
        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<BetaDistribution>(a, b);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_population_sizes();
        tree.fix_mutation_rate();
        tree.estimate_state_frequencies();
        tree.set_freq_1_prior(prior);
        tree.set_freq_1(0.5);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(! tree.state_frequencies_are_constrained());
        REQUIRE(! tree.state_frequencies_are_fixed());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_freq_1();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.005));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
    }

    SECTION("testing beta(1.0, 1.0) prior with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(2945720);
        std::shared_ptr<OperatorInterface> op = std::make_shared<FreqMover>(1.0, 0.1);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        double a = 1.0;
        double b = 1.0;
        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<BetaDistribution>(a, b);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_population_sizes();
        tree.fix_mutation_rate();
        tree.estimate_state_frequencies();
        tree.set_freq_1_prior(prior);
        tree.set_freq_1(0.5);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(! tree.state_frequencies_are_constrained());
        REQUIRE(! tree.state_frequencies_are_fixed());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_freq_1();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.005));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
    }

    SECTION("testing beta(5.0, 1.0) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(841984264);
        std::shared_ptr<OperatorInterface> op = std::make_shared<FreqMover>(1.0, 0.1);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        double a = 5.0;
        double b = 1.0;
        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<BetaDistribution>(a, b);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_population_sizes();
        tree.fix_mutation_rate();
        tree.estimate_state_frequencies();
        tree.set_freq_1_prior(prior);
        tree.set_freq_1(0.5);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(! tree.state_frequencies_are_constrained());
        REQUIRE(! tree.state_frequencies_are_fixed());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_freq_1();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.005));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
    }

    SECTION("testing beta(5.0, 1.0) prior with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(25456657);
        std::shared_ptr<OperatorInterface> op = std::make_shared<FreqMover>(1.0, 0.1);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        double a = 5.0;
        double b = 1.0;
        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<BetaDistribution>(a, b);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_population_sizes();
        tree.fix_mutation_rate();
        tree.estimate_state_frequencies();
        tree.set_freq_1_prior(prior);
        tree.set_freq_1(0.5);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(! tree.state_frequencies_are_constrained());
        REQUIRE(! tree.state_frequencies_are_fixed());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_freq_1();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.005));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
    }

    SECTION("testing beta(1.0, 5.0) prior with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(14458);
        std::shared_ptr<OperatorInterface> op = std::make_shared<FreqMover>(1.0, 0.1);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        double a = 1.0;
        double b = 5.0;
        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<BetaDistribution>(a, b);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_population_sizes();
        tree.fix_mutation_rate();
        tree.estimate_state_frequencies();
        tree.set_freq_1_prior(prior);
        tree.set_freq_1(0.5);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(! tree.state_frequencies_are_constrained());
        REQUIRE(! tree.state_frequencies_are_fixed());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_freq_1();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.005));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
    }
}
