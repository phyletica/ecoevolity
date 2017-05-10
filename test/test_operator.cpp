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
        os << "    equal_population_sizes: false\n";
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
        std::string test_path = "data/tmp-config-heightsizescaler-test2-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test2-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
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
        std::string test_path = "data/tmp-config-heightsizescaler-test3-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test3" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
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
        std::string test_path = "data/tmp-config-heightsizescaler-test4-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test14" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::string test_path = "data/tmp-config-heightsizescaler-test5-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test5-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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

TEST_CASE("Testing HeightSizeScaler with 4 pairs and shared event",
        "[HeightSizeScaler]") {

    SECTION("Testing 4 pairs with shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test6-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test6-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
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
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
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
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
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
        
        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
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

TEST_CASE("Testing HeightSizeScaler with mix of pairs and singletons and shared event",
        "[HeightSizeScaler]") {

    SECTION("Testing mix of pairs and singletons with shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test7-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test7-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
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
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
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
                REQUIRE(comparisons.get_number_of_events() == 1);

                height_summary.add_sample(comparisons.get_tree(0).get_height());

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
        
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

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

TEST_CASE("Testing HeightSizeScaler with 4 pairs with constrained sizes and shared event",
        "[HeightSizeScaler]") {

    SECTION("Testing 4 pairs with constrained sizes and shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizescaler-test8-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizescaler-test8-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: true\n";
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
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
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
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
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

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        }
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
        std::string test_path = "data/tmp-config-compheightsizescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizescaler-test1-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::string test_path = "data/tmp-config-compheightsizescaler-test2-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizescaler-test2-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
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
        std::string test_path = "data/tmp-config-compheightsizescaler-test3-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizescaler-test3-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
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
        std::string test_path = "data/tmp-config-compheightsizescaler-test4-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizescaler-test4-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::string test_path = "data/tmp-config-compheightsizescaler-test5-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizescaler-test5-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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

TEST_CASE("Testing CompositeHeightSizeScaler with 4 pairs and shared event",
        "[CompositeHeightSizeScaler]") {

    SECTION("Testing 4 pairs with shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizescaler-test6-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizescaler-test6-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
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
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
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
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
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
        
        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
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

TEST_CASE("Testing CompositeHeightSizeScaler with mix of pairs and singletons and shared event",
        "[CompositeHeightSizeScaler]") {

    SECTION("Testing mix of pairs and singletons with shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizescaler-test7-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizescaler-test7-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
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
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
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
                REQUIRE(comparisons.get_number_of_events() == 1);

                height_summary.add_sample(comparisons.get_tree(0).get_height());

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
        
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

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

TEST_CASE("Testing CompositeHeightSizeScaler with 4 pairs with constrained sizes and shared event",
        "[CompositeHeightSizeScaler]") {

    SECTION("Testing 4 pairs with constrained sizes and shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizescaler-test8-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizescaler-test8-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: true\n";
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
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
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
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
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

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing HeightSizeMixer with 4 pairs",
        "[HeightSizeMixer]") {

    SECTION("Testing 4 pairs with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizemixer-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizemixer-test1-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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

        RandomNumberGenerator rng = RandomNumberGenerator(1234567);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeMixer>(1.0, 0.5);
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
        unsigned int niterations = 900000;
        unsigned int sample_freq = 3;
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

TEST_CASE("Testing HeightSizeMixer with 4 pairs with constrained sizes",
        "[HeightSizeMixer]") {

    SECTION("Testing 4 pairs with constrained sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizemixer-test2-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizemixer-test2-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeMixer>(1.0, 0.5);
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

TEST_CASE("Testing HeightSizeMixer with 4 pairs with fixed sizes",
        "[HeightSizeMixer]") {

    SECTION("Testing 4 pairs with fixed sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizemixer-test3-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizemixer-test3-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeMixer>(1.0, 0.5);
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

TEST_CASE("Testing HeightSizeMixer with 4 singletons",
        "[HeightSizeMixer]") {

    SECTION("Testing 4 singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizemixer-test4-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizemixer-test4-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeMixer>(1.0, 0.5);
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

TEST_CASE("Testing HeightSizeMixer with mix of pairs and singletons",
        "[HeightSizeMixer]") {

    SECTION("Testing mix of pairs and singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizemixer-test5-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizemixer-test5-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeMixer>(1.0, 0.5);
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

TEST_CASE("Testing HeightSizeMixer with 4 pairs and shared event",
        "[HeightSizeMixer]") {

    SECTION("Testing 4 pairs with shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizemixer-test6-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizemixer-test6-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
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

        RandomNumberGenerator rng = RandomNumberGenerator(1234567);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
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
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
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
        
        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
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

TEST_CASE("Testing HeightSizeMixer with mix of pairs and singletons and shared event",
        "[HeightSizeMixer]") {

    SECTION("Testing mix of pairs and singletons with shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizemixer-test7-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizemixer-test7-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 5);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
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
                REQUIRE(comparisons.get_number_of_events() == 1);

                height_summary.add_sample(comparisons.get_tree(0).get_height());

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
        
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

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

TEST_CASE("Testing HeightSizeMixer with 4 pairs with constrained sizes and shared event",
        "[HeightSizeMixer]") {

    SECTION("Testing 4 pairs with constrained sizes and shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizemixer-test8-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizemixer-test8-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: true\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
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
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
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

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeMixer with 4 pairs",
        "[CompositeHeightSizeMixer]") {

    SECTION("Testing 4 pairs with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizemixer-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizemixer-test1-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeMixer>(1.0, 0.5);
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

TEST_CASE("Testing CompositeHeightSizeMixer with 4 pairs with constrained sizes",
        "[CompositeHeightSizeMixer]") {

    SECTION("Testing 4 pairs with constrained sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizemixer-test2-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizemixer-test2-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeMixer>(1.0, 0.5);
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

TEST_CASE("Testing CompositeHeightSizeMixer with 4 pairs with fixed sizes",
        "[CompositeHeightSizeMixer]") {

    SECTION("Testing 4 pairs with fixed sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizemixer-test3-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizemixer-test3-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeMixer>(1.0, 0.5);
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

TEST_CASE("Testing CompositeHeightSizeMixer with 4 singletons",
        "[CompositeHeightSizeMixer]") {

    SECTION("Testing 4 singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizemixer-test4-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizemixer-test4-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeMixer>(1.0, 0.5);
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

TEST_CASE("Testing CompositeHeightSizeMixer with mix of pairs and singletons",
        "[CompositeHeightSizeMixer]") {

    SECTION("Testing mix of pairs and singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizemixer-test5-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizemixer-test5-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeMixer>(1.0, 0.5);
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

TEST_CASE("Testing CompositeHeightSizeMixer with 4 pairs and shared event",
        "[CompositeHeightSizeMixer]") {

    SECTION("Testing 4 pairs with shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizemixer-test6-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizemixer-test6-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
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
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
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
        
        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
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

TEST_CASE("Testing CompositeHeightSizeMixer with mix of pairs and singletons and shared event",
        "[CompositeHeightSizeMixer]") {

    SECTION("Testing mix of pairs and singletons with shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizemixer-test7-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizemixer-test7-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 5);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
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
                REQUIRE(comparisons.get_number_of_events() == 1);

                height_summary.add_sample(comparisons.get_tree(0).get_height());

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
        
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

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

TEST_CASE("Testing CompositeHeightSizeMixer with 4 pairs with constrained sizes and shared event",
        "[CompositeHeightSizeMixer]") {

    SECTION("Testing 4 pairs with constrained sizes and shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizemixer-test8-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizemixer-test8-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: true\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
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
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
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

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing HeightSizeRateScaler with 4 pairs",
        "[HeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratescaler-test1-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateScaler>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
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
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing HeightSizeRateScaler with 4 pairs with constrained sizes",
        "[HeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with constrained sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratescaler-test2-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratescaler-test2-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateScaler>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 500000;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing HeightSizeRateScaler with 4 pairs with fixed sizes",
        "[HeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with fixed sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratescaler-test3-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratescaler-test3-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateScaler>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                        REQUIRE(size_leaf1 == 0.01);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double mu_sh;
        double mu_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(0.01));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(0.0));
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing HeightSizeRateScaler with 4 singletons",
        "[HeightSizeRateScaler]") {

    SECTION("Testing 4 singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratescaler-test4-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratescaler-test4-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateScaler>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(tree.get_leaf_node_count() == 1);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
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
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing HeightSizeRateScaler with mix of pairs and singletons",
        "[HeightSizeRateScaler]") {

    SECTION("Testing mix of pairs and singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratescaler-test5-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratescaler-test5-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname4-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateScaler>(1.0, 0.5);
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
        SampleSummarizer<double> mutation_rate_summary_pair1;
        SampleSummarizer<double> mutation_rate_summary_pair2;
        SampleSummarizer<double> mutation_rate_summary_pair3;
        SampleSummarizer<double> mutation_rate_summary_single1;
        SampleSummarizer<double> mutation_rate_summary_single2;

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

                mutation_rate_summary_pair1.add_sample(comparisons.get_tree(0).get_mutation_rate());
                mutation_rate_summary_pair2.add_sample(comparisons.get_tree(1).get_mutation_rate());
                mutation_rate_summary_pair3.add_sample(comparisons.get_tree(2).get_mutation_rate());
                mutation_rate_summary_single1.add_sample(comparisons.get_tree(3).get_mutation_rate());
                mutation_rate_summary_single2.add_sample(comparisons.get_tree(4).get_mutation_rate());

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

        REQUIRE(mutation_rate_summary_pair1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair1.variance() == Approx(0.0));

        double mu_sh = mu_shapes.at(0);
        double mu_sc = mu_scales.at(0);
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_pair2.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(1);
        mu_sc = mu_scales.at(1);
        REQUIRE(mutation_rate_summary_pair3.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(2);
        mu_sc = mu_scales.at(2);
        REQUIRE(mutation_rate_summary_single1.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_single1.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(3);
        mu_sc = mu_scales.at(3);
        REQUIRE(mutation_rate_summary_single2.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_single2.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
    }
}

TEST_CASE("Testing HeightSizeRateScaler with 4 pairs and shared event",
        "[HeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratescaler-test6-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratescaler-test6-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf2_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_root_summaries(ntrees);
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
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
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing HeightSizeRateScaler with mix of pairs and singletons and shared event",
        "[HeightSizeRateScaler]") {

    SECTION("Testing mix of pairs and singletons with shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratescaler-test7-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratescaler-test7-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname4-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 5);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        SampleSummarizer<double> size_root_summary_pair1;
        SampleSummarizer<double> size_root_summary_pair2;
        SampleSummarizer<double> size_root_summary_pair3;
        SampleSummarizer<double> size_root_summary_single1;
        SampleSummarizer<double> size_root_summary_single2;
        SampleSummarizer<double> size_leaf_summary_pair1;
        SampleSummarizer<double> size_leaf_summary_single1;
        SampleSummarizer<double> size_leaf_summary_single2;
        SampleSummarizer<double> mutation_rate_summary_pair1;
        SampleSummarizer<double> mutation_rate_summary_pair2;
        SampleSummarizer<double> mutation_rate_summary_pair3;
        SampleSummarizer<double> mutation_rate_summary_single1;
        SampleSummarizer<double> mutation_rate_summary_single2;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                REQUIRE(comparisons.get_number_of_events() == 1);

                height_summary.add_sample(comparisons.get_tree(0).get_height());

                mutation_rate_summary_pair1.add_sample(comparisons.get_tree(0).get_mutation_rate());
                mutation_rate_summary_pair2.add_sample(comparisons.get_tree(1).get_mutation_rate());
                mutation_rate_summary_pair3.add_sample(comparisons.get_tree(2).get_mutation_rate());
                mutation_rate_summary_single1.add_sample(comparisons.get_tree(3).get_mutation_rate());
                mutation_rate_summary_single2.add_sample(comparisons.get_tree(4).get_mutation_rate());

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
        
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

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

        REQUIRE(mutation_rate_summary_pair1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair1.variance() == Approx(0.0));

        double mu_sh = mu_shapes.at(0);
        double mu_sc = mu_scales.at(0);
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_pair2.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(1);
        mu_sc = mu_scales.at(1);
        REQUIRE(mutation_rate_summary_pair3.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(2);
        mu_sc = mu_scales.at(2);
        REQUIRE(mutation_rate_summary_single1.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_single1.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(3);
        mu_sc = mu_scales.at(3);
        REQUIRE(mutation_rate_summary_single2.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_single2.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
    }
}

TEST_CASE("Testing HeightSizeRateScaler with 4 pairs with constrained sizes and shared event",
        "[HeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with constrained sizes and shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratescaler-test8-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratescaler-test8-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 500000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing HeightSizeRateScaler with 4 pairs and fixed rates",
        "[HeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratescaler-test9-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratescaler-test9-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateScaler>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
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
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing HeightSizeRateScaler with 4 pairs with constrained sizes and fixed rates",
        "[HeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with constrained sizes and fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratescaler-test10-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratescaler-test10-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
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

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateScaler>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 500000;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
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
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing HeightSizeRateScaler with 4 pairs with fixed sizes and rates",
        "[HeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with fixed sizes and rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratescaler-test11-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratescaler-test11-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateScaler>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                        REQUIRE(size_leaf1 == 0.01);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
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

            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing HeightSizeRateScaler with 4 singletons and fixed rates",
        "[HeightSizeRateScaler]") {

    SECTION("Testing 4 singletons and fixed rates with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratescaler-test12-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratescaler-test12-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateScaler>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(tree.get_leaf_node_count() == 1);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
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
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing HeightSizeRateScaler with mix of pairs and singletons and fixed rates",
        "[HeightSizeRateScaler]") {

    SECTION("Testing mix of pairs and singletons with fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratescaler-test13-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratescaler-test13-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateScaler>(1.0, 0.5);
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
        SampleSummarizer<double> mutation_rate_summary_pair1;
        SampleSummarizer<double> mutation_rate_summary_pair2;
        SampleSummarizer<double> mutation_rate_summary_pair3;
        SampleSummarizer<double> mutation_rate_summary_single1;
        SampleSummarizer<double> mutation_rate_summary_single2;

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

                mutation_rate_summary_pair1.add_sample(comparisons.get_tree(0).get_mutation_rate());
                mutation_rate_summary_pair2.add_sample(comparisons.get_tree(1).get_mutation_rate());
                mutation_rate_summary_pair3.add_sample(comparisons.get_tree(2).get_mutation_rate());
                mutation_rate_summary_single1.add_sample(comparisons.get_tree(3).get_mutation_rate());
                mutation_rate_summary_single2.add_sample(comparisons.get_tree(4).get_mutation_rate());

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

        REQUIRE(mutation_rate_summary_pair1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair1.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_single1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_single1.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_single2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_single2.variance() == Approx(0.0));
    }
}

TEST_CASE("Testing HeightSizeRateScaler with 4 pairs and shared event and fixed rates",
        "[HeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with shared event and fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratescaler-test14-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratescaler-test14-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf2_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_root_summaries(ntrees);
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
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
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing HeightSizeRateScaler with mix of pairs and singletons and shared event and fixed rates",
        "[HeightSizeRateScaler]") {

    SECTION("Testing mix of pairs and singletons with shared event and fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratescaler-test15-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratescaler-test15-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 5);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        SampleSummarizer<double> size_root_summary_pair1;
        SampleSummarizer<double> size_root_summary_pair2;
        SampleSummarizer<double> size_root_summary_pair3;
        SampleSummarizer<double> size_root_summary_single1;
        SampleSummarizer<double> size_root_summary_single2;
        SampleSummarizer<double> size_leaf_summary_pair1;
        SampleSummarizer<double> size_leaf_summary_single1;
        SampleSummarizer<double> size_leaf_summary_single2;
        SampleSummarizer<double> mutation_rate_summary_pair1;
        SampleSummarizer<double> mutation_rate_summary_pair2;
        SampleSummarizer<double> mutation_rate_summary_pair3;
        SampleSummarizer<double> mutation_rate_summary_single1;
        SampleSummarizer<double> mutation_rate_summary_single2;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                REQUIRE(comparisons.get_number_of_events() == 1);

                height_summary.add_sample(comparisons.get_tree(0).get_height());

                mutation_rate_summary_pair1.add_sample(comparisons.get_tree(0).get_mutation_rate());
                mutation_rate_summary_pair2.add_sample(comparisons.get_tree(1).get_mutation_rate());
                mutation_rate_summary_pair3.add_sample(comparisons.get_tree(2).get_mutation_rate());
                mutation_rate_summary_single1.add_sample(comparisons.get_tree(3).get_mutation_rate());
                mutation_rate_summary_single2.add_sample(comparisons.get_tree(4).get_mutation_rate());

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
        
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

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

        REQUIRE(mutation_rate_summary_pair1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair1.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair2.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_pair3.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_single1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_single1.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_single2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_single2.variance() == Approx(0.0));
    }
}

TEST_CASE("Testing HeightSizeRateScaler with 4 pairs with constrained sizes and shared event and fixed rates",
        "[HeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with constrained sizes and shared event and fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratescaler-test16-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratescaler-test16-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: true\n";
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

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 500000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateScaler with 4 pairs",
        "[CompositeHeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratescaler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratescaler-test1-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateScaler>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
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
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateScaler with 4 pairs with constrained sizes",
        "[CompositeHeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with constrained sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratescaler-test2-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratescaler-test2-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateScaler>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 500000;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateScaler with 4 pairs with fixed sizes",
        "[CompositeHeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with fixed sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratescaler-test3-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratescaler-test3-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateScaler>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                        REQUIRE(size_leaf1 == 0.01);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double mu_sh;
        double mu_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(0.01));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(0.0));
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateScaler with 4 singletons",
        "[CompositeHeightSizeRateScaler]") {

    SECTION("Testing 4 singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratescaler-test4-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratescaler-test4-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateScaler>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(tree.get_leaf_node_count() == 1);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
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
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateScaler with mix of pairs and singletons",
        "[CompositeHeightSizeRateScaler]") {

    SECTION("Testing mix of pairs and singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratescaler-test5-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratescaler-test5-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname4-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateScaler>(1.0, 0.5);
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
        SampleSummarizer<double> mutation_rate_summary_pair1;
        SampleSummarizer<double> mutation_rate_summary_pair2;
        SampleSummarizer<double> mutation_rate_summary_pair3;
        SampleSummarizer<double> mutation_rate_summary_single1;
        SampleSummarizer<double> mutation_rate_summary_single2;

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

                mutation_rate_summary_pair1.add_sample(comparisons.get_tree(0).get_mutation_rate());
                mutation_rate_summary_pair2.add_sample(comparisons.get_tree(1).get_mutation_rate());
                mutation_rate_summary_pair3.add_sample(comparisons.get_tree(2).get_mutation_rate());
                mutation_rate_summary_single1.add_sample(comparisons.get_tree(3).get_mutation_rate());
                mutation_rate_summary_single2.add_sample(comparisons.get_tree(4).get_mutation_rate());

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

        REQUIRE(mutation_rate_summary_pair1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair1.variance() == Approx(0.0));

        double mu_sh = mu_shapes.at(0);
        double mu_sc = mu_scales.at(0);
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_pair2.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(1);
        mu_sc = mu_scales.at(1);
        REQUIRE(mutation_rate_summary_pair3.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(2);
        mu_sc = mu_scales.at(2);
        REQUIRE(mutation_rate_summary_single1.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_single1.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(3);
        mu_sc = mu_scales.at(3);
        REQUIRE(mutation_rate_summary_single2.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_single2.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
    }
}

TEST_CASE("Testing CompositeHeightSizeRateScaler with 4 pairs and shared event",
        "[CompositeHeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratescaler-test6-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratescaler-test6-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf2_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_root_summaries(ntrees);
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
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
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateScaler with mix of pairs and singletons and shared event",
        "[CompositeHeightSizeRateScaler]") {

    SECTION("Testing mix of pairs and singletons with shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratescaler-test7-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratescaler-test7-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname4-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 5);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        SampleSummarizer<double> size_root_summary_pair1;
        SampleSummarizer<double> size_root_summary_pair2;
        SampleSummarizer<double> size_root_summary_pair3;
        SampleSummarizer<double> size_root_summary_single1;
        SampleSummarizer<double> size_root_summary_single2;
        SampleSummarizer<double> size_leaf_summary_pair1;
        SampleSummarizer<double> size_leaf_summary_single1;
        SampleSummarizer<double> size_leaf_summary_single2;
        SampleSummarizer<double> mutation_rate_summary_pair1;
        SampleSummarizer<double> mutation_rate_summary_pair2;
        SampleSummarizer<double> mutation_rate_summary_pair3;
        SampleSummarizer<double> mutation_rate_summary_single1;
        SampleSummarizer<double> mutation_rate_summary_single2;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                REQUIRE(comparisons.get_number_of_events() == 1);

                height_summary.add_sample(comparisons.get_tree(0).get_height());

                mutation_rate_summary_pair1.add_sample(comparisons.get_tree(0).get_mutation_rate());
                mutation_rate_summary_pair2.add_sample(comparisons.get_tree(1).get_mutation_rate());
                mutation_rate_summary_pair3.add_sample(comparisons.get_tree(2).get_mutation_rate());
                mutation_rate_summary_single1.add_sample(comparisons.get_tree(3).get_mutation_rate());
                mutation_rate_summary_single2.add_sample(comparisons.get_tree(4).get_mutation_rate());

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
        
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

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

        REQUIRE(mutation_rate_summary_pair1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair1.variance() == Approx(0.0));

        double mu_sh = mu_shapes.at(0);
        double mu_sc = mu_scales.at(0);
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_pair2.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(1);
        mu_sc = mu_scales.at(1);
        REQUIRE(mutation_rate_summary_pair3.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(2);
        mu_sc = mu_scales.at(2);
        REQUIRE(mutation_rate_summary_single1.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_single1.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(3);
        mu_sc = mu_scales.at(3);
        REQUIRE(mutation_rate_summary_single2.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_single2.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
    }
}

TEST_CASE("Testing CompositeHeightSizeRateScaler with 4 pairs with constrained sizes and shared event",
        "[CompositeHeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with constrained sizes and shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratescaler-test8-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratescaler-test8-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 500000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateScaler with 4 pairs and fixed rates",
        "[CompositeHeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratescaler-test9-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratescaler-test9-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateScaler>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
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
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateScaler with 4 pairs with constrained sizes and fixed rates",
        "[CompositeHeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with constrained sizes and fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratescaler-test10-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratescaler-test10-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
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

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateScaler>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 500000;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
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
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateScaler with 4 pairs with fixed sizes and rates",
        "[CompositeHeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with fixed sizes and rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratescaler-test11-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratescaler-test11-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateScaler>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                        REQUIRE(size_leaf1 == 0.01);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
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

            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateScaler with 4 singletons and fixed rates",
        "[CompositeHeightSizeRateScaler]") {

    SECTION("Testing 4 singletons and fixed rates with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratescaler-test12-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratescaler-test12-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateScaler>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(tree.get_leaf_node_count() == 1);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
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
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateScaler with mix of pairs and singletons and fixed rates",
        "[CompositeHeightSizeRateScaler]") {

    SECTION("Testing mix of pairs and singletons with fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratescaler-test13-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratescaler-test13-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateScaler>(1.0, 0.5);
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
        SampleSummarizer<double> mutation_rate_summary_pair1;
        SampleSummarizer<double> mutation_rate_summary_pair2;
        SampleSummarizer<double> mutation_rate_summary_pair3;
        SampleSummarizer<double> mutation_rate_summary_single1;
        SampleSummarizer<double> mutation_rate_summary_single2;

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

                mutation_rate_summary_pair1.add_sample(comparisons.get_tree(0).get_mutation_rate());
                mutation_rate_summary_pair2.add_sample(comparisons.get_tree(1).get_mutation_rate());
                mutation_rate_summary_pair3.add_sample(comparisons.get_tree(2).get_mutation_rate());
                mutation_rate_summary_single1.add_sample(comparisons.get_tree(3).get_mutation_rate());
                mutation_rate_summary_single2.add_sample(comparisons.get_tree(4).get_mutation_rate());

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

        REQUIRE(mutation_rate_summary_pair1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair1.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_single1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_single1.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_single2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_single2.variance() == Approx(0.0));
    }
}

TEST_CASE("Testing CompositeHeightSizeRateScaler with 4 pairs and shared event and fixed rates",
        "[CompositeHeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with shared event and fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratescaler-test14-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratescaler-test14-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf2_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_root_summaries(ntrees);
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
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
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateScaler with mix of pairs and singletons and shared event and fixed rates",
        "[CompositeHeightSizeRateScaler]") {

    SECTION("Testing mix of pairs and singletons with shared event and fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratescaler-test15-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratescaler-test15-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 5);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        SampleSummarizer<double> size_root_summary_pair1;
        SampleSummarizer<double> size_root_summary_pair2;
        SampleSummarizer<double> size_root_summary_pair3;
        SampleSummarizer<double> size_root_summary_single1;
        SampleSummarizer<double> size_root_summary_single2;
        SampleSummarizer<double> size_leaf_summary_pair1;
        SampleSummarizer<double> size_leaf_summary_single1;
        SampleSummarizer<double> size_leaf_summary_single2;
        SampleSummarizer<double> mutation_rate_summary_pair1;
        SampleSummarizer<double> mutation_rate_summary_pair2;
        SampleSummarizer<double> mutation_rate_summary_pair3;
        SampleSummarizer<double> mutation_rate_summary_single1;
        SampleSummarizer<double> mutation_rate_summary_single2;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                REQUIRE(comparisons.get_number_of_events() == 1);

                height_summary.add_sample(comparisons.get_tree(0).get_height());

                mutation_rate_summary_pair1.add_sample(comparisons.get_tree(0).get_mutation_rate());
                mutation_rate_summary_pair2.add_sample(comparisons.get_tree(1).get_mutation_rate());
                mutation_rate_summary_pair3.add_sample(comparisons.get_tree(2).get_mutation_rate());
                mutation_rate_summary_single1.add_sample(comparisons.get_tree(3).get_mutation_rate());
                mutation_rate_summary_single2.add_sample(comparisons.get_tree(4).get_mutation_rate());

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
        
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

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

        REQUIRE(mutation_rate_summary_pair1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair1.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair2.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_pair3.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_single1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_single1.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_single2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_single2.variance() == Approx(0.0));
    }
}

TEST_CASE("Testing CompositeHeightSizeRateScaler with 4 pairs with constrained sizes and shared event and fixed rates",
        "[CompositeHeightSizeRateScaler]") {

    SECTION("Testing 4 pairs with constrained sizes and shared event and fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratescaler-test16-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratescaler-test16-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: true\n";
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

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateScaler>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 500000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}


TEST_CASE("Testing HeightSizeRateMixer with 4 pairs",
        "[HeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratemixer-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratemixer-test1-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateMixer>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
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
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing HeightSizeRateMixer with 4 pairs with constrained sizes",
        "[HeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with constrained sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratemixer-test2-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratemixer-test2-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateMixer>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 500000;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing HeightSizeRateMixer with 4 pairs with fixed sizes",
        "[HeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with fixed sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratemixer-test3-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratemixer-test3-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateMixer>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                        REQUIRE(size_leaf1 == 0.01);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double mu_sh;
        double mu_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(0.01));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(0.0));
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing HeightSizeRateMixer with 4 singletons",
        "[HeightSizeRateMixer]") {

    SECTION("Testing 4 singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratemixer-test4-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratemixer-test4-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateMixer>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 900000;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(tree.get_leaf_node_count() == 1);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
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
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing HeightSizeRateMixer with mix of pairs and singletons",
        "[HeightSizeRateMixer]") {

    SECTION("Testing mix of pairs and singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratemixer-test5-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratemixer-test5-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname4-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateMixer>(1.0, 0.5);
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
        SampleSummarizer<double> mutation_rate_summary_pair1;
        SampleSummarizer<double> mutation_rate_summary_pair2;
        SampleSummarizer<double> mutation_rate_summary_pair3;
        SampleSummarizer<double> mutation_rate_summary_single1;
        SampleSummarizer<double> mutation_rate_summary_single2;

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

                mutation_rate_summary_pair1.add_sample(comparisons.get_tree(0).get_mutation_rate());
                mutation_rate_summary_pair2.add_sample(comparisons.get_tree(1).get_mutation_rate());
                mutation_rate_summary_pair3.add_sample(comparisons.get_tree(2).get_mutation_rate());
                mutation_rate_summary_single1.add_sample(comparisons.get_tree(3).get_mutation_rate());
                mutation_rate_summary_single2.add_sample(comparisons.get_tree(4).get_mutation_rate());

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

        REQUIRE(mutation_rate_summary_pair1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair1.variance() == Approx(0.0));

        double mu_sh = mu_shapes.at(0);
        double mu_sc = mu_scales.at(0);
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_pair2.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(1);
        mu_sc = mu_scales.at(1);
        REQUIRE(mutation_rate_summary_pair3.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(2);
        mu_sc = mu_scales.at(2);
        REQUIRE(mutation_rate_summary_single1.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_single1.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(3);
        mu_sc = mu_scales.at(3);
        REQUIRE(mutation_rate_summary_single2.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_single2.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
    }
}

TEST_CASE("Testing HeightSizeRateMixer with 4 pairs and shared event",
        "[HeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratemixer-test6-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratemixer-test6-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf2_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_root_summaries(ntrees);
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 900000;
        unsigned int sample_freq = 3;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
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
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing HeightSizeRateMixer with mix of pairs and singletons and shared event",
        "[HeightSizeRateMixer]") {

    SECTION("Testing mix of pairs and singletons with shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratemixer-test7-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratemixer-test7-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname4-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 5);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        SampleSummarizer<double> size_root_summary_pair1;
        SampleSummarizer<double> size_root_summary_pair2;
        SampleSummarizer<double> size_root_summary_pair3;
        SampleSummarizer<double> size_root_summary_single1;
        SampleSummarizer<double> size_root_summary_single2;
        SampleSummarizer<double> size_leaf_summary_pair1;
        SampleSummarizer<double> size_leaf_summary_single1;
        SampleSummarizer<double> size_leaf_summary_single2;
        SampleSummarizer<double> mutation_rate_summary_pair1;
        SampleSummarizer<double> mutation_rate_summary_pair2;
        SampleSummarizer<double> mutation_rate_summary_pair3;
        SampleSummarizer<double> mutation_rate_summary_single1;
        SampleSummarizer<double> mutation_rate_summary_single2;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                REQUIRE(comparisons.get_number_of_events() == 1);

                height_summary.add_sample(comparisons.get_tree(0).get_height());

                mutation_rate_summary_pair1.add_sample(comparisons.get_tree(0).get_mutation_rate());
                mutation_rate_summary_pair2.add_sample(comparisons.get_tree(1).get_mutation_rate());
                mutation_rate_summary_pair3.add_sample(comparisons.get_tree(2).get_mutation_rate());
                mutation_rate_summary_single1.add_sample(comparisons.get_tree(3).get_mutation_rate());
                mutation_rate_summary_single2.add_sample(comparisons.get_tree(4).get_mutation_rate());

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
        
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

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

        REQUIRE(mutation_rate_summary_pair1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair1.variance() == Approx(0.0));

        double mu_sh = mu_shapes.at(0);
        double mu_sc = mu_scales.at(0);
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_pair2.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(1);
        mu_sc = mu_scales.at(1);
        REQUIRE(mutation_rate_summary_pair3.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(2);
        mu_sc = mu_scales.at(2);
        REQUIRE(mutation_rate_summary_single1.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_single1.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(3);
        mu_sc = mu_scales.at(3);
        REQUIRE(mutation_rate_summary_single2.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_single2.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
    }
}

TEST_CASE("Testing HeightSizeRateMixer with 4 pairs with constrained sizes and shared event",
        "[HeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with constrained sizes and shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratemixer-test8-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratemixer-test8-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 500000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing HeightSizeRateMixer with 4 pairs and fixed rates",
        "[HeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratemixer-test9-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratemixer-test9-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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

        RandomNumberGenerator rng = RandomNumberGenerator(1234);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateMixer>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 1200000;
        unsigned int sample_freq = 3;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
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
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing HeightSizeRateMixer with 4 pairs with constrained sizes and fixed rates",
        "[HeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with constrained sizes and fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratemixer-test10-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratemixer-test10-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
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

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateMixer>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 500000;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
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
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing HeightSizeRateMixer with 4 pairs with fixed sizes and rates",
        "[HeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with fixed sizes and rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratemixer-test11-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratemixer-test11-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateMixer>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                        REQUIRE(size_leaf1 == 0.01);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
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

            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing HeightSizeRateMixer with 4 singletons and fixed rates",
        "[HeightSizeRateMixer]") {

    SECTION("Testing 4 singletons and fixed rates with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratemixer-test12-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratemixer-test12-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateMixer>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 900000;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(tree.get_leaf_node_count() == 1);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
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
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing HeightSizeRateMixer with mix of pairs and singletons and fixed rates",
        "[HeightSizeRateMixer]") {

    SECTION("Testing mix of pairs and singletons with fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratemixer-test13-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratemixer-test13-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateMixer>(1.0, 0.5);
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
        SampleSummarizer<double> mutation_rate_summary_pair1;
        SampleSummarizer<double> mutation_rate_summary_pair2;
        SampleSummarizer<double> mutation_rate_summary_pair3;
        SampleSummarizer<double> mutation_rate_summary_single1;
        SampleSummarizer<double> mutation_rate_summary_single2;

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

                mutation_rate_summary_pair1.add_sample(comparisons.get_tree(0).get_mutation_rate());
                mutation_rate_summary_pair2.add_sample(comparisons.get_tree(1).get_mutation_rate());
                mutation_rate_summary_pair3.add_sample(comparisons.get_tree(2).get_mutation_rate());
                mutation_rate_summary_single1.add_sample(comparisons.get_tree(3).get_mutation_rate());
                mutation_rate_summary_single2.add_sample(comparisons.get_tree(4).get_mutation_rate());

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

        REQUIRE(mutation_rate_summary_pair1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair1.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_single1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_single1.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_single2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_single2.variance() == Approx(0.0));
    }
}

TEST_CASE("Testing HeightSizeRateMixer with 4 pairs and shared event and fixed rates",
        "[HeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with shared event and fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratemixer-test14-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratemixer-test14-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
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

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf2_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_root_summaries(ntrees);
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 900000;
        unsigned int sample_freq = 3;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
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
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing HeightSizeRateMixer with mix of pairs and singletons and shared event and fixed rates",
        "[HeightSizeRateMixer]") {

    SECTION("Testing mix of pairs and singletons with shared event and fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratemixer-test15-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratemixer-test15-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 5);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        SampleSummarizer<double> size_root_summary_pair1;
        SampleSummarizer<double> size_root_summary_pair2;
        SampleSummarizer<double> size_root_summary_pair3;
        SampleSummarizer<double> size_root_summary_single1;
        SampleSummarizer<double> size_root_summary_single2;
        SampleSummarizer<double> size_leaf_summary_pair1;
        SampleSummarizer<double> size_leaf_summary_single1;
        SampleSummarizer<double> size_leaf_summary_single2;
        SampleSummarizer<double> mutation_rate_summary_pair1;
        SampleSummarizer<double> mutation_rate_summary_pair2;
        SampleSummarizer<double> mutation_rate_summary_pair3;
        SampleSummarizer<double> mutation_rate_summary_single1;
        SampleSummarizer<double> mutation_rate_summary_single2;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                REQUIRE(comparisons.get_number_of_events() == 1);

                height_summary.add_sample(comparisons.get_tree(0).get_height());

                mutation_rate_summary_pair1.add_sample(comparisons.get_tree(0).get_mutation_rate());
                mutation_rate_summary_pair2.add_sample(comparisons.get_tree(1).get_mutation_rate());
                mutation_rate_summary_pair3.add_sample(comparisons.get_tree(2).get_mutation_rate());
                mutation_rate_summary_single1.add_sample(comparisons.get_tree(3).get_mutation_rate());
                mutation_rate_summary_single2.add_sample(comparisons.get_tree(4).get_mutation_rate());

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
        
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

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

        REQUIRE(mutation_rate_summary_pair1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair1.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair2.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_pair3.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_single1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_single1.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_single2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_single2.variance() == Approx(0.0));
    }
}

TEST_CASE("Testing HeightSizeRateMixer with 4 pairs with constrained sizes and shared event and fixed rates",
        "[HeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with constrained sizes and shared event and fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-heightsizeratemixer-test16-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-heightsizeratemixer-test16-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: true\n";
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

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<HeightSizeRateMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 500000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateMixer with 4 pairs",
        "[CompositeHeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratemixer-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratemixer-test1-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateMixer>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
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
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateMixer with 4 pairs with constrained sizes",
        "[CompositeHeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with constrained sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratemixer-test2-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratemixer-test2-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateMixer>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 500000;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateMixer with 4 pairs with fixed sizes",
        "[CompositeHeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with fixed sizes and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratemixer-test3-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratemixer-test3-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateMixer>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                        REQUIRE(size_leaf1 == 0.01);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double mu_sh;
        double mu_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(0.01));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(0.0));
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateMixer with 4 singletons",
        "[CompositeHeightSizeRateMixer]") {

    SECTION("Testing 4 singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratemixer-test4-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratemixer-test4-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateMixer>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(tree.get_leaf_node_count() == 1);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
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
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateMixer with mix of pairs and singletons",
        "[CompositeHeightSizeRateMixer]") {

    SECTION("Testing mix of pairs and singletons with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratemixer-test5-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratemixer-test5-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname4-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateMixer>(1.0, 0.5);
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
        SampleSummarizer<double> mutation_rate_summary_pair1;
        SampleSummarizer<double> mutation_rate_summary_pair2;
        SampleSummarizer<double> mutation_rate_summary_pair3;
        SampleSummarizer<double> mutation_rate_summary_single1;
        SampleSummarizer<double> mutation_rate_summary_single2;

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

                mutation_rate_summary_pair1.add_sample(comparisons.get_tree(0).get_mutation_rate());
                mutation_rate_summary_pair2.add_sample(comparisons.get_tree(1).get_mutation_rate());
                mutation_rate_summary_pair3.add_sample(comparisons.get_tree(2).get_mutation_rate());
                mutation_rate_summary_single1.add_sample(comparisons.get_tree(3).get_mutation_rate());
                mutation_rate_summary_single2.add_sample(comparisons.get_tree(4).get_mutation_rate());

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

        REQUIRE(mutation_rate_summary_pair1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair1.variance() == Approx(0.0));

        double mu_sh = mu_shapes.at(0);
        double mu_sc = mu_scales.at(0);
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_pair2.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(1);
        mu_sc = mu_scales.at(1);
        REQUIRE(mutation_rate_summary_pair3.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(2);
        mu_sc = mu_scales.at(2);
        REQUIRE(mutation_rate_summary_single1.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_single1.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(3);
        mu_sc = mu_scales.at(3);
        REQUIRE(mutation_rate_summary_single2.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_single2.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
    }
}

TEST_CASE("Testing CompositeHeightSizeRateMixer with 4 pairs and shared event",
        "[CompositeHeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratemixer-test6-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratemixer-test6-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf2_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_root_summaries(ntrees);
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
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
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateMixer with mix of pairs and singletons and shared event",
        "[CompositeHeightSizeRateMixer]") {

    SECTION("Testing mix of pairs and singletons with shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratemixer-test7-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratemixer-test7-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.01\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname4-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 5);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        SampleSummarizer<double> size_root_summary_pair1;
        SampleSummarizer<double> size_root_summary_pair2;
        SampleSummarizer<double> size_root_summary_pair3;
        SampleSummarizer<double> size_root_summary_single1;
        SampleSummarizer<double> size_root_summary_single2;
        SampleSummarizer<double> size_leaf_summary_pair1;
        SampleSummarizer<double> size_leaf_summary_single1;
        SampleSummarizer<double> size_leaf_summary_single2;
        SampleSummarizer<double> mutation_rate_summary_pair1;
        SampleSummarizer<double> mutation_rate_summary_pair2;
        SampleSummarizer<double> mutation_rate_summary_pair3;
        SampleSummarizer<double> mutation_rate_summary_single1;
        SampleSummarizer<double> mutation_rate_summary_single2;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                REQUIRE(comparisons.get_number_of_events() == 1);

                height_summary.add_sample(comparisons.get_tree(0).get_height());

                mutation_rate_summary_pair1.add_sample(comparisons.get_tree(0).get_mutation_rate());
                mutation_rate_summary_pair2.add_sample(comparisons.get_tree(1).get_mutation_rate());
                mutation_rate_summary_pair3.add_sample(comparisons.get_tree(2).get_mutation_rate());
                mutation_rate_summary_single1.add_sample(comparisons.get_tree(3).get_mutation_rate());
                mutation_rate_summary_single2.add_sample(comparisons.get_tree(4).get_mutation_rate());

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
        
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

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

        REQUIRE(mutation_rate_summary_pair1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair1.variance() == Approx(0.0));

        double mu_sh = mu_shapes.at(0);
        double mu_sc = mu_scales.at(0);
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_pair2.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(1);
        mu_sc = mu_scales.at(1);
        REQUIRE(mutation_rate_summary_pair3.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(2);
        mu_sc = mu_scales.at(2);
        REQUIRE(mutation_rate_summary_single1.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_single1.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));

        mu_sh = mu_shapes.at(3);
        mu_sc = mu_scales.at(3);
        REQUIRE(mutation_rate_summary_single2.mean() == Approx(mu_sh * mu_sc).epsilon(0.01));
        REQUIRE(mutation_rate_summary_single2.variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
    }
}

TEST_CASE("Testing CompositeHeightSizeRateMixer with 4 pairs with constrained sizes and shared event",
        "[CompositeHeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with constrained sizes and shared event and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratemixer-test8-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratemixer-test8-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
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
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(0) << "\n";
        os << "                    scale: " << mu_scales.at(0) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(1) << "\n";
        os << "                    scale: " << size_scales.at(1) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(1) << "\n";
        os << "                    scale: " << mu_scales.at(1) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(2) << "\n";
        os << "                    scale: " << size_scales.at(2) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(2) << "\n";
        os << "                    scale: " << mu_scales.at(2) << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname3.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size_shapes.at(3) << "\n";
        os << "                    scale: " << size_scales.at(3) << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mu_shapes.at(3) << "\n";
        os << "                    scale: " << mu_scales.at(3) << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 500000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        double mu_sh;
        double mu_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            mu_sh = mu_shapes.at(tree_idx);
            mu_sc = mu_scales.at(tree_idx);
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(mu_sh * mu_sc).epsilon(0.005));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(mu_sh * mu_sc * mu_sc).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateMixer with 4 pairs and fixed rates",
        "[CompositeHeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratemixer-test9-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratemixer-test9-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateMixer>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
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
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateMixer with 4 pairs with constrained sizes and fixed rates",
        "[CompositeHeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with constrained sizes and fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratemixer-test10-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratemixer-test10-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
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

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateMixer>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 500000;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
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
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateMixer with 4 pairs with fixed sizes and rates",
        "[CompositeHeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with fixed sizes and rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> mu_shapes {5.0, 3.0, 4.0, 10.0};
        std::vector<double> mu_scales {0.15, 0.2, 0.25, 0.05};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratemixer-test11-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratemixer-test11-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: true\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateMixer>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                        REQUIRE(size_leaf1 == 0.01);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
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

            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateMixer with 4 singletons and fixed rates",
        "[CompositeHeightSizeRateMixer]") {

    SECTION("Testing 4 singletons and fixed rates with optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratemixer-test12-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratemixer-test12-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateMixer>(1.0, 0.5);
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
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_root;
        double mutation_rate;
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
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(tree.get_leaf_node_count() == 1);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
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
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateMixer with mix of pairs and singletons and fixed rates",
        "[CompositeHeightSizeRateMixer]") {

    SECTION("Testing mix of pairs and singletons with fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratemixer-test13-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratemixer-test13-" + tag + "-state-run-1.log";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateMixer>(1.0, 0.5);
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
        SampleSummarizer<double> mutation_rate_summary_pair1;
        SampleSummarizer<double> mutation_rate_summary_pair2;
        SampleSummarizer<double> mutation_rate_summary_pair3;
        SampleSummarizer<double> mutation_rate_summary_single1;
        SampleSummarizer<double> mutation_rate_summary_single2;

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

                mutation_rate_summary_pair1.add_sample(comparisons.get_tree(0).get_mutation_rate());
                mutation_rate_summary_pair2.add_sample(comparisons.get_tree(1).get_mutation_rate());
                mutation_rate_summary_pair3.add_sample(comparisons.get_tree(2).get_mutation_rate());
                mutation_rate_summary_single1.add_sample(comparisons.get_tree(3).get_mutation_rate());
                mutation_rate_summary_single2.add_sample(comparisons.get_tree(4).get_mutation_rate());

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

        REQUIRE(mutation_rate_summary_pair1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair1.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_single1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_single1.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_single2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_single2.variance() == Approx(0.0));
    }
}

TEST_CASE("Testing CompositeHeightSizeRateMixer with 4 pairs and shared event and fixed rates",
        "[CompositeHeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with shared event and fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratemixer-test14-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratemixer-test14-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_leaf2_summaries(ntrees);
        std::vector< SampleSummarizer<double> > size_root_summaries(ntrees);
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 != size_leaf2);
                        REQUIRE(size_leaf1 != size_root);
                        REQUIRE(size_leaf2 != size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    size_leaf2_summaries.at(tree_idx).add_sample(size_leaf2);
                    size_root_summaries.at(tree_idx).add_sample(size_root);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
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
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
    }
}

TEST_CASE("Testing CompositeHeightSizeRateMixer with mix of pairs and singletons and shared event and fixed rates",
        "[CompositeHeightSizeRateMixer]") {

    SECTION("Testing mix of pairs and singletons with shared event and fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratemixer-test15-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratemixer-test15-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: false\n";
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
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 5);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        SampleSummarizer<double> size_root_summary_pair1;
        SampleSummarizer<double> size_root_summary_pair2;
        SampleSummarizer<double> size_root_summary_pair3;
        SampleSummarizer<double> size_root_summary_single1;
        SampleSummarizer<double> size_root_summary_single2;
        SampleSummarizer<double> size_leaf_summary_pair1;
        SampleSummarizer<double> size_leaf_summary_single1;
        SampleSummarizer<double> size_leaf_summary_single2;
        SampleSummarizer<double> mutation_rate_summary_pair1;
        SampleSummarizer<double> mutation_rate_summary_pair2;
        SampleSummarizer<double> mutation_rate_summary_pair3;
        SampleSummarizer<double> mutation_rate_summary_single1;
        SampleSummarizer<double> mutation_rate_summary_single2;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 4;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                REQUIRE(comparisons.get_number_of_events() == 1);

                height_summary.add_sample(comparisons.get_tree(0).get_height());

                mutation_rate_summary_pair1.add_sample(comparisons.get_tree(0).get_mutation_rate());
                mutation_rate_summary_pair2.add_sample(comparisons.get_tree(1).get_mutation_rate());
                mutation_rate_summary_pair3.add_sample(comparisons.get_tree(2).get_mutation_rate());
                mutation_rate_summary_single1.add_sample(comparisons.get_tree(3).get_mutation_rate());
                mutation_rate_summary_single2.add_sample(comparisons.get_tree(4).get_mutation_rate());

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
        
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.01));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));

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

        REQUIRE(mutation_rate_summary_pair1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair1.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_pair2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair2.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_pair3.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_pair3.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_single1.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_single1.variance() == Approx(0.0));
        REQUIRE(mutation_rate_summary_single2.mean() == Approx(1.0));
        REQUIRE(mutation_rate_summary_single2.variance() == Approx(0.0));
    }
}

TEST_CASE("Testing CompositeHeightSizeRateMixer with 4 pairs with constrained sizes and shared event and fixed rates",
        "[CompositeHeightSizeRateMixer]") {

    SECTION("Testing 4 pairs with constrained sizes and shared event and fixed rates and optimizing") {
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::vector<double> size_shapes {10.0, 2.0, 5.0, 4.0};
        std::vector<double> size_scales {0.1, 0.2, 0.2, 0.5};
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-compheightsizeratemixer-test16-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-compheightsizeratemixer-test16-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_model_prior:\n";
        os << "    fixed: [0, 0, 0, 0]\n";
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
        os << "    equal_population_sizes: true\n";
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

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<CompositeHeightSizeRateMixer>(1.0, 0.5);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        REQUIRE(comparisons.get_number_of_events() == 1);
        SampleSummarizer<double> height_summary;
        std::vector< SampleSummarizer<double> > size_leaf1_summaries(ntrees);
        std::vector< SampleSummarizer<double> > mutation_rate_summaries(ntrees);

        double size_leaf1;
        double size_leaf2;
        double size_root;
        double mutation_rate;
        unsigned int niterations = 500000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                height_summary.add_sample(comparisons.get_tree(0).get_height());
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    REQUIRE(comparisons.get_number_of_events() == 1);
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    size_leaf1 = tree.get_child_population_size(0);
                    size_leaf2 = tree.get_child_population_size(1);
                    size_root = tree.get_root_population_size();
                    mutation_rate = tree.get_mutation_rate();
                    if (i > (niterations - (sample_freq * 5))) {
                        REQUIRE(size_leaf1 == size_leaf2);
                        REQUIRE(size_leaf1 == size_root);
                    }
                    size_leaf1_summaries.at(tree_idx).add_sample(size_leaf1);
                    mutation_rate_summaries.at(tree_idx).add_sample(mutation_rate);
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(height_summary.sample_size() == nsamples);
        REQUIRE(height_summary.mean() == Approx(height_shape * height_scale).epsilon(0.005));
        REQUIRE(height_summary.variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            size_sh = size_shapes.at(tree_idx);
            size_sc = size_scales.at(tree_idx);
            REQUIRE(size_leaf1_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(size_leaf1_summaries.at(tree_idx).mean() == Approx(size_sh * size_sc).epsilon(0.005));
            REQUIRE(size_leaf1_summaries.at(tree_idx).variance() == Approx(size_sh * size_sc * size_sc).epsilon(0.01));
            REQUIRE(mutation_rate_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(mutation_rate_summaries.at(tree_idx).mean() == Approx(1.0));
            REQUIRE(mutation_rate_summaries.at(tree_idx).variance() == Approx(0.0));
        }
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
        out << "    equal_population_sizes: true\n";
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
        out << "    equal_population_sizes: true\n";
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
        out << "    equal_population_sizes: true\n";
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
        out << "    equal_population_sizes: true\n";
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

        double shape = 10.0;
        double scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream out;
        out.open(test_path);
        out << "global_comparison_settings:\n";
        out << "    genotypes_are_diploid: true\n";
        out << "    markers_are_dominant: false\n";
        out << "    population_name_delimiter: \" \"\n";
        out << "    population_name_is_prefix: true\n";
        out << "    constant_sites_removed: true\n";
        out << "    equal_population_sizes: false\n";
        out << "    parameters:\n";
        out << "        population_size:\n";
        out << "            value: 0.005\n";
        out << "            estimate: true\n";
        out << "            prior:\n";
        out << "                gamma_distribution:\n";
        out << "                    shape: " << shape << "\n";
        out << "                    scale: " << scale << "\n";
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

        RandomNumberGenerator rng = RandomNumberGenerator(12933);
        std::shared_ptr<OperatorInterface> op = std::make_shared<RootPopulationSizeScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(shape, scale);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(os);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            o.perform_collection_move(rng, comparisons, 1);
            double x = comparisons.get_tree(0).get_root_population_size();
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
        double shape = 10.0;
        double scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream out;
        out.open(test_path);
        out << "global_comparison_settings:\n";
        out << "    genotypes_are_diploid: true\n";
        out << "    markers_are_dominant: false\n";
        out << "    population_name_delimiter: \" \"\n";
        out << "    population_name_is_prefix: true\n";
        out << "    constant_sites_removed: true\n";
        out << "    equal_population_sizes: false\n";
        out << "    parameters:\n";
        out << "        population_size:\n";
        out << "            value: 0.005\n";
        out << "            estimate: true\n";
        out << "            prior:\n";
        out << "                gamma_distribution:\n";
        out << "                    shape: " << shape << "\n";
        out << "                    scale: " << scale << "\n";
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

        RandomNumberGenerator rng = RandomNumberGenerator(12933);
        std::shared_ptr<OperatorInterface> op = std::make_shared<RootPopulationSizeScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        //os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(shape, scale);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(os);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            o.perform_collection_move(rng, comparisons, 1);
            double x = comparisons.get_tree(0).get_root_population_size();
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
        double shape = 10.0;
        double scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream out;
        out.open(test_path);
        out << "global_comparison_settings:\n";
        out << "    genotypes_are_diploid: true\n";
        out << "    markers_are_dominant: false\n";
        out << "    population_name_delimiter: \" \"\n";
        out << "    population_name_is_prefix: true\n";
        out << "    constant_sites_removed: true\n";
        out << "    equal_population_sizes: false\n";
        out << "    parameters:\n";
        out << "        population_size:\n";
        out << "            value: 0.005\n";
        out << "            estimate: true\n";
        out << "            prior:\n";
        out << "                gamma_distribution:\n";
        out << "                    shape: " << shape << "\n";
        out << "                    scale: " << scale << "\n";
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

        RandomNumberGenerator rng = RandomNumberGenerator(12933);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ChildPopulationSizeScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(shape, scale);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(os);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            o.perform_collection_move(rng, comparisons, 1);
            double x = comparisons.get_tree(0).get_child_population_size(0);
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
        double shape = 10.0;
        double scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream out;
        out.open(test_path);
        out << "global_comparison_settings:\n";
        out << "    genotypes_are_diploid: true\n";
        out << "    markers_are_dominant: false\n";
        out << "    population_name_delimiter: \" \"\n";
        out << "    population_name_is_prefix: true\n";
        out << "    constant_sites_removed: true\n";
        out << "    equal_population_sizes: false\n";
        out << "    parameters:\n";
        out << "        population_size:\n";
        out << "            value: 0.005\n";
        out << "            estimate: true\n";
        out << "            prior:\n";
        out << "                gamma_distribution:\n";
        out << "                    shape: " << shape << "\n";
        out << "                    scale: " << scale << "\n";
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

        RandomNumberGenerator rng = RandomNumberGenerator(12933);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ChildPopulationSizeScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        // os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(shape, scale);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(os);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            o.perform_collection_move(rng, comparisons, 1);
            double x = comparisons.get_tree(0).get_child_population_size(0);
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
        double shape = 10.0;
        double scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream out;
        out.open(test_path);
        out << "global_comparison_settings:\n";
        out << "    genotypes_are_diploid: true\n";
        out << "    markers_are_dominant: false\n";
        out << "    population_name_delimiter: \" \"\n";
        out << "    population_name_is_prefix: true\n";
        out << "    constant_sites_removed: true\n";
        out << "    equal_population_sizes: true\n";
        out << "    parameters:\n";
        out << "        population_size:\n";
        out << "            value: 0.005\n";
        out << "            estimate: false\n";
        out << "        freq_1:\n";
        out << "            value: 0.5\n";
        out << "            estimate: false\n";
        out << "        mutation_rate:\n";
        out << "            value: 1.0\n";
        out << "            estimate: true\n";
        out << "            prior:\n";
        out << "                gamma_distribution:\n";
        out << "                    shape: " << shape << "\n";
        out << "                    scale: " << scale << "\n";
        out << "comparisons:\n";
        out << "- comparison:\n";
        out << "    path: hemi129.nex\n";
        out.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(928374);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ComparisonMutationRateScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(shape, scale);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(os);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            o.perform_collection_move(rng, comparisons, 1);
            double x = comparisons.get_tree(0).get_mutation_rate();
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
        double shape = 10.0;
        double scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream out;
        out.open(test_path);
        out << "global_comparison_settings:\n";
        out << "    genotypes_are_diploid: true\n";
        out << "    markers_are_dominant: false\n";
        out << "    population_name_delimiter: \" \"\n";
        out << "    population_name_is_prefix: true\n";
        out << "    constant_sites_removed: true\n";
        out << "    equal_population_sizes: true\n";
        out << "    parameters:\n";
        out << "        population_size:\n";
        out << "            value: 0.005\n";
        out << "            estimate: false\n";
        out << "        freq_1:\n";
        out << "            value: 0.5\n";
        out << "            estimate: false\n";
        out << "        mutation_rate:\n";
        out << "            value: 1.0\n";
        out << "            estimate: true\n";
        out << "            prior:\n";
        out << "                gamma_distribution:\n";
        out << "                    shape: " << shape << "\n";
        out << "                    scale: " << scale << "\n";
        out << "comparisons:\n";
        out << "- comparison:\n";
        out << "    path: hemi129.nex\n";
        out.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(928374);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ComparisonMutationRateScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        // os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(shape, scale);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(os);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            o.perform_collection_move(rng, comparisons, 1);
            double x = comparisons.get_tree(0).get_mutation_rate();
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
        double alpha = 1.0;
        double beta = 1.0;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream out;
        out.open(test_path);
        out << "global_comparison_settings:\n";
        out << "    genotypes_are_diploid: true\n";
        out << "    markers_are_dominant: false\n";
        out << "    population_name_delimiter: \" \"\n";
        out << "    population_name_is_prefix: true\n";
        out << "    constant_sites_removed: true\n";
        out << "    equal_population_sizes: true\n";
        out << "    parameters:\n";
        out << "        population_size:\n";
        out << "            value: 0.005\n";
        out << "            estimate: false\n";
        out << "        freq_1:\n";
        out << "            value: 0.5\n";
        out << "            estimate: true\n";
        out << "            prior:\n";
        out << "                beta_distribution:\n";
        out << "                    alpha: " << alpha << "\n";
        out << "                    beta: " << beta << "\n";
        out << "        mutation_rate:\n";
        out << "            value: 1.0\n";
        out << "            estimate: false\n";
        out << "comparisons:\n";
        out << "- comparison:\n";
        out << "    path: hemi129.nex\n";
        out.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(3648);
        std::shared_ptr<OperatorInterface> op = std::make_shared<FreqMover>(1.0, 0.1);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<BetaDistribution>(alpha, beta);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(os);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            o.perform_collection_move(rng, comparisons, 1);
            double x = comparisons.get_tree(0).get_freq_1();
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
        double alpha = 1.0;
        double beta = 1.0;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream out;
        out.open(test_path);
        out << "global_comparison_settings:\n";
        out << "    genotypes_are_diploid: true\n";
        out << "    markers_are_dominant: false\n";
        out << "    population_name_delimiter: \" \"\n";
        out << "    population_name_is_prefix: true\n";
        out << "    constant_sites_removed: true\n";
        out << "    equal_population_sizes: true\n";
        out << "    parameters:\n";
        out << "        population_size:\n";
        out << "            value: 0.005\n";
        out << "            estimate: false\n";
        out << "        freq_1:\n";
        out << "            value: 0.5\n";
        out << "            estimate: true\n";
        out << "            prior:\n";
        out << "                beta_distribution:\n";
        out << "                    alpha: " << alpha << "\n";
        out << "                    beta: " << beta << "\n";
        out << "        mutation_rate:\n";
        out << "            value: 1.0\n";
        out << "            estimate: false\n";
        out << "comparisons:\n";
        out << "- comparison:\n";
        out << "    path: hemi129.nex\n";
        out.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(3648);
        std::shared_ptr<OperatorInterface> op = std::make_shared<FreqMover>(1.0, 0.1);
        OperatorSchedule os = OperatorSchedule();
        // os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<BetaDistribution>(alpha, beta);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(os);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            o.perform_collection_move(rng, comparisons, 1);
            double x = comparisons.get_tree(0).get_freq_1();
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
        double alpha = 5.0;
        double beta = 1.0;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream out;
        out.open(test_path);
        out << "global_comparison_settings:\n";
        out << "    genotypes_are_diploid: true\n";
        out << "    markers_are_dominant: false\n";
        out << "    population_name_delimiter: \" \"\n";
        out << "    population_name_is_prefix: true\n";
        out << "    constant_sites_removed: true\n";
        out << "    equal_population_sizes: true\n";
        out << "    parameters:\n";
        out << "        population_size:\n";
        out << "            value: 0.005\n";
        out << "            estimate: false\n";
        out << "        freq_1:\n";
        out << "            value: 0.5\n";
        out << "            estimate: true\n";
        out << "            prior:\n";
        out << "                beta_distribution:\n";
        out << "                    alpha: " << alpha << "\n";
        out << "                    beta: " << beta << "\n";
        out << "        mutation_rate:\n";
        out << "            value: 1.0\n";
        out << "            estimate: false\n";
        out << "comparisons:\n";
        out << "- comparison:\n";
        out << "    path: hemi129.nex\n";
        out.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(3648);
        std::shared_ptr<OperatorInterface> op = std::make_shared<FreqMover>(1.0, 0.1);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<BetaDistribution>(alpha, beta);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(os);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            o.perform_collection_move(rng, comparisons, 1);
            double x = comparisons.get_tree(0).get_freq_1();
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

    SECTION("testing beta(5.0, 1.0) prior with optimizing") {
        double alpha = 5.0;
        double beta = 1.0;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream out;
        out.open(test_path);
        out << "global_comparison_settings:\n";
        out << "    genotypes_are_diploid: true\n";
        out << "    markers_are_dominant: false\n";
        out << "    population_name_delimiter: \" \"\n";
        out << "    population_name_is_prefix: true\n";
        out << "    constant_sites_removed: true\n";
        out << "    equal_population_sizes: true\n";
        out << "    parameters:\n";
        out << "        population_size:\n";
        out << "            value: 0.005\n";
        out << "            estimate: false\n";
        out << "        freq_1:\n";
        out << "            value: 0.5\n";
        out << "            estimate: true\n";
        out << "            prior:\n";
        out << "                beta_distribution:\n";
        out << "                    alpha: " << alpha << "\n";
        out << "                    beta: " << beta << "\n";
        out << "        mutation_rate:\n";
        out << "            value: 1.0\n";
        out << "            estimate: false\n";
        out << "comparisons:\n";
        out << "- comparison:\n";
        out << "    path: hemi129.nex\n";
        out.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(3648);
        std::shared_ptr<OperatorInterface> op = std::make_shared<FreqMover>(1.0, 0.1);
        OperatorSchedule os = OperatorSchedule();
        // os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<BetaDistribution>(alpha, beta);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(os);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            o.perform_collection_move(rng, comparisons, 1);
            double x = comparisons.get_tree(0).get_freq_1();
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

    SECTION("testing beta(1.0, 5.0) prior with optimizing") {
        double alpha = 1.0;
        double beta = 5.0;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-" + tag + "-state-run-1.log";
        std::ofstream out;
        out.open(test_path);
        out << "global_comparison_settings:\n";
        out << "    genotypes_are_diploid: true\n";
        out << "    markers_are_dominant: false\n";
        out << "    population_name_delimiter: \" \"\n";
        out << "    population_name_is_prefix: true\n";
        out << "    constant_sites_removed: true\n";
        out << "    equal_population_sizes: true\n";
        out << "    parameters:\n";
        out << "        population_size:\n";
        out << "            value: 0.005\n";
        out << "            estimate: false\n";
        out << "        freq_1:\n";
        out << "            value: 0.5\n";
        out << "            estimate: true\n";
        out << "            prior:\n";
        out << "                beta_distribution:\n";
        out << "                    alpha: " << alpha << "\n";
        out << "                    beta: " << beta << "\n";
        out << "        mutation_rate:\n";
        out << "            value: 1.0\n";
        out << "            estimate: false\n";
        out << "comparisons:\n";
        out << "- comparison:\n";
        out << "    path: hemi129.nex\n";
        out.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(3648);
        std::shared_ptr<OperatorInterface> op = std::make_shared<FreqMover>(1.0, 0.1);
        OperatorSchedule os = OperatorSchedule();
        // os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<BetaDistribution>(alpha, beta);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(os);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            o.perform_collection_move(rng, comparisons, 1);
            double x = comparisons.get_tree(0).get_freq_1();
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
}

TEST_CASE("Testing DirichletProcessGibbsSampler with 3 pairs and concentration 1.4142",
        "[DirichletProcessGibbsSampler]") {

    SECTION("Testing 3 pairs, conc 1.4142, with optimizing") {
        double concentration = 1.4142;
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-dpgibbssamper-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-dpgibbssampler-test1-" + tag + "-state-run-1.log";
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
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        population_size:\n";
        os << "            value: 0.002\n";
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

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<DirichletProcessGibbsSampler>(1.0, 2);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 3);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        std::vector<unsigned int> height_indices(ntrees, 0);
        unsigned int nevents;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                nevents = comparisons.get_number_of_events();
                height_indices = comparisons.get_standardized_height_indices();
                std::ostringstream stream;
                for (auto h_idx : height_indices) {
                    stream << h_idx;
                }
                std::string model_str = stream.str();
                if (model_counts.count(model_str) < 1) {
                    model_counts[model_str] = 1;
                }
                else {
                    ++model_counts[model_str];
                }
                if (nevent_counts.count(nevents) < 1) {
                    nevent_counts[nevents] = 1;
                }
                else {
                    ++nevent_counts[nevents];
                }
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
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
        }

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));
        unsigned int tally = 0;
        for (auto const & kv: model_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);
        tally = 0;
        for (auto const & kv: nevent_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: model_counts) {
            REQUIRE((kv.second / (double)nsamples) == Approx(std::exp(
                    get_dpp_log_prior_probability(kv.first, concentration))).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing DirichletProcessGibbsSampler with 6 pairs and concentration 1.7",
        "[DirichletProcessGibbsSampler]") {

    SECTION("Testing 6 pairs, conc 1.7, with optimizing") {
        double concentration = 1.7;
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-dpgibbssamper-test2-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-dpgibbssampler-test2-" + tag + "-state-run-1.log";
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
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        population_size:\n";
        os << "            value: 0.002\n";
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

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<DirichletProcessGibbsSampler>(1.0, 3);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 6);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        std::vector<unsigned int> height_indices(ntrees, 0);
        unsigned int nevents;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                nevents = comparisons.get_number_of_events();
                height_indices = comparisons.get_standardized_height_indices();
                std::ostringstream stream;
                for (auto h_idx : height_indices) {
                    stream << h_idx;
                }
                std::string model_str = stream.str();
                if (model_counts.count(model_str) < 1) {
                    model_counts[model_str] = 1;
                }
                else {
                    ++model_counts[model_str];
                }
                if (nevent_counts.count(nevents) < 1) {
                    nevent_counts[nevents] = 1;
                }
                else {
                    ++nevent_counts[nevents];
                }
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
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
        }

        REQUIRE(model_counts.at("000000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012345") == nevent_counts.at(6));
        unsigned int tally = 0;
        for (auto const & kv: model_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);
        tally = 0;
        for (auto const & kv: nevent_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: model_counts) {
            REQUIRE((kv.second / (double)nsamples) == Approx(std::exp(
                    get_dpp_log_prior_probability(kv.first, concentration))).epsilon(0.005));
        }
    }
}

TEST_CASE("Testing ReversibleJumpSampler with 3 pairs and split weight 1.0",
        "[ReversibleJumpSampler]") {

    SECTION("Testing 3 pairs and split weight 1.0,") {
        double split_weight = 1.0;
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-rjsampler-test1-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-rjsampler-test1-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    uniform:\n";
        os << "        parameters:\n";
        os << "            split_weight:\n";
        os << "                value: " << split_weight << "\n";
        os << "                estimate: false\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        population_size:\n";
        os << "            value: 0.002\n";
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

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ReversibleJumpSampler>(1.0);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 3);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        SampleSummarizer<double> split_weight_summary;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        std::vector<unsigned int> height_indices(ntrees, 0);
        unsigned int nevents;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                split_weight_summary.add_sample(comparisons.get_concentration());
                nevents = comparisons.get_number_of_events();
                height_indices = comparisons.get_standardized_height_indices();
                std::ostringstream stream;
                for (auto h_idx : height_indices) {
                    stream << h_idx;
                }
                std::string model_str = stream.str();
                if (model_counts.count(model_str) < 1) {
                    model_counts[model_str] = 1;
                }
                else {
                    ++model_counts[model_str];
                }
                if (nevent_counts.count(nevents) < 1) {
                    nevent_counts[nevents] = 1;
                }
                else {
                    ++nevent_counts[nevents];
                }
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));
        unsigned int tally = 0;
        for (auto const & kv: model_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);
        tally = 0;
        for (auto const & kv: nevent_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: model_counts) {
            REQUIRE((kv.second / (double)nsamples) == Approx(1.0/model_counts.size()).epsilon(0.001));
        }
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        }
        REQUIRE(split_weight_summary.sample_size() == nsamples);
        REQUIRE(split_weight_summary.mean() == Approx(1.0));
        REQUIRE(split_weight_summary.variance() == Approx(0.0));
    }
}

TEST_CASE("Testing ReversibleJumpSampler with 3 pairs and split weight 2.0",
        "[ReversibleJumpSampler]") {

    SECTION("Testing 3 pairs and split weight 2.0,") {
        double split_weight = 2.0;
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-rjsampler-test2-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-rjsampler-test2-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    uniform:\n";
        os << "        parameters:\n";
        os << "            split_weight:\n";
        os << "                value: " << split_weight << "\n";
        os << "                estimate: false\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        population_size:\n";
        os << "            value: 0.002\n";
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

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(123456);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ReversibleJumpSampler>(1.0);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 3);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        SampleSummarizer<double> split_weight_summary;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        std::vector<unsigned int> height_indices(ntrees, 0);
        unsigned int nevents;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                split_weight_summary.add_sample(comparisons.get_concentration());
                nevents = comparisons.get_number_of_events();
                height_indices = comparisons.get_standardized_height_indices();
                std::ostringstream stream;
                for (auto h_idx : height_indices) {
                    stream << h_idx;
                }
                std::string model_str = stream.str();
                if (model_counts.count(model_str) < 1) {
                    model_counts[model_str] = 1;
                }
                else {
                    ++model_counts[model_str];
                }
                if (nevent_counts.count(nevents) < 1) {
                    nevent_counts[nevents] = 1;
                }
                else {
                    ++nevent_counts[nevents];
                }
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));
        unsigned int tally = 0;
        for (auto const & kv: model_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);
        tally = 0;
        for (auto const & kv: nevent_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        REQUIRE(split_weight_summary.sample_size() == nsamples);
        REQUIRE(split_weight_summary.mean() == Approx(2.0));
        REQUIRE(split_weight_summary.variance() == Approx(0.0));
        for (std::string m: {"001", "010", "011"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(2.0/11.0).epsilon(0.001));
        }
        REQUIRE((model_counts.at("000") / (double)nsamples) == Approx(1.0/11.0).epsilon(0.001));
        REQUIRE((model_counts.at("012") / (double)nsamples) == Approx(4.0/11.0).epsilon(0.001));
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing ReversibleJumpSampler with 3 pairs and split weight 0.5",
        "[ReversibleJumpSampler]") {

    SECTION("Testing 3 pairs and split weight 0.5,") {
        double split_weight = 0.5;
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-rjsampler-test3-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-rjsampler-test3-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    uniform:\n";
        os << "        parameters:\n";
        os << "            split_weight:\n";
        os << "                value: " << split_weight << "\n";
        os << "                estimate: false\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        population_size:\n";
        os << "            value: 0.002\n";
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

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ReversibleJumpSampler>(1.0);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 3);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        SampleSummarizer<double> split_weight_summary;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        std::vector<unsigned int> height_indices(ntrees, 0);
        unsigned int nevents;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                split_weight_summary.add_sample(comparisons.get_concentration());
                nevents = comparisons.get_number_of_events();
                height_indices = comparisons.get_standardized_height_indices();
                std::ostringstream stream;
                for (auto h_idx : height_indices) {
                    stream << h_idx;
                }
                std::string model_str = stream.str();
                if (model_counts.count(model_str) < 1) {
                    model_counts[model_str] = 1;
                }
                else {
                    ++model_counts[model_str];
                }
                if (nevent_counts.count(nevents) < 1) {
                    nevent_counts[nevents] = 1;
                }
                else {
                    ++nevent_counts[nevents];
                }
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));
        unsigned int tally = 0;
        for (auto const & kv: model_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);
        tally = 0;
        for (auto const & kv: nevent_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        REQUIRE(split_weight_summary.sample_size() == nsamples);
        REQUIRE(split_weight_summary.mean() == Approx(0.5));
        REQUIRE(split_weight_summary.variance() == Approx(0.0));
        for (std::string m: {"001", "010", "011"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(2.0/11.0).epsilon(0.001));
        }
        REQUIRE((model_counts.at("012") / (double)nsamples) == Approx(1.0/11.0).epsilon(0.001));
        REQUIRE((model_counts.at("000") / (double)nsamples) == Approx(4.0/11.0).epsilon(0.001));
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing ReversibleJumpSampler with 4 pairs and split weight 3.0",
        "[ReversibleJumpSampler]") {

    SECTION("Testing 4 pairs and split weight 3.0,") {
        double split_weight = 3.0;
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-rjsampler-test4-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-rjsampler-test4-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    uniform:\n";
        os << "        parameters:\n";
        os << "            split_weight:\n";
        os << "                value: " << split_weight << "\n";
        os << "                estimate: false\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        population_size:\n";
        os << "            value: 0.002\n";
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
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ReversibleJumpSampler>(1.0);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        SampleSummarizer<double> split_weight_summary;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        std::vector<unsigned int> height_indices(ntrees, 0);
        unsigned int nevents;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                split_weight_summary.add_sample(comparisons.get_concentration());
                nevents = comparisons.get_number_of_events();
                height_indices = comparisons.get_standardized_height_indices();
                std::ostringstream stream;
                for (auto h_idx : height_indices) {
                    stream << h_idx;
                }
                std::string model_str = stream.str();
                if (model_counts.count(model_str) < 1) {
                    model_counts[model_str] = 1;
                }
                else {
                    ++model_counts[model_str];
                }
                if (nevent_counts.count(nevents) < 1) {
                    nevent_counts[nevents] = 1;
                }
                else {
                    ++nevent_counts[nevents];
                }
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(model_counts.at("0000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("0123") == nevent_counts.at(4));

        unsigned int tally = 0;
        for (auto const & kv: model_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);
        tally = 0;
        for (auto const & kv: nevent_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        REQUIRE(split_weight_summary.sample_size() == nsamples);
        REQUIRE(split_weight_summary.mean() == Approx(3.0));
        REQUIRE(split_weight_summary.variance() == Approx(0.0));

        // Below are the expected results. 'k' is the number of categories,
        // 'S(n,k)' is the number of possible partitions of n elements into k
        // categories (Stirling number of second kind), 'w' is the relative split
        // weight for any one of those possible partitions, 'total_w' is the
        // overall weight of the k class (w * S(n,k)), and 'prob' is the
        // probability of any one possible partition in the k class.
        //
        // k    S(n,k)  w       total_w     prob
        // -------------------------------------
        // 1    1       1       1           1/103
        // 2    7       3       21          3/103
        // 3    6       9       54          9/103
        // 4    1       27      27          27/103
        //                      103

        REQUIRE((model_counts.at("0123") / (double)nsamples) == Approx(27.0/103.0).epsilon(0.001));
        REQUIRE((model_counts.at("0122") / (double)nsamples) == Approx(9.0/103.0).epsilon(0.001));
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        }
    }
}

TEST_CASE("Testing ReversibleJumpSampler with 4 pairs and split weight 1/3.0",
        "[ReversibleJumpSampler]") {

    SECTION("Testing 4 pairs and split weight 1/3.0,") {
        double split_weight = 1.0/3.0;
        double height_shape = 5.0;
        double height_scale = 0.1;
        std::string tag = _TEST_OPERATOR_RNG.random_string(10);
        std::string test_path = "data/tmp-config-rjsampler-test5-" + tag + ".cfg";
        std::string log_path = "data/tmp-config-rjsampler-test5-" + tag + "-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    uniform:\n";
        os << "        parameters:\n";
        os << "            split_weight:\n";
        os << "                value: " << split_weight << "\n";
        os << "                estimate: false\n";
        os << "global_comparison_settings:\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: false\n";
        os << "    equal_population_sizes: true\n";
        os << "    parameters:\n";
        os << "        freq_1:\n";
        os << "            value: 0.5\n";
        os << "            estimate: false\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "        population_size:\n";
        os << "            value: 0.002\n";
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
        os.close();
        REQUIRE(path::exists(test_path));

        CollectionSettings settings = CollectionSettings(test_path);

        RandomNumberGenerator rng = RandomNumberGenerator(12345);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ReversibleJumpSampler>(1.0);
        OperatorSchedule op_schedule = OperatorSchedule();
        op_schedule.turn_on_auto_optimize();
        op_schedule.set_auto_optimize_delay(100);
        op_schedule.add_operator(op);

        ComparisonPopulationTreeCollection comparisons = ComparisonPopulationTreeCollection(settings, rng);
        comparisons.ignore_data();
        comparisons.set_operator_schedule(op_schedule);

        unsigned int ntrees = comparisons.get_number_of_trees();
        REQUIRE(ntrees == 4);
        std::vector< SampleSummarizer<double> > height_summaries(ntrees);
        SampleSummarizer<double> split_weight_summary;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 2;
        unsigned int nsamples = niterations / sample_freq;
        std::vector<unsigned int> height_indices(ntrees, 0);
        unsigned int nevents;
        for (unsigned int i = 0; i < niterations; ++i) {
            OperatorInterface& o = op_schedule.draw_operator(rng);
            o.operate(rng, comparisons, 1);
            if ((i + 1) % sample_freq == 0) {
                split_weight_summary.add_sample(comparisons.get_concentration());
                nevents = comparisons.get_number_of_events();
                height_indices = comparisons.get_standardized_height_indices();
                std::ostringstream stream;
                for (auto h_idx : height_indices) {
                    stream << h_idx;
                }
                std::string model_str = stream.str();
                if (model_counts.count(model_str) < 1) {
                    model_counts[model_str] = 1;
                }
                else {
                    ++model_counts[model_str];
                }
                if (nevent_counts.count(nevents) < 1) {
                    nevent_counts[nevents] = 1;
                }
                else {
                    ++nevent_counts[nevents];
                }
                for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
                    ComparisonPopulationTree & tree = comparisons.get_tree(tree_idx);
                    height_summaries.at(tree_idx).add_sample(tree.get_height());
                }
            }
        }
        std::cout << op->header_string();
        std::cout << op->to_string(op_schedule);

        REQUIRE(model_counts.at("0000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("0123") == nevent_counts.at(4));

        unsigned int tally = 0;
        for (auto const & kv: model_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);
        tally = 0;
        for (auto const & kv: nevent_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        REQUIRE(split_weight_summary.sample_size() == nsamples);
        REQUIRE(split_weight_summary.mean() == Approx(1.0/3.0));
        REQUIRE(split_weight_summary.variance() == Approx(0.0));

        // Below are the expected results. 'k' is the number of categories,
        // 'S(n,k)' is the number of possible partitions of n elements into k
        // categories (Stirling number of second kind), 'w' is the relative split
        // weight for any one of those possible partitions, 'total_w' is the
        // overall weight of the k class (w * S(n,k)), and 'prob' is the
        // probability of any one possible partition in the k class.
        //
        // k    S(n,k)  w       total_w     prob
        // -------------------------------------
        // 1    1       27      27          27/109
        // 2    7       9       63          9/109
        // 3    6       3       18          3/109
        // 4    1       1       1           1/109
        //                      109

        REQUIRE((model_counts.at("0000") / (double)nsamples) == Approx(27.0/109.0).epsilon(0.001));
        REQUIRE((model_counts.at("0001") / (double)nsamples) == Approx(9.0/109.0).epsilon(0.001));
        
        double size_sh;
        double size_sc;
        for (unsigned int tree_idx = 0; tree_idx < ntrees; ++tree_idx) {
            REQUIRE(height_summaries.at(tree_idx).sample_size() == nsamples);
            REQUIRE(height_summaries.at(tree_idx).mean() == Approx(height_shape * height_scale).epsilon(0.005));
            REQUIRE(height_summaries.at(tree_idx).variance() == Approx(height_shape * height_scale * height_scale).epsilon(0.01));
        }
    }
}
