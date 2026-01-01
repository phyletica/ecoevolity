#include "catch.hpp"

#include "ecoevolity/simphycoeval.hpp"
#include "ecoevolity/rng.hpp"
#include "ecoevolity/path.hpp"
#include "ecoevolity/spreadsheet.hpp"

#include "utils_for_testing.hpp"


TEST_CASE("Testing simphycoeval cli tree rejecting",
        "[simphycoeval]") {

    SECTION("Testing tree rejection") {
        RandomNumberGenerator rng = RandomNumberGenerator(7497983753);

        double mu_rate_shape = 10.0;
        double mu_rate_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> mu_rate_prior = std::make_shared<GammaDistribution>(
                mu_rate_shape,
                mu_rate_scale);

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double pop_size_shape = 10.0;
        double pop_size_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        double freq_a = 3.0;
        double freq_b = 2.0;
        std::shared_ptr<ContinuousProbabilityDistribution> freq_prior = std::make_shared<BetaDistribution>(
                freq_a,
                freq_b);

        double height_alpha_shape = 20.0;
        double height_alpha_scale = 0.4;
        std::shared_ptr<ContinuousProbabilityDistribution> height_alpha_prior = std::make_shared<GammaDistribution>(
                height_alpha_shape,
                height_alpha_scale);

        unsigned int chain_length = 500000;
        unsigned int sample_frequency = 20;

        std::string tag = rng.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-phyco.cfg";
        std::ofstream os;
        os.open(test_path);
        os << "data:\n";
        os << "    ploidy: 2\n";
        os << "    constant_sites_removed: false\n";
        os << "    yaml_allele_counts:\n";
        os << "        path: species-5-genomes-4-chars-1000.yml\n";
        os << "tree_model:\n";
        os << "    tree_space: generalized\n";
        os << "    starting_tree: comb\n";
        os << "    tree_prior:\n";
        os << "        uniform_root_and_betas:\n";
        os << "            parameters:\n";
        os << "                root_height:\n";
        os << "                    estimate: true\n";
        os << "                    prior:\n";
        os << "                        gamma_distribution:\n";
        os << "                            shape: " << root_height_shape << "\n";
        os << "                            scale: " << root_height_scale << "\n";
        os << "                alpha_of_node_height_beta_prior:\n";
        os << "                    estimate: true\n";
        os << "                    prior:\n";
        os << "                        gamma_distribution:\n";
        os << "                            shape: " << height_alpha_shape << "\n";
        os << "                            scale: " << height_alpha_scale << "\n";
        os << "branch_parameters:\n";
        os << "    population_size:\n";
        os << "        estimate: true\n";
        os << "        equal_population_sizes: false\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << pop_size_shape << "\n";
        os << "                scale: " << pop_size_scale << "\n";
        os << "mutation_parameters:\n";
        os << "    freq_1:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            beta_distribution:\n";
        os << "                alpha: " << freq_a << "\n";
        os << "                beta: " << freq_b << "\n";
        os << "    mutation_rate:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << mu_rate_shape << "\n";
        os << "                scale: " << mu_rate_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: " << chain_length << "\n";
        os << "    sample_frequency: " << sample_frequency << "\n";
        os.close();
        REQUIRE(path::exists(test_path));


        char arg0[] = "simphycoeval";
        char arg1[] = "--seed";
        char arg2[] = "1234";
        char arg3[] = "-n";
        unsigned int nsamples = 200;
        char arg4[] = "200";
        char arg5[] = "--prefix";
        char arg6[] = "test-12656154-";
        char arg7[] = "-t";
        char arg8[] = "100";
        char arg9[] = "--min-div-diff";
        double min_div_diff = 0.001;
        char arg10[] = "0.001";
        char arg11[] = "--parameters-only";
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
            &arg6[0],
            &arg7[0],
            &arg8[0],
            &arg9[0],
            &arg10[0],
            &arg11[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        ret = simphycoeval_main<BasePopulationTree>(argc, argv);
        REQUIRE(ret == 0);

        std::string log_path = "data/test-12656154-simphycoeval-true-parameters.txt";
        std::string tree_path = "data/test-12656154-simphycoeval-true-trees.phy";
        std::string rejected_trees_path = "data/test-12656154-simphycoeval-rejected-trees.phy";

        REQUIRE(path::exists(log_path));
        REQUIRE(path::exists(tree_path));
        REQUIRE(path::exists(rejected_trees_path));

        std::vector< BasePopulationTree > trees;
        get_trees<BasePopulationTree>(
                tree_path,
                "relaxedphyliptree",
                trees,
                0,
                1e-6);

        REQUIRE(trees.size() == nsamples);

        for (auto t : trees) {
            double min_diff = 9999999999.0;
            double last_time = t.get_height(0);
            for (unsigned int i = 1; i < t.get_number_of_node_heights(); ++i) {
                double diff = t.get_height(i) - last_time;
                if (diff < min_diff) {
                    min_diff = diff;
                }
                last_time = t.get_height(i);
            }
            REQUIRE(min_diff > min_div_diff);
        }

        std::vector< BasePopulationTree > rejected_trees;
        get_trees<BasePopulationTree>(
                rejected_trees_path,
                "relaxedphyliptree",
                rejected_trees,
                0,
                1e-6);

        REQUIRE(rejected_trees.size() > 1);

        for (auto t : rejected_trees) {
            double min_diff = 9999999999.0;
            double last_time = t.get_height(0);
            for (unsigned int i = 1; i < t.get_number_of_node_heights(); ++i) {
                double diff = t.get_height(i) - last_time;
                if (diff < min_diff) {
                    min_diff = diff;
                }
                last_time = t.get_height(i);
            }
            REQUIRE(min_diff < min_div_diff);
        }

        delete[] cfg_path;
    }
}


TEST_CASE("Testing 9 species with 2 genomes",
        "[simphycoeval]") {

    // This dataset was causing simphycoeval to crash due to invalid internal
    // node indices borking the simulation of gene trees.  Internal node
    // indices are irrelevant elsewhere, but they are used by
    // BasePopulationTree::simulate_gene_tree.
    SECTION("Testing generalized model with 9 species with 2 genomes") {
        RandomNumberGenerator rng = RandomNumberGenerator(16878646);

        std::string tag = rng.random_string(10);
        std::string out_prefix_str = "test-out-" + tag + "-species-9-genomes-2-";
        std::string test_path = "data/species-9-genomes-2-generalized-tree-random.yml";
        REQUIRE(path::exists(test_path));


        char arg0[] = "simphycoeval";
        char arg1[] = "--seed";
        char arg2[] = "1234";
        char arg3[] = "-n";
        unsigned int nsamples = 10;
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "-t";
        char arg7[] = "10";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * out_prefix = new char[out_prefix_str.size() + 1];
        std::copy(out_prefix_str.begin(), out_prefix_str.end(), out_prefix);
        out_prefix[out_prefix_str.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            out_prefix,
            &arg6[0],
            &arg7[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        ret = simphycoeval_main<BasePopulationTree>(argc, argv);
        REQUIRE(ret == 0);

        std::string out_model_cfg_path = "data/" + out_prefix_str + "simphycoeval-model-used-for-sims.yml";

        REQUIRE(path::exists(out_model_cfg_path));

        delete[] cfg_path;
        delete[] out_prefix;
    }
}
