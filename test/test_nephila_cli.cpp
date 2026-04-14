#include "catch.hpp"

#include "ecoevolity/nephila.hpp"
#include "ecoevolity/rng.hpp"
#include "ecoevolity/path.hpp"
#include "ecoevolity/spreadsheet.hpp"

#include "utils_for_testing.hpp"


TEST_CASE("Testing nephila cli with rbcl data and comb starting tree",
        "[nephila]") {

    SECTION("Testing rbcl data with comb starting tree") {
        RandomNumberGenerator rng = RandomNumberGenerator(87264654);

        double root_height_shape = 1.0;
        double root_height_scale = 0.01;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        unsigned int chain_length = 1000;
        unsigned int sample_frequency = 10;
        unsigned int nsamples = (chain_length / sample_frequency) + 1;

        std::string tag = rng.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-nephila.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-nephila-state-run-1.log";
        std::string tree_path = "data/tmp-config-" + tag + "-nephila-trees-run-1.nex";
        std::ofstream os;
        os.open(test_path);
        os << "data:\n";
        os << "    path: rbcl.phy\n";
        os << "tree_model:\n";
        os << "    tree_space: generalized\n";
        // os << "    starting_tree: random\n";
        os << "    starting_tree: comb\n";
        // os << "    starting_tree: ./rbcl-ultrametric-root-poly.tre\n";
        // os << "    starting_tree: ./rbcl-ultrametric.tre\n";
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
        os << "                    value: 1.0\n";
        os << "                    estimate: false\n";
        os << "mutation_parameters:\n";
        os << "    state_frequencies:\n";
        os << "        value: [0.309769, 0.163380, 0.121023, 0.405828]\n";
        os << "        estimate: false\n";
        os << "    rate_matrix:\n";
        os << "        value: [0.08394222, 0.34116704, 0.03603322, 0.15737940, 0.30297095, 0.07850717]\n";
        os << "        estimate: false\n";
        os << "    among_site_rate_variation:\n";
        os << "        discrete_gamma:\n";
        os << "            number_of_categories: 4\n";
        os << "            parameters:\n";
        os << "                one_over_shape:\n";
        os << "                    value: 2.0\n";
        os << "                    estimate: false\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: " << chain_length << "\n";
        os << "    sample_frequency: " << sample_frequency << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "nephila";
        char arg1[] = "--seed";
        char arg2[] = "1234";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        ret = nephila_main< ecoevolity::SeqTree<Node> >(argc, argv);
        REQUIRE(ret == 0);

        REQUIRE(path::exists(tree_path));

        std::vector< ecoevolity::SeqTree<Node> > trees;
        get_trees< ecoevolity::SeqTree<Node> >(
                tree_path,
                "nexus",
                trees,
                0,
                1e-6);

        REQUIRE(trees.size() == nsamples);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        for (auto tree : trees) {
            std::set< std::set<Split> > splits = tree.get_splits(false);
            if (split_counts.count(splits) > 0) {
                ++split_counts[splits];
            }
            else {
                split_counts[splits] = 1;
            }
        }

        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";

        delete[] cfg_path;
    }
}
