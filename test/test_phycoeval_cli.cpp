#include "catch.hpp"

#include "ecoevolity/phycoeval.hpp"
#include "ecoevolity/rng.hpp"
#include "ecoevolity/path.hpp"
#include "ecoevolity/spreadsheet.hpp"

#include "utils_for_testing.hpp"


TEST_CASE("Testing phycoeval cli with 3 leaves, constrained model",
        "[phycoeval]") {

    SECTION("Testing 3 leaves with constrained model") {
        RandomNumberGenerator rng = RandomNumberGenerator(839798);

        double root_height = 0.2;

        double mu_rate = 5.0;

        double pop_size = 0.01;

        double freq_1 = 0.3;

        double height_alpha = 2.0;

        std::string starting_tree_str = "[&R]((sp1:0.1,sp2:0.1):0.1,sp3:0.2):0;";
        BaseTree<Node> starting_tree(starting_tree_str);
        std::set< std::set<Split> > starting_tree_splits = starting_tree.get_splits(true);

        unsigned int chain_length = 50000;
        unsigned int sample_frequency = 10;
        unsigned int nsamples = (chain_length / sample_frequency) + 1;

        std::string tag = rng.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-phyco.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-phyco-state-run-1.log";
        std::string tree_path = "data/tmp-config-" + tag + "-phyco-trees-run-1.nex";
        std::ofstream os;
        os.open(test_path);
        os << "data:\n";
        os << "    ploidy: 2\n";
        os << "    constant_sites_removed: false\n";
        os << "    yaml_allele_counts:\n";
        os << "        path: species-3-genomes-4-chars-1000.yml\n";
        os << "tree_model:\n";
        os << "    tree_space: fixed\n";
        os << "    starting_tree: \"" << starting_tree_str << "\"\n";
        os << "    tree_prior:\n";
        os << "        uniform_root_and_betas:\n";
        os << "            parameters:\n";
        os << "                root_height:\n";
        os << "                    estimate: false\n";
        os << "                alpha_of_node_height_beta_prior:\n";
        os << "                    value: " << height_alpha << "\n";
        os << "                    estimate: false\n";
        os << "branch_parameters:\n";
        os << "    population_size:\n";
        os << "        equal_population_sizes: true\n";
        os << "        value: " << pop_size << "\n";
        os << "        estimate: false\n";
        os << "mutation_parameters:\n";
        os << "    freq_1:\n";
        os << "        value: " << freq_1 << "\n";
        os << "        estimate: false\n";
        os << "    mutation_rate:\n";
        os << "        value: " << mu_rate << "\n";
        os << "        estimate: false\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: " << chain_length << "\n";
        os << "    sample_frequency: " << sample_frequency << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "phycoeval";
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

        ret = phycoeval_main<BasePopulationTree>(argc, argv);
        REQUIRE(ret == 0);

        REQUIRE(path::exists(log_path));
        REQUIRE(path::exists(tree_path));

        std::vector< BasePopulationTree > trees;
        get_trees<BasePopulationTree>(
                tree_path,
                "nexus",
                trees,
                0,
                1e-6);

        REQUIRE(trees.size() == nsamples);

        spreadsheet::Spreadsheet prior_sample;
        prior_sample.update(log_path);

        SampleSummarizer<double> l_root_pop_size_summary = prior_sample.summarize<double>("pop_size_root");
        REQUIRE(l_root_pop_size_summary.sample_size() == nsamples);
        SampleSummarizer<double> l_leaf0_pop_size_summary = prior_sample.summarize<double>("pop_size_sp1");
        REQUIRE(l_leaf0_pop_size_summary.sample_size() == nsamples);
        SampleSummarizer<double> l_leaf1_pop_size_summary = prior_sample.summarize<double>("pop_size_sp2");
        REQUIRE(l_leaf1_pop_size_summary.sample_size() == nsamples);
        SampleSummarizer<double> l_leaf2_pop_size_summary = prior_sample.summarize<double>("pop_size_sp3");
        REQUIRE(l_leaf2_pop_size_summary.sample_size() == nsamples);

        SampleSummarizer<double> l_root_height_summary = prior_sample.summarize<double>("root_height");
        REQUIRE(l_root_height_summary.sample_size() == nsamples);
        SampleSummarizer<double> l_height_alpha_summary = prior_sample.summarize<double>("alpha_of_height_beta_prior");
        REQUIRE(l_height_alpha_summary.sample_size() == nsamples);
        SampleSummarizer<double> l_height_beta_summary = prior_sample.summarize<double>("beta_of_height_beta_prior");
        REQUIRE(l_height_beta_summary.sample_size() == nsamples);

        SampleSummarizer<double> l_mu_rate_summary = prior_sample.summarize<double>("mutation_rate");
        REQUIRE(l_mu_rate_summary.sample_size() == nsamples);
        SampleSummarizer<double> l_freq_summary = prior_sample.summarize<double>("freq_1");
        REQUIRE(l_freq_summary.sample_size() == nsamples);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;
        split_counts[starting_tree_splits] = 0;

        unsigned int count_nheights_1 = 0;
        unsigned int count_nheights_2 = 0;

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> root_pop_size_summary;
        SampleSummarizer<double> leaf0_pop_size_summary;
        SampleSummarizer<double> leaf1_pop_size_summary;
        SampleSummarizer<double> leaf2_pop_size_summary;
        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> pop_size_summary;

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes;

        for (auto tree : trees) {
            pop_sizes = tree.get_pointers_to_population_sizes();
            for (auto pop_size : pop_sizes) {
                pop_size_summary.add_sample(pop_size->get_value());
            }
            root_pop_size_summary.add_sample(tree.get_root_population_size());
            leaf0_pop_size_summary.add_sample(tree.get_node("sp1")->get_population_size());
            leaf1_pop_size_summary.add_sample(tree.get_node("sp2")->get_population_size());
            leaf2_pop_size_summary.add_sample(tree.get_node("sp3")->get_population_size());
            root_height_summary.add_sample(tree.get_root_height());
            internal_0_height_summary.add_sample(tree.get_height(0) / tree.get_height_of_youngest_parent(0));
            if (tree.get_number_of_node_heights() == 1) {
                ++count_nheights_1;
            }
            else if (tree.get_number_of_node_heights() == 2) {
                ++count_nheights_2;
            }
            std::set< std::set<Split> > splits = tree.get_splits(false);
            ++split_counts[splits];
        }

        REQUIRE(count_nheights_2 == nsamples);
        REQUIRE(count_nheights_1 == 0);
        REQUIRE(split_counts[starting_tree_splits] == nsamples);
        REQUIRE(split_counts.size() == 1);

        REQUIRE(root_pop_size_summary.mean() ==     Approx(l_root_pop_size_summary.mean()));
        REQUIRE(root_pop_size_summary.variance() == Approx(l_root_pop_size_summary.variance()));
        REQUIRE(leaf0_pop_size_summary.mean() ==     Approx(l_leaf0_pop_size_summary.mean()));
        REQUIRE(leaf0_pop_size_summary.variance() == Approx(l_leaf0_pop_size_summary.variance()));
        REQUIRE(leaf1_pop_size_summary.mean() ==     Approx(l_leaf1_pop_size_summary.mean()));
        REQUIRE(leaf1_pop_size_summary.variance() == Approx(l_leaf1_pop_size_summary.variance()));
        REQUIRE(leaf2_pop_size_summary.mean() ==     Approx(l_leaf2_pop_size_summary.mean()));
        REQUIRE(leaf2_pop_size_summary.variance() == Approx(l_leaf2_pop_size_summary.variance()));
        REQUIRE(root_height_summary.mean() ==     Approx(l_root_height_summary.mean()));
        REQUIRE(root_height_summary.variance() == Approx(l_root_height_summary.variance()));

        double eps = 0.01;
        REQUIRE(pop_size_summary.mean() == Approx(pop_size));
        REQUIRE(pop_size_summary.variance() == Approx(0.0));
        REQUIRE(root_pop_size_summary.mean() == Approx(pop_size));
        REQUIRE(root_pop_size_summary.variance() == Approx(0.0));
        REQUIRE(leaf0_pop_size_summary.mean() == Approx(pop_size));
        REQUIRE(leaf0_pop_size_summary.variance() == Approx(0.0));
        REQUIRE(leaf1_pop_size_summary.mean() == Approx(pop_size));
        REQUIRE(leaf1_pop_size_summary.variance() == Approx(0.0));
        REQUIRE(leaf2_pop_size_summary.mean() == Approx(pop_size));
        REQUIRE(leaf2_pop_size_summary.variance() == Approx(0.0));

        REQUIRE(root_height_summary.mean() == Approx(root_height));
        REQUIRE(root_height_summary.variance() == Approx(0.0));
        REQUIRE(l_height_alpha_summary.mean() == Approx(height_alpha));
        REQUIRE(l_height_alpha_summary.variance() == Approx(0.0));

        REQUIRE(l_height_beta_summary.mean() == Approx(1.0));
        REQUIRE(l_height_beta_summary.variance() == Approx(0.0));

        REQUIRE(l_mu_rate_summary.mean() == Approx(mu_rate));
        REQUIRE(l_mu_rate_summary.variance() == Approx(0.0));
        REQUIRE(l_freq_summary.mean() == Approx(freq_1));
        REQUIRE(l_freq_summary.variance() == Approx(0.0));

        BetaDistribution beta_ht_prior(height_alpha, 1.0);
        REQUIRE(internal_0_height_summary.mean() == Approx(beta_ht_prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(beta_ht_prior.get_variance()).epsilon(eps));

        delete[] cfg_path;
    }
}
