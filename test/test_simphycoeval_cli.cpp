#include "catch.hpp"

#include "ecoevolity/simphycoeval.hpp"
#include "ecoevolity/rng.hpp"
#include "ecoevolity/path.hpp"
#include "ecoevolity/spreadsheet.hpp"

#include "utils_for_testing.hpp"


TEST_CASE("Testing simphycoeval cli with 5 leaves, full model, unconstrained sizes",
        "[simphycoeval]") {

    SECTION("Testing 5 leaves with full model, unconstrained sizes") {
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
        unsigned int nsamples = 20000;
        char arg4[] = "20000";
        char arg5[] = "--prefix";
        char arg6[] = "test-1-";
        char arg7[] = "-t";
        char arg8[] = "100";
        char arg9[] = "--parameters-only";
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
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        ret = simphycoeval_main<BasePopulationTree>(argc, argv);
        REQUIRE(ret == 0);

        std::string log_path = "data/test-1-simphycoeval-true-parameters.txt";
        std::string tree_path = "data/test-1-simphycoeval-true-trees.phy";

        REQUIRE(path::exists(log_path));
        REQUIRE(path::exists(tree_path));

        std::vector< BasePopulationTree > trees;
        get_trees<BasePopulationTree>(
                tree_path,
                "relaxedphyliptree",
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
        SampleSummarizer<double> l_leaf3_pop_size_summary = prior_sample.summarize<double>("pop_size_sp4");
        REQUIRE(l_leaf3_pop_size_summary.sample_size() == nsamples);
        SampleSummarizer<double> l_leaf4_pop_size_summary = prior_sample.summarize<double>("pop_size_sp5");
        REQUIRE(l_leaf4_pop_size_summary.sample_size() == nsamples);

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

        unsigned int count_nheights_1 = 0;
        unsigned int count_nheights_2 = 0;
        unsigned int count_nheights_3 = 0;
        unsigned int count_nheights_4 = 0;

        SampleSummarizer<double> internal_height_prior_sample;
        for (unsigned int i = 0; i < chain_length; ++i) {
            double a = height_alpha_prior->draw(rng);
            double v = BetaDistribution::get_draw(rng, a, 1.0);
            internal_height_prior_sample.add_sample(v);
        }

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_height_summary;
        SampleSummarizer<double> root_pop_size_summary;
        SampleSummarizer<double> leaf0_pop_size_summary;
        SampleSummarizer<double> leaf1_pop_size_summary;
        SampleSummarizer<double> leaf2_pop_size_summary;
        SampleSummarizer<double> leaf3_pop_size_summary;
        SampleSummarizer<double> leaf4_pop_size_summary;
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
            leaf3_pop_size_summary.add_sample(tree.get_node("sp4")->get_population_size());
            leaf4_pop_size_summary.add_sample(tree.get_node("sp5")->get_population_size());
            root_height_summary.add_sample(tree.get_root_height());
            for (unsigned int height_idx = 0;
                    height_idx < (tree.get_number_of_node_heights() - 1);
                    ++height_idx) {
                internal_height_summary.add_sample(tree.get_height(height_idx) / tree.get_height_of_youngest_parent(height_idx));
                if (height_idx == 0) {
                    internal_0_height_summary.add_sample(tree.get_height(height_idx) / tree.get_height_of_youngest_parent(height_idx));
                }
            }
            if (tree.get_number_of_node_heights() == 1) {
                ++count_nheights_1;
            }
            else if (tree.get_number_of_node_heights() == 2) {
                ++count_nheights_2;
            }
            else if (tree.get_number_of_node_heights() == 3) {
                ++count_nheights_3;
            }
            else if (tree.get_number_of_node_heights() == 4) {
                ++count_nheights_4;
            }
            std::set< std::set<Split> > splits = tree.get_splits(false);
            if (split_counts.count(splits) > 0) {
                ++split_counts[splits];
            }
            else {
                split_counts[splits] = 1;
            }
        }

        REQUIRE((count_nheights_1 + count_nheights_2 + count_nheights_3 + count_nheights_4) == nsamples);

        double freq_nheights_1 = count_nheights_1 / (double)nsamples;
        double freq_nheights_2 = count_nheights_2 / (double)nsamples;
        double freq_nheights_3 = count_nheights_3 / (double)nsamples;
        double freq_nheights_4 = count_nheights_4 / (double)nsamples;

        double exp_freq = 1.0/336.0;
        double exp_count = nsamples/336.0;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.2;
        unsigned int total_trees_sampled = 0;
        std::map< std::set< std::set<Split> >, double> split_freqs;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            split_freqs[s_c.first] = s_c.second / (double)nsamples;
            /* std::cout << "Tree:\n"; */
            /* for (auto splitset : s_c.first) { */
            /*     unsigned int s_count = 0; */
            /*     for (auto split : splitset) { */
            /*         if (s_count > 0) { */
            /*             // Indent shared splits */
            /*             std::cout << "  "; */
            /*         } */
            /*         std::cout << "  " << split.as_string() << "\n"; */
            /*         ++s_count; */
            /*     } */
            /* } */
            /* double prop_error = ((double)s_c.second - exp_count) / exp_count; */
            /* std::cout << "  nsamples: " << s_c.second << "\n"; */
            /* std::cout << "  prop error: " << prop_error << "\n"; */
            /* if (fabs(prop_error) > prop_error_threshold) { */
            /*     bad_splits[s_c.first] = prop_error; */
            /* } */
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }

        double quantile_chi_sq_335_10 = 368.6;
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        std::cout << "Chi-square(335) 0.9 quantile: " << quantile_chi_sq_335_10 << "\n";

        /* std::cout << "BAD SPLITS (proportional error > " << prop_error_threshold << ")\n"; */
        /* for (auto s_e : bad_splits) { */
        /*     std::cout << "\nTree:\n"; */
        /*     for (auto splitset : s_e.first) { */
        /*         unsigned int s_count = 0; */
        /*         for (auto split : splitset) { */
        /*             if (s_count > 0) { */
        /*                 // Indent shared splits */
        /*                 std::cout << "  "; */
        /*             } */
        /*             std::cout << "  " << split.as_string() << "\n"; */
        /*             ++s_count; */
        /*         } */
        /*     } */
        /*     std::cout << "  prop error: " << s_e.second << "\n"; */
        /* } */

        write_r_script(split_counts, "../simphyco-5-leaf-general-tree-test-full-model-free-pop-sizes.r");

        // We should sample every possible tree
        REQUIRE(split_counts.size() == 336);

        REQUIRE(root_pop_size_summary.mean() ==     Approx(l_root_pop_size_summary.mean()));
        REQUIRE(root_pop_size_summary.variance() == Approx(l_root_pop_size_summary.variance()));
        REQUIRE(leaf0_pop_size_summary.mean() ==     Approx(l_leaf0_pop_size_summary.mean()));
        REQUIRE(leaf0_pop_size_summary.variance() == Approx(l_leaf0_pop_size_summary.variance()));
        REQUIRE(leaf1_pop_size_summary.mean() ==     Approx(l_leaf1_pop_size_summary.mean()));
        REQUIRE(leaf1_pop_size_summary.variance() == Approx(l_leaf1_pop_size_summary.variance()));
        REQUIRE(leaf2_pop_size_summary.mean() ==     Approx(l_leaf2_pop_size_summary.mean()));
        REQUIRE(leaf2_pop_size_summary.variance() == Approx(l_leaf2_pop_size_summary.variance()));
        REQUIRE(leaf3_pop_size_summary.mean() ==     Approx(l_leaf3_pop_size_summary.mean()));
        REQUIRE(leaf3_pop_size_summary.variance() == Approx(l_leaf3_pop_size_summary.variance()));
        REQUIRE(leaf4_pop_size_summary.mean() ==     Approx(l_leaf4_pop_size_summary.mean()));
        REQUIRE(leaf4_pop_size_summary.variance() == Approx(l_leaf4_pop_size_summary.variance()));
        REQUIRE(root_height_summary.mean() ==     Approx(l_root_height_summary.mean()));
        REQUIRE(root_height_summary.variance() == Approx(l_root_height_summary.variance()));

        double eps = 0.005;
        REQUIRE(pop_size_summary.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
        REQUIRE(pop_size_summary.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        REQUIRE(root_pop_size_summary.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
        REQUIRE(root_pop_size_summary.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        REQUIRE(leaf0_pop_size_summary.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
        REQUIRE(leaf0_pop_size_summary.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        REQUIRE(leaf1_pop_size_summary.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
        REQUIRE(leaf1_pop_size_summary.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        REQUIRE(leaf2_pop_size_summary.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
        REQUIRE(leaf2_pop_size_summary.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        REQUIRE(leaf3_pop_size_summary.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
        REQUIRE(leaf3_pop_size_summary.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        REQUIRE(leaf4_pop_size_summary.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
        REQUIRE(leaf4_pop_size_summary.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));

        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(l_height_alpha_summary.mean() == Approx(height_alpha_prior->get_mean()).epsilon(eps * 2.0));
        REQUIRE(l_height_alpha_summary.variance() == Approx(height_alpha_prior->get_variance()).epsilon(eps * 4.0));

        REQUIRE(l_height_beta_summary.mean() == Approx(1.0));
        REQUIRE(l_height_beta_summary.variance() == Approx(0.0));

        REQUIRE(l_mu_rate_summary.mean() == Approx(mu_rate_prior->get_mean()).epsilon(eps));
        REQUIRE(l_mu_rate_summary.variance() == Approx(mu_rate_prior->get_variance()).epsilon(eps));
        REQUIRE(l_freq_summary.mean() == Approx(freq_prior->get_mean()).epsilon(eps));
        REQUIRE(l_freq_summary.variance() == Approx(freq_prior->get_variance()).epsilon(eps));

        // Node heights collected by height index will not match beta prior
        // when more than 3 tips
        // REQUIRE(internal_0_height_summary.mean() == Approx(internal_height_prior_sample.mean()).epsilon(eps));
        // REQUIRE(internal_0_height_summary.variance() == Approx(internal_height_prior_sample.variance()).epsilon(eps));
        // REQUIRE(internal_height_summary.mean() == Approx(internal_height_prior_sample.mean()).epsilon(eps));
        // REQUIRE(internal_height_summary.variance() == Approx(internal_height_prior_sample.variance()).epsilon(eps));

        for (auto s_f : split_freqs) {
            REQUIRE(s_f.second == Approx(exp_freq).epsilon(eps));
        }

        REQUIRE(chi_sq_test_statistic < quantile_chi_sq_335_10);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing simphycoeval cli with 3 leaves, constrained model",
        "[simphycoeval]") {

    SECTION("Testing 3 leaves with constrained model") {
        RandomNumberGenerator rng = RandomNumberGenerator(843579834);

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

        std::string tag = rng.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-phyco.cfg";
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

        char arg0[] = "simphycoeval";
        char arg1[] = "--seed";
        char arg2[] = "1234";
        char arg3[] = "-n";
        unsigned int nsamples = 10000;
        char arg4[] = "10000";
        char arg5[] = "--prefix";
        char arg6[] = "test-2-";
        char arg7[] = "--parameters-only";
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
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        ret = simphycoeval_main<BasePopulationTree>(argc, argv);
        REQUIRE(ret == 0);

        std::string log_path = "data/test-2-simphycoeval-true-parameters.txt";
        std::string tree_path = "data/test-2-simphycoeval-true-trees.phy";

        REQUIRE(path::exists(log_path));
        REQUIRE(path::exists(tree_path));

        std::vector< BasePopulationTree > trees;
        get_trees<BasePopulationTree>(
                tree_path,
                "relaxedphyliptree",
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

        double eps = 0.005;
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

TEST_CASE("Testing simphycoeval cli with 3 leaves, fully constrained model",
        "[simphycoeval]") {

    SECTION("Testing 3 leaves with fully constrained model") {
        RandomNumberGenerator rng = RandomNumberGenerator(8237493789);

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

        std::string tag = rng.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-phyco.cfg";
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

        char arg0[] = "simphycoeval";
        char arg1[] = "--seed";
        char arg2[] = "1234";
        char arg3[] = "-n";
        unsigned int nsamples = 1000;
        char arg4[] = "1000";
        char arg5[] = "--prefix";
        char arg6[] = "test-3-";
        char arg7[] = "--parameters-only";
        char arg8[] = "--fix-model";
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
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        ret = simphycoeval_main<BasePopulationTree>(argc, argv);
        REQUIRE(ret == 0);

        std::string log_path = "data/test-3-simphycoeval-true-parameters.txt";
        std::string tree_path = "data/test-3-simphycoeval-true-trees.phy";

        REQUIRE(path::exists(log_path));
        REQUIRE(path::exists(tree_path));

        std::vector< BasePopulationTree > trees;
        get_trees<BasePopulationTree>(
                tree_path,
                "relaxedphyliptree",
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

        double eps = 0.005;
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

        REQUIRE(internal_0_height_summary.mean() == Approx(0.5));
        REQUIRE(internal_0_height_summary.variance() == Approx(0.0));

        delete[] cfg_path;
    }
}

TEST_CASE("Testing 9 species with 2 genomes",
        "[simphycoeval]") {

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
