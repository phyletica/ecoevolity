#include "catch.hpp"
#include "ecoevolity/general_tree_operator.hpp"
#include "ecoevolity/general_tree_operator_schedule.hpp"

#include <limits>
#include <memory>

#include "ecoevolity/tree.hpp"
#include "ecoevolity/node.hpp"
#include "ecoevolity/probability.hpp"
#include "ecoevolity/parameter.hpp"
#include "ecoevolity/stats_util.hpp"
#include "ecoevolity/rng.hpp"

#include "utils_for_testing.hpp"


TEST_CASE("Testing NeighborHeightNodePermute with 6 leaves",
        "[NeighborHeightNodePermute]") {

    SECTION("Testing 6 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(797853717);

        double root_ht = 0.5;
        std::shared_ptr<Node> root = std::make_shared<Node>(10, "root", root_ht);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(6, "internal0", 0.1);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(7, "internal1", 0.2);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(8, "internal2", 0.3);
        std::shared_ptr<Node> internal3 = std::make_shared<Node>(9, "internal3", 0.4);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>(5, "leaf5", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(internal0);
        internal2->add_child(leaf3);
        internal2->add_child(internal1);
        internal3->add_child(leaf4);
        internal3->add_child(internal2);
        root->add_child(leaf5);
        root->add_child(internal3);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        NeighborHeightNodePermute< BaseTree<Node> > op;
        NodeHeightScaler< BaseTree<Node> > op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);
        unsigned int num_height_moves = 2;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int niterations = 2000000;
        unsigned int sample_freq = 20;
        unsigned int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 10000;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            op2.operate(rng, &tree, 1, num_height_moves);
            if ((i + 1) % sample_freq == 0) {
                std::set< std::set<Split> > splits = tree.get_splits(false);
                if (split_counts.count(splits) > 0) {
                    ++split_counts[splits];
                }
                else {
                    split_counts[splits] = 1;
                }
                ++sample_count;
                if (sample_count % report_freq == 0) {
                    std::cout << "Sampled " << sample_count << " of " << nsamples << std::endl;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        double exp_freq = 1.0/945.0;
        double exp_count = nsamples/945.0;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.4;
        unsigned int total_trees_sampled = 0;
        std::map< std::set< std::set<Split> >, double> split_freqs;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            split_freqs[s_c.first] = s_c.second / (double)nsamples;
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }

        double quantile_chi_sq_944_95 = 1016.59;
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        std::cout << "Chi-square(944) 0.95 quantile: " << quantile_chi_sq_944_95 << "\n";

        std::cout << "BAD SPLITS (proportional error > " << prop_error_threshold << ")\n";
        for (auto s_e : bad_splits) {
            std::cout << "\nTree:\n";
            for (auto splitset : s_e.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            std::cout << "  prop error: " << s_e.second << "\n";
        }

        write_r_script(split_counts, 6, "../6-leaf-NeighborHeightNodePermute-test.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        REQUIRE(split_counts.size() == 945);

        double eps = 0.001;

        for (auto s_f : split_freqs) {
            REQUIRE(s_f.second == Approx(exp_freq).epsilon(eps));
        }

        REQUIRE(chi_sq_test_statistic < quantile_chi_sq_944_95);
    }
}


TEST_CASE("Testing NeighborHeightNodeSwap with 6 leaves",
        "[NeighborHeightNodeSwap]") {

    SECTION("Testing 6 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(797853717);

        double root_ht = 0.5;
        std::shared_ptr<Node> root = std::make_shared<Node>(10, "root", root_ht);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(6, "internal0", 0.1);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(7, "internal1", 0.2);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(8, "internal2", 0.3);
        std::shared_ptr<Node> internal3 = std::make_shared<Node>(9, "internal3", 0.4);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>(5, "leaf5", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(internal0);
        internal2->add_child(leaf3);
        internal2->add_child(internal1);
        internal3->add_child(leaf4);
        internal3->add_child(internal2);
        root->add_child(leaf5);
        root->add_child(internal3);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        NeighborHeightNodeSwap< BaseTree<Node> > op;
        NodeHeightScaler< BaseTree<Node> > op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);
        unsigned int num_height_moves = 2;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int niterations = 2000000;
        unsigned int sample_freq = 20;
        unsigned int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 10000;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            op2.operate(rng, &tree, 1, num_height_moves);
            if ((i + 1) % sample_freq == 0) {
                std::set< std::set<Split> > splits = tree.get_splits(false);
                if (split_counts.count(splits) > 0) {
                    ++split_counts[splits];
                }
                else {
                    split_counts[splits] = 1;
                }
                ++sample_count;
                if (sample_count % report_freq == 0) {
                    std::cout << "Sampled " << sample_count << " of " << nsamples << std::endl;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        double exp_freq = 1.0/945.0;
        double exp_count = nsamples/945.0;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.4;
        unsigned int total_trees_sampled = 0;
        std::map< std::set< std::set<Split> >, double> split_freqs;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            split_freqs[s_c.first] = s_c.second / (double)nsamples;
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }

        double quantile_chi_sq_944_95 = 1016.59;
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        std::cout << "Chi-square(944) 0.95 quantile: " << quantile_chi_sq_944_95 << "\n";

        std::cout << "BAD SPLITS (proportional error > " << prop_error_threshold << ")\n";
        for (auto s_e : bad_splits) {
            std::cout << "\nTree:\n";
            for (auto splitset : s_e.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            std::cout << "  prop error: " << s_e.second << "\n";
        }

        write_r_script(split_counts, 6, "../6-leaf-NeighborHeightNodeSwap-test.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        REQUIRE(split_counts.size() == 945);

        double eps = 0.001;

        for (auto s_f : split_freqs) {
            REQUIRE(s_f.second == Approx(exp_freq).epsilon(eps));
        }

        REQUIRE(chi_sq_test_statistic < quantile_chi_sq_944_95);
    }
}

TEST_CASE("Testing NeighborHeightNodeSwapAll with 6 leaves",
        "[NeighborHeightNodeSwapAll]") {

    SECTION("Testing 6 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(797853717);

        double root_ht = 0.5;
        std::shared_ptr<Node> root = std::make_shared<Node>(10, "root", root_ht);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(6, "internal0", 0.1);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(7, "internal1", 0.2);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(8, "internal2", 0.3);
        std::shared_ptr<Node> internal3 = std::make_shared<Node>(9, "internal3", 0.4);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>(5, "leaf5", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(internal0);
        internal2->add_child(leaf3);
        internal2->add_child(internal1);
        internal3->add_child(leaf4);
        internal3->add_child(internal2);
        root->add_child(leaf5);
        root->add_child(internal3);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        NeighborHeightNodeSwapAll< BaseTree<Node> > op;
        NodeHeightScaler< BaseTree<Node> > op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);
        unsigned int num_height_moves = 2;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int niterations = 2000000;
        unsigned int sample_freq = 20;
        unsigned int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 10000;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            op2.operate(rng, &tree, 1, num_height_moves);
            if ((i + 1) % sample_freq == 0) {
                std::set< std::set<Split> > splits = tree.get_splits(false);
                if (split_counts.count(splits) > 0) {
                    ++split_counts[splits];
                }
                else {
                    split_counts[splits] = 1;
                }
                ++sample_count;
                if (sample_count % report_freq == 0) {
                    std::cout << "Sampled " << sample_count << " of " << nsamples << std::endl;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        double exp_freq = 1.0/945.0;
        double exp_count = nsamples/945.0;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.4;
        unsigned int total_trees_sampled = 0;
        std::map< std::set< std::set<Split> >, double> split_freqs;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            split_freqs[s_c.first] = s_c.second / (double)nsamples;
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }

        double quantile_chi_sq_944_95 = 1016.59;
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        std::cout << "Chi-square(944) 0.95 quantile: " << quantile_chi_sq_944_95 << "\n";

        std::cout << "BAD SPLITS (proportional error > " << prop_error_threshold << ")\n";
        for (auto s_e : bad_splits) {
            std::cout << "\nTree:\n";
            for (auto splitset : s_e.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            std::cout << "  prop error: " << s_e.second << "\n";
        }

        write_r_script(split_counts, 6, "../6-leaf-NeighborHeightNodeSwapAll-test.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        REQUIRE(split_counts.size() == 945);

        double eps = 0.001;

        for (auto s_f : split_freqs) {
            REQUIRE(s_f.second == Approx(exp_freq).epsilon(eps));
        }

        REQUIRE(chi_sq_test_statistic < quantile_chi_sq_944_95);
    }
}

TEST_CASE("Testing NodeHeightSlideBumpSwapScaler with 6 leaves",
        "[NodeHeightSlideBumpSwapScaler]") {

    SECTION("Testing 6 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(797853717);

        double root_height_shape = 1.0;
        double root_height_scale = 0.5;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double root_ht = 0.5;
        std::shared_ptr<Node> root = std::make_shared<Node>(10, "root", root_ht);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(6, "internal0", 0.1);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(7, "internal1", 0.2);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(8, "internal2", 0.3);
        std::shared_ptr<Node> internal3 = std::make_shared<Node>(9, "internal3", 0.4);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>(5, "leaf5", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(internal0);
        internal2->add_child(leaf3);
        internal2->add_child(internal1);
        internal3->add_child(leaf4);
        internal3->add_child(internal2);
        root->add_child(leaf5);
        root->add_child(internal3);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.set_root_node_height_prior(root_height_prior);
        tree.estimate_root_height();

        NodeHeightSlideBumpSwapScaler< BaseTree<Node> > op;
        op.set_operate_on_root(true);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int niterations = 2000000;
        unsigned int sample_freq = 20;
        unsigned int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 10000;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                std::set< std::set<Split> > splits = tree.get_splits(false);
                if (split_counts.count(splits) > 0) {
                    ++split_counts[splits];
                }
                else {
                    split_counts[splits] = 1;
                }
                ++sample_count;
                if (sample_count % report_freq == 0) {
                    std::cout << "Sampled " << sample_count << " of " << nsamples << std::endl;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        double exp_freq = 1.0/945.0;
        double exp_count = nsamples/945.0;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.4;
        unsigned int total_trees_sampled = 0;
        std::map< std::set< std::set<Split> >, double> split_freqs;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            split_freqs[s_c.first] = s_c.second / (double)nsamples;
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }

        double quantile_chi_sq_944_95 = 1016.59;
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        std::cout << "Chi-square(944) 0.95 quantile: " << quantile_chi_sq_944_95 << "\n";

        std::cout << "BAD SPLITS (proportional error > " << prop_error_threshold << ")\n";
        for (auto s_e : bad_splits) {
            std::cout << "\nTree:\n";
            for (auto splitset : s_e.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            std::cout << "  prop error: " << s_e.second << "\n";
        }

        write_r_script(split_counts, 6, "../6-leaf-NodeHeightSlideBumpSwapScaler-test.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        REQUIRE(split_counts.size() == 945);

        double eps = 0.001;

        for (auto s_f : split_freqs) {
            REQUIRE(s_f.second == Approx(exp_freq).epsilon(eps));
        }

        REQUIRE(chi_sq_test_statistic < quantile_chi_sq_944_95);
    }
}

TEST_CASE("Testing NodeHeightSlideBumpSwapAllScaler with 6 leaves",
        "[NodeHeightSlideBumpSwapAllScaler]") {

    SECTION("Testing 6 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(797853717);

        double root_height_shape = 1.0;
        double root_height_scale = 0.5;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double root_ht = 0.5;
        std::shared_ptr<Node> root = std::make_shared<Node>(10, "root", root_ht);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(6, "internal0", 0.1);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(7, "internal1", 0.2);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(8, "internal2", 0.3);
        std::shared_ptr<Node> internal3 = std::make_shared<Node>(9, "internal3", 0.4);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>(5, "leaf5", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(internal0);
        internal2->add_child(leaf3);
        internal2->add_child(internal1);
        internal3->add_child(leaf4);
        internal3->add_child(internal2);
        root->add_child(leaf5);
        root->add_child(internal3);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.set_root_node_height_prior(root_height_prior);
        tree.estimate_root_height();

        NodeHeightSlideBumpSwapAllScaler< BaseTree<Node> > op;
        op.set_operate_on_root(true);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int niterations = 2000000;
        unsigned int sample_freq = 20;
        unsigned int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 10000;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                std::set< std::set<Split> > splits = tree.get_splits(false);
                if (split_counts.count(splits) > 0) {
                    ++split_counts[splits];
                }
                else {
                    split_counts[splits] = 1;
                }
                ++sample_count;
                if (sample_count % report_freq == 0) {
                    std::cout << "Sampled " << sample_count << " of " << nsamples << std::endl;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        double exp_freq = 1.0/945.0;
        double exp_count = nsamples/945.0;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.4;
        unsigned int total_trees_sampled = 0;
        std::map< std::set< std::set<Split> >, double> split_freqs;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            split_freqs[s_c.first] = s_c.second / (double)nsamples;
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }

        double quantile_chi_sq_944_95 = 1016.59;
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        std::cout << "Chi-square(944) 0.95 quantile: " << quantile_chi_sq_944_95 << "\n";

        std::cout << "BAD SPLITS (proportional error > " << prop_error_threshold << ")\n";
        for (auto s_e : bad_splits) {
            std::cout << "\nTree:\n";
            for (auto splitset : s_e.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            std::cout << "  prop error: " << s_e.second << "\n";
        }

        write_r_script(split_counts, 6, "../6-leaf-NodeHeightSlideBumpSwapAllScaler-test.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        REQUIRE(split_counts.size() == 945);

        double eps = 0.001;

        for (auto s_f : split_freqs) {
            REQUIRE(s_f.second == Approx(exp_freq).epsilon(eps));
        }

        REQUIRE(chi_sq_test_statistic < quantile_chi_sq_944_95);
    }
}


TEST_CASE("Testing NodeHeightSlideBumpPermuteScaler with 6 leaves",
        "[NodeHeightSlideBumpPermuteScaler]") {

    SECTION("Testing 6 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(797853717);

        double root_height_shape = 1.0;
        double root_height_scale = 0.5;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double root_ht = 0.5;
        std::shared_ptr<Node> root = std::make_shared<Node>(10, "root", root_ht);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(6, "internal0", 0.1);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(7, "internal1", 0.2);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(8, "internal2", 0.3);
        std::shared_ptr<Node> internal3 = std::make_shared<Node>(9, "internal3", 0.4);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>(5, "leaf5", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(internal0);
        internal2->add_child(leaf3);
        internal2->add_child(internal1);
        internal3->add_child(leaf4);
        internal3->add_child(internal2);
        root->add_child(leaf5);
        root->add_child(internal3);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.set_root_node_height_prior(root_height_prior);
        tree.estimate_root_height();

        NodeHeightSlideBumpPermuteScaler< BaseTree<Node> > op;
        op.set_operate_on_root(true);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int niterations = 2000000;
        unsigned int sample_freq = 20;
        unsigned int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 10000;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                std::set< std::set<Split> > splits = tree.get_splits(false);
                if (split_counts.count(splits) > 0) {
                    ++split_counts[splits];
                }
                else {
                    split_counts[splits] = 1;
                }
                ++sample_count;
                if (sample_count % report_freq == 0) {
                    std::cout << "Sampled " << sample_count << " of " << nsamples << std::endl;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        double exp_freq = 1.0/945.0;
        double exp_count = nsamples/945.0;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.4;
        unsigned int total_trees_sampled = 0;
        std::map< std::set< std::set<Split> >, double> split_freqs;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            split_freqs[s_c.first] = s_c.second / (double)nsamples;
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }

        double quantile_chi_sq_944_95 = 1016.59;
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        std::cout << "Chi-square(944) 0.95 quantile: " << quantile_chi_sq_944_95 << "\n";

        std::cout << "BAD SPLITS (proportional error > " << prop_error_threshold << ")\n";
        for (auto s_e : bad_splits) {
            std::cout << "\nTree:\n";
            for (auto splitset : s_e.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            std::cout << "  prop error: " << s_e.second << "\n";
        }

        write_r_script(split_counts, 6, "../6-leaf-NodeHeightSlideBumpPermuteScaler-test.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        REQUIRE(split_counts.size() == 945);

        double eps = 0.001;

        for (auto s_f : split_freqs) {
            REQUIRE(s_f.second == Approx(exp_freq).epsilon(eps));
        }

        REQUIRE(chi_sq_test_statistic < quantile_chi_sq_944_95);
    }
}


TEST_CASE("Testing NodeHeightSlideBumpSwapMover with 6 leaves",
        "[NodeHeightSlideBumpSwapMover]") {

    SECTION("Testing 6 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(78745872);

        double root_height_shape = 100.0;
        double root_height_scale = 0.01;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double root_ht = 1.0;
        std::shared_ptr<Node> root = std::make_shared<Node>(10, "root", root_ht);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(6, "internal0", 0.1);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(7, "internal1", 0.2);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(8, "internal2", 0.3);
        std::shared_ptr<Node> internal3 = std::make_shared<Node>(9, "internal3", 0.4);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>(5, "leaf5", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(internal0);
        internal2->add_child(leaf3);
        internal2->add_child(internal1);
        internal3->add_child(leaf4);
        internal3->add_child(internal2);
        root->add_child(leaf5);
        root->add_child(internal3);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.set_root_node_height_prior(root_height_prior);
        tree.estimate_root_height();

        // These *Mover operators have a difficult time sampling from the prior
        // distribution on their own. Having one proposal window width for all
        // the node heights in the tree does not work well given the prior on
        // node heights.
        NodeHeightSlideBumpSwapAllMover< BaseTree<Node> > op;
        op.set_operate_on_root(true);
        op.turn_off_auto_optimize();
        op.set_coercable_parameter_value(1.2);
        NodeHeightSlideBumpSwapAllMover< BaseTree<Node> > op2;
        op2.set_operate_on_root(true);
        op2.turn_off_auto_optimize();
        op2.set_coercable_parameter_value(1.0);
        NodeHeightSlideBumpSwapAllMover< BaseTree<Node> > op3;
        op3.set_operate_on_root(true);
        op3.turn_off_auto_optimize();
        op3.set_coercable_parameter_value(0.75);
        NodeHeightSlideBumpSwapAllMover< BaseTree<Node> > op4;
        op4.set_operate_on_root(true);
        op4.turn_off_auto_optimize();
        op4.set_coercable_parameter_value(0.05);
        NodeHeightSlideBumpSwapAllMover< BaseTree<Node> > op5;
        op5.set_operate_on_root(true);
        op5.turn_off_auto_optimize();
        op5.set_coercable_parameter_value(0.02);
        NodeHeightSlideBumpSwapAllMover< BaseTree<Node> > op6;
        op6.set_operate_on_root(true);
        op6.turn_on_auto_optimize();
        op6.set_auto_optimize_delay(100);
        NodeHeightScaler< BaseTree<Node> > op7;
        op7.turn_on_auto_optimize();
        op7.set_auto_optimize_delay(100);
        RootHeightScaler < BaseTree<Node> > op8;
        op8.turn_on_auto_optimize();
        op8.set_auto_optimize_delay(100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int niterations = 2000000;
        unsigned int sample_freq = 20;
        unsigned int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 10000;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);

            op2.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);

            op3.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);

            op4.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);

            op5.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);

            op6.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                std::set< std::set<Split> > splits = tree.get_splits(false);
                if (split_counts.count(splits) > 0) {
                    ++split_counts[splits];
                }
                else {
                    split_counts[splits] = 1;
                }
                ++sample_count;
                if (sample_count % report_freq == 0) {
                    std::cout << "Sampled " << sample_count << " of " << nsamples << std::endl;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();
        std::cout << op3.header_string();
        std::cout << op3.to_string();
        std::cout << op4.header_string();
        std::cout << op4.to_string();
        std::cout << op5.header_string();
        std::cout << op5.to_string();
        std::cout << op6.header_string();
        std::cout << op6.to_string();
        std::cout << op7.header_string();
        std::cout << op7.to_string();
        std::cout << op8.header_string();
        std::cout << op8.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        double exp_freq = 1.0/945.0;
        double exp_count = nsamples/945.0;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.4;
        unsigned int total_trees_sampled = 0;
        std::map< std::set< std::set<Split> >, double> split_freqs;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            split_freqs[s_c.first] = s_c.second / (double)nsamples;
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }

        double quantile_chi_sq_944_95 = 1016.59;
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        std::cout << "Chi-square(944) 0.95 quantile: " << quantile_chi_sq_944_95 << "\n";

        std::cout << "BAD SPLITS (proportional error > " << prop_error_threshold << ")\n";
        for (auto s_e : bad_splits) {
            std::cout << "\nTree:\n";
            for (auto splitset : s_e.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            std::cout << "  prop error: " << s_e.second << "\n";
        }

        write_r_script(split_counts, 6, "../6-leaf-NodeHeightSlideBumpSwapMover-test.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        REQUIRE(split_counts.size() == 945);

        double eps = 0.001;

        for (auto s_f : split_freqs) {
            REQUIRE(s_f.second == Approx(exp_freq).epsilon(eps));
        }

        REQUIRE(chi_sq_test_statistic < quantile_chi_sq_944_95);
    }
}


TEST_CASE("Testing NodeHeightSlideBumpSwapAllMover with 6 leaves",
        "[NodeHeightSlideBumpSwapAllMover]") {

    SECTION("Testing 6 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(24792874);

        double root_height_shape = 100.0;
        double root_height_scale = 0.01;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double root_ht = 1.0;
        std::shared_ptr<Node> root = std::make_shared<Node>(10, "root", root_ht);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(6, "internal0", 0.1);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(7, "internal1", 0.2);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(8, "internal2", 0.3);
        std::shared_ptr<Node> internal3 = std::make_shared<Node>(9, "internal3", 0.4);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>(5, "leaf5", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(internal0);
        internal2->add_child(leaf3);
        internal2->add_child(internal1);
        internal3->add_child(leaf4);
        internal3->add_child(internal2);
        root->add_child(leaf5);
        root->add_child(internal3);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.set_root_node_height_prior(root_height_prior);
        tree.estimate_root_height();

        // These *Mover operators have a difficult time sampling from the prior
        // distribution on their own. Having one proposal window width for all
        // the node heights in the tree does not work well given the prior on
        // node heights.
        NodeHeightSlideBumpSwapAllMover< BaseTree<Node> > op;
        op.set_operate_on_root(true);
        op.turn_off_auto_optimize();
        op.set_coercable_parameter_value(1.2);
        NodeHeightSlideBumpSwapAllMover< BaseTree<Node> > op2;
        op2.set_operate_on_root(true);
        op2.turn_off_auto_optimize();
        op2.set_coercable_parameter_value(1.0);
        NodeHeightSlideBumpSwapAllMover< BaseTree<Node> > op3;
        op3.set_operate_on_root(true);
        op3.turn_off_auto_optimize();
        op3.set_coercable_parameter_value(0.75);
        NodeHeightSlideBumpSwapAllMover< BaseTree<Node> > op4;
        op4.set_operate_on_root(true);
        op4.turn_off_auto_optimize();
        op4.set_coercable_parameter_value(0.05);
        NodeHeightSlideBumpSwapAllMover< BaseTree<Node> > op5;
        op5.set_operate_on_root(true);
        op5.turn_off_auto_optimize();
        op5.set_coercable_parameter_value(0.02);
        NodeHeightSlideBumpSwapAllMover< BaseTree<Node> > op6;
        op6.set_operate_on_root(true);
        op6.turn_on_auto_optimize();
        op6.set_auto_optimize_delay(100);
        NodeHeightScaler< BaseTree<Node> > op7;
        op7.turn_on_auto_optimize();
        op7.set_auto_optimize_delay(100);
        RootHeightScaler < BaseTree<Node> > op8;
        op8.turn_on_auto_optimize();
        op8.set_auto_optimize_delay(100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int niterations = 2000000;
        unsigned int sample_freq = 20;
        unsigned int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 10000;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);

            op2.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);

            op3.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);

            op4.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);

            op5.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);

            op6.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                std::set< std::set<Split> > splits = tree.get_splits(false);
                if (split_counts.count(splits) > 0) {
                    ++split_counts[splits];
                }
                else {
                    split_counts[splits] = 1;
                }
                ++sample_count;
                if (sample_count % report_freq == 0) {
                    std::cout << "Sampled " << sample_count << " of " << nsamples << std::endl;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();
        std::cout << op3.header_string();
        std::cout << op3.to_string();
        std::cout << op4.header_string();
        std::cout << op4.to_string();
        std::cout << op5.header_string();
        std::cout << op5.to_string();
        std::cout << op6.header_string();
        std::cout << op6.to_string();
        std::cout << op7.header_string();
        std::cout << op7.to_string();
        std::cout << op8.header_string();
        std::cout << op8.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        double exp_freq = 1.0/945.0;
        double exp_count = nsamples/945.0;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.4;
        unsigned int total_trees_sampled = 0;
        std::map< std::set< std::set<Split> >, double> split_freqs;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            split_freqs[s_c.first] = s_c.second / (double)nsamples;
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }

        double quantile_chi_sq_944_95 = 1016.59;
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        std::cout << "Chi-square(944) 0.95 quantile: " << quantile_chi_sq_944_95 << "\n";

        std::cout << "BAD SPLITS (proportional error > " << prop_error_threshold << ")\n";
        for (auto s_e : bad_splits) {
            std::cout << "\nTree:\n";
            for (auto splitset : s_e.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            std::cout << "  prop error: " << s_e.second << "\n";
        }

        write_r_script(split_counts, 6, "../6-leaf-NodeHeightSlideBumpSwapAllMover-test.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        REQUIRE(split_counts.size() == 945);

        double eps = 0.001;

        for (auto s_f : split_freqs) {
            REQUIRE(s_f.second == Approx(exp_freq).epsilon(eps));
        }

        REQUIRE(chi_sq_test_statistic < quantile_chi_sq_944_95);
    }
}


TEST_CASE("Testing NodeHeightSlideBumpPermuteMover with 6 leaves",
        "[NodeHeightSlideBumpPermuteMover]") {

    SECTION("Testing 6 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(2378983448);

        double root_height_shape = 100.0;
        double root_height_scale = 0.01;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double root_ht = 1.0;
        std::shared_ptr<Node> root = std::make_shared<Node>(10, "root", root_ht);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(6, "internal0", 0.1);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(7, "internal1", 0.2);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(8, "internal2", 0.3);
        std::shared_ptr<Node> internal3 = std::make_shared<Node>(9, "internal3", 0.4);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>(5, "leaf5", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(internal0);
        internal2->add_child(leaf3);
        internal2->add_child(internal1);
        internal3->add_child(leaf4);
        internal3->add_child(internal2);
        root->add_child(leaf5);
        root->add_child(internal3);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.set_root_node_height_prior(root_height_prior);
        tree.estimate_root_height();

        // These *Mover operators have a difficult time sampling from the prior
        // distribution on their own. Having one proposal window width for all
        // the node heights in the tree does not work well given the prior on
        // node heights.
        NodeHeightSlideBumpPermuteMover< BaseTree<Node> > op;
        op.set_operate_on_root(true);
        op.turn_off_auto_optimize();
        op.set_coercable_parameter_value(1.2);
        NodeHeightSlideBumpPermuteMover< BaseTree<Node> > op2;
        op2.set_operate_on_root(true);
        op2.turn_off_auto_optimize();
        op2.set_coercable_parameter_value(1.0);
        NodeHeightSlideBumpPermuteMover< BaseTree<Node> > op3;
        op3.set_operate_on_root(true);
        op3.turn_off_auto_optimize();
        op3.set_coercable_parameter_value(0.75);
        NodeHeightSlideBumpPermuteMover< BaseTree<Node> > op4;
        op4.set_operate_on_root(true);
        op4.turn_off_auto_optimize();
        op4.set_coercable_parameter_value(0.05);
        NodeHeightSlideBumpPermuteMover< BaseTree<Node> > op5;
        op5.set_operate_on_root(true);
        op5.turn_off_auto_optimize();
        op5.set_coercable_parameter_value(0.02);
        NodeHeightSlideBumpPermuteMover< BaseTree<Node> > op6;
        op6.set_operate_on_root(true);
        op6.turn_on_auto_optimize();
        op6.set_auto_optimize_delay(100);
        NodeHeightScaler< BaseTree<Node> > op7;
        op7.turn_on_auto_optimize();
        op7.set_auto_optimize_delay(100);
        RootHeightScaler < BaseTree<Node> > op8;
        op8.turn_on_auto_optimize();
        op8.set_auto_optimize_delay(100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int niterations = 2000000;
        unsigned int sample_freq = 20;
        unsigned int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 10000;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);

            op2.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);

            op3.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);

            op4.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);

            op5.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);

            op6.operate(rng, &tree, 1);
            op7.operate(rng, &tree, 1);
            op8.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                std::set< std::set<Split> > splits = tree.get_splits(false);
                if (split_counts.count(splits) > 0) {
                    ++split_counts[splits];
                }
                else {
                    split_counts[splits] = 1;
                }
                ++sample_count;
                if (sample_count % report_freq == 0) {
                    std::cout << "Sampled " << sample_count << " of " << nsamples << std::endl;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();
        std::cout << op3.header_string();
        std::cout << op3.to_string();
        std::cout << op4.header_string();
        std::cout << op4.to_string();
        std::cout << op5.header_string();
        std::cout << op5.to_string();
        std::cout << op6.header_string();
        std::cout << op6.to_string();
        std::cout << op7.header_string();
        std::cout << op7.to_string();
        std::cout << op8.header_string();
        std::cout << op8.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        double exp_freq = 1.0/945.0;
        double exp_count = nsamples/945.0;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.4;
        unsigned int total_trees_sampled = 0;
        std::map< std::set< std::set<Split> >, double> split_freqs;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            split_freqs[s_c.first] = s_c.second / (double)nsamples;
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }

        double quantile_chi_sq_944_95 = 1016.59;
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        std::cout << "Chi-square(944) 0.95 quantile: " << quantile_chi_sq_944_95 << "\n";

        std::cout << "BAD SPLITS (proportional error > " << prop_error_threshold << ")\n";
        for (auto s_e : bad_splits) {
            std::cout << "\nTree:\n";
            for (auto splitset : s_e.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            std::cout << "  prop error: " << s_e.second << "\n";
        }

        write_r_script(split_counts, 6, "../6-leaf-NodeHeightSlideBumpPermuteMover-test.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        REQUIRE(split_counts.size() == 945);

        double eps = 0.001;

        for (auto s_f : split_freqs) {
            REQUIRE(s_f.second == Approx(exp_freq).epsilon(eps));
        }

        REQUIRE(chi_sq_test_statistic < quantile_chi_sq_944_95);
    }
}
