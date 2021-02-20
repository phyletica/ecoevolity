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


TEST_CASE("Testing SplitLumpNodesRevJumpSampler with 6 leaves and fixed root",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing 6 leaves with fixed root") {
        RandomNumberGenerator rng = RandomNumberGenerator(2570892978);

        double root_ht = 0.5;
        std::shared_ptr<Node> root = std::make_shared<Node>(6, "root", root_ht);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>(5, "leaf5", 0.0);

        root->add_child(leaf0);
        root->add_child(leaf1);
        root->add_child(leaf2);
        root->add_child(leaf3);
        root->add_child(leaf4);
        root->add_child(leaf5);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        GeneralTreeOperatorSchedule< BaseTree<Node> > op_schedule;
        std::shared_ptr< GeneralTreeOperatorTemplate< BaseTree<Node> > > op;

        op = std::make_shared< SplitLumpNodesRevJumpSampler< BaseTree<Node> > >(1.0);
        op->turn_off_auto_optimize();
        op->set_coercable_parameter_value(1.5);
        op_schedule.add_operator(op);

        op = std::make_shared< SplitLumpNodesRevJumpSampler< BaseTree<Node> > >(1.0);
        op->turn_off_auto_optimize();
        op->set_coercable_parameter_value(5.0);
        op_schedule.add_operator(op);

        op = std::make_shared< SplitLumpNodesRevJumpSampler< BaseTree<Node> > >(1.0);
        op->turn_off_auto_optimize();
        op->set_coercable_parameter_value(50.0);
        op_schedule.add_operator(op);

        std::shared_ptr< NodeHeightScaler< BaseTree<Node> > > node_height_op = std::make_shared<NodeHeightScaler< BaseTree<Node> > >();
        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate< BaseTree<Node> > > > other_ops;
        other_ops.push_back(node_height_op);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int niterations = 1000000000;
        unsigned int sample_freq = 100;
        unsigned int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 200000;
        for (unsigned int i = 0; i < niterations; ++i) {
            op = op_schedule.draw_operator(rng);
            op->operate_plus(rng, &tree, other_ops, 1, 3);
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
        op_schedule.write_operator_rates(std::cout);

        // TODO: Figure this out
        unsigned int num_tree_models = split_counts.size();
        double exp_freq = 1.0/num_tree_models;
        double exp_count = nsamples/num_tree_models;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.2;
        unsigned int total_trees_sampled = 0;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            // std::cout << "Tree:\n";
            // for (auto splitset : s_c.first) {
            //     unsigned int s_count = 0;
            //     for (auto split : splitset) {
            //         if (s_count > 0) {
            //             // Indent shared splits
            //             std::cout << "  ";
            //         }
            //         std::cout << "  " << split.as_string() << "\n";
            //         ++s_count;
            //     }
            // }
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            // std::cout << "  nsamples: " << s_c.second << "\n";
            // std::cout << "  prop error: " << prop_error << "\n";
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }


        double quantile_chi_sq_5627_10 = 5763.4;
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        std::cout << "Chi-square(5627) 0.9 quantile: " << quantile_chi_sq_5627_10 << "\n";

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

        write_r_script(split_counts, 6, "../6-leaf-general-tree-test-beta-suite.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        // REQUIRE(split_counts.size() == ???);

        // REQUIRE(chi_sq_test_statistic < quantile_chi_sq_5627_10);
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler with 7 leaves and fixed root",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing 7 leaves with fixed root") {
        RandomNumberGenerator rng = RandomNumberGenerator(7429459230);

        double root_ht = 0.5;
        std::shared_ptr<Node> root = std::make_shared<Node>(7, "root", root_ht);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>(5, "leaf5", 0.0);
        std::shared_ptr<Node> leaf6 = std::make_shared<Node>(6, "leaf6", 0.0);

        root->add_child(leaf0);
        root->add_child(leaf1);
        root->add_child(leaf2);
        root->add_child(leaf3);
        root->add_child(leaf4);
        root->add_child(leaf5);
        root->add_child(leaf6);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        GeneralTreeOperatorSchedule< BaseTree<Node> > op_schedule;
        std::shared_ptr< GeneralTreeOperatorTemplate< BaseTree<Node> > > op;

        op = std::make_shared< SplitLumpNodesRevJumpSampler< BaseTree<Node> > >(1.0);
        op->turn_off_auto_optimize();
        op->set_coercable_parameter_value(1.5);
        op_schedule.add_operator(op);

        op = std::make_shared< SplitLumpNodesRevJumpSampler< BaseTree<Node> > >(1.0);
        op->turn_off_auto_optimize();
        op->set_coercable_parameter_value(5.0);
        op_schedule.add_operator(op);

        op = std::make_shared< SplitLumpNodesRevJumpSampler< BaseTree<Node> > >(1.0);
        op->turn_off_auto_optimize();
        op->set_coercable_parameter_value(50.0);
        op_schedule.add_operator(op);

        std::shared_ptr< NodeHeightScaler< BaseTree<Node> > > node_height_op = std::make_shared<NodeHeightScaler< BaseTree<Node> > >();
        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate< BaseTree<Node> > > > other_ops;
        other_ops.push_back(node_height_op);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned long long int niterations = 5000000000;
        unsigned int sample_freq = 100;
        unsigned long long int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 200000;
        for (unsigned long long int i = 0; i < niterations; ++i) {
            op = op_schedule.draw_operator(rng);
            op->operate_plus(rng, &tree, other_ops, 1, 3);
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
        op_schedule.write_operator_rates(std::cout);

        // niterations is likely beyond limit of unsigned int
        // REQUIRE(op.get_number_of_attempts() == niterations);

        // TODO: Figure this out
        unsigned int num_tree_models = split_counts.size();
        double exp_freq = 1.0/num_tree_models;
        double exp_count = nsamples/num_tree_models;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.2;
        unsigned long long int total_trees_sampled = 0;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            // std::cout << "Tree:\n";
            // for (auto splitset : s_c.first) {
            //     unsigned int s_count = 0;
            //     for (auto split : splitset) {
            //         if (s_count > 0) {
            //             // Indent shared splits
            //             std::cout << "  ";
            //         }
            //         std::cout << "  " << split.as_string() << "\n";
            //         ++s_count;
            //     }
            // }
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            // std::cout << "  nsamples: " << s_c.second << "\n";
            // std::cout << "  prop error: " << prop_error << "\n";
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }

        /* double quantile_chi_sq_5627_10 = 5763.4; */
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        /* std::cout << "Chi-square(5627) 0.9 quantile: " << quantile_chi_sq_5627_10 << "\n"; */

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

        write_r_script(split_counts, 7, "../7-leaf-general-tree-test-beta-suite.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        // REQUIRE(split_counts.size() == ???);

        /* REQUIRE(chi_sq_test_statistic < quantile_chi_sq_5627_10); */
    }
}
