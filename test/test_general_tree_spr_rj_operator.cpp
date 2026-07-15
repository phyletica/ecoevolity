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


TEST_CASE("Testing SubtreePruneRegraftRevJumpSampler::get_branch_and_node_target_weights",
        "[SubtreePruneRegraftRevJumpSampler]") {
    SECTION("Testing get_branch_and_node_target_weights with 7 leaf tree") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> n1 = std::make_shared<Node>("node1", 0.3);
        std::shared_ptr<Node> n2 = std::make_shared<Node>("node2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.2);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.1);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(0, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(1, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(2, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(3, "leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>(4, "leaf5", 0.0);
        std::shared_ptr<Node> leaf6 = std::make_shared<Node>(5, "leaf6", 0.0);
        std::shared_ptr<Node> leaf7 = std::make_shared<Node>(6, "leaf7", 0.0);
        n2->add_child(leaf1);
        n2->add_child(leaf2);
        internal1->add_child(leaf3);
        internal1->add_child(leaf4);
        internal2->add_child(leaf5);
        internal2->add_child(leaf6);
        n1->add_child(internal1);
        n1->add_child(internal2);
        n1->add_child(leaf7);
        root->add_child(n1);
        root->add_child(n2);

        BaseTree<Node> tree(root);
        unsigned int num_nodes = tree.get_node_count();
        REQUIRE(num_nodes == 12);

        double neg_inf = -std::numeric_limits<double>::infinity();
        std::vector<double> expected_weights(num_nodes * 2, neg_inf);

        double ln_distance_mult;

        SubtreePruneRegraftRevJumpSampler< BaseTree<Node> > op;
        ln_distance_mult = 0.0;
        op.set_coercable_parameter_value(ln_distance_mult);

        std::vector<double> ln_weights;
        ln_weights = op.get_branch_and_node_target_weights(&tree, leaf7);
        REQUIRE(ln_weights.size() == expected_weights.size());
        expected_weights.at(n1->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(internal1->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(internal2->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(leaf3->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(leaf4->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(leaf5->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(leaf6->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(root->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(n2->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(leaf1->get_index()) = 0.0 + (ln_distance_mult * 2);
        expected_weights.at(leaf2->get_index()) = 0.0 + (ln_distance_mult * 2);

        expected_weights.at(num_nodes + internal1->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(num_nodes + internal2->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(num_nodes + root->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(num_nodes + n2->get_index()) = 0.0 + (ln_distance_mult * 2);

        unsigned int n_i;
        for (unsigned int i = 0; i < ln_weights.size(); ++i) {
            n_i = i;
            std::cout << "\nweight index: " << i << "\n";
            if (i < num_nodes) {
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "branch weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            else {
                n_i = i - num_nodes;
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "node weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            REQUIRE(ln_weights.at(i) == expected_weights.at(i));
        }

        ln_distance_mult = 1.0;
        op.set_coercable_parameter_value(ln_distance_mult);
        ln_weights = op.get_branch_and_node_target_weights(&tree, leaf7);
        REQUIRE(ln_weights.size() == expected_weights.size());

        std::fill(expected_weights.begin(), expected_weights.end(), neg_inf); 
        expected_weights.at(n1->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(internal1->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(internal2->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(leaf3->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(leaf4->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(leaf5->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(leaf6->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(root->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(n2->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(leaf1->get_index()) = 0.0 + (ln_distance_mult * 2);
        expected_weights.at(leaf2->get_index()) = 0.0 + (ln_distance_mult * 2);

        expected_weights.at(num_nodes + internal1->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(num_nodes + internal2->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(num_nodes + root->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(num_nodes + n2->get_index()) = 0.0 + (ln_distance_mult * 2);

        for (unsigned int i = 0; i < ln_weights.size(); ++i) {
            n_i = i;
            std::cout << "\nweight index: " << i << "\n";
            if (i < num_nodes) {
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "branch weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            else {
                n_i = i - num_nodes;
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "node weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            REQUIRE(ln_weights.at(i) == expected_weights.at(i));
        }

        ln_distance_mult = -1.0;
        op.set_coercable_parameter_value(ln_distance_mult);
        ln_weights = op.get_branch_and_node_target_weights(&tree, leaf7);
        REQUIRE(ln_weights.size() == expected_weights.size());

        std::fill(expected_weights.begin(), expected_weights.end(), neg_inf); 
        expected_weights.at(n1->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(internal1->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(internal2->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(leaf3->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(leaf4->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(leaf5->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(leaf6->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(root->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(n2->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(leaf1->get_index()) = 0.0 + (ln_distance_mult * 2);
        expected_weights.at(leaf2->get_index()) = 0.0 + (ln_distance_mult * 2);

        expected_weights.at(num_nodes + internal1->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(num_nodes + internal2->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(num_nodes + root->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(num_nodes + n2->get_index()) = 0.0 + (ln_distance_mult * 2);

        for (unsigned int i = 0; i < ln_weights.size(); ++i) {
            n_i = i;
            std::cout << "\nweight index: " << i << "\n";
            if (i < num_nodes) {
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "branch weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            else {
                n_i = i - num_nodes;
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "node weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            REQUIRE(ln_weights.at(i) == expected_weights.at(i));
        }


        ln_distance_mult = 0.0;
        op.set_coercable_parameter_value(ln_distance_mult);
        ln_weights = op.get_branch_and_node_target_weights(&tree, n2);
        REQUIRE(ln_weights.size() == expected_weights.size());

        std::fill(expected_weights.begin(), expected_weights.end(), neg_inf); 
        expected_weights.at(n1->get_index()) = 0.0 + (ln_distance_mult * 0);

        for (unsigned int i = 0; i < ln_weights.size(); ++i) {
            n_i = i;
            std::cout << "\nweight index: " << i << "\n";
            if (i < num_nodes) {
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "branch weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            else {
                n_i = i - num_nodes;
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "node weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            REQUIRE(ln_weights.at(i) == expected_weights.at(i));
        }

        ln_distance_mult = 1.5;
        op.set_coercable_parameter_value(ln_distance_mult);
        ln_weights = op.get_branch_and_node_target_weights(&tree, n2);
        REQUIRE(ln_weights.size() == expected_weights.size());

        std::fill(expected_weights.begin(), expected_weights.end(), neg_inf); 
        expected_weights.at(n1->get_index()) = 0.0 + (ln_distance_mult * 0);

        for (unsigned int i = 0; i < ln_weights.size(); ++i) {
            n_i = i;
            std::cout << "\nweight index: " << i << "\n";
            if (i < num_nodes) {
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "branch weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            else {
                n_i = i - num_nodes;
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "node weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            REQUIRE(ln_weights.at(i) == expected_weights.at(i));
        }

        ln_distance_mult = -1.5;
        op.set_coercable_parameter_value(ln_distance_mult);
        ln_weights = op.get_branch_and_node_target_weights(&tree, n2);
        REQUIRE(ln_weights.size() == expected_weights.size());

        std::fill(expected_weights.begin(), expected_weights.end(), neg_inf); 
        expected_weights.at(n1->get_index()) = 0.0 + (ln_distance_mult * 0);

        for (unsigned int i = 0; i < ln_weights.size(); ++i) {
            n_i = i;
            std::cout << "\nweight index: " << i << "\n";
            if (i < num_nodes) {
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "branch weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            else {
                n_i = i - num_nodes;
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "node weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            REQUIRE(ln_weights.at(i) == expected_weights.at(i));
        }


        ln_distance_mult = 0.0;
        op.set_coercable_parameter_value(ln_distance_mult);
        ln_weights = op.get_branch_and_node_target_weights(&tree, n1);
        REQUIRE(ln_weights.size() == expected_weights.size());

        std::fill(expected_weights.begin(), expected_weights.end(), neg_inf); 
        expected_weights.at(n2->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(leaf1->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(leaf2->get_index()) = 0.0 + (ln_distance_mult * 0);

        expected_weights.at(num_nodes + n2->get_index()) = 0.0 + (ln_distance_mult * 0);

        for (unsigned int i = 0; i < ln_weights.size(); ++i) {
            n_i = i;
            std::cout << "\nweight index: " << i << "\n";
            if (i < num_nodes) {
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "branch weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            else {
                n_i = i - num_nodes;
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "node weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            REQUIRE(ln_weights.at(i) == expected_weights.at(i));
        }

        ln_distance_mult = 0.5;
        op.set_coercable_parameter_value(ln_distance_mult);
        ln_weights = op.get_branch_and_node_target_weights(&tree, n1);
        REQUIRE(ln_weights.size() == expected_weights.size());

        std::fill(expected_weights.begin(), expected_weights.end(), neg_inf); 
        expected_weights.at(n2->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(leaf1->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(leaf2->get_index()) = 0.0 + (ln_distance_mult * 0);

        expected_weights.at(num_nodes + n2->get_index()) = 0.0 + (ln_distance_mult * 0);

        for (unsigned int i = 0; i < ln_weights.size(); ++i) {
            n_i = i;
            std::cout << "\nweight index: " << i << "\n";
            if (i < num_nodes) {
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "branch weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            else {
                n_i = i - num_nodes;
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "node weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            REQUIRE(ln_weights.at(i) == expected_weights.at(i));
        }

        ln_distance_mult = -0.5;
        op.set_coercable_parameter_value(ln_distance_mult);
        ln_weights = op.get_branch_and_node_target_weights(&tree, n1);
        REQUIRE(ln_weights.size() == expected_weights.size());

        std::fill(expected_weights.begin(), expected_weights.end(), neg_inf); 
        expected_weights.at(n2->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(leaf1->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(leaf2->get_index()) = 0.0 + (ln_distance_mult * 0);

        expected_weights.at(num_nodes + n2->get_index()) = 0.0 + (ln_distance_mult * 0);

        for (unsigned int i = 0; i < ln_weights.size(); ++i) {
            n_i = i;
            std::cout << "\nweight index: " << i << "\n";
            if (i < num_nodes) {
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "branch weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            else {
                n_i = i - num_nodes;
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "node weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            REQUIRE(ln_weights.at(i) == expected_weights.at(i));
        }

        ln_distance_mult = 0.0;
        op.set_coercable_parameter_value(ln_distance_mult);
        ln_weights = op.get_branch_and_node_target_weights(&tree, internal1);
        REQUIRE(ln_weights.size() == expected_weights.size());

        std::fill(expected_weights.begin(), expected_weights.end(), neg_inf); 
        expected_weights.at(n2->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(root->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(leaf1->get_index()) = 0.0 + (ln_distance_mult * 2);
        expected_weights.at(leaf2->get_index()) = 0.0 + (ln_distance_mult * 2);
        expected_weights.at(n1->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(internal2->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(leaf7->get_index()) = 0.0 + (ln_distance_mult * 0);

        expected_weights.at(num_nodes + root->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(num_nodes + n2->get_index()) = 0.0 + (ln_distance_mult * 2);

        for (unsigned int i = 0; i < ln_weights.size(); ++i) {
            n_i = i;
            std::cout << "\nweight index: " << i << "\n";
            if (i < num_nodes) {
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "branch weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            else {
                n_i = i - num_nodes;
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "node weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            REQUIRE(ln_weights.at(i) == expected_weights.at(i));
        }

        ln_distance_mult = -0.1;
        op.set_coercable_parameter_value(ln_distance_mult);
        ln_weights = op.get_branch_and_node_target_weights(&tree, internal1);
        REQUIRE(ln_weights.size() == expected_weights.size());

        std::fill(expected_weights.begin(), expected_weights.end(), neg_inf); 
        expected_weights.at(n2->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(root->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(leaf1->get_index()) = 0.0 + (ln_distance_mult * 2);
        expected_weights.at(leaf2->get_index()) = 0.0 + (ln_distance_mult * 2);
        expected_weights.at(n1->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(internal2->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(leaf7->get_index()) = 0.0 + (ln_distance_mult * 0);

        expected_weights.at(num_nodes + root->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(num_nodes + n2->get_index()) = 0.0 + (ln_distance_mult * 2);

        for (unsigned int i = 0; i < ln_weights.size(); ++i) {
            n_i = i;
            std::cout << "\nweight index: " << i << "\n";
            if (i < num_nodes) {
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "branch weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            else {
                n_i = i - num_nodes;
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "node weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            REQUIRE(ln_weights.at(i) == expected_weights.at(i));
        }

        ln_distance_mult = 0.1;
        op.set_coercable_parameter_value(ln_distance_mult);
        ln_weights = op.get_branch_and_node_target_weights(&tree, internal1);
        REQUIRE(ln_weights.size() == expected_weights.size());

        std::fill(expected_weights.begin(), expected_weights.end(), neg_inf); 
        expected_weights.at(n2->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(root->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(leaf1->get_index()) = 0.0 + (ln_distance_mult * 2);
        expected_weights.at(leaf2->get_index()) = 0.0 + (ln_distance_mult * 2);
        expected_weights.at(n1->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(internal2->get_index()) = 0.0 + (ln_distance_mult * 0);
        expected_weights.at(leaf7->get_index()) = 0.0 + (ln_distance_mult * 0);

        expected_weights.at(num_nodes + root->get_index()) = 0.0 + (ln_distance_mult * 1);
        expected_weights.at(num_nodes + n2->get_index()) = 0.0 + (ln_distance_mult * 2);

        for (unsigned int i = 0; i < ln_weights.size(); ++i) {
            n_i = i;
            std::cout << "\nweight index: " << i << "\n";
            if (i < num_nodes) {
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "branch weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            else {
                n_i = i - num_nodes;
                std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                std::cout << "node weight: " << ln_weights.at(i) << "\n";
                std::cout << "expected weight: " << expected_weights.at(i) << "\n";
            }
            REQUIRE(ln_weights.at(i) == expected_weights.at(i));
        }

        std::vector<double> ln_dist_multipliers = {0.0, -1.2, 1.2};

        for (const auto ln_distance_mult : ln_dist_multipliers) {
            op.set_coercable_parameter_value(ln_distance_mult);
            ln_weights = op.get_branch_and_node_target_weights(&tree, internal2);
            REQUIRE(ln_weights.size() == expected_weights.size());

            std::fill(expected_weights.begin(), expected_weights.end(), neg_inf); 
            expected_weights.at(n2->get_index()) = 0.0 + (ln_distance_mult * 1);
            expected_weights.at(root->get_index()) = 0.0 + (ln_distance_mult * 1);
            expected_weights.at(leaf1->get_index()) = 0.0 + (ln_distance_mult * 2);
            expected_weights.at(leaf2->get_index()) = 0.0 + (ln_distance_mult * 2);
            expected_weights.at(n1->get_index()) = 0.0 + (ln_distance_mult * 0);
            expected_weights.at(internal1->get_index()) = 0.0 + (ln_distance_mult * 0);
            expected_weights.at(leaf7->get_index()) = 0.0 + (ln_distance_mult * 0);
            expected_weights.at(leaf3->get_index()) = 0.0 + (ln_distance_mult * 1);
            expected_weights.at(leaf4->get_index()) = 0.0 + (ln_distance_mult * 1);

            expected_weights.at(num_nodes + root->get_index()) = 0.0 + (ln_distance_mult * 1);
            expected_weights.at(num_nodes + n2->get_index()) = 0.0 + (ln_distance_mult * 2);
            expected_weights.at(num_nodes + internal1->get_index()) = 0.0 + (ln_distance_mult * 1);

            for (unsigned int i = 0; i < ln_weights.size(); ++i) {
                n_i = i;
                std::cout << "\nweight index: " << i << "\n";
                if (i < num_nodes) {
                    std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                    std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                    std::cout << "branch weight: " << ln_weights.at(i) << "\n";
                    std::cout << "expected weight: " << expected_weights.at(i) << "\n";
                }
                else {
                    n_i = i - num_nodes;
                    std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                    std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                    std::cout << "node weight: " << ln_weights.at(i) << "\n";
                    std::cout << "expected weight: " << expected_weights.at(i) << "\n";
                }
                REQUIRE(ln_weights.at(i) == expected_weights.at(i));
            }
        }

        for (const auto ln_distance_mult : ln_dist_multipliers) {
            op.set_coercable_parameter_value(ln_distance_mult);
            ln_weights = op.get_branch_and_node_target_weights(&tree, leaf5);
            REQUIRE(ln_weights.size() == expected_weights.size());

            std::fill(expected_weights.begin(), expected_weights.end(), neg_inf); 
            expected_weights.at(n2->get_index()) = 0.0 + (ln_distance_mult * 1);
            expected_weights.at(root->get_index()) = 0.0 + (ln_distance_mult * 1);
            expected_weights.at(leaf1->get_index()) = 0.0 + (ln_distance_mult * 2);
            expected_weights.at(leaf2->get_index()) = 0.0 + (ln_distance_mult * 2);
            expected_weights.at(n1->get_index()) = 0.0 + (ln_distance_mult * 0);
            expected_weights.at(internal1->get_index()) = 0.0 + (ln_distance_mult * 0);
            expected_weights.at(leaf6->get_index()) = 0.0 + (ln_distance_mult * 0);
            expected_weights.at(leaf7->get_index()) = 0.0 + (ln_distance_mult * 0);
            expected_weights.at(leaf3->get_index()) = 0.0 + (ln_distance_mult * 1);
            expected_weights.at(leaf4->get_index()) = 0.0 + (ln_distance_mult * 1);

            expected_weights.at(num_nodes + root->get_index()) = 0.0 + (ln_distance_mult * 1);
            expected_weights.at(num_nodes + n2->get_index()) = 0.0 + (ln_distance_mult * 2);
            expected_weights.at(num_nodes + n1->get_index()) = 0.0 + (ln_distance_mult * 0);
            expected_weights.at(num_nodes + internal1->get_index()) = 0.0 + (ln_distance_mult * 1);

            for (unsigned int i = 0; i < ln_weights.size(); ++i) {
                n_i = i;
                std::cout << "\nweight index: " << i << "\n";
                if (i < num_nodes) {
                    std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                    std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                    std::cout << "branch weight: " << ln_weights.at(i) << "\n";
                    std::cout << "expected weight: " << expected_weights.at(i) << "\n";
                }
                else {
                    n_i = i - num_nodes;
                    std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                    std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                    std::cout << "node weight: " << ln_weights.at(i) << "\n";
                    std::cout << "expected weight: " << expected_weights.at(i) << "\n";
                }
                REQUIRE(ln_weights.at(i) == expected_weights.at(i));
            }
        }

        for (const auto ln_distance_mult : ln_dist_multipliers) {
            op.set_coercable_parameter_value(ln_distance_mult);
            ln_weights = op.get_branch_and_node_target_weights(&tree, leaf1);
            REQUIRE(ln_weights.size() == expected_weights.size());

            std::fill(expected_weights.begin(), expected_weights.end(), neg_inf); 
            expected_weights.at(root->get_index()) = 0.0 + (ln_distance_mult * 0);
            expected_weights.at(leaf2->get_index()) = 0.0 + (ln_distance_mult * 0);
            expected_weights.at(n1->get_index()) = 0.0 + (ln_distance_mult * 0);
            expected_weights.at(internal1->get_index()) = 0.0 + (ln_distance_mult * 1);
            expected_weights.at(internal2->get_index()) = 0.0 + (ln_distance_mult * 1);
            expected_weights.at(leaf7->get_index()) = 0.0 + (ln_distance_mult * 1);
            expected_weights.at(leaf6->get_index()) = 0.0 + (ln_distance_mult * 2);
            expected_weights.at(leaf5->get_index()) = 0.0 + (ln_distance_mult * 2);
            expected_weights.at(leaf3->get_index()) = 0.0 + (ln_distance_mult * 2);
            expected_weights.at(leaf4->get_index()) = 0.0 + (ln_distance_mult * 2);

            expected_weights.at(num_nodes + root->get_index()) = 0.0 + (ln_distance_mult * 0);
            expected_weights.at(num_nodes + n1->get_index()) = 0.0 + (ln_distance_mult * 1);
            expected_weights.at(num_nodes + internal1->get_index()) = 0.0 + (ln_distance_mult * 2);
            expected_weights.at(num_nodes + internal2->get_index()) = 0.0 + (ln_distance_mult * 2);

            for (unsigned int i = 0; i < ln_weights.size(); ++i) {
                n_i = i;
                std::cout << "\nweight index: " << i << "\n";
                if (i < num_nodes) {
                    std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                    std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                    std::cout << "branch weight: " << ln_weights.at(i) << "\n";
                    std::cout << "expected weight: " << expected_weights.at(i) << "\n";
                }
                else {
                    n_i = i - num_nodes;
                    std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                    std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                    std::cout << "node weight: " << ln_weights.at(i) << "\n";
                    std::cout << "expected weight: " << expected_weights.at(i) << "\n";
                }
                REQUIRE(ln_weights.at(i) == expected_weights.at(i));
            }
        }
    }

    SECTION("Testing get_branch_and_node_target_weights with 8 leaf tree") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> root_child = std::make_shared<Node>("root_child", 0.5);
        std::shared_ptr<Node> n1 = std::make_shared<Node>("node1", 0.3);
        std::shared_ptr<Node> n2 = std::make_shared<Node>("node2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.2);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.1);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(0, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(1, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(2, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(3, "leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>(4, "leaf5", 0.0);
        std::shared_ptr<Node> leaf6 = std::make_shared<Node>(5, "leaf6", 0.0);
        std::shared_ptr<Node> leaf7 = std::make_shared<Node>(6, "leaf7", 0.0);
        std::shared_ptr<Node> leaf8 = std::make_shared<Node>(7, "leaf8", 0.0);
        n2->add_child(leaf1);
        n2->add_child(leaf2);
        internal1->add_child(leaf3);
        internal1->add_child(leaf4);
        internal2->add_child(leaf5);
        internal2->add_child(leaf6);
        n1->add_child(internal1);
        n1->add_child(internal2);
        n1->add_child(leaf7);
        root_child->add_child(n1);
        root_child->add_child(n2);
        root->add_child(root_child);
        root->add_child(leaf8);

        BaseTree<Node> tree(root);
        unsigned int num_nodes = tree.get_node_count();
        REQUIRE(num_nodes == 14);

        double neg_inf = -std::numeric_limits<double>::infinity();
        std::vector<double> expected_weights(num_nodes * 2, neg_inf);
        std::vector<double> ln_weights;

        SubtreePruneRegraftRevJumpSampler< BaseTree<Node> > op;

        std::vector<double> ln_dist_multipliers = {0.0, -2.1, 2.1};
        unsigned int n_i;

        for (const auto ln_distance_mult : ln_dist_multipliers) {
            op.set_coercable_parameter_value(ln_distance_mult);
            ln_weights = op.get_branch_and_node_target_weights(&tree, root_child);
            REQUIRE(ln_weights.size() == expected_weights.size());

            std::fill(expected_weights.begin(), expected_weights.end(), neg_inf); 

            expected_weights.at(leaf8->get_index()) = 0.0 + (ln_distance_mult * 0);

            for (unsigned int i = 0; i < ln_weights.size(); ++i) {
                n_i = i;
                std::cout << "\nweight index: " << i << "\n";
                if (i < num_nodes) {
                    std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                    std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                    std::cout << "branch weight: " << ln_weights.at(i) << "\n";
                    std::cout << "expected weight: " << expected_weights.at(i) << "\n";
                }
                else {
                    n_i = i - num_nodes;
                    std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                    std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                    std::cout << "node weight: " << ln_weights.at(i) << "\n";
                    std::cout << "expected weight: " << expected_weights.at(i) << "\n";
                }
                REQUIRE(ln_weights.at(i) == expected_weights.at(i));
            }
        }
    }

    SECTION("Testing get_branch_and_node_target_weights with 9 leaf tree") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> n1 = std::make_shared<Node>("node1", 0.3);
        std::shared_ptr<Node> n2 = std::make_shared<Node>("node2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.2);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.1);
        std::shared_ptr<Node> internal3 = std::make_shared<Node>("internal3", 0.05);
        std::shared_ptr<Node> internal4 = std::make_shared<Node>("internal4", 0.01);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(0, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(1, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(2, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(3, "leaf4", 0.0);
        std::shared_ptr<Node> leaf4b = std::make_shared<Node>(7, "leaf4b", 0.0);
        std::shared_ptr<Node> leaf4c = std::make_shared<Node>(8, "leaf4c", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>(4, "leaf5", 0.0);
        std::shared_ptr<Node> leaf6 = std::make_shared<Node>(5, "leaf6", 0.0);
        std::shared_ptr<Node> leaf7 = std::make_shared<Node>(6, "leaf7", 0.0);
        n2->add_child(leaf1);
        n2->add_child(leaf2);
        internal1->add_child(leaf3);
        internal2->add_child(leaf5);
        internal2->add_child(leaf6);
        internal4->add_child(leaf4b);
        internal4->add_child(leaf4c);
        internal3->add_child(internal4);
        internal3->add_child(leaf4);
        internal1->add_child(internal3);
        n1->add_child(internal1);
        n1->add_child(internal2);
        n1->add_child(leaf7);
        root->add_child(n1);
        root->add_child(n2);

        BaseTree<Node> tree(root);
        unsigned int num_nodes = tree.get_node_count();
        REQUIRE(num_nodes == 16);

        double neg_inf = -std::numeric_limits<double>::infinity();
        std::vector<double> expected_weights(num_nodes * 2, neg_inf);
        std::vector<double> ln_weights;

        SubtreePruneRegraftRevJumpSampler< BaseTree<Node> > op;

        std::vector<double> ln_dist_multipliers = {0.0, -2.1, 2.1};
        unsigned int n_i;

        for (const auto ln_distance_mult : ln_dist_multipliers) {
            op.set_coercable_parameter_value(ln_distance_mult);
            ln_weights = op.get_branch_and_node_target_weights(&tree, leaf4b);
            REQUIRE(ln_weights.size() == expected_weights.size());

            std::fill(expected_weights.begin(), expected_weights.end(), neg_inf); 

            expected_weights.at(leaf4c->get_index()) = 0.0 + (ln_distance_mult * 0);
            expected_weights.at(leaf4->get_index()) = 0.0 + (ln_distance_mult * 0);
            expected_weights.at(internal3->get_index()) = 0.0 + (ln_distance_mult * 0);
            expected_weights.at(leaf3->get_index()) = 0.0 + (ln_distance_mult * 1);
            expected_weights.at(internal1->get_index()) = 0.0 + (ln_distance_mult * 1);
            expected_weights.at(leaf7->get_index()) = 0.0 + (ln_distance_mult * 2);
            expected_weights.at(internal2->get_index()) = 0.0 + (ln_distance_mult * 2);
            expected_weights.at(n1->get_index()) = 0.0 + (ln_distance_mult * 2);
            expected_weights.at(leaf5->get_index()) = 0.0 + (ln_distance_mult * 3);
            expected_weights.at(leaf6->get_index()) = 0.0 + (ln_distance_mult * 3);
            expected_weights.at(root->get_index()) = 0.0 + (ln_distance_mult * 3);
            expected_weights.at(n2->get_index()) = 0.0 + (ln_distance_mult * 3);
            expected_weights.at(leaf1->get_index()) = 0.0 + (ln_distance_mult * 4);
            expected_weights.at(leaf2->get_index()) = 0.0 + (ln_distance_mult * 4);

            expected_weights.at(num_nodes + internal3->get_index()) = 0.0 + (ln_distance_mult * 0);
            expected_weights.at(num_nodes + internal1->get_index()) = 0.0 + (ln_distance_mult * 1);
            expected_weights.at(num_nodes + n1->get_index()) = 0.0 + (ln_distance_mult * 2);
            expected_weights.at(num_nodes + root->get_index()) = 0.0 + (ln_distance_mult * 3);
            expected_weights.at(num_nodes + internal2->get_index()) = 0.0 + (ln_distance_mult * 3);
            expected_weights.at(num_nodes + n2->get_index()) = 0.0 + (ln_distance_mult * 4);

            for (unsigned int i = 0; i < ln_weights.size(); ++i) {
                n_i = i;
                std::cout << "\nweight index: " << i << "\n";
                if (i < num_nodes) {
                    std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                    std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                    std::cout << "branch weight: " << ln_weights.at(i) << "\n";
                    std::cout << "expected weight: " << expected_weights.at(i) << "\n";
                }
                else {
                    n_i = i - num_nodes;
                    std::cout << "node index: " << tree.get_node(n_i)->get_index() << "\n";
                    std::cout << "node label: " << tree.get_node(n_i)->get_label() << "\n";
                    std::cout << "node weight: " << ln_weights.at(i) << "\n";
                    std::cout << "expected weight: " << expected_weights.at(i) << "\n";
                }
                REQUIRE(ln_weights.at(i) == expected_weights.at(i));
            }
        }
    }
}

TEST_CASE("Testing SubtreePruneRegraftRevJumpSampler with 3 leaves and free root",
        "[SubtreePruneRegraftRevJumpSampler]") {

    SECTION("Testing 3 leaves with free root") {
        RandomNumberGenerator rng = RandomNumberGenerator(18);

        double root_ht = 0.2;
        std::shared_ptr<Node> root = std::make_shared<Node>("root", root_ht);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.1);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);

        tree.ignore_data();

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_root_height();

        SubtreePruneRegraftRevJumpSampler< BaseTree<Node> > op;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;

        unsigned int count_012 = 0;
        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_12 = 0;

        unsigned int niterations = 2000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            // std::cout << i << " " << tree.to_parentheses() << "\n";
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                if (tree.get_root_ptr()->get_number_of_children() == 3) {
                    ++count_012;
                }
                else {
                    if (tree.get_root_ptr()->is_child("leaf0")) {
                        ++count_12;
                    }
                    if (tree.get_root_ptr()->is_child("leaf1")) {
                        ++count_02;
                    }
                    if (tree.get_root_ptr()->is_child("leaf2")) {
                        ++count_01;
                    }
                    root_height_summary.add_sample(tree.get_root_height());
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        REQUIRE((count_01 + count_02 + count_12 + count_012) == nsamples);
        REQUIRE(root_height_summary.sample_size() == (count_01 + count_02 + count_12));

        double freq_012 = count_012 / (double)nsamples;
        double freq_01 = count_01 / (double)nsamples;
        double freq_02 = count_02 / (double)nsamples;
        double freq_12 = count_12 / (double)nsamples;
        std::cout << "Freq of (0,1,2): " << freq_012 << "\n";
        std::cout << "Freq of ((0,1),2): " << freq_01 << "\n";
        std::cout << "Freq of ((0,2),1): " << freq_02 << "\n";
        std::cout << "Freq of ((1,2),0): " << freq_12 << "\n";

        std::vector<unsigned int> counts {count_012, count_01, count_02, count_12};
        write_r_script(counts, "../3-leaf-general-tree-spr-rj-test.r");

        double eps = 0.002;

        REQUIRE(freq_012 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_01 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_02 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_12 == Approx(0.25).epsilon(eps));
        
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing SubtreePruneRegraftRevJumpSampler with 3 leaves, fixed root and operate_plus",
        "[SubtreePruneRegraftRevJumpSampler]") {

    SECTION("Testing 3 leaves, free root and operate_plus") {
        RandomNumberGenerator rng = RandomNumberGenerator(81692135456);

        double root_ht = 0.2;
        std::shared_ptr<Node> root = std::make_shared<Node>("root", root_ht);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.1);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);

        tree.ignore_data();

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_root_height();

        SubtreePruneRegraftRevJumpSampler< BaseTree<Node> > op;
        std::shared_ptr< NodeHeightScaler< BaseTree<Node> > > node_height_op = std::make_shared<NodeHeightScaler< BaseTree<Node> > >();
        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate< BaseTree<Node> > > > other_ops;
        other_ops.push_back(node_height_op);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;

        unsigned int count_012 = 0;
        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_12 = 0;

        unsigned int niterations = 2000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            // std::cout << tree.to_parentheses() << "\n";
            op.operate_plus(rng, &tree, other_ops, 1, 2);
            if ((i + 1) % sample_freq == 0) {
                if (tree.get_root_ptr()->get_number_of_children() == 3) {
                    ++count_012;
                }
                else {
                    if (tree.get_root_ptr()->is_child("leaf0")) {
                        ++count_12;
                    }
                    if (tree.get_root_ptr()->is_child("leaf1")) {
                        ++count_02;
                    }
                    if (tree.get_root_ptr()->is_child("leaf2")) {
                        ++count_01;
                    }
                    root_height_summary.add_sample(tree.get_root_height());
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << node_height_op->to_string();

        REQUIRE((count_01 + count_02 + count_12 + count_012) == nsamples);
        REQUIRE(root_height_summary.sample_size() == (count_01 + count_02 + count_12));

        double freq_012 = count_012 / (double)nsamples;
        double freq_01 = count_01 / (double)nsamples;
        double freq_02 = count_02 / (double)nsamples;
        double freq_12 = count_12 / (double)nsamples;
        std::cout << "Freq of (0,1,2): " << freq_012 << "\n";
        std::cout << "Freq of ((0,1),2): " << freq_01 << "\n";
        std::cout << "Freq of ((0,2),1): " << freq_02 << "\n";
        std::cout << "Freq of ((1,2),0): " << freq_12 << "\n";

        double eps = 0.002;

        REQUIRE(freq_012 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_01 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_02 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_12 == Approx(0.25).epsilon(eps));
        
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing SubtreePruneRegraftRevJumpSampler with 4 leaves and estimated root",
        "[SubtreePruneRegraftRevJumpSampler]") {

    SECTION("Testing 4 leaves with estimated root") {
        RandomNumberGenerator rng = RandomNumberGenerator(20);

        double root_ht = 0.5;
        std::shared_ptr<Node> root = std::make_shared<Node>(4, "root", root_ht);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);

        root->add_child(leaf0);
        root->add_child(leaf1);
        root->add_child(leaf2);
        root->add_child(leaf3);

        BaseTree<Node> tree(root);

        tree.ignore_data();

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_root_height();

        SubtreePruneRegraftRevJumpSampler< BaseTree<Node> > op;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        unsigned int count_0123 = 0;
        unsigned int count_0_ = 0;
        unsigned int count_1_ = 0;
        unsigned int count_2_ = 0;
        unsigned int count_3_ = 0;
        unsigned int count_01_ = 0;
        unsigned int count_02_ = 0;
        unsigned int count_03_ = 0;
        unsigned int count_12_ = 0;
        unsigned int count_13_ = 0;
        unsigned int count_23_ = 0;
        unsigned int count_01_23 = 0;
        unsigned int count_02_13 = 0;
        unsigned int count_03_12 = 0;
        unsigned int count_gen_01_23 = 0;
        unsigned int count_gen_02_13 = 0;
        unsigned int count_gen_03_12 = 0;
        unsigned int count_gen_01_2_3 = 0;
        unsigned int count_gen_01_3_2 = 0;
        unsigned int count_gen_02_1_3 = 0;
        unsigned int count_gen_02_3_1 = 0;
        unsigned int count_gen_03_1_2 = 0;
        unsigned int count_gen_03_2_1 = 0;
        unsigned int count_gen_12_0_3 = 0;
        unsigned int count_gen_12_3_0 = 0;
        unsigned int count_gen_13_0_2 = 0;
        unsigned int count_gen_13_2_0 = 0;
        unsigned int count_gen_23_0_1 = 0;
        unsigned int count_gen_23_1_0 = 0;
        unsigned int count_3_heights = 0;
        unsigned int count_2_heights = 0;
        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int niterations = 5000000;
        unsigned int sample_freq = 20;
        unsigned int nsamples = niterations / sample_freq;
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

                if (tree.get_number_of_node_heights() == 1) {
                    ++count_0123;
                    REQUIRE(tree.get_root_ptr()->get_number_of_children() == 4);
                }
                else if (tree.get_number_of_node_heights() == 2) {
                    ++count_2_heights;
                    if (tree.get_root_ptr()->get_number_of_children() == 2) {
                        if (tree.get_root_ptr()->get_child(0)->is_leaf() ||
                                tree.get_root_ptr()->get_child(1)->is_leaf()) {
                            if(tree.get_root_ptr()->is_child("leaf0")) {
                                ++count_0_;
                            }
                            else if(tree.get_root_ptr()->is_child("leaf1")) {
                                ++count_1_;
                            }
                            else if(tree.get_root_ptr()->is_child("leaf2")) {
                                ++count_2_;
                            }
                            else if(tree.get_root_ptr()->is_child("leaf3")) {
                                ++count_3_;
                            }
                            else {
                                REQUIRE(0 == 1);
                            }
                        }
                        else {
                            if (
                                    (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(0)->is_child("leaf1"))
                                    ||
                                    (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                               ) {
                                ++count_01_23;
                            }
                            else if (
                                    (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(0)->is_child("leaf2"))
                                    ||
                                    (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                               ) {
                                ++count_02_13;
                            }
                            else if (
                                    (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(0)->is_child("leaf3"))
                                    ||
                                    (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                               ) {
                                ++count_03_12;
                            }
                            else {
                                REQUIRE(0 == 1);
                            }
                        }
                    }
                    else if (tree.get_root_ptr()->get_number_of_children() == 3) {
                        if (
                                tree.get_root_ptr()->is_child("leaf0") &&
                                tree.get_root_ptr()->is_child("leaf1")
                            )
                        {
                            ++count_23_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf0") &&
                                tree.get_root_ptr()->is_child("leaf2")
                            )
                        {
                            ++count_13_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf0") &&
                                tree.get_root_ptr()->is_child("leaf3")
                            )
                        {
                            ++count_12_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf1") &&
                                tree.get_root_ptr()->is_child("leaf2")
                            )
                        {
                            ++count_03_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf1") &&
                                tree.get_root_ptr()->is_child("leaf3")
                            )
                        {
                            ++count_02_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf2") &&
                                tree.get_root_ptr()->is_child("leaf3")
                            )
                        {
                            ++count_01_;

                        }
                        else {
                            REQUIRE(0 == 1);
                        }
                    }
                    else {
                        REQUIRE(0 == 1);
                    }
                }
                else if (tree.get_number_of_node_heights() == 3) {
                    ++count_3_heights;
                    if ((! tree.get_root_ptr()->get_child(0)->is_leaf()) &&
                        (! tree.get_root_ptr()->get_child(1)->is_leaf())) {
                        if (
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(0)->is_child("leaf1"))
                                ||
                                (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                           ) {
                            ++count_gen_01_23;
                        }
                        else if (
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(0)->is_child("leaf2"))
                                ||
                                (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                           ) {
                            ++count_gen_02_13;
                        }
                        else if (
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(0)->is_child("leaf3"))
                                ||
                                (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                           ) {
                            ++count_gen_03_12;
                        }
                        else {
                            REQUIRE(0 == 1);
                        }
                    }
                    else {
                        // general ladderized topology
                        if (tree.get_root_ptr()->is_child("leaf3") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf2") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                            )
                        {
                            ++count_gen_01_2_3;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf2") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf3") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                            )
                        {
                            ++count_gen_01_3_2;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf3") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf1") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                            )
                        {
                            ++count_gen_02_1_3;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf1") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf3") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                            )
                        {
                            ++count_gen_02_3_1;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf2") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf1") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                            )
                        {
                            ++count_gen_03_1_2;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf1") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf2") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                            )
                        {
                            ++count_gen_03_2_1;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf3") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf0"))
                            )
                        {
                            ++count_gen_12_0_3;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf0") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf3") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                            )
                        {
                            ++count_gen_12_3_0;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf2") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf0"))
                            )
                        {
                            ++count_gen_13_0_2;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf0") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf2") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                            )
                        {
                            ++count_gen_13_2_0;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf1") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf0"))
                            )
                        {
                            ++count_gen_23_0_1;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf0") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf1") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                            )
                        {
                            ++count_gen_23_1_0;
                        }
                        else {
                            REQUIRE(0 == 1);
                        }
                    }
                }
                else {
                    REQUIRE(0 == 1);
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        REQUIRE((count_0123 + count_2_heights + count_3_heights) == nsamples);
        REQUIRE((count_01_23 +
                count_02_13 +
                count_03_12 +
                count_0_ +
                count_1_ +
                count_2_ +
                count_3_ +
                count_01_ +
                count_02_ +
                count_03_ +
                count_12_ +
                count_13_ +
                count_23_) == count_2_heights);
        REQUIRE((count_gen_01_23 +
                count_gen_02_13 +
                count_gen_03_12 +
                count_gen_01_2_3 +
                count_gen_01_3_2 +
                count_gen_02_1_3 +
                count_gen_02_3_1 +
                count_gen_03_1_2 +
                count_gen_03_2_1 +
                count_gen_12_0_3 +
                count_gen_12_3_0 +
                count_gen_13_0_2 +
                count_gen_13_2_0 +
                count_gen_23_0_1 +
                count_gen_23_1_0) == count_3_heights);

        double freq_2_heights = count_2_heights / (double)nsamples;
        double freq_3_heights = count_3_heights / (double)nsamples;
        double freq_0123 = count_0123 / (double)nsamples;
        double freq_01_23 = count_01_23 / (double)nsamples;
        double freq_02_13 = count_02_13 / (double)nsamples;
        double freq_03_12 = count_03_12 / (double)nsamples;
        double freq_gen_01_23 = count_gen_01_23 / (double)nsamples;
        double freq_gen_02_13 = count_gen_02_13 / (double)nsamples;
        double freq_gen_03_12 = count_gen_03_12 / (double)nsamples;
        double freq_0_ = count_0_ / (double)nsamples;
        double freq_1_ = count_1_ / (double)nsamples;
        double freq_2_ = count_2_ / (double)nsamples;
        double freq_3_ = count_3_ / (double)nsamples;
        double freq_01_ = count_01_ / (double)nsamples;
        double freq_02_ = count_02_ / (double)nsamples;
        double freq_03_ = count_03_ / (double)nsamples;
        double freq_12_ = count_12_ / (double)nsamples;
        double freq_13_ = count_13_ / (double)nsamples;
        double freq_23_ = count_23_ / (double)nsamples;
        double freq_gen_01_2_3 = count_gen_01_2_3 / (double)nsamples;
        double freq_gen_01_3_2 = count_gen_01_3_2 / (double)nsamples;
        double freq_gen_02_1_3 = count_gen_02_1_3 / (double)nsamples;
        double freq_gen_02_3_1 = count_gen_02_3_1 / (double)nsamples;
        double freq_gen_03_1_2 = count_gen_03_1_2 / (double)nsamples;
        double freq_gen_03_2_1 = count_gen_03_2_1 / (double)nsamples;
        double freq_gen_12_0_3 = count_gen_12_0_3 / (double)nsamples;
        double freq_gen_12_3_0 = count_gen_12_3_0 / (double)nsamples;
        double freq_gen_13_0_2 = count_gen_13_0_2 / (double)nsamples;
        double freq_gen_13_2_0 = count_gen_13_2_0 / (double)nsamples;
        double freq_gen_23_0_1 = count_gen_23_0_1 / (double)nsamples;
        double freq_gen_23_1_0 = count_gen_23_1_0 / (double)nsamples;

        double exp_freq = 1.0/29.0;

        std::cout << "Freq of (0,1,2,3): " << freq_0123 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((0,1),2,3): " << freq_01_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((0,2),1,3): " << freq_02_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((0,3),1,2): " << freq_03_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((1,2),0,3): " << freq_12_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((1,3),0,2): " << freq_13_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((2,3),0,1): " << freq_23_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of (0,(1,2,3)): " << freq_0_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of (1,(0,2,3)): " << freq_1_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of (2,(1,0,3)): " << freq_2_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of (3,(1,2,0)): " << freq_3_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of shared ((0,1),(2,3)): " << freq_01_23 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of shared ((0,2),(1,3)): " << freq_02_13 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of shared ((0,3),(1,2)): " << freq_03_12 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen ((0,1),(2,3)): " << freq_gen_01_23 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen ((0,2),(1,3)): " << freq_gen_02_13 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen ((0,3),(1,2)): " << freq_gen_03_12 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,1),2),3): " << freq_gen_01_2_3 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,1),3),2): " << freq_gen_01_3_2 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,2),1),3): " << freq_gen_02_1_3 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,2),3),1): " << freq_gen_02_3_1 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,3),1),2): " << freq_gen_03_1_2 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,3),2),1): " << freq_gen_03_2_1 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((1,2),0),3): " << freq_gen_12_0_3 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((1,2),3),0): " << freq_gen_12_3_0 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((1,3),0),2): " << freq_gen_13_0_2 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((1,3),2),0): " << freq_gen_13_2_0 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((2,3),0),1): " << freq_gen_23_0_1 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((2,3),1),0): " << freq_gen_23_1_0 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of 2 heights: " << freq_2_heights << " (expected " << 13 * exp_freq << ")\n";
        std::cout << "Freq of 3 heights: " << freq_3_heights << " (expected " << 15 * exp_freq << ")\n";

        write_r_script(split_counts, 4, "../4-leaf-general-tree-spr-rj-test.r");

        double eps = 0.001;

        REQUIRE(freq_2_heights == Approx(13 * exp_freq).epsilon(eps));
        REQUIRE(freq_3_heights == Approx(15 * exp_freq).epsilon(eps));

        REQUIRE(freq_0123 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_01_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_02_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_03_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_12_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_13_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_23_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_0_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_1_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_2_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_3_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_01_23 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_02_13 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_03_12 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_01_23 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_02_13 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_03_12 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_01_2_3 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_01_3_2 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_02_1_3 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_02_3_1 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_03_1_2 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_03_2_1 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_12_0_3 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_12_3_0 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_13_0_2 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_13_2_0 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_23_0_1 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_23_1_0 == Approx(exp_freq).epsilon(eps));
    }
}

TEST_CASE("Testing SubtreePruneRegraftRevJumpSampler with 4 leaves, estimated root, and operate_plus",
        "[SubtreePruneRegraftRevJumpSampler]") {

    SECTION("Testing 4 leaves with estimated root and operate_plus") {
        RandomNumberGenerator rng = RandomNumberGenerator(21);

        double root_ht = 0.5;
        std::shared_ptr<Node> root = std::make_shared<Node>("root", root_ht);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);

        root->add_child(leaf0);
        root->add_child(leaf1);
        root->add_child(leaf2);
        root->add_child(leaf3);

        BaseTree<Node> tree(root);

        tree.ignore_data();

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_root_height();

        SubtreePruneRegraftRevJumpSampler< BaseTree<Node> > op;
        std::shared_ptr< NodeHeightScaler< BaseTree<Node> > > node_height_op = std::make_shared<NodeHeightScaler< BaseTree<Node> > >();
        /* std::shared_ptr< NeighborHeightNodeSwap< BasTree<Node> > > node_swap_op = std::make_shared<NeighborHeightNodeSwap< BaseTree<Node> > >(); */
        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate< BaseTree<Node> > > > other_ops;
        other_ops.push_back(node_height_op);
        /* other_ops.push_back(node_swap_op); */

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        unsigned int count_0123 = 0;
        unsigned int count_0_ = 0;
        unsigned int count_1_ = 0;
        unsigned int count_2_ = 0;
        unsigned int count_3_ = 0;
        unsigned int count_01_ = 0;
        unsigned int count_02_ = 0;
        unsigned int count_03_ = 0;
        unsigned int count_12_ = 0;
        unsigned int count_13_ = 0;
        unsigned int count_23_ = 0;
        unsigned int count_01_23 = 0;
        unsigned int count_02_13 = 0;
        unsigned int count_03_12 = 0;
        unsigned int count_gen_01_23 = 0;
        unsigned int count_gen_02_13 = 0;
        unsigned int count_gen_03_12 = 0;
        unsigned int count_gen_01_2_3 = 0;
        unsigned int count_gen_01_3_2 = 0;
        unsigned int count_gen_02_1_3 = 0;
        unsigned int count_gen_02_3_1 = 0;
        unsigned int count_gen_03_1_2 = 0;
        unsigned int count_gen_03_2_1 = 0;
        unsigned int count_gen_12_0_3 = 0;
        unsigned int count_gen_12_3_0 = 0;
        unsigned int count_gen_13_0_2 = 0;
        unsigned int count_gen_13_2_0 = 0;
        unsigned int count_gen_23_0_1 = 0;
        unsigned int count_gen_23_1_0 = 0;
        unsigned int count_3_heights = 0;
        unsigned int count_2_heights = 0;

        unsigned int niterations = 5000000;
        unsigned int sample_freq = 20;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate_plus(rng, &tree, other_ops, 1, 2, 2);
            if ((i + 1) % sample_freq == 0) {
                /* std::cout << "prior: " << tree.get_log_prior_density_value() << "\n"; */
                if (tree.get_number_of_node_heights() == 1) {
                    ++count_0123;
                    REQUIRE(tree.get_root_ptr()->get_number_of_children() == 4);
                }
                else if (tree.get_number_of_node_heights() == 2) {
                    ++count_2_heights;
                    if (tree.get_root_ptr()->get_number_of_children() == 2) {
                        if (tree.get_root_ptr()->get_child(0)->is_leaf() ||
                                tree.get_root_ptr()->get_child(1)->is_leaf()) {
                            if(tree.get_root_ptr()->is_child("leaf0")) {
                                ++count_0_;
                            }
                            else if(tree.get_root_ptr()->is_child("leaf1")) {
                                ++count_1_;
                            }
                            else if(tree.get_root_ptr()->is_child("leaf2")) {
                                ++count_2_;
                            }
                            else if(tree.get_root_ptr()->is_child("leaf3")) {
                                ++count_3_;
                            }
                            else {
                                REQUIRE(0 == 1);
                            }
                        }
                        else {
                            if (
                                    (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(0)->is_child("leaf1"))
                                    ||
                                    (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                               ) {
                                ++count_01_23;
                            }
                            else if (
                                    (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(0)->is_child("leaf2"))
                                    ||
                                    (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                               ) {
                                ++count_02_13;
                            }
                            else if (
                                    (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(0)->is_child("leaf3"))
                                    ||
                                    (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                               ) {
                                ++count_03_12;
                            }
                            else {
                                REQUIRE(0 == 1);
                            }
                        }
                    }
                    else if (tree.get_root_ptr()->get_number_of_children() == 3) {
                        if (
                                tree.get_root_ptr()->is_child("leaf0") &&
                                tree.get_root_ptr()->is_child("leaf1")
                            )
                        {
                            ++count_23_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf0") &&
                                tree.get_root_ptr()->is_child("leaf2")
                            )
                        {
                            ++count_13_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf0") &&
                                tree.get_root_ptr()->is_child("leaf3")
                            )
                        {
                            ++count_12_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf1") &&
                                tree.get_root_ptr()->is_child("leaf2")
                            )
                        {
                            ++count_03_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf1") &&
                                tree.get_root_ptr()->is_child("leaf3")
                            )
                        {
                            ++count_02_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf2") &&
                                tree.get_root_ptr()->is_child("leaf3")
                            )
                        {
                            ++count_01_;

                        }
                        else {
                            REQUIRE(0 == 1);
                        }
                    }
                    else {
                        REQUIRE(0 == 1);
                    }
                }
                else if (tree.get_number_of_node_heights() == 3) {
                    ++count_3_heights;
                    if ((! tree.get_root_ptr()->get_child(0)->is_leaf()) &&
                        (! tree.get_root_ptr()->get_child(1)->is_leaf())) {
                        if (
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(0)->is_child("leaf1"))
                                ||
                                (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                           ) {
                            ++count_gen_01_23;
                        }
                        else if (
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(0)->is_child("leaf2"))
                                ||
                                (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                           ) {
                            ++count_gen_02_13;
                        }
                        else if (
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(0)->is_child("leaf3"))
                                ||
                                (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                           ) {
                            ++count_gen_03_12;
                        }
                        else {
                            REQUIRE(0 == 1);
                        }
                    }
                    else {
                        // general ladderized topology
                        if (tree.get_root_ptr()->is_child("leaf3") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf2") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                            )
                        {
                            ++count_gen_01_2_3;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf2") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf3") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                            )
                        {
                            ++count_gen_01_3_2;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf3") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf1") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                            )
                        {
                            ++count_gen_02_1_3;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf1") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf3") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                            )
                        {
                            ++count_gen_02_3_1;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf2") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf1") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                            )
                        {
                            ++count_gen_03_1_2;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf1") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf2") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                            )
                        {
                            ++count_gen_03_2_1;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf3") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf0"))
                            )
                        {
                            ++count_gen_12_0_3;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf0") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf3") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                            )
                        {
                            ++count_gen_12_3_0;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf2") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf0"))
                            )
                        {
                            ++count_gen_13_0_2;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf0") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf2") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                            )
                        {
                            ++count_gen_13_2_0;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf1") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf0"))
                            )
                        {
                            ++count_gen_23_0_1;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf0") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf1") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                            )
                        {
                            ++count_gen_23_1_0;
                        }
                        else {
                            REQUIRE(0 == 1);
                        }
                    }
                }
                else {
                    REQUIRE(0 == 1);
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << node_height_op->to_string();
        /* std::cout << node_swap_op->to_string(); */


        REQUIRE((count_0123 + count_2_heights + count_3_heights) == nsamples);
        REQUIRE((count_01_23 +
                count_02_13 +
                count_03_12 +
                count_0_ +
                count_1_ +
                count_2_ +
                count_3_ +
                count_01_ +
                count_02_ +
                count_03_ +
                count_12_ +
                count_13_ +
                count_23_) == count_2_heights);
        REQUIRE((count_gen_01_23 +
                count_gen_02_13 +
                count_gen_03_12 +
                count_gen_01_2_3 +
                count_gen_01_3_2 +
                count_gen_02_1_3 +
                count_gen_02_3_1 +
                count_gen_03_1_2 +
                count_gen_03_2_1 +
                count_gen_12_0_3 +
                count_gen_12_3_0 +
                count_gen_13_0_2 +
                count_gen_13_2_0 +
                count_gen_23_0_1 +
                count_gen_23_1_0) == count_3_heights);

        double freq_2_heights = count_2_heights / (double)nsamples;
        double freq_3_heights = count_3_heights / (double)nsamples;
        double freq_0123 = count_0123 / (double)nsamples;
        double freq_01_23 = count_01_23 / (double)nsamples;
        double freq_02_13 = count_02_13 / (double)nsamples;
        double freq_03_12 = count_03_12 / (double)nsamples;
        double freq_gen_01_23 = count_gen_01_23 / (double)nsamples;
        double freq_gen_02_13 = count_gen_02_13 / (double)nsamples;
        double freq_gen_03_12 = count_gen_03_12 / (double)nsamples;
        double freq_0_ = count_0_ / (double)nsamples;
        double freq_1_ = count_1_ / (double)nsamples;
        double freq_2_ = count_2_ / (double)nsamples;
        double freq_3_ = count_3_ / (double)nsamples;
        double freq_01_ = count_01_ / (double)nsamples;
        double freq_02_ = count_02_ / (double)nsamples;
        double freq_03_ = count_03_ / (double)nsamples;
        double freq_12_ = count_12_ / (double)nsamples;
        double freq_13_ = count_13_ / (double)nsamples;
        double freq_23_ = count_23_ / (double)nsamples;
        double freq_gen_01_2_3 = count_gen_01_2_3 / (double)nsamples;
        double freq_gen_01_3_2 = count_gen_01_3_2 / (double)nsamples;
        double freq_gen_02_1_3 = count_gen_02_1_3 / (double)nsamples;
        double freq_gen_02_3_1 = count_gen_02_3_1 / (double)nsamples;
        double freq_gen_03_1_2 = count_gen_03_1_2 / (double)nsamples;
        double freq_gen_03_2_1 = count_gen_03_2_1 / (double)nsamples;
        double freq_gen_12_0_3 = count_gen_12_0_3 / (double)nsamples;
        double freq_gen_12_3_0 = count_gen_12_3_0 / (double)nsamples;
        double freq_gen_13_0_2 = count_gen_13_0_2 / (double)nsamples;
        double freq_gen_13_2_0 = count_gen_13_2_0 / (double)nsamples;
        double freq_gen_23_0_1 = count_gen_23_0_1 / (double)nsamples;
        double freq_gen_23_1_0 = count_gen_23_1_0 / (double)nsamples;

        double exp_freq = 1.0/29.0;

        std::cout << "Freq of (0,1,2,3): " << freq_0123 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((0,1),2,3): " << freq_01_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((0,2),1,3): " << freq_02_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((0,3),1,2): " << freq_03_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((1,2),0,3): " << freq_12_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((1,3),0,2): " << freq_13_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((2,3),0,1): " << freq_23_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of (0,(1,2,3)): " << freq_0_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of (1,(0,2,3)): " << freq_1_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of (2,(1,0,3)): " << freq_2_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of (3,(1,2,0)): " << freq_3_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of shared ((0,1),(2,3)): " << freq_01_23 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of shared ((0,2),(1,3)): " << freq_02_13 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of shared ((0,3),(1,2)): " << freq_03_12 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen ((0,1),(2,3)): " << freq_gen_01_23 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen ((0,2),(1,3)): " << freq_gen_02_13 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen ((0,3),(1,2)): " << freq_gen_03_12 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,1),2),3): " << freq_gen_01_2_3 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,1),3),2): " << freq_gen_01_3_2 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,2),1),3): " << freq_gen_02_1_3 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,2),3),1): " << freq_gen_02_3_1 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,3),1),2): " << freq_gen_03_1_2 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,3),2),1): " << freq_gen_03_2_1 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((1,2),0),3): " << freq_gen_12_0_3 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((1,2),3),0): " << freq_gen_12_3_0 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((1,3),0),2): " << freq_gen_13_0_2 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((1,3),2),0): " << freq_gen_13_2_0 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((2,3),0),1): " << freq_gen_23_0_1 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((2,3),1),0): " << freq_gen_23_1_0 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of 2 heights: " << freq_2_heights << " (expected " << 13 * exp_freq << ")\n";
        std::cout << "Freq of 3 heights: " << freq_3_heights << " (expected " << 15 * exp_freq << ")\n";

        double eps = 0.001;

        REQUIRE(freq_2_heights == Approx(13 * exp_freq).epsilon(eps));
        REQUIRE(freq_3_heights == Approx(15 * exp_freq).epsilon(eps));

        REQUIRE(freq_0123 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_01_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_02_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_03_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_12_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_13_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_23_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_0_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_1_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_2_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_3_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_01_23 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_02_13 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_03_12 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_01_23 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_02_13 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_03_12 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_01_2_3 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_01_3_2 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_02_1_3 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_02_3_1 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_03_1_2 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_03_2_1 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_12_0_3 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_12_3_0 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_13_0_2 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_13_2_0 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_23_0_1 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_23_1_0 == Approx(exp_freq).epsilon(eps));
    }
}


// Expectations for tree with 5 leaves:
// ----------------------------------------------------------------------------
// # of unlabeled                            # of trees (labelings)      
//  topologies
// ----------------------------------------------------------------------------
// (a,b,c,d,e)                                                     = 1
// (a,(b,c,d,e))           = 5 choose 4                            = 5
// (a,b,(c,d,e))           = 5 choose 3                            = 10 
// (a,b,c,(d,e))           = 5 choose 2                            = 10
// !((a,b,c),(d,e))        = 5 choose 3                            = 10
// (a,(b,(c,d,e)))         = (5 choose 3) * 2!                     = 20 
// !((a,b),(c,d),e)        = ((5 choose 2) * (3 choose 2)) / 2     = 15
// !(((a,b),(c,d)),e       = ((5 choose 2) * (3 choose 2)) / 2     = 15
// !(((a,b),c),(d,e))      = (5 choose 2) * (3 choose 2)           = 30
// (a,(b,(c,(d,e))))       = (5 choose 2) * 3!                     = 60
// (a,b,(c,(d,e)))         = (5 choose 2) * (3 choose 2)           = 30
// (a,(b,c,(d,e)))         = (5 choose 2) * (3 choose 2)           = 30
// ----------------------------------------------------------------------------
// TOTAL                                                           = 236
//
// 236 matches Felsenstein 1978, but we need to account for shared node
// heights. The topologies above prefixed with '!' are topologies that have
// potentially shared node heights. For each shared node configuration of these
// topologies, we have to add that many additional trees, which we do below.
//
// ----------------------------------------------------------------------------
// # of unlabeled                            # of trees (labelings)      
//  topologies
// ----------------------------------------------------------------------------
// ((a,b,c)*,(d,e)*)       = 5 choose 3                            = 10
// ((a,b)*,(c,d)*,e)       = ((5 choose 2) * (3 choose 2)) / 2     = 15
// (((a,b)*,(c,d)*),e      = ((5 choose 2) * (3 choose 2)) / 2     = 15
// (((a,b),c)*,(d,e)*)     = (5 choose 2) * (3 choose 2)           = 30
// (((a,b)*,c),(d,e)*)     = (5 choose 2) * (3 choose 2)           = 30
// ----------------------------------------------------------------------------
// GRAND TOTAL # OF TREE MODELS                                    = 336
//
// The asterisks in the topologies above indicated shared node heights.
TEST_CASE("Testing SubtreePruneRegraftRevJumpSampler with 5 leaves and estimated root",
        "[SubtreePruneRegraftRevJumpSampler]") {

    SECTION("Testing 5 leaves with estimated root") {
        RandomNumberGenerator rng = RandomNumberGenerator(9734598374);

        double root_ht = 0.5;
        std::shared_ptr<Node> root = std::make_shared<Node>(5, "root", root_ht);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);

        root->add_child(leaf0);
        root->add_child(leaf1);
        root->add_child(leaf2);
        root->add_child(leaf3);
        root->add_child(leaf4);

        BaseTree<Node> tree(root);

        tree.ignore_data();

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_root_height();

        SubtreePruneRegraftRevJumpSampler< BaseTree<Node> > op;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int count_nheights_1 = 0;
        unsigned int count_nheights_2 = 0;
        unsigned int count_nheights_3 = 0;
        unsigned int count_nheights_4 = 0;

        unsigned int niterations = 50000000;
        unsigned int sample_freq = 50;
        unsigned int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 10000;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
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
                ++sample_count;
                if (sample_count % report_freq == 0) {
                    std::cout << "Sampled " << sample_count << " of " << nsamples << std::endl;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        REQUIRE((count_nheights_1 + count_nheights_2 + count_nheights_3 + count_nheights_4) == nsamples);

        double freq_nheights_1 = count_nheights_1 / (double)nsamples;
        double freq_nheights_2 = count_nheights_2 / (double)nsamples;
        double freq_nheights_3 = count_nheights_3 / (double)nsamples;
        double freq_nheights_4 = count_nheights_4 / (double)nsamples;

        double exp_freq = 1.0/336.0;
        double exp_count = nsamples/336.0;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.1;
        unsigned int total_trees_sampled = 0;
        std::map< std::set< std::set<Split> >, double> split_freqs;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            split_freqs[s_c.first] = s_c.second / (double)nsamples;
            std::cout << "Tree:\n";
            for (auto splitset : s_c.first) {
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
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            std::cout << "  nsamples: " << s_c.second << "\n";
            std::cout << "  prop error: " << prop_error << "\n";
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }

        double quantile_chi_sq_335_10 = 368.6;
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        std::cout << "Chi-square(335) 0.9 quantile: " << quantile_chi_sq_335_10 << "\n";

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

        write_r_script(split_counts, 5, "../5-leaf-general-tree-spr-rj-test.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        REQUIRE(split_counts.size() == 336);

        double eps = 0.001;

        for (auto s_f : split_freqs) {
            REQUIRE(s_f.second == Approx(exp_freq).epsilon(eps));
        }

        //REQUIRE(chi_sq_test_statistic < quantile_chi_sq_335_10);
    }
}

TEST_CASE("Testing SubtreePruneRegraftRevJumpSampler with 6 leaves and estimated root",
        "[SubtreePruneRegraftRevJumpSampler]") {

    SECTION("Testing 6 leaves with estimated root") {
        RandomNumberGenerator rng = RandomNumberGenerator(37495738947);

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

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_root_height();

        SubtreePruneRegraftRevJumpSampler< BaseTree<Node> > op;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int niterations = 1000000000;
        unsigned int sample_freq = 100;
        unsigned int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 100000;
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
            std::cout << "Tree:\n";
            for (auto splitset : s_c.first) {
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
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            std::cout << "  nsamples: " << s_c.second << "\n";
            std::cout << "  prop error: " << prop_error << "\n";
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

        write_r_script(split_counts, 6, "../6-leaf-general-tree-spr-rj-test.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        // REQUIRE(split_counts.size() == ???);

        // REQUIRE(chi_sq_test_statistic < quantile_chi_sq_5627_10);
    }
}

TEST_CASE("Testing SubtreePruneRegraftRevJumpSampler with 7 leaves and estimated root",
        "[SubtreePruneRegraftRevJumpSampler]") {

    SECTION("Testing 7 leaves with estimated root") {
        RandomNumberGenerator rng = RandomNumberGenerator(874668437);

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

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_root_height();

        SubtreePruneRegraftRevJumpSampler< BaseTree<Node> > op;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned long long int niterations = 5000000000;
        unsigned int sample_freq = 100;
        unsigned long long int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 100000;
        for (unsigned long long int i = 0; i < niterations; ++i) {
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
            std::cout << "Tree:\n";
            for (auto splitset : s_c.first) {
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
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            std::cout << "  nsamples: " << s_c.second << "\n";
            std::cout << "  prop error: " << prop_error << "\n";
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

        write_r_script(split_counts, 7, "../7-leaf-general-tree-spr-rj-test.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        // REQUIRE(split_counts.size() == ???);

        /* REQUIRE(chi_sq_test_statistic < quantile_chi_sq_5627_10); */
    }
}

TEST_CASE("Testing SubtreePruneRegraftRevJumpSampler with BasePopulationTree, 5 leaves, full model, unconstrained sizes",
        "[SubtreePruneRegraftRevJumpSampler]") {

    SECTION("Testing 5 leaves with BasePopulationTree, full model, unconstrained sizes") {
        RandomNumberGenerator rng = RandomNumberGenerator(214354584);

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

        double height_beta_shape = 20.0;
        double height_beta_scale = 0.3;
        std::shared_ptr<ContinuousProbabilityDistribution> height_beta_prior = std::make_shared<GammaDistribution>(
                height_beta_shape,
                height_beta_scale);

        double root_ht = 0.5;
        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(5, "root", root_ht);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);
        std::shared_ptr<PopulationNode> leaf3 = std::make_shared<PopulationNode>(3, "leaf3", 0.0);
        std::shared_ptr<PopulationNode> leaf4 = std::make_shared<PopulationNode>(4, "leaf4", 0.0);

        root->add_child(leaf0);
        root->add_child(leaf1);
        root->add_child(leaf2);
        root->add_child(leaf3);
        root->add_child(leaf4);

        BasePopulationTree tree(root);
        tree.ignore_data();

        tree.set_mutation_rate_prior(mu_rate_prior);
        tree.set_root_node_height_prior(root_height_prior);
        tree.set_population_size_prior(pop_size_prior);
        tree.set_freq_1_prior(freq_prior);
        tree.set_prior_on_alpha_of_node_height_beta_prior(height_alpha_prior);
        tree.set_prior_on_beta_of_node_height_beta_prior(height_beta_prior);

        tree.estimate_mutation_rate();
        tree.estimate_root_height();
        tree.estimate_state_frequencies();
        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();

        GeneralTreeOperatorSchedule< BasePopulationTree > op_schedule;
        std::shared_ptr< GeneralTreeOperatorTemplate< BasePopulationTree > > op;

        op = std::make_shared< SubtreePruneRegraftRevJumpSampler<BasePopulationTree> >(12.0);
        op_schedule.add_operator(op);
        op = std::make_shared< MuRateScaler >(1.0);
        op_schedule.add_operator(op);
        op = std::make_shared< RootHeightScaler<BasePopulationTree> >(2.0);
        op_schedule.add_operator(op);
        op = std::make_shared< RootHeightSizeMixer >(1.0);
        op_schedule.add_operator(op);
        op = std::make_shared< HeightSizeMixer >(1.0);
        op_schedule.add_operator(op);
        op = std::make_shared< HeightSizeSlideBumpMixer >(1.0);
        op_schedule.add_operator(op);
        op = std::make_shared< NodeHeightScaler<BasePopulationTree> >(2.0);
        op_schedule.add_operator(op);
        op = std::make_shared< NodeHeightSlideBumpScaler<BasePopulationTree> >(2.0);
        op_schedule.add_operator(op);
        op = std::make_shared< NodeHeightSlideBumpSwapScaler<BasePopulationTree> >(2.0);
        op_schedule.add_operator(op);
        op = std::make_shared< NeighborHeightNodeSwap<BasePopulationTree> >(1.0);
        op_schedule.add_operator(op);
        op = std::make_shared< NodeHeightSlideBumpPermuteScaler<BasePopulationTree> >(2.0);
        op_schedule.add_operator(op);
        op = std::make_shared< NeighborHeightNodePermute<BasePopulationTree> >(1.0);
        op_schedule.add_operator(op);
        op = std::make_shared< NodeHeightDirichletOperator<BasePopulationTree> >(2.0);
        op_schedule.add_operator(op);
        op = std::make_shared< GlobalPopSizeScaler >(1.0);
        op_schedule.add_operator(op);
        op = std::make_shared< StateFreqMover >(1.0);
        op_schedule.add_operator(op);
        op = std::make_shared< StateFreqDirichletOperator >(1.0);
        op_schedule.add_operator(op);
        op = std::make_shared< NodeHeightPriorAlphaScaler<BasePopulationTree> >(3.0);
        op_schedule.add_operator(op);
        op = std::make_shared< NodeHeightPriorBetaScaler<BasePopulationTree> >(3.0);
        op_schedule.add_operator(op);
        op = std::make_shared< GlobalHeightSizeRateScaler >(1.0);
        op_schedule.add_operator(op);
        op = std::make_shared< GlobalNodeHeightDirichletOperator<BasePopulationTree> >(1.0);
        op_schedule.add_operator(op);
        op = std::make_shared< GlobalHeightSizeMixer >(1.0);
        op_schedule.add_operator(op);

        op = std::make_shared< PopSizeScaler >(2.0);
        op_schedule.add_operator(op);

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate< BasePopulationTree > > > time_ops = op_schedule.get_node_height_operators();
        time_ops.push_back(op);

        for (unsigned int i = 0; i < op_schedule.get_number_of_operators(); ++i) {
            op = op_schedule.get_operator(i);
            if (op->get_name() != "SubtreePruneRegraftRevJumpSampler") {
                op->turn_on_auto_optimize();
                op->set_auto_optimize_delay(1000);
            }
        }


        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int count_nheights_1 = 0;
        unsigned int count_nheights_2 = 0;
        unsigned int count_nheights_3 = 0;
        unsigned int count_nheights_4 = 0;

        SampleSummarizer<double> pop_size_summary;
        SampleSummarizer<double> root_pop_size_summary;
        SampleSummarizer<double> leaf0_pop_size_summary;
        SampleSummarizer<double> leaf1_pop_size_summary;
        SampleSummarizer<double> leaf2_pop_size_summary;
        SampleSummarizer<double> leaf3_pop_size_summary;
        SampleSummarizer<double> leaf4_pop_size_summary;
        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> height_alpha_summary;
        SampleSummarizer<double> height_beta_summary;

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes;

        unsigned int burnin = 10000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op = op_schedule.draw_operator(rng);
            if (op->get_type() == BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::topology_model_operator) {
                op->operate_plus(rng,
                        &tree,
                        time_ops,
                        1, 1, 1);
            }
            else {
                op->operate(rng, &tree, 1, 1);
            }
        }

        unsigned int niterations = 50000000;
        unsigned int sample_freq = 100;
        unsigned int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 10000;
        for (unsigned int i = 0; i < niterations; ++i) {
            op = op_schedule.draw_operator(rng);
            if (op->get_type() == BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::topology_model_operator) {
                op->operate_plus(rng,
                        &tree,
                        time_ops,
                        1, 1, 1);
            }
            else {
                op->operate(rng, &tree, 1, 1);
            }
            double a = height_alpha_prior->draw(rng);
            double b = height_beta_prior->draw(rng);
            double v = BetaDistribution::get_draw(rng, a, b);
            if ((i + 1) % sample_freq == 0) {
                pop_sizes = tree.get_pointers_to_population_sizes();
                for (auto pop_size : pop_sizes) {
                    pop_size_summary.add_sample(pop_size->get_value());
                }
                root_pop_size_summary.add_sample(tree.get_root_population_size());
                leaf0_pop_size_summary.add_sample(tree.get_node("leaf0")->get_population_size());
                leaf1_pop_size_summary.add_sample(tree.get_node("leaf1")->get_population_size());
                leaf2_pop_size_summary.add_sample(tree.get_node("leaf2")->get_population_size());
                leaf3_pop_size_summary.add_sample(tree.get_node("leaf3")->get_population_size());
                leaf4_pop_size_summary.add_sample(tree.get_node("leaf4")->get_population_size());
                root_height_summary.add_sample(tree.get_root_height());
                height_alpha_summary.add_sample(tree.get_alpha_of_node_height_beta_prior());
                height_beta_summary.add_sample(tree.get_beta_of_node_height_beta_prior());
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
                ++sample_count;
                if (sample_count % report_freq == 0) {
                    std::cout << "Sampled " << sample_count << " of " << nsamples << std::endl;
                }
            }
        }
        op_schedule.write_operator_rates(std::cout);

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
            std::cout << "Tree:\n";
            for (auto splitset : s_c.first) {
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
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            std::cout << "  nsamples: " << s_c.second << "\n";
            std::cout << "  prop error: " << prop_error << "\n";
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }

        double quantile_chi_sq_335_10 = 368.6;
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        std::cout << "Chi-square(335) 0.9 quantile: " << quantile_chi_sq_335_10 << "\n";

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

        write_r_script(split_counts, 5, "../5-leaf-general-tree-test-full-model-free-pop-sizes.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        REQUIRE(split_counts.size() == 336);

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
        REQUIRE(height_alpha_summary.mean() == Approx(height_alpha_prior->get_mean()).epsilon(eps * 2.0));
        REQUIRE(height_alpha_summary.variance() == Approx(height_alpha_prior->get_variance()).epsilon(eps * 2.0));
        REQUIRE(height_beta_summary.mean() == Approx(height_beta_prior->get_mean()).epsilon(eps * 2.0));
        REQUIRE(height_beta_summary.variance() == Approx(height_beta_prior->get_variance()).epsilon(eps * 2.0));

        for (auto s_f : split_freqs) {
            REQUIRE(s_f.second == Approx(exp_freq).epsilon(eps));
        }

        // REQUIRE(chi_sq_test_statistic < quantile_chi_sq_335_10);
    }
}
