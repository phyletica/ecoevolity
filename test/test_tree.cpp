#include "catch.hpp"
#include "ecoevolity/tree.hpp"
#include "ecoevolity/stats_util.hpp"


TEST_CASE("Testing BaseTree", "[BaseTree]") {
    SECTION("Testing three species") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.1);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.07);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.03);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        leaf0->fix_node_height();
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        leaf1->fix_node_height();
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        leaf2->fix_node_height();
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        leaf3->fix_node_height();

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(leaf3);
        root->add_child(internal0);
        root->add_child(internal1);
        BaseTree<Node> tree(root);

        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_leaf_node_count() == 4);
        REQUIRE(tree.get_node_count() == 7);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        std::vector<double> expected_heights {0.03, 0.07, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);
    }
}

TEST_CASE("Testing BaseTree::get_nearest_height_index", "[BaseTree]") {
    SECTION("Testing get_nearest_height_index") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.1);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.04);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.06);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.08);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        leaf0->fix_node_height();
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        leaf1->fix_node_height();
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        leaf2->fix_node_height();
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        leaf3->fix_node_height();
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
        leaf4->fix_node_height();

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(leaf3);
        internal2->add_child(internal0);
        internal2->add_child(internal1);
        root->add_child(internal2);
        root->add_child(leaf4);
        BaseTree<Node> tree(root);

        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_leaf_node_count() == 5);
        REQUIRE(tree.get_node_count() == 9);
        REQUIRE(tree.get_number_of_node_heights() == 4);
        std::vector<double> expected_heights {0.04, 0.06, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        REQUIRE(tree.get_nearest_height_index(0.0) == 0);
        REQUIRE(tree.get_nearest_height_index(2.0) == 3);

        // Tie break goes to larger node height
        REQUIRE(tree.get_nearest_height_index(0.0499) == 0);
        REQUIRE(tree.get_nearest_height_index(0.05) == 1);
        REQUIRE(tree.get_nearest_height_index(0.0699) == 1);
        REQUIRE(tree.get_nearest_height_index(0.07) == 2);
        REQUIRE(tree.get_nearest_height_index(0.0899) == 2);
        REQUIRE(tree.get_nearest_height_index(0.09000001) == 3);
    }
}

TEST_CASE("Testing BaseTree::get_intervening_height_indices", "[BaseTree]") {
    SECTION("Testing get_intervening_height_indices") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.1);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.04);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.06);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.08);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        leaf0->fix_node_height();
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        leaf1->fix_node_height();
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        leaf2->fix_node_height();
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        leaf3->fix_node_height();
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
        leaf4->fix_node_height();

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(leaf3);
        internal2->add_child(internal0);
        internal2->add_child(internal1);
        root->add_child(internal2);
        root->add_child(leaf4);
        BaseTree<Node> tree(root);

        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_leaf_node_count() == 5);
        REQUIRE(tree.get_node_count() == 9);
        REQUIRE(tree.get_number_of_node_heights() == 4);
        std::vector<double> expected_heights {0.04, 0.06, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        std::vector<unsigned int> expected_indices = {};
        REQUIRE(tree.get_intervening_height_indices(0, 0.0) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(0, 0.0599) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(1, 0.041) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(1, 0.0799) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(2, 0.061) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(2, 0.0999) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(3, 0.081) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(3, 1.1) == expected_indices);

        expected_indices = {0};
        // if value is equal, height should be included
        REQUIRE(tree.get_intervening_height_indices(1, 0.04) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(1, 0.039) == expected_indices);

        expected_indices = {1};
        REQUIRE(tree.get_intervening_height_indices(0, 0.06) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(0, 0.079) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(2, 0.06) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(2, 0.041) == expected_indices);

        expected_indices = {2};
        REQUIRE(tree.get_intervening_height_indices(1, 0.08) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(1, 0.099) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(3, 0.08) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(3, 0.061) == expected_indices);

        expected_indices = {3};
        REQUIRE(tree.get_intervening_height_indices(2, 1.0) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(2, 1.01) == expected_indices);

        expected_indices = {1, 0};
        REQUIRE(tree.get_intervening_height_indices(2, 0.04) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(2, 0.0) == expected_indices);

        expected_indices = {1, 2};
        REQUIRE(tree.get_intervening_height_indices(0, 0.08) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(0, 0.099) == expected_indices);
        expected_indices = {2, 1};
        REQUIRE(tree.get_intervening_height_indices(3, 0.06) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(3, 0.0401) == expected_indices);

        expected_indices = {2, 3};
        REQUIRE(tree.get_intervening_height_indices(1, 1.0) == expected_indices);
        REQUIRE(tree.get_intervening_height_indices(1, 2.0) == expected_indices);
    }
}

TEST_CASE("Testing BaseTree::height_is_splittable", "[BaseTree]") {
    SECTION("Testing height_is_splittable") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.1);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.04);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.06);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.08);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        leaf0->fix_node_height();
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        leaf1->fix_node_height();
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        leaf2->fix_node_height();
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        leaf3->fix_node_height();
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
        leaf4->fix_node_height();

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(leaf3);
        internal2->add_child(internal0);
        internal2->add_child(internal1);
        root->add_child(internal2);
        root->add_child(leaf4);
        BaseTree<Node> tree(root);

        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_leaf_node_count() == 5);
        REQUIRE(tree.get_node_count() == 9);
        REQUIRE(tree.get_number_of_node_heights() == 4);
        std::vector<double> expected_heights = {0.04, 0.06, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        for (unsigned int i = 0; i < tree.get_number_of_node_heights(); ++i) {
            REQUIRE(tree.height_is_splittable(i) == false);
        }

        internal0->set_height_parameter(internal1->get_height_parameter());
        tree.update_node_heights();

        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_leaf_node_count() == 5);
        REQUIRE(tree.get_node_count() == 9);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        expected_heights = {0.06, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        REQUIRE(tree.height_is_splittable(0) == true);
        REQUIRE(tree.height_is_splittable(1) == false);
        REQUIRE(tree.height_is_splittable(2) == false);

        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
        leaf5->fix_node_height();
        root->add_child(leaf5);
        tree.update_node_heights();

        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 3);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 10);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        expected_heights = {0.06, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        REQUIRE(tree.height_is_splittable(0) == true);
        REQUIRE(tree.height_is_splittable(1) == false);
        REQUIRE(tree.height_is_splittable(2) == true);
    }
}

TEST_CASE("Testing BaseTree::get_indices_of_splittable_heights", "[BaseTree]") {
    SECTION("Testing get_indices_of_splittable_heights") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.1);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.04);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.06);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.08);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        leaf0->fix_node_height();
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        leaf1->fix_node_height();
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        leaf2->fix_node_height();
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        leaf3->fix_node_height();
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
        leaf4->fix_node_height();

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(leaf3);
        internal2->add_child(internal0);
        internal2->add_child(internal1);
        root->add_child(internal2);
        root->add_child(leaf4);
        BaseTree<Node> tree(root);

        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_leaf_node_count() == 5);
        REQUIRE(tree.get_node_count() == 9);
        REQUIRE(tree.get_number_of_node_heights() == 4);
        std::vector<double> expected_heights = {0.04, 0.06, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        std::vector<unsigned int> expected_indices = {};
        REQUIRE(tree.get_indices_of_splittable_heights() == expected_indices);

        internal0->set_height_parameter(internal1->get_height_parameter());
        tree.update_node_heights();

        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_leaf_node_count() == 5);
        REQUIRE(tree.get_node_count() == 9);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        expected_heights = {0.06, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        expected_indices = {0};
        REQUIRE(tree.get_indices_of_splittable_heights() == expected_indices);

        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
        leaf5->fix_node_height();
        root->add_child(leaf5);
        tree.update_node_heights();

        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 3);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 10);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        expected_heights = {0.06, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        expected_indices = {0, 2};
        REQUIRE(tree.get_indices_of_splittable_heights() == expected_indices);
    }
}

TEST_CASE("Testing BaseTree::slide_bump_height", "[BaseTree]") {
    SECTION("Testing slide_bump_height") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.1);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.04);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.06);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.08);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        leaf0->fix_node_height();
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        leaf1->fix_node_height();
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        leaf2->fix_node_height();
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        leaf3->fix_node_height();
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
        leaf4->fix_node_height();
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
        leaf5->fix_node_height();

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(leaf3);
        internal2->add_child(internal0);
        internal2->add_child(internal1);
        root->add_child(internal2);
        root->add_child(leaf4);
        root->add_child(leaf5);

        internal0->set_height_parameter(internal1->get_height_parameter());

        BaseTree<Node> tree(root);

        REQUIRE(root->is_root());
        REQUIRE(root->is_child(internal2));
        REQUIRE(root->is_child(leaf4));
        REQUIRE(root->is_child(leaf5));
        REQUIRE(internal2->is_parent(root));
        REQUIRE(internal2->is_child(internal1));
        REQUIRE(internal2->is_child(internal0));
        REQUIRE(internal1->is_parent(internal2));
        REQUIRE(internal1->is_child(leaf2));
        REQUIRE(internal1->is_child(leaf3));
        REQUIRE(internal0->is_parent(internal2));
        REQUIRE(internal0->is_child(leaf0));
        REQUIRE(internal0->is_child(leaf1));
        REQUIRE(leaf0->is_leaf());
        REQUIRE(leaf0->is_parent(internal0));
        REQUIRE(leaf1->is_leaf());
        REQUIRE(leaf1->is_parent(internal0));
        REQUIRE(leaf2->is_leaf());
        REQUIRE(leaf2->is_parent(internal1));
        REQUIRE(leaf3->is_leaf());
        REQUIRE(leaf3->is_parent(internal1));
        REQUIRE(leaf4->is_leaf());
        REQUIRE(leaf4->is_parent(root));
        REQUIRE(leaf5->is_leaf());
        REQUIRE(leaf5->is_parent(root));

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 3);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 10);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        std::vector<double> expected_heights = {0.06, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        RandomNumberGenerator rng(111);

        tree.slide_bump_height(rng, 0, 0.07);

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 3);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 10);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        expected_heights = {0.07, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        tree.slide_bump_height(rng, 0, 0.01);

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 3);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 10);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        expected_heights = {0.01, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        tree.slide_bump_height(rng, 0, 0.2);

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_degree_of_root() == 3);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 10);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        expected_heights = {0.08, 0.1, 0.2};
        REQUIRE(tree.get_node_heights() == expected_heights);

        tree.slide_bump_height(rng, 2, 0.05);

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 3);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 10);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        expected_heights = {0.05, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        tree.slide_bump_height(rng, 2, 0.11);

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_root_height() == 0.11);
        REQUIRE(tree.get_degree_of_root() == 3);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 10);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        expected_heights = {0.05, 0.08, 0.11};
        REQUIRE(tree.get_node_heights() == expected_heights);

        tree.slide_bump_height(rng, 2, 0.1);

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 3);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 10);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        expected_heights = {0.05, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        tree.slide_bump_height(rng, 1, 0.07);

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 3);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 10);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        expected_heights = {0.05, 0.07, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        tree.slide_bump_height(rng, 1, 0.09);

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 3);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 10);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        expected_heights = {0.05, 0.09, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        tree.slide_bump_height(rng, 1, 0.11);

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_root_height() == 0.11);
        REQUIRE(tree.get_degree_of_root() == 3);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 10);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        expected_heights = {0.05, 0.1, 0.11};
        REQUIRE(tree.get_node_heights() == expected_heights);

        tree.slide_bump_height(rng, 1, 0.01);

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_root_height() == 0.11);
        REQUIRE(tree.get_degree_of_root() == 3);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 10);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        expected_heights = {0.01, 0.05, 0.11};
        REQUIRE(tree.get_node_heights() == expected_heights);

        REQUIRE(root->is_root());
        REQUIRE(root->is_child(internal2));
        REQUIRE(root->is_child(leaf4));
        REQUIRE(root->is_child(leaf5));
        REQUIRE(internal2->is_parent(root));
        REQUIRE(internal2->is_child(internal1));
        REQUIRE(internal2->is_child(internal0));
        REQUIRE(internal1->is_parent(internal2));
        REQUIRE(internal1->is_child(leaf2));
        REQUIRE(internal1->is_child(leaf3));
        REQUIRE(internal0->is_parent(internal2));
        REQUIRE(internal0->is_child(leaf0));
        REQUIRE(internal0->is_child(leaf1));
        REQUIRE(leaf0->is_leaf());
        REQUIRE(leaf0->is_parent(internal0));
        REQUIRE(leaf1->is_leaf());
        REQUIRE(leaf1->is_parent(internal0));
        REQUIRE(leaf2->is_leaf());
        REQUIRE(leaf2->is_parent(internal1));
        REQUIRE(leaf3->is_leaf());
        REQUIRE(leaf3->is_parent(internal1));
        REQUIRE(leaf4->is_leaf());
        REQUIRE(leaf4->is_parent(root));
        REQUIRE(leaf5->is_leaf());
        REQUIRE(leaf5->is_parent(root));
    }
}

TEST_CASE("Testing BaseTree::merge_node_height_up", "[BaseTree]") {
    SECTION("Testing merge_node_height_up") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.1);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.04);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.06);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.08);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        leaf0->fix_node_height();
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        leaf1->fix_node_height();
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        leaf2->fix_node_height();
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        leaf3->fix_node_height();
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
        leaf4->fix_node_height();
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
        leaf5->fix_node_height();

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(leaf3);
        internal2->add_child(internal0);
        internal2->add_child(internal1);
        root->add_child(internal2);
        root->add_child(leaf4);
        root->add_child(leaf5);

        internal0->set_height_parameter(internal1->get_height_parameter());

        BaseTree<Node> tree(root);

        REQUIRE(root->is_root());
        REQUIRE(root->is_child(internal2));
        REQUIRE(root->is_child(leaf4));
        REQUIRE(root->is_child(leaf5));
        REQUIRE(internal2->is_parent(root));
        REQUIRE(internal2->is_child(internal1));
        REQUIRE(internal2->is_child(internal0));
        REQUIRE(internal1->is_parent(internal2));
        REQUIRE(internal1->is_child(leaf2));
        REQUIRE(internal1->is_child(leaf3));
        REQUIRE(internal0->is_parent(internal2));
        REQUIRE(internal0->is_child(leaf0));
        REQUIRE(internal0->is_child(leaf1));
        REQUIRE(leaf0->is_leaf());
        REQUIRE(leaf0->is_parent(internal0));
        REQUIRE(leaf1->is_leaf());
        REQUIRE(leaf1->is_parent(internal0));
        REQUIRE(leaf2->is_leaf());
        REQUIRE(leaf2->is_parent(internal1));
        REQUIRE(leaf3->is_leaf());
        REQUIRE(leaf3->is_parent(internal1));
        REQUIRE(leaf4->is_leaf());
        REQUIRE(leaf4->is_parent(root));
        REQUIRE(leaf5->is_leaf());
        REQUIRE(leaf5->is_parent(root));

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 3);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 10);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        std::vector<double> expected_heights = {0.06, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        tree.merge_node_height_up(0);

        REQUIRE(root->is_root());
        REQUIRE(root->is_child(internal2));
        REQUIRE(root->is_child(leaf4));
        REQUIRE(root->is_child(leaf5));
        REQUIRE(internal2->is_parent(root));
        REQUIRE(internal2->is_child(leaf0));
        REQUIRE(internal2->is_child(leaf1));
        REQUIRE(internal2->is_child(leaf2));
        REQUIRE(internal2->is_child(leaf3));
        REQUIRE(leaf0->is_leaf());
        REQUIRE(leaf0->is_parent(internal2));
        REQUIRE(leaf1->is_leaf());
        REQUIRE(leaf1->is_parent(internal2));
        REQUIRE(leaf2->is_leaf());
        REQUIRE(leaf2->is_parent(internal2));
        REQUIRE(leaf3->is_leaf());
        REQUIRE(leaf3->is_parent(internal2));
        REQUIRE(leaf4->is_leaf());
        REQUIRE(leaf4->is_parent(root));
        REQUIRE(leaf5->is_leaf());
        REQUIRE(leaf5->is_parent(root));

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 3);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 8);
        REQUIRE(tree.get_number_of_node_heights() == 2);
        expected_heights = {0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        tree.merge_node_height_up(0);

        REQUIRE(root->is_root());
        REQUIRE(root->is_child(leaf0));
        REQUIRE(root->is_child(leaf1));
        REQUIRE(root->is_child(leaf2));
        REQUIRE(root->is_child(leaf3));
        REQUIRE(root->is_child(leaf4));
        REQUIRE(root->is_child(leaf5));
        REQUIRE(leaf0->is_leaf());
        REQUIRE(leaf0->is_parent(root));
        REQUIRE(leaf1->is_leaf());
        REQUIRE(leaf1->is_parent(root));
        REQUIRE(leaf2->is_leaf());
        REQUIRE(leaf2->is_parent(root));
        REQUIRE(leaf3->is_leaf());
        REQUIRE(leaf3->is_parent(root));
        REQUIRE(leaf4->is_leaf());
        REQUIRE(leaf4->is_parent(root));
        REQUIRE(leaf5->is_leaf());
        REQUIRE(leaf5->is_parent(root));

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 6);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 7);
        REQUIRE(tree.get_number_of_node_heights() == 1);
        expected_heights = {0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);
    }
}

TEST_CASE("Testing BaseTree::split_node_height_down", "[BaseTree]") {
    SECTION("Testing split_node_height_down") {
        double height_lower_bound;
        unsigned int number_of_mapped_nodes;
        std::vector<unsigned int> mapped_polytomy_sizes;

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.1);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        leaf0->fix_node_height();
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        leaf1->fix_node_height();
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        leaf2->fix_node_height();
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        leaf3->fix_node_height();
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
        leaf4->fix_node_height();
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
        leaf5->fix_node_height();

        root->add_child(leaf0);
        root->add_child(leaf1);
        root->add_child(leaf2);
        root->add_child(leaf3);
        root->add_child(leaf4);
        root->add_child(leaf5);

        std::unordered_set< std::shared_ptr<Node> > existing_nodes;
        existing_nodes.insert(root);
        existing_nodes.insert(leaf0);
        existing_nodes.insert(leaf1);
        existing_nodes.insert(leaf2);
        existing_nodes.insert(leaf3);
        existing_nodes.insert(leaf4);
        existing_nodes.insert(leaf5);

        BaseTree<Node> tree(root);

        REQUIRE(root->is_root());
        REQUIRE(root->is_child(leaf0));
        REQUIRE(root->is_child(leaf1));
        REQUIRE(root->is_child(leaf2));
        REQUIRE(root->is_child(leaf3));
        REQUIRE(root->is_child(leaf4));
        REQUIRE(root->is_child(leaf5));
        REQUIRE(leaf0->is_leaf());
        REQUIRE(leaf0->is_parent(root));
        REQUIRE(leaf1->is_leaf());
        REQUIRE(leaf1->is_parent(root));
        REQUIRE(leaf2->is_leaf());
        REQUIRE(leaf2->is_parent(root));
        REQUIRE(leaf3->is_leaf());
        REQUIRE(leaf3->is_parent(root));
        REQUIRE(leaf4->is_leaf());
        REQUIRE(leaf4->is_parent(root));
        REQUIRE(leaf5->is_leaf());
        REQUIRE(leaf5->is_parent(root));

        REQUIRE(tree.tree_is_valid());
        double expeced_root_height = 0.1;
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 6);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 7);
        REQUIRE(tree.get_number_of_node_heights() == 1);
        std::vector<double> expected_heights = {0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        RandomNumberGenerator rng = RandomNumberGenerator(111);

        std::vector<unsigned int> splittable_heights;
        splittable_heights = tree.get_indices_of_splittable_heights();
        REQUIRE(splittable_heights.size() == 1);
        REQUIRE(splittable_heights.at(0) == 0);

        tree.split_node_height_down(rng, 0,
                height_lower_bound,
                number_of_mapped_nodes,
                mapped_polytomy_sizes);

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() < 6);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() > 7);
        unsigned int expected_number_of_node_heights = 2;
        REQUIRE(tree.get_number_of_node_heights() == expected_number_of_node_heights);
        std::shared_ptr<PositiveRealParameter> new_height;
        bool new_height_set = false;
        for (unsigned int i = 0; i < root->get_number_of_children(); ++i) {
            if (existing_nodes.count(root->get_child(i)) < 1) {
                existing_nodes.insert(root->get_child(i));
                if (! new_height_set) {
                    new_height = root->get_child(i)->get_height_parameter();
                    new_height_set = true;
                    REQUIRE(new_height->get_value() < root->get_height());
                    expected_heights.push_back(new_height->get_value());
                    std::sort(expected_heights.begin(), expected_heights.end());
                }
                else {
                    REQUIRE(root->get_child(i)->get_height_parameter() == new_height);
                }
            }
        }
        REQUIRE(tree.get_node_heights() == expected_heights);

        std::vector< std::shared_ptr<Node> > mapped_nodes_new;
        std::vector< std::shared_ptr<Node> > mapped_nodes_orig;
        std::shared_ptr<PositiveRealParameter> orig_height;
        while (true) {
            splittable_heights = tree.get_indices_of_splittable_heights();
            if (splittable_heights.size() < 1) {
                break;
            }
            unsigned int splittable_ht = splittable_heights.at(0);
            orig_height = tree.get_height_parameter(splittable_ht);
            tree.split_node_height_down(rng, splittable_ht,
                    height_lower_bound,
                    number_of_mapped_nodes,
                    mapped_polytomy_sizes);
            ++expected_number_of_node_heights;
            mapped_nodes_orig = tree.get_mapped_nodes(splittable_ht + 1);
            mapped_nodes_new = tree.get_mapped_nodes(splittable_ht);

            REQUIRE(tree.tree_is_valid());
            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_leaf_node_count() == 6);
            REQUIRE(tree.get_number_of_node_heights() == expected_number_of_node_heights);
            bool new_height_set = false;
            for (auto mapped_nd : mapped_nodes_new) {
                if (existing_nodes.count(mapped_nd) < 1) {
                    existing_nodes.insert(mapped_nd);
                }
                if (! new_height_set) {
                    new_height = mapped_nd->get_height_parameter();
                    new_height_set = true;
                    REQUIRE(new_height->get_value() < orig_height->get_value());
                    expected_heights.push_back(new_height->get_value());
                    std::sort(expected_heights.begin(), expected_heights.end());
                }
                else {
                    REQUIRE(mapped_nd->get_height_parameter() == new_height);
                }
            }
            for (auto mapped_nd : mapped_nodes_orig) {
                REQUIRE(existing_nodes.count(mapped_nd) > 0);
                REQUIRE(mapped_nd->get_height_parameter() == orig_height);
            }
            REQUIRE(tree.get_node_heights() == expected_heights);
        }

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 11);
        REQUIRE(tree.get_number_of_node_heights() == 5);
        REQUIRE(tree.get_node_heights() == expected_heights);
    }
}

TEST_CASE("Testing BaseTree::get_collision_parents", "[BaseTree]") {
    SECTION("Testing get_collision_parents") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.1);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.04);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.06);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.08);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        leaf0->fix_node_height();
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        leaf1->fix_node_height();
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        leaf2->fix_node_height();
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        leaf3->fix_node_height();
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
        leaf4->fix_node_height();
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
        leaf5->fix_node_height();

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(leaf3);
        internal2->add_child(internal0);
        internal2->add_child(internal1);
        root->add_child(internal2);
        root->add_child(leaf4);
        root->add_child(leaf5);

        internal0->set_height_parameter(internal1->get_height_parameter());

        BaseTree<Node> tree(root);

        REQUIRE(root->is_root());
        REQUIRE(root->is_child(internal2));
        REQUIRE(root->is_child(leaf4));
        REQUIRE(root->is_child(leaf5));
        REQUIRE(internal2->is_parent(root));
        REQUIRE(internal2->is_child(internal1));
        REQUIRE(internal2->is_child(internal0));
        REQUIRE(internal1->is_parent(internal2));
        REQUIRE(internal1->is_child(leaf2));
        REQUIRE(internal1->is_child(leaf3));
        REQUIRE(internal0->is_parent(internal2));
        REQUIRE(internal0->is_child(leaf0));
        REQUIRE(internal0->is_child(leaf1));
        REQUIRE(leaf0->is_leaf());
        REQUIRE(leaf0->is_parent(internal0));
        REQUIRE(leaf1->is_leaf());
        REQUIRE(leaf1->is_parent(internal0));
        REQUIRE(leaf2->is_leaf());
        REQUIRE(leaf2->is_parent(internal1));
        REQUIRE(leaf3->is_leaf());
        REQUIRE(leaf3->is_parent(internal1));
        REQUIRE(leaf4->is_leaf());
        REQUIRE(leaf4->is_parent(root));
        REQUIRE(leaf5->is_leaf());
        REQUIRE(leaf5->is_parent(root));

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_degree_of_root() == 3);
        REQUIRE(tree.get_leaf_node_count() == 6);
        REQUIRE(tree.get_node_count() == 10);
        REQUIRE(tree.get_number_of_node_heights() == 3);
        std::vector<double> expected_heights = {0.06, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        std::vector< std::shared_ptr<Node> > expected_parents {root};
        std::vector< std::shared_ptr<Node> > parents = tree.get_collision_parents(2, 1);
        REQUIRE(parents == expected_parents);

        expected_parents = {internal2};
        parents = tree.get_collision_parents(1, 0);
        REQUIRE(parents == expected_parents);
    }
}

TEST_CASE("Testing BaseTree::collision_node_permute with 3 leaves", "[BaseTree]") {
    SECTION("Testing collision_node_permute") {
        unsigned int nsamples = 100000;
        RandomNumberGenerator rng = RandomNumberGenerator(111);
        unsigned int count_0_12 = 0;
        unsigned int count_1_02 = 0;
        unsigned int count_2_10 = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            leaf0->fix_node_height();
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            leaf1->fix_node_height();
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            leaf2->fix_node_height();

            internal0->add_child(leaf0);
            internal0->add_child(leaf1);
            root->add_child(internal0);
            root->add_child(leaf2);

            BaseTree<Node> tree(root);

            tree.collision_node_permute(rng, 1, 0);

            REQUIRE(tree.tree_is_valid());
            REQUIRE(root->get_number_of_children() == 2);
            REQUIRE(internal0->get_number_of_children() == 2);
            REQUIRE(leaf0->is_leaf());
            REQUIRE(leaf1->is_leaf());
            REQUIRE(leaf2->is_leaf());
            std::vector< std::shared_ptr<Node> > expected_leaves {leaf0, leaf1, leaf2};
            std::vector< std::shared_ptr<Node> > leaves = root->get_leaves();
            REQUIRE(std::is_permutation(
                        leaves.begin(), leaves.end(),
                        expected_leaves.begin()));
            if (root->is_child(leaf0)) {
                ++count_0_12;
            }
            if (root->is_child(leaf1)) {
                ++count_1_02;
            }
            if (root->is_child(leaf2)) {
                ++count_2_10;
            }
        }
        REQUIRE(count_0_12 + count_1_02 + count_2_10 == nsamples);
        double freq_0_12 = (count_0_12 / (double)nsamples);
        double freq_1_02 = (count_1_02 / (double)nsamples);
        double freq_2_10 = (count_2_10 / (double)nsamples);
        std::cout << "freq of (0,(1,2)): " << freq_0_12 << "\n";
        std::cout << "freq of (1,(0,2)): " << freq_1_02 << "\n";
        std::cout << "freq of (2,(0,1)): " << freq_2_10 << "\n";

        double eps = 0.005;
        // p(0 in pool) * p(0 drawn) = 0.5*0.5
        REQUIRE(freq_0_12 == Approx(0.25).epsilon(eps));
        // p(1 in pool) * p(1 drawn) = 0.5*0.5
        REQUIRE(freq_1_02 == Approx(0.25).epsilon(eps));
        // p(2 in pool) * p(1 drawn) = 1.0*0.5
        REQUIRE(freq_2_10 == Approx(0.5).epsilon(eps));
    }
}

TEST_CASE("Testing BaseTree::collision_node_permute 4 leaf polytomy", "[BaseTree]") {
    SECTION("Testing collision_node_permute") {
        unsigned int nsamples = 100000;
        RandomNumberGenerator rng = RandomNumberGenerator(111);
        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_03 = 0;
        unsigned int count_12 = 0;
        unsigned int count_13 = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            leaf0->fix_node_height();
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            leaf1->fix_node_height();
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            leaf2->fix_node_height();
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
            leaf3->fix_node_height();

            internal0->add_child(leaf0);
            internal0->add_child(leaf1);
            root->add_child(internal0);
            root->add_child(leaf2);
            root->add_child(leaf3);

            BaseTree<Node> tree(root);

            tree.collision_node_permute(rng, 1, 0);

            REQUIRE(tree.tree_is_valid());
            REQUIRE(root->get_number_of_children() == 3);
            REQUIRE(internal0->get_number_of_children() == 2);
            std::vector< std::shared_ptr<Node> > expected_leaves {leaf0, leaf1, leaf2, leaf3};
            std::vector< std::shared_ptr<Node> > leaves = root->get_leaves();
            REQUIRE(std::is_permutation(
                        leaves.begin(), leaves.end(),
                        expected_leaves.begin()));
            if (internal0->is_child(leaf0) && internal0->is_child(leaf1)) {
                ++count_01;
            }
            if (internal0->is_child(leaf0) && internal0->is_child(leaf2)) {
                ++count_02;
            }
            if (internal0->is_child(leaf0) && internal0->is_child(leaf3)) {
                ++count_03;
            }
            if (internal0->is_child(leaf1) && internal0->is_child(leaf2)) {
                ++count_12;
            }
            if (internal0->is_child(leaf1) && internal0->is_child(leaf3)) {
                ++count_13;
            }
        }
        REQUIRE(count_01 + count_02 + count_03 + count_12 + count_13 == nsamples);
        double freq_01 = (count_01 / (double)nsamples);
        double freq_02 = (count_02 / (double)nsamples);
        double freq_03 = (count_03 / (double)nsamples);
        double freq_12 = (count_12 / (double)nsamples);
        double freq_13 = (count_13 / (double)nsamples);
        std::cout << "freq of (2,3(0,1)): " << freq_01 << "\n";
        std::cout << "freq of (1,3(0,2)): " << freq_02 << "\n";
        std::cout << "freq of (2,1(0,3)): " << freq_03 << "\n";
        std::cout << "freq of (0,3(2,1)): " << freq_12 << "\n";
        std::cout << "freq of (2,0(3,1)): " << freq_13 << "\n";

        double eps = 0.005;
        // (p(0 and 2 in pool) * p(0 drawn)) = 0.5^2 * 0.5 = 0.125
        // (p(0 and 3 in pool) * p(0 drawn)) = 0.5^2 * 0.5 = 0.125
        // (p(1 and 2 in pool) * p(0 drawn)) = 0.5^2 * 0.5 = 0.125
        // (p(1 and 3 in pool) * p(0 drawn)) = 0.5^2 * 0.5 = 0.125
        // total = 0.5
        REQUIRE(freq_01 == Approx(0.5).epsilon(eps));
        // p(0 and 2 in pool) * p(2 drawn) = 0.5^2*0.5
        REQUIRE(freq_02 == Approx(0.125).epsilon(eps));
        // p(0 and 3 in pool) * p(3 drawn) = 0.5^2*0.5
        REQUIRE(freq_03 == Approx(0.125).epsilon(eps));
        // p(1 and 2 in pool) * p(2 drawn) = 0.5^2*0.5
        REQUIRE(freq_12 == Approx(0.125).epsilon(eps));
        // p(1 and 3 in pool) * p(3 drawn) = 0.5^5*0.5
        REQUIRE(freq_13 == Approx(0.125).epsilon(eps));
    }
}

TEST_CASE("Testing BaseTree::collision_node_permute 4 leaf 3 colliders", "[BaseTree]") {
    SECTION("Testing collision_node_permute") {
        unsigned int nsamples = 100000;
        RandomNumberGenerator rng = RandomNumberGenerator(111);
        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_03 = 0;
        unsigned int count_12 = 0;
        unsigned int count_13 = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
            std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.5);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            leaf0->fix_node_height();
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            leaf1->fix_node_height();
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            leaf2->fix_node_height();
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
            leaf3->fix_node_height();

            internal0->add_child(leaf0);
            internal0->add_child(leaf1);
            internal1->add_child(leaf2);
            internal1->add_child(leaf3);
            root->add_child(internal0);
            root->add_child(internal1);

            internal0->set_height_parameter(internal1->get_height_parameter());

            BaseTree<Node> tree(root);

            tree.collision_node_permute(rng, 1, 0);

            REQUIRE(tree.tree_is_valid());
            REQUIRE(root->get_number_of_children() == 2);
            REQUIRE(internal0->get_number_of_children() == 2);
            REQUIRE(internal1->get_number_of_children() == 2);
            std::vector< std::shared_ptr<Node> > expected_leaves {leaf0, leaf1, leaf2, leaf3};
            std::vector< std::shared_ptr<Node> > leaves = root->get_leaves();
            REQUIRE(std::is_permutation(
                        leaves.begin(), leaves.end(),
                        expected_leaves.begin()));
            if (internal0->is_child(leaf0) && internal0->is_child(leaf1)) {
                ++count_01;
            }
            if (internal0->is_child(leaf0) && internal0->is_child(leaf2)) {
                ++count_02;
            }
            if (internal0->is_child(leaf0) && internal0->is_child(leaf3)) {
                ++count_03;
            }
            if (internal0->is_child(leaf1) && internal0->is_child(leaf2)) {
                ++count_12;
            }
            if (internal0->is_child(leaf1) && internal0->is_child(leaf3)) {
                ++count_13;
            }
        }
        REQUIRE(count_01 + count_02 + count_03 + count_12 + count_13 == nsamples);
        double freq_01 = (count_01 / (double)nsamples);
        double freq_02 = (count_02 / (double)nsamples);
        double freq_03 = (count_03 / (double)nsamples);
        double freq_12 = (count_12 / (double)nsamples);
        double freq_13 = (count_13 / (double)nsamples);
        std::cout << "freq of ((2,3),(0,1)): " << freq_01 << "\n";
        std::cout << "freq of ((1,3),(0,2)): " << freq_02 << "\n";
        std::cout << "freq of ((2,1),(0,3)): " << freq_03 << "\n";
        std::cout << "freq of ((0,3),(2,1)): " << freq_12 << "\n";
        std::cout << "freq of ((2,0),(3,1)): " << freq_13 << "\n";

        double eps = 0.005;
        // (p(0 and 2 in pool) * p(0 drawn)) = 0.5^2 * 0.5 = 0.125
        // (p(0 and 3 in pool) * p(0 drawn)) = 0.5^2 * 0.5 = 0.125
        // (p(1 and 2 in pool) * p(1 drawn)) = 0.5^2 * 0.5 = 0.125
        // (p(1 and 3 in pool) * p(1 drawn)) = 0.5^2 * 0.5 = 0.125
        // total = 0.5
        REQUIRE(freq_01 == Approx(0.5).epsilon(eps));
        // p(0 and 2 in pool) * p(2 drawn) = 0.5^2*0.5
        REQUIRE(freq_02 == Approx(0.125).epsilon(eps));
        // p(0 and 3 in pool) * p(3 drawn) = 0.5^2*0.5
        REQUIRE(freq_03 == Approx(0.125).epsilon(eps));
        // p(1 and 2 in pool) * p(2 drawn) = 0.5^2*0.5
        REQUIRE(freq_12 == Approx(0.125).epsilon(eps));
        // p(1 and 3 in pool) * p(3 drawn) = 0.5^5*0.5
        REQUIRE(freq_13 == Approx(0.125).epsilon(eps));
    }
}

TEST_CASE("Testing BaseTree::collision_node_permute 6 leaf 4 colliders", "[BaseTree]") {
    SECTION("Testing collision_node_permute") {
        unsigned int nsamples = 100000;
        RandomNumberGenerator rng = RandomNumberGenerator(111);
        unsigned int count_01_23_45 = 0;
        unsigned int count_05_21_43 = 0;
        unsigned int count_05_23_41 = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
            std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.5);
            std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.5);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            leaf0->fix_node_height();
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            leaf1->fix_node_height();
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            leaf2->fix_node_height();
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
            leaf3->fix_node_height();
            std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
            leaf4->fix_node_height();
            std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
            leaf5->fix_node_height();

            internal0->add_child(leaf0);
            internal0->add_child(leaf1);
            internal1->add_child(leaf2);
            internal1->add_child(leaf3);
            internal2->add_child(leaf4);
            internal2->add_child(leaf5);
            root->add_child(internal0);
            root->add_child(internal1);
            root->add_child(internal2);

            internal0->set_height_parameter(internal1->get_height_parameter());
            internal2->set_height_parameter(internal1->get_height_parameter());

            BaseTree<Node> tree(root);

            tree.collision_node_permute(rng, 1, 0);

            REQUIRE(tree.tree_is_valid());
            REQUIRE(root->get_number_of_children() == 3);
            REQUIRE(internal0->get_number_of_children() == 2);
            REQUIRE(internal1->get_number_of_children() == 2);
            REQUIRE(internal2->get_number_of_children() == 2);
            std::vector< std::shared_ptr<Node> > expected_leaves {leaf0, leaf1, leaf2, leaf3, leaf4, leaf5};
            std::vector< std::shared_ptr<Node> > leaves = root->get_leaves();
            REQUIRE(std::is_permutation(
                        leaves.begin(), leaves.end(),
                        expected_leaves.begin()));
            if (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf1) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf3) &&
                    internal2->is_child(leaf4) &&
                    internal2->is_child(leaf5)) {
                ++count_01_23_45;
            }
            if (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf5) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf1) &&
                    internal2->is_child(leaf4) &&
                    internal2->is_child(leaf3)) {
                ++count_05_21_43;
            }
            if (
                    (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf5) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf3) &&
                    internal2->is_child(leaf4) &&
                    internal2->is_child(leaf1)) ||
                    (
                    internal2->is_child(leaf0) &&
                    internal2->is_child(leaf5) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf3) &&
                    internal0->is_child(leaf4) &&
                    internal0->is_child(leaf1))
                ) {
                ++count_05_23_41;
            }
        }
        double freq_01_23_45 = (count_01_23_45 / (double)nsamples);
        double freq_05_21_43 = (count_05_21_43 / (double)nsamples);
        double freq_05_23_41 = (count_05_23_41 / (double)nsamples);
        std::cout << "freq of ((0,1),(2,3),(4,5)): " << freq_01_23_45 << "\n";
        std::cout << "freq of ((0,5),(2,1),(4,3)): " << freq_05_21_43 << "\n";
        std::cout << "freq of ((0,5),(2,3),(4,1)): " << freq_05_23_41 << "\n";

        double eps = 0.005;
        ///////////////////////////////////////////////////////////////////////
        // How I like to think about the probability of any swap result.
        // First, once you have a pool of N nodes, the number of unique
        // ways to sample them without replacement is N!, and so each
        // way of sampling the nodes in the pool has probability 1 / N!
        //
        // What differs between different outcomes is the probability of
        // getting a pool of nodes that is consistent with the outcome.
        // For a node that ends up with a different child, X, from the pool
        // than it added to it, the probability of this is simpy
        // (1 / the number of possible children X's parent could have contributed to the pool)
        //
        // For a node that ends up unchanged (i.e., it contributes a child to
        // the pool, only to get it back), the probability of a pool consistent
        // with that outcome is 1.0.
        // This is because the parent node has to contribute a node to the
        // pool, and no matter which one it contributes, it is always possible
        // to sample it back out of the pool.
        // For example, a bump node with children A and B, will add A with prob
        // 1/2 and B with prob 1/2, but either way, it could sample A / B back
        // out of the pool again, so every possible pool is consistent with an
        // outcome of an unchanged bump node.
        ///////////////////////////////////////////////////////////////////////
        //
        // p(pool) * p(draws) = 1.0 * (1/3!) = 1.0 * 1/6
        // All possible pools can lead to getting the same state back, so
        // p(pool) = 1.0
        REQUIRE(freq_01_23_45 == Approx(1.0/6.0).epsilon(eps));
        // p(pool) * p(draws) = (1/2)^3 * (1/3!) = 1/8 * 1/6 = 1/48
        REQUIRE(freq_05_21_43 == Approx(1.0/48.0).epsilon(eps));
        // This swap is symmetric, so the 01 on clade can contribute either
        // leaf (prob = 1), and then the 45 has to contribute the corresponding
        // leaf for the swap (prob 1/2) [ or vice versa ]
        // p(pool) * p(draws) = (1)(1)(1/2) * (1/3!) = 
        // (1/2) * (1/3!) = 1/2 * 1/6 = 1/12
        REQUIRE(freq_05_23_41 == Approx(1.0/12.0).epsilon(eps));
    }
}

TEST_CASE("Testing BaseTree::collision_node_swap with 3 leaves", "[BaseTree]") {
    SECTION("Testing collision_node_swap") {
        unsigned int nsamples = 100000;
        RandomNumberGenerator rng = RandomNumberGenerator(111);
        unsigned int count_0_12 = 0;
        unsigned int count_1_02 = 0;
        unsigned int count_2_10 = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            leaf0->fix_node_height();
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            leaf1->fix_node_height();
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            leaf2->fix_node_height();

            internal0->add_child(leaf0);
            internal0->add_child(leaf1);
            root->add_child(internal0);
            root->add_child(leaf2);

            BaseTree<Node> tree(root);

            tree.collision_node_swap(rng, 1, 0);

            REQUIRE(tree.tree_is_valid());
            REQUIRE(root->get_number_of_children() == 2);
            REQUIRE(internal0->get_number_of_children() == 2);
            REQUIRE(leaf0->is_leaf());
            REQUIRE(leaf1->is_leaf());
            REQUIRE(leaf2->is_leaf());
            std::vector< std::shared_ptr<Node> > expected_leaves {leaf0, leaf1, leaf2};
            std::vector< std::shared_ptr<Node> > leaves = root->get_leaves();
            REQUIRE(std::is_permutation(
                        leaves.begin(), leaves.end(),
                        expected_leaves.begin()));
            if (root->is_child(leaf0)) {
                ++count_0_12;
            }
            if (root->is_child(leaf1)) {
                ++count_1_02;
            }
            if (root->is_child(leaf2)) {
                ++count_2_10;
            }
        }
        REQUIRE(count_2_10 == 0);
        REQUIRE(count_0_12 + count_1_02 == nsamples);
        double freq_0_12 = (count_0_12 / (double)nsamples);
        double freq_1_02 = (count_1_02 / (double)nsamples);
        std::cout << "freq of (0,(1,2)): " << freq_0_12 << "\n";
        std::cout << "freq of (1,(0,2)): " << freq_1_02 << "\n";

        double eps = 0.005;
        REQUIRE(freq_0_12 == Approx(0.5).epsilon(eps));
        REQUIRE(freq_1_02 == Approx(0.5).epsilon(eps));
    }
}

TEST_CASE("Testing BaseTree::collision_node_swap 4 leaf polytomy", "[BaseTree]") {
    SECTION("Testing collision_node_swap") {
        unsigned int nsamples = 100000;
        RandomNumberGenerator rng = RandomNumberGenerator(111);
        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_03 = 0;
        unsigned int count_12 = 0;
        unsigned int count_13 = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            leaf0->fix_node_height();
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            leaf1->fix_node_height();
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            leaf2->fix_node_height();
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
            leaf3->fix_node_height();

            internal0->add_child(leaf0);
            internal0->add_child(leaf1);
            root->add_child(internal0);
            root->add_child(leaf2);
            root->add_child(leaf3);

            BaseTree<Node> tree(root);

            tree.collision_node_swap(rng, 1, 0);

            REQUIRE(tree.tree_is_valid());
            REQUIRE(root->get_number_of_children() == 3);
            REQUIRE(internal0->get_number_of_children() == 2);
            std::vector< std::shared_ptr<Node> > expected_leaves {leaf0, leaf1, leaf2, leaf3};
            std::vector< std::shared_ptr<Node> > leaves = root->get_leaves();
            REQUIRE(std::is_permutation(
                        leaves.begin(), leaves.end(),
                        expected_leaves.begin()));
            if (internal0->is_child(leaf0) && internal0->is_child(leaf1)) {
                ++count_01;
            }
            if (internal0->is_child(leaf0) && internal0->is_child(leaf2)) {
                ++count_02;
            }
            if (internal0->is_child(leaf0) && internal0->is_child(leaf3)) {
                ++count_03;
            }
            if (internal0->is_child(leaf1) && internal0->is_child(leaf2)) {
                ++count_12;
            }
            if (internal0->is_child(leaf1) && internal0->is_child(leaf3)) {
                ++count_13;
            }
        }
        REQUIRE(count_01 == 0);
        REQUIRE(count_02 + count_03 + count_12 + count_13 == nsamples);
        double freq_02 = (count_02 / (double)nsamples);
        double freq_03 = (count_03 / (double)nsamples);
        double freq_12 = (count_12 / (double)nsamples);
        double freq_13 = (count_13 / (double)nsamples);
        std::cout << "freq of (1,3(0,2)): " << freq_02 << "\n";
        std::cout << "freq of (2,1(0,3)): " << freq_03 << "\n";
        std::cout << "freq of (0,3(2,1)): " << freq_12 << "\n";
        std::cout << "freq of (2,0(3,1)): " << freq_13 << "\n";

        double eps = 0.005;
        REQUIRE(freq_02 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_03 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_12 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_13 == Approx(0.25).epsilon(eps));
    }
}

TEST_CASE("Testing BaseTree::collision_node_swap 4 leaf 3 colliders", "[BaseTree]") {
    SECTION("Testing collision_node_swap") {
        unsigned int nsamples = 100000;
        RandomNumberGenerator rng = RandomNumberGenerator(111);
        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_03 = 0;
        unsigned int count_12 = 0;
        unsigned int count_13 = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
            std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.5);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            leaf0->fix_node_height();
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            leaf1->fix_node_height();
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            leaf2->fix_node_height();
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
            leaf3->fix_node_height();

            internal0->add_child(leaf0);
            internal0->add_child(leaf1);
            internal1->add_child(leaf2);
            internal1->add_child(leaf3);
            root->add_child(internal0);
            root->add_child(internal1);

            internal0->set_height_parameter(internal1->get_height_parameter());

            BaseTree<Node> tree(root);

            tree.collision_node_swap(rng, 1, 0);

            REQUIRE(tree.tree_is_valid());
            REQUIRE(root->get_number_of_children() == 2);
            REQUIRE(internal0->get_number_of_children() == 2);
            REQUIRE(internal1->get_number_of_children() == 2);
            std::vector< std::shared_ptr<Node> > expected_leaves {leaf0, leaf1, leaf2, leaf3};
            std::vector< std::shared_ptr<Node> > leaves = root->get_leaves();
            REQUIRE(std::is_permutation(
                        leaves.begin(), leaves.end(),
                        expected_leaves.begin()));
            if (internal0->is_child(leaf0) && internal0->is_child(leaf1)) {
                ++count_01;
            }
            if (internal0->is_child(leaf0) && internal0->is_child(leaf2)) {
                ++count_02;
            }
            if (internal0->is_child(leaf0) && internal0->is_child(leaf3)) {
                ++count_03;
            }
            if (internal0->is_child(leaf1) && internal0->is_child(leaf2)) {
                ++count_12;
            }
            if (internal0->is_child(leaf1) && internal0->is_child(leaf3)) {
                ++count_13;
            }
        }
        REQUIRE(count_01 == 0);
        REQUIRE(count_02 + count_03 + count_12 + count_13 == nsamples);
        double freq_02 = (count_02 / (double)nsamples);
        double freq_03 = (count_03 / (double)nsamples);
        double freq_12 = (count_12 / (double)nsamples);
        double freq_13 = (count_13 / (double)nsamples);
        std::cout << "freq of ((1,3),(0,2)): " << freq_02 << "\n";
        std::cout << "freq of ((2,1),(0,3)): " << freq_03 << "\n";
        std::cout << "freq of ((0,3),(2,1)): " << freq_12 << "\n";
        std::cout << "freq of ((2,0),(3,1)): " << freq_13 << "\n";

        double eps = 0.005;
        REQUIRE(freq_02 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_03 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_12 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_13 == Approx(0.25).epsilon(eps));
    }
}

TEST_CASE("Testing BaseTree::collision_node_swap 6 leaf 4 colliders", "[BaseTree]") {
    SECTION("Testing collision_node_swap") {
        unsigned int nsamples = 100000;
        RandomNumberGenerator rng = RandomNumberGenerator(111);
        unsigned int count_01_23_45 = 0;
        unsigned int count_05_21_43 = 0;
        unsigned int count_05_23_41 = 0;
        unsigned int count_04_23_15 = 0;
        unsigned int count_01_25_43 = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
            std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.5);
            std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.5);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            leaf0->fix_node_height();
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            leaf1->fix_node_height();
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            leaf2->fix_node_height();
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
            leaf3->fix_node_height();
            std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
            leaf4->fix_node_height();
            std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
            leaf5->fix_node_height();

            internal0->add_child(leaf0);
            internal0->add_child(leaf1);
            internal1->add_child(leaf2);
            internal1->add_child(leaf3);
            internal2->add_child(leaf4);
            internal2->add_child(leaf5);
            root->add_child(internal0);
            root->add_child(internal1);
            root->add_child(internal2);

            internal0->set_height_parameter(internal1->get_height_parameter());
            internal2->set_height_parameter(internal1->get_height_parameter());

            BaseTree<Node> tree(root);

            tree.collision_node_swap(rng, 1, 0);

            REQUIRE(tree.tree_is_valid());
            REQUIRE(root->get_number_of_children() == 3);
            REQUIRE(internal0->get_number_of_children() == 2);
            REQUIRE(internal1->get_number_of_children() == 2);
            REQUIRE(internal2->get_number_of_children() == 2);
            std::vector< std::shared_ptr<Node> > expected_leaves {leaf0, leaf1, leaf2, leaf3, leaf4, leaf5};
            std::vector< std::shared_ptr<Node> > leaves = root->get_leaves();
            REQUIRE(std::is_permutation(
                        leaves.begin(), leaves.end(),
                        expected_leaves.begin()));
            if (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf1) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf3) &&
                    internal2->is_child(leaf4) &&
                    internal2->is_child(leaf5)) {
                ++count_01_23_45;
            }
            if (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf5) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf1) &&
                    internal2->is_child(leaf4) &&
                    internal2->is_child(leaf3)) {
                ++count_05_21_43;
            }
            if (
                    (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf5) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf3) &&
                    internal2->is_child(leaf4) &&
                    internal2->is_child(leaf1)) ||
                    (
                    internal2->is_child(leaf0) &&
                    internal2->is_child(leaf5) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf3) &&
                    internal0->is_child(leaf4) &&
                    internal0->is_child(leaf1))
                ) {
                ++count_05_23_41;
            }
            if (
                    (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf4) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf3) &&
                    internal2->is_child(leaf1) &&
                    internal2->is_child(leaf5)) ||
                    (
                    internal2->is_child(leaf0) &&
                    internal2->is_child(leaf4) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf3) &&
                    internal0->is_child(leaf1) &&
                    internal0->is_child(leaf5))
                ) {
                ++count_04_23_15;
            }
            if (
                    (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf1) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf5) &&
                    internal2->is_child(leaf4) &&
                    internal2->is_child(leaf3)) ||
                    (
                    internal2->is_child(leaf2) &&
                    internal2->is_child(leaf5) &&
                    internal1->is_child(leaf4) &&
                    internal1->is_child(leaf3) &&
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf1))
                ) {
                ++count_01_25_43;
            }
        }
        REQUIRE(count_01_23_45 == 0);
        REQUIRE(count_05_21_43 == 0);
        double freq_05_23_41 = (count_05_23_41 / (double)nsamples);
        double freq_04_23_15 = (count_04_23_15 / (double)nsamples);
        double freq_01_25_43 = (count_01_25_43 / (double)nsamples);
        std::cout << "freq of ((0,5),(2,3),(4,1)): " << freq_05_23_41 << "\n";
        std::cout << "freq of ((0,4),(2,3),(1,5)): " << freq_04_23_15 << "\n";
        std::cout << "freq of ((0,1),(2,5),(4,3)): " << freq_01_25_43 << "\n";

        double eps = 0.005;
        ///////////////////////////////////////////////////////////////////////
        // Prob outcome = Prob(picking parents to swap) *
        //                Prob(picking children of selected parents) * 
        //                Number of ways to get outcome
        // Prob(picking parents to swap) = 1/(3 choose 2) = 1 / 3!/(2!1!) = 1/3
        // Prob(picking children of selected parents) = (1/2) * (1/2)
        // # of ways = 2
        // Prob = 1/3 * 1/4 * 2 = 1/6
        REQUIRE(freq_05_23_41 == Approx(1.0/6.0).epsilon(eps));
        REQUIRE(freq_04_23_15 == Approx(1.0/6.0).epsilon(eps));
        REQUIRE(freq_01_25_43 == Approx(1.0/6.0).epsilon(eps));
    }
}

TEST_CASE("Testing BaseTree::slide_bump_permute_height 9 leaf 4 colliders", "[BaseTree]") {
    SECTION("Testing slide_bump_permute_height") {
        unsigned int nsamples = 1000000;
        RandomNumberGenerator rng = RandomNumberGenerator(111);
        unsigned int count_01_234_5678_9_10 = 0; // no change
        unsigned int count_10_214_03_5_9678 = 0; // every clade changes on both bumps
        unsigned int count_03_5_9678_214_10 = 0; // 2nd bump symmetric swap
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.5);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
            std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.5);
            std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.5);
            std::shared_ptr<Node> internal3 = std::make_shared<Node>("internal3", 1.0);
            std::shared_ptr<Node> internal4 = std::make_shared<Node>("internal4", 1.0);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            leaf0->fix_node_height();
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            leaf1->fix_node_height();
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            leaf2->fix_node_height();
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
            leaf3->fix_node_height();
            std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
            leaf4->fix_node_height();
            std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
            leaf5->fix_node_height();
            std::shared_ptr<Node> leaf6 = std::make_shared<Node>("leaf6", 0.0);
            leaf6->fix_node_height();
            std::shared_ptr<Node> leaf7 = std::make_shared<Node>("leaf7", 0.0);
            leaf7->fix_node_height();
            std::shared_ptr<Node> leaf8 = std::make_shared<Node>("leaf8", 0.0);
            leaf8->fix_node_height();
            std::shared_ptr<Node> leaf9 = std::make_shared<Node>("leaf9", 0.0);
            leaf9->fix_node_height();
            std::shared_ptr<Node> leaf10 = std::make_shared<Node>("leaf10", 0.0);
            leaf10->fix_node_height();

            internal0->add_child(leaf0);
            internal0->add_child(leaf1);

            internal1->add_child(leaf2);
            internal1->add_child(leaf3);
            internal1->add_child(leaf4);

            internal2->add_child(leaf5);
            internal2->add_child(leaf6);
            internal2->add_child(leaf7);
            internal2->add_child(leaf8);

            internal3->add_child(internal0);
            internal3->add_child(internal1);

            internal4->add_child(internal2);
            internal4->add_child(leaf9);

            root->add_child(internal3);
            root->add_child(internal4);
            root->add_child(leaf10);

            internal0->set_height_parameter(internal1->get_height_parameter());
            internal2->set_height_parameter(internal1->get_height_parameter());

            internal3->set_height_parameter(internal4->get_height_parameter());

            BaseTree<Node> tree(root);

            REQUIRE(tree.get_number_of_node_heights() == 3);

            tree.slide_bump_permute_height(rng, 0, 2.0);

            REQUIRE(tree.get_number_of_node_heights() == 3);

            REQUIRE(tree.get_height(0) == 1.0);
            REQUIRE(tree.get_height(1) == 1.5);
            REQUIRE(tree.get_height(2) == 2.0);
            REQUIRE(tree.tree_is_valid());
            REQUIRE(root->get_number_of_children() == 3);
            REQUIRE(internal0->get_number_of_children() == 2);
            REQUIRE(internal1->get_number_of_children() == 3);
            REQUIRE(internal2->get_number_of_children() == 4);
            REQUIRE(internal3->get_number_of_children() == 2);
            REQUIRE(internal4->get_number_of_children() == 2);
            std::vector< std::shared_ptr<Node> > expected_leaves {leaf0, leaf1, leaf2, leaf3, leaf4, leaf5, leaf6, leaf7, leaf8, leaf9, leaf10};
            std::vector< std::shared_ptr<Node> > leaves = root->get_leaves();
            REQUIRE(std::is_permutation(
                        leaves.begin(), leaves.end(),
                        expected_leaves.begin()));
            if (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf1) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf3) &&
                    internal1->is_child(leaf4) &&
                    internal2->is_child(leaf5) &&
                    internal2->is_child(leaf6) &&
                    internal2->is_child(leaf7) &&
                    internal2->is_child(leaf8) &&
                    internal3->is_child(internal0) &&
                    internal3->is_child(internal1) &&
                    internal4->is_child(internal2) &&
                    internal4->is_child(leaf9) &&
                    root->is_child(internal3) &&
                    root->is_child(internal4) &&
                    root->is_child(leaf10)
                    ) {
                ++count_01_234_5678_9_10;
            }
            if (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf3) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf1) &&
                    internal1->is_child(leaf4) &&
                    internal2->is_child(leaf9) &&
                    internal2->is_child(leaf6) &&
                    internal2->is_child(leaf7) &&
                    internal2->is_child(leaf8) &&
                    internal3->is_child(leaf10) &&
                    internal3->is_child(internal1) &&
                    internal4->is_child(internal0) &&
                    internal4->is_child(leaf5) &&
                    root->is_child(internal3) &&
                    root->is_child(internal4) &&
                    root->is_child(internal2)
                    ) {
                ++count_10_214_03_5_9678;
            }
            if (
                    (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf3) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf1) &&
                    internal1->is_child(leaf4) &&
                    internal2->is_child(leaf9) &&
                    internal2->is_child(leaf6) &&
                    internal2->is_child(leaf7) &&
                    internal2->is_child(leaf8) &&
                    internal3->is_child(internal0) &&
                    internal3->is_child(leaf5) &&
                    internal4->is_child(internal2) &&
                    internal4->is_child(internal1) &&
                    root->is_child(internal3) &&
                    root->is_child(internal4) &&
                    root->is_child(leaf10)) ||
                    (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf3) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf1) &&
                    internal1->is_child(leaf4) &&
                    internal2->is_child(leaf9) &&
                    internal2->is_child(leaf6) &&
                    internal2->is_child(leaf7) &&
                    internal2->is_child(leaf8) &&
                    internal4->is_child(internal0) &&
                    internal4->is_child(leaf5) &&
                    internal3->is_child(internal2) &&
                    internal3->is_child(internal1) &&
                    root->is_child(internal3) &&
                    root->is_child(internal4) &&
                    root->is_child(leaf10))
                    ) {
                ++count_03_5_9678_214_10;
            }
        }
        double freq_01_234_5678_9_10 = (count_01_234_5678_9_10 / (double)nsamples);
        double freq_10_214_03_5_9678 = (count_10_214_03_5_9678 / (double)nsamples);
        double freq_03_5_9678_214_10 = (count_03_5_9678_214_10 / (double)nsamples);
        std::cout << "freq of (((0,1),(2,3,4)),((5,6,7,8),9),10): " << freq_01_234_5678_9_10 << "\n";
        std::cout << "freq of ((10,(2,1,4)),((0,3),5),(9,6,7,8)): " << freq_10_214_03_5_9678 << "\n";
        std::cout << "freq of (((0,3),5),((9,6,7,8),(2,1,4)),10): " << freq_03_5_9678_214_10 << "\n";

        double eps = 0.0001;
        ///////////////////////////////////////////////////////////////////////
        // How I like to think about the probability of any swap result.
        // First, once you have a pool of N nodes, the number of unique
        // ways to sample them without replacement is N!, and so each
        // way of sampling the nodes in the pool has probability 1 / N!
        //
        // What differs between different outcomes is the probability of
        // getting a pool of nodes that is consistent with the outcome.
        // For a node that ends up with a different child, X, from the pool
        // than it added to it, the probability of this is simpy
        // (1 / the number of possible children X's parent could have contributed to the pool)
        //
        // For a node that ends up unchanged (i.e., it contributes a child to
        // the pool, only to get it back), the probability of a pool consistent
        // with that outcome is 1.0.
        // This is because the parent node has to contribute a node to the
        // pool, and no matter which one it contributes, it is always possible
        // to sample it back out of the pool.
        // For example, a bump node with children A and B, will add A with prob
        // 1/2 and B with prob 1/2, but either way, it could sample A / B back
        // out of the pool again, so every possible pool is consistent with an
        // outcome of an unchanged bump node.
        ///////////////////////////////////////////////////////////////////////
        //
        // BUMP 1
        // p(pool internal 3) * p(draws internal 3) = 1.0 * (1/2!) = 1/2
        // p(pool internal 4) * p(draws internal 4) = 1.0 * (1/2!) = 1/2
        // prob = 1/4
        // All possible pools can lead to getting the same state back, so p(pool) = 1.0
        // BUMP 2
        // p(pool) * p(draws) = 1.0 * (1/3!) = 1.0 * 1/6
        // All possible pools can lead to getting the same state back, so
        // p(pool) = 1.0
        // PROB of BOTH = 1/4 * 1/6 = 1/24
        REQUIRE(freq_01_234_5678_9_10 == Approx(1.0/24.0).epsilon(eps*5));

        // BUMP 1
        // p(pool) = p(1 in int 3 pool) * p(3 in int 3 pool) * p(5 in int 4 pool) * p(9 in int 4 pool)
        //         = (1/2) * (1/3) * (1/4) * 1 = 1/24
        // p(outcome) = p(pool) * p(draw) = 1/24 * ( (1/2!) * (1/2!) )
        //            = 1/24 * 1/4 = 1/96
        // BUMP 2
        // p(pool) = p(clade 03 in pool) * p(clade 9678 in pool) * p(10 in pool)
        //         = 1/2 * 1/2 * 1 = 1/4
        // p(draws) = 1/3! = 1/6
        // p(pool) * p(draws) = 1/4 * 1/6 = 1/24
        // Prob of BOTH = 1/96 * 1/24 = 1/2304
        REQUIRE(freq_10_214_03_5_9678 == Approx(1.0/2304.0).epsilon(eps));

        // BUMP 1
        // p(pool) = p(1 in int 3 pool) * p(3 in int 3 pool) * p(5 in int 4 pool) * p(9 in int 4 pool)
        //         = (1/2) * (1/3) * (1/4) * 1 = 1/24
        // p(outcome) = p(pool) * p(draw) = 1/24 * ( (1/2!) * (1/2!) )
        //            = 1/24 * 1/4 = 1/96
        // BUMP 2
        // p(pool) = p(clade 03 or 214 in pool) * p(corresponding clade in pool) * p(10 in pool)
        //         = 1 * 1/2 * 1 = 1/2
        // p(draws) = 1/3! = 1/6
        // p(pool) * p(draws) = 1/2 * 1/6 = 1/12
        // Prob of BOTH = 1/96 * 1/12 = 1/1152
        REQUIRE(freq_03_5_9678_214_10 == Approx(1.0/1152.0).epsilon(eps));
    }

    SECTION("Testing reverse slide_bump_permute_height") {
        unsigned int nsamples = 1000000;
        RandomNumberGenerator rng = RandomNumberGenerator(111);
        unsigned int count_01_234_5678_9_10 = 0; // reverse move
        unsigned int count_03_5_9678_214_10 = 0; // no move
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::shared_ptr<Node> root = std::make_shared<Node>("root", 2.0);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 1.0);
            std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 1.0);
            std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 1.0);
            std::shared_ptr<Node> internal3 = std::make_shared<Node>("internal3", 1.5);
            std::shared_ptr<Node> internal4 = std::make_shared<Node>("internal4", 1.5);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            leaf0->fix_node_height();
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            leaf1->fix_node_height();
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            leaf2->fix_node_height();
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
            leaf3->fix_node_height();
            std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
            leaf4->fix_node_height();
            std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
            leaf5->fix_node_height();
            std::shared_ptr<Node> leaf6 = std::make_shared<Node>("leaf6", 0.0);
            leaf6->fix_node_height();
            std::shared_ptr<Node> leaf7 = std::make_shared<Node>("leaf7", 0.0);
            leaf7->fix_node_height();
            std::shared_ptr<Node> leaf8 = std::make_shared<Node>("leaf8", 0.0);
            leaf8->fix_node_height();
            std::shared_ptr<Node> leaf9 = std::make_shared<Node>("leaf9", 0.0);
            leaf9->fix_node_height();
            std::shared_ptr<Node> leaf10 = std::make_shared<Node>("leaf10", 0.0);
            leaf10->fix_node_height();

            internal0->add_child(leaf0);
            internal0->add_child(leaf3);

            internal1->add_child(leaf2);
            internal1->add_child(leaf1);
            internal1->add_child(leaf4);

            internal2->add_child(leaf9);
            internal2->add_child(leaf6);
            internal2->add_child(leaf7);
            internal2->add_child(leaf8);

            internal3->add_child(internal0);
            internal3->add_child(leaf5);

            internal4->add_child(internal1);
            internal4->add_child(internal2);

            root->add_child(internal3);
            root->add_child(internal4);
            root->add_child(leaf10);

            internal0->set_height_parameter(internal1->get_height_parameter());
            internal2->set_height_parameter(internal1->get_height_parameter());

            internal3->set_height_parameter(internal4->get_height_parameter());

            BaseTree<Node> tree(root);

            REQUIRE(tree.get_number_of_node_heights() == 3);

            tree.slide_bump_permute_height(rng, 2, 0.5);

            REQUIRE(tree.get_number_of_node_heights() == 3);

            REQUIRE(tree.get_height(0) == 0.5);
            REQUIRE(tree.get_height(1) == 1.0);
            REQUIRE(tree.get_height(2) == 1.5);

            REQUIRE(tree.tree_is_valid());
            REQUIRE(root->get_number_of_children() == 3);
            REQUIRE(internal0->get_number_of_children() == 2);
            REQUIRE(internal1->get_number_of_children() == 3);
            REQUIRE(internal2->get_number_of_children() == 4);
            REQUIRE(internal3->get_number_of_children() == 2);
            REQUIRE(internal4->get_number_of_children() == 2);
            std::vector< std::shared_ptr<Node> > expected_leaves {leaf0, leaf1, leaf2, leaf3, leaf4, leaf5, leaf6, leaf7, leaf8, leaf9, leaf10};
            std::vector< std::shared_ptr<Node> > leaves = root->get_leaves();
            REQUIRE(std::is_permutation(
                        leaves.begin(), leaves.end(),
                        expected_leaves.begin()));
            if (
                    (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf1) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf3) &&
                    internal1->is_child(leaf4) &&
                    internal2->is_child(leaf5) &&
                    internal2->is_child(leaf6) &&
                    internal2->is_child(leaf7) &&
                    internal2->is_child(leaf8) &&
                    internal3->is_child(internal0) &&
                    internal3->is_child(internal1) &&
                    internal4->is_child(internal2) &&
                    internal4->is_child(leaf9) &&
                    root->is_child(internal3) &&
                    root->is_child(internal4) &&
                    root->is_child(leaf10)
                    ) ||
                    (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf1) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf3) &&
                    internal1->is_child(leaf4) &&
                    internal2->is_child(leaf5) &&
                    internal2->is_child(leaf6) &&
                    internal2->is_child(leaf7) &&
                    internal2->is_child(leaf8) &&
                    internal4->is_child(internal0) &&
                    internal4->is_child(internal1) &&
                    internal3->is_child(internal2) &&
                    internal3->is_child(leaf9) &&
                    root->is_child(internal3) &&
                    root->is_child(internal4) &&
                    root->is_child(leaf10)
                    )
                ) {
                ++count_01_234_5678_9_10;
            }
            if (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf3) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf1) &&
                    internal1->is_child(leaf4) &&
                    internal2->is_child(leaf9) &&
                    internal2->is_child(leaf6) &&
                    internal2->is_child(leaf7) &&
                    internal2->is_child(leaf8) &&
                    internal3->is_child(internal0) &&
                    internal3->is_child(leaf5) &&
                    internal4->is_child(internal2) &&
                    internal4->is_child(internal1) &&
                    root->is_child(internal3) &&
                    root->is_child(internal4) &&
                    root->is_child(leaf10)
                    ) {
                ++count_03_5_9678_214_10;
            }
        }
        double freq_01_234_5678_9_10 = (count_01_234_5678_9_10 / (double)nsamples);
        double freq_03_5_9678_214_10 = (count_03_5_9678_214_10 / (double)nsamples);
        std::cout << "freq of (((0,1),(2,3,4)),((5,6,7,8),9),10): " << freq_01_234_5678_9_10 << "\n";
        std::cout << "freq of (((0,3),5),((9,6,7,8),(2,1,4)),10): " << freq_03_5_9678_214_10 << "\n";

        double eps = 0.0001;
        ///////////////////////////////////////////////////////////////////////
        // How I like to think about the probability of any swap result.
        // First, once you have a pool of N nodes, the number of unique
        // ways to sample them without replacement is N!, and so each
        // way of sampling the nodes in the pool has probability 1 / N!
        //
        // What differs between different outcomes is the probability of
        // getting a pool of nodes that is consistent with the outcome.
        // For a node that ends up with a different child, X, from the pool
        // than it added to it, the probability of this is simpy
        // (1 / the number of possible children X's parent could have contributed to the pool)
        //
        // For a node that ends up unchanged (i.e., it contributes a child to
        // the pool, only to get it back), the probability of a pool consistent
        // with that outcome is 1.0.
        // This is because the parent node has to contribute a node to the
        // pool, and no matter which one it contributes, it is always possible
        // to sample it back out of the pool.
        // For example, a bump node with children A and B, will add A with prob
        // 1/2 and B with prob 1/2, but either way, it could sample A / B back
        // out of the pool again, so every possible pool is consistent with an
        // outcome of an unchanged bump node.
        ///////////////////////////////////////////////////////////////////////
        //
        // BUMP 1
        // p(pool) * p(draws) = 1.0 * (1/3!) = 1.0 * 1/6
        // All possible pools can lead to getting the same state back, so
        // p(pool) = 1.0
        // BUMP 2
        // p(pool internal 3) * p(draws internal 3) = 1.0 * (1/2!) = 1/2
        // p(pool internal 4) * p(draws internal 4) = 1.0 * (1/2!) = 1/2
        // prob = 1/4
        // All possible pools can lead to getting the same state back, so p(pool) = 1.0
        // PROB of BOTH = 1/4 * 1/6 = 1/24
        REQUIRE(freq_03_5_9678_214_10 == Approx(1.0/24.0).epsilon(eps*5));

        // BUMP 1
        // p(pool) = p(clade 03 or 5 in pool) * p(corresponding clade in pool) * p(10 in pool)
        //         = 1 * 1/2 * 1 = 1/2
        // p(draws) = 1/3! = 1/6
        // p(pool) * p(draws) = 1/2 * 1/6 = 1/12
        // BUMP 2
        // p(pool) = p(1 in pool) * p(3 in pool) * p(5 in pool) * p(9 in pool)
        //         = (1/3) * (1/2) * 1 * (1/4) = 1/24
        // p(draw) = 1/2! for both colliding nodes
        // p(outcome) = p(pool) * p(draw) = 1/24 * ( (1/2!) * (1/2!) )
        //            = 1/24 * 1/4 = 1/96
        // Prob of BOTH = 1/96 * 1/12 = 1/1152
        REQUIRE(freq_01_234_5678_9_10 == Approx(1.0/1152.0).epsilon(eps));
    }
}

TEST_CASE("Testing BaseTree::slide_bump_swap_height 9 leaf 4 colliders", "[BaseTree]") {
    SECTION("Testing slide_bump_swap_height") {
        unsigned int nsamples = 200000;
        RandomNumberGenerator rng = RandomNumberGenerator(111);
        unsigned int count_01_234_5678_9_10 = 0; // no change
        unsigned int count_10_214_03_5_9678 = 0; // every clade changes on both bumps
        unsigned int count_03_5_9678_214_10 = 0; // 2nd bump symmetric swap
        unsigned int count_02_10_5679_8_134 = 0; // 2nd bump asymmetric swap
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.5);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
            std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.5);
            std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.5);
            std::shared_ptr<Node> internal3 = std::make_shared<Node>("internal3", 1.0);
            std::shared_ptr<Node> internal4 = std::make_shared<Node>("internal4", 1.0);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            leaf0->fix_node_height();
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            leaf1->fix_node_height();
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            leaf2->fix_node_height();
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
            leaf3->fix_node_height();
            std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
            leaf4->fix_node_height();
            std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
            leaf5->fix_node_height();
            std::shared_ptr<Node> leaf6 = std::make_shared<Node>("leaf6", 0.0);
            leaf6->fix_node_height();
            std::shared_ptr<Node> leaf7 = std::make_shared<Node>("leaf7", 0.0);
            leaf7->fix_node_height();
            std::shared_ptr<Node> leaf8 = std::make_shared<Node>("leaf8", 0.0);
            leaf8->fix_node_height();
            std::shared_ptr<Node> leaf9 = std::make_shared<Node>("leaf9", 0.0);
            leaf9->fix_node_height();
            std::shared_ptr<Node> leaf10 = std::make_shared<Node>("leaf10", 0.0);
            leaf10->fix_node_height();

            internal0->add_child(leaf0);
            internal0->add_child(leaf1);

            internal1->add_child(leaf2);
            internal1->add_child(leaf3);
            internal1->add_child(leaf4);

            internal2->add_child(leaf5);
            internal2->add_child(leaf6);
            internal2->add_child(leaf7);
            internal2->add_child(leaf8);

            internal3->add_child(internal0);
            internal3->add_child(internal1);

            internal4->add_child(internal2);
            internal4->add_child(leaf9);

            root->add_child(internal3);
            root->add_child(internal4);
            root->add_child(leaf10);

            internal0->set_height_parameter(internal1->get_height_parameter());
            internal2->set_height_parameter(internal1->get_height_parameter());

            internal3->set_height_parameter(internal4->get_height_parameter());

            BaseTree<Node> tree(root);

            REQUIRE(tree.get_number_of_node_heights() == 3);

            tree.slide_bump_swap_height(rng, 0, 2.0);

            REQUIRE(tree.get_number_of_node_heights() == 3);

            REQUIRE(tree.get_height(0) == 1.0);
            REQUIRE(tree.get_height(1) == 1.5);
            REQUIRE(tree.get_height(2) == 2.0);
            REQUIRE(tree.tree_is_valid());
            REQUIRE(root->get_number_of_children() == 3);
            REQUIRE(internal0->get_number_of_children() == 2);
            REQUIRE(internal1->get_number_of_children() == 3);
            REQUIRE(internal2->get_number_of_children() == 4);
            REQUIRE(internal3->get_number_of_children() == 2);
            REQUIRE(internal4->get_number_of_children() == 2);
            std::vector< std::shared_ptr<Node> > expected_leaves {leaf0, leaf1, leaf2, leaf3, leaf4, leaf5, leaf6, leaf7, leaf8, leaf9, leaf10};
            std::vector< std::shared_ptr<Node> > leaves = root->get_leaves();
            REQUIRE(std::is_permutation(
                        leaves.begin(), leaves.end(),
                        expected_leaves.begin()));
            if (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf1) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf3) &&
                    internal1->is_child(leaf4) &&
                    internal2->is_child(leaf5) &&
                    internal2->is_child(leaf6) &&
                    internal2->is_child(leaf7) &&
                    internal2->is_child(leaf8) &&
                    internal3->is_child(internal0) &&
                    internal3->is_child(internal1) &&
                    internal4->is_child(internal2) &&
                    internal4->is_child(leaf9) &&
                    root->is_child(internal3) &&
                    root->is_child(internal4) &&
                    root->is_child(leaf10)
                    ) {
                ++count_01_234_5678_9_10;
            }
            if (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf3) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf1) &&
                    internal1->is_child(leaf4) &&
                    internal2->is_child(leaf9) &&
                    internal2->is_child(leaf6) &&
                    internal2->is_child(leaf7) &&
                    internal2->is_child(leaf8) &&
                    internal3->is_child(leaf10) &&
                    internal3->is_child(internal1) &&
                    internal4->is_child(internal0) &&
                    internal4->is_child(leaf5) &&
                    root->is_child(internal3) &&
                    root->is_child(internal4) &&
                    root->is_child(internal2)
                    ) {
                ++count_10_214_03_5_9678;
            }
            if (
                    (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf3) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf1) &&
                    internal1->is_child(leaf4) &&
                    internal2->is_child(leaf9) &&
                    internal2->is_child(leaf6) &&
                    internal2->is_child(leaf7) &&
                    internal2->is_child(leaf8) &&
                    internal3->is_child(internal0) &&
                    internal3->is_child(leaf5) &&
                    internal4->is_child(internal2) &&
                    internal4->is_child(internal1) &&
                    root->is_child(internal3) &&
                    root->is_child(internal4) &&
                    root->is_child(leaf10)) ||
                    (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf3) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf1) &&
                    internal1->is_child(leaf4) &&
                    internal2->is_child(leaf9) &&
                    internal2->is_child(leaf6) &&
                    internal2->is_child(leaf7) &&
                    internal2->is_child(leaf8) &&
                    internal4->is_child(internal0) &&
                    internal4->is_child(leaf5) &&
                    internal3->is_child(internal2) &&
                    internal3->is_child(internal1) &&
                    root->is_child(internal3) &&
                    root->is_child(internal4) &&
                    root->is_child(leaf10))
                    ) {
                ++count_03_5_9678_214_10;
            }
            if (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf2) &&
                    internal1->is_child(leaf1) &&
                    internal1->is_child(leaf3) &&
                    internal1->is_child(leaf4) &&
                    internal2->is_child(leaf5) &&
                    internal2->is_child(leaf6) &&
                    internal2->is_child(leaf7) &&
                    internal2->is_child(leaf9) &&
                    internal3->is_child(internal0) &&
                    internal3->is_child(leaf10) &&
                    internal4->is_child(internal2) &&
                    internal4->is_child(leaf8) &&
                    root->is_child(internal3) &&
                    root->is_child(internal4) &&
                    root->is_child(internal1)
                    ) {
                ++count_02_10_5679_8_134;
            }
        }
        REQUIRE(count_01_234_5678_9_10 == 0);
        REQUIRE(count_10_214_03_5_9678 == 0);
        double freq_03_5_9678_214_10 = (count_03_5_9678_214_10 / (double)nsamples);
        double freq_02_10_5679_8_134 = (count_02_10_5679_8_134 / (double)nsamples);
        std::cout << "freq of (((0,3),5),((9,6,7,8),(2,1,4)),10): " << freq_03_5_9678_214_10 << "\n";
        std::cout << "freq of (((0,2),10),((5,6,7,9),8),(1,3,4)): " << freq_02_10_5679_8_134 << "\n";

        double eps = 0.0002;
        ///////////////////////////////////////////////////////////////////////
        // BUMP 1
        // p(int 3) = p(parents selected) * p(children of parents selected) * # of ways
        //          = 1 * (1/2) * (1/3) * 1 = 1/6
        // p(int 4) = p(parents selected) * p(children of parents selected) * # of ways
        //          = 1 * 1 * (1/4) * 1 = 1/4
        // p(outcome) = 1/6 * 1/4 = 1/24
        // BUMP 2
        // p(root) = p(parents selected) * p(children of parents selected) * # of ways
        //         =  (1 / (3 choose 2)) *  ( 1/2 * 1/2 ) * 2 [[ internals 3 and 4 can swap in 2 ways ]]
        //         = 1/3 * 1/2 = 1/6
        // Prob of BOTH = 1/24 * 1/6 = 1/144
        REQUIRE(freq_03_5_9678_214_10 == Approx(1.0/144.0).epsilon(eps));

        ///////////////////////////////////////////////////////////////////////
        // BUMP 1
        // p(int 3) = p(parents selected) * p(children of parents selected) * # of ways
        //          = 1 * (1/2) * (1/3) * 1 = 1/6
        // p(int 4) = p(parents selected) * p(children of parents selected) * # of ways
        //          = 1 * 1 * (1/4) * 1 = 1/4
        // p(outcome) = 1/6 * 1/4 = 1/24
        // BUMP 2
        // p(root) = p(parents selected) * p(children of parents selected) * # of ways
        //         =  (1 / (3 choose 2)) *  ( 1/2 * 1 ) * 1
        //         = 1/3 * 1/2 = 1/6
        // Prob of BOTH = 1/24 * 1/6 = 1/144
        REQUIRE(freq_02_10_5679_8_134 == Approx(1.0/144.0).epsilon(eps));
    }

    SECTION("Testing reverse slide_bump_swap_height") {
        unsigned int nsamples = 200000;
        RandomNumberGenerator rng = RandomNumberGenerator(111);
        unsigned int count_01_234_5678_9_10 = 0; // reverse move
        unsigned int count_03_5_9678_214_10 = 0; // no move
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::shared_ptr<Node> root = std::make_shared<Node>("root", 2.0);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 1.0);
            std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 1.0);
            std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 1.0);
            std::shared_ptr<Node> internal3 = std::make_shared<Node>("internal3", 1.5);
            std::shared_ptr<Node> internal4 = std::make_shared<Node>("internal4", 1.5);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            leaf0->fix_node_height();
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            leaf1->fix_node_height();
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            leaf2->fix_node_height();
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
            leaf3->fix_node_height();
            std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
            leaf4->fix_node_height();
            std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
            leaf5->fix_node_height();
            std::shared_ptr<Node> leaf6 = std::make_shared<Node>("leaf6", 0.0);
            leaf6->fix_node_height();
            std::shared_ptr<Node> leaf7 = std::make_shared<Node>("leaf7", 0.0);
            leaf7->fix_node_height();
            std::shared_ptr<Node> leaf8 = std::make_shared<Node>("leaf8", 0.0);
            leaf8->fix_node_height();
            std::shared_ptr<Node> leaf9 = std::make_shared<Node>("leaf9", 0.0);
            leaf9->fix_node_height();
            std::shared_ptr<Node> leaf10 = std::make_shared<Node>("leaf10", 0.0);
            leaf10->fix_node_height();

            internal0->add_child(leaf0);
            internal0->add_child(leaf3);

            internal1->add_child(leaf2);
            internal1->add_child(leaf1);
            internal1->add_child(leaf4);

            internal2->add_child(leaf9);
            internal2->add_child(leaf6);
            internal2->add_child(leaf7);
            internal2->add_child(leaf8);

            internal3->add_child(internal0);
            internal3->add_child(leaf5);

            internal4->add_child(internal1);
            internal4->add_child(internal2);

            root->add_child(internal3);
            root->add_child(internal4);
            root->add_child(leaf10);

            internal0->set_height_parameter(internal1->get_height_parameter());
            internal2->set_height_parameter(internal1->get_height_parameter());

            internal3->set_height_parameter(internal4->get_height_parameter());

            BaseTree<Node> tree(root);

            REQUIRE(tree.get_number_of_node_heights() == 3);

            tree.slide_bump_swap_height(rng, 2, 0.5);

            REQUIRE(tree.get_number_of_node_heights() == 3);

            REQUIRE(tree.get_height(0) == 0.5);
            REQUIRE(tree.get_height(1) == 1.0);
            REQUIRE(tree.get_height(2) == 1.5);

            REQUIRE(tree.tree_is_valid());
            REQUIRE(root->get_number_of_children() == 3);
            REQUIRE(internal0->get_number_of_children() == 2);
            REQUIRE(internal1->get_number_of_children() == 3);
            REQUIRE(internal2->get_number_of_children() == 4);
            REQUIRE(internal3->get_number_of_children() == 2);
            REQUIRE(internal4->get_number_of_children() == 2);
            std::vector< std::shared_ptr<Node> > expected_leaves {leaf0, leaf1, leaf2, leaf3, leaf4, leaf5, leaf6, leaf7, leaf8, leaf9, leaf10};
            std::vector< std::shared_ptr<Node> > leaves = root->get_leaves();
            REQUIRE(std::is_permutation(
                        leaves.begin(), leaves.end(),
                        expected_leaves.begin()));
            if (
                    (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf1) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf3) &&
                    internal1->is_child(leaf4) &&
                    internal2->is_child(leaf5) &&
                    internal2->is_child(leaf6) &&
                    internal2->is_child(leaf7) &&
                    internal2->is_child(leaf8) &&
                    internal3->is_child(internal0) &&
                    internal3->is_child(internal1) &&
                    internal4->is_child(internal2) &&
                    internal4->is_child(leaf9) &&
                    root->is_child(internal3) &&
                    root->is_child(internal4) &&
                    root->is_child(leaf10)
                    ) ||
                    (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf1) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf3) &&
                    internal1->is_child(leaf4) &&
                    internal2->is_child(leaf5) &&
                    internal2->is_child(leaf6) &&
                    internal2->is_child(leaf7) &&
                    internal2->is_child(leaf8) &&
                    internal4->is_child(internal0) &&
                    internal4->is_child(internal1) &&
                    internal3->is_child(internal2) &&
                    internal3->is_child(leaf9) &&
                    root->is_child(internal3) &&
                    root->is_child(internal4) &&
                    root->is_child(leaf10)
                    )
                ) {
                ++count_01_234_5678_9_10;
            }
            if (
                    internal0->is_child(leaf0) &&
                    internal0->is_child(leaf3) &&
                    internal1->is_child(leaf2) &&
                    internal1->is_child(leaf1) &&
                    internal1->is_child(leaf4) &&
                    internal2->is_child(leaf9) &&
                    internal2->is_child(leaf6) &&
                    internal2->is_child(leaf7) &&
                    internal2->is_child(leaf8) &&
                    internal3->is_child(internal0) &&
                    internal3->is_child(leaf5) &&
                    internal4->is_child(internal2) &&
                    internal4->is_child(internal1) &&
                    root->is_child(internal3) &&
                    root->is_child(internal4) &&
                    root->is_child(leaf10)
                    ) {
                ++count_03_5_9678_214_10;
            }
        }
        REQUIRE(count_03_5_9678_214_10 == 0);
        double freq_01_234_5678_9_10 = (count_01_234_5678_9_10 / (double)nsamples);
        std::cout << "freq of (((0,1),(2,3,4)),((5,6,7,8),9),10): " << freq_01_234_5678_9_10 << "\n";

        double eps = 0.0002;
        ///////////////////////////////////////////////////////////////////////
        // BUMP 1
        // p(int 3) = p(parents selected) * p(children of parents selected) * # of ways
        //          = 1 * (1/2) * 1 * 1 = 1/2
        // p(int 4) = p(parents selected) * p(children of parents selected) * # of ways
        //          = 1 * (1/3) * (1/4) * 1 = 1/12
        // p(outcome) = 1/2 * 1/12 = 1/24
        // BUMP 2
        // p(root) = p(parents selected) * p(children of parents selected) * # of ways
        //         =  (1 / (3 choose 2)) *  ( 1/2 * 1/2 ) * 2 [[ internals 3 and 4 can swap in 2 ways ]]
        //         = 1/3 * 1/2 = 1/6
        // Prob of BOTH = 1/24 * 1/6 = 1/144
        REQUIRE(freq_01_234_5678_9_10 == Approx(1.0/144.0).epsilon(eps));
    }
}

TEST_CASE("Testing BaseTree store and restore", "[BaseTree]") {
    SECTION("Testing store-restore of state") {
        double height_lower_bound;
        unsigned int number_of_mapped_nodes;
        std::vector<unsigned int> mapped_polytomy_sizes;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.5);
        std::shared_ptr<ContinuousProbabilityDistribution> root_node_height_prior = std::make_shared<ExponentialDistribution>(10.0);
        root->set_node_height_prior(root_node_height_prior);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.5);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.5);
        std::shared_ptr<Node> internal3 = std::make_shared<Node>("internal3", 1.0);
        std::shared_ptr<Node> internal4 = std::make_shared<Node>("internal4", 1.0);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        leaf0->fix_node_height();
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        leaf1->fix_node_height();
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        leaf2->fix_node_height();
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        leaf3->fix_node_height();
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
        leaf4->fix_node_height();
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
        leaf5->fix_node_height();
        std::shared_ptr<Node> leaf6 = std::make_shared<Node>("leaf6", 0.0);
        leaf6->fix_node_height();
        std::shared_ptr<Node> leaf7 = std::make_shared<Node>("leaf7", 0.0);
        leaf7->fix_node_height();
        std::shared_ptr<Node> leaf8 = std::make_shared<Node>("leaf8", 0.0);
        leaf8->fix_node_height();
        std::shared_ptr<Node> leaf9 = std::make_shared<Node>("leaf9", 0.0);
        leaf9->fix_node_height();
        std::shared_ptr<Node> leaf10 = std::make_shared<Node>("leaf10", 0.0);
        leaf10->fix_node_height();

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);

        internal1->add_child(leaf2);
        internal1->add_child(leaf3);
        internal1->add_child(leaf4);

        internal2->add_child(leaf5);
        internal2->add_child(leaf6);
        internal2->add_child(leaf7);
        internal2->add_child(leaf8);

        internal3->add_child(internal0);
        internal3->add_child(internal1);

        internal4->add_child(internal2);
        internal4->add_child(leaf9);

        root->add_child(internal3);
        root->add_child(internal4);
        root->add_child(leaf10);

        internal0->set_height_parameter(internal1->get_height_parameter());
        internal2->set_height_parameter(internal1->get_height_parameter());

        internal3->set_height_parameter(internal4->get_height_parameter());

        std::vector<double> expected_node_heights {0.5, 1.0, 1.5};
        std::vector<double> expected_move_node_heights;

        BaseTree<Node> tree(root);

        std::cout << root.get() << "\n";
        std::cout << &tree.get_root() << "\n";
        REQUIRE(root.get() == &tree.get_root());

        std::string expected_tree_str = tree.to_parentheses();

        std::cout << "Starting tree:\n";
        std::cout << tree.to_parentheses() << "\n";

        tree.compute_log_likelihood_and_prior();
        double expected_lnl = tree.get_log_likelihood_value();
        double expected_ln_prior = tree.get_log_prior_density_value();

        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_number_of_node_heights() == 3);
        REQUIRE(tree.get_node_heights() == expected_node_heights);
        REQUIRE(tree.get_root().get_number_of_children() == 3);
        REQUIRE(tree.get_log_likelihood_value() == expected_lnl);
        REQUIRE(tree.get_log_prior_density_value() == expected_ln_prior);

        tree.store_state();
        tree.slide_bump_permute_height(rng, 0, 2.0);
        std::cout << "Tree after slide_bump_permute_height(rng, 0, 2.0):\n";
        std::cout << tree.to_parentheses() << "\n";
        REQUIRE(tree.to_parentheses() != expected_tree_str);
        REQUIRE(tree.tree_is_valid());
        expected_move_node_heights = {1.0, 1.5, 2.0};
        REQUIRE(tree.get_node_heights() == expected_move_node_heights);
        tree.compute_log_likelihood_and_prior();
        tree.restore_state();

        std::cout << "Tree after restore:\n";
        std::cout << tree.to_parentheses() << "\n";

        std::cout << root.get() << "\n";
        std::cout << &tree.get_root() << "\n";
        REQUIRE(root.get() != &tree.get_root());

        REQUIRE(tree.to_parentheses() == expected_tree_str);
        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_number_of_node_heights() == 3);
        REQUIRE(tree.get_node_heights() == expected_node_heights);
        REQUIRE(tree.get_root().get_number_of_children() == 3);
        REQUIRE(tree.get_log_likelihood_value() == expected_lnl);
        REQUIRE(tree.get_log_prior_density_value() == expected_ln_prior);

        tree.store_state();
        tree.slide_bump_permute_height(rng, 2, 0.1);
        std::cout << "Tree after slide_bump_permute_height(rng, 2, 0.1):\n";
        std::cout << tree.to_parentheses() << "\n";
        REQUIRE(tree.to_parentheses() != expected_tree_str);
        REQUIRE(tree.tree_is_valid());
        expected_move_node_heights = {0.1, 0.5, 1.0};
        REQUIRE(tree.get_node_heights() == expected_move_node_heights);
        tree.compute_log_likelihood_and_prior();
        tree.restore_state();

        std::cout << "Tree after restore:\n";
        std::cout << tree.to_parentheses() << "\n";

        std::cout << root.get() << "\n";
        std::cout << &tree.get_root() << "\n";
        REQUIRE(root.get() != &tree.get_root());

        REQUIRE(tree.to_parentheses() == expected_tree_str);
        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_number_of_node_heights() == 3);
        REQUIRE(tree.get_node_heights() == expected_node_heights);
        REQUIRE(tree.get_root().get_number_of_children() == 3);
        REQUIRE(tree.get_log_likelihood_value() == expected_lnl);
        REQUIRE(tree.get_log_prior_density_value() == expected_ln_prior);

        tree.store_state();
        tree.merge_node_height_up(0);
        std::cout << "Tree after merge_node_height_up(0):\n";
        std::cout << tree.to_parentheses() << "\n";
        REQUIRE(tree.to_parentheses() != expected_tree_str);
        REQUIRE(tree.tree_is_valid());
        expected_move_node_heights = {1.0, 1.5};
        REQUIRE(tree.get_node_heights() == expected_move_node_heights);
        tree.compute_log_likelihood_and_prior();
        tree.restore_state();

        std::cout << "Tree after restore:\n";
        std::cout << tree.to_parentheses() << "\n";

        std::cout << root.get() << "\n";
        std::cout << &tree.get_root() << "\n";
        REQUIRE(root.get() != &tree.get_root());

        REQUIRE(tree.to_parentheses() == expected_tree_str);
        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_number_of_node_heights() == 3);
        REQUIRE(tree.get_node_heights() == expected_node_heights);
        REQUIRE(tree.get_root().get_number_of_children() == 3);
        REQUIRE(tree.get_log_likelihood_value() == expected_lnl);
        REQUIRE(tree.get_log_prior_density_value() == expected_ln_prior);

        tree.store_state();
        tree.merge_node_height_up(1);
        std::cout << "Tree after merge_node_height_up(1):\n";
        std::cout << tree.to_parentheses() << "\n";
        REQUIRE(tree.to_parentheses() != expected_tree_str);
        REQUIRE(tree.tree_is_valid());
        expected_move_node_heights = {0.5, 1.5};
        REQUIRE(tree.get_node_heights() == expected_move_node_heights);
        tree.compute_log_likelihood_and_prior();
        tree.restore_state();

        std::cout << "Tree after restore:\n";
        std::cout << tree.to_parentheses() << "\n";

        std::cout << root.get() << "\n";
        std::cout << &tree.get_root() << "\n";
        REQUIRE(root.get() != &tree.get_root());

        REQUIRE(tree.to_parentheses() == expected_tree_str);
        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_number_of_node_heights() == 3);
        REQUIRE(tree.get_node_heights() == expected_node_heights);
        REQUIRE(tree.get_root().get_number_of_children() == 3);
        REQUIRE(tree.get_log_likelihood_value() == expected_lnl);
        REQUIRE(tree.get_log_prior_density_value() == expected_ln_prior);

        tree.store_state();
        tree.split_node_height_down(rng, 2,
                height_lower_bound,
                number_of_mapped_nodes,
                mapped_polytomy_sizes);
        std::cout << "Tree after split_node_height_down(rng, 2):\n";
        std::cout << tree.to_parentheses() << "\n";
        REQUIRE(tree.to_parentheses() != expected_tree_str);
        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_node_heights().at(0) == 0.5);
        REQUIRE(tree.get_node_heights().at(1) == 1.0);
        REQUIRE(tree.get_node_heights().at(2) > 1.0);
        REQUIRE(tree.get_node_heights().at(2) < 1.5);
        REQUIRE(tree.get_node_heights().at(3) == 1.5);
        tree.compute_log_likelihood_and_prior();
        tree.restore_state();

        std::cout << "Tree after restore:\n";
        std::cout << tree.to_parentheses() << "\n";

        std::cout << root.get() << "\n";
        std::cout << &tree.get_root() << "\n";
        REQUIRE(root.get() != &tree.get_root());

        REQUIRE(tree.to_parentheses() == expected_tree_str);
        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_number_of_node_heights() == 3);
        REQUIRE(tree.get_node_heights() == expected_node_heights);
        REQUIRE(tree.get_root().get_number_of_children() == 3);
        REQUIRE(tree.get_log_likelihood_value() == expected_lnl);
        REQUIRE(tree.get_log_prior_density_value() == expected_ln_prior);

        tree.store_state();
        tree.split_node_height_down(rng, 1,
                height_lower_bound,
                number_of_mapped_nodes,
                mapped_polytomy_sizes);
        std::cout << "Tree after split_node_height_down(rng, 1):\n";
        std::cout << tree.to_parentheses() << "\n";
        REQUIRE(tree.to_parentheses() != expected_tree_str);
        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_node_heights().at(0) == 0.5);
        REQUIRE(tree.get_node_heights().at(1) > 0.5);
        REQUIRE(tree.get_node_heights().at(1) < 1.0);
        REQUIRE(tree.get_node_heights().at(2) == 1.0);
        REQUIRE(tree.get_node_heights().at(3) == 1.5);
        tree.compute_log_likelihood_and_prior();
        tree.restore_state();

        std::cout << "Tree after restore:\n";
        std::cout << tree.to_parentheses() << "\n";

        std::cout << root.get() << "\n";
        std::cout << &tree.get_root() << "\n";
        REQUIRE(root.get() != &tree.get_root());

        REQUIRE(tree.to_parentheses() == expected_tree_str);
        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_number_of_node_heights() == 3);
        REQUIRE(tree.get_node_heights() == expected_node_heights);
        REQUIRE(tree.get_root().get_number_of_children() == 3);
        REQUIRE(tree.get_log_likelihood_value() == expected_lnl);
        REQUIRE(tree.get_log_prior_density_value() == expected_ln_prior);

        tree.store_state();
        tree.split_node_height_down(rng, 0,
                height_lower_bound,
                number_of_mapped_nodes,
                mapped_polytomy_sizes);
        std::cout << "Tree after split_node_height_down(rng, 0):\n";
        std::cout << tree.to_parentheses() << "\n";
        REQUIRE(tree.to_parentheses() != expected_tree_str);
        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_node_heights().at(0) < 0.5);
        REQUIRE(tree.get_node_heights().at(1) == 0.5);
        REQUIRE(tree.get_node_heights().at(2) == 1.0);
        REQUIRE(tree.get_node_heights().at(3) == 1.5);
        tree.compute_log_likelihood_and_prior();
        tree.restore_state();

        std::cout << "Tree after restore:\n";
        std::cout << tree.to_parentheses() << "\n";

        std::cout << root.get() << "\n";
        std::cout << &tree.get_root() << "\n";
        REQUIRE(root.get() != &tree.get_root());

        REQUIRE(tree.to_parentheses() == expected_tree_str);
        REQUIRE(tree.tree_is_valid());
        REQUIRE(tree.get_number_of_node_heights() == 3);
        REQUIRE(tree.get_node_heights() == expected_node_heights);
        REQUIRE(tree.get_root().get_number_of_children() == 3);
        REQUIRE(tree.get_log_likelihood_value() == expected_lnl);
        REQUIRE(tree.get_log_prior_density_value() == expected_ln_prior);
    }
}


TEST_CASE("Testing scaling of simulate_gene_tree for three species",
        "[PopulationTree]") {

    SECTION("Testing three species") {
        std::vector<unsigned int> Nes = {
                50000,
                10000,
                25000,
                100000,
                200000};
        unsigned int ngenomes = 10;
        std::vector<unsigned int> nlineages = {
                ngenomes,
                ngenomes,
                ngenomes,
                2,
                2};
        double mu = 1e-8;
        std::vector<double> thetas(Nes.size());
        for (unsigned int i = 0; i < Nes.size(); ++i) {
            thetas.at(i) = 4 * Nes.at(i) * mu;
        }
        double time_internal = 2e7;
        double time_root = 4e7;

        std::vector<double> expected_means(Nes.size());
        std::vector<double> expected_variances(Nes.size());
        for (unsigned int i = 0; i < Nes.size(); ++i) {
            expected_means.at(i) = thetas.at(i)* (1.0 - (1.0 / nlineages.at(i)));
            expected_variances.at(i) = expected_means.at(i) * expected_means.at(i);
        }
        expected_means.at(4) += (time_root * mu);
        expected_means.at(3) += (time_internal * mu);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.1);
        std::shared_ptr<PopulationNode> internal = std::make_shared<PopulationNode>(3, "internal 0", 0.05);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf 0", 0.0, ngenomes);
        leaf0->fix_node_height();
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf 1", 0.0, ngenomes);
        leaf1->fix_node_height();
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf 2", 0.0, ngenomes);
        leaf2->fix_node_height();

        internal->add_child(leaf0);
        internal->add_child(leaf1);
        root->add_child(internal);
        root->add_child(leaf2);


        PopulationTree tree(root,
                100,   // number of loci
                1,     // length of loci
                true); // validate data

        leaf0->set_population_size(Nes.at(0));
        leaf1->set_population_size(Nes.at(1));
        leaf2->set_population_size(Nes.at(2));
        internal->set_population_size(Nes.at(3));
        root->set_population_size(Nes.at(4));

        internal->set_height(time_internal);
        root->set_height(time_root);

        tree.estimate_mutation_rate();

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(mu);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        std::shared_ptr<GeneTreeSimNode> gtree;
        std::vector<SampleSummarizer<double> > height_summaries(Nes.size());
        SampleSummarizer<double> tree_length;

        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            REQUIRE(gtree->is_root());
            REQUIRE(gtree->get_number_of_children() == 2);
            REQUIRE(gtree->get_leaf_node_count() == ngenomes * 3);
            std::vector<double> oldest_coals(5, -1.0);
            const std::vector< std::shared_ptr<GeneTreeSimNode> >& gnodes = gtree->get_nodes();
            for (auto const & node_iter: gnodes) {
                if (node_iter->get_height() > oldest_coals.at(node_iter->get_population_index())) {
                    oldest_coals.at(node_iter->get_population_index()) = node_iter->get_height();
                }
            }
            for (unsigned int branch_idx = 0; branch_idx < oldest_coals.size(); ++branch_idx) {
                height_summaries.at(branch_idx).add_sample(oldest_coals.at(branch_idx));
            }
            tree_length.add_sample(gtree->get_clade_length());
        }

        for (unsigned int i = 0; i < Nes.size(); ++i) {
            std::cout << "\nnode: " << i << "\n";
            std::cout << "expected height: " << expected_means.at(i) << "\n";
            std::cout << "mean height: " << height_summaries.at(i).mean() << "\n";
            std::cout << "nsamples: " << height_summaries.at(i).sample_size() << std::endl;
        }
        for (unsigned int i = 0; i < Nes.size(); ++i) {
            REQUIRE(height_summaries.at(i).mean() == Approx(expected_means.at(i)).epsilon(0.0005));
            REQUIRE(height_summaries.at(i).variance() == Approx(expected_variances.at(i)).epsilon(0.0005));
        }



        tree.set_mutation_rate(mu * 2.0);

        std::vector<SampleSummarizer<double> > scale_height_summaries(Nes.size());
        SampleSummarizer<double> scale_tree_length;

        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            gtree->scale(0.5);
            REQUIRE(gtree->is_root());
            REQUIRE(gtree->get_number_of_children() == 2);
            REQUIRE(gtree->get_leaf_node_count() == 30);
            std::vector<double> oldest_coals(5, -1.0);
            const std::vector< std::shared_ptr<GeneTreeSimNode> >& scale_gnodes = gtree->get_nodes();
            for (auto const & node_iter: scale_gnodes) {
                if (node_iter->get_height() > oldest_coals.at(node_iter->get_population_index())) {
                    oldest_coals.at(node_iter->get_population_index()) = node_iter->get_height();
                }
            }
            for (unsigned int branch_idx = 0; branch_idx < oldest_coals.size(); ++branch_idx) {
                scale_height_summaries.at(branch_idx).add_sample(oldest_coals.at(branch_idx));
            }
            scale_tree_length.add_sample(gtree->get_clade_length());
        }

        for (unsigned int i = 0; i < Nes.size(); ++i) {
            REQUIRE(scale_height_summaries.at(i).mean() == Approx(expected_means.at(i)).epsilon(0.0005));
            REQUIRE(scale_height_summaries.at(i).variance() == Approx(expected_variances.at(i)).epsilon(0.0005));
        }

        REQUIRE(scale_tree_length.mean() == Approx(tree_length.mean()).epsilon(0.0001));
        REQUIRE(scale_tree_length.variance() == Approx(tree_length.variance()).epsilon(0.0001));
    }
}

TEST_CASE("Testing likelihood of PopulationTree with three-way polytomy at root", "[PopulationTree]") {

    SECTION("Testing constructor and likelihood calc") {
        std::shared_ptr<PopulationNode> root0 = std::make_shared<PopulationNode>(3, "root", 0.1);
        std::shared_ptr<PopulationNode> leaf00 = std::make_shared<PopulationNode>(0, "leaf 0", 0.0, 4);
        leaf00->fix_node_height();
        std::shared_ptr<PopulationNode> leaf01 = std::make_shared<PopulationNode>(1, "leaf 1", 0.0, 4);
        leaf01->fix_node_height();
        std::shared_ptr<PopulationNode> leaf02 = std::make_shared<PopulationNode>(2, "leaf 2", 0.0, 4);
        leaf02->fix_node_height();

        root0->add_child(leaf00);
        root0->add_child(leaf01);
        root0->add_child(leaf02);

        std::shared_ptr<PopulationNode> internal1 = std::make_shared<PopulationNode>(3, "internal 0", 0.1);
        std::shared_ptr<PopulationNode> root1 = std::make_shared<PopulationNode>(4, "root", 0.1);
        std::shared_ptr<PopulationNode> leaf10 = std::make_shared<PopulationNode>(0, "leaf 0", 0.0, 4);
        leaf10->fix_node_height();
        std::shared_ptr<PopulationNode> leaf11 = std::make_shared<PopulationNode>(1, "leaf 1", 0.0, 4);
        leaf11->fix_node_height();
        std::shared_ptr<PopulationNode> leaf12 = std::make_shared<PopulationNode>(2, "leaf 2", 0.0, 4);
        leaf12->fix_node_height();

        internal1->add_child(leaf10);
        internal1->add_child(leaf11);
        root1->add_child(internal1);
        root1->add_child(leaf12);


        PopulationTree tree0(root0,
                100,   // number of loci
                1,     // length of loci
                true); // validate data
        PopulationTree tree1(root1,
                100,   // number of loci
                1,     // length of loci
                true); // validate data

        tree0.set_all_population_sizes(0.005);
        tree1.set_all_population_sizes(0.005);

        double l0 = tree0.compute_log_likelihood();
        double l1 = tree1.compute_log_likelihood();
        std::cout << "\nPolytomy lnL: " << l0;
        std::cout << "\nBifurcating lnL: " << l1 << "\n\n";
        REQUIRE(l0 == l1);

        RandomNumberGenerator rng = RandomNumberGenerator(1234);
        BiallelicData bd = tree0.simulate_linked_biallelic_data_set(rng,
                1.0,    // singleton sample probability
                false,  // max one variable site per locus
                true);  // validate data set

        tree0.set_data(bd, false);
        tree1.set_data(bd, false);

        l0 = tree0.compute_log_likelihood();
        l1 = tree1.compute_log_likelihood();
        std::cout << "\nPolytomy lnL: " << l0;
        std::cout << "\nBifurcating lnL: " << l1 << "\n\n";
        REQUIRE(l0 == l1);
    }
}

TEST_CASE("Testing likelihood of PopulationTree with four-way polytomy at root", "[PopulationTree]") {

    SECTION("Testing constructor and likelihood calc") {
        std::shared_ptr<PopulationNode> root0 = std::make_shared<PopulationNode>(4, "root", 0.1);
        std::shared_ptr<PopulationNode> leaf00 = std::make_shared<PopulationNode>(0, "leaf 0", 0.0, 4);
        leaf00->fix_node_height();
        std::shared_ptr<PopulationNode> leaf01 = std::make_shared<PopulationNode>(1, "leaf 1", 0.0, 4);
        leaf01->fix_node_height();
        std::shared_ptr<PopulationNode> leaf02 = std::make_shared<PopulationNode>(2, "leaf 2", 0.0, 4);
        leaf02->fix_node_height();
        std::shared_ptr<PopulationNode> leaf03 = std::make_shared<PopulationNode>(3, "leaf 3", 0.0, 4);
        leaf03->fix_node_height();

        root0->add_child(leaf00);
        root0->add_child(leaf01);
        root0->add_child(leaf02);
        root0->add_child(leaf03);

        std::shared_ptr<PopulationNode> internal10 = std::make_shared<PopulationNode>(4, "internal 0", 0.1);
        std::shared_ptr<PopulationNode> internal11 = std::make_shared<PopulationNode>(5, "internal 1", 0.1);
        std::shared_ptr<PopulationNode> root1 = std::make_shared<PopulationNode>(6, "root", 0.1);
        std::shared_ptr<PopulationNode> leaf10 = std::make_shared<PopulationNode>(0, "leaf 0", 0.0, 4);
        leaf10->fix_node_height();
        std::shared_ptr<PopulationNode> leaf11 = std::make_shared<PopulationNode>(1, "leaf 1", 0.0, 4);
        leaf11->fix_node_height();
        std::shared_ptr<PopulationNode> leaf12 = std::make_shared<PopulationNode>(2, "leaf 2", 0.0, 4);
        leaf12->fix_node_height();
        std::shared_ptr<PopulationNode> leaf13 = std::make_shared<PopulationNode>(3, "leaf 3", 0.0, 4);
        leaf13->fix_node_height();

        internal10->add_child(leaf10);
        internal10->add_child(leaf11);
        internal11->add_child(internal10);
        internal11->add_child(leaf12);
        root1->add_child(internal11);
        root1->add_child(leaf13);

        PopulationTree tree0(root0,
                100,   // number of loci
                1,     // length of loci
                true); // validate data
        PopulationTree tree1(root1,
                100,   // number of loci
                1,     // length of loci
                true); // validate data

        tree0.set_all_population_sizes(0.005);
        tree1.set_all_population_sizes(0.005);

        double l0 = tree0.compute_log_likelihood();
        double l1 = tree1.compute_log_likelihood();
        std::cout << "\nPolytomy lnL: " << l0;
        std::cout << "\nBifurcating lnL: " << l1 << "\n\n";
        REQUIRE(l0 == l1);

        RandomNumberGenerator rng = RandomNumberGenerator(1234);
        BiallelicData bd = tree0.simulate_linked_biallelic_data_set(rng,
                1.0,    // singleton sample probability
                false,  // max one variable site per locus
                true);  // validate data set

        tree0.set_data(bd, false);
        tree1.set_data(bd, false);

        l0 = tree0.compute_log_likelihood();
        l1 = tree1.compute_log_likelihood();
        std::cout << "\nPolytomy lnL: " << l0;
        std::cout << "\nBifurcating lnL: " << l1 << "\n\n";
        REQUIRE(l0 == l1);
    }
}

TEST_CASE("Testing likelihood of PopulationTree with three-way polytomy at internal", "[PopulationTree]") {

    SECTION("Testing constructor and likelihood calc") {
        std::shared_ptr<PopulationNode> root0 = std::make_shared<PopulationNode>(5, "root", 0.1);
        std::shared_ptr<PopulationNode> internal00 = std::make_shared<PopulationNode>(4, "internal 0", 0.05);
        std::shared_ptr<PopulationNode> leaf00 = std::make_shared<PopulationNode>(0, "leaf 0", 0.0, 4);
        leaf00->fix_node_height();
        std::shared_ptr<PopulationNode> leaf01 = std::make_shared<PopulationNode>(1, "leaf 1", 0.0, 4);
        leaf01->fix_node_height();
        std::shared_ptr<PopulationNode> leaf02 = std::make_shared<PopulationNode>(2, "leaf 2", 0.0, 4);
        leaf02->fix_node_height();
        std::shared_ptr<PopulationNode> leaf03 = std::make_shared<PopulationNode>(3, "leaf 3", 0.0, 4);
        leaf03->fix_node_height();

        internal00->add_child(leaf00);
        internal00->add_child(leaf01);
        internal00->add_child(leaf02);
        root0->add_child(internal00);
        root0->add_child(leaf03);

        std::shared_ptr<PopulationNode> root1 = std::make_shared<PopulationNode>(6, "root", 0.1);
        std::shared_ptr<PopulationNode> internal10 = std::make_shared<PopulationNode>(4, "internal 0", 0.05);
        std::shared_ptr<PopulationNode> internal11 = std::make_shared<PopulationNode>(5, "internal 1", 0.05);
        std::shared_ptr<PopulationNode> leaf10 = std::make_shared<PopulationNode>(0, "leaf 0", 0.0, 4);
        leaf10->fix_node_height();
        std::shared_ptr<PopulationNode> leaf11 = std::make_shared<PopulationNode>(1, "leaf 1", 0.0, 4);
        leaf11->fix_node_height();
        std::shared_ptr<PopulationNode> leaf12 = std::make_shared<PopulationNode>(2, "leaf 2", 0.0, 4);
        leaf12->fix_node_height();
        std::shared_ptr<PopulationNode> leaf13 = std::make_shared<PopulationNode>(3, "leaf 3", 0.0, 4);
        leaf13->fix_node_height();

        internal10->add_child(leaf10);
        internal10->add_child(leaf11);
        internal11->add_child(internal10);
        internal11->add_child(leaf12);
        root1->add_child(internal11);
        root1->add_child(leaf13);

        PopulationTree tree0(root0,
                100,   // number of loci
                1,     // length of loci
                true); // validate data
        PopulationTree tree1(root1,
                100,   // number of loci
                1,     // length of loci
                true); // validate data

        tree0.set_all_population_sizes(0.005);
        tree1.set_all_population_sizes(0.005);

        double l0 = tree0.compute_log_likelihood();
        double l1 = tree1.compute_log_likelihood();
        std::cout << "\nPolytomy lnL: " << l0;
        std::cout << "\nBifurcating lnL: " << l1 << "\n\n";
        REQUIRE(l0 == l1);

        RandomNumberGenerator rng = RandomNumberGenerator(1234);
        BiallelicData bd = tree0.simulate_linked_biallelic_data_set(rng,
                1.0,    // singleton sample probability
                false,  // max one variable site per locus
                true);  // validate data set

        tree0.set_data(bd, false);
        tree1.set_data(bd, false);

        l0 = tree0.compute_log_likelihood();
        l1 = tree1.compute_log_likelihood();
        std::cout << "\nPolytomy lnL: " << l0;
        std::cout << "\nBifurcating lnL: " << l1 << "\n\n";
        REQUIRE(l0 == l1);
    }
}

TEST_CASE("Testing likelihood of PopulationTree with four-way polytomy at internal", "[PopulationTree]") {

    SECTION("Testing constructor and likelihood calc") {
        std::shared_ptr<PopulationNode> root0 = std::make_shared<PopulationNode>(6, "root", 0.1);
        std::shared_ptr<PopulationNode> internal00 = std::make_shared<PopulationNode>(5, "internal 0", 0.05);
        std::shared_ptr<PopulationNode> leaf00 = std::make_shared<PopulationNode>(0, "leaf 0", 0.0, 4);
        leaf00->fix_node_height();
        std::shared_ptr<PopulationNode> leaf01 = std::make_shared<PopulationNode>(1, "leaf 1", 0.0, 6);
        leaf01->fix_node_height();
        std::shared_ptr<PopulationNode> leaf02 = std::make_shared<PopulationNode>(2, "leaf 2", 0.0, 8);
        leaf02->fix_node_height();
        std::shared_ptr<PopulationNode> leaf03 = std::make_shared<PopulationNode>(3, "leaf 3", 0.0, 10);
        leaf03->fix_node_height();
        std::shared_ptr<PopulationNode> leaf04 = std::make_shared<PopulationNode>(4, "leaf 4", 0.0, 12);
        leaf04->fix_node_height();

        internal00->add_child(leaf00);
        internal00->add_child(leaf01);
        internal00->add_child(leaf02);
        internal00->add_child(leaf03);
        root0->add_child(internal00);
        root0->add_child(leaf04);

        std::shared_ptr<PopulationNode> root1 = std::make_shared<PopulationNode>(8, "root", 0.1);
        std::shared_ptr<PopulationNode> internal10 = std::make_shared<PopulationNode>(5, "internal 0", 0.05);
        std::shared_ptr<PopulationNode> internal11 = std::make_shared<PopulationNode>(6, "internal 1", 0.05);
        std::shared_ptr<PopulationNode> internal12 = std::make_shared<PopulationNode>(7, "internal 2", 0.05);
        std::shared_ptr<PopulationNode> leaf10 = std::make_shared<PopulationNode>(0, "leaf 0", 0.0, 4);
        leaf10->fix_node_height();
        std::shared_ptr<PopulationNode> leaf11 = std::make_shared<PopulationNode>(1, "leaf 1", 0.0, 6);
        leaf11->fix_node_height();
        std::shared_ptr<PopulationNode> leaf12 = std::make_shared<PopulationNode>(2, "leaf 2", 0.0, 8);
        leaf12->fix_node_height();
        std::shared_ptr<PopulationNode> leaf13 = std::make_shared<PopulationNode>(3, "leaf 3", 0.0, 10);
        leaf13->fix_node_height();
        std::shared_ptr<PopulationNode> leaf14 = std::make_shared<PopulationNode>(4, "leaf 4", 0.0, 12);
        leaf14->fix_node_height();

        internal10->add_child(leaf10);
        internal10->add_child(leaf11);
        internal11->add_child(internal10);
        internal11->add_child(leaf12);
        internal12->add_child(internal11);
        internal12->add_child(leaf13);
        root1->add_child(internal12);
        root1->add_child(leaf14);

        PopulationTree tree0(root0,
                100,   // number of loci
                1,     // length of loci
                true); // validate data
        PopulationTree tree1(root1,
                100,   // number of loci
                1,     // length of loci
                true); // validate data

        tree0.set_all_population_sizes(0.005);
        tree1.set_all_population_sizes(0.005);

        double l0 = tree0.compute_log_likelihood();
        double l1 = tree1.compute_log_likelihood();
        std::cout << "\nPolytomy lnL: " << l0;
        std::cout << "\nBifurcating lnL: " << l1 << "\n\n";
        REQUIRE(l0 == l1);

        RandomNumberGenerator rng = RandomNumberGenerator(1234);
        BiallelicData bd = tree0.simulate_linked_biallelic_data_set(rng,
                1.0,    // singleton sample probability
                false,  // max one variable site per locus
                true);  // validate data set

        tree0.set_data(bd, false);
        tree1.set_data(bd, false);

        l0 = tree0.compute_log_likelihood();
        l1 = tree1.compute_log_likelihood();
        std::cout << "\nPolytomy lnL: " << l0;
        std::cout << "\nBifurcating lnL: " << l1 << "\n\n";
        REQUIRE(l0 == l1);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood = -31.77866581319647
TEST_CASE("Testing simple likelihood of PopulationTree", "[PopulationTree]") {

    SECTION("Testing constructor and likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        PopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing threaded likelihood of PopulationTree", "[PopulationTree]") {

    SECTION("Testing constructor and threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        PopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing over-threaded likelihood of PopulationTree", "[PopulationTree]") {

    SECTION("Testing constructor and over-threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        PopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(10);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(10);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood             = -248.93254688526213
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (0.01, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded likelihood (0.01, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7099.716015109998
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.01, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.01, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -6986.120524781545
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.01, 10.0, 1.0, 1.0, dominant)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6986.120524781545));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.01, 10.0, 1.0, 1.0, dominant)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(5);
        REQUIRE(l == Approx(-6986.120524781545));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}



// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -328.39238828878365
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (0.0, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded likelihood (0.0, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(7);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(7);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7256.501742344454
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.0, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.0, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7223.362711937651
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.0, 10.0, 1.0, 1.0, dominant)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7223.362711937651));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.0, 10.0, 1.0, 1.0, dominant)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7223.362711937651));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -227.41048391087554
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (0.2, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded likelihood (0.2, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7304.180743441677
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.2, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.2, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7405.145951634711
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.2, 10.0, 1.0, 1.0, dominant)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7405.145951634711));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.2, 10.0, 1.0, 1.0, dominant)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(6);
        REQUIRE(l == Approx(-7405.145951634711));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}




// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -327.7437811413033
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -6472.856486972301
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(5);
        REQUIRE(l == Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(3);
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -6494.774924871097
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.03, 10.0, 10.0, 10.0/19.0, dominant)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6494.774924871097));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.03, 10.0, 10.0, 10.0/19.0, dominant)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-6494.774924871097));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
  
  
// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -265.0023534261969
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(1);
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -10163.468886613919
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(5);
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -10999.288193543642
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.03, 10.0, 10.0/19.0, 10.0, dominant)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-10999.288193543642));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.03, 10.0, 10.0/19.0, 10.0, dominant)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-10999.288193543642));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}



// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -224.40177558289847
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-224.40177558289847));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(100000);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-224.40177558289847));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -8158.88094671241
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-8158.88094671241));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(1000000);
        REQUIRE(l == Approx(-8158.88094671241));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -8034.250341980543
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.03, 111.1, 10.0/19.0, 10.0, dominant)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-8034.250341980543));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.03, 111.1, 10.0/19.0, 10.0, dominant)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-8034.250341980543));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood = -31.77866581319647
TEST_CASE("Testing simple likelihood of ComparisonPopulationTree", "[ComparisonPopulationTree]") {

    SECTION("Testing constructor and likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_root().get_label() == "root-pop1");
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing simple threaded likelihood of ComparisonPopulationTree", "[ComparisonPopulationTree]") {

    SECTION("Testing constructor and threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_root().get_label() == "root-pop1");
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

TEST_CASE("Testing simple prior of PopulationTree", "[PopulationTree]") {

    SECTION("Testing constructor and prior calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        PopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);

        tree.set_root_height(0.1);
        tree.set_all_population_sizes(2.0/100.0);
        tree.set_freq_1(0.95);

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(100.0));
        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));
        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(2.0, 1.0));

        // height   -5.3948298140119091
        // sizes     3 * -155.90663080917298
        // f1       0.64185388617239469 
        // total    -472.47286835535851
        
        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-467.71989242751897));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == Approx(0.64185388617239469));

        double d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-472.47286835535851));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-472.47286835535851));


        tree.store_prior_density();
        tree.constrain_state_frequencies();

        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-467.71989242751897));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == 0.0);

        d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-473.1147222415309));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-473.1147222415309));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-472.47286835535851));


        tree.store_prior_density();
        tree.constrain_population_sizes();

        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-155.90663080917298));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == 0.0);

        d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-161.30146062318488));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-161.30146062318488));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-473.1147222415309));

        tree.restore_prior_density();
        REQUIRE(tree.get_log_prior_density_value() == Approx(-473.1147222415309));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-473.1147222415309));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}


TEST_CASE("Testing hemi129.nex state manipulation", "[ComparisonPopulationTree]") {

    SECTION("Testing state manipulation") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(10.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(1.5, 3.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_height(0.01);
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_child_population_size(0, 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_child_population_size(1, 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1(0.5);
        REQUIRE(tree.is_dirty());

        REQUIRE(tree.get_number_of_likelihood_calculations() == 0);

        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 1);

        tree.store_state();
        tree.set_root_height(0.0);
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 2);

        tree.restore_state();
        // ComparisonPopulationTree does not store/restore heights
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.compute_log_likelihood_and_prior();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.store_state();
        tree.set_root_height(0.0);
        tree.compute_log_likelihood_and_prior();
        REQUIRE(tree.get_number_of_likelihood_calculations() == 4);


        tree.store_state();
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 5);


        tree.store_state();
        tree.set_root_height(0.03);
        tree.set_freq_1(0.05);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.05);
        REQUIRE(tree.get_freq_0() == 0.95);
        REQUIRE(tree.get_u() == Approx(10.0));
        REQUIRE(tree.get_v() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.6991237645509));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 6);


        tree.store_state();
        tree.set_freq_1(0.95);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_freq_1() == 0.95);
        REQUIRE(tree.get_freq_0() == Approx(0.05));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1347.1157822333005));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.6991237645509));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 7);


        tree.store_state();
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-47.141114882027253));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1347.1157822333005));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 8);

        tree.store_state();
        tree.constrain_state_frequencies();
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.9800426669165));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-47.141114882027253));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 9);


        tree.store_state();
        tree.fold_patterns();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.9800426669165));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.9800426669165));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 10);


        tree.store_state();
        tree.constrain_population_sizes();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.66001422230551));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.9800426669165));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 11);

        tree.store_state();
        tree.estimate_mutation_rate();
        REQUIRE(! tree.mutation_rate_is_fixed());
        tree.set_mutation_rate(1.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());
        tree.set_mutation_rate_prior(std::make_shared<GammaDistribution>(10.0, 0.1));
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.66001422230551));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 12);

        tree.store_state();
        tree.set_mutation_rate(2.0);
        tree.set_root_height(0.1);
        tree.set_root_population_size(2.0/(20.0 * 2 * tree.get_ploidy()));
        tree.set_child_population_size(0, 2.0/(20.0 * 2 * tree.get_ploidy()));
        tree.set_child_population_size(1, 2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 2.0);
        REQUIRE(tree.get_root_population_size() == 2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-207.43599077244659));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        tree.restore_state();
        // ComparisonPopulationTree does not store/restore node heights
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

TEST_CASE("Testing hemi129.nex state manipulation for PopulationTree", "[PopulationTree]") {

    SECTION("Testing state manipulation") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(10.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(1.5, 3.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_height(0.01);
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1(0.5);
        REQUIRE(tree.is_dirty());

        REQUIRE(tree.get_number_of_likelihood_calculations() == 0);

        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1340.6289538974045));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 1);

        tree.store_state();
        tree.set_root_height(0.0);
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1340.5289538974043));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1340.6289538974045));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 2);

        tree.restore_state();
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1340.6289538974045));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1340.6289538974045));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.compute_log_likelihood_and_prior();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1340.6289538974045));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1340.6289538974045));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.store_state();
        tree.set_root_height(0.0);
        tree.compute_log_likelihood_and_prior();
        REQUIRE(tree.get_number_of_likelihood_calculations() == 4);


        tree.store_state();
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.5289538974043));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1340.5289538974043));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 5);


        tree.store_state();
        tree.set_root_height(0.03);
        tree.set_freq_1(0.05);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.05);
        REQUIRE(tree.get_freq_0() == Approx(0.95));
        REQUIRE(tree.get_u() == Approx(10.0));
        REQUIRE(tree.get_v() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1340.6965386715569));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.5289538974043));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 6);


        tree.store_state();
        tree.set_freq_1(0.95);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.95);
        REQUIRE(tree.get_freq_0() == Approx(0.05));
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1345.1131971403065));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1340.6965386715569));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 7);


        tree.store_state();
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-45.138529789033207));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1345.1131971403065));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 8);

        tree.store_state();
        tree.constrain_state_frequencies();
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.6774575739223));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-45.138529789033207));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 9);


        tree.store_state();
        tree.fold_patterns();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.6774575739223));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.6774575739223));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 10);


        tree.store_state();
        tree.constrain_population_sizes();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.35742912931147));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.6774575739223));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 11);

        tree.store_state();
        tree.estimate_mutation_rate();
        REQUIRE(! tree.mutation_rate_is_fixed());
        tree.set_mutation_rate(1.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());
        tree.set_mutation_rate_prior(std::make_shared<GammaDistribution>(10.0, 0.1));
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.35742912931147));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 12);

        tree.store_state();
        tree.set_mutation_rate(2.0);
        tree.set_root_height(0.1);
        tree.set_all_population_sizes(2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 2.0);
        REQUIRE(tree.get_root_population_size() == 2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-206.13340567945255));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        tree.restore_state();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 9fb5b3b7817a0bd4a21e3f90132f132cca72ce4e)
// SNAPP v1.3.0 (master 4f3f0f7366798f4fb38b766c15f6426a75ddf71e)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0/19.0
// Log likelihood            = -23.81984255023975
// Log likelihood correction = -6.87935580446044
//
// With constant sites inclucded and m_bUseNonPolymorphic = true
// Log likelihood            = -55.01646493341547
// Log likelihood correction = -6.87935580446044
TEST_CASE("Testing affect of constant sites on likelihood of PopulationTree", "[PopulationTree]") {

    SECTION("Testing haploid-standard-full-constant.nex") {
        std::string nex_path = "data/haploid-standard-full-constant.nex";
        PopulationTree t_included(
                nex_path, // path
                ' ',      // pop name delimiter
                true,     // pop name is prefix
                false,    // genotypes are diploid
                false,    // markers are dominant
                false,    // constant sites removed
                true);    // validate
        std::string nex_path2 = "data/haploid-standard-full-constant-removed.nex";
        PopulationTree t(
                nex_path2, // path
                ' ',       // pop name delimiter
                true,      // pop name is prefix
                false,     // genotypes are diploid
                false,     // markers are dominant
                true,      // constant sites removed
                true);     // validate

        t.set_root_height(0.03);
        t.set_freq_1(0.05);
        t.set_all_population_sizes(2.0/(10.0 * 2 * t.get_ploidy()));

        t_included.set_root_height(0.03);
        t_included.set_freq_1(0.05);
        t_included.set_all_population_sizes(2.0/(10.0 * 2 * t_included.get_ploidy()));

        double l = t.compute_log_likelihood();
        REQUIRE(l == Approx(-23.81984255023975));
        REQUIRE(t.get_likelihood_correction() == Approx(-6.87935580446044));

        double l_included = t_included.compute_log_likelihood();

        REQUIRE(l_included == Approx(-55.01646493341547));

        REQUIRE(t.get_likelihood_correction() == t_included.get_likelihood_correction());
    }
}

TEST_CASE("Testing affect of constant sites on threaded likelihood of PopulationTree", "[PopulationTree]") {

    SECTION("Testing haploid-standard-full-constant.nex") {
        std::string nex_path = "data/haploid-standard-full-constant.nex";
        PopulationTree t_included(
                nex_path, // path
                ' ',      // pop name delimiter
                true,     // pop name is prefix
                false,    // genotypes are diploid
                false,    // markers are dominant
                false,    // constant sites removed
                true);    // validate
        std::string nex_path2 = "data/haploid-standard-full-constant-removed.nex";
        PopulationTree t(
                nex_path2, // path
                ' ',       // pop name delimiter
                true,      // pop name is prefix
                false,     // genotypes are diploid
                false,     // markers are dominant
                true,      // constant sites removed
                true);     // validate

        t.set_root_height(0.03);
        t.set_freq_1(0.05);
        t.set_all_population_sizes(2.0/(10.0 * 2 * t.get_ploidy()));

        t_included.set_root_height(0.03);
        t_included.set_freq_1(0.05);
        t_included.set_all_population_sizes(2.0/(10.0 * 2 * t_included.get_ploidy()));

        double l = t.compute_log_likelihood(4);
        REQUIRE(l == Approx(-23.81984255023975));
        REQUIRE(t.get_likelihood_correction() == Approx(-6.87935580446044));

        double l_included = t_included.compute_log_likelihood(2);

        REQUIRE(l_included == Approx(-55.01646493341547));

        REQUIRE(t.get_likelihood_correction() == t_included.get_likelihood_correction());
    }
}


// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.01
// coalescent_rate = 100.0, 200.0, 500.0
// u = 1.0
// v = 1.0
// Log likelihood            = -226.11914854623677
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (0.01, 100.0, 200.0, 500.0)", "[ComparisonPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_child_population_size(0, 2.0/(100.0 * 2 * tree.get_ploidy()));
        tree.set_child_population_size(1, 2.0/(200.0 * 2 * tree.get_ploidy()));
        tree.set_root_population_size(2.0/(500.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-226.11914854623677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}


// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.00506843962151613554
// coalescent_rate = 2.0 / 0.00018955324120485613
// u = 1.0
// v = 1.0
// Log likelihood            = -277.06960543551577
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (0.00506843962151613554, 2.0 / 0.00018955324120485613)", "[ComparisonPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.00506843962151613554);
        tree.set_all_population_sizes(0.00018955324120485613 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-277.06960543551577));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 9.08323190033687971e-09
// coalescent_rate = 2.0 / 2.47975039926886321e-08
// u = 1.0
// v = 1.0
// Log likelihood            = -221.69627648370943
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (9.08323190033687971e-09, 2.0 / 2.47975039926886321e-08)", "[ComparisonPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(9.08323190033687971e-09);
        tree.set_all_population_sizes(2.47975039926886321e-08 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-221.69627648370943));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.04921319733994759e-08
// coalescent_rate = 2.0 / 2.75977168733651178e-10
// u = 1.0
// v = 1.0
// Log likelihood            = -324.2737564069293
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (1.04921319733994759e-08, 2.0 / 2.75977168733651178e-10)", "[ComparisonPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.04921319733994759e-08);
        tree.set_all_population_sizes(2.75977168733651178e-10 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-324.2737564069293));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.012386610001351e-08
// coalescent_rate = 2.0 / 5.75048645855884647e-30
// u = 1.0
// v = 1.0
// Log likelihood            = -1364.1427530000253
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (1.012386610001351e-08, 2.0 / 5.75048645855884647e-30)", "[ComparisonPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.012386610001351e-08);
        tree.set_all_population_sizes(5.75048645855884647e-30 / 4.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-1364.1427530000253));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.04856228318474786e-08
// coalescent_rate = 2.0 / 4.43934332792563837e-305
// u = 1.0
// v = 1.0
// Log likelihood            = NaN
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (1.04856228318474786e-08, 2.0 / 4.43934332792563837e-305)", "[ComparisonPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.04856228318474786e-08);
        tree.set_all_population_sizes(4.43934332792563837e-305 / 4.0);
        double l = tree.compute_log_likelihood();
        // REQUIRE(std::isnan(l));
        REQUIRE(std::isinf(l));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.012386610001351e-08
// coalescent_rate = 2.0 / 5.46641122085615013e-09,
//                   2.0 / 7.39871781998828579e-08,
//                   2.0 / 2.71077053326002069e-13
// u = 1.0
// v = 1.0
// Log likelihood            = -96.34394008351177
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing hemi129.nex weirdness", "[ComparisonPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.012386610001351e-08);
        tree.set_child_population_size(0, 5.46641122085615013e-09 / 4.0);
        tree.set_child_population_size(1, 7.39871781998828579e-08 / 4.0);
        tree.set_root_population_size(2.71077053326002069e-13 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-96.34394008351177));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.036374107244057e-08
// coalescent_rate = 2.0 / 4.57999694763258361e-09,
//                   2.0 / 6.70991782555376588e-08,
//                   2.0 / 1.33514111020266258e-08
// u = 1.0
// v = 1.0
// Log likelihood            = -44.95791900747736
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing hemi129.nex weirdness 2", "[ComparisonPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.036374107244057e-08);
        tree.set_child_population_size(0, 4.57999694763258361e-09 / 4.0);
        tree.set_child_population_size(1, 6.70991782555376588e-08 / 4.0);
        tree.set_root_population_size(1.33514111020266258e-08 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-44.95791900747736));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

TEST_CASE("Testing coalesce_in_branch for 2 lineages and theta of 1.0",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing coalesce_in_branch for 2 lineages and theta of 3.7",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 3.7;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing coalesce_in_branch for 2 lineages and theta of 0.17",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing coalesce_in_branch for 3 lineages and theta of 1.0",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing coalesce_in_branch for 3 lineages and theta of 1.47",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 1.47;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing coalesce_in_branch for 3 lineages and theta of 0.17",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing coalesce_in_branch for 10 lineages and theta of 1.0",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing coalesce_in_branch for 10 lineages and theta of 1.47",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 1.47;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing coalesce_in_branch for 10 lineages and theta of 0.17",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing scaling of simulate_gene_tree for singleton",
        "[ComparisonPopulationTree]") {

    SECTION("Testing singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        std::string nex_path = "data/singleton-n2.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, false, false);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.set_child_population_size(0, (Ne * mu));

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        unsigned int n;
        double mean;
        double variance;
        double sum_devs;
        double d;
        double d_n;
        double mn;
        double mx;
        double x;
        std::shared_ptr<GeneTreeSimNode> gtree;

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            x = gtree->get_height();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);

        REQUIRE(gtree->is_root());
        REQUIRE(gtree->get_number_of_children() == 2);
        REQUIRE(gtree->get_leaf_node_count() == 2);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.0005));
        REQUIRE(mn >= 0.0);

        tree.set_mutation_rate(mu);
        tree.set_root_population_size(Ne);
        tree.set_child_population_size(0, Ne);
        tree.set_root_height(0.1 / mu);

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            x = gtree->get_height();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);

        REQUIRE(gtree->is_root());
        REQUIRE(gtree->get_number_of_children() == 2);
        REQUIRE(gtree->get_leaf_node_count() == 2);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing scaling of simulate_gene_tree for pair",
        "[ComparisonPopulationTree]") {

    SECTION("Testing pair") {
        unsigned int Ne_root = 100000;
        unsigned int Ne_0 = 200000;
        unsigned int Ne_1 = 10000;
        double mu = 1e-8;
        double theta_root = 4 * Ne_root * mu;
        double theta_0 = 4 * Ne_0 * mu;
        double theta_1 = 4 * Ne_1 * mu;
        double time = 2e7;

        double expected_mean_root = theta_root * (1.0 - (1.0 / 2));
        double expected_variance_root = expected_mean_root * expected_mean_root;
        expected_mean_root += (time * mu);

        double expected_mean_0 = theta_0 * (1.0 - (1.0 / 10));
        double expected_mean_1 = theta_1 * (1.0 - (1.0 / 10));

        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne_root);
        tree.set_child_population_size(0, (Ne_0));
        tree.set_child_population_size(1, (Ne_1));

        tree.set_root_height(time);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(mu);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        std::shared_ptr<GeneTreeSimNode> gtree;
        SampleSummarizer<double> height_root;
        SampleSummarizer<double> height_0;
        SampleSummarizer<double> height_1;
        SampleSummarizer<double> tree_length;

        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            REQUIRE(gtree->is_root());
            REQUIRE(gtree->get_number_of_children() == 2);
            REQUIRE(gtree->get_leaf_node_count() == 20);
            height_root.add_sample(gtree->get_height());
            double h0 = gtree->get_child(0)->get_height();
            double h1 = gtree->get_child(1)->get_height();
            if (h0 < h1) {
                height_0.add_sample(h1);
                height_1.add_sample(h0);
            }
            else {
                height_0.add_sample(h0);
                height_1.add_sample(h1);
            }
            tree_length.add_sample(gtree->get_clade_length());
        }

        REQUIRE(height_root.mean() == Approx(expected_mean_root).epsilon(0.0005));
        REQUIRE(height_root.variance() == Approx(expected_variance_root).epsilon(0.0005));

        REQUIRE(height_0.mean() == Approx(expected_mean_0).epsilon(0.00003));
        REQUIRE(height_1.mean() == Approx(expected_mean_1).epsilon(0.00003));


        tree.set_mutation_rate(mu * 2.0);

        SampleSummarizer<double> scale_height_root;
        SampleSummarizer<double> scale_height_0;
        SampleSummarizer<double> scale_height_1;
        SampleSummarizer<double> scale_tree_length;

        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            gtree->scale(0.5);
            REQUIRE(gtree->is_root());
            REQUIRE(gtree->get_number_of_children() == 2);
            REQUIRE(gtree->get_leaf_node_count() == 20);
            scale_height_root.add_sample(gtree->get_height());
            double h0 = gtree->get_child(0)->get_height();
            double h1 = gtree->get_child(1)->get_height();
            if (h0 < h1) {
                scale_height_0.add_sample(h1);
                scale_height_1.add_sample(h0);
            }
            else {
                scale_height_0.add_sample(h0);
                scale_height_1.add_sample(h1);
            }
            scale_tree_length.add_sample(gtree->get_clade_length());
        }

        REQUIRE(scale_height_root.mean() == Approx(expected_mean_root).epsilon(0.0005));
        REQUIRE(scale_height_root.variance() == Approx(expected_variance_root).epsilon(0.0005));

        REQUIRE(scale_height_0.mean() == Approx(expected_mean_0).epsilon(0.00003));
        REQUIRE(scale_height_1.mean() == Approx(expected_mean_1).epsilon(0.00003));

        REQUIRE(scale_tree_length.mean() == Approx(tree_length.mean()).epsilon(0.0001));
        REQUIRE(scale_tree_length.variance() == Approx(tree_length.variance()).epsilon(0.0001));
    }
}

TEST_CASE("Testing dataset simulation", "[ComparisonPopulationTree]") {
    SECTION("Testing simulate_biallelic_data_set for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        ComparisonPopulationTree tree(nex_path, ' ', true, false, false, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        BiallelicData data = tree.simulate_biallelic_data_set(rng);

        REQUIRE(data.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data.markers_are_dominant() == false);
        REQUIRE(tree.get_data().markers_are_dominant() == false);

        data.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8).epsilon(0.01));
        REQUIRE(data.get_proportion_1() == Approx(0.625).epsilon(0.01));

        std::string io_nex_path = "data/tmp-data-test1.nex";
        std::ofstream out;
        out.open(io_nex_path);
        data.write_nexus(out, '-');
        out.close();
        REQUIRE(path::exists(io_nex_path));
        BiallelicData io_data(io_nex_path,
                '-',   // pop name delimiter
                true,  // pop name is prefix
                false, // genotypes are diploid
                false, // markers are dominant
                true,  // validate
                false   // store seq loci info
                );
        REQUIRE(io_data.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(io_data.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(io_data.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(io_data.markers_are_dominant() == false);

        io_data.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8).epsilon(0.01));
        REQUIRE(io_data.get_proportion_1() == Approx(0.625).epsilon(0.01));
    }
}

TEST_CASE("Testing complete linked dataset simulation", "[ComparisonPopulationTree]") {
    SECTION("Testing simulate_complete_biallelic_data_set for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        ComparisonPopulationTree tree(nex_path, ' ', true, false, false, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        unsigned int locus_size = 100;
        std::pair<BiallelicData, unsigned int> data_nloci= tree.simulate_complete_biallelic_data_set(
                rng,
                locus_size,
                1.0,
                true);
        BiallelicData data = data_nloci.first;
        unsigned int nloci = data_nloci.second;
        std::vector<unsigned int> expected_locus_end_indices;
        unsigned int end_idx = (locus_size - 1);
        while(end_idx < tree.get_data().get_number_of_sites()) {
            expected_locus_end_indices.push_back(end_idx);
            end_idx += locus_size;
        }
        expected_locus_end_indices.push_back(tree.get_data().get_number_of_sites() - 1);
        
        REQUIRE(data.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data.markers_are_dominant() == false);
        REQUIRE(tree.get_data().markers_are_dominant() == false);
        REQUIRE(data.get_locus_end_indices() == expected_locus_end_indices);
        REQUIRE(nloci == data.get_locus_end_indices().size());

        data.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8).epsilon(0.01));
        REQUIRE(data.get_proportion_1() == Approx(0.625).epsilon(0.01));

        std::string io_nex_path = "data/tmp-data-test-complete1.nex";
        std::ofstream out;
        out.open(io_nex_path);
        data.write_nexus(out, '-');
        out.close();
        REQUIRE(path::exists(io_nex_path));
        BiallelicData io_data(io_nex_path,
                '-',   // pop name delimiter
                true,  // pop name is prefix
                false, // genotypes are diploid
                false, // markers are dominant
                true,  // validate
                true   // store seq loci info
                );
        REQUIRE(io_data.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(io_data.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(io_data.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(io_data.markers_are_dominant() == false);
        REQUIRE(io_data.get_locus_end_indices() == expected_locus_end_indices);

        io_data.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8).epsilon(0.01));
        REQUIRE(io_data.get_proportion_1() == Approx(0.625).epsilon(0.01));
    }
}

TEST_CASE("Testing complete linked dataset simulation one locus", "[ComparisonPopulationTree]") {
    SECTION("Testing simulate_complete_biallelic_data_set for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        ComparisonPopulationTree tree(nex_path, ' ', true, false, false, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        unsigned int locus_size = 10000;
        std::pair<BiallelicData, unsigned int> data_nloci= tree.simulate_complete_biallelic_data_set(
                rng,
                locus_size,
                1.0,
                true);
        BiallelicData data = data_nloci.first;
        unsigned int nloci = data_nloci.second;
        REQUIRE(nloci == 1);
        std::vector<unsigned int> expected_locus_end_indices;
        expected_locus_end_indices.push_back(tree.get_data().get_number_of_sites() - 1);
        
        REQUIRE(data.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data.markers_are_dominant() == false);
        REQUIRE(tree.get_data().markers_are_dominant() == false);
        REQUIRE(data.get_locus_end_indices() == expected_locus_end_indices);

        data.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8).epsilon(0.01));
        REQUIRE(data.get_proportion_1() == Approx(0.625).epsilon(0.01));

        std::string io_nex_path = "data/tmp-data-test-complete1.nex";
        std::ofstream out;
        out.open(io_nex_path);
        data.write_nexus(out, '-');
        out.close();
        REQUIRE(path::exists(io_nex_path));
        BiallelicData io_data(io_nex_path,
                '-',   // pop name delimiter
                true,  // pop name is prefix
                false, // genotypes are diploid
                false, // markers are dominant
                true,  // validate
                true   // store seq loci info
                );
        REQUIRE(io_data.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(io_data.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(io_data.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(io_data.markers_are_dominant() == false);
        REQUIRE(io_data.get_locus_end_indices() == expected_locus_end_indices);

        io_data.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8).epsilon(0.01));
        REQUIRE(io_data.get_proportion_1() == Approx(0.625).epsilon(0.01));
    }
}

TEST_CASE("Testing linked dataset simulation", "[ComparisonPopulationTree]") {
    SECTION("Testing simulate_linked_biallelic_data_set for fully fixed pair") {
        std::string nex_path = "data/hemi129-with-missing.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        ComparisonPopulationTree tree(nex_path, ' ',
                true,   // population_name_is_prefix
                true,   // diploid
                false,  // dominant
                false,  // constant sites removed
                true,   // validate
                false,  // strict on constant
                false,  // strict on missing
                false,  // strict on triallelic
                2.0,    // ploidy
                true    // store charset info
                );

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        double u_sum = 0.0;
        double prop_sum = 0.0;
        double io_u_sum = 0.0;
        double io_prop_sum = 0.0;
        unsigned int nreps = 100;
        BiallelicData data;
        BiallelicData io_data;
        for (unsigned int i = 0; i < nreps; ++i) {
            data = tree.simulate_linked_biallelic_data_set(rng,
                    1.0,    // singleton sample prob
                    false,  // max one var site per locus
                    true    // validate
                    );

            REQUIRE(data.get_number_of_sites() == tree.get_data().get_number_of_sites());
            REQUIRE(data.get_number_of_populations() == tree.get_data().get_number_of_populations());
            REQUIRE(data.get_population_labels() == tree.get_data().get_population_labels());
            REQUIRE(data.markers_are_dominant() == false);
            REQUIRE(tree.get_data().markers_are_dominant() == false);

            REQUIRE(tree.get_data().get_locus_end_indices() == data.get_locus_end_indices());
            REQUIRE(tree.get_data().has_seq_loci_info() == data.has_seq_loci_info());

            data.get_empirical_u_v_rates(u, v);
            u_sum += u;
            prop_sum += data.get_proportion_1();

            std::string io_nex_path = "data/tmp-data-test2-" + std::to_string(i) + ".nex";
            std::ofstream out;
            out.open(io_nex_path);
            data.write_nexus(out, '-');
            out.close();
            REQUIRE(path::exists(io_nex_path));
            io_data = BiallelicData(io_nex_path,
                    '-',   // pop name delimiter
                    true,  // pop name is prefix
                    false, // genotypes are diploid
                    false, // markers are dominant
                    true,  // validate
                    true   // store seq loci info
                    );
            REQUIRE(io_data.get_number_of_sites() == tree.get_data().get_number_of_sites());
            REQUIRE(io_data.get_number_of_populations() == tree.get_data().get_number_of_populations());
            REQUIRE(io_data.get_population_labels() == tree.get_data().get_population_labels());
            REQUIRE(io_data.markers_are_dominant() == false);

            REQUIRE(tree.get_data().get_locus_end_indices() == io_data.get_locus_end_indices());
            REQUIRE(tree.get_data().has_seq_loci_info() == io_data.has_seq_loci_info());

            io_data.get_empirical_u_v_rates(u, v);
            io_u_sum += u;
            io_prop_sum += data.get_proportion_1();
        }
        REQUIRE(u_sum / nreps == Approx(0.8).epsilon(0.01));
        REQUIRE(prop_sum / nreps == Approx(0.625).epsilon(0.01));
        REQUIRE(io_u_sum == Approx(u_sum));
        REQUIRE(io_prop_sum == Approx(prop_sum));

        REQUIRE(data.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {3, 8, 13, 18};
        REQUIRE(data.get_contiguous_pattern_indices().size() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_locus_end_indices() == expected_locus_ends);

        REQUIRE(io_data.has_seq_loci_info() == true);
        REQUIRE(io_data.get_contiguous_pattern_indices().size() == tree.get_data().get_number_of_sites());
        REQUIRE(io_data.get_locus_end_indices() == expected_locus_ends);
    }
}

TEST_CASE("Testing linked dataset simulation with aflp dataset", "[ComparisonPopulationTree]") {
    SECTION("Testing simulate_linked_biallelic_data_set for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        ComparisonPopulationTree tree(nex_path, ' ',
                true,   // population_name_is_prefix
                false,  // diploid
                false,  // dominant
                false,  // constant sites removed
                true,   // validate
                false,  // strict on constant
                false,  // strict on missing
                false,  // strict on triallelic
                2.0,    // ploidy
                true    // store charset info
                );

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        double u_sum = 0.0;
        double prop_sum = 0.0;
        double io_u_sum = 0.0;
        double io_prop_sum = 0.0;
        unsigned int nreps = 10;
        BiallelicData data;
        BiallelicData io_data;
        for (unsigned int i = 0; i < nreps; ++i) {
            data = tree.simulate_linked_biallelic_data_set(rng,
                    1.0,    // singleton sample prob
                    false,  // max one var site per locus
                    true    // validate
                    );

            REQUIRE(data.get_number_of_sites() == tree.get_data().get_number_of_sites());
            REQUIRE(data.get_number_of_populations() == tree.get_data().get_number_of_populations());
            REQUIRE(data.get_population_labels() == tree.get_data().get_population_labels());
            REQUIRE(data.markers_are_dominant() == false);
            REQUIRE(tree.get_data().markers_are_dominant() == false);

            REQUIRE(tree.get_data().get_locus_end_indices() == data.get_locus_end_indices());
            REQUIRE(tree.get_data().has_seq_loci_info() == data.has_seq_loci_info());

            data.get_empirical_u_v_rates(u, v);
            u_sum += u;
            prop_sum += data.get_proportion_1();

            std::string io_nex_path = "data/tmp-data-test3-" + std::to_string(i) + ".nex";
            std::ofstream out;
            out.open(io_nex_path);
            data.write_nexus(out, '-');
            out.close();
            REQUIRE(path::exists(io_nex_path));
            io_data = BiallelicData(io_nex_path,
                    '-',   // pop name delimiter
                    true,  // pop name is prefix
                    false, // genotypes are diploid
                    false, // markers are dominant
                    true,  // validate
                    true   // store seq loci info
                    );
            REQUIRE(io_data.get_number_of_sites() == tree.get_data().get_number_of_sites());
            REQUIRE(io_data.get_number_of_populations() == tree.get_data().get_number_of_populations());
            REQUIRE(io_data.get_population_labels() == tree.get_data().get_population_labels());
            REQUIRE(io_data.markers_are_dominant() == false);

            REQUIRE(tree.get_data().get_locus_end_indices() == io_data.get_locus_end_indices());
            REQUIRE(tree.get_data().has_seq_loci_info() == io_data.has_seq_loci_info());

            io_data.get_empirical_u_v_rates(u, v);
            io_u_sum += u;
            io_prop_sum += data.get_proportion_1();
        }
        REQUIRE(u_sum / nreps == Approx(0.8).epsilon(0.01));
        REQUIRE(prop_sum / nreps == Approx(0.625).epsilon(0.01));
        REQUIRE(io_u_sum == Approx(u_sum));
        REQUIRE(io_prop_sum == Approx(prop_sum));

        REQUIRE(data.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {99, 199, 299, 399, 499, 599, 699, 799, 899, 999, 1099, 1199, 1216};
        REQUIRE(data.get_contiguous_pattern_indices().size() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_locus_end_indices() == expected_locus_ends);

        REQUIRE(io_data.has_seq_loci_info() == true);
        REQUIRE(io_data.get_contiguous_pattern_indices().size() == tree.get_data().get_number_of_sites());
        REQUIRE(io_data.get_locus_end_indices() == expected_locus_ends);
    }
}

TEST_CASE("Testing singleton acquisition bias", "[ComparisonPopulationTree]") {
    SECTION("Testing for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        ComparisonPopulationTree tree(nex_path, ' ', true, false, false, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        BiallelicData data1 = tree.simulate_biallelic_data_set(rng, 1.0);
        BiallelicData data0 = tree.simulate_biallelic_data_set(rng, 0.0);

        REQUIRE(data1.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data1.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data1.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data1.markers_are_dominant() == false);
        REQUIRE(data0.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data0.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data0.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data0.markers_are_dominant() == false);
        REQUIRE(tree.get_data().markers_are_dominant() == false);

        unsigned int singleton_count10 = 0;
        unsigned int singleton_count05 = 0;
        unsigned int singleton_count00 = 0;
        RandomNumberGenerator rng10 = RandomNumberGenerator(123);
        RandomNumberGenerator rng05 = RandomNumberGenerator(123);
        RandomNumberGenerator rng00 = RandomNumberGenerator(123);
        for (unsigned int rep = 0; rep < 100; ++rep) {
            BiallelicData data00 = tree.simulate_biallelic_data_set(rng00, 0.0);
            BiallelicData data05 = tree.simulate_biallelic_data_set(rng05, 0.5);
            BiallelicData data10 = tree.simulate_biallelic_data_set(rng10, 1.0);
            for (unsigned int i = 0; i < data00.get_number_of_patterns(); ++i) {
                const std::vector<unsigned int>& red_allele_counts = data00.get_red_allele_counts(i);
                const std::vector<unsigned int>& allele_counts = data00.get_allele_counts(i);
                unsigned int nreds = 0;
                unsigned int nalleles = 0;
                for (unsigned int j = 0; j < allele_counts.size(); ++j) {
                    nreds += red_allele_counts.at(j);
                    nalleles += allele_counts.at(j);
                }
                if ((nreds == 1) || (nreds == (nalleles - 1))) {
                    singleton_count00 += data00.get_pattern_weight(i);
                }
            }
            for (unsigned int i = 0; i < data05.get_number_of_patterns(); ++i) {
                const std::vector<unsigned int>& red_allele_counts = data05.get_red_allele_counts(i);
                const std::vector<unsigned int>& allele_counts = data05.get_allele_counts(i);
                unsigned int nreds = 0;
                unsigned int nalleles = 0;
                for (unsigned int j = 0; j < allele_counts.size(); ++j) {
                    nreds += red_allele_counts.at(j);
                    nalleles += allele_counts.at(j);
                }
                if ((nreds == 1) || (nreds == (nalleles - 1))) {
                    singleton_count05 += data05.get_pattern_weight(i);
                }
            }
            for (unsigned int i = 0; i < data10.get_number_of_patterns(); ++i) {
                const std::vector<unsigned int>& red_allele_counts = data10.get_red_allele_counts(i);
                const std::vector<unsigned int>& allele_counts = data10.get_allele_counts(i);
                unsigned int nreds = 0;
                unsigned int nalleles = 0;
                for (unsigned int j = 0; j < allele_counts.size(); ++j) {
                    nreds += red_allele_counts.at(j);
                    nalleles += allele_counts.at(j);
                }
                if ((nreds == 1) || (nreds == (nalleles - 1))) {
                    singleton_count10 += data10.get_pattern_weight(i);
                }
            }
        }
        REQUIRE(singleton_count10 > 0);
        REQUIRE(singleton_count05 > 0);
        REQUIRE(singleton_count00 == 0);
        REQUIRE((double)singleton_count05 == Approx(singleton_count10 * 0.5).epsilon(0.05));
    }
}

TEST_CASE("Testing singleton acquisition bias with charsets", "[ComparisonPopulationTree]") {
    SECTION("Testing for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        ComparisonPopulationTree tree(nex_path, ' ',
                true,   // population_name_is_prefix
                false,  // diploid
                false,  // dominant
                false,  // constant sites removed
                true,   // validate
                false,  // strict on constant
                false,  // strict on missing
                false,  // strict on triallelic
                2.0,    // ploidy
                true    // store charset info
                );

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        BiallelicData data1 = tree.simulate_linked_biallelic_data_set(rng, 1.0, false, true);
        BiallelicData data0 = tree.simulate_linked_biallelic_data_set(rng, 0.0, false, true);

        BiallelicData io_data1;
        BiallelicData io_data0;

        std::string io_nex_path1 = "data/tmp-data-test4.nex";
        std::string io_nex_path0 = "data/tmp-data-test5.nex";
        std::ofstream out1;
        std::ofstream out0;
        out1.open(io_nex_path1);
        out0.open(io_nex_path0);
        data1.write_nexus(out1, '-');
        data0.write_nexus(out0, '-');
        out1.close();
        out0.close();
        REQUIRE(path::exists(io_nex_path1));
        REQUIRE(path::exists(io_nex_path0));
        io_data1 = BiallelicData(io_nex_path1,
                '-',   // pop name delimiter
                true,  // pop name is prefix
                false, // genotypes are diploid
                false, // markers are dominant
                true,  // validate
                true   // store seq loci info
                );
        io_data0 = BiallelicData(io_nex_path0,
                '-',   // pop name delimiter
                true,  // pop name is prefix
                false, // genotypes are diploid
                false, // markers are dominant
                true,  // validate
                true   // store seq loci info
                );

        REQUIRE(data1.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data1.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data1.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data1.markers_are_dominant() == false);
        REQUIRE(data1.get_locus_end_indices() == tree.get_data().get_locus_end_indices());
        REQUIRE(data0.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data0.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data0.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data0.markers_are_dominant() == false);
        REQUIRE(data0.get_locus_end_indices() == tree.get_data().get_locus_end_indices());
        REQUIRE(tree.get_data().markers_are_dominant() == false);

        REQUIRE(io_data1.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(io_data1.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(io_data1.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(io_data1.markers_are_dominant() == false);
        REQUIRE(io_data1.get_locus_end_indices() == tree.get_data().get_locus_end_indices());
        REQUIRE(io_data0.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(io_data0.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(io_data0.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(io_data0.markers_are_dominant() == false);
        REQUIRE(io_data0.get_locus_end_indices() == tree.get_data().get_locus_end_indices());

        unsigned int singleton_count10 = 0;
        unsigned int singleton_count05 = 0;
        unsigned int singleton_count00 = 0;
        unsigned int io_singleton_count10 = 0;
        unsigned int io_singleton_count05 = 0;
        unsigned int io_singleton_count00 = 0;
        RandomNumberGenerator rng10 = RandomNumberGenerator(123);
        RandomNumberGenerator rng05 = RandomNumberGenerator(123);
        RandomNumberGenerator rng00 = RandomNumberGenerator(123);
        for (unsigned int rep = 0; rep < 100; ++rep) {
            BiallelicData data00 = tree.simulate_linked_biallelic_data_set(rng00, 0.0, false, true);
            BiallelicData data05 = tree.simulate_linked_biallelic_data_set(rng05, 0.5, false, true);
            BiallelicData data10 = tree.simulate_linked_biallelic_data_set(rng10, 1.0, false, true);
            BiallelicData io_data00;
            BiallelicData io_data05;
            BiallelicData io_data10;

            std::string io_nex_path00 = "data/tmp-data-test6-" + std::to_string(rep) + ".nex";
            std::string io_nex_path05 = "data/tmp-data-test7-" + std::to_string(rep) + ".nex";
            std::string io_nex_path10 = "data/tmp-data-test8-" + std::to_string(rep) + ".nex";
            std::ofstream out00;
            std::ofstream out05;
            std::ofstream out10;
            out00.open(io_nex_path00);
            out05.open(io_nex_path05);
            out10.open(io_nex_path10);
            data00.write_nexus(out00, '-');
            data05.write_nexus(out05, '-');
            data10.write_nexus(out10, '-');
            out00.close();
            out05.close();
            out10.close();
            REQUIRE(path::exists(io_nex_path00));
            REQUIRE(path::exists(io_nex_path05));
            REQUIRE(path::exists(io_nex_path10));
            io_data00 = BiallelicData(io_nex_path00,
                    '-',   // pop name delimiter
                    true,  // pop name is prefix
                    false, // genotypes are diploid
                    false, // markers are dominant
                    true,  // validate
                    true   // store seq loci info
                    );
            io_data05 = BiallelicData(io_nex_path05,
                    '-',   // pop name delimiter
                    true,  // pop name is prefix
                    false, // genotypes are diploid
                    false, // markers are dominant
                    true,  // validate
                    true   // store seq loci info
                    );
            io_data10 = BiallelicData(io_nex_path10,
                    '-',   // pop name delimiter
                    true,  // pop name is prefix
                    false, // genotypes are diploid
                    false, // markers are dominant
                    true,  // validate
                    true   // store seq loci info
                    );
            REQUIRE(data00.get_number_of_patterns() == io_data00.get_number_of_patterns());
            REQUIRE(data05.get_number_of_patterns() == io_data05.get_number_of_patterns());
            REQUIRE(data10.get_number_of_patterns() == io_data10.get_number_of_patterns());

            for (unsigned int i = 0; i < data00.get_number_of_patterns(); ++i) {
                const std::vector<unsigned int>& red_allele_counts = data00.get_red_allele_counts(i);
                const std::vector<unsigned int>& allele_counts = data00.get_allele_counts(i);
                unsigned int nreds = 0;
                unsigned int nalleles = 0;
                for (unsigned int j = 0; j < allele_counts.size(); ++j) {
                    nreds += red_allele_counts.at(j);
                    nalleles += allele_counts.at(j);
                }
                if ((nreds == 1) || (nreds == (nalleles - 1))) {
                    singleton_count00 += data00.get_pattern_weight(i);
                }

                const std::vector<unsigned int>& io_red_allele_counts = io_data00.get_red_allele_counts(i);
                const std::vector<unsigned int>& io_allele_counts = io_data00.get_allele_counts(i);
                unsigned int io_nreds = 0;
                unsigned int io_nalleles = 0;
                for (unsigned int j = 0; j < io_allele_counts.size(); ++j) {
                    io_nreds += io_red_allele_counts.at(j);
                    io_nalleles += io_allele_counts.at(j);
                }
                if ((io_nreds == 1) || (io_nreds == (io_nalleles - 1))) {
                    io_singleton_count00 += io_data00.get_pattern_weight(i);
                }
            }
            for (unsigned int i = 0; i < data05.get_number_of_patterns(); ++i) {
                const std::vector<unsigned int>& red_allele_counts = data05.get_red_allele_counts(i);
                const std::vector<unsigned int>& allele_counts = data05.get_allele_counts(i);
                unsigned int nreds = 0;
                unsigned int nalleles = 0;
                for (unsigned int j = 0; j < allele_counts.size(); ++j) {
                    nreds += red_allele_counts.at(j);
                    nalleles += allele_counts.at(j);
                }
                if ((nreds == 1) || (nreds == (nalleles - 1))) {
                    singleton_count05 += data05.get_pattern_weight(i);
                }

                const std::vector<unsigned int>& io_red_allele_counts = io_data05.get_red_allele_counts(i);
                const std::vector<unsigned int>& io_allele_counts = io_data05.get_allele_counts(i);
                unsigned int io_nreds = 0;
                unsigned int io_nalleles = 0;
                for (unsigned int j = 0; j < io_allele_counts.size(); ++j) {
                    io_nreds += io_red_allele_counts.at(j);
                    io_nalleles += io_allele_counts.at(j);
                }
                if ((io_nreds == 1) || (io_nreds == (io_nalleles - 1))) {
                    io_singleton_count05 += io_data05.get_pattern_weight(i);
                }
            }
            for (unsigned int i = 0; i < data10.get_number_of_patterns(); ++i) {
                const std::vector<unsigned int>& red_allele_counts = data10.get_red_allele_counts(i);
                const std::vector<unsigned int>& allele_counts = data10.get_allele_counts(i);
                unsigned int nreds = 0;
                unsigned int nalleles = 0;
                for (unsigned int j = 0; j < allele_counts.size(); ++j) {
                    nreds += red_allele_counts.at(j);
                    nalleles += allele_counts.at(j);
                }
                if ((nreds == 1) || (nreds == (nalleles - 1))) {
                    singleton_count10 += data10.get_pattern_weight(i);
                }

                const std::vector<unsigned int>& io_red_allele_counts = io_data10.get_red_allele_counts(i);
                const std::vector<unsigned int>& io_allele_counts = io_data10.get_allele_counts(i);
                unsigned int io_nreds = 0;
                unsigned int io_nalleles = 0;
                for (unsigned int j = 0; j < io_allele_counts.size(); ++j) {
                    io_nreds += io_red_allele_counts.at(j);
                    io_nalleles += io_allele_counts.at(j);
                }
                if ((io_nreds == 1) || (io_nreds == (io_nalleles - 1))) {
                    io_singleton_count10 += io_data10.get_pattern_weight(i);
                }
            }
        }
        REQUIRE(singleton_count10 > 0);
        REQUIRE(singleton_count05 > 0);
        REQUIRE(singleton_count00 == 0);
        REQUIRE((double)singleton_count05 == Approx(singleton_count10 * 0.5).epsilon(0.05));

        REQUIRE(io_singleton_count10 > 0);
        REQUIRE(io_singleton_count05 > 0);
        REQUIRE(io_singleton_count00 == 0);
        REQUIRE((double)io_singleton_count05 == Approx(io_singleton_count10 * 0.5).epsilon(0.05));
    }
}

TEST_CASE("Testing scaling of dataset simulation for singleton",
        "[ComparisonPopulationTree]") {

    SECTION("Testing simulate_biallelic_data_set for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        // Two times the expected gene tree root height
        double expected_mean = 2.0 * (theta * (1.0 - (1.0 / nlineages)));

        std::string nex_path = "data/dummy-singleton-n2-100k.nex";
        // Need to keep constant characters
        ComparisonPopulationTree tree(nex_path, '-', true, false, false, false);
        REQUIRE(tree.get_leaf_node_count() == 1);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.set_child_population_size(0, (Ne * mu));

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        BiallelicData data;
        SampleSummarizer<double> divergence;
        unsigned int nsamples = 100;

        for (unsigned int i = 0; i < nsamples; ++i) {
            data = tree.simulate_biallelic_data_set(rng, 1.0, false);
            REQUIRE(data.get_number_of_sites() == 100000);
            double x = (double)data.get_number_of_variable_sites() / (double)data.get_number_of_sites();
            divergence.add_sample(x);
        }

        std::cout << "Expected mean divergence: " << expected_mean << "\n";
        std::cout << "Simulated mean divergence: " << divergence.mean() << "\n";
        REQUIRE(divergence.mean() == Approx(expected_mean).epsilon(0.0001));

        tree.set_mutation_rate(mu);
        tree.set_root_population_size(Ne);
        tree.set_child_population_size(0, Ne);
        tree.set_root_height(0.1 / mu);

        SampleSummarizer<double> divergence2;

        for (unsigned int i = 0; i < nsamples; ++i) {
            data = tree.simulate_biallelic_data_set(rng, 1.0, false);
            REQUIRE(data.get_number_of_sites() == 100000);
            double x = (double)data.get_number_of_variable_sites() / (double)data.get_number_of_sites();
            divergence2.add_sample(x);
        }

        std::cout << "Simulated mean divergence: " << divergence2.mean() << "\n";
        REQUIRE(divergence2.mean() == Approx(expected_mean).epsilon(0.0001));
    }
}

TEST_CASE("Testing scaling of dataset simulation for singleton with charsets",
        "[ComparisonPopulationTree]") {

    SECTION("Testing simulate_biallelic_data_set for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        // Two times the expected gene tree root height
        double expected_mean = 2.0 * (theta * (1.0 - (1.0 / nlineages)));

        std::string nex_path = "data/dummy-singleton-n2-100k.nex";
        // Need to keep constant characters
        ComparisonPopulationTree tree(nex_path, '-',
                true,   // population_name_is_prefix
                false,  // diploid
                false,  // dominant
                false,  // constant sites removed
                true,   // validate
                false,  // strict on constant
                false,  // strict on missing
                false,  // strict on triallelic
                2.0,    // ploidy
                true    // store charset info
                );
        REQUIRE(tree.get_leaf_node_count() == 1);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.set_child_population_size(0, (Ne * mu));

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(12345678);

        BiallelicData data;
        SampleSummarizer<double> divergence;
        unsigned int nsamples = 100;

        for (unsigned int i = 0; i < nsamples; ++i) {
            data = tree.simulate_linked_biallelic_data_set(rng, 1.0, false);
            REQUIRE(data.get_number_of_sites() == 100000);
            double x = (double)data.get_number_of_variable_sites() / (double)data.get_number_of_sites();
            divergence.add_sample(x);
        }

        std::cout << "Expected mean divergence: " << expected_mean << "\n";
        std::cout << "Simulated mean divergence: " << divergence.mean() << "\n";
        REQUIRE(divergence.mean() == Approx(expected_mean).epsilon(0.0002));

        tree.set_mutation_rate(mu);
        tree.set_root_population_size(Ne);
        tree.set_child_population_size(0, Ne);
        tree.set_root_height(0.1 / mu);

        SampleSummarizer<double> divergence2;

        for (unsigned int i = 0; i < nsamples; ++i) {
            data = tree.simulate_linked_biallelic_data_set(rng, 1.0, false);
            REQUIRE(data.get_number_of_sites() == 100000);
            double x = (double)data.get_number_of_variable_sites() / (double)data.get_number_of_sites();
            divergence2.add_sample(x);
        }

        std::cout << "Simulated mean divergence: " << divergence2.mean() << "\n";
        REQUIRE(divergence2.mean() == Approx(expected_mean).epsilon(0.0002));

        REQUIRE(data.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends;
        unsigned int end_idx = 0;
        for (unsigned int i = 0; i < 100; ++i) {
            end_idx += 1000;
            expected_locus_ends.push_back(end_idx - 1);
        }
        REQUIRE(data.get_contiguous_pattern_indices().size() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_locus_end_indices() == expected_locus_ends);
    }
}

TEST_CASE("Testing scaling of simulation of loci for singleton",
        "[ComparisonPopulationTree]") {

    SECTION("Testing simulate_complete_biallelic_data_set for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        // Two times the expected gene tree root height
        double expected_mean = 2.0 * (theta * (1.0 - (1.0 / nlineages)));

        std::string nex_path = "data/dummy-singleton-n2-100k.nex";
        // Need to keep constant characters
        ComparisonPopulationTree tree(nex_path, '-', true, false, false, false);
        REQUIRE(tree.get_leaf_node_count() == 1);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.set_child_population_size(0, (Ne * mu));

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        SampleSummarizer<double> divergence;
        unsigned int nsamples = 100;

        for (unsigned int i = 0; i < nsamples; ++i) {
            auto data_nloci = tree.simulate_complete_biallelic_data_set(rng, 1000, 1.0, false);
            REQUIRE(data_nloci.second == 100);
            REQUIRE(data_nloci.first.get_number_of_sites() == 100000);
            double x = (double)data_nloci.first.get_number_of_variable_sites() / (double)data_nloci.first.get_number_of_sites();
            divergence.add_sample(x);
        }

        std::cout << "Expected mean divergence: " << expected_mean << "\n";
        std::cout << "Simulated mean divergence: " << divergence.mean() << "\n";
        REQUIRE(divergence.mean() == Approx(expected_mean).epsilon(0.0001));

        tree.set_mutation_rate(mu);
        tree.set_root_population_size(Ne);
        tree.set_child_population_size(0, Ne);
        tree.set_root_height(0.1 / mu);

        SampleSummarizer<double> divergence2;

        for (unsigned int i = 0; i < nsamples; ++i) {
            auto data_nloci = tree.simulate_complete_biallelic_data_set(rng, 1000, 1.0, false);
            REQUIRE(data_nloci.second == 100);
            REQUIRE(data_nloci.first.get_number_of_sites() == 100000);
            double x = (double)data_nloci.first.get_number_of_variable_sites() / (double)data_nloci.first.get_number_of_sites();
            divergence2.add_sample(x);
        }

        std::cout << "Simulated mean divergence: " << divergence2.mean() << "\n";
        REQUIRE(divergence2.mean() == Approx(expected_mean).epsilon(0.0001));
    }
}

TEST_CASE("Testing scaling of simulation of one variable site per locus for singleton",
        "[ComparisonPopulationTree]") {

    SECTION("Testing simulate_data_set_max_one_variable_site_per_locus for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        // Two times the expected gene tree root height
        double expected_mean = 2.0 * (theta * (1.0 - (1.0 / nlineages)));
        double epsilon = 0.0001;

        std::string nex_path = "data/dummy-singleton-n2-100k.nex";
        // Need to keep constant characters
        ComparisonPopulationTree tree(nex_path, '-', true, false, false, false);
        REQUIRE(tree.get_leaf_node_count() == 1);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.set_child_population_size(0, (Ne * mu));

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        /* RandomNumberGenerator rng = RandomNumberGenerator(54321); */
        RandomNumberGenerator rng = RandomNumberGenerator(12345678);

        SampleSummarizer<double> divergence;
        unsigned int nsamples = 100;

        for (unsigned int i = 0; i < nsamples; ++i) {
            auto data_nloci = tree.simulate_data_set_max_one_variable_site_per_locus(rng, 1000, 1.0, false);
            REQUIRE(data_nloci.second == 100);
            REQUIRE(data_nloci.first.get_number_of_sites() <= 100);
            double x = (double)data_nloci.first.get_number_of_variable_sites() / (data_nloci.first.get_number_of_sites() + data_nloci.first.get_number_of_constant_sites_removed());
            divergence.add_sample(x);
        }

        std::cout << "Expected mean divergence: " << expected_mean << "\n";
        std::cout << "Simulated mean divergence: " << divergence.mean() << "\n";
        // By stopping the simulation of sites for each locus when we get a
        // variable site, we are under-sampling sites from loci with older
        // coalescence times (over-sampling sites from loci with younger
        // coalescence times), so we should under estimate diversity:
        REQUIRE(divergence.mean() < (expected_mean - epsilon));

        tree.set_mutation_rate(mu);
        tree.set_root_population_size(Ne);
        tree.set_child_population_size(0, Ne);
        tree.set_root_height(0.1 / mu);

        SampleSummarizer<double> divergence2;

        for (unsigned int i = 0; i < nsamples; ++i) {
            auto data_nloci = tree.simulate_data_set_max_one_variable_site_per_locus(rng, 1000, 1.0, false);
            REQUIRE(data_nloci.second == 100);
            REQUIRE(data_nloci.first.get_number_of_sites() <= 100);
            double x = (double)data_nloci.first.get_number_of_variable_sites() / (data_nloci.first.get_number_of_sites() + data_nloci.first.get_number_of_constant_sites_removed());
            divergence2.add_sample(x);
        }

        std::cout << "Simulated mean divergence: " << divergence2.mean() << "\n";
        // By stopping the simulation of sites for each locus when we get a
        // variable site, we are under-sampling sites from loci with older
        // coalescence times (over-sampling sites from loci with younger
        // coalescence times), so we should under estimate diversity:
        REQUIRE(divergence2.mean() < (expected_mean - epsilon));
    }
}

TEST_CASE("Testing scaling of simulation of one variable site per locus for singleton with charsets",
        "[ComparisonPopulationTree]") {

    SECTION("Testing one SNP per locus for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        // Two times the expected gene tree root height
        double expected_mean = 2.0 * (theta * (1.0 - (1.0 / nlineages)));
        double epsilon = 0.0001;

        std::string nex_path = "data/dummy-singleton-n2-100k.nex";
        // Need to keep constant characters
        ComparisonPopulationTree tree(nex_path, '-',
                true,   // population_name_is_prefix
                false,  // diploid
                false,  // dominant
                false,  // constant sites removed
                true,   // validate
                true,   // strict on constant
                true,   // strict on missing
                true,   // strict on triallelic
                2.0,    // ploidy
                true    // store charset info
                );
        REQUIRE(tree.get_leaf_node_count() == 1);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.set_child_population_size(0, (Ne * mu));

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        /* RandomNumberGenerator rng = RandomNumberGenerator(54321); */
        RandomNumberGenerator rng = RandomNumberGenerator(12345678);

        SampleSummarizer<double> divergence;
        unsigned int nsamples = 100;

        BiallelicData data;
        for (unsigned int i = 0; i < nsamples; ++i) {
            data = tree.simulate_linked_biallelic_data_set(rng, 1.0, true, true);
            REQUIRE(data.get_number_of_sites() <= 100);
            double x = (double)data.get_number_of_variable_sites() / (data.get_number_of_sites() + data.get_number_of_constant_sites_removed());
            divergence.add_sample(x);
        }

        std::cout << "Expected mean divergence: " << expected_mean << "\n";
        std::cout << "Simulated mean divergence: " << divergence.mean() << "\n";
        // By stopping the simulation of sites for each locus when we get a
        // variable site, we are under-sampling sites from loci with older
        // coalescence times (over-sampling sites from loci with younger
        // coalescence times), so we should under estimate diversity:
        REQUIRE(divergence.mean() < (expected_mean - epsilon));

        tree.set_mutation_rate(mu);
        tree.set_root_population_size(Ne);
        tree.set_child_population_size(0, Ne);
        tree.set_root_height(0.1 / mu);

        SampleSummarizer<double> divergence2;

        for (unsigned int i = 0; i < nsamples; ++i) {
            data = tree.simulate_linked_biallelic_data_set(rng, 1.0, true, true);
            REQUIRE(data.get_number_of_sites() <= 100);
            double x = (double)data.get_number_of_variable_sites() / (data.get_number_of_sites() + data.get_number_of_constant_sites_removed());
            divergence2.add_sample(x);
        }

        std::cout << "Simulated mean divergence: " << divergence2.mean() << "\n";
        // By stopping the simulation of sites for each locus when we get a
        // variable site, we are under-sampling sites from loci with older
        // coalescence times (over-sampling sites from loci with younger
        // coalescence times), so we should under estimate diversity:
        REQUIRE(divergence2.mean() < (expected_mean - epsilon));

        REQUIRE(data.has_seq_loci_info() == false);
        REQUIRE(data.get_contiguous_pattern_indices().size() == 0);
        REQUIRE(data.get_locus_end_indices().size() == 0);


        data = tree.simulate_linked_biallelic_data_set(rng, 1.0, false, true);

        std::vector<unsigned int> expected_locus_ends;
        unsigned int end_idx = 0;
        for (unsigned int i = 0; i < 100; ++i) {
            end_idx += 1000;
            expected_locus_ends.push_back(end_idx - 1);
        }
        REQUIRE(data.get_contiguous_pattern_indices().size() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_locus_end_indices().size() == tree.get_data().get_locus_end_indices().size());
        REQUIRE(data.get_locus_end_indices() == expected_locus_ends);
    }
}


TEST_CASE("Testing errors when trying to sim loci with a template with constant characters removed",
        "[ComparisonPopulationTree]") {

    SECTION("Testing for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        // Two times the expected gene tree root height
        double expected_mean = 2.0 * (theta * (1.0 - (1.0 / nlineages)));
        double epsilon = 0.0001;

        std::string nex_path = "data/dummy-singleton-n2-100k.nex";
        // Need to keep constant characters
        ComparisonPopulationTree tree(nex_path, '-',
                true,   // population_name_is_prefix
                false,  // diploid
                false,  // dominant
                true,   // constant sites removed
                true,   // validate
                true,   // strict on constant
                true,   // strict on missing
                true,   // strict on triallelic
                2.0,    // ploidy
                true    // store charset info
                );
        REQUIRE(tree.get_leaf_node_count() == 1);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.set_child_population_size(0, (Ne * mu));
        tree.set_root_height(0.1);
        tree.set_freq_1(0.5);
        tree.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(12345678);

        // constant sites were removed, so all of these should fail
        REQUIRE_THROWS_AS(tree.simulate_linked_biallelic_data_set(rng, 1.0, false, false), EcoevolityBiallelicDataError &);
        REQUIRE_THROWS_AS(tree.simulate_linked_biallelic_data_set(rng, 1.0, true, true), EcoevolityBiallelicDataError &);
        REQUIRE_THROWS_AS(tree.simulate_linked_biallelic_data_set(rng, 1.0, true, false), EcoevolityBiallelicDataError &);
        REQUIRE_THROWS_AS(tree.simulate_linked_biallelic_data_set(rng, 1.0, false, true), EcoevolityBiallelicDataError &);

        REQUIRE_THROWS_AS(tree.simulate_complete_biallelic_data_set(rng, 1, 1.0, true), EcoevolityBiallelicDataError &);
        REQUIRE_THROWS_AS(tree.simulate_complete_biallelic_data_set(rng, 1, 1.0, false), EcoevolityBiallelicDataError &);
        REQUIRE_THROWS_AS(tree.simulate_complete_biallelic_data_set(rng, 10, 1.0, true), EcoevolityBiallelicDataError &);
        REQUIRE_THROWS_AS(tree.simulate_complete_biallelic_data_set(rng, 10, 1.0, false), EcoevolityBiallelicDataError &);
    }
}

TEST_CASE("Testing draw_from_prior for fully fixed", "[ComparisonPopulationTree]") {

    SECTION("Testing draw_from_prior for fully fixed pair") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.001);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_root_population_size() == 0.001);
            REQUIRE(tree.get_child_population_size(0) == 0.001);
            REQUIRE(tree.get_child_population_size(1) == 0.001);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);
        }
    }
}

TEST_CASE("Testing draw_from_prior for constrained sizes", "[ComparisonPopulationTree]") {

    SECTION("Testing draw_from_prior for constrained sizes") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        SampleSummarizer<double> pop_size;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() == tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() == tree.get_child_population_size(1));

            pop_size.add_sample(tree.get_root_population_size());
        }

        REQUIRE(pop_size.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
    }
}

TEST_CASE("Testing draw_from_prior for unconstrained sizes", "[ComparisonPopulationTree]") {

    SECTION("Testing draw_from_prior for unconstrained sizes") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.set_all_population_sizes(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            pop_size_root.add_sample(tree.get_root_population_size());
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
        }

        REQUIRE(pop_size_root.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_root.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_0.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_0.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_1.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_1.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
    }
}

TEST_CASE("Testing draw_from_prior for fully parameterized", "[ComparisonPopulationTree]") {

    SECTION("Testing draw_from_prior for fully parameterized") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(0.5, 0.8);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.set_all_population_sizes(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.estimate_state_frequencies();
        tree.set_freq_1(0.5);

        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;
        SampleSummarizer<double> rate;
        SampleSummarizer<double> f;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            pop_size_root.add_sample(tree.get_root_population_size());
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
            rate.add_sample(tree.get_mutation_rate());
            f.add_sample(tree.get_freq_1());
        }

        REQUIRE(pop_size_root.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_root.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_0.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_0.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_1.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_1.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(f.mean() == Approx(f_prior->get_mean()).epsilon(0.01));
        REQUIRE(f.variance() == Approx(f_prior->get_variance()).epsilon(0.01));
        REQUIRE(rate.mean() == Approx(rate_prior->get_mean()).epsilon(0.1));
        REQUIRE(rate.variance() == Approx(rate_prior->get_variance()).epsilon(0.1));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood = -31.77866581319647
TEST_CASE("Testing simple likelihood of DirichletPopulationTree", "[DirichletPopulationTree]") {

    SECTION("Testing constructor and likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing threaded likelihood of DirichletPopulationTree", "[DirichletPopulationTree]") {

    SECTION("Testing constructor and threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing over-threaded likelihood of DirichletPopulationTree", "[DirichletPopulationTree]") {

    SECTION("Testing constructor and over-threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(10);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(10);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood             = -248.93254688526213
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing hemi129.nex DPT likelihood (0.01, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex DPT threaded likelihood (0.01, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7099.716015109998
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.01, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.01, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -6986.120524781545
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.01, 10.0, 1.0, 1.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6986.120524781545));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.01, 10.0, 1.0, 1.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(5);
        REQUIRE(l == Approx(-6986.120524781545));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}



// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -328.39238828878365
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex DPT likelihood (0.0, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex DPT threaded likelihood (0.0, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(7);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(7);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7256.501742344454
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.0, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.0, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7223.362711937651
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.0, 10.0, 1.0, 1.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7223.362711937651));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.0, 10.0, 1.0, 1.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7223.362711937651));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -227.41048391087554
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex DPT likelihood (0.2, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex DPT threaded likelihood (0.2, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7304.180743441677
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.2, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.2, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7405.145951634711
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.2, 10.0, 1.0, 1.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7405.145951634711));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.2, 10.0, 1.0, 1.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(6);
        REQUIRE(l == Approx(-7405.145951634711));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}




// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -327.7437811413033
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex DPT likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex DPT threaded likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -6472.856486972301
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(5);
        REQUIRE(l == Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(3);
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -6494.774924871097
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.03, 10.0, 10.0, 10.0/19.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6494.774924871097));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.03, 10.0, 10.0, 10.0/19.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-6494.774924871097));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
  
  
// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -265.0023534261969
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex DPT likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex DPT threaded likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(1);
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -10163.468886613919
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(5);
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -10999.288193543642
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.03, 10.0, 10.0/19.0, 10.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-10999.288193543642));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.03, 10.0, 10.0/19.0, 10.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-10999.288193543642));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}



// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -224.40177558289847
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex DPT likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_mean_population_size(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-224.40177558289847));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex DPT threaded likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_mean_population_size(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(100000);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-224.40177558289847));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -8158.88094671241
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_mean_population_size(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-8158.88094671241));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_mean_population_size(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(1000000);
        REQUIRE(l == Approx(-8158.88094671241));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -8034.250341980543
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.03, 111.1, 10.0/19.0, 10.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_mean_population_size(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-8034.250341980543));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.03, 111.1, 10.0/19.0, 10.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_mean_population_size(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-8034.250341980543));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

TEST_CASE("Testing simple prior of DirichletPopulationTree", "[DirichletPopulationTree]") {

    SECTION("Testing constructor and prior calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);

        tree.set_root_height(0.1);
        tree.set_mean_population_size(2.0/100.0);
        tree.set_freq_1(0.95);

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(100.0));
        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));
        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(2.0, 1.0));

        // height           -5.3948298140119091
        // size             -155.90663080917298
        // size multipliers 0.69314718055994529 (relative is 0.0)
        // f1               0.64185388617239469 
        // total            -160.6596067370125
        
        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-155.90663080917298));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == Approx(0.64185388617239469));

        double d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-160.6596067370125));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-160.6596067370125));


        tree.store_prior_density();
        tree.constrain_state_frequencies();

        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-155.90663080917298));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == 0.0);

        d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-161.30146062318488));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-161.30146062318488));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-160.6596067370125));

        tree.restore_prior_density();
        REQUIRE(tree.get_log_prior_density_value() == Approx(-160.6596067370125));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-160.6596067370125));
        REQUIRE(tree.get_degree_of_root() == 2);

        tree.fix_population_size_multipliers();

        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-155.90663080917298));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == 0.0);

        d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-161.30146062318488));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-161.30146062318488));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-160.6596067370125));

        tree.store_prior_density();

        REQUIRE(tree.get_log_prior_density_value() == Approx(-161.30146062318488));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-161.30146062318488));
    }
}

TEST_CASE("Testing hemi129.nex state manipulation for DirichletPopulationTree", "[DirichletPopulationTree]") {

    SECTION("Testing state manipulation") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(10.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(1.5, 3.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        std::vector<double> parameters {10.0, 10.0, 10.0};
        tree.set_population_size_multiplier_prior(std::make_shared<DirichletDistribution>(parameters));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_height(0.01);
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1(0.5);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        std::vector<double> proportions {1.0/3.0, 1.0/3.0, 1.0/3.0};
        tree.set_population_sizes_as_proportions(proportions);

        REQUIRE(tree.is_dirty());

        REQUIRE(tree.get_number_of_likelihood_calculations() == 0);

        tree.compute_log_likelihood_and_prior();

        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-474.97145724683253));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 1);

        tree.store_state();
        tree.set_root_height(0.0);
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-474.87145724683256));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-474.97145724683253));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 2);

        tree.restore_state();
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-474.97145724683253));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-474.97145724683253));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.compute_log_likelihood_and_prior();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-474.97145724683253));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-474.97145724683253));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.store_state();
        tree.set_root_height(0.0);
        tree.compute_log_likelihood_and_prior();
        REQUIRE(tree.get_number_of_likelihood_calculations() == 4);


        tree.store_state();
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-476.87145724683256));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-474.87145724683256));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 5);


        tree.store_state();
        tree.set_root_height(0.03);
        tree.set_freq_1(0.05);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.05);
        REQUIRE(tree.get_freq_0() == Approx(0.95));
        REQUIRE(tree.get_u() == Approx(10.0));
        REQUIRE(tree.get_v() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-475.03904202098477));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-476.87145724683256));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 6);


        tree.store_state();
        tree.set_freq_1(0.95);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.95);
        REQUIRE(tree.get_freq_0() == Approx(0.05));
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-479.45570048973445));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-475.03904202098477));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 7);


        tree.store_state();
        tree.set_mean_population_size(2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(111.1 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-46.130811372643343));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-479.45570048973445));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 8);

        tree.store_state();
        tree.constrain_state_frequencies();
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-477.01996092335042));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-46.130811372643343));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 9);


        tree.store_state();
        tree.fold_patterns();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-477.01996092335042));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-477.01996092335042));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 10);


        tree.store_state();
        tree.fix_population_size_multipliers();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.35742912931147));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-477.01996092335042));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 11);

        tree.store_state();
        tree.estimate_mutation_rate();
        REQUIRE(! tree.mutation_rate_is_fixed());
        tree.set_mutation_rate(1.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());
        tree.set_mutation_rate_prior(std::make_shared<GammaDistribution>(10.0, 0.1));
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.35742912931147));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 12);

        tree.store_state();
        tree.set_mutation_rate(2.0);
        tree.set_root_height(0.1);
        tree.set_mean_population_size(2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 2.0);
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(20.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-206.13340567945255));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        tree.restore_state();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 9fb5b3b7817a0bd4a21e3f90132f132cca72ce4e)
// SNAPP v1.3.0 (master 4f3f0f7366798f4fb38b766c15f6426a75ddf71e)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0/19.0
// Log likelihood            = -23.81984255023975
// Log likelihood correction = -6.87935580446044
//
// With constant sites inclucded and m_bUseNonPolymorphic = true
// Log likelihood            = -55.01646493341547
// Log likelihood correction = -6.87935580446044
TEST_CASE("Testing affect of constant sites on likelihood of DirichletPopulationTree", "[DirichletPopulationTree]") {

    SECTION("Testing haploid-standard-full-constant.nex") {
        std::string nex_path = "data/haploid-standard-full-constant.nex";
        DirichletPopulationTree t_included(
                nex_path, // path
                ' ',      // pop name delimiter
                true,     // pop name is prefix
                false,    // genotypes are diploid
                false,    // markers are dominant
                false,    // constant sites removed
                true);    // validate
        std::string nex_path2 = "data/haploid-standard-full-constant-removed.nex";
        DirichletPopulationTree t(
                nex_path2, // path
                ' ',       // pop name delimiter
                true,      // pop name is prefix
                false,     // genotypes are diploid
                false,     // markers are dominant
                true,      // constant sites removed
                true);     // validate

        t.set_root_height(0.03);
        t.set_freq_1(0.05);
        t.set_mean_population_size(2.0/(10.0 * 2 * t.get_ploidy()));

        t_included.set_root_height(0.03);
        t_included.set_freq_1(0.05);
        t_included.set_mean_population_size(2.0/(10.0 * 2 * t_included.get_ploidy()));

        double l = t.compute_log_likelihood();
        REQUIRE(l == Approx(-23.81984255023975));
        REQUIRE(t.get_likelihood_correction() == Approx(-6.87935580446044));

        double l_included = t_included.compute_log_likelihood();

        REQUIRE(l_included == Approx(-55.01646493341547));

        REQUIRE(t.get_likelihood_correction() == t_included.get_likelihood_correction());
    }
}

TEST_CASE("Testing affect of constant sites on threaded likelihood of DirichletPopulationTree", "[DirichletPopulationTree]") {

    SECTION("Testing haploid-standard-full-constant.nex") {
        std::string nex_path = "data/haploid-standard-full-constant.nex";
        DirichletPopulationTree t_included(
                nex_path, // path
                ' ',      // pop name delimiter
                true,     // pop name is prefix
                false,    // genotypes are diploid
                false,    // markers are dominant
                false,    // constant sites removed
                true);    // validate
        std::string nex_path2 = "data/haploid-standard-full-constant-removed.nex";
        DirichletPopulationTree t(
                nex_path2, // path
                ' ',       // pop name delimiter
                true,      // pop name is prefix
                false,     // genotypes are diploid
                false,     // markers are dominant
                true,      // constant sites removed
                true);     // validate

        t.set_root_height(0.03);
        t.set_freq_1(0.05);
        t.set_mean_population_size(2.0/(10.0 * 2 * t.get_ploidy()));

        t_included.set_root_height(0.03);
        t_included.set_freq_1(0.05);
        t_included.set_mean_population_size(2.0/(10.0 * 2 * t_included.get_ploidy()));

        double l = t.compute_log_likelihood(4);
        REQUIRE(l == Approx(-23.81984255023975));
        REQUIRE(t.get_likelihood_correction() == Approx(-6.87935580446044));

        double l_included = t_included.compute_log_likelihood(2);

        REQUIRE(l_included == Approx(-55.01646493341547));

        REQUIRE(t.get_likelihood_correction() == t_included.get_likelihood_correction());
    }
}








// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood = -31.77866581319647
TEST_CASE("Testing simple likelihood of ComparisonDirichletPopulationTree", "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing constructor and likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_root().get_label() == "root-pop1");
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing simple threaded likelihood of ComparisonDirichletPopulationTree", "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing constructor and threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_root().get_label() == "root-pop1");
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree state manipulation",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing state manipulation") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(10.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        std::vector<double> alphas {10.0, 20.0, 30.0};
        tree.set_population_size_multiplier_prior(std::make_shared<DirichletDistribution>(alphas));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(1.5, 3.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_height(0.01);
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        std::vector<double> m = {1.0/3.0, 1.0/3.0, 1.0/3.0};
        tree.set_population_sizes_as_proportions(m);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1(0.5);
        REQUIRE(tree.is_dirty());

        REQUIRE(tree.get_number_of_likelihood_calculations() == 0);

        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 1);

        tree.store_state();
        tree.set_root_height(0.0);
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 2);

        tree.restore_state();
        // ComparisonDirichletPopulationTree does not store/restore heights
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.compute_log_likelihood_and_prior();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.store_state();
        tree.set_root_height(0.0);
        tree.compute_log_likelihood_and_prior();
        REQUIRE(tree.get_number_of_likelihood_calculations() == 4);


        tree.store_state();
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 5);


        tree.store_state();
        tree.set_root_height(0.03);
        tree.set_freq_1(0.05);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.05);
        REQUIRE(tree.get_freq_0() == 0.95);
        REQUIRE(tree.get_u() == Approx(10.0));
        REQUIRE(tree.get_v() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-509.99999577402207));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 6);


        tree.store_state();
        tree.set_freq_1(0.95);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_freq_1() == 0.95);
        REQUIRE(tree.get_freq_0() == Approx(0.05));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-514.41665424277176));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-509.99999577402207));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 7);


        tree.store_state();
        tree.set_mean_population_size(2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(111.1 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(111.1 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(111.1 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-81.091765125680695));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-514.41665424277176));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 8);

        tree.store_state();
        tree.constrain_state_frequencies();
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-510.28091467638774));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-81.091765125680695));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 9);


        tree.store_state();
        tree.fold_patterns();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-510.28091467638774));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-510.28091467638774));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 10);


        tree.store_state();
        tree.fix_population_size_multipliers();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.66001422230551));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-510.28091467638774));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 11);

        tree.store_state();
        tree.estimate_mutation_rate();
        REQUIRE(! tree.mutation_rate_is_fixed());
        tree.set_mutation_rate(1.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());
        tree.set_mutation_rate_prior(std::make_shared<GammaDistribution>(10.0, 0.1));
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.66001422230551));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 12);

        tree.store_state();
        tree.set_mutation_rate(2.0);
        tree.set_root_height(0.1);
        tree.set_mean_population_size(2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 2.0);
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(20.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(20.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(20.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-207.43599077244659));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        tree.restore_state();
        // ComparisonDirichletPopulationTree does not store/restore node heights
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        REQUIRE(tree.get_population_sizes_as_proportions() == m);

        tree.store_state();
        tree.estimate_population_size_multipliers();
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        m = {0.5/3.0, 1.0/3.0, 1.5/3.0};
        tree.set_population_sizes_as_proportions(m);
        REQUIRE(tree.get_population_sizes_as_proportions() == m);
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_mean_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_root_population_size() ==   Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(2) * 3.0));
        REQUIRE(tree.get_child_population_size(0) == Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(0) * 3.0));
        REQUIRE(tree.get_child_population_size(1) == Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(1) * 3.0));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(tree.get_log_prior_density_value() == Approx(-504.53672771643153));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 14);

        tree.store_state();
        tree.fix_mean_population_size();
        REQUIRE(tree.get_population_sizes_as_proportions() == m);
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_mean_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_root_population_size() == Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(2) * 3.0));
        REQUIRE(tree.get_child_population_size(0) == Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(0) * 3.0));
        REQUIRE(tree.get_child_population_size(1) == Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(1) * 3.0));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(tree.get_log_prior_density_value() == Approx(-56.876713494126008));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-504.53672771643153));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 15);

        tree.store_state();
        tree.fix_population_size_multipliers();
        REQUIRE(tree.get_population_sizes_as_proportions() == m);
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_mean_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_root_population_size() == Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(2) * 3.0));
        REQUIRE(tree.get_child_population_size(0) == Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(0) * 3.0));
        REQUIRE(tree.get_child_population_size(1) == Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(1) * 3.0));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(tree.get_log_prior_density_value() == Approx(0.22402344985899036));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-56.876713494126008));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 16);

        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.01
// coalescent_rate = 100.0, 200.0, 500.0
// u = 1.0
// v = 1.0
// Log likelihood            = -226.11914854623677
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing ComparisonDirichletPopulationTree likelihood (0.01, 100.0, 200.0, 500.0)",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        // tree.set_child_population_size(0, 2.0/(100.0 * 2 * tree.get_ploidy()));
        // tree.set_child_population_size(1, 2.0/(200.0 * 2 * tree.get_ploidy()));
        // tree.set_root_population_size(2.0/(500.0 * 2 * tree.get_ploidy()));
        tree.set_mean_population_size(0.0085 / 3.0);
        std::vector<double> proportions = {
                (0.005 / 0.0085), 
                (0.0025 / 0.0085), 
                (0.001 / 0.0085)}; 
        tree.set_population_sizes_as_proportions(proportions);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-226.11914854623677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}


// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.00506843962151613554
// coalescent_rate = 2.0 / 0.00018955324120485613
// u = 1.0
// v = 1.0
// Log likelihood            = -277.06960543551577
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonDirichletPopulationTree likelihood (0.00506843962151613554, 2.0 / 0.00018955324120485613)",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.00506843962151613554);
        tree.set_mean_population_size(0.00018955324120485613 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-277.06960543551577));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 9.08323190033687971e-09
// coalescent_rate = 2.0 / 2.47975039926886321e-08
// u = 1.0
// v = 1.0
// Log likelihood            = -221.69627648370943
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonDirichletPopulationTree likelihood (9.08323190033687971e-09, 2.0 / 2.47975039926886321e-08)", "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(9.08323190033687971e-09);
        tree.set_mean_population_size(2.47975039926886321e-08 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-221.69627648370943));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.04921319733994759e-08
// coalescent_rate = 2.0 / 2.75977168733651178e-10
// u = 1.0
// v = 1.0
// Log likelihood            = -324.2737564069293
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonDirichletPopulationTree likelihood (1.04921319733994759e-08, 2.0 / 2.75977168733651178e-10)", "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.04921319733994759e-08);
        tree.set_mean_population_size(2.75977168733651178e-10 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-324.2737564069293));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.012386610001351e-08
// coalescent_rate = 2.0 / 5.75048645855884647e-30
// u = 1.0
// v = 1.0
// Log likelihood            = -1364.1427530000253
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonDirichletPopulationTree likelihood (1.012386610001351e-08, 2.0 / 5.75048645855884647e-30)", "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.012386610001351e-08);
        tree.set_mean_population_size(5.75048645855884647e-30 / 4.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-1364.1427530000253));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.04856228318474786e-08
// coalescent_rate = 2.0 / 4.43934332792563837e-305
// u = 1.0
// v = 1.0
// Log likelihood            = NaN
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonDirichletPopulationTree likelihood (1.04856228318474786e-08, 2.0 / 4.43934332792563837e-305)", "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.04856228318474786e-08);
        tree.set_mean_population_size(4.43934332792563837e-305 / 4.0);
        double l = tree.compute_log_likelihood();
        // REQUIRE(std::isnan(l));
        REQUIRE(std::isinf(l));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.012386610001351e-08
// coalescent_rate = 2.0 / 5.46641122085615013e-09,
//                   2.0 / 7.39871781998828579e-08,
//                   2.0 / 2.71077053326002069e-13
// u = 1.0
// v = 1.0
// Log likelihood            = -96.34394008351177
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing ComparisonDirichletPopulationTree weirdness", "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.012386610001351e-08);
        tree.set_child_population_size(0, 5.46641122085615013e-09 / 4.0);
        tree.set_child_population_size(1, 7.39871781998828579e-08 / 4.0);
        tree.set_root_population_size(2.71077053326002069e-13 / 4.0);
        /* tree.set_mean_population_size(6.621155041482694e-09); */
        /* std::vector<double> proportions = {0.20639945699081685, 2.7935903077461655, 1.0235263017844203e-05}; */
        /* for (unsigned int i = 0; i < proportions.size(); ++i) { */
        /*     proportions.at(i) /= 3.0; */
        /* } */
        /* tree.set_population_sizes_as_proportions(proportions); */
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-96.34394008351177));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 2 lineages and theta of 1.0",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 2 lineages and theta of 3.7",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 3.7;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 2 lineages and theta of 0.17",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 3 lineages and theta of 1.0",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 3 lineages and theta of 1.47",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 1.47;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 3 lineages and theta of 0.17",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 10 lineages and theta of 1.0",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 10 lineages and theta of 1.47",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 1.47;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 10 lineages and theta of 0.17",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing scaling of ComparisonDirichletPopulationTree.simulate_gene_tree for singleton",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        std::string nex_path = "data/singleton-n2.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, false, false);
        tree.estimate_mutation_rate();

        tree.set_mean_population_size(Ne * mu);

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        unsigned int n;
        double mean;
        double variance;
        double sum_devs;
        double d;
        double d_n;
        double mn;
        double mx;
        double x;
        std::shared_ptr<GeneTreeSimNode> gtree;

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            x = gtree->get_height();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);

        REQUIRE(gtree->is_root());
        REQUIRE(gtree->get_number_of_children() == 2);
        REQUIRE(gtree->get_leaf_node_count() == 2);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.0005));
        REQUIRE(mn >= 0.0);

        tree.set_mutation_rate(mu);
        tree.set_mean_population_size(Ne);
        tree.set_root_height(0.1 / mu);

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            x = gtree->get_height();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);

        REQUIRE(gtree->is_root());
        REQUIRE(gtree->get_number_of_children() == 2);
        REQUIRE(gtree->get_leaf_node_count() == 2);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing scaling of ComparisonDirichletPopulationTree.simulate_gene_tree for pair",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing pair") {
        unsigned int Ne_root = 100000;
        unsigned int Ne_0 = 200000;
        unsigned int Ne_1 = 10000;
        double mu = 1e-8;
        double theta_root = 4 * Ne_root * mu;
        double theta_0 = 4 * Ne_0 * mu;
        double theta_1 = 4 * Ne_1 * mu;
        double time = 2e7;

        double expected_mean_root = theta_root * (1.0 - (1.0 / 2));
        double expected_variance_root = expected_mean_root * expected_mean_root;
        expected_mean_root += (time * mu);

        double expected_mean_0 = theta_0 * (1.0 - (1.0 / 10));
        double expected_mean_1 = theta_1 * (1.0 - (1.0 / 10));

        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        tree.estimate_mutation_rate();

        double sum_Ne = Ne_0 + Ne_1 + Ne_root;
        double ave_Ne =  sum_Ne / 3.0;
        std::vector<double> proportions = {
                (Ne_0 / sum_Ne),
                (Ne_1 / sum_Ne),
                (Ne_root / sum_Ne)};
        tree.set_mean_population_size(ave_Ne);
        tree.set_population_sizes_as_proportions(proportions);

        // tree.set_root_population_size(Ne_root);
        // tree.set_child_population_size(0, (Ne_0));
        // tree.set_child_population_size(1, (Ne_1));

        tree.set_root_height(time);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(mu);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        std::shared_ptr<GeneTreeSimNode> gtree;
        SampleSummarizer<double> height_root;
        SampleSummarizer<double> height_0;
        SampleSummarizer<double> height_1;
        SampleSummarizer<double> tree_length;

        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            REQUIRE(gtree->is_root());
            REQUIRE(gtree->get_number_of_children() == 2);
            REQUIRE(gtree->get_leaf_node_count() == 20);
            height_root.add_sample(gtree->get_height());
            double h0 = gtree->get_child(0)->get_height();
            double h1 = gtree->get_child(1)->get_height();
            if (h0 < h1) {
                height_0.add_sample(h1);
                height_1.add_sample(h0);
            }
            else {
                height_0.add_sample(h0);
                height_1.add_sample(h1);
            }
            tree_length.add_sample(gtree->get_clade_length());
        }

        REQUIRE(height_root.mean() == Approx(expected_mean_root).epsilon(0.0005));
        REQUIRE(height_root.variance() == Approx(expected_variance_root).epsilon(0.0005));

        REQUIRE(height_0.mean() == Approx(expected_mean_0).epsilon(0.00003));
        REQUIRE(height_1.mean() == Approx(expected_mean_1).epsilon(0.00003));


        tree.set_mutation_rate(mu * 2.0);

        SampleSummarizer<double> scale_height_root;
        SampleSummarizer<double> scale_height_0;
        SampleSummarizer<double> scale_height_1;
        SampleSummarizer<double> scale_tree_length;

        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            gtree->scale(0.5);
            REQUIRE(gtree->is_root());
            REQUIRE(gtree->get_number_of_children() == 2);
            REQUIRE(gtree->get_leaf_node_count() == 20);
            scale_height_root.add_sample(gtree->get_height());
            double h0 = gtree->get_child(0)->get_height();
            double h1 = gtree->get_child(1)->get_height();
            if (h0 < h1) {
                scale_height_0.add_sample(h1);
                scale_height_1.add_sample(h0);
            }
            else {
                scale_height_0.add_sample(h0);
                scale_height_1.add_sample(h1);
            }
            scale_tree_length.add_sample(gtree->get_clade_length());
        }

        REQUIRE(scale_height_root.mean() == Approx(expected_mean_root).epsilon(0.0005));
        REQUIRE(scale_height_root.variance() == Approx(expected_variance_root).epsilon(0.0005));

        REQUIRE(scale_height_0.mean() == Approx(expected_mean_0).epsilon(0.00003));
        REQUIRE(scale_height_1.mean() == Approx(expected_mean_1).epsilon(0.00003));

        REQUIRE(scale_tree_length.mean() == Approx(tree_length.mean()).epsilon(0.0001));
        REQUIRE(scale_tree_length.variance() == Approx(tree_length.variance()).epsilon(0.0001));
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree dataset simulation",
        "[ComparisonDirichletPopulationTree]") {
    SECTION("Testing simulate_biallelic_data_set for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, false, false, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.fix_population_size_multipliers();
        tree.set_mean_population_size(0.005);
        tree.fix_mean_population_size();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        BiallelicData data = tree.simulate_biallelic_data_set(rng);

        REQUIRE(data.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data.markers_are_dominant() == false);
        REQUIRE(tree.get_data().markers_are_dominant() == false);

        data.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8).epsilon(0.01));
        REQUIRE(data.get_proportion_1() == Approx(0.625).epsilon(0.01));
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.draw_from_prior for fully fixed",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing draw_from_prior for fully fixed pair") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.fix_population_size_multipliers();
        tree.set_mean_population_size(0.001);
        tree.fix_mean_population_size();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_root_population_size() == 0.001);
            REQUIRE(tree.get_child_population_size(0) == 0.001);
            REQUIRE(tree.get_child_population_size(1) == 0.001);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);
        }
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.draw_from_prior for constrained sizes",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing draw_from_prior for constrained sizes") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.fix_population_size_multipliers();
        tree.set_mean_population_size(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        SampleSummarizer<double> pop_size;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() == tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() == tree.get_child_population_size(1));

            pop_size.add_sample(tree.get_root_population_size());
        }

        REQUIRE(pop_size.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.draw_from_prior for unconstrained sizes",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing draw_from_prior for unconstrained sizes") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);
        std::vector<double> alphas {2.0, 1.0, 4.0};
        std::shared_ptr<DirichletDistribution> multiplier_prior = std::make_shared<DirichletDistribution>(alphas);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);
        tree.set_population_size_multiplier_prior(multiplier_prior);

        tree.set_root_height(0.1);
        tree.set_mean_population_size(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        SampleSummarizer<double> pop_size;
        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 100000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            pop_size.add_sample(tree.get_mean_population_size());
            std::vector<double> multipliers = tree.get_population_sizes_as_multipliers();
            pop_size_root.add_sample(multipliers.at(2));
            pop_size_0.add_sample(multipliers.at(0));
            pop_size_1.add_sample(multipliers.at(1));
        }

        REQUIRE(pop_size.mean() == Approx(size_prior->get_mean()).epsilon(0.01));
        REQUIRE(pop_size.variance() == Approx(size_prior->get_variance()).epsilon(0.01));
        REQUIRE(pop_size_root.mean() == Approx(multiplier_prior->get_mean(2) * 3.0).epsilon(0.01));
        REQUIRE(pop_size_root.variance() == Approx(multiplier_prior->get_variance(2) * 9.0).epsilon(0.01));
        REQUIRE(pop_size_0.mean() == Approx(multiplier_prior->get_mean(0) * 3.0).epsilon(0.01));
        REQUIRE(pop_size_0.variance() == Approx(multiplier_prior->get_variance(0) * 9.0).epsilon(0.01));
        REQUIRE(pop_size_1.mean() == Approx(multiplier_prior->get_mean(1) * 3.0).epsilon(0.01));
        REQUIRE(pop_size_1.variance() == Approx(multiplier_prior->get_variance(1) * 9.0).epsilon(0.01));
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.draw_from_prior for fully parameterized",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing draw_from_prior for fully parameterized") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(0.5, 0.8);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);
        std::vector<double> alphas {2.0, 1.0, 4.0};
        std::shared_ptr<DirichletDistribution> multiplier_prior = std::make_shared<DirichletDistribution>(alphas);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);
        tree.set_population_size_multiplier_prior(multiplier_prior);

        tree.set_root_height(0.1);
        tree.set_mean_population_size(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.estimate_state_frequencies();
        tree.set_freq_1(0.5);

        SampleSummarizer<double> pop_size;
        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;
        SampleSummarizer<double> rate;
        SampleSummarizer<double> f;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 100000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            pop_size.add_sample(tree.get_mean_population_size());
            std::vector<double> multipliers = tree.get_population_sizes_as_multipliers();
            pop_size_root.add_sample(multipliers.at(2));
            pop_size_0.add_sample(multipliers.at(0));
            pop_size_1.add_sample(multipliers.at(1));
            rate.add_sample(tree.get_mutation_rate());
            f.add_sample(tree.get_freq_1());
        }

        REQUIRE(pop_size.mean() == Approx(size_prior->get_mean()).epsilon(0.01));
        REQUIRE(pop_size.variance() == Approx(size_prior->get_variance()).epsilon(0.01));
        REQUIRE(pop_size_root.mean() == Approx(multiplier_prior->get_mean(2) * 3.0).epsilon(0.01));
        REQUIRE(pop_size_root.variance() == Approx(multiplier_prior->get_variance(2) * 9.0).epsilon(0.01));
        REQUIRE(pop_size_0.mean() == Approx(multiplier_prior->get_mean(0) * 3.0).epsilon(0.01));
        REQUIRE(pop_size_0.variance() == Approx(multiplier_prior->get_variance(0) * 9.0).epsilon(0.01));
        REQUIRE(pop_size_1.mean() == Approx(multiplier_prior->get_mean(1) * 3.0).epsilon(0.01));
        REQUIRE(pop_size_1.variance() == Approx(multiplier_prior->get_variance(1) * 9.0).epsilon(0.01));
        REQUIRE(f.mean() == Approx(f_prior->get_mean()).epsilon(0.01));
        REQUIRE(f.variance() == Approx(f_prior->get_variance()).epsilon(0.01));
        REQUIRE(rate.mean() == Approx(rate_prior->get_mean()).epsilon(0.1));
        REQUIRE(rate.variance() == Approx(rate_prior->get_variance()).epsilon(0.1));
    }
}



// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood = -31.77866581319647
TEST_CASE("Testing simple likelihood of RelativeRootPopulationTree", "[RelativeRootPopulationTree]") {

    SECTION("Testing constructor and likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing threaded likelihood of RelativeRootPopulationTree", "[RelativeRootPopulationTree]") {

    SECTION("Testing constructor and threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing over-threaded likelihood of RelativeRootPopulationTree", "[RelativeRootPopulationTree]") {

    SECTION("Testing constructor and over-threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(10);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(10);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood             = -248.93254688526213
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing hemi129.nex RelativeRootPopulationTree likelihood (0.01, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded RelativeRootPopulationTree likelihood (0.01, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7099.716015109998
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.01, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.01, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -6986.120524781545
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.01, 10.0, 1.0, 1.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6986.120524781545));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.01, 10.0, 1.0, 1.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(5);
        REQUIRE(l == Approx(-6986.120524781545));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}



// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -328.39238828878365
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex RelativeRootPopulationTree likelihood (0.0, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded RelativeRootPopulationTree likelihood (0.0, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(7);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(7);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7256.501742344454
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.0, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.0, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7223.362711937651
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.0, 10.0, 1.0, 1.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7223.362711937651));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.0, 10.0, 1.0, 1.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7223.362711937651));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -227.41048391087554
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex RelativeRootPopulationTree likelihood (0.2, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded RelativeRootPopulationTree likelihood (0.2, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7304.180743441677
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.2, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.2, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7405.145951634711
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.2, 10.0, 1.0, 1.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7405.145951634711));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.2, 10.0, 1.0, 1.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(6);
        REQUIRE(l == Approx(-7405.145951634711));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}




// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -327.7437811413033
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -6472.856486972301
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(5);
        REQUIRE(l == Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(3);
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -6494.774924871097
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0, 10.0/19.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6494.774924871097));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0, 10.0/19.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-6494.774924871097));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
  
  
// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -265.0023534261969
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(1);
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -10163.468886613919
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(5);
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -10999.288193543642
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0/19.0, 10.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-10999.288193543642));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0/19.0, 10.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-10999.288193543642));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}



// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -224.40177558289847
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex RelativeRootPopulationTree likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-224.40177558289847));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded RelativeRootPopulationTree likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(100000);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-224.40177558289847));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -8158.88094671241
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-8158.88094671241));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(1000000);
        REQUIRE(l == Approx(-8158.88094671241));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -8034.250341980543
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.03, 111.1, 10.0/19.0, 10.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-8034.250341980543));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.03, 111.1, 10.0/19.0, 10.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-8034.250341980543));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

TEST_CASE("Testing simple prior of RelativeRootPopulationTree", "[RelativeRootPopulationTree]") {

    SECTION("Testing constructor and prior calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);

        tree.set_root_height(0.1);
        tree.set_all_population_sizes(2.0/100.0);
        tree.set_freq_1(0.95);

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(100.0));
        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));
        tree.set_relative_root_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.1));
        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(2.0, 1.0));

        // height       -5.3948298140119091
        // root size    0.22402344985899036
        // leaf sizes   2 * -155.90663080917298
        // f1           0.64185388617239469 
        // total        -472.47286835535851
        
        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-311.58923816848699));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == Approx(0.64185388617239469));

        double d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-316.34221409632647));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-316.34221409632647));


        tree.store_prior_density();
        tree.constrain_state_frequencies();

        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-311.58923816848699));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == 0.0);

        d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-316.98406798249886));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-316.98406798249886));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-316.34221409632647));


        tree.store_prior_density();
        tree.constrain_population_sizes();

        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-155.90663080917298));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == 0.0);

        d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-161.30146062318488));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-161.30146062318488));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-316.98406798249886));

        tree.restore_prior_density();
        REQUIRE(tree.get_log_prior_density_value() == Approx(-316.98406798249886));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-316.98406798249886));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

TEST_CASE("Testing hemi129.nex state manipulation for RelativeRootPopulationTree", "[RelativeRootPopulationTree]") {

    SECTION("Testing state manipulation") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(10.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_relative_root_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.2));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(1.5, 3.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_height(0.01);
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1(0.5);
        REQUIRE(tree.is_dirty());

        REQUIRE(tree.get_number_of_likelihood_calculations() == 0);

        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-894.67638803083958));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 1);

        tree.store_state();
        tree.set_root_height(0.0);
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-894.57638803083955));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-894.67638803083958));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 2);

        tree.restore_state();
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-894.67638803083958));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-894.67638803083958));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.compute_log_likelihood_and_prior();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-894.67638803083958));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-894.67638803083958));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.store_state();
        tree.set_root_height(0.0);
        tree.compute_log_likelihood_and_prior();
        REQUIRE(tree.get_number_of_likelihood_calculations() == 4);


        tree.store_state();
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.57638803083955));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-894.57638803083955));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 5);


        tree.store_state();
        tree.set_root_height(0.03);
        tree.set_freq_1(0.05);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.05);
        REQUIRE(tree.get_freq_0() == Approx(0.95));
        REQUIRE(tree.get_u() == Approx(10.0));
        REQUIRE(tree.get_v() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-894.74397280499181));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.57638803083955));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 6);


        tree.store_state();
        tree.set_freq_1(0.95);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.95);
        REQUIRE(tree.get_freq_0() == Approx(0.05));
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-899.1606312737415));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-894.74397280499181));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 7);


        tree.store_state();
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(111.1 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-32.510853039559265));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-899.1606312737415));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 8);

        tree.store_state();
        tree.constrain_state_frequencies();
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.72489170735741));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-32.510853039559265));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 9);


        tree.store_state();
        tree.fold_patterns();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.72489170735741));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.72489170735741));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 10);


        tree.store_state();
        tree.constrain_population_sizes();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.35742912931147));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.72489170735741));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 11);

        tree.store_state();
        tree.estimate_mutation_rate();
        REQUIRE(! tree.mutation_rate_is_fixed());
        tree.set_mutation_rate(1.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());
        tree.set_mutation_rate_prior(std::make_shared<GammaDistribution>(10.0, 0.1));
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.35742912931147));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 12);

        tree.store_state();
        tree.set_mutation_rate(2.0);
        tree.set_root_height(0.1);
        tree.set_all_population_sizes(2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 2.0);
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(20.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-206.13340567945255));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        tree.restore_state();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 9fb5b3b7817a0bd4a21e3f90132f132cca72ce4e)
// SNAPP v1.3.0 (master 4f3f0f7366798f4fb38b766c15f6426a75ddf71e)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0/19.0
// Log likelihood            = -23.81984255023975
// Log likelihood correction = -6.87935580446044
//
// With constant sites inclucded and m_bUseNonPolymorphic = true
// Log likelihood            = -55.01646493341547
// Log likelihood correction = -6.87935580446044
TEST_CASE("Testing affect of constant sites on likelihood of RelativeRootPopulationTree", "[RelativeRootPopulationTree]") {

    SECTION("Testing haploid-standard-full-constant.nex") {
        std::string nex_path = "data/haploid-standard-full-constant.nex";
        RelativeRootPopulationTree t_included(
                nex_path, // path
                ' ',      // pop name delimiter
                true,     // pop name is prefix
                false,    // genotypes are diploid
                false,    // markers are dominant
                false,    // constant sites removed
                true);    // validate
        std::string nex_path2 = "data/haploid-standard-full-constant-removed.nex";
        RelativeRootPopulationTree t(
                nex_path2, // path
                ' ',       // pop name delimiter
                true,      // pop name is prefix
                false,     // genotypes are diploid
                false,     // markers are dominant
                true,      // constant sites removed
                true);     // validate

        t.set_root_height(0.03);
        t.set_freq_1(0.05);
        t.set_all_population_sizes(2.0/(10.0 * 2 * t.get_ploidy()));

        t_included.set_root_height(0.03);
        t_included.set_freq_1(0.05);
        t_included.set_all_population_sizes(2.0/(10.0 * 2 * t_included.get_ploidy()));

        double l = t.compute_log_likelihood();
        REQUIRE(l == Approx(-23.81984255023975));
        REQUIRE(t.get_likelihood_correction() == Approx(-6.87935580446044));

        double l_included = t_included.compute_log_likelihood();

        REQUIRE(l_included == Approx(-55.01646493341547));

        REQUIRE(t.get_likelihood_correction() == t_included.get_likelihood_correction());
    }
}

TEST_CASE("Testing affect of constant sites on threaded likelihood of RelativeRootPopulationTree", "[RelativeRootPopulationTree]") {

    SECTION("Testing haploid-standard-full-constant.nex") {
        std::string nex_path = "data/haploid-standard-full-constant.nex";
        RelativeRootPopulationTree t_included(
                nex_path, // path
                ' ',      // pop name delimiter
                true,     // pop name is prefix
                false,    // genotypes are diploid
                false,    // markers are dominant
                false,    // constant sites removed
                true);    // validate
        std::string nex_path2 = "data/haploid-standard-full-constant-removed.nex";
        RelativeRootPopulationTree t(
                nex_path2, // path
                ' ',       // pop name delimiter
                true,      // pop name is prefix
                false,     // genotypes are diploid
                false,     // markers are dominant
                true,      // constant sites removed
                true);     // validate

        t.set_root_height(0.03);
        t.set_freq_1(0.05);
        t.set_all_population_sizes(2.0/(10.0 * 2 * t.get_ploidy()));

        t_included.set_root_height(0.03);
        t_included.set_freq_1(0.05);
        t_included.set_all_population_sizes(2.0/(10.0 * 2 * t_included.get_ploidy()));

        double l = t.compute_log_likelihood(4);
        REQUIRE(l == Approx(-23.81984255023975));
        REQUIRE(t.get_likelihood_correction() == Approx(-6.87935580446044));

        double l_included = t_included.compute_log_likelihood(2);

        REQUIRE(l_included == Approx(-55.01646493341547));

        REQUIRE(t.get_likelihood_correction() == t_included.get_likelihood_correction());
    }
}


// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood = -31.77866581319647
TEST_CASE("Testing simple likelihood of ComparisonRelativeRootPopulationTree", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing constructor and likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_root().get_label() == "root-pop1");
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing simple threaded likelihood of ComparisonRelativeRootPopulationTree", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing constructor and threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_root().get_label() == "root-pop1");
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex state manipulation", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing state manipulation") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(10.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_relative_root_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.2));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(1.5, 3.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_height(0.01);
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_child_population_size(0, 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_child_population_size(1, 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));

        tree.set_freq_1(0.5);
        REQUIRE(tree.is_dirty());

        REQUIRE(tree.get_number_of_likelihood_calculations() == 0);

        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 1);

        tree.store_state();
        tree.set_root_height(0.0);
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 2);

        tree.restore_state();
        // ComparisonRelativeRootPopulationTree does not store/restore heights
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.compute_log_likelihood_and_prior();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.store_state();
        tree.set_root_height(0.0);
        tree.compute_log_likelihood_and_prior();
        REQUIRE(tree.get_number_of_likelihood_calculations() == 4);


        tree.store_state();
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 5);


        tree.store_state();
        tree.set_root_height(0.03);
        tree.set_freq_1(0.05);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.05);
        REQUIRE(tree.get_freq_0() == 0.95);
        REQUIRE(tree.get_u() == Approx(10.0));
        REQUIRE(tree.get_v() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.74655789798578));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 6);


        tree.store_state();
        tree.set_freq_1(0.95);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_freq_1() == 0.95);
        REQUIRE(tree.get_freq_0() == Approx(0.05));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-901.16321636673547));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.74655789798578));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 7);


        tree.store_state();
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(111.1 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(111.1 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(111.1 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-34.513438132553311));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-901.16321636673547));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 8);

        tree.store_state();
        tree.constrain_state_frequencies();
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-897.02747680035145));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-34.513438132553311));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 9);


        tree.store_state();
        tree.fold_patterns();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-897.02747680035145));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-897.02747680035145));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 10);


        tree.store_state();
        tree.constrain_population_sizes();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.66001422230551));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-897.02747680035145));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 11);

        tree.store_state();
        tree.estimate_mutation_rate();
        REQUIRE(! tree.mutation_rate_is_fixed());
        tree.set_mutation_rate(1.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());
        tree.set_mutation_rate_prior(std::make_shared<GammaDistribution>(10.0, 0.1));
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.66001422230551));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 12);

        tree.store_state();
        tree.set_mutation_rate(2.0);
        tree.set_root_height(0.1);
        tree.set_root_population_size(2.0/(20.0 * 2 * tree.get_ploidy()));
        tree.set_child_population_size(0, 2.0/(20.0 * 2 * tree.get_ploidy()));
        tree.set_child_population_size(1, 2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 2.0);
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(20.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(20.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(20.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-207.43599077244659));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        tree.restore_state();
        // ComparisonRelativeRootPopulationTree does not store/restore node heights
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.01
// coalescent_rate = 100.0, 200.0, 500.0
// u = 1.0
// v = 1.0
// Log likelihood            = -226.11914854623677
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex likelihood (0.01, 100.0, 200.0, 500.0)", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_child_population_size(0, 2.0/(100.0 * 2 * tree.get_ploidy()));
        tree.set_child_population_size(1, 2.0/(200.0 * 2 * tree.get_ploidy()));
        tree.set_root_population_size(2.0/(500.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-226.11914854623677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}


// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.00506843962151613554
// coalescent_rate = 2.0 / 0.00018955324120485613
// u = 1.0
// v = 1.0
// Log likelihood            = -277.06960543551577
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex likelihood (0.00506843962151613554, 2.0 / 0.00018955324120485613)", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.00506843962151613554);
        tree.set_all_population_sizes(0.00018955324120485613 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-277.06960543551577));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 9.08323190033687971e-09
// coalescent_rate = 2.0 / 2.47975039926886321e-08
// u = 1.0
// v = 1.0
// Log likelihood            = -221.69627648370943
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex likelihood (9.08323190033687971e-09, 2.0 / 2.47975039926886321e-08)", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(9.08323190033687971e-09);
        tree.set_all_population_sizes(2.47975039926886321e-08 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-221.69627648370943));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.04921319733994759e-08
// coalescent_rate = 2.0 / 2.75977168733651178e-10
// u = 1.0
// v = 1.0
// Log likelihood            = -324.2737564069293
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex likelihood (1.04921319733994759e-08, 2.0 / 2.75977168733651178e-10)", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.04921319733994759e-08);
        tree.set_all_population_sizes(2.75977168733651178e-10 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-324.2737564069293));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.012386610001351e-08
// coalescent_rate = 2.0 / 5.75048645855884647e-30
// u = 1.0
// v = 1.0
// Log likelihood            = -1364.1427530000253
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex likelihood (1.012386610001351e-08, 2.0 / 5.75048645855884647e-30)", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.012386610001351e-08);
        tree.set_all_population_sizes(5.75048645855884647e-30 / 4.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-1364.1427530000253));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.04856228318474786e-08
// coalescent_rate = 2.0 / 4.43934332792563837e-305
// u = 1.0
// v = 1.0
// Log likelihood            = NaN
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex likelihood (1.04856228318474786e-08, 2.0 / 4.43934332792563837e-305)", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.04856228318474786e-08);
        tree.set_all_population_sizes(4.43934332792563837e-305 / 4.0);
        double l = tree.compute_log_likelihood();
        // REQUIRE(std::isnan(l));
        REQUIRE(std::isinf(l));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.012386610001351e-08
// coalescent_rate = 2.0 / 5.46641122085615013e-09,
//                   2.0 / 7.39871781998828579e-08,
//                   2.0 / 2.71077053326002069e-13
// u = 1.0
// v = 1.0
// Log likelihood            = -96.34394008351177
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex weirdness", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.012386610001351e-08);
        tree.set_child_population_size(0, 5.46641122085615013e-09 / 4.0);
        tree.set_child_population_size(1, 7.39871781998828579e-08 / 4.0);
        tree.set_root_population_size(2.71077053326002069e-13 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-96.34394008351177));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.036374107244057e-08
// coalescent_rate = 2.0 / 4.57999694763258361e-09,
//                   2.0 / 6.70991782555376588e-08,
//                   2.0 / 1.33514111020266258e-08
// u = 1.0
// v = 1.0
// Log likelihood            = -44.95791900747736
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex weirdness 2", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.036374107244057e-08);
        tree.set_child_population_size(0, 4.57999694763258361e-09 / 4.0);
        tree.set_child_population_size(1, 6.70991782555376588e-08 / 4.0);
        tree.set_root_population_size(1.33514111020266258e-08 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-44.95791900747736));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 2 lineages and theta of 1.0",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 2 lineages and theta of 3.7",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 3.7;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 2 lineages and theta of 0.17",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 3 lineages and theta of 1.0",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 3 lineages and theta of 1.47",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 1.47;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 3 lineages and theta of 0.17",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 10 lineages and theta of 1.0",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 10 lineages and theta of 1.47",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 1.47;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 10 lineages and theta of 0.17",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree scaling of simulate_gene_tree for singleton",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        std::string nex_path = "data/singleton-n2.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.set_child_population_size(0, (Ne * mu));

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        unsigned int n;
        double mean;
        double variance;
        double sum_devs;
        double d;
        double d_n;
        double mn;
        double mx;
        double x;
        std::shared_ptr<GeneTreeSimNode> gtree;

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            x = gtree->get_height();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);

        REQUIRE(gtree->is_root());
        REQUIRE(gtree->get_number_of_children() == 2);
        REQUIRE(gtree->get_leaf_node_count() == 2);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.0005));
        REQUIRE(mn >= 0.0);

        tree.set_mutation_rate(mu);
        tree.set_root_population_size(Ne);
        tree.set_child_population_size(0, Ne);
        tree.set_root_height(0.1 / mu);

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            x = gtree->get_height();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);

        REQUIRE(gtree->is_root());
        REQUIRE(gtree->get_number_of_children() == 2);
        REQUIRE(gtree->get_leaf_node_count() == 2);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree scaling of simulate_gene_tree for pair",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing pair") {
        unsigned int Ne_root = 100000;
        unsigned int Ne_0 = 200000;
        unsigned int Ne_1 = 10000;
        double mu = 1e-8;
        double theta_root = 4 * Ne_root * mu;
        double theta_0 = 4 * Ne_0 * mu;
        double theta_1 = 4 * Ne_1 * mu;
        double time = 2e7;

        double expected_mean_root = theta_root * (1.0 - (1.0 / 2));
        double expected_variance_root = expected_mean_root * expected_mean_root;
        expected_mean_root += (time * mu);

        double expected_mean_0 = theta_0 * (1.0 - (1.0 / 10));
        double expected_mean_1 = theta_1 * (1.0 - (1.0 / 10));

        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        tree.estimate_mutation_rate();

        tree.set_child_population_size(0, (Ne_0));
        tree.set_child_population_size(1, (Ne_1));
        tree.set_root_population_size(Ne_root);

        tree.set_root_height(time);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(mu);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        std::shared_ptr<GeneTreeSimNode> gtree;
        SampleSummarizer<double> height_root;
        SampleSummarizer<double> height_0;
        SampleSummarizer<double> height_1;
        SampleSummarizer<double> tree_length;

        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            REQUIRE(gtree->is_root());
            REQUIRE(gtree->get_number_of_children() == 2);
            REQUIRE(gtree->get_leaf_node_count() == 20);
            height_root.add_sample(gtree->get_height());
            double h0 = gtree->get_child(0)->get_height();
            double h1 = gtree->get_child(1)->get_height();
            if (h0 < h1) {
                height_0.add_sample(h1);
                height_1.add_sample(h0);
            }
            else {
                height_0.add_sample(h0);
                height_1.add_sample(h1);
            }
            tree_length.add_sample(gtree->get_clade_length());
        }

        REQUIRE(height_root.mean() == Approx(expected_mean_root).epsilon(0.0005));
        REQUIRE(height_root.variance() == Approx(expected_variance_root).epsilon(0.0005));

        REQUIRE(height_0.mean() == Approx(expected_mean_0).epsilon(0.00003));
        REQUIRE(height_1.mean() == Approx(expected_mean_1).epsilon(0.00003));


        tree.set_mutation_rate(mu * 2.0);

        SampleSummarizer<double> scale_height_root;
        SampleSummarizer<double> scale_height_0;
        SampleSummarizer<double> scale_height_1;
        SampleSummarizer<double> scale_tree_length;

        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            gtree->scale(0.5);
            REQUIRE(gtree->is_root());
            REQUIRE(gtree->get_number_of_children() == 2);
            REQUIRE(gtree->get_leaf_node_count() == 20);
            scale_height_root.add_sample(gtree->get_height());
            double h0 = gtree->get_child(0)->get_height();
            double h1 = gtree->get_child(1)->get_height();
            if (h0 < h1) {
                scale_height_0.add_sample(h1);
                scale_height_1.add_sample(h0);
            }
            else {
                scale_height_0.add_sample(h0);
                scale_height_1.add_sample(h1);
            }
            scale_tree_length.add_sample(gtree->get_clade_length());
        }

        REQUIRE(scale_height_root.mean() == Approx(expected_mean_root).epsilon(0.0005));
        REQUIRE(scale_height_root.variance() == Approx(expected_variance_root).epsilon(0.0005));

        REQUIRE(scale_height_0.mean() == Approx(expected_mean_0).epsilon(0.00003));
        REQUIRE(scale_height_1.mean() == Approx(expected_mean_1).epsilon(0.00003));

        REQUIRE(scale_tree_length.mean() == Approx(tree_length.mean()).epsilon(0.0001));
        REQUIRE(scale_tree_length.variance() == Approx(tree_length.variance()).epsilon(0.0001));
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree dataset simulation", "[ComparisonRelativeRootPopulationTree]") {
    SECTION("Testing simulate_biallelic_data_set for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, false, false, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        BiallelicData data = tree.simulate_biallelic_data_set(rng);

        REQUIRE(data.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data.markers_are_dominant() == false);
        REQUIRE(tree.get_data().markers_are_dominant() == false);

        data.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8).epsilon(0.01));
        REQUIRE(data.get_proportion_1() == Approx(0.625).epsilon(0.01));
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree draw_from_prior for fully fixed", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fully fixed pair") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.001);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_root_population_size() == 0.001);
            REQUIRE(tree.get_child_population_size(0) == 0.001);
            REQUIRE(tree.get_child_population_size(1) == 0.001);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);
        }
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree draw_from_prior for constrained sizes", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for constrained sizes") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        SampleSummarizer<double> pop_size;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() == tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() == tree.get_child_population_size(1));

            pop_size.add_sample(tree.get_root_population_size());
        }

        REQUIRE(pop_size.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree draw_from_prior for unconstrained sizes", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for unconstrained sizes") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);
        std::shared_ptr<GammaDistribution> rel_root_size_prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);
        tree.set_relative_root_population_size_prior(rel_root_size_prior);

        tree.set_root_height(0.1);
        tree.set_all_population_sizes(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            pop_size_root.add_sample(tree.get_relative_root_population_size());
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
        }

        REQUIRE(pop_size_root.mean() == Approx(rel_root_size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_root.variance() == Approx(rel_root_size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_0.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_0.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_1.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_1.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree draw_from_prior for fully parameterized", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fully parameterized") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(0.5, 0.8);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);
        std::shared_ptr<GammaDistribution> rel_root_size_prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);
        tree.set_relative_root_population_size_prior(rel_root_size_prior);

        tree.set_root_height(0.1);
        tree.set_all_population_sizes(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.estimate_state_frequencies();
        tree.set_freq_1(0.5);

        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;
        SampleSummarizer<double> rate;
        SampleSummarizer<double> f;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            pop_size_root.add_sample(tree.get_relative_root_population_size());
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
            rate.add_sample(tree.get_mutation_rate());
            f.add_sample(tree.get_freq_1());
        }

        REQUIRE(pop_size_root.mean() == Approx(rel_root_size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_root.variance() == Approx(rel_root_size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_0.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_0.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_1.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_1.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(f.mean() == Approx(f_prior->get_mean()).epsilon(0.01));
        REQUIRE(f.variance() == Approx(f_prior->get_variance()).epsilon(0.01));
        REQUIRE(rate.mean() == Approx(rate_prior->get_mean()).epsilon(0.1));
        REQUIRE(rate.variance() == Approx(rate_prior->get_variance()).epsilon(0.1));
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree draw_from_prior for fixed relative root size",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fixed relative root size") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);
        std::shared_ptr<GammaDistribution> rel_root_size_prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);
        tree.set_relative_root_population_size_prior(rel_root_size_prior);

        tree.set_root_height(0.1);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();
        tree.set_all_population_sizes(0.001);
        tree.set_root_population_size(0.002);
        tree.fix_relative_root_population_size();

        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            REQUIRE(tree.get_relative_root_population_size() == Approx(2.0));
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
        }

        REQUIRE(pop_size_0.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_0.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_1.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_1.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree draw_from_prior for fixed leaf sizes",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fixed leaf sizes") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);
        std::shared_ptr<GammaDistribution> rel_root_size_prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);
        tree.set_relative_root_population_size_prior(rel_root_size_prior);

        tree.set_root_height(0.1);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();
        tree.set_all_population_sizes(0.001);
        tree.fix_population_sizes();

        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) == tree.get_child_population_size(1));

            pop_size_root.add_sample(tree.get_relative_root_population_size());
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
        }

        REQUIRE(pop_size_root.mean() == Approx(rel_root_size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_root.variance() == Approx(rel_root_size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_0.mean() == Approx(0.001));
        REQUIRE(pop_size_0.variance() == Approx(0.0));
        REQUIRE(pop_size_1.mean() == Approx(0.001));
        REQUIRE(pop_size_1.variance() == Approx(0.0));
    }
}

TEST_CASE("Testing config ComparisonRelativeRootPopulationTree draw_from_prior for fully fixed", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fully fixed pair") {

        std::stringstream os;
        os << "path: data/diploid-standard-data-ntax5-nchar5.nex\n";
        os << "genotypes_are_diploid: true\n";
        os << "markers_are_dominant: false\n";
        os << "population_name_delimiter: \" \"\n";
        os << "population_name_is_prefix: true\n";
        os << "constant_sites_removed: true\n";
        os << "equal_population_sizes: false\n";
        os << "parameters:\n";
        os << "    mutation_rate:\n";
        os << "        value: 0.8\n";
        os << "        estimate: false\n";
        os << "    population_size:\n";
        os << "        value: 0.001\n";
        os << "        estimate: false\n";
        os << "    root_relative_population_size:\n";
        os << "        value: 1.0\n";
        os << "        estimate: false\n";
        os << "    freq_1:\n";
        os << "        value: 0.5\n";
        os << "        estimate: false\n";

        YAML::Node n;
        n = YAML::Load(os);
        RelativeRootComparisonSettings settings = RelativeRootComparisonSettings(
                n,
                "dummy-config-path.yml");

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        ComparisonRelativeRootPopulationTree tree = ComparisonRelativeRootPopulationTree(settings, rng);

        for (unsigned int i = 0; i < 10; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_population_size() == 0.001);
            REQUIRE(tree.get_child_population_size(0) == 0.001);
            REQUIRE(tree.get_child_population_size(1) == 0.001);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);
        }
    }
}

TEST_CASE("Testing config ComparisonRelativeRootPopulationTree draw_from_prior for fully fixed and constrained", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fully fixed pair") {

        std::stringstream os;
        os << "path: data/diploid-standard-data-ntax5-nchar5.nex\n";
        os << "genotypes_are_diploid: true\n";
        os << "markers_are_dominant: false\n";
        os << "population_name_delimiter: \" \"\n";
        os << "population_name_is_prefix: true\n";
        os << "constant_sites_removed: true\n";
        os << "equal_population_sizes: true\n";
        os << "parameters:\n";
        os << "    mutation_rate:\n";
        os << "        value: 0.8\n";
        os << "        estimate: false\n";
        os << "    population_size:\n";
        os << "        value: 0.001\n";
        os << "        estimate: false\n";
        os << "    root_relative_population_size:\n";
        os << "        value: 1.0\n";
        os << "        estimate: false\n";
        os << "    freq_1:\n";
        os << "        value: 0.5\n";
        os << "        estimate: false\n";

        YAML::Node n;
        n = YAML::Load(os);
        RelativeRootComparisonSettings settings = RelativeRootComparisonSettings(
                n,
                "dummy-config-path.yml");

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        ComparisonRelativeRootPopulationTree tree = ComparisonRelativeRootPopulationTree(settings, rng);

        for (unsigned int i = 0; i < 10; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_population_size() == 0.001);
            REQUIRE(tree.get_child_population_size(0) == 0.001);
            REQUIRE(tree.get_child_population_size(1) == 0.001);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);
        }
    }
}

TEST_CASE("Testing config ComparisonRelativeRootPopulationTree draw_from_prior for constrained sizes", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for constrained sizes") {

        double size_shape = 2.0;
        double size_scale = 1.2;

        std::stringstream os;
        os << "path: data/diploid-standard-data-ntax5-nchar5.nex\n";
        os << "genotypes_are_diploid: true\n";
        os << "markers_are_dominant: false\n";
        os << "population_name_delimiter: \" \"\n";
        os << "population_name_is_prefix: true\n";
        os << "constant_sites_removed: true\n";
        os << "equal_population_sizes: true\n";
        os << "parameters:\n";
        os << "    mutation_rate:\n";
        os << "        value: 0.8\n";
        os << "        estimate: false\n";
        os << "    population_size:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << size_shape << "\n";
        os << "                scale: " << size_scale << "\n";
        os << "    root_relative_population_size:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: 20.0\n";
        os << "                scale: 0.05\n";
        os << "    freq_1:\n";
        os << "        value: 0.5\n";
        os << "        estimate: false\n";

        YAML::Node n;
        n = YAML::Load(os);
        RelativeRootComparisonSettings settings = RelativeRootComparisonSettings(
                n,
                "dummy-config-path.yml");

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        ComparisonRelativeRootPopulationTree tree = ComparisonRelativeRootPopulationTree(settings, rng);

        SampleSummarizer<double> pop_size;

        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_relative_root_population_size() == Approx(1.0));
            REQUIRE(tree.get_root_population_size() == tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() == tree.get_child_population_size(1));

            pop_size.add_sample(tree.get_root_population_size());
        }

        REQUIRE(pop_size.mean() == Approx(size_shape * size_scale).epsilon(0.1));
        REQUIRE(pop_size.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.1));
    }
}

TEST_CASE("Testing config ComparisonRelativeRootPopulationTree draw_from_prior for unconstrained sizes", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for unconstrained sizes") {
        double size_shape = 2.0;
        double size_scale = 1.2;
        double rel_size_shape = 10.0;
        double rel_size_scale = 0.1;

        std::stringstream os;
        os << "path: data/diploid-standard-data-ntax5-nchar5.nex\n";
        os << "genotypes_are_diploid: true\n";
        os << "markers_are_dominant: false\n";
        os << "population_name_delimiter: \" \"\n";
        os << "population_name_is_prefix: true\n";
        os << "constant_sites_removed: true\n";
        os << "equal_population_sizes: false\n";
        os << "parameters:\n";
        os << "    mutation_rate:\n";
        os << "        value: 0.8\n";
        os << "        estimate: false\n";
        os << "    population_size:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << size_shape << "\n";
        os << "                scale: " << size_scale << "\n";
        os << "    root_relative_population_size:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << rel_size_shape << "\n";
        os << "                scale: " << rel_size_scale << "\n";
        os << "    freq_1:\n";
        os << "        value: 0.5\n";
        os << "        estimate: false\n";

        YAML::Node n;
        n = YAML::Load(os);
        RelativeRootComparisonSettings settings = RelativeRootComparisonSettings(
                n,
                "dummy-config-path.yml");

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        ComparisonRelativeRootPopulationTree tree = ComparisonRelativeRootPopulationTree(settings, rng);

        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;

        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            pop_size_root.add_sample(tree.get_relative_root_population_size());
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
        }

        REQUIRE(pop_size_root.mean() == Approx(rel_size_shape * rel_size_scale).epsilon(0.1));
        REQUIRE(pop_size_root.variance() == Approx(rel_size_shape * rel_size_scale * rel_size_scale).epsilon(0.1));
        REQUIRE(pop_size_0.mean() == Approx(size_shape * size_scale).epsilon(0.1));
        REQUIRE(pop_size_0.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.1));
        REQUIRE(pop_size_1.mean() == Approx(size_shape * size_scale).epsilon(0.1));
        REQUIRE(pop_size_1.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.1));
    }
}

TEST_CASE("Testing config ComparisonRelativeRootPopulationTree draw_from_prior for fully parameterized", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fully parameterized") {
        double size_shape = 2.0;
        double size_scale = 1.2;
        double rel_size_shape = 10.0;
        double rel_size_scale = 0.1;
        double mu_shape = 3.0;
        double mu_scale = 1.1;
        double f_a = 0.5;
        double f_b = 0.8;

        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(f_a, f_b);

        std::stringstream os;
        os << "path: data/diploid-standard-data-ntax5-nchar5.nex\n";
        os << "genotypes_are_diploid: true\n";
        os << "markers_are_dominant: false\n";
        os << "population_name_delimiter: \" \"\n";
        os << "population_name_is_prefix: true\n";
        os << "constant_sites_removed: true\n";
        os << "equal_population_sizes: false\n";
        os << "parameters:\n";
        os << "    mutation_rate:\n";
        os << "        value: 0.8\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << mu_shape << "\n";
        os << "                scale: " << mu_scale << "\n";
        os << "    population_size:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << size_shape << "\n";
        os << "                scale: " << size_scale << "\n";
        os << "    root_relative_population_size:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << rel_size_shape << "\n";
        os << "                scale: " << rel_size_scale << "\n";
        os << "    freq_1:\n";
        os << "        value: 0.5\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            beta_distribution:\n";
        os << "                alpha: " << f_a << "\n";
        os << "                beta: " << f_b << "\n";

        YAML::Node n;
        n = YAML::Load(os);
        RelativeRootComparisonSettings settings = RelativeRootComparisonSettings(
                n,
                "dummy-config-path.yml");

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        ComparisonRelativeRootPopulationTree tree = ComparisonRelativeRootPopulationTree(settings, rng);

        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;
        SampleSummarizer<double> rate;
        SampleSummarizer<double> f;

        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            pop_size_root.add_sample(tree.get_relative_root_population_size());
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
            rate.add_sample(tree.get_mutation_rate());
            f.add_sample(tree.get_freq_1());
        }

        REQUIRE(pop_size_root.mean() == Approx(rel_size_shape * rel_size_scale).epsilon(0.1));
        REQUIRE(pop_size_root.variance() == Approx(rel_size_shape * rel_size_scale * rel_size_scale).epsilon(0.1));
        REQUIRE(pop_size_0.mean() == Approx(size_shape * size_scale).epsilon(0.1));
        REQUIRE(pop_size_0.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.1));
        REQUIRE(pop_size_1.mean() == Approx(size_shape * size_scale).epsilon(0.1));
        REQUIRE(pop_size_1.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.1));
        REQUIRE(f.mean() == Approx(f_prior->get_mean()).epsilon(0.01));
        REQUIRE(f.variance() == Approx(f_prior->get_variance()).epsilon(0.01));
        REQUIRE(rate.mean() == Approx(mu_shape * mu_scale).epsilon(0.1));
        REQUIRE(rate.variance() == Approx(mu_shape * mu_scale * mu_scale).epsilon(0.1));
    }
}

TEST_CASE("Testing config ComparisonRelativeRootPopulationTree draw_from_prior for fixed relative root size",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fixed relative root size") {
        double size_shape = 2.0;
        double size_scale = 1.2;

        std::stringstream os;
        os << "path: data/diploid-standard-data-ntax5-nchar5.nex\n";
        os << "genotypes_are_diploid: true\n";
        os << "markers_are_dominant: false\n";
        os << "population_name_delimiter: \" \"\n";
        os << "population_name_is_prefix: true\n";
        os << "constant_sites_removed: true\n";
        os << "equal_population_sizes: false\n";
        os << "parameters:\n";
        os << "    mutation_rate:\n";
        os << "        value: 0.8\n";
        os << "        estimate: false\n";
        os << "    population_size:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << size_shape << "\n";
        os << "                scale: " << size_scale << "\n";
        os << "    root_relative_population_size:\n";
        os << "        value: 2.0\n";
        os << "        estimate: false\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: 10.0\n";
        os << "                scale: 0.1\n";
        os << "    freq_1:\n";
        os << "        value: 0.5\n";
        os << "        estimate: false\n";

        YAML::Node n;
        n = YAML::Load(os);
        RelativeRootComparisonSettings settings = RelativeRootComparisonSettings(
                n,
                "dummy-config-path.yml");

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        ComparisonRelativeRootPopulationTree tree = ComparisonRelativeRootPopulationTree(settings, rng);

        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;

        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            REQUIRE(tree.get_relative_root_population_size() == Approx(2.0));
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
        }

        REQUIRE(pop_size_0.mean() == Approx(size_shape * size_scale).epsilon(0.1));
        REQUIRE(pop_size_0.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.1));
        REQUIRE(pop_size_1.mean() == Approx(size_shape * size_scale).epsilon(0.1));
        REQUIRE(pop_size_1.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.1));
    }
}

TEST_CASE("Testing config ComparisonRelativeRootPopulationTree draw_from_prior for fixed leaf sizes",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fixed leaf sizes") {
        double rel_size_shape = 10.0;
        double rel_size_scale = 0.1;

        std::stringstream os;
        os << "path: data/diploid-standard-data-ntax5-nchar5.nex\n";
        os << "genotypes_are_diploid: true\n";
        os << "markers_are_dominant: false\n";
        os << "population_name_delimiter: \" \"\n";
        os << "population_name_is_prefix: true\n";
        os << "constant_sites_removed: true\n";
        os << "equal_population_sizes: false\n";
        os << "parameters:\n";
        os << "    mutation_rate:\n";
        os << "        value: 0.8\n";
        os << "        estimate: false\n";
        os << "    population_size:\n";
        os << "        value: 0.001\n";
        os << "        estimate: false\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: 2.0\n";
        os << "                scale: 1.2\n";
        os << "    root_relative_population_size:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << rel_size_shape << "\n";
        os << "                scale: " << rel_size_scale << "\n";
        os << "    freq_1:\n";
        os << "        value: 0.5\n";
        os << "        estimate: false\n";

        YAML::Node n;
        n = YAML::Load(os);
        RelativeRootComparisonSettings settings = RelativeRootComparisonSettings(
                n,
                "dummy-config-path.yml");

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        ComparisonRelativeRootPopulationTree tree = ComparisonRelativeRootPopulationTree(settings, rng);


        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;

        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) == tree.get_child_population_size(1));

            pop_size_root.add_sample(tree.get_relative_root_population_size());
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
        }

        REQUIRE(pop_size_root.mean() == Approx(rel_size_shape * rel_size_scale).epsilon(0.1));
        REQUIRE(pop_size_root.variance() == Approx(rel_size_shape * rel_size_scale * rel_size_scale).epsilon(0.1));
        REQUIRE(pop_size_0.mean() == Approx(0.001));
        REQUIRE(pop_size_0.variance() == Approx(0.0));
        REQUIRE(pop_size_1.mean() == Approx(0.001));
        REQUIRE(pop_size_1.variance() == Approx(0.0));
    }
}



TEST_CASE("Testing missing gene copies",
        "[ComparisonPopulationTree]") {

    SECTION("Testing missing gene copies for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        std::string nex_path2 = "data/singleton-n2.nex";
        std::string nex_path3 = "data/singleton-n3.nex";
        // Need to keep constant characters
        ComparisonPopulationTree tree2(nex_path2, ' ', true, false, false, false);
        ComparisonPopulationTree tree3(nex_path3, ' ', true, false, false, false);
        tree2.estimate_mutation_rate();
        tree3.estimate_mutation_rate();

        tree2.set_root_population_size(Ne * mu);
        tree3.set_root_population_size(Ne * mu);

        tree2.set_child_population_size(0, (Ne * mu));
        tree3.set_child_population_size(0, (Ne * mu));

        tree2.set_root_height(0.1);
        tree3.set_root_height(0.1);

        tree2.set_freq_1(0.5);
        tree3.set_freq_1(0.5);

        tree2.set_mutation_rate(1.0);
        tree3.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        unsigned int nsamples = 1000000;

        std::shared_ptr<GeneTreeSimNode> gtree2;
        std::shared_ptr<GeneTreeSimNode> gtree3;
        std::vector<unsigned int> down_sampled_counts = {2};

        std::vector<unsigned int> red_allele_counts2 = {0, 0, 0};
        std::vector<unsigned int> red_allele_counts3 = {0, 0, 0};

        for (unsigned int i = 0; i < nsamples; ++i) {
            gtree2 = tree2.simulate_gene_tree(0, rng, true);
            gtree3 = tree3.simulate_gene_tree(0, rng, true);

            auto pattern2 = tree2.simulate_biallelic_site(gtree2, rng);
            auto pattern3 = tree3.simulate_biallelic_site_sans_missing(gtree3, down_sampled_counts, rng);

            ++red_allele_counts2.at(pattern2.first.at(0));
            ++red_allele_counts3.at(pattern3.first.at(0));
        }

        unsigned int total2 = 0;
        unsigned int total3 = 0;
        for (unsigned int i = 0; i < red_allele_counts2.size(); ++i) {
            total2 += red_allele_counts2.at(i);
            total3 += red_allele_counts3.at(i);
        }

        REQUIRE(total2 == nsamples);
        REQUIRE(total3 == nsamples);

        std::vector<double> red_allele_freqs2 = {0.0, 0.0, 0.0};
        std::vector<double> red_allele_freqs3 = {0.0, 0.0, 0.0};
        for (unsigned int i = 0; i < red_allele_counts2.size(); ++i) {
            red_allele_freqs2.at(i) = red_allele_counts2.at(i) / (double)nsamples;
            red_allele_freqs3.at(i) = red_allele_counts3.at(i) / (double)nsamples;
        }

        std::cout << "Freqs of red allele counts for 2 gene copies:\n";
        for (unsigned int i = 0; i < red_allele_freqs2.size(); ++i) {
            std::cout << i << ": " << red_allele_freqs2.at(i) << "\n";
        }
        std::cout << "Freqs of red allele counts for 3 gene copies with 1 dropped:\n";
        for (unsigned int i = 0; i < red_allele_freqs3.size(); ++i) {
            std::cout << i << ": " << red_allele_freqs3.at(i) << "\n";
        }

        for (unsigned int i = 0; i < red_allele_freqs2.size(); ++i) {
            REQUIRE(red_allele_freqs2.at(i) == Approx(red_allele_freqs3.at(i)).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing missing gene copies (1 of 3) with relative root",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing missing gene copies (1 of 3) for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        std::string nex_path2 = "data/singleton-n2.nex";
        std::string nex_path3 = "data/singleton-n3.nex";
        // Need to keep constant characters
        ComparisonRelativeRootPopulationTree tree2(nex_path2, ' ', true, false, false, false);
        ComparisonRelativeRootPopulationTree tree3(nex_path3, ' ', true, false, false, false);
        tree2.estimate_mutation_rate();
        tree3.estimate_mutation_rate();

        tree2.set_root_population_size(Ne * mu);
        tree3.set_root_population_size(Ne * mu);

        tree2.set_child_population_size(0, (Ne * mu));
        tree3.set_child_population_size(0, (Ne * mu));

        tree2.set_root_height(0.1);
        tree3.set_root_height(0.1);

        tree2.set_freq_1(0.5);
        tree3.set_freq_1(0.5);

        tree2.set_mutation_rate(1.0);
        tree3.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        unsigned int nsamples = 1000000;

        std::shared_ptr<GeneTreeSimNode> gtree2;
        std::shared_ptr<GeneTreeSimNode> gtree3;
        std::vector<unsigned int> down_sampled_counts = {2};

        std::vector<unsigned int> red_allele_counts2 = {0, 0, 0};
        std::vector<unsigned int> red_allele_counts3 = {0, 0, 0};

        for (unsigned int i = 0; i < nsamples; ++i) {
            gtree2 = tree2.simulate_gene_tree(0, rng, true);
            gtree3 = tree3.simulate_gene_tree(0, rng, true);

            auto pattern2 = tree2.simulate_biallelic_site(gtree2, rng);
            auto pattern3 = tree3.simulate_biallelic_site_sans_missing(gtree3, down_sampled_counts, rng);

            ++red_allele_counts2.at(pattern2.first.at(0));
            ++red_allele_counts3.at(pattern3.first.at(0));
        }

        unsigned int total2 = 0;
        unsigned int total3 = 0;
        for (unsigned int i = 0; i < red_allele_counts2.size(); ++i) {
            total2 += red_allele_counts2.at(i);
            total3 += red_allele_counts3.at(i);
        }

        REQUIRE(total2 == nsamples);
        REQUIRE(total3 == nsamples);

        std::vector<double> red_allele_freqs2 = {0.0, 0.0, 0.0};
        std::vector<double> red_allele_freqs3 = {0.0, 0.0, 0.0};
        for (unsigned int i = 0; i < red_allele_counts2.size(); ++i) {
            red_allele_freqs2.at(i) = red_allele_counts2.at(i) / (double)nsamples;
            red_allele_freqs3.at(i) = red_allele_counts3.at(i) / (double)nsamples;
        }

        std::cout << "Freqs of red allele counts for 2 gene copies:\n";
        for (unsigned int i = 0; i < red_allele_freqs2.size(); ++i) {
            std::cout << i << ": " << red_allele_freqs2.at(i) << "\n";
        }
        std::cout << "Freqs of red allele counts for 3 gene copies with 1 dropped:\n";
        for (unsigned int i = 0; i < red_allele_freqs3.size(); ++i) {
            std::cout << i << ": " << red_allele_freqs3.at(i) << "\n";
        }

        for (unsigned int i = 0; i < red_allele_freqs2.size(); ++i) {
            REQUIRE(red_allele_freqs2.at(i) == Approx(red_allele_freqs3.at(i)).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing missing gene copies (2 of 4) with relative root",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing missing gene copies (2 of 4) for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        std::string nex_path2 = "data/singleton-n2.nex";
        std::string nex_path3 = "data/singleton-n4.nex";
        // Need to keep constant characters
        ComparisonRelativeRootPopulationTree tree2(nex_path2, ' ', true, false, false, false);
        ComparisonRelativeRootPopulationTree tree3(nex_path3, ' ', true, false, false, false);
        tree2.estimate_mutation_rate();
        tree3.estimate_mutation_rate();

        tree2.set_root_population_size(Ne * mu);
        tree3.set_root_population_size(Ne * mu);

        tree2.set_child_population_size(0, (Ne * mu));
        tree3.set_child_population_size(0, (Ne * mu));

        tree2.set_root_height(0.1);
        tree3.set_root_height(0.1);

        tree2.set_freq_1(0.5);
        tree3.set_freq_1(0.5);

        tree2.set_mutation_rate(1.0);
        tree3.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        unsigned int nsamples = 1000000;

        std::shared_ptr<GeneTreeSimNode> gtree2;
        std::shared_ptr<GeneTreeSimNode> gtree3;
        std::vector<unsigned int> down_sampled_counts = {2};

        std::vector<unsigned int> red_allele_counts2 = {0, 0, 0};
        std::vector<unsigned int> red_allele_counts3 = {0, 0, 0};

        for (unsigned int i = 0; i < nsamples; ++i) {
            gtree2 = tree2.simulate_gene_tree(0, rng, true);
            gtree3 = tree3.simulate_gene_tree(0, rng, true);

            auto pattern2 = tree2.simulate_biallelic_site(gtree2, rng);
            auto pattern3 = tree3.simulate_biallelic_site_sans_missing(gtree3, down_sampled_counts, rng);

            ++red_allele_counts2.at(pattern2.first.at(0));
            ++red_allele_counts3.at(pattern3.first.at(0));
        }

        unsigned int total2 = 0;
        unsigned int total3 = 0;
        for (unsigned int i = 0; i < red_allele_counts2.size(); ++i) {
            total2 += red_allele_counts2.at(i);
            total3 += red_allele_counts3.at(i);
        }

        REQUIRE(total2 == nsamples);
        REQUIRE(total3 == nsamples);

        std::vector<double> red_allele_freqs2 = {0.0, 0.0, 0.0};
        std::vector<double> red_allele_freqs3 = {0.0, 0.0, 0.0};
        for (unsigned int i = 0; i < red_allele_counts2.size(); ++i) {
            red_allele_freqs2.at(i) = red_allele_counts2.at(i) / (double)nsamples;
            red_allele_freqs3.at(i) = red_allele_counts3.at(i) / (double)nsamples;
        }

        std::cout << "Freqs of red allele counts for 2 gene copies:\n";
        for (unsigned int i = 0; i < red_allele_freqs2.size(); ++i) {
            std::cout << i << ": " << red_allele_freqs2.at(i) << "\n";
        }
        std::cout << "Freqs of red allele counts for 3 gene copies with 1 dropped:\n";
        for (unsigned int i = 0; i < red_allele_freqs3.size(); ++i) {
            std::cout << i << ": " << red_allele_freqs3.at(i) << "\n";
        }

        for (unsigned int i = 0; i < red_allele_freqs2.size(); ++i) {
            REQUIRE(red_allele_freqs2.at(i) == Approx(red_allele_freqs3.at(i)).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing missing gene copies (6 of 10 and 8 of 10) with relative root",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing missing gene copies (6 of 10 and 8 of 10) for fully fixed pair") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        std::string nex_path2 = "data/hemi129-2-1.nex";
        std::string nex_path3 = "data/hemi129-5-5.nex";
        // Need to keep constant characters
        ComparisonRelativeRootPopulationTree tree2(nex_path2, ' ', true, true, false, false);
        ComparisonRelativeRootPopulationTree tree3(nex_path3, ' ', true, true, false, false);
        tree2.estimate_mutation_rate();
        tree3.estimate_mutation_rate();

        tree2.set_root_population_size(Ne * mu);
        tree3.set_root_population_size(Ne * mu);

        tree2.set_child_population_size(0, (Ne * mu));
        tree3.set_child_population_size(0, (Ne * mu));
        tree2.set_child_population_size(1, (Ne * mu));
        tree3.set_child_population_size(1, (Ne * mu));

        tree2.set_root_height(2.0 * Ne * mu);
        tree3.set_root_height(2.0 * Ne * mu);

        tree2.set_freq_1(0.5);
        tree3.set_freq_1(0.5);

        tree2.set_mutation_rate(1.0);
        tree3.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        unsigned int nsamples = 1000000;

        std::shared_ptr<GeneTreeSimNode> gtree2;
        std::shared_ptr<GeneTreeSimNode> gtree3;
        std::vector<unsigned int> down_sampled_counts = {4, 2};

        std::vector<unsigned int> pop1_red_allele_counts2 = {0, 0, 0, 0, 0};
        std::vector<unsigned int> pop1_red_allele_counts3 = {0, 0, 0, 0, 0};
        std::vector<unsigned int> pop2_red_allele_counts2 = {0, 0, 0};
        std::vector<unsigned int> pop2_red_allele_counts3 = {0, 0, 0};

        for (unsigned int i = 0; i < nsamples; ++i) {
            gtree2 = tree2.simulate_gene_tree(0, rng, true);
            gtree3 = tree3.simulate_gene_tree(0, rng, true);

            auto pattern2 = tree2.simulate_biallelic_site(gtree2, rng);
            auto pattern3 = tree3.simulate_biallelic_site_sans_missing(gtree3, down_sampled_counts, rng);

            ++pop1_red_allele_counts2.at(pattern2.first.at(0));
            ++pop1_red_allele_counts3.at(pattern3.first.at(0));
            ++pop2_red_allele_counts2.at(pattern2.first.at(1));
            ++pop2_red_allele_counts3.at(pattern3.first.at(1));
        }

        unsigned int pop1_total2 = 0;
        unsigned int pop1_total3 = 0;
        for (unsigned int i = 0; i < pop1_red_allele_counts2.size(); ++i) {
            pop1_total2 += pop1_red_allele_counts2.at(i);
            pop1_total3 += pop1_red_allele_counts3.at(i);
        }
        unsigned int pop2_total2 = 0;
        unsigned int pop2_total3 = 0;
        for (unsigned int i = 0; i < pop2_red_allele_counts2.size(); ++i) {
            pop2_total2 += pop2_red_allele_counts2.at(i);
            pop2_total3 += pop2_red_allele_counts3.at(i);
        }

        REQUIRE(pop1_total2 == nsamples);
        REQUIRE(pop1_total3 == nsamples);
        REQUIRE(pop2_total2 == nsamples);
        REQUIRE(pop2_total3 == nsamples);

        std::vector<double> pop1_red_allele_freqs2 = {0.0, 0.0, 0.0, 0.0, 0.0};
        std::vector<double> pop1_red_allele_freqs3 = {0.0, 0.0, 0.0, 0.0, 0.0};
        for (unsigned int i = 0; i < pop1_red_allele_counts2.size(); ++i) {
            pop1_red_allele_freqs2.at(i) = pop1_red_allele_counts2.at(i) / (double)nsamples;
            pop1_red_allele_freqs3.at(i) = pop1_red_allele_counts3.at(i) / (double)nsamples;
        }
        std::vector<double> pop2_red_allele_freqs2 = {0.0, 0.0, 0.0};
        std::vector<double> pop2_red_allele_freqs3 = {0.0, 0.0, 0.0};
        for (unsigned int i = 0; i < pop2_red_allele_counts2.size(); ++i) {
            pop2_red_allele_freqs2.at(i) = pop2_red_allele_counts2.at(i) / (double)nsamples;
            pop2_red_allele_freqs3.at(i) = pop2_red_allele_counts3.at(i) / (double)nsamples;
        }

        std::cout << "Pop 1 freqs of red allele counts for 4 gene copies:\n";
        for (unsigned int i = 0; i < pop1_red_allele_freqs2.size(); ++i) {
            std::cout << i << ": " << pop1_red_allele_freqs2.at(i) << "\n";
        }
        std::cout << "Pop 1 freqs of red allele counts for 10 gene copies with 6 dropped:\n";
        for (unsigned int i = 0; i < pop1_red_allele_freqs3.size(); ++i) {
            std::cout << i << ": " << pop1_red_allele_freqs3.at(i) << "\n";
        }
        std::cout << "Pop 2 freqs of red allele counts for 2 gene copies:\n";
        for (unsigned int i = 0; i < pop2_red_allele_freqs2.size(); ++i) {
            std::cout << i << ": " << pop2_red_allele_freqs2.at(i) << "\n";
        }
        std::cout << "Pop 2 freqs of red allele counts for 10 gene copies with 8 dropped:\n";
        for (unsigned int i = 0; i < pop2_red_allele_freqs3.size(); ++i) {
            std::cout << i << ": " << pop2_red_allele_freqs3.at(i) << "\n";
        }

        for (unsigned int i = 0; i < pop1_red_allele_freqs2.size(); ++i) {
            REQUIRE(pop1_red_allele_freqs2.at(i) == Approx(pop1_red_allele_freqs3.at(i)).epsilon(0.001));
        }
        for (unsigned int i = 0; i < pop2_red_allele_freqs2.size(); ++i) {
            REQUIRE(pop2_red_allele_freqs2.at(i) == Approx(pop2_red_allele_freqs3.at(i)).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing BaseTree::get_height_of_youngest_parent()", "[BaseTree]") {
    SECTION("Testing get_height_of_youngest_parent") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.1);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.08);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.06);
        std::shared_ptr<Node> internal3 = std::make_shared<Node>("internal3", 0.04);
        std::shared_ptr<Node> internal4 = std::make_shared<Node>("internal4", 0.04);
        internal4->set_height_parameter(internal3->get_height_parameter());
        std::shared_ptr<Node> internal5 = std::make_shared<Node>("internal5", 0.02);
        std::shared_ptr<Node> internal6 = std::make_shared<Node>("internal6", 0.02);
        internal6->set_height_parameter(internal5->get_height_parameter());
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        leaf1->fix_node_height();
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        leaf2->fix_node_height();
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        leaf3->fix_node_height();
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
        leaf4->fix_node_height();
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
        leaf5->fix_node_height();
        std::shared_ptr<Node> leaf6 = std::make_shared<Node>("leaf6", 0.0);
        leaf6->fix_node_height();
        std::shared_ptr<Node> leaf7 = std::make_shared<Node>("leaf7", 0.0);
        leaf7->fix_node_height();
        std::shared_ptr<Node> leaf8 = std::make_shared<Node>("leaf8", 0.0);
        leaf8->fix_node_height();

        internal3->add_child(leaf1);
        internal3->add_child(leaf2);
        internal4->add_child(leaf3);
        internal4->add_child(leaf4);
        internal5->add_child(leaf5);
        internal5->add_child(leaf6);
        internal6->add_child(leaf7);
        internal6->add_child(leaf8);

        internal1->add_child(internal5);
        internal1->add_child(internal3);
        internal2->add_child(internal6);
        internal2->add_child(internal4);

        root->add_child(internal1);
        root->add_child(internal2);
        BaseTree<Node> tree(root);

        REQUIRE(tree.get_leaf_node_count() == 8);
        REQUIRE(tree.get_node_count() == 15);
        REQUIRE(tree.get_number_of_node_heights() == 5);
        std::vector<double> expected_heights {0.02, 0.04, 0.06, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        REQUIRE_THROWS_AS(tree.get_youngest_parent(4), EcoevolityError &);
        REQUIRE_THROWS_AS(tree.get_height_of_youngest_parent(4), EcoevolityError &);
        REQUIRE_THROWS_AS(tree.get_index_of_youngest_parent(4), EcoevolityError &);
        REQUIRE(tree.get_youngest_parent(0)->get_height() == 0.06);
        REQUIRE(tree.get_height_of_youngest_parent(0) == 0.06);
        REQUIRE(tree.get_index_of_youngest_parent(0) == 2);
        REQUIRE(tree.get_youngest_parent(1)->get_height() == 0.06);
        REQUIRE(tree.get_height_of_youngest_parent(1) == 0.06);
        REQUIRE(tree.get_index_of_youngest_parent(1) == 2);
        REQUIRE(tree.get_youngest_parent(2)->get_height() == 0.1);
        REQUIRE(tree.get_height_of_youngest_parent(2) == 0.1);
        REQUIRE(tree.get_index_of_youngest_parent(2) == 4);
        REQUIRE(tree.get_youngest_parent(3)->get_height() == 0.1);
        REQUIRE(tree.get_height_of_youngest_parent(3) == 0.1);
        REQUIRE(tree.get_index_of_youngest_parent(3) == 4);
    }
}

TEST_CASE("Testing BaseTree::get_height_of_oldest_child()", "[BaseTree]") {
    SECTION("Testing get_height_of_oldest_child") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.1);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.08);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.06);
        std::shared_ptr<Node> internal3 = std::make_shared<Node>("internal3", 0.04);
        std::shared_ptr<Node> internal4 = std::make_shared<Node>("internal4", 0.04);
        internal4->set_height_parameter(internal3->get_height_parameter());
        std::shared_ptr<Node> internal5 = std::make_shared<Node>("internal5", 0.02);
        std::shared_ptr<Node> internal6 = std::make_shared<Node>("internal6", 0.02);
        internal6->set_height_parameter(internal5->get_height_parameter());
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        leaf1->fix_node_height();
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        leaf2->fix_node_height();
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        leaf3->fix_node_height();
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
        leaf4->fix_node_height();
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
        leaf5->fix_node_height();
        std::shared_ptr<Node> leaf6 = std::make_shared<Node>("leaf6", 0.0);
        leaf6->fix_node_height();
        std::shared_ptr<Node> leaf7 = std::make_shared<Node>("leaf7", 0.0);
        leaf7->fix_node_height();
        std::shared_ptr<Node> leaf8 = std::make_shared<Node>("leaf8", 0.0);
        leaf8->fix_node_height();

        internal3->add_child(leaf1);
        internal3->add_child(leaf2);
        internal4->add_child(leaf3);
        internal4->add_child(leaf4);
        internal5->add_child(leaf5);
        internal5->add_child(leaf6);
        internal6->add_child(leaf7);
        internal6->add_child(leaf8);

        internal1->add_child(internal5);
        internal1->add_child(internal3);
        internal2->add_child(internal6);
        internal2->add_child(internal4);

        root->add_child(internal1);
        root->add_child(internal2);
        BaseTree<Node> tree(root);

        REQUIRE(tree.get_leaf_node_count() == 8);
        REQUIRE(tree.get_node_count() == 15);
        REQUIRE(tree.get_number_of_node_heights() == 5);
        std::vector<double> expected_heights {0.02, 0.04, 0.06, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        REQUIRE(tree.get_oldest_child(0)->get_height() == 0.0);
        REQUIRE(tree.get_height_of_oldest_child(0) == 0.0);
        REQUIRE(tree.get_oldest_child(1)->get_height() == 0.0);
        REQUIRE(tree.get_height_of_oldest_child(1) == 0.0);
        REQUIRE(tree.get_oldest_child(2)->get_height() == 0.04);
        REQUIRE(tree.get_height_of_oldest_child(2) == 0.04);
        REQUIRE(tree.get_oldest_child(3)->get_height() == 0.04);
        REQUIRE(tree.get_height_of_oldest_child(3) == 0.04);
        REQUIRE(tree.get_oldest_child(4)->get_height() == 0.08);
        REQUIRE(tree.get_height_of_oldest_child(4) == 0.08);
    }
}
