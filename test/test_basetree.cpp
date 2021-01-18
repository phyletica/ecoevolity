#include "catch.hpp"
#include "ecoevolity/node.hpp"
#include "ecoevolity/basetree.hpp"


// NOTE: most tests of basetree.hpp are in test_tree.cpp. Starting this file
// for new tests of basetree.hpp and the goal of migrating old tests here too.

TEST_CASE("Testing BaseTree::get_min_height_diff with comb tree", "[BaseTree]") {
    SECTION("Testing BaseTree::draw_from_prior") {

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.25);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
        std::shared_ptr<Node> leaf6 = std::make_shared<Node>("leaf6", 0.0);
        std::shared_ptr<Node> leaf7 = std::make_shared<Node>("leaf7", 0.0);
        std::shared_ptr<Node> leaf8 = std::make_shared<Node>("leaf8", 0.0);
        std::shared_ptr<Node> leaf9 = std::make_shared<Node>("leaf9", 0.0);
        std::shared_ptr<Node> leaf10 = std::make_shared<Node>("leaf10", 0.0);
        std::shared_ptr<Node> leaf11 = std::make_shared<Node>("leaf11", 0.0);

        root->add_child(leaf1);
        root->add_child(leaf2);
        root->add_child(leaf3);
        root->add_child(leaf4);
        root->add_child(leaf5);
        root->add_child(leaf6);
        root->add_child(leaf7);
        root->add_child(leaf8);
        root->add_child(leaf9);
        root->add_child(leaf10);
        root->add_child(leaf11);

        BaseTree<Node> tree(root);

        REQUIRE(tree.get_min_height_diff() > 0.0);
        REQUIRE(std::isinf(tree.get_min_height_diff()));
    }
}

TEST_CASE("Testing BaseTree::get_min_height_diff", "[BaseTree]") {
    SECTION("Testing BaseTree::draw_from_prior") {

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.25);
        std::shared_ptr<Node> internal0a = std::make_shared<Node>("internal0a", 0.05);
        std::shared_ptr<Node> internal0b = std::make_shared<Node>("internal0b", 0.05);
        internal0a->set_height_parameter(internal0b->get_height_parameter());
        std::shared_ptr<Node> internal1a = std::make_shared<Node>("internal1a", 0.1);
        std::shared_ptr<Node> internal1b = std::make_shared<Node>("internal1b", 0.1);
        internal1a->set_height_parameter(internal1b->get_height_parameter());
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.15);
        std::shared_ptr<Node> internal3 = std::make_shared<Node>("internal3", 0.17);
        std::shared_ptr<Node> internal4 = std::make_shared<Node>("internal4", 0.2);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
        std::shared_ptr<Node> leaf6 = std::make_shared<Node>("leaf6", 0.0);
        std::shared_ptr<Node> leaf7 = std::make_shared<Node>("leaf7", 0.0);
        std::shared_ptr<Node> leaf8 = std::make_shared<Node>("leaf8", 0.0);
        std::shared_ptr<Node> leaf9 = std::make_shared<Node>("leaf9", 0.0);
        std::shared_ptr<Node> leaf10 = std::make_shared<Node>("leaf10", 0.0);
        std::shared_ptr<Node> leaf11 = std::make_shared<Node>("leaf11", 0.0);

        internal0a->add_child(leaf1);
        internal0a->add_child(leaf2);
        internal1a->add_child(leaf3);
        internal1a->add_child(leaf4);
        internal1a->add_child(leaf5);

        internal2->add_child(internal0a);
        internal2->add_child(internal1a);

        internal0b->add_child(leaf6);
        internal0b->add_child(leaf7);
        internal1b->add_child(leaf8);
        internal1b->add_child(leaf9);
        internal1b->add_child(leaf10);

        internal3->add_child(internal0b);
        internal3->add_child(internal1b);

        internal4->add_child(internal2);
        internal4->add_child(internal3);

        root->add_child(leaf11);
        root->add_child(internal4);

        BaseTree<Node> tree(root);

        REQUIRE(tree.get_min_height_diff() == Approx(0.02));
    }
}
