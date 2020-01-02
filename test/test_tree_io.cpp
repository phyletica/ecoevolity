#include "catch.hpp"
#include "ecoevolity/basetree.hpp"
#include "ecoevolity/node.hpp"

TEST_CASE("Testing 3 leaves", "[split]") {
    SECTION("Testing 3 leaves") {
        std::string newick_tree_str = "((spa[&length=0.1,height=0.0,pop_size=0.001]:0.1,spb[&length=0.1,height=0.0,pop_size=0.002]:0.1)[&length=0.2,support=1.0,height=0.1,pop_size=0.003]:0.2,spc[&length=0.3,height=0.0,pop_size=0.004]:0.3)[&height=0.3,support=1.0,pop_size=0.005];";
        BaseTree<Node> tree(newick_tree_str);
        std::set< std::pair< unsigned int, Split> > splits = tree.get_splits(false);


        std::shared_ptr<Node> root = std::make_shared<Node>(4, "root", 0.3);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(3, "internal1", 0.1);
        std::shared_ptr<Node> spa = std::make_shared<Node>(0, "spa", 0.0);
        spa->fix_node_height();
        std::shared_ptr<Node> spb = std::make_shared<Node>(1, "spb", 0.0);
        spb->fix_node_height();
        std::shared_ptr<Node> spc = std::make_shared<Node>(2, "spc", 0.0);
        spc->fix_node_height();

        internal1->add_child(spa);
        internal1->add_child(spb);

        root->add_child(internal1);
        root->add_child(spc);
        BaseTree<Node> expected_tree(root);

        std::set< std::pair< unsigned int, Split> > expected_splits = expected_tree.get_splits(false);

        REQUIRE(splits == expected_splits);

        // Changing positions of a and b should not matter
        newick_tree_str = "((spb[&length=0.1,height=0.0,pop_size=0.001]:0.1,spa[&length=0.1,height=0.0,pop_size=0.002]:0.1)[&length=0.2,support=1.0,height=0.1,pop_size=0.003]:0.2,spc[&length=0.3,height=0.0,pop_size=0.004]:0.3)[&height=0.3,support=1.0,pop_size=0.005];";
        BaseTree<Node> tree2(newick_tree_str);
        std::set< std::pair< unsigned int, Split> > splits2 = tree2.get_splits(false);
        REQUIRE(splits2 == expected_splits);

        // Changing positions of c and a + b should not matter
        newick_tree_str = "(spc[&length=0.3,height=0.0,pop_size=0.004]:0.3,(spb[&length=0.1,height=0.0,pop_size=0.001]:0.1,spa[&length=0.1,height=0.0,pop_size=0.002]:0.1)[&length=0.2,support=1.0,height=0.1,pop_size=0.003]:0.2)[&height=0.3,support=1.0,pop_size=0.005];";
        BaseTree<Node> tree3(newick_tree_str);
        std::set< std::pair< unsigned int, Split> > splits3 = tree3.get_splits(false);
        REQUIRE(splits3 == expected_splits);

        // Changing positions of a and c should matter
        newick_tree_str = "((spc[&length=0.1,height=0.0,pop_size=0.001]:0.1,spb[&length=0.1,height=0.0,pop_size=0.002]:0.1)[&length=0.2,support=1.0,height=0.1,pop_size=0.003]:0.2,spa[&length=0.3,height=0.0,pop_size=0.004]:0.3)[&height=0.3,support=1.0,pop_size=0.005];";
        BaseTree<Node> tree4(newick_tree_str);
        std::set< std::pair< unsigned int, Split> > splits4 = tree4.get_splits(false);
        REQUIRE(splits4 != expected_splits);
    }
}


TEST_CASE("Testing 3 leaves with no comments", "[split]") {
    SECTION("Testing 3 leaves with no comments") {
        std::string newick_tree_str = "((spa:0.1,spb:0.1):0.2,spc:0.3);";
        BaseTree<Node> tree(newick_tree_str);
        std::set< std::pair< unsigned int, Split> > splits = tree.get_splits(false);


        std::shared_ptr<Node> root = std::make_shared<Node>(4, "root", 0.3);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(3, "internal1", 0.1);
        std::shared_ptr<Node> spa = std::make_shared<Node>(0, "spa", 0.0);
        spa->fix_node_height();
        std::shared_ptr<Node> spb = std::make_shared<Node>(1, "spb", 0.0);
        spb->fix_node_height();
        std::shared_ptr<Node> spc = std::make_shared<Node>(2, "spc", 0.0);
        spc->fix_node_height();

        internal1->add_child(spa);
        internal1->add_child(spb);

        root->add_child(internal1);
        root->add_child(spc);
        BaseTree<Node> expected_tree(root);

        std::set< std::pair< unsigned int, Split> > expected_splits = expected_tree.get_splits(false);

        REQUIRE(splits == expected_splits);

        // Changing positions of a and b should not matter
        newick_tree_str = "((spb:0.1,spa:0.1):0.2,spc:0.3);";
        BaseTree<Node> tree2(newick_tree_str);
        std::set< std::pair< unsigned int, Split> > splits2 = tree2.get_splits(false);
        REQUIRE(splits2 == expected_splits);

        // Changing positions of c and a + b should not matter
        newick_tree_str = "(spc:0.3,(spb:0.1,spa:0.1):0.2);";
        BaseTree<Node> tree3(newick_tree_str);
        std::set< std::pair< unsigned int, Split> > splits3 = tree3.get_splits(false);
        REQUIRE(splits3 == expected_splits);

        // Changing positions of a and c should matter
        newick_tree_str = "((spc:0.1,spb:0.1):0.2,spa:0.3);";
        BaseTree<Node> tree4(newick_tree_str);
        std::set< std::pair< unsigned int, Split> > splits4 = tree4.get_splits(false);
        REQUIRE(splits4 != expected_splits);
    }
}

TEST_CASE("Testing 3 leaves with different expected leaf indices that should not matter", "[split]") {
    SECTION("Testing 3 leaves with different expected tip indices that should not matter") {
        const std::string newick_tree_str = "((spa[&length=0.1,height=0.0,pop_size=0.001]:0.1,spb[&length=0.1,height=0.0,pop_size=0.002]:0.1)[&length=0.2,support=1.0,height=0.1,pop_size=0.003]:0.2,spc[&length=0.3,height=0.0,pop_size=0.004]:0.3)[&height=0.3,support=1.0,pop_size=0.005];";
        BaseTree<Node> tree(newick_tree_str);
        std::set< std::pair< unsigned int, Split> > splits = tree.get_splits(false);


        std::shared_ptr<Node> root = std::make_shared<Node>(4, "root", 0.3);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(3, "internal1", 0.1);
        std::shared_ptr<Node> spa = std::make_shared<Node>(1, "spa", 0.0);
        spa->fix_node_height();
        std::shared_ptr<Node> spb = std::make_shared<Node>(0, "spb", 0.0);
        spb->fix_node_height();
        std::shared_ptr<Node> spc = std::make_shared<Node>(2, "spc", 0.0);
        spc->fix_node_height();

        internal1->add_child(spa);
        internal1->add_child(spb);

        root->add_child(internal1);
        root->add_child(spc);
        BaseTree<Node> expected_tree(root);

        std::set< std::pair< unsigned int, Split> > expected_splits = expected_tree.get_splits(false);

        // The BaseTree constructor indexes leaves after sorting them by their
        // labels, so these should NOT be equal
        REQUIRE(splits == expected_splits);
    }
}

TEST_CASE("Testing 3 leaves with different expected leaf indices that should matter", "[split]") {
    SECTION("Testing 3 leaves with different expected tip indices that should matter") {
        const std::string newick_tree_str = "((spa[&length=0.1,height=0.0,pop_size=0.001]:0.1,spb[&length=0.1,height=0.0,pop_size=0.002]:0.1)[&length=0.2,support=1.0,height=0.1,pop_size=0.003]:0.2,spc[&length=0.3,height=0.0,pop_size=0.004]:0.3)[&height=0.3,support=1.0,pop_size=0.005];";
        BaseTree<Node> tree(newick_tree_str);
        std::set< std::pair< unsigned int, Split> > splits = tree.get_splits(false);


        std::shared_ptr<Node> root = std::make_shared<Node>(4, "root", 0.3);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(3, "internal1", 0.1);
        std::shared_ptr<Node> spa = std::make_shared<Node>(0, "spa", 0.0);
        spa->fix_node_height();
        std::shared_ptr<Node> spb = std::make_shared<Node>(2, "spb", 0.0);
        spb->fix_node_height();
        std::shared_ptr<Node> spc = std::make_shared<Node>(1, "spc", 0.0);
        spc->fix_node_height();

        internal1->add_child(spa);
        internal1->add_child(spb);

        root->add_child(internal1);
        root->add_child(spc);
        BaseTree<Node> expected_tree(root);

        std::set< std::pair< unsigned int, Split> > expected_splits = expected_tree.get_splits(false);

        // The BaseTree constructor indexes leaves after sorting them by their
        // labels, so these should NOT be equal
        REQUIRE(splits != expected_splits);
    }
}
