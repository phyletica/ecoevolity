#include "catch.hpp"
#include "ecoevolity/basetree.hpp"
#include "ecoevolity/tree.hpp"
#include "ecoevolity/node.hpp"

TEST_CASE("Testing 3 leaves", "[treeio]") {
    SECTION("Testing 3 leaves") {
        std::string newick_tree_str = "((spa[&length=0.1,height=0.0,pop_size=0.001]:0.1,spb[&length=0.1,height=0.0,pop_size=0.002]:0.1)[&length=0.2,support=1.0,height=0.1,height_index=0,pop_size=0.003]:0.2,spc[&length=0.3,height=0.0,pop_size=0.004]:0.3)[&height=0.3,height_index=1,support=1.0,pop_size=0.005];";
        BaseTree<Node> tree(newick_tree_str);
        std::vector<double> expected_heights {0.1, 0.3};
        std::vector<double> heights = tree.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
        }

        std::map< int, std::set<Split> > splits = tree.get_splits_by_height_index(false);


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

        std::map< int, std::set<Split> > expected_splits = expected_tree.get_splits_by_height_index(false);

        REQUIRE(splits == expected_splits);

        // Changing positions of a and b should not matter
        newick_tree_str = "((spb[&length=0.1,height=0.0,pop_size=0.001]:0.1,spa[&length=0.1,height=0.0,pop_size=0.002]:0.1)[&length=0.2,support=1.0,height=0.1,height_index=0,pop_size=0.003]:0.2,spc[&length=0.3,height=0.0,pop_size=0.004]:0.3)[&height=0.3,height_index=1,support=1.0,pop_size=0.005];";
        BaseTree<Node> tree2(newick_tree_str);
        std::map< int, std::set<Split> > splits2 = tree2.get_splits_by_height_index(false);
        REQUIRE(splits2 == expected_splits);
        heights = tree2.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
        }

        // Changing positions of c and a + b should not matter
        newick_tree_str = "(spc[&length=0.3,height=0.0,pop_size=0.004]:0.3,(spb[&length=0.1,height=0.0,pop_size=0.001]:0.1,spa[&length=0.1,height=0.0,pop_size=0.002]:0.1)[&length=0.2,support=1.0,height=0.1,height_index=0,pop_size=0.003]:0.2)[&height=0.3,height_index=1,support=1.0,pop_size=0.005];";
        BaseTree<Node> tree3(newick_tree_str);
        std::map< int, std::set<Split> > splits3 = tree3.get_splits_by_height_index(false);
        REQUIRE(splits3 == expected_splits);
        heights = tree3.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
        }

        // Changing positions of a and c should matter
        newick_tree_str = "((spc[&length=0.1,height=0.0,pop_size=0.001]:0.1,spb[&length=0.1,height=0.0,pop_size=0.002]:0.1)[&length=0.2,support=1.0,height=0.1,height_index=0,pop_size=0.003]:0.2,spa[&length=0.3,height=0.0,pop_size=0.004]:0.3)[&height=0.3,height_index=1,support=1.0,pop_size=0.005];";
        BaseTree<Node> tree4(newick_tree_str);
        std::map< int, std::set<Split> > splits4 = tree4.get_splits_by_height_index(false);
        REQUIRE(splits4 != expected_splits);
        heights = tree4.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
        }
    }
}

TEST_CASE("Testing 3 leaves round trip", "[treeio]") {
    SECTION("Testing 3 leaves round trip") {
        std::string newick_tree_str = "((spa[&length=0.1,height=0.0,pop_size=0.001]:0.1,spb[&length=0.1,height=0.0,pop_size=0.002]:0.1)[&length=0.2,support=1.0,height=0.1,height_index=0,pop_size=0.003]:0.2,spc[&length=0.3,height=0.0,pop_size=0.004]:0.3)[&height=0.3,height_index=1,support=1.0,pop_size=0.005];";
        BaseTree<Node> tree(newick_tree_str);
        std::string written_newick_tree_str = tree.to_parentheses(true) + ";";
        BaseTree<Node> round_trip_tree(written_newick_tree_str);
        std::vector<double> expected_heights {0.1, 0.3};
        std::vector<double> heights = tree.get_node_heights();
        std::vector<double> round_trip_heights = round_trip_tree.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        REQUIRE(expected_heights.size() == round_trip_heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
            REQUIRE(round_trip_heights.at(i) == Approx(expected_heights.at(i)));
        }

        std::map< int, std::set<Split> > splits = tree.get_splits_by_height_index(false);
        std::map< int, std::set<Split> > round_trip_splits = round_trip_tree.get_splits_by_height_index(false);


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

        std::map< int, std::set<Split> > expected_splits = expected_tree.get_splits_by_height_index(false);

        REQUIRE(splits == expected_splits);
        REQUIRE(round_trip_splits == expected_splits);
    }
}


TEST_CASE("Testing 3 leaves with no comments", "[treeio]") {
    SECTION("Testing 3 leaves with no comments") {
        std::string newick_tree_str = "((spa:0.1,spb:0.1):0.2,spc:0.3);";
        BaseTree<Node> tree(newick_tree_str);
        std::vector<double> expected_heights {0.1, 0.3};
        std::vector<double> heights = tree.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
        }
        std::map< int, std::set<Split> > splits = tree.get_splits_by_height_index(false);


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

        std::map< int, std::set<Split> > expected_splits = expected_tree.get_splits_by_height_index(false);

        REQUIRE(splits == expected_splits);

        // Changing positions of a and b should not matter
        newick_tree_str = "((spb:0.1,spa:0.1):0.2,spc:0.3);";
        BaseTree<Node> tree2(newick_tree_str);
        std::map< int, std::set<Split> > splits2 = tree2.get_splits_by_height_index(false);
        REQUIRE(splits2 == expected_splits);
        heights = tree2.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
        }

        // Changing positions of c and a + b should not matter
        newick_tree_str = "(spc:0.3,(spb:0.1,spa:0.1):0.2);";
        BaseTree<Node> tree3(newick_tree_str);
        std::map< int, std::set<Split> > splits3 = tree3.get_splits_by_height_index(false);
        REQUIRE(splits3 == expected_splits);
        heights = tree3.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
        }

        // Changing positions of a and c should matter
        newick_tree_str = "((spc:0.1,spb:0.1):0.2,spa:0.3);";
        BaseTree<Node> tree4(newick_tree_str);
        std::map< int, std::set<Split> > splits4 = tree4.get_splits_by_height_index(false);
        REQUIRE(splits4 != expected_splits);
        heights = tree4.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
        }
    }
}

TEST_CASE("Testing 3 leaves with no comments round trip", "[treeio]") {
    SECTION("Testing 3 leaves with no comments round trip") {
        std::string newick_tree_str = "((spa:0.1,spb:0.1):0.2,spc:0.3);";
        BaseTree<Node> tree(newick_tree_str);
        std::string written_newick_tree_str = tree.to_parentheses(false) + ";";
        BaseTree<Node> round_trip_tree(written_newick_tree_str);
        std::vector<double> expected_heights {0.1, 0.3};
        std::vector<double> heights = tree.get_node_heights();
        std::vector<double> round_trip_heights = round_trip_tree.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        REQUIRE(expected_heights.size() == round_trip_heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
            REQUIRE(round_trip_heights.at(i) == Approx(expected_heights.at(i)));
        }
        std::map< int, std::set<Split> > splits = tree.get_splits_by_height_index(false);
        std::map< int, std::set<Split> > round_trip_splits = round_trip_tree.get_splits_by_height_index(false);


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

        std::map< int, std::set<Split> > expected_splits = expected_tree.get_splits_by_height_index(false);

        REQUIRE(splits == expected_splits);
        REQUIRE(round_trip_splits == expected_splits);
    }
}

TEST_CASE("Testing 3 leaves with different expected leaf indices that should not matter", "[treeio]") {
    SECTION("Testing 3 leaves with different expected tip indices that should not matter") {
        const std::string newick_tree_str = "((spa[&length=0.1,height=0.0,pop_size=0.001]:0.1,spb[&length=0.1,height=0.0,pop_size=0.002]:0.1)[&length=0.2,support=1.0,height=0.1,height_index=0,pop_size=0.003]:0.2,spc[&length=0.3,height=0.0,pop_size=0.004]:0.3)[&height=0.3,height_index=1,support=1.0,pop_size=0.005];";
        BaseTree<Node> tree(newick_tree_str);
        std::vector<double> expected_heights {0.1, 0.3};
        std::vector<double> heights = tree.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
        }
        std::map< int, std::set<Split> > splits = tree.get_splits_by_height_index(false);


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

        std::map< int, std::set<Split> > expected_splits = expected_tree.get_splits_by_height_index(false);

        // The BaseTree constructor indexes leaves after sorting them by their
        // labels, so these should NOT be equal
        REQUIRE(splits == expected_splits);
    }
}

TEST_CASE("Testing 3 leaves with different expected leaf indices that should matter", "[treeio]") {
    SECTION("Testing 3 leaves with different expected tip indices that should matter") {
        const std::string newick_tree_str = "((spa[&length=0.1,height=0.0,pop_size=0.001]:0.1,spb[&length=0.1,height=0.0,pop_size=0.002]:0.1)[&length=0.2,support=1.0,height=0.1,height_index=0,pop_size=0.003]:0.2,spc[&length=0.3,height=0.0,pop_size=0.004]:0.3)[&height=0.3,height_index=1,support=1.0,pop_size=0.005];";
        BaseTree<Node> tree(newick_tree_str);
        std::vector<double> expected_heights {0.1, 0.3};
        std::vector<double> heights = tree.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
        }
        std::map< int, std::set<Split> > splits = tree.get_splits_by_height_index(false);


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

        std::map< int, std::set<Split> > expected_splits = expected_tree.get_splits_by_height_index(false);

        // The BaseTree constructor indexes leaves after sorting them by their
        // labels, so these should NOT be equal
        REQUIRE(splits != expected_splits);
    }
}

TEST_CASE("Testing 5 leaves with polytomy and shared div", "[treeio]") {
    SECTION("Testing 5 leaves with polytomy and shared div") {
        std::string newick_tree_str = "((a:0.1,e:0.1)[&height=0.1,height_index=0]:0.2,(b:0.1,d:0.1)[&height=0.1,height_index=0]:0.2,c:0.3)[&height=0.3,height_index=1];";
        BaseTree<Node> tree(newick_tree_str);
        std::vector<double> expected_heights {0.1, 0.3};
        std::vector<double> heights = tree.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
        }

        std::map< int, std::set<Split> > splits = tree.get_splits_by_height_index(false);


        std::shared_ptr<Node> root = std::make_shared<Node>(4, "root", 0.3);
        std::shared_ptr<Node> internal_ae = std::make_shared<Node>(3, "internal_ae", 0.1);
        std::shared_ptr<Node> internal_bd = std::make_shared<Node>(3, "internal_bd", 0.1);
        internal_ae->set_height_parameter(internal_bd->get_height_parameter());
        std::shared_ptr<Node> a = std::make_shared<Node>(0, "a", 0.0);
        a->fix_node_height();
        std::shared_ptr<Node> b = std::make_shared<Node>(1, "b", 0.0);
        b->fix_node_height();
        std::shared_ptr<Node> c = std::make_shared<Node>(2, "c", 0.0);
        c->fix_node_height();
        std::shared_ptr<Node> d = std::make_shared<Node>(3, "d", 0.0);
        d->fix_node_height();
        std::shared_ptr<Node> e = std::make_shared<Node>(4, "e", 0.0);
        e->fix_node_height();

        internal_ae->add_child(a);
        internal_ae->add_child(e);

        internal_bd->add_child(b);
        internal_bd->add_child(d);

        root->add_child(internal_ae);
        root->add_child(internal_bd);
        root->add_child(c);
        BaseTree<Node> expected_tree(root);

        std::map< int, std::set<Split> > expected_splits = expected_tree.get_splits_by_height_index(false);

        REQUIRE(splits == expected_splits);

        // Changing positions of sisters should not matter
        newick_tree_str = "((d:0.1,b:0.1)[&height=0.1,height_index=0]:0.2,(e:0.1,a:0.1)[&height=0.1,height_index=0]:0.2,c:0.3)[&height=0.3,height_index=1];";
        BaseTree<Node> tree2(newick_tree_str);
        std::map< int, std::set<Split> > splits2 = tree2.get_splits_by_height_index(false);
        REQUIRE(splits2 == expected_splits);
        heights = tree2.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
        }

        // swapping sisters should matter
        newick_tree_str = "((a:0.1,d:0.1)[&height=0.1,height_index=0]:0.2,(b:0.1,e:0.1)[&height=0.1,height_index=0]:0.2,c:0.3)[&height=0.3,height_index=1];";
        BaseTree<Node> tree3(newick_tree_str);
        std::map< int, std::set<Split> > splits3 = tree3.get_splits_by_height_index(false);
        REQUIRE(splits3 != expected_splits);
        heights = tree3.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
        }
    }
}

TEST_CASE("Testing 5 leaves with polytomy and shared div, round trip", "[treeio]") {
    SECTION("Testing 5 leaves with polytomy and shared div, round trip") {
        std::string newick_tree_str = "((a:0.1,e:0.1)[&height=0.1,height_index=0]:0.2,(b:0.1,d:0.1)[&height=0.1,height_index=0]:0.2,c:0.3)[&height=0.3,height_index=1];";
        BaseTree<Node> tree(newick_tree_str);
        std::string written_newick_tree_str = tree.to_parentheses(true) + ";";
        BaseTree<Node> round_trip_tree(written_newick_tree_str);
        std::vector<double> expected_heights {0.1, 0.3};
        std::vector<double> heights = tree.get_node_heights();
        std::vector<double> round_trip_heights = round_trip_tree.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        REQUIRE(expected_heights.size() == round_trip_heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
            REQUIRE(round_trip_heights.at(i) == Approx(expected_heights.at(i)));
        }

        std::map< int, std::set<Split> > splits = tree.get_splits_by_height_index(false);
        std::map< int, std::set<Split> > round_trip_splits = round_trip_tree.get_splits_by_height_index(false);


        std::shared_ptr<Node> root = std::make_shared<Node>(4, "root", 0.3);
        std::shared_ptr<Node> internal_ae = std::make_shared<Node>(3, "internal_ae", 0.1);
        std::shared_ptr<Node> internal_bd = std::make_shared<Node>(3, "internal_bd", 0.1);
        internal_ae->set_height_parameter(internal_bd->get_height_parameter());
        std::shared_ptr<Node> a = std::make_shared<Node>(0, "a", 0.0);
        a->fix_node_height();
        std::shared_ptr<Node> b = std::make_shared<Node>(1, "b", 0.0);
        b->fix_node_height();
        std::shared_ptr<Node> c = std::make_shared<Node>(2, "c", 0.0);
        c->fix_node_height();
        std::shared_ptr<Node> d = std::make_shared<Node>(3, "d", 0.0);
        d->fix_node_height();
        std::shared_ptr<Node> e = std::make_shared<Node>(4, "e", 0.0);
        e->fix_node_height();

        internal_ae->add_child(a);
        internal_ae->add_child(e);

        internal_bd->add_child(b);
        internal_bd->add_child(d);

        root->add_child(internal_ae);
        root->add_child(internal_bd);
        root->add_child(c);
        BaseTree<Node> expected_tree(root);

        std::map< int, std::set<Split> > expected_splits = expected_tree.get_splits_by_height_index(false);

        REQUIRE(splits == expected_splits);
        REQUIRE(round_trip_splits == expected_splits);
    }
}

TEST_CASE("Testing 5 leaves with polytomy and shared div and no comments", "[treeio]") {
    SECTION("Testing 5 leaves with polytomy and no comments") {
        std::string newick_tree_str = "((a:0.1,e:0.1):0.2,(b:0.1,d:0.1):0.2,c:0.3);";
        BaseTree<Node> tree(newick_tree_str);
        std::vector<double> expected_heights {0.1, 0.1, 0.3};
        std::vector<double> heights = tree.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
        }

        // Cannot use height indices, because those must be parsed from
        // commments
        std::set<Split> splits;
        for (auto idx_splitset : tree.get_splits_by_height_index(false)) {
            for (auto split : idx_splitset.second) {
                splits.insert(split);
            }
        }

        std::shared_ptr<Node> root = std::make_shared<Node>(4, "root", 0.3);
        std::shared_ptr<Node> internal_ae = std::make_shared<Node>(3, "internal_ae", 0.1);
        std::shared_ptr<Node> internal_bd = std::make_shared<Node>(3, "internal_bd", 0.1);
        // With no comments, the shared div cannot be parsed
        // internal_ae->set_height_parameter(internal_bd->get_height_parameter());
        std::shared_ptr<Node> a = std::make_shared<Node>(0, "a", 0.0);
        a->fix_node_height();
        std::shared_ptr<Node> b = std::make_shared<Node>(1, "b", 0.0);
        b->fix_node_height();
        std::shared_ptr<Node> c = std::make_shared<Node>(2, "c", 0.0);
        c->fix_node_height();
        std::shared_ptr<Node> d = std::make_shared<Node>(3, "d", 0.0);
        d->fix_node_height();
        std::shared_ptr<Node> e = std::make_shared<Node>(4, "e", 0.0);
        e->fix_node_height();

        internal_ae->add_child(a);
        internal_ae->add_child(e);

        internal_bd->add_child(b);
        internal_bd->add_child(d);

        root->add_child(internal_ae);
        root->add_child(internal_bd);
        root->add_child(c);
        BaseTree<Node> expected_tree(root);

        std::set<Split> expected_splits;
        for (auto idx_splitset : expected_tree.get_splits_by_height_index(false)) {
            for (auto split : idx_splitset.second) {
                expected_splits.insert(split);
            }
        }

        REQUIRE(splits == expected_splits);

        // Changing positions of sisters should not matter
        newick_tree_str = "((d:0.1,b:0.1):0.2,(e:0.1,a:0.1):0.2,c:0.3);";
        BaseTree<Node> tree2(newick_tree_str);
        std::set<Split> splits2;
        for (auto idx_splitset : tree2.get_splits_by_height_index(false)) {
            for (auto split : idx_splitset.second) {
                splits2.insert(split);
            }
        }
        REQUIRE(splits2 == expected_splits);
        heights = tree2.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
        }

        // swapping sisters should matter
        newick_tree_str = "((a:0.1,d:0.1):0.2,(b:0.1,e:0.1):0.2,c:0.3);";
        BaseTree<Node> tree3(newick_tree_str);
        std::set<Split> splits3;
        for (auto idx_splitset : tree3.get_splits_by_height_index(false)) {
            for (auto split : idx_splitset.second) {
                splits3.insert(split);
            }
        }
        REQUIRE(splits3 != expected_splits);
        heights = tree3.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
        }
    }
}

TEST_CASE("Testing 5 leaves with polytomy and shared div and no comments, round trip", "[treeio]") {
    SECTION("Testing 5 leaves with polytomy and no comments, round trip") {
        std::string newick_tree_str = "((a:0.1,e:0.1):0.2,(b:0.1,d:0.1):0.2,c:0.3);";
        BaseTree<Node> tree(newick_tree_str);

        std::string written_newick_tree_str = tree.to_parentheses(false) + ";";
        std::string comments_written_newick_tree_str = tree.to_parentheses(true) + ";";
        BaseTree<Node> rt_tree(written_newick_tree_str);
        BaseTree<Node> rtc_tree(comments_written_newick_tree_str);

        std::vector<double> expected_heights {0.1, 0.1, 0.3};
        std::vector<double> heights = tree.get_node_heights();
        std::vector<double> rt_heights = rt_tree.get_node_heights();
        std::vector<double> rtc_heights = rtc_tree.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        REQUIRE(expected_heights.size() == rt_heights.size());
        REQUIRE(expected_heights.size() == rtc_heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
            REQUIRE(rt_heights.at(i) == Approx(expected_heights.at(i)));
            REQUIRE(rtc_heights.at(i) == Approx(expected_heights.at(i)));
        }

        // Cannot use height indices, because those must be parsed from
        // commments
        std::set<Split> splits;
        for (auto idx_splitset : tree.get_splits_by_height_index(false)) {
            for (auto split : idx_splitset.second) {
                splits.insert(split);
            }
        }
        std::set<Split> rt_splits;
        for (auto idx_splitset : rt_tree.get_splits_by_height_index(false)) {
            for (auto split : idx_splitset.second) {
                rt_splits.insert(split);
            }
        }
        std::set<Split> rtc_splits;
        for (auto idx_splitset : rtc_tree.get_splits_by_height_index(false)) {
            for (auto split : idx_splitset.second) {
                rtc_splits.insert(split);
            }
        }

        std::shared_ptr<Node> root = std::make_shared<Node>(4, "root", 0.3);
        std::shared_ptr<Node> internal_ae = std::make_shared<Node>(3, "internal_ae", 0.1);
        std::shared_ptr<Node> internal_bd = std::make_shared<Node>(3, "internal_bd", 0.1);
        // With no comments, the shared div cannot be parsed
        // internal_ae->set_height_parameter(internal_bd->get_height_parameter());
        std::shared_ptr<Node> a = std::make_shared<Node>(0, "a", 0.0);
        a->fix_node_height();
        std::shared_ptr<Node> b = std::make_shared<Node>(1, "b", 0.0);
        b->fix_node_height();
        std::shared_ptr<Node> c = std::make_shared<Node>(2, "c", 0.0);
        c->fix_node_height();
        std::shared_ptr<Node> d = std::make_shared<Node>(3, "d", 0.0);
        d->fix_node_height();
        std::shared_ptr<Node> e = std::make_shared<Node>(4, "e", 0.0);
        e->fix_node_height();

        internal_ae->add_child(a);
        internal_ae->add_child(e);

        internal_bd->add_child(b);
        internal_bd->add_child(d);

        root->add_child(internal_ae);
        root->add_child(internal_bd);
        root->add_child(c);
        BaseTree<Node> expected_tree(root);

        std::set<Split> expected_splits;
        for (auto idx_splitset : expected_tree.get_splits_by_height_index(false)) {
            for (auto split : idx_splitset.second) {
                expected_splits.insert(split);
            }
        }

        REQUIRE(splits == expected_splits);
        REQUIRE(rt_splits == expected_splits);
        REQUIRE(rtc_splits == expected_splits);
    }
}

TEST_CASE("Testing BasePopulationTree with 3 leaves", "[treeio]") {
    SECTION("Testing 3 leaves") {
        std::string newick_tree_str = "((spa[&length=0.1,height=0.0,pop_size=1.0]:0.1,spb[&length=0.1,height=0.0,pop_size=2.0]:0.1)[&length=0.2,support=1.0,height=0.1,height_index=0,pop_size=3.0]:0.2,spc[&length=0.3,height=0.0,pop_size=4.0]:0.3)[&height=0.3,height_index=1,support=1.0,pop_size=5.0];";
        std::istringstream newick_tree_stream(newick_tree_str);
        BasePopulationTree first_tree(newick_tree_stream, "relaxedphyliptree");
        std::string written_newick_tree_str = first_tree.to_parentheses(true) + ";";
        std::istringstream written_newick_tree_stream(written_newick_tree_str);
        BasePopulationTree tree(written_newick_tree_stream, "relaxedphyliptree");
        std::vector<double> expected_heights {0.1, 0.3};
        std::vector<double> heights = tree.get_node_heights();
        REQUIRE(expected_heights.size() == heights.size());
        for (unsigned int i = 0; i < heights.size(); ++ i) {
            REQUIRE(heights.at(i) == Approx(expected_heights.at(i)));
        }

        std::vector<double> expected_pop_sizes {1.0, 2.0, 3.0, 4.0, 5.0};
        double expected_root_pop_size = 5.0;
        std::vector<double> pop_sizes = tree.get_population_sizes();
        std::sort(pop_sizes.begin(), pop_sizes.end());
        REQUIRE(tree.get_root_population_size() == expected_root_pop_size);
        REQUIRE(pop_sizes == expected_pop_sizes);

        std::map< int, std::set<Split> > splits = tree.get_splits_by_height_index(false);


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

        std::map< int, std::set<Split> > expected_splits = expected_tree.get_splits_by_height_index(false);

        REQUIRE(splits == expected_splits);
    }
}
