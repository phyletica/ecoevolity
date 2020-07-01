#include "catch.hpp"
#include "ecoevolity/split.hpp"
#include "ecoevolity/basetree.hpp"
#include "ecoevolity/node.hpp"

TEST_CASE("Testing 2 leaves", "[split]") {
    SECTION("Testing 2 leaves") {
        Split s10;
        s10.resize(2);
        s10.set_leaf_bit(0);
        REQUIRE(s10.get_leaf_bit(0) == 1);
        REQUIRE(s10.get_leaf_bit(1) == 0);
        REQUIRE(s10.as_string() == "10");
        std::vector<unsigned int> expected_leaf_indices = {0};
        REQUIRE(s10.get_leaf_indices() == expected_leaf_indices);

        Split s01;
        s01.resize(2);
        s01.set_leaf_bit(1);
        REQUIRE(s01.get_leaf_bit(0) == 0);
        REQUIRE(s01.get_leaf_bit(1) == 1);
        REQUIRE(s01.as_string() == "01");
        expected_leaf_indices = {1};
        REQUIRE(s01.get_leaf_indices() == expected_leaf_indices);

        Split s00;
        s00.resize(2);
        REQUIRE(s00.get_leaf_bit(0) == 0);
        REQUIRE(s00.get_leaf_bit(1) == 0);
        REQUIRE(s00.as_string() == "00");
        expected_leaf_indices = {};
        REQUIRE(s00.get_leaf_indices() == expected_leaf_indices);

        Split s11;
        s11.resize(2);
        s11.set_leaf_bit(0);
        s11.set_leaf_bit(1);
        REQUIRE(s11.get_leaf_bit(0) == 1);
        REQUIRE(s11.get_leaf_bit(1) == 1);
        REQUIRE(s11.as_string() == "11");
        expected_leaf_indices = {0, 1};
        REQUIRE(s11.get_leaf_indices() == expected_leaf_indices);

        REQUIRE(s10 != s01);
        REQUIRE(! (s10 == s01));
        REQUIRE(! s10.is_equivalent(s01, true));
        REQUIRE(s10.is_equivalent(s01, false));
        REQUIRE(s10.is_compatible(s01));
        REQUIRE(! s10.conflicts_with(s01));

        REQUIRE(s00 != s11);
        REQUIRE(! (s00 == s11));
        REQUIRE(! s00.is_equivalent(s11, true));
        REQUIRE(s00.is_equivalent(s11, false));
        REQUIRE(s00.is_compatible(s11));
        REQUIRE(! s00.conflicts_with(s11));

        REQUIRE(s01 != s11);
        REQUIRE(! (s01 == s11));
        REQUIRE(! s01.is_equivalent(s11, true));
        REQUIRE(! s01.is_equivalent(s11, false));
        REQUIRE(s01.is_compatible(s11));
        REQUIRE(! s01.conflicts_with(s11));

        Split ss10(s10);
        REQUIRE(ss10.as_string() == "10");
        REQUIRE(ss10.get_leaf_bit(0) == 1);
        REQUIRE(ss10.get_leaf_bit(1) == 0);
        REQUIRE(ss10 == s10);
        REQUIRE(! (ss10 != s10));
        REQUIRE(ss10.is_equivalent(s10, true));
        REQUIRE(ss10.is_equivalent(s10, false));
        REQUIRE(ss10.is_compatible(s10));
        REQUIRE(! ss10.conflicts_with(s10));

        REQUIRE(ss10.as_string('-', '*') == "*-");

        REQUIRE(! s10.is_proper_subset_of(s01));
        REQUIRE(! s10.is_proper_superset_of(s01));
        REQUIRE(! s01.is_proper_subset_of(s10));
        REQUIRE(! s01.is_proper_superset_of(s10));

        REQUIRE(s10.is_proper_subset_of(s11));
        REQUIRE(s11.is_proper_superset_of(s10));
        REQUIRE(s01.is_proper_subset_of(s11));
        REQUIRE(s11.is_proper_superset_of(s01));

        REQUIRE(! s10.is_proper_superset_of(s11));
        REQUIRE(! s11.is_proper_subset_of(s10));
        REQUIRE(! s01.is_proper_superset_of(s11));
        REQUIRE(! s11.is_proper_subset_of(s01));

        // Wolfram tells me this is so; an empty set is proper subset of any
        // nonempty set
        REQUIRE(s00.is_proper_subset_of(s11));
        REQUIRE(s11.is_proper_superset_of(s00));

        // "proper" disqualifies equivalent sets from being considered
        // sub/supersets of each other
        REQUIRE(! ss10.is_proper_subset_of(s10));
        REQUIRE(! ss10.is_proper_superset_of(s10));
        REQUIRE(! s10.is_proper_subset_of(ss10));
        REQUIRE(! s10.is_proper_superset_of(ss10));
    }
}

TEST_CASE("Testing 5 leaves", "[split]") {
    SECTION("Testing 5 leaves") {
        Split s10100;
        s10100.resize(5);
        s10100.set_leaf_bit(0);
        s10100.set_leaf_bit(2);
        REQUIRE(s10100.as_string() == "10100");
        std::vector<unsigned int> expected_leaf_indices = {0, 2};
        REQUIRE(s10100.get_leaf_indices() == expected_leaf_indices);

        Split s01010;
        s01010.resize(5);
        s01010.set_leaf_bit(1);
        s01010.set_leaf_bit(3);
        REQUIRE(s01010.as_string() == "01010");
        expected_leaf_indices = {1, 3};
        REQUIRE(s01010.get_leaf_indices() == expected_leaf_indices);

        REQUIRE(! (s10100 == s01010));
        REQUIRE(s10100 != s01010);
        REQUIRE(! s10100.is_equivalent(s01010, false));
        REQUIRE(! s10100.is_equivalent(s01010, true));
        REQUIRE(s10100.is_compatible(s01010));
        REQUIRE(! s10100.conflicts_with(s01010));

        Split s01011;
        s01011.resize(5);
        s01011.set_leaf_bit(1);
        s01011.set_leaf_bit(3);
        s01011.set_leaf_bit(4);
        REQUIRE(s01011.as_string() == "01011");
        expected_leaf_indices = {1, 3, 4};
        REQUIRE(s01011.get_leaf_indices() == expected_leaf_indices);

        REQUIRE(! (s10100 == s01011));
        REQUIRE(s10100 != s01011);
        REQUIRE(s10100.is_equivalent(s01011, false));
        REQUIRE(! s10100.is_equivalent(s01011, true));
        REQUIRE(s10100.is_compatible(s01011));
        REQUIRE(! s10100.conflicts_with(s01011));

        REQUIRE(! (s01010 == s01011));
        REQUIRE(s01010 != s01011);
        REQUIRE(! s01010.is_equivalent(s01011, false));
        REQUIRE(! s01010.is_equivalent(s01011, true));
        REQUIRE(s01010.is_compatible(s01011));
        REQUIRE(! s01010.conflicts_with(s01011));

        Split s11000;
        s11000.resize(5);
        s11000.set_leaf_bit(0);
        s11000.set_leaf_bit(1);
        REQUIRE(s11000.as_string() == "11000");

        REQUIRE(! (s11000 == s10100));
        REQUIRE(s11000 != s10100);
        REQUIRE(! s11000.is_equivalent(s10100, false));
        REQUIRE(! s11000.is_equivalent(s10100, true));
        REQUIRE(! s11000.is_compatible(s10100));
        REQUIRE(s11000.conflicts_with(s10100));

        REQUIRE(! (s11000 == s01010));
        REQUIRE(s11000 != s01010);
        REQUIRE(! s11000.is_equivalent(s01010, false));
        REQUIRE(! s11000.is_equivalent(s01010, true));
        REQUIRE(! s11000.is_compatible(s01010));
        REQUIRE(s11000.conflicts_with(s01010));

        REQUIRE(! (s11000 == s01011));
        REQUIRE(s11000 != s01011);
        REQUIRE(! s11000.is_equivalent(s01011, false));
        REQUIRE(! s11000.is_equivalent(s01011, true));
        REQUIRE(! s11000.is_compatible(s01011));
        REQUIRE(s11000.conflicts_with(s01011));
    }
}

TEST_CASE("Testing 100000 leaves", "[split]") {
    SECTION("Testing 100000 leaves") {
        unsigned int nleaves = 100000;
        Split s;
        s.resize(nleaves);
        s.set_leaf_bit(0);
        s.set_leaf_bit(nleaves - 1);

        std::ostringstream expected;
        expected << "1";
        for (unsigned int i = 0; i < nleaves - 2; ++i) {
            expected << "0";
        }
        expected << "1";
        REQUIRE(s.as_string() == expected.str());
    }
}

TEST_CASE("Testing tree with 5 leaves", "[split]") {
    SECTION("Testing tree with 5 leaves") {
        std::string newick_tree_str = "((a:0.1,b:0.1)[&height=0.1,height_index=0]:0.3,(c:0.3,(d:0.2,e:0.2)[&height=0.2,height_index=1]:0.1)[&height=0.3,height_index=2]:0.1)[&height=0.4,height_index=3];";
        BaseTree<Node> tree_order_1(newick_tree_str);

        newick_tree_str = "((a:0.2,b:0.2)[&height=0.2,height_index=1]:0.2,(c:0.3,(d:0.1,e:0.1)[&height=0.1,height_index=0]:0.2)[&height=0.3,height_index=2]:0.1)[&height=0.4,height_index=3];";
        BaseTree<Node> tree_order_2(newick_tree_str);

        newick_tree_str = "((a:0.3,b:0.3)[&height=0.3,height_index=2]:0.1,(c:0.2,(d:0.1,e:0.1)[&height=0.1,height_index=0]:0.1)[&height=0.2,height_index=1]:0.2)[&height=0.4,height_index=3];";
        BaseTree<Node> tree_order_3(newick_tree_str);

        REQUIRE(tree_order_1.get_splits() == tree_order_2.get_splits());
        REQUIRE(tree_order_1.get_splits() == tree_order_3.get_splits());
        REQUIRE(tree_order_2.get_splits() == tree_order_3.get_splits());

        REQUIRE(tree_order_1.get_splits_by_height_index() != tree_order_2.get_splits_by_height_index());
        REQUIRE(tree_order_1.get_splits_by_height_index() != tree_order_3.get_splits_by_height_index());
        REQUIRE(tree_order_2.get_splits_by_height_index() != tree_order_3.get_splits_by_height_index());

        newick_tree_str = "((a:0.1,b:0.1)[&height=0.1,height_index=0]:0.3,(c:0.3,(d:0.1,e:0.1)[&height=0.1,height_index=0]:0.2)[&height=0.3,height_index=1]:0.1)[&height=0.4,height_index=2];";
        BaseTree<Node> shared_tree_1(newick_tree_str);

        newick_tree_str = "((a:0.3,b:0.3)[&height=0.3,height_index=1]:0.1,(c:0.3,(d:0.1,e:0.1)[&height=0.1,height_index=0]:0.2)[&height=0.3,height_index=1]:0.1)[&height=0.4,height_index=2];";
        BaseTree<Node> shared_tree_2(newick_tree_str);

        REQUIRE(shared_tree_1.get_splits() != tree_order_1.get_splits());
        REQUIRE(shared_tree_1.get_splits() != tree_order_2.get_splits());
        REQUIRE(shared_tree_1.get_splits() != tree_order_3.get_splits());
        REQUIRE(shared_tree_1.get_splits() != shared_tree_2.get_splits());

        REQUIRE(shared_tree_2.get_splits() != tree_order_1.get_splits());
        REQUIRE(shared_tree_2.get_splits() != tree_order_2.get_splits());
        REQUIRE(shared_tree_2.get_splits() != tree_order_3.get_splits());
        REQUIRE(shared_tree_2.get_splits() != shared_tree_1.get_splits());

        REQUIRE(shared_tree_1.get_splits_by_height_index() != tree_order_1.get_splits_by_height_index());
        REQUIRE(shared_tree_1.get_splits_by_height_index() != tree_order_2.get_splits_by_height_index());
        REQUIRE(shared_tree_1.get_splits_by_height_index() != tree_order_3.get_splits_by_height_index());
        REQUIRE(shared_tree_1.get_splits_by_height_index() != shared_tree_2.get_splits_by_height_index());

        REQUIRE(shared_tree_2.get_splits_by_height_index() != tree_order_1.get_splits_by_height_index());
        REQUIRE(shared_tree_2.get_splits_by_height_index() != tree_order_2.get_splits_by_height_index());
        REQUIRE(shared_tree_2.get_splits_by_height_index() != tree_order_3.get_splits_by_height_index());
        REQUIRE(shared_tree_2.get_splits_by_height_index() != shared_tree_1.get_splits_by_height_index());
    }
}

TEST_CASE("Testing trees with different shared orders" , "[split]") {
    SECTION("Testing tree with different shared orders") {
        std::string newick_tree_str = "(((a:0.1,b:0.1)[&height=0.1,height_index=0]:0.2,(c:0.2,d:0.2)[&height=0.2,height_index=1]:0.1)[&height=0.3,height_index=2]:0.1,((e:0.1,f:0.1)[&height=0.1,height_index=0]:0.2,(g:0.2,h:0.2)[&height=0.2,height_index=1]:0.1)[&height=0.3,height_index=2]:0.1)[&height=0.4,height_index=3];";
        BaseTree<Node> tree_order_1(newick_tree_str);

        newick_tree_str = "(((c:0.1,d:0.1)[&height=0.1,height_index=0]:0.2,(a:0.2,b:0.2)[&height=0.2,height_index=1]:0.1)[&height=0.3,height_index=2]:0.1,((g:0.1,h:0.1)[&height=0.1,height_index=0]:0.2,(e:0.2,f:0.2)[&height=0.2,height_index=1]:0.1)[&height=0.3,height_index=2]:0.1)[&height=0.4,height_index=3];";
        BaseTree<Node> tree_order_2(newick_tree_str);

        REQUIRE(tree_order_1.get_splits() == tree_order_2.get_splits());
        REQUIRE(tree_order_1.get_splits_by_height_index() != tree_order_2.get_splits_by_height_index());
    }
}
