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
        REQUIRE(s10.get_leaf_node_count() == 1);
        REQUIRE(! s10.is_empty());

        Split s01;
        s01.resize(2);
        s01.set_leaf_bit(1);
        REQUIRE(s01.get_leaf_bit(0) == 0);
        REQUIRE(s01.get_leaf_bit(1) == 1);
        REQUIRE(s01.as_string() == "01");
        expected_leaf_indices = {1};
        REQUIRE(s01.get_leaf_indices() == expected_leaf_indices);
        REQUIRE(s01.get_leaf_node_count() == 1);
        REQUIRE(! s01.is_empty());

        Split s00;
        s00.resize(2);
        REQUIRE(s00.get_leaf_bit(0) == 0);
        REQUIRE(s00.get_leaf_bit(1) == 0);
        REQUIRE(s00.as_string() == "00");
        expected_leaf_indices = {};
        REQUIRE(s00.get_leaf_indices() == expected_leaf_indices);
        REQUIRE(s00.get_leaf_node_count() == 0);
        REQUIRE(s00.is_empty());

        Split s11;
        s11.resize(2);
        s11.set_leaf_bit(0);
        s11.set_leaf_bit(1);
        REQUIRE(s11.get_leaf_bit(0) == 1);
        REQUIRE(s11.get_leaf_bit(1) == 1);
        REQUIRE(s11.as_string() == "11");
        expected_leaf_indices = {0, 1};
        REQUIRE(s11.get_leaf_indices() == expected_leaf_indices);
        REQUIRE(s11.get_leaf_node_count() == 2);
        REQUIRE(! s11.is_empty());


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
        REQUIRE(s00.is_proper_subset_of(s01));
        REQUIRE(s01.is_proper_superset_of(s00));
        REQUIRE(s00.is_proper_subset_of(s10));
        REQUIRE(s10.is_proper_superset_of(s00));

        // "proper" disqualifies equivalent sets from being considered
        // sub/supersets of each other
        REQUIRE(! ss10.is_proper_subset_of(s10));
        REQUIRE(! ss10.is_proper_superset_of(s10));
        REQUIRE(! s10.is_proper_subset_of(ss10));
        REQUIRE(! s10.is_proper_superset_of(ss10));
    }
}

TEST_CASE("Testing split sorting", "[split]") {
    SECTION("Testing sorting") {
        Split
        s1000;
        s1000.resize(4);
        s1000.set_leaf_bit(0);
        Split
        s0100;
        s0100.resize(4);
        s0100.set_leaf_bit(1);
        Split
        s0010;
        s0010.resize(4);
        s0010.set_leaf_bit(2);
        Split
        s0001;
        s0001.resize(4);
        s0001.set_leaf_bit(3);

        Split
        s1100;
        s1100.resize(4);
        s1100.set_leaf_bit(0);
        s1100.set_leaf_bit(1);
        Split
        s1110;
        s1110.resize(4);
        s1110.set_leaf_bit(0);
        s1110.set_leaf_bit(1);
        s1110.set_leaf_bit(2);
        Split
        s1111;
        s1111.resize(4);
        s1111.set_leaf_bit(0);
        s1111.set_leaf_bit(1);
        s1111.set_leaf_bit(2);
        s1111.set_leaf_bit(3);

        Split
        s0011;
        s0011.resize(4);
        s0011.set_leaf_bit(2);
        s0011.set_leaf_bit(3);
        Split
        s0111;
        s0111.resize(4);
        s0111.set_leaf_bit(1);
        s0111.set_leaf_bit(2);
        s0111.set_leaf_bit(3);

        std::vector<Split> splits {
            s1000,
            s0100,
            s1100,
            s0010,
            s1110,
            s0001,
            s1111
        };
        std::vector<Split> sorted_splits = splits;
        std::sort(std::begin(sorted_splits), std::end(sorted_splits));
        std::cout << "Sorted:\n";
        for (auto s : sorted_splits) {
            std::cout << s.as_string() << "\n";
        }
        REQUIRE(splits == sorted_splits);

        splits = {
            s1000,
            s0100,
            s0010,
            s0001,
            s0011,
            s0111,
            s1111
        };
        sorted_splits = splits;
        std::sort(std::begin(sorted_splits), std::end(sorted_splits));
        std::cout << "Sorted:\n";
        for (auto s : sorted_splits) {
            std::cout << s.as_string() << "\n";
        }
        REQUIRE(splits == sorted_splits);
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
        REQUIRE(s10100.get_leaf_node_count() == 2);
        REQUIRE(! s10100.is_empty());

        Split s01010;
        s01010.resize(5);
        s01010.set_leaf_bit(1);
        s01010.set_leaf_bit(3);
        REQUIRE(s01010.as_string() == "01010");
        expected_leaf_indices = {1, 3};
        REQUIRE(s01010.get_leaf_indices() == expected_leaf_indices);
        REQUIRE(s01010.get_leaf_node_count() == 2);
        REQUIRE(! s01010.is_empty());

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
        REQUIRE(s01011.get_leaf_node_count() == 3);

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

        s11000.clear();
        REQUIRE(s11000.is_empty());
    }
}

TEST_CASE("Testing 100000 leaves", "[split]") {
    SECTION("Testing 100000 leaves") {
        unsigned int nleaves = 100000;
        Split s;
        s.resize(nleaves);
        REQUIRE(s.is_empty());
        s.set_leaf_bit(0);
        REQUIRE(! s.is_empty());
        s.set_leaf_bit(nleaves - 1);
        REQUIRE(! s.is_empty());

        std::ostringstream expected;
        expected << "1";
        for (unsigned int i = 0; i < nleaves - 2; ++i) {
            expected << "0";
        }
        expected << "1";
        REQUIRE(s.as_string() == expected.str());

        s.clear();
        REQUIRE(s.is_empty());
    }
}

TEST_CASE("Testing is_parent_of", "[split]") {
    SECTION("Testing is_parent_of") {
        std::set<Split> splits;
        Split s00000;
        s00000.resize(5);

        Split s11000;
        s11000.resize(5);
        s11000.set_leaf_bit(0);
        s11000.set_leaf_bit(1);
        Split s00101;
        s00101.resize(5);
        s00101.set_leaf_bit(2);
        s00101.set_leaf_bit(4);
        Split s11101;
        s11101.resize(5);
        s11101.set_leaf_bit(0);
        s11101.set_leaf_bit(1);
        s11101.set_leaf_bit(2);
        s11101.set_leaf_bit(4);

        splits = {
            s11000,
            s00101};
        REQUIRE(s11101.is_parent_of(splits));
        splits = {
                  s11101,
                  s11000};
        REQUIRE(! s11101.is_parent_of(splits));
        splits = {
                  s11101,
                  s00101};
        REQUIRE(! s11000.is_parent_of(splits));

        Split s10000;
        s10000.resize(5);
        s10000.set_leaf_bit(0);
        Split s01000;
        s01000.resize(5);
        s01000.set_leaf_bit(1);
        Split s00100;
        s00100.resize(5);
        s00100.set_leaf_bit(2);
        Split s00010;
        s00010.resize(5);
        s00010.set_leaf_bit(3);
        Split s00001;
        s00001.resize(5);
        s00001.set_leaf_bit(4);

        splits = {
                s10000,
                s01000};
        REQUIRE(s11000.is_parent_of(splits));
        splits = {
                s00100,
                s00001};
        REQUIRE(s00101.is_parent_of(splits));
        splits = {
                s10000,
                s01000,
                s00100,
                s00001};
        REQUIRE(s11101.is_parent_of(splits));
        // Can't have empty split
        splits = {
                  s00000,
                  s10000,
                  s01000,
                  s00100,
                  s00001};
        REQUIRE(! s11101.is_parent_of(splits));
        // Can't have overlap
        splits = {
                  s10000,
                  s11000,
                  s01000,
                  s00100,
                  s00001};
        REQUIRE(! s11101.is_parent_of(splits));

        Split s01101;
        s01101.resize(5);
        s01101.set_leaf_bit(1);
        s01101.set_leaf_bit(2);
        s01101.set_leaf_bit(4);
        // Can't have overlap
        splits = {
                  s11000,
                  s01101};
        REQUIRE(! s11101.is_parent_of(splits));
        // Can't miss elements
        splits = {
                  s10000,
                  s00101};
        REQUIRE(! s11101.is_parent_of(splits));
        splits = {
                s10000,
                s01000,
                s00101};
        REQUIRE(s11101.is_parent_of(splits));
    }
}

TEST_CASE("Testing get_parent_of", "[split]") {
    SECTION("Testing get_parent_of") {
        std::set<Split> splits;

        Split s00000;
        s00000.resize(5);

        Split s11000;
        s11000.resize(5);
        s11000.set_leaf_bit(0);
        s11000.set_leaf_bit(1);
        Split s00101;
        s00101.resize(5);
        s00101.set_leaf_bit(2);
        s00101.set_leaf_bit(4);
        Split s11101;
        s11101.resize(5);
        s11101.set_leaf_bit(0);
        s11101.set_leaf_bit(1);
        s11101.set_leaf_bit(2);
        s11101.set_leaf_bit(4);

        splits = {
                s11000,
                s00101};
        REQUIRE(s11101 == Split::get_parent_of(splits));
        splits = {
                  s11101,
                  s11000};
        REQUIRE_THROWS_AS(Split::get_parent_of(splits),
                  EcoevolityError &);
        splits = {
                  s11101,
                  s00101};
        REQUIRE_THROWS_AS(Split::get_parent_of(splits),
                  EcoevolityError &);

        Split s10000;
        s10000.resize(5);
        s10000.set_leaf_bit(0);
        Split s01000;
        s01000.resize(5);
        s01000.set_leaf_bit(1);
        Split s00100;
        s00100.resize(5);
        s00100.set_leaf_bit(2);
        Split s00010;
        s00010.resize(5);
        s00010.set_leaf_bit(3);
        Split s00001;
        s00001.resize(5);
        s00001.set_leaf_bit(4);

        splits = {
                s10000,
                s01000};
        REQUIRE(s11000 == Split::get_parent_of(splits));
        splits = {
                s00100,
                s00001};
        REQUIRE(s00101 == Split::get_parent_of(splits));
        splits = {
                s10000,
                s01000,
                s00100,
                s00001};
        REQUIRE(s11101 == Split::get_parent_of(splits));
        // Can't have overlap
        splits = {
                  s10000,
                  s11000,
                  s01000,
                  s00100,
                  s00001};
        REQUIRE_THROWS_AS(Split::get_parent_of(splits),
                  EcoevolityError &);

        Split s01101;
        s01101.resize(5);
        s01101.set_leaf_bit(1);
        s01101.set_leaf_bit(2);
        s01101.set_leaf_bit(4);
        // Can't have overlap
        splits = {
                  s11000,
                  s01101};
        REQUIRE_THROWS_AS(Split::get_parent_of(splits),
                  EcoevolityError &);
        // Can't miss elements
        splits = {
                  s10000,
                  s00101};
        REQUIRE(s11101 != Split::get_parent_of(splits));
        splits = {
                s10000,
                s01000,
                s00101};
        REQUIRE(s11101 == Split::get_parent_of(splits));
        Split s10101;
        s10101.resize(5);
        s10101.set_leaf_bit(0);
        s10101.set_leaf_bit(2);
        s10101.set_leaf_bit(4);
        splits = {
                  s10000,
                  s00101};
        REQUIRE(s10101 == Split::get_parent_of(splits));
    }
}

TEST_CASE("Testing overlaps_with", "[split]") {
    SECTION("Testing overlaps_with") {
        Split s11000;
        s11000.resize(5);
        s11000.set_leaf_bit(0);
        s11000.set_leaf_bit(1);
        Split s00101;
        s00101.resize(5);
        s00101.set_leaf_bit(2);
        s00101.set_leaf_bit(4);
        Split s01001;
        s01001.resize(5);
        s01001.set_leaf_bit(1);
        s01001.set_leaf_bit(4);

        REQUIRE(s11000.overlaps_with(s11000));
        REQUIRE(! s11000.overlaps_with(s00101));
        REQUIRE(s11000.overlaps_with(s01001));

        REQUIRE(! s00101.overlaps_with(s11000));
        REQUIRE(s00101.overlaps_with(s00101));
        REQUIRE(s00101.overlaps_with(s01001));

        REQUIRE(s01001.overlaps_with(s11000));
        REQUIRE(s01001.overlaps_with(s00101));
        REQUIRE(s01001.overlaps_with(s01001));
    }
}

TEST_CASE("Testing can_be_siblings", "[split]") {
    SECTION("Testing can_be_siblings") {
        std::set<Split> splits;

        Split s00000;
        s00000.resize(5);

        Split s11000;
        s11000.resize(5);
        s11000.set_leaf_bit(0);
        s11000.set_leaf_bit(1);
        Split s00101;
        s00101.resize(5);
        s00101.set_leaf_bit(2);
        s00101.set_leaf_bit(4);
        Split s11101;
        s11101.resize(5);
        s11101.set_leaf_bit(0);
        s11101.set_leaf_bit(1);
        s11101.set_leaf_bit(2);
        s11101.set_leaf_bit(4);

        splits = {
                s11000,
                s00101};
        REQUIRE(Split::can_be_siblings(splits));
        splits = {
                  s11101,
                  s11000};
        REQUIRE(! Split::can_be_siblings(splits));
        splits = {
                  s11101,
                  s00101};
        REQUIRE(! Split::can_be_siblings(splits));

        Split s10000;
        s10000.resize(5);
        s10000.set_leaf_bit(0);
        Split s01000;
        s01000.resize(5);
        s01000.set_leaf_bit(1);
        Split s00100;
        s00100.resize(5);
        s00100.set_leaf_bit(2);
        Split s00010;
        s00010.resize(5);
        s00010.set_leaf_bit(3);
        Split s00001;
        s00001.resize(5);
        s00001.set_leaf_bit(4);

        splits = {
                s10000,
                s01000};
        REQUIRE(Split::can_be_siblings(splits));
        splits = {
                s00100,
                s00001};
        REQUIRE(Split::can_be_siblings(splits));
        splits = {
                s10000,
                s01000,
                s00100,
                s00001};
        REQUIRE(Split::can_be_siblings(splits));
        // Can't have empty split 
        splits = {
                  s00000,
                  s10000,
                  s01000,
                  s00100,
                  s00001};
        REQUIRE(! Split::can_be_siblings(splits));
        // Can't have overlap
        splits = {
                  s10000,
                  s11000,
                  s01000,
                  s00100,
                  s00001};
        REQUIRE(! Split::can_be_siblings(splits));

        Split s01101;
        s01101.resize(5);
        s01101.set_leaf_bit(1);
        s01101.set_leaf_bit(2);
        s01101.set_leaf_bit(4);
        // Can't have overlap
        splits = {
                  s11000,
                  s01101};
        REQUIRE(! Split::can_be_siblings(splits));
        splits = {
                  s10000,
                  s00101};
        REQUIRE(Split::can_be_siblings(splits));
        splits = {
                s10000,
                s01000,
                s00101};
        REQUIRE(Split::can_be_siblings(splits));
        Split s10101;
        s11101.resize(5);
        s11101.set_leaf_bit(0);
        s11101.set_leaf_bit(2);
        s11101.set_leaf_bit(4);
        splits = {
                  s10000,
                  s00101};
        REQUIRE(Split::can_be_siblings(splits));
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
