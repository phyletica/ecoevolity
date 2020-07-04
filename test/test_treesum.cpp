#include "catch.hpp"
#include "ecoevolity/treesum.hpp"

#include <limits>
#include <memory>

#include "ecoevolity/tree.hpp"
#include "ecoevolity/node.hpp"

TEST_CASE("Testing bare BaseSample", "[treesum]") {
    SECTION("Testing bare BaseSample") {
        treesum::BaseSamples bs = treesum::BaseSamples();
        REQUIRE(bs.get_sample_size() == 0);
    }
}

TEST_CASE("Bad label in source tree", "[treesum]") {
    SECTION("Bad label in source tree") {
        std::vector<std::string> source_tree_paths {
                "data/4-tip-trees-12-34.nex",
                "data/4-tip-trees-12.nex",
                "data/4-tip-trees-13-24.nex",
                "data/4-tip-trees-14-23-shared.nex",
                "data/4-tip-trees-34.nex",
                "data/4-tip-trees-comb-alt-tip-label.nex",
                "data/4-tip-trees-ladder-1234.nex",
                "data/4-tip-trees-ladder-4321.nex"
        };
        std::string target_tree_path = "data/4-tip-target-tree-14-23-shared.nex";

        treesum::TreeSample<PopulationNode> ts;
        ts.add_trees(source_tree_paths.at(0), "nexus");
        ts.add_trees(source_tree_paths.at(1), "nexus");
        ts.add_trees(source_tree_paths.at(2), "nexus");
        ts.add_trees(source_tree_paths.at(3), "nexus");
        ts.add_trees(source_tree_paths.at(4), "nexus");
        ts.add_trees(source_tree_paths.at(6), "nexus");
        ts.add_trees(source_tree_paths.at(7), "nexus");
        // Should get error for mismatching tip labels
        REQUIRE_THROWS_AS(ts.add_trees(source_tree_paths.at(5), "nexus"), EcoevolityError &);

        REQUIRE_THROWS_AS(
                treesum::TreeSample<PopulationNode>(
                        source_tree_paths,
                        "nexus",
                        0),
                EcoevolityError &);
        REQUIRE_THROWS_AS(
                treesum::TreeSample<PopulationNode>(
                        target_tree_path,
                        source_tree_paths,
                        "nexus",
                        0),
                EcoevolityError &);
    }
}

TEST_CASE("Extra label in source tree", "[treesum]") {
    SECTION("Extra label in source tree") {
        std::vector<std::string> source_tree_paths {
                "data/4-tip-trees-12-34.nex",
                "data/4-tip-trees-12.nex",
                "data/4-tip-trees-13-24.nex",
                "data/4-tip-trees-14-23-shared.nex",
                "data/4-tip-trees-34.nex",
                "data/4-tip-trees-comb.nex",
                "data/4-tip-trees-ladder-1234.nex",
                "data/4-tip-trees-ladder-4321.nex",
                "data/5-tip-trees-comb.nex"
        };
        std::string target_tree_path = "data/4-tip-target-tree-14-23-shared.nex";

        treesum::TreeSample<PopulationNode> ts;
        ts.add_trees(source_tree_paths.at(0), "nexus");
        ts.add_trees(source_tree_paths.at(1), "nexus");
        ts.add_trees(source_tree_paths.at(2), "nexus");
        ts.add_trees(source_tree_paths.at(3), "nexus");
        ts.add_trees(source_tree_paths.at(4), "nexus");
        ts.add_trees(source_tree_paths.at(5), "nexus");
        ts.add_trees(source_tree_paths.at(6), "nexus");
        ts.add_trees(source_tree_paths.at(7), "nexus");
        // Should get error for mismatching tip labels
        REQUIRE_THROWS_AS(ts.add_trees(source_tree_paths.at(8), "nexus"), EcoevolityError &);

        REQUIRE_THROWS_AS(
                treesum::TreeSample<PopulationNode>(
                        source_tree_paths,
                        "nexus",
                        0),
                EcoevolityError &);
        REQUIRE_THROWS_AS(
                treesum::TreeSample<PopulationNode>(
                        target_tree_path,
                        source_tree_paths,
                        "nexus",
                        0),
                EcoevolityError &);
    }
}


TEST_CASE("Missing label in source tree", "[treesum]") {
    SECTION("Missing label in source tree") {
        std::vector<std::string> source_tree_paths {
                "data/4-tip-trees-12-34.nex",
                "data/4-tip-trees-12.nex",
                "data/4-tip-trees-13-24.nex",
                "data/4-tip-trees-14-23-shared.nex",
                "data/4-tip-trees-34.nex",
                "data/4-tip-trees-comb.nex",
                "data/4-tip-trees-ladder-1234.nex",
                "data/4-tip-trees-ladder-4321.nex",
                "data/3-tip-trees-comb.nex"
        };
        std::string target_tree_path = "data/4-tip-target-tree-14-23-shared.nex";

        treesum::TreeSample<PopulationNode> ts;
        ts.add_trees(source_tree_paths.at(0), "nexus");
        ts.add_trees(source_tree_paths.at(1), "nexus");
        ts.add_trees(source_tree_paths.at(2), "nexus");
        ts.add_trees(source_tree_paths.at(3), "nexus");
        ts.add_trees(source_tree_paths.at(4), "nexus");
        ts.add_trees(source_tree_paths.at(5), "nexus");
        ts.add_trees(source_tree_paths.at(6), "nexus");
        ts.add_trees(source_tree_paths.at(7), "nexus");
        // Should get error for mismatching tip labels
        REQUIRE_THROWS_AS(ts.add_trees(source_tree_paths.at(8), "nexus"), EcoevolityError &);

        REQUIRE_THROWS_AS(
                treesum::TreeSample<PopulationNode>(
                        source_tree_paths,
                        "nexus",
                        0),
                EcoevolityError &);
        REQUIRE_THROWS_AS(
                treesum::TreeSample<PopulationNode>(
                        target_tree_path,
                        source_tree_paths,
                        "nexus",
                        0),
                EcoevolityError &);
    }
}

TEST_CASE("Bad label in target tree", "[treesum]") {
    SECTION("Bad label in target tree") {
        std::vector<std::string> source_tree_paths {
                "data/4-tip-trees-12-34.nex",
                "data/4-tip-trees-12.nex",
                "data/4-tip-trees-13-24.nex",
                "data/4-tip-trees-14-23-shared.nex",
                "data/4-tip-trees-34.nex",
                "data/4-tip-trees-comb.nex",
                "data/4-tip-trees-ladder-1234.nex",
                "data/4-tip-trees-ladder-4321.nex",
        };
        std::string target_tree_path = "data/4-tip-target-tree-14-23-shared-alt-tip-name.nex";

        treesum::TreeSample<PopulationNode> ts;
        ts.add_trees(source_tree_paths.at(0), "nexus");
        ts.add_trees(source_tree_paths.at(1), "nexus");
        ts.add_trees(source_tree_paths.at(2), "nexus");
        ts.add_trees(source_tree_paths.at(3), "nexus");
        ts.add_trees(source_tree_paths.at(4), "nexus");
        ts.add_trees(source_tree_paths.at(5), "nexus");
        ts.add_trees(source_tree_paths.at(6), "nexus");
        ts.add_trees(source_tree_paths.at(7), "nexus");

        // Should get error for mismatching tip labels
        REQUIRE_THROWS_AS(
                ts.set_target_tree(target_tree_path, "nexus"),
                EcoevolityError &);

        // Source trees should not cause error
        ts = treesum::TreeSample<PopulationNode>(
                source_tree_paths,
                "nexus",
                0);
        REQUIRE_THROWS_AS(
                treesum::TreeSample<PopulationNode>(
                        target_tree_path,
                        source_tree_paths,
                        "nexus",
                        0),
                EcoevolityError &);
    }
}

TEST_CASE("Extra label in target tree", "[treesum]") {
    SECTION("Extra label in target tree") {
        std::vector<std::string> source_tree_paths {
                "data/4-tip-trees-12-34.nex",
                "data/4-tip-trees-12.nex",
                "data/4-tip-trees-13-24.nex",
                "data/4-tip-trees-14-23-shared.nex",
                "data/4-tip-trees-34.nex",
                "data/4-tip-trees-comb.nex",
                "data/4-tip-trees-ladder-1234.nex",
                "data/4-tip-trees-ladder-4321.nex",
        };
        std::string target_tree_path = "data/5-tip-target-tree-comb.nex";

        treesum::TreeSample<PopulationNode> ts;
        ts.add_trees(source_tree_paths.at(0), "nexus");
        ts.add_trees(source_tree_paths.at(1), "nexus");
        ts.add_trees(source_tree_paths.at(2), "nexus");
        ts.add_trees(source_tree_paths.at(3), "nexus");
        ts.add_trees(source_tree_paths.at(4), "nexus");
        ts.add_trees(source_tree_paths.at(5), "nexus");
        ts.add_trees(source_tree_paths.at(6), "nexus");
        ts.add_trees(source_tree_paths.at(7), "nexus");

        // Should get error for mismatching tip labels
        REQUIRE_THROWS_AS(
                ts.set_target_tree(target_tree_path, "nexus"),
                EcoevolityError &);

        // Source trees should not cause error
        ts = treesum::TreeSample<PopulationNode>(
                source_tree_paths,
                "nexus",
                0);
        REQUIRE_THROWS_AS(
                treesum::TreeSample<PopulationNode>(
                        target_tree_path,
                        source_tree_paths,
                        "nexus",
                        0),
                EcoevolityError &);
    }
}

TEST_CASE("Missing label in target tree", "[treesum]") {
    SECTION("Missing label in target tree") {
        std::vector<std::string> source_tree_paths {
                "data/4-tip-trees-12-34.nex",
                "data/4-tip-trees-12.nex",
                "data/4-tip-trees-13-24.nex",
                "data/4-tip-trees-14-23-shared.nex",
                "data/4-tip-trees-34.nex",
                "data/4-tip-trees-comb.nex",
                "data/4-tip-trees-ladder-1234.nex",
                "data/4-tip-trees-ladder-4321.nex",
        };
        std::string target_tree_path = "data/3-tip-target-tree-comb.nex";

        treesum::TreeSample<PopulationNode> ts;
        ts.add_trees(source_tree_paths.at(0), "nexus");
        ts.add_trees(source_tree_paths.at(1), "nexus");
        ts.add_trees(source_tree_paths.at(2), "nexus");
        ts.add_trees(source_tree_paths.at(3), "nexus");
        ts.add_trees(source_tree_paths.at(4), "nexus");
        ts.add_trees(source_tree_paths.at(5), "nexus");
        ts.add_trees(source_tree_paths.at(6), "nexus");
        ts.add_trees(source_tree_paths.at(7), "nexus");

        // Should get error for mismatching tip labels
        REQUIRE_THROWS_AS(
                ts.set_target_tree(target_tree_path, "nexus"),
                EcoevolityError &);

        // Source trees should not cause error
        ts = treesum::TreeSample<PopulationNode>(
                source_tree_paths,
                "nexus",
                0);
        REQUIRE_THROWS_AS(
                treesum::TreeSample<PopulationNode>(
                        target_tree_path,
                        source_tree_paths,
                        "nexus",
                        0),
                EcoevolityError &);
    }
}


        std::vector<double> expected_12_lengths = {
            0.2,
            0.1,
            0.1,
            0.1,
            0.2,
            0.1,
            0.1,
        };
// 12-34
// TREE t3 = [&R] ((sp1[&height=0,pop_size=0.1]:0.1,sp2[&height=0,pop_size=0.2]:0.1)[&height_index=0,height=0.1,pop_size=0.3]:0.2,(sp3[&height=0,pop_size=0.1]:0.2,sp4[&height=0,pop_size=0.2]:0.2)[&height_index=1,height=0.2,pop_size=0.3]:0.1)[&height_index=2,height=0.3,pop_size=0.3]:0.0;
// TREE t4 = [&R] ((sp4[&height=0,pop_size=0.1]:0.1,sp3[&height=0,pop_size=0.2]:0.1)[&height_index=0,height=0.1,pop_size=0.3]:0.2,(sp2[&height=0,pop_size=0.1]:0.2,sp1[&height=0,pop_size=0.2]:0.2)[&height_index=1,height=0.2,pop_size=0.3]:0.1)[&height_index=2,height=0.3,pop_size=0.3]:0.0;
// 12
// TREE t3 = [&R] ((sp1[&height=0,pop_size=0.4]:0.4,sp2[&height=0,pop_size=0.5]:0.4)[&height_index=0,height=0.4,pop_size=0.6]:0.1,sp3[&height=0,pop_size=0.4]:0.5,sp4[&height=0,pop_size=0.5]:0.5)[&height_index=1,height=0.5,pop_size=0.6]:0.0;
// TREE t4 = [&R] ((sp2[&height=0,pop_size=0.4]:0.3,sp1[&height=0,pop_size=0.5]:0.3)[&height_index=0,height=0.3,pop_size=0.6]:0.1,sp4[&height=0,pop_size=0.4]:0.4,sp3[&height=0,pop_size=0.5]:0.4)[&height_index=1,height=0.4,pop_size=0.4]:0.0;
// TREE t5 = [&R] ((sp2[&height=0,pop_size=0.1]:0.1,sp1[&height=0,pop_size=0.2]:0.1)[&height_index=0,height=0.1,pop_size=0.2]:0.2,sp4[&height=0,pop_size=0.1]:0.3,sp3[&height=0,pop_size=0.3]:0.3)[&height_index=1,height=0.3,pop_size=0.1]:0.0;
// 13-24
// TREE t3 = [&R] ((sp3[&height=0,pop_size=0.1]:0.1,sp1[&height=0,pop_size=0.2]:0.1)[&height_index=0,height=0.1,pop_size=0.5]:0.2,(sp4[&height=0,pop_size=0.1]:0.2,sp2[&height=0,pop_size=0.2]:0.2)[&height_index=1,height=0.2,pop_size=0.2]:0.1)[&height_index=2,height=0.3,pop_size=0.1]:0.0;
// TREE t4 = [&R] ((sp1[&height=0,pop_size=0.1]:0.1,sp3[&height=0,pop_size=0.2]:0.1)[&height_index=0,height=0.1,pop_size=0.3]:0.2,(sp2[&height=0,pop_size=0.1]:0.2,sp4[&height=0,pop_size=0.2]:0.2)[&height_index=1,height=0.2,pop_size=0.3]:0.1)[&height_index=2,height=0.3,pop_size=0.3]:0.0;
// 14-23-shared
// TREE t3 = [&R] ((sp1[&height=0,pop_size=0.2]:0.3,sp4[&height=0,pop_size=0.1]:0.3)[&height_index=0,height=0.3,pop_size=0.1]:0.1,(sp2[&height=0,pop_size=0.2]:0.3,sp3[&height=0,pop_size=0.1]:0.3)[&height_index=0,height=0.3,pop_size=0.2]:0.1)[&height_index=1,height=0.4,pop_size=0.1]:0.0;
// TREE t4 = [&R] ((sp2[&height=0,pop_size=0.1]:0.4,sp3[&height=0,pop_size=0.2]:0.4)[&height_index=0,height=0.4,pop_size=0.1]:0.2,(sp1[&height=0,pop_size=0.1]:0.4,sp4[&height=0,pop_size=0.2]:0.4)[&height_index=0,height=0.4,pop_size=0.2]:0.2)[&height_index=1,height=0.6,pop_size=0.2]:0.0;
// TREE t5 = [&R] ((sp2[&height=0,pop_size=0.4]:0.2,sp3[&height=0,pop_size=0.2]:0.2)[&height_index=0,height=0.2,pop_size=0.3]:0.4,(sp1[&height=0,pop_size=0.1]:0.2,sp4[&height=0,pop_size=0.3]:0.2)[&height_index=0,height=0.2,pop_size=0.3]:0.4)[&height_index=1,height=0.6,pop_size=0.1]:0.0;
// TREE t6 = [&R] ((sp2[&height=0,pop_size=0.1]:0.5,sp3[&height=0,pop_size=0.5]:0.5)[&height_index=0,height=0.5,pop_size=0.3]:0.2,(sp1[&height=0,pop_size=0.1]:0.5,sp4[&height=0,pop_size=0.1]:0.5)[&height_index=0,height=0.5,pop_size=0.1]:0.2)[&height_index=1,height=0.7,pop_size=0.6]:0.0;
// 34
// TREE t3 = [&R] ((sp3[&height=0,pop_size=0.4]:0.4,sp4[&height=0,pop_size=0.5]:0.4)[&height_index=0,height=0.4,pop_size=0.6]:0.1,sp1[&height=0,pop_size=0.4]:0.5,sp2[&height=0,pop_size=0.5]:0.5)[&height_index=1,height=0.5,pop_size=0.6]:0.0;
// TREE t4 = [&R] ((sp4[&height=0,pop_size=0.4]:0.3,sp3[&height=0,pop_size=0.5]:0.3)[&height_index=0,height=0.3,pop_size=0.6]:0.1,sp2[&height=0,pop_size=0.4]:0.4,sp1[&height=0,pop_size=0.5]:0.4)[&height_index=1,height=0.4,pop_size=0.4]:0.0;
// comb
// TREE t3 = [&R] (sp1[&height=0,pop_size=0.5]:0.1,sp2[&height=0,pop_size=0.4]:0.1,sp3[&height=0,pop_size=0.5]:0.1,sp4[&height=0,pop_size=0.3]:0.1)[&height_index=1,height=0.1,pop_size=0.4]:0.0;
// TREE t4 = [&R] (sp3[&height=0,pop_size=0.5]:0.1,sp1[&height=0,pop_size=0.5]:0.1,sp2[&height=0,pop_size=0.4]:0.1,sp4[&height=0,pop_size=0.3]:0.1)[&height_index=1,height=0.1,pop_size=0.4]:0.0;
// ladder-1234
// TREE t3 = [&R] (((sp1[&height=0,pop_size=0.5]:0.1,sp2[&height=0,pop_size=0.2]:0.1)[&height_index=0,height=0.1,pop_size=0.1]:0.1,sp3[&height=0,pop_size=0.3]:0.2)[&height_index=1,height=0.2,pop_size=0.1]:0.1,sp4[&height=0,pop_size=0.1]:0.3)[&height_index=2,height=0.3,pop_size=0.2]:0.0;
// TREE t4 = [&R] (((sp2[&height=0,pop_size=0.5]:0.2,sp1[&height=0,pop_size=0.2]:0.2)[&height_index=0,height=0.2,pop_size=0.3]:0.1,sp3[&height=0,pop_size=0.1]:0.3)[&height_index=1,height=0.3,pop_size=0.2]:0.1,sp4[&height=0,pop_size=0.2]:0.4)[&height_index=2,height=0.4,pop_size=0.5]:0.0;
// ladder-4321
// TREE t3 = [&R] (((sp4[&height=0,pop_size=0.5]:0.1,sp3[&height=0,pop_size=0.2]:0.1)[&height_index=0,height=0.1,pop_size=0.1]:0.1,sp2[&height=0,pop_size=0.3]:0.2)[&height_index=1,height=0.2,pop_size=0.1]:0.1,sp1[&height=0,pop_size=0.1]:0.3)[&height_index=2,height=0.3,pop_size=0.2]:0.0;
TEST_CASE("Basic testing", "[treesum]") {
    SECTION("Basic testing") {
        std::vector<std::string> source_tree_paths {
                "data/4-tip-trees-12-34.nex",
                "data/4-tip-trees-12.nex",
                "data/4-tip-trees-13-24.nex",
                "data/4-tip-trees-14-23-shared.nex",
                "data/4-tip-trees-34.nex",
                "data/4-tip-trees-comb.nex",
                "data/4-tip-trees-ladder-1234.nex",
                "data/4-tip-trees-ladder-4321.nex",
        };
        std::string target_tree_path = "data/4-tip-target-tree-14-23-shared.nex";

        treesum::TreeSample<PopulationNode> ts(
                target_tree_path,
                source_tree_paths,
                "nexus",
                2);
        REQUIRE(ts.get_number_of_leaves() == 4);
        REQUIRE(ts.get_number_of_sources() == 8);
        REQUIRE(ts.get_sample_size() == 18);

        std::vector<unsigned int> expected_source_sample_sizes {
                2, // "data/4-tip-trees-12-34.nex"
                3, // "data/4-tip-trees-12.nex"
                2, // "data/4-tip-trees-13-24.nex"
                4, // "data/4-tip-trees-14-23-shared.nex"
                2, // "data/4-tip-trees-34.nex"
                2, // "data/4-tip-trees-comb.nex"
                2, // "data/4-tip-trees-ladder-1234.nex"
                1, // "data/4-tip-trees-ladder-4321.nex"
        };
        REQUIRE(ts.get_source_sample_sizes() == expected_source_sample_sizes);


        Split split_1;
        split_1.resize(4);
        split_1.set_leaf_bit(0);

        Split split_2;
        split_2.resize(4);
        split_2.set_leaf_bit(1);

        Split split_3;
        split_3.resize(4);
        split_3.set_leaf_bit(2);

        Split split_4;
        split_4.resize(4);
        split_4.set_leaf_bit(3);

        Split split_12;
        split_12.resize(4);
        split_12.set_leaf_bit(0);
        split_12.set_leaf_bit(1);

        Split split_34;
        split_34.resize(4);
        split_34.set_leaf_bit(2);
        split_34.set_leaf_bit(3);

        Split split_13;
        split_13.resize(4);
        split_13.set_leaf_bit(0);
        split_13.set_leaf_bit(2);

        Split split_24;
        split_24.resize(4);
        split_24.set_leaf_bit(1);
        split_24.set_leaf_bit(3);

        Split split_14;
        split_14.resize(4);
        split_14.set_leaf_bit(0);
        split_14.set_leaf_bit(3);

        Split split_23;
        split_23.resize(4);
        split_23.set_leaf_bit(1);
        split_23.set_leaf_bit(2);

        Split split_123;
        split_123.resize(4);
        split_123.set_leaf_bit(0);
        split_123.set_leaf_bit(1);
        split_123.set_leaf_bit(2);

        Split split_234;
        split_234.resize(4);
        split_234.set_leaf_bit(1);
        split_234.set_leaf_bit(2);
        split_234.set_leaf_bit(3);

        Split split_r;
        split_r.resize(4);
        split_r.set_leaf_bit(0);
        split_r.set_leaf_bit(1);
        split_r.set_leaf_bit(2);
        split_r.set_leaf_bit(3);

        std::set<Split> s_12;
        // s_12.insert("1100");
        s_12.insert(split_12);
        std::set<Split> s_34;
        // s_34.insert("0011");
        s_34.insert(split_34);
        std::set<Split> s_13;
        // s_13.insert("1010");
        s_13.insert(split_13);
        std::set<Split> s_24;
        // s_24.insert("0101");
        s_24.insert(split_24);
        std::set<Split> s_14;
        // s_14.insert("1001");
        s_14.insert(split_14);
        std::set<Split> s_23;
        // s_23.insert("0110");
        s_23.insert(split_23);
        std::set<Split> s_123;
        // s_123.insert("1110");
        s_123.insert(split_123);
        std::set<Split> s_234;
        // s_234.insert("0111");
        s_234.insert(split_234);
        std::set<Split> s_14_23;
        // s_14_23.insert("1001");
        // s_14_23.insert("0110");
        s_14_23.insert(split_14);
        s_14_23.insert(split_23);
        std::set<Split> s_r;
        // s_r.insert("1111");
        s_r.insert(split_r);
        
        std::set< std::set<Split> > t_12_34;
        t_12_34.insert(s_12);
        t_12_34.insert(s_34);
        t_12_34.insert(s_r);
        
        std::set< std::set<Split> > t_12;
        t_12.insert(s_12);
        t_12.insert(s_r);
        
        std::set< std::set<Split> > t_13_24;
        t_13_24.insert(s_13);
        t_13_24.insert(s_24);
        t_13_24.insert(s_r);
        
        std::set< std::set<Split> > t_14_23_shared;
        t_14_23_shared.insert(s_14_23);
        t_14_23_shared.insert(s_r);
        
        std::set< std::set<Split> > t_34;
        t_34.insert(s_34);
        t_34.insert(s_r);
        
        std::set< std::set<Split> > t_comb;
        t_comb.insert(s_r);
        
        std::set< std::set<Split> > t_ladder_1234;
        t_ladder_1234.insert(s_12);
        t_ladder_1234.insert(s_123);
        t_ladder_1234.insert(s_r);
        
        std::set< std::set<Split> > t_ladder_4321;
        t_ladder_4321.insert(s_34);
        t_ladder_4321.insert(s_234);
        t_ladder_4321.insert(s_r);

        std::set< std::set<Split> > t_13_missing;
        t_13_24.insert(s_13);
        t_13_24.insert(s_r);

        // Split counts
        std::map< Split, unsigned int > expected_split_count_map {
            {split_1,   18},
            {split_2,   18},
            {split_3,   18},
            {split_4,   18},
            {split_r,   18},
            {split_12,   7},
            {split_34,   5},
            {split_14,   4},
            {split_23,   4},
            {split_13,   2},
            {split_24,   2},
            {split_123,  2},
            {split_234,  1}
        };
        std::vector<unsigned int> expected_split_counts {
            18,
            18,
            18,
            18,
            18,
             7,
             5,
             4,
             4,
             2,
             2,
             2,
             1,
        };

        for (auto key_count : expected_split_count_map) {
            REQUIRE(ts.get_split_count(key_count.first) == key_count.second);
            REQUIRE(ts.get_split_frequency(key_count.first) == key_count.second / 18.0);
        }

        std::map< Split, unsigned int > split_count_map;
        std::vector<unsigned int > split_counts;
        for (auto s : ts.get_splits()) {
            REQUIRE(split_count_map.count(s->get_split()) == 0);
            split_count_map[s->get_split()] = s->get_sample_size();
            split_counts.push_back(s->get_sample_size());
        }
        REQUIRE(split_counts == expected_split_counts);
        REQUIRE(split_count_map == expected_split_count_map);

        std::map< Split, unsigned int > expected_non_trivial_split_count_map {
            {split_12,   7},
            {split_34,   5},
            {split_14,   4},
            {split_23,   4},
            {split_13,   2},
            {split_24,   2},
            {split_123,  2},
            {split_234,  1}
        };
        std::vector<unsigned int> expected_non_trivial_split_counts {
             7,
             5,
             4,
             4,
             2,
             2,
             2,
             1,
        };

        std::map< Split, unsigned int > non_trivial_split_count_map;
        std::vector<unsigned int > non_trivial_split_counts;
        for (auto s : ts.get_non_trivial_splits()) {
            REQUIRE(non_trivial_split_count_map.count(s->get_split()) == 0);
            non_trivial_split_count_map[s->get_split()] = s->get_sample_size();
            non_trivial_split_counts.push_back(s->get_sample_size());
        }
        REQUIRE(non_trivial_split_counts == expected_non_trivial_split_counts);
        REQUIRE(non_trivial_split_count_map == expected_non_trivial_split_count_map);

        // num heights counts
        std::map< unsigned int, unsigned int> expected_nheights_count_map;
        expected_nheights_count_map[2] = 9;
        expected_nheights_count_map[3] = 7;
        expected_nheights_count_map[1] = 2;

        for (auto key_count : expected_nheights_count_map) {
            // std::cout << "key: " << key_count.first << "\n";
            // std::cout << "count: " << ts.get_number_of_heights_count(key_count.first) << "\n";
            // std::cout << "freq: " << ts.get_number_of_heights_frequency(key_count.first) << "\n";
            REQUIRE(ts.get_number_of_heights_count(key_count.first) == key_count.second);
            REQUIRE(ts.get_number_of_heights_frequency(key_count.first) == Approx(key_count.second / 18.0));
        }
        
        std::vector<unsigned int> expected_nheights_counts = {9, 7, 2};

        REQUIRE(ts.get_number_of_heights_credibility_level(2) == 1.0);
        REQUIRE(ts.get_number_of_heights_credibility_level(3) == Approx(9.0/18.0));
        REQUIRE(ts.get_number_of_heights_credibility_level(1) == Approx(2.0/18.0));
        REQUIRE(ts.get_number_of_heights_credibility_level(4) == 0.0);

        std::map< unsigned int, unsigned int> nheights_count_map;
        std::vector<unsigned int> nheights_counts;
        for (auto nh : ts.get_all_numbers_of_heights()) {
            unsigned int number_of_hts = nh->get_number_of_heights();
            REQUIRE(nheights_count_map.count(number_of_hts) == 0);
            nheights_count_map[number_of_hts] = nh->get_sample_size();
            nheights_counts.push_back(nh->get_sample_size());
        }
        REQUIRE(nheights_count_map == expected_nheights_count_map);
        REQUIRE(nheights_counts == expected_nheights_counts);
        
        // topology counts
        std::map< std::set< std::set<Split> >, unsigned int> expected_topo_count_map;
        expected_topo_count_map[t_14_23_shared] = 4;
        expected_topo_count_map[t_12] = 3;
        expected_topo_count_map[t_34] = 2;
        expected_topo_count_map[t_12_34] = 2;
        expected_topo_count_map[t_13_24] = 2;
        expected_topo_count_map[t_comb] = 2;
        expected_topo_count_map[t_ladder_1234] = 2;
        expected_topo_count_map[t_ladder_4321] = 1;

        for (auto key_count : expected_topo_count_map) {
            REQUIRE(ts.get_topology_count(key_count.first) == key_count.second);
            REQUIRE(ts.get_topology_frequency(key_count.first) == Approx(key_count.second / 18.0));
        }
        
        std::vector<unsigned int> expected_topo_counts = {4,3,2,2,2,2,2,1};

        REQUIRE(ts.get_topology_credibility_level(t_14_23_shared) == 1.0);
        REQUIRE(ts.get_topology_credibility_level(t_12) == Approx(14.0/18.0));
        REQUIRE(ts.get_topology_credibility_level(t_34) == Approx(11.0/18.0));
        REQUIRE(ts.get_topology_credibility_level(t_12_34) == Approx(11.0/18.0));
        REQUIRE(ts.get_topology_credibility_level(t_13_24) == Approx(11.0/18.0));
        REQUIRE(ts.get_topology_credibility_level(t_comb) == Approx(11.0/18.0));
        REQUIRE(ts.get_topology_credibility_level(t_ladder_1234) == Approx(11.0/18.0));
        REQUIRE(ts.get_topology_credibility_level(t_ladder_4321) == Approx(1.0/18.0));
        REQUIRE(ts.get_topology_credibility_level(t_13_missing) == 0.0);

        std::map< std::set< std::set<Split> >, unsigned int> topo_count_map;
        std::vector<unsigned int> topo_counts;
        for (auto t : ts.get_topologies()) {
            std::set< std::set<Split> > tree;
            for (auto s_set : t->get_split_set()) {
                std::set<Split> splits;
                for (auto split : s_set) {
                    splits.insert(split);
                }
                tree.insert(splits);
            }
            REQUIRE(topo_count_map.count(tree) == 0);
            topo_count_map[tree] = t->get_sample_size();
            topo_counts.push_back(t->get_sample_size());
        }
        REQUIRE(topo_count_map == expected_topo_count_map);
        REQUIRE(topo_counts == expected_topo_counts);

        // Height counts
        std::map< std::set<Split> , unsigned int> expected_height_count_map;
        expected_height_count_map[s_r] = 18;
        expected_height_count_map[s_12] = 7;
        expected_height_count_map[s_34] = 5;
        expected_height_count_map[s_14_23] = 4;
        expected_height_count_map[s_13] = 2;
        expected_height_count_map[s_24] = 2;
        expected_height_count_map[s_123] = 2;
        expected_height_count_map[s_234] = 1;

        for (auto key_count : expected_height_count_map) {
            REQUIRE(ts.get_height_count(key_count.first) == key_count.second);
            REQUIRE(ts.get_height_frequency(key_count.first) == Approx(key_count.second / 18.0));
        }

        std::vector<unsigned int> expected_height_counts = {18,7,5,4,2,2,2,1};

        std::map< std::set<Split>, unsigned int> height_count_map;
        std::vector<unsigned int> height_counts;
        for (auto h : ts.get_heights()) {
            std::set< Split > split_set;
            for (auto s : h->get_split_set()) {
                split_set.insert(s);
            }
            REQUIRE(height_count_map.count(split_set) == 0);
            height_count_map[split_set] = h->get_sample_size();
            height_counts.push_back(h->get_sample_size());
        }
        REQUIRE(height_count_map == expected_height_count_map);
        REQUIRE(height_counts == expected_height_counts);

        std::vector<double> expected_root_heights = {
            0.3,
            0.3,
            0.5,
            0.4,
            0.3,
            0.3,
            0.3,
            0.4,
            0.6,
            0.6,
            0.7,
            0.5,
            0.4,
            0.1,
            0.1,
            0.3,
            0.4,
            0.3
        };
        std::vector<double> root_heights = ts.get_height(s_r)->get_heights();
        REQUIRE(expected_root_heights.size() == root_heights.size());
        REQUIRE(std::is_permutation(
                    std::begin(root_heights),
                    std::end(root_heights),
                    std::begin(expected_root_heights)));


        std::vector<double> expected_shared_heights = {
            0.3,
            0.4,
            0.2,
            0.5,
        };
        std::vector<double> shared_heights = ts.get_height(s_14_23)->get_heights();
        REQUIRE(expected_shared_heights.size() == shared_heights.size());
        REQUIRE(std::is_permutation(
                    std::begin(shared_heights),
                    std::end(shared_heights),
                    std::begin(expected_shared_heights)));

        std::vector<double> expected_12_heights = {
            0.1,
            0.2,
            0.4,
            0.3,
            0.1,
            0.1,
            0.2
        };
        std::vector<double> s_12_heights = ts.get_height(s_12)->get_heights();
        REQUIRE(expected_12_heights.size() == s_12_heights.size());
        REQUIRE(std::is_permutation(
                    std::begin(s_12_heights),
                    std::end(s_12_heights),
                    std::begin(expected_12_heights)));

        std::vector<double> t_14_23_shared_expected_root_heights = {
            0.4,
            0.6,
            0.6,
            0.7,
        };
        std::vector<double> t_14_23_shared_expected_shared_heights = {
            0.3,
            0.4,
            0.2,
            0.5,
        };

        std::vector<double> t_14_23_shared_root_heights = ts.get_topology(t_14_23_shared)->get_heights(s_r);
        std::vector<double> t_14_23_shared_shared_heights = ts.get_topology(t_14_23_shared)->get_heights(s_14_23);

        REQUIRE(t_14_23_shared_expected_root_heights.size() == t_14_23_shared_root_heights.size());
        REQUIRE(std::is_permutation(
                    std::begin(t_14_23_shared_root_heights),
                    std::end(t_14_23_shared_root_heights),
                    std::begin(t_14_23_shared_expected_root_heights)));
        REQUIRE(t_14_23_shared_expected_shared_heights.size() == t_14_23_shared_shared_heights.size());
        REQUIRE(std::is_permutation(
                    std::begin(t_14_23_shared_shared_heights),
                    std::end(t_14_23_shared_shared_heights),
                    std::begin(t_14_23_shared_expected_shared_heights)));


        // Split parameters
        root_heights = ts.get_split(split_r)->get_values("height");
        REQUIRE(expected_root_heights.size() == root_heights.size());
        REQUIRE(std::is_permutation(
                    std::begin(root_heights),
                    std::end(root_heights),
                    std::begin(expected_root_heights)));

        std::vector<double> expected_root_lengths = {
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
        };

        std::vector<double> root_lengths = ts.get_split(split_r)->get_values("length");
        REQUIRE(expected_root_lengths.size() == root_lengths.size());
        REQUIRE(std::is_permutation(
                    std::begin(root_lengths),
                    std::end(root_lengths),
                    std::begin(expected_root_lengths)));

        std::vector<double> expected_root_sizes = {
            0.3,
            0.3,
            0.6,
            0.4,
            0.1,
            0.1,
            0.3,
            0.1,
            0.2,
            0.1,
            0.6,
            0.6,
            0.4,
            0.4,
            0.4,
            0.2,
            0.5,
            0.2
        };

        std::vector<double> root_sizes = ts.get_split(split_r)->get_values("pop_size");
        REQUIRE(expected_root_sizes.size() == root_sizes.size());
        REQUIRE(std::is_permutation(
                    std::begin(root_sizes),
                    std::end(root_sizes),
                    std::begin(expected_root_sizes)));


        s_12_heights = ts.get_split(split_12)->get_values("height");
        REQUIRE(expected_12_heights.size() == s_12_heights.size());
        REQUIRE(std::is_permutation(
                    std::begin(s_12_heights),
                    std::end(s_12_heights),
                    std::begin(expected_12_heights)));

        std::vector<double> expected_12_lengths = {
            0.2,
            0.1,
            0.1,
            0.1,
            0.2,
            0.1,
            0.1,
        };
        std::vector<double> s_12_lengths = ts.get_split(split_12)->get_values("length");
        REQUIRE(expected_12_lengths.size() == s_12_lengths.size());
        // Lengths require calcs (not just parsing), so floats won't exactly
        // match
        for (unsigned int i = 0; i < s_12_lengths.size(); ++i) {
            REQUIRE(s_12_lengths.at(i) == Approx(expected_12_lengths.at(i)));
        }

        std::vector<double> expected_12_sizes = {
            0.3,
            0.3,
            0.6,
            0.6,
            0.2,
            0.1,
            0.3,
        };
        std::vector<double> s_12_sizes = ts.get_split(split_12)->get_values("pop_size");
        REQUIRE(expected_12_sizes.size() == s_12_sizes.size());
        REQUIRE(std::is_permutation(
                    std::begin(s_12_sizes),
                    std::end(s_12_sizes),
                    std::begin(expected_12_sizes)));

        // split_12  : 7
        // split_34  : 5
        // split_14  : 4
        // split_23  : 4
        // split_13  : 2
        // split_24  : 2
        // split_123 : 2
        // split_234 : 1
        std::vector<double> split_12_source_freqs {
            1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0}; // 7/18
        std::vector<double> split_34_source_freqs {
            1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0}; // 5/18
        std::vector<double> split_14_source_freqs {
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0}; // 4/18
        std::vector<double> split_23_source_freqs {
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0}; // 4/18
        SampleSummarizer<double> split_12_source_freq_summary;
        for (auto x : split_12_source_freqs) {
            split_12_source_freq_summary.add_sample(x);
        }
        SampleSummarizer<double> split_34_source_freq_summary;
        for (auto x : split_34_source_freqs) {
            split_34_source_freq_summary.add_sample(x);
        }
        SampleSummarizer<double> split_14_source_freq_summary;
        for (auto x : split_14_source_freqs) {
            split_14_source_freq_summary.add_sample(x);
        }
        SampleSummarizer<double> split_23_source_freq_summary;
        for (auto x : split_23_source_freqs) {
            split_23_source_freq_summary.add_sample(x);
        }
        SampleSummarizer<double> split_freq_stdevs;

        double asdsf = ts.get_average_std_dev_of_split_freqs(6.0/18.0);
        split_freq_stdevs.add_sample(split_12_source_freq_summary.std_dev());
        REQUIRE(asdsf == Approx(split_freq_stdevs.mean()));

        asdsf = ts.get_average_std_dev_of_split_freqs(5.0/18.0);
        split_freq_stdevs.add_sample(split_34_source_freq_summary.std_dev());
        REQUIRE(asdsf == Approx(split_freq_stdevs.mean()));

        asdsf = ts.get_average_std_dev_of_split_freqs(3.0/18.0);
        split_freq_stdevs.add_sample(split_14_source_freq_summary.std_dev());
        split_freq_stdevs.add_sample(split_23_source_freq_summary.std_dev());
        REQUIRE(asdsf == Approx(split_freq_stdevs.mean()));

        std::stringstream ss;
        ts.write_summary_of_nontrivial_split(split_12,
                ss,
                true,
                "",
                2);
        std::stringstream ess;
        ess << "leaf_indices: [0, 1]\n"
            << "count: 7\n"
            << "frequency: 0.39\n"
            << "height_n: 7\n"
            << "height_ess: 2.7\n"
            << "height_mean: 0.2\n"
            << "height_median: 0.2\n"
            << "height_std_dev: 0.12\n"
            << "height_range: [0.1, 0.4]\n"
            << "height_eti_95: [0.1, 0.38]\n"
            << "height_hpdi_95: [0.1, 0.4]\n"
            << "length_n: 7\n"
            << "length_ess: 7\n"
            << "length_mean: 0.13\n"
            << "length_median: 0.1\n"
            << "length_std_dev: 0.049\n"
            << "length_range: [0.1, 0.2]\n"
            << "length_eti_95: [0.1, 0.2]\n"
            << "length_hpdi_95: [0.1, 0.2]\n"
            << "pop_size_n: 7\n"
            << "pop_size_ess: 2.4\n"
            << "pop_size_mean: 0.34\n"
            << "pop_size_median: 0.3\n"
            << "pop_size_std_dev: 0.19\n"
            << "pop_size_range: [0.1, 0.6]\n"
            << "pop_size_eti_95: [0.12, 0.6]\n"
            << "pop_size_hpdi_95: [0.1, 0.6]\n";
        REQUIRE(ss.str() == ess.str());

        std::stringstream hs;
        ts.write_summary_of_height(s_14_23,
                hs,
                "",
                2);
        std::stringstream ehs;
        ehs << "number_of_nodes: 2\n"
            << "clades:\n"
            << "    - leaf_indices: [1, 2]\n"
            << "    - leaf_indices: [0, 3]\n"
            << "count: 4\n"
            << "frequency: 0.22\n"
            << "n: 4\n"
            << "ess: 0\n"
            << "mean: 0.35\n"
            << "median: 0.35\n"
            << "std_dev: 0.13\n"
            << "range: [0.2, 0.5]\n"
            << "eti_95: [0.21, 0.49]\n"
            << "hpdi_95: [0.2, 0.5]\n";
        REQUIRE(hs.str() == ehs.str());

        std::stringstream rhs;
        ts.write_summary_of_height(s_r,
                rhs,
                "",
                2);
        std::stringstream erhs;
        erhs << "number_of_nodes: 1\n"
            << "clades:\n"
            << "    - leaf_indices: [0, 1, 2, 3]\n"
            << "count: 18\n"
            << "frequency: 1\n"
            << "n: 18\n"
            << "ess: 4.5\n"
            << "mean: 0.38\n"
            << "median: 0.35\n"
            << "std_dev: 0.16\n"
            << "range: [0.1, 0.7]\n"
            << "eti_95: [0.1, 0.66]\n"
            << "hpdi_95: [0.1, 0.7]\n";
        REQUIRE(rhs.str() == erhs.str());

        std::stringstream tree_stream;
        ts.write_summary_of_topology(t_12,
                tree_stream,
                0.0,
                "",
                2);
        std::stringstream ets;
        ets << "count: 3\n"
            << "frequency: 0.17\n"
            << "cumulative_frequency: 0.17\n"
            << "number_of_heights: 2\n"
            << "heights:\n"
            << "    - number_of_nodes: 1\n"
            << "      clades:\n"
            << "          - leaf_indices: [0, 1]\n"
            << "      n: 3\n"
            << "      ess: 3\n"
            << "      mean: 0.27\n"
            << "      median: 0.3\n"
            << "      std_dev: 0.15\n"
            << "      range: [0.1, 0.4]\n"
            << "      eti_95: [0.11, 0.4]\n"
            << "      hpdi_95: [0.1, 0.4]\n"
            << "    - number_of_nodes: 1\n"
            << "      clades:\n"
            << "          - leaf_indices: [0, 1, 2, 3]\n"
            << "      n: 3\n"
            << "      ess: 3\n"
            << "      mean: 0.4\n"
            << "      median: 0.4\n"
            << "      std_dev: 0.1\n"
            << "      range: [0.3, 0.5]\n"
            << "      eti_95: [0.3, 0.49]\n"
            << "      hpdi_95: [0.3, 0.5]\n";
        REQUIRE(tree_stream.str() == ets.str());

        std::stringstream nohs;
        ts.write_summary_of_all_numbers_of_heights(
                nohs,
                "",
                2);
        std::stringstream enohs;
        enohs << "numbers_of_heights:\n"
              << "    -\n"
              << "      number_of_heights: 2\n"
              << "      count: 9\n"
              << "      frequency: 0.5\n"
              << "      cumulative_frequency: 0.5\n"
              << "    -\n"
              << "      number_of_heights: 3\n"
              << "      count: 7\n"
              << "      frequency: 0.39\n"
              << "      cumulative_frequency: 0.89\n"
              << "    -\n"
              << "      number_of_heights: 1\n"
              << "      count: 2\n"
              << "      frequency: 0.11\n"
              << "      cumulative_frequency: 1\n";
        REQUIRE(nohs.str() == enohs.str());

        // ts.write_summary(std::cout,
        //         false, // use_median_heights
        //         0.1,   // min_freq_for_asdsf
        //         "",    // margin
        //         2);    // precision

        std::stringstream sumstream;
        ts.write_summary(sumstream,
                false, // use_median_heights
                0.1,   // min_freq_for_asdsf
                "",    // margin
                2);    // precision
        YAML::Node summary;
        try {
            summary = YAML::Load(sumstream);
        }
        catch (...) {
            REQUIRE(0 == 1);
        }
        REQUIRE(summary["clades"]["root"]["count"].as<int>() == 18);
        REQUIRE(summary["clades"]["leaves"][0]["count"].as<int>() == 18);
        REQUIRE(summary["clades"]["leaves"][3]["count"].as<int>() == 18);
        REQUIRE(summary["clades"]["nontrivial_clades"][0]["count"].as<int>() == 7);
        REQUIRE(summary["clades"]["nontrivial_clades"][0]["frequency"].as<double>() == 0.39);
        REQUIRE(summary["clades"]["nontrivial_clades"][0]["height_mean"].as<double>() == 0.2);
        REQUIRE(summary["clades"]["nontrivial_clades"][0]["height_std_dev"].as<double>() == 0.12);
        REQUIRE(summary["clades"]["nontrivial_clades"][0]["length_mean"].as<double>() == 0.13);
        REQUIRE(summary["clades"]["nontrivial_clades"][0]["length_std_dev"].as<double>() == 0.049);
        REQUIRE(summary["clades"]["nontrivial_clades"][0]["height_eti_95"][0].as<double>() == 0.1);
        REQUIRE(summary["clades"]["nontrivial_clades"][0]["height_eti_95"][1].as<double>() == 0.38);
        REQUIRE(summary["clades"]["nontrivial_clades"][0]["pop_size_mean"].as<double>() == 0.34);
        REQUIRE(summary["clades"]["nontrivial_clades"][0]["pop_size_eti_95"][0].as<double>() == 0.12);
        REQUIRE(summary["clades"]["nontrivial_clades"][0]["pop_size_eti_95"][1].as<double>() == 0.6);
        REQUIRE(summary["heights"][2]["count"].as<int>() == 4);
        REQUIRE(summary["heights"][2]["number_of_nodes"].as<int>() == 2);
        REQUIRE(summary["heights"][2]["frequency"].as<double>() == 0.22);
        REQUIRE(summary["heights"][2]["mean"].as<double>() == 0.35);
        REQUIRE(summary["heights"][2]["std_dev"].as<double>() == 0.13);
        REQUIRE(summary["heights"][2]["eti_95"][0].as<double>() == 0.21);
        REQUIRE(summary["heights"][2]["eti_95"][1].as<double>() == 0.49);
        REQUIRE(summary["heights"][2]["clades"][0]["leaf_indices"][0].as<int>() == 1);
        REQUIRE(summary["heights"][2]["clades"][0]["leaf_indices"][1].as<int>() == 2);
        REQUIRE(summary["heights"][2]["clades"][1]["leaf_indices"][0].as<int>() == 0);
        REQUIRE(summary["heights"][2]["clades"][1]["leaf_indices"][1].as<int>() == 3);
        REQUIRE(summary["topologies"][1]["count"].as<int>() == 3);
        REQUIRE(summary["topologies"][1]["number_of_heights"].as<int>() == 2);
        REQUIRE(summary["topologies"][1]["frequency"].as<double>() == 0.17);
        REQUIRE(summary["topologies"][1]["cumulative_frequency"].as<double>() == 0.39);
        REQUIRE(summary["topologies"][1]["heights"][0]["number_of_nodes"].as<int>() == 1);
        REQUIRE(summary["topologies"][1]["heights"][1]["number_of_nodes"].as<int>() == 1);
        REQUIRE(summary["topologies"][1]["heights"][0]["clades"][0]["leaf_indices"][0].as<int>() == 0);
        REQUIRE(summary["topologies"][1]["heights"][0]["clades"][0]["leaf_indices"][1].as<int>() == 1);
        REQUIRE(summary["topologies"][1]["heights"][1]["clades"][0]["leaf_indices"][0].as<int>() == 0);
        REQUIRE(summary["topologies"][1]["heights"][1]["clades"][0]["leaf_indices"][1].as<int>() == 1);
        REQUIRE(summary["topologies"][1]["heights"][1]["clades"][0]["leaf_indices"][2].as<int>() == 2);
        REQUIRE(summary["topologies"][1]["heights"][1]["clades"][0]["leaf_indices"][3].as<int>() == 3);
        REQUIRE(summary["topologies"][1]["heights"][0]["mean"].as<double>() == 0.27);
        REQUIRE(summary["topologies"][1]["heights"][1]["mean"].as<double>() == 0.4);
        REQUIRE(summary["topologies"][1]["heights"][0]["std_dev"].as<double>() == 0.15);
        REQUIRE(summary["topologies"][1]["heights"][1]["std_dev"].as<double>() == 0.1);
        REQUIRE(summary["topologies"][1]["heights"][0]["eti_95"][0].as<double>() == 0.11);
        REQUIRE(summary["topologies"][1]["heights"][1]["eti_95"][0].as<double>() == 0.3);
        REQUIRE(summary["topologies"][1]["heights"][0]["eti_95"][1].as<double>() == 0.4);
        REQUIRE(summary["topologies"][1]["heights"][1]["eti_95"][1].as<double>() == 0.49);
        REQUIRE(summary["topologies"][1]["heights"][0]["hpdi_95"][0].as<double>() == 0.1);
        REQUIRE(summary["topologies"][1]["heights"][1]["hpdi_95"][0].as<double>() == 0.3);
        REQUIRE(summary["topologies"][1]["heights"][0]["hpdi_95"][1].as<double>() == 0.4);
        REQUIRE(summary["topologies"][1]["heights"][1]["hpdi_95"][1].as<double>() == 0.5);
        REQUIRE(summary["summary_of_map_trees"][0]["count"].as<int>() == 4);
        REQUIRE(summary["summary_of_map_trees"][0]["frequency"].as<double>() == 0.22);

        // ts.write_map_trees_to_nexus(std::cout);
        // ts.write_target_tree_to_nexus(std::cout);

        source_tree_paths = {
                "data/4-tip-trees-12-34.nex",
                "data/4-tip-trees-13-24.nex",
                "data/4-tip-trees-34.nex",
                "data/4-tip-trees-comb.nex",
                "data/4-tip-trees-ladder-1234.nex",
        };
        target_tree_path = "data/4-tip-target-tree-comb.nex";

        treesum::TreeSample<PopulationNode> ts2(
                target_tree_path,
                source_tree_paths,
                "nexus",
                1);
        REQUIRE(ts2.get_number_of_leaves() == 4);
        REQUIRE(ts2.get_number_of_sources() == 5);
        REQUIRE(ts2.get_sample_size() == 15);

        expected_source_sample_sizes = {
            3, 3, 3, 3, 3
        };
        REQUIRE(ts2.get_source_sample_sizes() == expected_source_sample_sizes);
        REQUIRE(ts2.has_multiple_sources_with_equal_n() == true);

        ts2.write_summary(std::cout,
                false, // use_median_heights
                0.1,   // min_freq_for_asdsf
                "",    // margin
                2);    // precision
    }
}
