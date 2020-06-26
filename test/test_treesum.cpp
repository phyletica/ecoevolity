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

        std::map< std::string, unsigned int > expected_split_count_map {
            {"1000" , 18},
            {"0100" , 18},
            {"0010" , 18},
            {"0001" , 18},
            {"1111" , 18},
            {"1100" ,  7},
            {"0011" ,  5},
            {"1001" ,  4},
            {"0110" ,  4},
            {"1010" ,  2},
            {"0101" ,  2},
            {"1110" ,  2},
            {"0111" ,  1}
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

        std::map< std::string, unsigned int > split_count_map;
        std::vector<unsigned int > split_counts;
        for (auto s : ts.get_splits()) {
            REQUIRE(split_count_map.count(s->get_split_as_string()) == 0);
            split_count_map[s->get_split_as_string()] = s->get_sample_size();
            split_counts.push_back(s->get_sample_size());
        }
        REQUIRE(split_counts == expected_split_counts);
        REQUIRE(split_count_map == expected_split_count_map);

        std::map< std::string, unsigned int > expected_non_trivial_split_count_map {
            {"1100" ,  7},
            {"0011" ,  5},
            {"1001" ,  4},
            {"0110" ,  4},
            {"1010" ,  2},
            {"0101" ,  2},
            {"1110" ,  2},
            {"0111" ,  1}
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

        std::map< std::string, unsigned int > non_trivial_split_count_map;
        std::vector<unsigned int > non_trivial_split_counts;
        for (auto s : ts.get_non_trivial_splits()) {
            REQUIRE(non_trivial_split_count_map.count(s->get_split_as_string()) == 0);
            non_trivial_split_count_map[s->get_split_as_string()] = s->get_sample_size();
            non_trivial_split_counts.push_back(s->get_sample_size());
        }
        REQUIRE(non_trivial_split_counts == expected_non_trivial_split_counts);
        REQUIRE(non_trivial_split_count_map == expected_non_trivial_split_count_map);

        std::set<std::string> s_12;
        s_12.insert("1100");
        std::set<std::string> s_34;
        s_34.insert("0011");
        std::set<std::string> s_13;
        s_13.insert("1010");
        std::set<std::string> s_24;
        s_24.insert("0101");
        std::set<std::string> s_14;
        s_14.insert("1001");
        std::set<std::string> s_23;
        s_23.insert("0110");
        std::set<std::string> s_123;
        s_123.insert("1110");
        std::set<std::string> s_234;
        s_234.insert("0111");
        std::set<std::string> s_14_23;
        s_14_23.insert("1001");
        s_14_23.insert("0110");
        std::set<std::string> s_r;
        s_r.insert("1111");
        
        std::set< std::set<std::string> > t_12_34;
        t_12_34.insert(s_12);
        t_12_34.insert(s_34);
        t_12_34.insert(s_r);
        
        std::set< std::set<std::string> > t_12;
        t_12.insert(s_12);
        t_12.insert(s_r);
        
        std::set< std::set<std::string> > t_13_24;
        t_13_24.insert(s_13);
        t_13_24.insert(s_24);
        t_13_24.insert(s_r);
        
        std::set< std::set<std::string> > t_14_23_shared;
        t_14_23_shared.insert(s_14_23);
        t_14_23_shared.insert(s_r);
        
        std::set< std::set<std::string> > t_34;
        t_34.insert(s_34);
        t_34.insert(s_r);
        
        std::set< std::set<std::string> > t_comb;
        t_comb.insert(s_r);
        
        std::set< std::set<std::string> > t_ladder_1234;
        t_ladder_1234.insert(s_12);
        t_ladder_1234.insert(s_123);
        t_ladder_1234.insert(s_r);
        
        std::set< std::set<std::string> > t_ladder_4321;
        t_ladder_4321.insert(s_34);
        t_ladder_4321.insert(s_234);
        t_ladder_4321.insert(s_r);
        
        // topology counts
        std::map< std::set< std::set<std::string> >, unsigned int> expected_topo_count_map;
        expected_topo_count_map[t_14_23_shared] = 4;
        expected_topo_count_map[t_12] = 3;
        expected_topo_count_map[t_34] = 2;
        expected_topo_count_map[t_12_34] = 2;
        expected_topo_count_map[t_13_24] = 2;
        expected_topo_count_map[t_comb] = 2;
        expected_topo_count_map[t_ladder_1234] = 2;
        expected_topo_count_map[t_ladder_4321] = 1;
        
        std::vector<unsigned int> expected_topo_counts = {4,3,2,2,2,2,2,1};

        std::map< std::set< std::set<std::string> >, unsigned int> topo_count_map;
        std::vector<unsigned int> topo_counts;
        for (auto t : ts.get_topologies()) {
            std::set< std::set<std::string> > tree;
            for (auto s_set : t->get_split_set()) {
                std::set<std::string> splits;
                for (auto split : s_set) {
                    splits.insert(split.as_string());
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
        std::map< std::set<std::string> , unsigned int> expected_height_count_map;
        expected_height_count_map[s_r] = 18;
        expected_height_count_map[s_12] = 7;
        expected_height_count_map[s_34] = 5;
        expected_height_count_map[s_14_23] = 4;
        expected_height_count_map[s_13] = 2;
        expected_height_count_map[s_24] = 2;
        expected_height_count_map[s_123] = 2;
        expected_height_count_map[s_234] = 1;

        std::vector<unsigned int> expected_height_counts = {18,7,5,4,2,2,2,1};

        std::map< std::set<std::string>, unsigned int> height_count_map;
        std::vector<unsigned int> height_counts;
        for (auto h : ts.get_heights()) {
            std::set< std::string > split_set;
            for (auto s : h->get_split_set()) {
                split_set.insert(s.as_string());
            }
            REQUIRE(height_count_map.count(split_set) == 0);
            height_count_map[split_set] = h->get_sample_size();
            height_counts.push_back(h->get_sample_size());
        }
        REQUIRE(height_count_map == expected_height_count_map);
        REQUIRE(height_counts == expected_height_counts);
    }
}
