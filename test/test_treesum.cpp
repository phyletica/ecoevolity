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
        // Should get error for mismatching tip label
        /* ts.add_trees(source_tree_paths.at(5), "nexus"); */
    }
}

// Test extra tip in source tree
// Test missing tip in source tree
// Test bad label in target tree
// Test extra tip in target tree
// Test mssing tip in target tree
