#include "catch.hpp"
#include "ecoevolity/treecomp.hpp"

#include <limits>
#include <memory>

#include "ecoevolity/tree.hpp"
#include "ecoevolity/node.hpp"


TEST_CASE("Testing euclidean_distance with comb BaseTree<Node> trees", "[treecomp]") {
    SECTION("Testing comb trees") {
        std::string t1 = "(sp1[&height=0.0,pop_size=0.001]:2.0,sp2[&height=0.0,pop_size=0.002]:2.0,sp3[&height=0.0,pop_size=0.002]:2.0,sp4[&height=0.0,pop_size=0.002]:2.0)[&height_index=0,height=2.0,pop_size=0.003];";
        std::string t2 = "(sp1[&height=0.0,pop_size=0.001]:4.0,sp2[&height=0.0,pop_size=0.002]:4.0,sp3[&height=0.0,pop_size=0.002]:4.0,sp4[&height=0.0,pop_size=0.002]:4.0)[&height_index=0,height=4.0,pop_size=0.003];";
        BaseTree<Node> tree1(t1);
        BaseTree<Node> tree2(t2);

        REQUIRE(treecomp::euclidean_distance(tree1, tree2) == Approx(4.0));
    }
}

TEST_CASE("Testing euclidean_distance with comb BaseTree<PopulationNode> trees", "[treecomp]") {
    SECTION("Testing comb trees") {
        std::string t1 = "(sp1[&height=0.0,pop_size=0.001]:2.0,sp2[&height=0.0,pop_size=0.002]:2.0,sp3[&height=0.0,pop_size=0.002]:2.0,sp4[&height=0.0,pop_size=0.002]:2.0)[&height_index=0,height=2.0,pop_size=0.003];";
        std::string t2 = "(sp1[&height=0.0,pop_size=0.001]:4.0,sp2[&height=0.0,pop_size=0.002]:4.0,sp3[&height=0.0,pop_size=0.002]:4.0,sp4[&height=0.0,pop_size=0.002]:4.0)[&height_index=0,height=4.0,pop_size=0.003];";
        BaseTree<PopulationNode> tree1(t1);
        BaseTree<PopulationNode> tree2(t2);

        REQUIRE(treecomp::euclidean_distance(tree1, tree2) == Approx(4.0));
    }
}

TEST_CASE("Testing euclidean_distance with identical comb BaseTree<Node> trees", "[treecomp]") {
    SECTION("Testing comb trees") {
        std::string t1 = "(sp1[&height=0.0,pop_size=0.001]:2.0,sp2[&height=0.0,pop_size=0.002]:2.0,sp3[&height=0.0,pop_size=0.002]:2.0,sp4[&height=0.0,pop_size=0.002]:2.0)[&height_index=0,height=2.0,pop_size=0.003];";
        std::string t2 = "(sp4[&height=0.0,pop_size=0.001]:2.0,sp3[&height=0.0,pop_size=0.002]:2.0,sp2[&height=0.0,pop_size=0.002]:2.0,sp1[&height=0.0,pop_size=0.002]:2.0)[&height_index=0,height=2.0,pop_size=0.003];";
        BaseTree<Node> tree1(t1);
        BaseTree<Node> tree2(t2);

        REQUIRE(treecomp::euclidean_distance(tree1, tree2) == Approx(0.0));
    }
}

TEST_CASE("Testing euclidean_distance with comb vs 1-clade", "[treecomp]") {
    SECTION("Testing comb trees") {
        std::string t1 = "(sp1[&height=0.0,pop_size=0.001]:4.0,sp2[&height=0.0,pop_size=0.002]:4.0,sp3[&height=0.0,pop_size=0.002]:4.0,sp4[&height=0.0,pop_size=0.002]:4.0)[&height_index=0,height=4.0,pop_size=0.003];";
        std::string t2 = "((sp1[&height=0.0,pop_size=0.001]:2.0,sp2[&height=0.0,pop_size=0.002]:2.0)[&height_index=0,height=2.0,pop_size=0.002]:2.0,sp3[&height=0.0,pop_size=0.002]:4.0,sp4[&height=0.0,pop_size=0.002]:4.0)[&height_index=1,height=4.0,pop_size=0.003];";
        BaseTree<Node> tree1(t1);
        BaseTree<Node> tree2(t2);

        REQUIRE(treecomp::euclidean_distance(tree1, tree2, true) == Approx(std::sqrt(12.0)));
    }
}
