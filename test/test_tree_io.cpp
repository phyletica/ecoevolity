#include "catch.hpp"
#include "ecoevolity/basetree.hpp"
#include "ecoevolity/node.hpp"

TEST_CASE("Testing 2 leaves", "[split]") {
    SECTION("Testing 2 leaves") {
        const std::string newick_tree_str = "((spa[&length=0.1,height=0.0,pop_size=0.001]:0.1,spb[&length=0.1,height=0.0,pop_size=0.002]:0.1)[&length=0.2,support=1.0,height=0.1,pop_size=0.003]:0.2,spc[&length=0.3,height=0.0,pop_size=0.004]:0.3)[&height=0.3,support=1.0,pop_size=0.005];";
        BaseTree<Node> tree(newick_tree_str);
    }
}
