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
