#include "catch.hpp"
#include "ecoevolity/tree.hpp"

TEST_CASE("Testing constructor of PopulationTree", "[PopulationTree]") {

    SECTION("Testing constructor") {
        //std::string nex_path = "data/simple.nex";
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        PopulationTree tree(nex_path, '_', true, true);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-31.77866581319647));
    }
}
