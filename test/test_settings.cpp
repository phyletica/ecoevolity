#include "catch.hpp"
#include "ecoevolity/settings.hpp"

TEST_CASE("Testing bare constructor of ContinuousDistributionSettings",
        "[ContinuousDistributionSettings]") {

    SECTION("Testing bare") {
        ContinuousDistributionSettings cds;
        REQUIRE(cds.get_name() == "none");
        REQUIRE(cds.to_string() == "");
        REQUIRE_THROWS_AS(cds.get_instance(),
                ContinuousDistributionSettingError);

        cds = ContinuousDistributionSettings();
        REQUIRE(cds.get_name() == "none");
        REQUIRE(cds.to_string() == "");
        REQUIRE_THROWS_AS(cds.get_instance(),
                ContinuousDistributionSettingError);
    }
}

TEST_CASE("Testing invalid distribution in ContinuousDistributionSettings constructor",
        "[ContinuousDistributionSettings]") {

    SECTION("Testing invalid distribution") {
        std::string name = "blah_distribution";
        std::unordered_map<std::string, double> parameters;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                ContinuousDistributionSettingError);
    }
}

TEST_CASE("Testing gamma settings", "[ContinuousDistributionSettings]") {

    SECTION("Testing missing gamma parameters") {
        std::string name = "gamma_distribution";
        std::unordered_map<std::string, double> parameters;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                ContinuousDistributionSettingError);
    }
}
