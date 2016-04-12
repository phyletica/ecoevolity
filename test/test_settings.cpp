#include "catch.hpp"
#include "ecoevolity/settings.hpp"

TEST_CASE("Testing bare constructor of ContinuousDistributionSettings",
        "[ContinuousDistributionSettings]") {

    SECTION("Testing bare") {
        ContinuousDistributionSettings cds;
        REQUIRE(cds.get_name() == "none");
        REQUIRE(cds.to_string() == "");
        REQUIRE_THROWS_AS(cds.get_instance(),
                EcoevolityContinuousDistributionSettingError);

        cds = ContinuousDistributionSettings();
        REQUIRE(cds.get_name() == "none");
        REQUIRE(cds.to_string() == "");
        REQUIRE_THROWS_AS(cds.get_instance(),
                EcoevolityContinuousDistributionSettingError);
    }
}

TEST_CASE("Testing invalid distribution in ContinuousDistributionSettings constructor",
        "[ContinuousDistributionSettings]") {

    SECTION("Testing invalid distribution") {
        std::string name = "blah_distribution";
        std::unordered_map<std::string, double> parameters;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }
}

TEST_CASE("Testing missing gamma parameter settings", "[ContinuousDistributionSettings]") {

    SECTION("Testing no parameters") {
        std::string name = "gamma_distribution";
        std::unordered_map<std::string, double> parameters;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }

    SECTION("Testing missing parameters") {
        std::string name = "gamma_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["shape"] = 1.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
        parameters["offset"] = 0.01;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);

        parameters.clear();
        parameters["scale"] = 1.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
        parameters["offset"] = 0.01;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }

    SECTION("Testing rate parameter") {
        std::string name = "gamma_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["shape"] = 10.0;
        parameters["rate"] = 1.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }
}

TEST_CASE("Testing gamma parameter settings errors", "[ContinuousDistributionSettings]") {
        std::string name = "gamma_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["shape"] = 0.0;
        parameters["scale"] = 1.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);

        parameters.clear();
        parameters["shape"] = 10.0;
        parameters["scale"] = 0.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
}

TEST_CASE("Testing gamma settings to string", "[ContinuousDistributionSettings]") {
    SECTION("Testing gamma 10, 1") {
        std::string name = "gamma_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["shape"] = 10.0;
        parameters["scale"] = 1.0;
        ContinuousDistributionSettings settings = ContinuousDistributionSettings(name, parameters);

        REQUIRE(settings.get_name() == "gamma_distribution");
        std::string s = settings.to_string(0);
        std::string e = "gamma_distribution:\n    shape: 10\n    scale: 1\n";
        REQUIRE(s == e);

        parameters.clear();
        parameters["shape"] = 10.0;
        parameters["scale"] = 1.0;
        parameters["offset"] = 0.01;
        settings = ContinuousDistributionSettings(name, parameters);

        REQUIRE(settings.get_name() == "gamma_distribution");
        s = settings.to_string(0);
        e = "gamma_distribution:\n    shape: 10\n    scale: 1\n    offset: 0.01\n";
        REQUIRE(s == e);

        parameters.clear();
        parameters["shape"] = 0.12345678910;
        parameters["scale"] = 1.00001;
        parameters["offset"] = 0.001;
        settings = ContinuousDistributionSettings(name, parameters);

        REQUIRE(settings.get_name() == "gamma_distribution");
        s = settings.to_string(0);
        e = "gamma_distribution:\n    shape: 0.123457\n    scale: 1.00001\n    offset: 0.001\n";
        REQUIRE(s == e);
    }
}
