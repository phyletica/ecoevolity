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

TEST_CASE("Testing missing exponential parameter settings", "[ContinuousDistributionSettings]") {

    SECTION("Testing no parameters") {
        std::string name = "exponential_distribution";
        std::unordered_map<std::string, double> parameters;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }

    SECTION("Testing missing parameters") {
        std::string name = "exponential_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["offset"] = 0.01;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }

    SECTION("Testing lambda parameter") {
        std::string name = "exponential_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["lambda"] = 1.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }
}

TEST_CASE("Testing gamma parameter settings errors", "[ContinuousDistributionSettings]") {
    SECTION("Testing gamma errors") {
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
}

TEST_CASE("Testing exponential parameter settings errors", "[ContinuousDistributionSettings]") {
    SECTION("Testing exponential errors") {
        std::string name = "exponential_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["rate"] = 0.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }
}

TEST_CASE("Testing gamma settings to string", "[ContinuousDistributionSettings]") {
    SECTION("Testing gamma to string") {
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

TEST_CASE("Testing exponential settings to string", "[ContinuousDistributionSettings]") {
    SECTION("Testing exponential to string") {
        std::string name = "exponential_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["rate"] = 1.0;
        ContinuousDistributionSettings settings = ContinuousDistributionSettings(name, parameters);

        REQUIRE(settings.get_name() == "exponential_distribution");
        std::string s = settings.to_string(0);
        std::string e = "exponential_distribution:\n    rate: 1\n";
        REQUIRE(s == e);

        parameters.clear();
        parameters["rate"] = 1.0;
        parameters["offset"] = 0.01;
        settings = ContinuousDistributionSettings(name, parameters);

        REQUIRE(settings.get_name() == "exponential_distribution");
        s = settings.to_string(0);
        e = "exponential_distribution:\n    rate: 1\n    offset: 0.01\n";
        REQUIRE(s == e);

        parameters.clear();
        parameters["rate"] = 0.12345678910;
        parameters["offset"] = 0.001;
        settings = ContinuousDistributionSettings(name, parameters);

        REQUIRE(settings.get_name() == "exponential_distribution");
        s = settings.to_string(0);
        e = "exponential_distribution:\n    rate: 0.123457\n    offset: 0.001\n";
        REQUIRE(s == e);
    }
}

TEST_CASE("Testing gamma settings to instance", "[ContinuousDistributionSettings]") {
    SECTION("Testing gamma to instance") {
        std::string name = "gamma_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["shape"] = 10.0;
        parameters["scale"] = 1.0;
        ContinuousDistributionSettings settings = ContinuousDistributionSettings(name, parameters);

        std::shared_ptr<ContinuousProbabilityDistribution> g;
        g = settings.get_instance();

        REQUIRE(g->get_min() == 0.0);
        REQUIRE(g->to_string() == "gamma(shape = 10, scale = 1)");


        parameters.clear();
        parameters["shape"] = 10.0;
        parameters["scale"] = 1.0;
        parameters["offset"] = 0.01;
        settings = ContinuousDistributionSettings(name, parameters);

        g = settings.get_instance();

        REQUIRE(g->to_string() == "gamma(shape = 10, scale = 1, offset = 0.01)");
    }
}

TEST_CASE("Testing exponential settings to instance", "[ContinuousDistributionSettings]") {
    SECTION("Testing exponential to instance") {
        std::string name = "exponential_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["rate"] = 1.0;
        ContinuousDistributionSettings settings = ContinuousDistributionSettings(name, parameters);

        std::shared_ptr<ContinuousProbabilityDistribution> g;
        g = settings.get_instance();

        REQUIRE(g->get_min() == 0.0);
        REQUIRE(g->to_string() == "exp(lambda = 1)");


        parameters["rate"] = 10.0;
        parameters["offset"] = 0.01;
        settings = ContinuousDistributionSettings(name, parameters);

        g = settings.get_instance();

        REQUIRE(g->to_string() == "exp(lambda = 10, offset = 0.01)");
    }
}

TEST_CASE("Testing uniform parameter settings errors", "[ContinuousDistributionSettings]") {
    SECTION("Testing uniform errors") {
        std::string name = "uniform_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["min"] = 1.0;
        parameters["max"] = 0.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }
}

TEST_CASE("Testing uniform settings to string", "[ContinuousDistributionSettings]") {
    SECTION("Testing uniform to string") {
        std::string name = "uniform_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["min"] = 1.0;
        parameters["max"] = 2.0;
        ContinuousDistributionSettings settings = ContinuousDistributionSettings(name, parameters);

        REQUIRE(settings.get_name() == "uniform_distribution");
        std::string s = settings.to_string(0);
        std::string e = "uniform_distribution:\n    min: 1\n    max: 2\n";
        REQUIRE(s == e);
    }
}

TEST_CASE("Testing uniform settings to instance", "[ContinuousDistributionSettings]") {
    SECTION("Testing uniform to instance") {
        std::string name = "uniform_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["min"] = 1.0;
        parameters["max"] = 2.0;
        ContinuousDistributionSettings settings = ContinuousDistributionSettings(name, parameters);

        std::shared_ptr<ContinuousProbabilityDistribution> g;
        g = settings.get_instance();

        REQUIRE(g->get_min() == 1.0);
        REQUIRE(g->to_string() == "uniform(1, 2)");
    }
}

TEST_CASE("Testing missing uniform parameter settings", "[ContinuousDistributionSettings]") {

    SECTION("Testing no parameters") {
        std::string name = "uniform_distribution";
        std::unordered_map<std::string, double> parameters;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }

    SECTION("Testing missing parameters") {
        std::string name = "uniform_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["min"] = 0.01;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);

        parameters.clear();
        parameters["max"] = 1.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }
}
