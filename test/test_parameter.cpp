#include "catch.hpp"
#include "ecoevolity/parameter.hpp"

#include <limits>

TEST_CASE("Testing RealParameter constructors", "[RealParameter]") {

    SECTION("Testing bare") {
        RealParameter p = RealParameter();
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());
    }

    SECTION("Testing value") {
        RealParameter p = RealParameter(1.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());
        REQUIRE(p.get_value() == Approx(1.1));
        p.store();
        REQUIRE(p.get_stored_value() == Approx(1.1));
    }
}

TEST_CASE("Testing PositiveRealParameter constructors", "[PositiveRealParameter]") {

    SECTION("Testing bare") {
        PositiveRealParameter p = PositiveRealParameter();
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == 0.0);
    }

    SECTION("Testing value") {
        PositiveRealParameter p = PositiveRealParameter(1.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == 0.0);
        REQUIRE(p.get_value() == Approx(1.1));
        p.store();
        REQUIRE(p.get_stored_value() == Approx(1.1));
    }
}

TEST_CASE("Testing IntParameter constructors", "[IntParameter]") {

    SECTION("Testing bare") {
        IntParameter p = IntParameter();
        REQUIRE(p.get_max() == std::numeric_limits<int>::max());
        REQUIRE(p.get_min() == -std::numeric_limits<int>::max());
        REQUIRE(p.get_upper() == std::numeric_limits<int>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<int>::lowest());
    }

    SECTION("Testing value") {
        IntParameter p = IntParameter(2);
        REQUIRE(p.get_max() == std::numeric_limits<int>::max());
        REQUIRE(p.get_min() == -std::numeric_limits<int>::max());
        REQUIRE(p.get_upper() == std::numeric_limits<int>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<int>::lowest());
        REQUIRE(p.get_value() == 2);
        p.store();
        REQUIRE(p.get_stored_value() == 2);
    }
}

TEST_CASE("Testing Probability constructors", "[Probability]") {

    SECTION("Testing bare") {
        Probability p = Probability();
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == 1.0);
        REQUIRE(p.get_lower() == 0.0);
    }

    SECTION("Testing value") {
        Probability p = Probability(0.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == 1.0);
        REQUIRE(p.get_lower() == 0.0);
        REQUIRE(p.get_value() == 0.1);
        p.store();
        REQUIRE(p.get_stored_value() == 0.1);
    }
}

TEST_CASE("Testing ProbabilityDensity constructors", "[ProbabilityDensity]") {

    SECTION("Testing bare") {
        ProbabilityDensity p = ProbabilityDensity();
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == 0.0);
    }

    SECTION("Testing value") {
        ProbabilityDensity p = ProbabilityDensity(0.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == 0.0);
        REQUIRE(p.get_value() == 0.1);
        p.store();
        REQUIRE(p.get_stored_value() == 0.1);
    }
}

TEST_CASE("Testing LogProbability constructors", "[LogProbability]") {

    SECTION("Testing bare") {
        LogProbability p = LogProbability();
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == 0.0);
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());
    }

    SECTION("Testing value") {
        LogProbability p = LogProbability(-0.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == 0.0);
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());
        REQUIRE(p.get_value() == -0.1);
        p.store();
        REQUIRE(p.get_stored_value() == -0.1);
    }
}

TEST_CASE("Testing LogProbabilityDensity constructors", "[LogProbabilityDensity]") {

    SECTION("Testing bare") {
        LogProbabilityDensity p = LogProbabilityDensity();
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());
    }

    SECTION("Testing value") {
        LogProbabilityDensity p = LogProbabilityDensity(-0.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());
        REQUIRE(p.get_value() == -0.1);
        p.store();
        REQUIRE(p.get_stored_value() == -0.1);
    }
}

TEST_CASE("Testing RealParameter value methods", "[RealParameter]") {

    SECTION("Testing value methods") {
        RealParameter p = RealParameter(1.1);
        REQUIRE(p.get_value() == 1.1);
        p.store();
        p.set_value(-2.2);
        REQUIRE(p.get_value() == -2.2);
        REQUIRE(p.get_stored_value() == 1.1);

        p.restore();
        REQUIRE(p.get_value() == 1.1);
        REQUIRE(p.get_stored_value() == 1.1);

        p.update_value(3.0);
        REQUIRE(p.get_value() == 3.0);
        REQUIRE(p.get_stored_value() == 1.1);

        p.update_value(4.0);
        REQUIRE(p.get_value() == 4.0);
        REQUIRE(p.get_stored_value() == 3.0);

        p.restore();
        REQUIRE(p.get_value() == 3.0);
        REQUIRE(p.get_stored_value() == 3.0);
    }
}

TEST_CASE("Testing IntParameter value methods", "[IntParameter]") {

    SECTION("Testing value methods") {
        IntParameter p = IntParameter(1);
        REQUIRE(p.get_value() == 1);
        p.store();
        p.set_value(-2);
        REQUIRE(p.get_value() == -2);
        REQUIRE(p.get_stored_value() == 1);

        p.restore();
        REQUIRE(p.get_value() == 1);
        REQUIRE(p.get_stored_value() == 1);

        p.update_value(3);
        REQUIRE(p.get_value() == 3);
        REQUIRE(p.get_stored_value() == 1);

        p.update_value(4);
        REQUIRE(p.get_value() == 4);
        REQUIRE(p.get_stored_value() == 3);

        p.restore();
        REQUIRE(p.get_value() == 3);
        REQUIRE(p.get_stored_value() == 3);
    }
}

TEST_CASE("Testing Probability value methods", "[Probability]") {

    SECTION("Testing value methods") {
        REQUIRE_THROWS_AS(Probability(1.1), EcoevolityParameterValueError);
        REQUIRE_THROWS_AS(Probability(-0.1), EcoevolityParameterValueError);
        Probability p = Probability(0.1);
        REQUIRE(p.get_value() == 0.1);
        p.store();
        p.set_value(0.2);
        REQUIRE(p.get_value() == 0.2);
        REQUIRE(p.get_stored_value() == 0.1);

        p.restore();
        REQUIRE(p.get_value() == 0.1);
        REQUIRE(p.get_stored_value() == 0.1);

        p.update_value(0.3);
        REQUIRE(p.get_value() == 0.3);
        REQUIRE(p.get_stored_value() == 0.1);

        p.update_value(0.4);
        REQUIRE(p.get_value() == 0.4);
        REQUIRE(p.get_stored_value() == 0.3);

        p.restore();
        REQUIRE(p.get_value() == 0.3);
        REQUIRE(p.get_stored_value() == 0.3);

        REQUIRE_THROWS_AS(p.set_value(1.01), EcoevolityParameterValueError);
    }
}

TEST_CASE("Testing LogProbability value methods", "[LogProbability]") {

    SECTION("Testing value methods") {
        REQUIRE_THROWS_AS(LogProbability(0.1), EcoevolityParameterValueError);
        LogProbability p = LogProbability(0.0);
        REQUIRE(p.get_value() == 0.0);
        p.store();
        p.set_value(-0.2);
        REQUIRE(p.get_value() == -0.2);
        REQUIRE(p.get_stored_value() == 0.0);

        p.restore();
        REQUIRE(p.get_value() == 0.0);
        REQUIRE(p.get_stored_value() == 0.0);

        p.update_value(-0.3);
        REQUIRE(p.get_value() == -0.3);
        REQUIRE(p.get_stored_value() == 0.0);

        p.update_value(-0.4);
        REQUIRE(p.get_value() == -0.4);
        REQUIRE(p.get_stored_value() == -0.3);

        p.restore();
        REQUIRE(p.get_value() == -0.3);
        REQUIRE(p.get_stored_value() == -0.3);

        REQUIRE_THROWS_AS(p.set_value(0.01), EcoevolityParameterValueError);
    }
}

TEST_CASE("Testing PositiveRealParameter value methods", "[PositiveRealParameter]") {

    SECTION("Testing value methods") {
        REQUIRE_THROWS_AS(PositiveRealParameter(-0.1), EcoevolityParameterValueError);
        PositiveRealParameter p = PositiveRealParameter(0.1);
        REQUIRE(p.get_value() == 0.1);
        p.store();
        p.set_value(0.2);
        REQUIRE(p.get_value() == 0.2);
        REQUIRE(p.get_stored_value() == 0.1);

        p.restore();
        REQUIRE(p.get_value() == 0.1);
        REQUIRE(p.get_stored_value() == 0.1);

        p.update_value(0.3);
        REQUIRE(p.get_value() == 0.3);
        REQUIRE(p.get_stored_value() == 0.1);

        p.update_value(0.4);
        REQUIRE(p.get_value() == 0.4);
        REQUIRE(p.get_stored_value() == 0.3);

        p.restore();
        REQUIRE(p.get_value() == 0.3);
        REQUIRE(p.get_stored_value() == 0.3);

        REQUIRE_THROWS_AS(p.set_value(-0.01), EcoevolityParameterValueError);
    }
}
