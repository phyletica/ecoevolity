#include "catch.hpp"
#include "ecoevolity/probability.hpp"

#include <limits>

TEST_CASE("Testing ImproperUniformDist", "[ImproperUniformDist]") {

    SECTION("Testing ImproperUniformDist") {
        ImproperUniformDist u = ImproperUniformDist();
        REQUIRE(u.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(u.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(u.to_string() == "uniform(-inf, +inf)");
        REQUIRE_THROWS_AS(u.get_mean(), EcoevolityProbabilityDistError);
        REQUIRE_THROWS_AS(u.get_variance(), EcoevolityProbabilityDistError);
        REQUIRE_THROWS_AS(u.ln_pdf(1.0), EcoevolityProbabilityDistError);
        REQUIRE(u.relative_ln_pdf(1.0) == 0.0);
        REQUIRE(u.relative_ln_pdf(100.0) == 0.0);
        REQUIRE(u.relative_ln_pdf(-1.0) == 0.0);
        REQUIRE(u.relative_ln_pdf(0.0) == 0.0);
    }
}

TEST_CASE("Testing ImproperPositiveUniformDist", "[ImproperPositiveUniformDist]") {

    SECTION("Testing ImproperPositiveUniformDist") {
        ImproperPositiveUniformDist u = ImproperPositiveUniformDist();
        REQUIRE(u.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(u.get_min() == 0.0);
        REQUIRE(u.to_string() == "uniform(0, +inf)");
        REQUIRE_THROWS_AS(u.get_mean(), EcoevolityProbabilityDistError);
        REQUIRE_THROWS_AS(u.get_variance(), EcoevolityProbabilityDistError);
        REQUIRE_THROWS_AS(u.ln_pdf(1.0), EcoevolityProbabilityDistError);
        REQUIRE(u.relative_ln_pdf(1.0) == 0.0);
        REQUIRE(u.relative_ln_pdf(100.0) == 0.0);
        REQUIRE(u.relative_ln_pdf(-1.0) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.relative_ln_pdf(0.0) == 0.0);
    }
}

TEST_CASE("Testing UniformDist", "[UniformDist]") {

    SECTION("Testing UniformDist bare constructor") {
        UniformDist u = UniformDist();
        REQUIRE(u.get_max() == 1.0);
        REQUIRE(u.get_min() == 0.0);
        REQUIRE(u.to_string() == "uniform(0, 1)");
        REQUIRE(u.get_mean() == 0.5);
        REQUIRE(u.get_variance() == Approx(1.0/12.0));
        REQUIRE(u.relative_ln_pdf(1.0) == 0.0);
        REQUIRE(u.relative_ln_pdf(0.0) == 0.0);
        REQUIRE(u.ln_pdf(1.0) == 0.0);
        REQUIRE(u.ln_pdf(0.0) == 0.0);
        REQUIRE(u.relative_ln_pdf(-1.0) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.relative_ln_pdf(1.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.ln_pdf(-1.0) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.ln_pdf(1.1) == -std::numeric_limits<double>::infinity());
    }

    SECTION("Testing UniformDist bare constructor") {
        UniformDist u = UniformDist(10.0, 20.0);
        REQUIRE(u.get_max() == 20.0);
        REQUIRE(u.get_min() == 10.0);
        REQUIRE(u.to_string() == "uniform(10, 20)");
        REQUIRE(u.get_mean() == 15.0);
        REQUIRE(u.get_variance() == Approx(25.0/3.0));
        REQUIRE(u.relative_ln_pdf(10.0) == Approx(std::log(0.1)));
        REQUIRE(u.relative_ln_pdf(20.0) == Approx(std::log(0.1)));
        REQUIRE(u.relative_ln_pdf(11.0) == Approx(std::log(0.1)));
        REQUIRE(u.ln_pdf(10.0) == Approx(std::log(0.1)));
        REQUIRE(u.ln_pdf(20.0) == Approx(std::log(0.1)));
        REQUIRE(u.ln_pdf(11.0) == Approx(std::log(0.1)));
        REQUIRE(u.relative_ln_pdf(9.9) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.relative_ln_pdf(20.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.relative_ln_pdf(-15.0) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.ln_pdf(9.9) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.ln_pdf(20.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.ln_pdf(-15.0) == -std::numeric_limits<double>::infinity());
    }
}

TEST_CASE("Testing OffsetGammaDist", "[OffsetGammaDist]") {
    SECTION("Testing bare constructor") {
        OffsetGammaDist f = OffsetGammaDist();
        REQUIRE(f.to_string() == "gamma(shape = 1, scale = 1, offset = 0)");
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_mean() == 1.0);
        REQUIRE(f.get_variance() == 1.0);

        REQUIRE(f.get_offset() == 0.0);
        REQUIRE(f.get_shape() == 1.0);
        REQUIRE(f.get_scale() == 1.0);

        REQUIRE(f.ln_pdf(-0.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(0.0) == 0.0);
        REQUIRE(f.relative_ln_pdf(0.0) == 0.0);
        REQUIRE(f.ln_pdf(0.01) == -0.01);
        REQUIRE(f.relative_ln_pdf(0.01) == -0.01);
        REQUIRE(f.ln_pdf(1.0) == -1.0);
        REQUIRE(f.relative_ln_pdf(1.0) == -1.0);
        REQUIRE(f.ln_pdf(100.0) == -100.0);
        REQUIRE(f.relative_ln_pdf(100.0) == -100.0);
    }

    SECTION("Testing constructor errors") {
        REQUIRE_THROWS_AS(OffsetGammaDist(0.0, 1.0, 1.0), EcoevolityProbabilityDistError);
        REQUIRE_THROWS_AS(OffsetGammaDist(1.0, 0.0, 1.0), EcoevolityProbabilityDistError);
        REQUIRE_THROWS_AS(OffsetGammaDist(-1.0, 1.0, 1.0), EcoevolityProbabilityDistError);
        REQUIRE_THROWS_AS(OffsetGammaDist(1.0, -1.0, 1.0), EcoevolityProbabilityDistError);
    }

    SECTION("Testing OffsetGammaDist(2, 4, 5)") {
        OffsetGammaDist f = OffsetGammaDist(2.0, 4.0, 5.0);
        REQUIRE(f.to_string() == "gamma(shape = 2, scale = 4, offset = 5)");
        REQUIRE(f.get_min() == 5.0);
        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_mean() == 13.0);
        REQUIRE(f.get_variance() == 32.0);
        REQUIRE(f.ln_pdf(5.0) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(4.9) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(6.0) == Approx(-3.0225887222397811));
        REQUIRE(f.ln_pdf(13.0) == Approx(-2.6931471805599454));
        REQUIRE(f.ln_pdf(105.0) == Approx(-23.167418536251692));

        REQUIRE(f.get_offset() == 5.0);
        REQUIRE(f.get_shape() == 2.0);
        REQUIRE(f.get_scale() == 4.0);
    }

    SECTION("Testing OffsetGammaDist(0.1, 100, -5)") {
        OffsetGammaDist f = OffsetGammaDist(0.1, 100.0, -5.0);
        REQUIRE(f.to_string() == "gamma(shape = 0.1, scale = 100, offset = -5)");
        REQUIRE(f.get_min() == -5.0);
        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_mean() == 5.0);
        REQUIRE(f.get_variance() == 1000.0);
        REQUIRE(f.ln_pdf(-5.0) == std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(-5.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(-4.0) == Approx(-2.7232296703330157));
        REQUIRE(f.ln_pdf(0.0) == Approx(-4.2117237915237062));
        REQUIRE(f.ln_pdf(5.0) == Approx(-4.8855562540276569));
        REQUIRE(f.ln_pdf(45.0) == Approx(-6.7340503752183469));

        REQUIRE(f.get_offset() == -5.0);
        REQUIRE(f.get_shape() == 0.1);
        REQUIRE(f.get_scale() == 100.0);
    }
}

TEST_CASE("Testing GammaDist", "[GammaDist]") {
    SECTION("Testing bare constructor") {
        GammaDist f = GammaDist();
        REQUIRE(f.to_string() == "gamma(shape = 1, scale = 1)");
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_mean() == 1.0);
        REQUIRE(f.get_variance() == 1.0);

        REQUIRE(f.ln_pdf(-0.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(0.0) == 0.0);
        REQUIRE(f.relative_ln_pdf(0.0) == 0.0);
        REQUIRE(f.ln_pdf(0.01) == -0.01);
        REQUIRE(f.relative_ln_pdf(0.01) == -0.01);
        REQUIRE(f.ln_pdf(1.0) == -1.0);
        REQUIRE(f.relative_ln_pdf(1.0) == -1.0);
        REQUIRE(f.ln_pdf(100.0) == -100.0);
        REQUIRE(f.relative_ln_pdf(100.0) == -100.0);

        REQUIRE(f.get_offset() == 0.0);
        REQUIRE(f.get_shape() == 1.0);
        REQUIRE(f.get_scale() == 1.0);
    }
    SECTION("Testing constructor errors") {
        REQUIRE_THROWS_AS(GammaDist(0.0, 1.0), EcoevolityProbabilityDistError);
        REQUIRE_THROWS_AS(GammaDist(1.0, 0.0), EcoevolityProbabilityDistError);
        REQUIRE_THROWS_AS(GammaDist(-1.0, 1.0), EcoevolityProbabilityDistError);
        REQUIRE_THROWS_AS(GammaDist(1.0, -1.0), EcoevolityProbabilityDistError);
    }

    SECTION("Testing GammaDist(2, 4)") {
        GammaDist f = GammaDist(2.0, 4.0);
        REQUIRE(f.to_string() == "gamma(shape = 2, scale = 4)");
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_mean() == 8.0);
        REQUIRE(f.get_variance() == 32.0);
        REQUIRE(f.ln_pdf(0.0) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(1.0) == Approx(-3.0225887222397811));
        REQUIRE(f.ln_pdf(8.0) == Approx(-2.6931471805599454));
        REQUIRE(f.ln_pdf(100.0) == Approx(-23.167418536251692));

        REQUIRE(f.get_offset() == 0.0);
        REQUIRE(f.get_shape() == 2.0);
        REQUIRE(f.get_scale() == 4.0);
    }

    SECTION("Testing GammaDist(0.1, 100)") {
        GammaDist f = GammaDist(0.1, 100.0);
        REQUIRE(f.to_string() == "gamma(shape = 0.1, scale = 100)");
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_mean() == 10.0);
        REQUIRE(f.get_variance() == 1000.0);
        REQUIRE(f.ln_pdf(0.0) == std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(1.0) == Approx(-2.7232296703330157));
        REQUIRE(f.ln_pdf(5.0) == Approx(-4.2117237915237062));
        REQUIRE(f.ln_pdf(10.0) == Approx(-4.8855562540276569));
        REQUIRE(f.ln_pdf(50.0) == Approx(-6.7340503752183469));

        REQUIRE(f.get_offset() == 0.0);
        REQUIRE(f.get_shape() == 0.1);
        REQUIRE(f.get_scale() == 100.0);
    }
}

TEST_CASE("Testing OffsetExponentialDist", "[OffsetExponentialDist]") {
    SECTION("Testing bare constructor") {
        OffsetExponentialDist f = OffsetExponentialDist();

        REQUIRE(f.to_string() == "exp(lambda = 1, offset = 0)");

        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_offset() == 0.0);
        REQUIRE(f.get_shape() == 1.0);
        REQUIRE(f.get_scale() == 1.0);
        REQUIRE(f.get_lambda() == 1.0);

        REQUIRE(f.get_mean() == 1.0);
        REQUIRE(f.get_variance() == 1.0);

        REQUIRE(f.ln_pdf(-0.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(0.0) == 0.0);
        REQUIRE(f.ln_pdf(0.01) == -0.01);
        REQUIRE(f.ln_pdf(1.0) == -1.0);
        REQUIRE(f.ln_pdf(100.0) == -100.0);
    }

    SECTION("Testing constructor errors") {
        REQUIRE_THROWS_AS(OffsetExponentialDist(0.0, 1.0), EcoevolityProbabilityDistError);
        REQUIRE_THROWS_AS(OffsetExponentialDist(-0.01, 1.0), EcoevolityProbabilityDistError);
    }

    SECTION("Testing OffsetExponentialDist(5, -5)") {
        OffsetExponentialDist f = OffsetExponentialDist(5.0, -5.0);

        REQUIRE(f.to_string() == "exp(lambda = 5, offset = -5)");

        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_min() == -5.0);
        REQUIRE(f.get_offset() == -5.0);
        REQUIRE(f.get_shape() == 1.0);
        REQUIRE(f.get_scale() == Approx(1.0/5.0));
        REQUIRE(f.get_lambda() == 5.0);

        REQUIRE(f.get_mean() == Approx(-4.8));
        REQUIRE(f.get_variance() == Approx(1.0/25.0));

        REQUIRE(f.ln_pdf(-5.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(-5.0) == Approx(1.6094379124341003));
        REQUIRE(f.ln_pdf(-4.99) == Approx(1.5594379124341002));
        REQUIRE(f.ln_pdf(-4.0) == Approx(-3.3905620875658995));
        REQUIRE(f.ln_pdf(95.0) == Approx(-498.39056208756591));
    }
}

TEST_CASE("Testing ExponentialDist", "[ExponentialDist]") {
    SECTION("Testing bare constructor") {
        ExponentialDist f = ExponentialDist();

        REQUIRE(f.to_string() == "exp(lambda = 1)");

        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_offset() == 0.0);
        REQUIRE(f.get_shape() == 1.0);
        REQUIRE(f.get_scale() == 1.0);
        REQUIRE(f.get_lambda() == 1.0);

        REQUIRE(f.get_mean() == 1.0);
        REQUIRE(f.get_variance() == 1.0);

        REQUIRE(f.ln_pdf(-0.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(0.0) == 0.0);
        REQUIRE(f.ln_pdf(0.01) == -0.01);
        REQUIRE(f.ln_pdf(1.0) == -1.0);
        REQUIRE(f.ln_pdf(100.0) == -100.0);
    }

    SECTION("Testing constructor errors") {
        REQUIRE_THROWS_AS(ExponentialDist(0.0), EcoevolityProbabilityDistError);
        REQUIRE_THROWS_AS(ExponentialDist(-0.01), EcoevolityProbabilityDistError);
    }

    SECTION("Testing ExponentialDist(5)") {
        ExponentialDist f = ExponentialDist(5.0);

        REQUIRE(f.to_string() == "exp(lambda = 5)");

        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_offset() == 0.0);
        REQUIRE(f.get_shape() == 1.0);
        REQUIRE(f.get_scale() == Approx(1.0/5.0));
        REQUIRE(f.get_lambda() == 5.0);

        REQUIRE(f.get_mean() == Approx(0.2));
        REQUIRE(f.get_variance() == Approx(1.0/25.0));

        REQUIRE(f.ln_pdf(-0.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(0.0) == Approx(1.6094379124341003));
        REQUIRE(f.ln_pdf(0.01) == Approx(1.5594379124341002));
        REQUIRE(f.ln_pdf(1.0) == Approx(-3.3905620875658995));
        REQUIRE(f.ln_pdf(100.0) == Approx(-498.39056208756591));
    }
}
