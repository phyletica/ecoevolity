#include "catch.hpp"
#include "ecoevolity/probability.hpp"

#include <limits>

TEST_CASE("Testing ImproperUniformDistribution", "[ImproperUniformDistribution]") {

    SECTION("Testing ImproperUniformDistribution") {
        ImproperUniformDistribution u = ImproperUniformDistribution();
        REQUIRE(u.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(u.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(u.to_string() == "uniform(-inf, +inf)");
        REQUIRE_THROWS_AS(u.get_mean(), EcoevolityProbabilityDistributionError);
        REQUIRE_THROWS_AS(u.get_variance(), EcoevolityProbabilityDistributionError);
        REQUIRE_THROWS_AS(u.ln_pdf(1.0), EcoevolityProbabilityDistributionError);
        REQUIRE(u.relative_ln_pdf(1.0) == 0.0);
        REQUIRE(u.relative_ln_pdf(-1.0) == 0.0);
        REQUIRE(u.relative_ln_pdf(0.0) == 0.0);
    }
}

TEST_CASE("Testing UniformDistribution", "[UniformDistribution]") {

    SECTION("Testing UniformDistribution bare constructor") {
        UniformDistribution u = UniformDistribution();
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

    SECTION("Testing UniformDistribution bare constructor") {
        UniformDistribution u = UniformDistribution(10.0, 20.0);
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

TEST_CASE("Testing OffsetGammaDistribution", "[OffsetGammaDistribution]") {
    SECTION("Testing bare constructor") {
        OffsetGammaDistribution f = OffsetGammaDistribution();
        REQUIRE(f.to_string() == "gamma(shape = 1, scale = 1, offset = 0)");
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
    }
}
