#include "catch.hpp"
#include "ecoevolity/rng.hpp"

TEST_CASE("Testing seed constructor", "[RandomNumberGenerator]") {

    SECTION("Testing seed constructor") {
        RandomNumberGenerator rng = RandomNumberGenerator(123);
        REQUIRE(rng.get_seed() == 123);
    }
}

TEST_CASE("Testing bare constructor", "[RandomNumberGenerator]") {

    SECTION("Testing seed constructor") {
        RandomNumberGenerator rng = RandomNumberGenerator();
        REQUIRE(rng.get_seed() > 0);
    }
}

TEST_CASE("Testing uniform_real", "[RandomNumberGenerator]") {

    SECTION("Testing uniform_real") {
        RandomNumberGenerator rng(123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = rng.uniform_real();
            REQUIRE(x >= 0.0);
            REQUIRE(x < 1.0);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(0.5).epsilon(0.001));
        REQUIRE(variance == Approx(1.0/12.0).epsilon(0.001));
    }
}
