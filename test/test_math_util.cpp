#include "catch.hpp"
#include "ecoevolity/math_util.hpp"

TEST_CASE("Testing normalize_log_likelihoods", "[math_util]") {

    SECTION("Testing round trip") {
        std::vector<double> v = {0.2, 0.2, 0.2, 0.2, 0.2};
        std::vector<double> v2(v.size());
        for (unsigned int i = 0; i < v.size(); ++i) {
            v2.at(i) = std::log(v.at(i));
        }
        normalize_log_likelihoods(v2);
        for (unsigned int i = 0; i < v.size(); ++i) {
            REQUIRE(v.at(i) == Approx(v2.at(i)));
        }
    }

    SECTION("Testing uneven round trip") {
        std::vector<double> v = {0.1, 0.2, 0.3, 0.4};
        std::vector<double> v2(v.size());
        for (unsigned int i = 0; i < v.size(); ++i) {
            v2.at(i) = std::log(v.at(i));
        }
        normalize_log_likelihoods(v2);
        for (unsigned int i = 0; i < v.size(); ++i) {
            REQUIRE(v.at(i) == Approx(v2.at(i)));
        }
    }

    SECTION("Testing uneven unnormalized round trip") {
        std::vector<double> v = {0.00001, 0.00002, 0.00003, 0.00004};
        std::vector<double> v2(v.size());
        for (unsigned int i = 0; i < v.size(); ++i) {
            v2.at(i) = std::log(v.at(i));
        }
        normalize_log_likelihoods(v2);
        for (unsigned int i = 0; i < v.size(); ++i) {
            REQUIRE(v.at(i) * 10000.0 == Approx(v2.at(i)));
        }
    }
}

TEST_CASE("Testing normalize_log_likelihoods of copy", "[math_util]") {
    SECTION("Testing copy normalization") {
        std::vector<double> v = {-23.45, -12.24, -13.46, -45.12};
        std::vector<double> v2(v);
        normalize_log_likelihoods(v2);
        REQUIRE(v.at(0) == -23.45);
        REQUIRE(v.at(1) == -12.24);
        REQUIRE(v.at(2) == -13.46);
        REQUIRE(v.at(3) == -45.12);
    }
}
