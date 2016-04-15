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

TEST_CASE("Testing get_dpp_expected_number_of_categories 0.218 7282", "[math_util]") {
    double e = get_dpp_expected_number_of_categories(0.218, 7282);
    REQUIRE(e == Approx(3.0).epsilon(0.001));
}

TEST_CASE("Testing get_dpp_expected_number_of_categories 0.449 7282", "[math_util]") {
    double e = get_dpp_expected_number_of_categories(0.449, 7282);
    REQUIRE(e == Approx(5.0).epsilon(0.001));
}

TEST_CASE("Testing get_dpp_expected_number_of_categories 0.814 7282", "[math_util]") {
    double e = get_dpp_expected_number_of_categories(0.814, 7282);
    REQUIRE(e == Approx(8.0).epsilon(0.001));
}

TEST_CASE("Testing get_dpp_expected_number_of_categories 1.068 7282", "[math_util]") {
    double e = get_dpp_expected_number_of_categories(1.068, 7282);
    REQUIRE(e == Approx(10.0).epsilon(0.0001));
}

TEST_CASE("Testing get_dpp_concentration 3.0 7282", "[math_util]") {
    double e = get_dpp_concentration(3.0, 7282);
    REQUIRE(e == Approx(0.218).epsilon(0.001));
}
TEST_CASE("Testing get_dpp_concentration 5.0 7282", "[math_util]") {
    double e = get_dpp_concentration(5.0, 7282);
    REQUIRE(e == Approx(0.449).epsilon(0.001));
}
TEST_CASE("Testing get_dpp_concentration 8.0 7282", "[math_util]") {
    double e = get_dpp_concentration(8.0, 7282);
    REQUIRE(e == Approx(0.814).epsilon(0.0001));
}
TEST_CASE("Testing get_dpp_concentration 10.0 7282", "[math_util]") {
    double e = get_dpp_concentration(10.0, 7282);
    REQUIRE(e == Approx(1.068).epsilon(0.0001));
}

TEST_CASE("Testing get_dpp_gamma_scale", "[math_util]") {
    unsigned int n = 7282;
    double ncats = 10.0;
    double concentration = 1.068;
    std::vector<double> shapes {0.001, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0, 10000.0};
    for (auto shape : shapes) {
        double scale = get_dpp_gamma_scale(ncats, n, shape);
        REQUIRE(concentration == Approx(shape * scale).epsilon(0.001));
    }
}
