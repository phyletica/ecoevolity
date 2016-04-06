#include "catch.hpp"
#include "ecoevolity/util.hpp"

#include <map>
#include <unordered_map>

SCENARIO("map_at provides easy map element access", "[util]") {

    GIVEN("A (ordered) map with some key-value pairs") {
        std::map<std::string, int> test_map;
        test_map["a"] = 1;
        test_map["b"] = 2;
        test_map["c"] = 3;

        REQUIRE(test_map.size() == 3);

        WHEN("map_at is called with an existing key") {
            int value = map_at(test_map, "b");
            
            THEN("the corresponding element is returned") {
                REQUIRE(value == 2);
            }
        }

        WHEN("map_at is called with a key that is not in the map") {

            THEN("an out_of_range error is thrown") {

                REQUIRE_THROWS_AS(map_at(test_map, "d"), std::out_of_range);
            }
        }
    }

    GIVEN("An unordered map with some key-value pairs") {
        std::unordered_map<std::string, int> test_map;
        test_map["a"] = 1;
        test_map["b"] = 2;
        test_map["c"] = 3;

        REQUIRE(test_map.size() == 3);

        WHEN("map_at is called with an existing key") {
            int value = map_at(test_map, "b");
            
            THEN("the corresponding element is returned") {
                REQUIRE(value == 2);
            }
        }

        WHEN("map_at is called with a key that is not in the map") {

            THEN("an out_of_range error is thrown") {

                REQUIRE_THROWS_AS(map_at(test_map, "d"), std::out_of_range);
            }
        }
    }
}
            
SCENARIO("split provides Python-like splitting of strings", "[util]") {
    std::vector<std::string> elements;

    REQUIRE(elements.size() == 0);

    SECTION("splitting string into pre-allocated vector") {
        split("Test_string__TO_Split", '_', elements);
        REQUIRE(elements.size() == 5);
        REQUIRE(elements[0] == "Test");
        REQUIRE(elements[1] == "string");
        REQUIRE(elements[2] == "");
        REQUIRE(elements[3] == "TO");
        REQUIRE(elements[4] == "Split");
    }

    SECTION("splitting string and returning elements in a new vector") {
        std::vector<std::string> words = split("Test_string__TO_Split", '_');
        REQUIRE(words.size() == 5);
        REQUIRE(words[0] == "Test");
        REQUIRE(words[1] == "string");
        REQUIRE(words[2] == "");
        REQUIRE(words[3] == "TO");
        REQUIRE(words[4] == "Split");
    }
}

TEST_CASE("Testing normalize_log_likelihoods", "[util]") {

    SECTION("Testing round trip", "[util]") {
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

TEST_CASE("Testing normalize_log_likelihoods of copy", "[util]") {
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
