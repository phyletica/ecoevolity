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

SCENARIO("sort_pairs provides easy sorting", "[util]") {

    GIVEN("An unsorted vector of pairs of ints and unsigned ints") {
        std::vector< std::pair<int, unsigned int> > pairs;
        pairs.push_back(std::make_pair(1, 3));
        pairs.push_back(std::make_pair(-1, 5));
        pairs.push_back(std::make_pair(2, 2));
        pairs.push_back(std::make_pair(-3, 9));
        pairs.push_back(std::make_pair(0, 0));

        WHEN("sort_pairs(pairs, true, false) is called") {
            sort_pairs(pairs, true, false);
            
            THEN("the pairs are sorted ascending by first") {
                std::vector< std::pair<int, unsigned int> > expected;
                expected.push_back(std::make_pair(-3, 9));
                expected.push_back(std::make_pair(-1, 5));
                expected.push_back(std::make_pair(0, 0));
                expected.push_back(std::make_pair(1, 3));
                expected.push_back(std::make_pair(2, 2));

                REQUIRE(pairs == expected);
            }
        }

        WHEN("sort_pairs(pairs, true, true) is called") {
            sort_pairs(pairs, true, true);
            
            THEN("the pairs are sorted descending by first") {
                std::vector< std::pair<int, unsigned int> > expected;
                expected.push_back(std::make_pair(2, 2));
                expected.push_back(std::make_pair(1, 3));
                expected.push_back(std::make_pair(0, 0));
                expected.push_back(std::make_pair(-1, 5));
                expected.push_back(std::make_pair(-3, 9));

                REQUIRE(pairs == expected);
            }
        }

        WHEN("sort_pairs(pairs, false, false) is called") {
            sort_pairs(pairs, false, false);
            
            THEN("the pairs are sorted ascending by second") {
                std::vector< std::pair<int, unsigned int> > expected;
                expected.push_back(std::make_pair(0, 0));
                expected.push_back(std::make_pair(2, 2));
                expected.push_back(std::make_pair(1, 3));
                expected.push_back(std::make_pair(-1, 5));
                expected.push_back(std::make_pair(-3, 9));

                REQUIRE(pairs == expected);
            }
        }

        WHEN("sort_pairs(pairs, false, true) is called") {
            sort_pairs(pairs, false, true);
            
            THEN("the pairs are sorted descending by second") {
                std::vector< std::pair<int, unsigned int> > expected;
                expected.push_back(std::make_pair(-3, 9));
                expected.push_back(std::make_pair(-1, 5));
                expected.push_back(std::make_pair(1, 3));
                expected.push_back(std::make_pair(2, 2));
                expected.push_back(std::make_pair(0, 0));

                REQUIRE(pairs == expected);
            }
        }
    }
}
