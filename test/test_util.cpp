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
