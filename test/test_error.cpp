#include "catch.hpp"
#include "coevolity/error.hpp"

TEST_CASE("derived error classes can be thrown", "[error]") {

    std::string message = "Testing Error";
    std::string file_name = "dummy-file-name.txt";
    size_t line_number = 1;

    SECTION("throwing CoevolityError") {
        REQUIRE_THROWS_AS(throw CoevolityError(message), CoevolityError);
    }
    SECTION("throwing CoevolityParsingError") {
        REQUIRE_THROWS_AS(throw CoevolityParsingError(message, file_name, line_number),
                CoevolityParsingError);
    }
}
