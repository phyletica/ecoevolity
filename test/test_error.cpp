#include "catch.hpp"
#include "ecoevolity/error.hpp"

TEST_CASE("derived error classes can be thrown", "[error]") {

    std::string message = "Testing Error";
    std::string file_name = "dummy-file-name.txt";
    size_t line_number = 1;

    SECTION("throwing EcoevolityError") {
        REQUIRE_THROWS_AS(throw EcoevolityError(message), EcoevolityError);
    }
    SECTION("throwing EcoevolityBiallelicDataError") {
        REQUIRE_THROWS_AS(throw EcoevolityBiallelicDataError(message, file_name),
                EcoevolityBiallelicDataError);
    }
    SECTION("throwing EcoevolityParsingError") {
        REQUIRE_THROWS_AS(throw EcoevolityParsingError(message, file_name, line_number),
                EcoevolityParsingError);
    }
    SECTION("throwing EcoevolityInvalidCharacterError") {
        REQUIRE_THROWS_AS(throw EcoevolityInvalidCharacterError(
                    message,
                    file_name,
                    "test-taxon-name",
                    1),
                EcoevolityInvalidCharacterError);
    }
}
