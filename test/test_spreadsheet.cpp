#include "catch.hpp"
#include "ecoevolity/spreadsheet.hpp"

TEST_CASE("Testing parse_header", "[spreadsheet]") {

    SECTION("Testing simple header with three columns") {
        std::stringstream stream;
        stream << "col1\tcol2\tcol3\n";
        stream << "1.0\t2.0\t3.0\n";
        std::vector<std::string> expected = {"col1", "col2", "col3"};

        std::vector<std::string> header;
        header = spreadsheet::parse_header(stream);
        REQUIRE(header == expected);
    }
}
