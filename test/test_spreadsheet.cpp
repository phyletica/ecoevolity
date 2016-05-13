#include "catch.hpp"
#include "ecoevolity/spreadsheet.hpp"

TEST_CASE("Testing simple header parsing",
        "[spreadsheet]") {

    SECTION("Testing bare") {
        std::stringstream stream;
        stream << "col1\tcol2\tcol3\n";

        std::vector<std::string> header;
        header = spreadsheet::parse_header(stream);
    }
}
