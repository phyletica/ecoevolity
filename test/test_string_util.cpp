#include "catch.hpp"
#include "ecoevolity/string_util.hpp"

SCENARIO("split provides Python-like splitting of strings", "[string_util]") {
    std::vector<std::string> elements;

    REQUIRE(elements.size() == 0);

    SECTION("splitting string into pre-allocated vector") {
        string_util::split("Test_string__TO_Split", '_', elements);
        REQUIRE(elements.size() == 5);
        REQUIRE(elements[0] == "Test");
        REQUIRE(elements[1] == "string");
        REQUIRE(elements[2] == "");
        REQUIRE(elements[3] == "TO");
        REQUIRE(elements[4] == "Split");
    }

    SECTION("splitting string and returning elements in a new vector") {
        std::vector<std::string> words = string_util::split("Test_string__TO_Split", '_');
        REQUIRE(words.size() == 5);
        REQUIRE(words[0] == "Test");
        REQUIRE(words[1] == "string");
        REQUIRE(words[2] == "");
        REQUIRE(words[3] == "TO");
        REQUIRE(words[4] == "Split");
    }
}

TEST_CASE("Testing get_indent", "[util]") {
    SECTION("Testing get_indent") {
        REQUIRE(string_util::get_indent() == "    ");
        REQUIRE(string_util::get_indent(0) == "");
        REQUIRE(string_util::get_indent(1) == "    ");
        REQUIRE(string_util::get_indent(2) == "        ");
        REQUIRE(string_util::get_indent(3) == "            ");
    }
}
