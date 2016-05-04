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

    SECTION("splitting string by tabs") {
        string_util::split("Test\tstring\t\tTO\tSplit", '\t', elements);
        REQUIRE(elements.size() == 5);
        REQUIRE(elements[0] == "Test");
        REQUIRE(elements[1] == "string");
        REQUIRE(elements[2] == "");
        REQUIRE(elements[3] == "TO");
        REQUIRE(elements[4] == "Split");
    }
}

TEST_CASE("Testing get_indent", "[string_util]") {
    SECTION("Testing get_indent") {
        REQUIRE(string_util::get_indent() == "    ");
        REQUIRE(string_util::get_indent(0) == "");
        REQUIRE(string_util::get_indent(1) == "    ");
        REQUIRE(string_util::get_indent(2) == "        ");
        REQUIRE(string_util::get_indent(3) == "            ");
    }
}

TEST_CASE("Testing lstrip", "[string_util]") {
    SECTION("Testing all white space") {
        std::string s = "\n\t   \t\n";
        std::string t = string_util::lstrip(s);
        REQUIRE(s == "\n\t   \t\n");
        REQUIRE(t == "");
    }

    SECTION("Testing internal white space") {
        std::string s = "hello \t\n world";
        std::string t = string_util::lstrip(s);
        REQUIRE(s == t);
    }
    
    SECTION("Testing left white space") {
        std::string s = " \t \nhello";
        std::string t = string_util::lstrip(s);
        REQUIRE(s == " \t \nhello");
        REQUIRE(t == "hello");
    }

    SECTION("Testing right white space") {
        std::string s = "hello \t \n";
        std::string t = string_util::lstrip(s);
        REQUIRE(s == t);
    }

    SECTION("Testing empty string") {
        std::string s = "";
        std::string t = string_util::lstrip(s);
        REQUIRE(s == t);
    }
}

TEST_CASE("Testing rstrip", "[string_util]") {
    SECTION("Testing all white space") {
        std::string s = "\n\t   \t\n";
        std::string t = string_util::rstrip(s);
        REQUIRE(s == "\n\t   \t\n");
        REQUIRE(t == "");
    }

    SECTION("Testing internal white space") {
        std::string s = "hello \t\n world";
        std::string t = string_util::rstrip(s);
        REQUIRE(s == t);
    }
    
    SECTION("Testing left white space") {
        std::string s = " \t \nhello";
        std::string t = string_util::rstrip(s);
        REQUIRE(s == t);
    }

    SECTION("Testing right white space") {
        std::string s = "hello \t \n";
        std::string t = string_util::rstrip(s);
        REQUIRE(s == "hello \t \n");
        REQUIRE(t == "hello");
    }

    SECTION("Testing empty string") {
        std::string s = "";
        std::string t = string_util::rstrip(s);
        REQUIRE(s == t);
    }
}

TEST_CASE("Testing strip", "[string_util]") {
    SECTION("Testing all white space") {
        std::string s = "\n\t   \t\n";
        std::string t = string_util::strip(s);
        REQUIRE(s == "\n\t   \t\n");
        REQUIRE(t == "");
    }

    SECTION("Testing internal white space") {
        std::string s = "hello \t\n world";
        std::string t = string_util::strip(s);
        REQUIRE(s == t);
    }

    SECTION("Testing left and right white space") {
        std::string s = "\n  \t \t hello \t \n";
        std::string t = string_util::strip(s);
        REQUIRE(s == "\n  \t \t hello \t \n");
        REQUIRE(t == "hello");
    }

    SECTION("Testing empty string") {
        std::string s = "";
        std::string t = string_util::strip(s);
        REQUIRE(s == t);
    }
}

TEST_CASE("Testing join", "[string_util]") {
    std::vector<std::string> strings {"1", "2", "3"};
    std::string s = string_util::join(strings, "---");
    REQUIRE(s == "1---2---3");
}

TEST_CASE("Testing join on empty vector", "[string_util]") {
    std::vector<std::string> strings;
    std::string s = string_util::join(strings, "---");
    REQUIRE(s == "");
}

TEST_CASE("Testing join on singleton vector", "[string_util]") {
    std::vector<std::string> strings {"1"};
    std::string s = string_util::join(strings, "---");
    REQUIRE(s == "1");
}
