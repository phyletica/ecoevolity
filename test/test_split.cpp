#include "catch.hpp"
#include "ecoevolity/split.hpp"

TEST_CASE("Testing 2 leaves", "[split]") {
    SECTION("Testing 2 leaves") {
        Split s10;
        s10.resize(2);
        s10.set_leaf_bit(0);
        REQUIRE(s10.get_leaf_bit(0) == 1);
        REQUIRE(s10.get_leaf_bit(1) == 0);
        REQUIRE(s10.as_string() == "10");

        Split s01;
        s01.resize(2);
        s01.set_leaf_bit(1);
        REQUIRE(s01.get_leaf_bit(0) == 0);
        REQUIRE(s01.get_leaf_bit(1) == 1);
        REQUIRE(s01.as_string() == "01");

        Split s00;
        s00.resize(2);
        REQUIRE(s00.get_leaf_bit(0) == 0);
        REQUIRE(s00.get_leaf_bit(1) == 0);
        REQUIRE(s00.as_string() == "00");

        Split s11;
        s11.resize(2);
        s11.set_leaf_bit(0);
        s11.set_leaf_bit(1);
        REQUIRE(s11.get_leaf_bit(0) == 1);
        REQUIRE(s11.get_leaf_bit(1) == 1);
        REQUIRE(s11.as_string() == "11");

        REQUIRE(s10 != s01);
        REQUIRE(! (s10 == s01));
        REQUIRE(! s10.is_equivalent(s01, true));
        REQUIRE(s10.is_equivalent(s01, false));
        REQUIRE(s10.is_compatible(s01));
        REQUIRE(! s10.conflicts_with(s01));

        REQUIRE(s00 != s11);
        REQUIRE(! (s00 == s11));
        REQUIRE(! s00.is_equivalent(s11, true));
        REQUIRE(s00.is_equivalent(s11, false));
        REQUIRE(s00.is_compatible(s11));
        REQUIRE(! s00.conflicts_with(s11));

        REQUIRE(s01 != s11);
        REQUIRE(! (s01 == s11));
        REQUIRE(! s01.is_equivalent(s11, true));
        REQUIRE(! s01.is_equivalent(s11, false));
        REQUIRE(s01.is_compatible(s11));
        REQUIRE(! s01.conflicts_with(s11));

        Split ss10(s10);
        REQUIRE(ss10.as_string() == "10");
        REQUIRE(ss10.get_leaf_bit(0) == 1);
        REQUIRE(ss10.get_leaf_bit(1) == 0);
        REQUIRE(ss10 == s10);
        REQUIRE(! (ss10 != s10));
        REQUIRE(ss10.is_equivalent(s10, true));
        REQUIRE(ss10.is_equivalent(s10, false));
        REQUIRE(ss10.is_compatible(s10));
        REQUIRE(! ss10.conflicts_with(s10));

        REQUIRE(ss10.as_string('-', '*') == "*-");
    }
}

TEST_CASE("Testing 5 leaves", "[split]") {
    SECTION("Testing 5 leaves") {
        Split s10100;
        s10100.resize(5);
        s10100.set_leaf_bit(0);
        s10100.set_leaf_bit(2);
        REQUIRE(s10100.as_string() == "10100");

        Split s01010;
        s01010.resize(5);
        s01010.set_leaf_bit(1);
        s01010.set_leaf_bit(3);
        REQUIRE(s01010.as_string() == "01010");

        REQUIRE(! (s10100 == s01010));
        REQUIRE(s10100 != s01010);
        REQUIRE(! s10100.is_equivalent(s01010, false));
        REQUIRE(! s10100.is_equivalent(s01010, true));
        REQUIRE(s10100.is_compatible(s01010));
        REQUIRE(! s10100.conflicts_with(s01010));

        Split s01011;
        s01011.resize(5);
        s01011.set_leaf_bit(1);
        s01011.set_leaf_bit(3);
        s01011.set_leaf_bit(4);
        REQUIRE(s01011.as_string() == "01011");

        REQUIRE(! (s10100 == s01011));
        REQUIRE(s10100 != s01011);
        REQUIRE(s10100.is_equivalent(s01011, false));
        REQUIRE(! s10100.is_equivalent(s01011, true));
        REQUIRE(s10100.is_compatible(s01011));
        REQUIRE(! s10100.conflicts_with(s01011));

        REQUIRE(! (s01010 == s01011));
        REQUIRE(s01010 != s01011);
        REQUIRE(! s01010.is_equivalent(s01011, false));
        REQUIRE(! s01010.is_equivalent(s01011, true));
        REQUIRE(s01010.is_compatible(s01011));
        REQUIRE(! s01010.conflicts_with(s01011));

        Split s11000;
        s11000.resize(5);
        s11000.set_leaf_bit(0);
        s11000.set_leaf_bit(1);
        REQUIRE(s11000.as_string() == "11000");

        REQUIRE(! (s11000 == s10100));
        REQUIRE(s11000 != s10100);
        REQUIRE(! s11000.is_equivalent(s10100, false));
        REQUIRE(! s11000.is_equivalent(s10100, true));
        REQUIRE(! s11000.is_compatible(s10100));
        REQUIRE(s11000.conflicts_with(s10100));

        REQUIRE(! (s11000 == s01010));
        REQUIRE(s11000 != s01010);
        REQUIRE(! s11000.is_equivalent(s01010, false));
        REQUIRE(! s11000.is_equivalent(s01010, true));
        REQUIRE(! s11000.is_compatible(s01010));
        REQUIRE(s11000.conflicts_with(s01010));

        REQUIRE(! (s11000 == s01011));
        REQUIRE(s11000 != s01011);
        REQUIRE(! s11000.is_equivalent(s01011, false));
        REQUIRE(! s11000.is_equivalent(s01011, true));
        REQUIRE(! s11000.is_compatible(s01011));
        REQUIRE(s11000.conflicts_with(s01011));
    }
}

TEST_CASE("Testing 100000 leaves", "[split]") {
    SECTION("Testing 100000 leaves") {
        unsigned int nleaves = 100000;
        Split s;
        s.resize(nleaves);
        s.set_leaf_bit(0);
        s.set_leaf_bit(nleaves - 1);

        std::ostringstream expected;
        expected << "1";
        for (unsigned int i = 0; i < nleaves - 2; ++i) {
            expected << "0";
        }
        expected << "1";
        REQUIRE(s.as_string() == expected.str());
    }
}
