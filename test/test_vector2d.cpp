#include "catch.hpp"
#include "ecoevolity/vector2d.hpp"


TEST_CASE("Testing nrows, ncols constructor of Vector2d", "[Vector2d]") {

    SECTION("Testing 0x0 constructor") {
        Vector2d m(0,0);
        REQUIRE(m.get_nrows() == 0);
        REQUIRE(m.get_ncols() == 0);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,1), std::out_of_range);
    }

    SECTION("Testing 1x0 constructor") {
        Vector2d m(1,0);
        REQUIRE(m.get_nrows() == 1);
        REQUIRE(m.get_ncols() == 0);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,1), std::out_of_range);
    }

    SECTION("Testing 0x1 constructor") {
        Vector2d m(0,1);
        REQUIRE(m.get_nrows() == 0);
        REQUIRE(m.get_ncols() == 1);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,1), std::out_of_range);
    }

    SECTION("Testing 1x1 constructor") {
        Vector2d m(1,1);
        REQUIRE(m.get_nrows() == 1);
        REQUIRE(m.get_ncols() == 1);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);
        REQUIRE(m.at(1,1) == Approx(0.0));
        REQUIRE_THROWS_AS(m.at(2,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,2), std::out_of_range);
    }

    SECTION("Testing 2x4 constructor") {
        Vector2d m(2,4);
        REQUIRE(m.get_nrows() == 2);
        REQUIRE(m.get_ncols() == 4);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);

        REQUIRE(m.at(1,1) == Approx(0.0));
        REQUIRE(m.at(1,2) == Approx(0.0));
        REQUIRE(m.at(1,3) == Approx(0.0));
        REQUIRE(m.at(1,4) == Approx(0.0));
        REQUIRE(m.at(2,1) == Approx(0.0));
        REQUIRE(m.at(2,2) == Approx(0.0));
        REQUIRE(m.at(2,3) == Approx(0.0));
        REQUIRE(m.at(2,4) == Approx(0.0));
        REQUIRE_THROWS_AS(m.at(3,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,5), std::out_of_range);
    }

    SECTION("Testing 4x2 constructor") {
        Vector2d m(4,2);
        REQUIRE(m.get_nrows() == 4);
        REQUIRE(m.get_ncols() == 2);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);

        REQUIRE(m.at(1,1) == Approx(0.0));
        REQUIRE(m.at(1,2) == Approx(0.0));
        REQUIRE(m.at(2,1) == Approx(0.0));
        REQUIRE(m.at(2,2) == Approx(0.0));
        REQUIRE(m.at(3,1) == Approx(0.0));
        REQUIRE(m.at(3,2) == Approx(0.0));
        REQUIRE(m.at(4,1) == Approx(0.0));
        REQUIRE(m.at(4,2) == Approx(0.0));
        REQUIRE_THROWS_AS(m.at(5,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,3), std::out_of_range);
    }
}

TEST_CASE("Testing square identity constructor of Vector2d", "[Vector2d]") {

    SECTION("Testing 0x0 square constructor") {
        Vector2d m(0);
        REQUIRE(m.get_nrows() == 0);
        REQUIRE(m.get_ncols() == 0);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,1), std::out_of_range);
    }

    SECTION("Testing 1x1 square constructor") {
        Vector2d m(1);
        REQUIRE(m.get_nrows() == 1);
        REQUIRE(m.get_ncols() == 1);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);
        REQUIRE(m.at(1,1) == Approx(1.0));
        REQUIRE_THROWS_AS(m.at(2,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,2), std::out_of_range);
    }

    SECTION("Testing 2x2 square constructor") {
        Vector2d m(2);
        REQUIRE(m.get_nrows() == 2);
        REQUIRE(m.get_ncols() == 2);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);

        REQUIRE(m.at(1,1) == Approx(1.0));
        REQUIRE(m.at(1,2) == Approx(0.0));
        REQUIRE(m.at(2,1) == Approx(0.0));
        REQUIRE(m.at(2,2) == Approx(1.0));
        REQUIRE_THROWS_AS(m.at(3,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,3), std::out_of_range);
    }

    SECTION("Testing 3x3 square constructor") {
        Vector2d m(3);
        REQUIRE(m.get_nrows() == 3);
        REQUIRE(m.get_ncols() == 3);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);

        REQUIRE(m.at(1,1) == Approx(1.0));
        REQUIRE(m.at(1,2) == Approx(0.0));
        REQUIRE(m.at(1,3) == Approx(0.0));
        REQUIRE(m.at(2,1) == Approx(0.0));
        REQUIRE(m.at(2,2) == Approx(1.0));
        REQUIRE(m.at(2,3) == Approx(0.0));
        REQUIRE(m.at(3,1) == Approx(0.0));
        REQUIRE(m.at(3,2) == Approx(0.0));
        REQUIRE(m.at(3,3) == Approx(1.0));
        REQUIRE_THROWS_AS(m.at(4,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,4), std::out_of_range);
    }
}

TEST_CASE("Testing value constructor of Vector2d", "[Vector2d]") {

    SECTION("Testing 2x3 value constructor") {
        std::vector<double> values = {0.11, 0.12, 0.13, 0.21, 0.22, 0.23};
        Vector2d m(2, 3, values);
        REQUIRE(m.get_nrows() == 2);
        REQUIRE(m.get_ncols() == 3);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);

        REQUIRE(m.at(1,1) == Approx(0.11));
        REQUIRE(m.at(1,2) == Approx(0.12));
        REQUIRE(m.at(1,3) == Approx(0.13));
        REQUIRE(m.at(2,1) == Approx(0.21));
        REQUIRE(m.at(2,2) == Approx(0.22));
        REQUIRE(m.at(2,3) == Approx(0.23));
        REQUIRE_THROWS_AS(m.at(3,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,4), std::out_of_range);
    }

    SECTION("Testing 3x2 value constructor") {
        std::vector<double> values = {0.11, 0.12, 0.21, 0.22, 0.31, 0.32};
        Vector2d m(3, 2, values);
        REQUIRE(m.get_nrows() == 3);
        REQUIRE(m.get_ncols() == 2);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);

        REQUIRE(m.at(1,1) == Approx(0.11));
        REQUIRE(m.at(1,2) == Approx(0.12));
        REQUIRE(m.at(2,1) == Approx(0.21));
        REQUIRE(m.at(2,2) == Approx(0.22));
        REQUIRE(m.at(3,1) == Approx(0.31));
        REQUIRE(m.at(3,2) == Approx(0.32));
        REQUIRE_THROWS_AS(m.at(4,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,3), std::out_of_range);
    }
}

TEST_CASE("Testing setter of Vector2d", "[Vector2d]") {

    SECTION("Testing setter") {
        Vector2d m(2, 3);
        REQUIRE(m.get_nrows() == 2);
        REQUIRE(m.get_ncols() == 3);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(3,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,4), std::out_of_range);

        REQUIRE(m.at(1,1) == Approx(0.0));
        REQUIRE(m.at(1,2) == Approx(0.0));
        REQUIRE(m.at(1,3) == Approx(0.0));
        REQUIRE(m.at(2,1) == Approx(0.0));
        REQUIRE(m.at(2,2) == Approx(0.0));
        REQUIRE(m.at(2,3) == Approx(0.0));

        m.set(1, 1, 0.11);
        m.set(1, 2, 0.12);
        m.set(1, 3, 0.13);
        m.set(2, 1, 0.21);
        m.set(2, 2, 0.22);
        m.set(2, 3, 0.23);

        REQUIRE(m.at(1,1) == Approx(0.11));
        REQUIRE(m.at(1,2) == Approx(0.12));
        REQUIRE(m.at(1,3) == Approx(0.13));
        REQUIRE(m.at(2,1) == Approx(0.21));
        REQUIRE(m.at(2,2) == Approx(0.22));
        REQUIRE(m.at(2,3) == Approx(0.23));
    }
}

TEST_CASE("Testing copy operator of Vector2d", "[Vector2d]") {

    SECTION("Testing copy operator") {
        Vector2d m(2, 3);
        m.set(2, 2, -0.3);
        Vector2d m2 = m;
        m2.set(2, 2, 1.5);

        REQUIRE(m.at(2,2) == Approx(-0.3));
        REQUIRE(m2.at(2,2) == Approx(1.5));

        m.resize(1,1);
        REQUIRE(m.get_nrows() == 1);
        REQUIRE(m.get_ncols() == 1);
        REQUIRE(m2.get_nrows() == 2);
        REQUIRE(m2.get_ncols() == 3);
    }
}

TEST_CASE("Testing transpose of Vector2d", "[Vector2d]") {
    SECTION("Testing transpose") {
        std::vector<double> values = {0.11, 0.12, 0.13, 0.21, 0.22, 0.23};
        Vector2d m(2, 3, values);
        REQUIRE(m.get_nrows() == 2);
        REQUIRE(m.get_ncols() == 3);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);

        REQUIRE(m.at(1,1) == Approx(0.11));
        REQUIRE(m.at(1,2) == Approx(0.12));
        REQUIRE(m.at(1,3) == Approx(0.13));
        REQUIRE(m.at(2,1) == Approx(0.21));
        REQUIRE(m.at(2,2) == Approx(0.22));
        REQUIRE(m.at(2,3) == Approx(0.23));
        REQUIRE_THROWS_AS(m.at(3,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,4), std::out_of_range);

        Vector2d m2 = m.transpose();
        REQUIRE(m2.get_nrows() == 3);
        REQUIRE(m2.get_ncols() == 2);
        REQUIRE_THROWS_AS(m2.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m2.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m2.at(0,1), std::out_of_range);

        REQUIRE(m2.at(1,1) == Approx(0.11));
        REQUIRE(m2.at(1,2) == Approx(0.21));
        REQUIRE(m2.at(2,1) == Approx(0.12));
        REQUIRE(m2.at(2,2) == Approx(0.22));
        REQUIRE(m2.at(3,1) == Approx(0.13));
        REQUIRE(m2.at(3,2) == Approx(0.23));
        REQUIRE_THROWS_AS(m2.at(4,1), std::out_of_range);
        REQUIRE_THROWS_AS(m2.at(1,3), std::out_of_range);

        // Make sure original was unaffected.
        REQUIRE(m.get_nrows() == 2);
        REQUIRE(m.get_ncols() == 3);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);

        REQUIRE(m.at(1,1) == Approx(0.11));
        REQUIRE(m.at(1,2) == Approx(0.12));
        REQUIRE(m.at(1,3) == Approx(0.13));
        REQUIRE(m.at(2,1) == Approx(0.21));
        REQUIRE(m.at(2,2) == Approx(0.22));
        REQUIRE(m.at(2,3) == Approx(0.23));
        REQUIRE_THROWS_AS(m.at(3,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,4), std::out_of_range);
    }
}

TEST_CASE("Testing cell math methods Vector2d", "[Vector2d]") {
    SECTION("Testing math methods") {
        Vector2d m(1, 2);
        REQUIRE(m.at(1,1) == Approx(0.0));
        REQUIRE(m.at(1,2) == Approx(0.0));

        m.add(1, 1, 1.0);
        REQUIRE(m.at(1,1) == Approx(1.0));
        REQUIRE(m.at(1,2) == Approx(0.0));

        m.multiply(1, 1, 4.0);
        REQUIRE(m.at(1,1) == Approx(4.0));
        REQUIRE(m.at(1,2) == Approx(0.0));
        
        m.divide(1, 1, 0.5);
        REQUIRE(m.at(1,1) == Approx(8.0));
        REQUIRE(m.at(1,2) == Approx(0.0));

        m.subtract(1, 1, 6.0);
        REQUIRE(m.at(1,1) == Approx(2.0));
        REQUIRE(m.at(1,2) == Approx(0.0));

        std::vector<double> addends = {2.0, 3.0};
        m.add(0.0, addends);

        REQUIRE(m.at(1,1) == Approx(2.0));
        REQUIRE(m.at(1,2) == Approx(0.0));

        m.add(1.0, addends);
        REQUIRE(m.at(1,1) == Approx(4.0));
        REQUIRE(m.at(1,2) == Approx(3.0));

        m.add(-0.5, addends);
        REQUIRE(m.at(1,1) == Approx(3.0));
        REQUIRE(m.at(1,2) == Approx(1.5));
    }
}

TEST_CASE("Testing scale method of Vector2d", "[Vector2d]") {
    SECTION("Testing scale method") {
        Vector2d m(1, 2);
        m.set(1,1, 2.0);
        m.set(1,2, 4.0);
        REQUIRE(m.at(1,1) == Approx(2.0));
        REQUIRE(m.at(1,2) == Approx(4.0));

        m.scale(2.0);
        REQUIRE(m.at(1,1) == Approx(4.0));
        REQUIRE(m.at(1,2) == Approx(8.0));

        m.scale(-0.25);
        REQUIRE(m.at(1,1) == Approx(-1.0));
        REQUIRE(m.at(1,2) == Approx(-2.0));
    }
}

TEST_CASE("Testing get column methods Vector2d", "[Vector2d]") {
    SECTION("Testing get column copy") {
        std::vector<double> values = {0.11, 0.12, 0.13, 0.21, 0.22, 0.23};
        Vector2d m(2, 3, values);
        REQUIRE_THROWS_AS(m.get_column(0), std::out_of_range);
        REQUIRE_THROWS_AS(m.get_column(4), std::out_of_range);

        std::vector<double> c = m.get_column(3);
        REQUIRE(c.size() == 3);
        REQUIRE(c.at(0) == Approx(0.0));
        REQUIRE(c.at(1) == Approx(0.13));
        REQUIRE(c.at(2) == Approx(0.23));

        c = m.get_column(2);
        REQUIRE(c.size() == 3);
        REQUIRE(c.at(0) == Approx(0.0));
        REQUIRE(c.at(1) == Approx(0.12));
        REQUIRE(c.at(2) == Approx(0.22));

        c = m.get_column(1);
        REQUIRE(c.size() == 3);
        REQUIRE(c.at(0) == Approx(0.0));
        REQUIRE(c.at(1) == Approx(0.11));
        REQUIRE(c.at(2) == Approx(0.21));
    }

    SECTION("Testing get column fill") {
        std::vector<double> values = {0.11, 0.12, 0.13, 0.21, 0.22, 0.23};
        Vector2d m(2, 3, values);
        REQUIRE_THROWS_AS(m.get_column(0), std::out_of_range);
        REQUIRE_THROWS_AS(m.get_column(4), std::out_of_range);

        std::vector<double> c = {-1.0};
        m.get_column(3, c);
        REQUIRE(c.size() == 3);
        REQUIRE(c.at(0) == Approx(-1.0));
        REQUIRE(c.at(1) == Approx(0.13));
        REQUIRE(c.at(2) == Approx(0.23));

        m.get_column(2, c);
        REQUIRE(c.size() == 3);
        REQUIRE(c.at(0) == Approx(-1.0));
        REQUIRE(c.at(1) == Approx(0.12));
        REQUIRE(c.at(2) == Approx(0.22));

        m.get_column(1, c);
        REQUIRE(c.size() == 3);
        REQUIRE(c.at(0) == Approx(-1.0));
        REQUIRE(c.at(1) == Approx(0.11));
        REQUIRE(c.at(2) == Approx(0.21));
    }
}

TEST_CASE("Testing resize of Vector2d", "[Vector2d]") {
    SECTION("Testing resize") {
        std::vector<double> values = {0.11, 0.12, 0.13, 0.21, 0.22, 0.23};
        Vector2d m(2, 3, values);
        REQUIRE(m.get_nrows() == 2);
        REQUIRE(m.get_ncols() == 3);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);

        REQUIRE(m.at(1,1) == Approx(0.11));
        REQUIRE(m.at(1,2) == Approx(0.12));
        REQUIRE(m.at(1,3) == Approx(0.13));
        REQUIRE(m.at(2,1) == Approx(0.21));
        REQUIRE(m.at(2,2) == Approx(0.22));
        REQUIRE(m.at(2,3) == Approx(0.23));
        REQUIRE_THROWS_AS(m.at(3,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,4), std::out_of_range);

        m.resize(3,4);
        REQUIRE(m.get_nrows() == 3);
        REQUIRE(m.get_ncols() == 4);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);

        REQUIRE(m.at(1,1) == Approx(0.11));
        REQUIRE(m.at(1,2) == Approx(0.12));
        REQUIRE(m.at(1,3) == Approx(0.13));
        REQUIRE(m.at(1,4) == Approx(0.0));
        REQUIRE(m.at(2,1) == Approx(0.21));
        REQUIRE(m.at(2,2) == Approx(0.22));
        REQUIRE(m.at(2,3) == Approx(0.23));
        REQUIRE(m.at(2,4) == Approx(0.0));
        REQUIRE(m.at(3,1) == Approx(0.0));
        REQUIRE(m.at(3,2) == Approx(0.0));
        REQUIRE(m.at(3,3) == Approx(0.0));
        REQUIRE(m.at(3,4) == Approx(0.0));
        REQUIRE_THROWS_AS(m.at(4,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,5), std::out_of_range);

        m.resize(3, 2);
        REQUIRE(m.get_nrows() == 3);
        REQUIRE(m.get_ncols() == 2);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);

        REQUIRE(m.at(1,1) == Approx(0.11));
        REQUIRE(m.at(1,2) == Approx(0.12));
        REQUIRE(m.at(2,1) == Approx(0.21));
        REQUIRE(m.at(2,2) == Approx(0.22));
        REQUIRE(m.at(3,1) == Approx(0.0));
        REQUIRE(m.at(3,2) == Approx(0.0));
        REQUIRE_THROWS_AS(m.at(4,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,3), std::out_of_range);

        m.resize(2, 1);
        REQUIRE(m.get_nrows() == 2);
        REQUIRE(m.get_ncols() == 1);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);

        REQUIRE(m.at(1,1) == Approx(0.11));
        REQUIRE(m.at(2,1) == Approx(0.21));
        REQUIRE_THROWS_AS(m.at(3,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,2), std::out_of_range);
    }
}

TEST_CASE("Testing get zero based vector of Vector2d", "[Vector2d]") {
    SECTION("Testing get_zero_based_vector") {
        std::vector<double> values = {0.11, 0.12, 0.13, 0.21, 0.22, 0.23};
        Vector2d m(2, 3, values);
        REQUIRE(m.get_nrows() == 2);
        REQUIRE(m.get_ncols() == 3);
        REQUIRE_THROWS_AS(m.at(0,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,0), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(0,1), std::out_of_range);

        REQUIRE(m.at(1,1) == Approx(0.11));
        REQUIRE(m.at(1,2) == Approx(0.12));
        REQUIRE(m.at(1,3) == Approx(0.13));
        REQUIRE(m.at(2,1) == Approx(0.21));
        REQUIRE(m.at(2,2) == Approx(0.22));
        REQUIRE(m.at(2,3) == Approx(0.23));
        REQUIRE_THROWS_AS(m.at(3,1), std::out_of_range);
        REQUIRE_THROWS_AS(m.at(1,4), std::out_of_range);

        std::vector<double> v = m.get_zero_based_vector();
        REQUIRE(v.size() == values.size());
        for (unsigned int i = 0; i < v.size(); ++i) {
            REQUIRE(v.at(i) == Approx(values.at(i)));
        }

        v.at(0) = 3.1;
        REQUIRE(v.at(0) == 3.1);
        REQUIRE(values.at(0) == 0.11);
    }
}

