#include "catch.hpp"
#include "ecoevolity/matrix.hpp"


TEST_CASE("Testing constructors of BiallelicPatternProbabilityMatrix",
        "[BiallelicPatternProbabilityMatrix]") {

    SECTION("Testing bare constructor") {
        BiallelicPatternProbabilityMatrix m;
        REQUIRE(m.get_allele_count() == 0);
        std::vector<double> expected_prob_matrix;
        REQUIRE(m.get_pattern_prob_matrix() == expected_prob_matrix);
        REQUIRE_THROWS_AS(m.get_pattern_probability(0, 0), std::out_of_range);
        REQUIRE_THROWS_AS(m.set_pattern_probability(0, 0, 1.0), std::out_of_range);
        REQUIRE_THROWS_AS(m.get_pattern_probability(1, 0), std::out_of_range);
        REQUIRE_THROWS_AS(m.set_pattern_probability(1, 0, 1.0), std::out_of_range);

        m.resize(1);
        REQUIRE(m.get_allele_count() == 1);
        REQUIRE(m.get_pattern_probability(1, 0) == 0.0);
        REQUIRE(m.get_pattern_probability(1, 1) == 0.0);
        expected_prob_matrix.push_back(0.0);
        expected_prob_matrix.push_back(0.0);
        REQUIRE(m.get_pattern_prob_matrix() == expected_prob_matrix);

        m.set_pattern_probability(1, 0, 1.0);
        REQUIRE(m.get_pattern_probability(1, 0) == 1.0);
        REQUIRE(m.get_pattern_probability(1, 1) == 0.0);
        expected_prob_matrix.at(0) = 1.0;
        REQUIRE(m.get_pattern_prob_matrix() == expected_prob_matrix);

        m.set_pattern_probability(1, 1, 1.0);
        REQUIRE(m.get_pattern_probability(1, 0) == 1.0);
        REQUIRE(m.get_pattern_probability(1, 1) == 1.0);
        expected_prob_matrix.at(1) = 1.0;
        REQUIRE(m.get_pattern_prob_matrix() == expected_prob_matrix);

        m.resize(1);
        REQUIRE(m.get_allele_count() == 1);
        REQUIRE(m.get_pattern_prob_matrix() == expected_prob_matrix);

        REQUIRE_THROWS_AS(m.get_pattern_probability(0, 0), std::out_of_range);
        REQUIRE_THROWS_AS(m.set_pattern_probability(0, 0, 1.0), std::out_of_range);
        REQUIRE_THROWS_AS(m.get_pattern_probability(2, 0), std::out_of_range);
        REQUIRE_THROWS_AS(m.set_pattern_probability(2, 0, 1.0), std::out_of_range);

        m.reset(2);
        REQUIRE(m.get_allele_count() == 2);
        expected_prob_matrix.clear();
        expected_prob_matrix = {0.0, 0.0, 0.0, 0.0, 0.0};
        REQUIRE(m.get_pattern_prob_matrix() == expected_prob_matrix);
        REQUIRE(m.get_pattern_probability(1, 0) == 0.0);
        REQUIRE(m.get_pattern_probability(1, 1) == 0.0);
        REQUIRE(m.get_pattern_probability(2, 0) == 0.0);
        REQUIRE(m.get_pattern_probability(2, 1) == 0.0);
        REQUIRE(m.get_pattern_probability(2, 2) == 0.0);

        REQUIRE_THROWS_AS(m.get_pattern_probability(0, 0), std::out_of_range);
        REQUIRE_THROWS_AS(m.set_pattern_probability(0, 0, 1.0), std::out_of_range);
        REQUIRE_THROWS_AS(m.get_pattern_probability(3, 0), std::out_of_range);
        REQUIRE_THROWS_AS(m.set_pattern_probability(3, 0, 1.0), std::out_of_range);

        BiallelicPatternProbabilityMatrix m2;
        m2.copy(m);

        REQUIRE(m2.get_allele_count() == 2);
        expected_prob_matrix.clear();
        expected_prob_matrix = {0.0, 0.0, 0.0, 0.0, 0.0};
        REQUIRE(m2.get_pattern_prob_matrix() == expected_prob_matrix);
        REQUIRE(m2.get_pattern_probability(1, 0) == 0.0);
        REQUIRE(m2.get_pattern_probability(1, 1) == 0.0);
        REQUIRE(m2.get_pattern_probability(2, 0) == 0.0);
        REQUIRE(m2.get_pattern_probability(2, 1) == 0.0);
        REQUIRE(m2.get_pattern_probability(2, 2) == 0.0);

        REQUIRE_THROWS_AS(m2.get_pattern_probability(0, 0), std::out_of_range);
        REQUIRE_THROWS_AS(m2.set_pattern_probability(0, 0, 1.0), std::out_of_range);
        REQUIRE_THROWS_AS(m2.get_pattern_probability(3, 0), std::out_of_range);
        REQUIRE_THROWS_AS(m2.set_pattern_probability(3, 0, 1.0), std::out_of_range);

        m2.set_pattern_probability(2, 2, 1.0);
        REQUIRE(m2.get_pattern_probability(2, 2) == 1.0);
        REQUIRE(m.get_pattern_probability(2, 2) == 0.0);

        m2.resize(1);
        REQUIRE(m2.get_allele_count() == 1);
        REQUIRE(m.get_allele_count() == 2);

    }

    SECTION("Testing allele count constructor") {
        BiallelicPatternProbabilityMatrix m(1);
        REQUIRE(m.get_allele_count() == 1);
        std::vector<double> expected_prob_matrix = {0.0, 0.0};
        REQUIRE(m.get_pattern_prob_matrix() == expected_prob_matrix);
        REQUIRE(m.get_pattern_probability(1, 0) == 0.0);
        REQUIRE(m.get_pattern_probability(1, 1) == 0.0);

        REQUIRE_THROWS_AS(m.get_pattern_probability(0, 0), std::out_of_range);
        REQUIRE_THROWS_AS(m.set_pattern_probability(0, 0, 1.0), std::out_of_range);
        REQUIRE_THROWS_AS(m.get_pattern_probability(2, 0), std::out_of_range);
        REQUIRE_THROWS_AS(m.set_pattern_probability(2, 0, 1.0), std::out_of_range);
    }

    SECTION("Testing copy constructor") {
        BiallelicPatternProbabilityMatrix m(2);
        REQUIRE(m.get_allele_count() == 2);
        m.set_pattern_probability(2, 2, 1.0);
        std::vector<double> expected_prob_matrix = {0.0, 0.0, 0.0, 0.0, 1.0};
        REQUIRE(m.get_pattern_prob_matrix() == expected_prob_matrix);

        BiallelicPatternProbabilityMatrix m2;
        m2.copy(m);
        REQUIRE(m2.get_allele_count() == 2);
        REQUIRE(m2.get_pattern_prob_matrix() == expected_prob_matrix);

        m.set_pattern_probability(2, 1, 1.0);
        REQUIRE(m2.get_pattern_prob_matrix() == expected_prob_matrix);
        expected_prob_matrix.clear();
        expected_prob_matrix = {0.0, 0.0, 0.0, 1.0, 1.0};
        REQUIRE(m.get_pattern_prob_matrix() == expected_prob_matrix);
    }

    SECTION("Testing allele count and vector constructor") {
        std::vector<double> expected_prob_matrix = {0.0, 0.5, 0.0, 0.0, 0.5};
        BiallelicPatternProbabilityMatrix m(2, expected_prob_matrix);
        REQUIRE(m.get_allele_count() == 2);
        REQUIRE(m.get_pattern_prob_matrix() == expected_prob_matrix);
        expected_prob_matrix.at(1) = 0.0;
        expected_prob_matrix.at(4) = 1.0;
        REQUIRE(m.get_pattern_prob_matrix() != expected_prob_matrix);

        REQUIRE_THROWS_AS(BiallelicPatternProbabilityMatrix m2(1, expected_prob_matrix),
                EcoevolityError);
        REQUIRE_THROWS_AS(BiallelicPatternProbabilityMatrix m3(3, expected_prob_matrix),
                EcoevolityError);
    }
}