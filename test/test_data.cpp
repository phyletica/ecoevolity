#include "catch.hpp"
#include "ecoevolity/data.hpp"

/**
 *
 * #NEXUS
 * Begin data;
 *     Dimensions ntax=5 nchar=5;
 *     Format datatype=standard symbols="012" gap=-;
 *     Matrix
 * pop1_a  00000
 * pop1_b  00202
 * pop1_c  01101
 * pop2_c  2-122
 * pop2_d  02202
 *     ;
 * End;
 */
TEST_CASE("Testing small, diploid, standard data set", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-data-ntax5-nchar5.nex") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        BiallelicData bd(nex_path);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_site_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);

        std::vector<unsigned int> expected_wts = {2,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(4);
        expected_allele_counts[0] = {6, 4};
        expected_allele_counts[1] = {6, 2};
        expected_allele_counts[2] = {6, 4};
        expected_allele_counts[3] = {6, 4};

        std::vector< std::vector<unsigned int> > expected_red_counts(4);
        expected_red_counts[0] = {0, 2};
        expected_red_counts[1] = {1, 2};
        expected_red_counts[2] = {3, 3};
        expected_red_counts[3] = {3, 4};

        for (unsigned int pattern_idx = 0; pattern_idx < bd.get_number_of_patterns(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);
    }
}

TEST_CASE("Testing dominant data het error", "[BialleleicData]") {

    SECTION("Testing for error if dominant data has het genotype codes") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, '_', true, true), EcoevolityInvalidCharacterError);
    }
}
