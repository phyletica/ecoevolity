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
        REQUIRE(bd.has_constant_patterns() == false);
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

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
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
TEST_CASE("Testing diploid standard data with 012 as dominant", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-data-ntax5-nchar5.nex as dominant") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, '_', true, true, true), EcoevolityInvalidCharacterError);
    }
}

/**
 * #NEXUS
 * Begin data;
 *     Dimensions ntax=5 nchar=5;
 *     Format datatype=standard symbols="012" gap=-;
 *     Matrix
 * pop1_a  11100
 * pop1_b  00001
 * pop2_c  11101
 * pop2_d  10011
 * pop2_e  1-?01

 *     ;
 * End;
 */
TEST_CASE("Testing standard diploid with only 0/1 genotypes", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-only-01.nex") {
        std::string nex_path = "data/diploid-standard-only-01.nex";
        BiallelicData bd(nex_path);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);

        std::vector<unsigned int> expected_wts = {2,2,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(3);
        expected_allele_counts[0] = {4, 6};
        expected_allele_counts[1] = {4, 4};
        expected_allele_counts[2] = {4, 6};

        std::vector< std::vector<unsigned int> > expected_red_counts(3);
        expected_red_counts[0] = {1, 3};
        expected_red_counts[1] = {1, 1};
        expected_red_counts[2] = {0, 1};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);
    }
}

TEST_CASE("Testing standard haploid with a 2 genotype", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard-012.nex") {
        std::string nex_path = "data/haploid-standard-012.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, '_', true, false), NxsException);
    }

    SECTION("Testing data/haploid-standard-012.nex as diploid") {
        std::string nex_path = "data/haploid-standard-012.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, '_', true, true), NxsException);
    }

    SECTION("Testing data/haploid-standard-012.nex as dominant") {
        std::string nex_path = "data/haploid-standard-012.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, '_', true, false, true), EcoevolityBiallelicDataError);
    }

    SECTION("Testing data/haploid-standard-012.nex as diploid and dominant") {
        std::string nex_path = "data/haploid-standard-012.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, '_', true, true, true), NxsException);
    }
}

/**
 * #NEXUS
 * Begin data;
 *     Dimensions ntax=5 nchar=5;
 *     Format datatype=standard symbols="01" gap=-;
 *     Matrix
 * pop1_a  11100
 * pop1_b  00001
 * pop2_c  11101
 * pop2_d  10011
 * pop2_e  1-?00
 *     ;
 * End;
 */
TEST_CASE("Testing standard haploid", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard.nex") {
        std::string nex_path = "data/haploid-standard.nex";
        BiallelicData bd(nex_path, '_', true, false);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);

        std::vector<unsigned int> expected_wts = {1,2,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(4);
        expected_allele_counts[0] = {2, 3};
        expected_allele_counts[1] = {2, 2};
        expected_allele_counts[2] = {2, 3};
        expected_allele_counts[3] = {2, 3};

        std::vector< std::vector<unsigned int> > expected_red_counts(4);
        expected_red_counts[0] = {1, 3};
        expected_red_counts[1] = {1, 1};
        expected_red_counts[2] = {0, 1};
        expected_red_counts[3] = {1, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
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

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);
    }
}

/**
 * #NEXUS
 * Begin data;
 *     Dimensions ntax=5 nchar=5;
 *     Format datatype=standard symbols="01" gap=-;
 *     Matrix
 * pop1_a  11100
 * pop1_b  00001
 * pop2_c  11101
 * pop2_d  10011
 * pop2_e  1-?00
 *     ;
 * End;
 */
TEST_CASE("Testing standard haploid dominant", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard.nex as dominant") {
        std::string nex_path = "data/haploid-standard.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, '_', true, false, true), EcoevolityBiallelicDataError);
    }

    SECTION("Testing data/haploid-standard.nex as diploid and dominant") {
        std::string nex_path = "data/haploid-standard.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, '_', true, true, true), EcoevolityBiallelicDataError);
    }

    SECTION("Testing data/haploid-standard.nex as diploid") {
        std::string nex_path = "data/haploid-standard.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, '_', true, true, false), EcoevolityBiallelicDataError);
    }
}

/**
 * #NEXUS
 * Begin data;
 *     Dimensions ntax=5 nchar=5;
 *     Format datatype=standard symbols="012" gap=-;
 *     Matrix
 * pop1_a  2220000
 * pop1_b  0000200
 * pop2_c  2220202
 * pop2_d  2002200
 * pop2_e  2-?0220
 *     ;
 * End;
 */
TEST_CASE("Testing standard diploid dominant", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-dominant.nex as dominant") {
        std::string nex_path = "data/diploid-standard-dominant.nex";
        BiallelicData bd(nex_path, '_', true, true, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);

        std::vector<unsigned int> expected_wts = {2,2,3};

        std::vector< std::vector<unsigned int> > expected_allele_counts(3);
        expected_allele_counts[0] = {4, 6};
        expected_allele_counts[1] = {4, 4};
        expected_allele_counts[2] = {4, 6};

        std::vector< std::vector<unsigned int> > expected_red_counts(3);
        expected_red_counts[0] = {2, 6};
        expected_red_counts[1] = {2, 2};
        expected_red_counts[2] = {0, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);
    }
}

/**
 * #NEXUS
 * Begin data;
 *     Dimensions ntax=5 nchar=5;
 *     Format datatype=standard symbols="012" gap=-;
 *     Matrix
 * pop1_a  2220000
 * pop1_b  0000200
 * pop2_c  2220202
 * pop2_d  2002200
 * pop2_e  2-?0220
 *     ;
 * End;
 */
TEST_CASE("Testing standard diploid dominant as NOT dominant", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-dominant.nex as NOT dominant") {
        std::string nex_path = "data/diploid-standard-dominant.nex";
        BiallelicData bd(nex_path, '_', true, true, false);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);

        std::vector<unsigned int> expected_wts = {2,2,3};

        std::vector< std::vector<unsigned int> > expected_allele_counts(3);
        expected_allele_counts[0] = {4, 6};
        expected_allele_counts[1] = {4, 4};
        expected_allele_counts[2] = {4, 6};

        std::vector< std::vector<unsigned int> > expected_red_counts(3);
        expected_red_counts[0] = {2, 6};
        expected_red_counts[1] = {2, 2};
        expected_red_counts[2] = {0, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);
    }
}

TEST_CASE("Testing for constant diploid site patterns", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant0.nex") {
        std::string nex_path = "data/diploid-standard-constant0.nex";
        BiallelicData bd(nex_path);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);

        std::vector<unsigned int> expected_wts = {1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(6);
        expected_allele_counts[0] = {6, 4};
        expected_allele_counts[1] = {6, 2};
        expected_allele_counts[2] = {6, 4};
        expected_allele_counts[3] = {6, 4};
        expected_allele_counts[4] = {6, 4};
        expected_allele_counts[5] = {4, 2};

        std::vector< std::vector<unsigned int> > expected_red_counts(6);
        expected_red_counts[0] = {0, 0};
        expected_red_counts[1] = {1, 2};
        expected_red_counts[2] = {3, 3};
        expected_red_counts[3] = {0, 2};
        expected_red_counts[4] = {3, 4};
        expected_red_counts[5] = {0, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);

        // Removing patterns
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);

        std::vector<unsigned int> rm_expected_wts = {1,1,1,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(4);
        rm_expected_allele_counts[0] = {6, 2};
        rm_expected_allele_counts[1] = {6, 4};
        rm_expected_allele_counts[2] = {6, 4};
        rm_expected_allele_counts[3] = {6, 4};

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(4);
        rm_expected_red_counts[0] = {1, 2};
        rm_expected_red_counts[1] = {3, 3};
        rm_expected_red_counts[2] = {0, 2};
        rm_expected_red_counts[3] = {3, 4};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
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

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);
    }


    SECTION("Testing data/diploid-standard-constant2.nex") {
        std::string nex_path = "data/diploid-standard-constant2.nex";
        BiallelicData bd(nex_path);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);

        std::vector<unsigned int> expected_wts = {1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(6);
        expected_allele_counts[0] = {6, 4};
        expected_allele_counts[1] = {6, 2};
        expected_allele_counts[2] = {6, 4};
        expected_allele_counts[3] = {6, 4};
        expected_allele_counts[4] = {6, 4};
        expected_allele_counts[5] = {4, 2};

        std::vector< std::vector<unsigned int> > expected_red_counts(6);
        expected_red_counts[0] = {6, 4};
        expected_red_counts[1] = {1, 2};
        expected_red_counts[2] = {3, 3};
        expected_red_counts[3] = {0, 2};
        expected_red_counts[4] = {3, 4};
        expected_red_counts[5] = {4, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);

        // Removing patterns
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);

        std::vector<unsigned int> rm_expected_wts = {1,1,1,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(4);
        rm_expected_allele_counts[0] = {6, 2};
        rm_expected_allele_counts[1] = {6, 4};
        rm_expected_allele_counts[2] = {6, 4};
        rm_expected_allele_counts[3] = {6, 4};

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(4);
        rm_expected_red_counts[0] = {1, 2};
        rm_expected_red_counts[1] = {3, 3};
        rm_expected_red_counts[2] = {0, 2};
        rm_expected_red_counts[3] = {3, 4};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
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

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);
    }
}


TEST_CASE("Testing for constant haploid site patterns", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard-constant.nex") {
        std::string nex_path = "data/haploid-standard-constant.nex";
        BiallelicData bd(nex_path, '_', true, false);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);

        std::vector<unsigned int> expected_wts = {1,1,2,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {3, 2};
        expected_allele_counts[1] = {3, 1};
        expected_allele_counts[2] = {3, 2};
        expected_allele_counts[3] = {3, 2};
        expected_allele_counts[4] = {2, 1};

        std::vector< std::vector<unsigned int> > expected_red_counts(5);
        expected_red_counts[0] = {0, 0};
        expected_red_counts[1] = {1, 1};
        expected_red_counts[2] = {2, 2};
        expected_red_counts[3] = {0, 1};
        expected_red_counts[4] = {2, 1};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);

        // Removing patterns
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);

        std::vector<unsigned int> rm_expected_wts = {1,2,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(3);
        rm_expected_allele_counts[0] = {3, 1};
        rm_expected_allele_counts[1] = {3, 2};
        rm_expected_allele_counts[2] = {3, 2};

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(3);
        rm_expected_red_counts[0] = {1, 1};
        rm_expected_red_counts[1] = {2, 2};
        rm_expected_red_counts[2] = {0, 1};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);
    }

    SECTION("Testing data/haploid-standard-constant.nex as dominant") {
        std::string nex_path = "data/haploid-standard-constant.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, '_', true, false, true), EcoevolityBiallelicDataError);
    }

    SECTION("Testing data/haploid-standard-constant.nex as diploid") {
        std::string nex_path = "data/haploid-standard-constant.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, '_', true, true, false), EcoevolityBiallelicDataError);
    }

    SECTION("Testing data/haploid-standard-constant.nex as diploid and dominant") {
        std::string nex_path = "data/haploid-standard-constant.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, '_', true, true, true), EcoevolityBiallelicDataError);
    }
}



TEST_CASE("Testing for constant dominant diploid site patterns", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant-dominant.nex as dominant") {
        std::string nex_path = "data/diploid-standard-constant-dominant.nex";
        BiallelicData bd(nex_path, '_', true, true, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);

        std::vector<unsigned int> expected_wts = {1,1,2,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {6, 4};
        expected_allele_counts[1] = {6, 2};
        expected_allele_counts[2] = {6, 4};
        expected_allele_counts[3] = {6, 4};
        expected_allele_counts[4] = {4, 2};

        std::vector< std::vector<unsigned int> > expected_red_counts(5);
        expected_red_counts[0] = {6, 4};
        expected_red_counts[1] = {2, 2};
        expected_red_counts[2] = {4, 4};
        expected_red_counts[3] = {0, 2};
        expected_red_counts[4] = {0, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);

        // Removing patterns
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);

        std::vector<unsigned int> rm_expected_wts = {1,2,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(3);
        rm_expected_allele_counts[0] = {6, 2};
        rm_expected_allele_counts[1] = {6, 4};
        rm_expected_allele_counts[2] = {6, 4};

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(3);
        rm_expected_red_counts[0] = {2, 2};
        rm_expected_red_counts[1] = {4, 4};
        rm_expected_red_counts[2] = {0, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);
    }
}

TEST_CASE("Testing for missing haploid site patterns", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard-missing.nex") {
        std::string nex_path = "data/haploid-standard-missing.nex";
        BiallelicData bd(nex_path, '_', true, false);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);

        std::vector<unsigned int> expected_wts = {2,2,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {1, 0};
        expected_allele_counts[1] = {2, 2};
        expected_allele_counts[2] = {2, 3};
        expected_allele_counts[3] = {0, 0};
        expected_allele_counts[4] = {0, 3};

        std::vector< std::vector<unsigned int> > expected_red_counts(5);
        expected_red_counts[0] = {1, 0};
        expected_red_counts[1] = {1, 1};
        expected_red_counts[2] = {1, 2};
        expected_red_counts[3] = {0, 0};
        expected_red_counts[4] = {0, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);

        // Removing patterns
        unsigned int number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 3);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 2);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);

        std::vector<unsigned int> rm_expected_wts = {2,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(2);
        rm_expected_allele_counts[0] = {2, 2};
        rm_expected_allele_counts[1] = {2, 3};

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(2);
        rm_expected_red_counts[0] = {1, 1};
        rm_expected_red_counts[1] = {1, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(2), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(2), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(2), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);
    }
}

TEST_CASE("Testing for constant AND missing haploid site patterns", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard-missing.nex") {
        std::string nex_path = "data/haploid-standard-missing.nex";
        BiallelicData bd(nex_path, '_', true, false);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);

        std::vector<unsigned int> expected_wts = {2,2,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {1, 0};
        expected_allele_counts[1] = {2, 2};
        expected_allele_counts[2] = {2, 3};
        expected_allele_counts[3] = {0, 0};
        expected_allele_counts[4] = {0, 3};

        std::vector< std::vector<unsigned int> > expected_red_counts(5);
        expected_red_counts[0] = {1, 0};
        expected_red_counts[1] = {1, 1};
        expected_red_counts[2] = {1, 2};
        expected_red_counts[3] = {0, 0};
        expected_red_counts[4] = {0, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);

        // Removing patterns
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);

        std::vector<unsigned int> rm_expected_wts = {2,1,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(3);
        rm_expected_allele_counts[0] = {2, 2};
        rm_expected_allele_counts[1] = {2, 3};
        rm_expected_allele_counts[2] = {0, 3};

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(3);
        rm_expected_red_counts[0] = {1, 1};
        rm_expected_red_counts[1] = {1, 2};
        rm_expected_red_counts[2] = {0, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);

        // Removing patterns
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 2);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);

        std::vector<unsigned int> rm_rm_expected_wts = {2,1};

        std::vector< std::vector<unsigned int> > rm_rm_expected_allele_counts(2);
        rm_rm_expected_allele_counts[0] = {2, 2};
        rm_rm_expected_allele_counts[1] = {2, 3};

        std::vector< std::vector<unsigned int> > rm_rm_expected_red_counts(2);
        rm_rm_expected_red_counts[0] = {1, 1};
        rm_rm_expected_red_counts[1] = {1, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(2), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(2), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(2), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> rm_rm_expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == rm_rm_expected_labels);
        rm_rm_expected_labels.clear();
        rm_rm_expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);
    }
}

TEST_CASE("Testing for constant AND missing dominant haploid site patterns", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard-missing.nex as dominant") {
        std::string nex_path = "data/haploid-standard-missing.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, '_', true, false, true), EcoevolityBiallelicDataError);
    }

    SECTION("Testing data/haploid-standard-missing.nex as diploid") {
        std::string nex_path = "data/haploid-standard-missing.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, '_', true, true, false), EcoevolityBiallelicDataError);
    }

    SECTION("Testing data/haploid-standard-missing.nex as diploid and dominant") {
        std::string nex_path = "data/haploid-standard-missing.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, '_', true, true, true), EcoevolityBiallelicDataError);
    }
}

TEST_CASE("Testing for constant AND missing diploid site patterns", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-missing.nex") {
        std::string nex_path = "data/diploid-standard-missing.nex";
        BiallelicData bd(nex_path);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);

        std::vector<unsigned int> expected_wts = {3,1,2,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {6, 0};
        expected_allele_counts[1] = {0, 4};
        expected_allele_counts[2] = {6, 4};
        expected_allele_counts[3] = {6, 2};
        expected_allele_counts[4] = {6, 2};

        std::vector< std::vector<unsigned int> > expected_red_counts(5);
        expected_red_counts[0] = {3, 0};
        expected_red_counts[1] = {0, 3};
        expected_red_counts[2] = {3, 3};
        expected_red_counts[3] = {0, 0};
        expected_red_counts[4] = {6, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);

        // Removing patterns
        unsigned int number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);

        std::vector<unsigned int> rm_expected_wts = {2,1,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(3);
        rm_expected_allele_counts[0] = {6, 4};
        rm_expected_allele_counts[1] = {6, 2};
        rm_expected_allele_counts[2] = {6, 2};

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(3);
        rm_expected_red_counts[0] = {3, 3};
        rm_expected_red_counts[1] = {0, 0};
        rm_expected_red_counts[2] = {6, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);

        // Removing patterns
        number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 1);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);

        std::vector<unsigned int> rm_rm_expected_wts = {2};

        std::vector< std::vector<unsigned int> > rm_rm_expected_allele_counts(1);
        rm_rm_expected_allele_counts[0] = {6, 4};

        std::vector< std::vector<unsigned int> > rm_rm_expected_red_counts(1);
        rm_rm_expected_red_counts[0] = {3, 3};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(1), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_allele_counts(1), std::out_of_range);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(1), std::out_of_range);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range);

        std::vector<std::string> rm_rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_rm_expected_labels);
        rm_rm_expected_labels.clear();
        rm_rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range);
    }
}

TEST_CASE("Testing for removing all site patterns", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-all-remove.nex") {
        std::string nex_path = "data/diploid-standard-all-remove.nex";
        BiallelicData bd(nex_path);
        int number_removed = bd.remove_constant_patterns();
        REQUIRE_THROWS_AS(bd.remove_missing_population_patterns(), EcoevolityBiallelicDataError);
    }
}
