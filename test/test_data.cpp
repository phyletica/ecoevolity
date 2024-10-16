#include "catch.hpp"
#include "ecoevolity/data.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include "ecoevolity/rng.hpp"

RandomNumberGenerator _ECOEVOLITY_DATA_RNG = RandomNumberGenerator();

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
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(bd.get_number_of_variable_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);
        REQUIRE(bd.has_seq_loci_info() == false);

        std::vector<unsigned int> expected_wts = {2,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(4);
        expected_allele_counts[0] = {6, 4};
        expected_allele_counts[1] = {6, 2};
        expected_allele_counts[2] = {6, 4};
        expected_allele_counts[3] = {6, 4};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 4;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(4);
        expected_red_counts[0] = {0, 2};
        expected_red_counts[1] = {1, 2};
        expected_red_counts[2] = {3, 3};
        expected_red_counts[3] = {3, 4};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
            for (unsigned int pop_idx = 0; pop_idx < bd.get_number_of_populations(); ++pop_idx) {
                REQUIRE(bd.get_allele_count(pattern_idx, pop_idx) ==
                        expected_allele_counts.at(pattern_idx).at(pop_idx));
                REQUIRE(bd.get_red_allele_count(pattern_idx, pop_idx) ==
                        expected_red_counts.at(pattern_idx).at(pop_idx));
            }
        }

        std::vector<unsigned int> expected_max_cts = {6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        //fold patterns
        unsigned int number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(bd.get_number_of_variable_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_seq_loci_info() == false);

        expected_red_counts[2] = {3, 1};
        expected_red_counts[3] = {3, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }
    }
}

TEST_CASE("Testing small, diploid, standard data set with charsets", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-data-ntax5-nchar5.nex") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        BiallelicData bd(nex_path, ' ', true, true, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(bd.get_number_of_variable_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {2,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(4);
        expected_allele_counts[0] = {6, 4};
        expected_allele_counts[1] = {6, 2};
        expected_allele_counts[2] = {6, 4};
        expected_allele_counts[3] = {6, 4};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 4;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(4);
        expected_red_counts[0] = {0, 2};
        expected_red_counts[1] = {1, 2};
        expected_red_counts[2] = {3, 3};
        expected_red_counts[3] = {3, 4};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
            for (unsigned int pop_idx = 0; pop_idx < bd.get_number_of_populations(); ++pop_idx) {
                REQUIRE(bd.get_allele_count(pattern_idx, pop_idx) ==
                        expected_allele_counts.at(pattern_idx).at(pop_idx));
                REQUIRE(bd.get_red_allele_count(pattern_idx, pop_idx) ==
                        expected_red_counts.at(pattern_idx).at(pop_idx));
            }
        }

        std::vector<unsigned int> expected_max_cts = {6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {2, 4};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 2, 0, 3};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        //fold patterns
        unsigned int number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(bd.get_number_of_variable_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.get_path() == nex_path);

        expected_red_counts[2] = {3, 1};
        expected_red_counts[3] = {3, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 4};
        expected_pattern_indices = {0, 1, 2, 0, 3};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
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
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, true, true), EcoevolityBiallelicDataError &);
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
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(bd.get_number_of_variable_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {2,2,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(3);
        expected_allele_counts[0] = {4, 6};
        expected_allele_counts[1] = {4, 4};
        expected_allele_counts[2] = {4, 6};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 3;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(3);
        expected_red_counts[0] = {1, 3};
        expected_red_counts[1] = {1, 1};
        expected_red_counts[2] = {0, 1};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);
    }
}

TEST_CASE("Testing standard diploid with only 0/1 genotypes and charsets", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-only-01.nex") {
        std::string nex_path = "data/diploid-standard-only-01.nex";
        BiallelicData bd(nex_path, ' ', true, true, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(bd.get_number_of_variable_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {2,2,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(3);
        expected_allele_counts[0] = {4, 6};
        expected_allele_counts[1] = {4, 4};
        expected_allele_counts[2] = {4, 6};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 3;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(3);
        expected_red_counts[0] = {1, 3};
        expected_red_counts[1] = {1, 1};
        expected_red_counts[2] = {0, 1};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {1, 4};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 1, 2, 0};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }
}

TEST_CASE("Testing standard haploid with a 2 genotype", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard-012.nex") {
        std::string nex_path = "data/haploid-standard-012.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, false), NxsException &);
    }

    SECTION("Testing data/haploid-standard-012.nex as diploid") {
        std::string nex_path = "data/haploid-standard-012.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, true), NxsException &);
    }

    SECTION("Testing data/haploid-standard-012.nex as dominant") {
        std::string nex_path = "data/haploid-standard-012.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, false, true), NxsException &);
    }

    SECTION("Testing data/haploid-standard-012.nex as diploid and dominant") {
        std::string nex_path = "data/haploid-standard-012.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, true, true), EcoevolityBiallelicDataError &);
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
        BiallelicData bd(nex_path, ' ', true, false);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(bd.get_number_of_variable_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {1,2,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(4);
        expected_allele_counts[0] = {2, 3};
        expected_allele_counts[1] = {2, 2};
        expected_allele_counts[2] = {2, 3};
        expected_allele_counts[3] = {2, 3};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 3;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Folding
        unsigned int number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(bd.get_number_of_variable_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_red_counts[0] = {1, 0};
        expected_red_counts[3] = {1, 1};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }
    }
}

TEST_CASE("Testing standard haploid with charsets", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard.nex") {
        std::string nex_path = "data/haploid-standard.nex";
        BiallelicData bd(nex_path, ' ', true, false, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(bd.get_number_of_variable_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {1,2,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(4);
        expected_allele_counts[0] = {2, 3};
        expected_allele_counts[1] = {2, 2};
        expected_allele_counts[2] = {2, 3};
        expected_allele_counts[3] = {2, 3};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 3;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {3, 4};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 1, 2, 3};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        unsigned int number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(bd.get_number_of_variable_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_red_counts[0] = {1, 0};
        expected_red_counts[3] = {1, 1};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {3, 4};
        expected_pattern_indices = {0, 1, 1, 2, 3};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
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
        BiallelicData bd(nex_path, ' ', true, false, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(bd.get_number_of_variable_sites() == 5);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {1,2,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(4);
        expected_allele_counts[0] = {2, 3};
        expected_allele_counts[1] = {2, 2};
        expected_allele_counts[2] = {2, 3};
        expected_allele_counts[3] = {2, 3};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 3;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Folding
        REQUIRE_THROWS_AS(bd.fold_patterns(), EcoevolityBiallelicDataError &);
    }

    SECTION("Testing data/haploid-standard.nex as diploid and dominant") {
        std::string nex_path = "data/haploid-standard.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, true, true), EcoevolityBiallelicDataError &);
    }

    SECTION("Testing data/haploid-standard.nex as diploid") {
        std::string nex_path = "data/haploid-standard.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, true, false), EcoevolityBiallelicDataError &);
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
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, true, true), EcoevolityBiallelicDataError &);
    }

    SECTION("Testing data/diploid-standard-dominant.nex as dominant and haploid") {
        std::string nex_path = "data/diploid-standard-dominant.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, false, true), EcoevolityBiallelicDataError &);
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
        BiallelicData bd(nex_path, ' ', true, true, false);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(bd.get_number_of_variable_sites() == 7);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {2,2,3};

        std::vector< std::vector<unsigned int> > expected_allele_counts(3);
        expected_allele_counts[0] = {4, 6};
        expected_allele_counts[1] = {4, 4};
        expected_allele_counts[2] = {4, 6};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(3);
        expected_red_counts[0] = {2, 6};
        expected_red_counts[1] = {2, 2};
        expected_red_counts[2] = {0, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Folding
        unsigned int number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(bd.get_number_of_variable_sites() == 7);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_red_counts[0] = {2, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }
    }
}

TEST_CASE("Testing for constant diploid site patterns", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant0.nex") {
        std::string nex_path = "data/diploid-standard-constant0.nex";
        BiallelicData bd(nex_path);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(6);
        expected_allele_counts[0] = {6, 4};
        expected_allele_counts[1] = {6, 2};
        expected_allele_counts[2] = {6, 4};
        expected_allele_counts[3] = {6, 4};
        expected_allele_counts[4] = {6, 4};
        expected_allele_counts[5] = {4, 2};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 4;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Removing patterns
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_expected_wts = {1,1,1,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(4);
        rm_expected_allele_counts[0] = {6, 2};
        rm_expected_allele_counts[1] = {6, 4};
        rm_expected_allele_counts[2] = {6, 4};
        rm_expected_allele_counts[3] = {6, 4};

        std::map<std::vector<unsigned int>, unsigned int> rm_expected_unique_allele_counts;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(0)] = 1;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(1)] = 3;
        REQUIRE(bd.get_unique_allele_counts() == rm_expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);
        
        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        rm_expected_red_counts[1] = {3, 1};
        rm_expected_red_counts[3] = {3, 0};
        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }
    }

    SECTION("Testing data/diploid-standard-constant0.nex with charsets") {
        std::string nex_path = "data/diploid-standard-constant0.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, true, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(6);
        expected_allele_counts[0] = {6, 4};
        expected_allele_counts[1] = {6, 2};
        expected_allele_counts[2] = {6, 4};
        expected_allele_counts[3] = {6, 4};
        expected_allele_counts[4] = {6, 4};
        expected_allele_counts[5] = {4, 2};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 4;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {1, 3, 5};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 2, 3, 4, 5};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Removing patterns
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_expected_wts = {1,1,1,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(4);
        rm_expected_allele_counts[0] = {6, 2};
        rm_expected_allele_counts[1] = {6, 4};
        rm_expected_allele_counts[2] = {6, 4};
        rm_expected_allele_counts[3] = {6, 4};

        std::map<std::vector<unsigned int>, unsigned int> rm_expected_unique_allele_counts;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(0)] = 1;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(1)] = 3;
        REQUIRE(bd.get_unique_allele_counts() == rm_expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {0, 2, 3};
        expected_pattern_indices = {0, 1, 2, 3};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
        
        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        rm_expected_red_counts[1] = {3, 1};
        rm_expected_red_counts[3] = {3, 0};
        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {0, 2, 3};
        expected_pattern_indices = {0, 1, 2, 3};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }

    SECTION("Testing data/diploid-standard-constant2.nex with charsets") {
        std::string nex_path = "data/diploid-standard-constant2.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, true, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {1, 3, 5};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 2, 3, 4, 5};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Removing patterns
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {0, 2, 3};
        expected_pattern_indices = {0, 1, 2, 3};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        rm_expected_red_counts[1] = {3, 1};
        rm_expected_red_counts[3] = {3, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {0, 2, 3};
        expected_pattern_indices = {0, 1, 2, 3};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }


    SECTION("Testing data/diploid-standard-constant2.nex") {
        std::string nex_path = "data/diploid-standard-constant2.nex";
        BiallelicData bd(nex_path);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Removing patterns
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        rm_expected_red_counts[1] = {3, 1};
        rm_expected_red_counts[3] = {3, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }
    }
}


TEST_CASE("Testing for constant haploid site patterns", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard-constant.nex") {
        std::string nex_path = "data/haploid-standard-constant.nex";
        BiallelicData bd(nex_path, ' ', true, false);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {1,1,2,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {3, 2};
        expected_allele_counts[1] = {3, 1};
        expected_allele_counts[2] = {3, 2};
        expected_allele_counts[3] = {3, 2};
        expected_allele_counts[4] = {2, 1};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 4;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(4)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Removing patterns
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_expected_wts = {1,2,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(3);
        rm_expected_allele_counts[0] = {3, 1};
        rm_expected_allele_counts[1] = {3, 2};
        rm_expected_allele_counts[2] = {3, 2};

        std::map<std::vector<unsigned int>, unsigned int> rm_expected_unique_allele_counts;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(0)] = 1;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(1)] = 3;
        REQUIRE(bd.get_unique_allele_counts() == rm_expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(3);
        rm_expected_red_counts[0] = {1, 1};
        rm_expected_red_counts[1] = {2, 2};
        rm_expected_red_counts[2] = {0, 1};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        rm_expected_red_counts[1] = {1, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }
    }

    SECTION("Testing data/haploid-standard-constant.nex with charsets") {
        std::string nex_path = "data/haploid-standard-constant.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, false, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {1,1,2,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {3, 2};
        expected_allele_counts[1] = {3, 1};
        expected_allele_counts[2] = {3, 2};
        expected_allele_counts[3] = {3, 2};
        expected_allele_counts[4] = {2, 1};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 4;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(4)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {1, 3, 5};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 2, 3, 2, 4};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Removing patterns
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_expected_wts = {1,2,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(3);
        rm_expected_allele_counts[0] = {3, 1};
        rm_expected_allele_counts[1] = {3, 2};
        rm_expected_allele_counts[2] = {3, 2};

        std::map<std::vector<unsigned int>, unsigned int> rm_expected_unique_allele_counts;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(0)] = 1;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(1)] = 3;
        REQUIRE(bd.get_unique_allele_counts() == rm_expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(3);
        rm_expected_red_counts[0] = {1, 1};
        rm_expected_red_counts[1] = {2, 2};
        rm_expected_red_counts[2] = {0, 1};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {0, 2, 3};
        expected_pattern_indices = {0, 1, 2, 1};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        rm_expected_red_counts[1] = {1, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {0, 2, 3};
        expected_pattern_indices = {0, 1, 2, 1};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }

    SECTION("Testing data/haploid-standard-constant.nex as dominant") {
        std::string nex_path = "data/haploid-standard-constant.nex";
        BiallelicData bd(nex_path, ' ', true, false, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Removing patterns
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Folding
        REQUIRE_THROWS_AS(bd.fold_patterns(), EcoevolityBiallelicDataError &);
    }

    SECTION("Testing data/haploid-standard-constant.nex as dominant with charsets") {
        std::string nex_path = "data/haploid-standard-constant.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, false, true, true, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {1, 3, 5};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 2, 3, 2, 4};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Removing patterns
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {0, 2, 3};
        expected_pattern_indices = {0, 1, 2, 1};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        REQUIRE_THROWS_AS(bd.fold_patterns(), EcoevolityBiallelicDataError &);
    }

    SECTION("Testing data/haploid-standard-constant.nex as diploid") {
        std::string nex_path = "data/haploid-standard-constant.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, true, false), EcoevolityBiallelicDataError &);
    }

    SECTION("Testing data/haploid-standard-constant.nex as diploid and dominant") {
        std::string nex_path = "data/haploid-standard-constant.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, true, true), EcoevolityBiallelicDataError &);
    }
}



TEST_CASE("Testing for constant dominant diploid site patterns", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant-dominant.nex as dominant") {
        std::string nex_path = "data/diploid-standard-constant-dominant.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, true, true), EcoevolityBiallelicDataError &);
    }
}

TEST_CASE("Testing for missing haploid site patterns", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard-missing.nex") {
        std::string nex_path = "data/haploid-standard-missing.nex";
        BiallelicData bd(nex_path, ' ', true, false);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {2,2,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {1, 0};
        expected_allele_counts[1] = {2, 2};
        expected_allele_counts[2] = {2, 3};
        expected_allele_counts[3] = {0, 0};
        expected_allele_counts[4] = {0, 3};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(4)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Removing patterns
        REQUIRE(bd.get_number_of_missing_sites_removed() == 0);
        unsigned int number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 3);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 4);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 2);
        REQUIRE(bd.get_number_of_sites() == 3);
        REQUIRE(bd.get_number_of_variable_sites() == 3);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_expected_wts = {2,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(2);
        rm_expected_allele_counts[0] = {2, 2};
        rm_expected_allele_counts[1] = {2, 3};

        std::map<std::vector<unsigned int>, unsigned int> rm_expected_unique_allele_counts;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(0)] = 2;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(1)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == rm_expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(2);
        rm_expected_red_counts[0] = {1, 1};
        rm_expected_red_counts[1] = {1, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(2), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(2), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(2), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 2);
        REQUIRE(bd.get_number_of_sites() == 3);
        REQUIRE(bd.get_number_of_variable_sites() == 3);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);

        rm_expected_red_counts[1] = {1, 1};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }
    }

    SECTION("Testing data/haploid-standard-missing.nex with charsets") {
        std::string nex_path = "data/haploid-standard-missing.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, false, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {2,2,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {1, 0};
        expected_allele_counts[1] = {2, 2};
        expected_allele_counts[2] = {2, 3};
        expected_allele_counts[3] = {0, 0};
        expected_allele_counts[4] = {0, 3};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(4)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {2, 6};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 1, 2, 3, 4, 0};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Removing patterns
        REQUIRE(bd.get_number_of_missing_sites_removed() == 0);
        unsigned int number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 3);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 4);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 2);
        REQUIRE(bd.get_number_of_sites() == 3);
        REQUIRE(bd.get_number_of_variable_sites() == 3);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_expected_wts = {2,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(2);
        rm_expected_allele_counts[0] = {2, 2};
        rm_expected_allele_counts[1] = {2, 3};

        std::map<std::vector<unsigned int>, unsigned int> rm_expected_unique_allele_counts;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(0)] = 2;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(1)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == rm_expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(2);
        rm_expected_red_counts[0] = {1, 1};
        rm_expected_red_counts[1] = {1, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(2), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(2), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(2), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {1, 2};
        expected_pattern_indices = {0, 0, 1};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 2);
        REQUIRE(bd.get_number_of_sites() == 3);
        REQUIRE(bd.get_number_of_variable_sites() == 3);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);

        rm_expected_red_counts[1] = {1, 1};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {1, 2};
        expected_pattern_indices = {0, 0, 1};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }
}

TEST_CASE("Testing for missing haploid site patterns as dominant", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard-missing.nex") {
        std::string nex_path = "data/haploid-standard-missing.nex";
        BiallelicData bd(nex_path, ' ', true, false, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {2,2,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {1, 0};
        expected_allele_counts[1] = {2, 2};
        expected_allele_counts[2] = {2, 3};
        expected_allele_counts[3] = {0, 0};
        expected_allele_counts[4] = {0, 3};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(4)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Removing patterns
        REQUIRE(bd.get_number_of_missing_sites_removed() == 0);
        unsigned int number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 3);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 4);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 2);
        REQUIRE(bd.get_number_of_sites() == 3);
        REQUIRE(bd.get_number_of_variable_sites() == 3);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_expected_wts = {2,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(2);
        rm_expected_allele_counts[0] = {2, 2};
        rm_expected_allele_counts[1] = {2, 3};

        std::map<std::vector<unsigned int>, unsigned int> rm_expected_unique_allele_counts;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(0)] = 2;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(1)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == rm_expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(2);
        rm_expected_red_counts[0] = {1, 1};
        rm_expected_red_counts[1] = {1, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(2), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(2), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(2), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Folding
        REQUIRE_THROWS_AS(bd.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}

TEST_CASE("Testing for missing haploid site patterns as dominant with charsets", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard-missing.nex") {
        std::string nex_path = "data/haploid-standard-missing.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, false, true, true, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {2,2,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {1, 0};
        expected_allele_counts[1] = {2, 2};
        expected_allele_counts[2] = {2, 3};
        expected_allele_counts[3] = {0, 0};
        expected_allele_counts[4] = {0, 3};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(4)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {2, 6};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 1, 2, 3, 4, 0};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Removing patterns
        REQUIRE(bd.get_number_of_missing_sites_removed() == 0);
        unsigned int number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 3);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 4);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 2);
        REQUIRE(bd.get_number_of_sites() == 3);
        REQUIRE(bd.get_number_of_variable_sites() == 3);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_expected_wts = {2,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(2);
        rm_expected_allele_counts[0] = {2, 2};
        rm_expected_allele_counts[1] = {2, 3};

        std::map<std::vector<unsigned int>, unsigned int> rm_expected_unique_allele_counts;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(0)] = 2;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(1)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == rm_expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(2);
        rm_expected_red_counts[0] = {1, 1};
        rm_expected_red_counts[1] = {1, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(2), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(2), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(2), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {1, 2};
        expected_pattern_indices = {0, 0, 1};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        REQUIRE_THROWS_AS(bd.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}

TEST_CASE("Testing for constant AND missing haploid site patterns", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard-missing.nex") {
        std::string nex_path = "data/haploid-standard-missing.nex";
        BiallelicData bd(nex_path, ' ', true, false);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {2,2,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {1, 0};
        expected_allele_counts[1] = {2, 2};
        expected_allele_counts[2] = {2, 3};
        expected_allele_counts[3] = {0, 0};
        expected_allele_counts[4] = {0, 3};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(4)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Removing patterns
        REQUIRE(bd.get_number_of_constant_sites_removed() == 0);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 0);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 0);
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_constant_sites_removed() == 3);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 1);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_expected_wts = {2,1,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(3);
        rm_expected_allele_counts[0] = {2, 2};
        rm_expected_allele_counts[1] = {2, 3};
        rm_expected_allele_counts[2] = {0, 3};

        std::map<std::vector<unsigned int>, unsigned int> rm_expected_unique_allele_counts;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(0)] = 2;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(1)] = 1;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(2)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == rm_expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(3);
        rm_expected_red_counts[0] = {1, 1};
        rm_expected_red_counts[1] = {1, 2};
        rm_expected_red_counts[2] = {0, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Removing patterns
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 2);
        REQUIRE(bd.get_number_of_sites() == 3);
        REQUIRE(bd.get_number_of_variable_sites() == 3);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_rm_expected_wts = {2,1};

        std::vector< std::vector<unsigned int> > rm_rm_expected_allele_counts(2);
        rm_rm_expected_allele_counts[0] = {2, 2};
        rm_rm_expected_allele_counts[1] = {2, 3};

        std::map<std::vector<unsigned int>, unsigned int> rm_rm_expected_unique_allele_counts;
        rm_rm_expected_unique_allele_counts[rm_rm_expected_allele_counts.at(0)] = 2;
        rm_rm_expected_unique_allele_counts[rm_rm_expected_allele_counts.at(1)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == rm_rm_expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > rm_rm_expected_red_counts(2);
        rm_rm_expected_red_counts[0] = {1, 1};
        rm_rm_expected_red_counts[1] = {1, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(2), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(2), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(2), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_rm_expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == rm_rm_expected_labels);
        rm_rm_expected_labels.clear();
        rm_rm_expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 2);
        REQUIRE(bd.get_number_of_sites() == 3);
        REQUIRE(bd.get_number_of_variable_sites() == 3);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        rm_rm_expected_red_counts[1] = {1, 1};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_rm_expected_red_counts.at(pattern_idx));
        }
    }
}

TEST_CASE("Testing for constant AND missing haploid site patterns with charsets", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard-missing.nex") {
        std::string nex_path = "data/haploid-standard-missing.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, false, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {2,2,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {1, 0};
        expected_allele_counts[1] = {2, 2};
        expected_allele_counts[2] = {2, 3};
        expected_allele_counts[3] = {0, 0};
        expected_allele_counts[4] = {0, 3};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(4)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {2, 6};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 1, 2, 3, 4, 0};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Removing patterns
        REQUIRE(bd.get_number_of_constant_sites_removed() == 0);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 0);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 0);
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_constant_sites_removed() == 3);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 1);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_expected_wts = {2,1,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(3);
        rm_expected_allele_counts[0] = {2, 2};
        rm_expected_allele_counts[1] = {2, 3};
        rm_expected_allele_counts[2] = {0, 3};

        std::map<std::vector<unsigned int>, unsigned int> rm_expected_unique_allele_counts;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(0)] = 2;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(1)] = 1;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(2)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == rm_expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(3);
        rm_expected_red_counts[0] = {1, 1};
        rm_expected_red_counts[1] = {1, 2};
        rm_expected_red_counts[2] = {0, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {1, 3};
        expected_pattern_indices = {0, 0, 1, 2};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Removing patterns
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 2);
        REQUIRE(bd.get_number_of_sites() == 3);
        REQUIRE(bd.get_number_of_variable_sites() == 3);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_rm_expected_wts = {2,1};

        std::vector< std::vector<unsigned int> > rm_rm_expected_allele_counts(2);
        rm_rm_expected_allele_counts[0] = {2, 2};
        rm_rm_expected_allele_counts[1] = {2, 3};

        std::map<std::vector<unsigned int>, unsigned int> rm_rm_expected_unique_allele_counts;
        rm_rm_expected_unique_allele_counts[rm_rm_expected_allele_counts.at(0)] = 2;
        rm_rm_expected_unique_allele_counts[rm_rm_expected_allele_counts.at(1)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == rm_rm_expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > rm_rm_expected_red_counts(2);
        rm_rm_expected_red_counts[0] = {1, 1};
        rm_rm_expected_red_counts[1] = {1, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(2), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(2), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(2), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_rm_expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == rm_rm_expected_labels);
        rm_rm_expected_labels.clear();
        rm_rm_expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {1, 2};
        expected_pattern_indices = {0, 0, 1};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 2);
        REQUIRE(bd.get_number_of_sites() == 3);
        REQUIRE(bd.get_number_of_variable_sites() == 3);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        rm_rm_expected_red_counts[1] = {1, 1};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {1, 2};
        expected_pattern_indices = {0, 0, 1};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }
}

TEST_CASE("Testing for constant AND missing haploid site patterns as dominant", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard-missing.nex") {
        std::string nex_path = "data/haploid-standard-missing.nex";
        BiallelicData bd(nex_path, ' ', true, false, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Removing patterns
        REQUIRE(bd.get_number_of_constant_sites_removed() == 0);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 0);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 0);
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_constant_sites_removed() == 3);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 1);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Removing patterns
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 2);
        REQUIRE(bd.get_number_of_sites() == 3);
        REQUIRE(bd.get_number_of_variable_sites() == 3);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(2), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(2), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(2), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_rm_expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == rm_rm_expected_labels);
        rm_rm_expected_labels.clear();
        rm_rm_expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Folding
        REQUIRE_THROWS_AS(bd.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}

TEST_CASE("Testing for constant AND missing haploid site patterns as dominant with charsets", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard-missing.nex") {
        std::string nex_path = "data/haploid-standard-missing.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, false, true, true, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {2, 6};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 1, 2, 3, 4, 0};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Removing patterns
        REQUIRE(bd.get_number_of_constant_sites_removed() == 0);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 0);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 0);
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_constant_sites_removed() == 3);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 1);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 4);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {1, 3};
        expected_pattern_indices = {0, 0, 1, 2};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Removing patterns
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 2);
        REQUIRE(bd.get_number_of_sites() == 3);
        REQUIRE(bd.get_number_of_variable_sites() == 3);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(2), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(2), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(2), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_rm_expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == rm_rm_expected_labels);
        rm_rm_expected_labels.clear();
        rm_rm_expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {1, 2};
        expected_pattern_indices = {0, 0, 1};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        REQUIRE_THROWS_AS(bd.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}

TEST_CASE("Testing for constant AND missing dominant haploid site patterns", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard-missing.nex as diploid") {
        std::string nex_path = "data/haploid-standard-missing.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, true, false), EcoevolityBiallelicDataError &);
    }

    SECTION("Testing data/haploid-standard-missing.nex as diploid and dominant") {
        std::string nex_path = "data/haploid-standard-missing.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, true, true), EcoevolityBiallelicDataError &);
    }
}

TEST_CASE("Testing for constant AND missing diploid site patterns", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-missing.nex") {
        std::string nex_path = "data/diploid-standard-missing.nex";
        BiallelicData bd(nex_path);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 8);
        REQUIRE(bd.get_number_of_variable_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {3,1,2,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {6, 0};
        expected_allele_counts[1] = {0, 4};
        expected_allele_counts[2] = {6, 4};
        expected_allele_counts[3] = {6, 2};
        expected_allele_counts[4] = {6, 2};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 3;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Removing patterns
        unsigned int number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 2);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_expected_wts = {2,1,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(3);
        rm_expected_allele_counts[0] = {6, 4};
        rm_expected_allele_counts[1] = {6, 2};
        rm_expected_allele_counts[2] = {6, 2};

        std::map<std::vector<unsigned int>, unsigned int> rm_expected_unique_allele_counts;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(0)] = 2;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(1)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == rm_expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(3);
        rm_expected_red_counts[0] = {3, 3};
        rm_expected_red_counts[1] = {0, 0};
        rm_expected_red_counts[2] = {6, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Removing patterns
        number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 1);
        REQUIRE(bd.get_number_of_sites() == 2);
        REQUIRE(bd.get_number_of_variable_sites() == 2);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_rm_expected_wts = {2};

        std::vector< std::vector<unsigned int> > rm_rm_expected_allele_counts(1);
        rm_rm_expected_allele_counts[0] = {6, 4};

        std::map<std::vector<unsigned int>, unsigned int> rm_rm_expected_unique_allele_counts;
        rm_rm_expected_unique_allele_counts[rm_rm_expected_allele_counts.at(0)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == rm_rm_expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > rm_rm_expected_red_counts(1);
        rm_rm_expected_red_counts[0] = {3, 3};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(1), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(1), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(1), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_rm_expected_labels);
        rm_rm_expected_labels.clear();
        rm_rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 1);
        REQUIRE(bd.get_number_of_sites() == 2);
        REQUIRE(bd.get_number_of_variable_sites() == 2);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        rm_rm_expected_red_counts[0] = {3, 1};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_rm_expected_red_counts.at(pattern_idx));
        }
    }
}

TEST_CASE("Testing for constant AND missing diploid site patterns with charsets", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-missing.nex") {
        std::string nex_path = "data/diploid-standard-missing.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, true, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 8);
        REQUIRE(bd.get_number_of_variable_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {3,1,2,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {6, 0};
        expected_allele_counts[1] = {0, 4};
        expected_allele_counts[2] = {6, 4};
        expected_allele_counts[3] = {6, 2};
        expected_allele_counts[4] = {6, 2};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 3;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {2, 7};
        std::vector<unsigned int> expected_pattern_indices = {0, 0, 1, 2, 3, 4, 2, 0};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Removing patterns
        unsigned int number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.get_number_of_sites() == 4);
        REQUIRE(bd.get_number_of_variable_sites() == 2);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_expected_wts = {2,1,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(3);
        rm_expected_allele_counts[0] = {6, 4};
        rm_expected_allele_counts[1] = {6, 2};
        rm_expected_allele_counts[2] = {6, 2};

        std::map<std::vector<unsigned int>, unsigned int> rm_expected_unique_allele_counts;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(0)] = 2;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(1)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == rm_expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(3);
        rm_expected_red_counts[0] = {3, 3};
        rm_expected_red_counts[1] = {0, 0};
        rm_expected_red_counts[2] = {6, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(3), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(3), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {3};
        expected_pattern_indices = {0, 1, 2, 0};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Removing patterns
        number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 1);
        REQUIRE(bd.get_number_of_sites() == 2);
        REQUIRE(bd.get_number_of_variable_sites() == 2);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_rm_expected_wts = {2};

        std::vector< std::vector<unsigned int> > rm_rm_expected_allele_counts(1);
        rm_rm_expected_allele_counts[0] = {6, 4};

        std::map<std::vector<unsigned int>, unsigned int> rm_rm_expected_unique_allele_counts;
        rm_rm_expected_unique_allele_counts[rm_rm_expected_allele_counts.at(0)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == rm_rm_expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > rm_rm_expected_red_counts(1);
        rm_rm_expected_red_counts[0] = {3, 3};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(1), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(1), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(1), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_rm_expected_labels);
        rm_rm_expected_labels.clear();
        rm_rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {1};
        expected_pattern_indices = {0, 0};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 1);
        REQUIRE(bd.get_number_of_sites() == 2);
        REQUIRE(bd.get_number_of_variable_sites() == 2);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        rm_rm_expected_red_counts[0] = {3, 1};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {1};
        expected_pattern_indices = {0, 0};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }
}

TEST_CASE("Testing for removing all site patterns", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-all-remove.nex") {
        std::string nex_path = "data/diploid-standard-all-remove.nex";
        BiallelicData bd(nex_path);
        int number_removed = bd.remove_constant_patterns();
        REQUIRE_THROWS_AS(bd.remove_missing_population_patterns(), EcoevolityBiallelicDataError &);
    }
}

TEST_CASE("Testing for mirrored diploid site patterns", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-missing.nex") {
        std::string nex_path = "data/diploid-standard-missing.nex";
        BiallelicData bd(nex_path);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 8);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {3,1,2,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {6, 0};
        expected_allele_counts[1] = {0, 4};
        expected_allele_counts[2] = {6, 4};
        expected_allele_counts[3] = {6, 2};
        expected_allele_counts[4] = {6, 2};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 3;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Removing patterns
        unsigned int number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 8);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_expected_wts = {3,1,2,2};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(4);
        rm_expected_allele_counts[0] = {6, 0};
        rm_expected_allele_counts[1] = {0, 4};
        rm_expected_allele_counts[2] = {6, 4};
        rm_expected_allele_counts[3] = {6, 2};

        std::map<std::vector<unsigned int>, unsigned int> rm_expected_unique_allele_counts;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(0)] = 3;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(1)] = 1;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(2)] = 2;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(3)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == rm_expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(4);
        rm_expected_red_counts[0] = {3, 0};
        rm_expected_red_counts[1] = {0, 1};
        rm_expected_red_counts[2] = {3, 1};
        rm_expected_red_counts[3] = {0, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);
    }
}

TEST_CASE("Testing for mirrored diploid site patterns with charsets", "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-missing.nex") {
        std::string nex_path = "data/diploid-standard-missing.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, true, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 8);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {3,1,2,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {6, 0};
        expected_allele_counts[1] = {0, 4};
        expected_allele_counts[2] = {6, 4};
        expected_allele_counts[3] = {6, 2};
        expected_allele_counts[4] = {6, 2};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 3;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {2, 7};
        std::vector<unsigned int> expected_pattern_indices = {0, 0, 1, 2, 3, 4, 2, 0};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Removing patterns
        unsigned int number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 8);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_expected_wts = {3,1,2,2};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(4);
        rm_expected_allele_counts[0] = {6, 0};
        rm_expected_allele_counts[1] = {0, 4};
        rm_expected_allele_counts[2] = {6, 4};
        rm_expected_allele_counts[3] = {6, 2};

        std::map<std::vector<unsigned int>, unsigned int> rm_expected_unique_allele_counts;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(0)] = 3;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(1)] = 1;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(2)] = 2;
        rm_expected_unique_allele_counts[rm_expected_allele_counts.at(3)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == rm_expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(4);
        rm_expected_red_counts[0] = {3, 0};
        rm_expected_red_counts[1] = {0, 1};
        rm_expected_red_counts[2] = {3, 1};
        rm_expected_red_counts[3] = {0, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 7};
        expected_pattern_indices = {0, 0, 1, 2, 3, 3, 2, 0};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }
}

/**
 *
 * #NEXUS
 * Begin data;
 *     Dimensions ntax=5 nchar=6;
 *     Format datatype=standard symbols="012" gap=-;
 *     Matrix
 * pop1_a  ATC-CR
 * pop1_b  GCCCYA
 * pop1_c  ACTGTA
 * pop2_c  ATTGCA
 * pop2_d  GCC?YA
 *     ;
 * End;
 */
TEST_CASE("Testing small, diploid, dna data set", "[BiallelicData]") {

    SECTION("Testing data/diploid-dna.nex") {
        std::string nex_path = "data/diploid-dna.nex";
        BiallelicData bd(nex_path);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {2,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {6, 4};
        expected_allele_counts[1] = {6, 4};
        expected_allele_counts[2] = {4, 2};
        expected_allele_counts[3] = {6, 4};
        expected_allele_counts[4] = {6, 4};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(5);
        expected_red_counts[0] = {2, 2};
        expected_red_counts[1] = {4, 2};
        expected_red_counts[2] = {2, 2};
        expected_red_counts[3] = {3, 1};
        expected_red_counts[4] = {1, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);
    }

    SECTION("Testing data/diploid-dna.nex as haploid") {
        std::string nex_path = "data/diploid-dna.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, false), EcoevolityInvalidCharacterError &);
    }

    SECTION("Testing data/diploid-dna.nex as dominant") {
        std::string nex_path = "data/diploid-dna.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, true, true), EcoevolityBiallelicDataError &);
    }

    SECTION("Testing folding of data/diploid-dna.nex") {
        std::string nex_path = "data/diploid-dna.nex";
        BiallelicData bd(nex_path);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {2,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {6, 4};
        expected_allele_counts[1] = {6, 4};
        expected_allele_counts[2] = {4, 2};
        expected_allele_counts[3] = {6, 4};
        expected_allele_counts[4] = {6, 4};

        std::vector< std::vector<unsigned int> > expected_red_counts(5);
        expected_red_counts[0] = {2, 2};
        expected_red_counts[1] = {4, 2};
        expected_red_counts[2] = {2, 2};
        expected_red_counts[3] = {3, 1};
        expected_red_counts[4] = {1, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Folding
        unsigned int number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_expected_wts = {3,1,1,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(4);
        rm_expected_allele_counts[0] = {6, 4};
        rm_expected_allele_counts[1] = {4, 2};
        rm_expected_allele_counts[2] = {6, 4};
        rm_expected_allele_counts[3] = {6, 4};

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(4);
        rm_expected_red_counts[0] = {2, 2};
        rm_expected_red_counts[1] = {2, 0};
        rm_expected_red_counts[2] = {3, 1};
        rm_expected_red_counts[3] = {1, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);
    }
}

TEST_CASE("Testing small, diploid, dna data set with charsets", "[BiallelicData]") {

    SECTION("Testing data/diploid-dna.nex") {
        std::string nex_path = "data/diploid-dna.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, true, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {2,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {6, 4};
        expected_allele_counts[1] = {6, 4};
        expected_allele_counts[2] = {4, 2};
        expected_allele_counts[3] = {6, 4};
        expected_allele_counts[4] = {6, 4};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(5);
        expected_red_counts[0] = {2, 2};
        expected_red_counts[1] = {4, 2};
        expected_red_counts[2] = {2, 2};
        expected_red_counts[3] = {3, 1};
        expected_red_counts[4] = {1, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {2, 5};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 0, 2, 3, 4};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }

    SECTION("Testing data/diploid-dna.nex as haploid") {
        std::string nex_path = "data/diploid-dna.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, false), EcoevolityInvalidCharacterError &);
    }

    SECTION("Testing data/diploid-dna.nex as dominant") {
        std::string nex_path = "data/diploid-dna.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, true, true), EcoevolityBiallelicDataError &);
    }

    SECTION("Testing folding of data/diploid-dna.nex") {
        std::string nex_path = "data/diploid-dna.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, true, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {2,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(5);
        expected_allele_counts[0] = {6, 4};
        expected_allele_counts[1] = {6, 4};
        expected_allele_counts[2] = {4, 2};
        expected_allele_counts[3] = {6, 4};
        expected_allele_counts[4] = {6, 4};

        std::vector< std::vector<unsigned int> > expected_red_counts(5);
        expected_red_counts[0] = {2, 2};
        expected_red_counts[1] = {4, 2};
        expected_red_counts[2] = {2, 2};
        expected_red_counts[3] = {3, 1};
        expected_red_counts[4] = {1, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {2, 5};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 0, 2, 3, 4};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        unsigned int number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> rm_expected_wts = {3,1,1,1};

        std::vector< std::vector<unsigned int> > rm_expected_allele_counts(4);
        rm_expected_allele_counts[0] = {6, 4};
        rm_expected_allele_counts[1] = {4, 2};
        rm_expected_allele_counts[2] = {6, 4};
        rm_expected_allele_counts[3] = {6, 4};

        std::vector< std::vector<unsigned int> > rm_expected_red_counts(4);
        rm_expected_red_counts[0] = {2, 2};
        rm_expected_red_counts[1] = {2, 0};
        rm_expected_red_counts[2] = {3, 1};
        rm_expected_red_counts[3] = {1, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < rm_expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == rm_expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == rm_expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == rm_expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> rm_expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == rm_expected_labels);
        rm_expected_labels.clear();
        rm_expected_labels = {"pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == rm_expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 5};
        expected_pattern_indices = {0, 0, 0, 1, 2, 3};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }
}

TEST_CASE("Testing diploid dna with missing, mirrored, and constant sites", "[BiallelicData]") {

    SECTION("Testing data/diploid-dna-constant-missing.nex") {
        std::string nex_path = "data/diploid-dna-constant-missing.nex";
        BiallelicData bd(nex_path);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 8);
        REQUIRE(bd.get_number_of_sites() == 9);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {2,1,1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(8);
        expected_allele_counts[0] = {6, 6, 4};
        expected_allele_counts[1] = {6, 6, 4};
        expected_allele_counts[2] = {6, 6, 4};
        expected_allele_counts[3] = {4, 2, 2};
        expected_allele_counts[4] = {6, 6, 4};
        expected_allele_counts[5] = {6, 6, 0};
        expected_allele_counts[6] = {4, 2, 2};
        expected_allele_counts[7] = {6, 6, 4};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 6;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(8);
        expected_red_counts[0] = {2, 4, 2};
        expected_red_counts[1] = {4, 4, 2};
        expected_red_counts[2] = {2, 2, 2};
        expected_red_counts[3] = {2, 2, 0};
        expected_red_counts[4] = {3, 2, 2};
        expected_red_counts[5] = {1, 0, 0};
        expected_red_counts[6] = {0, 0, 0};
        expected_red_counts[7] = {0, 0, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(8), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(8), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(8), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        // Folding
        unsigned int number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 7);
        REQUIRE(bd.get_number_of_sites() == 9);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {2,2,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 0});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 6;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(4)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({3, 2, 2});
        expected_red_counts.push_back({1, 0, 0});
        expected_red_counts.push_back({0, 0, 0});
        expected_red_counts.push_back({0, 0, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(7), std::out_of_range &);
    
        // Remove missing
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 8);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {2,2,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 6;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({3, 2, 2});
        expected_red_counts.push_back({0, 0, 0});
        expected_red_counts.push_back({0, 0, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        // Remove missing
        number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {2,2,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({3, 2, 2});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        expected_labels.clear();
        expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        std::vector<unsigned int> expected_max_cts = {6,6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);
    }
}

TEST_CASE("Testing diploid dna with missing, mirrored, and constant sites with charsets", "[BiallelicData]") {

    SECTION("Testing data/diploid-dna-constant-missing.nex") {
        std::string nex_path = "data/diploid-dna-constant-missing.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, true, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 8);
        REQUIRE(bd.get_number_of_sites() == 9);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {2,1,1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(8);
        expected_allele_counts[0] = {6, 6, 4};
        expected_allele_counts[1] = {6, 6, 4};
        expected_allele_counts[2] = {6, 6, 4};
        expected_allele_counts[3] = {4, 2, 2};
        expected_allele_counts[4] = {6, 6, 4};
        expected_allele_counts[5] = {6, 6, 0};
        expected_allele_counts[6] = {4, 2, 2};
        expected_allele_counts[7] = {6, 6, 4};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 6;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(8);
        expected_red_counts[0] = {2, 4, 2};
        expected_red_counts[1] = {4, 4, 2};
        expected_red_counts[2] = {2, 2, 2};
        expected_red_counts[3] = {2, 2, 0};
        expected_red_counts[4] = {3, 2, 2};
        expected_red_counts[5] = {1, 0, 0};
        expected_red_counts[6] = {0, 0, 0};
        expected_red_counts[7] = {0, 0, 0};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(8), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(8), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(8), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {2, 5, 8};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 2, 3, 4, 5, 6, 7, 0};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        unsigned int number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 7);
        REQUIRE(bd.get_number_of_sites() == 9);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {2,2,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 0});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 6;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(4)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({3, 2, 2});
        expected_red_counts.push_back({1, 0, 0});
        expected_red_counts.push_back({0, 0, 0});
        expected_red_counts.push_back({0, 0, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(7), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 5, 8};
        expected_pattern_indices = {0, 1, 1, 2, 3, 4, 5, 6, 0};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    
        // Remove missing
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 8);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {2,2,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 6;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({3, 2, 2});
        expected_red_counts.push_back({0, 0, 0});
        expected_red_counts.push_back({0, 0, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 4, 7};
        expected_pattern_indices = {0, 1, 1, 2, 3, 4, 5, 0};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Remove constant
        number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {2,2,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({3, 2, 2});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        expected_labels.clear();
        expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        std::vector<unsigned int> expected_max_cts = {6,6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 4, 5};
        expected_pattern_indices = {0, 1, 1, 2, 3, 0};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }
}

TEST_CASE("Testing diploid dna with missing, mirrored, constant sites, and no hets", "[BiallelicData]") {

    SECTION("Testing data/diploid-dna-constant-missing-nohets.nex") {
        std::string nex_path = "data/diploid-dna-constant-missing-nohets.nex";
        BiallelicData bd(nex_path);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 9);
        REQUIRE(bd.get_number_of_sites() == 9);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        REQUIRE(bd.get_number_of_constant_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 0);

        std::vector<unsigned int> expected_wts = {1,1,1,1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(9);
        expected_allele_counts[0] = {6, 6, 4};
        expected_allele_counts[1] = {6, 6, 4};
        expected_allele_counts[2] = {6, 6, 4};
        expected_allele_counts[3] = {4, 2, 2};
        expected_allele_counts[4] = {6, 6, 4};
        expected_allele_counts[5] = {6, 6, 0};
        expected_allele_counts[6] = {4, 2, 2};
        expected_allele_counts[7] = {6, 6, 4};
        expected_allele_counts[8] = {6, 6, 4};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 6;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(9);
        expected_red_counts[0] = {2, 4, 2};
        expected_red_counts[1] = {4, 4, 2};
        expected_red_counts[2] = {2, 2, 2};
        expected_red_counts[3] = {2, 2, 0};
        expected_red_counts[4] = {2, 0, 0};
        expected_red_counts[5] = {4, 6, 0};
        expected_red_counts[6] = {0, 0, 0};
        expected_red_counts[7] = {0, 0, 0};
        expected_red_counts[8] = {4, 2, 4};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        std::vector<unsigned int> expected_max_cts = {6,6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE_THROWS_AS(bd.get_pattern_weight(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(9), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        // Remove constant
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 7);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        REQUIRE(bd.get_number_of_constant_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 0);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 0});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({4, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({4, 6, 0});
        expected_red_counts.push_back({4, 2, 4});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(7), std::out_of_range &);
    
        // Remove missing
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        REQUIRE(bd.get_number_of_constant_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 1);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({4, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({4, 2, 4});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        REQUIRE(bd.get_number_of_constant_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 1);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({2, 4, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        expected_labels.clear();
        expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);
    }

    SECTION("Testing data/diploid-dna-constant-missing-nohets.nex as dominant") {
        std::string nex_path = "data/diploid-dna-constant-missing-nohets.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, true, true), EcoevolityBiallelicDataError &);
    }

    SECTION("Testing data/diploid-dna-constant-missing-nohets.nex as dominant and haploid") {
        std::string nex_path = "data/diploid-dna-constant-missing-nohets.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, false, true), EcoevolityBiallelicDataError &);
    }

    SECTION("Testing data/diploid-dna-constant-missing-nohets.nex as haploid") {
        std::string nex_path = "data/diploid-dna-constant-missing-nohets.nex";
        BiallelicData bd(nex_path, ' ', true, false);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 9);
        REQUIRE(bd.get_number_of_sites() == 9);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {1,1,1,1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(9);
        expected_allele_counts[0] = {3, 3, 2};
        expected_allele_counts[1] = {3, 3, 2};
        expected_allele_counts[2] = {3, 3, 2};
        expected_allele_counts[3] = {2, 1, 1};
        expected_allele_counts[4] = {3, 3, 2};
        expected_allele_counts[5] = {3, 3, 0};
        expected_allele_counts[6] = {2, 1, 1};
        expected_allele_counts[7] = {3, 3, 2};
        expected_allele_counts[8] = {3, 3, 2};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 6;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(9);
        expected_red_counts[0] = {1, 2, 1};
        expected_red_counts[1] = {2, 2, 1};
        expected_red_counts[2] = {1, 1, 1};
        expected_red_counts[3] = {1, 1, 0};
        expected_red_counts[4] = {1, 0, 0};
        expected_red_counts[5] = {2, 3, 0};
        expected_red_counts[6] = {0, 0, 0};
        expected_red_counts[7] = {0, 0, 0};
        expected_red_counts[8] = {2, 1, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(9), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        // Remove constant
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 7);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({2, 1, 1});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 0});
        expected_allele_counts.push_back({3, 3, 2});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({1, 2, 1});
        expected_red_counts.push_back({2, 2, 1});
        expected_red_counts.push_back({1, 1, 1});
        expected_red_counts.push_back({1, 1, 0});
        expected_red_counts.push_back({1, 0, 0});
        expected_red_counts.push_back({2, 3, 0});
        expected_red_counts.push_back({2, 1, 2});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(7), std::out_of_range &);
    
        // Remove missing
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({2, 1, 1});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({1, 2, 1});
        expected_red_counts.push_back({2, 2, 1});
        expected_red_counts.push_back({1, 1, 1});
        expected_red_counts.push_back({1, 1, 0});
        expected_red_counts.push_back({1, 0, 0});
        expected_red_counts.push_back({2, 1, 2});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({2, 1, 1});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({1, 2, 1});
        expected_red_counts.push_back({1, 1, 1});
        expected_red_counts.push_back({1, 1, 0});
        expected_red_counts.push_back({1, 0, 0});
        expected_red_counts.push_back({1, 2, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        expected_labels.clear();
        expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        std::vector<unsigned int> expected_max_cts = {3,3,2};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);
    }
}

TEST_CASE("Testing diploid dna with missing, mirrored, constant sites, no hets, and charsets", "[BiallelicData]") {

    SECTION("Testing data/diploid-dna-constant-missing-nohets.nex") {
        std::string nex_path = "data/diploid-dna-constant-missing-nohets.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, true, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 9);
        REQUIRE(bd.get_number_of_sites() == 9);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        REQUIRE(bd.get_number_of_constant_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 0);

        std::vector<unsigned int> expected_wts = {1,1,1,1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(9);
        expected_allele_counts[0] = {6, 6, 4};
        expected_allele_counts[1] = {6, 6, 4};
        expected_allele_counts[2] = {6, 6, 4};
        expected_allele_counts[3] = {4, 2, 2};
        expected_allele_counts[4] = {6, 6, 4};
        expected_allele_counts[5] = {6, 6, 0};
        expected_allele_counts[6] = {4, 2, 2};
        expected_allele_counts[7] = {6, 6, 4};
        expected_allele_counts[8] = {6, 6, 4};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 6;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(9);
        expected_red_counts[0] = {2, 4, 2};
        expected_red_counts[1] = {4, 4, 2};
        expected_red_counts[2] = {2, 2, 2};
        expected_red_counts[3] = {2, 2, 0};
        expected_red_counts[4] = {2, 0, 0};
        expected_red_counts[5] = {4, 6, 0};
        expected_red_counts[6] = {0, 0, 0};
        expected_red_counts[7] = {0, 0, 0};
        expected_red_counts[8] = {4, 2, 4};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        std::vector<unsigned int> expected_max_cts = {6,6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE_THROWS_AS(bd.get_pattern_weight(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(9), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {2, 5, 8};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Remove constant
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 7);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        REQUIRE(bd.get_number_of_constant_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 0);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 0});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({4, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({4, 6, 0});
        expected_red_counts.push_back({4, 2, 4});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(7), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 5, 6};
        expected_pattern_indices = {0, 1, 2, 3, 4, 5, 6};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    
        // Remove missing
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        REQUIRE(bd.get_number_of_constant_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 1);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({4, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({4, 2, 4});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 4, 5};
        expected_pattern_indices = {0, 1, 2, 3, 4, 5};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        REQUIRE(bd.get_number_of_constant_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 1);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({2, 4, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        expected_labels.clear();
        expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 4, 5};
        expected_pattern_indices = {0, 1, 1, 2, 3, 4};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }

    SECTION("Testing data/diploid-dna-constant-missing-nohets.nex as haploid") {
        std::string nex_path = "data/diploid-dna-constant-missing-nohets.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, false, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 9);
        REQUIRE(bd.get_number_of_sites() == 9);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {1,1,1,1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(9);
        expected_allele_counts[0] = {3, 3, 2};
        expected_allele_counts[1] = {3, 3, 2};
        expected_allele_counts[2] = {3, 3, 2};
        expected_allele_counts[3] = {2, 1, 1};
        expected_allele_counts[4] = {3, 3, 2};
        expected_allele_counts[5] = {3, 3, 0};
        expected_allele_counts[6] = {2, 1, 1};
        expected_allele_counts[7] = {3, 3, 2};
        expected_allele_counts[8] = {3, 3, 2};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 6;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(9);
        expected_red_counts[0] = {1, 2, 1};
        expected_red_counts[1] = {2, 2, 1};
        expected_red_counts[2] = {1, 1, 1};
        expected_red_counts[3] = {1, 1, 0};
        expected_red_counts[4] = {1, 0, 0};
        expected_red_counts[5] = {2, 3, 0};
        expected_red_counts[6] = {0, 0, 0};
        expected_red_counts[7] = {0, 0, 0};
        expected_red_counts[8] = {2, 1, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(9), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {2, 5, 8};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Remove constant
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 7);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({2, 1, 1});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 0});
        expected_allele_counts.push_back({3, 3, 2});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({1, 2, 1});
        expected_red_counts.push_back({2, 2, 1});
        expected_red_counts.push_back({1, 1, 1});
        expected_red_counts.push_back({1, 1, 0});
        expected_red_counts.push_back({1, 0, 0});
        expected_red_counts.push_back({2, 3, 0});
        expected_red_counts.push_back({2, 1, 2});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(7), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 5, 6};
        expected_pattern_indices = {0, 1, 2, 3, 4, 5, 6};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    
        // Remove missing
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({2, 1, 1});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({1, 2, 1});
        expected_red_counts.push_back({2, 2, 1});
        expected_red_counts.push_back({1, 1, 1});
        expected_red_counts.push_back({1, 1, 0});
        expected_red_counts.push_back({1, 0, 0});
        expected_red_counts.push_back({2, 1, 2});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 4, 5};
        expected_pattern_indices = {0, 1, 2, 3, 4, 5};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({2, 1, 1});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({1, 2, 1});
        expected_red_counts.push_back({1, 1, 1});
        expected_red_counts.push_back({1, 1, 0});
        expected_red_counts.push_back({1, 0, 0});
        expected_red_counts.push_back({1, 2, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        expected_labels.clear();
        expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        std::vector<unsigned int> expected_max_cts = {3,3,2};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 4, 5};
        expected_pattern_indices = {0, 1, 1, 2, 3, 4};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }
}

/*
 * #NEXUS
 * Begin data;
 *     Dimensions ntax=8 nchar=9;
 *     Format datatype=dna gap=-;
 *     Matrix
 * pop1_a  ATC-CGGAT
 * pop1_b  ---C-AGA-
 * pop1_c  ACTGTA?AC
 * pop2_d  ATCGCAGAC
 * pop2_e  GCT?CA-AT
 * pop2_f  GCC?CA-AT
 * pop3_g  GCT?C?GAC
 * pop3_h  ATCCC??AC
 *     ;
 * End;
 */
TEST_CASE("Testing change in max sample size", "[BiallelicData]") {
    SECTION("Testing data/diploid-dna-sample-size-changes-nohets.nex") {
        std::string nex_path = "data/diploid-dna-sample-size-changes-nohets.nex";
        BiallelicData bd(nex_path);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 9);
        REQUIRE(bd.get_number_of_sites() == 9);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {1,1,1,1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(9);
        expected_allele_counts[0] = {4, 6, 4};
        expected_allele_counts[1] = {4, 6, 4};
        expected_allele_counts[2] = {4, 6, 4};
        expected_allele_counts[3] = {4, 2, 2};
        expected_allele_counts[4] = {4, 6, 4};
        expected_allele_counts[5] = {6, 6, 0};
        expected_allele_counts[6] = {4, 2, 2};
        expected_allele_counts[7] = {6, 6, 4};
        expected_allele_counts[8] = {4, 6, 4};

        std::vector< std::vector<unsigned int> > expected_red_counts(9);
        expected_red_counts[0] = {0, 4, 2};
        expected_red_counts[1] = {2, 4, 2};
        expected_red_counts[2] = {2, 2, 2};
        expected_red_counts[3] = {2, 2, 0};
        expected_red_counts[4] = {2, 0, 0};
        expected_red_counts[5] = {4, 6, 0};
        expected_red_counts[6] = {0, 0, 0};
        expected_red_counts[7] = {0, 0, 0};
        expected_red_counts[8] = {2, 2, 4};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        std::vector<unsigned int> expected_max_cts = {6,6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE_THROWS_AS(bd.get_pattern_weight(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(9), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        // Remove constant
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 7);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({6, 6, 0});
        expected_allele_counts.push_back({4, 6, 4});

        expected_red_counts.clear();
        expected_red_counts.push_back({0, 4, 2});
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({4, 6, 0});
        expected_red_counts.push_back({2, 2, 4});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        expected_max_cts = {6,6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE_THROWS_AS(bd.get_pattern_weight(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(7), std::out_of_range &);
    
        // Remove missing
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 6, 4});

        expected_red_counts.clear();
        expected_red_counts.push_back({0, 4, 2});
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({2, 2, 4});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        expected_max_cts = {4,6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 6, 4});

        expected_red_counts.clear();
        expected_red_counts.push_back({0, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({2, 4, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
            for (unsigned int pop_idx = 0; pop_idx < bd.get_number_of_populations(); ++pop_idx) {
                REQUIRE(bd.get_allele_count(pattern_idx, pop_idx) ==
                        expected_allele_counts.at(pattern_idx).at(pop_idx));
                REQUIRE(bd.get_red_allele_count(pattern_idx, pop_idx) ==
                        expected_red_counts.at(pattern_idx).at(pop_idx));
            }
        }

        expected_max_cts = {4,6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        expected_labels.clear();
        expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);
    }
}

TEST_CASE("Testing change in max sample size with charsets", "[BiallelicData]") {
    SECTION("Testing data/diploid-dna-sample-size-changes-nohets.nex") {
        std::string nex_path = "data/diploid-dna-sample-size-changes-nohets.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, true, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 9);
        REQUIRE(bd.get_number_of_sites() == 9);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        std::vector<unsigned int> expected_wts = {1,1,1,1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(9);
        expected_allele_counts[0] = {4, 6, 4};
        expected_allele_counts[1] = {4, 6, 4};
        expected_allele_counts[2] = {4, 6, 4};
        expected_allele_counts[3] = {4, 2, 2};
        expected_allele_counts[4] = {4, 6, 4};
        expected_allele_counts[5] = {6, 6, 0};
        expected_allele_counts[6] = {4, 2, 2};
        expected_allele_counts[7] = {6, 6, 4};
        expected_allele_counts[8] = {4, 6, 4};

        std::vector< std::vector<unsigned int> > expected_red_counts(9);
        expected_red_counts[0] = {0, 4, 2};
        expected_red_counts[1] = {2, 4, 2};
        expected_red_counts[2] = {2, 2, 2};
        expected_red_counts[3] = {2, 2, 0};
        expected_red_counts[4] = {2, 0, 0};
        expected_red_counts[5] = {4, 6, 0};
        expected_red_counts[6] = {0, 0, 0};
        expected_red_counts[7] = {0, 0, 0};
        expected_red_counts[8] = {2, 2, 4};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        std::vector<unsigned int> expected_max_cts = {6,6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE_THROWS_AS(bd.get_pattern_weight(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(9), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {2, 5, 8};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Remove constant
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 7);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({6, 6, 0});
        expected_allele_counts.push_back({4, 6, 4});

        expected_red_counts.clear();
        expected_red_counts.push_back({0, 4, 2});
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({4, 6, 0});
        expected_red_counts.push_back({2, 2, 4});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        expected_max_cts = {6,6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE_THROWS_AS(bd.get_pattern_weight(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(7), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 5, 6};
        expected_pattern_indices = {0, 1, 2, 3, 4, 5, 6};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    
        // Remove missing
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 6, 4});

        expected_red_counts.clear();
        expected_red_counts.push_back({0, 4, 2});
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({2, 2, 4});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        expected_max_cts = {4,6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 4, 5};
        expected_pattern_indices = {0, 1, 2, 3, 4, 5};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({4, 6, 4});
        expected_allele_counts.push_back({4, 6, 4});

        expected_red_counts.clear();
        expected_red_counts.push_back({0, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({2, 4, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
            for (unsigned int pop_idx = 0; pop_idx < bd.get_number_of_populations(); ++pop_idx) {
                REQUIRE(bd.get_allele_count(pattern_idx, pop_idx) ==
                        expected_allele_counts.at(pattern_idx).at(pop_idx));
                REQUIRE(bd.get_red_allele_count(pattern_idx, pop_idx) ==
                        expected_red_counts.at(pattern_idx).at(pop_idx));
            }
        }

        expected_max_cts = {4,6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        expected_labels.clear();
        expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 4, 5};
        expected_pattern_indices = {0, 1, 1, 2, 3, 4};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }
}


TEST_CASE("Testing empirical u/v rates", "[BiallelicData]") {
    SECTION("Testing data/diploid-standard-data-ntax5-nchar5.nex") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        BiallelicData bd(nex_path);
        bd.remove_constant_patterns();
        bd.remove_missing_population_patterns();
        double u;
        double v;
        REQUIRE(bd.get_proportion_of_red_alleles() == Approx(5.0/12.0));
        bd.get_empirical_u_v_rates(u, v);
        // Verified these values in SNAPP
        REQUIRE(u == Approx(1.2));
        REQUIRE(v == Approx(0.8571429));
    }

    SECTION("Testing data/haploid-standard.nex") {
        std::string nex_path = "data/haploid-standard.nex";
        BiallelicData bd(nex_path, ' ', true, false);
        bd.remove_constant_patterns();
        bd.remove_missing_population_patterns();
        double u;
        double v;
        REQUIRE(bd.get_proportion_of_red_alleles() == Approx(12.0/23.0));
        bd.get_empirical_u_v_rates(u, v);
        // Verified these values in SNAPP
        REQUIRE(u == Approx(0.9583333));
        REQUIRE(v == Approx(1.045455));
    }

    SECTION("Testing data/haploid-standard.nex as dominant") {
        std::string nex_path = "data/haploid-standard.nex";
        BiallelicData bd(nex_path, ' ', true, false, true);
        bd.remove_constant_patterns();
        bd.remove_missing_population_patterns();
        double u;
        double v;
        REQUIRE(bd.get_proportion_of_red_alleles() == Approx(12.0/23.0));
        bd.get_empirical_u_v_rates(u, v);
        // Verified these values in SNAPP
        REQUIRE(u == Approx(0.9583333));
        REQUIRE(v == Approx(1.045455));
    }

    SECTION("Testing data/haploid-standard-constant.nex") {
        std::string nex_path = "data/haploid-standard-constant.nex";
        BiallelicData bd(nex_path, ' ', true, false);
        REQUIRE(bd.get_proportion_of_red_alleles() == Approx(14.0/27.0));
        double u;
        double v;
        // SNAPP: u = 0.90625 v = 1.1153846153846154
        // SNAPP seems to be getting the proportion of zeros (green) by doing:
        //     # of zeros / all cells except ?
        // the denominator is the count of all cells in the alignement except
        // "?"; this includes "-".
        bd.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.9642857142857143));
        REQUIRE(v == Approx(1.0384615384615383));

        bd.remove_constant_patterns();
        bd.remove_missing_population_patterns();
        REQUIRE(bd.get_proportion_of_red_alleles() == Approx(11.0/19.0));
        bd.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8636364));
        REQUIRE(v == Approx(1.1875));
    }

    SECTION("Testing data/haploid-standard-missing.nex") {
        std::string nex_path = "data/haploid-standard-missing.nex";
        BiallelicData bd(nex_path, ' ', true, false);
        REQUIRE(bd.get_proportion_of_red_alleles() == Approx(0.6111111111111112));
        double u;
        double v;
        // SNAPP: u = 0.8235294117647058 v = 1.2727272727272727
        // SNAPP seems to be getting the proportion of zeros (green) by doing:
        //     # of zeros / all cells except ?
        // the denominator is the count of all cells in the alignement except
        // "?"; this includes "-".
        bd.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8181818181818181));
        REQUIRE(v == Approx(1.2857142857142858));

        bd.remove_constant_patterns();
        bd.remove_missing_population_patterns();
        REQUIRE(bd.get_proportion_of_red_alleles() == Approx(7.0/13.0));
        bd.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.9285714));
        REQUIRE(v == Approx(1.083333));
    }

    SECTION("Testing data/diploid-dna-constant-missing.nex") {
        std::string nex_path = "data/diploid-dna-constant-missing.nex";
        BiallelicData bd(nex_path);
        REQUIRE(bd.get_proportion_of_red_alleles() == Approx(0.3548387096774194));
        double u;
        double v;
        bd.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(1.409090909090909));
        REQUIRE(v == Approx(0.775));

        bd.remove_constant_patterns();
        bd.remove_missing_population_patterns();
        REQUIRE(bd.get_proportion_of_red_alleles() == Approx(0.48863636363636365));
        bd.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(1.0232558139534884));
        REQUIRE(v == Approx(0.9777777777777777));
    }
}

TEST_CASE("Testing quoted underscores", "[BiallelicData]") {

    SECTION("Testing quoted underscores") {
        std::string nex_path = "data/haploid-standard-quoted-underscores.nex";
        BiallelicData bd(nex_path, '_', true, false);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1_a", "pop1_b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2_c", "pop2_d", "pop2_e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Folding
        unsigned int number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);

        expected_red_counts[0] = {1, 0};
        expected_red_counts[3] = {1, 1};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }
    }
}

TEST_CASE("Testing quoted spaces", "[BiallelicData]") {

    SECTION("Testing quoted spaces") {
        std::string nex_path = "data/haploid-standard-quoted-spaces.nex";
        BiallelicData bd(nex_path, ' ', true, false);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);

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

        REQUIRE_THROWS_AS(bd.get_pattern_weight(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(4), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(4), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE_THROWS_AS(bd.get_population_label(2), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 c", "pop2 d", "pop2 e"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(2), std::out_of_range &);

        // Folding
        unsigned int number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);

        expected_red_counts[0] = {1, 0};
        expected_red_counts[3] = {1, 1};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }
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
TEST_CASE("Testing writing methods for diploid standard data set",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-data-ntax5-nchar5.nex") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                true,  // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.get_path() == nex_path);

        char delim = ' ';
        std::stringstream expected_alignment;
        expected_alignment << "\'pop1" << delim << "0000\'  00000\n"
                           << "\'pop1" << delim << "0001\'  00000\n"
                           << "\'pop1" << delim << "0002\'  00000\n"
                           << "\'pop1" << delim << "0003\'  00011\n"
                           << "\'pop1" << delim << "0004\'  00011\n"
                           << "\'pop1" << delim << "0005\'  00111\n"
                           << "\'pop2" << delim << "0000\'  00101\n"
                           << "\'pop2" << delim << "0001\'  00111\n"
                           << "\'pop2" << delim << "0002\'  11?11\n"
                           << "\'pop2" << delim << "0003\'  11?11\n";
        std::stringstream alignment;
        bd.write_alignment(alignment, delim);

        REQUIRE(alignment.str() == expected_alignment.str());

        std::stringstream expected_nexus;
        expected_nexus << "#NEXUS\n"
                       << "Begin data;\n"
                       << "    Dimensions ntax=10 nchar=5;\n"
                       << "    Format datatype=standard symbols=\"01\" missing=? gap=-;\n"
                       << "    Matrix\n"
                       << expected_alignment.str()
                       << "    ;\n"
                       << "End;\n";

        std::stringstream nexus;
        bd.write_nexus(nexus, delim);

        REQUIRE(nexus.str() == expected_nexus.str());

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t25.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
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
TEST_CASE("Testing writing methods for standard diploid with only 0/1 genotypes",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-only-01.nex") {
        std::string nex_path = "data/diploid-standard-only-01.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                true,  // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 3);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);

        char delim = ' ';
        std::stringstream expected_alignment;
        expected_alignment << "\'pop1" << delim << "0000\'  00000\n"
                           << "\'pop1" << delim << "0001\'  00000\n"
                           << "\'pop1" << delim << "0002\'  00000\n"
                           << "\'pop1" << delim << "0003\'  11110\n"
                           << "\'pop2" << delim << "0000\'  00000\n"
                           << "\'pop2" << delim << "0001\'  00000\n"
                           << "\'pop2" << delim << "0002\'  00000\n"
                           << "\'pop2" << delim << "0003\'  11110\n"
                           << "\'pop2" << delim << "0004\'  11??0\n"
                           << "\'pop2" << delim << "0005\'  11??1\n";
        std::stringstream alignment;
        bd.write_alignment(alignment, delim);

        REQUIRE(alignment.str() == expected_alignment.str());

        std::stringstream expected_nexus;
        expected_nexus << "#NEXUS\n"
                       << "Begin data;\n"
                       << "    Dimensions ntax=10 nchar=5;\n"
                       << "    Format datatype=standard symbols=\"01\" missing=? gap=-;\n"
                       << "    Matrix\n"
                       << expected_alignment.str()
                       << "    ;\n"
                       << "End;\n";

        std::stringstream nexus;
        bd.write_nexus(nexus, delim);

        REQUIRE(nexus.str() == expected_nexus.str());

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t26.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
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
TEST_CASE("Testing writing methods for standard haploid", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard.nex") {
        std::string nex_path = "data/haploid-standard.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);

        char delim = ' ';
        std::stringstream expected_alignment;
        expected_alignment << "\'pop1" << delim << "0000\'  00000\n"
                           << "\'pop1" << delim << "0001\'  11101\n"
                           << "\'pop2" << delim << "0000\'  10000\n"
                           << "\'pop2" << delim << "0001\'  11101\n"
                           << "\'pop2" << delim << "0002\'  1??11\n";

        std::stringstream alignment;
        bd.write_alignment(alignment, delim);

        REQUIRE(alignment.str() == expected_alignment.str());

        std::stringstream expected_nexus;
        expected_nexus << "#NEXUS\n"
                       << "Begin data;\n"
                       << "    Dimensions ntax=5 nchar=5;\n"
                       << "    Format datatype=standard symbols=\"01\" missing=? gap=-;\n"
                       << "    Matrix\n"
                       << expected_alignment.str()
                       << "    ;\n"
                       << "End;\n";

        std::stringstream nexus;
        bd.write_nexus(nexus, delim);

        REQUIRE(nexus.str() == expected_nexus.str());

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t27.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
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
TEST_CASE("Testing writing methods for standard haploid dominant", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard.nex as dominant") {
        std::string nex_path = "data/haploid-standard.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                true,  // markers are dominant (false)
                true); // validate (true)
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);

        char delim = ' ';
        std::stringstream expected_alignment;
        expected_alignment << "\'pop1" << delim << "0000\'  00000\n"
                           << "\'pop1" << delim << "0001\'  11101\n"
                           << "\'pop2" << delim << "0000\'  10000\n"
                           << "\'pop2" << delim << "0001\'  11101\n"
                           << "\'pop2" << delim << "0002\'  1??11\n";

        std::stringstream alignment;
        bd.write_alignment(alignment, delim);

        REQUIRE(alignment.str() == expected_alignment.str());

        std::stringstream expected_nexus;
        expected_nexus << "#NEXUS\n"
                       << "Begin data;\n"
                       << "    Dimensions ntax=5 nchar=5;\n"
                       << "    Format datatype=standard symbols=\"01\" missing=? gap=-;\n"
                       << "    Matrix\n"
                       << expected_alignment.str()
                       << "    ;\n"
                       << "End;\n";

        std::stringstream nexus;
        bd.write_nexus(nexus, delim);

        REQUIRE(nexus.str() == expected_nexus.str());

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t28.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                true, // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing writing methods for constant diploid site patterns",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant0.nex") {
        std::string nex_path = "data/diploid-standard-constant0.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                true,  // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);

        char delim = ' ';
        std::stringstream expected_alignment;
        expected_alignment << "\'pop1" << delim << "0000\'  000000\n"
                           << "\'pop1" << delim << "0001\'  000000\n"
                           << "\'pop1" << delim << "0002\'  000000\n"
                           << "\'pop1" << delim << "0003\'  001010\n"
                           << "\'pop1" << delim << "0004\'  00101?\n"
                           << "\'pop1" << delim << "0005\'  01101?\n"
                           << "\'pop2" << delim << "0000\'  010010\n"
                           << "\'pop2" << delim << "0001\'  011010\n"
                           << "\'pop2" << delim << "0002\'  0?111?\n"
                           << "\'pop2" << delim << "0003\'  0?111?\n";

        std::stringstream alignment;
        bd.write_alignment(alignment, delim);

        REQUIRE(alignment.str() == expected_alignment.str());

        std::stringstream expected_nexus;
        expected_nexus << "#NEXUS\n"
                       << "Begin data;\n"
                       << "    Dimensions ntax=10 nchar=6;\n"
                       << "    Format datatype=standard symbols=\"01\" missing=? gap=-;\n"
                       << "    Matrix\n"
                       << expected_alignment.str()
                       << "    ;\n"
                       << "End;\n";

        std::stringstream nexus;
        bd.write_nexus(nexus, delim);

        REQUIRE(nexus.str() == expected_nexus.str());

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t29.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }

        bd.fold_patterns();

        tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        test_path = "data/tmp-data-" + tag + "-t30.cfg";
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d2(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d2.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d2.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d2.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d2.get_number_of_sites());

        for (unsigned int i = 0; i < d2.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d2.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d2.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d2.get_pattern_weight(i));
        }

        bd.remove_constant_patterns();

        tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        test_path = "data/tmp-data-" + tag + "-t31.cfg";
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d3(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d3.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d3.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d3.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d3.get_number_of_sites());

        for (unsigned int i = 0; i < d3.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d3.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d3.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d3.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing writing methods for aflp data",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant0.nex") {
        std::string nex_path = "data/aflp_25-with-constant.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                false,  // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        char delim = ' ';
        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t32.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing writing methods for aflp data with folding",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant0.nex") {
        std::string nex_path = "data/aflp_25-with-constant.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        bd.fold_patterns();

        char delim = ' ';
        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t33.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing writing methods for aflp data with removing",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant0.nex") {
        std::string nex_path = "data/aflp_25-with-constant.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        bd.remove_constant_patterns();

        char delim = ' ';
        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t34.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing writing methods for aflp data with folding and removing",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant0.nex") {
        std::string nex_path = "data/aflp_25-with-constant.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        bd.fold_patterns();
        bd.remove_constant_patterns();

        char delim = ' ';
        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t35.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing writing methods for dominant aflp data",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant0.nex") {
        std::string nex_path = "data/aflp_25-with-constant.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                false,  // genotypes are diploid (true)
                true,  // markers are dominant (false)
                true); // validate (true)

        char delim = ' ';
        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t36.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                true,  // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing writing methods for dominant aflp data with removing",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant0.nex") {
        std::string nex_path = "data/aflp_25-with-constant.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                true,  // markers are dominant (false)
                true); // validate (true)

        bd.remove_constant_patterns();

        char delim = ' ';
        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t37.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                true,  // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing writing methods for hemi data",
        "[BiallelicData]") {

    SECTION("Testing data/hemi129.nex") {
        std::string nex_path = "data/hemi129.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                true,  // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        char delim = ' ';
        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t38.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing writing methods for hemi data with folding",
        "[BiallelicData]") {

    SECTION("Testing data/hemi129.nex") {
        std::string nex_path = "data/hemi129.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                true,  // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        bd.fold_patterns();

        char delim = ' ';
        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t39.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing writing methods for hemi data with removing",
        "[BiallelicData]") {

    SECTION("Testing data/hemi129.nex") {
        std::string nex_path = "data/hemi129.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                true,  // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        bd.remove_constant_patterns();

        char delim = ' ';
        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t40.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing writing methods for hemi data with folding and removing",
        "[BiallelicData]") {

    SECTION("Testing data/hemi129.nex") {
        std::string nex_path = "data/hemi129.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                true,  // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        bd.fold_patterns();
        bd.remove_constant_patterns();

        char delim = ' ';
        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t41.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing writing methods for quoted underscores",
        "[BiallelicData]") {

    SECTION("Testing quoted underscores") {
        std::string nex_path = "data/haploid-standard-quoted-underscores.nex";
        BiallelicData bd(
                nex_path,
                '_',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);

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

        char delim = '_';
        std::stringstream expected_alignment;
        expected_alignment << "\'pop1" << delim << "0000\'  00000\n"
                           << "\'pop1" << delim << "0001\'  11101\n"
                           << "\'pop2" << delim << "0000\'  10000\n"
                           << "\'pop2" << delim << "0001\'  11101\n"
                           << "\'pop2" << delim << "0002\'  1??11\n";
        std::stringstream alignment;
        bd.write_alignment(alignment, delim);

        REQUIRE(alignment.str() == expected_alignment.str());

        std::stringstream expected_nexus;
        expected_nexus << "#NEXUS\n"
                       << "Begin data;\n"
                       << "    Dimensions ntax=5 nchar=5;\n"
                       << "    Format datatype=standard symbols=\"01\" missing=? gap=-;\n"
                       << "    Matrix\n"
                       << expected_alignment.str()
                       << "    ;\n"
                       << "End;\n";

        std::stringstream nexus;
        bd.write_nexus(nexus, delim);

        REQUIRE(nexus.str() == expected_nexus.str());

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t42.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d(
                test_path,
                '_',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing writing methods for quoted spaces",
        "[BiallelicData]") {

    SECTION("Testing quoted spaces") {
        std::string nex_path = "data/haploid-standard-quoted-spaces.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)
        REQUIRE(bd.get_number_of_populations() == 2);
        REQUIRE(bd.get_number_of_patterns() == 4);
        REQUIRE(bd.get_number_of_sites() == 5);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == false);

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

        char delim = ' ';
        std::stringstream expected_alignment;
        expected_alignment << "\'pop1" << delim << "0000\'  00000\n"
                           << "\'pop1" << delim << "0001\'  11101\n"
                           << "\'pop2" << delim << "0000\'  10000\n"
                           << "\'pop2" << delim << "0001\'  11101\n"
                           << "\'pop2" << delim << "0002\'  1??11\n";
        std::stringstream alignment;
        bd.write_alignment(alignment, delim);

        REQUIRE(alignment.str() == expected_alignment.str());

        std::stringstream expected_nexus;
        expected_nexus << "#NEXUS\n"
                       << "Begin data;\n"
                       << "    Dimensions ntax=5 nchar=5;\n"
                       << "    Format datatype=standard symbols=\"01\" missing=? gap=-;\n"
                       << "    Matrix\n"
                       << expected_alignment.str()
                       << "    ;\n"
                       << "End;\n";

        std::stringstream nexus;
        bd.write_nexus(nexus, delim);

        REQUIRE(nexus.str() == expected_nexus.str());

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t43.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_nexus(os, delim);
        os.close();

        BiallelicData d(
                test_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

//
// Test alignment:
//
// #NEXUS
// Begin data;
// 	Dimensions ntax=8 nchar=10;
// 	Format datatype=dna gap=-;
// 	Matrix
// pop1_a  ATC-CRGAYA
// pop1_b  GCCCYAGAYG
// pop1_c  ACTGTA?ACA
// pop2_d  ATCGCAGATR
// pop2_e  GCT?YA-ACC
// pop2_f  GCC?YA-ARA
// pop3_g  GCT?Y?GAAT
// pop3_h  ATCCY??AGY
// 	;
// End;
//
// Recoded as biallelic:
//
// pop1_a  000-010010
// pop1_b  2200100012
// pop1_c  022220?000
// pop2_d  0002000021
// pop2_e  222?10-002
// pop2_f  220?10-020
// pop3_g  222?1?0022
// pop3_h  00001??022
//
// Distilled to allele counts | mirror allele counts:
//
// 2/6 4/6 2/4  | 2 4 2 
// 4/6 4/6 2/4  | 2 2 2*
// 2/6 2/6 2/4  | 2 2 2
// 2/4 2/2 0/2  | 2 2 0
// 3/6 2/6 2/4  | 3 2 4
// 1/6 0/6 0/0  | 1 0 0
// 0/4 0/2 0/2  | 0 0 0
// 0/6 0/6 0/4  | 0 0 0
// 2/6 4/6 4/4  | 4 2 0*
// 2/6 3/6 4/4  | 4 3 0*

TEST_CASE("Testing diploid dna with triallelic, missing, mirrored, and constant sites", "[BiallelicData]") {

    SECTION("Testing data/diploid-dna-constant-missing-triallelic.nex") {
        std::string nex_path = "data/diploid-dna-constant-missing-triallelic.nex";
        BiallelicData bd(nex_path);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 10);
        REQUIRE(bd.get_number_of_sites() == 10);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        std::vector<unsigned int> expected_wts = {1,1,1,1,1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(10);
        expected_allele_counts[0] = {6, 6, 4};
        expected_allele_counts[1] = {6, 6, 4};
        expected_allele_counts[2] = {6, 6, 4};
        expected_allele_counts[3] = {4, 2, 2};
        expected_allele_counts[4] = {6, 6, 4};
        expected_allele_counts[5] = {6, 6, 0};
        expected_allele_counts[6] = {4, 2, 2};
        expected_allele_counts[7] = {6, 6, 4};
        expected_allele_counts[8] = {6, 6, 4};
        expected_allele_counts[9] = {6, 6, 4};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 7;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(10);
        expected_red_counts[0] = {2, 4, 2};
        expected_red_counts[1] = {4, 4, 2};
        expected_red_counts[2] = {2, 2, 2};
        expected_red_counts[3] = {2, 2, 0};
        expected_red_counts[4] = {3, 2, 2};
        expected_red_counts[5] = {1, 0, 0};
        expected_red_counts[6] = {0, 0, 0};
        expected_red_counts[7] = {0, 0, 0};
        expected_red_counts[8] = {2, 4, 4};
        expected_red_counts[9] = {2, 3, 4};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(10), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(10), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(10), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        // Folding
        unsigned int number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 9);
        REQUIRE(bd.get_number_of_sites() == 10);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 0});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 7;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(4)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({3, 2, 2});
        expected_red_counts.push_back({1, 0, 0});
        expected_red_counts.push_back({0, 0, 0});
        expected_red_counts.push_back({0, 0, 0});
        expected_red_counts.push_back({4, 2, 0});
        expected_red_counts.push_back({4, 3, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(9), std::out_of_range &);
    
        // Remove missing
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 8);
        REQUIRE(bd.get_number_of_sites() == 9);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 7;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({3, 2, 2});
        expected_red_counts.push_back({0, 0, 0});
        expected_red_counts.push_back({0, 0, 0});
        expected_red_counts.push_back({4, 2, 0});
        expected_red_counts.push_back({4, 3, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(8), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(8), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(8), std::out_of_range &);

        // Remove constant
        number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 6;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({3, 2, 2});
        expected_red_counts.push_back({4, 2, 0});
        expected_red_counts.push_back({4, 3, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        expected_labels.clear();
        expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        std::vector<unsigned int> expected_max_cts = {6,6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);
    }
}

TEST_CASE("Testing diploid dna with triallelic, missing, mirrored, constant sites, and charsets", "[BiallelicData]") {

    SECTION("Testing data/diploid-dna-constant-missing-triallelic.nex") {
        std::string nex_path = "data/diploid-dna-constant-missing-triallelic.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, true, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 10);
        REQUIRE(bd.get_number_of_sites() == 10);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        std::vector<unsigned int> expected_wts = {1,1,1,1,1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(10);
        expected_allele_counts[0] = {6, 6, 4};
        expected_allele_counts[1] = {6, 6, 4};
        expected_allele_counts[2] = {6, 6, 4};
        expected_allele_counts[3] = {4, 2, 2};
        expected_allele_counts[4] = {6, 6, 4};
        expected_allele_counts[5] = {6, 6, 0};
        expected_allele_counts[6] = {4, 2, 2};
        expected_allele_counts[7] = {6, 6, 4};
        expected_allele_counts[8] = {6, 6, 4};
        expected_allele_counts[9] = {6, 6, 4};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 7;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(10);
        expected_red_counts[0] = {2, 4, 2};
        expected_red_counts[1] = {4, 4, 2};
        expected_red_counts[2] = {2, 2, 2};
        expected_red_counts[3] = {2, 2, 0};
        expected_red_counts[4] = {3, 2, 2};
        expected_red_counts[5] = {1, 0, 0};
        expected_red_counts[6] = {0, 0, 0};
        expected_red_counts[7] = {0, 0, 0};
        expected_red_counts[8] = {2, 4, 4};
        expected_red_counts[9] = {2, 3, 4};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(10), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(10), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(10), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {2, 5, 8, 9};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        unsigned int number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 9);
        REQUIRE(bd.get_number_of_sites() == 10);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 0});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 7;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(4)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({3, 2, 2});
        expected_red_counts.push_back({1, 0, 0});
        expected_red_counts.push_back({0, 0, 0});
        expected_red_counts.push_back({0, 0, 0});
        expected_red_counts.push_back({4, 2, 0});
        expected_red_counts.push_back({4, 3, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(9), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 5, 8, 9};
        expected_pattern_indices = {0, 1, 1, 2, 3, 4, 5, 6, 7, 8};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    
        // Remove missing
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 8);
        REQUIRE(bd.get_number_of_sites() == 9);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 7;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({3, 2, 2});
        expected_red_counts.push_back({0, 0, 0});
        expected_red_counts.push_back({0, 0, 0});
        expected_red_counts.push_back({4, 2, 0});
        expected_red_counts.push_back({4, 3, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(8), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(8), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(8), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 4, 7, 8};
        expected_pattern_indices = {0, 1, 1, 2, 3, 4, 5, 6, 7};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Remove constant
        number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 6;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({3, 2, 2});
        expected_red_counts.push_back({4, 2, 0});
        expected_red_counts.push_back({4, 3, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        expected_labels.clear();
        expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        std::vector<unsigned int> expected_max_cts = {6,6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 4, 5, 6};
        expected_pattern_indices = {0, 1, 1, 2, 3, 4, 5};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }
}

TEST_CASE("Testing diploid dna with missing, mirrored, constant, triallelic sites, and no hets",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-dna-constant-missing-triallelic-nohets.nex") {
        std::string nex_path = "data/diploid-dna-constant-missing-triallelic-nohets.nex";
        BiallelicData bd(nex_path);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 9);
        REQUIRE(bd.get_number_of_sites() == 9);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        REQUIRE(bd.get_number_of_constant_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 0);

        std::vector<unsigned int> expected_wts = {1,1,1,1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(9);
        expected_allele_counts[0] = {6, 6, 4};
        expected_allele_counts[1] = {6, 6, 4};
        expected_allele_counts[2] = {6, 6, 4};
        expected_allele_counts[3] = {4, 2, 2};
        expected_allele_counts[4] = {6, 6, 4};
        expected_allele_counts[5] = {6, 6, 0};
        expected_allele_counts[6] = {4, 2, 2};
        expected_allele_counts[7] = {6, 6, 4};
        expected_allele_counts[8] = {6, 6, 4};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 6;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(9);
        expected_red_counts[0] = {2, 4, 2};
        expected_red_counts[1] = {4, 4, 2};
        expected_red_counts[2] = {2, 2, 2};
        expected_red_counts[3] = {2, 2, 0};
        expected_red_counts[4] = {2, 0, 0};
        expected_red_counts[5] = {4, 6, 0};
        expected_red_counts[6] = {0, 0, 0};
        expected_red_counts[7] = {0, 0, 0};
        expected_red_counts[8] = {4, 2, 4};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        std::vector<unsigned int> expected_max_cts = {6,6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE_THROWS_AS(bd.get_pattern_weight(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(9), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        // Remove constant
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 7);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        REQUIRE(bd.get_number_of_constant_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 0);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 0});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({4, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({4, 6, 0});
        expected_red_counts.push_back({4, 2, 4});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(7), std::out_of_range &);
    
        // Remove missing
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        REQUIRE(bd.get_number_of_constant_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 1);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({4, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({4, 2, 4});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        REQUIRE(bd.get_number_of_constant_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 1);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({2, 4, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        expected_labels.clear();
        expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);
    }

    SECTION("Testing data/diploid-dna-constant-missing-triallelic-nohets.nex as dominant") {
        std::string nex_path = "data/diploid-dna-constant-missing-triallelic-nohets.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, true, true), EcoevolityBiallelicDataError &);
    }

    SECTION("Testing data/diploid-dna-constant-missing-triallelic-nohets.nex as dominant and haploid") {
        std::string nex_path = "data/diploid-dna-constant-missing-triallelic-nohets.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, false, true), EcoevolityBiallelicDataError &);
    }

    SECTION("Testing data/diploid-dna-constant-missing-triallelic-nohets.nex as haploid") {
        std::string nex_path = "data/diploid-dna-constant-missing-triallelic-nohets.nex";
        BiallelicData bd(nex_path, ' ', true, false);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 9);
        REQUIRE(bd.get_number_of_sites() == 9);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        std::vector<unsigned int> expected_wts = {1,1,1,1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(9);
        expected_allele_counts[0] = {3, 3, 2};
        expected_allele_counts[1] = {3, 3, 2};
        expected_allele_counts[2] = {3, 3, 2};
        expected_allele_counts[3] = {2, 1, 1};
        expected_allele_counts[4] = {3, 3, 2};
        expected_allele_counts[5] = {3, 3, 0};
        expected_allele_counts[6] = {2, 1, 1};
        expected_allele_counts[7] = {3, 3, 2};
        expected_allele_counts[8] = {3, 3, 2};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 6;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(9);
        expected_red_counts[0] = {1, 2, 1};
        expected_red_counts[1] = {2, 2, 1};
        expected_red_counts[2] = {1, 1, 1};
        expected_red_counts[3] = {1, 1, 0};
        expected_red_counts[4] = {1, 0, 0};
        expected_red_counts[5] = {2, 3, 0};
        expected_red_counts[6] = {0, 0, 0};
        expected_red_counts[7] = {0, 0, 0};
        expected_red_counts[8] = {2, 1, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(9), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        // Remove constant
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 7);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({2, 1, 1});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 0});
        expected_allele_counts.push_back({3, 3, 2});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({1, 2, 1});
        expected_red_counts.push_back({2, 2, 1});
        expected_red_counts.push_back({1, 1, 1});
        expected_red_counts.push_back({1, 1, 0});
        expected_red_counts.push_back({1, 0, 0});
        expected_red_counts.push_back({2, 3, 0});
        expected_red_counts.push_back({2, 1, 2});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(7), std::out_of_range &);
    
        // Remove missing
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({2, 1, 1});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({1, 2, 1});
        expected_red_counts.push_back({2, 2, 1});
        expected_red_counts.push_back({1, 1, 1});
        expected_red_counts.push_back({1, 1, 0});
        expected_red_counts.push_back({1, 0, 0});
        expected_red_counts.push_back({2, 1, 2});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({2, 1, 1});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({1, 2, 1});
        expected_red_counts.push_back({1, 1, 1});
        expected_red_counts.push_back({1, 1, 0});
        expected_red_counts.push_back({1, 0, 0});
        expected_red_counts.push_back({1, 2, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        expected_labels.clear();
        expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        std::vector<unsigned int> expected_max_cts = {3,3,2};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);
    }
}

TEST_CASE("Testing diploid dna with missing, mirrored, constant, triallelic sites, no hets, and charsets",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-dna-constant-missing-triallelic-nohets.nex") {
        std::string nex_path = "data/diploid-dna-constant-missing-triallelic-nohets.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, true, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 9);
        REQUIRE(bd.get_number_of_sites() == 9);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        REQUIRE(bd.get_number_of_constant_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 0);

        std::vector<unsigned int> expected_wts = {1,1,1,1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(9);
        expected_allele_counts[0] = {6, 6, 4};
        expected_allele_counts[1] = {6, 6, 4};
        expected_allele_counts[2] = {6, 6, 4};
        expected_allele_counts[3] = {4, 2, 2};
        expected_allele_counts[4] = {6, 6, 4};
        expected_allele_counts[5] = {6, 6, 0};
        expected_allele_counts[6] = {4, 2, 2};
        expected_allele_counts[7] = {6, 6, 4};
        expected_allele_counts[8] = {6, 6, 4};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 6;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(9);
        expected_red_counts[0] = {2, 4, 2};
        expected_red_counts[1] = {4, 4, 2};
        expected_red_counts[2] = {2, 2, 2};
        expected_red_counts[3] = {2, 2, 0};
        expected_red_counts[4] = {2, 0, 0};
        expected_red_counts[5] = {4, 6, 0};
        expected_red_counts[6] = {0, 0, 0};
        expected_red_counts[7] = {0, 0, 0};
        expected_red_counts[8] = {4, 2, 4};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        std::vector<unsigned int> expected_max_cts = {6,6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE_THROWS_AS(bd.get_pattern_weight(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(9), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {1, 3, 5, 7, 8};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Remove constant
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 7);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        REQUIRE(bd.get_number_of_constant_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 0);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 0});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({4, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({4, 6, 0});
        expected_red_counts.push_back({4, 2, 4});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(7), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {1, 3, 5, 6};
        expected_pattern_indices = {0, 1, 2, 3, 4, 5, 6};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    
        // Remove missing
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        REQUIRE(bd.get_number_of_constant_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 1);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({4, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({4, 2, 4});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {1, 3, 4, 5};
        expected_pattern_indices = {0, 1, 2, 3, 4, 5};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        REQUIRE(bd.get_number_of_constant_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_green_sites_removed() == 2);
        REQUIRE(bd.get_number_of_constant_red_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 1);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({2, 0, 0});
        expected_red_counts.push_back({2, 4, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        expected_labels.clear();
        expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {1, 3, 4, 5};
        expected_pattern_indices = {0, 1, 1, 2, 3, 4};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }

    SECTION("Testing data/diploid-dna-constant-missing-triallelic-nohets.nex as dominant") {
        std::string nex_path = "data/diploid-dna-constant-missing-triallelic-nohets.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, true, true), EcoevolityBiallelicDataError &);
    }

    SECTION("Testing data/diploid-dna-constant-missing-triallelic-nohets.nex as dominant and haploid") {
        std::string nex_path = "data/diploid-dna-constant-missing-triallelic-nohets.nex";
        REQUIRE_THROWS_AS(BiallelicData bd(nex_path, ' ', true, false, true), EcoevolityBiallelicDataError &);
    }

    SECTION("Testing data/diploid-dna-constant-missing-triallelic-nohets.nex as haploid") {
        std::string nex_path = "data/diploid-dna-constant-missing-triallelic-nohets.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, false, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 9);
        REQUIRE(bd.get_number_of_sites() == 9);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        std::vector<unsigned int> expected_wts = {1,1,1,1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(9);
        expected_allele_counts[0] = {3, 3, 2};
        expected_allele_counts[1] = {3, 3, 2};
        expected_allele_counts[2] = {3, 3, 2};
        expected_allele_counts[3] = {2, 1, 1};
        expected_allele_counts[4] = {3, 3, 2};
        expected_allele_counts[5] = {3, 3, 0};
        expected_allele_counts[6] = {2, 1, 1};
        expected_allele_counts[7] = {3, 3, 2};
        expected_allele_counts[8] = {3, 3, 2};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 6;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(9);
        expected_red_counts[0] = {1, 2, 1};
        expected_red_counts[1] = {2, 2, 1};
        expected_red_counts[2] = {1, 1, 1};
        expected_red_counts[3] = {1, 1, 0};
        expected_red_counts[4] = {1, 0, 0};
        expected_red_counts[5] = {2, 3, 0};
        expected_red_counts[6] = {0, 0, 0};
        expected_red_counts[7] = {0, 0, 0};
        expected_red_counts[8] = {2, 1, 2};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(9), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {1, 3, 5, 7, 8};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Remove constant
        unsigned int number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 7);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({2, 1, 1});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 0});
        expected_allele_counts.push_back({3, 3, 2});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({1, 2, 1});
        expected_red_counts.push_back({2, 2, 1});
        expected_red_counts.push_back({1, 1, 1});
        expected_red_counts.push_back({1, 1, 0});
        expected_red_counts.push_back({1, 0, 0});
        expected_red_counts.push_back({2, 3, 0});
        expected_red_counts.push_back({2, 1, 2});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(7), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(7), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {2, 5, 6};
        expected_locus_ends = {1, 3, 5, 6};
        expected_pattern_indices = {0, 1, 2, 3, 4, 5, 6};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    
        // Remove missing
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        expected_wts.clear();
        expected_wts = {1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({2, 1, 1});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({1, 2, 1});
        expected_red_counts.push_back({2, 2, 1});
        expected_red_counts.push_back({1, 1, 1});
        expected_red_counts.push_back({1, 1, 0});
        expected_red_counts.push_back({1, 0, 0});
        expected_red_counts.push_back({2, 1, 2});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {1, 3, 4, 5};
        expected_pattern_indices = {0, 1, 2, 3, 4, 5};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 5);
        REQUIRE(bd.get_number_of_sites() == 6);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({2, 1, 1});
        expected_allele_counts.push_back({3, 3, 2});
        expected_allele_counts.push_back({3, 3, 2});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({1, 2, 1});
        expected_red_counts.push_back({1, 1, 1});
        expected_red_counts.push_back({1, 1, 0});
        expected_red_counts.push_back({1, 0, 0});
        expected_red_counts.push_back({1, 2, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(5), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(5), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        expected_labels.clear();
        expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        std::vector<unsigned int> expected_max_cts = {3,3,2};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {1, 3, 4, 5};
        expected_pattern_indices = {0, 1, 1, 2, 3, 4};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }
}

TEST_CASE("Testing charsets", "[BiallelicData]") {

    SECTION("Testing data/diploid-dna-constant-missing-triallelic.nex") {
        std::string nex_path = "data/diploid-dna-constant-missing-triallelic-cs10.nex";
        // file, delim, pop_is_prefix, diploid, dominant, validate, seq_loci
        BiallelicData bd(nex_path, ' ', true, true, false, true, true);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 10);
        REQUIRE(bd.get_number_of_sites() == 10);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == true);
        REQUIRE(bd.patterns_are_folded() == false);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        std::vector<unsigned int> expected_wts = {1,1,1,1,1,1,1,1,1,1};

        std::vector< std::vector<unsigned int> > expected_allele_counts(10);
        expected_allele_counts[0] = {6, 6, 4};
        expected_allele_counts[1] = {6, 6, 4};
        expected_allele_counts[2] = {6, 6, 4};
        expected_allele_counts[3] = {4, 2, 2};
        expected_allele_counts[4] = {6, 6, 4};
        expected_allele_counts[5] = {6, 6, 0};
        expected_allele_counts[6] = {4, 2, 2};
        expected_allele_counts[7] = {6, 6, 4};
        expected_allele_counts[8] = {6, 6, 4};
        expected_allele_counts[9] = {6, 6, 4};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 7;
        expected_unique_allele_counts[expected_allele_counts.at(3)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(5)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(10);
        expected_red_counts[0] = {2, 4, 2};
        expected_red_counts[1] = {4, 4, 2};
        expected_red_counts[2] = {2, 2, 2};
        expected_red_counts[3] = {2, 2, 0};
        expected_red_counts[4] = {3, 2, 2};
        expected_red_counts[5] = {1, 0, 0};
        expected_red_counts[6] = {0, 0, 0};
        expected_red_counts[7] = {0, 0, 0};
        expected_red_counts[8] = {2, 4, 4};
        expected_red_counts[9] = {2, 3, 4};

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(10), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(10), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(10), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        std::vector<std::string> expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        std::vector<unsigned int> expected_pattern_indices = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Folding
        unsigned int number_removed = bd.fold_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 9);
        REQUIRE(bd.get_number_of_sites() == 10);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 0});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 7;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 2;
        expected_unique_allele_counts[expected_allele_counts.at(4)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({3, 2, 2});
        expected_red_counts.push_back({1, 0, 0});
        expected_red_counts.push_back({0, 0, 0});
        expected_red_counts.push_back({0, 0, 0});
        expected_red_counts.push_back({4, 2, 0});
        expected_red_counts.push_back({4, 3, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(9), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(9), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        expected_pattern_indices = {0, 1, 1, 2, 3, 4, 5, 6, 7, 8};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    
        // Remove missing
        number_removed = bd.remove_missing_population_patterns();
        REQUIRE(number_removed == 1);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 8);
        REQUIRE(bd.get_number_of_sites() == 9);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 7;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 2;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({3, 2, 2});
        expected_red_counts.push_back({0, 0, 0});
        expected_red_counts.push_back({0, 0, 0});
        expected_red_counts.push_back({4, 2, 0});
        expected_red_counts.push_back({4, 3, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(8), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(8), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(8), std::out_of_range &);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        expected_pattern_indices = {0, 1, 1, 2, 3, 4, 5, 6, 7};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);

        // Remove constant
        number_removed = bd.remove_constant_patterns();
        REQUIRE(number_removed == 2);
        REQUIRE(bd.get_number_of_populations() == 3);
        REQUIRE(bd.get_number_of_patterns() == 6);
        REQUIRE(bd.get_number_of_sites() == 7);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.has_recoded_triallelic_sites() == true);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 2);

        expected_wts.clear();
        expected_wts = {1,2,1,1,1,1};

        expected_allele_counts.clear();
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({4, 2, 2});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});
        expected_allele_counts.push_back({6, 6, 4});

        expected_unique_allele_counts.clear();
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 6;
        expected_unique_allele_counts[expected_allele_counts.at(2)] = 1;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        expected_red_counts.clear();
        expected_red_counts.push_back({2, 4, 2});
        expected_red_counts.push_back({2, 2, 2});
        expected_red_counts.push_back({2, 2, 0});
        expected_red_counts.push_back({3, 2, 2});
        expected_red_counts.push_back({4, 2, 0});
        expected_red_counts.push_back({4, 3, 0});

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
        }

        REQUIRE_THROWS_AS(bd.get_pattern_weight(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_allele_counts(6), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_red_allele_counts(6), std::out_of_range &);

        REQUIRE(bd.get_population_index("pop1") == 0);
        REQUIRE(bd.get_population_index("pop2") == 1);
        REQUIRE(bd.get_population_index("pop3") == 2);
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE(bd.get_population_label(0) == "pop1");
        REQUIRE(bd.get_population_label(1) == "pop2");
        REQUIRE(bd.get_population_label(2) == "pop3");
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);

        expected_labels.clear();
        expected_labels = {"pop1 a", "pop1 b", "pop1 c"};
        REQUIRE(bd.get_sequence_labels(0) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop2 d", "pop2 e", "pop2 f"};
        REQUIRE(bd.get_sequence_labels(1) == expected_labels);
        expected_labels.clear();
        expected_labels = {"pop3 g", "pop3 h"};
        REQUIRE(bd.get_sequence_labels(2) == expected_labels);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(3), std::out_of_range &);

        std::vector<unsigned int> expected_max_cts = {6,6,4};
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        REQUIRE(bd.has_seq_loci_info() == true);
        expected_locus_ends = {0, 1, 2, 3, 4, 5, 6};
        expected_pattern_indices = {0, 1, 1, 2, 3, 4, 5};
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
    }
}

TEST_CASE("Testing pop label constructor with 2 pops, 4 genomes, locus length 1", "[BiallelicData]") {

    SECTION("Testing label constructor") {
        unsigned int npops = 2;
        unsigned int ngenomes = 4;
        unsigned int nloci = 10;
        unsigned int locus_length = 1;
        unsigned int nsites = nloci * locus_length;
        std::vector<std::string> pop_labels;
        std::vector< std::vector<std::string> > expected_seq_labels;
        for (unsigned int i = 0; i < npops; ++i) {
            std::ostringstream p_label;
            p_label << "pop" << i;
            pop_labels.push_back(p_label.str());
            std::vector<std::string> seq_labels;
            for (unsigned int j = 0; j < ngenomes; ++j) {
                std::ostringstream s_label;
                s_label << p_label.str() << "-genome" << j;
                seq_labels.push_back(s_label.str());
            }
            expected_seq_labels.push_back(seq_labels);
        }
        BiallelicData bd(pop_labels, ngenomes, nloci, locus_length, true);
        REQUIRE(bd.get_number_of_populations() == npops);
        REQUIRE(bd.get_number_of_patterns() == 1);
        REQUIRE(bd.get_number_of_sites() == nsites);
        REQUIRE(bd.get_number_of_variable_sites() == 0);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.get_path() == "");
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends(nloci);
        unsigned int locus_end = locus_length - 1;
        expected_locus_ends.at(0) = locus_end;
        for (unsigned int i = 1; i < nloci; ++i) {
            locus_end += locus_length;
            expected_locus_ends.at(i) = locus_end;
        }
        std::vector<unsigned int> expected_pattern_indices(nsites, 0);
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_contiguous_pattern_indices().size() == nsites);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
        REQUIRE(bd.get_locus_end_indices().size() == nloci);
        REQUIRE(bd.get_locus_end_indices().back() == nsites - 1);

        for (unsigned int i = 0; i < npops; ++i) {
            REQUIRE(bd.get_population_index(pop_labels.at(i)) == i);
            REQUIRE(bd.get_population_label(i) == pop_labels.at(i));
            REQUIRE(bd.get_sequence_labels(i) == expected_seq_labels.at(i));
        }
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_population_label(npops), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(npops), std::out_of_range &);

        std::vector<unsigned int> expected_wts (1, nsites);

        std::vector< std::vector<unsigned int> > expected_allele_counts(1);
        expected_allele_counts[0] = std::vector<unsigned int>(npops, ngenomes);
        std::vector<unsigned int> expected_max_cts = expected_allele_counts.at(0);
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = nsites;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(1);
        expected_red_counts[0] = std::vector<unsigned int>(npops, 0);

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
            for (unsigned int pop_idx = 0; pop_idx < npops; ++pop_idx) {
                REQUIRE(bd.get_allele_count(pattern_idx, pop_idx) ==
                        expected_allele_counts.at(pattern_idx).at(pop_idx));
                REQUIRE(bd.get_red_allele_count(pattern_idx, pop_idx) ==
                        expected_red_counts.at(pattern_idx).at(pop_idx));
            }
        }
    }
}

TEST_CASE("Testing pop label constructor with 3 pops, 6 genomes, locus length 2", "[BiallelicData]") {

    SECTION("Testing label constructor") {
        unsigned int npops = 3;
        unsigned int ngenomes = 6;
        unsigned int nloci = 10;
        unsigned int locus_length = 2;
        unsigned int nsites = nloci * locus_length;
        std::vector<std::string> pop_labels;
        std::vector< std::vector<std::string> > expected_seq_labels;
        for (unsigned int i = 0; i < npops; ++i) {
            std::ostringstream p_label;
            p_label << "pop" << i;
            pop_labels.push_back(p_label.str());
            std::vector<std::string> seq_labels;
            for (unsigned int j = 0; j < ngenomes; ++j) {
                std::ostringstream s_label;
                s_label << p_label.str() << "-genome" << j;
                seq_labels.push_back(s_label.str());
            }
            expected_seq_labels.push_back(seq_labels);
        }
        BiallelicData bd(pop_labels, ngenomes, nloci, locus_length, true);
        REQUIRE(bd.get_number_of_populations() == npops);
        REQUIRE(bd.get_number_of_patterns() == 1);
        REQUIRE(bd.get_number_of_sites() == nsites);
        REQUIRE(bd.get_number_of_variable_sites() == 0);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.get_path() == "");
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends(nloci);
        unsigned int locus_end = locus_length - 1;
        expected_locus_ends.at(0) = locus_end;
        for (unsigned int i = 1; i < nloci; ++i) {
            locus_end += locus_length;
            expected_locus_ends.at(i) = locus_end;
        }
        std::vector<unsigned int> expected_pattern_indices(nsites, 0);
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_contiguous_pattern_indices().size() == nsites);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
        REQUIRE(bd.get_locus_end_indices().size() == nloci);
        REQUIRE(bd.get_locus_end_indices().back() == nsites - 1);

        for (unsigned int i = 0; i < npops; ++i) {
            REQUIRE(bd.get_population_index(pop_labels.at(i)) == i);
            REQUIRE(bd.get_population_label(i) == pop_labels.at(i));
            REQUIRE(bd.get_sequence_labels(i) == expected_seq_labels.at(i));
        }
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_population_label(npops), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(npops), std::out_of_range &);

        std::vector<unsigned int> expected_wts (1, nsites);

        std::vector< std::vector<unsigned int> > expected_allele_counts(1);
        expected_allele_counts[0] = std::vector<unsigned int>(npops, ngenomes);
        std::vector<unsigned int> expected_max_cts = expected_allele_counts.at(0);
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = nsites;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(1);
        expected_red_counts[0] = std::vector<unsigned int>(npops, 0);

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
            for (unsigned int pop_idx = 0; pop_idx < npops; ++pop_idx) {
                REQUIRE(bd.get_allele_count(pattern_idx, pop_idx) ==
                        expected_allele_counts.at(pattern_idx).at(pop_idx));
                REQUIRE(bd.get_red_allele_count(pattern_idx, pop_idx) ==
                        expected_red_counts.at(pattern_idx).at(pop_idx));
            }
        }
    }
}

TEST_CASE("Testing pop label constructor with 11 pops, 13 genomes, locus length 87, 53 loci", "[BiallelicData]") {

    SECTION("Testing label constructor") {
        unsigned int npops = 11;
        unsigned int ngenomes = 13;
        unsigned int nloci = 53;
        unsigned int locus_length = 87;
        unsigned int nsites = nloci * locus_length;
        std::vector<std::string> pop_labels;
        std::vector< std::vector<std::string> > expected_seq_labels;
        for (unsigned int i = 0; i < npops; ++i) {
            std::ostringstream p_label;
            p_label << "pop" << i;
            pop_labels.push_back(p_label.str());
            std::vector<std::string> seq_labels;
            for (unsigned int j = 0; j < ngenomes; ++j) {
                std::ostringstream s_label;
                s_label << p_label.str() << "-genome" << j;
                seq_labels.push_back(s_label.str());
            }
            expected_seq_labels.push_back(seq_labels);
        }
        BiallelicData bd(pop_labels, ngenomes, nloci, locus_length, true);
        REQUIRE(bd.get_number_of_populations() == npops);
        REQUIRE(bd.get_number_of_patterns() == 1);
        REQUIRE(bd.get_number_of_sites() == nsites);
        REQUIRE(bd.get_number_of_variable_sites() == 0);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.get_path() == "");
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends(nloci);
        unsigned int locus_end = locus_length - 1;
        expected_locus_ends.at(0) = locus_end;
        for (unsigned int i = 1; i < nloci; ++i) {
            locus_end += locus_length;
            expected_locus_ends.at(i) = locus_end;
        }
        std::vector<unsigned int> expected_pattern_indices(nsites, 0);
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_contiguous_pattern_indices().size() == nsites);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
        REQUIRE(bd.get_locus_end_indices().size() == nloci);
        REQUIRE(bd.get_locus_end_indices().back() == nsites - 1);

        for (unsigned int i = 0; i < npops; ++i) {
            REQUIRE(bd.get_population_index(pop_labels.at(i)) == i);
            REQUIRE(bd.get_population_label(i) == pop_labels.at(i));
            REQUIRE(bd.get_sequence_labels(i) == expected_seq_labels.at(i));
        }
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_population_label(npops), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(npops), std::out_of_range &);

        std::vector<unsigned int> expected_wts (1, nsites);

        std::vector< std::vector<unsigned int> > expected_allele_counts(1);
        expected_allele_counts[0] = std::vector<unsigned int>(npops, ngenomes);
        std::vector<unsigned int> expected_max_cts = expected_allele_counts.at(0);
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = nsites;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(1);
        expected_red_counts[0] = std::vector<unsigned int>(npops, 0);

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
            for (unsigned int pop_idx = 0; pop_idx < npops; ++pop_idx) {
                REQUIRE(bd.get_allele_count(pattern_idx, pop_idx) ==
                        expected_allele_counts.at(pattern_idx).at(pop_idx));
                REQUIRE(bd.get_red_allele_count(pattern_idx, pop_idx) ==
                        expected_red_counts.at(pattern_idx).at(pop_idx));
            }
        }
    }
}

TEST_CASE("Testing pop label constructor with 2 pops, 4-6 genomes, locus length 1", "[BiallelicData]") {

    SECTION("Testing label constructor") {
        unsigned int npops = 2;
        std::vector<unsigned int> ngenomes = {4, 6};
        unsigned int nloci = 10;
        unsigned int locus_length = 1;
        unsigned int nsites = nloci * locus_length;
        std::vector<std::string> pop_labels;
        std::vector< std::vector<std::string> > expected_seq_labels;
        for (unsigned int i = 0; i < npops; ++i) {
            std::ostringstream p_label;
            p_label << "pop" << i;
            pop_labels.push_back(p_label.str());
            std::vector<std::string> seq_labels;
            for (unsigned int j = 0; j < ngenomes.at(i); ++j) {
                std::ostringstream s_label;
                s_label << p_label.str() << "-genome" << j;
                seq_labels.push_back(s_label.str());
            }
            expected_seq_labels.push_back(seq_labels);
        }
        BiallelicData bd(pop_labels, ngenomes, nloci, locus_length, true);
        REQUIRE(bd.get_number_of_populations() == npops);
        REQUIRE(bd.get_number_of_patterns() == 1);
        REQUIRE(bd.get_number_of_sites() == nsites);
        REQUIRE(bd.get_number_of_variable_sites() == 0);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.get_path() == "");
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends(nloci);
        unsigned int locus_end = locus_length - 1;
        expected_locus_ends.at(0) = locus_end;
        for (unsigned int i = 1; i < nloci; ++i) {
            locus_end += locus_length;
            expected_locus_ends.at(i) = locus_end;
        }
        std::vector<unsigned int> expected_pattern_indices(nsites, 0);
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_contiguous_pattern_indices().size() == nsites);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
        REQUIRE(bd.get_locus_end_indices().size() == nloci);
        REQUIRE(bd.get_locus_end_indices().back() == nsites - 1);

        for (unsigned int i = 0; i < npops; ++i) {
            REQUIRE(bd.get_population_index(pop_labels.at(i)) == i);
            REQUIRE(bd.get_population_label(i) == pop_labels.at(i));
            REQUIRE(bd.get_sequence_labels(i) == expected_seq_labels.at(i));
        }
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_population_label(npops), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(npops), std::out_of_range &);

        std::vector<unsigned int> expected_wts (1, nsites);

        std::vector< std::vector<unsigned int> > expected_allele_counts(1);
        expected_allele_counts[0] = ngenomes;
        std::vector<unsigned int> expected_max_cts = expected_allele_counts.at(0);
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = nsites;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(1);
        expected_red_counts[0] = std::vector<unsigned int>(npops, 0);

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
            for (unsigned int pop_idx = 0; pop_idx < npops; ++pop_idx) {
                REQUIRE(bd.get_allele_count(pattern_idx, pop_idx) ==
                        expected_allele_counts.at(pattern_idx).at(pop_idx));
                REQUIRE(bd.get_red_allele_count(pattern_idx, pop_idx) ==
                        expected_red_counts.at(pattern_idx).at(pop_idx));
            }
        }
    }
}

TEST_CASE("Testing pop label constructor with 3 pops, 11-4-6 genomes, locus length 5", "[BiallelicData]") {

    SECTION("Testing label constructor") {
        unsigned int npops = 3;
        std::vector<unsigned int> ngenomes = {11, 4, 6};
        unsigned int nloci = 10;
        unsigned int locus_length = 5;
        unsigned int nsites = nloci * locus_length;
        std::vector<std::string> pop_labels;
        std::vector< std::vector<std::string> > expected_seq_labels;
        for (unsigned int i = 0; i < npops; ++i) {
            std::ostringstream p_label;
            p_label << "pop" << i;
            pop_labels.push_back(p_label.str());
            std::vector<std::string> seq_labels;
            for (unsigned int j = 0; j < ngenomes.at(i); ++j) {
                std::ostringstream s_label;
                s_label << p_label.str() << "-genome" << j;
                seq_labels.push_back(s_label.str());
            }
            expected_seq_labels.push_back(seq_labels);
        }
        BiallelicData bd(pop_labels, ngenomes, nloci, locus_length, true);
        REQUIRE(bd.get_number_of_populations() == npops);
        REQUIRE(bd.get_number_of_patterns() == 1);
        REQUIRE(bd.get_number_of_sites() == nsites);
        REQUIRE(bd.get_number_of_variable_sites() == 0);
        REQUIRE(! bd.markers_are_dominant());
        REQUIRE(! bd.genotypes_are_diploid());
        REQUIRE(bd.has_constant_patterns() == true);
        REQUIRE(bd.has_missing_population_patterns() == false);
        REQUIRE(bd.has_mirrored_patterns() == false);
        REQUIRE(bd.patterns_are_folded() == true);
        REQUIRE(bd.get_path() == "");
        REQUIRE(bd.has_recoded_triallelic_sites() == false);
        REQUIRE(bd.get_number_of_triallelic_sites_recoded() == 0);

        REQUIRE(bd.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends(nloci);
        unsigned int locus_end = locus_length - 1;
        expected_locus_ends.at(0) = locus_end;
        for (unsigned int i = 1; i < nloci; ++i) {
            locus_end += locus_length;
            expected_locus_ends.at(i) = locus_end;
        }
        std::vector<unsigned int> expected_pattern_indices(nsites, 0);
        REQUIRE(bd.get_contiguous_pattern_indices() == expected_pattern_indices);
        REQUIRE(bd.get_contiguous_pattern_indices().size() == nsites);
        REQUIRE(bd.get_locus_end_indices() == expected_locus_ends);
        REQUIRE(bd.get_locus_end_indices().size() == nloci);
        REQUIRE(bd.get_locus_end_indices().back() == nsites - 1);

        for (unsigned int i = 0; i < npops; ++i) {
            REQUIRE(bd.get_population_index(pop_labels.at(i)) == i);
            REQUIRE(bd.get_population_label(i) == pop_labels.at(i));
            REQUIRE(bd.get_sequence_labels(i) == expected_seq_labels.at(i));
        }
        REQUIRE_THROWS_AS(bd.get_population_index("bogus_label"), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_population_label(npops), std::out_of_range &);
        REQUIRE_THROWS_AS(bd.get_sequence_labels(npops), std::out_of_range &);

        std::vector<unsigned int> expected_wts (1, nsites);

        std::vector< std::vector<unsigned int> > expected_allele_counts(1);
        expected_allele_counts[0] = ngenomes;
        std::vector<unsigned int> expected_max_cts = expected_allele_counts.at(0);
        REQUIRE(bd.get_max_allele_counts() == expected_max_cts);

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = nsites;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::vector< std::vector<unsigned int> > expected_red_counts(1);
        expected_red_counts[0] = std::vector<unsigned int>(npops, 0);

        for (unsigned int pattern_idx = 0; pattern_idx < expected_wts.size(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == expected_wts.at(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == expected_allele_counts.at(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == expected_red_counts.at(pattern_idx));
            for (unsigned int pop_idx = 0; pop_idx < npops; ++pop_idx) {
                REQUIRE(bd.get_allele_count(pattern_idx, pop_idx) ==
                        expected_allele_counts.at(pattern_idx).at(pop_idx));
                REQUIRE(bd.get_red_allele_count(pattern_idx, pop_idx) ==
                        expected_red_counts.at(pattern_idx).at(pop_idx));
            }
        }
    }
}

TEST_CASE("Yaml parsing diploid dna with missing, mirrored, and constant sites", "[BiallelicData]") {

    SECTION("Testing data/diploid-dna-constant-missing.yml") {
        std::string nex_path = "data/diploid-dna-constant-missing.nex";
        std::string yml_path = "data/diploid-dna-constant-missing.yml";
        BiallelicData bd(nex_path);
        BiallelicData ybd;
        ybd.init_from_yaml_path(yml_path);
        REQUIRE(bd.get_number_of_populations() == ybd.get_number_of_populations());
        REQUIRE(bd.get_number_of_patterns() == ybd.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == ybd.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == ybd.markers_are_dominant());
        REQUIRE(bd.has_constant_patterns() == ybd.has_constant_patterns());
        REQUIRE(bd.has_missing_population_patterns() == ybd.has_missing_population_patterns());
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(ybd.get_path() == yml_path);
        REQUIRE(bd.has_mirrored_patterns() == ybd.has_mirrored_patterns());
        REQUIRE(bd.patterns_are_folded() == ybd.patterns_are_folded());

        REQUIRE(bd.get_unique_allele_counts() == ybd.get_unique_allele_counts());
        REQUIRE(bd.get_max_allele_counts() == ybd.get_max_allele_counts());

        for (unsigned int pattern_idx = 0; pattern_idx < bd.get_number_of_patterns(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == ybd.get_pattern_weight(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == ybd.get_allele_counts(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == ybd.get_red_allele_counts(pattern_idx));
        }

        REQUIRE(bd.get_population_index("pop1") == ybd.get_population_index("pop1"));
        REQUIRE(bd.get_population_index("pop2") == ybd.get_population_index("pop2"));
        REQUIRE(bd.get_population_index("pop3") == ybd.get_population_index("pop3"));
        REQUIRE(bd.get_population_label(0) == ybd.get_population_label(0));
        REQUIRE(bd.get_population_label(1) == ybd.get_population_label(1));
        REQUIRE(bd.get_population_label(2) == ybd.get_population_label(2));
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);
        REQUIRE_THROWS_AS(ybd.get_population_label(3), std::out_of_range &);

        // Folding
        unsigned int number_removed = bd.fold_patterns();
        unsigned int ynumber_removed = ybd.fold_patterns();
        REQUIRE(number_removed == ynumber_removed);

        REQUIRE(bd.get_number_of_populations() == ybd.get_number_of_populations());
        REQUIRE(bd.get_number_of_patterns() == ybd.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == ybd.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == ybd.markers_are_dominant());
        REQUIRE(bd.has_constant_patterns() == ybd.has_constant_patterns());
        REQUIRE(bd.has_missing_population_patterns() == ybd.has_missing_population_patterns());
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(ybd.get_path() == yml_path);
        REQUIRE(bd.has_mirrored_patterns() == ybd.has_mirrored_patterns());
        REQUIRE(bd.patterns_are_folded() == ybd.patterns_are_folded());

        REQUIRE(bd.get_unique_allele_counts() == ybd.get_unique_allele_counts());
        REQUIRE(bd.get_max_allele_counts() == ybd.get_max_allele_counts());

        for (unsigned int pattern_idx = 0; pattern_idx < bd.get_number_of_patterns(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == ybd.get_pattern_weight(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == ybd.get_allele_counts(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == ybd.get_red_allele_counts(pattern_idx));
        }

        REQUIRE(bd.get_population_index("pop1") == ybd.get_population_index("pop1"));
        REQUIRE(bd.get_population_index("pop2") == ybd.get_population_index("pop2"));
        REQUIRE(bd.get_population_index("pop3") == ybd.get_population_index("pop3"));
        REQUIRE(bd.get_population_label(0) == ybd.get_population_label(0));
        REQUIRE(bd.get_population_label(1) == ybd.get_population_label(1));
        REQUIRE(bd.get_population_label(2) == ybd.get_population_label(2));
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);
        REQUIRE_THROWS_AS(ybd.get_population_label(3), std::out_of_range &);

        // Remove missing
        number_removed = bd.remove_missing_population_patterns();
        ynumber_removed = ybd.remove_missing_population_patterns();
        REQUIRE(number_removed == ynumber_removed);

        REQUIRE(bd.get_number_of_populations() == ybd.get_number_of_populations());
        REQUIRE(bd.get_number_of_patterns() == ybd.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == ybd.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == ybd.markers_are_dominant());
        REQUIRE(bd.has_constant_patterns() == ybd.has_constant_patterns());
        REQUIRE(bd.has_missing_population_patterns() == ybd.has_missing_population_patterns());
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(ybd.get_path() == yml_path);
        REQUIRE(bd.has_mirrored_patterns() == ybd.has_mirrored_patterns());
        REQUIRE(bd.patterns_are_folded() == ybd.patterns_are_folded());

        REQUIRE(bd.get_unique_allele_counts() == ybd.get_unique_allele_counts());
        REQUIRE(bd.get_max_allele_counts() == ybd.get_max_allele_counts());

        for (unsigned int pattern_idx = 0; pattern_idx < bd.get_number_of_patterns(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == ybd.get_pattern_weight(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == ybd.get_allele_counts(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == ybd.get_red_allele_counts(pattern_idx));
        }

        REQUIRE(bd.get_population_index("pop1") == ybd.get_population_index("pop1"));
        REQUIRE(bd.get_population_index("pop2") == ybd.get_population_index("pop2"));
        REQUIRE(bd.get_population_index("pop3") == ybd.get_population_index("pop3"));
        REQUIRE(bd.get_population_label(0) == ybd.get_population_label(0));
        REQUIRE(bd.get_population_label(1) == ybd.get_population_label(1));
        REQUIRE(bd.get_population_label(2) == ybd.get_population_label(2));
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);
        REQUIRE_THROWS_AS(ybd.get_population_label(3), std::out_of_range &);

        // Remove constant
        number_removed = bd.remove_constant_patterns();
        ynumber_removed = ybd.remove_constant_patterns();
        REQUIRE(number_removed == ynumber_removed);

        REQUIRE(bd.get_number_of_populations() == ybd.get_number_of_populations());
        REQUIRE(bd.get_number_of_patterns() == ybd.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == ybd.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == ybd.markers_are_dominant());
        REQUIRE(bd.has_constant_patterns() == ybd.has_constant_patterns());
        REQUIRE(bd.has_missing_population_patterns() == ybd.has_missing_population_patterns());
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(ybd.get_path() == yml_path);
        REQUIRE(bd.has_mirrored_patterns() == ybd.has_mirrored_patterns());
        REQUIRE(bd.patterns_are_folded() == ybd.patterns_are_folded());

        REQUIRE(bd.get_unique_allele_counts() == ybd.get_unique_allele_counts());
        REQUIRE(bd.get_max_allele_counts() == ybd.get_max_allele_counts());

        for (unsigned int pattern_idx = 0; pattern_idx < bd.get_number_of_patterns(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == ybd.get_pattern_weight(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == ybd.get_allele_counts(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == ybd.get_red_allele_counts(pattern_idx));
        }

        REQUIRE(bd.get_population_index("pop1") == ybd.get_population_index("pop1"));
        REQUIRE(bd.get_population_index("pop2") == ybd.get_population_index("pop2"));
        REQUIRE(bd.get_population_index("pop3") == ybd.get_population_index("pop3"));
        REQUIRE(bd.get_population_label(0) == ybd.get_population_label(0));
        REQUIRE(bd.get_population_label(1) == ybd.get_population_label(1));
        REQUIRE(bd.get_population_label(2) == ybd.get_population_label(2));
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);
        REQUIRE_THROWS_AS(ybd.get_population_label(3), std::out_of_range &);

        REQUIRE(bd.get_max_allele_counts() == ybd.get_max_allele_counts());
    }
}

TEST_CASE("Yaml parsing default dominance diploid dna with missing, mirrored, and constant sites", "[BiallelicData]") {

    SECTION("Testing data/diploid-dna-constant-missing-default-dominance.yml") {
        std::string nex_path = "data/diploid-dna-constant-missing.nex";
        std::string yml_path = "data/diploid-dna-constant-missing-default-dominance.yml";
        BiallelicData bd(nex_path);
        BiallelicData ybd;
        ybd.init_from_yaml_path(yml_path);
        REQUIRE(bd.get_number_of_populations() == ybd.get_number_of_populations());
        REQUIRE(bd.get_number_of_patterns() == ybd.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == ybd.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == ybd.markers_are_dominant());
        REQUIRE(bd.has_constant_patterns() == ybd.has_constant_patterns());
        REQUIRE(bd.has_missing_population_patterns() == ybd.has_missing_population_patterns());
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(ybd.get_path() == yml_path);
        REQUIRE(bd.has_mirrored_patterns() == ybd.has_mirrored_patterns());
        REQUIRE(bd.patterns_are_folded() == ybd.patterns_are_folded());

        REQUIRE(bd.get_unique_allele_counts() == ybd.get_unique_allele_counts());
        REQUIRE(bd.get_max_allele_counts() == ybd.get_max_allele_counts());

        for (unsigned int pattern_idx = 0; pattern_idx < bd.get_number_of_patterns(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == ybd.get_pattern_weight(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == ybd.get_allele_counts(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == ybd.get_red_allele_counts(pattern_idx));
        }

        REQUIRE(bd.get_population_index("pop1") == ybd.get_population_index("pop1"));
        REQUIRE(bd.get_population_index("pop2") == ybd.get_population_index("pop2"));
        REQUIRE(bd.get_population_index("pop3") == ybd.get_population_index("pop3"));
        REQUIRE(bd.get_population_label(0) == ybd.get_population_label(0));
        REQUIRE(bd.get_population_label(1) == ybd.get_population_label(1));
        REQUIRE(bd.get_population_label(2) == ybd.get_population_label(2));
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);
        REQUIRE_THROWS_AS(ybd.get_population_label(3), std::out_of_range &);

        // Folding
        unsigned int number_removed = bd.fold_patterns();
        unsigned int ynumber_removed = ybd.fold_patterns();
        REQUIRE(number_removed == ynumber_removed);

        REQUIRE(bd.get_number_of_populations() == ybd.get_number_of_populations());
        REQUIRE(bd.get_number_of_patterns() == ybd.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == ybd.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == ybd.markers_are_dominant());
        REQUIRE(bd.has_constant_patterns() == ybd.has_constant_patterns());
        REQUIRE(bd.has_missing_population_patterns() == ybd.has_missing_population_patterns());
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(ybd.get_path() == yml_path);
        REQUIRE(bd.has_mirrored_patterns() == ybd.has_mirrored_patterns());
        REQUIRE(bd.patterns_are_folded() == ybd.patterns_are_folded());

        REQUIRE(bd.get_unique_allele_counts() == ybd.get_unique_allele_counts());
        REQUIRE(bd.get_max_allele_counts() == ybd.get_max_allele_counts());

        for (unsigned int pattern_idx = 0; pattern_idx < bd.get_number_of_patterns(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == ybd.get_pattern_weight(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == ybd.get_allele_counts(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == ybd.get_red_allele_counts(pattern_idx));
        }

        REQUIRE(bd.get_population_index("pop1") == ybd.get_population_index("pop1"));
        REQUIRE(bd.get_population_index("pop2") == ybd.get_population_index("pop2"));
        REQUIRE(bd.get_population_index("pop3") == ybd.get_population_index("pop3"));
        REQUIRE(bd.get_population_label(0) == ybd.get_population_label(0));
        REQUIRE(bd.get_population_label(1) == ybd.get_population_label(1));
        REQUIRE(bd.get_population_label(2) == ybd.get_population_label(2));
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);
        REQUIRE_THROWS_AS(ybd.get_population_label(3), std::out_of_range &);

        // Remove missing
        number_removed = bd.remove_missing_population_patterns();
        ynumber_removed = ybd.remove_missing_population_patterns();
        REQUIRE(number_removed == ynumber_removed);

        REQUIRE(bd.get_number_of_populations() == ybd.get_number_of_populations());
        REQUIRE(bd.get_number_of_patterns() == ybd.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == ybd.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == ybd.markers_are_dominant());
        REQUIRE(bd.has_constant_patterns() == ybd.has_constant_patterns());
        REQUIRE(bd.has_missing_population_patterns() == ybd.has_missing_population_patterns());
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(ybd.get_path() == yml_path);
        REQUIRE(bd.has_mirrored_patterns() == ybd.has_mirrored_patterns());
        REQUIRE(bd.patterns_are_folded() == ybd.patterns_are_folded());

        REQUIRE(bd.get_unique_allele_counts() == ybd.get_unique_allele_counts());
        REQUIRE(bd.get_max_allele_counts() == ybd.get_max_allele_counts());

        for (unsigned int pattern_idx = 0; pattern_idx < bd.get_number_of_patterns(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == ybd.get_pattern_weight(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == ybd.get_allele_counts(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == ybd.get_red_allele_counts(pattern_idx));
        }

        REQUIRE(bd.get_population_index("pop1") == ybd.get_population_index("pop1"));
        REQUIRE(bd.get_population_index("pop2") == ybd.get_population_index("pop2"));
        REQUIRE(bd.get_population_index("pop3") == ybd.get_population_index("pop3"));
        REQUIRE(bd.get_population_label(0) == ybd.get_population_label(0));
        REQUIRE(bd.get_population_label(1) == ybd.get_population_label(1));
        REQUIRE(bd.get_population_label(2) == ybd.get_population_label(2));
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);
        REQUIRE_THROWS_AS(ybd.get_population_label(3), std::out_of_range &);

        // Remove constant
        number_removed = bd.remove_constant_patterns();
        ynumber_removed = ybd.remove_constant_patterns();
        REQUIRE(number_removed == ynumber_removed);

        REQUIRE(bd.get_number_of_populations() == ybd.get_number_of_populations());
        REQUIRE(bd.get_number_of_patterns() == ybd.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == ybd.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == ybd.markers_are_dominant());
        REQUIRE(bd.has_constant_patterns() == ybd.has_constant_patterns());
        REQUIRE(bd.has_missing_population_patterns() == ybd.has_missing_population_patterns());
        REQUIRE(bd.get_path() == nex_path);
        REQUIRE(ybd.get_path() == yml_path);
        REQUIRE(bd.has_mirrored_patterns() == ybd.has_mirrored_patterns());
        REQUIRE(bd.patterns_are_folded() == ybd.patterns_are_folded());

        REQUIRE(bd.get_unique_allele_counts() == ybd.get_unique_allele_counts());
        REQUIRE(bd.get_max_allele_counts() == ybd.get_max_allele_counts());

        for (unsigned int pattern_idx = 0; pattern_idx < bd.get_number_of_patterns(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == ybd.get_pattern_weight(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == ybd.get_allele_counts(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == ybd.get_red_allele_counts(pattern_idx));
        }

        REQUIRE(bd.get_population_index("pop1") == ybd.get_population_index("pop1"));
        REQUIRE(bd.get_population_index("pop2") == ybd.get_population_index("pop2"));
        REQUIRE(bd.get_population_index("pop3") == ybd.get_population_index("pop3"));
        REQUIRE(bd.get_population_label(0) == ybd.get_population_label(0));
        REQUIRE(bd.get_population_label(1) == ybd.get_population_label(1));
        REQUIRE(bd.get_population_label(2) == ybd.get_population_label(2));
        REQUIRE_THROWS_AS(bd.get_population_label(3), std::out_of_range &);
        REQUIRE_THROWS_AS(ybd.get_population_label(3), std::out_of_range &);

        REQUIRE(bd.get_max_allele_counts() == ybd.get_max_allele_counts());
    }
}

TEST_CASE("Testing yaml parsing standard haploid dominant", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard-dominant.yml") {
        std::string nex_path = "data/haploid-standard.nex";
        std::string yml_path = "data/haploid-standard-dominant.yml";
        BiallelicData bd(nex_path, ' ', true, false, true);
        BiallelicData ybd;
        ybd.init_from_yaml_path(yml_path);
        REQUIRE(bd.get_number_of_populations() == ybd.get_number_of_populations());
        REQUIRE(bd.get_number_of_patterns() == ybd.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == ybd.get_number_of_sites());
        REQUIRE(bd.get_number_of_variable_sites() == ybd.get_number_of_variable_sites());
        REQUIRE(bd.markers_are_dominant() == ybd.markers_are_dominant());
        REQUIRE(bd.has_constant_patterns() == ybd.has_constant_patterns());
        REQUIRE(bd.has_missing_population_patterns() == ybd.has_missing_population_patterns());
        REQUIRE(bd.has_mirrored_patterns() == ybd.has_mirrored_patterns());
        REQUIRE(bd.patterns_are_folded() == ybd.patterns_are_folded());

        REQUIRE(bd.get_unique_allele_counts() == ybd.get_unique_allele_counts());

        for (unsigned int pattern_idx = 0; pattern_idx < bd.get_number_of_patterns(); ++pattern_idx) {
            REQUIRE(bd.get_pattern_weight(pattern_idx) == ybd.get_pattern_weight(pattern_idx));
            REQUIRE(bd.get_allele_counts(pattern_idx) == ybd.get_allele_counts(pattern_idx));
            REQUIRE(bd.get_red_allele_counts(pattern_idx) == ybd.get_red_allele_counts(pattern_idx));
        }

        REQUIRE(bd.get_population_index("pop1") == ybd.get_population_index("pop1"));
        REQUIRE(bd.get_population_index("pop2") == ybd.get_population_index("pop2"));
        REQUIRE(bd.get_population_label(0) == ybd.get_population_label(0));
        REQUIRE(bd.get_population_label(1) == ybd.get_population_label(1));


        // Folding
        REQUIRE_THROWS_AS(bd.fold_patterns(), EcoevolityBiallelicDataError &);
        REQUIRE_THROWS_AS(ybd.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}

TEST_CASE("Testing yaml parsing with duplicated pattern", "[BiallelicData]") {

    SECTION("Testing data/duplicated-pattern.yml") {
        std::string yml_path = "data/duplicated-pattern.yml";
        BiallelicData bd;
        REQUIRE_THROWS_AS(bd.init_from_yaml_path(yml_path), EcoevolityYamlDataError &);
    }
}

TEST_CASE("Testing yaml writing methods for diploid standard data set",
        "[BiallelicData]") {

    SECTION("Testing yaml writing of data/diploid-standard-data-ntax5-nchar5.nex") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                true,  // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        std::stringstream yaml;
        bd.write_yaml(yaml);
        std::cout << "\n" << yaml.str() << "\n";

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t927475.yml";
        std::ofstream os;
        os.open(test_path);
        bd.write_yaml(os);
        os.close();

        BiallelicData d;
        d.init_from_yaml_path(test_path);

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == d.markers_are_dominant());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing yaml writing methods for standard diploid with only 0/1 genotypes",
        "[BiallelicData]") {

    SECTION("Testing yaml writing of data/diploid-standard-only-01.nex") {
        std::string nex_path = "data/diploid-standard-only-01.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                true,  // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        std::stringstream yaml;
        bd.write_yaml(yaml);
        std::cout << "\n" << yaml.str() << "\n";

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t658972387.yml";
        std::ofstream os;
        os.open(test_path);
        bd.write_yaml(os);
        os.close();

        BiallelicData d;
        d.init_from_yaml_path(test_path);

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == d.markers_are_dominant());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing yaml writing methods for standard haploid", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard.nex") {
        std::string nex_path = "data/haploid-standard.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        std::stringstream yaml;
        bd.write_yaml(yaml);
        std::cout << "\n" << yaml.str() << "\n";

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t140879536.yml";
        std::ofstream os;
        os.open(test_path);
        bd.write_yaml(os);
        os.close();

        BiallelicData d;
        d.init_from_yaml_path(test_path);

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == d.markers_are_dominant());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing yaml writing methods for standard haploid dominant", "[BiallelicData]") {

    SECTION("Testing data/haploid-standard.nex as dominant") {
        std::string nex_path = "data/haploid-standard.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                true,  // markers are dominant (false)
                true); // validate (true)

        std::stringstream yaml;
        bd.write_yaml(yaml);
        std::cout << "\n" << yaml.str() << "\n";

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t47132890876.yml";
        std::ofstream os;
        os.open(test_path);
        bd.write_yaml(os);
        os.close();

        BiallelicData d;
        d.init_from_yaml_path(test_path);

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == d.markers_are_dominant());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing yaml writing methods for constant diploid site patterns",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant0.nex") {
        std::string nex_path = "data/diploid-standard-constant0.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                true,  // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t593874364.yml";
        std::ofstream os;
        os.open(test_path);
        bd.write_yaml(os);
        os.close();

        BiallelicData d;
        d.init_from_yaml_path(test_path);

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == d.markers_are_dominant());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }

        bd.fold_patterns();

        tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        test_path = "data/tmp-data-" + tag + "-t92745983274.cfg";
        os.open(test_path);
        bd.write_yaml(os);
        os.close();

        BiallelicData d2;
        d2.init_from_yaml_path(test_path);

        REQUIRE(bd.get_number_of_populations() == d2.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d2.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d2.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d2.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == d2.markers_are_dominant());

        for (unsigned int i = 0; i < d2.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d2.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d2.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d2.get_pattern_weight(i));
        }

        bd.remove_constant_patterns();

        tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        test_path = "data/tmp-data-" + tag + "-t987987234.cfg";
        os.open(test_path);
        bd.write_yaml(os);
        os.close();

        BiallelicData d3;
        d3.init_from_yaml_path(test_path);

        REQUIRE(bd.get_number_of_populations() == d3.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d3.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d3.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d3.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == d3.markers_are_dominant());

        for (unsigned int i = 0; i < d3.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d3.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d3.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d3.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing yaml writing methods for aflp data",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant0.nex") {
        std::string nex_path = "data/aflp_25-with-constant.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                false,  // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        std::stringstream yaml;
        bd.write_yaml(yaml);
        std::cout << "\n" << yaml.str() << "\n";

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t6876864877435.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_yaml(os);
        os.close();

        BiallelicData d;
        d.init_from_yaml_path(test_path);

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == d.markers_are_dominant());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing yaml writing methods for aflp data with folding",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant0.nex") {
        std::string nex_path = "data/aflp_25-with-constant.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        bd.fold_patterns();

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t8892787349.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_yaml(os);
        os.close();

        BiallelicData d;
        d.init_from_yaml_path(test_path);

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == d.markers_are_dominant());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing yaml writing methods for aflp data with removing",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant0.nex") {
        std::string nex_path = "data/aflp_25-with-constant.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        bd.remove_constant_patterns();

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t7676865398.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_yaml(os);
        os.close();

        BiallelicData d;
        d.init_from_yaml_path(test_path);

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == d.markers_are_dominant());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing yaml writing methods for aflp data with folding and removing",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant0.nex") {
        std::string nex_path = "data/aflp_25-with-constant.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        bd.fold_patterns();
        bd.remove_constant_patterns();

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t342398987961.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_yaml(os);
        os.close();

        BiallelicData d;
        d.init_from_yaml_path(test_path);

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == d.markers_are_dominant());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing yaml writing methods for dominant aflp data",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant0.nex") {
        std::string nex_path = "data/aflp_25-with-constant.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                false,  // genotypes are diploid (true)
                true,  // markers are dominant (false)
                true); // validate (true)

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t92734987232174.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_yaml(os);
        os.close();

        BiallelicData d;
        d.init_from_yaml_path(test_path);

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == d.markers_are_dominant());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing yaml writing methods for dominant aflp data with removing",
        "[BiallelicData]") {

    SECTION("Testing data/diploid-standard-constant0.nex") {
        std::string nex_path = "data/aflp_25-with-constant.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                false, // genotypes are diploid (true)
                true,  // markers are dominant (false)
                true); // validate (true)

        bd.remove_constant_patterns();

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t3749829093845.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_yaml(os);
        os.close();

        BiallelicData d;
        d.init_from_yaml_path(test_path);

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == d.markers_are_dominant());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing yaml writing methods for hemi data",
        "[BiallelicData]") {

    SECTION("Testing data/hemi129.nex") {
        std::string nex_path = "data/hemi129.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                true,  // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        std::stringstream yaml;
        bd.write_yaml(yaml);
        std::cout << "\n" << yaml.str() << "\n";

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t43298785456.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_yaml(os);
        os.close();

        BiallelicData d;
        d.init_from_yaml_path(test_path);

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == d.markers_are_dominant());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing yaml writing methods for hemi data with folding",
        "[BiallelicData]") {

    SECTION("Testing data/hemi129.nex") {
        std::string nex_path = "data/hemi129.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                true,  // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        bd.fold_patterns();

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t194378258667.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_yaml(os);
        os.close();

        BiallelicData d;
        d.init_from_yaml_path(test_path);

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == d.markers_are_dominant());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing yaml writing methods for hemi data with removing",
        "[BiallelicData]") {

    SECTION("Testing data/hemi129.nex") {
        std::string nex_path = "data/hemi129.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                true,  // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        bd.remove_constant_patterns();

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t87236470234.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_yaml(os);
        os.close();

        BiallelicData d;
        d.init_from_yaml_path(test_path);

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == d.markers_are_dominant());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}

TEST_CASE("Testing yaml writing methods for hemi data with folding and removing",
        "[BiallelicData]") {

    SECTION("Testing data/hemi129.nex") {
        std::string nex_path = "data/hemi129.nex";
        BiallelicData bd(
                nex_path,
                ' ',
                true,  // pop name is prefix (true)
                true,  // genotypes are diploid (true)
                false, // markers are dominant (false)
                true); // validate (true)

        bd.fold_patterns();
        bd.remove_constant_patterns();

        std::string tag = _ECOEVOLITY_DATA_RNG.random_string(10);
        std::string test_path = "data/tmp-data-" + tag + "-t73827648765.cfg";
        std::ofstream os;
        os.open(test_path);
        bd.write_yaml(os);
        os.close();

        BiallelicData d;
        d.init_from_yaml_path(test_path);

        REQUIRE(bd.get_number_of_populations() == d.get_number_of_populations());
        REQUIRE(bd.get_population_labels() == d.get_population_labels());
        REQUIRE(bd.get_number_of_patterns() == d.get_number_of_patterns());
        REQUIRE(bd.get_number_of_sites() == d.get_number_of_sites());
        REQUIRE(bd.markers_are_dominant() == d.markers_are_dominant());

        for (unsigned int i = 0; i < d.get_number_of_patterns(); ++i) {
            REQUIRE(bd.get_allele_counts(i) == d.get_allele_counts(i));
            REQUIRE(bd.get_red_allele_counts(i) == d.get_red_allele_counts(i));
            REQUIRE(bd.get_pattern_weight(i) == d.get_pattern_weight(i));
        }
    }
}
