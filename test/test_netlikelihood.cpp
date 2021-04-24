#include "catch.hpp"
#include "ecoevolity/netlikelihood.hpp"


/**
 * The probability that an allele came from the left / right parent is 4/12 /
 * 8/12, respectively, and we have the following conditional probs at the top
 * of the daugher branch:
 *
 *   Daughter prob  Ways to split                       Prob (over 1728)
 *   ngreen,nred
 *   0,0 = 1/12
 *                  0,0 | 0,0 = 1/12                    144
 *   1,0 = 1/12
 *                  1,0 | 0,0 = 1/12 * 4/12             48
 *                  0,0 | 1,0 = 1/12 * 8/12             96
 *   0,1 = 1/12
 *                  0,1 | 0,0 = 1/12 * 4/12             48
 *                  0,0 | 0,1 = 1/12 * 8/12             96
 *   2,0 = 2/12
 *                  2,0 | 0,0 = 2/12 * 4/12 * 4/12      32
 *                  0,0 | 2,0 = 2/12 * 8/12 * 8/12      128
 *                  1,0 | 1,0 = 2(2/12 * 4/12 * 8/12)   128
 *                              x 2 above because 2
 *                              ways for red alleles
 *                              to end up in each
 *                              parent
 *   1,1 = 3/12
 *                  1,1 | 0,0 = 3/12 * 4/12 * 4/12      48
 *                  0,0 | 1,1 = 3/12 * 8/12 * 8/12      192
 *                  1,0 | 0,1 = 3/12 * 4/12 * 8/12      96
 *                  0,1 | 1,0 = 3/12 * 4/12 * 8/12      96
 *   0,2 = 4/12
 *                  0,2 | 0,0 = 4/12 * 4/12 * 4/12      64
 *                  0,0 | 0,2 = 4/12 * 8/12 * 8/12      256
 *                  0,1 | 0,1 = 2(4/12 * 4/12 * 8/12)   256
 *                              x 2 above because 2
 *                              ways for green alleles
 *                              to end up in each
 *                              parent
 *                                                      Total = 1728/1728
 *
 * From above, we can calculate the conditional probability of all possible
 * allele patterns at the bottom of each parent branch (which is the goal of
 * this function)
 *
 *   Left parent bottom probs (over 1728)
 *   0,0 = 144+96+96+128+192+256 = 912
 *   1,0 = 48+128+96             = 272
 *   0,1 = 48+96+256             = 400
 *   2,0 = 32                    = 32
 *   1,1 = 48                    = 48
 *   0,2 = 64                    = 64
 *   total                       = 1728 / 1728
 *
 *   Right parent bottom probs (over 1728)
 *   0,0 = 144+48+48+32+48+64    = 384
 *   1,0 = 96+128+96             = 320
 *   0,1 = 96+96+256             = 448
 *   2,0 = 128                   = 128
 *   1,1 = 192                   = 192
 *   0,2 = 256                   = 256
 *   total                       = 1728 / 1728
 */
TEST_CASE("Testing split_top_of_branch_partials",
        "[SplitTopOfBranchPartials]") {

    SECTION("Testing split_top_of_branch_partials") {
        unsigned int max_num_alleles = 2;
        BiallelicPatternProbabilityMatrix top_child_partials;
        top_child_partials.resize(max_num_alleles);
        top_child_partials.set_pattern_probability(0, 0, 1/12.0);
        top_child_partials.set_pattern_probability(1, 0, 1/12.0);
        top_child_partials.set_pattern_probability(1, 1, 1/12.0);
        top_child_partials.set_pattern_probability(2, 0, 2/12.0);
        top_child_partials.set_pattern_probability(2, 1, 3/12.0);
        top_child_partials.set_pattern_probability(2, 2, 4/12.0);
        BiallelicPatternProbabilityMatrix bottom_parent1_partials;
        BiallelicPatternProbabilityMatrix bottom_parent2_partials;
        double prob_to_parent1 = 4/12;
        double prob_to_parent2 = 8/12;
        netlikelihood::split_top_of_branch_partials(
                max_num_alleles,
                top_child_partials,
                prob_to_parent1,
                prob_to_parent2,
                bottom_parent1_partials,
                bottom_parent2_partials);

        REQUIRE(bottom_parent1_partials.get_pattern_probability(0,0) == Approx(912/1728.0).epsilon(1e-8));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(1,0) == Approx(272/1728.0).epsilon(1e-8));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(1,1) == Approx(400/1728.0).epsilon(1e-8));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(2,0) == Approx(32/1728.0).epsilon(1e-8));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(2,1) == Approx(48/1728.0).epsilon(1e-8));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(2,2) == Approx(64/1728.0).epsilon(1e-8));

        REQUIRE(bottom_parent2_partials.get_pattern_probability(0,0) == Approx(384/1728.0).epsilon(1e-8));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(1,0) == Approx(320/1728.0).epsilon(1e-8));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(1,1) == Approx(448/1728.0).epsilon(1e-8));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(2,0) == Approx(128/1728.0).epsilon(1e-8));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(2,1) == Approx(192/1728.0).epsilon(1e-8));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(2,2) == Approx(256/1728.0).epsilon(1e-8));
    }
}

/**
 * The probability that an allele came from the left / right parent is 4/12 /
 * 8/12, respectively, and we have the following conditional probs at the top
 * of the daugher branch:
 *
 *   Child 1
 *   ngreen,nred
 *   0,0 = 0
 *   1,0 = 1/12
 *   0,1 = 2/12
 *   2,0 = 2/12
 *   1,1 = 3/12
 *   0,2 = 4/12
 *
 *   Child 2
 *   ngreen,nred
 *   0,0 = 0
 *   1,0 = 4/12
 *   0,1 = 3/12
 *   2,0 = 2/12
 *   1,1 = 2/12
 *   0,2 = 1/12
 *
 *   Merged
 *   0,0 = 0
 *   1,0 = 0
 *   0,1 = 0
 *   2,0 = (1/12 * 4/12)                                                  = 4
 *   1,1 = (1/12 * 3/12) + (2/12 * 4/12)                                  = 11
 *   0,2 = (2/12 * 3/12)                                                  = 6
 *   3,0 = (1/12 * 2/12) + (2/12 * 4/12)                                  = 10
 *   2,1 = (2/12 * 3/12) + (2/12 * 2/12) + (3/12 * 4/12) + (1/12 * 2/12)  = 24
 *   1,2 = (1/12 * 1/12) + (4/12 + 4/12) + (3/12 * 3/12) + (2/12 * 2/12)  = 30
 *   0,3 = (4/12 * 3/12) + (1/12 + 2/12)                                  = 14
 *   4,0 = (2/12 * 2/12)                                                  = 4
 *   3,1 = (2/12 * 2/12) + (2/12 * 3/12)                                  = 10
 *   2,2 = (2/12 * 1/12) + (2/12 * 4/12) + (3/12 * 2/12)                  = 16
 *   1,3 = (4/12 * 2/12) + (1/12 * 3/12)                                  = 11
 *   0,4 = (4/12 * 1/12)                                                  = 4
 *   TOTAL                                                                = 144/144
 */
TEST_CASE("Testing merge_top_of_branch_partials with no missing prob",
        "[MergeTopOfBranchPartials]") {

    SECTION("Testing merge_top_of_branch_partials") {
        unsigned int max_num_alleles = 2;
        BiallelicPatternProbabilityMatrix top_child1_partials;
        top_child1_partials.resize(max_num_alleles);
        top_child1_partials.set_pattern_probability(0, 0, 0.0);
        top_child1_partials.set_pattern_probability(1, 0, 1/12.0);
        top_child1_partials.set_pattern_probability(1, 1, 2/12.0);
        top_child1_partials.set_pattern_probability(2, 0, 2/12.0);
        top_child1_partials.set_pattern_probability(2, 1, 3/12.0);
        top_child1_partials.set_pattern_probability(2, 2, 4/12.0);
        BiallelicPatternProbabilityMatrix top_child2_partials;
        top_child2_partials.resize(max_num_alleles);
        top_child2_partials.set_pattern_probability(0, 0, 0.0);
        top_child2_partials.set_pattern_probability(1, 0, 4/12.0);
        top_child2_partials.set_pattern_probability(1, 1, 3/12.0);
        top_child2_partials.set_pattern_probability(2, 0, 2/12.0);
        top_child2_partials.set_pattern_probability(2, 1, 2/12.0);
        top_child2_partials.set_pattern_probability(2, 2, 1/12.0);

        std::vector<double> top_ch1_partials = top_child1_partials.get_pattern_prob_matrix();
        std::vector<double> top_ch2_partials = top_child2_partials.get_pattern_prob_matrix();

        std::vector<double> merged_probs;

        unsigned int max_n_alleles_merged;
        netlikelihood::merge_top_of_branch_partials(
                max_num_alleles,
                max_num_alleles,
                top_ch1_partials,
                top_ch2_partials,
                max_n_alleles_merged,
                merged_probs);
        REQUIRE(max_n_alleles_merged == 4);
        BiallelicPatternProbabilityMatrix merged(max_n_alleles_merged, merged_probs);
        REQUIRE(merged.get_pattern_probability(0,0) == 0.0);
        REQUIRE(merged.get_pattern_probability(1,0) == 0.0);
        REQUIRE(merged.get_pattern_probability(0,1) == 0.0);
        REQUIRE(merged.get_pattern_probability(2,0) == Approx(4/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(2,1) == Approx(11/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(2,2) == Approx(6/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(3,0) == Approx(10/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(3,1) == Approx(24/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(3,2) == Approx(30/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(3,3) == Approx(14/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,0) == Approx(4/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,1) == Approx(10/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,2) == Approx(16/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,3) == Approx(11/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,4) == Approx(4/144.0).epsilon(1e-8));
    }
}

/**
 * The probability that an allele came from the left / right parent is 4/12 /
 * 8/12, respectively, and we have the following conditional probs at the top
 * of the daugher branch:
 *
 *   Child 1
 *   ngreen,nred
 *   0,0 = 1/12
 *   1,0 = 1/12
 *   0,1 = 1/12
 *   2,0 = 2/12
 *   1,1 = 3/12
 *   0,2 = 4/12
 *
 *   Child 2
 *   ngreen,nred
 *   0,0 = 2/12
 *   1,0 = 4/12
 *   0,1 = 3/12
 *   2,0 = 1/12
 *   1,1 = 1/12
 *   0,2 = 1/12
 *
 *   Merged
 *   0,0 = (1 * 2)                                  = 2
 *   1,0 = (1 * 2) + (4 * 1)                        = 6
 *   0,1 = (1 * 2) + (3 * 1)                        = 5
 *   2,0 = (2 * 2) + (1 * 1) + (1 * 4)              = 9
 *   1,1 = (3 * 2) + (1 * 1) + (1 * 3) + (1 * 4)    = 14
 *   0,2 = (4 * 2) + (1 * 1) + (1 * 3)              = 12
 *   3,0 = (2 * 4) + (1 * 1)                        = 9
 *   2,1 = (2 * 3) + (1 * 1) + (3 * 4) + (1 * 1)    = 20
 *   1,2 = (1 * 1) + (4 * 4) + (3 * 3) + (1 * 1)    = 27
 *   0,3 = (4 * 3) + (1 * 1)                        = 13
 *   4,0 = (2 * 1)                                  = 2
 *   3,1 = (2 * 1) + (1 * 3)                        = 5
 *   2,2 = (2 * 1) + (4 * 1) + (3 * 1)              = 9
 *   1,3 = (3 * 1) + (1 * 4)                        = 7
 *   0,4 = (4 * 1)                                  = 4
 *   TOTAL                                          = 144/144
 */
TEST_CASE("Testing merge_top_of_branch_partials with missing probs",
        "[MergeTopOfBranchPartials]") {

    SECTION("Testing merge_top_of_branch_partials") {
        unsigned int max_num_alleles = 2;
        BiallelicPatternProbabilityMatrix top_child1_partials;
        top_child1_partials.resize(max_num_alleles);
        top_child1_partials.set_pattern_probability(0, 0, 1/12.0);
        top_child1_partials.set_pattern_probability(1, 0, 1/12.0);
        top_child1_partials.set_pattern_probability(1, 1, 1/12.0);
        top_child1_partials.set_pattern_probability(2, 0, 2/12.0);
        top_child1_partials.set_pattern_probability(2, 1, 3/12.0);
        top_child1_partials.set_pattern_probability(2, 2, 4/12.0);
        BiallelicPatternProbabilityMatrix top_child2_partials;
        top_child2_partials.resize(max_num_alleles);
        top_child2_partials.set_pattern_probability(0, 0, 2/12.0);
        top_child2_partials.set_pattern_probability(1, 0, 4/12.0);
        top_child2_partials.set_pattern_probability(1, 1, 3/12.0);
        top_child2_partials.set_pattern_probability(2, 0, 1/12.0);
        top_child2_partials.set_pattern_probability(2, 1, 1/12.0);
        top_child2_partials.set_pattern_probability(2, 2, 1/12.0);

        std::vector<double> top_ch1_partials = top_child1_partials.get_pattern_prob_matrix();
        std::vector<double> top_ch2_partials = top_child2_partials.get_pattern_prob_matrix();

        std::vector<double> merged_probs;

        unsigned int max_n_alleles_merged;
        netlikelihood::merge_top_of_branch_partials(
                max_num_alleles,
                max_num_alleles,
                top_ch1_partials,
                top_ch2_partials,
                max_n_alleles_merged,
                merged_probs);
        REQUIRE(max_n_alleles_merged == 4);
        BiallelicPatternProbabilityMatrix merged(max_n_alleles_merged, merged_probs);
        REQUIRE(merged.get_pattern_probability(0,0) == Approx(2/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(1,0) == Approx(6/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(1,1) == Approx(5/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(2,0) == Approx(9/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(2,1) == Approx(14/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(2,2) == Approx(12/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(3,0) == Approx(9/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(3,1) == Approx(20/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(3,2) == Approx(27/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(3,3) == Approx(13/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,0) == Approx(2/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,1) == Approx(5/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,2) == Approx(9/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,3) == Approx(7/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,4) == Approx(4/144.0).epsilon(1e-8));
    }
}
