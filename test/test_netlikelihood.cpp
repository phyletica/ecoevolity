#include "catch.hpp"
#include "ecoevolity/netlikelihood.hpp"
#include "utils_for_testing.hpp"


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
        double prob_to_parent1 = 4/12.0;
        double prob_to_parent2 = 8/12.0;
        netlikelihood::split_top_of_branch_partials(
                max_num_alleles,
                top_child_partials,
                prob_to_parent1,
                prob_to_parent2,
                bottom_parent1_partials,
                bottom_parent2_partials);

        REQUIRE(bottom_parent1_partials.get_pattern_probability(0,0) == Approx(912/1728.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(1,0) == Approx(272/1728.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(1,1) == Approx(400/1728.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(2,0) == Approx(32/1728.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(2,1) == Approx(48/1728.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(2,2) == Approx(64/1728.0).epsilon(1e-10));

        REQUIRE(bottom_parent2_partials.get_pattern_probability(0,0) == Approx(384/1728.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(1,0) == Approx(320/1728.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(1,1) == Approx(448/1728.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(2,0) == Approx(128/1728.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(2,1) == Approx(192/1728.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(2,2) == Approx(256/1728.0).epsilon(1e-10));
    }
}


/**
 * The probability that an allele came from the left / right parent is 0 /
 * 1, respectively, and we have the following conditional probs at the top
 * of the daugher branch:
 *
 *   Daughter prob  Ways to split                       Prob (over 1728)
 *   ngreen,nred
 *   0,0 = 1/12
 *                  0,0 | 0,0 = 1/12                    144
 *   1,0 = 1/12
 *                  1,0 | 0,0 = 1/12 * 0                0
 *                  0,0 | 1,0 = 1/12 * 1                144
 *   0,1 = 1/12
 *                  0,1 | 0,0 = 1/12 * 0                0
 *                  0,0 | 0,1 = 1/12 * 1                144
 *   2,0 = 2/12
 *                  2,0 | 0,0 = 2/12 * 0 * 0            0
 *                  0,0 | 2,0 = 2/12 * 1 * 1            288
 *                  1,0 | 1,0 = 2(2/12 * 0 * 1)         0
 *                              x 2 above because 2
 *                              ways for red alleles
 *                              to end up in each
 *                              parent
 *   1,1 = 3/12
 *                  1,1 | 0,0 = 3/12 * 0 * 0            0
 *                  0,0 | 1,1 = 3/12 * 1 * 1            432
 *                  1,0 | 0,1 = 3/12 * 0 * 1            0
 *                  0,1 | 1,0 = 3/12 * 0 * 1            0
 *   0,2 = 4/12
 *                  0,2 | 0,0 = 4/12 * 0 * 0            0
 *                  0,0 | 0,2 = 4/12 * 1 * 1            576
 *                  0,1 | 0,1 = 2(4/12 * 0 * 1)         0
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
 *   0,0 = 144+144+144+288+432+576 = 1728
 *   1,0 =                         = 0
 *   0,1 =                         = 0
 *   2,0 =                         = 0
 *   1,1 =                         = 0
 *   0,2 =                         = 0
 *   total                         = 1728 / 1728
 *
 *   Right parent bottom probs (over 1728)
 *   0,0 = 144                   = 144
 *   1,0 = 144                   = 144
 *   0,1 = 144                   = 144
 *   2,0 = 288                   = 288
 *   1,1 = 432                   = 432
 *   0,2 = 576                   = 576
 *   total                       = 1728 / 1728
 */
TEST_CASE("Testing split_top_of_branch_partials with all to one parent",
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
        double prob_to_parent1 = 0.0;
        double prob_to_parent2 = 1.0;
        netlikelihood::split_top_of_branch_partials(
                max_num_alleles,
                top_child_partials,
                prob_to_parent1,
                prob_to_parent2,
                bottom_parent1_partials,
                bottom_parent2_partials);

        REQUIRE(bottom_parent1_partials.get_pattern_probability(0,0) == Approx(1.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(1,0) == Approx(0.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(1,1) == Approx(0.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(2,0) == Approx(0.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(2,1) == Approx(0.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(2,2) == Approx(0.0).epsilon(1e-10));

        REQUIRE(bottom_parent2_partials.get_pattern_probability(0,0) == Approx(144/1728.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(1,0) == Approx(144/1728.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(1,1) == Approx(144/1728.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(2,0) == Approx(288/1728.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(2,1) == Approx(432/1728.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(2,2) == Approx(576/1728.0).epsilon(1e-10));
    }
}


/**
 * The probability that an allele came from the left / right parent is 9/15 /
 * 6/15, respectively, and we have the following conditional probs at the top
 * of the daugher branch:
 *
 *   Daughter prob  Ways to split               Prob (over 759375)
 *   ngreen,nred
 *   0,0 = 1/15
 *                  0,0 | 0,0 = 1               50625
 *   1,0 = 1/15
 *                  1,0 | 0,0 = 1 * 9           30375
 *                  0,0 | 1,0 = 1 * 6           20250
 *   0,1 = 1/15
 *                  0,1 | 0,0 = 1 * 9           30375
 *                  0,0 | 0,1 = 1 * 6           20250
 *   2,0 = 1/15
 *                  2,0 | 0,0 = 1 * 9^2         18225
 *                  0,0 | 2,0 = 1 * 6^2         8100
 *                  1,0 | 1,0 = 2(1 * 9 * 6)    24300
 *   1,1 = 1/15
 *                  1,1 | 0,0 = 1 * 9^2         18225
 *                  0,0 | 1,1 = 1 * 6^2         8100
 *                  1,0 | 0,1 = 1 * 9 * 6       12150
 *                  0,1 | 1,0 = 1 * 9 * 6       12150
 *   0,2 = 1/15
 *                  0,2 | 0,0 = 1 * 9^2         18225
 *                  0,0 | 0,2 = 1 * 6^2         8100
 *                  0,1 | 0,1 = 2(1 * 9 * 6)    24300
 *   3,0 = 1/15
 *                  3,0 | 0,0 = 1 * 9^3         10935
 *                  2,0 | 1,0 = 3(1 * 9^2 * 6)  21870
 *                  1,0 | 2,0 = 3(1 * 9 * 6^2)  14580
 *                  0,0 | 3,0 = 1 * 6^3         3240
 *   2,1 = 1/15
 *                  2,1 | 0,0 = 1 * 9^3         10935
 *                  1,1 | 1,0 = 2(1 * 9^2 * 6)  14580
 *                  0,1 | 2,0 = 1 * 9 * 6^2     4860
 *                  2,0 | 0,1 = 1 * 9^2 * 6     7290
 *                  0,0 | 2,1 = 1 * 6^3         3240
 *                  1,0 | 1,1 = 2(1 * 9 & 6^2)  9720
 *   1,2 = 1/15
 *                  1,2 | 0,0 = 1 * 9^3         10935
 *                  1,1 | 0,1 = 2(1 * 9^2 * 6)  14580
 *                  1,0 | 0,2 = 1 * 9 * 6^2     4860
 *                  0,2 | 1,0 = 1 * 9^2 * 6     7290
 *                  0,0 | 1,2 = 1 * 6^3         3240
 *                  0,1 | 1,1 = 2(1 * 9 * 6^2)  9720
 *   0,3 = 1/15
 *                  0,3 | 0,0 = 9^3             10935
 *                  0,2 | 0,1 = 3(9^2 * 6)      21870
 *                  0,1 | 0,2 = 3(9 * 6^2)      14580
 *                  0,0 | 0,3 = 6^3             3240
 *   4,0 = 1/15
 *                  4,0 | 0,0 = 9^4             6561
 *                  3,0 | 1,0 = 4(9^3 * 6)      17496
 *                  2,0 | 2,0 = 6(9^2 * 6^2)    17496
 *                  1,0 | 3,0 = 4(9 * 6^3)      7776
 *                  0,0 | 4,0 = 6^4             1296
 *   3,1 = 1/15
 *                  3,1 | 0,0 = 9^4             6561
 *                  0,0 | 3,1 = 6^4             1296
 *                  3,0 | 0,1 = 9^3 * 6         4374
 *                  0,1 | 3,0 = 9 * 6^3         1944
 *                  2,1 | 1,0 = 3(9^3 * 6)      13122
 *                  1,0 | 2,1 = 3(9 * 6^3)      5832
 *                  1,1 | 2,0 = 3(9^2 * 6^2)    8748
 *                  2,0 | 1,1 = 3(9^2 * 6^2)    8748
 *   2,2 = 1/15
 *                  2,2 | 0,0 = 9^4             6561
 *                  0,0 | 2,2 = 6^4             1296
 *                  2,1 | 0,1 = 2(9^3 * 6)      8748
 *                  0,1 | 2,1 = 2(9 * 6^3)      3888
 *                  1,2 | 1,0 = 2(9^3 * 6)      8748
 *                  1,0 | 1,2 = 2(9 * 6^3)      3888
 *                  2,0 | 0,2 = 9^2 * 6^2       2916
 *                  0,2 | 2,0 = 9^2 * 6^2       2916
 *                  1,1 | 1,1 = 2*2(9^2 * 6^2)  11664
 *   1,3 = 1/15
 *                  1,3 | 0,0 = 9^4             6561
 *                  0,0 | 1,3 = 6^4             1296
 *                  0,3 | 1,0 = 9^3 * 6         4374
 *                  1,0 | 0,3 = 9 * 6^3         1944
 *                  1,2 | 0,1 = 3(9^3 * 6)      13122
 *                  0,1 | 1,2 = 3(9 * 6^3)      5832
 *                  1,1 | 0,2 = 3(9^2 * 6^2)    8748
 *                  0,2 | 1,1 = 3(9^2 * 6^2)    8748
 *   0,4 = 1/15
 *                  0,4 | 0,0 = 9^4             6561
 *                  0,3 | 0,1 = 4(9^3 * 6)      17496
 *                  0,2 | 0,2 = 6(9^2 * 6^2)    17496
 *                  0,1 | 0,3 = 4(9 * 6^3)      7776
 *                  0,0 | 0,4 = 6^4             1296
 *                                              Total = 759375
 *
 * From above, we can calculate the conditional probability of all possible
 * allele patterns at the bottom of each parent branch (which is the goal of
 * this function)
 *
 *   Left parent:
 *     0,0 134865
 *     0,1 115425
 *     0,2 76545
 *     0,3 32805
 *     0,4 6561
 *     1,0 115425
 *     1,1 76545
 *     1,2 32805
 *     1,3 6561
 *     2,0 76545
 *     2,1 32805
 *     2,2 6561
 *     3,0 32805
 *     3,1 6561
 *     4,0 6561
 *     Total = 759375
 *
 *   Rigth parent:
 *     0,0 242595
 *     0,1 144180
 *     0,2 56700
 *     0,3 12960
 *     0,4 1296
 *     1,0 144180
 *     1,1 56700
 *     1,2 12960
 *     1,3 1296
 *     2,0 56700
 *     2,1 12960
 *     2,2 1296
 *     3,0 12960
 *     3,1 1296
 *     4,0 1296
 *     Total = 759375
 */
TEST_CASE("Testing split_top_of_branch_partials with more alleles",
        "[SplitTopOfBranchPartials]") {

    SECTION("Testing split_top_of_branch_partials") {
        unsigned int max_num_alleles = 4;
        BiallelicPatternProbabilityMatrix top_child_partials;
        top_child_partials.resize(max_num_alleles);
        top_child_partials.set_pattern_probability(0, 0, 1/15.0);
        top_child_partials.set_pattern_probability(1, 0, 1/15.0);
        top_child_partials.set_pattern_probability(1, 1, 1/15.0);
        top_child_partials.set_pattern_probability(2, 0, 1/15.0);
        top_child_partials.set_pattern_probability(2, 1, 1/15.0);
        top_child_partials.set_pattern_probability(2, 2, 1/15.0);
        top_child_partials.set_pattern_probability(3, 0, 1/15.0);
        top_child_partials.set_pattern_probability(3, 1, 1/15.0);
        top_child_partials.set_pattern_probability(3, 2, 1/15.0);
        top_child_partials.set_pattern_probability(3, 3, 1/15.0);
        top_child_partials.set_pattern_probability(4, 0, 1/15.0);
        top_child_partials.set_pattern_probability(4, 1, 1/15.0);
        top_child_partials.set_pattern_probability(4, 2, 1/15.0);
        top_child_partials.set_pattern_probability(4, 3, 1/15.0);
        top_child_partials.set_pattern_probability(4, 4, 1/15.0);
        BiallelicPatternProbabilityMatrix bottom_parent1_partials;
        BiallelicPatternProbabilityMatrix bottom_parent2_partials;
        double prob_to_parent1 = 9/15.0;
        double prob_to_parent2 = 6/15.0;
        netlikelihood::split_top_of_branch_partials(
                max_num_alleles,
                top_child_partials,
                prob_to_parent1,
                prob_to_parent2,
                bottom_parent1_partials,
                bottom_parent2_partials);

        REQUIRE(bottom_parent1_partials.get_pattern_probability(0,0) == Approx(134865/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(1,1) == Approx(115425/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(2,2) == Approx(76545/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(3,3) == Approx(32805/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(4,4) == Approx(6561/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(1,0) == Approx(115425/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(2,1) == Approx(76545/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(3,2) == Approx(32805/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(4,3) == Approx(6561/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(2,0) == Approx(76545/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(3,1) == Approx(32805/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(4,2) == Approx(6561/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(3,0) == Approx(32805/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(4,1) == Approx(6561/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent1_partials.get_pattern_probability(4,0) == Approx(6561/759375.0).epsilon(1e-10));

        REQUIRE(bottom_parent2_partials.get_pattern_probability(0,0) == Approx(242595/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(1,1) == Approx(144180/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(2,2) == Approx(56700/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(3,3) == Approx(12960/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(4,4) == Approx(1296/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(1,0) == Approx(144180/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(2,1) == Approx(56700/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(3,2) == Approx(12960/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(4,3) == Approx(1296/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(2,0) == Approx(56700/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(3,1) == Approx(12960/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(4,2) == Approx(1296/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(3,0) == Approx(12960/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(4,1) == Approx(1296/759375.0).epsilon(1e-10));
        REQUIRE(bottom_parent2_partials.get_pattern_probability(4,0) == Approx(1296/759375.0).epsilon(1e-10));
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
 *   1,2 = (1/12 * 1/12) + (4/12 * 4/12) + (3/12 * 3/12) + (2/12 * 2/12)  = 30
 *   0,3 = (4/12 * 3/12) + (1/12 * 2/12)                                  = 14
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
        double merged_prob_no_alleles;

        unsigned int max_n_alleles_merged;
        netlikelihood::merge_top_of_branch_partials(
                max_num_alleles,
                max_num_alleles,
                top_child1_partials.get_pattern_probability(0,0),
                top_child2_partials.get_pattern_probability(0,0),
                top_ch1_partials,
                top_ch2_partials,
                max_n_alleles_merged,
                merged_probs,
                merged_prob_no_alleles,
                false);
        REQUIRE(max_n_alleles_merged == 4);
        BiallelicPatternProbabilityMatrix merged(max_n_alleles_merged, merged_probs);

        std::cout << "0,0 : " << merged.get_pattern_probability(0,0) * 144 << "\n";
        std::cout << "1,0 : " << merged.get_pattern_probability(1,0) * 144 << "\n";
        std::cout << "1,1 : " << merged.get_pattern_probability(1,1) * 144 << "\n";
        std::cout << "2,0 : " << merged.get_pattern_probability(2,0) * 144 << "\n";
        std::cout << "2,1 : " << merged.get_pattern_probability(2,1) * 144 << "\n";
        std::cout << "2,2 : " << merged.get_pattern_probability(2,2) * 144 << "\n";
        std::cout << "3,0 : " << merged.get_pattern_probability(3,0) * 144 << "\n";
        std::cout << "3,1 : " << merged.get_pattern_probability(3,1) * 144 << "\n";
        std::cout << "3,2 : " << merged.get_pattern_probability(3,2) * 144 << "\n";
        std::cout << "3,3 : " << merged.get_pattern_probability(3,3) * 144 << "\n";
        std::cout << "4,0 : " << merged.get_pattern_probability(4,0) * 144 << "\n";
        std::cout << "4,1 : " << merged.get_pattern_probability(4,1) * 144 << "\n";
        std::cout << "4,2 : " << merged.get_pattern_probability(4,2) * 144 << "\n";
        std::cout << "4,3 : " << merged.get_pattern_probability(4,3) * 144 << "\n";
        std::cout << "4,4 : " << merged.get_pattern_probability(4,4) * 144 << "\n";

        REQUIRE(merged.get_pattern_probability(0,0) == 0.0);
        REQUIRE(merged.get_pattern_probability(1,0) == 0.0);
        REQUIRE(merged.get_pattern_probability(1,1) == 0.0);
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
 * Below, the factor under terms is the hypergeometric probability:
 *
 *   ((n1 choose r1) * (n2 choose r2)) / (n choose r)
 *
 * See Eq. 19 in Bryant et al. 2012
 *   
 *
 *   Merged
 *   0,0 = 0
 *   1,0 = 0
 *   0,1 = 0
 *   2,0 = (1/12 * 4/12)                                                  = 4
 *   1,1 = (1/12 * 3/12) + (2/12 * 4/12)                                  = 5.5
 *            * 1/2            * 1/2
 *   0,2 = (2/12 * 3/12)                                                  = 6
 *   3,0 = (1/12 * 2/12) + (2/12 * 4/12)                                  = 10
 *   2,1 = (2/12 * 3/12) + (2/12 * 2/12) + (3/12 * 4/12) + (1/12 * 2/12)  = 38/3
 *            * 1/3            * 1/3           * 2/3           * 2/3
 *   1,2 = (1/12 * 1/12) + (4/12 * 4/12) + (3/12 * 3/12) + (2/12 * 2/12)  = 43/3
 *            * 1/3           * 1/3            * 2/3           * 2/3
 *   0,3 = (4/12 * 3/12) + (1/12 * 2/12)                                  = 14
 *   4,0 = (2/12 * 2/12)                                                  = 4
 *   3,1 = (2/12 * 2/12) + (2/12 * 3/12)                                  = 5
 *            * 1/2           * 1/2
 *   2,2 = (2/12 * 1/12) + (2/12 * 4/12) + (3/12 * 2/12)                  = 34/6
 *            * 1/6           * 1/6            * 4/6
 *   1,3 = (4/12 * 2/12) + (1/12 * 3/12)                                  = 5.5
 *            * 1/2           * 1/2
 *   0,4 = (4/12 * 1/12)                                                  = 4
 *   TOTAL                                                                = 144/144
 */
TEST_CASE("Testing merge_top_of_branch_partials with no missing prob and rescaling",
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
        double merged_prob_no_alleles;

        unsigned int max_n_alleles_merged;
        netlikelihood::merge_top_of_branch_partials(
                max_num_alleles,
                max_num_alleles,
                top_child1_partials.get_pattern_probability(0,0),
                top_child2_partials.get_pattern_probability(0,0),
                top_ch1_partials,
                top_ch2_partials,
                max_n_alleles_merged,
                merged_probs,
                merged_prob_no_alleles,
                true);
        REQUIRE(max_n_alleles_merged == 4);
        BiallelicPatternProbabilityMatrix merged(max_n_alleles_merged, merged_probs);

        std::cout << "0,0 : " << merged.get_pattern_probability(0,0) * 144 << "\n";
        std::cout << "1,0 : " << merged.get_pattern_probability(1,0) * 144 << "\n";
        std::cout << "1,1 : " << merged.get_pattern_probability(1,1) * 144 << "\n";
        std::cout << "2,0 : " << merged.get_pattern_probability(2,0) * 144 << "\n";
        std::cout << "2,1 : " << merged.get_pattern_probability(2,1) * 144 << "\n";
        std::cout << "2,2 : " << merged.get_pattern_probability(2,2) * 144 << "\n";
        std::cout << "3,0 : " << merged.get_pattern_probability(3,0) * 144 << "\n";
        std::cout << "3,1 : " << merged.get_pattern_probability(3,1) * 144 << "\n";
        std::cout << "3,2 : " << merged.get_pattern_probability(3,2) * 144 << "\n";
        std::cout << "3,3 : " << merged.get_pattern_probability(3,3) * 144 << "\n";
        std::cout << "4,0 : " << merged.get_pattern_probability(4,0) * 144 << "\n";
        std::cout << "4,1 : " << merged.get_pattern_probability(4,1) * 144 << "\n";
        std::cout << "4,2 : " << merged.get_pattern_probability(4,2) * 144 << "\n";
        std::cout << "4,3 : " << merged.get_pattern_probability(4,3) * 144 << "\n";
        std::cout << "4,4 : " << merged.get_pattern_probability(4,4) * 144 << "\n";

        REQUIRE(merged.get_pattern_probability(0,0) == 0.0);
        REQUIRE(merged.get_pattern_probability(1,0) == 0.0);
        REQUIRE(merged.get_pattern_probability(1,1) == 0.0);
        REQUIRE(merged.get_pattern_probability(2,0) == Approx(4/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(2,1) == Approx(5.5/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(2,2) == Approx(6/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(3,0) == Approx(10/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(3,1) == Approx(38.0/3.0/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(3,2) == Approx(43.0/3.0/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(3,3) == Approx(14/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,0) == Approx(4/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,1) == Approx(5/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,2) == Approx(34.0/6.0/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,3) == Approx(5.5/144.0).epsilon(1e-8));
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
        double merged_prob_no_alleles;

        unsigned int max_n_alleles_merged;
        netlikelihood::merge_top_of_branch_partials(
                max_num_alleles,
                max_num_alleles,
                top_child1_partials.get_pattern_probability(0,0),
                top_child2_partials.get_pattern_probability(0,0),
                top_ch1_partials,
                top_ch2_partials,
                max_n_alleles_merged,
                merged_probs,
                merged_prob_no_alleles,
                false);
        REQUIRE(max_n_alleles_merged == 4);
        BiallelicPatternProbabilityMatrix merged(max_n_alleles_merged, merged_probs);
        merged.set_pattern_probability(0, 0, merged_prob_no_alleles);
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
 * Below, the factor under terms is the hypergeometric probability:
 *
 *   ((n1 choose r1) * (n2 choose r2)) / (n choose r)
 *
 * See Eq. 19 in Bryant et al. 2012
 *
 *   Merged
 *   0,0 = (1 * 2)                                  = 2
 *   1,0 = (1 * 2) + (4 * 1)                        = 6
 *   0,1 = (1 * 2) + (3 * 1)                        = 5
 *   2,0 = (2 * 2) + (1 * 1) + (1 * 4)              = 9
 *   1,1 = (3 * 2) + (1 * 1) + (1 * 3) + (1 * 4)    = 10.5
 *                              * 1/2     * 1/2
 *   0,2 = (4 * 2) + (1 * 1) + (1 * 3)              = 12
 *   3,0 = (2 * 4) + (1 * 1)                        = 9
 *   2,1 = (2 * 3) + (1 * 1) + (3 * 4) + (1 * 1)    = 33/3
 *          * 1/3     * 1/3     * 2/3     * 2/3
 *   1,2 = (1 * 1) + (4 * 4) + (3 * 3) + (1 * 1)    = 37/3
 *          * 1/3     * 1/3     * 2/3     * 2/3
 *   0,3 = (4 * 3) + (1 * 1)                        = 13
 *   4,0 = (2 * 1)                                  = 2
 *   3,1 = (2 * 1) + (1 * 3)                        = 2.5
 *          * 1/2     * 1/2
 *   2,2 = (2 * 1) + (4 * 1) + (3 * 1)              = 3
 *          * 1/6     * 1/6     * 4/6
 *   1,3 = (3 * 1) + (1 * 4)                        = 3.5
 *          * 1/2     * 1/2
 *   0,4 = (4 * 1)                                  = 4
 *   TOTAL                                          = 144/144
 */
TEST_CASE("Testing merge_top_of_branch_partials with missing probs and scaling",
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
        double merged_prob_no_alleles;

        unsigned int max_n_alleles_merged;
        netlikelihood::merge_top_of_branch_partials(
                max_num_alleles,
                max_num_alleles,
                top_child1_partials.get_pattern_probability(0,0),
                top_child2_partials.get_pattern_probability(0,0),
                top_ch1_partials,
                top_ch2_partials,
                max_n_alleles_merged,
                merged_probs,
                merged_prob_no_alleles,
                true);
        REQUIRE(max_n_alleles_merged == 4);
        BiallelicPatternProbabilityMatrix merged(max_n_alleles_merged, merged_probs);
        merged.set_pattern_probability(0, 0, merged_prob_no_alleles);
        REQUIRE(merged.get_pattern_probability(0,0) == Approx(2/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(1,0) == Approx(6/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(1,1) == Approx(5/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(2,0) == Approx(9/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(2,1) == Approx(10.5/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(2,2) == Approx(12/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(3,0) == Approx(9/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(3,1) == Approx(11/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(3,2) == Approx(37.0/3.0/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(3,3) == Approx(13/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,0) == Approx(2/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,1) == Approx(2.5/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,2) == Approx(3/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,3) == Approx(3.5/144.0).epsilon(1e-8));
        REQUIRE(merged.get_pattern_probability(4,4) == Approx(4/144.0).epsilon(1e-8));
    }
}
