#include "catch.hpp"
#include "ecoevolity/tree.hpp"
#include "ecoevolity/stats_util.hpp"

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood = -31.77866581319647
TEST_CASE("Testing simple likelihood of PopulationTree", "[PopulationTree]") {

    SECTION("Testing constructor and likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        PopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing threaded likelihood of PopulationTree", "[PopulationTree]") {

    SECTION("Testing constructor and threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        PopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing over-threaded likelihood of PopulationTree", "[PopulationTree]") {

    SECTION("Testing constructor and over-threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        PopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(10);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(10);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood             = -248.93254688526213
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (0.01, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded likelihood (0.01, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7099.716015109998
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.01, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.01, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -6986.120524781545
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.01, 10.0, 1.0, 1.0, dominant)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6986.120524781545));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.01, 10.0, 1.0, 1.0, dominant)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(5);
        REQUIRE(l == Approx(-6986.120524781545));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}



// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -328.39238828878365
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (0.0, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded likelihood (0.0, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(7);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(7);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7256.501742344454
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.0, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.0, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7223.362711937651
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.0, 10.0, 1.0, 1.0, dominant)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7223.362711937651));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.0, 10.0, 1.0, 1.0, dominant)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7223.362711937651));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -227.41048391087554
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (0.2, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded likelihood (0.2, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7304.180743441677
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.2, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.2, 10.0, 1.0, 1.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7405.145951634711
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.2, 10.0, 1.0, 1.0, dominant)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7405.145951634711));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.2, 10.0, 1.0, 1.0, dominant)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(6);
        REQUIRE(l == Approx(-7405.145951634711));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}




// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -327.7437811413033
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -6472.856486972301
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(5);
        REQUIRE(l == Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(3);
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -6494.774924871097
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.03, 10.0, 10.0, 10.0/19.0, dominant)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6494.774924871097));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.03, 10.0, 10.0, 10.0/19.0, dominant)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-6494.774924871097));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
  
  
// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -265.0023534261969
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(1);
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -10163.468886613919
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(5);
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -10999.288193543642
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.03, 10.0, 10.0/19.0, 10.0, dominant)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-10999.288193543642));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.03, 10.0, 10.0/19.0, 10.0, dominant)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-10999.288193543642));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}



// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -224.40177558289847
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-224.40177558289847));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(100000);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-224.40177558289847));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -8158.88094671241
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-8158.88094671241));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(1000000);
        REQUIRE(l == Approx(-8158.88094671241));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -8034.250341980543
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex likelihood (0.03, 111.1, 10.0/19.0, 10.0, dominant)", "[PopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-8034.250341980543));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded likelihood (0.03, 111.1, 10.0/19.0, 10.0, dominant)", "[PopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        PopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-8034.250341980543));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood = -31.77866581319647
TEST_CASE("Testing simple likelihood of ComparisonPopulationTree", "[ComparisonPopulationTree]") {

    SECTION("Testing constructor and likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_root().get_label() == "root-pop1");
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing simple threaded likelihood of ComparisonPopulationTree", "[ComparisonPopulationTree]") {

    SECTION("Testing constructor and threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_root().get_label() == "root-pop1");
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

TEST_CASE("Testing simple prior of PopulationTree", "[PopulationTree]") {

    SECTION("Testing constructor and prior calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        PopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);

        tree.set_root_height(0.1);
        tree.set_all_population_sizes(2.0/100.0);
        tree.set_freq_1(0.95);

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(100.0));
        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));
        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(2.0, 1.0));

        // height   -5.3948298140119091
        // sizes     3 * -155.90663080917298
        // f1       0.64185388617239469 
        // total    -472.47286835535851
        
        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-467.71989242751897));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == Approx(0.64185388617239469));

        double d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-472.47286835535851));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-472.47286835535851));


        tree.store_prior_density();
        tree.constrain_state_frequencies();

        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-467.71989242751897));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == 0.0);

        d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-473.1147222415309));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-473.1147222415309));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-472.47286835535851));


        tree.store_prior_density();
        tree.constrain_population_sizes();

        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-155.90663080917298));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == 0.0);

        d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-161.30146062318488));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-161.30146062318488));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-473.1147222415309));

        tree.restore_prior_density();
        REQUIRE(tree.get_log_prior_density_value() == Approx(-473.1147222415309));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-473.1147222415309));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}


TEST_CASE("Testing hemi129.nex state manipulation", "[ComparisonPopulationTree]") {

    SECTION("Testing state manipulation") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(10.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(1.5, 3.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_height(0.01);
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_child_population_size(0, 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_child_population_size(1, 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1(0.5);
        REQUIRE(tree.is_dirty());

        REQUIRE(tree.get_number_of_likelihood_calculations() == 0);

        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 1);

        tree.store_state();
        tree.set_root_height(0.0);
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 2);

        tree.restore_state();
        // ComparisonPopulationTree does not store/restore heights
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.compute_log_likelihood_and_prior();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.store_state();
        tree.set_root_height(0.0);
        tree.compute_log_likelihood_and_prior();
        REQUIRE(tree.get_number_of_likelihood_calculations() == 4);


        tree.store_state();
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 5);


        tree.store_state();
        tree.set_root_height(0.03);
        tree.set_freq_1(0.05);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.05);
        REQUIRE(tree.get_freq_0() == 0.95);
        REQUIRE(tree.get_u() == Approx(10.0));
        REQUIRE(tree.get_v() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.6991237645509));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.8315389903985));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 6);


        tree.store_state();
        tree.set_freq_1(0.95);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_freq_1() == 0.95);
        REQUIRE(tree.get_freq_0() == Approx(0.05));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1347.1157822333005));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.6991237645509));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 7);


        tree.store_state();
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-47.141114882027253));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1347.1157822333005));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 8);

        tree.store_state();
        tree.constrain_state_frequencies();
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.9800426669165));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-47.141114882027253));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 9);


        tree.store_state();
        tree.fold_patterns();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.9800426669165));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.9800426669165));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 10);


        tree.store_state();
        tree.constrain_population_sizes();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.66001422230551));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.9800426669165));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 11);

        tree.store_state();
        tree.estimate_mutation_rate();
        REQUIRE(! tree.mutation_rate_is_fixed());
        tree.set_mutation_rate(1.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());
        tree.set_mutation_rate_prior(std::make_shared<GammaDistribution>(10.0, 0.1));
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.66001422230551));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 12);

        tree.store_state();
        tree.set_mutation_rate(2.0);
        tree.set_root_height(0.1);
        tree.set_root_population_size(2.0/(20.0 * 2 * tree.get_ploidy()));
        tree.set_child_population_size(0, 2.0/(20.0 * 2 * tree.get_ploidy()));
        tree.set_child_population_size(1, 2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 2.0);
        REQUIRE(tree.get_root_population_size() == 2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-207.43599077244659));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        tree.restore_state();
        // ComparisonPopulationTree does not store/restore node heights
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

TEST_CASE("Testing hemi129.nex state manipulation for PopulationTree", "[PopulationTree]") {

    SECTION("Testing state manipulation") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(10.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(1.5, 3.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_height(0.01);
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1(0.5);
        REQUIRE(tree.is_dirty());

        REQUIRE(tree.get_number_of_likelihood_calculations() == 0);

        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1340.6289538974045));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 1);

        tree.store_state();
        tree.set_root_height(0.0);
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1340.5289538974043));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1340.6289538974045));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 2);

        tree.restore_state();
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1340.6289538974045));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1340.6289538974045));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.compute_log_likelihood_and_prior();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1340.6289538974045));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1340.6289538974045));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.store_state();
        tree.set_root_height(0.0);
        tree.compute_log_likelihood_and_prior();
        REQUIRE(tree.get_number_of_likelihood_calculations() == 4);


        tree.store_state();
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.5289538974043));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1340.5289538974043));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 5);


        tree.store_state();
        tree.set_root_height(0.03);
        tree.set_freq_1(0.05);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.05);
        REQUIRE(tree.get_freq_0() == Approx(0.95));
        REQUIRE(tree.get_u() == Approx(10.0));
        REQUIRE(tree.get_v() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1340.6965386715569));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.5289538974043));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 6);


        tree.store_state();
        tree.set_freq_1(0.95);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.95);
        REQUIRE(tree.get_freq_0() == Approx(0.05));
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1345.1131971403065));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1340.6965386715569));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 7);


        tree.store_state();
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-45.138529789033207));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1345.1131971403065));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 8);

        tree.store_state();
        tree.constrain_state_frequencies();
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.6774575739223));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-45.138529789033207));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 9);


        tree.store_state();
        tree.fold_patterns();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1342.6774575739223));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.6774575739223));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 10);


        tree.store_state();
        tree.constrain_population_sizes();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.35742912931147));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1342.6774575739223));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 11);

        tree.store_state();
        tree.estimate_mutation_rate();
        REQUIRE(! tree.mutation_rate_is_fixed());
        tree.set_mutation_rate(1.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());
        tree.set_mutation_rate_prior(std::make_shared<GammaDistribution>(10.0, 0.1));
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.35742912931147));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 12);

        tree.store_state();
        tree.set_mutation_rate(2.0);
        tree.set_root_height(0.1);
        tree.set_all_population_sizes(2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 2.0);
        REQUIRE(tree.get_root_population_size() == 2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-206.13340567945255));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        tree.restore_state();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 9fb5b3b7817a0bd4a21e3f90132f132cca72ce4e)
// SNAPP v1.3.0 (master 4f3f0f7366798f4fb38b766c15f6426a75ddf71e)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0/19.0
// Log likelihood            = -23.81984255023975
// Log likelihood correction = -6.87935580446044
//
// With constant sites inclucded and m_bUseNonPolymorphic = true
// Log likelihood            = -55.01646493341547
// Log likelihood correction = -6.87935580446044
TEST_CASE("Testing affect of constant sites on likelihood of PopulationTree", "[PopulationTree]") {

    SECTION("Testing haploid-standard-full-constant.nex") {
        std::string nex_path = "data/haploid-standard-full-constant.nex";
        PopulationTree t_included(
                nex_path, // path
                ' ',      // pop name delimiter
                true,     // pop name is prefix
                false,    // genotypes are diploid
                false,    // markers are dominant
                false,    // constant sites removed
                true);    // validate
        std::string nex_path2 = "data/haploid-standard-full-constant-removed.nex";
        PopulationTree t(
                nex_path2, // path
                ' ',       // pop name delimiter
                true,      // pop name is prefix
                false,     // genotypes are diploid
                false,     // markers are dominant
                true,      // constant sites removed
                true);     // validate

        t.set_root_height(0.03);
        t.set_freq_1(0.05);
        t.set_all_population_sizes(2.0/(10.0 * 2 * t.get_ploidy()));

        t_included.set_root_height(0.03);
        t_included.set_freq_1(0.05);
        t_included.set_all_population_sizes(2.0/(10.0 * 2 * t_included.get_ploidy()));

        double l = t.compute_log_likelihood();
        REQUIRE(l == Approx(-23.81984255023975));
        REQUIRE(t.get_likelihood_correction() == Approx(-6.87935580446044));

        double l_included = t_included.compute_log_likelihood();

        REQUIRE(l_included == Approx(-55.01646493341547));

        REQUIRE(t.get_likelihood_correction() == t_included.get_likelihood_correction());
    }
}

TEST_CASE("Testing affect of constant sites on threaded likelihood of PopulationTree", "[PopulationTree]") {

    SECTION("Testing haploid-standard-full-constant.nex") {
        std::string nex_path = "data/haploid-standard-full-constant.nex";
        PopulationTree t_included(
                nex_path, // path
                ' ',      // pop name delimiter
                true,     // pop name is prefix
                false,    // genotypes are diploid
                false,    // markers are dominant
                false,    // constant sites removed
                true);    // validate
        std::string nex_path2 = "data/haploid-standard-full-constant-removed.nex";
        PopulationTree t(
                nex_path2, // path
                ' ',       // pop name delimiter
                true,      // pop name is prefix
                false,     // genotypes are diploid
                false,     // markers are dominant
                true,      // constant sites removed
                true);     // validate

        t.set_root_height(0.03);
        t.set_freq_1(0.05);
        t.set_all_population_sizes(2.0/(10.0 * 2 * t.get_ploidy()));

        t_included.set_root_height(0.03);
        t_included.set_freq_1(0.05);
        t_included.set_all_population_sizes(2.0/(10.0 * 2 * t_included.get_ploidy()));

        double l = t.compute_log_likelihood(4);
        REQUIRE(l == Approx(-23.81984255023975));
        REQUIRE(t.get_likelihood_correction() == Approx(-6.87935580446044));

        double l_included = t_included.compute_log_likelihood(2);

        REQUIRE(l_included == Approx(-55.01646493341547));

        REQUIRE(t.get_likelihood_correction() == t_included.get_likelihood_correction());
    }
}


// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.01
// coalescent_rate = 100.0, 200.0, 500.0
// u = 1.0
// v = 1.0
// Log likelihood            = -226.11914854623677
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (0.01, 100.0, 200.0, 500.0)", "[ComparisonPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_child_population_size(0, 2.0/(100.0 * 2 * tree.get_ploidy()));
        tree.set_child_population_size(1, 2.0/(200.0 * 2 * tree.get_ploidy()));
        tree.set_root_population_size(2.0/(500.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-226.11914854623677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}


// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.00506843962151613554
// coalescent_rate = 2.0 / 0.00018955324120485613
// u = 1.0
// v = 1.0
// Log likelihood            = -277.06960543551577
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (0.00506843962151613554, 2.0 / 0.00018955324120485613)", "[ComparisonPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.00506843962151613554);
        tree.set_all_population_sizes(0.00018955324120485613 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-277.06960543551577));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 9.08323190033687971e-09
// coalescent_rate = 2.0 / 2.47975039926886321e-08
// u = 1.0
// v = 1.0
// Log likelihood            = -221.69627648370943
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (9.08323190033687971e-09, 2.0 / 2.47975039926886321e-08)", "[ComparisonPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(9.08323190033687971e-09);
        tree.set_all_population_sizes(2.47975039926886321e-08 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-221.69627648370943));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.04921319733994759e-08
// coalescent_rate = 2.0 / 2.75977168733651178e-10
// u = 1.0
// v = 1.0
// Log likelihood            = -324.2737564069293
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (1.04921319733994759e-08, 2.0 / 2.75977168733651178e-10)", "[ComparisonPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.04921319733994759e-08);
        tree.set_all_population_sizes(2.75977168733651178e-10 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-324.2737564069293));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.012386610001351e-08
// coalescent_rate = 2.0 / 5.75048645855884647e-30
// u = 1.0
// v = 1.0
// Log likelihood            = -1364.1427530000253
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (1.012386610001351e-08, 2.0 / 5.75048645855884647e-30)", "[ComparisonPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.012386610001351e-08);
        tree.set_all_population_sizes(5.75048645855884647e-30 / 4.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-1364.1427530000253));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.04856228318474786e-08
// coalescent_rate = 2.0 / 4.43934332792563837e-305
// u = 1.0
// v = 1.0
// Log likelihood            = NaN
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex likelihood (1.04856228318474786e-08, 2.0 / 4.43934332792563837e-305)", "[ComparisonPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.04856228318474786e-08);
        tree.set_all_population_sizes(4.43934332792563837e-305 / 4.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(std::isnan(l));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.012386610001351e-08
// coalescent_rate = 2.0 / 5.46641122085615013e-09,
//                   2.0 / 7.39871781998828579e-08,
//                   2.0 / 2.71077053326002069e-13
// u = 1.0
// v = 1.0
// Log likelihood            = -96.34394008351177
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing hemi129.nex weirdness", "[ComparisonPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.012386610001351e-08);
        tree.set_child_population_size(0, 5.46641122085615013e-09 / 4.0);
        tree.set_child_population_size(1, 7.39871781998828579e-08 / 4.0);
        tree.set_root_population_size(2.71077053326002069e-13 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-96.34394008351177));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.036374107244057e-08
// coalescent_rate = 2.0 / 4.57999694763258361e-09,
//                   2.0 / 6.70991782555376588e-08,
//                   2.0 / 1.33514111020266258e-08
// u = 1.0
// v = 1.0
// Log likelihood            = -44.95791900747736
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing hemi129.nex weirdness 2", "[ComparisonPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.036374107244057e-08);
        tree.set_child_population_size(0, 4.57999694763258361e-09 / 4.0);
        tree.set_child_population_size(1, 6.70991782555376588e-08 / 4.0);
        tree.set_root_population_size(1.33514111020266258e-08 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-44.95791900747736));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

TEST_CASE("Testing coalesce_in_branch for 2 lineages and theta of 1.0",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing coalesce_in_branch for 2 lineages and theta of 3.7",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 3.7;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing coalesce_in_branch for 2 lineages and theta of 0.17",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing coalesce_in_branch for 3 lineages and theta of 1.0",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing coalesce_in_branch for 3 lineages and theta of 1.47",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 1.47;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing coalesce_in_branch for 3 lineages and theta of 0.17",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing coalesce_in_branch for 10 lineages and theta of 1.0",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing coalesce_in_branch for 10 lineages and theta of 1.47",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 1.47;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing coalesce_in_branch for 10 lineages and theta of 0.17",
        "[ComparisonPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing scaling of simulate_gene_tree for singleton",
        "[ComparisonPopulationTree]") {

    SECTION("Testing singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        std::string nex_path = "data/singleton-n2.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, false, false);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.set_child_population_size(0, (Ne * mu));

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        unsigned int n;
        double mean;
        double variance;
        double sum_devs;
        double d;
        double d_n;
        double mn;
        double mx;
        double x;
        std::shared_ptr<GeneTreeSimNode> gtree;

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            x = gtree->get_height();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);

        REQUIRE(gtree->is_root());
        REQUIRE(gtree->get_number_of_children() == 2);
        REQUIRE(gtree->get_leaf_node_count() == 2);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.0005));
        REQUIRE(mn >= 0.0);

        tree.set_mutation_rate(mu);
        tree.set_root_population_size(Ne);
        tree.set_child_population_size(0, Ne);
        tree.set_root_height(0.1 / mu);

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            x = gtree->get_height();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);

        REQUIRE(gtree->is_root());
        REQUIRE(gtree->get_number_of_children() == 2);
        REQUIRE(gtree->get_leaf_node_count() == 2);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing scaling of simulate_gene_tree for pair",
        "[ComparisonPopulationTree]") {

    SECTION("Testing pair") {
        unsigned int Ne_root = 100000;
        unsigned int Ne_0 = 200000;
        unsigned int Ne_1 = 10000;
        double mu = 1e-8;
        double theta_root = 4 * Ne_root * mu;
        double theta_0 = 4 * Ne_0 * mu;
        double theta_1 = 4 * Ne_1 * mu;
        double time = 2e7;

        double expected_mean_root = theta_root * (1.0 - (1.0 / 2));
        double expected_variance_root = expected_mean_root * expected_mean_root;
        expected_mean_root += (time * mu);

        double expected_mean_0 = theta_0 * (1.0 - (1.0 / 10));
        double expected_mean_1 = theta_1 * (1.0 - (1.0 / 10));

        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne_root);
        tree.set_child_population_size(0, (Ne_0));
        tree.set_child_population_size(1, (Ne_1));

        tree.set_root_height(time);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(mu);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        std::shared_ptr<GeneTreeSimNode> gtree;
        SampleSummarizer<double> height_root;
        SampleSummarizer<double> height_0;
        SampleSummarizer<double> height_1;
        SampleSummarizer<double> tree_length;

        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            REQUIRE(gtree->is_root());
            REQUIRE(gtree->get_number_of_children() == 2);
            REQUIRE(gtree->get_leaf_node_count() == 20);
            height_root.add_sample(gtree->get_height());
            double h0 = gtree->get_child(0)->get_height();
            double h1 = gtree->get_child(1)->get_height();
            if (h0 < h1) {
                height_0.add_sample(h1);
                height_1.add_sample(h0);
            }
            else {
                height_0.add_sample(h0);
                height_1.add_sample(h1);
            }
            tree_length.add_sample(gtree->get_clade_length());
        }

        REQUIRE(height_root.mean() == Approx(expected_mean_root).epsilon(0.0005));
        REQUIRE(height_root.variance() == Approx(expected_variance_root).epsilon(0.0005));

        REQUIRE(height_0.mean() == Approx(expected_mean_0).epsilon(0.00003));
        REQUIRE(height_1.mean() == Approx(expected_mean_1).epsilon(0.00003));


        tree.set_mutation_rate(mu * 2.0);

        SampleSummarizer<double> scale_height_root;
        SampleSummarizer<double> scale_height_0;
        SampleSummarizer<double> scale_height_1;
        SampleSummarizer<double> scale_tree_length;

        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            gtree->scale(0.5);
            REQUIRE(gtree->is_root());
            REQUIRE(gtree->get_number_of_children() == 2);
            REQUIRE(gtree->get_leaf_node_count() == 20);
            scale_height_root.add_sample(gtree->get_height());
            double h0 = gtree->get_child(0)->get_height();
            double h1 = gtree->get_child(1)->get_height();
            if (h0 < h1) {
                scale_height_0.add_sample(h1);
                scale_height_1.add_sample(h0);
            }
            else {
                scale_height_0.add_sample(h0);
                scale_height_1.add_sample(h1);
            }
            scale_tree_length.add_sample(gtree->get_clade_length());
        }

        REQUIRE(scale_height_root.mean() == Approx(expected_mean_root).epsilon(0.0005));
        REQUIRE(scale_height_root.variance() == Approx(expected_variance_root).epsilon(0.0005));

        REQUIRE(scale_height_0.mean() == Approx(expected_mean_0).epsilon(0.00003));
        REQUIRE(scale_height_1.mean() == Approx(expected_mean_1).epsilon(0.00003));

        REQUIRE(scale_tree_length.mean() == Approx(tree_length.mean()).epsilon(0.0001));
        REQUIRE(scale_tree_length.variance() == Approx(tree_length.variance()).epsilon(0.0001));
    }
}

TEST_CASE("Testing dataset simulation", "[ComparisonPopulationTree]") {
    SECTION("Testing simulate_biallelic_data_set for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        ComparisonPopulationTree tree(nex_path, ' ', true, false, false, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        BiallelicData data = tree.simulate_biallelic_data_set(rng);

        REQUIRE(data.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data.markers_are_dominant() == false);
        REQUIRE(tree.get_data().markers_are_dominant() == false);

        data.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8).epsilon(0.01));
        REQUIRE(data.get_proportion_1() == Approx(0.625).epsilon(0.01));
    }
}

TEST_CASE("Testing singleton acquisition bioas", "[ComparisonPopulationTree]") {
    SECTION("Testing for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        ComparisonPopulationTree tree(nex_path, ' ', true, false, false, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        BiallelicData data1 = tree.simulate_biallelic_data_set(rng, 1.0);
        BiallelicData data0 = tree.simulate_biallelic_data_set(rng, 0.0);

        REQUIRE(data1.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data1.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data1.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data1.markers_are_dominant() == false);
        REQUIRE(data0.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data0.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data0.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data0.markers_are_dominant() == false);
        REQUIRE(tree.get_data().markers_are_dominant() == false);

        unsigned int singleton_count10 = 0;
        unsigned int singleton_count05 = 0;
        unsigned int singleton_count00 = 0;
        RandomNumberGenerator rng10 = RandomNumberGenerator(123);
        RandomNumberGenerator rng05 = RandomNumberGenerator(123);
        RandomNumberGenerator rng00 = RandomNumberGenerator(123);
        for (unsigned int rep = 0; rep < 100; ++rep) {
            BiallelicData data00 = tree.simulate_biallelic_data_set(rng00, 0.0);
            BiallelicData data05 = tree.simulate_biallelic_data_set(rng05, 0.5);
            BiallelicData data10 = tree.simulate_biallelic_data_set(rng10, 1.0);
            for (unsigned int i = 0; i < data00.get_number_of_patterns(); ++i) {
                const std::vector<unsigned int>& red_allele_counts = data00.get_red_allele_counts(i);
                const std::vector<unsigned int>& allele_counts = data00.get_allele_counts(i);
                unsigned int nreds = 0;
                unsigned int nalleles = 0;
                for (unsigned int j = 0; j < allele_counts.size(); ++j) {
                    nreds += red_allele_counts.at(j);
                    nalleles += allele_counts.at(j);
                }
                if ((nreds == 1) || (nreds == (nalleles - 1))) {
                    singleton_count00 += data00.get_pattern_weight(i);
                }
            }
            for (unsigned int i = 0; i < data05.get_number_of_patterns(); ++i) {
                const std::vector<unsigned int>& red_allele_counts = data05.get_red_allele_counts(i);
                const std::vector<unsigned int>& allele_counts = data05.get_allele_counts(i);
                unsigned int nreds = 0;
                unsigned int nalleles = 0;
                for (unsigned int j = 0; j < allele_counts.size(); ++j) {
                    nreds += red_allele_counts.at(j);
                    nalleles += allele_counts.at(j);
                }
                if ((nreds == 1) || (nreds == (nalleles - 1))) {
                    singleton_count05 += data05.get_pattern_weight(i);
                }
            }
            for (unsigned int i = 0; i < data10.get_number_of_patterns(); ++i) {
                const std::vector<unsigned int>& red_allele_counts = data10.get_red_allele_counts(i);
                const std::vector<unsigned int>& allele_counts = data10.get_allele_counts(i);
                unsigned int nreds = 0;
                unsigned int nalleles = 0;
                for (unsigned int j = 0; j < allele_counts.size(); ++j) {
                    nreds += red_allele_counts.at(j);
                    nalleles += allele_counts.at(j);
                }
                if ((nreds == 1) || (nreds == (nalleles - 1))) {
                    singleton_count10 += data10.get_pattern_weight(i);
                }
            }
        }
        REQUIRE(singleton_count10 > 0);
        REQUIRE(singleton_count05 > 0);
        REQUIRE(singleton_count00 == 0);
        REQUIRE((double)singleton_count05 == Approx(singleton_count10 * 0.5).epsilon(0.05));
    }
}

TEST_CASE("Testing scaling of dataset simulation for singleton",
        "[ComparisonPopulationTree]") {

    SECTION("Testing simulate_biallelic_data_set for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        // Two times the expected gene tree root height
        double expected_mean = 2.0 * (theta * (1.0 - (1.0 / nlineages)));

        std::string nex_path = "data/dummy-singleton-n2-100k.nex";
        // Need to keep constant characters
        ComparisonPopulationTree tree(nex_path, '-', true, false, false, false);
        REQUIRE(tree.get_leaf_node_count() == 1);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.set_child_population_size(0, (Ne * mu));

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        BiallelicData data;
        SampleSummarizer<double> divergence;
        unsigned int nsamples = 100;

        for (unsigned int i = 0; i < nsamples; ++i) {
            data = tree.simulate_biallelic_data_set(rng, 1.0, false);
            REQUIRE(data.get_number_of_sites() == 100000);
            double x = (double)data.get_number_of_variable_sites() / (double)data.get_number_of_sites();
            divergence.add_sample(x);
        }

        std::cout << "Expected mean divergence: " << expected_mean << "\n";
        std::cout << "Simulated mean divergence: " << divergence.mean() << "\n";
        REQUIRE(divergence.mean() == Approx(expected_mean).epsilon(0.0001));

        tree.set_mutation_rate(mu);
        tree.set_root_population_size(Ne);
        tree.set_child_population_size(0, Ne);
        tree.set_root_height(0.1 / mu);

        SampleSummarizer<double> divergence2;

        for (unsigned int i = 0; i < nsamples; ++i) {
            data = tree.simulate_biallelic_data_set(rng, 1.0, false);
            REQUIRE(data.get_number_of_sites() == 100000);
            double x = (double)data.get_number_of_variable_sites() / (double)data.get_number_of_sites();
            divergence2.add_sample(x);
        }

        std::cout << "Simulated mean divergence: " << divergence2.mean() << "\n";
        REQUIRE(divergence2.mean() == Approx(expected_mean).epsilon(0.0001));
    }
}

TEST_CASE("Testing scaling of simulation of loci for singleton",
        "[ComparisonPopulationTree]") {

    SECTION("Testing simulate_complete_biallelic_data_set for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        // Two times the expected gene tree root height
        double expected_mean = 2.0 * (theta * (1.0 - (1.0 / nlineages)));

        std::string nex_path = "data/dummy-singleton-n2-100k.nex";
        // Need to keep constant characters
        ComparisonPopulationTree tree(nex_path, '-', true, false, false, false);
        REQUIRE(tree.get_leaf_node_count() == 1);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.set_child_population_size(0, (Ne * mu));

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        SampleSummarizer<double> divergence;
        unsigned int nsamples = 100;

        for (unsigned int i = 0; i < nsamples; ++i) {
            auto data_nloci = tree.simulate_complete_biallelic_data_set(rng, 1000, 1.0, false);
            REQUIRE(data_nloci.second == 100);
            REQUIRE(data_nloci.first.get_number_of_sites() == 100000);
            double x = (double)data_nloci.first.get_number_of_variable_sites() / (double)data_nloci.first.get_number_of_sites();
            divergence.add_sample(x);
        }

        std::cout << "Expected mean divergence: " << expected_mean << "\n";
        std::cout << "Simulated mean divergence: " << divergence.mean() << "\n";
        REQUIRE(divergence.mean() == Approx(expected_mean).epsilon(0.0001));

        tree.set_mutation_rate(mu);
        tree.set_root_population_size(Ne);
        tree.set_child_population_size(0, Ne);
        tree.set_root_height(0.1 / mu);

        SampleSummarizer<double> divergence2;

        for (unsigned int i = 0; i < nsamples; ++i) {
            auto data_nloci = tree.simulate_complete_biallelic_data_set(rng, 1000, 1.0, false);
            REQUIRE(data_nloci.second == 100);
            REQUIRE(data_nloci.first.get_number_of_sites() == 100000);
            double x = (double)data_nloci.first.get_number_of_variable_sites() / (double)data_nloci.first.get_number_of_sites();
            divergence2.add_sample(x);
        }

        std::cout << "Simulated mean divergence: " << divergence2.mean() << "\n";
        REQUIRE(divergence2.mean() == Approx(expected_mean).epsilon(0.0001));
    }
}

TEST_CASE("Testing scaling of simulation of one variable site per locus for singleton",
        "[ComparisonPopulationTree]") {

    SECTION("Testing simulate_data_set_max_one_variable_site_per_locus for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        // Two times the expected gene tree root height
        double expected_mean = 2.0 * (theta * (1.0 - (1.0 / nlineages)));
        double epsilon = 0.0001;

        std::string nex_path = "data/dummy-singleton-n2-100k.nex";
        // Need to keep constant characters
        ComparisonPopulationTree tree(nex_path, '-', true, false, false, false);
        REQUIRE(tree.get_leaf_node_count() == 1);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.set_child_population_size(0, (Ne * mu));

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        /* RandomNumberGenerator rng = RandomNumberGenerator(54321); */
        RandomNumberGenerator rng = RandomNumberGenerator(12345678);

        SampleSummarizer<double> divergence;
        unsigned int nsamples = 100;

        for (unsigned int i = 0; i < nsamples; ++i) {
            auto data_nloci = tree.simulate_data_set_max_one_variable_site_per_locus(rng, 1000, 1.0, false);
            REQUIRE(data_nloci.second == 100);
            REQUIRE(data_nloci.first.get_number_of_sites() <= 100);
            double x = (double)data_nloci.first.get_number_of_variable_sites() / (data_nloci.first.get_number_of_sites() + data_nloci.first.get_number_of_constant_sites_removed());
            divergence.add_sample(x);
        }

        std::cout << "Expected mean divergence: " << expected_mean << "\n";
        std::cout << "Simulated mean divergence: " << divergence.mean() << "\n";
        // By stopping the simulation of sites for each locus when we get a
        // variable site, we are under-sampling sites from loci with older
        // coalescence times (over-sampling sites from loci with younger
        // coalescence times), so we should under estimate diversity:
        REQUIRE(divergence.mean() < (expected_mean - epsilon));

        tree.set_mutation_rate(mu);
        tree.set_root_population_size(Ne);
        tree.set_child_population_size(0, Ne);
        tree.set_root_height(0.1 / mu);

        SampleSummarizer<double> divergence2;

        for (unsigned int i = 0; i < nsamples; ++i) {
            auto data_nloci = tree.simulate_data_set_max_one_variable_site_per_locus(rng, 1000, 1.0, false);
            REQUIRE(data_nloci.second == 100);
            REQUIRE(data_nloci.first.get_number_of_sites() <= 100);
            double x = (double)data_nloci.first.get_number_of_variable_sites() / (data_nloci.first.get_number_of_sites() + data_nloci.first.get_number_of_constant_sites_removed());
            divergence2.add_sample(x);
        }

        std::cout << "Simulated mean divergence: " << divergence2.mean() << "\n";
        // By stopping the simulation of sites for each locus when we get a
        // variable site, we are under-sampling sites from loci with older
        // coalescence times (over-sampling sites from loci with younger
        // coalescence times), so we should under estimate diversity:
        REQUIRE(divergence2.mean() < (expected_mean - epsilon));
    }
}

TEST_CASE("Testing draw_from_prior for fully fixed", "[ComparisonPopulationTree]") {

    SECTION("Testing draw_from_prior for fully fixed pair") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.001);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_root_population_size() == 0.001);
            REQUIRE(tree.get_child_population_size(0) == 0.001);
            REQUIRE(tree.get_child_population_size(1) == 0.001);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);
        }
    }
}

TEST_CASE("Testing draw_from_prior for constrained sizes", "[ComparisonPopulationTree]") {

    SECTION("Testing draw_from_prior for constrained sizes") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        SampleSummarizer<double> pop_size;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() == tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() == tree.get_child_population_size(1));

            pop_size.add_sample(tree.get_root_population_size());
        }

        REQUIRE(pop_size.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
    }
}

TEST_CASE("Testing draw_from_prior for unconstrained sizes", "[ComparisonPopulationTree]") {

    SECTION("Testing draw_from_prior for unconstrained sizes") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.set_all_population_sizes(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            pop_size_root.add_sample(tree.get_root_population_size());
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
        }

        REQUIRE(pop_size_root.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_root.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_0.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_0.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_1.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_1.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
    }
}

TEST_CASE("Testing draw_from_prior for fully parameterized", "[ComparisonPopulationTree]") {

    SECTION("Testing draw_from_prior for fully parameterized") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(0.5, 0.8);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.set_all_population_sizes(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.estimate_state_frequencies();
        tree.set_freq_1(0.5);

        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;
        SampleSummarizer<double> rate;
        SampleSummarizer<double> f;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            pop_size_root.add_sample(tree.get_root_population_size());
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
            rate.add_sample(tree.get_mutation_rate());
            f.add_sample(tree.get_freq_1());
        }

        REQUIRE(pop_size_root.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_root.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_0.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_0.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_1.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_1.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(f.mean() == Approx(f_prior->get_mean()).epsilon(0.01));
        REQUIRE(f.variance() == Approx(f_prior->get_variance()).epsilon(0.01));
        REQUIRE(rate.mean() == Approx(rate_prior->get_mean()).epsilon(0.1));
        REQUIRE(rate.variance() == Approx(rate_prior->get_variance()).epsilon(0.1));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood = -31.77866581319647
TEST_CASE("Testing simple likelihood of DirichletPopulationTree", "[DirichletPopulationTree]") {

    SECTION("Testing constructor and likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing threaded likelihood of DirichletPopulationTree", "[DirichletPopulationTree]") {

    SECTION("Testing constructor and threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing over-threaded likelihood of DirichletPopulationTree", "[DirichletPopulationTree]") {

    SECTION("Testing constructor and over-threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(10);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(10);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood             = -248.93254688526213
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing hemi129.nex DPT likelihood (0.01, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex DPT threaded likelihood (0.01, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7099.716015109998
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.01, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.01, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -6986.120524781545
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.01, 10.0, 1.0, 1.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6986.120524781545));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.01, 10.0, 1.0, 1.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(5);
        REQUIRE(l == Approx(-6986.120524781545));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}



// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -328.39238828878365
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex DPT likelihood (0.0, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex DPT threaded likelihood (0.0, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(7);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(7);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7256.501742344454
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.0, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.0, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7223.362711937651
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.0, 10.0, 1.0, 1.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7223.362711937651));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.0, 10.0, 1.0, 1.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7223.362711937651));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -227.41048391087554
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex DPT likelihood (0.2, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex DPT threaded likelihood (0.2, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7304.180743441677
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.2, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.2, 10.0, 1.0, 1.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7405.145951634711
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.2, 10.0, 1.0, 1.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7405.145951634711));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.2, 10.0, 1.0, 1.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(6);
        REQUIRE(l == Approx(-7405.145951634711));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}




// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -327.7437811413033
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex DPT likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex DPT threaded likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -6472.856486972301
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(5);
        REQUIRE(l == Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(3);
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -6494.774924871097
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.03, 10.0, 10.0, 10.0/19.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6494.774924871097));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.03, 10.0, 10.0, 10.0/19.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-6494.774924871097));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
  
  
// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -265.0023534261969
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex DPT likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex DPT threaded likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(1);
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -10163.468886613919
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(5);
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -10999.288193543642
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.03, 10.0, 10.0/19.0, 10.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-10999.288193543642));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.03, 10.0, 10.0/19.0, 10.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-10999.288193543642));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}



// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -224.40177558289847
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex DPT likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_mean_population_size(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-224.40177558289847));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex DPT threaded likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_mean_population_size(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(100000);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-224.40177558289847));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -8158.88094671241
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_mean_population_size(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-8158.88094671241));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_mean_population_size(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(1000000);
        REQUIRE(l == Approx(-8158.88094671241));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -8034.250341980543
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex DPT likelihood (0.03, 111.1, 10.0/19.0, 10.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_mean_population_size(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-8034.250341980543));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex DPT threaded likelihood (0.03, 111.1, 10.0/19.0, 10.0, dominant)", "[DirichletPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_mean_population_size(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-8034.250341980543));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

TEST_CASE("Testing simple prior of DirichletPopulationTree", "[DirichletPopulationTree]") {

    SECTION("Testing constructor and prior calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);

        tree.set_root_height(0.1);
        tree.set_mean_population_size(2.0/100.0);
        tree.set_freq_1(0.95);

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(100.0));
        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));
        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(2.0, 1.0));

        // height           -5.3948298140119091
        // size             -155.90663080917298
        // size multipliers 0.69314718055994529 (relative is 0.0)
        // f1               0.64185388617239469 
        // total            -160.6596067370125
        
        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-155.90663080917298));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == Approx(0.64185388617239469));

        double d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-160.6596067370125));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-160.6596067370125));


        tree.store_prior_density();
        tree.constrain_state_frequencies();

        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-155.90663080917298));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == 0.0);

        d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-161.30146062318488));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-161.30146062318488));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-160.6596067370125));

        tree.restore_prior_density();
        REQUIRE(tree.get_log_prior_density_value() == Approx(-160.6596067370125));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-160.6596067370125));
        REQUIRE(tree.get_degree_of_root() == 2);

        tree.fix_population_size_multipliers();

        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-155.90663080917298));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == 0.0);

        d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-161.30146062318488));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-161.30146062318488));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-160.6596067370125));

        tree.store_prior_density();

        REQUIRE(tree.get_log_prior_density_value() == Approx(-161.30146062318488));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-161.30146062318488));
    }
}

TEST_CASE("Testing hemi129.nex state manipulation for DirichletPopulationTree", "[DirichletPopulationTree]") {

    SECTION("Testing state manipulation") {
        std::string nex_path = "data/hemi129.nex";
        DirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(10.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(1.5, 3.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        std::vector<double> parameters {10.0, 10.0, 10.0};
        tree.set_population_size_multiplier_prior(std::make_shared<DirichletDistribution>(parameters));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_height(0.01);
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1(0.5);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        std::vector<double> proportions {1.0/3.0, 1.0/3.0, 1.0/3.0};
        tree.set_population_sizes_as_proportions(proportions);

        REQUIRE(tree.is_dirty());

        REQUIRE(tree.get_number_of_likelihood_calculations() == 0);

        tree.compute_log_likelihood_and_prior();

        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-474.97145724683253));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 1);

        tree.store_state();
        tree.set_root_height(0.0);
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-474.87145724683256));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-474.97145724683253));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 2);

        tree.restore_state();
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-474.97145724683253));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-474.97145724683253));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.compute_log_likelihood_and_prior();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-474.97145724683253));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-474.97145724683253));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.store_state();
        tree.set_root_height(0.0);
        tree.compute_log_likelihood_and_prior();
        REQUIRE(tree.get_number_of_likelihood_calculations() == 4);


        tree.store_state();
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-476.87145724683256));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-474.87145724683256));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 5);


        tree.store_state();
        tree.set_root_height(0.03);
        tree.set_freq_1(0.05);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.05);
        REQUIRE(tree.get_freq_0() == Approx(0.95));
        REQUIRE(tree.get_u() == Approx(10.0));
        REQUIRE(tree.get_v() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-475.03904202098477));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-476.87145724683256));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 6);


        tree.store_state();
        tree.set_freq_1(0.95);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.95);
        REQUIRE(tree.get_freq_0() == Approx(0.05));
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-479.45570048973445));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-475.03904202098477));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 7);


        tree.store_state();
        tree.set_mean_population_size(2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(111.1 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-46.130811372643343));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-479.45570048973445));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 8);

        tree.store_state();
        tree.constrain_state_frequencies();
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-477.01996092335042));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-46.130811372643343));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 9);


        tree.store_state();
        tree.fold_patterns();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-477.01996092335042));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-477.01996092335042));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 10);


        tree.store_state();
        tree.fix_population_size_multipliers();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.35742912931147));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-477.01996092335042));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 11);

        tree.store_state();
        tree.estimate_mutation_rate();
        REQUIRE(! tree.mutation_rate_is_fixed());
        tree.set_mutation_rate(1.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());
        tree.set_mutation_rate_prior(std::make_shared<GammaDistribution>(10.0, 0.1));
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.35742912931147));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 12);

        tree.store_state();
        tree.set_mutation_rate(2.0);
        tree.set_root_height(0.1);
        tree.set_mean_population_size(2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 2.0);
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(20.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-206.13340567945255));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        tree.restore_state();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 9fb5b3b7817a0bd4a21e3f90132f132cca72ce4e)
// SNAPP v1.3.0 (master 4f3f0f7366798f4fb38b766c15f6426a75ddf71e)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0/19.0
// Log likelihood            = -23.81984255023975
// Log likelihood correction = -6.87935580446044
//
// With constant sites inclucded and m_bUseNonPolymorphic = true
// Log likelihood            = -55.01646493341547
// Log likelihood correction = -6.87935580446044
TEST_CASE("Testing affect of constant sites on likelihood of DirichletPopulationTree", "[DirichletPopulationTree]") {

    SECTION("Testing haploid-standard-full-constant.nex") {
        std::string nex_path = "data/haploid-standard-full-constant.nex";
        DirichletPopulationTree t_included(
                nex_path, // path
                ' ',      // pop name delimiter
                true,     // pop name is prefix
                false,    // genotypes are diploid
                false,    // markers are dominant
                false,    // constant sites removed
                true);    // validate
        std::string nex_path2 = "data/haploid-standard-full-constant-removed.nex";
        DirichletPopulationTree t(
                nex_path2, // path
                ' ',       // pop name delimiter
                true,      // pop name is prefix
                false,     // genotypes are diploid
                false,     // markers are dominant
                true,      // constant sites removed
                true);     // validate

        t.set_root_height(0.03);
        t.set_freq_1(0.05);
        t.set_mean_population_size(2.0/(10.0 * 2 * t.get_ploidy()));

        t_included.set_root_height(0.03);
        t_included.set_freq_1(0.05);
        t_included.set_mean_population_size(2.0/(10.0 * 2 * t_included.get_ploidy()));

        double l = t.compute_log_likelihood();
        REQUIRE(l == Approx(-23.81984255023975));
        REQUIRE(t.get_likelihood_correction() == Approx(-6.87935580446044));

        double l_included = t_included.compute_log_likelihood();

        REQUIRE(l_included == Approx(-55.01646493341547));

        REQUIRE(t.get_likelihood_correction() == t_included.get_likelihood_correction());
    }
}

TEST_CASE("Testing affect of constant sites on threaded likelihood of DirichletPopulationTree", "[DirichletPopulationTree]") {

    SECTION("Testing haploid-standard-full-constant.nex") {
        std::string nex_path = "data/haploid-standard-full-constant.nex";
        DirichletPopulationTree t_included(
                nex_path, // path
                ' ',      // pop name delimiter
                true,     // pop name is prefix
                false,    // genotypes are diploid
                false,    // markers are dominant
                false,    // constant sites removed
                true);    // validate
        std::string nex_path2 = "data/haploid-standard-full-constant-removed.nex";
        DirichletPopulationTree t(
                nex_path2, // path
                ' ',       // pop name delimiter
                true,      // pop name is prefix
                false,     // genotypes are diploid
                false,     // markers are dominant
                true,      // constant sites removed
                true);     // validate

        t.set_root_height(0.03);
        t.set_freq_1(0.05);
        t.set_mean_population_size(2.0/(10.0 * 2 * t.get_ploidy()));

        t_included.set_root_height(0.03);
        t_included.set_freq_1(0.05);
        t_included.set_mean_population_size(2.0/(10.0 * 2 * t_included.get_ploidy()));

        double l = t.compute_log_likelihood(4);
        REQUIRE(l == Approx(-23.81984255023975));
        REQUIRE(t.get_likelihood_correction() == Approx(-6.87935580446044));

        double l_included = t_included.compute_log_likelihood(2);

        REQUIRE(l_included == Approx(-55.01646493341547));

        REQUIRE(t.get_likelihood_correction() == t_included.get_likelihood_correction());
    }
}








// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood = -31.77866581319647
TEST_CASE("Testing simple likelihood of ComparisonDirichletPopulationTree", "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing constructor and likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_root().get_label() == "root-pop1");
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing simple threaded likelihood of ComparisonDirichletPopulationTree", "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing constructor and threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_root().get_label() == "root-pop1");
        tree.set_root_height(0.01);
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree state manipulation",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing state manipulation") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(10.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        std::vector<double> alphas {10.0, 20.0, 30.0};
        tree.set_population_size_multiplier_prior(std::make_shared<DirichletDistribution>(alphas));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(1.5, 3.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_height(0.01);
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        std::vector<double> m = {1.0/3.0, 1.0/3.0, 1.0/3.0};
        tree.set_population_sizes_as_proportions(m);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1(0.5);
        REQUIRE(tree.is_dirty());

        REQUIRE(tree.get_number_of_likelihood_calculations() == 0);

        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 1);

        tree.store_state();
        tree.set_root_height(0.0);
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 2);

        tree.restore_state();
        // ComparisonDirichletPopulationTree does not store/restore heights
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.compute_log_likelihood_and_prior();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.store_state();
        tree.set_root_height(0.0);
        tree.compute_log_likelihood_and_prior();
        REQUIRE(tree.get_number_of_likelihood_calculations() == 4);


        tree.store_state();
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 5);


        tree.store_state();
        tree.set_root_height(0.03);
        tree.set_freq_1(0.05);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.05);
        REQUIRE(tree.get_freq_0() == 0.95);
        REQUIRE(tree.get_u() == Approx(10.0));
        REQUIRE(tree.get_v() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-509.99999577402207));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-510.13241099986988));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 6);


        tree.store_state();
        tree.set_freq_1(0.95);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_freq_1() == 0.95);
        REQUIRE(tree.get_freq_0() == Approx(0.05));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-514.41665424277176));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-509.99999577402207));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 7);


        tree.store_state();
        tree.set_mean_population_size(2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(111.1 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(111.1 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(111.1 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-81.091765125680695));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-514.41665424277176));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 8);

        tree.store_state();
        tree.constrain_state_frequencies();
        tree.set_mean_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-510.28091467638774));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-81.091765125680695));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 9);


        tree.store_state();
        tree.fold_patterns();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-510.28091467638774));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-510.28091467638774));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 10);


        tree.store_state();
        tree.fix_population_size_multipliers();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.66001422230551));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-510.28091467638774));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 11);

        tree.store_state();
        tree.estimate_mutation_rate();
        REQUIRE(! tree.mutation_rate_is_fixed());
        tree.set_mutation_rate(1.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());
        tree.set_mutation_rate_prior(std::make_shared<GammaDistribution>(10.0, 0.1));
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.66001422230551));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 12);

        tree.store_state();
        tree.set_mutation_rate(2.0);
        tree.set_root_height(0.1);
        tree.set_mean_population_size(2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 2.0);
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(20.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(20.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(20.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-207.43599077244659));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        tree.restore_state();
        // ComparisonDirichletPopulationTree does not store/restore node heights
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        REQUIRE(tree.get_population_sizes_as_proportions() == m);

        tree.store_state();
        tree.estimate_population_size_multipliers();
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        m = {0.5/3.0, 1.0/3.0, 1.5/3.0};
        tree.set_population_sizes_as_proportions(m);
        REQUIRE(tree.get_population_sizes_as_proportions() == m);
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_mean_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_root_population_size() ==   Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(2) * 3.0));
        REQUIRE(tree.get_child_population_size(0) == Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(0) * 3.0));
        REQUIRE(tree.get_child_population_size(1) == Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(1) * 3.0));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(tree.get_log_prior_density_value() == Approx(-504.53672771643153));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 14);

        tree.store_state();
        tree.fix_mean_population_size();
        REQUIRE(tree.get_population_sizes_as_proportions() == m);
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_mean_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_root_population_size() == Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(2) * 3.0));
        REQUIRE(tree.get_child_population_size(0) == Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(0) * 3.0));
        REQUIRE(tree.get_child_population_size(1) == Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(1) * 3.0));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(tree.get_log_prior_density_value() == Approx(-56.876713494126008));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-504.53672771643153));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 15);

        tree.store_state();
        tree.fix_population_size_multipliers();
        REQUIRE(tree.get_population_sizes_as_proportions() == m);
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_mean_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_root_population_size() == Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(2) * 3.0));
        REQUIRE(tree.get_child_population_size(0) == Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(0) * 3.0));
        REQUIRE(tree.get_child_population_size(1) == Approx((2.0/(10.0 * 2 * tree.get_ploidy())) * m.at(1) * 3.0));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(tree.get_log_prior_density_value() == Approx(0.22402344985899036));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-56.876713494126008));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 16);

        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.01
// coalescent_rate = 100.0, 200.0, 500.0
// u = 1.0
// v = 1.0
// Log likelihood            = -226.11914854623677
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing ComparisonDirichletPopulationTree likelihood (0.01, 100.0, 200.0, 500.0)",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        // tree.set_child_population_size(0, 2.0/(100.0 * 2 * tree.get_ploidy()));
        // tree.set_child_population_size(1, 2.0/(200.0 * 2 * tree.get_ploidy()));
        // tree.set_root_population_size(2.0/(500.0 * 2 * tree.get_ploidy()));
        tree.set_mean_population_size(0.0085 / 3.0);
        std::vector<double> proportions = {
                (0.005 / 0.0085), 
                (0.0025 / 0.0085), 
                (0.001 / 0.0085)}; 
        tree.set_population_sizes_as_proportions(proportions);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-226.11914854623677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}


// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.00506843962151613554
// coalescent_rate = 2.0 / 0.00018955324120485613
// u = 1.0
// v = 1.0
// Log likelihood            = -277.06960543551577
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonDirichletPopulationTree likelihood (0.00506843962151613554, 2.0 / 0.00018955324120485613)",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.00506843962151613554);
        tree.set_mean_population_size(0.00018955324120485613 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-277.06960543551577));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 9.08323190033687971e-09
// coalescent_rate = 2.0 / 2.47975039926886321e-08
// u = 1.0
// v = 1.0
// Log likelihood            = -221.69627648370943
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonDirichletPopulationTree likelihood (9.08323190033687971e-09, 2.0 / 2.47975039926886321e-08)", "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(9.08323190033687971e-09);
        tree.set_mean_population_size(2.47975039926886321e-08 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-221.69627648370943));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.04921319733994759e-08
// coalescent_rate = 2.0 / 2.75977168733651178e-10
// u = 1.0
// v = 1.0
// Log likelihood            = -324.2737564069293
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonDirichletPopulationTree likelihood (1.04921319733994759e-08, 2.0 / 2.75977168733651178e-10)", "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.04921319733994759e-08);
        tree.set_mean_population_size(2.75977168733651178e-10 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-324.2737564069293));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.012386610001351e-08
// coalescent_rate = 2.0 / 5.75048645855884647e-30
// u = 1.0
// v = 1.0
// Log likelihood            = -1364.1427530000253
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonDirichletPopulationTree likelihood (1.012386610001351e-08, 2.0 / 5.75048645855884647e-30)", "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.012386610001351e-08);
        tree.set_mean_population_size(5.75048645855884647e-30 / 4.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-1364.1427530000253));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.04856228318474786e-08
// coalescent_rate = 2.0 / 4.43934332792563837e-305
// u = 1.0
// v = 1.0
// Log likelihood            = NaN
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonDirichletPopulationTree likelihood (1.04856228318474786e-08, 2.0 / 4.43934332792563837e-305)", "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.04856228318474786e-08);
        tree.set_mean_population_size(4.43934332792563837e-305 / 4.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(std::isnan(l));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.012386610001351e-08
// coalescent_rate = 2.0 / 5.46641122085615013e-09,
//                   2.0 / 7.39871781998828579e-08,
//                   2.0 / 2.71077053326002069e-13
// u = 1.0
// v = 1.0
// Log likelihood            = -96.34394008351177
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing ComparisonDirichletPopulationTree weirdness", "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.012386610001351e-08);
        tree.set_child_population_size(0, 5.46641122085615013e-09 / 4.0);
        tree.set_child_population_size(1, 7.39871781998828579e-08 / 4.0);
        tree.set_root_population_size(2.71077053326002069e-13 / 4.0);
        /* tree.set_mean_population_size(6.621155041482694e-09); */
        /* std::vector<double> proportions = {0.20639945699081685, 2.7935903077461655, 1.0235263017844203e-05}; */
        /* for (unsigned int i = 0; i < proportions.size(); ++i) { */
        /*     proportions.at(i) /= 3.0; */
        /* } */
        /* tree.set_population_sizes_as_proportions(proportions); */
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-96.34394008351177));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 2 lineages and theta of 1.0",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 2 lineages and theta of 3.7",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 3.7;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 2 lineages and theta of 0.17",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 3 lineages and theta of 1.0",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 3 lineages and theta of 1.47",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 1.47;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 3 lineages and theta of 0.17",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 10 lineages and theta of 1.0",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 10 lineages and theta of 1.47",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 1.47;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.coalesce_in_branch for 10 lineages and theta of 0.17",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonDirichletPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing scaling of ComparisonDirichletPopulationTree.simulate_gene_tree for singleton",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        std::string nex_path = "data/singleton-n2.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, false, false);
        tree.estimate_mutation_rate();

        tree.set_mean_population_size(Ne * mu);

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        unsigned int n;
        double mean;
        double variance;
        double sum_devs;
        double d;
        double d_n;
        double mn;
        double mx;
        double x;
        std::shared_ptr<GeneTreeSimNode> gtree;

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            x = gtree->get_height();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);

        REQUIRE(gtree->is_root());
        REQUIRE(gtree->get_number_of_children() == 2);
        REQUIRE(gtree->get_leaf_node_count() == 2);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.0005));
        REQUIRE(mn >= 0.0);

        tree.set_mutation_rate(mu);
        tree.set_mean_population_size(Ne);
        tree.set_root_height(0.1 / mu);

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            x = gtree->get_height();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);

        REQUIRE(gtree->is_root());
        REQUIRE(gtree->get_number_of_children() == 2);
        REQUIRE(gtree->get_leaf_node_count() == 2);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing scaling of ComparisonDirichletPopulationTree.simulate_gene_tree for pair",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing pair") {
        unsigned int Ne_root = 100000;
        unsigned int Ne_0 = 200000;
        unsigned int Ne_1 = 10000;
        double mu = 1e-8;
        double theta_root = 4 * Ne_root * mu;
        double theta_0 = 4 * Ne_0 * mu;
        double theta_1 = 4 * Ne_1 * mu;
        double time = 2e7;

        double expected_mean_root = theta_root * (1.0 - (1.0 / 2));
        double expected_variance_root = expected_mean_root * expected_mean_root;
        expected_mean_root += (time * mu);

        double expected_mean_0 = theta_0 * (1.0 - (1.0 / 10));
        double expected_mean_1 = theta_1 * (1.0 - (1.0 / 10));

        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);
        tree.estimate_mutation_rate();

        double sum_Ne = Ne_0 + Ne_1 + Ne_root;
        double ave_Ne =  sum_Ne / 3.0;
        std::vector<double> proportions = {
                (Ne_0 / sum_Ne),
                (Ne_1 / sum_Ne),
                (Ne_root / sum_Ne)};
        tree.set_mean_population_size(ave_Ne);
        tree.set_population_sizes_as_proportions(proportions);

        // tree.set_root_population_size(Ne_root);
        // tree.set_child_population_size(0, (Ne_0));
        // tree.set_child_population_size(1, (Ne_1));

        tree.set_root_height(time);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(mu);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        std::shared_ptr<GeneTreeSimNode> gtree;
        SampleSummarizer<double> height_root;
        SampleSummarizer<double> height_0;
        SampleSummarizer<double> height_1;
        SampleSummarizer<double> tree_length;

        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            REQUIRE(gtree->is_root());
            REQUIRE(gtree->get_number_of_children() == 2);
            REQUIRE(gtree->get_leaf_node_count() == 20);
            height_root.add_sample(gtree->get_height());
            double h0 = gtree->get_child(0)->get_height();
            double h1 = gtree->get_child(1)->get_height();
            if (h0 < h1) {
                height_0.add_sample(h1);
                height_1.add_sample(h0);
            }
            else {
                height_0.add_sample(h0);
                height_1.add_sample(h1);
            }
            tree_length.add_sample(gtree->get_clade_length());
        }

        REQUIRE(height_root.mean() == Approx(expected_mean_root).epsilon(0.0005));
        REQUIRE(height_root.variance() == Approx(expected_variance_root).epsilon(0.0005));

        REQUIRE(height_0.mean() == Approx(expected_mean_0).epsilon(0.00003));
        REQUIRE(height_1.mean() == Approx(expected_mean_1).epsilon(0.00003));


        tree.set_mutation_rate(mu * 2.0);

        SampleSummarizer<double> scale_height_root;
        SampleSummarizer<double> scale_height_0;
        SampleSummarizer<double> scale_height_1;
        SampleSummarizer<double> scale_tree_length;

        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            gtree->scale(0.5);
            REQUIRE(gtree->is_root());
            REQUIRE(gtree->get_number_of_children() == 2);
            REQUIRE(gtree->get_leaf_node_count() == 20);
            scale_height_root.add_sample(gtree->get_height());
            double h0 = gtree->get_child(0)->get_height();
            double h1 = gtree->get_child(1)->get_height();
            if (h0 < h1) {
                scale_height_0.add_sample(h1);
                scale_height_1.add_sample(h0);
            }
            else {
                scale_height_0.add_sample(h0);
                scale_height_1.add_sample(h1);
            }
            scale_tree_length.add_sample(gtree->get_clade_length());
        }

        REQUIRE(scale_height_root.mean() == Approx(expected_mean_root).epsilon(0.0005));
        REQUIRE(scale_height_root.variance() == Approx(expected_variance_root).epsilon(0.0005));

        REQUIRE(scale_height_0.mean() == Approx(expected_mean_0).epsilon(0.00003));
        REQUIRE(scale_height_1.mean() == Approx(expected_mean_1).epsilon(0.00003));

        REQUIRE(scale_tree_length.mean() == Approx(tree_length.mean()).epsilon(0.0001));
        REQUIRE(scale_tree_length.variance() == Approx(tree_length.variance()).epsilon(0.0001));
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree dataset simulation",
        "[ComparisonDirichletPopulationTree]") {
    SECTION("Testing simulate_biallelic_data_set for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, false, false, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.fix_population_size_multipliers();
        tree.set_mean_population_size(0.005);
        tree.fix_mean_population_size();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        BiallelicData data = tree.simulate_biallelic_data_set(rng);

        REQUIRE(data.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data.markers_are_dominant() == false);
        REQUIRE(tree.get_data().markers_are_dominant() == false);

        data.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8).epsilon(0.01));
        REQUIRE(data.get_proportion_1() == Approx(0.625).epsilon(0.01));
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.draw_from_prior for fully fixed",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing draw_from_prior for fully fixed pair") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.fix_population_size_multipliers();
        tree.set_mean_population_size(0.001);
        tree.fix_mean_population_size();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_root_population_size() == 0.001);
            REQUIRE(tree.get_child_population_size(0) == 0.001);
            REQUIRE(tree.get_child_population_size(1) == 0.001);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);
        }
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.draw_from_prior for constrained sizes",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing draw_from_prior for constrained sizes") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.fix_population_size_multipliers();
        tree.set_mean_population_size(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        SampleSummarizer<double> pop_size;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() == tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() == tree.get_child_population_size(1));

            pop_size.add_sample(tree.get_root_population_size());
        }

        REQUIRE(pop_size.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.draw_from_prior for unconstrained sizes",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing draw_from_prior for unconstrained sizes") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);
        std::vector<double> alphas {2.0, 1.0, 4.0};
        std::shared_ptr<DirichletDistribution> multiplier_prior = std::make_shared<DirichletDistribution>(alphas);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);
        tree.set_population_size_multiplier_prior(multiplier_prior);

        tree.set_root_height(0.1);
        tree.set_mean_population_size(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        SampleSummarizer<double> pop_size;
        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 100000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            pop_size.add_sample(tree.get_mean_population_size());
            std::vector<double> multipliers = tree.get_population_sizes_as_multipliers();
            pop_size_root.add_sample(multipliers.at(2));
            pop_size_0.add_sample(multipliers.at(0));
            pop_size_1.add_sample(multipliers.at(1));
        }

        REQUIRE(pop_size.mean() == Approx(size_prior->get_mean()).epsilon(0.01));
        REQUIRE(pop_size.variance() == Approx(size_prior->get_variance()).epsilon(0.01));
        REQUIRE(pop_size_root.mean() == Approx(multiplier_prior->get_mean(2) * 3.0).epsilon(0.01));
        REQUIRE(pop_size_root.variance() == Approx(multiplier_prior->get_variance(2) * 9.0).epsilon(0.01));
        REQUIRE(pop_size_0.mean() == Approx(multiplier_prior->get_mean(0) * 3.0).epsilon(0.01));
        REQUIRE(pop_size_0.variance() == Approx(multiplier_prior->get_variance(0) * 9.0).epsilon(0.01));
        REQUIRE(pop_size_1.mean() == Approx(multiplier_prior->get_mean(1) * 3.0).epsilon(0.01));
        REQUIRE(pop_size_1.variance() == Approx(multiplier_prior->get_variance(1) * 9.0).epsilon(0.01));
    }
}

TEST_CASE("Testing ComparisonDirichletPopulationTree.draw_from_prior for fully parameterized",
        "[ComparisonDirichletPopulationTree]") {

    SECTION("Testing draw_from_prior for fully parameterized") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonDirichletPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(0.5, 0.8);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);
        std::vector<double> alphas {2.0, 1.0, 4.0};
        std::shared_ptr<DirichletDistribution> multiplier_prior = std::make_shared<DirichletDistribution>(alphas);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);
        tree.set_population_size_multiplier_prior(multiplier_prior);

        tree.set_root_height(0.1);
        tree.set_mean_population_size(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.estimate_state_frequencies();
        tree.set_freq_1(0.5);

        SampleSummarizer<double> pop_size;
        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;
        SampleSummarizer<double> rate;
        SampleSummarizer<double> f;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 100000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            pop_size.add_sample(tree.get_mean_population_size());
            std::vector<double> multipliers = tree.get_population_sizes_as_multipliers();
            pop_size_root.add_sample(multipliers.at(2));
            pop_size_0.add_sample(multipliers.at(0));
            pop_size_1.add_sample(multipliers.at(1));
            rate.add_sample(tree.get_mutation_rate());
            f.add_sample(tree.get_freq_1());
        }

        REQUIRE(pop_size.mean() == Approx(size_prior->get_mean()).epsilon(0.01));
        REQUIRE(pop_size.variance() == Approx(size_prior->get_variance()).epsilon(0.01));
        REQUIRE(pop_size_root.mean() == Approx(multiplier_prior->get_mean(2) * 3.0).epsilon(0.01));
        REQUIRE(pop_size_root.variance() == Approx(multiplier_prior->get_variance(2) * 9.0).epsilon(0.01));
        REQUIRE(pop_size_0.mean() == Approx(multiplier_prior->get_mean(0) * 3.0).epsilon(0.01));
        REQUIRE(pop_size_0.variance() == Approx(multiplier_prior->get_variance(0) * 9.0).epsilon(0.01));
        REQUIRE(pop_size_1.mean() == Approx(multiplier_prior->get_mean(1) * 3.0).epsilon(0.01));
        REQUIRE(pop_size_1.variance() == Approx(multiplier_prior->get_variance(1) * 9.0).epsilon(0.01));
        REQUIRE(f.mean() == Approx(f_prior->get_mean()).epsilon(0.01));
        REQUIRE(f.variance() == Approx(f_prior->get_variance()).epsilon(0.01));
        REQUIRE(rate.mean() == Approx(rate_prior->get_mean()).epsilon(0.1));
        REQUIRE(rate.variance() == Approx(rate_prior->get_variance()).epsilon(0.1));
    }
}



// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood = -31.77866581319647
TEST_CASE("Testing simple likelihood of RelativeRootPopulationTree", "[RelativeRootPopulationTree]") {

    SECTION("Testing constructor and likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing threaded likelihood of RelativeRootPopulationTree", "[RelativeRootPopulationTree]") {

    SECTION("Testing constructor and threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing over-threaded likelihood of RelativeRootPopulationTree", "[RelativeRootPopulationTree]") {

    SECTION("Testing constructor and over-threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(10);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(10);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood             = -248.93254688526213
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing hemi129.nex RelativeRootPopulationTree likelihood (0.01, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded RelativeRootPopulationTree likelihood (0.01, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7099.716015109998
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.01, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.01, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -6986.120524781545
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.01, 10.0, 1.0, 1.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6986.120524781545));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.01, 10.0, 1.0, 1.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(5);
        REQUIRE(l == Approx(-6986.120524781545));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}



// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -328.39238828878365
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex RelativeRootPopulationTree likelihood (0.0, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded RelativeRootPopulationTree likelihood (0.0, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(7);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(7);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7256.501742344454
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.0, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.0, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.0
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7223.362711937651
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.0, 10.0, 1.0, 1.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7223.362711937651));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.0, 10.0, 1.0, 1.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7223.362711937651));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -227.41048391087554
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex RelativeRootPopulationTree likelihood (0.2, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded RelativeRootPopulationTree likelihood (0.2, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7304.180743441677
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.2, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.2, 10.0, 1.0, 1.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.2
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// Log likelihood            = -7405.145951634711
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.2, 10.0, 1.0, 1.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7405.145951634711));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.2, 10.0, 1.0, 1.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(6);
        REQUIRE(l == Approx(-7405.145951634711));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}




// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -327.7437811413033
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -6472.856486972301
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(5);
        REQUIRE(l == Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(3);
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0 / 19.0
// Log likelihood            = -6494.774924871097
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0, 10.0/19.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6494.774924871097));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0, 10.0/19.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-6494.774924871097));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
  
  
// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -265.0023534261969
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(1);
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -10163.468886613919
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(5);
        REQUIRE(! std::isnan(l));
        REQUIRE(l != Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -10999.288193543642
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0/19.0, 10.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-10999.288193543642));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.03, 10.0, 10.0/19.0, 10.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-10999.288193543642));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError);
    }
}



// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -224.40177558289847
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing hemi129.nex RelativeRootPopulationTree likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-224.40177558289847));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing hemi129.nex threaded RelativeRootPopulationTree likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(100000);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-224.40177558289847));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -8158.88094671241
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-8158.88094671241));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(1000000);
        REQUIRE(l == Approx(-8158.88094671241));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// aflp_25.nex 
// dominant
// height = 0.03
// coalescent_rate = 111.1, 111.1, 111.1
// u = 10.0 / 19.0
// v = 10.0
// Log likelihood            = -8034.250341980543
// Log likelihood correction = -3317.567573476714
TEST_CASE("Testing aflp_25.nex RelativeRootPopulationTree likelihood (0.03, 111.1, 10.0/19.0, 10.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-8034.250341980543));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing aflp_25.nex threaded RelativeRootPopulationTree likelihood (0.03, 111.1, 10.0/19.0, 10.0, dominant)", "[RelativeRootPopulationTree]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-8034.250341980543));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

TEST_CASE("Testing simple prior of RelativeRootPopulationTree", "[RelativeRootPopulationTree]") {

    SECTION("Testing constructor and prior calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);

        tree.set_root_height(0.1);
        tree.set_all_population_sizes(2.0/100.0);
        tree.set_freq_1(0.95);

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(100.0));
        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));
        tree.set_relative_root_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.1));
        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(2.0, 1.0));

        // height       -5.3948298140119091
        // root size    0.22402344985899036
        // leaf sizes   2 * -155.90663080917298
        // f1           0.64185388617239469 
        // total        -472.47286835535851
        
        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-311.58923816848699));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == Approx(0.64185388617239469));

        double d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-316.34221409632647));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-316.34221409632647));


        tree.store_prior_density();
        tree.constrain_state_frequencies();

        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-311.58923816848699));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == 0.0);

        d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-316.98406798249886));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-316.98406798249886));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-316.34221409632647));


        tree.store_prior_density();
        tree.constrain_population_sizes();

        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_population_sizes() == Approx(-155.90663080917298));
        REQUIRE(tree.compute_log_prior_density_of_state_frequencies() == 0.0);

        d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-161.30146062318488));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-161.30146062318488));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-316.98406798249886));

        tree.restore_prior_density();
        REQUIRE(tree.get_log_prior_density_value() == Approx(-316.98406798249886));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-316.98406798249886));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

TEST_CASE("Testing hemi129.nex state manipulation for RelativeRootPopulationTree", "[RelativeRootPopulationTree]") {

    SECTION("Testing state manipulation") {
        std::string nex_path = "data/hemi129.nex";
        RelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(10.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_relative_root_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.2));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(1.5, 3.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_height(0.01);
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1(0.5);
        REQUIRE(tree.is_dirty());

        REQUIRE(tree.get_number_of_likelihood_calculations() == 0);

        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-894.67638803083958));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 1);

        tree.store_state();
        tree.set_root_height(0.0);
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-894.57638803083955));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-894.67638803083958));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 2);

        tree.restore_state();
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-894.67638803083958));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-894.67638803083958));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.compute_log_likelihood_and_prior();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-894.67638803083958));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-894.67638803083958));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.store_state();
        tree.set_root_height(0.0);
        tree.compute_log_likelihood_and_prior();
        REQUIRE(tree.get_number_of_likelihood_calculations() == 4);


        tree.store_state();
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.57638803083955));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-894.57638803083955));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 5);


        tree.store_state();
        tree.set_root_height(0.03);
        tree.set_freq_1(0.05);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.05);
        REQUIRE(tree.get_freq_0() == Approx(0.95));
        REQUIRE(tree.get_u() == Approx(10.0));
        REQUIRE(tree.get_v() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-894.74397280499181));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.57638803083955));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 6);


        tree.store_state();
        tree.set_freq_1(0.95);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.95);
        REQUIRE(tree.get_freq_0() == Approx(0.05));
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-899.1606312737415));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-894.74397280499181));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 7);


        tree.store_state();
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(111.1 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-32.510853039559265));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-899.1606312737415));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 8);

        tree.store_state();
        tree.constrain_state_frequencies();
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.72489170735741));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-32.510853039559265));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 9);


        tree.store_state();
        tree.fold_patterns();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.72489170735741));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.72489170735741));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 10);


        tree.store_state();
        tree.constrain_population_sizes();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.35742912931147));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.72489170735741));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 11);

        tree.store_state();
        tree.estimate_mutation_rate();
        REQUIRE(! tree.mutation_rate_is_fixed());
        tree.set_mutation_rate(1.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());
        tree.set_mutation_rate_prior(std::make_shared<GammaDistribution>(10.0, 0.1));
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.35742912931147));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 12);

        tree.store_state();
        tree.set_mutation_rate(2.0);
        tree.set_root_height(0.1);
        tree.set_all_population_sizes(2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 2.0);
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(20.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-206.13340567945255));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        tree.restore_state();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.13340567945249));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 9fb5b3b7817a0bd4a21e3f90132f132cca72ce4e)
// SNAPP v1.3.0 (master 4f3f0f7366798f4fb38b766c15f6426a75ddf71e)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0/19.0
// Log likelihood            = -23.81984255023975
// Log likelihood correction = -6.87935580446044
//
// With constant sites inclucded and m_bUseNonPolymorphic = true
// Log likelihood            = -55.01646493341547
// Log likelihood correction = -6.87935580446044
TEST_CASE("Testing affect of constant sites on likelihood of RelativeRootPopulationTree", "[RelativeRootPopulationTree]") {

    SECTION("Testing haploid-standard-full-constant.nex") {
        std::string nex_path = "data/haploid-standard-full-constant.nex";
        RelativeRootPopulationTree t_included(
                nex_path, // path
                ' ',      // pop name delimiter
                true,     // pop name is prefix
                false,    // genotypes are diploid
                false,    // markers are dominant
                false,    // constant sites removed
                true);    // validate
        std::string nex_path2 = "data/haploid-standard-full-constant-removed.nex";
        RelativeRootPopulationTree t(
                nex_path2, // path
                ' ',       // pop name delimiter
                true,      // pop name is prefix
                false,     // genotypes are diploid
                false,     // markers are dominant
                true,      // constant sites removed
                true);     // validate

        t.set_root_height(0.03);
        t.set_freq_1(0.05);
        t.set_all_population_sizes(2.0/(10.0 * 2 * t.get_ploidy()));

        t_included.set_root_height(0.03);
        t_included.set_freq_1(0.05);
        t_included.set_all_population_sizes(2.0/(10.0 * 2 * t_included.get_ploidy()));

        double l = t.compute_log_likelihood();
        REQUIRE(l == Approx(-23.81984255023975));
        REQUIRE(t.get_likelihood_correction() == Approx(-6.87935580446044));

        double l_included = t_included.compute_log_likelihood();

        REQUIRE(l_included == Approx(-55.01646493341547));

        REQUIRE(t.get_likelihood_correction() == t_included.get_likelihood_correction());
    }
}

TEST_CASE("Testing affect of constant sites on threaded likelihood of RelativeRootPopulationTree", "[RelativeRootPopulationTree]") {

    SECTION("Testing haploid-standard-full-constant.nex") {
        std::string nex_path = "data/haploid-standard-full-constant.nex";
        RelativeRootPopulationTree t_included(
                nex_path, // path
                ' ',      // pop name delimiter
                true,     // pop name is prefix
                false,    // genotypes are diploid
                false,    // markers are dominant
                false,    // constant sites removed
                true);    // validate
        std::string nex_path2 = "data/haploid-standard-full-constant-removed.nex";
        RelativeRootPopulationTree t(
                nex_path2, // path
                ' ',       // pop name delimiter
                true,      // pop name is prefix
                false,     // genotypes are diploid
                false,     // markers are dominant
                true,      // constant sites removed
                true);     // validate

        t.set_root_height(0.03);
        t.set_freq_1(0.05);
        t.set_all_population_sizes(2.0/(10.0 * 2 * t.get_ploidy()));

        t_included.set_root_height(0.03);
        t_included.set_freq_1(0.05);
        t_included.set_all_population_sizes(2.0/(10.0 * 2 * t_included.get_ploidy()));

        double l = t.compute_log_likelihood(4);
        REQUIRE(l == Approx(-23.81984255023975));
        REQUIRE(t.get_likelihood_correction() == Approx(-6.87935580446044));

        double l_included = t_included.compute_log_likelihood(2);

        REQUIRE(l_included == Approx(-55.01646493341547));

        REQUIRE(t.get_likelihood_correction() == t_included.get_likelihood_correction());
    }
}


// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.01
// coalescent_rate = 10.0, 10.0, 10.0
// u = 1.0
// v = 1.0
// log likelihood = -31.77866581319647
TEST_CASE("Testing simple likelihood of ComparisonRelativeRootPopulationTree", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing constructor and likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_root().get_label() == "root-pop1");
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("Testing simple threaded likelihood of ComparisonRelativeRootPopulationTree", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing constructor and threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        REQUIRE(tree.get_root().get_label() == "root-pop1");
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex state manipulation", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing state manipulation") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(10.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_relative_root_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.2));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_freq_1_prior(std::make_shared<BetaDistribution>(1.5, 3.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_height(0.01);
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_child_population_size(0, 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_child_population_size(1, 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_population_size(2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));

        tree.set_freq_1(0.5);
        REQUIRE(tree.is_dirty());

        REQUIRE(tree.get_number_of_likelihood_calculations() == 0);

        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 1);

        tree.store_state();
        tree.set_root_height(0.0);
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 2);

        tree.restore_state();
        // ComparisonRelativeRootPopulationTree does not store/restore heights
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.compute_log_likelihood_and_prior();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.store_state();
        tree.set_root_height(0.0);
        tree.compute_log_likelihood_and_prior();
        REQUIRE(tree.get_number_of_likelihood_calculations() == 4);


        tree.store_state();
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-328.39238828878365));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 5);


        tree.store_state();
        tree.set_root_height(0.03);
        tree.set_freq_1(0.05);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_freq_1() == 0.05);
        REQUIRE(tree.get_freq_0() == 0.95);
        REQUIRE(tree.get_u() == Approx(10.0));
        REQUIRE(tree.get_v() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-896.74655789798578));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.87897312383359));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 6);


        tree.store_state();
        tree.set_freq_1(0.95);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_freq_1() == 0.95);
        REQUIRE(tree.get_freq_0() == Approx(0.05));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-327.7437811413033));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-901.16321636673547));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-896.74655789798578));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 7);


        tree.store_state();
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == Approx(10.0));
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(111.1 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(111.1 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(111.1 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-265.0023534261969));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-34.513438132553311));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-901.16321636673547));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 8);

        tree.store_state();
        tree.constrain_state_frequencies();
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-224.40177558289847));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-897.02747680035145));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-34.513438132553311));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 9);


        tree.store_state();
        tree.fold_patterns();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_root_population_size() == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(0) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_child_population_size(1) == 2.0/(10.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-897.02747680035145));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-897.02747680035145));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 10);


        tree.store_state();
        tree.constrain_population_sizes();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.66001422230551));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-897.02747680035145));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 11);

        tree.store_state();
        tree.estimate_mutation_rate();
        REQUIRE(! tree.mutation_rate_is_fixed());
        tree.set_mutation_rate(1.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());
        tree.set_mutation_rate_prior(std::make_shared<GammaDistribution>(10.0, 0.1));
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.66001422230551));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 12);

        tree.store_state();
        tree.set_mutation_rate(2.0);
        tree.set_root_height(0.1);
        tree.set_root_population_size(2.0/(20.0 * 2 * tree.get_ploidy()));
        tree.set_child_population_size(0, 2.0/(20.0 * 2 * tree.get_ploidy()));
        tree.set_child_population_size(1, 2.0/(20.0 * 2 * tree.get_ploidy()));
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 2.0);
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(20.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(20.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(20.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-207.43599077244659));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        tree.restore_state();
        // ComparisonRelativeRootPopulationTree does not store/restore node heights
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == Approx(1.0));
        REQUIRE(tree.get_u() == Approx(1.0));
        REQUIRE(tree.get_mutation_rate() == 1.0);
        REQUIRE(tree.get_root_population_size() ==   Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(0) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.get_child_population_size(1) == Approx(2.0/(10.0 * 2 * tree.get_ploidy())));
        REQUIRE(tree.is_dirty());

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(tree.get_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() != Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-447.43599077244653));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.01
// coalescent_rate = 100.0, 200.0, 500.0
// u = 1.0
// v = 1.0
// Log likelihood            = -226.11914854623677
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex likelihood (0.01, 100.0, 200.0, 500.0)", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_child_population_size(0, 2.0/(100.0 * 2 * tree.get_ploidy()));
        tree.set_child_population_size(1, 2.0/(200.0 * 2 * tree.get_ploidy()));
        tree.set_root_population_size(2.0/(500.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-226.11914854623677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}


// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 0.00506843962151613554
// coalescent_rate = 2.0 / 0.00018955324120485613
// u = 1.0
// v = 1.0
// Log likelihood            = -277.06960543551577
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex likelihood (0.00506843962151613554, 2.0 / 0.00018955324120485613)", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.00506843962151613554);
        tree.set_all_population_sizes(0.00018955324120485613 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-277.06960543551577));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 9.08323190033687971e-09
// coalescent_rate = 2.0 / 2.47975039926886321e-08
// u = 1.0
// v = 1.0
// Log likelihood            = -221.69627648370943
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex likelihood (9.08323190033687971e-09, 2.0 / 2.47975039926886321e-08)", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(9.08323190033687971e-09);
        tree.set_all_population_sizes(2.47975039926886321e-08 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-221.69627648370943));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.04921319733994759e-08
// coalescent_rate = 2.0 / 2.75977168733651178e-10
// u = 1.0
// v = 1.0
// Log likelihood            = -324.2737564069293
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex likelihood (1.04921319733994759e-08, 2.0 / 2.75977168733651178e-10)", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.04921319733994759e-08);
        tree.set_all_population_sizes(2.75977168733651178e-10 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-324.2737564069293));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.012386610001351e-08
// coalescent_rate = 2.0 / 5.75048645855884647e-30
// u = 1.0
// v = 1.0
// Log likelihood            = -1364.1427530000253
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex likelihood (1.012386610001351e-08, 2.0 / 5.75048645855884647e-30)", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.012386610001351e-08);
        tree.set_all_population_sizes(5.75048645855884647e-30 / 4.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-1364.1427530000253));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.04856228318474786e-08
// coalescent_rate = 2.0 / 4.43934332792563837e-305
// u = 1.0
// v = 1.0
// Log likelihood            = NaN
// Log likelihood correction = -135.97095011239867
TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex likelihood (1.04856228318474786e-08, 2.0 / 4.43934332792563837e-305)", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.04856228318474786e-08);
        tree.set_all_population_sizes(4.43934332792563837e-305 / 4.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(std::isnan(l));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.012386610001351e-08
// coalescent_rate = 2.0 / 5.46641122085615013e-09,
//                   2.0 / 7.39871781998828579e-08,
//                   2.0 / 2.71077053326002069e-13
// u = 1.0
// v = 1.0
// Log likelihood            = -96.34394008351177
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex weirdness", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.012386610001351e-08);
        tree.set_child_population_size(0, 5.46641122085615013e-09 / 4.0);
        tree.set_child_population_size(1, 7.39871781998828579e-08 / 4.0);
        tree.set_root_population_size(2.71077053326002069e-13 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-96.34394008351177));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

// BEAST v2.4.0 (master 3731dff6884f7dd27b288099027dc1d500a3a9d8)
// SNAPP v1.3.0 (master 24d18026c774b10f2e79de16100d96e1a5df1b96)
// hemi129.nex 
// height = 1.036374107244057e-08
// coalescent_rate = 2.0 / 4.57999694763258361e-09,
//                   2.0 / 6.70991782555376588e-08,
//                   2.0 / 1.33514111020266258e-08
// u = 1.0
// v = 1.0
// Log likelihood            = -44.95791900747736
// log likelihood correction  = -135.97095011239867
TEST_CASE("Testing ComparisonRelativeRootPopulationTree hemi129.nex weirdness 2", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(1.036374107244057e-08);
        tree.set_child_population_size(0, 4.57999694763258361e-09 / 4.0);
        tree.set_child_population_size(1, 6.70991782555376588e-08 / 4.0);
        tree.set_root_population_size(1.33514111020266258e-08 / 4.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-44.95791900747736));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 2 lineages and theta of 1.0",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 2 lineages and theta of 3.7",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 3.7;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 2 lineages and theta of 0.17",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 3 lineages and theta of 1.0",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 3 lineages and theta of 1.47",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 1.47;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 3 lineages and theta of 0.17",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 10 lineages and theta of 1.0",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 10 lineages and theta of 1.47",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 1.47;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree coalesce_in_branch for 10 lineages and theta of 0.17",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = ComparisonRelativeRootPopulationTree::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree scaling of simulate_gene_tree for singleton",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        std::string nex_path = "data/singleton-n2.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, false, false);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.set_child_population_size(0, (Ne * mu));

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        unsigned int n;
        double mean;
        double variance;
        double sum_devs;
        double d;
        double d_n;
        double mn;
        double mx;
        double x;
        std::shared_ptr<GeneTreeSimNode> gtree;

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            x = gtree->get_height();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);

        REQUIRE(gtree->is_root());
        REQUIRE(gtree->get_number_of_children() == 2);
        REQUIRE(gtree->get_leaf_node_count() == 2);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.0005));
        REQUIRE(mn >= 0.0);

        tree.set_mutation_rate(mu);
        tree.set_root_population_size(Ne);
        tree.set_child_population_size(0, Ne);
        tree.set_root_height(0.1 / mu);

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            x = gtree->get_height();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);

        REQUIRE(gtree->is_root());
        REQUIRE(gtree->get_number_of_children() == 2);
        REQUIRE(gtree->get_leaf_node_count() == 2);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree scaling of simulate_gene_tree for pair",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing pair") {
        unsigned int Ne_root = 100000;
        unsigned int Ne_0 = 200000;
        unsigned int Ne_1 = 10000;
        double mu = 1e-8;
        double theta_root = 4 * Ne_root * mu;
        double theta_0 = 4 * Ne_0 * mu;
        double theta_1 = 4 * Ne_1 * mu;
        double time = 2e7;

        double expected_mean_root = theta_root * (1.0 - (1.0 / 2));
        double expected_variance_root = expected_mean_root * expected_mean_root;
        expected_mean_root += (time * mu);

        double expected_mean_0 = theta_0 * (1.0 - (1.0 / 10));
        double expected_mean_1 = theta_1 * (1.0 - (1.0 / 10));

        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);
        tree.estimate_mutation_rate();

        tree.set_child_population_size(0, (Ne_0));
        tree.set_child_population_size(1, (Ne_1));
        tree.set_root_population_size(Ne_root);

        tree.set_root_height(time);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(mu);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        std::shared_ptr<GeneTreeSimNode> gtree;
        SampleSummarizer<double> height_root;
        SampleSummarizer<double> height_0;
        SampleSummarizer<double> height_1;
        SampleSummarizer<double> tree_length;

        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            REQUIRE(gtree->is_root());
            REQUIRE(gtree->get_number_of_children() == 2);
            REQUIRE(gtree->get_leaf_node_count() == 20);
            height_root.add_sample(gtree->get_height());
            double h0 = gtree->get_child(0)->get_height();
            double h1 = gtree->get_child(1)->get_height();
            if (h0 < h1) {
                height_0.add_sample(h1);
                height_1.add_sample(h0);
            }
            else {
                height_0.add_sample(h0);
                height_1.add_sample(h1);
            }
            tree_length.add_sample(gtree->get_clade_length());
        }

        REQUIRE(height_root.mean() == Approx(expected_mean_root).epsilon(0.0005));
        REQUIRE(height_root.variance() == Approx(expected_variance_root).epsilon(0.0005));

        REQUIRE(height_0.mean() == Approx(expected_mean_0).epsilon(0.00003));
        REQUIRE(height_1.mean() == Approx(expected_mean_1).epsilon(0.00003));


        tree.set_mutation_rate(mu * 2.0);

        SampleSummarizer<double> scale_height_root;
        SampleSummarizer<double> scale_height_0;
        SampleSummarizer<double> scale_height_1;
        SampleSummarizer<double> scale_tree_length;

        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            gtree->scale(0.5);
            REQUIRE(gtree->is_root());
            REQUIRE(gtree->get_number_of_children() == 2);
            REQUIRE(gtree->get_leaf_node_count() == 20);
            scale_height_root.add_sample(gtree->get_height());
            double h0 = gtree->get_child(0)->get_height();
            double h1 = gtree->get_child(1)->get_height();
            if (h0 < h1) {
                scale_height_0.add_sample(h1);
                scale_height_1.add_sample(h0);
            }
            else {
                scale_height_0.add_sample(h0);
                scale_height_1.add_sample(h1);
            }
            scale_tree_length.add_sample(gtree->get_clade_length());
        }

        REQUIRE(scale_height_root.mean() == Approx(expected_mean_root).epsilon(0.0005));
        REQUIRE(scale_height_root.variance() == Approx(expected_variance_root).epsilon(0.0005));

        REQUIRE(scale_height_0.mean() == Approx(expected_mean_0).epsilon(0.00003));
        REQUIRE(scale_height_1.mean() == Approx(expected_mean_1).epsilon(0.00003));

        REQUIRE(scale_tree_length.mean() == Approx(tree_length.mean()).epsilon(0.0001));
        REQUIRE(scale_tree_length.variance() == Approx(tree_length.variance()).epsilon(0.0001));
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree dataset simulation", "[ComparisonRelativeRootPopulationTree]") {
    SECTION("Testing simulate_biallelic_data_set for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, false, false, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        BiallelicData data = tree.simulate_biallelic_data_set(rng);

        REQUIRE(data.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data.markers_are_dominant() == false);
        REQUIRE(tree.get_data().markers_are_dominant() == false);

        data.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8).epsilon(0.01));
        REQUIRE(data.get_proportion_1() == Approx(0.625).epsilon(0.01));
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree draw_from_prior for fully fixed", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fully fixed pair") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.001);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_root_population_size() == 0.001);
            REQUIRE(tree.get_child_population_size(0) == 0.001);
            REQUIRE(tree.get_child_population_size(1) == 0.001);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);
        }
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree draw_from_prior for constrained sizes", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for constrained sizes") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        SampleSummarizer<double> pop_size;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() == tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() == tree.get_child_population_size(1));

            pop_size.add_sample(tree.get_root_population_size());
        }

        REQUIRE(pop_size.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree draw_from_prior for unconstrained sizes", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for unconstrained sizes") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);
        std::shared_ptr<GammaDistribution> rel_root_size_prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);
        tree.set_relative_root_population_size_prior(rel_root_size_prior);

        tree.set_root_height(0.1);
        tree.set_all_population_sizes(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();

        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            pop_size_root.add_sample(tree.get_relative_root_population_size());
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
        }

        REQUIRE(pop_size_root.mean() == Approx(rel_root_size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_root.variance() == Approx(rel_root_size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_0.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_0.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_1.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_1.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree draw_from_prior for fully parameterized", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fully parameterized") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(0.5, 0.8);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);
        std::shared_ptr<GammaDistribution> rel_root_size_prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);
        tree.set_relative_root_population_size_prior(rel_root_size_prior);

        tree.set_root_height(0.1);
        tree.set_all_population_sizes(0.001);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.estimate_state_frequencies();
        tree.set_freq_1(0.5);

        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;
        SampleSummarizer<double> rate;
        SampleSummarizer<double> f;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            pop_size_root.add_sample(tree.get_relative_root_population_size());
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
            rate.add_sample(tree.get_mutation_rate());
            f.add_sample(tree.get_freq_1());
        }

        REQUIRE(pop_size_root.mean() == Approx(rel_root_size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_root.variance() == Approx(rel_root_size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_0.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_0.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_1.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_1.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(f.mean() == Approx(f_prior->get_mean()).epsilon(0.01));
        REQUIRE(f.variance() == Approx(f_prior->get_variance()).epsilon(0.01));
        REQUIRE(rate.mean() == Approx(rate_prior->get_mean()).epsilon(0.1));
        REQUIRE(rate.variance() == Approx(rate_prior->get_variance()).epsilon(0.1));
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree draw_from_prior for fixed relative root size",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fixed relative root size") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);
        std::shared_ptr<GammaDistribution> rel_root_size_prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);
        tree.set_relative_root_population_size_prior(rel_root_size_prior);

        tree.set_root_height(0.1);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();
        tree.set_all_population_sizes(0.001);
        tree.set_root_population_size(0.002);
        tree.fix_relative_root_population_size();

        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            REQUIRE(tree.get_relative_root_population_size() == Approx(2.0));
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
        }

        REQUIRE(pop_size_0.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_0.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_1.mean() == Approx(size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_1.variance() == Approx(size_prior->get_variance()).epsilon(0.1));
    }
}

TEST_CASE("Testing ComparisonRelativeRootPopulationTree draw_from_prior for fixed leaf sizes",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fixed leaf sizes") {
        std::string nex_path = "data/hemi129-5-5.nex";
        ComparisonRelativeRootPopulationTree tree(nex_path, ' ', true, true, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);
        std::shared_ptr<GammaDistribution> rel_root_size_prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        tree.set_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);
        tree.set_relative_root_population_size_prior(rel_root_size_prior);

        tree.set_root_height(0.1);
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(0.8);
        tree.fix_mutation_rate();
        tree.constrain_state_frequencies();
        tree.set_all_population_sizes(0.001);
        tree.fix_population_sizes();

        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_height() == 0.1);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) == tree.get_child_population_size(1));

            pop_size_root.add_sample(tree.get_relative_root_population_size());
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
        }

        REQUIRE(pop_size_root.mean() == Approx(rel_root_size_prior->get_mean()).epsilon(0.1));
        REQUIRE(pop_size_root.variance() == Approx(rel_root_size_prior->get_variance()).epsilon(0.1));
        REQUIRE(pop_size_0.mean() == Approx(0.001));
        REQUIRE(pop_size_0.variance() == Approx(0.0));
        REQUIRE(pop_size_1.mean() == Approx(0.001));
        REQUIRE(pop_size_1.variance() == Approx(0.0));
    }
}

TEST_CASE("Testing config ComparisonRelativeRootPopulationTree draw_from_prior for fully fixed", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fully fixed pair") {

        std::stringstream os;
        os << "path: data/diploid-standard-data-ntax5-nchar5.nex\n";
        os << "genotypes_are_diploid: true\n";
        os << "markers_are_dominant: false\n";
        os << "population_name_delimiter: \" \"\n";
        os << "population_name_is_prefix: true\n";
        os << "constant_sites_removed: true\n";
        os << "equal_population_sizes: false\n";
        os << "parameters:\n";
        os << "    mutation_rate:\n";
        os << "        value: 0.8\n";
        os << "        estimate: false\n";
        os << "    population_size:\n";
        os << "        value: 0.001\n";
        os << "        estimate: false\n";
        os << "    root_relative_population_size:\n";
        os << "        value: 1.0\n";
        os << "        estimate: false\n";
        os << "    freq_1:\n";
        os << "        value: 0.5\n";
        os << "        estimate: false\n";

        YAML::Node n;
        n = YAML::Load(os);
        RelativeRootComparisonSettings settings = RelativeRootComparisonSettings(
                n,
                "dummy-config-path.yml");

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        ComparisonRelativeRootPopulationTree tree = ComparisonRelativeRootPopulationTree(settings, rng);

        for (unsigned int i = 0; i < 10; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_population_size() == 0.001);
            REQUIRE(tree.get_child_population_size(0) == 0.001);
            REQUIRE(tree.get_child_population_size(1) == 0.001);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);
        }
    }
}

TEST_CASE("Testing config ComparisonRelativeRootPopulationTree draw_from_prior for fully fixed and constrained", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fully fixed pair") {

        std::stringstream os;
        os << "path: data/diploid-standard-data-ntax5-nchar5.nex\n";
        os << "genotypes_are_diploid: true\n";
        os << "markers_are_dominant: false\n";
        os << "population_name_delimiter: \" \"\n";
        os << "population_name_is_prefix: true\n";
        os << "constant_sites_removed: true\n";
        os << "equal_population_sizes: true\n";
        os << "parameters:\n";
        os << "    mutation_rate:\n";
        os << "        value: 0.8\n";
        os << "        estimate: false\n";
        os << "    population_size:\n";
        os << "        value: 0.001\n";
        os << "        estimate: false\n";
        os << "    root_relative_population_size:\n";
        os << "        value: 1.0\n";
        os << "        estimate: false\n";
        os << "    freq_1:\n";
        os << "        value: 0.5\n";
        os << "        estimate: false\n";

        YAML::Node n;
        n = YAML::Load(os);
        RelativeRootComparisonSettings settings = RelativeRootComparisonSettings(
                n,
                "dummy-config-path.yml");

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        ComparisonRelativeRootPopulationTree tree = ComparisonRelativeRootPopulationTree(settings, rng);

        for (unsigned int i = 0; i < 10; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_population_size() == 0.001);
            REQUIRE(tree.get_child_population_size(0) == 0.001);
            REQUIRE(tree.get_child_population_size(1) == 0.001);
            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);
        }
    }
}

TEST_CASE("Testing config ComparisonRelativeRootPopulationTree draw_from_prior for constrained sizes", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for constrained sizes") {

        double size_shape = 2.0;
        double size_scale = 1.2;

        std::stringstream os;
        os << "path: data/diploid-standard-data-ntax5-nchar5.nex\n";
        os << "genotypes_are_diploid: true\n";
        os << "markers_are_dominant: false\n";
        os << "population_name_delimiter: \" \"\n";
        os << "population_name_is_prefix: true\n";
        os << "constant_sites_removed: true\n";
        os << "equal_population_sizes: true\n";
        os << "parameters:\n";
        os << "    mutation_rate:\n";
        os << "        value: 0.8\n";
        os << "        estimate: false\n";
        os << "    population_size:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << size_shape << "\n";
        os << "                scale: " << size_scale << "\n";
        os << "    root_relative_population_size:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: 20.0\n";
        os << "                scale: 0.05\n";
        os << "    freq_1:\n";
        os << "        value: 0.5\n";
        os << "        estimate: false\n";

        YAML::Node n;
        n = YAML::Load(os);
        RelativeRootComparisonSettings settings = RelativeRootComparisonSettings(
                n,
                "dummy-config-path.yml");

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        ComparisonRelativeRootPopulationTree tree = ComparisonRelativeRootPopulationTree(settings, rng);

        SampleSummarizer<double> pop_size;

        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_relative_root_population_size() == Approx(1.0));
            REQUIRE(tree.get_root_population_size() == tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() == tree.get_child_population_size(1));

            pop_size.add_sample(tree.get_root_population_size());
        }

        REQUIRE(pop_size.mean() == Approx(size_shape * size_scale).epsilon(0.1));
        REQUIRE(pop_size.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.1));
    }
}

TEST_CASE("Testing config ComparisonRelativeRootPopulationTree draw_from_prior for unconstrained sizes", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for unconstrained sizes") {
        double size_shape = 2.0;
        double size_scale = 1.2;
        double rel_size_shape = 10.0;
        double rel_size_scale = 0.1;

        std::stringstream os;
        os << "path: data/diploid-standard-data-ntax5-nchar5.nex\n";
        os << "genotypes_are_diploid: true\n";
        os << "markers_are_dominant: false\n";
        os << "population_name_delimiter: \" \"\n";
        os << "population_name_is_prefix: true\n";
        os << "constant_sites_removed: true\n";
        os << "equal_population_sizes: false\n";
        os << "parameters:\n";
        os << "    mutation_rate:\n";
        os << "        value: 0.8\n";
        os << "        estimate: false\n";
        os << "    population_size:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << size_shape << "\n";
        os << "                scale: " << size_scale << "\n";
        os << "    root_relative_population_size:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << rel_size_shape << "\n";
        os << "                scale: " << rel_size_scale << "\n";
        os << "    freq_1:\n";
        os << "        value: 0.5\n";
        os << "        estimate: false\n";

        YAML::Node n;
        n = YAML::Load(os);
        RelativeRootComparisonSettings settings = RelativeRootComparisonSettings(
                n,
                "dummy-config-path.yml");

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        ComparisonRelativeRootPopulationTree tree = ComparisonRelativeRootPopulationTree(settings, rng);

        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;

        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            pop_size_root.add_sample(tree.get_relative_root_population_size());
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
        }

        REQUIRE(pop_size_root.mean() == Approx(rel_size_shape * rel_size_scale).epsilon(0.1));
        REQUIRE(pop_size_root.variance() == Approx(rel_size_shape * rel_size_scale * rel_size_scale).epsilon(0.1));
        REQUIRE(pop_size_0.mean() == Approx(size_shape * size_scale).epsilon(0.1));
        REQUIRE(pop_size_0.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.1));
        REQUIRE(pop_size_1.mean() == Approx(size_shape * size_scale).epsilon(0.1));
        REQUIRE(pop_size_1.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.1));
    }
}

TEST_CASE("Testing config ComparisonRelativeRootPopulationTree draw_from_prior for fully parameterized", "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fully parameterized") {
        double size_shape = 2.0;
        double size_scale = 1.2;
        double rel_size_shape = 10.0;
        double rel_size_scale = 0.1;
        double mu_shape = 3.0;
        double mu_scale = 1.1;
        double f_a = 0.5;
        double f_b = 0.8;

        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(f_a, f_b);

        std::stringstream os;
        os << "path: data/diploid-standard-data-ntax5-nchar5.nex\n";
        os << "genotypes_are_diploid: true\n";
        os << "markers_are_dominant: false\n";
        os << "population_name_delimiter: \" \"\n";
        os << "population_name_is_prefix: true\n";
        os << "constant_sites_removed: true\n";
        os << "equal_population_sizes: false\n";
        os << "parameters:\n";
        os << "    mutation_rate:\n";
        os << "        value: 0.8\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << mu_shape << "\n";
        os << "                scale: " << mu_scale << "\n";
        os << "    population_size:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << size_shape << "\n";
        os << "                scale: " << size_scale << "\n";
        os << "    root_relative_population_size:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << rel_size_shape << "\n";
        os << "                scale: " << rel_size_scale << "\n";
        os << "    freq_1:\n";
        os << "        value: 0.5\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            beta_distribution:\n";
        os << "                alpha: " << f_a << "\n";
        os << "                beta: " << f_b << "\n";

        YAML::Node n;
        n = YAML::Load(os);
        RelativeRootComparisonSettings settings = RelativeRootComparisonSettings(
                n,
                "dummy-config-path.yml");

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        ComparisonRelativeRootPopulationTree tree = ComparisonRelativeRootPopulationTree(settings, rng);

        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;
        SampleSummarizer<double> rate;
        SampleSummarizer<double> f;

        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            pop_size_root.add_sample(tree.get_relative_root_population_size());
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
            rate.add_sample(tree.get_mutation_rate());
            f.add_sample(tree.get_freq_1());
        }

        REQUIRE(pop_size_root.mean() == Approx(rel_size_shape * rel_size_scale).epsilon(0.1));
        REQUIRE(pop_size_root.variance() == Approx(rel_size_shape * rel_size_scale * rel_size_scale).epsilon(0.1));
        REQUIRE(pop_size_0.mean() == Approx(size_shape * size_scale).epsilon(0.1));
        REQUIRE(pop_size_0.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.1));
        REQUIRE(pop_size_1.mean() == Approx(size_shape * size_scale).epsilon(0.1));
        REQUIRE(pop_size_1.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.1));
        REQUIRE(f.mean() == Approx(f_prior->get_mean()).epsilon(0.01));
        REQUIRE(f.variance() == Approx(f_prior->get_variance()).epsilon(0.01));
        REQUIRE(rate.mean() == Approx(mu_shape * mu_scale).epsilon(0.1));
        REQUIRE(rate.variance() == Approx(mu_shape * mu_scale * mu_scale).epsilon(0.1));
    }
}

TEST_CASE("Testing config ComparisonRelativeRootPopulationTree draw_from_prior for fixed relative root size",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fixed relative root size") {
        double size_shape = 2.0;
        double size_scale = 1.2;

        std::stringstream os;
        os << "path: data/diploid-standard-data-ntax5-nchar5.nex\n";
        os << "genotypes_are_diploid: true\n";
        os << "markers_are_dominant: false\n";
        os << "population_name_delimiter: \" \"\n";
        os << "population_name_is_prefix: true\n";
        os << "constant_sites_removed: true\n";
        os << "equal_population_sizes: false\n";
        os << "parameters:\n";
        os << "    mutation_rate:\n";
        os << "        value: 0.8\n";
        os << "        estimate: false\n";
        os << "    population_size:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << size_shape << "\n";
        os << "                scale: " << size_scale << "\n";
        os << "    root_relative_population_size:\n";
        os << "        value: 2.0\n";
        os << "        estimate: false\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: 10.0\n";
        os << "                scale: 0.1\n";
        os << "    freq_1:\n";
        os << "        value: 0.5\n";
        os << "        estimate: false\n";

        YAML::Node n;
        n = YAML::Load(os);
        RelativeRootComparisonSettings settings = RelativeRootComparisonSettings(
                n,
                "dummy-config-path.yml");

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        ComparisonRelativeRootPopulationTree tree = ComparisonRelativeRootPopulationTree(settings, rng);

        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;

        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) != tree.get_child_population_size(1));

            REQUIRE(tree.get_relative_root_population_size() == Approx(2.0));
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
        }

        REQUIRE(pop_size_0.mean() == Approx(size_shape * size_scale).epsilon(0.1));
        REQUIRE(pop_size_0.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.1));
        REQUIRE(pop_size_1.mean() == Approx(size_shape * size_scale).epsilon(0.1));
        REQUIRE(pop_size_1.variance() == Approx(size_shape * size_scale * size_scale).epsilon(0.1));
    }
}

TEST_CASE("Testing config ComparisonRelativeRootPopulationTree draw_from_prior for fixed leaf sizes",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing draw_from_prior for fixed leaf sizes") {
        double rel_size_shape = 10.0;
        double rel_size_scale = 0.1;

        std::stringstream os;
        os << "path: data/diploid-standard-data-ntax5-nchar5.nex\n";
        os << "genotypes_are_diploid: true\n";
        os << "markers_are_dominant: false\n";
        os << "population_name_delimiter: \" \"\n";
        os << "population_name_is_prefix: true\n";
        os << "constant_sites_removed: true\n";
        os << "equal_population_sizes: false\n";
        os << "parameters:\n";
        os << "    mutation_rate:\n";
        os << "        value: 0.8\n";
        os << "        estimate: false\n";
        os << "    population_size:\n";
        os << "        value: 0.001\n";
        os << "        estimate: false\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: 2.0\n";
        os << "                scale: 1.2\n";
        os << "    root_relative_population_size:\n";
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            gamma_distribution:\n";
        os << "                shape: " << rel_size_shape << "\n";
        os << "                scale: " << rel_size_scale << "\n";
        os << "    freq_1:\n";
        os << "        value: 0.5\n";
        os << "        estimate: false\n";

        YAML::Node n;
        n = YAML::Load(os);
        RelativeRootComparisonSettings settings = RelativeRootComparisonSettings(
                n,
                "dummy-config-path.yml");

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        ComparisonRelativeRootPopulationTree tree = ComparisonRelativeRootPopulationTree(settings, rng);


        SampleSummarizer<double> pop_size_root;
        SampleSummarizer<double> pop_size_0;
        SampleSummarizer<double> pop_size_1;

        for (unsigned int i = 0; i < 10000; ++i) {
            tree.draw_from_prior(rng);

            REQUIRE(tree.get_u() == Approx(1.0));
            REQUIRE(tree.get_v() == Approx(1.0));
            REQUIRE(tree.get_freq_1() == 0.5);
            REQUIRE(tree.get_mutation_rate() == 0.8);

            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(0));
            REQUIRE(tree.get_root_population_size() != tree.get_child_population_size(1));
            REQUIRE(tree.get_child_population_size(0) == tree.get_child_population_size(1));

            pop_size_root.add_sample(tree.get_relative_root_population_size());
            pop_size_0.add_sample(tree.get_child_population_size(0));
            pop_size_1.add_sample(tree.get_child_population_size(1));
        }

        REQUIRE(pop_size_root.mean() == Approx(rel_size_shape * rel_size_scale).epsilon(0.1));
        REQUIRE(pop_size_root.variance() == Approx(rel_size_shape * rel_size_scale * rel_size_scale).epsilon(0.1));
        REQUIRE(pop_size_0.mean() == Approx(0.001));
        REQUIRE(pop_size_0.variance() == Approx(0.0));
        REQUIRE(pop_size_1.mean() == Approx(0.001));
        REQUIRE(pop_size_1.variance() == Approx(0.0));
    }
}



TEST_CASE("Testing missing gene copies",
        "[ComparisonPopulationTree]") {

    SECTION("Testing missing gene copies for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        std::string nex_path2 = "data/singleton-n2.nex";
        std::string nex_path3 = "data/singleton-n3.nex";
        // Need to keep constant characters
        ComparisonPopulationTree tree2(nex_path2, ' ', true, false, false, false);
        ComparisonPopulationTree tree3(nex_path3, ' ', true, false, false, false);
        tree2.estimate_mutation_rate();
        tree3.estimate_mutation_rate();

        tree2.set_root_population_size(Ne * mu);
        tree3.set_root_population_size(Ne * mu);

        tree2.set_child_population_size(0, (Ne * mu));
        tree3.set_child_population_size(0, (Ne * mu));

        tree2.set_root_height(0.1);
        tree3.set_root_height(0.1);

        tree2.set_freq_1(0.5);
        tree3.set_freq_1(0.5);

        tree2.set_mutation_rate(1.0);
        tree3.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        unsigned int nsamples = 1000000;

        std::shared_ptr<GeneTreeSimNode> gtree2;
        std::shared_ptr<GeneTreeSimNode> gtree3;
        std::vector<unsigned int> down_sampled_counts = {2};

        std::vector<unsigned int> red_allele_counts2 = {0, 0, 0};
        std::vector<unsigned int> red_allele_counts3 = {0, 0, 0};

        for (unsigned int i = 0; i < nsamples; ++i) {
            gtree2 = tree2.simulate_gene_tree(0, rng, true);
            gtree3 = tree3.simulate_gene_tree(0, rng, true);

            auto pattern2 = tree2.simulate_biallelic_site(gtree2, rng);
            auto pattern3 = tree3.simulate_biallelic_site_sans_missing(gtree3, down_sampled_counts, rng);

            ++red_allele_counts2.at(pattern2.first.at(0));
            ++red_allele_counts3.at(pattern3.first.at(0));
        }

        unsigned int total2 = 0;
        unsigned int total3 = 0;
        for (unsigned int i = 0; i < red_allele_counts2.size(); ++i) {
            total2 += red_allele_counts2.at(i);
            total3 += red_allele_counts3.at(i);
        }

        REQUIRE(total2 == nsamples);
        REQUIRE(total3 == nsamples);

        std::vector<double> red_allele_freqs2 = {0.0, 0.0, 0.0};
        std::vector<double> red_allele_freqs3 = {0.0, 0.0, 0.0};
        for (unsigned int i = 0; i < red_allele_counts2.size(); ++i) {
            red_allele_freqs2.at(i) = red_allele_counts2.at(i) / (double)nsamples;
            red_allele_freqs3.at(i) = red_allele_counts3.at(i) / (double)nsamples;
        }

        std::cout << "Freqs of red allele counts for 2 gene copies:\n";
        for (unsigned int i = 0; i < red_allele_freqs2.size(); ++i) {
            std::cout << i << ": " << red_allele_freqs2.at(i) << "\n";
        }
        std::cout << "Freqs of red allele counts for 3 gene copies with 1 dropped:\n";
        for (unsigned int i = 0; i < red_allele_freqs3.size(); ++i) {
            std::cout << i << ": " << red_allele_freqs3.at(i) << "\n";
        }

        for (unsigned int i = 0; i < red_allele_freqs2.size(); ++i) {
            REQUIRE(red_allele_freqs2.at(i) == Approx(red_allele_freqs3.at(i)).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing missing gene copies (1 of 3) with relative root",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing missing gene copies (1 of 3) for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        std::string nex_path2 = "data/singleton-n2.nex";
        std::string nex_path3 = "data/singleton-n3.nex";
        // Need to keep constant characters
        ComparisonRelativeRootPopulationTree tree2(nex_path2, ' ', true, false, false, false);
        ComparisonRelativeRootPopulationTree tree3(nex_path3, ' ', true, false, false, false);
        tree2.estimate_mutation_rate();
        tree3.estimate_mutation_rate();

        tree2.set_root_population_size(Ne * mu);
        tree3.set_root_population_size(Ne * mu);

        tree2.set_child_population_size(0, (Ne * mu));
        tree3.set_child_population_size(0, (Ne * mu));

        tree2.set_root_height(0.1);
        tree3.set_root_height(0.1);

        tree2.set_freq_1(0.5);
        tree3.set_freq_1(0.5);

        tree2.set_mutation_rate(1.0);
        tree3.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        unsigned int nsamples = 1000000;

        std::shared_ptr<GeneTreeSimNode> gtree2;
        std::shared_ptr<GeneTreeSimNode> gtree3;
        std::vector<unsigned int> down_sampled_counts = {2};

        std::vector<unsigned int> red_allele_counts2 = {0, 0, 0};
        std::vector<unsigned int> red_allele_counts3 = {0, 0, 0};

        for (unsigned int i = 0; i < nsamples; ++i) {
            gtree2 = tree2.simulate_gene_tree(0, rng, true);
            gtree3 = tree3.simulate_gene_tree(0, rng, true);

            auto pattern2 = tree2.simulate_biallelic_site(gtree2, rng);
            auto pattern3 = tree3.simulate_biallelic_site_sans_missing(gtree3, down_sampled_counts, rng);

            ++red_allele_counts2.at(pattern2.first.at(0));
            ++red_allele_counts3.at(pattern3.first.at(0));
        }

        unsigned int total2 = 0;
        unsigned int total3 = 0;
        for (unsigned int i = 0; i < red_allele_counts2.size(); ++i) {
            total2 += red_allele_counts2.at(i);
            total3 += red_allele_counts3.at(i);
        }

        REQUIRE(total2 == nsamples);
        REQUIRE(total3 == nsamples);

        std::vector<double> red_allele_freqs2 = {0.0, 0.0, 0.0};
        std::vector<double> red_allele_freqs3 = {0.0, 0.0, 0.0};
        for (unsigned int i = 0; i < red_allele_counts2.size(); ++i) {
            red_allele_freqs2.at(i) = red_allele_counts2.at(i) / (double)nsamples;
            red_allele_freqs3.at(i) = red_allele_counts3.at(i) / (double)nsamples;
        }

        std::cout << "Freqs of red allele counts for 2 gene copies:\n";
        for (unsigned int i = 0; i < red_allele_freqs2.size(); ++i) {
            std::cout << i << ": " << red_allele_freqs2.at(i) << "\n";
        }
        std::cout << "Freqs of red allele counts for 3 gene copies with 1 dropped:\n";
        for (unsigned int i = 0; i < red_allele_freqs3.size(); ++i) {
            std::cout << i << ": " << red_allele_freqs3.at(i) << "\n";
        }

        for (unsigned int i = 0; i < red_allele_freqs2.size(); ++i) {
            REQUIRE(red_allele_freqs2.at(i) == Approx(red_allele_freqs3.at(i)).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing missing gene copies (2 of 4) with relative root",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing missing gene copies (2 of 4) for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        std::string nex_path2 = "data/singleton-n2.nex";
        std::string nex_path3 = "data/singleton-n4.nex";
        // Need to keep constant characters
        ComparisonRelativeRootPopulationTree tree2(nex_path2, ' ', true, false, false, false);
        ComparisonRelativeRootPopulationTree tree3(nex_path3, ' ', true, false, false, false);
        tree2.estimate_mutation_rate();
        tree3.estimate_mutation_rate();

        tree2.set_root_population_size(Ne * mu);
        tree3.set_root_population_size(Ne * mu);

        tree2.set_child_population_size(0, (Ne * mu));
        tree3.set_child_population_size(0, (Ne * mu));

        tree2.set_root_height(0.1);
        tree3.set_root_height(0.1);

        tree2.set_freq_1(0.5);
        tree3.set_freq_1(0.5);

        tree2.set_mutation_rate(1.0);
        tree3.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        unsigned int nsamples = 1000000;

        std::shared_ptr<GeneTreeSimNode> gtree2;
        std::shared_ptr<GeneTreeSimNode> gtree3;
        std::vector<unsigned int> down_sampled_counts = {2};

        std::vector<unsigned int> red_allele_counts2 = {0, 0, 0};
        std::vector<unsigned int> red_allele_counts3 = {0, 0, 0};

        for (unsigned int i = 0; i < nsamples; ++i) {
            gtree2 = tree2.simulate_gene_tree(0, rng, true);
            gtree3 = tree3.simulate_gene_tree(0, rng, true);

            auto pattern2 = tree2.simulate_biallelic_site(gtree2, rng);
            auto pattern3 = tree3.simulate_biallelic_site_sans_missing(gtree3, down_sampled_counts, rng);

            ++red_allele_counts2.at(pattern2.first.at(0));
            ++red_allele_counts3.at(pattern3.first.at(0));
        }

        unsigned int total2 = 0;
        unsigned int total3 = 0;
        for (unsigned int i = 0; i < red_allele_counts2.size(); ++i) {
            total2 += red_allele_counts2.at(i);
            total3 += red_allele_counts3.at(i);
        }

        REQUIRE(total2 == nsamples);
        REQUIRE(total3 == nsamples);

        std::vector<double> red_allele_freqs2 = {0.0, 0.0, 0.0};
        std::vector<double> red_allele_freqs3 = {0.0, 0.0, 0.0};
        for (unsigned int i = 0; i < red_allele_counts2.size(); ++i) {
            red_allele_freqs2.at(i) = red_allele_counts2.at(i) / (double)nsamples;
            red_allele_freqs3.at(i) = red_allele_counts3.at(i) / (double)nsamples;
        }

        std::cout << "Freqs of red allele counts for 2 gene copies:\n";
        for (unsigned int i = 0; i < red_allele_freqs2.size(); ++i) {
            std::cout << i << ": " << red_allele_freqs2.at(i) << "\n";
        }
        std::cout << "Freqs of red allele counts for 3 gene copies with 1 dropped:\n";
        for (unsigned int i = 0; i < red_allele_freqs3.size(); ++i) {
            std::cout << i << ": " << red_allele_freqs3.at(i) << "\n";
        }

        for (unsigned int i = 0; i < red_allele_freqs2.size(); ++i) {
            REQUIRE(red_allele_freqs2.at(i) == Approx(red_allele_freqs3.at(i)).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing missing gene copies (6 of 10 and 8 of 10) with relative root",
        "[ComparisonRelativeRootPopulationTree]") {

    SECTION("Testing missing gene copies (6 of 10 and 8 of 10) for fully fixed pair") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        std::string nex_path2 = "data/hemi129-2-1.nex";
        std::string nex_path3 = "data/hemi129-5-5.nex";
        // Need to keep constant characters
        ComparisonRelativeRootPopulationTree tree2(nex_path2, ' ', true, true, false, false);
        ComparisonRelativeRootPopulationTree tree3(nex_path3, ' ', true, true, false, false);
        tree2.estimate_mutation_rate();
        tree3.estimate_mutation_rate();

        tree2.set_root_population_size(Ne * mu);
        tree3.set_root_population_size(Ne * mu);

        tree2.set_child_population_size(0, (Ne * mu));
        tree3.set_child_population_size(0, (Ne * mu));
        tree2.set_child_population_size(1, (Ne * mu));
        tree3.set_child_population_size(1, (Ne * mu));

        tree2.set_root_height(2.0 * Ne * mu);
        tree3.set_root_height(2.0 * Ne * mu);

        tree2.set_freq_1(0.5);
        tree3.set_freq_1(0.5);

        tree2.set_mutation_rate(1.0);
        tree3.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        unsigned int nsamples = 1000000;

        std::shared_ptr<GeneTreeSimNode> gtree2;
        std::shared_ptr<GeneTreeSimNode> gtree3;
        std::vector<unsigned int> down_sampled_counts = {4, 2};

        std::vector<unsigned int> pop1_red_allele_counts2 = {0, 0, 0, 0, 0};
        std::vector<unsigned int> pop1_red_allele_counts3 = {0, 0, 0, 0, 0};
        std::vector<unsigned int> pop2_red_allele_counts2 = {0, 0, 0};
        std::vector<unsigned int> pop2_red_allele_counts3 = {0, 0, 0};

        for (unsigned int i = 0; i < nsamples; ++i) {
            gtree2 = tree2.simulate_gene_tree(0, rng, true);
            gtree3 = tree3.simulate_gene_tree(0, rng, true);

            auto pattern2 = tree2.simulate_biallelic_site(gtree2, rng);
            auto pattern3 = tree3.simulate_biallelic_site_sans_missing(gtree3, down_sampled_counts, rng);

            ++pop1_red_allele_counts2.at(pattern2.first.at(0));
            ++pop1_red_allele_counts3.at(pattern3.first.at(0));
            ++pop2_red_allele_counts2.at(pattern2.first.at(1));
            ++pop2_red_allele_counts3.at(pattern3.first.at(1));
        }

        unsigned int pop1_total2 = 0;
        unsigned int pop1_total3 = 0;
        for (unsigned int i = 0; i < pop1_red_allele_counts2.size(); ++i) {
            pop1_total2 += pop1_red_allele_counts2.at(i);
            pop1_total3 += pop1_red_allele_counts3.at(i);
        }
        unsigned int pop2_total2 = 0;
        unsigned int pop2_total3 = 0;
        for (unsigned int i = 0; i < pop2_red_allele_counts2.size(); ++i) {
            pop2_total2 += pop2_red_allele_counts2.at(i);
            pop2_total3 += pop2_red_allele_counts3.at(i);
        }

        REQUIRE(pop1_total2 == nsamples);
        REQUIRE(pop1_total3 == nsamples);
        REQUIRE(pop2_total2 == nsamples);
        REQUIRE(pop2_total3 == nsamples);

        std::vector<double> pop1_red_allele_freqs2 = {0.0, 0.0, 0.0, 0.0, 0.0};
        std::vector<double> pop1_red_allele_freqs3 = {0.0, 0.0, 0.0, 0.0, 0.0};
        for (unsigned int i = 0; i < pop1_red_allele_counts2.size(); ++i) {
            pop1_red_allele_freqs2.at(i) = pop1_red_allele_counts2.at(i) / (double)nsamples;
            pop1_red_allele_freqs3.at(i) = pop1_red_allele_counts3.at(i) / (double)nsamples;
        }
        std::vector<double> pop2_red_allele_freqs2 = {0.0, 0.0, 0.0};
        std::vector<double> pop2_red_allele_freqs3 = {0.0, 0.0, 0.0};
        for (unsigned int i = 0; i < pop2_red_allele_counts2.size(); ++i) {
            pop2_red_allele_freqs2.at(i) = pop2_red_allele_counts2.at(i) / (double)nsamples;
            pop2_red_allele_freqs3.at(i) = pop2_red_allele_counts3.at(i) / (double)nsamples;
        }

        std::cout << "Pop 1 freqs of red allele counts for 4 gene copies:\n";
        for (unsigned int i = 0; i < pop1_red_allele_freqs2.size(); ++i) {
            std::cout << i << ": " << pop1_red_allele_freqs2.at(i) << "\n";
        }
        std::cout << "Pop 1 freqs of red allele counts for 10 gene copies with 6 dropped:\n";
        for (unsigned int i = 0; i < pop1_red_allele_freqs3.size(); ++i) {
            std::cout << i << ": " << pop1_red_allele_freqs3.at(i) << "\n";
        }
        std::cout << "Pop 2 freqs of red allele counts for 2 gene copies:\n";
        for (unsigned int i = 0; i < pop2_red_allele_freqs2.size(); ++i) {
            std::cout << i << ": " << pop2_red_allele_freqs2.at(i) << "\n";
        }
        std::cout << "Pop 2 freqs of red allele counts for 10 gene copies with 8 dropped:\n";
        for (unsigned int i = 0; i < pop2_red_allele_freqs3.size(); ++i) {
            std::cout << i << ": " << pop2_red_allele_freqs3.at(i) << "\n";
        }

        for (unsigned int i = 0; i < pop1_red_allele_freqs2.size(); ++i) {
            REQUIRE(pop1_red_allele_freqs2.at(i) == Approx(pop1_red_allele_freqs3.at(i)).epsilon(0.001));
        }
        for (unsigned int i = 0; i < pop2_red_allele_freqs2.size(); ++i) {
            REQUIRE(pop2_red_allele_freqs2.at(i) == Approx(pop2_red_allele_freqs3.at(i)).epsilon(0.001));
        }
    }
}
