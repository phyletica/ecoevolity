#include "catch.hpp"
#include "ecoevolity/tree.hpp"

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
        PopulationTree tree(nex_path, '_', true, true);
        tree.set_root_height(0.01);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-31.77866581319647));
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
        PopulationTree tree(nex_path, '_', true, true, false);
        tree.set_root_height(0.01);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
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
        PopulationTree tree(nex_path, '_', true, false, false);
        tree.set_root_height(0.01);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
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
        PopulationTree tree(nex_path, '_', true, false, true);
        tree.set_root_height(0.01);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6986.120524781545));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        // Can't fold patterns for dominant markers even if u==v, so this
        // should NOT be equal
        REQUIRE(l != Approx(-6986.120524781545));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
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
        PopulationTree tree(nex_path, '_', true, true, false);
        tree.set_root_height(0.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
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
        PopulationTree tree(nex_path, '_', true, false, false);
        tree.set_root_height(0.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
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
        PopulationTree tree(nex_path, '_', true, false, true);
        tree.set_root_height(0.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7223.362711937651));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        // Can't fold patterns for dominant markers even if u==v, so this
        // should NOT be equal
        REQUIRE(l != Approx(-7223.362711937651));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
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
        PopulationTree tree(nex_path, '_', true, true, false);
        tree.set_root_height(0.2);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
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
        PopulationTree tree(nex_path, '_', true, false, false);
        tree.set_root_height(0.2);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
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
        PopulationTree tree(nex_path, '_', true, false, true);
        tree.set_root_height(0.2);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7405.145951634711));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        // Can't fold patterns for dominant markers even if u==v, so this
        // should NOT be equal
        REQUIRE(l != Approx(-7405.145951634711));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
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
        PopulationTree tree(nex_path, '_', true, true, false);
        tree.set_root_height(0.03);
        tree.set_u(10.0);
        tree.set_v(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
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
        PopulationTree tree(nex_path, '_', true, false, false);
        tree.set_root_height(0.03);
        tree.set_u(10.0);
        tree.set_v(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l != Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
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
        PopulationTree tree(nex_path, '_', true, false, true);
        tree.set_root_height(0.03);
        tree.set_u(10.0);
        tree.set_v(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6494.774924871097));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l != Approx(-6494.774924871097));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
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
        PopulationTree tree(nex_path, '_', true, true, false);
        tree.set_root_height(0.03);
        tree.set_u(10.0/19.0);
        tree.set_v(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
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
        PopulationTree tree(nex_path, '_', true, false, false);
        tree.set_root_height(0.03);
        tree.set_u(10.0/19.0);
        tree.set_v(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l != Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
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
        PopulationTree tree(nex_path, '_', true, false, true);
        tree.set_root_height(0.03);
        tree.set_u(10.0/19.0);
        tree.set_v(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-10999.288193543642));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l != Approx(-10999.288193543642));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
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
        PopulationTree tree(nex_path, '_', true, true, false);
        tree.set_root_height(0.03);
        tree.set_u(10.0/19.0);
        tree.set_v(10.0);
        tree.set_coalescence_rate(111.1);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-224.40177558289847));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
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
        PopulationTree tree(nex_path, '_', true, false, false);
        tree.set_root_height(0.03);
        tree.set_u(10.0/19.0);
        tree.set_v(10.0);
        tree.set_coalescence_rate(111.1);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-8158.88094671241));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
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
        PopulationTree tree(nex_path, '_', true, false, true);
        tree.set_root_height(0.03);
        tree.set_u(10.0/19.0);
        tree.set_v(10.0);
        tree.set_coalescence_rate(111.1);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-8034.250341980543));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
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
        ComparisonPopulationTree tree(nex_path, '_', true, true);
        REQUIRE(tree.get_root()->get_label() == "root-pop1");
        tree.set_height(0.01);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-31.77866581319647));
    }
}
