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
        PopulationTree tree(nex_path, '_', true, true, false);
        tree.set_root_height(0.03);
        tree.set_u(10.0);
        tree.set_v(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
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
        REQUIRE(! std::isnan(l));
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
        PopulationTree tree(nex_path, '_', true, true, false);
        tree.set_root_height(0.03);
        tree.set_u(10.0/19.0);
        tree.set_v(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
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
        REQUIRE(! std::isnan(l));
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

TEST_CASE("Testing simple prior of PopulationTree", "[PopulationTree]") {

    SECTION("Testing constructor and prior calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        PopulationTree tree(nex_path, '_', true, true);

        tree.set_root_height(0.1);
        tree.set_coalescence_rate(100.0);
        tree.set_u(10.0/19.0);
        tree.set_v(10.0);

        tree.set_node_height_prior(std::make_shared<ExponentialDistribution>(100.0));
        tree.set_population_size_prior(std::make_shared<GammaDistribution>(10.0, 0.0001));
        tree.set_u_prior(std::make_shared<ExponentialDistribution>(10.0));
        tree.set_v_prior(std::make_shared<ExponentialDistribution>(0.1));

        // height   -5.3948298140119091
        // sizes     3 * -155.90663080917298
        // u        -2.9605728017427961
        // v        -3.3025850929940459
        // total    -479.37788013626772
        
        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_coalescence_rates() == Approx(-467.71989242751897));
        REQUIRE(tree.compute_log_prior_density_of_mutation_rates() == Approx(-6.2631578947368425));

        double d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-479.37788013626772));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-479.37788013626772));


        tree.store_prior_density();
        tree.constrain_mutation_rates();

        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_coalescence_rates() == Approx(-467.71989242751897));
        REQUIRE(tree.compute_log_prior_density_of_mutation_rates() == 0.0);

        d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-473.1147222415309));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-473.1147222415309));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-479.37788013626772));


        tree.store_prior_density();
        tree.constrain_coalescence_rates();

        REQUIRE(tree.compute_log_prior_density_of_node_heights() == Approx(-5.3948298140119091));
        REQUIRE(tree.compute_log_prior_density_of_coalescence_rates() == Approx(-155.90663080917298));
        REQUIRE(tree.compute_log_prior_density_of_mutation_rates() == 0.0);

        d = tree.compute_log_prior_density();

        REQUIRE(d == Approx(-161.30146062318488));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-161.30146062318488));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-473.1147222415309));

        tree.restore_prior_density();
        REQUIRE(tree.get_log_prior_density_value() == Approx(-473.1147222415309));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-473.1147222415309));
    }
}




TEST_CASE("Testing hemi129.nex state manipulation", "[ComparisonPopulationTree]") {

    SECTION("Testing state manipulation") {
        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, '_', true, true, false);

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

        tree.set_u_prior(std::make_shared<ExponentialDistribution>(1.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_v_prior(std::make_shared<ExponentialDistribution>(1.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_height(0.01);
        REQUIRE(tree.get_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_coalescence_rate(10.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_child_coalescence_rate(0, 10.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_child_coalescence_rate(1, 10.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_u(1.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_v(1.0);
        REQUIRE(tree.is_dirty());

        REQUIRE(tree.get_number_of_likelihood_calculations() == 0);

        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5807.5500949166799));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 1);

        tree.store_state();
        tree.set_height(0.0);
        REQUIRE(tree.get_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5807.5500949166799));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5807.5500949166799));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 2);

        tree.restore_state();
        REQUIRE(tree.get_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-248.93254688526213));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5807.5500949166799));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5807.5500949166799));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.compute_log_likelihood_and_prior();
        REQUIRE(tree.get_log_likelihood_value() == Approx(-248.93254688526213));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5807.5500949166799));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5807.5500949166799));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.store_state();
        tree.set_height(0.0);
        tree.compute_log_likelihood_and_prior();
        REQUIRE(tree.get_number_of_likelihood_calculations() == 4);


        tree.store_state();
        tree.set_height(0.2);
        REQUIRE(tree.get_height() == 0.2);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-328.39238828878365));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5807.5500949166799));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5807.5500949166799));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 5);


        tree.store_state();
        tree.set_height(0.03);
        tree.set_u(10.0);
        tree.set_v(10.0/19.0);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_u() == 10.0);
        REQUIRE(tree.get_v() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_coalescence_rate() == 10.0);
        REQUIRE(tree.get_child_coalescence_rate(0) == 10.0);
        REQUIRE(tree.get_child_coalescence_rate(1) == 10.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-327.7437811413033));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5816.0764107061532));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5807.5500949166799));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 6);


        tree.store_state();
        tree.set_u(10.0/19.0);
        tree.set_v(10.0);
        REQUIRE(tree.get_height() == 0.03);
        REQUIRE(tree.get_v() == 10.0);
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-265.0023534261969));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-327.7437811413033));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5816.0764107061532));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5816.0764107061532));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 7);


        tree.store_state();
        tree.set_coalescence_rate(111.1);
        REQUIRE(tree.get_height() == 0.03);
        REQUIRE(tree.get_v() == 10.0);
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_coalescence_rate() == 111.1);
        REQUIRE(tree.get_child_coalescence_rate(0) == 111.1);
        REQUIRE(tree.get_child_coalescence_rate(1) == 111.1);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-224.40177558289847));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-265.0023534261969));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-421.14224740528493));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5816.0764107061532));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 8);

        tree.store_state();
        tree.constrain_mutation_rates();
        tree.set_coalescence_rate(10.0);
        tree.set_height(0.2);
        REQUIRE(tree.get_height() == 0.2);
        REQUIRE(tree.get_v() == 1.0);
        REQUIRE(tree.get_u() == 1.0);
        REQUIRE(tree.get_root_coalescence_rate() == 10.0);
        REQUIRE(tree.get_child_coalescence_rate(0) == 10.0);
        REQUIRE(tree.get_child_coalescence_rate(1) == 10.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-224.40177558289847));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5805.5500949166799));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-421.14224740528493));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 9);


        tree.store_state();
        tree.fold_patterns();
        REQUIRE(tree.get_height() == 0.2);
        REQUIRE(tree.get_v() == 1.0);
        REQUIRE(tree.get_u() == 1.0);
        REQUIRE(tree.get_root_coalescence_rate() == 10.0);
        REQUIRE(tree.get_child_coalescence_rate(0) == 10.0);
        REQUIRE(tree.get_child_coalescence_rate(1) == 10.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5805.5500949166799));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5805.5500949166799));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 10);


        tree.store_state();
        tree.constrain_coalescence_rates();
        REQUIRE(tree.get_height() == 0.2);
        REQUIRE(tree.get_v() == 1.0);
        REQUIRE(tree.get_u() == 1.0);
        REQUIRE(tree.get_node_height_multiplier() == 1.0);
        REQUIRE(tree.get_root_coalescence_rate() == 10.0);
        REQUIRE(tree.get_child_coalescence_rate(0) == 10.0);
        REQUIRE(tree.get_child_coalescence_rate(1) == 10.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1935.1833649722266));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5805.5500949166799));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 11);

        tree.store_state();
        tree.estimate_node_height_multiplier();
        REQUIRE(! tree.node_height_multiplier_is_fixed());
        tree.set_node_height_multiplier(1.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());
        tree.set_node_height_multiplier_prior(std::make_shared<GammaDistribution>(10.0, 0.1));
        REQUIRE(tree.get_height() == 0.2);
        REQUIRE(tree.get_v() == 1.0);
        REQUIRE(tree.get_u() == 1.0);
        REQUIRE(tree.get_node_height_multiplier() == 1.0);
        REQUIRE(tree.get_root_coalescence_rate() == 10.0);
        REQUIRE(tree.get_child_coalescence_rate(0) == 10.0);
        REQUIRE(tree.get_child_coalescence_rate(1) == 10.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1934.9593415223676));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1935.1833649722266));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 12);

        tree.store_state();
        tree.set_node_height_multiplier(2.0);
        tree.set_height(0.1);
        REQUIRE(tree.get_height() == 0.1);
        REQUIRE(tree.get_v() == 1.0);
        REQUIRE(tree.get_u() == 1.0);
        REQUIRE(tree.get_node_height_multiplier() == 2.0);
        REQUIRE(tree.get_root_coalescence_rate() == 10.0);
        REQUIRE(tree.get_child_coalescence_rate(0) == 10.0);
        REQUIRE(tree.get_child_coalescence_rate(1) == 10.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1938.7210168973281));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1934.9593415223676));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        tree.restore_state();
        REQUIRE(tree.get_height() == 0.2);
        REQUIRE(tree.get_v() == 1.0);
        REQUIRE(tree.get_u() == 1.0);
        REQUIRE(tree.get_node_height_multiplier() == 1.0);
        REQUIRE(tree.get_root_coalescence_rate() == 10.0);
        REQUIRE(tree.get_child_coalescence_rate(0) == 10.0);
        REQUIRE(tree.get_child_coalescence_rate(1) == 10.0);
        REQUIRE(tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1934.9593415223676));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1934.9593415223676));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);
    }
}

TEST_CASE("Testing hemi129.nex state manipulation for PopulationTree", "[PopulationTree]") {

    SECTION("Testing state manipulation") {
        std::string nex_path = "data/hemi129.nex";
        PopulationTree tree(nex_path, '_', true, true, false);

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

        tree.set_u_prior(std::make_shared<ExponentialDistribution>(1.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_v_prior(std::make_shared<ExponentialDistribution>(1.0));

        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_root_height(0.01);
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_coalescence_rate(10.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_u(1.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());

        tree.set_v(1.0);
        REQUIRE(tree.is_dirty());

        REQUIRE(tree.get_number_of_likelihood_calculations() == 0);

        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5805.3475098236859));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 1);

        tree.store_state();
        tree.set_root_height(0.0);
        REQUIRE(tree.get_root_height() == 0.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-328.39238828878365));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5805.2475098236855));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5805.3475098236859));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 2);

        tree.restore_state();
        REQUIRE(tree.get_root_height() == 0.01);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-248.93254688526213));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5805.3475098236859));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5805.3475098236859));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 3);

        tree.compute_log_likelihood_and_prior();
        REQUIRE(tree.get_log_likelihood_value() == Approx(-248.93254688526213));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-248.93254688526213));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5805.3475098236859));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5805.3475098236859));
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
        REQUIRE(tree.get_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-328.39238828878365));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5807.2475098236855));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5805.2475098236855));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 5);


        tree.store_state();
        tree.set_root_height(0.03);
        tree.set_u(10.0);
        tree.set_v(10.0/19.0);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_u() == 10.0);
        REQUIRE(tree.get_v() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_coalescence_rate() == 10.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-327.7437811413033));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5814.073825613159));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5807.2475098236855));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 6);


        tree.store_state();
        tree.set_u(10.0/19.0);
        tree.set_v(10.0);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == 10.0);
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-265.0023534261969));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-327.7437811413033));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5814.073825613159));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5814.073825613159));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 7);


        tree.store_state();
        tree.set_coalescence_rate(111.1);
        REQUIRE(tree.get_root_height() == 0.03);
        REQUIRE(tree.get_v() == 10.0);
        REQUIRE(tree.get_u() == Approx(10.0/19.0));
        REQUIRE(tree.get_root_coalescence_rate() == 111.1);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-224.40177558289847));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-265.0023534261969));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-419.13966231229091));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5814.073825613159));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 8);

        tree.store_state();
        tree.constrain_mutation_rates();
        tree.set_coalescence_rate(10.0);
        tree.set_root_height(0.2);
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == 1.0);
        REQUIRE(tree.get_u() == 1.0);
        REQUIRE(tree.get_root_coalescence_rate() == 10.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-224.40177558289847));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5805.2475098236855));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-419.13966231229091));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 9);


        tree.store_state();
        tree.fold_patterns();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == 1.0);
        REQUIRE(tree.get_u() == 1.0);
        REQUIRE(tree.get_root_coalescence_rate() == 10.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-5805.2475098236855));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5805.2475098236855));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 10);


        tree.store_state();
        tree.constrain_coalescence_rates();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == 1.0);
        REQUIRE(tree.get_u() == 1.0);
        REQUIRE(tree.get_node_height_multiplier() == 1.0);
        REQUIRE(tree.get_root_coalescence_rate() == 10.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1934.8807798792325));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-5805.2475098236855));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 11);

        tree.store_state();
        tree.estimate_node_height_multiplier();
        REQUIRE(! tree.node_height_multiplier_is_fixed());
        tree.set_node_height_multiplier(1.0);
        REQUIRE(tree.is_dirty());
        tree.make_clean();
        REQUIRE(! tree.is_dirty());
        tree.set_node_height_multiplier_prior(std::make_shared<GammaDistribution>(10.0, 0.1));
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == 1.0);
        REQUIRE(tree.get_u() == 1.0);
        REQUIRE(tree.get_node_height_multiplier() == 1.0);
        REQUIRE(tree.get_root_coalescence_rate() == 10.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1934.6567564293734));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1934.8807798792325));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 12);

        tree.store_state();
        tree.set_node_height_multiplier(2.0);
        tree.set_root_height(0.1);
        REQUIRE(tree.get_root_height() == 0.1);
        REQUIRE(tree.get_v() == 1.0);
        REQUIRE(tree.get_u() == 1.0);
        REQUIRE(tree.get_node_height_multiplier() == 2.0);
        REQUIRE(tree.get_root_coalescence_rate() == 10.0);
        REQUIRE(tree.is_dirty());
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1937.418431804334));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1934.6567564293734));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);

        tree.restore_state();
        REQUIRE(tree.get_root_height() == 0.2);
        REQUIRE(tree.get_v() == 1.0);
        REQUIRE(tree.get_u() == 1.0);
        REQUIRE(tree.get_node_height_multiplier() == 1.0);
        REQUIRE(tree.get_root_coalescence_rate() == 10.0);
        REQUIRE(tree.is_dirty());
        REQUIRE(tree.get_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_stored_log_likelihood_value() == Approx(-227.41048391087554));
        REQUIRE(tree.get_log_prior_density_value() == Approx(-1934.6567564293734));
        REQUIRE(tree.get_stored_log_prior_density_value() == Approx(-1934.6567564293734));
        REQUIRE(tree.get_number_of_likelihood_calculations() == 13);
    }
}
