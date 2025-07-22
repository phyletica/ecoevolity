#include "catch.hpp"
#include <cmath>
#include "ecoevolity/general_tree_settings.hpp"
#include "ecoevolity/rng.hpp"
#include "ecoevolity/seqtree.hpp"

TEST_CASE("Testing no partitioning JC model on ML tree", "[pll]") {
    SECTION("Testing simple JC model on ML tree") {
        std::string config_path = "data/rbcl-ml-jc-config.yml";
        NucTreeAnalysisSettings settings(config_path);
        RandomNumberGenerator rng;
        rng.set_seed(1);

        ecoevolity::SeqTree<Node> tree(settings, rng);

        double lnl = tree.compute_log_likelihood(1);
        double expected_lnl = -278.838;
        double diff = std::fabs(lnl - expected_lnl);
        REQUIRE(diff == Approx(0.0).epsilon(0.0005));
    }
}

TEST_CASE("Testing no partitioning JC model on rooted bifurcating ML tree", "[pll]") {
    SECTION("Testing simple JC model on rooted bifurcating ML tree") {
        std::string config_path = "data/rbcl-ml-rooted-bifurcating-jc-config.yml";
        NucTreeAnalysisSettings settings(config_path);
        RandomNumberGenerator rng;
        rng.set_seed(1);

        ecoevolity::SeqTree<Node> tree(settings, rng);

        double lnl = tree.compute_log_likelihood(1);
        double expected_lnl = -278.838;
        double diff = std::fabs(lnl - expected_lnl);
        REQUIRE(diff == Approx(0.0).epsilon(0.0005));
    }
}

TEST_CASE("Testing no partitioning JC model on ultrametric tree", "[pll]") {
    SECTION("Testing simple JC model on ultrametric tree") {
        std::string config_path = "data/rbcl-ultrametric-jc-config.yml";
        NucTreeAnalysisSettings settings(config_path);
        RandomNumberGenerator rng;
        rng.set_seed(1);

        ecoevolity::SeqTree<Node> tree(settings, rng);

        double lnl = tree.compute_log_likelihood(1);
        double expected_lnl = -340.493;
        double diff = std::fabs(lnl - expected_lnl);
        REQUIRE(diff == Approx(0.0).epsilon(0.0005));
    }
}

TEST_CASE("Testing no partitioning JC model on ultrametric tree with polytomy root", "[pll]") {
    SECTION("Testing simple JC model on ultrametric tree with polytomy root") {
        std::string config_path = "data/rbcl-ultrametric-root-poly-jc-config.yml";
        NucTreeAnalysisSettings settings(config_path);
        RandomNumberGenerator rng;
        rng.set_seed(1);

        ecoevolity::SeqTree<Node> tree(settings, rng);

        double lnl = tree.compute_log_likelihood(1);
        double expected_lnl = -340.493;
        double diff = std::fabs(lnl - expected_lnl);
        REQUIRE(diff == Approx(0.0).epsilon(0.0005));
    }
}

TEST_CASE("Testing no partitioning GTR+G model on ML tree", "[pll]") {
    SECTION("Testing simple GTR+G model on ML tree") {
        std::string config_path = "data/rbcl-ml-gtrg-config.yml";
        NucTreeAnalysisSettings settings(config_path);
        RandomNumberGenerator rng;
        rng.set_seed(1);

        ecoevolity::SeqTree<Node> tree(settings, rng);

        double lnl = tree.compute_log_likelihood(1);
        double expected_lnl = -274.681;
        double diff = std::fabs(lnl - expected_lnl);
        REQUIRE(diff == Approx(0.0).epsilon(0.0005));
    }
}
