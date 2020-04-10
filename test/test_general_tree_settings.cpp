#include "catch.hpp"
#include "ecoevolity/general_tree_settings.hpp"
#include "ecoevolity/parameter.hpp"


TEST_CASE("Testing basic settings", "[PopSizeSettings]") {
    SECTION("Testing basic settings") {
        std::stringstream ss;
        ss << "value: 1.0\n";
        ss << "estimate: true\n";
        ss << "prior:\n";
        ss << "    gamma_distribution:\n";
        ss << "        shape: 2.0\n";
        ss << "        mean: 0.1\n";

        YAML::Node n;
        n = YAML::Load(ss);
        PopSizeSettings settings(n);

        REQUIRE(settings.population_sizes_are_constrained() == true);
        REQUIRE(settings.get_value() == 1.0);
        REQUIRE(settings.is_fixed() == false);
        std::string s = settings.to_string();
        std::string e = (
                "equal_population_sizes: true\n"
                "value: 1\n"
                "estimate: true\n"
                "prior:\n"
                "    gamma_distribution:\n"
                "        shape: 2\n"
                "        scale: 0.05\n"
                );
        REQUIRE(s == e);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        PositiveRealParameter p = PositiveRealParameter(settings, rng);
        REQUIRE(p.get_value() == 1.0);
        REQUIRE(p.is_fixed() == false);
        REQUIRE(p.get_prior_string() == "gamma(shape = 2, scale = 0.05)");
    }
}

TEST_CASE("Testing explicit pop size constraint", "[PopSizeSettings]") {
    SECTION("Testing explicit constraint") {
        std::stringstream ss;
        ss << "value: 1.0\n";
        ss << "estimate: true\n";
        ss << "equal_population_sizes: true\n";
        ss << "prior:\n";
        ss << "    gamma_distribution:\n";
        ss << "        shape: 2.0\n";
        ss << "        mean: 0.1\n";

        YAML::Node n;
        n = YAML::Load(ss);
        PopSizeSettings settings(n);

        REQUIRE(settings.population_sizes_are_constrained() == true);
        REQUIRE(settings.get_value() == 1.0);
        REQUIRE(settings.is_fixed() == false);
        std::string s = settings.to_string();
        std::string e = (
                "equal_population_sizes: true\n"
                "value: 1\n"
                "estimate: true\n"
                "prior:\n"
                "    gamma_distribution:\n"
                "        shape: 2\n"
                "        scale: 0.05\n"
                );
        REQUIRE(s == e);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        PositiveRealParameter p = PositiveRealParameter(settings, rng);
        REQUIRE(p.get_value() == 1.0);
        REQUIRE(p.is_fixed() == false);
        REQUIRE(p.get_prior_string() == "gamma(shape = 2, scale = 0.05)");
    }
}

TEST_CASE("Testing unconstrained", "[PopSizeSettings]") {
    SECTION("Testing no constraint") {
        std::stringstream ss;
        ss << "value: 1.0\n";
        ss << "estimate: true\n";
        ss << "equal_population_sizes: false\n";
        ss << "prior:\n";
        ss << "    gamma_distribution:\n";
        ss << "        shape: 2.0\n";
        ss << "        mean: 0.1\n";

        YAML::Node n;
        n = YAML::Load(ss);
        PopSizeSettings settings(n);

        REQUIRE(settings.population_sizes_are_constrained() == false);
        REQUIRE(settings.get_value() == 1.0);
        REQUIRE(settings.is_fixed() == false);
        std::string s = settings.to_string();
        std::string e = (
                "equal_population_sizes: false\n"
                "value: 1\n"
                "estimate: true\n"
                "prior:\n"
                "    gamma_distribution:\n"
                "        shape: 2\n"
                "        scale: 0.05\n"
                );
        REQUIRE(s == e);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        PositiveRealParameter p = PositiveRealParameter(settings, rng);
        REQUIRE(p.get_value() == 1.0);
        REQUIRE(p.is_fixed() == false);
        REQUIRE(p.get_prior_string() == "gamma(shape = 2, scale = 0.05)");
    }
}

TEST_CASE("Testing PopulationTreeSettings with uniform_root_and_betas tree prior",
        "[PopulationTreeSettings]") {
    SECTION("Testing with uniform_root_and_betas tree prior") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "---\n"
                   << "tree_model:\n"
                   << "    tree_space: generalized\n"
                   // << "    tree_space:\n"
                   // << "        generalized:\n"
                   // << "    tree_prior:\n"
                   // << "        uniform_root_and_betas:\n"
                   // << "            parameters:\n"
                   // << "                root_height:\n"
                   // << "                    value: 0.3\n"
                   // << "                    estimate: true\n"
                   // << "                    prior:\n"
                   // << "                        gamma_distribution:\n"
                   // << "                            shape: 8.0\n"
                   // << "                            mean: 0.3\n"
                   // << "                alpha_of_node_height_beta_prior:\n"
                   // << "                    value: 1.0\n"
                   // << "                    estimate: true\n"
                   // << "                    prior:\n"
                   // << "                        gamma_distribution:\n"
                   // << "                            shape: 4.0\n"
                   // << "                            mean: 1.0\n"
                   << "data:\n"
                   << "    ploidy: 1\n"
                   << "    constant_sites_removed: false\n"
                   << "    yaml_allele_counts:\n"
                   << "        path: \"diploid-dna-constant-missing.yml\"\n";

        std::cout << "Testing the following config:\n";
        std::cout << cfg_stream.str() << "\n";

        PopulationTreeSettings settings = PopulationTreeSettings(cfg_stream, cfg_path);

        REQUIRE(settings.tree_model_settings.tree_prior->get_type() ==
                EcoevolityOptions::TreePrior::uniform_root_and_betas);
        REQUIRE(settings.tree_model_settings.get_tree_space() ==
                EcoevolityOptions::TreeSpace::generalized);

        settings.write_settings(std::cout);
    }
}
