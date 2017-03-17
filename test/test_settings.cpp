#include "catch.hpp"
#include "ecoevolity/settings.hpp"
#include "ecoevolity/parameter.hpp"
#include "ecoevolity/tree.hpp"

TEST_CASE("Testing bare constructor of ContinuousDistributionSettings",
        "[ContinuousDistributionSettings]") {

    SECTION("Testing bare") {
        ContinuousDistributionSettings cds;
        REQUIRE(cds.get_name() == "none");
        REQUIRE(cds.to_string() == "");
        std::shared_ptr<ContinuousProbabilityDistribution> p = cds.get_instance();
        REQUIRE(! p);

        cds = ContinuousDistributionSettings();
        REQUIRE(cds.get_name() == "none");
        REQUIRE(cds.to_string() == "");
        std::shared_ptr<ContinuousProbabilityDistribution> p2 = cds.get_instance();
        REQUIRE(! p2);
    }
}

TEST_CASE("Testing invalid distribution in ContinuousDistributionSettings constructor",
        "[ContinuousDistributionSettings]") {

    SECTION("Testing invalid distribution") {
        std::string name = "blah_distribution";
        std::unordered_map<std::string, double> parameters;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }
}

TEST_CASE("Testing missing gamma parameter settings", "[ContinuousDistributionSettings]") {

    SECTION("Testing no parameters") {
        std::string name = "gamma_distribution";
        std::unordered_map<std::string, double> parameters;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }

    SECTION("Testing missing parameters") {
        std::string name = "gamma_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["shape"] = 1.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
        parameters["offset"] = 0.01;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);

        parameters.clear();
        parameters["scale"] = 1.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
        parameters["offset"] = 0.01;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }

    SECTION("Testing rate parameter") {
        std::string name = "gamma_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["shape"] = 10.0;
        parameters["rate"] = 1.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }
}

TEST_CASE("Testing missing exponential parameter settings", "[ContinuousDistributionSettings]") {

    SECTION("Testing no parameters") {
        std::string name = "exponential_distribution";
        std::unordered_map<std::string, double> parameters;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }

    SECTION("Testing missing parameters") {
        std::string name = "exponential_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["offset"] = 0.01;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }

    SECTION("Testing lambda parameter") {
        std::string name = "exponential_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["lambda"] = 1.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }
}

TEST_CASE("Testing gamma parameter settings errors", "[ContinuousDistributionSettings]") {
    SECTION("Testing gamma errors") {
        std::string name = "gamma_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["shape"] = 0.0;
        parameters["scale"] = 1.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);

        parameters.clear();
        parameters["shape"] = 10.0;
        parameters["scale"] = 0.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }
}

TEST_CASE("Testing exponential parameter settings errors", "[ContinuousDistributionSettings]") {
    SECTION("Testing exponential errors") {
        std::string name = "exponential_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["rate"] = 0.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }
}

TEST_CASE("Testing gamma settings to string", "[ContinuousDistributionSettings]") {
    SECTION("Testing gamma to string") {
        std::string name = "gamma_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["shape"] = 10.0;
        parameters["scale"] = 1.0;
        ContinuousDistributionSettings settings = ContinuousDistributionSettings(name, parameters);

        REQUIRE(settings.get_name() == "gamma_distribution");
        std::string s = settings.to_string(0);
        std::string e = "gamma_distribution:\n    shape: 10\n    scale: 1\n";
        REQUIRE(s == e);

        parameters.clear();
        parameters["shape"] = 10.0;
        parameters["scale"] = 1.0;
        parameters["offset"] = 0.01;
        settings = ContinuousDistributionSettings(name, parameters);

        REQUIRE(settings.get_name() == "gamma_distribution");
        s = settings.to_string(0);
        e = "gamma_distribution:\n    shape: 10\n    scale: 1\n    offset: 0.01\n";
        REQUIRE(s == e);

        parameters.clear();
        parameters["shape"] = 0.12345678910;
        parameters["scale"] = 1.00001;
        parameters["offset"] = 0.001;
        settings = ContinuousDistributionSettings(name, parameters);

        REQUIRE(settings.get_name() == "gamma_distribution");
        s = settings.to_string(0);
        e = "gamma_distribution:\n    shape: 0.123457\n    scale: 1.00001\n    offset: 0.001\n";
        REQUIRE(s == e);
    }
}

TEST_CASE("Testing exponential settings to string", "[ContinuousDistributionSettings]") {
    SECTION("Testing exponential to string") {
        std::string name = "exponential_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["rate"] = 1.0;
        ContinuousDistributionSettings settings = ContinuousDistributionSettings(name, parameters);

        REQUIRE(settings.get_name() == "exponential_distribution");
        std::string s = settings.to_string(0);
        std::string e = "exponential_distribution:\n    rate: 1\n";
        REQUIRE(s == e);

        parameters.clear();
        parameters["rate"] = 1.0;
        parameters["offset"] = 0.01;
        settings = ContinuousDistributionSettings(name, parameters);

        REQUIRE(settings.get_name() == "exponential_distribution");
        s = settings.to_string(0);
        e = "exponential_distribution:\n    rate: 1\n    offset: 0.01\n";
        REQUIRE(s == e);

        parameters.clear();
        parameters["rate"] = 0.12345678910;
        parameters["offset"] = 0.001;
        settings = ContinuousDistributionSettings(name, parameters);

        REQUIRE(settings.get_name() == "exponential_distribution");
        s = settings.to_string(0);
        e = "exponential_distribution:\n    rate: 0.123457\n    offset: 0.001\n";
        REQUIRE(s == e);
    }
}

TEST_CASE("Testing gamma settings to instance", "[ContinuousDistributionSettings]") {
    SECTION("Testing gamma to instance") {
        std::string name = "gamma_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["shape"] = 10.0;
        parameters["scale"] = 1.0;
        ContinuousDistributionSettings settings = ContinuousDistributionSettings(name, parameters);

        std::shared_ptr<ContinuousProbabilityDistribution> g;
        g = settings.get_instance();

        REQUIRE(g->get_min() == 0.0);
        REQUIRE(g->to_string() == "gamma(shape = 10, scale = 1)");


        parameters.clear();
        parameters["shape"] = 10.0;
        parameters["scale"] = 1.0;
        parameters["offset"] = 0.01;
        settings = ContinuousDistributionSettings(name, parameters);

        g = settings.get_instance();

        REQUIRE(g->to_string() == "gamma(shape = 10, scale = 1, offset = 0.01)");
    }
}

TEST_CASE("Testing exponential settings to instance", "[ContinuousDistributionSettings]") {
    SECTION("Testing exponential to instance") {
        std::string name = "exponential_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["rate"] = 1.0;
        ContinuousDistributionSettings settings = ContinuousDistributionSettings(name, parameters);

        std::shared_ptr<ContinuousProbabilityDistribution> g;
        g = settings.get_instance();

        REQUIRE(g->get_min() == 0.0);
        REQUIRE(g->to_string() == "exp(lambda = 1)");


        parameters["rate"] = 10.0;
        parameters["offset"] = 0.01;
        settings = ContinuousDistributionSettings(name, parameters);

        g = settings.get_instance();

        REQUIRE(g->to_string() == "exp(lambda = 10, offset = 0.01)");
    }
}

TEST_CASE("Testing uniform parameter settings errors", "[ContinuousDistributionSettings]") {
    SECTION("Testing uniform errors") {
        std::string name = "uniform_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["min"] = 1.0;
        parameters["max"] = 0.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }
}

TEST_CASE("Testing uniform settings to string", "[ContinuousDistributionSettings]") {
    SECTION("Testing uniform to string") {
        std::string name = "uniform_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["min"] = 1.0;
        parameters["max"] = 2.0;
        ContinuousDistributionSettings settings = ContinuousDistributionSettings(name, parameters);

        REQUIRE(settings.get_name() == "uniform_distribution");
        std::string s = settings.to_string(0);
        std::string e = "uniform_distribution:\n    min: 1\n    max: 2\n";
        REQUIRE(s == e);
    }
}

TEST_CASE("Testing uniform settings to instance", "[ContinuousDistributionSettings]") {
    SECTION("Testing uniform to instance") {
        std::string name = "uniform_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["min"] = 1.0;
        parameters["max"] = 2.0;
        ContinuousDistributionSettings settings = ContinuousDistributionSettings(name, parameters);

        std::shared_ptr<ContinuousProbabilityDistribution> g;
        g = settings.get_instance();

        REQUIRE(g->get_min() == 1.0);
        REQUIRE(g->to_string() == "uniform(1, 2)");
    }
}

TEST_CASE("Testing missing uniform parameter settings", "[ContinuousDistributionSettings]") {

    SECTION("Testing no parameters") {
        std::string name = "uniform_distribution";
        std::unordered_map<std::string, double> parameters;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }

    SECTION("Testing missing parameters") {
        std::string name = "uniform_distribution";
        std::unordered_map<std::string, double> parameters;
        parameters["min"] = 0.01;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);

        parameters.clear();
        parameters["max"] = 1.0;
        REQUIRE_THROWS_AS(ContinuousDistributionSettings(name, parameters),
                EcoevolityContinuousDistributionSettingError);
    }
}

TEST_CASE("Testing error parameter settings", "[PositiveRealParameterSettings]") {
    SECTION("Testing negative parameter") {
        std::unordered_map<std::string, double> prior_parameters;
        prior_parameters["shape"] = 1.0;
        prior_parameters["scale"] = 1.0;

        REQUIRE_THROWS_AS(PositiveRealParameterSettings(
                -0.001,
                true,
                "gamma_distribution",
                prior_parameters),
                EcoevolityPositiveRealParameterSettingError);
    }
}

TEST_CASE("Testing fixed parameter settings", "[PositiveRealParameterSettings]") {
    SECTION("Testing fixed parameter") {
        std::unordered_map<std::string, double> prior_parameters;
        prior_parameters["shape"] = 1.0;
        prior_parameters["scale"] = 1.0;

        PositiveRealParameterSettings settings = PositiveRealParameterSettings(
                0.001,
                true,
                "gamma_distribution",
                prior_parameters);

        REQUIRE(settings.get_value() == 0.001);
        REQUIRE(settings.is_fixed() == true);
        std::string s = settings.to_string();
        std::string e = "value: 0.001\nestimate: false\n";
        REQUIRE(s == e);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        PositiveRealParameter p = PositiveRealParameter(settings, rng);
        REQUIRE_THROWS_AS(p.check_prior(), EcoevolityNullPointerError);
        REQUIRE(p.get_value() == 0.001);
        REQUIRE(p.is_fixed() == true);
    }
}

TEST_CASE("Testing free parameter settings", "[PositiveRealParameterSettings]") {
    SECTION("Testing free parameter") {
        std::unordered_map<std::string, double> prior_parameters;
        prior_parameters["shape"] = 1.0;
        prior_parameters["scale"] = 1.0;

        PositiveRealParameterSettings settings = PositiveRealParameterSettings(
                0.001,
                false,
                "gamma_distribution",
                prior_parameters);

        REQUIRE(settings.get_value() == 0.001);
        REQUIRE(settings.is_fixed() == false);
        std::string s = settings.to_string();
        std::string e = "value: 0.001\nestimate: true\nprior:\n    gamma_distribution:\n        shape: 1\n        scale: 1\n";
        REQUIRE(s == e);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        PositiveRealParameter p = PositiveRealParameter(settings, rng);
        REQUIRE(p.get_value() == 0.001);
        REQUIRE(p.is_fixed() == false);
        REQUIRE(p.get_prior_min() == 0.0);
        REQUIRE(p.get_prior_mean() == 1.0);
        REQUIRE(p.get_prior_variance() == 1.0);
        REQUIRE(p.get_prior_string() == "gamma(shape = 1, scale = 1)");
    }
}

TEST_CASE("Testing nan parameter settings", "[PositiveRealParameterSettings]") {
    SECTION("Testing nan parameter") {
        std::unordered_map<std::string, double> prior_parameters;
        prior_parameters["shape"] = 1.0;
        prior_parameters["scale"] = 1.0;

        PositiveRealParameterSettings settings = PositiveRealParameterSettings(
                std::numeric_limits<double>::quiet_NaN(),
                false,
                "gamma_distribution",
                prior_parameters);

        REQUIRE(std::isnan(settings.get_value()));
        REQUIRE(settings.is_fixed() == false);
        std::string s = settings.to_string();
        std::string e = "estimate: true\nprior:\n    gamma_distribution:\n        shape: 1\n        scale: 1\n";
        REQUIRE(s == e);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        PositiveRealParameter p = PositiveRealParameter(settings, rng);
        REQUIRE(p.is_fixed() == false);
        REQUIRE(p.get_prior_min() == 0.0);
        REQUIRE(p.get_prior_mean() == 1.0);
        REQUIRE(p.get_prior_variance() == 1.0);
        REQUIRE(p.get_prior_string() == "gamma(shape = 1, scale = 1)");
    }
}

TEST_CASE("Testing fixed nan parameter settings", "[PositiveRealParameterSettings]") {
    SECTION("Testing fixed nan parameter") {
        std::unordered_map<std::string, double> prior_parameters;
        prior_parameters["shape"] = 1.0;
        prior_parameters["scale"] = 1.0;

        REQUIRE_THROWS_AS(PositiveRealParameterSettings(
                std::numeric_limits<double>::quiet_NaN(),
                true,
                "gamma_distribution",
                prior_parameters),
                EcoevolityPositiveRealParameterSettingError);
    }
}

TEST_CASE("Testing comparison setting constructor", "[ComparisonSettings]") {
    SECTION("Testing data/diploid-standard-data-ntax5-nchar5.nex") {

        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";

        std::unordered_map<std::string, double> no_prior_parameters;
        std::unordered_map<std::string, double> pop_prior_parameters;
        pop_prior_parameters["shape"] = 10.0;
        pop_prior_parameters["scale"] = 0.0001;
        PositiveRealParameterSettings pop_size = PositiveRealParameterSettings(
                std::numeric_limits<double>::quiet_NaN(),
                false,
                "gamma_distribution",
                pop_prior_parameters);
        PositiveRealParameterSettings u = PositiveRealParameterSettings(
                1.0,
                true,
                "none",
                no_prior_parameters);
        PositiveRealParameterSettings mutation_rate = PositiveRealParameterSettings(
                1.0,
                true,
                "none",
                no_prior_parameters);

        ComparisonSettings settings = ComparisonSettings(
                nex_path,
                pop_size,
                u,
                mutation_rate,
                ' ',
                true,
                true,
                false,
                true,
                true,
                false,
                false);

        std::string s = settings.to_string(0);
        std::string e =  "";
        e += "path: data/diploid-standard-data-ntax5-nchar5.nex\n";
        e += "ploidy: 2\n";
        e += "genotypes_are_diploid: true\n";
        e += "markers_are_dominant: false\n";
        e += "population_name_delimiter: ' '\n";
        e += "population_name_is_prefix: true\n";
        e += "constant_sites_removed: true\n";
        e += "use_empirical_starting_value_for_freq_1: true\n";
        e += "equal_population_sizes: false\n";
        e += "equal_state_frequencies: false\n";
        e += "parameters:\n";
        e += "    population_size:\n";
        e += "        estimate: true\n";
        e += "        prior:\n";
        e += "            gamma_distribution:\n";
        e += "                shape: 10\n";
        e += "                scale: 0.0001\n";
        e += "    mutation_rate:\n";
        e += "        value: 1\n";
        e += "        estimate: false\n";
        e += "    freq_1:\n";
        e += "        value: 0.416667\n";
        e += "        estimate: false\n";
        REQUIRE(s == e);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        ComparisonPopulationTree t = ComparisonPopulationTree(settings, rng);
        REQUIRE(t.get_degree_of_root() == 2);

        REQUIRE(! std::isnan(t.get_root_population_size()));
        REQUIRE(! std::isnan(t.get_child_population_size(0)));
        REQUIRE(! std::isnan(t.get_child_population_size(1)));
        REQUIRE(t.get_u() == Approx(1.2));
        REQUIRE(t.get_v() == Approx(0.8571429));
        REQUIRE(t.get_freq_1() == Approx(0.4166666666666667));
        REQUIRE(t.get_freq_0() == Approx(1.0 - 0.4166666666666667));
        REQUIRE(t.get_freq_1_parameter()->is_fixed() == true);
        REQUIRE(t.get_root_population_size_parameter()->is_fixed() == false);
        REQUIRE(t.get_child_population_size_parameter(0)->is_fixed() == false);
        REQUIRE(t.get_child_population_size_parameter(1)->is_fixed() == false);
        REQUIRE(t.get_root_population_size_parameter() != t.get_child_population_size_parameter(0));
        REQUIRE(t.get_root_population_size_parameter() != t.get_child_population_size_parameter(1));
        REQUIRE(t.get_child_population_size_parameter(0) != t.get_child_population_size_parameter(1));
        REQUIRE(t.get_mutation_rate_parameter()->is_fixed() == true);
        REQUIRE(t.state_frequencies_are_fixed() == true);
        REQUIRE(t.population_sizes_are_fixed() == false);
        REQUIRE(t.mutation_rate_is_fixed() == true);

        REQUIRE(! t.get_freq_1_prior());
        REQUIRE(t.get_population_size_prior()->to_string() == "gamma(shape = 10, scale = 0.0001)");
        REQUIRE(! t.get_mutation_rate_prior());
    }

    SECTION("Testing data/diploid-standard-data-ntax5-nchar5.nex likelihood") {

        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";

        std::unordered_map<std::string, double> no_prior_parameters;
        std::unordered_map<std::string, double> freq_1_prior_parameters;
        std::unordered_map<std::string, double> size_prior_parameters;
        freq_1_prior_parameters["alpha"] = 2.0;
        freq_1_prior_parameters["beta"] = 2.0;
        size_prior_parameters["shape"] = 10.0;
        size_prior_parameters["scale"] = 0.1;
        PositiveRealParameterSettings pop_size = PositiveRealParameterSettings(
                2.0 / 10.0,
                true,
                "gamma_distribution",
                size_prior_parameters);
        PositiveRealParameterSettings freq_1 = PositiveRealParameterSettings(
                0.5,
                false,
                "beta_distribution",
                freq_1_prior_parameters);
        PositiveRealParameterSettings mutation_rate = PositiveRealParameterSettings(
                1.0,
                true,
                "none",
                no_prior_parameters);

        ComparisonSettings settings = ComparisonSettings(
                nex_path,
                pop_size,
                freq_1,
                mutation_rate,
                ' ',
                true,
                true,
                false,
                true,
                true,
                false,
                true);

        std::string s = settings.to_string(0);
        std::string e =  "";
        e += "path: data/diploid-standard-data-ntax5-nchar5.nex\n";
        e += "ploidy: 2\n";
        e += "genotypes_are_diploid: true\n";
        e += "markers_are_dominant: false\n";
        e += "population_name_delimiter: ' '\n";
        e += "population_name_is_prefix: true\n";
        e += "constant_sites_removed: true\n";
        e += "use_empirical_starting_value_for_freq_1: false\n";
        e += "equal_population_sizes: false\n";
        e += "equal_state_frequencies: true\n";
        e += "parameters:\n";
        e += "    population_size:\n";
        e += "        value: 0.2\n";
        e += "        estimate: false\n";
        e += "    mutation_rate:\n";
        e += "        value: 1\n";
        e += "        estimate: false\n";
        e += "    freq_1:\n";
        e += "        value: 0.5\n";
        e += "        estimate: false\n";
        REQUIRE(s == e);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        ComparisonPopulationTree t = ComparisonPopulationTree(settings, rng);
        REQUIRE(t.get_degree_of_root() == 2);

        REQUIRE(t.get_root_population_size() == Approx(2.0/10.0));
        REQUIRE(t.get_child_population_size(0) == Approx(2.0/10.0));
        REQUIRE(t.get_child_population_size(1) == Approx(2.0/10.0));
        REQUIRE(t.get_u() == Approx(1.0));
        REQUIRE(t.get_v() == Approx(1.0));
        REQUIRE(t.get_freq_1() == Approx(0.5));
        REQUIRE(t.get_freq_0() == Approx(0.5));
        REQUIRE(t.get_freq_1_parameter()->is_fixed() == true);
        REQUIRE(t.get_root_population_size_parameter()->is_fixed() == true);
        REQUIRE(t.get_child_population_size_parameter(0)->is_fixed() == true);
        REQUIRE(t.get_child_population_size_parameter(1)->is_fixed() == true);
        REQUIRE(t.get_root_population_size_parameter() != t.get_child_population_size_parameter(0));
        REQUIRE(t.get_root_population_size_parameter() != t.get_child_population_size_parameter(1));
        REQUIRE(t.get_child_population_size_parameter(0) != t.get_child_population_size_parameter(1));
        REQUIRE(t.get_mutation_rate_parameter()->is_fixed() == true);
        REQUIRE(t.state_frequencies_are_fixed() == true);
        REQUIRE(t.population_sizes_are_fixed() == true);
        REQUIRE(t.mutation_rate_is_fixed() == true);

        REQUIRE(! t.get_freq_1_prior());
        REQUIRE(! t.get_population_size_prior());
        REQUIRE(! t.get_mutation_rate_prior());

        t.set_root_height(0.01);
        t.estimate_population_sizes();
        t.set_population_size(2.0 / (10.0 * 4.0));
        double l = t.compute_log_likelihood();
        REQUIRE(l == Approx(-31.77866581319647));
        REQUIRE(t.get_degree_of_root() == 2);
    }

    SECTION("Testing data/diploid-standard-data-ntax5-nchar5.nex likelihood with constrained sizes") {

        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";

        std::unordered_map<std::string, double> no_prior_parameters;
        std::unordered_map<std::string, double> pop_prior_parameters;
        pop_prior_parameters["rate"] = 1.0;
        std::unordered_map<std::string, double> freq_1_prior_parameters;
        freq_1_prior_parameters["alpha"] = 2.0;
        freq_1_prior_parameters["beta"] = 4.0;
        std::unordered_map<std::string, double> m_prior_parameters;
        m_prior_parameters["min"] = 0.9;
        m_prior_parameters["max"] = 1.1;
        PositiveRealParameterSettings pop_size = PositiveRealParameterSettings(
                2.0 / 10.0,
                false,
                "exponential_distribution",
                pop_prior_parameters);
        PositiveRealParameterSettings freq_1 = PositiveRealParameterSettings(
                0.5,
                false,
                "beta_distribution",
                freq_1_prior_parameters);
        PositiveRealParameterSettings mutation_rate = PositiveRealParameterSettings(
                1.0,
                false,
                "uniform_distribution",
                m_prior_parameters);

        ComparisonSettings settings = ComparisonSettings(
                nex_path,
                pop_size,
                freq_1,
                mutation_rate,
                ' ',
                true,
                true,
                false,
                true,
                false,
                true,
                false);

        std::string s = settings.to_string(0);
        std::string e =  "";
        e += "path: data/diploid-standard-data-ntax5-nchar5.nex\n";
        e += "ploidy: 2\n";
        e += "genotypes_are_diploid: true\n";
        e += "markers_are_dominant: false\n";
        e += "population_name_delimiter: ' '\n";
        e += "population_name_is_prefix: true\n";
        e += "constant_sites_removed: true\n";
        e += "use_empirical_starting_value_for_freq_1: false\n";
        e += "equal_population_sizes: true\n";
        e += "equal_state_frequencies: false\n";
        e += "parameters:\n";
        e += "    population_size:\n";
        e += "        value: 0.2\n";
        e += "        estimate: true\n";
        e += "        prior:\n";
        e += "            exponential_distribution:\n";
        e += "                rate: 1\n";
        e += "    mutation_rate:\n";
        e += "        value: 1\n";
        e += "        estimate: true\n";
        e += "        prior:\n";
        e += "            uniform_distribution:\n";
        e += "                min: 0.9\n";
        e += "                max: 1.1\n";
        e += "    freq_1:\n";
        e += "        value: 0.5\n";
        e += "        estimate: true\n";
        e += "        prior:\n";
        e += "            beta_distribution:\n";
        e += "                alpha: 2\n";
        e += "                beta: 4\n";
        REQUIRE(s == e);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        ComparisonPopulationTree t = ComparisonPopulationTree(settings, rng);
        REQUIRE(t.get_degree_of_root() == 2);

        REQUIRE(t.get_root_population_size() == Approx(2.0/10.0));
        REQUIRE(t.get_child_population_size(0) == Approx(2.0/10.0));
        REQUIRE(t.get_child_population_size(1) == Approx(2.0/10.0));
        REQUIRE(t.get_u() == Approx(1.0));
        REQUIRE(t.get_v() == Approx(1.0));
        REQUIRE(t.get_freq_1() == Approx(0.5));
        REQUIRE(t.get_freq_0() == Approx(0.5));
        REQUIRE(t.get_freq_1_parameter()->is_fixed() == false);
        REQUIRE(t.get_root_population_size_parameter()->is_fixed() == false);
        REQUIRE(t.get_child_population_size_parameter(0)->is_fixed() == false);
        REQUIRE(t.get_child_population_size_parameter(1)->is_fixed() == false);
        REQUIRE(t.get_root_population_size_parameter() == t.get_child_population_size_parameter(0));
        REQUIRE(t.get_root_population_size_parameter() == t.get_child_population_size_parameter(1));
        REQUIRE(t.get_child_population_size_parameter(0) == t.get_child_population_size_parameter(1));
        REQUIRE(t.get_mutation_rate_parameter()->is_fixed() == false);
        REQUIRE(t.state_frequencies_are_fixed() == false);
        REQUIRE(t.population_sizes_are_fixed() == false);
        REQUIRE(t.mutation_rate_is_fixed() == false);

        REQUIRE(t.get_freq_1_prior()->to_string() == "beta(2, 4)");
        REQUIRE(t.get_population_size_prior()->to_string() == "exp(lambda = 1)");
        REQUIRE(t.get_mutation_rate_prior()->to_string() == "uniform(0.9, 1.1)");

        t.set_root_height(0.01);
        t.set_population_size(2.0 / (10.0 * 4.0));
        double l = t.compute_log_likelihood();
        REQUIRE(l == Approx(-31.77866581319647));
        REQUIRE(t.get_degree_of_root() == 2);
    }
}

TEST_CASE("Testing collection settings from minimal config", "[CollectionSettings]") {
    SECTION("Testing data/minimal-config.yml") {
        std::string cfg_path = "data/minimal-config.yml";
        CollectionSettings settings = CollectionSettings(cfg_path);
        std::string e =  "";
        e += "---\n";
        e += "event_model_prior:\n";
        e += "    uniform:\n";
        e += "event_time_prior:\n";
        e += "    exponential_distribution:\n";
        e += "        rate: 100\n";
        e += "mcmc_settings:\n";
        e += "    chain_length: 100000\n";
        e += "    sample_frequency: 100\n";
        e += "comparisons:\n";
        e += "- comparison:\n";
        e += "    path: data/hemi129.nex\n";
        e += "    ploidy: 2\n";
        e += "    genotypes_are_diploid: true\n";
        e += "    markers_are_dominant: false\n";
        e += "    population_name_delimiter: ' '\n";
        e += "    population_name_is_prefix: true\n";
        e += "    constant_sites_removed: true\n";
        e += "    use_empirical_starting_value_for_freq_1: false\n";
        e += "    equal_population_sizes: false\n";
        e += "    equal_state_frequencies: true\n";
        e += "    parameters:\n";
        e += "        population_size:\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                exponential_distribution:\n";
        e += "                    rate: 1000\n";
        e += "        mutation_rate:\n";
        e += "            value: 1\n";
        e += "            estimate: false\n";
        e += "        freq_1:\n";
        e += "            value: 0.5\n";
        e += "            estimate: false\n";
        e += "operator_settings:\n";
        e += "    auto_optimize: true\n";
        e += "    auto_optimize_delay: 10000\n";
        e += "    operators:\n";
        e += "        ModelOperator:\n";
        e += "            weight: 0\n";
        e += "            number_of_auxiliary_categories: 4\n";
        e += "        ConcentrationScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.5\n";
        e += "        CompositeHeightSizeRateMixer:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        CompositeHeightSizeRateScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.5\n";
        e += "        ComparisonHeightScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        ComparisonMutationRateScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.3\n";
        e += "        RootPopulationSizeScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        ChildPopulationSizeScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        FreqMover:\n";
        e += "            weight: 0\n";
        e += "            window: 0.1\n";

        REQUIRE(settings.to_string() == e);
        REQUIRE(settings.get_path() == "data/minimal-config.yml");
        REQUIRE(settings.using_dpp() == false);
        REQUIRE(settings.get_chain_length() == 100000);
        REQUIRE(settings.get_sample_frequency() == 100);
        REQUIRE(settings.get_number_of_comparisons() == 1);
        REQUIRE(settings.get_number_of_comparisons_with_free_mutation_rate() == 0);
        REQUIRE(settings.get_number_of_comparisons_with_free_state_frequencies() == 0);
        REQUIRE(settings.get_number_of_comparisons_with_free_population_size() == 1);
    }
}

TEST_CASE("Testing collection settings from minimal config with two comparisons", "[CollectionSettings]") {
    SECTION("Testing minimal two comparisons") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "comparisons:\n";
        cfg_stream << "- comparison:\n";
        cfg_stream << "    path: \"hemi129.nex\"\n";
        cfg_stream << "- comparison:\n";
        cfg_stream << "    path: \"diploid-dna.nex\"\n";

        CollectionSettings settings = CollectionSettings(cfg_stream, cfg_path);

        std::string e =  "";
        e += "---\n";
        e += "event_model_prior:\n";
        e += "    dirichlet_process:\n";
        e += "        parameters:\n";
        e += "            concentration:\n";
        e += "                estimate: true\n";
        e += "                prior:\n";
        e += "                    gamma_distribution:\n";
        e += "                        shape: 2\n";
        e += "                        scale: 0.5\n";
        e += "event_time_prior:\n";
        e += "    exponential_distribution:\n";
        e += "        rate: 100\n";
        e += "mcmc_settings:\n";
        e += "    chain_length: 100000\n";
        e += "    sample_frequency: 100\n";
        e += "comparisons:\n";
        e += "- comparison:\n";
        e += "    path: data/hemi129.nex\n";
        e += "    ploidy: 2\n";
        e += "    genotypes_are_diploid: true\n";
        e += "    markers_are_dominant: false\n";
        e += "    population_name_delimiter: ' '\n";
        e += "    population_name_is_prefix: true\n";
        e += "    constant_sites_removed: true\n";
        e += "    use_empirical_starting_value_for_freq_1: false\n";
        e += "    equal_population_sizes: false\n";
        e += "    equal_state_frequencies: true\n";
        e += "    parameters:\n";
        e += "        population_size:\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                exponential_distribution:\n";
        e += "                    rate: 1000\n";
        e += "        mutation_rate:\n";
        e += "            value: 1\n";
        e += "            estimate: false\n";
        e += "        freq_1:\n";
        e += "            value: 0.5\n";
        e += "            estimate: false\n";
        e += "- comparison:\n";
        e += "    path: data/diploid-dna.nex\n";
        e += "    ploidy: 2\n";
        e += "    genotypes_are_diploid: true\n";
        e += "    markers_are_dominant: false\n";
        e += "    population_name_delimiter: ' '\n";
        e += "    population_name_is_prefix: true\n";
        e += "    constant_sites_removed: true\n";
        e += "    use_empirical_starting_value_for_freq_1: false\n";
        e += "    equal_population_sizes: false\n";
        e += "    equal_state_frequencies: true\n";
        e += "    parameters:\n";
        e += "        population_size:\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                exponential_distribution:\n";
        e += "                    rate: 1000\n";
        e += "        mutation_rate:\n";
        e += "            value: 1\n";
        e += "            estimate: false\n";
        e += "        freq_1:\n";
        e += "            value: 0.5\n";
        e += "            estimate: false\n";
        e += "operator_settings:\n";
        e += "    auto_optimize: true\n";
        e += "    auto_optimize_delay: 10000\n";
        e += "    operators:\n";
        e += "        ModelOperator:\n";
        e += "            weight: 3\n";
        e += "            number_of_auxiliary_categories: 4\n";
        e += "        ConcentrationScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        CompositeHeightSizeRateMixer:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        CompositeHeightSizeRateScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.5\n";
        e += "        ComparisonHeightScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        ComparisonMutationRateScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.3\n";
        e += "        RootPopulationSizeScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        ChildPopulationSizeScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        FreqMover:\n";
        e += "            weight: 0\n";
        e += "            window: 0.1\n";

        REQUIRE(settings.to_string() == e);
        REQUIRE(settings.get_path() == "data/dummy.yml");
        REQUIRE(settings.using_dpp() == true);
        REQUIRE(settings.get_chain_length() == 100000);
        REQUIRE(settings.get_sample_frequency() == 100);
        REQUIRE(settings.get_number_of_comparisons() == 2);
        REQUIRE(settings.get_number_of_comparisons_with_free_mutation_rate() == 0);
        REQUIRE(settings.get_number_of_comparisons_with_free_state_frequencies() == 0);
        REQUIRE(settings.get_number_of_comparisons_with_free_population_size() == 2);
    }
}

TEST_CASE("Testing override model prior", "[CollectionSettings]") {
    SECTION("Testing model prior override") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "event_model_prior:\n";
        cfg_stream << "    uniform:\n";
        cfg_stream << "comparisons:\n";
        cfg_stream << "- comparison:\n";
        cfg_stream << "    path: \"hemi129.nex\"\n";
        cfg_stream << "- comparison:\n";
        cfg_stream << "    path: \"diploid-dna.nex\"\n";

        CollectionSettings settings = CollectionSettings(cfg_stream, cfg_path);

        std::string e =  "";
        e += "---\n";
        e += "event_model_prior:\n";
        e += "    uniform:\n";
        e += "event_time_prior:\n";
        e += "    exponential_distribution:\n";
        e += "        rate: 100\n";
        e += "mcmc_settings:\n";
        e += "    chain_length: 100000\n";
        e += "    sample_frequency: 100\n";
        e += "comparisons:\n";
        e += "- comparison:\n";
        e += "    path: data/hemi129.nex\n";
        e += "    ploidy: 2\n";
        e += "    genotypes_are_diploid: true\n";
        e += "    markers_are_dominant: false\n";
        e += "    population_name_delimiter: ' '\n";
        e += "    population_name_is_prefix: true\n";
        e += "    constant_sites_removed: true\n";
        e += "    use_empirical_starting_value_for_freq_1: false\n";
        e += "    equal_population_sizes: false\n";
        e += "    equal_state_frequencies: true\n";
        e += "    parameters:\n";
        e += "        population_size:\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                exponential_distribution:\n";
        e += "                    rate: 1000\n";
        e += "        mutation_rate:\n";
        e += "            value: 1\n";
        e += "            estimate: false\n";
        e += "        freq_1:\n";
        e += "            value: 0.5\n";
        e += "            estimate: false\n";
        e += "- comparison:\n";
        e += "    path: data/diploid-dna.nex\n";
        e += "    ploidy: 2\n";
        e += "    genotypes_are_diploid: true\n";
        e += "    markers_are_dominant: false\n";
        e += "    population_name_delimiter: ' '\n";
        e += "    population_name_is_prefix: true\n";
        e += "    constant_sites_removed: true\n";
        e += "    use_empirical_starting_value_for_freq_1: false\n";
        e += "    equal_population_sizes: false\n";
        e += "    equal_state_frequencies: true\n";
        e += "    parameters:\n";
        e += "        population_size:\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                exponential_distribution:\n";
        e += "                    rate: 1000\n";
        e += "        mutation_rate:\n";
        e += "            value: 1\n";
        e += "            estimate: false\n";
        e += "        freq_1:\n";
        e += "            value: 0.5\n";
        e += "            estimate: false\n";
        e += "operator_settings:\n";
        e += "    auto_optimize: true\n";
        e += "    auto_optimize_delay: 10000\n";
        e += "    operators:\n";
        e += "        ModelOperator:\n";
        e += "            weight: 3\n";
        e += "            number_of_auxiliary_categories: 4\n";
        e += "        ConcentrationScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.5\n";
        e += "        CompositeHeightSizeRateMixer:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        CompositeHeightSizeRateScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.5\n";
        e += "        ComparisonHeightScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        ComparisonMutationRateScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.3\n";
        e += "        RootPopulationSizeScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        ChildPopulationSizeScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        FreqMover:\n";
        e += "            weight: 0\n";
        e += "            window: 0.1\n";

        REQUIRE(settings.to_string() == e);
        REQUIRE(settings.get_path() == "data/dummy.yml");
        REQUIRE(settings.using_dpp() == false);
        REQUIRE(settings.get_chain_length() == 100000);
        REQUIRE(settings.get_sample_frequency() == 100);
        REQUIRE(settings.get_number_of_comparisons() == 2);
        REQUIRE(settings.get_number_of_comparisons_with_free_mutation_rate() == 0);
        REQUIRE(settings.get_number_of_comparisons_with_free_state_frequencies() == 0);
        REQUIRE(settings.get_number_of_comparisons_with_free_population_size() == 2);
    }
}

TEST_CASE("Testing override model prior with DPP", "[CollectionSettings]") {
    SECTION("Testing DPP override") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "event_model_prior:\n";
        cfg_stream << "    dirichlet_process:\n";
        cfg_stream << "        parameters:\n";
        cfg_stream << "            concentration:\n";
        cfg_stream << "                value: 10.0\n";
        cfg_stream << "                estimate: true\n";
        cfg_stream << "                prior:\n";
        cfg_stream << "                    gamma_distribution:\n";
        cfg_stream << "                        shape: 100.0\n";
        cfg_stream << "                        prior_mean_number_of_events: 1.9\n";
        cfg_stream << "comparisons:\n";
        cfg_stream << "- comparison:\n";
        cfg_stream << "    path: \"hemi129.nex\"\n";
        cfg_stream << "- comparison:\n";
        cfg_stream << "    path: \"diploid-dna.nex\"\n";

        CollectionSettings settings = CollectionSettings(cfg_stream, cfg_path);

        std::string e =  "";
        e += "---\n";
        e += "event_model_prior:\n";
        e += "    dirichlet_process:\n";
        e += "        parameters:\n";
        e += "            concentration:\n";
        e += "                value: 10\n";
        e += "                estimate: true\n";
        e += "                prior:\n";
        e += "                    gamma_distribution:\n";
        e += "                        shape: 100\n";
        e += "                        scale: 0.09\n";
        e += "event_time_prior:\n";
        e += "    exponential_distribution:\n";
        e += "        rate: 100\n";
        e += "mcmc_settings:\n";
        e += "    chain_length: 100000\n";
        e += "    sample_frequency: 100\n";
        e += "comparisons:\n";
        e += "- comparison:\n";
        e += "    path: data/hemi129.nex\n";
        e += "    ploidy: 2\n";
        e += "    genotypes_are_diploid: true\n";
        e += "    markers_are_dominant: false\n";
        e += "    population_name_delimiter: ' '\n";
        e += "    population_name_is_prefix: true\n";
        e += "    constant_sites_removed: true\n";
        e += "    use_empirical_starting_value_for_freq_1: false\n";
        e += "    equal_population_sizes: false\n";
        e += "    equal_state_frequencies: true\n";
        e += "    parameters:\n";
        e += "        population_size:\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                exponential_distribution:\n";
        e += "                    rate: 1000\n";
        e += "        mutation_rate:\n";
        e += "            value: 1\n";
        e += "            estimate: false\n";
        e += "        freq_1:\n";
        e += "            value: 0.5\n";
        e += "            estimate: false\n";
        e += "- comparison:\n";
        e += "    path: data/diploid-dna.nex\n";
        e += "    ploidy: 2\n";
        e += "    genotypes_are_diploid: true\n";
        e += "    markers_are_dominant: false\n";
        e += "    population_name_delimiter: ' '\n";
        e += "    population_name_is_prefix: true\n";
        e += "    constant_sites_removed: true\n";
        e += "    use_empirical_starting_value_for_freq_1: false\n";
        e += "    equal_population_sizes: false\n";
        e += "    equal_state_frequencies: true\n";
        e += "    parameters:\n";
        e += "        population_size:\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                exponential_distribution:\n";
        e += "                    rate: 1000\n";
        e += "        mutation_rate:\n";
        e += "            value: 1\n";
        e += "            estimate: false\n";
        e += "        freq_1:\n";
        e += "            value: 0.5\n";
        e += "            estimate: false\n";
        e += "operator_settings:\n";
        e += "    auto_optimize: true\n";
        e += "    auto_optimize_delay: 10000\n";
        e += "    operators:\n";
        e += "        ModelOperator:\n";
        e += "            weight: 3\n";
        e += "            number_of_auxiliary_categories: 4\n";
        e += "        ConcentrationScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        CompositeHeightSizeRateMixer:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        CompositeHeightSizeRateScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.5\n";
        e += "        ComparisonHeightScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        ComparisonMutationRateScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.3\n";
        e += "        RootPopulationSizeScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        ChildPopulationSizeScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        FreqMover:\n";
        e += "            weight: 0\n";
        e += "            window: 0.1\n";

        REQUIRE(settings.to_string() == e);
        REQUIRE(settings.get_path() == "data/dummy.yml");
        REQUIRE(settings.using_dpp() == true);
        REQUIRE(settings.get_chain_length() == 100000);
        REQUIRE(settings.get_sample_frequency() == 100);
        REQUIRE(settings.get_number_of_comparisons() == 2);
        REQUIRE(settings.get_number_of_comparisons_with_free_mutation_rate() == 0);
        REQUIRE(settings.get_number_of_comparisons_with_free_state_frequencies() == 0);
        REQUIRE(settings.get_number_of_comparisons_with_free_population_size() == 2);
    }
}

TEST_CASE("Testing no comparisons", "[CollectionSettings]") {
    SECTION("Testing no comparisons") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "event_model_prior:\n";
        cfg_stream << "    dirichlet_process:\n";
        cfg_stream << "        parameters:\n";
        cfg_stream << "            concentration:\n";
        cfg_stream << "                value: 10.0\n";
        cfg_stream << "                estimate: true\n";
        cfg_stream << "                prior:\n";
        cfg_stream << "                    gamma_distribution:\n";
        cfg_stream << "                        shape: 100.0\n";
        cfg_stream << "                        prior_mean_number_of_events: 1.9\n";

        REQUIRE_THROWS_AS(CollectionSettings(cfg_stream, cfg_path), EcoevolityYamlConfigError);
    }
}

TEST_CASE("Testing empty comparisons", "[CollectionSettings]") {
    SECTION("Testing empty comparisons") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "event_model_prior:\n";
        cfg_stream << "    dirichlet_process:\n";
        cfg_stream << "        parameters:\n";
        cfg_stream << "            concentration:\n";
        cfg_stream << "                value: 10.0\n";
        cfg_stream << "                estimate: true\n";
        cfg_stream << "                prior:\n";
        cfg_stream << "                    gamma_distribution:\n";
        cfg_stream << "                        shape: 100.0\n";
        cfg_stream << "                        prior_mean_number_of_events: 1.9\n";
        cfg_stream << "comparisons:\n";

        REQUIRE_THROWS_AS(CollectionSettings(cfg_stream, cfg_path), EcoevolityYamlConfigError);
    }
}

TEST_CASE("Testing empty config", "[CollectionSettings]") {
    SECTION("Testing empty config") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "";

        REQUIRE_THROWS_AS(CollectionSettings(cfg_stream, cfg_path), EcoevolityYamlConfigError);
    }
}

TEST_CASE("Testing bad YAML formatting", "[CollectionSettings]") {
    SECTION("Testing bad YAML") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "comparisons: :\n";
        cfg_stream << "- comparison:\n";
        cfg_stream << "    path: \"hemi129.nex\"\n";
        cfg_stream << "- comparison:\n";
        cfg_stream << "    path: \"diploid-dna.nex\"\n";

        REQUIRE_THROWS_AS(CollectionSettings(cfg_stream, cfg_path), std::exception);
    }
}

TEST_CASE("Testing fixing with no value", "[CollectionSettings]") {
    SECTION("Testing valueless fix") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "comparisons:\n";
        cfg_stream << "- comparison:\n";
        cfg_stream << "    path: \"hemi129.nex\"\n";
        cfg_stream << "- comparison:\n";
        cfg_stream << "    path: \"diploid-dna.nex\"\n";
        cfg_stream << "    parameters:\n";
        cfg_stream << "        population_size:\n";
        cfg_stream << "            estimate: false\n";

        REQUIRE_THROWS_AS(CollectionSettings(cfg_stream, cfg_path), EcoevolityPositiveRealParameterSettingError);
    }
}

TEST_CASE("Testing override with global settings", "[CollectionSettings]") {
    SECTION("Testing global settings") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "global_comparison_settings:\n";
        cfg_stream << "    genotypes_are_diploid: false\n";
        cfg_stream << "    markers_are_dominant: true\n";
        cfg_stream << "    population_name_delimiter: '-'\n";
        cfg_stream << "    population_name_is_prefix: false\n";
        cfg_stream << "    constant_sites_removed: false\n";
        cfg_stream << "    use_empirical_starting_value_for_freq_1: true\n";
        cfg_stream << "    equal_population_sizes: true\n";
        cfg_stream << "    equal_state_frequencies: false\n";
        cfg_stream << "comparisons:\n";
        cfg_stream << "- comparison:\n";
        cfg_stream << "    path: haploid-standard.nex\n";
        cfg_stream << "- comparison:\n";
        cfg_stream << "    path: haploid-standard-missing.nex\n";

        CollectionSettings settings = CollectionSettings(cfg_stream, cfg_path);

        std::string e =  "";
        e += "---\n";
        e += "event_model_prior:\n";
        e += "    dirichlet_process:\n";
        e += "        parameters:\n";
        e += "            concentration:\n";
        e += "                estimate: true\n";
        e += "                prior:\n";
        e += "                    gamma_distribution:\n";
        e += "                        shape: 2\n";
        e += "                        scale: 0.5\n";
        e += "event_time_prior:\n";
        e += "    exponential_distribution:\n";
        e += "        rate: 100\n";
        e += "mcmc_settings:\n";
        e += "    chain_length: 100000\n";
        e += "    sample_frequency: 100\n";
        e += "comparisons:\n";
        e += "- comparison:\n";
        e += "    path: data/haploid-standard.nex\n";
        e += "    ploidy: 2\n";
        e += "    genotypes_are_diploid: false\n";
        e += "    markers_are_dominant: true\n";
        e += "    population_name_delimiter: '-'\n";
        e += "    population_name_is_prefix: false\n";
        e += "    constant_sites_removed: false\n";
        e += "    use_empirical_starting_value_for_freq_1: true\n";
        e += "    equal_population_sizes: true\n";
        e += "    equal_state_frequencies: false\n";
        e += "    parameters:\n";
        e += "        population_size:\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                exponential_distribution:\n";
        e += "                    rate: 1000\n";
        e += "        mutation_rate:\n";
        e += "            value: 1\n";
        e += "            estimate: false\n";
        e += "        freq_1:\n";
        e += "            value: 0.521739\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                beta_distribution:\n";
        e += "                    alpha: 1\n";
        e += "                    beta: 1\n";
        e += "- comparison:\n";
        e += "    path: data/haploid-standard-missing.nex\n";
        e += "    ploidy: 2\n";
        e += "    genotypes_are_diploid: false\n";
        e += "    markers_are_dominant: true\n";
        e += "    population_name_delimiter: '-'\n";
        e += "    population_name_is_prefix: false\n";
        e += "    constant_sites_removed: false\n";
        e += "    use_empirical_starting_value_for_freq_1: true\n";
        e += "    equal_population_sizes: true\n";
        e += "    equal_state_frequencies: false\n";
        e += "    parameters:\n";
        e += "        population_size:\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                exponential_distribution:\n";
        e += "                    rate: 1000\n";
        e += "        mutation_rate:\n";
        e += "            value: 1\n";
        e += "            estimate: false\n";
        e += "        freq_1:\n";
        e += "            value: 0.611111\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                beta_distribution:\n";
        e += "                    alpha: 1\n";
        e += "                    beta: 1\n";
        e += "operator_settings:\n";
        e += "    auto_optimize: true\n";
        e += "    auto_optimize_delay: 10000\n";
        e += "    operators:\n";
        e += "        ModelOperator:\n";
        e += "            weight: 3\n";
        e += "            number_of_auxiliary_categories: 4\n";
        e += "        ConcentrationScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        CompositeHeightSizeRateMixer:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        CompositeHeightSizeRateScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.5\n";
        e += "        ComparisonHeightScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        ComparisonMutationRateScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.3\n";
        e += "        RootPopulationSizeScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        ChildPopulationSizeScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        FreqMover:\n";
        e += "            weight: 1\n";
        e += "            window: 0.1\n";

        REQUIRE(settings.to_string() == e);
        REQUIRE(settings.get_path() == "data/dummy.yml");
        REQUIRE(settings.using_dpp() == true);
        REQUIRE(settings.get_chain_length() == 100000);
        REQUIRE(settings.get_sample_frequency() == 100);
        REQUIRE(settings.get_number_of_comparisons() == 2);
        REQUIRE(settings.get_number_of_comparisons_with_free_mutation_rate() == 0);
        REQUIRE(settings.get_number_of_comparisons_with_free_state_frequencies() == 2);
        REQUIRE(settings.get_number_of_comparisons_with_free_population_size() == 2);
    }
}

TEST_CASE("Testing override with global settings with parameters", "[CollectionSettings]") {
    SECTION("Testing global settings") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "global_comparison_settings:\n";
        cfg_stream << "    genotypes_are_diploid: false\n";
        cfg_stream << "    markers_are_dominant: true\n";
        cfg_stream << "    population_name_delimiter: '-'\n";
        cfg_stream << "    population_name_is_prefix: false\n";
        cfg_stream << "    constant_sites_removed: false\n";
        cfg_stream << "    use_empirical_starting_value_for_freq_1: true\n";
        cfg_stream << "    equal_population_sizes: true\n";
        cfg_stream << "    equal_state_frequencies: false\n";
        cfg_stream << "    parameters:\n";
        cfg_stream << "        mutation_rate:\n";
        cfg_stream << "            estimate: true\n";
        cfg_stream << "        freq_1:\n";
        cfg_stream << "            value: 0.56\n";
        cfg_stream << "            estimate: true\n";
        cfg_stream << "            prior:\n";
        cfg_stream << "                beta_distribution:\n";
        cfg_stream << "                    alpha: 3\n";
        cfg_stream << "                    beta: 2\n";
        cfg_stream << "        population_size:\n";
        cfg_stream << "            value: 0.01\n";
        cfg_stream << "            estimate: false\n";
        cfg_stream << "            prior:\n";
        cfg_stream << "                gamma_distribution:\n";
        cfg_stream << "                    shape: 2.0\n";
        cfg_stream << "                    scale: 0.001\n";
        cfg_stream << "comparisons:\n";
        cfg_stream << "- comparison:\n";
        cfg_stream << "    path: haploid-standard.nex\n";
        cfg_stream << "- comparison:\n";
        cfg_stream << "    path: haploid-standard-missing.nex\n";

        CollectionSettings settings = CollectionSettings(cfg_stream, cfg_path);

        std::string e =  "";
        e += "---\n";
        e += "event_model_prior:\n";
        e += "    dirichlet_process:\n";
        e += "        parameters:\n";
        e += "            concentration:\n";
        e += "                estimate: true\n";
        e += "                prior:\n";
        e += "                    gamma_distribution:\n";
        e += "                        shape: 2\n";
        e += "                        scale: 0.5\n";
        e += "event_time_prior:\n";
        e += "    exponential_distribution:\n";
        e += "        rate: 100\n";
        e += "mcmc_settings:\n";
        e += "    chain_length: 100000\n";
        e += "    sample_frequency: 100\n";
        e += "comparisons:\n";
        e += "- comparison:\n";
        e += "    path: data/haploid-standard.nex\n";
        e += "    ploidy: 2\n";
        e += "    genotypes_are_diploid: false\n";
        e += "    markers_are_dominant: true\n";
        e += "    population_name_delimiter: '-'\n";
        e += "    population_name_is_prefix: false\n";
        e += "    constant_sites_removed: false\n";
        e += "    use_empirical_starting_value_for_freq_1: true\n";
        e += "    equal_population_sizes: true\n";
        e += "    equal_state_frequencies: false\n";
        e += "    parameters:\n";
        e += "        population_size:\n";
        e += "            value: 0.01\n";
        e += "            estimate: false\n";
        e += "        mutation_rate:\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                gamma_distribution:\n";
        e += "                    shape: 1000\n";
        e += "                    scale: 0.001\n";
        e += "        freq_1:\n";
        e += "            value: 0.521739\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                beta_distribution:\n";
        e += "                    alpha: 3\n";
        e += "                    beta: 2\n";
        e += "- comparison:\n";
        e += "    path: data/haploid-standard-missing.nex\n";
        e += "    ploidy: 2\n";
        e += "    genotypes_are_diploid: false\n";
        e += "    markers_are_dominant: true\n";
        e += "    population_name_delimiter: '-'\n";
        e += "    population_name_is_prefix: false\n";
        e += "    constant_sites_removed: false\n";
        e += "    use_empirical_starting_value_for_freq_1: true\n";
        e += "    equal_population_sizes: true\n";
        e += "    equal_state_frequencies: false\n";
        e += "    parameters:\n";
        e += "        population_size:\n";
        e += "            value: 0.01\n";
        e += "            estimate: false\n";
        e += "        mutation_rate:\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                gamma_distribution:\n";
        e += "                    shape: 1000\n";
        e += "                    scale: 0.001\n";
        e += "        freq_1:\n";
        e += "            value: 0.611111\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                beta_distribution:\n";
        e += "                    alpha: 3\n";
        e += "                    beta: 2\n";
        e += "operator_settings:\n";
        e += "    auto_optimize: true\n";
        e += "    auto_optimize_delay: 10000\n";
        e += "    operators:\n";
        e += "        ModelOperator:\n";
        e += "            weight: 3\n";
        e += "            number_of_auxiliary_categories: 4\n";
        e += "        ConcentrationScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        CompositeHeightSizeRateMixer:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        CompositeHeightSizeRateScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.5\n";
        e += "        ComparisonHeightScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        ComparisonMutationRateScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.3\n";
        e += "        RootPopulationSizeScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.5\n";
        e += "        ChildPopulationSizeScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.5\n";
        e += "        FreqMover:\n";
        e += "            weight: 1\n";
        e += "            window: 0.1\n";

        REQUIRE(settings.to_string() == e);
        REQUIRE(settings.get_path() == "data/dummy.yml");
        REQUIRE(settings.using_dpp() == true);
        REQUIRE(settings.get_chain_length() == 100000);
        REQUIRE(settings.get_sample_frequency() == 100);
        REQUIRE(settings.get_number_of_comparisons() == 2);
        REQUIRE(settings.get_number_of_comparisons_with_free_mutation_rate() == 2);
        REQUIRE(settings.get_number_of_comparisons_with_free_state_frequencies() == 2);
        REQUIRE(settings.get_number_of_comparisons_with_free_population_size() == 0);
    }
}

TEST_CASE("Testing collection settings from full config", "[CollectionSettings]") {
    SECTION("Testing data/config.yml") {
        std::string cfg_path = "data/config.yml";
        CollectionSettings settings = CollectionSettings(cfg_path);
        std::string e =  "";
        e += "---\n";
        e += "event_model_prior:\n";
        e += "    dirichlet_process:\n";
        e += "        parameters:\n";
        e += "            concentration:\n";
        e += "                value: 5\n";
        e += "                estimate: true\n";
        e += "                prior:\n";
        e += "                    gamma_distribution:\n";
        e += "                        shape: 10\n";
        e += "                        scale: 0.141422\n";
        e += "event_time_prior:\n";
        e += "    gamma_distribution:\n";
        e += "        shape: 2\n";
        e += "        scale: 0.001\n";
        e += "        offset: 0\n";
        e += "mcmc_settings:\n";
        e += "    chain_length: 2000000\n";
        e += "    sample_frequency: 2000\n";
        e += "comparisons:\n";
        e += "- comparison:\n";
        e += "    path: data/data3.nex\n";
        e += "    ploidy: 2\n";
        e += "    genotypes_are_diploid: true\n";
        e += "    markers_are_dominant: false\n";
        e += "    population_name_delimiter: ' '\n";
        e += "    population_name_is_prefix: true\n";
        e += "    constant_sites_removed: true\n";
        e += "    use_empirical_starting_value_for_freq_1: false\n";
        e += "    equal_population_sizes: false\n";
        e += "    equal_state_frequencies: false\n";
        e += "    parameters:\n";
        e += "        population_size:\n";
        e += "            value: 0.005\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                gamma_distribution:\n";
        e += "                    shape: 10\n";
        e += "                    scale: 0.0001\n";
        e += "                    offset: 0\n";
        e += "        mutation_rate:\n";
        e += "            value: 1\n";
        e += "            estimate: false\n";
        e += "        freq_1:\n";
        e += "            value: 0.5\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                beta_distribution:\n";
        e += "                    alpha: 2.5\n";
        e += "                    beta: 5.3\n";
        e += "- comparison:\n";
        e += "    path: data/data1.nex\n";
        e += "    ploidy: 2\n";
        e += "    genotypes_are_diploid: true\n";
        e += "    markers_are_dominant: false\n";
        e += "    population_name_delimiter: ' '\n";
        e += "    population_name_is_prefix: true\n";
        e += "    constant_sites_removed: true\n";
        e += "    use_empirical_starting_value_for_freq_1: false\n";
        e += "    equal_population_sizes: true\n";
        e += "    equal_state_frequencies: false\n";
        e += "    parameters:\n";
        e += "        population_size:\n";
        e += "            value: 0.01\n";
        e += "            estimate: false\n";
        e += "        mutation_rate:\n";
        e += "            value: 1\n";
        e += "            estimate: false\n";
        e += "        freq_1:\n";
        e += "            value: 0.5\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                beta_distribution:\n";
        e += "                    alpha: 2.5\n";
        e += "                    beta: 5.3\n";
        e += "- comparison:\n";
        e += "    path: data/data2.nex\n";
        e += "    ploidy: 2\n";
        e += "    genotypes_are_diploid: true\n";
        e += "    markers_are_dominant: false\n";
        e += "    population_name_delimiter: '-'\n";
        e += "    population_name_is_prefix: false\n";
        e += "    constant_sites_removed: false\n";
        e += "    use_empirical_starting_value_for_freq_1: false\n";
        e += "    equal_population_sizes: false\n";
        e += "    equal_state_frequencies: true\n";
        e += "    parameters:\n";
        e += "        population_size:\n";
        e += "            value: 0.005\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                gamma_distribution:\n";
        e += "                    shape: 10\n";
        e += "                    scale: 0.0001\n";
        e += "                    offset: 0\n";
        e += "        mutation_rate:\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                gamma_distribution:\n";
        e += "                    shape: 100\n";
        e += "                    scale: 0.01\n";
        e += "        freq_1:\n";
        e += "            value: 0.5\n";
        e += "            estimate: false\n";
        e += "operator_settings:\n";
        e += "    auto_optimize: true\n";
        e += "    auto_optimize_delay: 20000\n";
        e += "    operators:\n";
        e += "        ModelOperator:\n";
        e += "            weight: 5\n";
        e += "            number_of_auxiliary_categories: 5\n";
        e += "        ConcentrationScaler:\n";
        e += "            weight: 3\n";
        e += "            scale: 0.2\n";
        e += "        CompositeHeightSizeRateMixer:\n";
        e += "            weight: 2\n";
        e += "            scale: 0.4\n";
        e += "        CompositeHeightSizeRateScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.5\n";
        e += "        ComparisonHeightScaler:\n";
        e += "            weight: 3\n";
        e += "            scale: 0.3\n";
        e += "        ComparisonMutationRateScaler:\n";
        e += "            weight: 2\n";
        e += "            scale: 0.5\n";
        e += "        RootPopulationSizeScaler:\n";
        e += "            weight: 2\n";
        e += "            scale: 0.2\n";
        e += "        ChildPopulationSizeScaler:\n";
        e += "            weight: 2\n";
        e += "            scale: 0.2\n";
        e += "        FreqMover:\n";
        e += "            weight: 0.5\n";
        e += "            window: 0.2\n";

        REQUIRE(settings.to_string() == e);
        REQUIRE(settings.get_path() == "data/config.yml");
        REQUIRE(settings.using_dpp() == true);
        REQUIRE(settings.get_chain_length() == 2000000);
        REQUIRE(settings.get_sample_frequency() == 2000);
        REQUIRE(settings.get_number_of_comparisons() == 3);
        REQUIRE(settings.get_number_of_comparisons_with_free_mutation_rate() == 1);
        REQUIRE(settings.get_number_of_comparisons_with_free_state_frequencies() == 2);
        REQUIRE(settings.get_number_of_comparisons_with_free_population_size() == 2);
    }
}

TEST_CASE("Testing path comparison between collection settings",
        "[CollectionSettings]") {
    SECTION("Testing path comparisons between collection settings") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream1;
        cfg_stream1 << "comparisons:\n";
        cfg_stream1 << "- comparison:\n";
        cfg_stream1 << "    path: \"hemi129.nex\"\n";
        cfg_stream1 << "- comparison:\n";
        cfg_stream1 << "    path: \"hemi129-altname1.nex\"\n";
        cfg_stream1 << "- comparison:\n";
        cfg_stream1 << "    path: \"hemi129-altname2.nex\"\n";

        CollectionSettings settings1 = CollectionSettings(cfg_stream1, cfg_path);

        std::stringstream cfg_stream2;
        cfg_stream2 << "comparisons:\n";
        cfg_stream2 << "- comparison:\n";
        cfg_stream2 << "    path: \"hemi129-altname2.nex\"\n";
        cfg_stream2 << "- comparison:\n";
        cfg_stream2 << "    path: \"hemi129.nex\"\n";
        cfg_stream2 << "- comparison:\n";
        cfg_stream2 << "    path: \"hemi129-altname1.nex\"\n";

        CollectionSettings settings2 = CollectionSettings(cfg_stream2, cfg_path);

        std::stringstream cfg_stream3;
        cfg_stream3 << "comparisons:\n";
        cfg_stream3 << "- comparison:\n";
        cfg_stream3 << "    path: \"hemi129.nex\"\n";
        cfg_stream3 << "- comparison:\n";
        cfg_stream3 << "    path: \"hemi129-altname1.nex\"\n";
        cfg_stream3 << "- comparison:\n";
        cfg_stream3 << "    path: \"hemi129-altname2.nex\"\n";
        cfg_stream3 << "- comparison:\n";
        cfg_stream3 << "    path: \"hemi129-altname3.nex\"\n";

        CollectionSettings settings3 = CollectionSettings(cfg_stream3, cfg_path);

        std::stringstream cfg_stream4;
        cfg_stream4 << "comparisons:\n";
        cfg_stream4 << "- comparison:\n";
        cfg_stream4 << "    path: \"hemi129.nex\"\n";
        cfg_stream4 << "- comparison:\n";
        cfg_stream4 << "    path: \"hemi129-altname1.nex\"\n";

        CollectionSettings settings4 = CollectionSettings(cfg_stream4, cfg_path);

        REQUIRE(settings1.same_comparison_paths(settings2) == true);
        REQUIRE(settings2.same_comparison_paths(settings1) == true);
        REQUIRE(settings1.same_comparison_paths(settings3) == false);
        REQUIRE(settings1.same_comparison_paths(settings4) == false);
        REQUIRE(settings2.same_comparison_paths(settings3) == false);
        REQUIRE(settings2.same_comparison_paths(settings4) == false);
        REQUIRE(settings3.same_comparison_paths(settings1) == false);
        REQUIRE(settings3.same_comparison_paths(settings2) == false);
        REQUIRE(settings3.same_comparison_paths(settings4) == false);
        REQUIRE(settings4.same_comparison_paths(settings1) == false);
        REQUIRE(settings4.same_comparison_paths(settings2) == false);
        REQUIRE(settings4.same_comparison_paths(settings3) == false);
    }
}

TEST_CASE("Testing fixed model prior", "[CollectionSettings]") {
    SECTION("Testing fixed model") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "event_model_prior:\n";
        cfg_stream << "    fixed: [0, 1]\n";
        cfg_stream << "comparisons:\n";
        cfg_stream << "- comparison:\n";
        cfg_stream << "    path: \"hemi129.nex\"\n";
        cfg_stream << "- comparison:\n";
        cfg_stream << "    path: \"diploid-dna.nex\"\n";

        CollectionSettings settings = CollectionSettings(cfg_stream, cfg_path);

        std::string e =  "";
        e += "---\n";
        e += "event_model_prior:\n";
        e += "    fixed: [0, 1]\n";
        e += "event_time_prior:\n";
        e += "    exponential_distribution:\n";
        e += "        rate: 100\n";
        e += "mcmc_settings:\n";
        e += "    chain_length: 100000\n";
        e += "    sample_frequency: 100\n";
        e += "comparisons:\n";
        e += "- comparison:\n";
        e += "    path: data/hemi129.nex\n";
        e += "    ploidy: 2\n";
        e += "    genotypes_are_diploid: true\n";
        e += "    markers_are_dominant: false\n";
        e += "    population_name_delimiter: ' '\n";
        e += "    population_name_is_prefix: true\n";
        e += "    constant_sites_removed: true\n";
        e += "    use_empirical_starting_value_for_freq_1: false\n";
        e += "    equal_population_sizes: false\n";
        e += "    equal_state_frequencies: true\n";
        e += "    parameters:\n";
        e += "        population_size:\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                exponential_distribution:\n";
        e += "                    rate: 1000\n";
        e += "        mutation_rate:\n";
        e += "            value: 1\n";
        e += "            estimate: false\n";
        e += "        freq_1:\n";
        e += "            value: 0.5\n";
        e += "            estimate: false\n";
        e += "- comparison:\n";
        e += "    path: data/diploid-dna.nex\n";
        e += "    ploidy: 2\n";
        e += "    genotypes_are_diploid: true\n";
        e += "    markers_are_dominant: false\n";
        e += "    population_name_delimiter: ' '\n";
        e += "    population_name_is_prefix: true\n";
        e += "    constant_sites_removed: true\n";
        e += "    use_empirical_starting_value_for_freq_1: false\n";
        e += "    equal_population_sizes: false\n";
        e += "    equal_state_frequencies: true\n";
        e += "    parameters:\n";
        e += "        population_size:\n";
        e += "            estimate: true\n";
        e += "            prior:\n";
        e += "                exponential_distribution:\n";
        e += "                    rate: 1000\n";
        e += "        mutation_rate:\n";
        e += "            value: 1\n";
        e += "            estimate: false\n";
        e += "        freq_1:\n";
        e += "            value: 0.5\n";
        e += "            estimate: false\n";
        e += "operator_settings:\n";
        e += "    auto_optimize: true\n";
        e += "    auto_optimize_delay: 10000\n";
        e += "    operators:\n";
        e += "        ModelOperator:\n";
        e += "            weight: 0\n";
        e += "            number_of_auxiliary_categories: 4\n";
        e += "        ConcentrationScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.5\n";
        e += "        CompositeHeightSizeRateMixer:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        CompositeHeightSizeRateScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.5\n";
        e += "        ComparisonHeightScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        ComparisonMutationRateScaler:\n";
        e += "            weight: 0\n";
        e += "            scale: 0.3\n";
        e += "        RootPopulationSizeScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        ChildPopulationSizeScaler:\n";
        e += "            weight: 1\n";
        e += "            scale: 0.5\n";
        e += "        FreqMover:\n";
        e += "            weight: 0\n";
        e += "            window: 0.1\n";

        REQUIRE(settings.to_string() == e);
        REQUIRE(settings.get_path() == "data/dummy.yml");
        REQUIRE(settings.using_dpp() == false);
        REQUIRE(settings.get_chain_length() == 100000);
        REQUIRE(settings.get_sample_frequency() == 100);
        REQUIRE(settings.get_number_of_comparisons() == 2);
        REQUIRE(settings.get_number_of_comparisons_with_free_mutation_rate() == 0);
        REQUIRE(settings.get_number_of_comparisons_with_free_state_frequencies() == 0);
        REQUIRE(settings.get_number_of_comparisons_with_free_population_size() == 2);
        REQUIRE(settings.event_model_is_fixed() == true);
    }
}

