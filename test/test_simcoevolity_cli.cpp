#include "catch.hpp"
#include "ecoevolity/simcoevolity.hpp"

#include "ecoevolity/rng.hpp"
#include "ecoevolity/path.hpp"

RandomNumberGenerator _SIMCOEVOLITY_CLI_RNG = RandomNumberGenerator();

TEST_CASE("Testing simcoevolity constant sites error", "[SimcoevolityCLI]") {

    SECTION("Testing constant sites error") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_alpha = 2.0;
        double f1_beta = 1.0;
        double f2_alpha = 1.0;
        double f2_beta = 0.5;
        double f3_alpha = 1.0;
        double f3_beta = 2.0;
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _SIMCOEVOLITY_CLI_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-t10.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-t10-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 10\n";
        os << "    sample_frequency: 1\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    operators:\n";
        os << "        EventTimeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        LeafPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    equal_population_sizes: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129-with-constant.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_alpha << "\n";
        os << "                    beta: " << f1_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_alpha << "\n";
        os << "                    beta: " << f2_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_alpha<< "\n";
        os << "                    beta: " << f3_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "simcoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "-n";
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "test1-";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            &arg6[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        REQUIRE_THROWS_AS(
                (simcoevolity_main<CollectionSettings, ComparisonPopulationTreeCollection>(argc, argv)),
                EcoevolityConstantSitesError);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing simcoevolity constant sites error for dirichlet trees", "[SimcoevolityCLI]") {

    SECTION("Testing constant sites error for dirichlet trees") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_alpha = 2.0;
        double f1_beta = 1.0;
        double f2_alpha = 1.0;
        double f2_beta = 0.5;
        double f3_alpha = 1.0;
        double f3_beta = 2.0;
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _SIMCOEVOLITY_CLI_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-t11.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-t11-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 10\n";
        os << "    sample_frequency: 1\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    operators:\n";
        os << "        EventTimeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MeanPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RelativePopulationSizeMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    parameters:\n";
        os << "        population_size_multipliers:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                dirichlet_distribution:\n";
        os << "                    alpha: [10.0, 10.0, 10.0]\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129-with-constant.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_alpha << "\n";
        os << "                    beta: " << f1_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_alpha << "\n";
        os << "                    beta: " << f2_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_alpha<< "\n";
        os << "                    beta: " << f3_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "simcoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "-n";
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "test2-";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            &arg6[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        REQUIRE_THROWS_AS(
                (simcoevolity_main<DirichletCollectionSettings, ComparisonDirichletPopulationTreeCollection>(argc, argv)),
                EcoevolityConstantSitesError);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing simcoevolity relaxed constant sites setting", "[SimcoevolityCLI]") {

    SECTION("Testing constant sites relaxed") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_alpha = 2.0;
        double f1_beta = 1.0;
        double f2_alpha = 1.0;
        double f2_beta = 0.5;
        double f3_alpha = 1.0;
        double f3_beta = 2.0;
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _SIMCOEVOLITY_CLI_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-t12.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-t12-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 10\n";
        os << "    sample_frequency: 1\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    operators:\n";
        os << "        EventTimeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        LeafPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    equal_population_sizes: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129-with-constant.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_alpha << "\n";
        os << "                    beta: " << f1_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_alpha << "\n";
        os << "                    beta: " << f2_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_alpha << "\n";
        os << "                    beta: " << f3_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "simcoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "-n";
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "test3-";
        char arg7[] = "--relax-constant-sites";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            &arg6[0],
            &arg7[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        ret = simcoevolity_main<CollectionSettings, ComparisonPopulationTreeCollection>(argc, argv);
        REQUIRE(ret == 0);

        REQUIRE(path::exists("data/test3-simcoevolity-model-used-for-sims.yml"));

        for (unsigned int i = 0; i < 10; ++i) {
            std::string sim_rep = string_util::pad_int(i, 2);
            std::string sim_prefix = "data/test3-simcoevolity-sim-" + sim_rep + "-";
            std::string expected_true_path = sim_prefix + "true-values.txt";
            std::string expected_config_path = sim_prefix + "config.yml";
            std::string expected_align_path1 = sim_prefix + "hemi129-with-constant.nex";
            std::string expected_align_path2 = sim_prefix + "hemi129-altname1.nex";
            std::string expected_align_path3 = sim_prefix + "hemi129-altname2.nex";
            REQUIRE(path::exists(expected_true_path));
            REQUIRE(path::exists(expected_config_path));
            REQUIRE(path::exists(expected_align_path1));
            REQUIRE(path::exists(expected_align_path2));
            REQUIRE(path::exists(expected_align_path3));
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing simcoevolity relaxed constant sites setting for dirichlet trees",
        "[SimcoevolityCLI]") {

    SECTION("Testing constant sites relaxed for dirichlet trees") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_alpha = 2.0;
        double f1_beta = 1.0;
        double f2_alpha = 1.0;
        double f2_beta = 0.5;
        double f3_alpha = 1.0;
        double f3_beta = 2.0;
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _SIMCOEVOLITY_CLI_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-t13.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-t13-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 10\n";
        os << "    sample_frequency: 1\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 3\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    operators:\n";
        os << "        EventTimeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MeanPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RelativePopulationSizeMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    parameters:\n";
        os << "        population_size_multipliers:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                dirichlet_distribution:\n";
        os << "                    alpha: [10.0, 10.0, 10.0]\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129-with-constant.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_alpha << "\n";
        os << "                    beta: " << f1_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_alpha << "\n";
        os << "                    beta: " << f2_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_alpha << "\n";
        os << "                    beta: " << f3_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "simcoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "-n";
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "test4-";
        char arg7[] = "--relax-constant-sites";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            &arg6[0],
            &arg7[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        ret = simcoevolity_main<DirichletCollectionSettings, ComparisonDirichletPopulationTreeCollection>(argc, argv);
        REQUIRE(ret == 0);

        REQUIRE(path::exists("data/test4-simcoevolity-model-used-for-sims.yml"));

        for (unsigned int i = 0; i < 10; ++i) {
            std::string sim_rep = string_util::pad_int(i, 2);
            std::string sim_prefix = "data/test4-simcoevolity-sim-" + sim_rep + "-";
            std::string expected_true_path = sim_prefix + "true-values.txt";
            std::string expected_config_path = sim_prefix + "config.yml";
            std::string expected_align_path1 = sim_prefix + "hemi129-with-constant.nex";
            std::string expected_align_path2 = sim_prefix + "hemi129-altname1.nex";
            std::string expected_align_path3 = sim_prefix + "hemi129-altname2.nex";
            REQUIRE(path::exists(expected_true_path));
            REQUIRE(path::exists(expected_config_path));
            REQUIRE(path::exists(expected_align_path1));
            REQUIRE(path::exists(expected_align_path2));
            REQUIRE(path::exists(expected_align_path3));
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing simcoevolity missing sites error", "[SimcoevolityCLI]") {

    SECTION("Testing missing sites error") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_alpha = 2.0;
        double f1_beta = 1.0;
        double f2_alpha = 1.0;
        double f2_beta = 0.5;
        double f3_alpha = 1.0;
        double f3_beta = 2.0;
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _SIMCOEVOLITY_CLI_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-t14.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-t14-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 10\n";
        os << "    sample_frequency: 1\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    operators:\n";
        os << "        EventTimeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        LeafPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    equal_population_sizes: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129-with-missing.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_alpha << "\n";
        os << "                    beta: " << f1_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_alpha << "\n";
        os << "                    beta: " << f2_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_alpha << "\n";
        os << "                    beta: " << f3_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "-n";
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "test5-";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            &arg6[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        REQUIRE_THROWS_AS(
                (simcoevolity_main<CollectionSettings, ComparisonPopulationTreeCollection>(argc, argv)),
                EcoevolityMissingDataError);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing simcoevolity missing sites error for dirichlet trees", "[SimcoevolityCLI]") {

    SECTION("Testing missing sites error for dirichlet trees") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_alpha = 2.0;
        double f1_beta = 1.0;
        double f2_alpha = 1.0;
        double f2_beta = 0.5;
        double f3_alpha = 1.0;
        double f3_beta = 2.0;
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _SIMCOEVOLITY_CLI_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-t15.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-t15-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 10\n";
        os << "    sample_frequency: 1\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    operators:\n";
        os << "        EventTimeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MeanPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RelativePopulationSizeMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    parameters:\n";
        os << "        population_size_multipliers:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                dirichlet_distribution:\n";
        os << "                    alpha: [10.0, 10.0, 10.0]\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129-with-missing.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_alpha << "\n";
        os << "                    beta: " << f1_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_alpha << "\n";
        os << "                    beta: " << f2_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_alpha << "\n";
        os << "                    beta: " << f3_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "-n";
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "test6-";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            &arg6[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        REQUIRE_THROWS_AS(
                (simcoevolity_main<DirichletCollectionSettings, ComparisonDirichletPopulationTreeCollection>(argc, argv)),
                EcoevolityMissingDataError);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing simcoevolity relaxed missing sites setting", "[SimcoevolityCLI]") {

    SECTION("Testing missing sites relaxed") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_alpha = 2.0;
        double f1_beta = 1.0;
        double f2_alpha = 1.0;
        double f2_beta = 0.5;
        double f3_alpha = 1.0;
        double f3_beta = 2.0;
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _SIMCOEVOLITY_CLI_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-t16.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-t16-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 10\n";
        os << "    sample_frequency: 1\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    operators:\n";
        os << "        EventTimeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        LeafPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    equal_population_sizes: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129-with-missing.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_alpha << "\n";
        os << "                    beta: " << f1_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_alpha << "\n";
        os << "                    beta: " << f2_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_alpha << "\n";
        os << "                    beta: " << f3_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "-n";
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "test7-";
        char arg7[] = "--relax-missing-sites";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            &arg6[0],
            &arg7[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        ret = simcoevolity_main<CollectionSettings, ComparisonPopulationTreeCollection>(argc, argv);
        REQUIRE(ret == 0);

        REQUIRE(path::exists("data/test7-simcoevolity-model-used-for-sims.yml"));

        for (unsigned int i = 0; i < 10; ++i) {
            std::string sim_rep = string_util::pad_int(i, 2);
            std::string sim_prefix = "data/test7-simcoevolity-sim-" + sim_rep + "-";
            std::string expected_true_path = sim_prefix + "true-values.txt";
            std::string expected_config_path = sim_prefix + "config.yml";
            std::string expected_align_path1 = sim_prefix + "hemi129-with-missing.nex";
            std::string expected_align_path2 = sim_prefix + "hemi129-altname1.nex";
            std::string expected_align_path3 = sim_prefix + "hemi129-altname2.nex";
            REQUIRE(path::exists(expected_true_path));
            REQUIRE(path::exists(expected_config_path));
            REQUIRE(path::exists(expected_align_path1));
            REQUIRE(path::exists(expected_align_path2));
            REQUIRE(path::exists(expected_align_path3));
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing simcoevolity relaxed missing sites setting for dirichlet trees",
        "[SimcoevolityCLI]") {

    SECTION("Testing missing sites relaxed for dirichlet trees") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_alpha = 2.0;
        double f1_beta = 1.0;
        double f2_alpha = 1.0;
        double f2_beta = 0.5;
        double f3_alpha = 1.0;
        double f3_beta = 2.0;
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _SIMCOEVOLITY_CLI_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-t17.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-t17-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 10\n";
        os << "    sample_frequency: 1\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    operators:\n";
        os << "        EventTimeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MeanPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RelativePopulationSizeMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    parameters:\n";
        os << "        population_size_multipliers:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                dirichlet_distribution:\n";
        os << "                    alpha: [10.0, 10.0, 10.0]\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129-with-missing.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_alpha << "\n";
        os << "                    beta: " << f1_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_alpha << "\n";
        os << "                    beta: " << f2_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_alpha << "\n";
        os << "                    beta: " << f3_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "-n";
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "test8-";
        char arg7[] = "--relax-missing-sites";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            &arg6[0],
            &arg7[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        ret = simcoevolity_main<DirichletCollectionSettings, ComparisonDirichletPopulationTreeCollection>(argc, argv);
        REQUIRE(ret == 0);

        REQUIRE(path::exists("data/test8-simcoevolity-model-used-for-sims.yml"));

        for (unsigned int i = 0; i < 10; ++i) {
            std::string sim_rep = string_util::pad_int(i, 2);
            std::string sim_prefix = "data/test8-simcoevolity-sim-" + sim_rep + "-";
            std::string expected_true_path = sim_prefix + "true-values.txt";
            std::string expected_config_path = sim_prefix + "config.yml";
            std::string expected_align_path1 = sim_prefix + "hemi129-with-missing.nex";
            std::string expected_align_path2 = sim_prefix + "hemi129-altname1.nex";
            std::string expected_align_path3 = sim_prefix + "hemi129-altname2.nex";
            REQUIRE(path::exists(expected_true_path));
            REQUIRE(path::exists(expected_config_path));
            REQUIRE(path::exists(expected_align_path1));
            REQUIRE(path::exists(expected_align_path2));
            REQUIRE(path::exists(expected_align_path3));
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing simcoevolity constrained singleton error",
        "[SimcoevolityCLI]") {

    SECTION("Testing error for singleton with constrained pop size") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_alpha = 2.0;
        double f1_beta = 1.0;
        double f2_alpha = 1.0;
        double f2_beta = 0.5;
        double f3_alpha = 1.0;
        double f3_beta = 2.0;
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _SIMCOEVOLITY_CLI_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-t18.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-t18-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 2000000\n";
        os << "    sample_frequency: 100\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    operators:\n";
        os << "        EventTimeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        LeafPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    equal_population_sizes: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    equal_population_sizes: true\n";
        os << "    path: hemi129-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_alpha << "\n";
        os << "                    beta: " << f1_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_alpha << "\n";
        os << "                    beta: " << f2_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_alpha << "\n";
        os << "                    beta: " << f3_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "simcoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "-n";
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "test9-";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            &arg6[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        REQUIRE_THROWS_AS(
                (simcoevolity_main<CollectionSettings, ComparisonPopulationTreeCollection>(argc, argv)),
                EcoevolityComparisonSettingError);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing simcoevolity fixed singleton error",
        "[SimcoevolityCLI]") {

    SECTION("Testing error for singleton with fixed pop size") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_alpha = 2.0;
        double f1_beta = 1.0;
        double f2_alpha = 1.0;
        double f2_beta = 0.5;
        double f3_alpha = 1.0;
        double f3_beta = 2.0;
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _SIMCOEVOLITY_CLI_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-t19.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-t19-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 2000000\n";
        os << "    sample_frequency: 100\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    operators:\n";
        os << "        EventTimeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        LeafPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    equal_population_sizes: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            value: 0.001\n";
        os << "            estimate: false\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_alpha << "\n";
        os << "                    beta: " << f1_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_alpha << "\n";
        os << "                    beta: " << f2_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_alpha << "\n";
        os << "                    beta: " << f3_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "simcoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "-n";
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "test10-";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            &arg6[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        REQUIRE_THROWS_AS(
                (simcoevolity_main<CollectionSettings, ComparisonPopulationTreeCollection>(argc, argv)),
                EcoevolityComparisonSettingError);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing simcoevolity fixed singleton error for dirichlet tree",
        "[SimcoevolityCLI]") {

    SECTION("Testing error for dirichlet singleton with fixed pop size") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_alpha = 2.0;
        double f1_beta = 1.0;
        double f2_alpha = 1.0;
        double f2_beta = 0.5;
        double f3_alpha = 1.0;
        double f3_beta = 2.0;
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _SIMCOEVOLITY_CLI_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-t20.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-t20-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 2000000\n";
        os << "    sample_frequency: 100\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    operators:\n";
        os << "        EventTimeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MeanPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RelativePopulationSizeMixer:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        population_size_multipliers:\n";
        os << "            value: [1.0, 1.0]\n";
        os << "            estimate: false\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_alpha << "\n";
        os << "                    beta: " << f1_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1-singleton.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        population_size_multipliers:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                dirichlet_distribution:\n";
        os << "                    alpha: [10.0, 10.0]\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_alpha << "\n";
        os << "                    beta: " << f2_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        population_size_multipliers:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                dirichlet_distribution:\n";
        os << "                    alpha: [10.0, 10.0, 10.0]\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_alpha << "\n";
        os << "                    beta: " << f3_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "simcoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "-n";
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "test11-";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            &arg6[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;
        REQUIRE_THROWS_AS(
                (simcoevolity_main<DirichletCollectionSettings, ComparisonDirichletPopulationTreeCollection>(argc, argv)),
                EcoevolityComparisonSettingError);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing simcoevolity population label conflict", "[SimcoevolityCLI]") {

    SECTION("Testing duplicate pop labels") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_alpha = 2.0;
        double f1_beta = 1.0;
        double f2_alpha = 1.0;
        double f2_beta = 0.5;
        double f3_alpha = 1.0;
        double f3_beta = 2.0;
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _SIMCOEVOLITY_CLI_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-t21.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-t21-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 10\n";
        os << "    sample_frequency: 1\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    operators:\n";
        os << "        EventTimeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        LeafPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    equal_population_sizes: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_alpha << "\n";
        os << "                    beta: " << f1_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-name-conflict.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_alpha << "\n";
        os << "                    beta: " << f2_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_alpha << "\n";
        os << "                    beta: " << f3_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "simcoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "-n";
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "test12-";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            &arg6[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        REQUIRE_THROWS_AS(
                (simcoevolity_main<CollectionSettings, ComparisonPopulationTreeCollection>(argc, argv)),
                EcoevolityCollectionSettingError);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing simcoevolity parameters only setting", "[SimcoevolityCLI]") {

    SECTION("Testing parameters only.") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_alpha = 2.0;
        double f1_beta = 1.0;
        double f2_alpha = 1.0;
        double f2_beta = 0.5;
        double f3_alpha = 1.0;
        double f3_beta = 2.0;
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _SIMCOEVOLITY_CLI_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-t22.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-t22-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 10\n";
        os << "    sample_frequency: 1\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    operators:\n";
        os << "        EventTimeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        LeafPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: false\n";
        os << "    equal_population_sizes: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129-with-constant.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_alpha << "\n";
        os << "                    beta: " << f1_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_alpha << "\n";
        os << "                    beta: " << f2_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_alpha << "\n";
        os << "                    beta: " << f3_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "simcoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "-n";
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "test-params-only-";
        char arg7[] = "--parameters-only";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            &arg6[0],
            &arg7[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        ret = simcoevolity_main<CollectionSettings, ComparisonPopulationTreeCollection>(argc, argv);
        REQUIRE(ret == 0);

        REQUIRE(! path::exists("data/test-params-only-simcoevolity-model-used-for-sims.yml"));
        REQUIRE(! path::exists("data/test-params-only-simcoevolity-parameter-values.txt"));

        for (unsigned int i = 0; i < 10; ++i) {
            std::string sim_rep = string_util::pad_int(i, 2);
            std::string sim_prefix = "data/test-params-only-simcoevolity-sim-" + sim_rep + "-";
            std::string expected_true_path = sim_prefix + "true-values.txt";
            std::string expected_config_path = sim_prefix + "config.yml";
            std::string expected_align_path1 = sim_prefix + "hemi129-with-constant.nex";
            std::string expected_align_path2 = sim_prefix + "hemi129-altname1.nex";
            std::string expected_align_path3 = sim_prefix + "hemi129-altname2.nex";
            REQUIRE(! path::exists(expected_true_path));
            REQUIRE(! path::exists(expected_config_path));
            REQUIRE(! path::exists(expected_align_path1));
            REQUIRE(! path::exists(expected_align_path2));
            REQUIRE(! path::exists(expected_align_path3));
        }

        delete[] cfg_path;
    }
}


TEST_CASE("Testing simcoevolity triallelic sites error", "[SimcoevolityCLI]") {

    SECTION("Testing triallelic sites error") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_alpha = 2.0;
        double f1_beta = 1.0;
        double f2_alpha = 1.0;
        double f2_beta = 0.5;
        double f3_alpha = 1.0;
        double f3_beta = 2.0;
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _SIMCOEVOLITY_CLI_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-t23.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-t23-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 10\n";
        os << "    sample_frequency: 1\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    operators:\n";
        os << "        EventTimeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        LeafPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    equal_population_sizes: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: diploid-dna-triallelic.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_alpha << "\n";
        os << "                    beta: " << f1_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_alpha << "\n";
        os << "                    beta: " << f2_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_alpha << "\n";
        os << "                    beta: " << f3_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "-n";
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "test13-";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            &arg6[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        REQUIRE_THROWS_AS(
                (simcoevolity_main<CollectionSettings, ComparisonPopulationTreeCollection>(argc, argv)),
                EcoevolityTriallelicDataError);

        delete[] cfg_path;
    }
}


TEST_CASE("Testing simcoevolity relax triallelic sites option", "[SimcoevolityCLI]") {

    SECTION("Testing relax triallelic sites option") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_alpha = 2.0;
        double f1_beta = 1.0;
        double f2_alpha = 1.0;
        double f2_beta = 0.5;
        double f3_alpha = 1.0;
        double f3_beta = 2.0;
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _SIMCOEVOLITY_CLI_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-t24.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-t24-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 10\n";
        os << "    sample_frequency: 1\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    operators:\n";
        os << "        EventTimeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        LeafPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    equal_population_sizes: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: diploid-dna-triallelic.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_alpha << "\n";
        os << "                    beta: " << f1_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_alpha << "\n";
        os << "                    beta: " << f2_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_alpha << "\n";
        os << "                    beta: " << f3_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "ecoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "-n";
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "test14-";
        char arg7[] = "--relax-triallelic-sites";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            &arg6[0],
            &arg7[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        ret = simcoevolity_main<CollectionSettings, ComparisonPopulationTreeCollection>(argc, argv);
        REQUIRE(ret == 0);

        REQUIRE(path::exists("data/test14-simcoevolity-model-used-for-sims.yml"));

        for (unsigned int i = 0; i < 10; ++i) {
            std::string sim_rep = string_util::pad_int(i, 2);
            std::string sim_prefix = "data/test14-simcoevolity-sim-" + sim_rep + "-";
            std::string expected_true_path = sim_prefix + "true-values.txt";
            std::string expected_config_path = sim_prefix + "config.yml";
            std::string expected_align_path1 = sim_prefix + "diploid-dna-triallelic.nex";
            std::string expected_align_path2 = sim_prefix + "hemi129-altname1.nex";
            std::string expected_align_path3 = sim_prefix + "hemi129-altname2.nex";
            REQUIRE(path::exists(expected_true_path));
            REQUIRE(path::exists(expected_config_path));
            REQUIRE(path::exists(expected_align_path1));
            REQUIRE(path::exists(expected_align_path2));
            REQUIRE(path::exists(expected_align_path3));
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing simcoevolity missing charsets error", "[SimcoevolityCLI]") {

    SECTION("Testing missing charsets error") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_alpha = 2.0;
        double f1_beta = 1.0;
        double f2_alpha = 1.0;
        double f2_beta = 0.5;
        double f3_alpha = 1.0;
        double f3_beta = 2.0;
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _SIMCOEVOLITY_CLI_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-t10.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-t10-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 10\n";
        os << "    sample_frequency: 1\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    operators:\n";
        os << "        EventTimeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        LeafPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: true\n";
        os << "    equal_population_sizes: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_alpha << "\n";
        os << "                    beta: " << f1_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_alpha << "\n";
        os << "                    beta: " << f2_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_alpha<< "\n";
        os << "                    beta: " << f3_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "simcoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "-n";
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "test15-";
        char arg7[] = "--charsets";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            &arg6[0],
            &arg7[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        REQUIRE_THROWS_AS(
                (simcoevolity_main<CollectionSettings, ComparisonPopulationTreeCollection>(argc, argv)),
                EcoevolityError);

        delete[] cfg_path;
    }
}

TEST_CASE("Testing simcoevolity charsets setting", "[SimcoevolityCLI]") {

    SECTION("Testing charsets") {
        double height_shape = 10.0;
        double height_scale = 0.001;
        double size1_shape = 10.0;
        double size1_scale = 0.0001;
        double size2_shape = 2.0;
        double size2_scale = 0.001;
        double size3_shape = 5.0;
        double size3_scale = 0.0005;
        double f1_alpha = 2.0;
        double f1_beta = 1.0;
        double f2_alpha = 1.0;
        double f2_beta = 0.5;
        double f3_alpha = 1.0;
        double f3_beta = 2.0;
        double mult2_shape = 100.0;
        double mult2_scale = 0.005;
        double mult3_shape = 100.0;
        double mult3_scale = 0.02;
        double concentration_shape = 5.0;
        double concentration_scale = 0.2;
        std::string auto_optimize = "true";
        std::string tag = _SIMCOEVOLITY_CLI_RNG.random_string(10);
        std::string test_path = "data/tmp-config-" + tag + "-t12.cfg";
        std::string log_path = "data/tmp-config-" + tag + "-t12-state-run-1.log";
        std::ofstream os;
        os.open(test_path);
        os << "event_time_prior:\n";
        os << "    gamma_distribution:\n";
        os << "        shape: " << height_shape << "\n";
        os << "        scale: " << height_scale << "\n";
        os << "event_model_prior:\n";
        os << "    dirichlet_process:\n";
        os << "        parameters:\n";
        os << "            concentration:\n";
        os << "                estimate: true\n";
        os << "                prior:\n";
        os << "                    gamma_distribution:\n";
        os << "                        shape: " << concentration_shape << "\n";
        os << "                        scale: " << concentration_scale << "\n";
        os << "mcmc_settings:\n";
        os << "    chain_length: 10\n";
        os << "    sample_frequency: 1\n";
        os << "operator_settings:\n";
        os << "    auto_optimize: " << auto_optimize << "\n";
        os << "    auto_optimize_delay: 10000\n";
        os << "    operators:\n";
        os << "        ModelOperator:\n";
        os << "            number_of_auxiliary_categories: 5\n";
        os << "            weight: 1.0\n";
        os << "        ConcentrationScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "global_comparison_settings:\n";
        os << "    operators:\n";
        os << "        EventTimeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        MutationRateScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        RootPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        LeafPopulationSizeScaler:\n";
        os << "            scale: 0.5\n";
        os << "            weight: 1.0\n";
        os << "        FreqMover:\n";
        os << "            window: 0.1\n";
        os << "            weight: 1.0\n";
        os << "    genotypes_are_diploid: true\n";
        os << "    markers_are_dominant: false\n";
        os << "    population_name_delimiter: \" \"\n";
        os << "    population_name_is_prefix: true\n";
        os << "    constant_sites_removed: false\n";
        os << "    equal_population_sizes: false\n";
        os << "comparisons:\n";
        os << "- comparison:\n";
        os << "    path: hemi129-with-constant.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size1_shape << "\n";
        os << "                    scale: " << size1_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f1_alpha << "\n";
        os << "                    beta: " << f1_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            value: 1.0\n";
        os << "            estimate: false\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname1.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size2_shape << "\n";
        os << "                    scale: " << size2_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f2_alpha << "\n";
        os << "                    beta: " << f2_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult2_shape << "\n";
        os << "                    scale: " << mult2_scale << "\n";
        os << "- comparison:\n";
        os << "    path: hemi129-altname2.nex\n";
        os << "    parameters:\n";
        os << "        population_size:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << size3_shape << "\n";
        os << "                    scale: " << size3_scale << "\n";
        os << "        freq_1:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                beta_distribution:\n";
        os << "                    alpha: " << f3_alpha << "\n";
        os << "                    beta: " << f3_beta << "\n";
        os << "        mutation_rate:\n";
        os << "            estimate: true\n";
        os << "            prior:\n";
        os << "                gamma_distribution:\n";
        os << "                    shape: " << mult3_shape << "\n";
        os << "                    scale: " << mult3_scale << "\n";
        os.close();
        REQUIRE(path::exists(test_path));

        char arg0[] = "simcoevolity";
        char arg1[] = "--seed";
        char arg2[] = "283402";
        char arg3[] = "-n";
        char arg4[] = "10";
        char arg5[] = "--prefix";
        char arg6[] = "test16-";
        char arg7[] = "--relax-constant-sites";
        char arg8[] = "--charsets";
        char * cfg_path = new char[test_path.size() + 1];
        std::copy(test_path.begin(), test_path.end(), cfg_path);
        cfg_path[test_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            &arg2[0],
            &arg3[0],
            &arg4[0],
            &arg5[0],
            &arg6[0],
            &arg7[0],
            &arg8[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        ret = simcoevolity_main<CollectionSettings, ComparisonPopulationTreeCollection>(argc, argv);
        REQUIRE(ret == 0);

        REQUIRE(path::exists("data/test16-simcoevolity-model-used-for-sims.yml"));

        for (unsigned int i = 0; i < 10; ++i) {
            std::string sim_rep = string_util::pad_int(i, 2);
            std::string sim_prefix = "data/test16-simcoevolity-sim-" + sim_rep + "-";
            std::string expected_true_path = sim_prefix + "true-values.txt";
            std::string expected_config_path = sim_prefix + "config.yml";
            std::string expected_align_path1 = sim_prefix + "hemi129-with-constant.nex";
            std::string expected_align_path2 = sim_prefix + "hemi129-altname1.nex";
            std::string expected_align_path3 = sim_prefix + "hemi129-altname2.nex";
            REQUIRE(path::exists(expected_true_path));
            REQUIRE(path::exists(expected_config_path));
            REQUIRE(path::exists(expected_align_path1));
            REQUIRE(path::exists(expected_align_path2));
            REQUIRE(path::exists(expected_align_path3));
        }

        delete[] cfg_path;
    }
}
