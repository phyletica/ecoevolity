#include "catch.hpp"
#include "ecoevolity/general_tree_operator_schedule.hpp"
#include "ecoevolity/general_tree_settings.hpp"
#include "ecoevolity/settings_io.hpp"
#include "ecoevolity/parameter.hpp"
#include "ecoevolity/tree.hpp"


TEST_CASE("Testing generalized BasePopulationTree config",
        "[GeneralTreeOperatorSchedule]") {
    SECTION("Testing with uniform_root_and_betas tree prior") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "---\n";
        cfg_stream << "tree_model:\n";
        cfg_stream << "    tree_space: \"generalized\"\n";
        cfg_stream << "    tree_prior:\n";
        cfg_stream << "        uniform_root_and_betas:\n";
        cfg_stream << "            parameters:\n";
        cfg_stream << "                root_height:\n";
        cfg_stream << "                    value: 0.3\n";
        cfg_stream << "                    estimate: true\n";
        cfg_stream << "                    prior:\n";
        cfg_stream << "                        gamma_distribution:\n";
        cfg_stream << "                            shape: 8.0\n";
        cfg_stream << "                            mean: 0.3\n";
        cfg_stream << "                alpha_of_node_height_beta_prior:\n";
        cfg_stream << "                    value: 1.0\n";
        cfg_stream << "                    estimate: true\n";
        cfg_stream << "                    prior:\n";
        cfg_stream << "                        gamma_distribution:\n";
        cfg_stream << "                            shape: 4.0\n";
        cfg_stream << "                            mean: 1.0\n";
        cfg_stream << "data:\n";
        cfg_stream << "    ploidy: 1\n";
        cfg_stream << "    constant_sites_removed: false\n";
        cfg_stream << "    yaml_allele_counts:\n";
        cfg_stream << "        path: \"diploid-dna-constant-missing.yml\"\n";

        PopulationTreeSettings settings = PopulationTreeSettings(cfg_stream, cfg_path);

        GeneralTreeOperatorSchedule<BasePopulationTree> op_schedule(
                settings.operator_settings,
                5); // number of leaves

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<BasePopulationTree> > > ops;

        ops = op_schedule.get_split_lump_rj_operators();
        REQUIRE(! ops.empty());

        for (auto oper : ops) {
            REQUIRE(oper->helper_ops.size() > 0);
        }

        ops.clear();
        ops = op_schedule.get_operators(BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::topology);
        REQUIRE(ops.size() > 0);
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler operator",
        "[GeneralTreeOperatorSchedule]") {
    SECTION("Testing SplitLumpNodesRevJumpSampler") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "---\n";
        cfg_stream << "tree_model:\n";
        cfg_stream << "    tree_space: \"generalized\"\n";
        cfg_stream << "    tree_prior:\n";
        cfg_stream << "        uniform_root_and_betas:\n";
        cfg_stream << "            parameters:\n";
        cfg_stream << "                root_height:\n";
        cfg_stream << "                    value: 0.3\n";
        cfg_stream << "                    estimate: true\n";
        cfg_stream << "                    prior:\n";
        cfg_stream << "                        gamma_distribution:\n";
        cfg_stream << "                            shape: 8.0\n";
        cfg_stream << "                            mean: 0.3\n";
        cfg_stream << "                alpha_of_node_height_beta_prior:\n";
        cfg_stream << "                    value: 1.0\n";
        cfg_stream << "                    estimate: true\n";
        cfg_stream << "                    prior:\n";
        cfg_stream << "                        gamma_distribution:\n";
        cfg_stream << "                            shape: 4.0\n";
        cfg_stream << "                            mean: 1.0\n";
        cfg_stream << "data:\n";
        cfg_stream << "    ploidy: 1\n";
        cfg_stream << "    constant_sites_removed: false\n";
        cfg_stream << "    yaml_allele_counts:\n";
        cfg_stream << "        path: \"diploid-dna-constant-missing.yml\"\n";
        cfg_stream << "mcmc_settings:\n";
        cfg_stream << "    chain_length: 15000\n";
        cfg_stream << "    sample_frequency: 10\n";
        cfg_stream << "    operators:\n";
        cfg_stream << "        SplitLumpNodesRevJumpSampler:\n";
        cfg_stream << "            weight: [2, 3, 5]\n";
        cfg_stream << "            tuning_parameter: [2.0, 1.0, 10.0]\n";
        cfg_stream << "            auto_optimize: [true, false, false]\n";
        cfg_stream << "            auto_optimize_delay: [100, 50, 20]\n";

        PopulationTreeSettings settings = PopulationTreeSettings(cfg_stream, cfg_path);

        auto rj_op = settings.operator_settings->tunable_operators["SplitLumpNodesRevJumpSampler"];
        for (unsigned int i = 0; i < rj_op.get_number_of_operators(); ++i) {
            std::cout << "weight " << i << ": " << rj_op.get_weight(i) << "\n";
            std::cout << "tuning param " << i << ": " << rj_op.get_tuning_parameter(i) << "\n";
            std::cout << "auto opt " << i << ": " << rj_op.auto_optimizing(i) << "\n";
            std::cout << "auto opt delay " << i << ": " << rj_op.get_auto_optimize_delay(i) << "\n";
        }

        GeneralTreeOperatorSchedule<BasePopulationTree> op_schedule(
                settings.operator_settings,
                5); // number of leaves

        op_schedule.write_op_settings(std::cout);

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<BasePopulationTree> > > ops;

        ops = op_schedule.get_split_lump_rj_operators();
        REQUIRE(ops.size() == 3);

        std::set<double> expected_weights;
        expected_weights.insert(2);
        expected_weights.insert(3);
        expected_weights.insert(5);


        std::set<double> weights;
        for (auto oper : ops) {
            REQUIRE(oper->helper_ops.size() > 0);
            weights.insert(oper->get_weight());
            if (oper->get_weight() == 2.0) {
                REQUIRE(oper->get_coercable_parameter_value() == 2.0);
                REQUIRE(oper->auto_optimizing() == true);
                REQUIRE(oper->get_auto_optimize_delay() == 100);
            }
            else if (oper->get_weight() == 3.0) {
                REQUIRE(oper->get_coercable_parameter_value() == 1.0);
                REQUIRE(oper->auto_optimizing() == false);
                REQUIRE(oper->get_auto_optimize_delay() == 50);
            }
            else if (oper->get_weight() == 5.0) {
                REQUIRE(oper->get_coercable_parameter_value() == 10.0);
                REQUIRE(oper->auto_optimizing() == false);
                REQUIRE(oper->get_auto_optimize_delay() == 20);
            }
            else {
                REQUIRE(0 == 1);
            }
        }
        REQUIRE(weights == expected_weights);
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler operator with some defaults",
        "[GeneralTreeOperatorSchedule]") {
    SECTION("Testing SplitLumpNodesRevJumpSampler") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "---\n";
        cfg_stream << "tree_model:\n";
        cfg_stream << "    tree_space: \"generalized\"\n";
        cfg_stream << "    tree_prior:\n";
        cfg_stream << "        uniform_root_and_betas:\n";
        cfg_stream << "            parameters:\n";
        cfg_stream << "                root_height:\n";
        cfg_stream << "                    value: 0.3\n";
        cfg_stream << "                    estimate: true\n";
        cfg_stream << "                    prior:\n";
        cfg_stream << "                        gamma_distribution:\n";
        cfg_stream << "                            shape: 8.0\n";
        cfg_stream << "                            mean: 0.3\n";
        cfg_stream << "                alpha_of_node_height_beta_prior:\n";
        cfg_stream << "                    value: 1.0\n";
        cfg_stream << "                    estimate: true\n";
        cfg_stream << "                    prior:\n";
        cfg_stream << "                        gamma_distribution:\n";
        cfg_stream << "                            shape: 4.0\n";
        cfg_stream << "                            mean: 1.0\n";
        cfg_stream << "data:\n";
        cfg_stream << "    ploidy: 1\n";
        cfg_stream << "    constant_sites_removed: false\n";
        cfg_stream << "    yaml_allele_counts:\n";
        cfg_stream << "        path: \"diploid-dna-constant-missing.yml\"\n";
        cfg_stream << "mcmc_settings:\n";
        cfg_stream << "    chain_length: 15000\n";
        cfg_stream << "    sample_frequency: 10\n";
        cfg_stream << "    operators:\n";
        cfg_stream << "        SplitLumpNodesRevJumpSampler:\n";
        cfg_stream << "            tuning_parameter: [2.0, 5.0, 10.0]\n";

        PopulationTreeSettings settings = PopulationTreeSettings(cfg_stream, cfg_path);

        auto rj_op = settings.operator_settings->tunable_operators["SplitLumpNodesRevJumpSampler"];
        for (unsigned int i = 0; i < rj_op.get_number_of_operators(); ++i) {
            std::cout << "weight " << i << ": " << rj_op.get_weight(i) << "\n";
            std::cout << "tuning param " << i << ": " << rj_op.get_tuning_parameter(i) << "\n";
            std::cout << "auto opt " << i << ": " << rj_op.auto_optimizing(i) << "\n";
            std::cout << "auto opt delay " << i << ": " << rj_op.get_auto_optimize_delay(i) << "\n";
        }

        GeneralTreeOperatorSchedule<BasePopulationTree> op_schedule(
                settings.operator_settings,
                9); // number of leaves

        op_schedule.write_op_settings(std::cout);

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<BasePopulationTree> > > ops;

        ops = op_schedule.get_split_lump_rj_operators();
        REQUIRE(ops.size() == 3);

        std::set<double> expected_tuning_parameters;
        expected_tuning_parameters.insert(2);
        expected_tuning_parameters.insert(5);
        expected_tuning_parameters.insert(10);


        std::set<double> tuning_parameters;
        for (auto oper : ops) {
            REQUIRE(oper->helper_ops.size() > 0);
            tuning_parameters.insert(oper->get_coercable_parameter_value());
            if (oper->get_coercable_parameter_value() == 2.0) {
                REQUIRE(oper->get_weight() == 3.0);
                REQUIRE(oper->auto_optimizing() == true);
                REQUIRE(oper->get_auto_optimize_delay() == 50);
            }
            else if (oper->get_coercable_parameter_value() == 5.0) {
                REQUIRE(oper->get_weight() == 3.0);
                REQUIRE(oper->auto_optimizing() == true);
                REQUIRE(oper->get_auto_optimize_delay() == 50);
            }
            else if (oper->get_coercable_parameter_value() == 10.0) {
                REQUIRE(oper->get_weight() == 3.0);
                REQUIRE(oper->auto_optimizing() == true);
                REQUIRE(oper->get_auto_optimize_delay() == 50);
            }
            else {
                REQUIRE(0 == 1);
            }
        }
        REQUIRE(tuning_parameters == expected_tuning_parameters);
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler operator with zero weight",
        "[GeneralTreeOperatorSchedule]") {
    SECTION("Testing SplitLumpNodesRevJumpSampler with zero weight") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "---\n";
        cfg_stream << "tree_model:\n";
        cfg_stream << "    tree_space: \"generalized\"\n";
        cfg_stream << "    tree_prior:\n";
        cfg_stream << "        uniform_root_and_betas:\n";
        cfg_stream << "            parameters:\n";
        cfg_stream << "                root_height:\n";
        cfg_stream << "                    value: 0.3\n";
        cfg_stream << "                    estimate: true\n";
        cfg_stream << "                    prior:\n";
        cfg_stream << "                        gamma_distribution:\n";
        cfg_stream << "                            shape: 8.0\n";
        cfg_stream << "                            mean: 0.3\n";
        cfg_stream << "                alpha_of_node_height_beta_prior:\n";
        cfg_stream << "                    value: 1.0\n";
        cfg_stream << "                    estimate: true\n";
        cfg_stream << "                    prior:\n";
        cfg_stream << "                        gamma_distribution:\n";
        cfg_stream << "                            shape: 4.0\n";
        cfg_stream << "                            mean: 1.0\n";
        cfg_stream << "data:\n";
        cfg_stream << "    ploidy: 1\n";
        cfg_stream << "    constant_sites_removed: false\n";
        cfg_stream << "    yaml_allele_counts:\n";
        cfg_stream << "        path: \"diploid-dna-constant-missing.yml\"\n";
        cfg_stream << "mcmc_settings:\n";
        cfg_stream << "    chain_length: 15000\n";
        cfg_stream << "    sample_frequency: 10\n";
        cfg_stream << "    operators:\n";
        cfg_stream << "        SplitLumpNodesRevJumpSampler:\n";
        cfg_stream << "            weight: 0\n";

        PopulationTreeSettings settings = PopulationTreeSettings(cfg_stream, cfg_path);

        auto rj_op = settings.operator_settings->tunable_operators["SplitLumpNodesRevJumpSampler"];
        for (unsigned int i = 0; i < rj_op.get_number_of_operators(); ++i) {
            std::cout << "weight " << i << ": " << rj_op.get_weight(i) << "\n";
            std::cout << "tuning param " << i << ": " << rj_op.get_tuning_parameter(i) << "\n";
            std::cout << "auto opt " << i << ": " << rj_op.auto_optimizing(i) << "\n";
            std::cout << "auto opt delay " << i << ": " << rj_op.get_auto_optimize_delay(i) << "\n";
        }

        GeneralTreeOperatorSchedule<BasePopulationTree> op_schedule(
                settings.operator_settings,
                9); // number of leaves

        op_schedule.write_op_settings(std::cout);

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<BasePopulationTree> > > ops;

        ops = op_schedule.get_split_lump_rj_operators();
        REQUIRE(ops.empty());
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler operator turn off",
        "[GeneralTreeOperatorSchedule]") {
    SECTION("Testing SplitLumpNodesRevJumpSampler turn off") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "---\n";
        cfg_stream << "tree_model:\n";
        cfg_stream << "    tree_space: \"bifurcating\"\n";
        cfg_stream << "    tree_prior:\n";
        cfg_stream << "        uniform_root_and_betas:\n";
        cfg_stream << "            parameters:\n";
        cfg_stream << "                root_height:\n";
        cfg_stream << "                    value: 0.3\n";
        cfg_stream << "                    estimate: true\n";
        cfg_stream << "                    prior:\n";
        cfg_stream << "                        gamma_distribution:\n";
        cfg_stream << "                            shape: 8.0\n";
        cfg_stream << "                            mean: 0.3\n";
        cfg_stream << "                alpha_of_node_height_beta_prior:\n";
        cfg_stream << "                    value: 1.0\n";
        cfg_stream << "                    estimate: true\n";
        cfg_stream << "                    prior:\n";
        cfg_stream << "                        gamma_distribution:\n";
        cfg_stream << "                            shape: 4.0\n";
        cfg_stream << "                            mean: 1.0\n";
        cfg_stream << "data:\n";
        cfg_stream << "    ploidy: 1\n";
        cfg_stream << "    constant_sites_removed: false\n";
        cfg_stream << "    yaml_allele_counts:\n";
        cfg_stream << "        path: \"diploid-dna-constant-missing.yml\"\n";
        cfg_stream << "mcmc_settings:\n";
        cfg_stream << "    chain_length: 15000\n";
        cfg_stream << "    sample_frequency: 10\n";
        cfg_stream << "    operators:\n";
        cfg_stream << "        SplitLumpNodesRevJumpSampler:\n";
        cfg_stream << "            weight: [2, 3, 5]\n";
        cfg_stream << "            tuning_parameter: [2.0, 1.0, 10.0]\n";
        cfg_stream << "            auto_optimize: [true, false, false]\n";
        cfg_stream << "            auto_optimize_delay: [100, 50, 20]\n";

        PopulationTreeSettings settings = PopulationTreeSettings(cfg_stream, cfg_path);

        auto rj_op = settings.operator_settings->tunable_operators["SplitLumpNodesRevJumpSampler"];
        for (unsigned int i = 0; i < rj_op.get_number_of_operators(); ++i) {
            std::cout << "weight " << i << ": " << rj_op.get_weight(i) << "\n";
            std::cout << "tuning param " << i << ": " << rj_op.get_tuning_parameter(i) << "\n";
            std::cout << "auto opt " << i << ": " << rj_op.auto_optimizing(i) << "\n";
            std::cout << "auto opt delay " << i << ": " << rj_op.get_auto_optimize_delay(i) << "\n";
        }

        GeneralTreeOperatorSchedule<BasePopulationTree> op_schedule(
                settings.operator_settings,
                5); // number of leaves

        op_schedule.write_op_settings(std::cout);

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<BasePopulationTree> > > ops;

        ops = op_schedule.get_split_lump_rj_operators();
        REQUIRE(ops.empty());
    }
}

TEST_CASE("Testing bifurcating BasePopulationTree config",
        "[GeneralTreeOperatorSchedule]") {
    SECTION("Testing with uniform_root_and_betas tree prior") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "---\n";
        cfg_stream << "tree_model:\n";
        cfg_stream << "    tree_space: \"bifurcating\"\n";
        cfg_stream << "    tree_prior:\n";
        cfg_stream << "        uniform_root_and_betas:\n";
        cfg_stream << "            parameters:\n";
        cfg_stream << "                root_height:\n";
        cfg_stream << "                    value: 0.3\n";
        cfg_stream << "                    estimate: true\n";
        cfg_stream << "                    prior:\n";
        cfg_stream << "                        gamma_distribution:\n";
        cfg_stream << "                            shape: 8.0\n";
        cfg_stream << "                            mean: 0.3\n";
        cfg_stream << "                alpha_of_node_height_beta_prior:\n";
        cfg_stream << "                    value: 1.0\n";
        cfg_stream << "                    estimate: true\n";
        cfg_stream << "                    prior:\n";
        cfg_stream << "                        gamma_distribution:\n";
        cfg_stream << "                            shape: 4.0\n";
        cfg_stream << "                            mean: 1.0\n";
        cfg_stream << "data:\n";
        cfg_stream << "    ploidy: 1\n";
        cfg_stream << "    constant_sites_removed: false\n";
        cfg_stream << "    yaml_allele_counts:\n";
        cfg_stream << "        path: \"diploid-dna-constant-missing.yml\"\n";

        PopulationTreeSettings settings = PopulationTreeSettings(cfg_stream, cfg_path);

        GeneralTreeOperatorSchedule<BasePopulationTree> op_schedule(
                settings.operator_settings,
                5); // number of leaves

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<BasePopulationTree> > > ops;

        ops = op_schedule.get_split_lump_rj_operators();
        REQUIRE(ops.empty());

        ops.clear();

        ops = op_schedule.get_operators(BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::topology);
        REQUIRE(ops.size() > 0);
    }
}

TEST_CASE("Testing fixed tree BasePopulationTree config",
        "[GeneralTreeOperatorSchedule]") {
    SECTION("Testing with uniform_root_and_betas tree prior") {
        std::string cfg_path = "data/dummy.yml";

        std::stringstream cfg_stream;
        cfg_stream << "---\n";
        cfg_stream << "tree_model:\n";
        cfg_stream << "    tree_space: \"fixed\"\n";
        cfg_stream << "    starting_tree: \"[&R]((A:0.1,B:0.1):0.1,C:0.2):0.0\"\n";
        cfg_stream << "    tree_prior:\n";
        cfg_stream << "        uniform_root_and_betas:\n";
        cfg_stream << "            parameters:\n";
        cfg_stream << "                root_height:\n";
        cfg_stream << "                    estimate: true\n";
        cfg_stream << "                    prior:\n";
        cfg_stream << "                        gamma_distribution:\n";
        cfg_stream << "                            shape: 8.0\n";
        cfg_stream << "                            mean: 0.3\n";
        cfg_stream << "                alpha_of_node_height_beta_prior:\n";
        cfg_stream << "                    value: 1.0\n";
        cfg_stream << "                    estimate: true\n";
        cfg_stream << "                    prior:\n";
        cfg_stream << "                        gamma_distribution:\n";
        cfg_stream << "                            shape: 4.0\n";
        cfg_stream << "                            mean: 1.0\n";
        cfg_stream << "data:\n";
        cfg_stream << "    ploidy: 1\n";
        cfg_stream << "    constant_sites_removed: false\n";
        cfg_stream << "    yaml_allele_counts:\n";
        cfg_stream << "        path: \"diploid-dna-constant-missing.yml\"\n";

        PopulationTreeSettings settings = PopulationTreeSettings(cfg_stream, cfg_path);

        GeneralTreeOperatorSchedule<BasePopulationTree> op_schedule(
                settings.operator_settings,
                5); // number of leaves

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<BasePopulationTree> > > ops;

        ops = op_schedule.get_split_lump_rj_operators();
        REQUIRE(ops.empty());

        ops.clear();

        ops = op_schedule.get_operators(BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::topology);
        REQUIRE(ops.size() == 0);
    }
}
