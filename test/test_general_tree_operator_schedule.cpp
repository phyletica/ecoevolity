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
        std::shared_ptr< GeneralTreeOperatorTemplate<BasePopulationTree> > op;

        op = op_schedule.get_split_lump_rj_operator();
        REQUIRE(op);

        REQUIRE(op->helper_ops.size() > 0);

        ops = op_schedule.get_operators(BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::topology);
        REQUIRE(ops.size() > 0);
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
        std::shared_ptr< GeneralTreeOperatorTemplate<BasePopulationTree> > op;

        op = op_schedule.get_split_lump_rj_operator();
        REQUIRE(! op);

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
        std::shared_ptr< GeneralTreeOperatorTemplate<BasePopulationTree> > op;

        op = op_schedule.get_split_lump_rj_operator();
        REQUIRE(! op);

        ops = op_schedule.get_operators(BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::topology);
        REQUIRE(ops.size() == 0);
    }
}
