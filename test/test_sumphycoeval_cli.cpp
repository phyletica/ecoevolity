#include "catch.hpp"

#include "ecoevolity/sumphycoeval.hpp"
#include "ecoevolity/node.hpp"
#include "ecoevolity/treesum.hpp"
#include "ecoevolity/path.hpp"


TEST_CASE("Testing missing target tree", "[sumphycoeval]") {

    SECTION("Testing missing target tree") {
        char source1[] = "data/4-tip-trees-12-34.nex";
        char source2[] = "data/4-tip-trees-12.nex";
        char source3[] = "data/4-tip-trees-13-24.nex";
        char source4[] = "data/4-tip-trees-14-23-shared.nex";
        char source5[] = "data/4-tip-trees-34.nex";
        char source6[] = "data/4-tip-trees-comb.nex";
        char source7[] = "data/4-tip-trees-ladder-1234.nex";
        char source8[] = "data/4-tip-trees-ladder-4321.nex";

        char t_out_path[] = "data/tmp-target-tree-out-phyco.nex";
        char m_out_path[] = "data/tmp-map-tree-out-phyco.nex";

        char exe[] = "sumphycoeval";
        char burnin_flag[] = "--burnin";
        char burnin[] = "2";
        char t_flag[] = "--target-tree";
        char t_out_flag[] = "--target-tree-out";
        char m_out_flag[] = "--map-tree-out";
        char force[] = "--force";

        char * argv[] = {
            &exe[0],
            &force[0],
            &burnin_flag[0],
            &burnin[0],
            &t_out_flag[0],
            &t_out_path[0],
            &source1[0],
            &source2[0],
            &source3[0],
            &source4[0],
            &source5[0],
            &source6[0],
            &source7[0],
            &source8[0],
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        // int ret;

        // ret = sumphycoeval_main<PopulationNode>(argc, argv);
        // REQUIRE(ret == 0);

        REQUIRE_THROWS_AS(
                (sumphycoeval_main<PopulationNode>(argc, argv)),
                EcoevolityError &);
    }
}

TEST_CASE("Testing bad freq argument", "[sumphycoeval]") {

    SECTION("Testing bad freq argument") {
        char source1[] = "data/4-tip-trees-12-34.nex";
        char source2[] = "data/4-tip-trees-12.nex";
        char source3[] = "data/4-tip-trees-13-24.nex";
        char source4[] = "data/4-tip-trees-14-23-shared.nex";
        char source5[] = "data/4-tip-trees-34.nex";
        char source6[] = "data/4-tip-trees-comb.nex";
        char source7[] = "data/4-tip-trees-ladder-1234.nex";
        char source8[] = "data/4-tip-trees-ladder-4321.nex";

        char t_out_path[] = "data/tmp-target-tree-out-phyco.nex";
        char m_out_path[] = "data/tmp-map-tree-out-phyco.nex";

        char exe[] = "sumphycoeval";
        char burnin_flag[] = "--burnin";
        char burnin[] = "2";
        char t_flag[] = "--target-tree";
        char t_out_flag[] = "--target-tree-out";
        char m_out_flag[] = "--map-tree-out";
        char sf_flag[] = "--min-split-freq";
        char sf[] = "1.1";
        char force[] = "--force";

        char * argv[] = {
            &exe[0],
            &force[0],
            &burnin_flag[0],
            &burnin[0],
            &sf_flag[0],
            &sf[0],
            &source1[0],
            &source2[0],
            &source3[0],
            &source4[0],
            &source5[0],
            &source6[0],
            &source7[0],
            &source8[0],
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        // int ret;

        // ret = sumphycoeval_main<PopulationNode>(argc, argv);
        // REQUIRE(ret == 0);

        REQUIRE_THROWS_AS(
                (sumphycoeval_main<PopulationNode>(argc, argv)),
                EcoevolityError &);
    }
}


TEST_CASE("Testing working example", "[sumphycoeval]") {

    SECTION("Testing working example") {
        char source1[] = "data/4-tip-trees-12-34.nex";
        char source2[] = "data/4-tip-trees-12.nex";
        char source3[] = "data/4-tip-trees-13-24.nex";
        char source4[] = "data/4-tip-trees-14-23-shared.nex";
        char source5[] = "data/4-tip-trees-34.nex";
        char source6[] = "data/4-tip-trees-comb.nex";
        char source7[] = "data/4-tip-trees-ladder-1234.nex";
        char source8[] = "data/4-tip-trees-ladder-4321.nex";
        char t_path[] = "data/4-tip-target-tree-14-23-shared.nex";

        char t_out_path[] = "data/tmp-target-tree-out-sumphyco1.nex";
        char m_out_path[] = "data/tmp-map-tree-out-sumphyco1.nex";

        char exe[] = "sumphycoeval";
        char burnin_flag[] = "--burnin";
        char burnin[] = "2";
        char t_flag[] = "--target-tree";
        char t_out_flag[] = "--target-tree-out";
        char m_out_flag[] = "--map-tree-out";
        char sf_flag[] = "--min-split-freq";
        char sf[] = "0.1";
        char force[] = "--force";

        char * argv[] = {
            &exe[0],
            &force[0],
            &burnin_flag[0],
            &burnin[0],
            &sf_flag[0],
            &sf[0],
            &t_flag[0],
            &t_path[0],
            &t_out_flag[0],
            &t_out_path[0],
            &m_out_flag[0],
            &m_out_path[0],
            &source1[0],
            &source2[0],
            &source3[0],
            &source4[0],
            &source5[0],
            &source6[0],
            &source7[0],
            &source8[0],
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        ret = sumphycoeval_main<PopulationNode>(argc, argv);
        REQUIRE(ret == 0);

        REQUIRE(path::exists(t_out_path));
        REQUIRE(path::exists(m_out_path));

        std::vector<std::string> source_tree_paths {
                "data/4-tip-trees-12-34.nex",
                "data/4-tip-trees-12.nex",
                "data/4-tip-trees-13-24.nex",
                "data/4-tip-trees-14-23-shared.nex",
                "data/4-tip-trees-34.nex",
                "data/4-tip-trees-comb.nex",
                "data/4-tip-trees-ladder-1234.nex",
                "data/4-tip-trees-ladder-4321.nex",
        };
        std::string target_tree_path = "data/4-tip-target-tree-14-23-shared.nex";

        treesum::TreeSample<PopulationNode> ts(
                target_tree_path,
                source_tree_paths,
                "nexus",
                2);
        std::stringstream expected_target_stream;
        std::stringstream expected_map_stream;
        ts.write_map_trees_to_nexus(expected_target_stream, false, 18);
        ts.write_target_tree_to_nexus(expected_map_stream, false, 18);

        std::stringstream target_stream;
        std::ifstream target_file_in_stream(t_out_path);
        if (target_file_in_stream) {
            target_stream << target_file_in_stream.rdbuf();
            target_file_in_stream.close();
        }

        std::stringstream map_stream;
        std::ifstream map_file_in_stream(m_out_path);
        if (map_file_in_stream) {
            map_stream << map_file_in_stream.rdbuf();
            map_file_in_stream.close();
        }

        REQUIRE(expected_target_stream.str() == target_stream.str());
        REQUIRE(expected_map_stream.str() == map_stream.str());

    }
}
