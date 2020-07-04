#include "catch.hpp"

#include "ecoevolity/sumphycoeval.hpp"
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

        // ret = simphycoeval_main<BasePopulationTree>(argc, argv);
        // REQUIRE(ret == 0);

        REQUIRE_THROWS_AS(
                (sumphycoeval_main<BasePopulationTree>(argc, argv)),
                EcoevolityError &);
    }
}
