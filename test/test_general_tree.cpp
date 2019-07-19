#include "catch.hpp"
#include "ecoevolity/general_tree.hpp"

TEST_CASE("Testing likelihood of GeneralPopulationTree with three-way polytomy at root", "[GeneralPopulationTree]") {

    SECTION("Testing constructor and likelihood calc") {
        std::shared_ptr<PopulationNode> root0 = std::make_shared<PopulationNode>(3, "root", 0.1);
        std::shared_ptr<PopulationNode> leaf00 = std::make_shared<PopulationNode>(0, "leaf 0", 0.0, 4);
        leaf00->fix_node_height();
        std::shared_ptr<PopulationNode> leaf01 = std::make_shared<PopulationNode>(1, "leaf 1", 0.0, 4);
        leaf01->fix_node_height();
        std::shared_ptr<PopulationNode> leaf02 = std::make_shared<PopulationNode>(2, "leaf 2", 0.0, 4);
        leaf02->fix_node_height();

        root0->add_child(leaf00);
        root0->add_child(leaf01);
        root0->add_child(leaf02);

        std::shared_ptr<PopulationNode> internal1 = std::make_shared<PopulationNode>(3, "internal 0", 0.1);
        std::shared_ptr<PopulationNode> root1 = std::make_shared<PopulationNode>(4, "root", 0.1);
        std::shared_ptr<PopulationNode> leaf10 = std::make_shared<PopulationNode>(0, "leaf 0", 0.0, 4);
        leaf10->fix_node_height();
        std::shared_ptr<PopulationNode> leaf11 = std::make_shared<PopulationNode>(1, "leaf 1", 0.0, 4);
        leaf11->fix_node_height();
        std::shared_ptr<PopulationNode> leaf12 = std::make_shared<PopulationNode>(2, "leaf 2", 0.0, 4);
        leaf12->fix_node_height();

        internal1->add_child(leaf10);
        internal1->add_child(leaf11);
        root1->add_child(internal1);
        root1->add_child(leaf12);


        GeneralPopulationTree tree0(root0,
                100,   // number of loci
                1,     // length of loci
                true); // validate data
        GeneralPopulationTree tree1(root1,
                100,   // number of loci
                1,     // length of loci
                true); // validate data

        REQUIRE(tree0.get_number_of_node_heights() == 1);
        REQUIRE(tree1.get_number_of_node_heights() == 2);

        const std::vector< std::shared_ptr<PositiveRealParameter> >& tree0_heights = tree0.get_node_height_parameters();
        REQUIRE(tree0_heights.at(0) == root0->get_height_parameter());

        const std::vector< std::shared_ptr<PositiveRealParameter> >& tree1_heights = tree1.get_node_height_parameters();
        std::shared_ptr<PositiveRealParameter> ht = internal1->get_height_parameter();
        bool found = false;
        for (unsigned int i = 0; i < tree1_heights.size(); ++i) {
            if (tree1_heights.at(i) == ht) {
                found = true;
            }
        }
        REQUIRE(found);
        

        tree0.set_all_population_sizes(0.005);
        tree1.set_all_population_sizes(0.005);

        double l0 = tree0.compute_log_likelihood();
        double l1 = tree1.compute_log_likelihood();
        std::cout << "\nPolytomy lnL: " << l0;
        std::cout << "\nBifurcating lnL: " << l1 << "\n\n";
        REQUIRE(l0 == l1);

        RandomNumberGenerator rng = RandomNumberGenerator(1234);
        BiallelicData bd = tree0.simulate_linked_biallelic_data_set(rng,
                1.0,    // singleton sample probability
                false,  // max one variable site per locus
                true);  // validate data set

        tree0.set_data(bd, false);
        tree1.set_data(bd, false);

        l0 = tree0.compute_log_likelihood();
        l1 = tree1.compute_log_likelihood();
        std::cout << "\nPolytomy lnL: " << l0;
        std::cout << "\nBifurcating lnL: " << l1 << "\n\n";
        REQUIRE(l0 == l1);
    }
}
