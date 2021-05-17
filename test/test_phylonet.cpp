#include "catch.hpp"
#include "ecoevolity/phylonet.hpp"
#include "ecoevolity/stats_util.hpp"


TEST_CASE("phylonet; Testing simulations against likelihood for one pop with 2,1 pattern",
        "[BasePopulationNetwork]") {

    SECTION("Testing sims v likelihood for one pop with 2 alleles") {
        double pop_size = 0.1;

        std::string nex_path = "data/singleton-2-1.nex";
        // Need to keep constant characters
        BasePopulationNetwork tree(nex_path, ' ',
                true,   // population_name_is_prefix
                false,  // diploid
                false,  // dominant
                false,  // constant sites removed
                true,   // validate
                true,  // strict on constant
                true,  // strict on missing
                true,  // strict on triallelic
                2.0,    // ploidy
                false    // store charset info
                );
        REQUIRE(tree.get_leaf_node_count() == 1);
        tree.estimate_mutation_rate();

        tree.set_all_population_sizes(pop_size);

        tree.set_root_height(0.01);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        double ln_l = tree.compute_log_likelihood();
        double ln_l_correction = tree.get_likelihood_correction();
        double raw_ln_l = ln_l - ln_l_correction;
        double l = std::exp(ln_l);
        double raw_l = std::exp(raw_ln_l);
        double l_correction = std::exp(ln_l_correction);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        unsigned int nsamples = 1000000;

        unsigned int het_count = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            auto pattern_tree = tree.simulate_biallelic_site(
                    0,
                    rng,
                    false);
            auto pattern = pattern_tree.first;
            /* std::vector<unsigned int> red_allele_counts = pattern.first; */
            /* std::vector<unsigned int> allele_counts = pattern.second; */
            /* if (red_allele_counts.at(0) == 1) { */
            if (pattern.first.at(0) == 1) {
                ++het_count;
            }
        }
        double approx_l = het_count / (double)nsamples;
        std::cout << "approx like: " << approx_l << "\n";
        std::cout << "like: " << l << "\n";
        std::cout << "like correction: " << l_correction << "\n";
        std::cout << "raw like: " << raw_l << "\n";
        REQUIRE(approx_l == Approx(raw_l).epsilon(0.001));
    }
}

TEST_CASE("phylonet; Testing simulations against likelihood for one pop with 3,1 pattern",
        "[BasePopulationNetwork]") {

    SECTION("Testing sims v likelihood for one pop with 3 alleles") {
        double pop_size = 0.1;

        std::string nex_path = "data/singleton-3-1.nex";
        // Need to keep constant characters
        BasePopulationNetwork tree(nex_path, ' ',
                true,   // population_name_is_prefix
                false,  // diploid
                false,  // dominant
                false,  // constant sites removed
                true,   // validate
                true,  // strict on constant
                true,  // strict on missing
                true,  // strict on triallelic
                2.0,    // ploidy
                false    // store charset info
                );
        REQUIRE(tree.get_leaf_node_count() == 1);
        tree.estimate_mutation_rate();

        tree.set_all_population_sizes(pop_size);

        tree.set_root_height(0.01);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        double ln_l = tree.compute_log_likelihood();
        double ln_l_correction = tree.get_likelihood_correction();
        double raw_ln_l = ln_l - ln_l_correction;
        double l = std::exp(ln_l);
        double raw_l = std::exp(raw_ln_l);
        double l_correction = std::exp(ln_l_correction);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        unsigned int nsamples = 1000000;

        unsigned int het_count = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            auto pattern_tree = tree.simulate_biallelic_site(
                    0,
                    rng,
                    false);
            auto pattern = pattern_tree.first;
            /* std::vector<unsigned int> red_allele_counts = pattern.first; */
            /* std::vector<unsigned int> allele_counts = pattern.second; */
            /* if (red_allele_counts.at(0) == 1) { */
            if (pattern.first.at(0) == 1) {
                ++het_count;
            }
        }
        double approx_l = het_count / (double)nsamples;
        std::cout << "approx like: " << approx_l << "\n";
        std::cout << "like: " << l << "\n";
        std::cout << "like correction: " << l_correction << "\n";
        std::cout << "raw like: " << raw_l << "\n";
        REQUIRE(approx_l == Approx(raw_l).epsilon(0.001));
    }
}

TEST_CASE("phylonet; Testing simulations against likelihood for 2,1;2,2 pattern",
        "[BasePopulationNetwork]") {

    SECTION("Testing sims v likelihood for 2 pops with 2,1;2,2") {
        double pop_size = 0.1;

        std::string nex_path = "data/s2_2-1_2-2.nex";
        // Need to keep constant characters
        BasePopulationNetwork tree(nex_path, ' ',
                true,   // population_name_is_prefix
                false,  // diploid
                false,  // dominant
                false,  // constant sites removed
                true,   // validate
                true,  // strict on constant
                true,  // strict on missing
                true,  // strict on triallelic
                2.0,    // ploidy
                false    // store charset info
                );
        REQUIRE(tree.get_leaf_node_count() == 2);
        REQUIRE(tree.get_node_count() == 3);
        tree.estimate_mutation_rate();

        tree.set_all_population_sizes(pop_size);

        tree.set_root_height(0.02);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        double ln_l = tree.compute_log_likelihood();
        double ln_l_correction = tree.get_likelihood_correction();
        double raw_ln_l = ln_l - ln_l_correction;
        double l = std::exp(ln_l);
        double raw_l = std::exp(raw_ln_l);
        double l_correction = std::exp(ln_l_correction);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        unsigned int nsamples = 1000000;

        std::vector<unsigned int> allele_counts = tree.get_data().get_allele_counts(0);
        std::vector<unsigned int> red_allele_counts = tree.get_data().get_red_allele_counts(0);

        unsigned int pattern_count = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            auto sim_pattern_tree = tree.simulate_biallelic_site(
                    0,
                    rng,
                    false);
            auto sim_pattern = sim_pattern_tree.first;
            REQUIRE(sim_pattern.second == allele_counts);
            /* std::vector<unsigned int> red_allele_counts = pattern.first; */
            /* std::vector<unsigned int> allele_counts = pattern.second; */
            /* if (red_allele_counts.at(0) == 1) { */
            if (sim_pattern.first == red_allele_counts) {
                ++pattern_count;
            }
        }
        double approx_l = pattern_count / (double)nsamples;
        std::cout << "approx like: " << approx_l << "\n";
        std::cout << "like: " << l << "\n";
        std::cout << "like correction: " << l_correction << "\n";
        std::cout << "raw like: " << raw_l << "\n";
        REQUIRE(approx_l == Approx(raw_l).epsilon(0.001));
    }
}

TEST_CASE("phylonet; Testing simulations against likelihood for network with 2,1;2,2 pattern",
        /* "[BasePopulationNetwork]") { */
        "[xx]") {

    SECTION("Testing sims v likelihood for 2 retic pops with 2,1;2,2") {
        double pop_size = 0.1;

        std::string nex_path = "data/s2_2-1_2-2.nex";
        // Need to keep constant characters
        BasePopulationNetwork default_tree(nex_path, ' ',
                true,   // population_name_is_prefix
                false,  // diploid
                false,  // dominant
                false,  // constant sites removed
                true,   // validate
                true,  // strict on constant
                true,  // strict on missing
                true,  // strict on triallelic
                2.0,    // ploidy
                false    // store charset info
                );

        BiallelicData data = default_tree.get_data();

        std::shared_ptr<PopulationNetNode> root = default_tree.get_root_ptr()->get_clade_clone();
        std::shared_ptr<PopulationNetNode> troot = default_tree.get_root_ptr()->get_clade_clone();
        REQUIRE(root->get_leaf_node_count() == 2);
        REQUIRE(root->get_node_count() == 3);
        std::shared_ptr<PopulationNetNode> child0 = root->get_child(0);        
        std::shared_ptr<PopulationNetNode> child1 = root->get_child(1);        
        std::shared_ptr<PopulationNetNode> tchild0 = troot->get_child(0);        
        std::shared_ptr<PopulationNetNode> tchild1 = troot->get_child(1);        

        root->set_height(0.1);
        root->set_label("root");
        troot->set_height(0.1);
        troot->set_label("troot");
        std::shared_ptr<PopulationNetNode> internal0 = std::make_shared<PopulationNetNode>(3, "internal0", 0.02);
        std::shared_ptr<PopulationNetNode> internal1 = std::make_shared<PopulationNetNode>(4, "internal1", 0.04);
        std::shared_ptr<PopulationNetNode> tinternal0 = std::make_shared<PopulationNetNode>(3, "tinternal0", 0.02);
        std::shared_ptr<PopulationNetNode> tinternal1 = std::make_shared<PopulationNetNode>(4, "tinternal1", 0.04);
        internal0->add_child(child0);
        tinternal0->add_child(tchild0);
        internal0->add_parent(root);
        tinternal0->add_parent(troot);
        internal0->add_parent(internal1);
        internal1->add_child(child1);
        tinternal1->add_child(tchild1);
        internal1->add_parent(root);
        tinternal1->add_parent(troot);
        child0->remove_parent(root);
        child1->remove_parent(root);
        tchild0->remove_parent(troot);
        tchild1->remove_parent(troot);

        REQUIRE(root->get_leaf_node_count() == 2);
        REQUIRE(root->get_node_count() == 5);
        REQUIRE(troot->get_leaf_node_count() == 2);
        REQUIRE(troot->get_node_count() == 5);


        REQUIRE(internal0->get_number_of_children() == 1);
        REQUIRE(internal0->get_number_of_parents() == 2);
        REQUIRE(internal1->get_number_of_children() == 2);
        REQUIRE(internal1->get_number_of_parents() == 1);
        REQUIRE(root->get_number_of_children() == 2);
        REQUIRE(root->get_number_of_parents() == 0);
        REQUIRE(child0->get_number_of_children() == 0);
        REQUIRE(child0->get_number_of_parents() == 1);
        REQUIRE(child1->get_number_of_children() == 0);
        REQUIRE(child1->get_number_of_parents() == 1);

        REQUIRE(root->is_root());
        REQUIRE(! root->is_leaf());
        REQUIRE(! child0->is_root());
        REQUIRE(child0->is_leaf());
        REQUIRE(! child1->is_root());
        REQUIRE(child1->is_leaf());
        REQUIRE(! internal0->is_root());
        REQUIRE(! internal0->is_leaf());
        REQUIRE(! internal1->is_root());
        REQUIRE(! internal1->is_leaf());

        REQUIRE(root->is_child(internal0));
        REQUIRE(root->is_child(internal1));
        REQUIRE(! root->is_child(child0));
        REQUIRE(! root->is_child(child1));
        REQUIRE(internal0->is_child(child0));
        REQUIRE(! internal0->is_child(child1));
        REQUIRE(! internal0->is_child(root));
        REQUIRE(internal1->is_child(child1));
        REQUIRE(internal1->is_child(internal0));
        REQUIRE(! internal1->is_child(child0));
        REQUIRE(! internal1->is_child(root));

        REQUIRE(internal0->is_parent(root));
        REQUIRE(internal0->is_parent(internal1));
        REQUIRE(internal1->is_parent(root));
        REQUIRE(child0->is_parent(internal0));
        REQUIRE(! child0->is_parent(root));
        REQUIRE(child1->is_parent(internal1));
        REQUIRE(! child1->is_parent(root));

        root->set_all_population_sizes(pop_size);
        troot->set_all_population_sizes(pop_size);

        std::cout << root->to_parentheses(2, true) << "\n";

        BasePopulationNetwork tree(root);
        tree.set_data(data, false);

        BasePopulationNetwork ttree(troot);
        ttree.set_data(data, false);

        std::vector<double> expected_heights = {0.02, 0.04, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        /* std::vector<double> texpected_heights = {0.1}; */
        REQUIRE(ttree.get_node_heights() == expected_heights);

        // Turning off gene flow
        internal0->set_inheritance_proportion(0, 1.0);


        tree.estimate_mutation_rate();
        tree.set_all_population_sizes(pop_size);
        tree.set_freq_1(0.5);
        tree.set_mutation_rate(1.0);

        ttree.estimate_mutation_rate();
        ttree.set_all_population_sizes(pop_size);
        ttree.set_freq_1(0.5);
        ttree.set_mutation_rate(1.0);

        std::cout << "Net likelihood:\n";
        double ln_l = tree.compute_log_likelihood();
        double ln_l_correction = tree.get_likelihood_correction();
        double raw_ln_l = ln_l - ln_l_correction;
        double l = std::exp(ln_l);
        double raw_l = std::exp(raw_ln_l);
        double l_correction = std::exp(ln_l_correction);

        std::cout << "Tree likelihood:\n";
        double tln_l = ttree.compute_log_likelihood();
        double tln_l_correction = ttree.get_likelihood_correction();
        double traw_ln_l = tln_l - tln_l_correction;
        double tl = std::exp(tln_l);
        double traw_l = std::exp(traw_ln_l);
        double tl_correction = std::exp(tln_l_correction);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        unsigned int nsamples = 1000000;

        std::vector<unsigned int> allele_counts = tree.get_data().get_allele_counts(0);
        std::vector<unsigned int> red_allele_counts = tree.get_data().get_red_allele_counts(0);

        std::vector<unsigned int> tallele_counts = ttree.get_data().get_allele_counts(0);
        std::vector<unsigned int> tred_allele_counts = ttree.get_data().get_red_allele_counts(0);

        REQUIRE(allele_counts == tallele_counts);
        REQUIRE(red_allele_counts == tred_allele_counts);

        unsigned int pattern_count = 0;
        unsigned int tpattern_count = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            auto sim_pattern_tree = tree.simulate_biallelic_site(
                    0,
                    rng,
                    false);
            auto sim_pattern = sim_pattern_tree.first;
            auto tsim_pattern_tree = ttree.simulate_biallelic_site(
                    0,
                    rng,
                    false);
            auto tsim_pattern = tsim_pattern_tree.first;
            REQUIRE(sim_pattern.second == allele_counts);
            REQUIRE(tsim_pattern.second == allele_counts);
            /* std::vector<unsigned int> red_allele_counts = pattern.first; */
            /* std::vector<unsigned int> allele_counts = pattern.second; */
            /* if (red_allele_counts.at(0) == 1) { */
            if (sim_pattern.first == red_allele_counts) {
                ++pattern_count;
            }
            if (tsim_pattern.first == red_allele_counts) {
                ++tpattern_count;
            }
        }
        double approx_l = pattern_count / (double)nsamples;
        double tapprox_l = tpattern_count / (double)nsamples;
        std::cout << "approx like: " << approx_l << "\n";
        std::cout << "approx tree like: " << tapprox_l << "\n";
        std::cout << "like: " << l << "\n";
        std::cout << "tree like: " << tl << "\n";
        std::cout << "like correction: " << l_correction << "\n";
        std::cout << "tree like correction: " << tl_correction << "\n";
        std::cout << "raw like: " << raw_l << "\n";
        std::cout << "tree raw like: " << traw_l << "\n";
        REQUIRE(approx_l == Approx(raw_l).epsilon(0.001));
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
TEST_CASE("phylonet; Testing simple likelihood of BasePopulationNetwork", "[BasePopulationNetwork]") {

    SECTION("Testing constructor and likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("phylonet; Testing threaded likelihood of BasePopulationNetwork", "[BasePopulationNetwork]") {

    SECTION("Testing constructor and threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("phylonet; Testing over-threaded likelihood of BasePopulationNetwork", "[BasePopulationNetwork]") {

    SECTION("Testing constructor and over-threaded likelihood calc") {
        std::string nex_path = "data/diploid-standard-data-ntax5-nchar5.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(10);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(10);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-31.77866581319647));
        REQUIRE(tree.get_degree_of_root() == 2);
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
TEST_CASE("phylonet; Testing hemi129.nex likelihood (0.01, 10.0, 1.0, 1.0)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("phylonet; Testing hemi129.nex threaded likelihood (0.01, 10.0, 1.0, 1.0)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-248.93254688526213));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
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
TEST_CASE("phylonet; Testing aflp_25.nex likelihood (0.01, 10.0, 1.0, 1.0)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("phylonet; Testing aflp_25.nex threaded likelihood (0.01, 10.0, 1.0, 1.0)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7099.716015109998));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
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
TEST_CASE("phylonet; Testing aflp_25.nex likelihood (0.01, 10.0, 1.0, 1.0, dominant)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6986.120524781545));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("phylonet; Testing aflp_25.nex threaded likelihood (0.01, 10.0, 1.0, 1.0, dominant)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(5);
        REQUIRE(l == Approx(-6986.120524781545));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
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
TEST_CASE("phylonet; Testing hemi129.nex likelihood (0.0, 10.0, 1.0, 1.0)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("phylonet; Testing hemi129.nex threaded likelihood (0.0, 10.0, 1.0, 1.0)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(7);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(7);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-328.39238828878365));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
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
TEST_CASE("phylonet; Testing aflp_25.nex likelihood (0.0, 10.0, 1.0, 1.0)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("phylonet; Testing aflp_25.nex threaded likelihood (0.0, 10.0, 1.0, 1.0)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-7256.501742344454));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
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
TEST_CASE("phylonet; Testing aflp_25.nex likelihood (0.0, 10.0, 1.0, 1.0, dominant)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7223.362711937651));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("phylonet; Testing aflp_25.nex threaded likelihood (0.0, 10.0, 1.0, 1.0, dominant)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.0);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7223.362711937651));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
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
TEST_CASE("phylonet; Testing hemi129.nex likelihood (0.2, 10.0, 1.0, 1.0)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("phylonet; Testing hemi129.nex threaded likelihood (0.2, 10.0, 1.0, 1.0)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-227.41048391087554));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
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
TEST_CASE("phylonet; Testing aflp_25.nex likelihood (0.2, 10.0, 1.0, 1.0)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("phylonet; Testing aflp_25.nex threaded likelihood (0.2, 10.0, 1.0, 1.0)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-7304.180743441677));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
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
TEST_CASE("phylonet; Testing aflp_25.nex likelihood (0.2, 10.0, 1.0, 1.0, dominant)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-7405.145951634711));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("phylonet; Testing aflp_25.nex threaded likelihood (0.2, 10.0, 1.0, 1.0, dominant)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.2);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(6);
        REQUIRE(l == Approx(-7405.145951634711));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
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
TEST_CASE("phylonet; Testing hemi129.nex likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("phylonet; Testing hemi129.nex threaded likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(2);
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-327.7437811413033));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
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
TEST_CASE("phylonet; Testing aflp_25.nex likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("phylonet; Testing aflp_25.nex threaded likelihood (0.03, 10.0, 10.0, 10.0/19.0)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(5);
        REQUIRE(l == Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(3);
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-6472.856486972301));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
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
TEST_CASE("phylonet; Testing aflp_25.nex likelihood (0.03, 10.0, 10.0, 10.0/19.0, dominant)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-6494.774924871097));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("phylonet; Testing aflp_25.nex threaded likelihood (0.03, 10.0, 10.0, 10.0/19.0, dominant)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.05); // tree.set_u(10.0);
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-6494.774924871097));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
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
TEST_CASE("phylonet; Testing hemi129.nex likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("phylonet; Testing hemi129.nex threaded likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(4);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(1);
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-265.0023534261969));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
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
TEST_CASE("phylonet; Testing aflp_25.nex likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood();
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("phylonet; Testing aflp_25.nex threaded likelihood (0.03, 10.0, 10.0/19.0, 10.0)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(3);
        REQUIRE(l == Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));

        tree.fold_patterns();
        l = tree.compute_log_likelihood(5);
        REQUIRE(! std::isnan(l));
        REQUIRE(! std::isinf(l));
        REQUIRE(l != Approx(-10163.468886613919));
        REQUIRE(tree.get_likelihood_correction(true) == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
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
TEST_CASE("phylonet; Testing aflp_25.nex likelihood (0.03, 10.0, 10.0/19.0, 10.0, dominant)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-10999.288193543642));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
    }
}
TEST_CASE("phylonet; Testing aflp_25.nex threaded likelihood (0.03, 10.0, 10.0/19.0, 10.0, dominant)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_all_population_sizes(2.0/(10.0 * 2 * tree.get_ploidy()));
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-10999.288193543642));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);

        REQUIRE_THROWS_AS(tree.fold_patterns(), EcoevolityBiallelicDataError &);
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
TEST_CASE("phylonet; Testing hemi129.nex likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-224.40177558289847));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("phylonet; Testing hemi129.nex threaded likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/hemi129.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(100000);

        // Now handling correction for constant patterns differently than SNAPP
        // when there are allele count patterns with missing data, so are no
        // longer expected to match
        REQUIRE(l != Approx(-224.40177558289847));
        REQUIRE(tree.get_likelihood_correction() == Approx(-135.97095011239867));
        REQUIRE(tree.get_degree_of_root() == 2);
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
TEST_CASE("phylonet; Testing aflp_25.nex likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-8158.88094671241));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("phylonet; Testing aflp_25.nex threaded likelihood (0.03, 111.1, 10.0/19.0, 10.0)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(1000000);
        REQUIRE(l == Approx(-8158.88094671241));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
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
TEST_CASE("phylonet; Testing aflp_25.nex likelihood (0.03, 111.1, 10.0/19.0, 10.0, dominant)", "[BasePopulationNetwork]") {

    SECTION("Testing likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood();
        REQUIRE(l == Approx(-8034.250341980543));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}
TEST_CASE("phylonet; Testing aflp_25.nex threaded likelihood (0.03, 111.1, 10.0/19.0, 10.0, dominant)", "[BasePopulationNetwork]") {

    SECTION("Testing threaded likelihood calc") {
        std::string nex_path = "data/aflp_25.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, false, true);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.03);
        tree.set_freq_1(0.95); // tree.set_u(10.0/19.0);
        tree.set_all_population_sizes(2.0/(111.1 * 2 * tree.get_ploidy()));
        double l = tree.compute_log_likelihood(4);
        REQUIRE(l == Approx(-8034.250341980543));
        REQUIRE(tree.get_likelihood_correction() == Approx(-3317.567573476714));
        REQUIRE(tree.get_degree_of_root() == 2);
    }
}

// BEAST v2.4.0 (master 9fb5b3b7817a0bd4a21e3f90132f132cca72ce4e)
// SNAPP v1.3.0 (master 4f3f0f7366798f4fb38b766c15f6426a75ddf71e)
// diploid-standard-data-ntax5-nchar5.nex
// height = 0.03
// coalescent_rate = 10.0, 10.0, 10.0
// u = 10.0
// v = 10.0/19.0
// Log likelihood            = -23.81984255023975
// Log likelihood correction = -6.87935580446044
//
// With constant sites inclucded and m_bUseNonPolymorphic = true
// Log likelihood            = -55.01646493341547
// Log likelihood correction = -6.87935580446044
TEST_CASE("phylonet; Testing affect of constant sites on likelihood of BasePopulationNetwork", "[BasePopulationNetwork]") {

    SECTION("Testing haploid-standard-full-constant.nex") {
        std::string nex_path = "data/haploid-standard-full-constant.nex";
        BasePopulationNetwork t_included(
                nex_path, // path
                ' ',      // pop name delimiter
                true,     // pop name is prefix
                false,    // genotypes are diploid
                false,    // markers are dominant
                false,    // constant sites removed
                true);    // validate
        std::string nex_path2 = "data/haploid-standard-full-constant-removed.nex";
        BasePopulationNetwork t(
                nex_path2, // path
                ' ',       // pop name delimiter
                true,      // pop name is prefix
                false,     // genotypes are diploid
                false,     // markers are dominant
                true,      // constant sites removed
                true);     // validate

        t.set_root_height(0.03);
        t.set_freq_1(0.05);
        t.set_all_population_sizes(2.0/(10.0 * 2 * t.get_ploidy()));

        t_included.set_root_height(0.03);
        t_included.set_freq_1(0.05);
        t_included.set_all_population_sizes(2.0/(10.0 * 2 * t_included.get_ploidy()));

        double l = t.compute_log_likelihood();
        REQUIRE(l == Approx(-23.81984255023975));
        REQUIRE(t.get_likelihood_correction() == Approx(-6.87935580446044));

        double l_included = t_included.compute_log_likelihood();

        REQUIRE(l_included == Approx(-55.01646493341547));

        REQUIRE(t.get_likelihood_correction() == t_included.get_likelihood_correction());
    }
}

TEST_CASE("phylonet; Testing affect of constant sites on threaded likelihood of BasePopulationNetwork", "[BasePopulationNetwork]") {

    SECTION("Testing haploid-standard-full-constant.nex") {
        std::string nex_path = "data/haploid-standard-full-constant.nex";
        BasePopulationNetwork t_included(
                nex_path, // path
                ' ',      // pop name delimiter
                true,     // pop name is prefix
                false,    // genotypes are diploid
                false,    // markers are dominant
                false,    // constant sites removed
                true);    // validate
        std::string nex_path2 = "data/haploid-standard-full-constant-removed.nex";
        BasePopulationNetwork t(
                nex_path2, // path
                ' ',       // pop name delimiter
                true,      // pop name is prefix
                false,     // genotypes are diploid
                false,     // markers are dominant
                true,      // constant sites removed
                true);     // validate

        t.set_root_height(0.03);
        t.set_freq_1(0.05);
        t.set_all_population_sizes(2.0/(10.0 * 2 * t.get_ploidy()));

        t_included.set_root_height(0.03);
        t_included.set_freq_1(0.05);
        t_included.set_all_population_sizes(2.0/(10.0 * 2 * t_included.get_ploidy()));

        double l = t.compute_log_likelihood(4);
        REQUIRE(l == Approx(-23.81984255023975));
        REQUIRE(t.get_likelihood_correction() == Approx(-6.87935580446044));

        double l_included = t_included.compute_log_likelihood(2);

        REQUIRE(l_included == Approx(-55.01646493341547));

        REQUIRE(t.get_likelihood_correction() == t_included.get_likelihood_correction());
    }
}


TEST_CASE("phylonet; Testing coalesce_in_branch for 2 lineages and theta of 1.0",
        "[BasePopulationNetwork]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = BasePopulationNetwork::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("phylonet; Testing coalesce_in_branch for 2 lineages and theta of 3.7",
        "[BasePopulationNetwork]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 3.7;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = BasePopulationNetwork::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("phylonet; Testing coalesce_in_branch for 2 lineages and theta of 0.17",
        "[BasePopulationNetwork]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 2;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = BasePopulationNetwork::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("phylonet; Testing coalesce_in_branch for 3 lineages and theta of 1.0",
        "[BasePopulationNetwork]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = BasePopulationNetwork::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("phylonet; Testing coalesce_in_branch for 3 lineages and theta of 1.47",
        "[BasePopulationNetwork]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 1.47;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = BasePopulationNetwork::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("phylonet; Testing coalesce_in_branch for 3 lineages and theta of 0.17",
        "[BasePopulationNetwork]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 3;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = BasePopulationNetwork::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("phylonet; Testing coalesce_in_branch for 10 lineages and theta of 1.0",
        "[BasePopulationNetwork]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 1.0;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = BasePopulationNetwork::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("phylonet; Testing coalesce_in_branch for 10 lineages and theta of 1.47",
        "[BasePopulationNetwork]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 1.47;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = BasePopulationNetwork::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.005));
        REQUIRE(mn >= 0.0);
    }
}

TEST_CASE("phylonet; Testing coalesce_in_branch for 10 lineages and theta of 0.17",
        "[BasePopulationNetwork]") {

    SECTION("Testing coalesce_in_branch") {

        unsigned int nlineages = 10;
        double theta = 0.17;
        double population_size = theta;

        std::vector<std::shared_ptr<GeneTreeSimNode> > lineages;
        lineages.reserve(nlineages);

        double expected_mean = theta * (1.0 - (1.0 / nlineages));
        double expected_variance = expected_mean * expected_mean;

        RandomNumberGenerator rng = RandomNumberGenerator(5431);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            lineages.clear();
            for (unsigned int k = 0; k < nlineages; ++k) {
                lineages.push_back(std::make_shared<GeneTreeSimNode>(0.0));
            }
            double x = BasePopulationNetwork::coalesce_in_branch(
                    lineages,
                    population_size,
                    rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);

        REQUIRE(lineages.size() == 1);
        REQUIRE(mean == Approx(expected_mean).epsilon(0.0005));
        REQUIRE(mn >= 0.0);
    }
}


TEST_CASE("phylonet; Testing scaling of simulate_gene_tree for pair",
        "[BasePopulationNetwork]") {

    SECTION("Testing pair") {
        unsigned int Ne_root = 100000;
        unsigned int Ne_0 = 200000;
        unsigned int Ne_1 = 10000;
        double mu = 1e-8;
        double theta_root = 4 * Ne_root * mu;
        double theta_0 = 4 * Ne_0 * mu;
        double theta_1 = 4 * Ne_1 * mu;
        double time = 2e7;

        double expected_mean_root = theta_root * (1.0 - (1.0 / 2));
        double expected_variance_root = expected_mean_root * expected_mean_root;
        expected_mean_root += (time * mu);

        double expected_mean_0 = theta_0 * (1.0 - (1.0 / 10));
        double expected_mean_1 = theta_1 * (1.0 - (1.0 / 10));

        std::string nex_path = "data/hemi129-5-5.nex";
        BasePopulationNetwork tree(nex_path, ' ', true, true, false);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne_root);
        tree.get_root_ptr()->get_child(0)->set_population_size(Ne_0);
        tree.get_root_ptr()->get_child(1)->set_population_size(Ne_1);

        tree.set_root_height(time);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(mu);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        std::shared_ptr<GeneTreeSimNode> gtree;
        SampleSummarizer<double> height_root;
        SampleSummarizer<double> height_0;
        SampleSummarizer<double> height_1;
        SampleSummarizer<double> tree_length;

        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            REQUIRE(gtree->is_root());
            REQUIRE(gtree->get_number_of_children() == 2);
            REQUIRE(gtree->get_leaf_node_count() == 20);
            height_root.add_sample(gtree->get_height());
            double h0 = gtree->get_child(0)->get_height();
            double h1 = gtree->get_child(1)->get_height();
            if (h0 < h1) {
                height_0.add_sample(h1);
                height_1.add_sample(h0);
            }
            else {
                height_0.add_sample(h0);
                height_1.add_sample(h1);
            }
            tree_length.add_sample(gtree->get_clade_length());
        }

        REQUIRE(height_root.mean() == Approx(expected_mean_root).epsilon(0.0005));
        REQUIRE(height_root.variance() == Approx(expected_variance_root).epsilon(0.0005));

        REQUIRE(height_0.mean() == Approx(expected_mean_0).epsilon(0.00003));
        REQUIRE(height_1.mean() == Approx(expected_mean_1).epsilon(0.00003));


        tree.set_mutation_rate(mu * 2.0);

        SampleSummarizer<double> scale_height_root;
        SampleSummarizer<double> scale_height_0;
        SampleSummarizer<double> scale_height_1;
        SampleSummarizer<double> scale_tree_length;

        for (unsigned int i = 0; i < 100000; ++i) {
            gtree = tree.simulate_gene_tree(0, rng);
            gtree->scale(0.5);
            REQUIRE(gtree->is_root());
            REQUIRE(gtree->get_number_of_children() == 2);
            REQUIRE(gtree->get_leaf_node_count() == 20);
            scale_height_root.add_sample(gtree->get_height());
            double h0 = gtree->get_child(0)->get_height();
            double h1 = gtree->get_child(1)->get_height();
            if (h0 < h1) {
                scale_height_0.add_sample(h1);
                scale_height_1.add_sample(h0);
            }
            else {
                scale_height_0.add_sample(h0);
                scale_height_1.add_sample(h1);
            }
            scale_tree_length.add_sample(gtree->get_clade_length());
        }

        REQUIRE(scale_height_root.mean() == Approx(expected_mean_root).epsilon(0.0005));
        REQUIRE(scale_height_root.variance() == Approx(expected_variance_root).epsilon(0.0005));

        REQUIRE(scale_height_0.mean() == Approx(expected_mean_0).epsilon(0.00003));
        REQUIRE(scale_height_1.mean() == Approx(expected_mean_1).epsilon(0.00003));

        REQUIRE(scale_tree_length.mean() == Approx(tree_length.mean()).epsilon(0.0001));
        REQUIRE(scale_tree_length.variance() == Approx(tree_length.variance()).epsilon(0.0001));
    }
}

TEST_CASE("phylonet; Testing dataset simulation", "[BasePopulationNetwork]") {
    SECTION("Testing simulate_biallelic_data_set for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        BasePopulationNetwork tree(nex_path, ' ', true, false, false, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_root_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        BiallelicData data = tree.simulate_biallelic_data_set(rng);

        REQUIRE(data.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data.markers_are_dominant() == false);
        REQUIRE(tree.get_data().markers_are_dominant() == false);

        data.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8).epsilon(0.01));
        REQUIRE(data.get_proportion_1() == Approx(0.625).epsilon(0.01));

        std::string io_nex_path = "data/tmp-data-test1.nex";
        std::ofstream out;
        out.open(io_nex_path);
        data.write_nexus(out, '-');
        out.close();
        REQUIRE(path::exists(io_nex_path));
        BiallelicData io_data(io_nex_path,
                '-',   // pop name delimiter
                true,  // pop name is prefix
                false, // genotypes are diploid
                false, // markers are dominant
                true,  // validate
                false   // store seq loci info
                );
        REQUIRE(io_data.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(io_data.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(io_data.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(io_data.markers_are_dominant() == false);

        io_data.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8).epsilon(0.01));
        REQUIRE(io_data.get_proportion_1() == Approx(0.625).epsilon(0.01));
    }
}

TEST_CASE("phylonet; Testing complete linked dataset simulation", "[BasePopulationNetwork]") {
    SECTION("Testing simulate_complete_biallelic_data_set for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        BasePopulationNetwork tree(nex_path, ' ', true, false, false, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_root_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        unsigned int locus_size = 100;
        std::pair<BiallelicData, unsigned int> data_nloci= tree.simulate_complete_biallelic_data_set(
                rng,
                locus_size,
                1.0,
                true);
        BiallelicData data = data_nloci.first;
        unsigned int nloci = data_nloci.second;
        std::vector<unsigned int> expected_locus_end_indices;
        unsigned int end_idx = (locus_size - 1);
        while(end_idx < tree.get_data().get_number_of_sites()) {
            expected_locus_end_indices.push_back(end_idx);
            end_idx += locus_size;
        }
        expected_locus_end_indices.push_back(tree.get_data().get_number_of_sites() - 1);
        
        REQUIRE(data.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data.markers_are_dominant() == false);
        REQUIRE(tree.get_data().markers_are_dominant() == false);
        REQUIRE(data.get_locus_end_indices() == expected_locus_end_indices);
        REQUIRE(nloci == data.get_locus_end_indices().size());

        data.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8).epsilon(0.01));
        REQUIRE(data.get_proportion_1() == Approx(0.625).epsilon(0.01));

        std::string io_nex_path = "data/tmp-data-test-complete1.nex";
        std::ofstream out;
        out.open(io_nex_path);
        data.write_nexus(out, '-');
        out.close();
        REQUIRE(path::exists(io_nex_path));
        BiallelicData io_data(io_nex_path,
                '-',   // pop name delimiter
                true,  // pop name is prefix
                false, // genotypes are diploid
                false, // markers are dominant
                true,  // validate
                true   // store seq loci info
                );
        REQUIRE(io_data.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(io_data.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(io_data.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(io_data.markers_are_dominant() == false);
        REQUIRE(io_data.get_locus_end_indices() == expected_locus_end_indices);

        io_data.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8).epsilon(0.01));
        REQUIRE(io_data.get_proportion_1() == Approx(0.625).epsilon(0.01));
    }
}

TEST_CASE("phylonet; Testing complete linked dataset simulation one locus", "[BasePopulationNetwork]") {
    SECTION("Testing simulate_complete_biallelic_data_set for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        BasePopulationNetwork tree(nex_path, ' ', true, false, false, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_root_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        unsigned int locus_size = 10000;
        std::pair<BiallelicData, unsigned int> data_nloci= tree.simulate_complete_biallelic_data_set(
                rng,
                locus_size,
                1.0,
                true);
        BiallelicData data = data_nloci.first;
        unsigned int nloci = data_nloci.second;
        REQUIRE(nloci == 1);
        std::vector<unsigned int> expected_locus_end_indices;
        expected_locus_end_indices.push_back(tree.get_data().get_number_of_sites() - 1);
        
        REQUIRE(data.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data.markers_are_dominant() == false);
        REQUIRE(tree.get_data().markers_are_dominant() == false);
        REQUIRE(data.get_locus_end_indices() == expected_locus_end_indices);

        data.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8).epsilon(0.01));
        REQUIRE(data.get_proportion_1() == Approx(0.625).epsilon(0.01));

        std::string io_nex_path = "data/tmp-data-test-complete1.nex";
        std::ofstream out;
        out.open(io_nex_path);
        data.write_nexus(out, '-');
        out.close();
        REQUIRE(path::exists(io_nex_path));
        BiallelicData io_data(io_nex_path,
                '-',   // pop name delimiter
                true,  // pop name is prefix
                false, // genotypes are diploid
                false, // markers are dominant
                true,  // validate
                true   // store seq loci info
                );
        REQUIRE(io_data.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(io_data.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(io_data.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(io_data.markers_are_dominant() == false);
        REQUIRE(io_data.get_locus_end_indices() == expected_locus_end_indices);

        io_data.get_empirical_u_v_rates(u, v);
        REQUIRE(u == Approx(0.8).epsilon(0.01));
        REQUIRE(io_data.get_proportion_1() == Approx(0.625).epsilon(0.01));
    }
}

TEST_CASE("phylonet; Testing linked dataset simulation", "[BasePopulationNetwork]") {
    SECTION("Testing simulate_linked_biallelic_data_set for fully fixed pair") {
        std::string nex_path = "data/hemi129-with-missing.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        BasePopulationNetwork tree(nex_path, ' ',
                true,   // population_name_is_prefix
                true,   // diploid
                false,  // dominant
                false,  // constant sites removed
                true,   // validate
                false,  // strict on constant
                false,  // strict on missing
                false,  // strict on triallelic
                2.0,    // ploidy
                true    // store charset info
                );

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_root_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        double u_sum = 0.0;
        double prop_sum = 0.0;
        double io_u_sum = 0.0;
        double io_prop_sum = 0.0;
        unsigned int nreps = 100;
        BiallelicData data;
        BiallelicData io_data;
        for (unsigned int i = 0; i < nreps; ++i) {
            data = tree.simulate_linked_biallelic_data_set(rng,
                    1.0,    // singleton sample prob
                    false,  // max one var site per locus
                    true    // validate
                    );

            REQUIRE(data.get_number_of_sites() == tree.get_data().get_number_of_sites());
            REQUIRE(data.get_number_of_populations() == tree.get_data().get_number_of_populations());
            REQUIRE(data.get_population_labels() == tree.get_data().get_population_labels());
            REQUIRE(data.markers_are_dominant() == false);
            REQUIRE(tree.get_data().markers_are_dominant() == false);

            REQUIRE(tree.get_data().get_locus_end_indices() == data.get_locus_end_indices());
            REQUIRE(tree.get_data().has_seq_loci_info() == data.has_seq_loci_info());

            data.get_empirical_u_v_rates(u, v);
            u_sum += u;
            prop_sum += data.get_proportion_1();

            std::string io_nex_path = "data/tmp-data-test2-" + std::to_string(i) + ".nex";
            std::ofstream out;
            out.open(io_nex_path);
            data.write_nexus(out, '-');
            out.close();
            REQUIRE(path::exists(io_nex_path));
            io_data = BiallelicData(io_nex_path,
                    '-',   // pop name delimiter
                    true,  // pop name is prefix
                    false, // genotypes are diploid
                    false, // markers are dominant
                    true,  // validate
                    true   // store seq loci info
                    );
            REQUIRE(io_data.get_number_of_sites() == tree.get_data().get_number_of_sites());
            REQUIRE(io_data.get_number_of_populations() == tree.get_data().get_number_of_populations());
            REQUIRE(io_data.get_population_labels() == tree.get_data().get_population_labels());
            REQUIRE(io_data.markers_are_dominant() == false);

            REQUIRE(tree.get_data().get_locus_end_indices() == io_data.get_locus_end_indices());
            REQUIRE(tree.get_data().has_seq_loci_info() == io_data.has_seq_loci_info());

            io_data.get_empirical_u_v_rates(u, v);
            io_u_sum += u;
            io_prop_sum += data.get_proportion_1();
        }
        REQUIRE(u_sum / nreps == Approx(0.8).epsilon(0.01));
        REQUIRE(prop_sum / nreps == Approx(0.625).epsilon(0.01));
        REQUIRE(io_u_sum == Approx(u_sum));
        REQUIRE(io_prop_sum == Approx(prop_sum));

        REQUIRE(data.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {3, 8, 13, 18};
        REQUIRE(data.get_contiguous_pattern_indices().size() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_locus_end_indices() == expected_locus_ends);

        REQUIRE(io_data.has_seq_loci_info() == true);
        REQUIRE(io_data.get_contiguous_pattern_indices().size() == tree.get_data().get_number_of_sites());
        REQUIRE(io_data.get_locus_end_indices() == expected_locus_ends);
    }
}

TEST_CASE("phylonet; Testing linked dataset simulation with aflp dataset", "[BasePopulationNetwork]") {
    SECTION("Testing simulate_linked_biallelic_data_set for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        BasePopulationNetwork tree(nex_path, ' ',
                true,   // population_name_is_prefix
                false,  // diploid
                false,  // dominant
                false,  // constant sites removed
                true,   // validate
                false,  // strict on constant
                false,  // strict on missing
                false,  // strict on triallelic
                2.0,    // ploidy
                true    // store charset info
                );

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_root_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        double u_sum = 0.0;
        double prop_sum = 0.0;
        double io_u_sum = 0.0;
        double io_prop_sum = 0.0;
        unsigned int nreps = 10;
        BiallelicData data;
        BiallelicData io_data;
        for (unsigned int i = 0; i < nreps; ++i) {
            data = tree.simulate_linked_biallelic_data_set(rng,
                    1.0,    // singleton sample prob
                    false,  // max one var site per locus
                    true    // validate
                    );

            REQUIRE(data.get_number_of_sites() == tree.get_data().get_number_of_sites());
            REQUIRE(data.get_number_of_populations() == tree.get_data().get_number_of_populations());
            REQUIRE(data.get_population_labels() == tree.get_data().get_population_labels());
            REQUIRE(data.markers_are_dominant() == false);
            REQUIRE(tree.get_data().markers_are_dominant() == false);

            REQUIRE(tree.get_data().get_locus_end_indices() == data.get_locus_end_indices());
            REQUIRE(tree.get_data().has_seq_loci_info() == data.has_seq_loci_info());

            data.get_empirical_u_v_rates(u, v);
            u_sum += u;
            prop_sum += data.get_proportion_1();

            std::string io_nex_path = "data/tmp-data-test3-" + std::to_string(i) + ".nex";
            std::ofstream out;
            out.open(io_nex_path);
            data.write_nexus(out, '-');
            out.close();
            REQUIRE(path::exists(io_nex_path));
            io_data = BiallelicData(io_nex_path,
                    '-',   // pop name delimiter
                    true,  // pop name is prefix
                    false, // genotypes are diploid
                    false, // markers are dominant
                    true,  // validate
                    true   // store seq loci info
                    );
            REQUIRE(io_data.get_number_of_sites() == tree.get_data().get_number_of_sites());
            REQUIRE(io_data.get_number_of_populations() == tree.get_data().get_number_of_populations());
            REQUIRE(io_data.get_population_labels() == tree.get_data().get_population_labels());
            REQUIRE(io_data.markers_are_dominant() == false);

            REQUIRE(tree.get_data().get_locus_end_indices() == io_data.get_locus_end_indices());
            REQUIRE(tree.get_data().has_seq_loci_info() == io_data.has_seq_loci_info());

            io_data.get_empirical_u_v_rates(u, v);
            io_u_sum += u;
            io_prop_sum += data.get_proportion_1();
        }
        REQUIRE(u_sum / nreps == Approx(0.8).epsilon(0.01));
        REQUIRE(prop_sum / nreps == Approx(0.625).epsilon(0.01));
        REQUIRE(io_u_sum == Approx(u_sum));
        REQUIRE(io_prop_sum == Approx(prop_sum));

        REQUIRE(data.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends = {99, 199, 299, 399, 499, 599, 699, 799, 899, 999, 1099, 1199, 1216};
        REQUIRE(data.get_contiguous_pattern_indices().size() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_locus_end_indices() == expected_locus_ends);

        REQUIRE(io_data.has_seq_loci_info() == true);
        REQUIRE(io_data.get_contiguous_pattern_indices().size() == tree.get_data().get_number_of_sites());
        REQUIRE(io_data.get_locus_end_indices() == expected_locus_ends);
    }
}

TEST_CASE("phylonet; Testing singleton acquisition bias", "[BasePopulationNetwork]") {
    SECTION("Testing for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        BasePopulationNetwork tree(nex_path, ' ', true, false, false, false);

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_root_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        BiallelicData data1 = tree.simulate_biallelic_data_set(rng, 1.0);
        BiallelicData data0 = tree.simulate_biallelic_data_set(rng, 0.0);

        REQUIRE(data1.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data1.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data1.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data1.markers_are_dominant() == false);
        REQUIRE(data0.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data0.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data0.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data0.markers_are_dominant() == false);
        REQUIRE(tree.get_data().markers_are_dominant() == false);

        unsigned int singleton_count10 = 0;
        unsigned int singleton_count05 = 0;
        unsigned int singleton_count00 = 0;
        RandomNumberGenerator rng10 = RandomNumberGenerator(123);
        RandomNumberGenerator rng05 = RandomNumberGenerator(123);
        RandomNumberGenerator rng00 = RandomNumberGenerator(123);
        for (unsigned int rep = 0; rep < 100; ++rep) {
            BiallelicData data00 = tree.simulate_biallelic_data_set(rng00, 0.0);
            BiallelicData data05 = tree.simulate_biallelic_data_set(rng05, 0.5);
            BiallelicData data10 = tree.simulate_biallelic_data_set(rng10, 1.0);
            for (unsigned int i = 0; i < data00.get_number_of_patterns(); ++i) {
                const std::vector<unsigned int>& red_allele_counts = data00.get_red_allele_counts(i);
                const std::vector<unsigned int>& allele_counts = data00.get_allele_counts(i);
                unsigned int nreds = 0;
                unsigned int nalleles = 0;
                for (unsigned int j = 0; j < allele_counts.size(); ++j) {
                    nreds += red_allele_counts.at(j);
                    nalleles += allele_counts.at(j);
                }
                if ((nreds == 1) || (nreds == (nalleles - 1))) {
                    singleton_count00 += data00.get_pattern_weight(i);
                }
            }
            for (unsigned int i = 0; i < data05.get_number_of_patterns(); ++i) {
                const std::vector<unsigned int>& red_allele_counts = data05.get_red_allele_counts(i);
                const std::vector<unsigned int>& allele_counts = data05.get_allele_counts(i);
                unsigned int nreds = 0;
                unsigned int nalleles = 0;
                for (unsigned int j = 0; j < allele_counts.size(); ++j) {
                    nreds += red_allele_counts.at(j);
                    nalleles += allele_counts.at(j);
                }
                if ((nreds == 1) || (nreds == (nalleles - 1))) {
                    singleton_count05 += data05.get_pattern_weight(i);
                }
            }
            for (unsigned int i = 0; i < data10.get_number_of_patterns(); ++i) {
                const std::vector<unsigned int>& red_allele_counts = data10.get_red_allele_counts(i);
                const std::vector<unsigned int>& allele_counts = data10.get_allele_counts(i);
                unsigned int nreds = 0;
                unsigned int nalleles = 0;
                for (unsigned int j = 0; j < allele_counts.size(); ++j) {
                    nreds += red_allele_counts.at(j);
                    nalleles += allele_counts.at(j);
                }
                if ((nreds == 1) || (nreds == (nalleles - 1))) {
                    singleton_count10 += data10.get_pattern_weight(i);
                }
            }
        }
        REQUIRE(singleton_count10 > 0);
        REQUIRE(singleton_count05 > 0);
        REQUIRE(singleton_count00 == 0);
        REQUIRE((double)singleton_count05 == Approx(singleton_count10 * 0.5).epsilon(0.05));
    }
}

TEST_CASE("phylonet; Testing singleton acquisition bias with charsets", "[BasePopulationNetwork]") {
    SECTION("Testing for fully fixed pair") {
        std::string nex_path = "data/aflp_25.nex";
        // Need to keep constant characters to get expected
        // character state frequencies
        BasePopulationNetwork tree(nex_path, ' ',
                true,   // population_name_is_prefix
                false,  // diploid
                false,  // dominant
                false,  // constant sites removed
                true,   // validate
                false,  // strict on constant
                false,  // strict on missing
                false,  // strict on triallelic
                2.0,    // ploidy
                true    // store charset info
                );

        std::shared_ptr<ExponentialDistribution> height_prior = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<GammaDistribution> size_prior = std::make_shared<GammaDistribution>(2.0, 1.2);
        std::shared_ptr<BetaDistribution> f_prior = std::make_shared<BetaDistribution>(2.0, 1.5);
        std::shared_ptr<GammaDistribution> rate_prior = std::make_shared<GammaDistribution>(3.0, 1.1);

        tree.set_root_node_height_prior(height_prior);
        tree.set_population_size_prior(size_prior);
        tree.set_freq_1_prior(f_prior);
        tree.set_mutation_rate_prior(rate_prior);

        tree.set_root_height(0.1);
        tree.constrain_population_sizes();
        tree.set_all_population_sizes(0.005);
        tree.fix_population_sizes();
        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);
        tree.fix_mutation_rate();

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.625);
        tree.fix_state_frequencies();

        double u;
        double v;

        tree.get_data().get_empirical_u_v_rates(u, v);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        BiallelicData data1 = tree.simulate_linked_biallelic_data_set(rng, 1.0, false, true);
        BiallelicData data0 = tree.simulate_linked_biallelic_data_set(rng, 0.0, false, true);

        BiallelicData io_data1;
        BiallelicData io_data0;

        std::string io_nex_path1 = "data/tmp-data-test4.nex";
        std::string io_nex_path0 = "data/tmp-data-test5.nex";
        std::ofstream out1;
        std::ofstream out0;
        out1.open(io_nex_path1);
        out0.open(io_nex_path0);
        data1.write_nexus(out1, '-');
        data0.write_nexus(out0, '-');
        out1.close();
        out0.close();
        REQUIRE(path::exists(io_nex_path1));
        REQUIRE(path::exists(io_nex_path0));
        io_data1 = BiallelicData(io_nex_path1,
                '-',   // pop name delimiter
                true,  // pop name is prefix
                false, // genotypes are diploid
                false, // markers are dominant
                true,  // validate
                true   // store seq loci info
                );
        io_data0 = BiallelicData(io_nex_path0,
                '-',   // pop name delimiter
                true,  // pop name is prefix
                false, // genotypes are diploid
                false, // markers are dominant
                true,  // validate
                true   // store seq loci info
                );

        REQUIRE(data1.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data1.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data1.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data1.markers_are_dominant() == false);
        REQUIRE(data1.get_locus_end_indices() == tree.get_data().get_locus_end_indices());
        REQUIRE(data0.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(data0.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(data0.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(data0.markers_are_dominant() == false);
        REQUIRE(data0.get_locus_end_indices() == tree.get_data().get_locus_end_indices());
        REQUIRE(tree.get_data().markers_are_dominant() == false);

        REQUIRE(io_data1.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(io_data1.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(io_data1.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(io_data1.markers_are_dominant() == false);
        REQUIRE(io_data1.get_locus_end_indices() == tree.get_data().get_locus_end_indices());
        REQUIRE(io_data0.get_number_of_sites() == tree.get_data().get_number_of_sites());
        REQUIRE(io_data0.get_number_of_populations() == tree.get_data().get_number_of_populations());
        REQUIRE(io_data0.get_population_labels() == tree.get_data().get_population_labels());
        REQUIRE(io_data0.markers_are_dominant() == false);
        REQUIRE(io_data0.get_locus_end_indices() == tree.get_data().get_locus_end_indices());

        unsigned int singleton_count10 = 0;
        unsigned int singleton_count05 = 0;
        unsigned int singleton_count00 = 0;
        unsigned int io_singleton_count10 = 0;
        unsigned int io_singleton_count05 = 0;
        unsigned int io_singleton_count00 = 0;
        RandomNumberGenerator rng10 = RandomNumberGenerator(123);
        RandomNumberGenerator rng05 = RandomNumberGenerator(123);
        RandomNumberGenerator rng00 = RandomNumberGenerator(123);
        for (unsigned int rep = 0; rep < 100; ++rep) {
            BiallelicData data00 = tree.simulate_linked_biallelic_data_set(rng00, 0.0, false, true);
            BiallelicData data05 = tree.simulate_linked_biallelic_data_set(rng05, 0.5, false, true);
            BiallelicData data10 = tree.simulate_linked_biallelic_data_set(rng10, 1.0, false, true);
            BiallelicData io_data00;
            BiallelicData io_data05;
            BiallelicData io_data10;

            std::string io_nex_path00 = "data/tmp-data-test6-" + std::to_string(rep) + ".nex";
            std::string io_nex_path05 = "data/tmp-data-test7-" + std::to_string(rep) + ".nex";
            std::string io_nex_path10 = "data/tmp-data-test8-" + std::to_string(rep) + ".nex";
            std::ofstream out00;
            std::ofstream out05;
            std::ofstream out10;
            out00.open(io_nex_path00);
            out05.open(io_nex_path05);
            out10.open(io_nex_path10);
            data00.write_nexus(out00, '-');
            data05.write_nexus(out05, '-');
            data10.write_nexus(out10, '-');
            out00.close();
            out05.close();
            out10.close();
            REQUIRE(path::exists(io_nex_path00));
            REQUIRE(path::exists(io_nex_path05));
            REQUIRE(path::exists(io_nex_path10));
            io_data00 = BiallelicData(io_nex_path00,
                    '-',   // pop name delimiter
                    true,  // pop name is prefix
                    false, // genotypes are diploid
                    false, // markers are dominant
                    true,  // validate
                    true   // store seq loci info
                    );
            io_data05 = BiallelicData(io_nex_path05,
                    '-',   // pop name delimiter
                    true,  // pop name is prefix
                    false, // genotypes are diploid
                    false, // markers are dominant
                    true,  // validate
                    true   // store seq loci info
                    );
            io_data10 = BiallelicData(io_nex_path10,
                    '-',   // pop name delimiter
                    true,  // pop name is prefix
                    false, // genotypes are diploid
                    false, // markers are dominant
                    true,  // validate
                    true   // store seq loci info
                    );
            REQUIRE(data00.get_number_of_patterns() == io_data00.get_number_of_patterns());
            REQUIRE(data05.get_number_of_patterns() == io_data05.get_number_of_patterns());
            REQUIRE(data10.get_number_of_patterns() == io_data10.get_number_of_patterns());

            for (unsigned int i = 0; i < data00.get_number_of_patterns(); ++i) {
                const std::vector<unsigned int>& red_allele_counts = data00.get_red_allele_counts(i);
                const std::vector<unsigned int>& allele_counts = data00.get_allele_counts(i);
                unsigned int nreds = 0;
                unsigned int nalleles = 0;
                for (unsigned int j = 0; j < allele_counts.size(); ++j) {
                    nreds += red_allele_counts.at(j);
                    nalleles += allele_counts.at(j);
                }
                if ((nreds == 1) || (nreds == (nalleles - 1))) {
                    singleton_count00 += data00.get_pattern_weight(i);
                }

                const std::vector<unsigned int>& io_red_allele_counts = io_data00.get_red_allele_counts(i);
                const std::vector<unsigned int>& io_allele_counts = io_data00.get_allele_counts(i);
                unsigned int io_nreds = 0;
                unsigned int io_nalleles = 0;
                for (unsigned int j = 0; j < io_allele_counts.size(); ++j) {
                    io_nreds += io_red_allele_counts.at(j);
                    io_nalleles += io_allele_counts.at(j);
                }
                if ((io_nreds == 1) || (io_nreds == (io_nalleles - 1))) {
                    io_singleton_count00 += io_data00.get_pattern_weight(i);
                }
            }
            for (unsigned int i = 0; i < data05.get_number_of_patterns(); ++i) {
                const std::vector<unsigned int>& red_allele_counts = data05.get_red_allele_counts(i);
                const std::vector<unsigned int>& allele_counts = data05.get_allele_counts(i);
                unsigned int nreds = 0;
                unsigned int nalleles = 0;
                for (unsigned int j = 0; j < allele_counts.size(); ++j) {
                    nreds += red_allele_counts.at(j);
                    nalleles += allele_counts.at(j);
                }
                if ((nreds == 1) || (nreds == (nalleles - 1))) {
                    singleton_count05 += data05.get_pattern_weight(i);
                }

                const std::vector<unsigned int>& io_red_allele_counts = io_data05.get_red_allele_counts(i);
                const std::vector<unsigned int>& io_allele_counts = io_data05.get_allele_counts(i);
                unsigned int io_nreds = 0;
                unsigned int io_nalleles = 0;
                for (unsigned int j = 0; j < io_allele_counts.size(); ++j) {
                    io_nreds += io_red_allele_counts.at(j);
                    io_nalleles += io_allele_counts.at(j);
                }
                if ((io_nreds == 1) || (io_nreds == (io_nalleles - 1))) {
                    io_singleton_count05 += io_data05.get_pattern_weight(i);
                }
            }
            for (unsigned int i = 0; i < data10.get_number_of_patterns(); ++i) {
                const std::vector<unsigned int>& red_allele_counts = data10.get_red_allele_counts(i);
                const std::vector<unsigned int>& allele_counts = data10.get_allele_counts(i);
                unsigned int nreds = 0;
                unsigned int nalleles = 0;
                for (unsigned int j = 0; j < allele_counts.size(); ++j) {
                    nreds += red_allele_counts.at(j);
                    nalleles += allele_counts.at(j);
                }
                if ((nreds == 1) || (nreds == (nalleles - 1))) {
                    singleton_count10 += data10.get_pattern_weight(i);
                }

                const std::vector<unsigned int>& io_red_allele_counts = io_data10.get_red_allele_counts(i);
                const std::vector<unsigned int>& io_allele_counts = io_data10.get_allele_counts(i);
                unsigned int io_nreds = 0;
                unsigned int io_nalleles = 0;
                for (unsigned int j = 0; j < io_allele_counts.size(); ++j) {
                    io_nreds += io_red_allele_counts.at(j);
                    io_nalleles += io_allele_counts.at(j);
                }
                if ((io_nreds == 1) || (io_nreds == (io_nalleles - 1))) {
                    io_singleton_count10 += io_data10.get_pattern_weight(i);
                }
            }
        }
        REQUIRE(singleton_count10 > 0);
        REQUIRE(singleton_count05 > 0);
        REQUIRE(singleton_count00 == 0);
        REQUIRE((double)singleton_count05 == Approx(singleton_count10 * 0.5).epsilon(0.05));

        REQUIRE(io_singleton_count10 > 0);
        REQUIRE(io_singleton_count05 > 0);
        REQUIRE(io_singleton_count00 == 0);
        REQUIRE((double)io_singleton_count05 == Approx(io_singleton_count10 * 0.5).epsilon(0.05));
    }
}

TEST_CASE("phylonet; Testing scaling of dataset simulation for singleton",
        "[BasePopulationNetwork]") {

    SECTION("Testing simulate_biallelic_data_set for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        // Two times the expected gene tree root height
        double expected_mean = 2.0 * (theta * (1.0 - (1.0 / nlineages)));

        std::string nex_path = "data/dummy-singleton-n2-100k.nex";
        // Need to keep constant characters
        BasePopulationNetwork tree(nex_path, '-', true, false, false, false);
        REQUIRE(tree.get_leaf_node_count() == 1);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.get_root_ptr()->get_child(0)->set_population_size(Ne * mu);

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        BiallelicData data;
        SampleSummarizer<double> divergence;
        unsigned int nsamples = 100;

        for (unsigned int i = 0; i < nsamples; ++i) {
            data = tree.simulate_biallelic_data_set(rng, 1.0, false);
            REQUIRE(data.get_number_of_sites() == 100000);
            double x = (double)data.get_number_of_variable_sites() / (double)data.get_number_of_sites();
            divergence.add_sample(x);
        }

        std::cout << "Expected mean divergence: " << expected_mean << "\n";
        std::cout << "Simulated mean divergence: " << divergence.mean() << "\n";
        REQUIRE(divergence.mean() == Approx(expected_mean).epsilon(0.0001));

        tree.set_mutation_rate(mu);
        tree.set_root_population_size(Ne);
        tree.get_root_ptr()->get_child(0)->set_population_size(Ne);
        tree.set_root_height(0.1 / mu);

        SampleSummarizer<double> divergence2;

        for (unsigned int i = 0; i < nsamples; ++i) {
            data = tree.simulate_biallelic_data_set(rng, 1.0, false);
            REQUIRE(data.get_number_of_sites() == 100000);
            double x = (double)data.get_number_of_variable_sites() / (double)data.get_number_of_sites();
            divergence2.add_sample(x);
        }

        std::cout << "Simulated mean divergence: " << divergence2.mean() << "\n";
        REQUIRE(divergence2.mean() == Approx(expected_mean).epsilon(0.0001));
    }
}

TEST_CASE("phylonet; Testing scaling of dataset simulation for singleton with charsets",
        "[BasePopulationNetwork]") {

    SECTION("Testing simulate_biallelic_data_set for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        // Two times the expected gene tree root height
        double expected_mean = 2.0 * (theta * (1.0 - (1.0 / nlineages)));

        std::string nex_path = "data/dummy-singleton-n2-100k.nex";
        // Need to keep constant characters
        BasePopulationNetwork tree(nex_path, '-',
                true,   // population_name_is_prefix
                false,  // diploid
                false,  // dominant
                false,  // constant sites removed
                true,   // validate
                false,  // strict on constant
                false,  // strict on missing
                false,  // strict on triallelic
                2.0,    // ploidy
                true    // store charset info
                );
        REQUIRE(tree.get_leaf_node_count() == 1);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.get_root_ptr()->get_child(0)->set_population_size(Ne * mu);

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(12345678);

        BiallelicData data;
        SampleSummarizer<double> divergence;
        unsigned int nsamples = 100;

        for (unsigned int i = 0; i < nsamples; ++i) {
            data = tree.simulate_linked_biallelic_data_set(rng, 1.0, false);
            REQUIRE(data.get_number_of_sites() == 100000);
            double x = (double)data.get_number_of_variable_sites() / (double)data.get_number_of_sites();
            divergence.add_sample(x);
        }

        std::cout << "Expected mean divergence: " << expected_mean << "\n";
        std::cout << "Simulated mean divergence: " << divergence.mean() << "\n";
        REQUIRE(divergence.mean() == Approx(expected_mean).epsilon(0.0002));

        tree.set_mutation_rate(mu);
        tree.set_root_population_size(Ne);
        tree.get_root_ptr()->get_child(0)->set_population_size(Ne);
        tree.set_root_height(0.1 / mu);

        SampleSummarizer<double> divergence2;

        for (unsigned int i = 0; i < nsamples; ++i) {
            data = tree.simulate_linked_biallelic_data_set(rng, 1.0, false);
            REQUIRE(data.get_number_of_sites() == 100000);
            double x = (double)data.get_number_of_variable_sites() / (double)data.get_number_of_sites();
            divergence2.add_sample(x);
        }

        std::cout << "Simulated mean divergence: " << divergence2.mean() << "\n";
        REQUIRE(divergence2.mean() == Approx(expected_mean).epsilon(0.0002));

        REQUIRE(data.has_seq_loci_info() == true);
        std::vector<unsigned int> expected_locus_ends;
        unsigned int end_idx = 0;
        for (unsigned int i = 0; i < 100; ++i) {
            end_idx += 1000;
            expected_locus_ends.push_back(end_idx - 1);
        }
        REQUIRE(data.get_contiguous_pattern_indices().size() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_locus_end_indices() == expected_locus_ends);
    }
}

TEST_CASE("phylonet; Testing scaling of simulation of loci for singleton",
        "[BasePopulationNetwork]") {

    SECTION("Testing simulate_complete_biallelic_data_set for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        // Two times the expected gene tree root height
        double expected_mean = 2.0 * (theta * (1.0 - (1.0 / nlineages)));

        std::string nex_path = "data/dummy-singleton-n2-100k.nex";
        // Need to keep constant characters
        BasePopulationNetwork tree(nex_path, '-', true, false, false, false);
        REQUIRE(tree.get_leaf_node_count() == 1);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.get_root_ptr()->get_child(0)->set_population_size(Ne * mu);

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(54321);

        SampleSummarizer<double> divergence;
        unsigned int nsamples = 100;

        for (unsigned int i = 0; i < nsamples; ++i) {
            auto data_nloci = tree.simulate_complete_biallelic_data_set(rng, 1000, 1.0, false);
            REQUIRE(data_nloci.second == 100);
            REQUIRE(data_nloci.first.get_number_of_sites() == 100000);
            double x = (double)data_nloci.first.get_number_of_variable_sites() / (double)data_nloci.first.get_number_of_sites();
            divergence.add_sample(x);
        }

        std::cout << "Expected mean divergence: " << expected_mean << "\n";
        std::cout << "Simulated mean divergence: " << divergence.mean() << "\n";
        REQUIRE(divergence.mean() == Approx(expected_mean).epsilon(0.0001));

        tree.set_mutation_rate(mu);
        tree.set_root_population_size(Ne);
        tree.get_root_ptr()->get_child(0)->set_population_size(Ne);
        tree.set_root_height(0.1 / mu);

        SampleSummarizer<double> divergence2;

        for (unsigned int i = 0; i < nsamples; ++i) {
            auto data_nloci = tree.simulate_complete_biallelic_data_set(rng, 1000, 1.0, false);
            REQUIRE(data_nloci.second == 100);
            REQUIRE(data_nloci.first.get_number_of_sites() == 100000);
            double x = (double)data_nloci.first.get_number_of_variable_sites() / (double)data_nloci.first.get_number_of_sites();
            divergence2.add_sample(x);
        }

        std::cout << "Simulated mean divergence: " << divergence2.mean() << "\n";
        REQUIRE(divergence2.mean() == Approx(expected_mean).epsilon(0.0001));
    }
}

TEST_CASE("phylonet; Testing scaling of simulation of one variable site per locus for singleton",
        "[BasePopulationNetwork]") {

    SECTION("Testing simulate_data_set_max_one_variable_site_per_locus for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        // Two times the expected gene tree root height
        double expected_mean = 2.0 * (theta * (1.0 - (1.0 / nlineages)));
        double epsilon = 0.0001;

        std::string nex_path = "data/dummy-singleton-n2-100k.nex";
        // Need to keep constant characters
        BasePopulationNetwork tree(nex_path, '-', true, false, false, false);
        REQUIRE(tree.get_leaf_node_count() == 1);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.get_root_ptr()->get_child(0)->set_population_size(Ne * mu);

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        /* RandomNumberGenerator rng = RandomNumberGenerator(54321); */
        RandomNumberGenerator rng = RandomNumberGenerator(12345678);

        SampleSummarizer<double> divergence;
        unsigned int nsamples = 100;

        for (unsigned int i = 0; i < nsamples; ++i) {
            auto data_nloci = tree.simulate_data_set_max_one_variable_site_per_locus(rng, 1000, 1.0, false);
            REQUIRE(data_nloci.second == 100);
            REQUIRE(data_nloci.first.get_number_of_sites() <= 100);
            double x = (double)data_nloci.first.get_number_of_variable_sites() / (data_nloci.first.get_number_of_sites() + data_nloci.first.get_number_of_constant_sites_removed());
            divergence.add_sample(x);
        }

        std::cout << "Expected mean divergence: " << expected_mean << "\n";
        std::cout << "Simulated mean divergence: " << divergence.mean() << "\n";
        // By stopping the simulation of sites for each locus when we get a
        // variable site, we are under-sampling sites from loci with older
        // coalescence times (over-sampling sites from loci with younger
        // coalescence times), so we should under estimate diversity:
        REQUIRE(divergence.mean() < (expected_mean - epsilon));

        tree.set_mutation_rate(mu);
        tree.set_root_population_size(Ne);
        tree.get_root_ptr()->get_child(0)->set_population_size(Ne);
        tree.set_root_height(0.1 / mu);

        SampleSummarizer<double> divergence2;

        for (unsigned int i = 0; i < nsamples; ++i) {
            auto data_nloci = tree.simulate_data_set_max_one_variable_site_per_locus(rng, 1000, 1.0, false);
            REQUIRE(data_nloci.second == 100);
            REQUIRE(data_nloci.first.get_number_of_sites() <= 100);
            double x = (double)data_nloci.first.get_number_of_variable_sites() / (data_nloci.first.get_number_of_sites() + data_nloci.first.get_number_of_constant_sites_removed());
            divergence2.add_sample(x);
        }

        std::cout << "Simulated mean divergence: " << divergence2.mean() << "\n";
        // By stopping the simulation of sites for each locus when we get a
        // variable site, we are under-sampling sites from loci with older
        // coalescence times (over-sampling sites from loci with younger
        // coalescence times), so we should under estimate diversity:
        REQUIRE(divergence2.mean() < (expected_mean - epsilon));
    }
}

TEST_CASE("phylonet; Testing scaling of simulation of one variable site per locus for singleton with charsets",
        "[BasePopulationNetwork]") {

    SECTION("Testing one SNP per locus for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        // Two times the expected gene tree root height
        double expected_mean = 2.0 * (theta * (1.0 - (1.0 / nlineages)));
        double epsilon = 0.0001;

        std::string nex_path = "data/dummy-singleton-n2-100k.nex";
        // Need to keep constant characters
        BasePopulationNetwork tree(nex_path, '-',
                true,   // population_name_is_prefix
                false,  // diploid
                false,  // dominant
                false,  // constant sites removed
                true,   // validate
                true,   // strict on constant
                true,   // strict on missing
                true,   // strict on triallelic
                2.0,    // ploidy
                true    // store charset info
                );
        REQUIRE(tree.get_leaf_node_count() == 1);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.get_root_ptr()->get_child(0)->set_population_size(Ne * mu);

        tree.set_root_height(0.1);

        tree.set_freq_1(0.5);

        tree.set_mutation_rate(1.0);

        /* RandomNumberGenerator rng = RandomNumberGenerator(54321); */
        RandomNumberGenerator rng = RandomNumberGenerator(12345678);

        SampleSummarizer<double> divergence;
        unsigned int nsamples = 100;

        BiallelicData data;
        for (unsigned int i = 0; i < nsamples; ++i) {
            data = tree.simulate_linked_biallelic_data_set(rng, 1.0, true, true);
            REQUIRE(data.get_number_of_sites() <= 100);
            double x = (double)data.get_number_of_variable_sites() / (data.get_number_of_sites() + data.get_number_of_constant_sites_removed());
            divergence.add_sample(x);
        }

        std::cout << "Expected mean divergence: " << expected_mean << "\n";
        std::cout << "Simulated mean divergence: " << divergence.mean() << "\n";
        // By stopping the simulation of sites for each locus when we get a
        // variable site, we are under-sampling sites from loci with older
        // coalescence times (over-sampling sites from loci with younger
        // coalescence times), so we should under estimate diversity:
        REQUIRE(divergence.mean() < (expected_mean - epsilon));

        tree.set_mutation_rate(mu);
        tree.set_root_population_size(Ne);
        tree.get_root_ptr()->get_child(0)->set_population_size(Ne);
        tree.set_root_height(0.1 / mu);

        SampleSummarizer<double> divergence2;

        for (unsigned int i = 0; i < nsamples; ++i) {
            data = tree.simulate_linked_biallelic_data_set(rng, 1.0, true, true);
            REQUIRE(data.get_number_of_sites() <= 100);
            double x = (double)data.get_number_of_variable_sites() / (data.get_number_of_sites() + data.get_number_of_constant_sites_removed());
            divergence2.add_sample(x);
        }

        std::cout << "Simulated mean divergence: " << divergence2.mean() << "\n";
        // By stopping the simulation of sites for each locus when we get a
        // variable site, we are under-sampling sites from loci with older
        // coalescence times (over-sampling sites from loci with younger
        // coalescence times), so we should under estimate diversity:
        REQUIRE(divergence2.mean() < (expected_mean - epsilon));

        REQUIRE(data.has_seq_loci_info() == false);
        REQUIRE(data.get_contiguous_pattern_indices().size() == 0);
        REQUIRE(data.get_locus_end_indices().size() == 0);


        data = tree.simulate_linked_biallelic_data_set(rng, 1.0, false, true);

        std::vector<unsigned int> expected_locus_ends;
        unsigned int end_idx = 0;
        for (unsigned int i = 0; i < 100; ++i) {
            end_idx += 1000;
            expected_locus_ends.push_back(end_idx - 1);
        }
        REQUIRE(data.get_contiguous_pattern_indices().size() == tree.get_data().get_number_of_sites());
        REQUIRE(data.get_locus_end_indices().size() == tree.get_data().get_locus_end_indices().size());
        REQUIRE(data.get_locus_end_indices() == expected_locus_ends);
    }
}


TEST_CASE("phylonet; Testing errors when trying to sim loci with a template with constant characters removed",
        "[BasePopulationNetwork]") {

    SECTION("Testing for fully fixed singleton") {
        unsigned int Ne = 100000;
        double mu = 1e-8;
        double theta = 4 * Ne * mu;

        unsigned int nlineages = 2;

        // Two times the expected gene tree root height
        double expected_mean = 2.0 * (theta * (1.0 - (1.0 / nlineages)));
        double epsilon = 0.0001;

        std::string nex_path = "data/dummy-singleton-n2-100k.nex";
        // Need to keep constant characters
        BasePopulationNetwork tree(nex_path, '-',
                true,   // population_name_is_prefix
                false,  // diploid
                false,  // dominant
                true,   // constant sites removed
                true,   // validate
                true,   // strict on constant
                true,   // strict on missing
                true,   // strict on triallelic
                2.0,    // ploidy
                true    // store charset info
                );
        REQUIRE(tree.get_leaf_node_count() == 1);
        tree.estimate_mutation_rate();

        tree.set_root_population_size(Ne * mu);
        tree.get_root_ptr()->get_child(0)->set_population_size(Ne * mu);
        tree.set_root_height(0.1);
        tree.set_freq_1(0.5);
        tree.set_mutation_rate(1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(12345678);

        // constant sites were removed, so all of these should fail
        REQUIRE_THROWS_AS(tree.simulate_linked_biallelic_data_set(rng, 1.0, false, false), EcoevolityBiallelicDataError &);
        REQUIRE_THROWS_AS(tree.simulate_linked_biallelic_data_set(rng, 1.0, true, true), EcoevolityBiallelicDataError &);
        REQUIRE_THROWS_AS(tree.simulate_linked_biallelic_data_set(rng, 1.0, true, false), EcoevolityBiallelicDataError &);
        REQUIRE_THROWS_AS(tree.simulate_linked_biallelic_data_set(rng, 1.0, false, true), EcoevolityBiallelicDataError &);

        REQUIRE_THROWS_AS(tree.simulate_complete_biallelic_data_set(rng, 1, 1.0, true), EcoevolityBiallelicDataError &);
        REQUIRE_THROWS_AS(tree.simulate_complete_biallelic_data_set(rng, 1, 1.0, false), EcoevolityBiallelicDataError &);
        REQUIRE_THROWS_AS(tree.simulate_complete_biallelic_data_set(rng, 10, 1.0, true), EcoevolityBiallelicDataError &);
        REQUIRE_THROWS_AS(tree.simulate_complete_biallelic_data_set(rng, 10, 1.0, false), EcoevolityBiallelicDataError &);
    }
}

TEST_CASE("phylonet; Testing BaseTree::get_height_of_youngest_parent()", "[BaseTree]") {
    SECTION("Testing get_height_of_youngest_parent") {
        std::shared_ptr<NetNode> root = std::make_shared<NetNode>("root", 0.1);
        std::shared_ptr<NetNode> internal1 = std::make_shared<NetNode>("internal1", 0.08);
        std::shared_ptr<NetNode> internal2 = std::make_shared<NetNode>("internal2", 0.06);
        std::shared_ptr<NetNode> internal3 = std::make_shared<NetNode>("internal3", 0.04);
        std::shared_ptr<NetNode> internal4 = std::make_shared<NetNode>("internal4", 0.04);
        internal4->set_height_parameter(internal3->get_height_parameter());
        std::shared_ptr<NetNode> internal5 = std::make_shared<NetNode>("internal5", 0.02);
        std::shared_ptr<NetNode> internal6 = std::make_shared<NetNode>("internal6", 0.02);
        internal6->set_height_parameter(internal5->get_height_parameter());
        std::shared_ptr<NetNode> leaf1 = std::make_shared<NetNode>("leaf1", 0.0);
        leaf1->fix_node_height();
        std::shared_ptr<NetNode> leaf2 = std::make_shared<NetNode>("leaf2", 0.0);
        leaf2->fix_node_height();
        std::shared_ptr<NetNode> leaf3 = std::make_shared<NetNode>("leaf3", 0.0);
        leaf3->fix_node_height();
        std::shared_ptr<NetNode> leaf4 = std::make_shared<NetNode>("leaf4", 0.0);
        leaf4->fix_node_height();
        std::shared_ptr<NetNode> leaf5 = std::make_shared<NetNode>("leaf5", 0.0);
        leaf5->fix_node_height();
        std::shared_ptr<NetNode> leaf6 = std::make_shared<NetNode>("leaf6", 0.0);
        leaf6->fix_node_height();
        std::shared_ptr<NetNode> leaf7 = std::make_shared<NetNode>("leaf7", 0.0);
        leaf7->fix_node_height();
        std::shared_ptr<NetNode> leaf8 = std::make_shared<NetNode>("leaf8", 0.0);
        leaf8->fix_node_height();

        internal3->add_child(leaf1);
        internal3->add_child(leaf2);
        internal4->add_child(leaf3);
        internal4->add_child(leaf4);
        internal5->add_child(leaf5);
        internal5->add_child(leaf6);
        internal6->add_child(leaf7);
        internal6->add_child(leaf8);

        internal1->add_child(internal5);
        internal1->add_child(internal3);
        internal2->add_child(internal6);
        internal2->add_child(internal4);

        root->add_child(internal1);
        root->add_child(internal2);
        BaseTree<NetNode> tree(root);

        REQUIRE(tree.get_leaf_node_count() == 8);
        REQUIRE(tree.get_node_count() == 15);
        REQUIRE(tree.get_number_of_node_heights() == 5);
        std::vector<double> expected_heights {0.02, 0.04, 0.06, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        REQUIRE_THROWS_AS(tree.get_youngest_parent(4), EcoevolityError &);
        REQUIRE_THROWS_AS(tree.get_height_of_youngest_parent(4), EcoevolityError &);
        REQUIRE_THROWS_AS(tree.get_index_of_youngest_parent(4), EcoevolityError &);
        REQUIRE(tree.get_youngest_parent(0)->get_height() == 0.06);
        REQUIRE(tree.get_height_of_youngest_parent(0) == 0.06);
        REQUIRE(tree.get_index_of_youngest_parent(0) == 2);
        REQUIRE(tree.get_youngest_parent(1)->get_height() == 0.06);
        REQUIRE(tree.get_height_of_youngest_parent(1) == 0.06);
        REQUIRE(tree.get_index_of_youngest_parent(1) == 2);
        REQUIRE(tree.get_youngest_parent(2)->get_height() == 0.1);
        REQUIRE(tree.get_height_of_youngest_parent(2) == 0.1);
        REQUIRE(tree.get_index_of_youngest_parent(2) == 4);
        REQUIRE(tree.get_youngest_parent(3)->get_height() == 0.1);
        REQUIRE(tree.get_height_of_youngest_parent(3) == 0.1);
        REQUIRE(tree.get_index_of_youngest_parent(3) == 4);
    }
}

TEST_CASE("phylonet; Testing BaseTree::get_height_of_oldest_child()", "[BaseTree]") {
    SECTION("Testing get_height_of_oldest_child") {
        std::shared_ptr<NetNode> root = std::make_shared<NetNode>("root", 0.1);
        std::shared_ptr<NetNode> internal1 = std::make_shared<NetNode>("internal1", 0.08);
        std::shared_ptr<NetNode> internal2 = std::make_shared<NetNode>("internal2", 0.06);
        std::shared_ptr<NetNode> internal3 = std::make_shared<NetNode>("internal3", 0.04);
        std::shared_ptr<NetNode> internal4 = std::make_shared<NetNode>("internal4", 0.04);
        internal4->set_height_parameter(internal3->get_height_parameter());
        std::shared_ptr<NetNode> internal5 = std::make_shared<NetNode>("internal5", 0.02);
        std::shared_ptr<NetNode> internal6 = std::make_shared<NetNode>("internal6", 0.02);
        internal6->set_height_parameter(internal5->get_height_parameter());
        std::shared_ptr<NetNode> leaf1 = std::make_shared<NetNode>("leaf1", 0.0);
        leaf1->fix_node_height();
        std::shared_ptr<NetNode> leaf2 = std::make_shared<NetNode>("leaf2", 0.0);
        leaf2->fix_node_height();
        std::shared_ptr<NetNode> leaf3 = std::make_shared<NetNode>("leaf3", 0.0);
        leaf3->fix_node_height();
        std::shared_ptr<NetNode> leaf4 = std::make_shared<NetNode>("leaf4", 0.0);
        leaf4->fix_node_height();
        std::shared_ptr<NetNode> leaf5 = std::make_shared<NetNode>("leaf5", 0.0);
        leaf5->fix_node_height();
        std::shared_ptr<NetNode> leaf6 = std::make_shared<NetNode>("leaf6", 0.0);
        leaf6->fix_node_height();
        std::shared_ptr<NetNode> leaf7 = std::make_shared<NetNode>("leaf7", 0.0);
        leaf7->fix_node_height();
        std::shared_ptr<NetNode> leaf8 = std::make_shared<NetNode>("leaf8", 0.0);
        leaf8->fix_node_height();

        internal3->add_child(leaf1);
        internal3->add_child(leaf2);
        internal4->add_child(leaf3);
        internal4->add_child(leaf4);
        internal5->add_child(leaf5);
        internal5->add_child(leaf6);
        internal6->add_child(leaf7);
        internal6->add_child(leaf8);

        internal1->add_child(internal5);
        internal1->add_child(internal3);
        internal2->add_child(internal6);
        internal2->add_child(internal4);

        root->add_child(internal1);
        root->add_child(internal2);
        BaseTree<NetNode> tree(root);

        REQUIRE(tree.get_leaf_node_count() == 8);
        REQUIRE(tree.get_node_count() == 15);
        REQUIRE(tree.get_number_of_node_heights() == 5);
        std::vector<double> expected_heights {0.02, 0.04, 0.06, 0.08, 0.1};
        REQUIRE(tree.get_node_heights() == expected_heights);

        REQUIRE(tree.get_oldest_child(0)->get_height() == 0.0);
        REQUIRE(tree.get_height_of_oldest_child(0) == 0.0);
        REQUIRE(tree.get_oldest_child(1)->get_height() == 0.0);
        REQUIRE(tree.get_height_of_oldest_child(1) == 0.0);
        REQUIRE(tree.get_oldest_child(2)->get_height() == 0.04);
        REQUIRE(tree.get_height_of_oldest_child(2) == 0.04);
        REQUIRE(tree.get_oldest_child(3)->get_height() == 0.04);
        REQUIRE(tree.get_height_of_oldest_child(3) == 0.04);
        REQUIRE(tree.get_oldest_child(4)->get_height() == 0.08);
        REQUIRE(tree.get_height_of_oldest_child(4) == 0.08);
    }
}

TEST_CASE("phylonet; Testing BaseTree::store_splits()", "[BaseTree]") {
    SECTION("Testing store_splits") {
        std::shared_ptr<NetNode> root = std::make_shared<NetNode>(14, "root", 0.1);
        std::shared_ptr<NetNode> internal1 = std::make_shared<NetNode>(12, "internal1", 0.08);
        std::shared_ptr<NetNode> internal2 = std::make_shared<NetNode>(13, "internal2", 0.06);
        std::shared_ptr<NetNode> internal3 = std::make_shared<NetNode>(8, "internal3", 0.04);
        std::shared_ptr<NetNode> internal4 = std::make_shared<NetNode>(9, "internal4", 0.04);
        internal4->set_height_parameter(internal3->get_height_parameter());
        std::shared_ptr<NetNode> internal5 = std::make_shared<NetNode>(10, "internal5", 0.02);
        std::shared_ptr<NetNode> internal6 = std::make_shared<NetNode>(11, "internal6", 0.02);
        internal6->set_height_parameter(internal5->get_height_parameter());
        std::shared_ptr<NetNode> leaf1 = std::make_shared<NetNode>(0, "leaf1", 0.0);
        leaf1->fix_node_height();
        std::shared_ptr<NetNode> leaf2 = std::make_shared<NetNode>(1, "leaf2", 0.0);
        leaf2->fix_node_height();
        std::shared_ptr<NetNode> leaf3 = std::make_shared<NetNode>(2, "leaf3", 0.0);
        leaf3->fix_node_height();
        std::shared_ptr<NetNode> leaf4 = std::make_shared<NetNode>(3, "leaf4", 0.0);
        leaf4->fix_node_height();
        std::shared_ptr<NetNode> leaf5 = std::make_shared<NetNode>(4, "leaf5", 0.0);
        leaf5->fix_node_height();
        std::shared_ptr<NetNode> leaf6 = std::make_shared<NetNode>(5, "leaf6", 0.0);
        leaf6->fix_node_height();
        std::shared_ptr<NetNode> leaf7 = std::make_shared<NetNode>(6, "leaf7", 0.0);
        leaf7->fix_node_height();
        std::shared_ptr<NetNode> leaf8 = std::make_shared<NetNode>(7, "leaf8", 0.0);
        leaf8->fix_node_height();

        internal3->add_child(leaf1);
        internal3->add_child(leaf2);
        internal4->add_child(leaf3);
        internal4->add_child(leaf4);
        internal5->add_child(leaf5);
        internal5->add_child(leaf6);
        internal6->add_child(leaf7);
        internal6->add_child(leaf8);

        internal1->add_child(internal5);
        internal1->add_child(internal3);
        internal2->add_child(internal6);
        internal2->add_child(internal4);

        root->add_child(internal1);
        root->add_child(internal2);
        BaseTree<NetNode> tree(root);

        std::map< int, std::set<Split> > split_set;
        tree.store_splits_by_height_index(split_set);
        std::cout << "split set:\n";
        for (auto height_splits : split_set) {
            std::cout << height_splits.first << ": ";
            unsigned int split_count = 0;
            for (auto split : height_splits.second) {
                if (split_count > 0) {
                    std::cout << "   ";
                }
                std::cout << split.as_string() << "\n";
                ++split_count;
            }
        }

        Split s;
        s.resize(8);
        std::map< int, std::set<Split> > expected_set;
        for (unsigned int i = 0; i < 8; ++i) {
            s.set_leaf_bit(i);
        }
        expected_set[4].insert(s);

        s.clear();
        for (auto i : {0, 1, 4, 5}) {
            s.set_leaf_bit(i);
        }
        expected_set[3].insert(s);

        s.clear();
        for (auto i : {2, 3, 6, 7}) {
            s.set_leaf_bit(i);
        }
        expected_set[2].insert(s);

        s.clear();
        for (auto i : {0, 1}) {
            s.set_leaf_bit(i);
        }
        expected_set[1].insert(s);

        s.clear();
        for (auto i : {2, 3}) {
            s.set_leaf_bit(i);
        }
        expected_set[1].insert(s);

        s.clear();
        for (auto i : {4, 5}) {
            s.set_leaf_bit(i);
        }
        expected_set[0].insert(s);

        s.clear();
        for (auto i : {6, 7}) {
            s.set_leaf_bit(i);
        }
        expected_set[0].insert(s);

        std::cout << "expected split set:\n";
        for (auto height_splits : expected_set) {
            std::cout << height_splits.first << ": ";
            unsigned int split_count = 0;
            for (auto split : height_splits.second) {
                if (split_count > 0) {
                    std::cout << "   ";
                }
                std::cout << split.as_string() << "\n";
                ++split_count;
            }
        }

        REQUIRE(split_set == expected_set);

        // What if a height index is off?
        std::map< int, std::set<Split> > index_off_set;
        for (unsigned int i = 0; i < 8; ++i) {
            s.set_leaf_bit(i);
        }
        index_off_set[4].insert(s);

        s.clear();
        for (auto i : {0, 1, 4, 5}) {
            s.set_leaf_bit(i);
        }
        index_off_set[3].insert(s);

        s.clear();
        for (auto i : {2, 3, 6, 7}) {
            s.set_leaf_bit(i);
        }
        index_off_set[2].insert(s);

        s.clear();
        for (auto i : {0, 1}) {
            s.set_leaf_bit(i);
        }
        index_off_set[1].insert(s);

        s.clear();
        for (auto i : {2, 3}) {
            s.set_leaf_bit(i);
        }
        index_off_set[1].insert(s);

        s.clear();
        for (auto i : {4, 5}) {
            s.set_leaf_bit(i);
        }
        index_off_set[1].insert(s);

        s.clear();
        for (auto i : {6, 7}) {
            s.set_leaf_bit(i);
        }
        index_off_set[0].insert(s);

        REQUIRE(split_set != index_off_set);

        // What if one leaf index is off?
        std::map< int, std::set<Split> > one_leaf_off_set;
        for (unsigned int i = 0; i < 8; ++i) {
            s.set_leaf_bit(i);
        }
        one_leaf_off_set[4].insert(s);

        s.clear();
        for (auto i : {0, 1, 4, 5}) {
            s.set_leaf_bit(i);
        }
        one_leaf_off_set[3].insert(s);

        s.clear();
        for (auto i : {2, 3, 6, 7}) {
            s.set_leaf_bit(i);
        }
        one_leaf_off_set[2].insert(s);

        s.clear();
        for (auto i : {0, 4}) {
            s.set_leaf_bit(i);
        }
        one_leaf_off_set[1].insert(s);

        s.clear();
        for (auto i : {2, 3}) {
            s.set_leaf_bit(i);
        }
        one_leaf_off_set[1].insert(s);

        s.clear();
        for (auto i : {1, 5}) {
            s.set_leaf_bit(i);
        }
        one_leaf_off_set[0].insert(s);

        s.clear();
        for (auto i : {6, 7}) {
            s.set_leaf_bit(i);
        }
        one_leaf_off_set[0].insert(s);

        REQUIRE(split_set != one_leaf_off_set);
    }
}
