#include "catch.hpp"
#include "ecoevolity/ecoevolity.hpp"

#include "ecoevolity/path.hpp"
#include "ecoevolity/general_tree_operator.hpp"
#include "ecoevolity/tree.hpp"
#include "ecoevolity/node.hpp"
#include "ecoevolity/probability.hpp"
#include "ecoevolity/parameter.hpp"
#include "ecoevolity/stats_util.hpp"
#include "ecoevolity/rng.hpp"

TEST_CASE("Testing Watterson's theta with one pop with missing data",
        "[BasePopulationTree]") {

    SECTION("Testing sim of 1 species with missing data") {
        RandomNumberGenerator rng = RandomNumberGenerator(123);

        std::string yml_path = "data/species-1-genomes-20-chars-10000-missing.yml";
        BasePopulationTree tree(yml_path,
                ' ',    // population_name_delimiter
                true,   // population_name_is_prefix
                false,  // genotypes_are_diploid
                false,  // markers_are_dominant
                false,  // constant_sites_removed
                true,   // validate
                true,   // strict_on_constant_sites
                true,   // strict_on_missing_sites
                true,   // strict_on_triallelic_sites
                2.0,    // ploidy
                false); //store_seq_loci_info
        REQUIRE(tree.constant_sites_removed() == false);

        double ne = 0.001;
        double expected_theta = 4.0 * ne;
        tree.get_node("pop1")->set_population_size(ne);
        tree.get_root_ptr()->set_population_size(ne);

        tree.estimate_mutation_rate();
        tree.set_mutation_rate(1.0);

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.5);

        tree.fix_mutation_rate();
        tree.fix_state_frequencies();

        tree.estimate_root_height();
        tree.get_root_ptr()->set_height(0.5);

        SampleSummarizer<double> theta_summary;

        for (unsigned int i = 0; i < 10; ++i) {
            BiallelicData sbd = tree.simulate_biallelic_data_set(
                    rng,  // RandomNumberGenerator
                    1.0,                // singleton_sample_probability
                    true);              // validate
            theta_summary.add_sample(sbd.get_wattersons_theta("pop1"));
        }
        std::cout << "Expected theta: " << expected_theta << "\n";
        std::cout << "Mean Watterson's theta: " << theta_summary.mean() << "\n";
        double eps = 0.001;
        REQUIRE(theta_summary.mean() == Approx(expected_theta).epsilon(eps));
    }
}
