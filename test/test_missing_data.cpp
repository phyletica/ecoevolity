#include "catch.hpp"
#include "ecoevolity/ecoevolity.hpp"

#include "ecoevolity/path.hpp"
#include "ecoevolity/general_tree_operator.hpp"
#include "ecoevolity/poptree.hpp"
#include "ecoevolity/node.hpp"
#include "ecoevolity/probability.hpp"
#include "ecoevolity/parameter.hpp"
#include "ecoevolity/stats_util.hpp"
#include "ecoevolity/rng.hpp"

RandomNumberGenerator _MISSING_DATA_RNG = RandomNumberGenerator();

TEST_CASE("Testing simulation and MCMC of 2 pops with no shared data",
        "[BasePopulationTree]") {

    SECTION("Testing sim and MCMC of 2 species with no shared data") {
        RandomNumberGenerator rng = RandomNumberGenerator(123);

        std::string yml_path = "data/species-2-genomes-20-chars-10000-no-overlap.yml";
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
        BiallelicData bd = tree.get_data();

        // Make sure the data from the file match our expections
        REQUIRE(bd.has_missing_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_number_of_missing_population_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 0);
        unsigned int number_removed = bd.remove_missing_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.has_missing_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_number_of_missing_population_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 0);

        std::vector< std::vector<unsigned int> > expected_allele_counts(2);
        expected_allele_counts[0] = {20, 0};
        expected_allele_counts[1] = {0, 20};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5001;
        expected_unique_allele_counts[expected_allele_counts.at(1)] = 4999;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::map<unsigned int, unsigned int> pop1_expected_unique_allele_counts;
        pop1_expected_unique_allele_counts[20] = 5001;
        pop1_expected_unique_allele_counts[0] = 4999;
        std::map<unsigned int, unsigned int> pop2_expected_unique_allele_counts;
        pop2_expected_unique_allele_counts[20] = 4999;
        pop2_expected_unique_allele_counts[0] = 5001;
        REQUIRE(bd.get_unique_allele_counts_for_population("pop1") == pop1_expected_unique_allele_counts);
        REQUIRE(bd.get_unique_allele_counts_for_population("pop2") == pop2_expected_unique_allele_counts);

        double pop1_ne = 0.001;
        double pop2_ne = 0.005;
        double root_ne = 0.0005;
        double mu = 1.0;
        double time_root = 0.25;

        tree.get_node("pop1")->set_population_size(pop1_ne);
        tree.get_node("pop2")->set_population_size(pop2_ne);
        tree.get_root_ptr()->set_population_size(root_ne);

        tree.estimate_mutation_rate();
        tree.set_mutation_rate(mu);

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.5);

        tree.fix_mutation_rate();
        tree.fix_state_frequencies();

        tree.estimate_root_height();
        tree.get_root_ptr()->set_height(time_root);

        // Simulate a dataset and make sure the allele counts match the
        // original data to confirm the simulation code is working
        BiallelicData sbd = tree.simulate_biallelic_data_set(
                rng,  // RandomNumberGenerator
                1.0,                // singleton_sample_probability
                true);              // validate
        REQUIRE(sbd.has_missing_patterns() == false);
        REQUIRE(sbd.has_missing_population_patterns() == true);
        REQUIRE(sbd.get_number_of_missing_population_sites_removed() == 0);
        REQUIRE(sbd.get_number_of_missing_sites_removed() == 0);
        number_removed = sbd.remove_missing_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(sbd.has_missing_patterns() == false);
        REQUIRE(sbd.has_missing_population_patterns() == true);
        REQUIRE(sbd.get_number_of_missing_population_sites_removed() == 0);
        REQUIRE(sbd.get_number_of_missing_sites_removed() == 0);
        REQUIRE(sbd.get_unique_allele_counts() == expected_unique_allele_counts);
        REQUIRE(sbd.get_unique_allele_counts_for_population("pop1") == pop1_expected_unique_allele_counts);
        REQUIRE(sbd.get_unique_allele_counts_for_population("pop2") == pop2_expected_unique_allele_counts);

        tree.set_data(sbd, false);

        REQUIRE(tree.constant_sites_removed() == false);

        // Now let's see if sampling from the posterior with MCMC behaves as we
        // expect given the missing data
        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double pop_size_shape = 5.0;
        double pop_size_scale = 0.002;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        tree.set_population_size_prior(pop_size_prior);
        tree.set_root_node_height_prior(root_height_prior);

        RootHeightSizeMixer op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        PopSizeScaler op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);
        RootHeightScaler< BasePopulationTree > op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);

        // Initialize prior probs
        tree.draw_from_prior(rng);
        tree.compute_log_likelihood_and_prior(true);

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 3);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> root_size_summary;
        SampleSummarizer<double> pop1_size_summary;
        SampleSummarizer<double> pop2_size_summary;

        // std::shared_ptr<PopulationNode> node;
        unsigned int nmoves_per_op = 1;
        unsigned int niterations = 10000;
        unsigned int nburnin = 100;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < nburnin; ++i) {
            op.operate(rng, &tree, 1);
            op2.operate(rng, &tree, 1, nmoves_per_op);
            op3.operate(rng, &tree, 1, nmoves_per_op);
        }
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            op2.operate(rng, &tree, 1, nmoves_per_op);
            op3.operate(rng, &tree, 1, nmoves_per_op);
            if ((i + 1) % sample_freq == 0) {
                root_height_summary.add_sample(tree.get_root_height());

                // node = tree.get_root_ptr();
                root_size_summary.add_sample(tree.get_root_ptr()->get_population_size());

                // node = tree.get_node("pop1");
                pop1_size_summary.add_sample(tree.get_node("pop1")->get_population_size());
                // node = tree.get_node("pop2");
                pop2_size_summary.add_sample(tree.get_node("pop2")->get_population_size());
                if ((i + 1) % (sample_freq * 100) == 0) {
                    std::cout << "MCMC iter " << (i + 1) << " of " << niterations << "\n";
                }
            }
        }
        std::cout << "\n";
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();
        std::cout << op3.header_string();
        std::cout << op3.to_string();

        std::cout << "\n";
        std::cout << "Pop1 Watterson's theta: " << sbd.get_wattersons_theta("pop1") << "\n";
        std::cout << "Pop2 Watterson's theta: " << sbd.get_wattersons_theta("pop2") << "\n";

        std::cout << "\n";
        std::cout << "Pop size prior mean: " << pop_size_prior->get_mean() << "\n";
        std::cout << "Pop size prior var: " << pop_size_prior->get_variance() << "\n";
        std::cout << "\n";
        std::cout << "True root pop size: " << root_ne << "\n";
        std::cout << "Root pop size mean: " << root_size_summary.mean() << "\n";
        std::cout << "Root pop size var: " << root_size_summary.variance() << "\n";
        std::cout << "\n";
        std::cout << "True pop1 pop size: " << pop1_ne << "\n";
        std::cout << "pop1 size mean: " << pop1_size_summary.mean() << "\n";
        std::cout << "pop1 size var: " << pop1_size_summary.variance() << "\n";
        std::cout << "\n";
        std::cout << "True pop2 pop size: " << pop2_ne << "\n";
        std::cout << "pop2 size mean: " << pop2_size_summary.mean() << "\n";
        std::cout << "pop2 size var: " << pop2_size_summary.variance() << "\n";

        std::cout << "\n";
        std::cout << "Root height prior mean: " << root_height_prior->get_mean() << "\n";
        std::cout << "Root height prior var: " << root_height_prior->get_variance() << "\n";
        std::cout << "True root height: " << time_root << "\n";
        std::cout << "Root height mean: " << root_height_summary.mean() << "\n";
        std::cout << "Root height var: " << root_height_summary.variance() << "\n";

        double eps = 0.01;

        // Given the true divergence time is much larger than the leaf pop
        // sizes, and no characters overlap between the 2 leaf populations, no
        // coalescence is observed in the root population, and thus there
        // should be no information in the data about the root height or pop
        // size and MCMC should sample the prior for both these parameters
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));

        REQUIRE(root_size_summary.mean() > root_ne);
        REQUIRE(root_size_summary.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
        REQUIRE(root_size_summary.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));

        // There should be plenty of info in the data to estimate the size of
        // the leaf populations, which are much smaller than expected under the
        // prior. Let's make sure our posterior means are in the ballpark of
        // the true parameter values.
        REQUIRE(pop1_size_summary.mean() == Approx(pop1_ne).epsilon(0.001));
        REQUIRE(pop2_size_summary.mean() == Approx(pop2_ne).epsilon(0.001));
    }
}

TEST_CASE("Testing simulation and MCMC of 2 pops with one missing data",
        "[BasePopulationTree]") {

    SECTION("Testing sim and MCMC of 2 species with one missing data") {
        RandomNumberGenerator rng = RandomNumberGenerator(123);

        std::string yml_path = "data/species-2-genomes-20-chars-5000-one-sp-missing.yml";
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
        BiallelicData bd = tree.get_data();

        // Make sure the data from the file match our expections
        REQUIRE(bd.has_missing_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_number_of_missing_population_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 0);
        unsigned int number_removed = bd.remove_missing_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(bd.has_missing_patterns() == false);
        REQUIRE(bd.has_missing_population_patterns() == true);
        REQUIRE(bd.get_number_of_missing_population_sites_removed() == 0);
        REQUIRE(bd.get_number_of_missing_sites_removed() == 0);

        std::vector< std::vector<unsigned int> > expected_allele_counts(1);
        expected_allele_counts[0] = {20, 0};

        std::map<std::vector<unsigned int>, unsigned int> expected_unique_allele_counts;
        expected_unique_allele_counts[expected_allele_counts.at(0)] = 5000;
        REQUIRE(bd.get_unique_allele_counts() == expected_unique_allele_counts);

        std::map<unsigned int, unsigned int> pop1_expected_unique_allele_counts;
        pop1_expected_unique_allele_counts[20] = 5000;
        std::map<unsigned int, unsigned int> pop2_expected_unique_allele_counts;
        pop2_expected_unique_allele_counts[0] = 5000;
        REQUIRE(bd.get_unique_allele_counts_for_population("pop1") == pop1_expected_unique_allele_counts);
        REQUIRE(bd.get_unique_allele_counts_for_population("pop2") == pop2_expected_unique_allele_counts);

        double pop1_ne = 0.001;
        double pop2_ne = 0.004;
        double root_ne = 0.001;
        double mu = 1.0;
        double time_root = 0.05;

        tree.get_node("pop1")->set_population_size(pop1_ne);
        tree.get_node("pop2")->set_population_size(pop2_ne);
        tree.get_root_ptr()->set_population_size(root_ne);

        tree.estimate_mutation_rate();
        tree.set_mutation_rate(mu);

        tree.estimate_state_frequencies();
        tree.set_freq_1(0.5);

        tree.fix_mutation_rate();
        tree.fix_state_frequencies();

        tree.estimate_root_height();
        tree.get_root_ptr()->set_height(time_root);

        // Simulate a dataset and make sure the allele counts match the
        // original data to confirm the simulation code is working
        BiallelicData sbd = tree.simulate_biallelic_data_set(
                rng,  // RandomNumberGenerator
                1.0,                // singleton_sample_probability
                true);              // validate
        REQUIRE(sbd.has_missing_patterns() == false);
        REQUIRE(sbd.has_missing_population_patterns() == true);
        REQUIRE(sbd.get_number_of_missing_population_sites_removed() == 0);
        REQUIRE(sbd.get_number_of_missing_sites_removed() == 0);
        number_removed = sbd.remove_missing_patterns();
        REQUIRE(number_removed == 0);
        REQUIRE(sbd.has_missing_patterns() == false);
        REQUIRE(sbd.has_missing_population_patterns() == true);
        REQUIRE(sbd.get_number_of_missing_population_sites_removed() == 0);
        REQUIRE(sbd.get_number_of_missing_sites_removed() == 0);
        REQUIRE(sbd.get_unique_allele_counts() == expected_unique_allele_counts);
        REQUIRE(sbd.get_unique_allele_counts_for_population("pop1") == pop1_expected_unique_allele_counts);
        REQUIRE(sbd.get_unique_allele_counts_for_population("pop2") == pop2_expected_unique_allele_counts);

        tree.set_data(sbd, false);

        REQUIRE(tree.constant_sites_removed() == false);

        // Now let's see if sampling from the posterior with MCMC behaves as we
        // expect given the missing data
        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double pop_size_shape = 5.0;
        double pop_size_scale = 0.002;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        tree.set_population_size_prior(pop_size_prior);
        tree.set_root_node_height_prior(root_height_prior);

        RootHeightSizeMixer op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        PopSizeScaler op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);
        RootHeightScaler< BasePopulationTree > op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);

        // Initialize prior probs
        tree.draw_from_prior(rng);
        tree.compute_log_likelihood_and_prior(true);

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 3);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> root_size_summary;
        SampleSummarizer<double> pop1_size_summary;
        SampleSummarizer<double> pop2_size_summary;

        // std::shared_ptr<PopulationNode> node;
        unsigned int nmoves_per_op = 1;
        unsigned int niterations = 10000;
        unsigned int nburnin = 100;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < nburnin; ++i) {
            op.operate(rng, &tree, 1);
            op2.operate(rng, &tree, 1, nmoves_per_op);
            op3.operate(rng, &tree, 1, nmoves_per_op);
        }
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            op2.operate(rng, &tree, 1, nmoves_per_op);
            op3.operate(rng, &tree, 1, nmoves_per_op);
            if ((i + 1) % sample_freq == 0) {
                root_height_summary.add_sample(tree.get_root_height());

                // node = tree.get_root_ptr();
                root_size_summary.add_sample(tree.get_root_ptr()->get_population_size());

                // node = tree.get_node("pop1");
                pop1_size_summary.add_sample(tree.get_node("pop1")->get_population_size());
                // node = tree.get_node("pop2");
                pop2_size_summary.add_sample(tree.get_node("pop2")->get_population_size());
                if ((i + 1) % (sample_freq * 100) == 0) {
                    std::cout << "MCMC iter " << (i + 1) << " of " << niterations << "\n";
                }
            }
        }
        std::cout << "\n";
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();
        std::cout << op3.header_string();
        std::cout << op3.to_string();

        std::cout << "\n";
        std::cout << "Pop size prior mean: " << pop_size_prior->get_mean() << "\n";
        std::cout << "Pop size prior var: " << pop_size_prior->get_variance() << "\n";
        std::cout << "\n";
        std::cout << "True root pop size: " << root_ne << "\n";
        std::cout << "Root pop size mean: " << root_size_summary.mean() << "\n";
        std::cout << "Root pop size var: " << root_size_summary.variance() << "\n";
        std::cout << "\n";
        std::cout << "True pop1 pop size: " << pop1_ne << "\n";
        std::cout << "pop1 size mean: " << pop1_size_summary.mean() << "\n";
        std::cout << "pop1 size var: " << pop1_size_summary.variance() << "\n";
        std::cout << "\n";
        std::cout << "True pop2 pop size: " << pop2_ne << "\n";
        std::cout << "pop2 size mean: " << pop2_size_summary.mean() << "\n";
        std::cout << "pop2 size var: " << pop2_size_summary.variance() << "\n";

        std::cout << "\n";
        std::cout << "Root height prior mean: " << root_height_prior->get_mean() << "\n";
        std::cout << "Root height prior var: " << root_height_prior->get_variance() << "\n";
        std::cout << "True root height: " << time_root << "\n";
        std::cout << "Root height mean: " << root_height_summary.mean() << "\n";
        std::cout << "Root height var: " << root_height_summary.variance() << "\n";

        double eps = 0.01;

        // The samples from Pop 1 should inform the size of the pop 1 but there
        // should be no information about the size of pop 2, the root, or the
        // divergence time, because there is no data from Pop2 and there was no
        // change in pop size between Pop 1 and the root
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));

        REQUIRE(pop2_size_summary.mean() > pop2_ne);
        REQUIRE(pop2_size_summary.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
        REQUIRE(pop2_size_summary.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));

        REQUIRE(root_size_summary.mean() > root_ne);
        REQUIRE(root_size_summary.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
        REQUIRE(root_size_summary.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));

        REQUIRE(pop1_size_summary.mean() == Approx(pop1_ne).epsilon(0.001));
    }
}
