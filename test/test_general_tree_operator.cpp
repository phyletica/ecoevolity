#include "catch.hpp"
#include "ecoevolity/general_tree_operator.hpp"

#include <limits>
#include <memory>

#include "ecoevolity/probability.hpp"
#include "ecoevolity/parameter.hpp"
#include "ecoevolity/stats_util.hpp"
#include "ecoevolity/rng.hpp"


TEST_CASE("Testing NodeHeightSlideBumpScaler with 3 leaves, fixed root, and optimizing",
        "[NodeHeightSlideBumpScaler]") {

    SECTION("Testing 3 leaves with fixed root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(1);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        NodeHeightSlideBumpScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 200000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0));
                REQUIRE(tree.get_height(1) == tree.get_root_height());
                REQUIRE(tree.get_height(1) == 1.0);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.001;
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing NodeHeightSlideBumpScaler with 3 leaves, fixed root, and no optimizing",
        "[NodeHeightSlideBumpScaler]") {

    SECTION("Testing 3 leaves with fixed root and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(2);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        NodeHeightSlideBumpScaler<Node> op;
        op.turn_off_auto_optimize();
        op.set_coercable_parameter_value(1.8);

        REQUIRE(! op.auto_optimizing());
        REQUIRE(op.get_coercable_parameter_value() == 1.8);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 200000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0));
                REQUIRE(tree.get_height(1) == tree.get_root_height());
                REQUIRE(tree.get_height(1) == 1.0);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - op.get_auto_optimize_delay()));
        REQUIRE(op.get_coercable_parameter_value() == 1.8);

        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.001;
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing NodeHeightSlideBumpScaler with 3 leaves, variable root, and optimizing",
        "[NodeHeightSlideBumpScaler]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(3);

        double root_height_lower = 0.999;
        double root_height_upper = 1.001;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<UniformDistribution>(
                root_height_lower,
                root_height_upper);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightSlideBumpScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 200000;
        unsigned int burnin = 10000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1);
        }
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0));
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(1) == tree.get_root_height());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        SampleSummarizer<double> expected_summary;
        for (unsigned int i = 0; i < niterations; ++i) {
            double root_ht = root_height_prior->draw(rng);
            double sample = rng.uniform_real(0.0, root_ht);
            expected_summary.add_sample(sample);
        }

        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(expected_summary.mean()).epsilon(eps * 5));
        REQUIRE(internal_height_summary.variance() == Approx(expected_summary.variance()).epsilon(eps * 5));
    }
}

TEST_CASE("Testing NodeHeightSlideBumpScaler with 3 leaves, gamma root, and optimizing",
        "[NodeHeightSlideBumpScaler]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(4);

        double root_height_shape = 100.0;
        double root_height_scale = 0.01;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightSlideBumpScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 200000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0));
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(1) == tree.get_root_height());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        SampleSummarizer<double> expected_summary;
        for (unsigned int i = 0; i < niterations; ++i) {
            double root_ht = root_height_prior->draw(rng);
            double sample = rng.uniform_real(0.0, root_ht);
            expected_summary.add_sample(sample);
        }

        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(expected_summary.mean()).epsilon(eps * 5));
        REQUIRE(internal_height_summary.variance() == Approx(expected_summary.variance()).epsilon(eps * 5));
    }
}



TEST_CASE("Testing NodeHeightSlideBumpMover with 3 leaves, fixed root, and optimizing",
        "[NodeHeightSlideBumpMover]") {

    SECTION("Testing 3 leaves with fixed root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(50000);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        NodeHeightSlideBumpMover<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 200000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0));
                REQUIRE(tree.get_height(1) == tree.get_root_height());
                REQUIRE(tree.get_height(1) == 1.0);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.001;
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing NodeHeightSlideBumpMover with 3 leaves, fixed root, and no optimizing",
        "[NodeHeightSlideBumpMover]") {

    SECTION("Testing 3 leaves with fixed root and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(6);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        NodeHeightSlideBumpMover<Node> op;
        op.turn_off_auto_optimize();
        op.set_coercable_parameter_value(1.8);

        REQUIRE(! op.auto_optimizing());
        REQUIRE(op.get_coercable_parameter_value() == 1.8);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 200000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0));
                REQUIRE(tree.get_height(1) == tree.get_root_height());
                REQUIRE(tree.get_height(1) == 1.0);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - op.get_auto_optimize_delay()));
        REQUIRE(op.get_coercable_parameter_value() == 1.8);

        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.001;
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing NodeHeightSlideBumpMover with 3 leaves, variable root, and optimizing",
        "[NodeHeightSlideBumpMover]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(7);

        double root_height_lower = 0.999;
        double root_height_upper = 1.001;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<UniformDistribution>(
                root_height_lower,
                root_height_upper);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightSlideBumpMover<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 400000;
        unsigned int burnin = 1000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1);
        }
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0));
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(1) == tree.get_root_height());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        SampleSummarizer<double> expected_summary;
        for (unsigned int i = 0; i < niterations; ++i) {
            double root_ht = root_height_prior->draw(rng);
            double sample = rng.uniform_real(0.0, root_ht);
            expected_summary.add_sample(sample);
        }

        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(expected_summary.mean()).epsilon(eps * 5));
        REQUIRE(internal_height_summary.variance() == Approx(expected_summary.variance()).epsilon(eps * 5));
    }
}

TEST_CASE("Testing NodeHeightSlideBumpMover with 3 leaves, gamma root, and optimizing",
        "[NodeHeightSlideBumpMover]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(8);

        double root_height_shape = 100.0;
        double root_height_scale = 0.01;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightSlideBumpMover<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 200000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0));
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(1) == tree.get_root_height());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        SampleSummarizer<double> expected_summary;
        for (unsigned int i = 0; i < niterations; ++i) {
            double root_ht = root_height_prior->draw(rng);
            double sample = rng.uniform_real(0.0, root_ht);
            expected_summary.add_sample(sample);
        }

        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(expected_summary.mean()).epsilon(eps * 5));
        REQUIRE(internal_height_summary.variance() == Approx(expected_summary.variance()).epsilon(eps * 5));
    }
}




TEST_CASE("Testing NodeHeightSlideBumpSwapScaler with 3 leaves, gamma root, and optimizing",
        "[NodeHeightSlideBumpSwapScaler]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(9);

        double root_height_shape = 100.0;
        double root_height_scale = 0.01;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightSlideBumpSwapScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_12 = 0;

        unsigned int niterations = 400000;
        unsigned int sample_freq = 20;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < 100000; ++i) {
            op.operate(rng, &tree, 1);
        }
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0));
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(1) == tree.get_root_height());
                if (tree.get_root_ptr()->is_child("leaf0")) {
                    ++count_12;
                }
                if (tree.get_root_ptr()->is_child("leaf1")) {
                    ++count_02;
                }
                if (tree.get_root_ptr()->is_child("leaf2")) {
                    ++count_01;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        SampleSummarizer<double> expected_summary;
        for (unsigned int i = 0; i < niterations; ++i) {
            double root_ht = root_height_prior->draw(rng);
            double sample = rng.uniform_real(0.0, root_ht);
            expected_summary.add_sample(sample);
        }

        double eps = 0.005;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(expected_summary.mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(expected_summary.variance()).epsilon(eps * 2));

        double freq_01 = count_01 / (double)nsamples;
        double freq_02 = count_02 / (double)nsamples;
        double freq_12 = count_12 / (double)nsamples;
        std::cout << "freq of ((0,1),2): " << freq_01 << "\n";
        std::cout << "freq of ((0,2),1): " << freq_02 << "\n";
        std::cout << "freq of ((1,2),0): " << freq_12 << "\n";
        REQUIRE(freq_01 == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(freq_02 == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(freq_12 == Approx(1.0/3.0).epsilon(eps));
    }
}

TEST_CASE("Testing NodeHeightSlideBumpSwapMover with 3 leaves, gamma root, and optimizing",
        "[NodeHeightSlideBumpSwapMover]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(10);

        double root_height_shape = 100.0;
        double root_height_scale = 0.01;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightSlideBumpSwapMover<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_12 = 0;

        unsigned int niterations = 400000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        /* for (unsigned int i = 0; i < 100000; ++i) { */
        /*     op.operate(rng, &tree, 1); */
        /* } */
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0));
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(1) == tree.get_root_height());
                if (tree.get_root_ptr()->is_child("leaf0")) {
                    ++count_12;
                }
                if (tree.get_root_ptr()->is_child("leaf1")) {
                    ++count_02;
                }
                if (tree.get_root_ptr()->is_child("leaf2")) {
                    ++count_01;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        SampleSummarizer<double> expected_summary;
        for (unsigned int i = 0; i < niterations; ++i) {
            double root_ht = root_height_prior->draw(rng);
            double sample = rng.uniform_real(0.0, root_ht);
            expected_summary.add_sample(sample);
        }

        double eps = 0.005;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(expected_summary.mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(expected_summary.variance()).epsilon(eps * 2));

        double freq_01 = count_01 / (double)nsamples;
        double freq_02 = count_02 / (double)nsamples;
        double freq_12 = count_12 / (double)nsamples;
        std::cout << "freq of ((0,1),2): " << freq_01 << "\n";
        std::cout << "freq of ((0,2),1): " << freq_02 << "\n";
        std::cout << "freq of ((1,2),0): " << freq_12 << "\n";
        REQUIRE(freq_01 == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(freq_02 == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(freq_12 == Approx(1.0/3.0).epsilon(eps));
    }
}



TEST_CASE("Testing NeighborHeightNodeSwap with 3 leaves, gamma root",
        "[NeighborHeightNodeSwap]") {

    SECTION("Testing 3 leaves with variable root") {
        RandomNumberGenerator rng = RandomNumberGenerator(11000);

        double root_height_shape = 100.0;
        double root_height_scale = 0.01;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.estimate_root_height();

        NeighborHeightNodeSwap<Node> op;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_12 = 0;

        unsigned int niterations = 200000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                REQUIRE(tree.get_height(1) == tree.get_root_height());
                if (tree.get_root_ptr()->is_child("leaf0")) {
                    ++count_12;
                }
                if (tree.get_root_ptr()->is_child("leaf1")) {
                    ++count_02;
                }
                if (tree.get_root_ptr()->is_child("leaf2")) {
                    ++count_01;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        
        double eps = 0.005;

        double freq_01 = count_01 / (double)nsamples;
        double freq_02 = count_02 / (double)nsamples;
        double freq_12 = count_12 / (double)nsamples;
        std::cout << "freq of ((0,1),2): " << freq_01 << "\n";
        std::cout << "freq of ((0,2),1): " << freq_02 << "\n";
        std::cout << "freq of ((1,2),0): " << freq_12 << "\n";
        REQUIRE(freq_01 == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(freq_02 == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(freq_12 == Approx(1.0/3.0).epsilon(eps));
    }
}

TEST_CASE("Testing NeighborHeightNodeSwap with 3 leaves, fixed root",
        "[NeighborHeightNodeSwap]") {

    SECTION("Testing 3 leaves with fixed root") {
        RandomNumberGenerator rng = RandomNumberGenerator(12);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        NeighborHeightNodeSwap<Node> op;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_12 = 0;

        unsigned int niterations = 200000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                REQUIRE(tree.get_height(1) == tree.get_root_height());
                if (tree.get_root_ptr()->is_child("leaf0")) {
                    ++count_12;
                }
                if (tree.get_root_ptr()->is_child("leaf1")) {
                    ++count_02;
                }
                if (tree.get_root_ptr()->is_child("leaf2")) {
                    ++count_01;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        
        double eps = 0.005;

        double freq_01 = count_01 / (double)nsamples;
        double freq_02 = count_02 / (double)nsamples;
        double freq_12 = count_12 / (double)nsamples;
        std::cout << "freq of ((0,1),2): " << freq_01 << "\n";
        std::cout << "freq of ((0,2),1): " << freq_02 << "\n";
        std::cout << "freq of ((1,2),0): " << freq_12 << "\n";
        REQUIRE(freq_01 == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(freq_02 == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(freq_12 == Approx(1.0/3.0).epsilon(eps));
    }
}


TEST_CASE("Testing NeighborHeightNodePermute with 3 leaves, gamma root",
        "[NeighborHeightNodePermute]") {

    SECTION("Testing 3 leaves with variable root") {
        RandomNumberGenerator rng = RandomNumberGenerator(13);

        double root_height_shape = 100.0;
        double root_height_scale = 0.01;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.estimate_root_height();

        NeighborHeightNodePermute<Node> op;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_12 = 0;

        unsigned int niterations = 200000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                REQUIRE(tree.get_height(1) == tree.get_root_height());
                if (tree.get_root_ptr()->is_child("leaf0")) {
                    ++count_12;
                }
                if (tree.get_root_ptr()->is_child("leaf1")) {
                    ++count_02;
                }
                if (tree.get_root_ptr()->is_child("leaf2")) {
                    ++count_01;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        
        double eps = 0.005;

        double freq_01 = count_01 / (double)nsamples;
        double freq_02 = count_02 / (double)nsamples;
        double freq_12 = count_12 / (double)nsamples;
        std::cout << "freq of ((0,1),2): " << freq_01 << "\n";
        std::cout << "freq of ((0,2),1): " << freq_02 << "\n";
        std::cout << "freq of ((1,2),0): " << freq_12 << "\n";
        REQUIRE(freq_01 == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(freq_02 == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(freq_12 == Approx(1.0/3.0).epsilon(eps));
    }
}

TEST_CASE("Testing NeighborHeightNodePermute with 3 leaves, fixed root",
        "[NeighborHeightNodePermute]") {

    SECTION("Testing 3 leaves with fixed root") {
        RandomNumberGenerator rng = RandomNumberGenerator(14);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        NeighborHeightNodePermute<Node> op;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_12 = 0;

        unsigned int niterations = 200000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                REQUIRE(tree.get_height(1) == tree.get_root_height());
                if (tree.get_root_ptr()->is_child("leaf0")) {
                    ++count_12;
                }
                if (tree.get_root_ptr()->is_child("leaf1")) {
                    ++count_02;
                }
                if (tree.get_root_ptr()->is_child("leaf2")) {
                    ++count_01;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        
        double eps = 0.005;

        double freq_01 = count_01 / (double)nsamples;
        double freq_02 = count_02 / (double)nsamples;
        double freq_12 = count_12 / (double)nsamples;
        std::cout << "freq of ((0,1),2): " << freq_01 << "\n";
        std::cout << "freq of ((0,2),1): " << freq_02 << "\n";
        std::cout << "freq of ((1,2),0): " << freq_12 << "\n";
        REQUIRE(freq_01 == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(freq_02 == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(freq_12 == Approx(1.0/3.0).epsilon(eps));
    }
}




TEST_CASE("Testing NodeHeightSlideBumpPermuteScaler with 3 leaves, gamma root, and optimizing",
        "[NodeHeightSlideBumpPermuteScaler]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(15);

        double root_height_shape = 100.0;
        double root_height_scale = 0.01;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightSlideBumpPermuteScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_12 = 0;

        unsigned int niterations = 400000;
        unsigned int sample_freq = 20;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < 100000; ++i) {
            op.operate(rng, &tree, 1);
        }
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0));
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(1) == tree.get_root_height());
                if (tree.get_root_ptr()->is_child("leaf0")) {
                    ++count_12;
                }
                if (tree.get_root_ptr()->is_child("leaf1")) {
                    ++count_02;
                }
                if (tree.get_root_ptr()->is_child("leaf2")) {
                    ++count_01;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        SampleSummarizer<double> expected_summary;
        for (unsigned int i = 0; i < niterations; ++i) {
            double root_ht = root_height_prior->draw(rng);
            double sample = rng.uniform_real(0.0, root_ht);
            expected_summary.add_sample(sample);
        }

        double eps = 0.005;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(expected_summary.mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(expected_summary.variance()).epsilon(eps * 2));

        double freq_01 = count_01 / (double)nsamples;
        double freq_02 = count_02 / (double)nsamples;
        double freq_12 = count_12 / (double)nsamples;
        std::cout << "freq of ((0,1),2): " << freq_01 << "\n";
        std::cout << "freq of ((0,2),1): " << freq_02 << "\n";
        std::cout << "freq of ((1,2),0): " << freq_12 << "\n";
        REQUIRE(freq_01 == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(freq_02 == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(freq_12 == Approx(1.0/3.0).epsilon(eps));
    }
}

TEST_CASE("Testing NodeHeightSlideBumpPermuteMover with 3 leaves, gamma root, and optimizing",
        "[NodeHeightSlideBumpPermuteMover]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(16);

        double root_height_shape = 100.0;
        double root_height_scale = 0.01;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightSlideBumpPermuteMover<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_12 = 0;

        unsigned int niterations = 400000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        /* for (unsigned int i = 0; i < 100000; ++i) { */
        /*     op.operate(rng, &tree, 1); */
        /* } */
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0));
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(1) == tree.get_root_height());
                if (tree.get_root_ptr()->is_child("leaf0")) {
                    ++count_12;
                }
                if (tree.get_root_ptr()->is_child("leaf1")) {
                    ++count_02;
                }
                if (tree.get_root_ptr()->is_child("leaf2")) {
                    ++count_01;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        SampleSummarizer<double> expected_summary;
        for (unsigned int i = 0; i < niterations; ++i) {
            double root_ht = root_height_prior->draw(rng);
            double sample = rng.uniform_real(0.0, root_ht);
            expected_summary.add_sample(sample);
        }

        double eps = 0.005;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(expected_summary.mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(expected_summary.variance()).epsilon(eps * 2));

        double freq_01 = count_01 / (double)nsamples;
        double freq_02 = count_02 / (double)nsamples;
        double freq_12 = count_12 / (double)nsamples;
        std::cout << "freq of ((0,1),2): " << freq_01 << "\n";
        std::cout << "freq of ((0,2),1): " << freq_02 << "\n";
        std::cout << "freq of ((1,2),0): " << freq_12 << "\n";
        REQUIRE(freq_01 == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(freq_02 == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(freq_12 == Approx(1.0/3.0).epsilon(eps));
    }
}


TEST_CASE("Testing NodeHeightSlideBumpMover with 2 nested internals, fixed root, and optimizing",
        "[NodeHeightSlideBumpMover]") {

    SECTION("Testing 2 nested internals with fixed root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(17);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.7);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.3);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        root->add_child(internal1);
        root->add_child(leaf3);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        NodeHeightSlideBumpMover<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal0_height_summary;
        SampleSummarizer<double> internal1_height_summary;

        unsigned int niterations = 200000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal0_height_summary.add_sample(tree.get_node("internal0")->get_height());
                internal1_height_summary.add_sample(tree.get_node("internal1")->get_height());
                REQUIRE(tree.get_height(2) == tree.get_root_height());
                REQUIRE(tree.get_height(2) == 1.0);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(internal0_height_summary.sample_size() == nsamples);
        REQUIRE(internal1_height_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        std::cout << "internal0 mean: " << internal0_height_summary.mean() << "\n";
        std::cout << "internal1 mean: " << internal1_height_summary.mean() << "\n";

        // Because the nodes are nested they both cannot sample from the
        // uniform prior
        REQUIRE(internal0_height_summary.mean() < prior.get_mean());
        REQUIRE(internal0_height_summary.variance() < prior.get_variance());
        REQUIRE(internal1_height_summary.mean() > prior.get_mean());
        REQUIRE(internal1_height_summary.variance() < prior.get_variance());
    }
}
