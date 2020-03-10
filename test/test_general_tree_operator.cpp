#include "catch.hpp"
#include "ecoevolity/general_tree_operator.hpp"

#include <limits>
#include <memory>

#include "ecoevolity/tree.hpp"
#include "ecoevolity/node.hpp"
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

        unsigned int niterations = 600000;
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

        double root_height_lower = 0.26;
        double root_height_upper = 0.34;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<UniformDistribution>(
                root_height_lower,
                root_height_upper);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.3);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.15);
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

        RootHeightScaler<Node> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 200000;
        unsigned int burnin = 100;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < burnin; ++i) {
            rop.operate(rng, &tree, 1);
            op.operate(rng, &tree, 1);
        }
        for (unsigned int i = 0; i < niterations; ++i) {
            rop.operate(rng, &tree, 1);
            op.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                root_height_summary.add_sample(tree.get_root_height());
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_root_height());
                REQUIRE(tree.get_height(1) == tree.get_root_height());
            }
        }
        std::cout << rop.header_string();
        std::cout << rop.to_string();
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));
        REQUIRE(rop.get_number_of_attempts() == niterations + burnin);
        REQUIRE(rop.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing NodeHeightSlideBumpScaler with 3 leaves, gamma root, and optimizing",
        "[NodeHeightSlideBumpScaler]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(4);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.25);
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

        RootHeightScaler<Node> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 400000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            rop.operate(rng, &tree, 1);
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_root_height());
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(1) == tree.get_root_height());
            }
        }
        std::cout << rop.header_string();
        std::cout << rop.to_string();
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(rop.get_number_of_attempts() == niterations);
        REQUIRE(rop.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing NodeHeightSlideBumpScaler with 3 leaves, beta(4, 2), fixed root, and optimizing",
        "[NodeHeightSlideBumpScaler]") {

    SECTION("Testing 3 leaves with fixed root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(1);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.3);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.2);
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
        
        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(4.0);
        tree.set_beta_of_node_height_beta_prior(2.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

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
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                REQUIRE(tree.get_height(1) == tree.get_root_height());
                REQUIRE(tree.get_height(1) == 0.3);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        BetaDistribution prior(4.0, 2.0);

        double eps = 0.001;
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing NodeHeightSlideBumpScaler with 4 leaves, fixed root, and optimizing",
        "[NodeHeightSlideBumpScaler]") {

    SECTION("Testing 4 leaves with fixed root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(143132134);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.25);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(internal0);
        root->add_child(internal1);
        root->add_child(leaf3);

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

        SampleSummarizer<double> internal_height_0_summary;
        SampleSummarizer<double> internal_height_1_summary;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_0_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                internal_height_1_summary.add_sample(tree.get_height(1) / tree.get_height(2));
                REQUIRE(tree.get_height(2) == tree.get_root_height());
                REQUIRE(tree.get_height(2) == 1.0);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(internal_height_0_summary.sample_size() == nsamples);
        REQUIRE(internal_height_1_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.001;
        REQUIRE(internal_height_0_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_0_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_height_1_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_1_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}


TEST_CASE("Testing NodeHeightSlideBumpScaler with 4 leaves, balanced, fixed root, and optimizing",
        "[NodeHeightSlideBumpScaler]") {

    SECTION("Testing balanced 4 leaves with fixed root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(54709754);
        /* RandomNumberGenerator rng = RandomNumberGenerator(111); */

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.25);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.15);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);

        internal1->add_child(leaf0);
        internal1->add_child(leaf1);
        internal0->add_child(leaf2);
        internal0->add_child(leaf3);
        root->add_child(internal0);
        root->add_child(internal1);

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

        SampleSummarizer<double> internal_height_0_summary;
        SampleSummarizer<double> internal_height_1_summary;

        unsigned int niterations = 800000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;

        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_0_summary.add_sample(tree.get_node("internal0")->get_height());
                internal_height_1_summary.add_sample(tree.get_node("internal1")->get_height());
                REQUIRE(tree.get_height(2) == tree.get_root_height());
                REQUIRE(tree.get_height(2) == 0.4);
                REQUIRE(tree.get_height_of_youngest_parent(0) == 0.4);
                REQUIRE(tree.get_height_of_youngest_parent(1) == 0.4);
                REQUIRE(tree.get_height(0) < tree.get_height(1));
                REQUIRE(tree.get_node("internal0")->get_height() != tree.get_node("internal1")->get_height());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(internal_height_0_summary.sample_size() == nsamples);
        REQUIRE(internal_height_1_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 0.4);

        double eps = 0.001;
        std::cout << "int0 mean: " << internal_height_0_summary.mean() << "\n";
        std::cout << "int1 mean: " << internal_height_1_summary.mean() << "\n";
        REQUIRE(internal_height_0_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_0_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_height_1_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_1_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}


TEST_CASE("Testing NodeHeightSlideBumpScaler with 4 leaves, gamma root, and optimizing",
        "[NodeHeightSlideBumpScaler]") {

    SECTION("Testing 4 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(250925098743);

        double root_height_shape = 5.0;
        double root_height_scale = 0.1;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.125);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.25);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(internal0);
        root->add_child(internal1);
        root->add_child(leaf3);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.estimate_root_height();

        RootHeightScaler<Node> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        NodeHeightSlideBumpScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);
        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_0_summary;
        SampleSummarizer<double> internal_height_1_summary;

        unsigned int niterations = 1400000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            rop.operate(rng, &tree, 1);
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_0_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                internal_height_1_summary.add_sample(tree.get_height(1) / tree.get_height(2));
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(2) == tree.get_root_height());
                /* std::cout << i + 1 << "\t" << tree.get_root_height() << "\n"; */
            }
        }
        std::cout << rop.header_string();
        std::cout << rop.to_string();
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(rop.get_number_of_attempts() == niterations);
        REQUIRE(rop.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_0_summary.sample_size() == nsamples);
        REQUIRE(internal_height_1_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_0_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_0_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_height_1_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_1_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}


TEST_CASE("Testing NodeHeightSlideBumpScaler with 4 leaves, balanced, gamma root, and optimizing",
        "[NodeHeightSlideBumpScaler]") {

    SECTION("Testing balanced 4 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(4723445);

        double root_height_shape = 10.0;
        double root_height_scale = 0.08;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.125);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.25);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(leaf3);
        root->add_child(internal1);
        root->add_child(internal0);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.estimate_root_height();

        RootHeightScaler<Node> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        NodeHeightSlideBumpScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);
        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_0_summary;
        SampleSummarizer<double> internal_height_1_summary;

        unsigned int niterations = 1200000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            rop.operate(rng, &tree, 1);
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_0_summary.add_sample(tree.get_node("internal0")->get_height() / tree.get_root_height());
                internal_height_1_summary.add_sample(tree.get_node("internal1")->get_height() / tree.get_root_height());
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(2) == tree.get_root_height());
                REQUIRE(tree.get_height_of_youngest_parent(0) == tree.get_root_height());
                REQUIRE(tree.get_height_of_youngest_parent(1) == tree.get_root_height());
                REQUIRE(tree.get_height(0) < tree.get_height(1));
                REQUIRE(tree.get_node("internal0")->get_height() != tree.get_node("internal1")->get_height());
            }
        }
        std::cout << rop.header_string();
        std::cout << rop.to_string();
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(rop.get_number_of_attempts() == niterations);
        REQUIRE(rop.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_0_summary.sample_size() == nsamples);
        REQUIRE(internal_height_1_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.001;
        std::cout << "int0 mean: " << internal_height_0_summary.mean() << "\n";
        std::cout << "int1 mean: " << internal_height_1_summary.mean() << "\n";
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_0_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_0_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_height_1_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_1_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
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

        double root_height_lower = 0.48;
        double root_height_upper = 0.52;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<UniformDistribution>(
                root_height_lower,
                root_height_upper);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.25);
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

        RootHeightScaler<Node> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 600000;
        unsigned int burnin = 1000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1);
            rop.operate(rng, &tree, 1);
        }
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            rop.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_root_height());
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(1) == tree.get_root_height());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << rop.header_string();
        std::cout << rop.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));
        REQUIRE(rop.get_number_of_attempts() == niterations + burnin);
        REQUIRE(rop.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing NodeHeightSlideBumpMover with 3 leaves, gamma root, and optimizing",
        "[NodeHeightSlideBumpMover]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(8);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.25);
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

        RootHeightScaler<Node> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 400000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            rop.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_root_height());
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(1) == tree.get_root_height());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << rop.header_string();
        std::cout << rop.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(rop.get_number_of_attempts() == niterations);
        REQUIRE(rop.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}



TEST_CASE("Testing NodeHeightSlideBumpSwapScaler with 3 leaves, gamma root, optimizing, no op root",
        "[NodeHeightSlideBumpSwapScaler]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(9);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.25);
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
        op.set_operate_on_root(false);

        RootHeightScaler<Node> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

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
            rop.operate(rng, &tree, 1);
        }
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            rop.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_root_height());
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
        std::cout << rop.header_string();
        std::cout << rop.to_string();

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.005;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));

        double freq_01 = count_01 / (double)nsamples;
        double freq_02 = count_02 / (double)nsamples;
        double freq_12 = count_12 / (double)nsamples;
        std::cout << "freq of ((0,1),2): " << freq_01 << "\n";
        std::cout << "freq of ((0,2),1): " << freq_02 << "\n";
        std::cout << "freq of ((1,2),0): " << freq_12 << "\n";
        REQUIRE(freq_01 == 1.0);
        REQUIRE(freq_02 == 0.0);
        REQUIRE(freq_12 == 0.0);
    }
}

TEST_CASE("Testing NodeHeightSlideBumpSwapScaler with 3 leaves, gamma root, optimizing, and op root",
        "[NodeHeightSlideBumpSwapScaler]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(9);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.25);
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
        op.set_operate_on_root(true);

        /* RootHeightScaler<Node> rop; */
        /* rop.turn_on_auto_optimize(); */
        /* rop.set_auto_optimize_delay(100); */

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        /* REQUIRE(rop.auto_optimizing()); */
        /* REQUIRE(rop.get_auto_optimize_delay() == 100); */

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_12 = 0;

        unsigned int niterations = 800000;
        unsigned int sample_freq = 20;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < 100000; ++i) {
            op.operate(rng, &tree, 1);
            /* rop.operate(rng, &tree, 1); */
        }
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            /* rop.operate(rng, &tree, 1); */
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_root_height());
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
        /* std::cout << rop.header_string(); */
        /* std::cout << rop.to_string(); */

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.005;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));

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

TEST_CASE("Testing NodeHeightSlideBumpSwapMover with 3 leaves, gamma root, optimizing, no root op",
        "[NodeHeightSlideBumpSwapMover]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(10);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.25);
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
        op.set_operate_on_root(false);

        RootHeightScaler<Node> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

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
            rop.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_root_height());
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
        std::cout << rop.header_string();
        std::cout << rop.to_string();

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.005;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));

        double freq_01 = count_01 / (double)nsamples;
        double freq_02 = count_02 / (double)nsamples;
        double freq_12 = count_12 / (double)nsamples;
        std::cout << "freq of ((0,1),2): " << freq_01 << "\n";
        std::cout << "freq of ((0,2),1): " << freq_02 << "\n";
        std::cout << "freq of ((1,2),0): " << freq_12 << "\n";
        REQUIRE(freq_01 == 1.0);
        REQUIRE(freq_02 == 0.0);
        REQUIRE(freq_12 == 0.0);
    }
}

TEST_CASE("Testing NodeHeightSlideBumpSwapMover with 3 leaves, gamma root, optimizing, and root op",
        "[NodeHeightSlideBumpSwapMover]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(10);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.25);
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
        op.set_operate_on_root(true);

        /* RootHeightScaler<Node> rop; */
        /* rop.turn_on_auto_optimize(); */
        /* rop.set_auto_optimize_delay(100); */

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        /* REQUIRE(rop.auto_optimizing()); */
        /* REQUIRE(rop.get_auto_optimize_delay() == 100); */

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
            /* rop.operate(rng, &tree, 1); */
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_root_height());
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
        /* std::cout << rop.header_string(); */
        /* std::cout << rop.to_string(); */

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.005;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));

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

TEST_CASE("Testing NodeHeightSlideBumpPermuteScaler with 3 leaves, gamma root, optimizing, no root op",
        "[NodeHeightSlideBumpPermuteScaler]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(15);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.25);
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
        op.set_operate_on_root(false);

        RootHeightScaler<Node> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

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
            rop.operate(rng, &tree, 1);
        }
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            rop.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_root_height());
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
        std::cout << rop.header_string();
        std::cout << rop.to_string();

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.005;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));

        double freq_01 = count_01 / (double)nsamples;
        double freq_02 = count_02 / (double)nsamples;
        double freq_12 = count_12 / (double)nsamples;
        std::cout << "freq of ((0,1),2): " << freq_01 << "\n";
        std::cout << "freq of ((0,2),1): " << freq_02 << "\n";
        std::cout << "freq of ((1,2),0): " << freq_12 << "\n";
        REQUIRE(freq_01 == 1.0);
        REQUIRE(freq_02 == 0.0);
        REQUIRE(freq_12 == 0.0);
    }
}

TEST_CASE("Testing NodeHeightSlideBumpPermuteScaler with 3 leaves, gamma root, optimizing, and root op",
        "[NodeHeightSlideBumpPermuteScaler]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(15);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.25);
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
        op.set_operate_on_root(true);

        /* RootHeightScaler<Node> rop; */
        /* rop.turn_on_auto_optimize(); */
        /* rop.set_auto_optimize_delay(100); */

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        /* REQUIRE(rop.auto_optimizing()); */
        /* REQUIRE(rop.get_auto_optimize_delay() == 100); */

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
            /* rop.operate(rng, &tree, 1); */
        }
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            /* rop.operate(rng, &tree, 1); */
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_root_height());
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
        /* std::cout << rop.header_string(); */
        /* std::cout << rop.to_string(); */

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.005;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));

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

TEST_CASE("Testing NodeHeightSlideBumpPermuteMover with 3 leaves, gamma root, optimizing, no root op",
        "[NodeHeightSlideBumpPermuteMover]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(16);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.25);
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
        op.set_operate_on_root(false);

        RootHeightScaler<Node> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

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
            rop.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_root_height());
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
        std::cout << rop.header_string();
        std::cout << rop.to_string();

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.005;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));

        double freq_01 = count_01 / (double)nsamples;
        double freq_02 = count_02 / (double)nsamples;
        double freq_12 = count_12 / (double)nsamples;
        std::cout << "freq of ((0,1),2): " << freq_01 << "\n";
        std::cout << "freq of ((0,2),1): " << freq_02 << "\n";
        std::cout << "freq of ((1,2),0): " << freq_12 << "\n";
        REQUIRE(freq_01 == 1.0);
        REQUIRE(freq_02 == 0.0);
        REQUIRE(freq_12 == 0.0);
    }
}

TEST_CASE("Testing NodeHeightSlideBumpPermuteMover with 3 leaves, gamma root, optimizing, and root op",
        "[NodeHeightSlideBumpPermuteMover]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(16);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.25);
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
        op.set_operate_on_root(true);

        /* RootHeightScaler<Node> rop; */
        /* rop.turn_on_auto_optimize(); */
        /* rop.set_auto_optimize_delay(100); */

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        /* REQUIRE(rop.auto_optimizing()); */
        /* REQUIRE(rop.get_auto_optimize_delay() == 100); */

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
            /* rop.operate(rng, &tree, 1); */
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_root_height());
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
        /* std::cout << rop.header_string(); */
        /* std::cout << rop.to_string(); */

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        double eps = 0.005;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));

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
        RandomNumberGenerator rng = RandomNumberGenerator(924735);

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

        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1);
        }

        unsigned int niterations = 1500000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal0_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                internal1_height_summary.add_sample(tree.get_height(1) / tree.get_height(2));
                REQUIRE(tree.get_height(2) == tree.get_root_height());
                REQUIRE(tree.get_height(2) == 1.0);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(internal0_height_summary.sample_size() == nsamples);
        REQUIRE(internal1_height_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);

        std::cout << "internal0 mean: " << internal0_height_summary.mean() << "\n";
        std::cout << "internal1 mean: " << internal1_height_summary.mean() << "\n";

        double eps = 0.001;
        REQUIRE(internal0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing HeightSizeSlideBumpMixer with 3 leaves, unconstrained sizes, and optimizing",
        "[HeightSizeSlideBumpMixer]") {

    SECTION("Testing 3 leaves, unconstrained sizes, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(7982729874);

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double pop_size_shape = 10.0;
        double pop_size_scale = 0.02;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.estimate_root_height();

        tree.set_population_size_prior(pop_size_prior);
        tree.set_root_node_height_prior(root_height_prior);

        HeightSizeSlideBumpMixer op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);
        RootHeightSizeMixer rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        PopSizeScaler op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);
        NodeHeightScaler<PopulationNode> op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 5);
        std::map<std::string, SampleSummarizer<double> > pop_size_summaries;

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        std::shared_ptr<PopulationNode> node;
        unsigned int nmoves_per_op = 1;
        unsigned int niterations = 500000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            rop.operate(rng, &tree, 1);
            op2.operate(rng, &tree, 1, nmoves_per_op);
            op3.operate(rng, &tree, 1, nmoves_per_op);
            if ((i + 1) % sample_freq == 0) {
                root_height_summary.add_sample(tree.get_root_height());
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                node = tree.get_node("internal0");
                pop_size_summaries["internal0"].add_sample(node->get_population_size());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");

                node = tree.get_node("root");
                pop_size_summaries["root"].add_sample(node->get_population_size());

                node = tree.get_node("leaf0");
                pop_size_summaries["leaf0"].add_sample(node->get_population_size());
                node = tree.get_node("leaf1");
                pop_size_summaries["leaf1"].add_sample(node->get_population_size());
                node = tree.get_node("leaf2");
                pop_size_summaries["leaf2"].add_sample(node->get_population_size());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << rop.header_string();
        std::cout << rop.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();
        std::cout << op3.header_string();
        std::cout << op3.to_string();

        std::cout << "Expected pop size mean: " << pop_size_prior->get_mean() << "\n";
        std::cout << "Expected pop size var: " << pop_size_prior->get_variance() << "\n";
        for (auto k_v : pop_size_summaries) {
            std::string l = k_v.first;
            std::cout << "Pop size " << l << " mean: " << k_v.second.mean() << "\n";
            std::cout << "Pop size " << l << " var: " << k_v.second.variance() << "\n";
        }

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(rop.get_number_of_attempts() == niterations);
        REQUIRE(rop.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(op2.get_number_of_attempts() == (niterations * nmoves_per_op));
        REQUIRE(op2.get_number_of_attempts_for_correction() == ((nmoves_per_op * niterations) - 100));
        REQUIRE(op3.get_number_of_attempts() == (niterations * nmoves_per_op));
        REQUIRE(op3.get_number_of_attempts_for_correction() == ((nmoves_per_op * niterations) - 100));

        double eps = 0.001;

        BetaDistribution prior(1.0, 1.0);

        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        for (auto k_v : pop_size_summaries) {
            REQUIRE(k_v.second.sample_size() == nsamples);
            REQUIRE(k_v.second.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
            REQUIRE(k_v.second.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        }
    }
}


TEST_CASE("Testing NodeHeightScaler with 4 leaves, fixed root, and optimizing",
        "[NodeHeightScaler]") {

    SECTION("Testing 4 leaves with fixed root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(9475297498);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>(6, "root", 0.5);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(5, "internal1", 0.3);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(4, "internal0", 0.1);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        root->add_child(internal1);
        root->add_child(leaf3);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.fix_root_height();

        NodeHeightScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal0_height_summary;
        SampleSummarizer<double> internal1_height_summary;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                root_height_summary.add_sample(tree.get_root_height());
                internal0_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                internal1_height_summary.add_sample(tree.get_height(1) / tree.get_height(2));
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal0_height_summary.sample_size() == nsamples);
        REQUIRE(internal1_height_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);
        
        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(0.5));
        REQUIRE(root_height_summary.variance() == Approx(0.0));
        REQUIRE(internal0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}


TEST_CASE("Testing RootHeightScaler with 2 leaves, gamma root, and optimizing",
        "[RootHeightScaler]") {

    SECTION("Testing 2 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(54168);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>(2, "root", 0.5);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);

        root->add_child(leaf0);
        root->add_child(leaf1);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.estimate_root_height();

        RootHeightScaler<Node> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;

        unsigned int niterations = 200000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            rop.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                root_height_summary.add_sample(tree.get_root_height());
            }
        }
        std::cout << rop.header_string();
        std::cout << rop.to_string();

        REQUIRE(rop.get_number_of_attempts() == niterations);
        REQUIRE(rop.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);
        
        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing RootHeightScaler with 3 leaves, gamma root, and optimizing",
        "[RootHeightScaler]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(54168);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.001);
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

        RootHeightScaler<Node> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        NodeHeightScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 1400000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            rop.operate(rng, &tree, 1);
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(1) == tree.get_root_height());
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
            }
        }
        std::cout << rop.header_string();
        std::cout << rop.to_string();
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(rop.get_number_of_attempts() == niterations);
        REQUIRE(rop.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        UniformDistribution prior(0.0, 1.0);
        
        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}


TEST_CASE("Testing RootHeightScaler with 4 leaves, gamma root, and optimizing",
        "[RootHeightScaler]") {

    SECTION("Testing 4 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(213579);

        double root_height_shape = 5.0;
        double root_height_scale = 0.1;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.125);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.25);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(internal0);
        root->add_child(internal1);
        root->add_child(leaf3);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.estimate_root_height();

        RootHeightScaler<Node> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        NodeHeightScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            rop.operate(rng, &tree, 1);
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(2) == tree.get_root_height());
                internal_0_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                internal_1_height_summary.add_sample(tree.get_height(1) / tree.get_height(2));
            }
        }
        std::cout << rop.header_string();
        std::cout << rop.to_string();
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(rop.get_number_of_attempts() == niterations);
        REQUIRE(rop.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);

        UniformDistribution prior(0.0, 1.0);

        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
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



TEST_CASE("Testing NodeHeightPriorAlphaScaler with optimizing",
        "[NodeHeightPriorAlphaOperator]") {

    SECTION("Testing NodeHeightPriorAlphaScaler with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(274589347597);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.4);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.2);
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

        tree.estimate_alpha_of_node_height_beta_prior();
        std::shared_ptr<ContinuousProbabilityDistribution> alpha_prior = std::make_shared<GammaDistribution>(
                20.0,
                0.1);
        tree.set_prior_on_alpha_of_node_height_beta_prior(alpha_prior);

        NodeHeightPriorAlphaScaler<Node> op;
        op.set_coercable_parameter_value(1.0);
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        NodeHeightScaler<Node> op2;
        op2.set_coercable_parameter_value(1.0);
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(op2.auto_optimizing());
        REQUIRE(op2.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> alpha_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 600000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            op2.operate(rng, &tree, 1, 2);
            if ((i + 1) % sample_freq == 0) {
                alpha_summary.add_sample(tree.get_alpha_of_node_height_beta_prior());
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(alpha_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        double eps = 0.001;
        REQUIRE(alpha_summary.mean() == Approx(alpha_prior->get_mean()).epsilon(eps));
        REQUIRE(alpha_summary.variance() == Approx(alpha_prior->get_variance()).epsilon(eps * 2));
        REQUIRE(internal_height_summary.mean() > 0.65);
        REQUIRE(internal_height_summary.mean() < 0.68);
    }
}

TEST_CASE("Testing NodeHeightPriorAlphaMover with optimizing",
        "[NodeHeightPriorAlphaOperator]") {

    SECTION("Testing NodeHeightPriorAlphaMover with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(583648364);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.4);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.2);
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

        tree.estimate_alpha_of_node_height_beta_prior();
        std::shared_ptr<ContinuousProbabilityDistribution> alpha_prior = std::make_shared<GammaDistribution>(
                20.0,
                0.1);
        tree.set_prior_on_alpha_of_node_height_beta_prior(alpha_prior);

        NodeHeightPriorAlphaMover<Node> op;
        op.set_coercable_parameter_value(1.0);
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        NodeHeightScaler<Node> op2;
        op2.set_coercable_parameter_value(1.0);
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(op2.auto_optimizing());
        REQUIRE(op2.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> alpha_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 600000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            op2.operate(rng, &tree, 1, 2);
            if ((i + 1) % sample_freq == 0) {
                alpha_summary.add_sample(tree.get_alpha_of_node_height_beta_prior());
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(alpha_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        double eps = 0.001;
        REQUIRE(alpha_summary.mean() == Approx(alpha_prior->get_mean()).epsilon(eps));
        REQUIRE(alpha_summary.variance() == Approx(alpha_prior->get_variance()).epsilon(eps * 2));
        REQUIRE(internal_height_summary.mean() > 0.65);
        REQUIRE(internal_height_summary.mean() < 0.68);
    }
}


TEST_CASE("Testing NodeHeightPriorBetaScaler with optimizing",
        "[NodeHeightPriorBetaOperator]") {

    SECTION("Testing NodeHeightPriorBetaScaler with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(274589347597);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.4);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.2);
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

        tree.estimate_beta_of_node_height_beta_prior();
        std::shared_ptr<ContinuousProbabilityDistribution> beta_prior = std::make_shared<GammaDistribution>(
                20.0,
                0.1);
        tree.set_prior_on_beta_of_node_height_beta_prior(beta_prior);

        NodeHeightPriorBetaScaler<Node> op;
        op.set_coercable_parameter_value(1.0);
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        NodeHeightScaler<Node> op2;
        op2.set_coercable_parameter_value(1.0);
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(op2.auto_optimizing());
        REQUIRE(op2.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> beta_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 600000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            op2.operate(rng, &tree, 1, 2);
            if ((i + 1) % sample_freq == 0) {
                beta_summary.add_sample(tree.get_beta_of_node_height_beta_prior());
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(beta_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        double eps = 0.001;
        REQUIRE(beta_summary.mean() == Approx(beta_prior->get_mean()).epsilon(eps));
        REQUIRE(beta_summary.variance() == Approx(beta_prior->get_variance()).epsilon(eps * 2));
        REQUIRE(internal_height_summary.mean() > 0.32);
        REQUIRE(internal_height_summary.mean() < 0.35);
    }
}

TEST_CASE("Testing NodeHeightPriorBetaMover with optimizing",
        "[NodeHeightPriorBetaOperator]") {

    SECTION("Testing NodeHeightPriorBetaMover with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(583648364);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.4);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.2);
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

        tree.estimate_beta_of_node_height_beta_prior();
        std::shared_ptr<ContinuousProbabilityDistribution> beta_prior = std::make_shared<GammaDistribution>(
                20.0,
                0.1);
        tree.set_prior_on_beta_of_node_height_beta_prior(beta_prior);

        NodeHeightPriorBetaMover<Node> op;
        op.set_coercable_parameter_value(1.0);
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        NodeHeightScaler<Node> op2;
        op2.set_coercable_parameter_value(1.0);
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(op2.auto_optimizing());
        REQUIRE(op2.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> beta_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 600000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            op2.operate(rng, &tree, 1, 2);
            if ((i + 1) % sample_freq == 0) {
                beta_summary.add_sample(tree.get_beta_of_node_height_beta_prior());
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(beta_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        double eps = 0.001;
        REQUIRE(beta_summary.mean() == Approx(beta_prior->get_mean()).epsilon(eps));
        REQUIRE(beta_summary.variance() == Approx(beta_prior->get_variance()).epsilon(eps * 2));
        REQUIRE(internal_height_summary.mean() > 0.32);
        REQUIRE(internal_height_summary.mean() < 0.35);
    }
}

TEST_CASE("Testing TreeScaler with 3 leaves, gamma root, and optimizing",
        "[TreeScaler]") {

    SECTION("Testing 3 leaves with variable root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(4);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.25);
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

        TreeScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 400000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_root_height());
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
        
        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing TreeScaler with 3 leaves, gamma root, internal free, and optimizing",
        "[TreeScaler]") {

    SECTION("Testing 3 leaves with variable root, internal free, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(1422315338);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.25);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(3.0);
        tree.set_beta_of_node_height_beta_prior(1.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.estimate_root_height();

        TreeScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        NodeHeightScaler<Node> op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(op2.auto_optimizing());
        REQUIRE(op2.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int niterations = 400000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            op2.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_root_height());
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(1) == tree.get_root_height());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(op2.get_number_of_attempts() == niterations);
        REQUIRE(op2.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        BetaDistribution prior(3.0, 1.0);
        
        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing TreeScaler with 4 leaves, gamma root, internals fixed, and optimizing",
        "[TreeScaler]") {

    SECTION("Testing 4 leaves with variable root, internals fixed, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(4);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.15);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.25);
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

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(3.0);
        tree.set_beta_of_node_height_beta_prior(1.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.estimate_root_height();

        TreeScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal0_height_summary;
        SampleSummarizer<double> internal1_height_summary;

        unsigned int niterations = 400000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal0_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                internal1_height_summary.add_sample(tree.get_height(1) / tree.get_height(2));
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(2) == tree.get_root_height());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal0_height_summary.sample_size() == nsamples);
        REQUIRE(internal1_height_summary.sample_size() == nsamples);
        
        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing TreeScaler with 4 leaves, gamma root, internals free, and optimizing",
        "[TreeScaler]") {

    SECTION("Testing 4 leaves with variable root, internal frees, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(4);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.15);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.25);
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

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(3.0);
        tree.set_beta_of_node_height_beta_prior(1.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.set_root_node_height_prior(root_height_prior);

        tree.ignore_data();
        tree.estimate_root_height();

        TreeScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        NodeHeightScaler<Node> op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(op2.auto_optimizing());
        REQUIRE(op2.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal0_height_summary;
        SampleSummarizer<double> internal1_height_summary;

        unsigned int niterations = 400000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            op2.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal0_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                internal1_height_summary.add_sample(tree.get_height(1) / tree.get_height(2));
                root_height_summary.add_sample(tree.get_root_height());
                REQUIRE(tree.get_height(2) == tree.get_root_height());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(op2.get_number_of_attempts() == niterations);
        REQUIRE(op2.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal0_height_summary.sample_size() == nsamples);
        REQUIRE(internal1_height_summary.sample_size() == nsamples);
        
        BetaDistribution prior(3.0, 1.0);
        
        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing GlobalPopSizeScaler with 3 leaves, constrained sizes, and optimizing",
        "[GlobalPopSizeScaler]") {

    SECTION("Testing 3 leaves, constrained sizes, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(4);

        double pop_size_shape = 10.0;
        double pop_size_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        tree.set_population_size_prior(pop_size_prior);
        tree.constrain_population_sizes();

        GlobalPopSizeScaler op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> pop_size_summary;

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes;

        unsigned int niterations = 400000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                pop_sizes = tree.get_pointers_to_population_sizes();
                REQUIRE(pop_sizes.size() == 1);
                pop_size_summary.add_sample(pop_sizes.at(0)->get_value());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(pop_size_summary.sample_size() == nsamples);
        
        double eps = 0.001;
        REQUIRE(pop_size_summary.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
        REQUIRE(pop_size_summary.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing PopSizeScaler with 3 leaves, constrained sizes, and optimizing",
        "[PopSizeScaler]") {

    SECTION("Testing 3 leaves, constrained sizes, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(4);

        double pop_size_shape = 10.0;
        double pop_size_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        tree.set_population_size_prior(pop_size_prior);
        tree.constrain_population_sizes();

        PopSizeScaler op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> pop_size_summary;

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes;

        unsigned int niterations = 400000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                pop_sizes = tree.get_pointers_to_population_sizes();
                REQUIRE(pop_sizes.size() == 1);
                pop_size_summary.add_sample(pop_sizes.at(0)->get_value());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(pop_size_summary.sample_size() == nsamples);
        
        double eps = 0.001;
        REQUIRE(pop_size_summary.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
        REQUIRE(pop_size_summary.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing PopSizeScaler with 3 leaves, unconstrained sizes, and optimizing",
        "[PopSizeScaler]") {

    SECTION("Testing 3 leaves, unconstrained sizes, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(4);

        double pop_size_shape = 10.0;
        double pop_size_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        tree.set_population_size_prior(pop_size_prior);

        PopSizeScaler op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 5);
        std::vector<SampleSummarizer<double> > pop_size_summaries(pop_sizes.size());


        unsigned int niterations = 400000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 5);
            if ((i + 1) % sample_freq == 0) {
                pop_sizes = tree.get_pointers_to_population_sizes();
                REQUIRE(pop_sizes.size() == 5);
                for (unsigned int i = 0; i < pop_sizes.size(); ++i) {
                    pop_size_summaries.at(i).add_sample(pop_sizes.at(i)->get_value());
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == (niterations * 5));
        REQUIRE(op.get_number_of_attempts_for_correction() == ((niterations * 5) - 100));

        
        double eps = 0.001;
        for (unsigned int i = 0; i < pop_size_summaries.size(); ++i) {
            REQUIRE(pop_size_summaries.at(i).sample_size() == nsamples);
            REQUIRE(pop_size_summaries.at(i).mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
            REQUIRE(pop_size_summaries.at(i).variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        }
    }
}

TEST_CASE("Testing GlobalPopSizeScaler with 3 leaves, unconstrained sizes, and optimizing",
        "[GlobalPopSizeScaler]") {

    SECTION("Testing 3 leaves, constrained sizes, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(4);

        double pop_size_shape = 10.0;
        double pop_size_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        tree.set_population_size_prior(pop_size_prior);

        GlobalPopSizeScaler op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        PopSizeScaler op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(op2.auto_optimizing());
        REQUIRE(op2.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 5);
        std::vector<SampleSummarizer<double> > pop_size_summaries(pop_sizes.size());

        unsigned int nmoves_per_op = 3;
        unsigned int niterations = 200000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            op2.operate(rng, &tree, 1, nmoves_per_op);
            if ((i + 1) % sample_freq == 0) {
                pop_sizes = tree.get_pointers_to_population_sizes();
                REQUIRE(pop_sizes.size() == 5);
                for (unsigned int i = 0; i < pop_sizes.size(); ++i) {
                    pop_size_summaries.at(i).add_sample(pop_sizes.at(i)->get_value());
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(op2.get_number_of_attempts() == (niterations * nmoves_per_op));
        REQUIRE(op2.get_number_of_attempts_for_correction() == ((nmoves_per_op * niterations) - 100));

        double eps = 0.001;
        for (unsigned int i = 0; i < pop_size_summaries.size(); ++i) {
            REQUIRE(pop_size_summaries.at(i).sample_size() == nsamples);
            REQUIRE(pop_size_summaries.at(i).mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
            REQUIRE(pop_size_summaries.at(i).variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        }
    }
}


TEST_CASE("Testing MuRateScaler with 3 leaves, constrained sizes, and optimizing",
        "[MuRateScaler]") {

    SECTION("Testing 3 leaves, constrained sizes, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(4);

        double rate_shape = 10.0;
        double rate_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> rate_prior = std::make_shared<GammaDistribution>(
                rate_shape,
                rate_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.fix_root_height();
        tree.fix_population_sizes();

        tree.set_mutation_rate_prior(rate_prior);
        tree.estimate_mutation_rate();

        MuRateScaler op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> rate_summary;

        unsigned int niterations = 400000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                rate_summary.add_sample(tree.get_mutation_rate());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(rate_summary.sample_size() == nsamples);
        
        double eps = 0.001;
        REQUIRE(rate_summary.mean() == Approx(rate_prior->get_mean()).epsilon(eps));
        REQUIRE(rate_summary.variance() == Approx(rate_prior->get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing GlobalHeightSizeMixer with 3 leaves, unconstrained sizes, and optimizing",
        "[GlobalHeightSizeMixer]") {

    SECTION("Testing 3 leaves, unconstrained sizes, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(9187243457);

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);
        double pop_size_shape = 10.0;
        double pop_size_scale = 0.02;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.estimate_root_height();

        tree.set_population_size_prior(pop_size_prior);
        tree.set_root_node_height_prior(root_height_prior);

        GlobalHeightSizeMixer op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        PopSizeScaler op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);
        NodeHeightScaler<PopulationNode> op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 5);
        std::vector<SampleSummarizer<double> > pop_size_summaries(pop_sizes.size());

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int nmoves_per_op = 2;
        unsigned int niterations = 400000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            op2.operate(rng, &tree, 1, nmoves_per_op);
            op3.operate(rng, &tree, 1, nmoves_per_op);
            if ((i + 1) % sample_freq == 0) {
                pop_sizes = tree.get_pointers_to_population_sizes();
                REQUIRE(pop_sizes.size() == 5);
                for (unsigned int i = 0; i < pop_sizes.size(); ++i) {
                    pop_size_summaries.at(i).add_sample(pop_sizes.at(i)->get_value());
                }
                root_height_summary.add_sample(tree.get_root_height());
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();
        std::cout << op3.header_string();
        std::cout << op3.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(op2.get_number_of_attempts() == (niterations * nmoves_per_op));
        REQUIRE(op2.get_number_of_attempts_for_correction() == ((nmoves_per_op * niterations) - 100));
        REQUIRE(op3.get_number_of_attempts() == (niterations * nmoves_per_op));
        REQUIRE(op3.get_number_of_attempts_for_correction() == ((nmoves_per_op * niterations) - 100));

        double eps = 0.001;
        for (unsigned int i = 0; i < pop_size_summaries.size(); ++i) {
            REQUIRE(pop_size_summaries.at(i).sample_size() == nsamples);
            REQUIRE(pop_size_summaries.at(i).mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
            REQUIRE(pop_size_summaries.at(i).variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        }

        BetaDistribution prior(1.0, 1.0);

        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing GlobalHeightSizeMixer with 3 leaves, constrained sizes, and optimizing",
        "[GlobalHeightSizeMixer]") {

    SECTION("Testing 3 leaves, constrained sizes, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(9872357426);

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);
        double pop_size_shape = 10.0;
        double pop_size_scale = 0.02;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.estimate_root_height();

        tree.set_population_size_prior(pop_size_prior);
        tree.set_root_node_height_prior(root_height_prior);

        tree.constrain_population_sizes();

        GlobalHeightSizeMixer op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        PopSizeScaler op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);
        NodeHeightScaler<PopulationNode> op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 1);
        SampleSummarizer<double> pop_size_summary;

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int nmoves_per_op = 3;
        unsigned int niterations = 200000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            op2.operate(rng, &tree, 1, 1);
            op3.operate(rng, &tree, 1, nmoves_per_op);
            if ((i + 1) % sample_freq == 0) {
                pop_sizes = tree.get_pointers_to_population_sizes();
                REQUIRE(pop_sizes.size() == 1);
                pop_size_summary.add_sample(pop_sizes.at(0)->get_value());
                root_height_summary.add_sample(tree.get_root_height());
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();
        std::cout << op3.header_string();
        std::cout << op3.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(op2.get_number_of_attempts() == niterations);
        REQUIRE(op2.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(op3.get_number_of_attempts() == (niterations * nmoves_per_op));
        REQUIRE(op3.get_number_of_attempts_for_correction() == ((nmoves_per_op * niterations) - 100));

        double eps = 0.001;
        REQUIRE(pop_size_summary.sample_size() == nsamples);
        REQUIRE(pop_size_summary.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
        REQUIRE(pop_size_summary.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));

        BetaDistribution prior(1.0, 1.0);

        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing HeightSizeMixer with 3 leaves, unconstrained sizes, and optimizing",
        "[HeightSizeMixer]") {

    SECTION("Testing 3 leaves, unconstrained sizes, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(6416461);

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double pop_size_shape = 10.0;
        double pop_size_scale = 0.02;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.estimate_root_height();

        tree.set_population_size_prior(pop_size_prior);
        tree.set_root_node_height_prior(root_height_prior);

        HeightSizeMixer op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);
        RootHeightSizeMixer rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        PopSizeScaler op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);
        NodeHeightScaler<PopulationNode> op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 5);
        std::vector<SampleSummarizer<double> > pop_size_summaries(pop_sizes.size());

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;

        unsigned int nmoves_per_op = 2;
        unsigned int niterations = 200000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            rop.operate(rng, &tree, 1);
            op2.operate(rng, &tree, 1, nmoves_per_op);
            op3.operate(rng, &tree, 1, nmoves_per_op);
            if ((i + 1) % sample_freq == 0) {
                pop_sizes = tree.get_pointers_to_population_sizes();
                REQUIRE(pop_sizes.size() == 5);
                for (unsigned int i = 0; i < pop_sizes.size(); ++i) {
                    pop_size_summaries.at(i).add_sample(pop_sizes.at(i)->get_value());
                }
                root_height_summary.add_sample(tree.get_root_height());
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << rop.header_string();
        std::cout << rop.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();
        std::cout << op3.header_string();
        std::cout << op3.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(rop.get_number_of_attempts() == niterations);
        REQUIRE(rop.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(op2.get_number_of_attempts() == (niterations * nmoves_per_op));
        REQUIRE(op2.get_number_of_attempts_for_correction() == ((nmoves_per_op * niterations) - 100));
        REQUIRE(op3.get_number_of_attempts() == (niterations * nmoves_per_op));
        REQUIRE(op3.get_number_of_attempts_for_correction() == ((nmoves_per_op * niterations) - 100));

        double eps = 0.001;
        for (unsigned int i = 0; i < pop_size_summaries.size(); ++i) {
            REQUIRE(pop_size_summaries.at(i).sample_size() == nsamples);
            REQUIRE(pop_size_summaries.at(i).mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
            REQUIRE(pop_size_summaries.at(i).variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        }

        BetaDistribution prior(1.0, 1.0);

        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing GlobalHeightSizeRateScaler with 3 leaves, unconstrained sizes, and optimizing",
        "[GlobalHeightSizeRateScaler]") {

    SECTION("Testing 3 leaves, unconstrained sizes, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(51464654);

        double mu_rate_shape = 10.0;
        double mu_rate_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> mu_rate_prior = std::make_shared<GammaDistribution>(
                mu_rate_shape,
                mu_rate_scale);

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double pop_size_shape = 10.0;
        double pop_size_scale = 0.02;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.estimate_root_height();
        tree.estimate_mutation_rate();

        tree.set_population_size_prior(pop_size_prior);
        tree.set_root_node_height_prior(root_height_prior);
        tree.set_mutation_rate_prior(mu_rate_prior);

        GlobalHeightSizeRateScaler op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);
        RootHeightScaler<PopulationNode> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        PopSizeScaler op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);
        NodeHeightScaler<PopulationNode> op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);
        MuRateScaler op4;
        op4.turn_on_auto_optimize();
        op4.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 5);
        std::map<std::string, SampleSummarizer<double> > pop_size_summaries;

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;
        SampleSummarizer<double> mu_rate_summary;

        unsigned int nmoves_per_op = 3;
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
            rop.operate(rng, &tree, 1, 1);
            op2.operate(rng, &tree, 1, nmoves_per_op);
            op3.operate(rng, &tree, 1, nmoves_per_op);
            op4.operate(rng, &tree, 1, 1);
        }

        std::shared_ptr<PopulationNode> node;
        unsigned int niterations = 500000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            rop.operate(rng, &tree, 1, 1);
            op2.operate(rng, &tree, 1, nmoves_per_op);
            op3.operate(rng, &tree, 1, nmoves_per_op);
            op4.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                root_height_summary.add_sample(tree.get_root_height());
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                mu_rate_summary.add_sample(tree.get_mutation_rate());
                node = tree.get_node("root");
                pop_size_summaries["root"].add_sample(node->get_population_size());
                node = tree.get_node("internal0");
                pop_size_summaries["internal0"].add_sample(node->get_population_size());
                node = tree.get_node("leaf0");
                pop_size_summaries["leaf0"].add_sample(node->get_population_size());
                node = tree.get_node("leaf1");
                pop_size_summaries["leaf1"].add_sample(node->get_population_size());
                node = tree.get_node("leaf2");
                pop_size_summaries["leaf2"].add_sample(node->get_population_size());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << rop.header_string();
        std::cout << rop.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();
        std::cout << op3.header_string();
        std::cout << op3.to_string();
        std::cout << op4.header_string();
        std::cout << op4.to_string();

        std::cout << "Expected pop size mean: " << pop_size_prior->get_mean() << "\n";
        std::cout << "Expected pop size var: " << pop_size_prior->get_variance() << "\n";
        for (auto k_v : pop_size_summaries) {
            std::string l = k_v.first;
            std::cout << "Pop size " << l << " mean: " << k_v.second.mean() << "\n";
            std::cout << "Pop size " << l << " var: " << k_v.second.variance() << "\n";
        }

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));
        REQUIRE(rop.get_number_of_attempts() == niterations + burnin);
        REQUIRE(rop.get_number_of_attempts_for_correction() == (niterations + burnin - 100));
        REQUIRE(op2.get_number_of_attempts() == ((niterations + burnin) * nmoves_per_op));
        REQUIRE(op2.get_number_of_attempts_for_correction() == ((nmoves_per_op * (burnin + niterations)) - 100));
        REQUIRE(op3.get_number_of_attempts() == ((niterations + burnin) * nmoves_per_op));
        REQUIRE(op3.get_number_of_attempts_for_correction() == ((nmoves_per_op * (burnin + niterations)) - 100));
        REQUIRE(op4.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op4.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        double eps = 0.001;

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(mu_rate_summary.sample_size() == nsamples);

        BetaDistribution prior(1.0, 1.0);

        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(mu_rate_summary.mean() == Approx(mu_rate_prior->get_mean()).epsilon(eps));
        REQUIRE(mu_rate_summary.variance() == Approx(mu_rate_prior->get_variance()).epsilon(eps));
        for (auto k_v : pop_size_summaries) {
            REQUIRE(k_v.second.sample_size() == nsamples);
            REQUIRE(k_v.second.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
            REQUIRE(k_v.second.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        }
    }
}

TEST_CASE("Testing GlobalHeightSizeRateScaler with 3 leaves, constrained sizes, and optimizing",
        "[GlobalHeightSizeRateScaler]") {

    SECTION("Testing 3 leaves, unconstrained sizes, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(472374912);

        double mu_rate_shape = 10.0;
        double mu_rate_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> mu_rate_prior = std::make_shared<GammaDistribution>(
                mu_rate_shape,
                mu_rate_scale);

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double pop_size_shape = 10.0;
        double pop_size_scale = 0.02;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.estimate_root_height();
        tree.estimate_mutation_rate();

        tree.set_population_size_prior(pop_size_prior);
        tree.set_root_node_height_prior(root_height_prior);
        tree.set_mutation_rate_prior(mu_rate_prior);
        tree.constrain_population_sizes();

        GlobalHeightSizeRateScaler op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);
        RootHeightScaler<PopulationNode> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        PopSizeScaler op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);
        NodeHeightScaler<PopulationNode> op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);
        MuRateScaler op4;
        op4.turn_on_auto_optimize();
        op4.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 1);
        std::vector<SampleSummarizer<double> > pop_size_summaries(pop_sizes.size());

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;
        SampleSummarizer<double> mu_rate_summary;

        unsigned int nmoves_per_op = 2;
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
            rop.operate(rng, &tree, 1, 1);
            op2.operate(rng, &tree, 1, nmoves_per_op);
            op3.operate(rng, &tree, 1, nmoves_per_op);
            op4.operate(rng, &tree, 1, 1);
        }

        unsigned int niterations = 600000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            rop.operate(rng, &tree, 1, 1);
            op2.operate(rng, &tree, 1, nmoves_per_op);
            op3.operate(rng, &tree, 1, nmoves_per_op);
            op4.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                pop_sizes = tree.get_pointers_to_population_sizes();
                REQUIRE(pop_sizes.size() == 1);
                for (unsigned int i = 0; i < pop_sizes.size(); ++i) {
                    pop_size_summaries.at(i).add_sample(pop_sizes.at(i)->get_value());
                }
                root_height_summary.add_sample(tree.get_root_height());
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                mu_rate_summary.add_sample(tree.get_mutation_rate());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << rop.header_string();
        std::cout << rop.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();
        std::cout << op3.header_string();
        std::cout << op3.to_string();
        std::cout << op4.header_string();
        std::cout << op4.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));
        REQUIRE(rop.get_number_of_attempts() == niterations + burnin);
        REQUIRE(rop.get_number_of_attempts_for_correction() == (niterations + burnin - 100));
        REQUIRE(op2.get_number_of_attempts() == ((niterations + burnin) * nmoves_per_op));
        REQUIRE(op2.get_number_of_attempts_for_correction() == ((nmoves_per_op * (burnin + niterations)) - 100));
        REQUIRE(op3.get_number_of_attempts() == ((niterations + burnin) * nmoves_per_op));
        REQUIRE(op3.get_number_of_attempts_for_correction() == ((nmoves_per_op * (burnin + niterations)) - 100));
        REQUIRE(op4.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op4.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        double eps = 0.001;
        for (unsigned int i = 0; i < pop_size_summaries.size(); ++i) {
            REQUIRE(pop_size_summaries.at(i).sample_size() == nsamples);
            REQUIRE(pop_size_summaries.at(i).mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
            REQUIRE(pop_size_summaries.at(i).variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        }

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(mu_rate_summary.sample_size() == nsamples);

        BetaDistribution prior(1.0, 1.0);

        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(mu_rate_summary.mean() == Approx(mu_rate_prior->get_mean()).epsilon(eps));
        REQUIRE(mu_rate_summary.variance() == Approx(mu_rate_prior->get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing GlobalHeightSizeScaler with 3 leaves, unconstrained sizes, and optimizing",
        "[GlobalHeightSizeScaler]") {

    SECTION("Testing 3 leaves, unconstrained sizes, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(98464610);

        double mu_rate_shape = 10.0;
        double mu_rate_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> mu_rate_prior = std::make_shared<GammaDistribution>(
                mu_rate_shape,
                mu_rate_scale);

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double pop_size_shape = 10.0;
        double pop_size_scale = 0.02;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.estimate_root_height();
        tree.estimate_mutation_rate();

        tree.set_population_size_prior(pop_size_prior);
        tree.set_root_node_height_prior(root_height_prior);
        tree.set_mutation_rate_prior(mu_rate_prior);

        GlobalHeightSizeScaler op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);
        RootHeightScaler<PopulationNode> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        PopSizeScaler op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);
        NodeHeightScaler<PopulationNode> op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);
        /* MuRateScaler op4; */
        /* op4.turn_on_auto_optimize(); */
        /* op4.set_auto_optimize_delay(100); */

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 5);
        std::vector<SampleSummarizer<double> > pop_size_summaries(pop_sizes.size());

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;
        SampleSummarizer<double> mu_rate_summary;

        unsigned int nmoves_per_op = 2;
        unsigned int niterations = 400000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            rop.operate(rng, &tree, 1, 1);
            op2.operate(rng, &tree, 1, nmoves_per_op);
            op3.operate(rng, &tree, 1, nmoves_per_op);
            /* op4.operate(rng, &tree, 1, 1); */
            if ((i + 1) % sample_freq == 0) {
                pop_sizes = tree.get_pointers_to_population_sizes();
                REQUIRE(pop_sizes.size() == 5);
                for (unsigned int i = 0; i < pop_sizes.size(); ++i) {
                    pop_size_summaries.at(i).add_sample(pop_sizes.at(i)->get_value());
                }
                root_height_summary.add_sample(tree.get_root_height());
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                mu_rate_summary.add_sample(tree.get_mutation_rate());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << rop.header_string();
        std::cout << rop.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();
        std::cout << op3.header_string();
        std::cout << op3.to_string();
        /* std::cout << op4.header_string(); */
        /* std::cout << op4.to_string(); */

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(rop.get_number_of_attempts() == niterations);
        REQUIRE(rop.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(op2.get_number_of_attempts() == (niterations * nmoves_per_op));
        REQUIRE(op2.get_number_of_attempts_for_correction() == ((nmoves_per_op * niterations) - 100));
        REQUIRE(op3.get_number_of_attempts() == (niterations * nmoves_per_op));
        REQUIRE(op3.get_number_of_attempts_for_correction() == ((nmoves_per_op * niterations) - 100));
        /* REQUIRE(op4.get_number_of_attempts() == niterations); */
        /* REQUIRE(op4.get_number_of_attempts_for_correction() == (niterations - 100)); */

        for (unsigned int i = 0; i < pop_size_summaries.size(); ++i) {
            std::cout << "pop size mean: " <<  pop_size_summaries.at(i).mean() << "\n";
            std::cout << "pop size var:  " <<  pop_size_summaries.at(i).variance() << "\n";
        }

        double eps = 0.001;
        for (unsigned int i = 0; i < pop_size_summaries.size(); ++i) {
            REQUIRE(pop_size_summaries.at(i).sample_size() == nsamples);
            REQUIRE(pop_size_summaries.at(i).mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
            REQUIRE(pop_size_summaries.at(i).variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        }

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(mu_rate_summary.sample_size() == nsamples);

        BetaDistribution prior(1.0, 1.0);

        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        /* REQUIRE(mu_rate_summary.mean() == Approx(mu_rate_prior->get_mean()).epsilon(eps)); */
        REQUIRE(mu_rate_summary.variance() == Approx(0.0).epsilon(eps));
    }
}

TEST_CASE("Testing GlobalHeightSizeScaler with 3 leaves, constrained sizes, and optimizing",
        "[GlobalHeightSizeScaler]") {

    SECTION("Testing 3 leaves, unconstrained sizes, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(98464610);

        double mu_rate_shape = 10.0;
        double mu_rate_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> mu_rate_prior = std::make_shared<GammaDistribution>(
                mu_rate_shape,
                mu_rate_scale);

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double pop_size_shape = 10.0;
        double pop_size_scale = 0.02;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.estimate_root_height();
        tree.estimate_mutation_rate();
        tree.constrain_population_sizes();

        tree.set_population_size_prior(pop_size_prior);
        tree.set_root_node_height_prior(root_height_prior);
        tree.set_mutation_rate_prior(mu_rate_prior);

        GlobalHeightSizeScaler op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);
        RootHeightScaler<PopulationNode> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        PopSizeScaler op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);
        NodeHeightScaler<PopulationNode> op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);
        /* MuRateScaler op4; */
        /* op4.turn_on_auto_optimize(); */
        /* op4.set_auto_optimize_delay(100); */

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 1);
        std::vector<SampleSummarizer<double> > pop_size_summaries(pop_sizes.size());

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;
        SampleSummarizer<double> mu_rate_summary;

        unsigned int nmoves_per_op = 1;
        unsigned int niterations = 400000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            rop.operate(rng, &tree, 1, 1);
            op2.operate(rng, &tree, 1, nmoves_per_op);
            op3.operate(rng, &tree, 1, nmoves_per_op);
            /* op4.operate(rng, &tree, 1, 1); */
            if ((i + 1) % sample_freq == 0) {
                pop_sizes = tree.get_pointers_to_population_sizes();
                REQUIRE(pop_sizes.size() == 1);
                for (unsigned int i = 0; i < pop_sizes.size(); ++i) {
                    pop_size_summaries.at(i).add_sample(pop_sizes.at(i)->get_value());
                }
                root_height_summary.add_sample(tree.get_root_height());
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                mu_rate_summary.add_sample(tree.get_mutation_rate());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << rop.header_string();
        std::cout << rop.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();
        std::cout << op3.header_string();
        std::cout << op3.to_string();
        /* std::cout << op4.header_string(); */
        /* std::cout << op4.to_string(); */

        REQUIRE(op.get_number_of_attempts() == niterations);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(rop.get_number_of_attempts() == niterations);
        REQUIRE(rop.get_number_of_attempts_for_correction() == (niterations - 100));
        REQUIRE(op2.get_number_of_attempts() == (niterations * nmoves_per_op));
        REQUIRE(op2.get_number_of_attempts_for_correction() == ((nmoves_per_op * niterations) - 100));
        REQUIRE(op3.get_number_of_attempts() == (niterations * nmoves_per_op));
        REQUIRE(op3.get_number_of_attempts_for_correction() == ((nmoves_per_op * niterations) - 100));
        /* REQUIRE(op4.get_number_of_attempts() == niterations); */
        /* REQUIRE(op4.get_number_of_attempts_for_correction() == (niterations - 100)); */

        for (unsigned int i = 0; i < pop_size_summaries.size(); ++i) {
            std::cout << "pop size mean: " <<  pop_size_summaries.at(i).mean() << "\n";
            std::cout << "pop size var:  " <<  pop_size_summaries.at(i).variance() << "\n";
        }

        double eps = 0.001;
        for (unsigned int i = 0; i < pop_size_summaries.size(); ++i) {
            REQUIRE(pop_size_summaries.at(i).sample_size() == nsamples);
            REQUIRE(pop_size_summaries.at(i).mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
            REQUIRE(pop_size_summaries.at(i).variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        }

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(mu_rate_summary.sample_size() == nsamples);

        BetaDistribution prior(1.0, 1.0);

        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        /* REQUIRE(mu_rate_summary.mean() == Approx(mu_rate_prior->get_mean()).epsilon(eps)); */
        REQUIRE(mu_rate_summary.variance() == Approx(0.0).epsilon(eps));
    }
}

TEST_CASE("Testing GlobalHeightRateScaler with 3 leaves, unconstrained sizes, and optimizing",
        "[GlobalHeightRateScaler]") {

    SECTION("Testing 3 leaves, unconstrained sizes, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(24652345);

        double mu_rate_shape = 10.0;
        double mu_rate_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> mu_rate_prior = std::make_shared<GammaDistribution>(
                mu_rate_shape,
                mu_rate_scale);

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double pop_size_shape = 10.0;
        double pop_size_scale = 0.02;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.estimate_root_height();
        tree.estimate_mutation_rate();

        tree.set_population_size_prior(pop_size_prior);
        tree.set_root_node_height_prior(root_height_prior);
        tree.set_mutation_rate_prior(mu_rate_prior);

        GlobalHeightRateScaler op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);
        RootHeightScaler<PopulationNode> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        /* PopSizeScaler op2; */
        /* op2.turn_on_auto_optimize(); */
        /* op2.set_auto_optimize_delay(100); */
        NodeHeightScaler<PopulationNode> op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);
        MuRateScaler op4;
        op4.turn_on_auto_optimize();
        op4.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 5);
        std::vector<SampleSummarizer<double> > pop_size_summaries(pop_sizes.size());

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;
        SampleSummarizer<double> mu_rate_summary;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
            rop.operate(rng, &tree, 1, 1);
            /* op2.operate(rng, &tree, 1, nmoves_per_op); */
            op3.operate(rng, &tree, 1, 1);
            op4.operate(rng, &tree, 1, 1);
        }

        unsigned int nmoves_per_op = 1;
        unsigned int niterations = 400000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            rop.operate(rng, &tree, 1, 1);
            /* op2.operate(rng, &tree, 1, nmoves_per_op); */
            op3.operate(rng, &tree, 1, 1);
            op4.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                pop_sizes = tree.get_pointers_to_population_sizes();
                REQUIRE(pop_sizes.size() == 5);
                for (unsigned int i = 0; i < pop_sizes.size(); ++i) {
                    pop_size_summaries.at(i).add_sample(pop_sizes.at(i)->get_value());
                }
                root_height_summary.add_sample(tree.get_root_height());
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                mu_rate_summary.add_sample(tree.get_mutation_rate());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << rop.header_string();
        std::cout << rop.to_string();
        /* std::cout << op2.header_string(); */
        /* std::cout << op2.to_string(); */
        std::cout << op3.header_string();
        std::cout << op3.to_string();
        std::cout << op4.header_string();
        std::cout << op4.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));
        REQUIRE(rop.get_number_of_attempts() == niterations + burnin);
        REQUIRE(rop.get_number_of_attempts_for_correction() == (niterations + burnin - 100));
        /* REQUIRE(op2.get_number_of_attempts() == ((niterations + burnin) * nmoves_per_op)); */
        /* REQUIRE(op2.get_number_of_attempts_for_correction() == ((nmoves_per_op * (burnin + niterations)) - 100)); */
        REQUIRE(op3.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op3.get_number_of_attempts_for_correction() == (niterations + burnin - 100));
        REQUIRE(op4.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op4.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        double eps = 0.001;
        for (unsigned int i = 0; i < pop_size_summaries.size(); ++i) {
            REQUIRE(pop_size_summaries.at(i).sample_size() == nsamples);
            /* REQUIRE(pop_size_summaries.at(i).mean() == Approx(pop_size_prior->get_mean()).epsilon(eps)); */
            REQUIRE(pop_size_summaries.at(i).variance() == Approx(0.0).epsilon(eps));
        }

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(mu_rate_summary.sample_size() == nsamples);

        BetaDistribution prior(1.0, 1.0);

        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(mu_rate_summary.mean() == Approx(mu_rate_prior->get_mean()).epsilon(eps));
        REQUIRE(mu_rate_summary.variance() == Approx(mu_rate_prior->get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing GlobalHeightRateScaler with 3 leaves, constrained sizes, and optimizing",
        "[GlobalHeightRateScaler]") {

    SECTION("Testing 3 leaves, unconstrained sizes, and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(24652345);

        double mu_rate_shape = 10.0;
        double mu_rate_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> mu_rate_prior = std::make_shared<GammaDistribution>(
                mu_rate_shape,
                mu_rate_scale);

        double root_height_shape = 20.0;
        double root_height_scale = 0.025;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double pop_size_shape = 10.0;
        double pop_size_scale = 0.02;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.estimate_root_height();
        tree.estimate_mutation_rate();
        tree.constrain_population_sizes();

        tree.set_population_size_prior(pop_size_prior);
        tree.set_root_node_height_prior(root_height_prior);
        tree.set_mutation_rate_prior(mu_rate_prior);

        GlobalHeightRateScaler op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);
        RootHeightScaler<PopulationNode> rop;
        rop.turn_on_auto_optimize();
        rop.set_auto_optimize_delay(100);

        /* PopSizeScaler op2; */
        /* op2.turn_on_auto_optimize(); */
        /* op2.set_auto_optimize_delay(100); */
        NodeHeightScaler<PopulationNode> op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);
        MuRateScaler op4;
        op4.turn_on_auto_optimize();
        op4.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(rop.auto_optimizing());
        REQUIRE(rop.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 1);
        std::vector<SampleSummarizer<double> > pop_size_summaries(pop_sizes.size());

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_height_summary;
        SampleSummarizer<double> mu_rate_summary;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
            rop.operate(rng, &tree, 1, 1);
            /* op2.operate(rng, &tree, 1, nmoves_per_op); */
            op3.operate(rng, &tree, 1, 1);
            op4.operate(rng, &tree, 1, 1);
        }

        unsigned int nmoves_per_op = 1;
        unsigned int niterations = 400000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            rop.operate(rng, &tree, 1, 1);
            /* op2.operate(rng, &tree, 1, nmoves_per_op); */
            op3.operate(rng, &tree, 1, 1);
            op4.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                pop_sizes = tree.get_pointers_to_population_sizes();
                REQUIRE(pop_sizes.size() == 1);
                for (unsigned int i = 0; i < pop_sizes.size(); ++i) {
                    pop_size_summaries.at(i).add_sample(pop_sizes.at(i)->get_value());
                }
                root_height_summary.add_sample(tree.get_root_height());
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                mu_rate_summary.add_sample(tree.get_mutation_rate());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << rop.header_string();
        std::cout << rop.to_string();
        /* std::cout << op2.header_string(); */
        /* std::cout << op2.to_string(); */
        std::cout << op3.header_string();
        std::cout << op3.to_string();
        std::cout << op4.header_string();
        std::cout << op4.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));
        REQUIRE(rop.get_number_of_attempts() == niterations + burnin);
        REQUIRE(rop.get_number_of_attempts_for_correction() == (niterations + burnin - 100));
        /* REQUIRE(op2.get_number_of_attempts() == ((niterations + burnin) * nmoves_per_op)); */
        /* REQUIRE(op2.get_number_of_attempts_for_correction() == ((nmoves_per_op * (burnin + niterations)) - 100)); */
        REQUIRE(op3.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op3.get_number_of_attempts_for_correction() == (niterations + burnin - 100));
        REQUIRE(op4.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op4.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        double eps = 0.001;
        for (unsigned int i = 0; i < pop_size_summaries.size(); ++i) {
            REQUIRE(pop_size_summaries.at(i).sample_size() == nsamples);
            /* REQUIRE(pop_size_summaries.at(i).mean() == Approx(pop_size_prior->get_mean()).epsilon(eps)); */
            REQUIRE(pop_size_summaries.at(i).variance() == Approx(0.0).epsilon(eps));
        }

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_height_summary.sample_size() == nsamples);
        REQUIRE(mu_rate_summary.sample_size() == nsamples);

        BetaDistribution prior(1.0, 1.0);

        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(mu_rate_summary.mean() == Approx(mu_rate_prior->get_mean()).epsilon(eps));
        REQUIRE(mu_rate_summary.variance() == Approx(mu_rate_prior->get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing StateFreqMover with 3 leaves and optimizing",
        "[StateFreqOperator]") {

    SECTION("Testing 3 leaves and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(286431321);

        double freq_a = 3.0;
        double freq_b = 2.0;
        std::shared_ptr<ContinuousProbabilityDistribution> freq_prior = std::make_shared<BetaDistribution>(
                freq_a,
                freq_b);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.fix_root_height();
        tree.fix_mutation_rate();
        tree.fix_population_sizes();
        tree.estimate_state_frequencies();

        tree.set_freq_1_prior(freq_prior);

        StateFreqMover op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> freq_summary;
        std::vector<double> freq_samples;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
        }

        unsigned int niterations = 200000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                freq_summary.add_sample(tree.get_freq_1());
                freq_samples.push_back(tree.get_freq_1());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess = effective_sample_size<double>(freq_samples);
        std::cout << "ESS: " << ess << "\n";

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(freq_summary.sample_size() == nsamples);

        double eps = 0.001;

        REQUIRE(freq_summary.mean() == Approx(freq_prior->get_mean()).epsilon(eps));
        REQUIRE(freq_summary.variance() == Approx(freq_prior->get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing StateFreqDirichletOperator with 3 leaves and optimizing",
        "[StateFreqOperator]") {

    SECTION("Testing 3 leaves and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(264123456);

        double freq_a = 3.0;
        double freq_b = 2.0;
        std::shared_ptr<ContinuousProbabilityDistribution> freq_prior = std::make_shared<BetaDistribution>(
                freq_a,
                freq_b);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.fix_root_height();
        tree.fix_mutation_rate();
        tree.fix_population_sizes();
        tree.estimate_state_frequencies();

        tree.set_freq_1_prior(freq_prior);

        StateFreqDirichletOperator op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> freq_summary;
        std::vector<double> freq_samples;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
        }

        unsigned int niterations = 200000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                freq_summary.add_sample(tree.get_freq_1());
                freq_samples.push_back(tree.get_freq_1());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess = effective_sample_size<double>(freq_samples);
        std::cout << "ESS: " << ess << "\n";

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(freq_summary.sample_size() == nsamples);

        double eps = 0.001;

        REQUIRE(freq_summary.mean() == Approx(freq_prior->get_mean()).epsilon(eps));
        REQUIRE(freq_summary.variance() == Approx(freq_prior->get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing StateFreqBetaOperator with 3 leaves and optimizing",
        "[StateFreqOperator]") {

    SECTION("Testing 3 leaves and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(251725179);

        double freq_a = 3.0;
        double freq_b = 2.0;
        std::shared_ptr<ContinuousProbabilityDistribution> freq_prior = std::make_shared<BetaDistribution>(
                freq_a,
                freq_b);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(4, "root", 0.5);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(3, "internal0", 0.25);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BasePopulationTree tree(root);

        tree.ignore_data();
        tree.fix_root_height();
        tree.fix_mutation_rate();
        tree.fix_population_sizes();
        tree.estimate_state_frequencies();

        tree.set_freq_1_prior(freq_prior);

        StateFreqBetaOperator op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> freq_summary;
        std::vector<double> freq_samples;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
        }

        unsigned int niterations = 200000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                freq_summary.add_sample(tree.get_freq_1());
                freq_samples.push_back(tree.get_freq_1());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess = effective_sample_size<double>(freq_samples);
        std::cout << "ESS: " << ess << "\n";

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(freq_summary.sample_size() == nsamples);

        double eps = 0.001;

        REQUIRE(freq_summary.mean() == Approx(freq_prior->get_mean()).epsilon(eps));
        REQUIRE(freq_summary.variance() == Approx(freq_prior->get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing GlobalNodeHeightDirichletOperator with 3 leaves, fixed root, and optimizing",
        "[GlobalNodeHeightDirichletOperator]") {

    SECTION("Testing 3 leaves with fixed root and optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(1);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.3);
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

        GlobalNodeHeightDirichletOperator<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_height_summary;
        std::vector<double> internal_height_samples;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
        }

        unsigned int niterations = 200000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                internal_height_samples.push_back(tree.get_height(0) / tree.get_height(1));
                REQUIRE(tree.get_root_height() == 0.5);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess = effective_sample_size<double>(internal_height_samples);
        std::cout << "ESS: " << ess << "\n";

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(internal_height_summary.sample_size() == nsamples);
        
        BetaDistribution prior(1.0, 1.0);

        double eps = 0.001;
        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing GlobalNodeHeightDirichletOperator with 5 leaf ladder and optimizing",
        "[GlobalNodeHeightDirichletOperator]") {

    SECTION("Testing 5 leaf ladder optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(9873746320);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>(8, "root", 0.5);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(7, "internal2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(6, "internal1", 0.3);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(5, "internal0", 0.2);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        internal2->add_child(internal1);
        internal2->add_child(leaf3);
        root->add_child(internal2);
        root->add_child(leaf4);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(4.0);
        tree.set_beta_of_node_height_beta_prior(3.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        GlobalNodeHeightDirichletOperator<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        SampleSummarizer<double> internal_2_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;
        std::vector<double> internal_2_height_samples;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
        }

        unsigned int niterations = 500000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_0_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                internal_1_height_summary.add_sample(tree.get_height(1) / tree.get_height(2));
                internal_2_height_summary.add_sample(tree.get_height(2) / tree.get_height(3));
                internal_0_height_samples.push_back(tree.get_height(0) / tree.get_height(1));
                internal_1_height_samples.push_back(tree.get_height(1) / tree.get_height(2));
                internal_2_height_samples.push_back(tree.get_height(2) / tree.get_height(3));
                REQUIRE(tree.get_root_height() == 0.5);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        double ess_2 = effective_sample_size<double>(internal_2_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        std::cout << "ESS 2: " << ess_2 << "\n";
        
        BetaDistribution prior(4.0, 3.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
        std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);
        REQUIRE(internal_2_height_summary.sample_size() == nsamples);

        double eps = 0.001;
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing GlobalNodeHeightDirichletOperator with 5 leaf balanced and optimizing",
        "[GlobalNodeHeightDirichletOperator]") {

    SECTION("Testing 5 leaf balanced optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(164651453);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>(8, "root", 0.5);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(7, "internal2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(6, "internal1", 0.3);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(5, "internal0", 0.2);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        internal2->add_child(leaf3);
        internal2->add_child(leaf4);
        root->add_child(internal1);
        root->add_child(internal2);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(4.0);
        tree.set_beta_of_node_height_beta_prior(3.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        GlobalNodeHeightDirichletOperator<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        SampleSummarizer<double> internal_2_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;
        std::vector<double> internal_2_height_samples;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
        }

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        std::shared_ptr<Node> node;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                node = tree.get_node("internal0");
                internal_0_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_0_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->has_parent());
                REQUIRE(node->get_parent()->get_label() == "internal1");
                node = tree.get_node("internal1");
                internal_1_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_1_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");
                node = tree.get_node("internal2");
                internal_2_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_2_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");
                REQUIRE(tree.get_root_height() == 0.5);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        double ess_2 = effective_sample_size<double>(internal_2_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        std::cout << "ESS 2: " << ess_2 << "\n";
        
        BetaDistribution prior(4.0, 3.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
        std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);
        REQUIRE(internal_2_height_summary.sample_size() == nsamples);

        double eps = 0.001;
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing NodeHeightScaler with 5 leaf ladder and optimizing",
        "[NodeHeightScaler]") {

    SECTION("Testing 5 leaf ladder optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(164651453);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>(8, "root", 0.5);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(7, "internal2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(6, "internal1", 0.3);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(5, "internal0", 0.2);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        internal2->add_child(internal1);
        internal2->add_child(leaf3);
        root->add_child(internal2);
        root->add_child(leaf4);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(4.0);
        tree.set_beta_of_node_height_beta_prior(3.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        SampleSummarizer<double> internal_2_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;
        std::vector<double> internal_2_height_samples;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
        }

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_0_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                internal_1_height_summary.add_sample(tree.get_height(1) / tree.get_height(2));
                internal_2_height_summary.add_sample(tree.get_height(2) / tree.get_height(3));
                internal_0_height_samples.push_back(tree.get_height(0) / tree.get_height(1));
                internal_1_height_samples.push_back(tree.get_height(1) / tree.get_height(2));
                internal_2_height_samples.push_back(tree.get_height(2) / tree.get_height(3));
                REQUIRE(tree.get_root_height() == 0.5);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        double ess_2 = effective_sample_size<double>(internal_2_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        std::cout << "ESS 2: " << ess_2 << "\n";
        
        BetaDistribution prior(4.0, 3.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
        std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);
        REQUIRE(internal_2_height_summary.sample_size() == nsamples);

        double eps = 0.001;
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing NodeHeightScaler with 5 leaf balanced and optimizing",
        "[NodeHeightScaler]") {

    SECTION("Testing 5 leaf balanced optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(164651453);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>(8, "root", 0.5);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(7, "internal2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(6, "internal1", 0.3);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(5, "internal0", 0.2);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        internal2->add_child(leaf3);
        internal2->add_child(leaf4);
        root->add_child(internal1);
        root->add_child(internal2);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(4.0);
        tree.set_beta_of_node_height_beta_prior(3.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        SampleSummarizer<double> internal_2_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;
        std::vector<double> internal_2_height_samples;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
        }

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        std::shared_ptr<Node> node;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                node = tree.get_node("internal0");
                internal_0_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_0_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->has_parent());
                REQUIRE(node->get_parent()->get_label() == "internal1");
                node = tree.get_node("internal1");
                internal_1_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_1_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");
                node = tree.get_node("internal2");
                internal_2_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_2_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");
                REQUIRE(tree.get_root_height() == 0.5);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        double ess_2 = effective_sample_size<double>(internal_2_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        std::cout << "ESS 2: " << ess_2 << "\n";
        
        BetaDistribution prior(4.0, 3.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
        std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);
        REQUIRE(internal_2_height_summary.sample_size() == nsamples);

        double eps = 0.001;
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}


TEST_CASE("Testing NodeHeightDirichletOperator with 5 leaf ladder and optimizing",
        "[NodeHeightDirichletOperator]") {

    SECTION("Testing 5 leaf ladder optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(164651453);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>(8, "root", 0.5);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(7, "internal2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(6, "internal1", 0.3);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(5, "internal0", 0.2);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        internal2->add_child(internal1);
        internal2->add_child(leaf3);
        root->add_child(internal2);
        root->add_child(leaf4);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(4.0);
        tree.set_beta_of_node_height_beta_prior(3.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightDirichletOperator<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        SampleSummarizer<double> internal_2_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;
        std::vector<double> internal_2_height_samples;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
        }

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_0_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                internal_1_height_summary.add_sample(tree.get_height(1) / tree.get_height(2));
                internal_2_height_summary.add_sample(tree.get_height(2) / tree.get_height(3));
                internal_0_height_samples.push_back(tree.get_height(0) / tree.get_height(1));
                internal_1_height_samples.push_back(tree.get_height(1) / tree.get_height(2));
                internal_2_height_samples.push_back(tree.get_height(2) / tree.get_height(3));
                REQUIRE(tree.get_root_height() == 0.5);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        double ess_2 = effective_sample_size<double>(internal_2_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        std::cout << "ESS 2: " << ess_2 << "\n";
        
        BetaDistribution prior(4.0, 3.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
        std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);
        REQUIRE(internal_2_height_summary.sample_size() == nsamples);

        double eps = 0.001;
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing NodeHeightDirichletOperator with 5 leaf balanced and optimizing",
        "[NodeHeightDirichletOperator]") {

    SECTION("Testing 5 leaf balanced optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(164651453);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>(8, "root", 0.5);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(7, "internal2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(6, "internal1", 0.3);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(5, "internal0", 0.2);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        internal2->add_child(leaf3);
        internal2->add_child(leaf4);
        root->add_child(internal1);
        root->add_child(internal2);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(4.0);
        tree.set_beta_of_node_height_beta_prior(3.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightDirichletOperator<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        SampleSummarizer<double> internal_2_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;
        std::vector<double> internal_2_height_samples;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
        }

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        std::shared_ptr<Node> node;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                node = tree.get_node("internal0");
                internal_0_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_0_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->has_parent());
                REQUIRE(node->get_parent()->get_label() == "internal1");
                node = tree.get_node("internal1");
                internal_1_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_1_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");
                node = tree.get_node("internal2");
                internal_2_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_2_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");
                REQUIRE(tree.get_root_height() == 0.5);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        double ess_2 = effective_sample_size<double>(internal_2_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        std::cout << "ESS 2: " << ess_2 << "\n";
        
        BetaDistribution prior(4.0, 3.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
        std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);
        REQUIRE(internal_2_height_summary.sample_size() == nsamples);

        double eps = 0.001;
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// This operator is note working correctly
///////////////////////////////////////////////////////////////////////////////
// TEST_CASE("Testing NodeHeightBetaOperator with 5 leaf ladder and optimizing",
//         "[NodeHeightBetaOperator]") {
// 
//     SECTION("Testing 5 leaf ladder optimizing") {
//         RandomNumberGenerator rng = RandomNumberGenerator(164651453);
// 
//         double root_height_shape = 10.0;
//         double root_height_scale = 0.05;
//         std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
//                 root_height_shape,
//                 root_height_scale);
// 
//         std::shared_ptr<Node> root = std::make_shared<Node>(8, "root", 0.5);
//         std::shared_ptr<Node> internal2 = std::make_shared<Node>(7, "internal2", 0.4);
//         std::shared_ptr<Node> internal1 = std::make_shared<Node>(6, "internal1", 0.3);
//         std::shared_ptr<Node> internal0 = std::make_shared<Node>(5, "internal0", 0.2);
//         std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
//         std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
//         std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
//         std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
//         std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);
// 
//         internal0->add_child(leaf0);
//         internal0->add_child(leaf1);
//         internal1->add_child(internal0);
//         internal1->add_child(leaf2);
//         internal2->add_child(internal1);
//         internal2->add_child(leaf3);
//         root->add_child(internal2);
//         root->add_child(leaf4);
// 
//         BaseTree<Node> tree(root);
//         tree.set_root_node_height_prior(root_height_prior);
// 
//         tree.estimate_alpha_of_node_height_beta_prior();
//         tree.estimate_beta_of_node_height_beta_prior();
//         tree.set_alpha_of_node_height_beta_prior(4.0);
//         tree.set_beta_of_node_height_beta_prior(3.0);
//         tree.fix_alpha_of_node_height_beta_prior();
//         tree.fix_beta_of_node_height_beta_prior();
// 
//         tree.ignore_data();
//         tree.estimate_root_height();
// 
//         NodeHeightBetaOperator<Node> op;
//         op.turn_on_auto_optimize();
//         op.set_auto_optimize_delay(100);
// 
//         REQUIRE(op.auto_optimizing());
//         REQUIRE(op.get_auto_optimize_delay() == 100);
// 
//         // Initialize prior probs
//         tree.compute_log_likelihood_and_prior(true);
// 
//         SampleSummarizer<double> internal_0_height_summary;
//         SampleSummarizer<double> internal_1_height_summary;
//         SampleSummarizer<double> internal_2_height_summary;
//         std::vector<double> internal_0_height_samples;
//         std::vector<double> internal_1_height_samples;
//         std::vector<double> internal_2_height_samples;
// 
//         // burnin
//         unsigned int burnin = 1000;
//         for (unsigned int i = 0; i < burnin; ++i) {
//             op.operate(rng, &tree, 1, 1);
//         }
// 
//         unsigned int niterations = 1000000;
//         unsigned int sample_freq = 10;
//         unsigned int nsamples = niterations / sample_freq;
//         for (unsigned int i = 0; i < niterations; ++i) {
//             op.operate(rng, &tree, 1, 1);
//             if ((i + 1) % sample_freq == 0) {
//                 internal_0_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
//                 internal_1_height_summary.add_sample(tree.get_height(1) / tree.get_height(2));
//                 internal_2_height_summary.add_sample(tree.get_height(2) / tree.get_height(3));
//                 internal_0_height_samples.push_back(tree.get_height(0) / tree.get_height(1));
//                 internal_1_height_samples.push_back(tree.get_height(1) / tree.get_height(2));
//                 internal_2_height_samples.push_back(tree.get_height(2) / tree.get_height(3));
//                 REQUIRE(tree.get_root_height() == 0.5);
//             }
//         }
//         std::cout << op.header_string();
//         std::cout << op.to_string();
//         double ess_0 = effective_sample_size<double>(internal_0_height_samples);
//         double ess_1 = effective_sample_size<double>(internal_1_height_samples);
//         double ess_2 = effective_sample_size<double>(internal_2_height_samples);
//         std::cout << "ESS 0: " << ess_0 << "\n";
//         std::cout << "ESS 1: " << ess_1 << "\n";
//         std::cout << "ESS 2: " << ess_2 << "\n";
//         
//         BetaDistribution prior(4.0, 3.0);
// 
//         std::cout << "Expected mean: " << prior.get_mean() << "\n";
//         std::cout << "Expected var: " << prior.get_variance() << "\n";
//         std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
//         std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
//         std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
//         std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
//         std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
//         std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";
// 
//         REQUIRE(op.get_number_of_attempts() == niterations + burnin);
//         REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));
// 
//         REQUIRE(internal_0_height_summary.sample_size() == nsamples);
//         REQUIRE(internal_1_height_summary.sample_size() == nsamples);
//         REQUIRE(internal_2_height_summary.sample_size() == nsamples);
// 
//         double eps = 0.001;
//         REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
//         REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
//         REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
//         REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
//         REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
//         REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
//     }
// }
// 
// TEST_CASE("Testing NodeHeightBetaOperator with 5 leaf balanced and optimizing",
//         "[NodeHeightBetaOperator]") {
// 
//     SECTION("Testing 5 leaf balanced optimizing") {
//         RandomNumberGenerator rng = RandomNumberGenerator(164651453);
// 
//         double root_height_shape = 10.0;
//         double root_height_scale = 0.05;
//         std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
//                 root_height_shape,
//                 root_height_scale);
// 
//         std::shared_ptr<Node> root = std::make_shared<Node>(8, "root", 0.5);
//         std::shared_ptr<Node> internal2 = std::make_shared<Node>(7, "internal2", 0.4);
//         std::shared_ptr<Node> internal1 = std::make_shared<Node>(6, "internal1", 0.3);
//         std::shared_ptr<Node> internal0 = std::make_shared<Node>(5, "internal0", 0.2);
//         std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
//         std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
//         std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
//         std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
//         std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);
// 
//         internal0->add_child(leaf0);
//         internal0->add_child(leaf1);
//         internal1->add_child(internal0);
//         internal1->add_child(leaf2);
//         internal2->add_child(leaf3);
//         internal2->add_child(leaf4);
//         root->add_child(internal1);
//         root->add_child(internal2);
// 
//         BaseTree<Node> tree(root);
//         tree.set_root_node_height_prior(root_height_prior);
// 
//         tree.estimate_alpha_of_node_height_beta_prior();
//         tree.estimate_beta_of_node_height_beta_prior();
//         tree.set_alpha_of_node_height_beta_prior(4.0);
//         tree.set_beta_of_node_height_beta_prior(3.0);
//         tree.fix_alpha_of_node_height_beta_prior();
//         tree.fix_beta_of_node_height_beta_prior();
// 
//         tree.ignore_data();
//         tree.estimate_root_height();
// 
//         NodeHeightBetaOperator<Node> op;
//         op.turn_on_auto_optimize();
//         op.set_auto_optimize_delay(100);
// 
//         REQUIRE(op.auto_optimizing());
//         REQUIRE(op.get_auto_optimize_delay() == 100);
// 
//         // Initialize prior probs
//         tree.compute_log_likelihood_and_prior(true);
// 
//         SampleSummarizer<double> internal_0_height_summary;
//         SampleSummarizer<double> internal_1_height_summary;
//         SampleSummarizer<double> internal_2_height_summary;
//         std::vector<double> internal_0_height_samples;
//         std::vector<double> internal_1_height_samples;
//         std::vector<double> internal_2_height_samples;
// 
//         // burnin
//         unsigned int burnin = 1000;
//         for (unsigned int i = 0; i < burnin; ++i) {
//             op.operate(rng, &tree, 1, 1);
//         }
// 
//         unsigned int niterations = 1000000;
//         unsigned int sample_freq = 10;
//         unsigned int nsamples = niterations / sample_freq;
//         std::shared_ptr<Node> node;
//         for (unsigned int i = 0; i < niterations; ++i) {
//             op.operate(rng, &tree, 1, 1);
//             if ((i + 1) % sample_freq == 0) {
//                 node = tree.get_node("internal0");
//                 internal_0_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
//                 internal_0_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
//                 REQUIRE(node->get_parent()->has_parent());
//                 REQUIRE(node->get_parent()->get_label() == "internal1");
//                 node = tree.get_node("internal1");
//                 internal_1_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
//                 internal_1_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
//                 REQUIRE(node->get_parent()->is_root());
//                 REQUIRE(node->get_parent()->get_label() == "root");
//                 node = tree.get_node("internal2");
//                 internal_2_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
//                 internal_2_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
//                 REQUIRE(node->get_parent()->is_root());
//                 REQUIRE(node->get_parent()->get_label() == "root");
//                 REQUIRE(tree.get_root_height() == 0.5);
//             }
//         }
//         std::cout << op.header_string();
//         std::cout << op.to_string();
//         double ess_0 = effective_sample_size<double>(internal_0_height_samples);
//         double ess_1 = effective_sample_size<double>(internal_1_height_samples);
//         double ess_2 = effective_sample_size<double>(internal_2_height_samples);
//         std::cout << "ESS 0: " << ess_0 << "\n";
//         std::cout << "ESS 1: " << ess_1 << "\n";
//         std::cout << "ESS 2: " << ess_2 << "\n";
//         
//         BetaDistribution prior(4.0, 3.0);
// 
//         std::cout << "Expected mean: " << prior.get_mean() << "\n";
//         std::cout << "Expected var: " << prior.get_variance() << "\n";
//         std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
//         std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
//         std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
//         std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
//         std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
//         std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";
// 
//         REQUIRE(op.get_number_of_attempts() == niterations + burnin);
//         REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));
// 
//         REQUIRE(internal_0_height_summary.sample_size() == nsamples);
//         REQUIRE(internal_1_height_summary.sample_size() == nsamples);
//         REQUIRE(internal_2_height_summary.sample_size() == nsamples);
// 
//         double eps = 0.001;
//         REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
//         REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
//         REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
//         REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
//         REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
//         REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
//     }
// }
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("Testing NodeHeightMover with 5 leaf ladder and optimizing",
        "[NodeHeightMover]") {

    SECTION("Testing 5 leaf ladder optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(164651453);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>(8, "root", 0.5);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(7, "internal2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(6, "internal1", 0.3);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(5, "internal0", 0.2);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        internal2->add_child(internal1);
        internal2->add_child(leaf3);
        root->add_child(internal2);
        root->add_child(leaf4);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(4.0);
        tree.set_beta_of_node_height_beta_prior(3.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightMover<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        SampleSummarizer<double> internal_2_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;
        std::vector<double> internal_2_height_samples;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
        }

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_0_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                internal_1_height_summary.add_sample(tree.get_height(1) / tree.get_height(2));
                internal_2_height_summary.add_sample(tree.get_height(2) / tree.get_height(3));
                internal_0_height_samples.push_back(tree.get_height(0) / tree.get_height(1));
                internal_1_height_samples.push_back(tree.get_height(1) / tree.get_height(2));
                internal_2_height_samples.push_back(tree.get_height(2) / tree.get_height(3));
                REQUIRE(tree.get_root_height() == 0.5);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        double ess_2 = effective_sample_size<double>(internal_2_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        std::cout << "ESS 2: " << ess_2 << "\n";
        
        BetaDistribution prior(4.0, 3.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
        std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);
        REQUIRE(internal_2_height_summary.sample_size() == nsamples);

        double eps = 0.001;
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing NodeHeightMover with 5 leaf balanced and optimizing",
        "[NodeHeightMover]") {

    SECTION("Testing 5 leaf balanced optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(164635246);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>(8, "root", 0.5);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(7, "internal2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(6, "internal1", 0.3);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(5, "internal0", 0.2);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        internal2->add_child(leaf3);
        internal2->add_child(leaf4);
        root->add_child(internal1);
        root->add_child(internal2);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(4.0);
        tree.set_beta_of_node_height_beta_prior(3.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightMover<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        SampleSummarizer<double> internal_2_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;
        std::vector<double> internal_2_height_samples;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
        }

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        std::shared_ptr<Node> node;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                node = tree.get_node("internal0");
                internal_0_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_0_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->has_parent());
                REQUIRE(node->get_parent()->get_label() == "internal1");
                node = tree.get_node("internal1");
                internal_1_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_1_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");
                node = tree.get_node("internal2");
                internal_2_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_2_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");
                REQUIRE(tree.get_root_height() == 0.5);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        double ess_2 = effective_sample_size<double>(internal_2_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        std::cout << "ESS 2: " << ess_2 << "\n";
        
        BetaDistribution prior(4.0, 3.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
        std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);
        REQUIRE(internal_2_height_summary.sample_size() == nsamples);

        double eps = 0.001;
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}


TEST_CASE("Testing NodeHeightSlideBumpScaler with 5 leaf ladder and optimizing",
        "[NodeHeightSlideBumpScaler]") {

    SECTION("Testing 5 leaf ladder optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(164651453);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>(8, "root", 0.5);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(7, "internal2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(6, "internal1", 0.3);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(5, "internal0", 0.2);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        internal2->add_child(internal1);
        internal2->add_child(leaf3);
        root->add_child(internal2);
        root->add_child(leaf4);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(4.0);
        tree.set_beta_of_node_height_beta_prior(3.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightSlideBumpScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        SampleSummarizer<double> internal_2_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;
        std::vector<double> internal_2_height_samples;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
        }

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_0_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                internal_1_height_summary.add_sample(tree.get_height(1) / tree.get_height(2));
                internal_2_height_summary.add_sample(tree.get_height(2) / tree.get_height(3));
                internal_0_height_samples.push_back(tree.get_height(0) / tree.get_height(1));
                internal_1_height_samples.push_back(tree.get_height(1) / tree.get_height(2));
                internal_2_height_samples.push_back(tree.get_height(2) / tree.get_height(3));
                REQUIRE(tree.get_root_height() == 0.5);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        double ess_2 = effective_sample_size<double>(internal_2_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        std::cout << "ESS 2: " << ess_2 << "\n";
        
        BetaDistribution prior(4.0, 3.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
        std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);
        REQUIRE(internal_2_height_summary.sample_size() == nsamples);

        double eps = 0.001;
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing NodeHeightSlideBumpScaler with 5 leaf balanced and optimizing",
        "[NodeHeightSlideBumpScaler]") {

    SECTION("Testing 5 leaf balanced optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(164651453);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>(8, "root", 0.5);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(7, "internal2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(6, "internal1", 0.3);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(5, "internal0", 0.2);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        internal2->add_child(leaf3);
        internal2->add_child(leaf4);
        root->add_child(internal1);
        root->add_child(internal2);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(4.0);
        tree.set_beta_of_node_height_beta_prior(3.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightSlideBumpScaler<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        SampleSummarizer<double> internal_2_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;
        std::vector<double> internal_2_height_samples;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
        }

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        std::shared_ptr<Node> node;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                node = tree.get_node("internal0");
                internal_0_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_0_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->has_parent());
                REQUIRE(node->get_parent()->get_label() == "internal1");
                node = tree.get_node("internal1");
                internal_1_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_1_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");
                node = tree.get_node("internal2");
                internal_2_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_2_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");
                REQUIRE(tree.get_root_height() == 0.5);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        double ess_2 = effective_sample_size<double>(internal_2_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        std::cout << "ESS 2: " << ess_2 << "\n";
        
        BetaDistribution prior(4.0, 3.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
        std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);
        REQUIRE(internal_2_height_summary.sample_size() == nsamples);

        double eps = 0.001;
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}


TEST_CASE("Testing NodeHeightSlideBumpMover with 5 leaf ladder and optimizing",
        "[NodeHeightSlideBumpMover]") {

    SECTION("Testing 5 leaf ladder optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(234524583);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>(8, "root", 0.5);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(7, "internal2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(6, "internal1", 0.3);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(5, "internal0", 0.2);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        internal2->add_child(internal1);
        internal2->add_child(leaf3);
        root->add_child(internal2);
        root->add_child(leaf4);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(4.0);
        tree.set_beta_of_node_height_beta_prior(3.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightSlideBumpMover<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        SampleSummarizer<double> internal_2_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;
        std::vector<double> internal_2_height_samples;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
        }

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                internal_0_height_summary.add_sample(tree.get_height(0) / tree.get_height(1));
                internal_1_height_summary.add_sample(tree.get_height(1) / tree.get_height(2));
                internal_2_height_summary.add_sample(tree.get_height(2) / tree.get_height(3));
                internal_0_height_samples.push_back(tree.get_height(0) / tree.get_height(1));
                internal_1_height_samples.push_back(tree.get_height(1) / tree.get_height(2));
                internal_2_height_samples.push_back(tree.get_height(2) / tree.get_height(3));
                REQUIRE(tree.get_root_height() == 0.5);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        double ess_2 = effective_sample_size<double>(internal_2_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        std::cout << "ESS 2: " << ess_2 << "\n";
        
        BetaDistribution prior(4.0, 3.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
        std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);
        REQUIRE(internal_2_height_summary.sample_size() == nsamples);

        double eps = 0.001;
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing NodeHeightSlideBumpMover with 5 leaf balanced and optimizing",
        "[NodeHeightSlideBumpMover]") {

    SECTION("Testing 5 leaf balanced optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(164651453);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        std::shared_ptr<Node> root = std::make_shared<Node>(8, "root", 0.5);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>(7, "internal2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>(6, "internal1", 0.3);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>(5, "internal0", 0.2);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        internal2->add_child(leaf3);
        internal2->add_child(leaf4);
        root->add_child(internal1);
        root->add_child(internal2);

        BaseTree<Node> tree(root);
        tree.set_root_node_height_prior(root_height_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(4.0);
        tree.set_beta_of_node_height_beta_prior(3.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        NodeHeightSlideBumpMover<Node> op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        SampleSummarizer<double> internal_2_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;
        std::vector<double> internal_2_height_samples;

        // burnin
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
        }

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        std::shared_ptr<Node> node;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                node = tree.get_node("internal0");
                internal_0_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_0_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->has_parent());
                REQUIRE(node->get_parent()->get_label() == "internal1");
                node = tree.get_node("internal1");
                internal_1_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_1_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");
                node = tree.get_node("internal2");
                internal_2_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_2_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");
                REQUIRE(tree.get_root_height() == 0.5);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        double ess_2 = effective_sample_size<double>(internal_2_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        std::cout << "ESS 2: " << ess_2 << "\n";
        
        BetaDistribution prior(4.0, 3.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
        std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));

        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);
        REQUIRE(internal_2_height_summary.sample_size() == nsamples);

        double eps = 0.001;
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing GlobalHeightSizeMixer with 5 leaf balanced and optimizing",
        "[GlobalHeightSizeMixer]") {

    SECTION("Testing 5 leaf balanced optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(164651453);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double pop_size_shape = 10.0;
        double pop_size_scale = 0.02;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(8, "root", 0.5);
        std::shared_ptr<PopulationNode> internal2 = std::make_shared<PopulationNode>(7, "internal2", 0.4);
        std::shared_ptr<PopulationNode> internal1 = std::make_shared<PopulationNode>(6, "internal1", 0.3);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(5, "internal0", 0.2);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);
        std::shared_ptr<PopulationNode> leaf3 = std::make_shared<PopulationNode>(3, "leaf3", 0.0);
        std::shared_ptr<PopulationNode> leaf4 = std::make_shared<PopulationNode>(4, "leaf4", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        internal2->add_child(leaf3);
        internal2->add_child(leaf4);
        root->add_child(internal1);
        root->add_child(internal2);

        BasePopulationTree tree(root);
        tree.set_root_node_height_prior(root_height_prior);
        tree.set_population_size_prior(pop_size_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(4.0);
        tree.set_beta_of_node_height_beta_prior(3.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        GlobalHeightSizeMixer op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);
        RootHeightScaler<PopulationNode> op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);
        NodeHeightScaler<PopulationNode> op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);
        PopSizeScaler op4;
        op4.turn_on_auto_optimize();
        op4.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(op2.auto_optimizing());
        REQUIRE(op2.get_auto_optimize_delay() == 100);
        REQUIRE(op3.auto_optimizing());
        REQUIRE(op3.get_auto_optimize_delay() == 100);
        REQUIRE(op4.auto_optimizing());
        REQUIRE(op4.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> root_height_summary;
        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        SampleSummarizer<double> internal_2_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;
        std::vector<double> internal_2_height_samples;

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 9);
        std::vector<SampleSummarizer<double> > pop_size_summaries(pop_sizes.size());

        // burnin
        unsigned int n_moves_per_op = 3;
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, 1);
            op2.operate(rng, &tree, 1, 1);
            op3.operate(rng, &tree, 1, n_moves_per_op);
            op4.operate(rng, &tree, 1, n_moves_per_op);
        }

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        std::shared_ptr<PopulationNode> node;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            op2.operate(rng, &tree, 1, 1);
            op3.operate(rng, &tree, 1, n_moves_per_op);
            op4.operate(rng, &tree, 1, n_moves_per_op);
            if ((i + 1) % sample_freq == 0) {
                node = tree.get_node("internal0");
                internal_0_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_0_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->has_parent());
                REQUIRE(node->get_parent()->get_label() == "internal1");
                node = tree.get_node("internal1");
                internal_1_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_1_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");
                node = tree.get_node("internal2");
                internal_2_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_2_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");

                root_height_summary.add_sample(tree.get_root_height());

                pop_sizes = tree.get_pointers_to_population_sizes();
                REQUIRE(pop_sizes.size() == 9);
                for (unsigned int i = 0; i < pop_sizes.size(); ++i) {
                    pop_size_summaries.at(i).add_sample(pop_sizes.at(i)->get_value());
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        double ess_2 = effective_sample_size<double>(internal_2_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        std::cout << "ESS 2: " << ess_2 << "\n";
        
        BetaDistribution prior(4.0, 3.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
        std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";

        REQUIRE(op.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op.get_number_of_attempts_for_correction() == (niterations + burnin - 100));
        REQUIRE(op2.get_number_of_attempts() == niterations + burnin);
        REQUIRE(op2.get_number_of_attempts_for_correction() == (niterations + burnin - 100));
        REQUIRE(op3.get_number_of_attempts() == (niterations + burnin) * n_moves_per_op);
        REQUIRE(op3.get_number_of_attempts_for_correction() == (((niterations + burnin) * n_moves_per_op) - 100));
        REQUIRE(op4.get_number_of_attempts() == (niterations + burnin) * n_moves_per_op);
        REQUIRE(op4.get_number_of_attempts_for_correction() == (((niterations + burnin) * n_moves_per_op) - 100));

        REQUIRE(root_height_summary.sample_size() == nsamples);
        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);
        REQUIRE(internal_2_height_summary.sample_size() == nsamples);

        double eps = 0.001;
        REQUIRE(root_height_summary.mean() == Approx(root_height_prior->get_mean()).epsilon(eps));
        REQUIRE(root_height_summary.variance() == Approx(root_height_prior->get_variance()).epsilon(eps));
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        for (unsigned int i = 0; i < pop_size_summaries.size(); ++i) {
            REQUIRE(pop_size_summaries.at(i).sample_size() == nsamples);
            REQUIRE(pop_size_summaries.at(i).mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
            REQUIRE(pop_size_summaries.at(i).variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        }
    }
}

TEST_CASE("Testing HeightSizeMixer with 5 leaf balanced and optimizing",
        "[HeightSizeMixer]") {

    SECTION("Testing 5 leaf balanced optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(164651453);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double pop_size_shape = 10.0;
        double pop_size_scale = 0.02;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(8, "root", 0.5);
        std::shared_ptr<PopulationNode> internal2 = std::make_shared<PopulationNode>(7, "internal2", 0.4);
        std::shared_ptr<PopulationNode> internal1 = std::make_shared<PopulationNode>(6, "internal1", 0.3);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(5, "internal0", 0.2);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);
        std::shared_ptr<PopulationNode> leaf3 = std::make_shared<PopulationNode>(3, "leaf3", 0.0);
        std::shared_ptr<PopulationNode> leaf4 = std::make_shared<PopulationNode>(4, "leaf4", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        internal2->add_child(leaf3);
        internal2->add_child(leaf4);
        root->add_child(internal1);
        root->add_child(internal2);

        BasePopulationTree tree(root);
        tree.set_root_node_height_prior(root_height_prior);
        tree.set_population_size_prior(pop_size_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(4.0);
        tree.set_beta_of_node_height_beta_prior(3.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        HeightSizeMixer op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);
        NodeHeightScaler<PopulationNode> op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);
        PopSizeScaler op4;
        op4.turn_on_auto_optimize();
        op4.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(op3.auto_optimizing());
        REQUIRE(op3.get_auto_optimize_delay() == 100);
        REQUIRE(op4.auto_optimizing());
        REQUIRE(op4.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        SampleSummarizer<double> internal_2_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;
        std::vector<double> internal_2_height_samples;

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 9);
        std::vector<SampleSummarizer<double> > pop_size_summaries(pop_sizes.size());

        // burnin
        unsigned int n_moves_per_op = 3;
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, n_moves_per_op);
            op3.operate(rng, &tree, 1, n_moves_per_op);
            op4.operate(rng, &tree, 1, n_moves_per_op);
        }

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        std::shared_ptr<PopulationNode> node;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, n_moves_per_op);
            op3.operate(rng, &tree, 1, n_moves_per_op);
            op4.operate(rng, &tree, 1, n_moves_per_op);
            if ((i + 1) % sample_freq == 0) {
                node = tree.get_node("internal0");
                internal_0_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_0_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->has_parent());
                REQUIRE(node->get_parent()->get_label() == "internal1");
                node = tree.get_node("internal1");
                internal_1_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_1_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");
                node = tree.get_node("internal2");
                internal_2_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_2_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");

                pop_sizes = tree.get_pointers_to_population_sizes();
                REQUIRE(pop_sizes.size() == 9);
                for (unsigned int i = 0; i < pop_sizes.size(); ++i) {
                    pop_size_summaries.at(i).add_sample(pop_sizes.at(i)->get_value());
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        double ess_2 = effective_sample_size<double>(internal_2_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        std::cout << "ESS 2: " << ess_2 << "\n";
        
        BetaDistribution prior(4.0, 3.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
        std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";

        REQUIRE(op.get_number_of_attempts() == (niterations + burnin) * n_moves_per_op);
        REQUIRE(op.get_number_of_attempts_for_correction() == (((niterations + burnin) * n_moves_per_op) - 100));
        REQUIRE(op3.get_number_of_attempts() == (niterations + burnin) * n_moves_per_op);
        REQUIRE(op3.get_number_of_attempts_for_correction() == (((niterations + burnin) * n_moves_per_op) - 100));
        REQUIRE(op4.get_number_of_attempts() == (niterations + burnin) * n_moves_per_op);
        REQUIRE(op4.get_number_of_attempts_for_correction() == (((niterations + burnin) * n_moves_per_op) - 100));

        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);
        REQUIRE(internal_2_height_summary.sample_size() == nsamples);

        double eps = 0.001;
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        for (unsigned int i = 0; i < pop_size_summaries.size(); ++i) {
            REQUIRE(pop_size_summaries.at(i).sample_size() == nsamples);
            REQUIRE(pop_size_summaries.at(i).mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
            REQUIRE(pop_size_summaries.at(i).variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        }
    }
}

TEST_CASE("Testing HeightSizeSlideBumpMixer with 5 leaf balanced and optimizing",
        "[HeightSizeSlideBumpMixer]") {

    SECTION("Testing 5 leaf balanced optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(97458389374);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double pop_size_shape = 20.0;
        double pop_size_scale = 0.02;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(8, "root", 0.5);
        std::shared_ptr<PopulationNode> internal2 = std::make_shared<PopulationNode>(7, "internal2", 0.4);
        std::shared_ptr<PopulationNode> internal1 = std::make_shared<PopulationNode>(6, "internal1", 0.3);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(5, "internal0", 0.2);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);
        std::shared_ptr<PopulationNode> leaf3 = std::make_shared<PopulationNode>(3, "leaf3", 0.0);
        std::shared_ptr<PopulationNode> leaf4 = std::make_shared<PopulationNode>(4, "leaf4", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        internal2->add_child(leaf3);
        internal2->add_child(leaf4);
        root->add_child(internal1);
        root->add_child(internal2);

        BasePopulationTree tree(root);
        tree.set_root_node_height_prior(root_height_prior);
        tree.set_population_size_prior(pop_size_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(7.0);
        tree.set_beta_of_node_height_beta_prior(5.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        HeightSizeSlideBumpMixer op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);
        NodeHeightScaler<PopulationNode> op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);
        PopSizeScaler op4;
        op4.turn_on_auto_optimize();
        op4.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(op3.auto_optimizing());
        REQUIRE(op3.get_auto_optimize_delay() == 100);
        REQUIRE(op4.auto_optimizing());
        REQUIRE(op4.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        SampleSummarizer<double> internal_2_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;
        std::vector<double> internal_2_height_samples;

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 9);
        std::map<std::string, SampleSummarizer<double> > pop_size_summaries;

        // burnin
        unsigned int n_moves_per_op = 4;
        unsigned int burnin = 10000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, n_moves_per_op);
            op3.operate(rng, &tree, 1, n_moves_per_op);
            op4.operate(rng, &tree, 1, n_moves_per_op);
        }

        unsigned int niterations = 500000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        std::shared_ptr<PopulationNode> node;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, n_moves_per_op);
            op3.operate(rng, &tree, 1, n_moves_per_op);
            op4.operate(rng, &tree, 1, n_moves_per_op);
            if ((i + 1) % sample_freq == 0) {
                node = tree.get_node("internal0");
                internal_0_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_0_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                pop_size_summaries["internal0"].add_sample(node->get_population_size());
                REQUIRE(node->get_parent()->has_parent());
                REQUIRE(node->get_parent()->get_label() == "internal1");
                node = tree.get_node("internal1");
                internal_1_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_1_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                pop_size_summaries["internal1"].add_sample(node->get_population_size());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");
                node = tree.get_node("internal2");
                internal_2_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_2_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                pop_size_summaries["internal2"].add_sample(node->get_population_size());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");

                node = tree.get_node("root");
                pop_size_summaries["root"].add_sample(node->get_population_size());
                REQUIRE(node->get_height() == 0.5);
                node = tree.get_node("leaf0");
                pop_size_summaries["leaf0"].add_sample(node->get_population_size());
                node = tree.get_node("leaf1");
                pop_size_summaries["leaf1"].add_sample(node->get_population_size());
                node = tree.get_node("leaf2");
                pop_size_summaries["leaf2"].add_sample(node->get_population_size());
                node = tree.get_node("leaf3");
                pop_size_summaries["leaf3"].add_sample(node->get_population_size());
                node = tree.get_node("leaf4");
                pop_size_summaries["leaf4"].add_sample(node->get_population_size());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        double ess_2 = effective_sample_size<double>(internal_2_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        std::cout << "ESS 2: " << ess_2 << "\n";
        
        BetaDistribution prior(7.0, 5.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
        std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";
        std::cout << "\n";
        std::cout << "Expected pop size mean: " << pop_size_prior->get_mean() << "\n";
        std::cout << "Expected pop size var: " << pop_size_prior->get_variance() << "\n";
        for (auto k_v : pop_size_summaries) {
            std::string l = k_v.first;
            std::cout << "Pop size " << l << " mean: " << k_v.second.mean() << "\n";
            std::cout << "Pop size " << l << " var: " << k_v.second.variance() << "\n";
        }

        REQUIRE(op.get_number_of_attempts() == (niterations + burnin) * n_moves_per_op);
        REQUIRE(op.get_number_of_attempts_for_correction() == (((niterations + burnin) * n_moves_per_op) - 100));
        REQUIRE(op3.get_number_of_attempts() == (niterations + burnin) * n_moves_per_op);
        REQUIRE(op3.get_number_of_attempts_for_correction() == (((niterations + burnin) * n_moves_per_op) - 100));
        REQUIRE(op4.get_number_of_attempts() == (niterations + burnin) * n_moves_per_op);
        REQUIRE(op4.get_number_of_attempts_for_correction() == (((niterations + burnin) * n_moves_per_op) - 100));

        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);
        REQUIRE(internal_2_height_summary.sample_size() == nsamples);

        double eps = 0.002;
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        for (auto k_v : pop_size_summaries) {
            REQUIRE(k_v.second.sample_size() == nsamples);
            REQUIRE(k_v.second.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
            REQUIRE(k_v.second.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        }
    }
}

TEST_CASE("Testing HeightSizeSlideBumpMixer with 5 leaf ladder and optimizing",
        "[HeightSizeSlideBumpMixer]") {

    SECTION("Testing 5 leaf ladder optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(4254233651);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double pop_size_shape = 20.0;
        double pop_size_scale = 0.02;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(8, "root", 0.5);
        std::shared_ptr<PopulationNode> internal2 = std::make_shared<PopulationNode>(7, "internal2", 0.4);
        std::shared_ptr<PopulationNode> internal1 = std::make_shared<PopulationNode>(6, "internal1", 0.3);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(5, "internal0", 0.2);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);
        std::shared_ptr<PopulationNode> leaf3 = std::make_shared<PopulationNode>(3, "leaf3", 0.0);
        std::shared_ptr<PopulationNode> leaf4 = std::make_shared<PopulationNode>(4, "leaf4", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        internal2->add_child(internal1);
        internal2->add_child(leaf3);
        root->add_child(internal2);
        root->add_child(leaf4);

        BasePopulationTree tree(root);
        tree.set_root_node_height_prior(root_height_prior);
        tree.set_population_size_prior(pop_size_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(7.0);
        tree.set_beta_of_node_height_beta_prior(5.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        HeightSizeSlideBumpMixer op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);
        NodeHeightScaler<PopulationNode> op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);
        PopSizeScaler op4;
        op4.turn_on_auto_optimize();
        op4.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(op3.auto_optimizing());
        REQUIRE(op3.get_auto_optimize_delay() == 100);
        REQUIRE(op4.auto_optimizing());
        REQUIRE(op4.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        SampleSummarizer<double> internal_2_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;
        std::vector<double> internal_2_height_samples;

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 9);
        std::map<std::string, SampleSummarizer<double> > pop_size_summaries;

        // burnin
        unsigned int n_moves_per_op = 4;
        unsigned int burnin = 10000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, n_moves_per_op);
            op3.operate(rng, &tree, 1, n_moves_per_op);
            op4.operate(rng, &tree, 1, n_moves_per_op);
        }

        unsigned int niterations = 500000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        std::shared_ptr<PopulationNode> node;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, n_moves_per_op);
            op3.operate(rng, &tree, 1, n_moves_per_op);
            op4.operate(rng, &tree, 1, n_moves_per_op);
            if ((i + 1) % sample_freq == 0) {
                node = tree.get_node("internal0");
                internal_0_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_0_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                pop_size_summaries["internal0"].add_sample(node->get_population_size());
                REQUIRE(node->get_parent()->has_parent());
                REQUIRE(node->get_parent()->get_label() == "internal1");
                node = tree.get_node("internal1");
                internal_1_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_1_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                pop_size_summaries["internal1"].add_sample(node->get_population_size());
                REQUIRE(node->get_parent()->has_parent());
                REQUIRE(node->get_parent()->get_label() == "internal2");
                node = tree.get_node("internal2");
                internal_2_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_2_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                pop_size_summaries["internal2"].add_sample(node->get_population_size());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");

                node = tree.get_node("root");
                pop_size_summaries["root"].add_sample(node->get_population_size());
                REQUIRE(node->get_height() == 0.5);
                node = tree.get_node("leaf0");
                pop_size_summaries["leaf0"].add_sample(node->get_population_size());
                node = tree.get_node("leaf1");
                pop_size_summaries["leaf1"].add_sample(node->get_population_size());
                node = tree.get_node("leaf2");
                pop_size_summaries["leaf2"].add_sample(node->get_population_size());
                node = tree.get_node("leaf3");
                pop_size_summaries["leaf3"].add_sample(node->get_population_size());
                node = tree.get_node("leaf4");
                pop_size_summaries["leaf4"].add_sample(node->get_population_size());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        double ess_2 = effective_sample_size<double>(internal_2_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        std::cout << "ESS 2: " << ess_2 << "\n";
        
        BetaDistribution prior(7.0, 5.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "Internal 2 mean: " << internal_2_height_summary.mean() << "\n";
        std::cout << "Internal 2 var: "  << internal_2_height_summary.variance() << "\n";
        std::cout << "\n";
        std::cout << "Expected pop size mean: " << pop_size_prior->get_mean() << "\n";
        std::cout << "Expected pop size var: " << pop_size_prior->get_variance() << "\n";
        for (auto k_v : pop_size_summaries) {
            std::string l = k_v.first;
            std::cout << "Pop size " << l << " mean: " << k_v.second.mean() << "\n";
            std::cout << "Pop size " << l << " var: " << k_v.second.variance() << "\n";
        }

        REQUIRE(op.get_number_of_attempts() == (niterations + burnin) * n_moves_per_op);
        REQUIRE(op.get_number_of_attempts_for_correction() == (((niterations + burnin) * n_moves_per_op) - 100));
        REQUIRE(op3.get_number_of_attempts() == (niterations + burnin) * n_moves_per_op);
        REQUIRE(op3.get_number_of_attempts_for_correction() == (((niterations + burnin) * n_moves_per_op) - 100));
        REQUIRE(op4.get_number_of_attempts() == (niterations + burnin) * n_moves_per_op);
        REQUIRE(op4.get_number_of_attempts_for_correction() == (((niterations + burnin) * n_moves_per_op) - 100));

        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);
        REQUIRE(internal_2_height_summary.sample_size() == nsamples);

        double eps = 0.002;
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_2_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_2_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        for (auto k_v : pop_size_summaries) {
            REQUIRE(k_v.second.sample_size() == nsamples);
            REQUIRE(k_v.second.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
            REQUIRE(k_v.second.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        }
    }
}

TEST_CASE("Testing HeightSizeSlideBumpMixer with 4 leaf ladder and optimizing",
        "[HeightSizeSlideBumpMixer]") {

    SECTION("Testing 4 leaf ladder optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(973984247);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double pop_size_shape = 20.0;
        double pop_size_scale = 0.02;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(6, "root", 0.5);
        std::shared_ptr<PopulationNode> internal1 = std::make_shared<PopulationNode>(5, "internal1", 0.3);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(4, "internal0", 0.2);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);
        std::shared_ptr<PopulationNode> leaf3 = std::make_shared<PopulationNode>(3, "leaf3", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(internal0);
        internal1->add_child(leaf2);
        root->add_child(internal1);
        root->add_child(leaf3);

        BasePopulationTree tree(root);
        tree.set_root_node_height_prior(root_height_prior);
        tree.set_population_size_prior(pop_size_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(7.0);
        tree.set_beta_of_node_height_beta_prior(5.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        HeightSizeSlideBumpMixer op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);
        NodeHeightScaler<PopulationNode> op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);
        PopSizeScaler op4;
        op4.turn_on_auto_optimize();
        op4.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(op3.auto_optimizing());
        REQUIRE(op3.get_auto_optimize_delay() == 100);
        REQUIRE(op4.auto_optimizing());
        REQUIRE(op4.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 7);
        std::map<std::string, SampleSummarizer<double> > pop_size_summaries;

        // burnin
        unsigned int n_moves_per_op = 4;
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, n_moves_per_op);
            op3.operate(rng, &tree, 1, n_moves_per_op);
            op4.operate(rng, &tree, 1, n_moves_per_op);
        }

        unsigned int niterations = 100000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        std::shared_ptr<PopulationNode> node;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, n_moves_per_op);
            op3.operate(rng, &tree, 1, n_moves_per_op);
            op4.operate(rng, &tree, 1, n_moves_per_op);
            if ((i + 1) % sample_freq == 0) {
                node = tree.get_node("internal0");
                internal_0_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_0_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                pop_size_summaries["internal0"].add_sample(node->get_population_size());
                REQUIRE(node->get_parent()->has_parent());
                REQUIRE(node->get_parent()->get_label() == "internal1");
                node = tree.get_node("internal1");
                internal_1_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_1_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                pop_size_summaries["internal1"].add_sample(node->get_population_size());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");

                node = tree.get_node("root");
                pop_size_summaries["root"].add_sample(node->get_population_size());
                REQUIRE(node->get_height() == 0.5);
                node = tree.get_node("leaf0");
                pop_size_summaries["leaf0"].add_sample(node->get_population_size());
                node = tree.get_node("leaf1");
                pop_size_summaries["leaf1"].add_sample(node->get_population_size());
                node = tree.get_node("leaf2");
                pop_size_summaries["leaf2"].add_sample(node->get_population_size());
                node = tree.get_node("leaf3");
                pop_size_summaries["leaf3"].add_sample(node->get_population_size());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        
        BetaDistribution prior(7.0, 5.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "\n";
        std::cout << "Expected pop size mean: " << pop_size_prior->get_mean() << "\n";
        std::cout << "Expected pop size var: " << pop_size_prior->get_variance() << "\n";
        for (auto k_v : pop_size_summaries) {
            std::string l = k_v.first;
            std::cout << "Pop size " << l << " mean: " << k_v.second.mean() << "\n";
            std::cout << "Pop size " << l << " var: " << k_v.second.variance() << "\n";
        }

        REQUIRE(op.get_number_of_attempts() == (niterations + burnin) * n_moves_per_op);
        REQUIRE(op.get_number_of_attempts_for_correction() == (((niterations + burnin) * n_moves_per_op) - 100));
        REQUIRE(op3.get_number_of_attempts() == (niterations + burnin) * n_moves_per_op);
        REQUIRE(op3.get_number_of_attempts_for_correction() == (((niterations + burnin) * n_moves_per_op) - 100));
        REQUIRE(op4.get_number_of_attempts() == (niterations + burnin) * n_moves_per_op);
        REQUIRE(op4.get_number_of_attempts_for_correction() == (((niterations + burnin) * n_moves_per_op) - 100));

        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);

        double eps = 0.002;
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        for (auto k_v : pop_size_summaries) {
            REQUIRE(k_v.second.sample_size() == nsamples);
            REQUIRE(k_v.second.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
            REQUIRE(k_v.second.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        }
    }
}

TEST_CASE("Testing HeightSizeSlideBumpMixer with 4 leaf balanced and optimizing",
        "[HeightSizeSlideBumpMixer]") {

    SECTION("Testing 4 leaf balanced optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(879587345);

        double root_height_shape = 10.0;
        double root_height_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<GammaDistribution>(
                root_height_shape,
                root_height_scale);

        double pop_size_shape = 20.0;
        double pop_size_scale = 0.02;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(6, "root", 0.5);
        std::shared_ptr<PopulationNode> internal1 = std::make_shared<PopulationNode>(5, "internal1", 0.3);
        std::shared_ptr<PopulationNode> internal0 = std::make_shared<PopulationNode>(4, "internal0", 0.2);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);
        std::shared_ptr<PopulationNode> leaf3 = std::make_shared<PopulationNode>(3, "leaf3", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(leaf3);
        root->add_child(internal0);
        root->add_child(internal1);

        BasePopulationTree tree(root);
        tree.set_root_node_height_prior(root_height_prior);
        tree.set_population_size_prior(pop_size_prior);

        tree.estimate_alpha_of_node_height_beta_prior();
        tree.estimate_beta_of_node_height_beta_prior();
        tree.set_alpha_of_node_height_beta_prior(7.0);
        tree.set_beta_of_node_height_beta_prior(5.0);
        tree.fix_alpha_of_node_height_beta_prior();
        tree.fix_beta_of_node_height_beta_prior();

        tree.ignore_data();
        tree.estimate_root_height();

        HeightSizeSlideBumpMixer op;
        op.turn_on_auto_optimize();
        op.set_auto_optimize_delay(100);
        NodeHeightScaler<PopulationNode> op3;
        op3.turn_on_auto_optimize();
        op3.set_auto_optimize_delay(100);
        PopSizeScaler op4;
        op4.turn_on_auto_optimize();
        op4.set_auto_optimize_delay(100);

        REQUIRE(op.auto_optimizing());
        REQUIRE(op.get_auto_optimize_delay() == 100);
        REQUIRE(op3.auto_optimizing());
        REQUIRE(op3.get_auto_optimize_delay() == 100);
        REQUIRE(op4.auto_optimizing());
        REQUIRE(op4.get_auto_optimize_delay() == 100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_0_height_summary;
        SampleSummarizer<double> internal_1_height_summary;
        std::vector<double> internal_0_height_samples;
        std::vector<double> internal_1_height_samples;

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes = tree.get_pointers_to_population_sizes();
        REQUIRE(pop_sizes.size() == 7);
        std::map<std::string, SampleSummarizer<double> > pop_size_summaries;

        // burnin
        unsigned int n_moves_per_op = 4;
        unsigned int burnin = 1000;
        for (unsigned int i = 0; i < burnin; ++i) {
            op.operate(rng, &tree, 1, n_moves_per_op);
            op3.operate(rng, &tree, 1, n_moves_per_op);
            op4.operate(rng, &tree, 1, n_moves_per_op);
        }

        unsigned int niterations = 100000;
        unsigned int sample_freq = 5;
        unsigned int nsamples = niterations / sample_freq;
        std::shared_ptr<PopulationNode> node;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, n_moves_per_op);
            op3.operate(rng, &tree, 1, n_moves_per_op);
            op4.operate(rng, &tree, 1, n_moves_per_op);
            if ((i + 1) % sample_freq == 0) {
                node = tree.get_node("internal0");
                internal_0_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_0_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                pop_size_summaries["internal0"].add_sample(node->get_population_size());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");
                node = tree.get_node("internal1");
                internal_1_height_summary.add_sample(node->get_height() / node->get_parent()->get_height());
                internal_1_height_samples.push_back(node->get_height() / node->get_parent()->get_height());
                pop_size_summaries["internal1"].add_sample(node->get_population_size());
                REQUIRE(node->get_parent()->is_root());
                REQUIRE(node->get_parent()->get_label() == "root");

                node = tree.get_node("root");
                pop_size_summaries["root"].add_sample(node->get_population_size());
                REQUIRE(node->get_height() == 0.5);
                node = tree.get_node("leaf0");
                pop_size_summaries["leaf0"].add_sample(node->get_population_size());
                node = tree.get_node("leaf1");
                pop_size_summaries["leaf1"].add_sample(node->get_population_size());
                node = tree.get_node("leaf2");
                pop_size_summaries["leaf2"].add_sample(node->get_population_size());
                node = tree.get_node("leaf3");
                pop_size_summaries["leaf3"].add_sample(node->get_population_size());
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        double ess_0 = effective_sample_size<double>(internal_0_height_samples);
        double ess_1 = effective_sample_size<double>(internal_1_height_samples);
        std::cout << "ESS 0: " << ess_0 << "\n";
        std::cout << "ESS 1: " << ess_1 << "\n";
        
        BetaDistribution prior(7.0, 5.0);

        std::cout << "Expected mean: " << prior.get_mean() << "\n";
        std::cout << "Expected var: " << prior.get_variance() << "\n";
        std::cout << "Internal 0 mean: " << internal_0_height_summary.mean() << "\n";
        std::cout << "Internal 0 var: "  << internal_0_height_summary.variance() << "\n";
        std::cout << "Internal 1 mean: " << internal_1_height_summary.mean() << "\n";
        std::cout << "Internal 1 var: "  << internal_1_height_summary.variance() << "\n";
        std::cout << "\n";
        std::cout << "Expected pop size mean: " << pop_size_prior->get_mean() << "\n";
        std::cout << "Expected pop size var: " << pop_size_prior->get_variance() << "\n";
        for (auto k_v : pop_size_summaries) {
            std::string l = k_v.first;
            std::cout << "Pop size " << l << " mean: " << k_v.second.mean() << "\n";
            std::cout << "Pop size " << l << " var: " << k_v.second.variance() << "\n";
        }

        REQUIRE(op.get_number_of_attempts() == (niterations + burnin) * n_moves_per_op);
        REQUIRE(op.get_number_of_attempts_for_correction() == (((niterations + burnin) * n_moves_per_op) - 100));
        REQUIRE(op3.get_number_of_attempts() == (niterations + burnin) * n_moves_per_op);
        REQUIRE(op3.get_number_of_attempts_for_correction() == (((niterations + burnin) * n_moves_per_op) - 100));
        REQUIRE(op4.get_number_of_attempts() == (niterations + burnin) * n_moves_per_op);
        REQUIRE(op4.get_number_of_attempts_for_correction() == (((niterations + burnin) * n_moves_per_op) - 100));

        REQUIRE(internal_0_height_summary.sample_size() == nsamples);
        REQUIRE(internal_1_height_summary.sample_size() == nsamples);

        double eps = 0.002;
        REQUIRE(internal_0_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_0_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        REQUIRE(internal_1_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_1_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
        for (auto k_v : pop_size_summaries) {
            REQUIRE(k_v.second.sample_size() == nsamples);
            REQUIRE(k_v.second.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
            REQUIRE(k_v.second.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));
        }
    }
}
