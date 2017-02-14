#include "catch.hpp"
#include "ecoevolity/operator.hpp"

#include <limits>
#include <memory>

#include "ecoevolity/probability.hpp"
#include "ecoevolity/parameter.hpp"
#include "ecoevolity/rng.hpp"

TEST_CASE("Testing ComparisonHeightScaler", "[ComparisonHeightScaler]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(123);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ComparisonHeightScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);
        PositiveRealParameter node_height = PositiveRealParameter(prior, 1.0, false);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            node_height.store();
            double hastings = o.propose(rng, node_height);
            double prior_ratio = node_height.relative_prior_ln_pdf() -
                node_height.relative_prior_ln_pdf(node_height.get_stored_value());
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                node_height.restore();
            }
            o.optimize(os, acceptance_prob);
            mn = std::min(mn, node_height.get_value());
            mx = std::max(mx, node_height.get_value());
            ++n;
            d = node_height.get_value() - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
        
    }

    SECTION("Testing gamma(10.0, 0.1) prior and with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(123);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ComparisonHeightScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);
        PositiveRealParameter node_height = PositiveRealParameter(prior, 1.0, false);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            node_height.store();
            double hastings = o.propose(rng, node_height);
            double prior_ratio = node_height.relative_prior_ln_pdf() -
                node_height.relative_prior_ln_pdf(node_height.get_stored_value());
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                node_height.restore();
            }
            o.optimize(os, acceptance_prob);
            mn = std::min(mn, node_height.get_value());
            mx = std::max(mx, node_height.get_value());
            ++n;
            d = node_height.get_value() - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
        
    }
}

TEST_CASE("Testing ComparisonHeightMover", "[ComparisonHeightMover]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(123);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ComparisonHeightMover>(1.0, 0.2);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);
        PositiveRealParameter node_height = PositiveRealParameter(prior, 1.0, false);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            node_height.store();
            double hastings = o.propose(rng, node_height);
            double prior_ratio = node_height.relative_prior_ln_pdf() -
                node_height.relative_prior_ln_pdf(node_height.get_stored_value());
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                node_height.restore();
            }
            o.optimize(os, acceptance_prob);
            mn = std::min(mn, node_height.get_value());
            mx = std::max(mx, node_height.get_value());
            ++n;
            d = node_height.get_value() - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
        
    }

    SECTION("Testing gamma(10.0, 0.1) prior and with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(123);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ComparisonHeightMover>(1.0, 0.2);
        OperatorSchedule os = OperatorSchedule();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);
        PositiveRealParameter node_height = PositiveRealParameter(prior, 1.0, false);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            node_height.store();
            double hastings = o.propose(rng, node_height);
            double prior_ratio = node_height.relative_prior_ln_pdf() -
                node_height.relative_prior_ln_pdf(node_height.get_stored_value());
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                node_height.restore();
            }
            o.optimize(os, acceptance_prob);
            mn = std::min(mn, node_height.get_value());
            mx = std::max(mx, node_height.get_value());
            ++n;
            d = node_height.get_value() - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
        
    }
}


TEST_CASE("Testing RootPopulationSizeScaler", "[RootPopulationSizeScaler]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(12933);
        std::shared_ptr<OperatorInterface> op = std::make_shared<RootPopulationSizeScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_root_population_size(1.0);
        tree.set_population_size_prior(prior);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            // REQUIRE(tree.is_dirty());
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_root_population_size();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
        
    }

    SECTION("Testing gamma(10.0, 0.1) prior with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(12933);
        std::shared_ptr<OperatorInterface> op = std::make_shared<RootPopulationSizeScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        //os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_root_population_size(1.0);
        tree.set_population_size_prior(prior);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            // REQUIRE(tree.is_dirty());
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_root_population_size();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
        
    }
}

TEST_CASE("Testing ChildPopulationSizeScaler", "[ChildPopulationSizeScaler]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(12933);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ChildPopulationSizeScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.constrain_population_sizes();
        tree.set_child_population_size(0, 1.0);
        tree.set_population_size_prior(prior);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            // REQUIRE(tree.is_dirty());
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_child_population_size(0);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
    }

    SECTION("Testing gamma(10.0, 0.1) prior with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(12933);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ChildPopulationSizeScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        // os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.constrain_population_sizes();
        tree.set_child_population_size(0, 1.0);
        tree.set_population_size_prior(prior);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 1000000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            // REQUIRE(tree.is_dirty());
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_child_population_size(0);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
    }
}

TEST_CASE("Testing ComparisonMutationRateScaler", "[ComparisonMutationRateScaler]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(928374);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ComparisonMutationRateScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_state_frequencies();
        tree.fix_population_sizes();
        tree.set_mutation_rate(1.0);
        tree.set_mutation_rate_prior(prior);
        tree.estimate_mutation_rate();
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double old_v = tree.get_mutation_rate();
            double hastings = o.propose(rng, tree);
            double new_v = tree.get_mutation_rate();
            REQUIRE(tree.is_dirty());
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            /* double prior_ratio = */
            /*     prior->relative_ln_pdf(new_v) - */
            /*     prior->relative_ln_pdf(old_v); */
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_mutation_rate();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
        
    }

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(928374);
        std::shared_ptr<OperatorInterface> op = std::make_shared<ComparisonMutationRateScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        // os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_state_frequencies();
        tree.fix_population_sizes();
        tree.set_mutation_rate(1.0);
        tree.set_mutation_rate_prior(prior);
        tree.estimate_mutation_rate();
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double old_v = tree.get_mutation_rate();
            double hastings = o.propose(rng, tree);
            double new_v = tree.get_mutation_rate();
            REQUIRE(tree.is_dirty());
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            /* double prior_ratio = */
            /*     prior->relative_ln_pdf(new_v) - */
            /*     prior->relative_ln_pdf(old_v); */
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_mutation_rate();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
        
    }
}

TEST_CASE("Testing FreqMover", "[FreqMover]") {

    SECTION("testing beta(1.0, 1.0) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(3648);
        std::shared_ptr<OperatorInterface> op = std::make_shared<FreqMover>(1.0, 0.1);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        double a = 1.0;
        double b = 1.0;
        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<BetaDistribution>(a, b);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_population_sizes();
        tree.fix_mutation_rate();
        tree.estimate_state_frequencies();
        tree.set_freq_1_prior(prior);
        tree.set_freq_1(0.5);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(! tree.state_frequencies_are_constrained());
        REQUIRE(! tree.state_frequencies_are_fixed());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_freq_1();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.005));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
    }

    SECTION("testing beta(1.0, 1.0) prior with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(2945720);
        std::shared_ptr<OperatorInterface> op = std::make_shared<FreqMover>(1.0, 0.1);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        double a = 1.0;
        double b = 1.0;
        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<BetaDistribution>(a, b);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_population_sizes();
        tree.fix_mutation_rate();
        tree.estimate_state_frequencies();
        tree.set_freq_1_prior(prior);
        tree.set_freq_1(0.5);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(! tree.state_frequencies_are_constrained());
        REQUIRE(! tree.state_frequencies_are_fixed());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_freq_1();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.005));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
    }

    SECTION("testing beta(5.0, 1.0) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(841984264);
        std::shared_ptr<OperatorInterface> op = std::make_shared<FreqMover>(1.0, 0.1);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        double a = 5.0;
        double b = 1.0;
        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<BetaDistribution>(a, b);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_population_sizes();
        tree.fix_mutation_rate();
        tree.estimate_state_frequencies();
        tree.set_freq_1_prior(prior);
        tree.set_freq_1(0.5);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(! tree.state_frequencies_are_constrained());
        REQUIRE(! tree.state_frequencies_are_fixed());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_freq_1();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.005));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
    }

    SECTION("testing beta(5.0, 1.0) prior with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(25456657);
        std::shared_ptr<OperatorInterface> op = std::make_shared<FreqMover>(1.0, 0.1);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        double a = 5.0;
        double b = 1.0;
        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<BetaDistribution>(a, b);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_population_sizes();
        tree.fix_mutation_rate();
        tree.estimate_state_frequencies();
        tree.set_freq_1_prior(prior);
        tree.set_freq_1(0.5);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(! tree.state_frequencies_are_constrained());
        REQUIRE(! tree.state_frequencies_are_fixed());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_freq_1();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.005));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
    }

    SECTION("testing beta(1.0, 5.0) prior with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(14458);
        std::shared_ptr<OperatorInterface> op = std::make_shared<FreqMover>(1.0, 0.1);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        double a = 1.0;
        double b = 5.0;
        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<BetaDistribution>(a, b);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, ' ', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_population_sizes();
        tree.fix_mutation_rate();
        tree.estimate_state_frequencies();
        tree.set_freq_1_prior(prior);
        tree.set_freq_1(0.5);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(! tree.state_frequencies_are_constrained());
        REQUIRE(! tree.state_frequencies_are_fixed());
    
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            OperatorInterface& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            tree.compute_log_likelihood_and_prior();
            double prior_ratio = tree.get_log_prior_density_value() -
                tree.get_stored_log_prior_density_value();
            double acceptance_prob = prior_ratio + hastings;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_prob)) {
                o.accept(os);
            }
            else {
                o.reject(os);
                tree.restore_state();
            }
            o.optimize(os, acceptance_prob);
            double x = tree.get_freq_1();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        std::cout << op->header_string();
        std::cout << op->to_string(os);
        std::cout << "prior mean: " << mean << "\n";
        std::cout << "prior variance: " << variance << "\n";
        std::cout << "expected prior mean: " << prior->get_mean() << "\n";
        std::cout << "expected prior variance: " << prior->get_variance() << "\n";
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.005));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
    }
}
