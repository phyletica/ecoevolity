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
        std::shared_ptr<Operator> op = std::make_shared<ComparisonHeightScaler>(1.0, 0.5);
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
            Operator& o = os.draw_operator(rng);
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
        std::shared_ptr<Operator> op = std::make_shared<ComparisonHeightScaler>(1.0, 0.5);
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
            Operator& o = os.draw_operator(rng);
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

    SECTION("Testing tree with gamma(10.0, 0.1) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(123);
        std::shared_ptr<Operator> op = std::make_shared<ComparisonHeightScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);
        std::shared_ptr<PositiveRealParameter> node_height = std::make_shared<PositiveRealParameter>(prior, 1.0, false);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, '_', true, true, false);
        tree.set_node_height_prior(prior);
        tree.set_population_size_prior(std::make_shared<GammaDistribution>(1.0, 0.001));
        tree.set_height_parameter(node_height);
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
            Operator& o = os.draw_operator(rng);
            node_height->store();
            tree.store_state();
            double hastings = o.propose(rng, *node_height);
            tree.make_dirty();
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
                node_height->restore();
                tree.restore_likelihood();
                tree.restore_prior_density();
            }
            o.optimize(os, acceptance_prob);

            //REQUIRE(node_height->get_value() == Approx(tree.get_root_height()));
            mn = std::min(mn, node_height->get_value());
            mx = std::max(mx, node_height->get_value());
            ++n;
            d = node_height->get_value() - mean;
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

TEST_CASE("Testing RootCoalescenceRateScaler", "[RootCoalescenceRateScaler]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(123);
        std::shared_ptr<Operator> op = std::make_shared<RootCoalescenceRateScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, '_', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.set_root_coalescence_rate(1.0);
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
            Operator& o = os.draw_operator(rng);
            tree.store_state();
            double hastings = o.propose(rng, tree);
            REQUIRE(tree.is_dirty());
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
            double x = CoalescenceRateParameter::get_population_size_from_rate(tree.get_root_coalescence_rate());
            //double x = tree.get_root_coalescence_rate();
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
        
        REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
        REQUIRE(mn >= prior->get_min());
        REQUIRE(mx < prior->get_max());
        //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
        //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
        
    }
}
