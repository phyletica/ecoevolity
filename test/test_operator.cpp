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
}

TEST_CASE("Testing ComparisonHeightMover", "[ComparisonHeightMover]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(123);
        std::shared_ptr<Operator> op = std::make_shared<ComparisonHeightMover>(1.0, 0.2);
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
        std::shared_ptr<Operator> op = std::make_shared<ComparisonHeightMover>(1.0, 0.2);
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
}

// The 'CoalescenceRateScaler' operators are not working, and this is
// documented in the code. Commenting-out these tests, because these operators
// are no longer being used in the package.
// TEST_CASE("Testing RootCoalescenceRateScaler", "[RootCoalescenceRateScaler]") {
// 
//     SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
//         RandomNumberGenerator rng = RandomNumberGenerator(1292733);
//         std::shared_ptr<Operator> op = std::make_shared<RootCoalescenceRateScaler>(1.0, 0.5);
//         OperatorSchedule os = OperatorSchedule();
//         os.turn_off_auto_optimize();
//         // os.turn_on_auto_optimize();
//         // os.set_auto_optimize_delay(10000);
//         os.add_operator(op);
// 
//         std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);
// 
//         std::string nex_path = "data/hemi129.nex";
//         ComparisonPopulationTree tree(nex_path, '_', true, true, false);
//         REQUIRE(tree.get_degree_of_root() == 2);
//         tree.set_root_height(0.01);
//         tree.set_root_coalescence_rate(1.0);
//         tree.set_population_size_prior(prior);
//         tree.ignore_data();
// 
//         tree.make_dirty();
//         tree.compute_log_likelihood_and_prior();
//         REQUIRE(! tree.is_dirty());
//     
//         unsigned int n = 0;
//         double mean = 0.0;
//         double sum_devs = 0.0;
//         double d;
//         double d_n;
//         double mn = std::numeric_limits<double>::max();
//         double mx = -std::numeric_limits<double>::max();
//         for (unsigned int i = 0; i < 100000; ++i) {
//             Operator& o = os.draw_operator(rng);
//             tree.store_state();
//             double old_v = tree.get_root_population_size();
//             double hastings = o.propose(rng, tree);
//             double new_v = tree.get_root_population_size();
//             REQUIRE(tree.is_dirty());
//             tree.compute_log_likelihood_and_prior();
//             double prior_ratio = tree.get_log_prior_density_value() -
//                 tree.get_stored_log_prior_density_value();
//             /* double prior_ratio = tree.get_root_coalescence_rate_parameter()->relative_prior_ln_pdf() - */
//             /*     tree.get_root_coalescence_rate_parameter()->relative_prior_ln_pdf( */
//             /*             tree.get_root_coalescence_rate_parameter()->get_stored_value()); */
//             /* double prior_ratio = */
//             /*     prior->relative_ln_pdf( */
//             /*             2.0/tree.get_root_coalescence_rate_parameter()->get_value()) - */
//             /*     prior->relative_ln_pdf( */
//             /*             2.0/tree.get_root_coalescence_rate_parameter()->get_stored_value()); */
//             /* double prior_ratio = */
//             /*     prior->relative_ln_pdf(new_v) - */
//             /*     prior->relative_ln_pdf(old_v); */
//             double acceptance_prob = prior_ratio + hastings;
//             double u = rng.uniform_real();
//             if (u < std::exp(acceptance_prob)) {
//                 o.accept(os);
//             }
//             else {
//                 o.reject(os);
//                 tree.restore_state();
//             }
//             o.optimize(os, acceptance_prob);
//             double x = tree.get_root_population_size();
//             mn = std::min(mn, x);
//             mx = std::max(mx, x);
//             ++n;
//             d = x - mean;
//             d_n = d / n;
//             mean += d_n;
//             sum_devs += d * d_n * (n - 1);
//         }
//         double variance = sum_devs / (n - 1);
//         std::cout << op->header_string();
//         std::cout << op->to_string(os);
//         std::cout << "prior mean: " << mean << "\n";
//         std::cout << "prior variance: " << variance << "\n";
//         std::cout << "expected prior mean: " << prior->get_mean() << "\n";
//         std::cout << "expected prior variance: " << prior->get_variance() << "\n";
//         
//         REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
//         REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
//         REQUIRE(mn >= prior->get_min());
//         REQUIRE(mx < prior->get_max());
//         //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
//         //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
//         
//     }
// 
//     SECTION("Testing gamma(10.0, 0.1) prior with optimizing") {
//         RandomNumberGenerator rng = RandomNumberGenerator(1292733);
//         std::shared_ptr<Operator> op = std::make_shared<RootCoalescenceRateScaler>(1.0, 0.5);
//         OperatorSchedule os = OperatorSchedule();
//         //os.turn_off_auto_optimize();
//         os.turn_on_auto_optimize();
//         os.set_auto_optimize_delay(10000);
//         os.add_operator(op);
// 
//         std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);
// 
//         std::string nex_path = "data/hemi129.nex";
//         ComparisonPopulationTree tree(nex_path, '_', true, true, false);
//         REQUIRE(tree.get_degree_of_root() == 2);
//         tree.set_root_height(0.01);
//         tree.set_root_coalescence_rate(1.0);
//         tree.set_population_size_prior(prior);
//         tree.ignore_data();
// 
//         tree.make_dirty();
//         tree.compute_log_likelihood_and_prior();
//         REQUIRE(! tree.is_dirty());
//     
//         unsigned int n = 0;
//         double mean = 0.0;
//         double sum_devs = 0.0;
//         double d;
//         double d_n;
//         double mn = std::numeric_limits<double>::max();
//         double mx = -std::numeric_limits<double>::max();
//         for (unsigned int i = 0; i < 100000; ++i) {
//             Operator& o = os.draw_operator(rng);
//             tree.store_state();
//             double old_v = tree.get_root_population_size();
//             double hastings = o.propose(rng, tree);
//             double new_v = tree.get_root_population_size();
//             REQUIRE(tree.is_dirty());
//             tree.compute_log_likelihood_and_prior();
//             double prior_ratio = tree.get_log_prior_density_value() -
//                 tree.get_stored_log_prior_density_value();
//             /* double prior_ratio = tree.get_root_coalescence_rate_parameter()->relative_prior_ln_pdf() - */
//             /*     tree.get_root_coalescence_rate_parameter()->relative_prior_ln_pdf( */
//             /*             tree.get_root_coalescence_rate_parameter()->get_stored_value()); */
//             /* double prior_ratio = */
//             /*     prior->relative_ln_pdf( */
//             /*             2.0/tree.get_root_coalescence_rate_parameter()->get_value()) - */
//             /*     prior->relative_ln_pdf( */
//             /*             2.0/tree.get_root_coalescence_rate_parameter()->get_stored_value()); */
//             /* double prior_ratio = */
//             /*     prior->relative_ln_pdf(new_v) - */
//             /*     prior->relative_ln_pdf(old_v); */
//             double acceptance_prob = prior_ratio + hastings;
//             double u = rng.uniform_real();
//             if (u < std::exp(acceptance_prob)) {
//                 o.accept(os);
//             }
//             else {
//                 o.reject(os);
//                 tree.restore_state();
//             }
//             o.optimize(os, acceptance_prob);
//             double x = tree.get_root_population_size();
//             mn = std::min(mn, x);
//             mx = std::max(mx, x);
//             ++n;
//             d = x - mean;
//             d_n = d / n;
//             mean += d_n;
//             sum_devs += d * d_n * (n - 1);
//         }
//         double variance = sum_devs / (n - 1);
//         std::cout << op->header_string();
//         std::cout << op->to_string(os);
//         std::cout << "prior mean: " << mean << "\n";
//         std::cout << "prior variance: " << variance << "\n";
//         std::cout << "expected prior mean: " << prior->get_mean() << "\n";
//         std::cout << "expected prior variance: " << prior->get_variance() << "\n";
//         
//         REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
//         REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
//         REQUIRE(mn >= prior->get_min());
//         REQUIRE(mx < prior->get_max());
//         //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
//         //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
//         
//     }
// }
// 
// TEST_CASE("Testing ChildCoalescenceRateScaler", "[ChildCoalescenceRateScaler]") {
// 
//     SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
//         RandomNumberGenerator rng = RandomNumberGenerator(928374);
//         std::shared_ptr<Operator> op = std::make_shared<ChildCoalescenceRateScaler>(1.0, 0.5);
//         OperatorSchedule os = OperatorSchedule();
//         os.turn_off_auto_optimize();
//         // os.turn_on_auto_optimize();
//         // os.set_auto_optimize_delay(10000);
//         os.add_operator(op);
// 
//         std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);
// 
//         std::string nex_path = "data/hemi129.nex";
//         ComparisonPopulationTree tree(nex_path, '_', true, true, false);
//         REQUIRE(tree.get_degree_of_root() == 2);
//         tree.set_root_height(0.01);
//         tree.set_root_coalescence_rate(1.0);
//         tree.set_population_size_prior(prior);
//         tree.constrain_coalescence_rates();
//         tree.ignore_data();
// 
//         tree.make_dirty();
//         tree.compute_log_likelihood_and_prior();
//         REQUIRE(! tree.is_dirty());
//     
//         unsigned int n = 0;
//         double mean = 0.0;
//         double sum_devs = 0.0;
//         double d;
//         double d_n;
//         double mn = std::numeric_limits<double>::max();
//         double mx = -std::numeric_limits<double>::max();
//         for (unsigned int i = 0; i < 100000; ++i) {
//             Operator& o = os.draw_operator(rng);
//             tree.store_state();
//             double old_v = tree.get_child_population_size(0);
//             double hastings = o.propose(rng, tree);
//             double new_v = tree.get_child_population_size(0);
//             REQUIRE(tree.is_dirty());
//             REQUIRE(new_v == tree.get_child_population_size(1));
//             REQUIRE(new_v == tree.get_root_population_size());
//             tree.compute_log_likelihood_and_prior();
//             double prior_ratio = tree.get_log_prior_density_value() -
//                 tree.get_stored_log_prior_density_value();
//             /* double prior_ratio = */
//             /*     prior->relative_ln_pdf(new_v) - */
//             /*     prior->relative_ln_pdf(old_v); */
//             double acceptance_prob = prior_ratio + hastings;
//             double u = rng.uniform_real();
//             if (u < std::exp(acceptance_prob)) {
//                 o.accept(os);
//             }
//             else {
//                 o.reject(os);
//                 tree.restore_state();
//                 REQUIRE(tree.get_child_population_size(0) == old_v);
//             }
//             o.optimize(os, acceptance_prob);
//             double x = tree.get_child_population_size(0);
//             mn = std::min(mn, x);
//             mx = std::max(mx, x);
//             ++n;
//             d = x - mean;
//             d_n = d / n;
//             mean += d_n;
//             sum_devs += d * d_n * (n - 1);
//         }
//         double variance = sum_devs / (n - 1);
//         std::cout << op->header_string();
//         std::cout << op->to_string(os);
//         std::cout << "prior mean: " << mean << "\n";
//         std::cout << "prior variance: " << variance << "\n";
//         std::cout << "expected prior mean: " << prior->get_mean() << "\n";
//         std::cout << "expected prior variance: " << prior->get_variance() << "\n";
//         
//         REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
//         REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
//         REQUIRE(mn >= prior->get_min());
//         REQUIRE(mx < prior->get_max());
//         //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
//         //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
//         
//     }
// 
//     SECTION("Testing gamma(10.0, 0.1) prior with optimizing") {
//         RandomNumberGenerator rng = RandomNumberGenerator(928374);
//         std::shared_ptr<Operator> op = std::make_shared<ChildCoalescenceRateScaler>(1.0, 0.5);
//         OperatorSchedule os = OperatorSchedule();
//         // os.turn_off_auto_optimize();
//         os.turn_on_auto_optimize();
//         os.set_auto_optimize_delay(10000);
//         os.add_operator(op);
// 
//         std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);
// 
//         std::string nex_path = "data/hemi129.nex";
//         ComparisonPopulationTree tree(nex_path, '_', true, true, false);
//         REQUIRE(tree.get_degree_of_root() == 2);
//         tree.set_root_height(0.01);
//         tree.set_root_coalescence_rate(1.0);
//         tree.set_population_size_prior(prior);
//         tree.constrain_coalescence_rates();
//         tree.ignore_data();
// 
//         tree.make_dirty();
//         tree.compute_log_likelihood_and_prior();
//         REQUIRE(! tree.is_dirty());
//     
//         unsigned int n = 0;
//         double mean = 0.0;
//         double sum_devs = 0.0;
//         double d;
//         double d_n;
//         double mn = std::numeric_limits<double>::max();
//         double mx = -std::numeric_limits<double>::max();
//         for (unsigned int i = 0; i < 100000; ++i) {
//             Operator& o = os.draw_operator(rng);
//             tree.store_state();
//             double old_v = tree.get_child_population_size(0);
//             double hastings = o.propose(rng, tree);
//             double new_v = tree.get_child_population_size(0);
//             REQUIRE(tree.is_dirty());
//             REQUIRE(new_v == tree.get_child_population_size(1));
//             REQUIRE(new_v == tree.get_root_population_size());
//             tree.compute_log_likelihood_and_prior();
//             double prior_ratio = tree.get_log_prior_density_value() -
//                 tree.get_stored_log_prior_density_value();
//             /* double prior_ratio = */
//             /*     prior->relative_ln_pdf(new_v) - */
//             /*     prior->relative_ln_pdf(old_v); */
//             double acceptance_prob = prior_ratio + hastings;
//             double u = rng.uniform_real();
//             if (u < std::exp(acceptance_prob)) {
//                 o.accept(os);
//             }
//             else {
//                 o.reject(os);
//                 tree.restore_state();
//                 REQUIRE(tree.get_child_population_size(0) == old_v);
//             }
//             o.optimize(os, acceptance_prob);
//             double x = tree.get_child_population_size(0);
//             mn = std::min(mn, x);
//             mx = std::max(mx, x);
//             ++n;
//             d = x - mean;
//             d_n = d / n;
//             mean += d_n;
//             sum_devs += d * d_n * (n - 1);
//         }
//         double variance = sum_devs / (n - 1);
//         std::cout << op->header_string();
//         std::cout << op->to_string(os);
//         std::cout << "prior mean: " << mean << "\n";
//         std::cout << "prior variance: " << variance << "\n";
//         std::cout << "expected prior mean: " << prior->get_mean() << "\n";
//         std::cout << "expected prior variance: " << prior->get_variance() << "\n";
//         
//         REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
//         REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
//         REQUIRE(mn >= prior->get_min());
//         REQUIRE(mx < prior->get_max());
//         //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
//         //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
//         
//     }
// }

TEST_CASE("Testing RootPopulationSizeScaler", "[RootPopulationSizeScaler]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(12933);
        std::shared_ptr<Operator> op = std::make_shared<RootPopulationSizeScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, '_', true, true, false);
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
            Operator& o = os.draw_operator(rng);
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
        std::shared_ptr<Operator> op = std::make_shared<RootPopulationSizeScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        //os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, '_', true, true, false);
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
            Operator& o = os.draw_operator(rng);
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
        std::shared_ptr<Operator> op = std::make_shared<ChildPopulationSizeScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, '_', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.constrain_coalescence_rates();
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
            Operator& o = os.draw_operator(rng);
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
        std::shared_ptr<Operator> op = std::make_shared<ChildPopulationSizeScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        // os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, '_', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.set_root_height(0.01);
        tree.constrain_coalescence_rates();
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
            Operator& o = os.draw_operator(rng);
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

TEST_CASE("Testing ComparisonRateMultiplierScaler", "[ComparisonRateMultiplierScaler]") {

    SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(928374);
        std::shared_ptr<Operator> op = std::make_shared<ComparisonRateMultiplierScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, '_', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_u_v_rates();
        tree.fix_coalescence_rates();
        tree.set_rate_multiplier(1.0);
        tree.set_rate_multiplier_prior(prior);
        tree.estimate_rate_multiplier();
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
            Operator& o = os.draw_operator(rng);
            tree.store_state();
            double old_v = tree.get_rate_multiplier();
            double hastings = o.propose(rng, tree);
            double new_v = tree.get_rate_multiplier();
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
            double x = tree.get_rate_multiplier();
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
        std::shared_ptr<Operator> op = std::make_shared<ComparisonRateMultiplierScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        // os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(10.0, 0.1);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, '_', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_u_v_rates();
        tree.fix_coalescence_rates();
        tree.set_rate_multiplier(1.0);
        tree.set_rate_multiplier_prior(prior);
        tree.estimate_rate_multiplier();
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
            Operator& o = os.draw_operator(rng);
            tree.store_state();
            double old_v = tree.get_rate_multiplier();
            double hastings = o.propose(rng, tree);
            double new_v = tree.get_rate_multiplier();
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
            double x = tree.get_rate_multiplier();
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

// TEST_CASE("Testing UMover", "[UMover]") {
// 
//     SECTION("Testing gamma(10.0, 0.001) prior and no optimizing") {
//         RandomNumberGenerator rng = RandomNumberGenerator(9284);
//         std::shared_ptr<Operator> op = std::make_shared<UMover>(1.0, 0.1);
//         OperatorSchedule os = OperatorSchedule();
//         os.turn_off_auto_optimize();
//         // os.turn_on_auto_optimize();
//         // os.set_auto_optimize_delay(10000);
//         os.add_operator(op);
// 
//         std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<OffsetGammaDistribution>(10.0, 0.05, 0.5);
// 
//         std::string nex_path = "data/hemi129.nex";
//         ComparisonPopulationTree tree(nex_path, '_', true, true, false);
//         REQUIRE(tree.get_degree_of_root() == 2);
//         tree.fix_coalescence_rates();
//         tree.fix_rate_multiplier();
//         tree.estimate_u_v_rates();
//         tree.set_u_prior(prior);
//         tree.set_u(1.0);
//         tree.ignore_data();
// 
//         tree.make_dirty();
//         tree.compute_log_likelihood_and_prior();
//         REQUIRE(! tree.is_dirty());
// 
//         REQUIRE(! tree.u_v_rates_are_constrained());
//         REQUIRE(! tree.u_v_rates_are_fixed());
//     
//         unsigned int n = 0;
//         double mean = 0.0;
//         double sum_devs = 0.0;
//         double d;
//         double d_n;
//         double mn = std::numeric_limits<double>::max();
//         double mx = -std::numeric_limits<double>::max();
//         for (unsigned int i = 0; i < 1000000; ++i) {
//             Operator& o = os.draw_operator(rng);
//             tree.store_state();
//             double old_v = tree.get_u();
//             double hastings = o.propose(rng, tree);
//             double new_v = tree.get_u();
//             tree.compute_log_likelihood_and_prior();
//             double prior_ratio = tree.get_log_prior_density_value() -
//                 tree.get_stored_log_prior_density_value();
//             /* double prior_ratio = */
//             /*     prior->relative_ln_pdf(new_v) - */
//             /*     prior->relative_ln_pdf(old_v); */
//             double acceptance_prob = prior_ratio + hastings;
//             double u = rng.uniform_real();
//             if (u < std::exp(acceptance_prob)) {
//                 o.accept(os);
//             }
//             else {
//                 o.reject(os);
//                 tree.restore_state();
//             }
//             o.optimize(os, acceptance_prob);
//             double x = tree.get_u();
//             mn = std::min(mn, x);
//             mx = std::max(mx, x);
//             ++n;
//             d = x - mean;
//             d_n = d / n;
//             mean += d_n;
//             sum_devs += d * d_n * (n - 1);
//         }
//         double variance = sum_devs / (n - 1);
//         std::cout << op->header_string();
//         std::cout << op->to_string(os);
//         std::cout << "prior mean: " << mean << "\n";
//         std::cout << "prior variance: " << variance << "\n";
//         std::cout << "expected prior mean: " << prior->get_mean() << "\n";
//         std::cout << "expected prior variance: " << prior->get_variance() << "\n";
//         
//         REQUIRE(mean == Approx(prior->get_mean()).epsilon(0.001));
//         REQUIRE(variance == Approx(prior->get_variance()).epsilon(0.001));
//         REQUIRE(mn >= prior->get_min());
//         REQUIRE(mx < prior->get_max());
//         //REQUIRE(mn == Approx(prior->get_min()).epsilon(0.001));
//         //REQUIRE(mx == Approx(prior->get_max()).epsilon(0.001));
//         
//     }
// }

TEST_CASE("Testing UScaler", "[UScaler]") {

    SECTION("Testing offset gamma(10.0, 0.05, 0.5) prior and no optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(9284);
        std::shared_ptr<Operator> op = std::make_shared<UScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        os.turn_off_auto_optimize();
        // os.turn_on_auto_optimize();
        // os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<OffsetGammaDistribution>(10.0, 0.05, 0.5);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, '_', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_coalescence_rates();
        tree.fix_rate_multiplier();
        tree.estimate_u_v_rates();
        tree.set_u_prior(prior);
        tree.set_u(1.0);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(! tree.u_v_rates_are_constrained());
        REQUIRE(! tree.u_v_rates_are_fixed());
    
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
            double old_v = tree.get_u();
            double hastings = o.propose(rng, tree);
            double new_v = tree.get_u();
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
            double x = tree.get_u();
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

    SECTION("Testing offset gamma(10.0, 0.05, 0.5) prior with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(9284);
        std::shared_ptr<Operator> op = std::make_shared<UScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        // os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<OffsetGammaDistribution>(10.0, 0.05, 0.5);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, '_', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_coalescence_rates();
        tree.fix_rate_multiplier();
        tree.estimate_u_v_rates();
        tree.set_u_prior(prior);
        tree.set_u(1.0);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(! tree.u_v_rates_are_constrained());
        REQUIRE(! tree.u_v_rates_are_fixed());
    
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
            double old_v = tree.get_u();
            double hastings = o.propose(rng, tree);
            double new_v = tree.get_u();
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
            double x = tree.get_u();
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

    SECTION("Testing offset exp(2.0, 0.5) prior with optimizing") {
        RandomNumberGenerator rng = RandomNumberGenerator(9284);
        std::shared_ptr<Operator> op = std::make_shared<UScaler>(1.0, 0.5);
        OperatorSchedule os = OperatorSchedule();
        // os.turn_off_auto_optimize();
        os.turn_on_auto_optimize();
        os.set_auto_optimize_delay(10000);
        os.add_operator(op);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<OffsetExponentialDistribution>(2.0, 0.5);

        std::string nex_path = "data/hemi129.nex";
        ComparisonPopulationTree tree(nex_path, '_', true, true, false);
        REQUIRE(tree.get_degree_of_root() == 2);
        tree.fix_coalescence_rates();
        tree.fix_rate_multiplier();
        tree.estimate_u_v_rates();
        tree.set_u_prior(prior);
        tree.set_u(1.0);
        tree.ignore_data();

        tree.make_dirty();
        tree.compute_log_likelihood_and_prior();
        REQUIRE(! tree.is_dirty());

        REQUIRE(! tree.u_v_rates_are_constrained());
        REQUIRE(! tree.u_v_rates_are_fixed());
    
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
            double old_v = tree.get_u();
            double hastings = o.propose(rng, tree);
            double new_v = tree.get_u();
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
            double x = tree.get_u();
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


/* TEST_CASE("Testing RootCoalescenceRateScaler", "[RootCoalescenceRateScaler]") { */
/*     SECTION("Testing gamma(10.0, 0.1) prior and no optimizing") { */
/*         std::string cfg_path = "data/dummy.yml"; */

/*         std::stringstream cfg_stream; */
/*         cfg_stream << "comparisons:\n"; */
/*         cfg_stream << "- comparison:\n"; */
/*         cfg_stream << "    path: haploid-standard-missing.nex\n"; */
/*         cfg_stream << "    genotypes_are_diploid: true\n"; */
/*         cfg_stream << "    markers_are_dominant: false\n"; */
/*         cfg_stream << "    population_name_delimiter: '_'\n"; */
/*         cfg_stream << "    population_name_is_prefix: true\n"; */
/*         cfg_stream << "    constant_sites_removed: true\n"; */
/*         cfg_stream << "    use_empirical_u_rate_starting_value: false\n"; */
/*         cfg_stream << "    constrain_population_sizes: true\n"; */
/*         cfg_stream << "    constrain_u_v_rates: true\n"; */
/*         cfg_stream << "    parameters:\n"; */
/*         cfg_stream << "        rate_multiplier:\n"; */
/*         cfg_stream << "            value: 1.0\n"; */
/*         cfg_stream << "            estimate: false\n"; */
/*         cfg_stream << "        u_rate:\n"; */
/*         cfg_stream << "            value: 1.0\n"; */
/*         cfg_stream << "            estimate: false\n"; */
/*         cfg_stream << "        population_size:\n"; */
/*         cfg_stream << "            value: 0.01\n"; */
/*         cfg_stream << "            estimate: false\n"; */
/*         cfg_stream << "            prior:\n"; */
/*         cfg_stream << "                gamma_distribution:\n"; */
/*         cfg_stream << "                    shape: 2.0\n"; */
/*         cfg_stream << "                    scale: 0.001\n"; */
/*         cfg_stream << "    path: haploid-standard.nex\n"; */

/*         CollectionSettings settings = CollectionSettings(cfg_stream, cfg_path); */
/*     } */
/* } */
