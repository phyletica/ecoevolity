#include "catch.hpp"
#include "ecoevolity/parameter.hpp"

#include "ecoevolity/probability.hpp"
#include "ecoevolity/rng.hpp"

#include <limits>

TEST_CASE("Testing RealVariable constructors", "[RealVariable]") {

    SECTION("Testing bare") {
        RealVariable p = RealVariable();
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
    }

    SECTION("Testing value") {
        RealVariable p = RealVariable(1.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_value() == Approx(1.1));
        p.store();
        REQUIRE(p.get_stored_value() == Approx(1.1));
    }
}

TEST_CASE("Testing RealParameter constructors", "[RealParameter]") {

    SECTION("Testing bare") {
        RealParameter p = RealParameter();
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE_THROWS_AS(p.fix(), EcoevolityParameterValueError &);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        REQUIRE_THROWS_AS(p.check_prior(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.draw_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.set_value_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.update_value_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_mean(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_variance(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_min(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_max(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_name(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_string(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(0.1), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(0.1), EcoevolityNullPointerError &);
    }

    SECTION("Testing value") {
        RealParameter p = RealParameter(1.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_value() == Approx(1.1));
        p.store();
        REQUIRE(p.get_stored_value() == Approx(1.1));

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        REQUIRE_THROWS_AS(p.check_prior(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.draw_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.set_value_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.update_value_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_mean(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_variance(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_min(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_max(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_name(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_string(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(0.1), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(0.1), EcoevolityNullPointerError &);
    }

    SECTION("Testing prior") {
        std::shared_ptr<UniformDistribution> u = std::make_shared<UniformDistribution>(10.0, 20.0);
        RealParameter p = RealParameter(u);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());

        REQUIRE(p.get_prior_mean() == 15.0);
        REQUIRE(p.get_prior_variance() == Approx(25.0/3.0));
        REQUIRE(p.get_prior_min() == 10.0);
        REQUIRE(p.get_prior_max() == 20.0);
        REQUIRE(p.get_prior_name() == "uniform");
        REQUIRE(p.get_prior_string() == "uniform(10, 20)");
        p.set_value(9.9);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(9.9) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(9.9) == -std::numeric_limits<double>::infinity());
        p.set_value(10.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(10.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(10.0) == Approx(std::log(0.1)));
        p.set_value(15.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(15.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(15.0) == Approx(std::log(0.1)));
        p.set_value(20.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(20.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(20.0) == Approx(std::log(0.1)));
        p.set_value(20.1);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(20.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(20.1) == -std::numeric_limits<double>::infinity());

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mx < p.get_prior_max());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
        REQUIRE(mx == Approx(p.get_prior_max()).epsilon(0.001));

        std::shared_ptr<UniformDistribution> u2 = std::make_shared<UniformDistribution>(0.0, 1.0);
        p.set_prior(u2);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());

        REQUIRE(p.get_prior_mean() == 0.5);
        REQUIRE(p.get_prior_variance() == Approx(1.0/12.0));
        REQUIRE(p.get_prior_min() == 0.0);
        REQUIRE(p.get_prior_max() == 1.0);
        REQUIRE(p.get_prior_name() == "uniform");
        REQUIRE(p.get_prior_string() == "uniform(0, 1)");
        p.set_value(-0.1);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        p.set_value(0.0);
        REQUIRE(p.prior_ln_pdf() == 0.0);
        REQUIRE(p.prior_ln_pdf(0.0) == 0.0);
        REQUIRE(p.relative_prior_ln_pdf() == 0.0);
        REQUIRE(p.relative_prior_ln_pdf(0.0) == 0.0);
        p.set_value(0.5);
        REQUIRE(p.prior_ln_pdf() == 0.0);
        REQUIRE(p.prior_ln_pdf(0.5) == 0.0);
        REQUIRE(p.relative_prior_ln_pdf() == 0.0);
        REQUIRE(p.relative_prior_ln_pdf(0.5) == 0.0);
        p.set_value(1.0);
        REQUIRE(p.prior_ln_pdf() == 0.0);
        REQUIRE(p.prior_ln_pdf(1.0) == 0.0);
        REQUIRE(p.relative_prior_ln_pdf() == 0.0);
        REQUIRE(p.relative_prior_ln_pdf(1.0) == 0.0);
        p.set_value(1.1);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(1.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(1.1) == -std::numeric_limits<double>::infinity());

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mx < p.get_prior_max());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
        REQUIRE(mx == Approx(p.get_prior_max()).epsilon(0.001));

        // Test more derived prior
        std::shared_ptr<ExponentialDistribution> f = std::make_shared<ExponentialDistribution>(5.0);
        p.set_prior(f);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());

        REQUIRE(p.get_prior_mean() == Approx(0.2));
        REQUIRE(p.get_prior_variance() == Approx(1.0/25.0));
        REQUIRE(p.get_prior_min() == 0.0);
        REQUIRE(p.get_prior_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_prior_name() == "exp");
        REQUIRE(p.get_prior_string() == "exp(lambda = 5)");
        p.set_value(-0.1);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());

        REQUIRE(p.prior_ln_pdf(0.0) == Approx(1.6094379124341003));
        REQUIRE(p.prior_ln_pdf(0.01) == Approx(1.5594379124341002));
        REQUIRE(p.prior_ln_pdf(1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.prior_ln_pdf(100.0) == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf(0.0) == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf(0.01) == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf(1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf(100.0) == Approx(-498.39056208756591));

        p.set_value(0.0);
        REQUIRE(p.prior_ln_pdf() == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.6094379124341003));

        p.set_value(0.01);
        REQUIRE(p.prior_ln_pdf() == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.5594379124341002));

        p.set_value(1.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-3.3905620875658995));

        p.set_value(100.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-498.39056208756591));

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
    }

    SECTION("Testing derived exp prior") {
        RandomNumberGenerator rng = RandomNumberGenerator(123);
        std::shared_ptr<ExponentialDistribution> f = std::make_shared<ExponentialDistribution>(5.0);
        RealParameter p = RealParameter(f);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());

        REQUIRE(p.get_prior_mean() == Approx(0.2));
        REQUIRE(p.get_prior_variance() == Approx(1.0/25.0));
        REQUIRE(p.get_prior_min() == 0.0);
        REQUIRE(p.get_prior_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_prior_name() == "exp");
        REQUIRE(p.get_prior_string() == "exp(lambda = 5)");
        p.set_value(-0.1);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());

        REQUIRE(p.prior_ln_pdf(0.0) == Approx(1.6094379124341003));
        REQUIRE(p.prior_ln_pdf(0.01) == Approx(1.5594379124341002));
        REQUIRE(p.prior_ln_pdf(1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.prior_ln_pdf(100.0) == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf(0.0) == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf(0.01) == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf(1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf(100.0) == Approx(-498.39056208756591));

        p.set_value(0.0);
        REQUIRE(p.prior_ln_pdf() == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.6094379124341003));

        p.set_value(0.01);
        REQUIRE(p.prior_ln_pdf() == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.5594379124341002));

        p.set_value(1.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-3.3905620875658995));

        p.set_value(100.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-498.39056208756591));

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
    }

    SECTION("Testing derived exp prior with base pointer") {
        RandomNumberGenerator rng = RandomNumberGenerator(123);
        std::shared_ptr<ContinuousProbabilityDistribution> f = std::make_shared<ExponentialDistribution>(5.0);
        RealParameter p = RealParameter(f);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());

        REQUIRE(p.get_prior_mean() == Approx(0.2));
        REQUIRE(p.get_prior_variance() == Approx(1.0/25.0));
        REQUIRE(p.get_prior_min() == 0.0);
        REQUIRE(p.get_prior_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_prior_name() == "exp");
        REQUIRE(p.get_prior_string() == "exp(lambda = 5)");
        p.set_value(-0.1);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());

        REQUIRE(p.prior_ln_pdf(0.0) == Approx(1.6094379124341003));
        REQUIRE(p.prior_ln_pdf(0.01) == Approx(1.5594379124341002));
        REQUIRE(p.prior_ln_pdf(1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.prior_ln_pdf(100.0) == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf(0.0) == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf(0.01) == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf(1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf(100.0) == Approx(-498.39056208756591));

        p.set_value(0.0);
        REQUIRE(p.prior_ln_pdf() == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.6094379124341003));

        p.set_value(0.01);
        REQUIRE(p.prior_ln_pdf() == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.5594379124341002));

        p.set_value(1.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-3.3905620875658995));

        p.set_value(100.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-498.39056208756591));

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
    }

    SECTION("Testing prior and value") {
        std::shared_ptr<UniformDistribution> u = std::make_shared<UniformDistribution>(10.0, 20.0);
        RealParameter p = RealParameter(u, 1.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_value() == Approx(1.1));
        p.store();
        REQUIRE(p.get_stored_value() == Approx(1.1));

        REQUIRE(p.get_prior_mean() == 15.0);
        REQUIRE(p.get_prior_variance() == Approx(25.0/3.0));
        REQUIRE(p.get_prior_min() == 10.0);
        REQUIRE(p.get_prior_max() == 20.0);
        REQUIRE(p.get_prior_name() == "uniform");
        REQUIRE(p.get_prior_string() == "uniform(10, 20)");
        p.set_value(9.9);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(9.9) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(9.9) == -std::numeric_limits<double>::infinity());
        p.set_value(10.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(10.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(10.0) == Approx(std::log(0.1)));
        p.set_value(15.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(15.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(15.0) == Approx(std::log(0.1)));
        p.set_value(20.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(20.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(20.0) == Approx(std::log(0.1)));
        p.set_value(20.1);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(20.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(20.1) == -std::numeric_limits<double>::infinity());

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mx < p.get_prior_max());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
        REQUIRE(mx == Approx(p.get_prior_max()).epsilon(0.001));

        std::shared_ptr<UniformDistribution> u2 = std::make_shared<UniformDistribution>(0.0, 1.0);
        p.set_prior(u2);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_value() == Approx(20.1));
        REQUIRE(p.get_stored_value() == Approx(1.1));
        p.store();
        REQUIRE(p.get_stored_value() == Approx(20.1));

        REQUIRE(p.get_prior_mean() == 0.5);
        REQUIRE(p.get_prior_variance() == Approx(1.0/12.0));
        REQUIRE(p.get_prior_min() == 0.0);
        REQUIRE(p.get_prior_max() == 1.0);
        REQUIRE(p.get_prior_name() == "uniform");
        REQUIRE(p.get_prior_string() == "uniform(0, 1)");
        p.set_value(-0.1);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        p.set_value(0.0);
        REQUIRE(p.prior_ln_pdf() == 0.0);
        REQUIRE(p.prior_ln_pdf(0.0) == 0.0);
        REQUIRE(p.relative_prior_ln_pdf() == 0.0);
        REQUIRE(p.relative_prior_ln_pdf(0.0) == 0.0);
        p.set_value(0.5);
        REQUIRE(p.prior_ln_pdf() == 0.0);
        REQUIRE(p.prior_ln_pdf(0.5) == 0.0);
        REQUIRE(p.relative_prior_ln_pdf() == 0.0);
        REQUIRE(p.relative_prior_ln_pdf(0.5) == 0.0);
        p.set_value(1.0);
        REQUIRE(p.prior_ln_pdf() == 0.0);
        REQUIRE(p.prior_ln_pdf(1.0) == 0.0);
        REQUIRE(p.relative_prior_ln_pdf() == 0.0);
        REQUIRE(p.relative_prior_ln_pdf(1.0) == 0.0);
        p.set_value(1.1);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(1.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(1.1) == -std::numeric_limits<double>::infinity());

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mx < p.get_prior_max());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
        REQUIRE(mx == Approx(p.get_prior_max()).epsilon(0.001));

        // Test more derived prior
        std::shared_ptr<ExponentialDistribution> f = std::make_shared<ExponentialDistribution>(5.0);
        p.set_prior(f);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());

        REQUIRE(p.get_prior_mean() == Approx(0.2));
        REQUIRE(p.get_prior_variance() == Approx(1.0/25.0));
        REQUIRE(p.get_prior_min() == 0.0);
        REQUIRE(p.get_prior_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_prior_name() == "exp");
        REQUIRE(p.get_prior_string() == "exp(lambda = 5)");
        p.set_value(-0.1);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());

        REQUIRE(p.prior_ln_pdf(0.0) == Approx(1.6094379124341003));
        REQUIRE(p.prior_ln_pdf(0.01) == Approx(1.5594379124341002));
        REQUIRE(p.prior_ln_pdf(1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.prior_ln_pdf(100.0) == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf(0.0) == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf(0.01) == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf(1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf(100.0) == Approx(-498.39056208756591));

        p.set_value(0.0);
        REQUIRE(p.prior_ln_pdf() == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.6094379124341003));

        p.set_value(0.01);
        REQUIRE(p.prior_ln_pdf() == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.5594379124341002));

        p.set_value(1.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-3.3905620875658995));

        p.set_value(100.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-498.39056208756591));

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
    }
}

TEST_CASE("Testing PositiveRealVariable constructors", "[PositiveRealVariable]") {

    SECTION("Testing bare") {
        PositiveRealVariable p = PositiveRealVariable();
        REQUIRE(std::isnan(p.get_value()));
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);
    }

    SECTION("Testing value") {
        PositiveRealVariable p = PositiveRealVariable(1.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);
        REQUIRE(p.get_value() == Approx(1.1));
        p.store();
        REQUIRE(p.get_stored_value() == Approx(1.1));
    }
}

TEST_CASE("Testing PositiveRealParameter constructors", "[PositiveRealParameter]") {

    SECTION("Testing bare") {
        PositiveRealParameter p = PositiveRealParameter();
        REQUIRE(std::isnan(p.get_value()));
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        REQUIRE_THROWS_AS(p.check_prior(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.draw_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.set_value_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.update_value_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_mean(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_variance(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_min(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_max(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_name(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_string(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(0.1), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(0.1), EcoevolityNullPointerError &);
    }

    SECTION("Testing value") {
        PositiveRealParameter p = PositiveRealParameter(1.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);
        REQUIRE(p.get_value() == Approx(1.1));
        p.store();
        REQUIRE(p.get_stored_value() == Approx(1.1));

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        REQUIRE_THROWS_AS(p.check_prior(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.draw_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.set_value_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.update_value_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_mean(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_variance(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_min(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_max(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_name(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_string(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(0.1), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(0.1), EcoevolityNullPointerError &);
    }

    SECTION("Testing value and prior") {
        std::shared_ptr<UniformDistribution> u = std::make_shared<UniformDistribution>(10.0, 20.0);
        PositiveRealParameter p = PositiveRealParameter(u, 11.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);
        REQUIRE(p.get_value() == Approx(11.1));
        p.store();
        REQUIRE(p.get_stored_value() == Approx(11.1));

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        p.initialize_value(rng);
        REQUIRE(p.get_value() == 11.1);
    }

    SECTION("Testing bad value and prior") {
        std::shared_ptr<UniformDistribution> u = std::make_shared<UniformDistribution>(10.0, 20.0);
        PositiveRealParameter p = PositiveRealParameter(u, 1.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);
        REQUIRE(p.get_value() == Approx(1.1));
        p.store();
        REQUIRE(p.get_stored_value() == Approx(1.1));

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        REQUIRE_THROWS_AS(p.initialize_value(rng), EcoevolityParameterValueError &);
    }

    SECTION("Testing prior") {
        std::shared_ptr<UniformDistribution> u = std::make_shared<UniformDistribution>(10.0, 20.0);
        PositiveRealParameter p = PositiveRealParameter(u);
        REQUIRE(std::isnan(p.get_value()));
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        p.initialize_value(rng);
        REQUIRE(! std::isnan(p.get_value()));

        REQUIRE(p.get_prior_mean() == 15.0);
        REQUIRE(p.get_prior_variance() == Approx(25.0/3.0));
        REQUIRE(p.get_prior_min() == 10.0);
        REQUIRE(p.get_prior_max() == 20.0);
        REQUIRE(p.get_prior_name() == "uniform");
        REQUIRE(p.get_prior_string() == "uniform(10, 20)");
        p.set_value(9.9);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(9.9) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(9.9) == -std::numeric_limits<double>::infinity());
        p.set_value(10.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(10.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(10.0) == Approx(std::log(0.1)));
        p.set_value(15.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(15.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(15.0) == Approx(std::log(0.1)));
        p.set_value(20.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(20.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(20.0) == Approx(std::log(0.1)));
        p.set_value(20.1);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(20.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(20.1) == -std::numeric_limits<double>::infinity());

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mx < p.get_prior_max());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
        REQUIRE(mx == Approx(p.get_prior_max()).epsilon(0.001));

        std::shared_ptr<UniformDistribution> u2 = std::make_shared<UniformDistribution>(0.0, 1.0);
        p.set_prior(u2);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);

        REQUIRE(p.get_prior_mean() == 0.5);
        REQUIRE(p.get_prior_variance() == Approx(1.0/12.0));
        REQUIRE(p.get_prior_min() == 0.0);
        REQUIRE(p.get_prior_max() == 1.0);
        REQUIRE(p.get_prior_name() == "uniform");
        REQUIRE(p.get_prior_string() == "uniform(0, 1)");
        REQUIRE_THROWS_AS(p.set_value(-0.1), EcoevolityParameterValueError &);
        REQUIRE(p.prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        p.set_value(0.0);
        REQUIRE(p.prior_ln_pdf() == 0.0);
        REQUIRE(p.prior_ln_pdf(0.0) == 0.0);
        REQUIRE(p.relative_prior_ln_pdf() == 0.0);
        REQUIRE(p.relative_prior_ln_pdf(0.0) == 0.0);
        p.set_value(0.5);
        REQUIRE(p.prior_ln_pdf() == 0.0);
        REQUIRE(p.prior_ln_pdf(0.5) == 0.0);
        REQUIRE(p.relative_prior_ln_pdf() == 0.0);
        REQUIRE(p.relative_prior_ln_pdf(0.5) == 0.0);
        p.set_value(1.0);
        REQUIRE(p.prior_ln_pdf() == 0.0);
        REQUIRE(p.prior_ln_pdf(1.0) == 0.0);
        REQUIRE(p.relative_prior_ln_pdf() == 0.0);
        REQUIRE(p.relative_prior_ln_pdf(1.0) == 0.0);
        p.set_value(1.1);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(1.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(1.1) == -std::numeric_limits<double>::infinity());

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mx < p.get_prior_max());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
        REQUIRE(mx == Approx(p.get_prior_max()).epsilon(0.001));

        // Test more derived prior
        std::shared_ptr<ExponentialDistribution> f = std::make_shared<ExponentialDistribution>(5.0);
        p.set_prior(f);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);

        REQUIRE(p.get_prior_mean() == Approx(0.2));
        REQUIRE(p.get_prior_variance() == Approx(1.0/25.0));
        REQUIRE(p.get_prior_min() == 0.0);
        REQUIRE(p.get_prior_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_prior_name() == "exp");
        REQUIRE(p.get_prior_string() == "exp(lambda = 5)");

        REQUIRE(p.prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());

        REQUIRE(p.prior_ln_pdf(0.0) == Approx(1.6094379124341003));
        REQUIRE(p.prior_ln_pdf(0.01) == Approx(1.5594379124341002));
        REQUIRE(p.prior_ln_pdf(1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.prior_ln_pdf(100.0) == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf(0.0) == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf(0.01) == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf(1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf(100.0) == Approx(-498.39056208756591));

        p.set_value(0.0);
        REQUIRE(p.prior_ln_pdf() == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.6094379124341003));

        p.set_value(0.01);
        REQUIRE(p.prior_ln_pdf() == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.5594379124341002));

        p.set_value(1.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-3.3905620875658995));

        p.set_value(100.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-498.39056208756591));

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
    }

    SECTION("Testing derived exp prior") {
        RandomNumberGenerator rng = RandomNumberGenerator(123);
        std::shared_ptr<ExponentialDistribution> f = std::make_shared<ExponentialDistribution>(5.0);
        PositiveRealParameter p = PositiveRealParameter(f);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);

        REQUIRE(p.get_prior_mean() == Approx(0.2));
        REQUIRE(p.get_prior_variance() == Approx(1.0/25.0));
        REQUIRE(p.get_prior_min() == 0.0);
        REQUIRE(p.get_prior_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_prior_name() == "exp");
        REQUIRE(p.get_prior_string() == "exp(lambda = 5)");

        REQUIRE(p.prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());

        REQUIRE(p.prior_ln_pdf(0.0) == Approx(1.6094379124341003));
        REQUIRE(p.prior_ln_pdf(0.01) == Approx(1.5594379124341002));
        REQUIRE(p.prior_ln_pdf(1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.prior_ln_pdf(100.0) == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf(0.0) == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf(0.01) == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf(1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf(100.0) == Approx(-498.39056208756591));

        p.set_value(0.0);
        REQUIRE(p.prior_ln_pdf() == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.6094379124341003));

        p.set_value(0.01);
        REQUIRE(p.prior_ln_pdf() == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.5594379124341002));

        p.set_value(1.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-3.3905620875658995));

        p.set_value(100.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-498.39056208756591));

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
    }

    SECTION("Testing prior and value") {
        std::shared_ptr<UniformDistribution> u = std::make_shared<UniformDistribution>(10.0, 20.0);
        PositiveRealParameter p = PositiveRealParameter(u, 1.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);
        REQUIRE(p.get_value() == Approx(1.1));
        p.store();
        REQUIRE(p.get_stored_value() == Approx(1.1));

        REQUIRE(p.get_prior_mean() == 15.0);
        REQUIRE(p.get_prior_variance() == Approx(25.0/3.0));
        REQUIRE(p.get_prior_min() == 10.0);
        REQUIRE(p.get_prior_max() == 20.0);
        REQUIRE(p.get_prior_name() == "uniform");
        REQUIRE(p.get_prior_string() == "uniform(10, 20)");
        p.set_value(9.9);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(9.9) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(9.9) == -std::numeric_limits<double>::infinity());
        p.set_value(10.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(10.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(10.0) == Approx(std::log(0.1)));
        p.set_value(15.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(15.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(15.0) == Approx(std::log(0.1)));
        p.set_value(20.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(20.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(20.0) == Approx(std::log(0.1)));
        p.set_value(20.1);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(20.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(20.1) == -std::numeric_limits<double>::infinity());

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mx < p.get_prior_max());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
        REQUIRE(mx == Approx(p.get_prior_max()).epsilon(0.001));

        std::shared_ptr<UniformDistribution> u2 = std::make_shared<UniformDistribution>(0.0, 1.0);
        p.set_prior(u2);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);
        REQUIRE(p.get_value() == Approx(20.1));
        REQUIRE(p.get_stored_value() == Approx(1.1));
        p.store();
        REQUIRE(p.get_stored_value() == Approx(20.1));

        REQUIRE(p.get_prior_mean() == 0.5);
        REQUIRE(p.get_prior_variance() == Approx(1.0/12.0));
        REQUIRE(p.get_prior_min() == 0.0);
        REQUIRE(p.get_prior_max() == 1.0);
        REQUIRE(p.get_prior_name() == "uniform");
        REQUIRE(p.get_prior_string() == "uniform(0, 1)");
        REQUIRE_THROWS_AS(p.set_value(-0.1), EcoevolityParameterValueError &);
        REQUIRE(p.prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        p.set_value(0.0);
        REQUIRE(p.prior_ln_pdf() == 0.0);
        REQUIRE(p.prior_ln_pdf(0.0) == 0.0);
        REQUIRE(p.relative_prior_ln_pdf() == 0.0);
        REQUIRE(p.relative_prior_ln_pdf(0.0) == 0.0);
        p.set_value(0.5);
        REQUIRE(p.prior_ln_pdf() == 0.0);
        REQUIRE(p.prior_ln_pdf(0.5) == 0.0);
        REQUIRE(p.relative_prior_ln_pdf() == 0.0);
        REQUIRE(p.relative_prior_ln_pdf(0.5) == 0.0);
        p.set_value(1.0);
        REQUIRE(p.prior_ln_pdf() == 0.0);
        REQUIRE(p.prior_ln_pdf(1.0) == 0.0);
        REQUIRE(p.relative_prior_ln_pdf() == 0.0);
        REQUIRE(p.relative_prior_ln_pdf(1.0) == 0.0);
        p.set_value(1.1);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(1.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(1.1) == -std::numeric_limits<double>::infinity());

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mx < p.get_prior_max());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
        REQUIRE(mx == Approx(p.get_prior_max()).epsilon(0.001));

        // Test more derived prior
        std::shared_ptr<ExponentialDistribution> f = std::make_shared<ExponentialDistribution>(5.0);
        p.set_prior(f);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);

        REQUIRE(p.get_prior_mean() == Approx(0.2));
        REQUIRE(p.get_prior_variance() == Approx(1.0/25.0));
        REQUIRE(p.get_prior_min() == 0.0);
        REQUIRE(p.get_prior_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_prior_name() == "exp");
        REQUIRE(p.get_prior_string() == "exp(lambda = 5)");

        REQUIRE(p.prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());

        REQUIRE(p.prior_ln_pdf(0.0) == Approx(1.6094379124341003));
        REQUIRE(p.prior_ln_pdf(0.01) == Approx(1.5594379124341002));
        REQUIRE(p.prior_ln_pdf(1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.prior_ln_pdf(100.0) == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf(0.0) == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf(0.01) == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf(1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf(100.0) == Approx(-498.39056208756591));

        p.set_value(0.0);
        REQUIRE(p.prior_ln_pdf() == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.6094379124341003));

        p.set_value(0.01);
        REQUIRE(p.prior_ln_pdf() == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.5594379124341002));

        p.set_value(1.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-3.3905620875658995));

        p.set_value(100.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-498.39056208756591));

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
    }
}

TEST_CASE("Testing IntVariable constructors", "[IntVariable]") {

    SECTION("Testing bare") {
        IntVariable p = IntVariable();
        REQUIRE(p.get_max() == std::numeric_limits<int>::max());
        REQUIRE(p.get_min() == -std::numeric_limits<int>::max());
    }

    SECTION("Testing value") {
        IntVariable p = IntVariable(2);
        REQUIRE(p.get_max() == std::numeric_limits<int>::max());
        REQUIRE(p.get_min() == -std::numeric_limits<int>::max());
        REQUIRE(p.get_value() == 2);
        p.store();
        REQUIRE(p.get_stored_value() == 2);
    }
}

TEST_CASE("Testing Probability constructors", "[Probability]") {

    SECTION("Testing bare") {
        Probability p = Probability();
        REQUIRE(p.get_max() == 1.0);
        REQUIRE(p.get_min() == 0.0);
    }

    SECTION("Testing value") {
        Probability p = Probability(0.1);
        REQUIRE(p.get_max() == 1.0);
        REQUIRE(p.get_min() == 0.0);
        REQUIRE(p.get_value() == 0.1);
        p.store();
        REQUIRE(p.get_stored_value() == 0.1);
    }
}

TEST_CASE("Testing ProbabilityDensity constructors", "[ProbabilityDensity]") {

    SECTION("Testing bare") {
        ProbabilityDensity p = ProbabilityDensity();
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);
    }

    SECTION("Testing value") {
        ProbabilityDensity p = ProbabilityDensity(0.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);
        REQUIRE(p.get_value() == 0.1);
        p.store();
        REQUIRE(p.get_stored_value() == 0.1);
    }
}

TEST_CASE("Testing LogProbability constructors", "[LogProbability]") {

    SECTION("Testing bare") {
        LogProbability p = LogProbability();
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_max() == 0.0);
    }

    SECTION("Testing value") {
        LogProbability p = LogProbability(-0.1);
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_max() == 0.0);
        REQUIRE(p.get_value() == -0.1);
        p.store();
        REQUIRE(p.get_stored_value() == -0.1);
    }
}

TEST_CASE("Testing LogProbabilityDensity constructors", "[LogProbabilityDensity]") {

    SECTION("Testing bare") {
        LogProbabilityDensity p = LogProbabilityDensity();
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
    }

    SECTION("Testing value") {
        LogProbabilityDensity p = LogProbabilityDensity(-0.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_value() == -0.1);
        p.store();
        REQUIRE(p.get_stored_value() == -0.1);
    }
}

TEST_CASE("Testing RealVariable value methods", "[RealVariable]") {

    SECTION("Testing value methods") {
        RealVariable p = RealVariable(1.1);
        REQUIRE(p.get_value() == 1.1);
        p.store();
        p.set_value(-2.2);
        REQUIRE(p.get_value() == -2.2);
        REQUIRE(p.get_stored_value() == 1.1);

        p.restore();
        REQUIRE(p.get_value() == 1.1);
        REQUIRE(p.get_stored_value() == 1.1);

        p.update_value(3.0);
        REQUIRE(p.get_value() == 3.0);
        REQUIRE(p.get_stored_value() == 1.1);

        p.update_value(4.0);
        REQUIRE(p.get_value() == 4.0);
        REQUIRE(p.get_stored_value() == 3.0);

        p.restore();
        REQUIRE(p.get_value() == 3.0);
        REQUIRE(p.get_stored_value() == 3.0);
    }
}

TEST_CASE("Testing IntVariable value methods", "[IntVariable]") {

    SECTION("Testing value methods") {
        IntVariable p = IntVariable(1);
        REQUIRE(p.get_value() == 1);
        p.store();
        p.set_value(-2);
        REQUIRE(p.get_value() == -2);
        REQUIRE(p.get_stored_value() == 1);

        p.restore();
        REQUIRE(p.get_value() == 1);
        REQUIRE(p.get_stored_value() == 1);

        p.update_value(3);
        REQUIRE(p.get_value() == 3);
        REQUIRE(p.get_stored_value() == 1);

        p.update_value(4);
        REQUIRE(p.get_value() == 4);
        REQUIRE(p.get_stored_value() == 3);

        p.restore();
        REQUIRE(p.get_value() == 3);
        REQUIRE(p.get_stored_value() == 3);
    }
}

TEST_CASE("Testing Probability value methods", "[Probability]") {

    SECTION("Testing value methods") {
        REQUIRE_THROWS_AS(Probability(1.1), EcoevolityParameterValueError &);
        REQUIRE_THROWS_AS(Probability(-0.1), EcoevolityParameterValueError &);
        Probability p = Probability(0.1);
        REQUIRE(p.get_value() == 0.1);
        p.store();
        p.set_value(0.2);
        REQUIRE(p.get_value() == 0.2);
        REQUIRE(p.get_stored_value() == 0.1);

        p.restore();
        REQUIRE(p.get_value() == 0.1);
        REQUIRE(p.get_stored_value() == 0.1);

        p.update_value(0.3);
        REQUIRE(p.get_value() == 0.3);
        REQUIRE(p.get_stored_value() == 0.1);

        p.update_value(0.4);
        REQUIRE(p.get_value() == 0.4);
        REQUIRE(p.get_stored_value() == 0.3);

        p.restore();
        REQUIRE(p.get_value() == 0.3);
        REQUIRE(p.get_stored_value() == 0.3);

        REQUIRE_THROWS_AS(p.set_value(1.01), EcoevolityParameterValueError &);
    }
}

TEST_CASE("Testing LogProbability value methods", "[LogProbability]") {

    SECTION("Testing value methods") {
        REQUIRE_THROWS_AS(LogProbability(0.1), EcoevolityParameterValueError &);
        LogProbability p = LogProbability(0.0);
        REQUIRE(p.get_value() == 0.0);
        p.store();
        p.set_value(-0.2);
        REQUIRE(p.get_value() == -0.2);
        REQUIRE(p.get_stored_value() == 0.0);

        p.restore();
        REQUIRE(p.get_value() == 0.0);
        REQUIRE(p.get_stored_value() == 0.0);

        p.update_value(-0.3);
        REQUIRE(p.get_value() == -0.3);
        REQUIRE(p.get_stored_value() == 0.0);

        p.update_value(-0.4);
        REQUIRE(p.get_value() == -0.4);
        REQUIRE(p.get_stored_value() == -0.3);

        p.restore();
        REQUIRE(p.get_value() == -0.3);
        REQUIRE(p.get_stored_value() == -0.3);

        REQUIRE_THROWS_AS(p.set_value(0.01), EcoevolityParameterValueError &);
    }
}

TEST_CASE("Testing PositiveRealVariable value methods", "[PositiveRealVariable]") {

    SECTION("Testing value methods") {
        REQUIRE_THROWS_AS(PositiveRealVariable(-0.1), EcoevolityParameterValueError &);
        PositiveRealVariable p = PositiveRealVariable(0.1);
        REQUIRE(p.get_value() == 0.1);
        p.store();
        p.set_value(0.2);
        REQUIRE(p.get_value() == 0.2);
        REQUIRE(p.get_stored_value() == 0.1);

        p.restore();
        REQUIRE(p.get_value() == 0.1);
        REQUIRE(p.get_stored_value() == 0.1);

        p.update_value(0.3);
        REQUIRE(p.get_value() == 0.3);
        REQUIRE(p.get_stored_value() == 0.1);

        p.update_value(0.4);
        REQUIRE(p.get_value() == 0.4);
        REQUIRE(p.get_stored_value() == 0.3);

        p.restore();
        REQUIRE(p.get_value() == 0.3);
        REQUIRE(p.get_stored_value() == 0.3);

        REQUIRE_THROWS_AS(p.set_value(-0.01), EcoevolityParameterValueError &);
    }
}

TEST_CASE("Testing CoalescenceRateParameter bare constructor", "[CoalescenceRateParameter]") {

    SECTION("Testing bare") {
        CoalescenceRateParameter p = CoalescenceRateParameter();
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        REQUIRE_THROWS_AS(p.check_prior(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.draw_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.set_value_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.update_value_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_mean(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_variance(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_min(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_max(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_name(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_string(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(0.1), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(0.1), EcoevolityNullPointerError &);

        REQUIRE(p.is_fixed() == false);

        REQUIRE_THROWS_AS(p.fix(), EcoevolityParameterValueError &);
        p.set_value(0.0);
        p.fix();
        REQUIRE(p.is_fixed() == true);
        REQUIRE(p.prior_ln_pdf() == 0.0);
        REQUIRE(p.prior_ln_pdf(0.1) == 0.0);
        REQUIRE(p.relative_prior_ln_pdf() == 0.0);
        REQUIRE(p.relative_prior_ln_pdf(0.1) == 0.0);

        p.estimate();
        REQUIRE(p.is_fixed() == false);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(0.1), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(0.1), EcoevolityNullPointerError &);
    }
}

TEST_CASE("Testing CoalescenceRateParameter value constructor", "[CoalescenceRateParameter]") {
    SECTION("Testing value") {
        CoalescenceRateParameter p = CoalescenceRateParameter(1.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);
        REQUIRE(p.get_value() == Approx(1.1));
        p.store();
        REQUIRE(p.get_stored_value() == Approx(1.1));

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        REQUIRE_THROWS_AS(p.check_prior(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.draw_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.set_value_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.update_value_from_prior(rng), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_mean(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_variance(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_min(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_max(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_name(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.get_prior_string(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(0.1), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(0.1), EcoevolityNullPointerError &);

        REQUIRE(p.is_fixed() == false);

        p.fix();
        REQUIRE(p.is_fixed() == true);
        REQUIRE(p.prior_ln_pdf() == 0.0);
        REQUIRE(p.prior_ln_pdf(0.1) == 0.0);
        REQUIRE(p.relative_prior_ln_pdf() == 0.0);
        REQUIRE(p.relative_prior_ln_pdf(0.1) == 0.0);

        p.estimate();
        REQUIRE(p.is_fixed() == false);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(0.1), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(), EcoevolityNullPointerError &);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(0.1), EcoevolityNullPointerError &);
    }
}

TEST_CASE("Testing CoalescenceRateParameter prior", "[CoalescenceRateParameter]") {
    SECTION("Testing prior") {
        std::shared_ptr<UniformDistribution> u = std::make_shared<UniformDistribution>(10.0, 20.0);
        CoalescenceRateParameter p = CoalescenceRateParameter(u);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);

        REQUIRE(p.get_prior_mean() == 15.0);
        REQUIRE(p.get_prior_variance() == Approx(25.0/3.0));
        REQUIRE(p.get_prior_min() == 10.0);
        REQUIRE(p.get_prior_max() == 20.0);
        REQUIRE(p.get_prior_name() == "uniform");
        REQUIRE(p.get_prior_string() == "uniform(10, 20)");
        p.set_value(2.0/9.9);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(2.0/9.9) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(2.0/9.9) == -std::numeric_limits<double>::infinity());
        p.set_value(2.0/10.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(2.0/10.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(2.0/10.0) == Approx(std::log(0.1)));
        p.set_value(2.0/15.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(2.0/15.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(2.0/15.0) == Approx(std::log(0.1)));
        p.set_value(2.0/20.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(2.0/20.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(2.0/20.0) == Approx(std::log(0.1)));
        p.set_value(2.0/20.1);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(2.0/20.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(2.0/20.1) == -std::numeric_limits<double>::infinity());

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = 2.0/p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mx < p.get_prior_max());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
        REQUIRE(mx == Approx(p.get_prior_max()).epsilon(0.001));

        std::shared_ptr<UniformDistribution> u2 = std::make_shared<UniformDistribution>(0.0, 1.0);
        p.set_prior(u2);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);

        REQUIRE(p.get_prior_mean() == 0.5);
        REQUIRE(p.get_prior_variance() == Approx(1.0/12.0));
        REQUIRE(p.get_prior_min() == 0.0);
        REQUIRE(p.get_prior_max() == 1.0);
        REQUIRE(p.get_prior_name() == "uniform");
        REQUIRE(p.get_prior_string() == "uniform(0, 1)");
        REQUIRE_THROWS_AS(p.set_value(-0.1), EcoevolityParameterValueError &);
        REQUIRE(p.prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        p.set_value(std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf() == 0.0);
        REQUIRE(p.prior_ln_pdf(std::numeric_limits<double>::infinity()) == 0.0);
        REQUIRE(p.relative_prior_ln_pdf() == 0.0);
        REQUIRE(p.relative_prior_ln_pdf(std::numeric_limits<double>::infinity()) == 0.0);
        p.set_value(2.0/0.5);
        REQUIRE(p.prior_ln_pdf() == 0.0);
        REQUIRE(p.prior_ln_pdf(2.0/0.5) == 0.0);
        REQUIRE(p.relative_prior_ln_pdf() == 0.0);
        REQUIRE(p.relative_prior_ln_pdf(2.0/0.5) == 0.0);
        p.set_value(2.0/1.0);
        REQUIRE(p.prior_ln_pdf() == 0.0);
        REQUIRE(p.prior_ln_pdf(2.0/1.0) == 0.0);
        REQUIRE(p.relative_prior_ln_pdf() == 0.0);
        REQUIRE(p.relative_prior_ln_pdf(2.0/1.0) == 0.0);
        p.set_value(2.0/1.1);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(2.0/1.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(2.0/1.1) == -std::numeric_limits<double>::infinity());

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = 2.0/p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mx < p.get_prior_max());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
        REQUIRE(mx == Approx(p.get_prior_max()).epsilon(0.001));

        // Test more derived prior
        std::shared_ptr<ExponentialDistribution> f = std::make_shared<ExponentialDistribution>(5.0);
        p.set_prior(f);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);

        REQUIRE(p.get_prior_mean() == Approx(0.2));
        REQUIRE(p.get_prior_variance() == Approx(1.0/25.0));
        REQUIRE(p.get_prior_min() == 0.0);
        REQUIRE(p.get_prior_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_prior_name() == "exp");
        REQUIRE(p.get_prior_string() == "exp(lambda = 5)");

        REQUIRE(p.prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());

        REQUIRE(p.prior_ln_pdf(std::numeric_limits<double>::infinity()) == Approx(1.6094379124341003));
        REQUIRE(p.prior_ln_pdf(2.0/0.01) == Approx(1.5594379124341002));
        REQUIRE(p.prior_ln_pdf(2.0/1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.prior_ln_pdf(2.0/100.0) == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf(std::numeric_limits<double>::infinity()) == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf(2.0/0.01) == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf(2.0/1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf(2.0/100.0) == Approx(-498.39056208756591));

        p.set_value(std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf() == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.6094379124341003));

        p.set_value(2.0/0.01);
        REQUIRE(p.prior_ln_pdf() == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.5594379124341002));

        p.set_value(2.0/1.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-3.3905620875658995));

        p.set_value(2.0/100.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-498.39056208756591));

        n = 0;
        mean = 0.0;
        sum_devs = 0.0;
        mn = std::numeric_limits<double>::max();
        mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = 2.0/p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
    }
}

TEST_CASE("Testing CoalescenceRateParameter derived prior", "[CoalescenceRateParameter]") {
    SECTION("Testing derived exp prior") {
        RandomNumberGenerator rng = RandomNumberGenerator(123);
        std::shared_ptr<ExponentialDistribution> f = std::make_shared<ExponentialDistribution>(5.0);
        CoalescenceRateParameter p = CoalescenceRateParameter(f);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);

        REQUIRE(p.get_prior_mean() == Approx(0.2));
        REQUIRE(p.get_prior_variance() == Approx(1.0/25.0));
        REQUIRE(p.get_prior_min() == 0.0);
        REQUIRE(p.get_prior_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_prior_name() == "exp");
        REQUIRE(p.get_prior_string() == "exp(lambda = 5)");

        REQUIRE(p.prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());

        REQUIRE(p.prior_ln_pdf(std::numeric_limits<double>::infinity()) == Approx(1.6094379124341003));
        REQUIRE(p.prior_ln_pdf(2.0/0.01) == Approx(1.5594379124341002));
        REQUIRE(p.prior_ln_pdf(2.0/1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.prior_ln_pdf(2.0/100.0) == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf(std::numeric_limits<double>::infinity()) == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf(2.0/0.01) == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf(2.0/1.0) == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf(2.0/100.0) == Approx(-498.39056208756591));

        p.set_value(std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf() == Approx(1.6094379124341003));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.6094379124341003));

        p.set_value(2.0/0.01);
        REQUIRE(p.prior_ln_pdf() == Approx(1.5594379124341002));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(1.5594379124341002));

        p.set_value(2.0/1.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-3.3905620875658995));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-3.3905620875658995));

        p.set_value(2.0/100.0);
        REQUIRE(p.prior_ln_pdf() == Approx(-498.39056208756591));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(-498.39056208756591));

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = 2.0/p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
    }
}


TEST_CASE("Testing CoalescenceRateParameter prior and value", "[CoalescenceRateParameter]") {
    SECTION("Testing prior and value") {
        std::shared_ptr<UniformDistribution> u = std::make_shared<UniformDistribution>(10.0, 20.0);
        CoalescenceRateParameter p = CoalescenceRateParameter(u, 1.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == 0.0);
        REQUIRE(p.get_value() == Approx(1.1));
        p.store();
        REQUIRE(p.get_stored_value() == Approx(1.1));

        REQUIRE(p.get_prior_mean() == 15.0);
        REQUIRE(p.get_prior_variance() == Approx(25.0/3.0));
        REQUIRE(p.get_prior_min() == 10.0);
        REQUIRE(p.get_prior_max() == 20.0);
        REQUIRE(p.get_prior_name() == "uniform");
        REQUIRE(p.get_prior_string() == "uniform(10, 20)");
        p.set_value(2.0/9.9);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(2.0/9.9) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(2.0/9.9) == -std::numeric_limits<double>::infinity());
        p.set_value(2.0/10.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(2.0/10.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(2.0/10.0) == Approx(std::log(0.1)));
        p.set_value(2.0/15.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(2.0/15.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(2.0/15.0) == Approx(std::log(0.1)));
        p.set_value(2.0/20.0);
        REQUIRE(p.prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.prior_ln_pdf(2.0/20.0) == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf() == Approx(std::log(0.1)));
        REQUIRE(p.relative_prior_ln_pdf(2.0/20.0) == Approx(std::log(0.1)));
        p.set_value(2.0/20.1);
        REQUIRE(p.prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.prior_ln_pdf(2.0/20.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.relative_prior_ln_pdf(2.0/20.1) == -std::numeric_limits<double>::infinity());

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = 2.0/p.draw_from_prior(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(p.get_prior_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(p.get_prior_variance()).epsilon(0.01));
        REQUIRE(mn >= p.get_prior_min());
        REQUIRE(mx < p.get_prior_max());
        REQUIRE(mn == Approx(p.get_prior_min()).epsilon(0.001));
        REQUIRE(mx == Approx(p.get_prior_max()).epsilon(0.001));
    }
}
