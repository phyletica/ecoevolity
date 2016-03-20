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
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());
    }

    SECTION("Testing value") {
        RealVariable p = RealVariable(1.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());
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
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        REQUIRE_THROWS_AS(p.check_prior(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.draw_from_prior(rng), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.set_value_from_prior(rng), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.update_value_from_prior(rng), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_mean(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_variance(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_min(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_max(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_name(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_string(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(0.1), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(0.1), EcoevolityNullPointerError);
    }

    SECTION("Testing value") {
        RealParameter p = RealParameter(1.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());
        REQUIRE(p.get_value() == Approx(1.1));
        p.store();
        REQUIRE(p.get_stored_value() == Approx(1.1));

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        REQUIRE_THROWS_AS(p.check_prior(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.draw_from_prior(rng), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.set_value_from_prior(rng), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.update_value_from_prior(rng), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_mean(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_variance(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_min(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_max(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_name(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_string(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(0.1), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(0.1), EcoevolityNullPointerError);
    }

    SECTION("Testing prior") {
        UniformDistribution * u = new UniformDistribution(10.0, 20.0);
        RealParameter p = RealParameter(u);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());

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

        UniformDistribution * u2 = new UniformDistribution(0.0, 1.0);
        p.set_prior(u2);
        delete u;
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());

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
        d;
        d_n;
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
    }

    SECTION("Testing prior and value") {
        UniformDistribution * u = new UniformDistribution(10.0, 20.0);
        RealParameter p = RealParameter(u, 1.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());
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

        UniformDistribution * u2 = new UniformDistribution(0.0, 1.0);
        p.set_prior(u2);
        delete u;
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());
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
        d;
        d_n;
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
    }
}

TEST_CASE("Testing PositiveRealVariable constructors", "[PositiveRealVariable]") {

    SECTION("Testing bare") {
        PositiveRealVariable p = PositiveRealVariable();
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == 0.0);
    }

    SECTION("Testing value") {
        PositiveRealVariable p = PositiveRealVariable(1.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == 0.0);
        REQUIRE(p.get_value() == Approx(1.1));
        p.store();
        REQUIRE(p.get_stored_value() == Approx(1.1));
    }
}

TEST_CASE("Testing PositiveRealParameter constructors", "[PositiveRealParameter]") {

    SECTION("Testing bare") {
        PositiveRealParameter p = PositiveRealParameter();
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == 0.0);

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        REQUIRE_THROWS_AS(p.check_prior(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.draw_from_prior(rng), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.set_value_from_prior(rng), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.update_value_from_prior(rng), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_mean(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_variance(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_min(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_max(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_name(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_string(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(0.1), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(0.1), EcoevolityNullPointerError);
    }

    SECTION("Testing value") {
        PositiveRealParameter p = PositiveRealParameter(1.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == 0.0);
        REQUIRE(p.get_value() == Approx(1.1));
        p.store();
        REQUIRE(p.get_stored_value() == Approx(1.1));

        RandomNumberGenerator rng = RandomNumberGenerator(123);

        REQUIRE_THROWS_AS(p.check_prior(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.draw_from_prior(rng), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.set_value_from_prior(rng), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.update_value_from_prior(rng), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_mean(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_variance(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_min(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_max(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_name(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.get_prior_string(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.prior_ln_pdf(0.1), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(), EcoevolityNullPointerError);
        REQUIRE_THROWS_AS(p.relative_prior_ln_pdf(0.1), EcoevolityNullPointerError);
    }

    SECTION("Testing prior") {
        UniformDistribution * u = new UniformDistribution(10.0, 20.0);
        PositiveRealParameter p = PositiveRealParameter(u);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == 0.0);

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

        UniformDistribution * u2 = new UniformDistribution(0.0, 1.0);
        p.set_prior(u2);
        delete u;
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == 0.0);

        REQUIRE(p.get_prior_mean() == 0.5);
        REQUIRE(p.get_prior_variance() == Approx(1.0/12.0));
        REQUIRE(p.get_prior_min() == 0.0);
        REQUIRE(p.get_prior_max() == 1.0);
        REQUIRE(p.get_prior_name() == "uniform");
        REQUIRE(p.get_prior_string() == "uniform(0, 1)");
        REQUIRE_THROWS_AS(p.set_value(-0.1), EcoevolityParameterValueError);
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
        d;
        d_n;
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
    }

    SECTION("Testing prior and value") {
        UniformDistribution * u = new UniformDistribution(10.0, 20.0);
        PositiveRealParameter p = PositiveRealParameter(u, 1.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == 0.0);
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

        UniformDistribution * u2 = new UniformDistribution(0.0, 1.0);
        p.set_prior(u2);
        delete u;
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == 0.0);
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
        REQUIRE_THROWS_AS(p.set_value(-0.1), EcoevolityParameterValueError);
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
        d;
        d_n;
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
    }
}

TEST_CASE("Testing IntVariable constructors", "[IntVariable]") {

    SECTION("Testing bare") {
        IntVariable p = IntVariable();
        REQUIRE(p.get_max() == std::numeric_limits<int>::max());
        REQUIRE(p.get_min() == -std::numeric_limits<int>::max());
        REQUIRE(p.get_upper() == std::numeric_limits<int>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<int>::lowest());
    }

    SECTION("Testing value") {
        IntVariable p = IntVariable(2);
        REQUIRE(p.get_max() == std::numeric_limits<int>::max());
        REQUIRE(p.get_min() == -std::numeric_limits<int>::max());
        REQUIRE(p.get_upper() == std::numeric_limits<int>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<int>::lowest());
        REQUIRE(p.get_value() == 2);
        p.store();
        REQUIRE(p.get_stored_value() == 2);
    }
}

TEST_CASE("Testing Probability constructors", "[Probability]") {

    SECTION("Testing bare") {
        Probability p = Probability();
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == 1.0);
        REQUIRE(p.get_lower() == 0.0);
    }

    SECTION("Testing value") {
        Probability p = Probability(0.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == 1.0);
        REQUIRE(p.get_lower() == 0.0);
        REQUIRE(p.get_value() == 0.1);
        p.store();
        REQUIRE(p.get_stored_value() == 0.1);
    }
}

TEST_CASE("Testing ProbabilityDensity constructors", "[ProbabilityDensity]") {

    SECTION("Testing bare") {
        ProbabilityDensity p = ProbabilityDensity();
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == 0.0);
    }

    SECTION("Testing value") {
        ProbabilityDensity p = ProbabilityDensity(0.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == 0.0);
        REQUIRE(p.get_value() == 0.1);
        p.store();
        REQUIRE(p.get_stored_value() == 0.1);
    }
}

TEST_CASE("Testing LogProbability constructors", "[LogProbability]") {

    SECTION("Testing bare") {
        LogProbability p = LogProbability();
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == 0.0);
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());
    }

    SECTION("Testing value") {
        LogProbability p = LogProbability(-0.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == 0.0);
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());
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
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());
    }

    SECTION("Testing value") {
        LogProbabilityDensity p = LogProbabilityDensity(-0.1);
        REQUIRE(p.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(p.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(p.get_upper() == std::numeric_limits<double>::max());
        REQUIRE(p.get_lower() == std::numeric_limits<double>::lowest());
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
        REQUIRE_THROWS_AS(Probability(1.1), EcoevolityParameterValueError);
        REQUIRE_THROWS_AS(Probability(-0.1), EcoevolityParameterValueError);
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

        REQUIRE_THROWS_AS(p.set_value(1.01), EcoevolityParameterValueError);
    }
}

TEST_CASE("Testing LogProbability value methods", "[LogProbability]") {

    SECTION("Testing value methods") {
        REQUIRE_THROWS_AS(LogProbability(0.1), EcoevolityParameterValueError);
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

        REQUIRE_THROWS_AS(p.set_value(0.01), EcoevolityParameterValueError);
    }
}

TEST_CASE("Testing PositiveRealVariable value methods", "[PositiveRealVariable]") {

    SECTION("Testing value methods") {
        REQUIRE_THROWS_AS(PositiveRealVariable(-0.1), EcoevolityParameterValueError);
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

        REQUIRE_THROWS_AS(p.set_value(-0.01), EcoevolityParameterValueError);
    }
}
