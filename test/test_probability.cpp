#include "catch.hpp"
#include "ecoevolity/probability.hpp"
#include "ecoevolity/stats_util.hpp"

#include <limits>
#include <memory>

TEST_CASE("Testing ImproperUniformDistribution", "[ImproperUniformDistribution]") {

    SECTION("Testing ImproperUniformDistribution") {
        ImproperUniformDistribution u = ImproperUniformDistribution();
        REQUIRE(u.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(u.get_min() == -std::numeric_limits<double>::infinity());
        REQUIRE(u.to_string() == "uniform(-inf, +inf)");
        REQUIRE_THROWS_AS(u.get_mean(), EcoevolityProbabilityDistributionError &);
        REQUIRE_THROWS_AS(u.get_variance(), EcoevolityProbabilityDistributionError &);
        REQUIRE_THROWS_AS(u.ln_pdf(1.0), EcoevolityProbabilityDistributionError &);
        REQUIRE(u.relative_ln_pdf(1.0) == 0.0);
        REQUIRE(u.relative_ln_pdf(100.0) == 0.0);
        REQUIRE(u.relative_ln_pdf(-1.0) == 0.0);
        REQUIRE(u.relative_ln_pdf(0.0) == 0.0);
    }
}

TEST_CASE("Testing ImproperPositiveUniformDistribution", "[ImproperPositiveUniformDistribution]") {

    SECTION("Testing ImproperPositiveUniformDistribution") {
        ImproperPositiveUniformDistribution u = ImproperPositiveUniformDistribution();
        REQUIRE(u.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(u.get_min() == 0.0);
        REQUIRE(u.to_string() == "uniform(0, +inf)");
        REQUIRE_THROWS_AS(u.get_mean(), EcoevolityProbabilityDistributionError &);
        REQUIRE_THROWS_AS(u.get_variance(), EcoevolityProbabilityDistributionError &);
        REQUIRE_THROWS_AS(u.ln_pdf(1.0), EcoevolityProbabilityDistributionError &);
        REQUIRE(u.relative_ln_pdf(1.0) == 0.0);
        REQUIRE(u.relative_ln_pdf(100.0) == 0.0);
        REQUIRE(u.relative_ln_pdf(-1.0) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.relative_ln_pdf(0.0) == 0.0);
    }
}

TEST_CASE("Testing UniformDistribution", "[UniformDistribution]") {

    SECTION("Testing UniformDistribution bare constructor") {
        UniformDistribution u = UniformDistribution();
        REQUIRE(u.get_max() == 1.0);
        REQUIRE(u.get_min() == 0.0);
        REQUIRE(u.to_string() == "uniform(0, 1)");
        REQUIRE(u.get_mean() == 0.5);
        REQUIRE(u.get_variance() == Approx(1.0/12.0));
        REQUIRE(u.relative_ln_pdf(1.0) == 0.0);
        REQUIRE(u.relative_ln_pdf(0.0) == 0.0);
        REQUIRE(u.ln_pdf(1.0) == 0.0);
        REQUIRE(u.ln_pdf(0.0) == 0.0);
        REQUIRE(u.relative_ln_pdf(-1.0) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.relative_ln_pdf(1.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.ln_pdf(-1.0) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.ln_pdf(1.1) == -std::numeric_limits<double>::infinity());

        RandomNumberGenerator rng = RandomNumberGenerator(123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = u.draw(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(u.get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(u.get_variance()).epsilon(0.001));
        REQUIRE(mn >= 0.0);
        REQUIRE(mx < 1.0);
        REQUIRE(mn == Approx(0.0).epsilon(0.001));
        REQUIRE(mx == Approx(1.0).epsilon(0.001));
    }

    SECTION("Testing UniformDistribution bare constructor") {
        UniformDistribution u = UniformDistribution(10.0, 20.0);
        REQUIRE(u.get_max() == 20.0);
        REQUIRE(u.get_min() == 10.0);
        REQUIRE(u.to_string() == "uniform(10, 20)");
        REQUIRE(u.get_mean() == 15.0);
        REQUIRE(u.get_variance() == Approx(25.0/3.0));
        REQUIRE(u.relative_ln_pdf(10.0) == Approx(std::log(0.1)));
        REQUIRE(u.relative_ln_pdf(20.0) == Approx(std::log(0.1)));
        REQUIRE(u.relative_ln_pdf(11.0) == Approx(std::log(0.1)));
        REQUIRE(u.ln_pdf(10.0) == Approx(std::log(0.1)));
        REQUIRE(u.ln_pdf(20.0) == Approx(std::log(0.1)));
        REQUIRE(u.ln_pdf(11.0) == Approx(std::log(0.1)));
        REQUIRE(u.relative_ln_pdf(9.9) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.relative_ln_pdf(20.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.relative_ln_pdf(-15.0) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.ln_pdf(9.9) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.ln_pdf(20.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(u.ln_pdf(-15.0) == -std::numeric_limits<double>::infinity());

        RandomNumberGenerator rng = RandomNumberGenerator(111);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = u.draw(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(u.get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(u.get_variance()).epsilon(0.01));
        REQUIRE(mn >= 10.0);
        REQUIRE(mx < 20.0);
        REQUIRE(mn == Approx(10.0).epsilon(0.001));
        REQUIRE(mx == Approx(20.0).epsilon(0.001));
    }
}

TEST_CASE("Testing BetaDistribution", "[BetaDistribution]") {

    SECTION("Testing constructor errors") {
        REQUIRE_THROWS_AS(BetaDistribution(0.0, 1.0), EcoevolityProbabilityDistributionError &);
        REQUIRE_THROWS_AS(BetaDistribution(1.0, 0.0), EcoevolityProbabilityDistributionError &);
        REQUIRE_THROWS_AS(BetaDistribution(-1.0, 1.0), EcoevolityProbabilityDistributionError &);
        REQUIRE_THROWS_AS(BetaDistribution(1.0, -1.0), EcoevolityProbabilityDistributionError &);
    }

    SECTION("Testing bare constructor") {
        double a = 1.0;
        double b = 1.0;
        double expected_mean = a / (a + b);
        double expected_variance = (a * b) / ((a + b) * (a + b) * (a + b + 1.0));

        BetaDistribution f = BetaDistribution();
        REQUIRE(f.to_string() == "beta(1, 1)");
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_max() == 1.0);
        REQUIRE(f.get_alpha() == a);
        REQUIRE(f.get_beta() == b);
        REQUIRE(f.get_mean() == Approx(expected_mean));
        REQUIRE(f.get_variance() == Approx(expected_variance));

        REQUIRE(f.ln_pdf(-0.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(0) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(1.0) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(1.1) == -std::numeric_limits<double>::infinity());

        REQUIRE(f.ln_pdf(0.1) == Approx(0.0));
        REQUIRE(f.ln_pdf(0.5) == Approx(0.0));
        REQUIRE(f.ln_pdf(0.9) == Approx(0.0));

        RandomNumberGenerator rng(321);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = f.draw(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn > 0.0);
        REQUIRE(mx < 1.0);
        REQUIRE(mn == Approx(0.0).epsilon(0.001));
        REQUIRE(mx == Approx(1.0).epsilon(0.001));
    }

    SECTION("Testing BetaDistribution(0.5, 0.5)") {
        double a = 0.5;
        double b = 0.5;
        double expected_mean = a / (a + b);
        double expected_variance = (a * b) / ((a + b) * (a + b) * (a + b + 1.0));

        BetaDistribution f = BetaDistribution(a, b);
        REQUIRE(f.to_string() == "beta(0.5, 0.5)");
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_max() == 1.0);
        REQUIRE(f.get_alpha() == a);
        REQUIRE(f.get_beta() == b);
        REQUIRE(f.get_mean() == Approx(expected_mean));
        REQUIRE(f.get_variance() == Approx(expected_variance));

        REQUIRE(f.ln_pdf(-0.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(0) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(1.0) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(1.1) == -std::numeric_limits<double>::infinity());

        // numbers from scipy.stats.beta.logpdf
        REQUIRE(f.ln_pdf(0.1) == Approx(0.059242918476535955));
        REQUIRE(f.ln_pdf(0.5) == Approx(-0.45158270528945466));
        REQUIRE(f.ln_pdf(0.9) == Approx(0.059242918476536177));

        RandomNumberGenerator rng(321);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = f.draw(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn > 0.0);
        REQUIRE(mx < 1.0);
        REQUIRE(mn == Approx(0.0).epsilon(0.001));
        REQUIRE(mx == Approx(1.0).epsilon(0.001));
    }

    SECTION("Testing BetaDistribution(5, 1)") {
        double a = 5.0;
        double b = 1.0;
        double expected_mean = a / (a + b);
        double expected_variance = (a * b) / ((a + b) * (a + b) * (a + b + 1.0));

        BetaDistribution f = BetaDistribution(a, b);
        REQUIRE(f.to_string() == "beta(5, 1)");
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_max() == 1.0);
        REQUIRE(f.get_alpha() == a);
        REQUIRE(f.get_beta() == b);
        REQUIRE(f.get_mean() == Approx(expected_mean));
        REQUIRE(f.get_variance() == Approx(expected_variance));

        REQUIRE(f.ln_pdf(-0.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(0) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(1.0) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(1.1) == -std::numeric_limits<double>::infinity());

        // numbers from scipy.stats.beta.logpdf
        REQUIRE(f.ln_pdf(0.1) == Approx(-7.6009024595420813));
        REQUIRE(f.ln_pdf(0.5) == Approx(-1.1631508098056809));
        REQUIRE(f.ln_pdf(0.9) == Approx(1.1879958498027952));

        RandomNumberGenerator rng(321);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = f.draw(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn > 0.0);
        REQUIRE(mx < 1.0);
        REQUIRE(mx == Approx(1.0).epsilon(0.001));
    }

    SECTION("Testing BetaDistribution(1, 5)") {
        double a = 1.0;
        double b = 5.0;
        double expected_mean = a / (a + b);
        double expected_variance = (a * b) / ((a + b) * (a + b) * (a + b + 1.0));

        BetaDistribution f = BetaDistribution(a, b);
        REQUIRE(f.to_string() == "beta(1, 5)");
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_max() == 1.0);
        REQUIRE(f.get_alpha() == a);
        REQUIRE(f.get_beta() == b);
        REQUIRE(f.get_mean() == Approx(expected_mean));
        REQUIRE(f.get_variance() == Approx(expected_variance));

        REQUIRE(f.ln_pdf(-0.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(0) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(1.0) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(1.1) == -std::numeric_limits<double>::infinity());

        // numbers from scipy.stats.beta.logpdf
        REQUIRE(f.ln_pdf(0.9) == Approx(-7.6009024595420813));
        REQUIRE(f.ln_pdf(0.5) == Approx(-1.1631508098056809));
        REQUIRE(f.ln_pdf(0.1) == Approx(1.1879958498027952));

        RandomNumberGenerator rng(321);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = f.draw(rng);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(expected_mean).epsilon(0.001));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.001));
        REQUIRE(mn > 0.0);
        REQUIRE(mx < 1.0);
        REQUIRE(mn == Approx(0.0).epsilon(0.001));
    }
}

TEST_CASE("Testing OffsetGammaDistribution", "[OffsetGammaDistribution]") {
    SECTION("Testing bare constructor") {
        OffsetGammaDistribution f = OffsetGammaDistribution();
        REQUIRE(f.to_string() == "gamma(shape = 1, scale = 1, offset = 0)");
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_mean() == 1.0);
        REQUIRE(f.get_variance() == 1.0);

        REQUIRE(f.get_offset() == 0.0);
        REQUIRE(f.get_shape() == 1.0);
        REQUIRE(f.get_scale() == 1.0);

        REQUIRE(f.ln_pdf(-0.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(0.0) == 0.0);
        REQUIRE(f.relative_ln_pdf(0.0) == 0.0);
        REQUIRE(f.ln_pdf(0.01) == -0.01);
        REQUIRE(f.relative_ln_pdf(0.01) == -0.01);
        REQUIRE(f.ln_pdf(1.0) == -1.0);
        REQUIRE(f.relative_ln_pdf(1.0) == -1.0);
        REQUIRE(f.ln_pdf(100.0) == -100.0);
        REQUIRE(f.relative_ln_pdf(100.0) == -100.0);

        RandomNumberGenerator rng = RandomNumberGenerator(123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = f.draw(rng);
            mn = std::min(mn, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(f.get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(f.get_variance()).epsilon(0.01));
        REQUIRE(mn >= 0.0);
        REQUIRE(mn == Approx(0.0).epsilon(0.001));
    }

    SECTION("Testing constructor errors") {
        REQUIRE_THROWS_AS(OffsetGammaDistribution(0.0, 1.0, 1.0), EcoevolityProbabilityDistributionError &);
        REQUIRE_THROWS_AS(OffsetGammaDistribution(1.0, 0.0, 1.0), EcoevolityProbabilityDistributionError &);
        REQUIRE_THROWS_AS(OffsetGammaDistribution(-1.0, 1.0, 1.0), EcoevolityProbabilityDistributionError &);
        REQUIRE_THROWS_AS(OffsetGammaDistribution(1.0, -1.0, 1.0), EcoevolityProbabilityDistributionError &);
    }

    SECTION("Testing OffsetGammaDistribution(2, 4, 5)") {
        OffsetGammaDistribution f = OffsetGammaDistribution(2.0, 4.0, 5.0);
        REQUIRE(f.to_string() == "gamma(shape = 2, scale = 4, offset = 5)");
        REQUIRE(f.get_min() == 5.0);
        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_mean() == 13.0);
        REQUIRE(f.get_variance() == 32.0);
        REQUIRE(f.ln_pdf(5.0) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(4.9) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(6.0) == Approx(-3.0225887222397811));
        REQUIRE(f.ln_pdf(13.0) == Approx(-2.6931471805599454));
        REQUIRE(f.ln_pdf(105.0) == Approx(-23.167418536251692));

        REQUIRE(f.get_offset() == 5.0);
        REQUIRE(f.get_shape() == 2.0);
        REQUIRE(f.get_scale() == 4.0);

        RandomNumberGenerator rng = RandomNumberGenerator(123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = f.draw(rng);
            mn = std::min(mn, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(f.get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(f.get_variance()).epsilon(0.01));
        REQUIRE(mn >= 5.0);
        REQUIRE(mn == Approx(5.0).epsilon(0.01));
    }

    SECTION("Testing OffsetGammaDistribution(0.1, 100, -5)") {
        OffsetGammaDistribution f = OffsetGammaDistribution(0.1, 100.0, -5.0);
        REQUIRE(f.to_string() == "gamma(shape = 0.1, scale = 100, offset = -5)");
        REQUIRE(f.get_min() == -5.0);
        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_mean() == 5.0);
        REQUIRE(f.get_variance() == 1000.0);
        REQUIRE(f.ln_pdf(-5.0) == std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(-5.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(-4.0) == Approx(-2.7232296703330157));
        REQUIRE(f.ln_pdf(0.0) == Approx(-4.2117237915237062));
        REQUIRE(f.ln_pdf(5.0) == Approx(-4.8855562540276569));
        REQUIRE(f.ln_pdf(45.0) == Approx(-6.7340503752183469));

        REQUIRE(f.get_offset() == -5.0);
        REQUIRE(f.get_shape() == 0.1);
        REQUIRE(f.get_scale() == 100.0);

        RandomNumberGenerator rng = RandomNumberGenerator(123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = f.draw(rng);
            mn = std::min(mn, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(f.get_mean()).epsilon(0.01));
        REQUIRE(variance == Approx(f.get_variance()).epsilon(1.0));
        REQUIRE(mn >= -5.0);
        REQUIRE(mn == Approx(-5.0).epsilon(0.001));
    }
}

TEST_CASE("Testing GammaDistribution", "[GammaDistribution]") {
    SECTION("Testing bare constructor") {
        GammaDistribution f = GammaDistribution();
        REQUIRE(f.to_string() == "gamma(shape = 1, scale = 1)");
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_mean() == 1.0);
        REQUIRE(f.get_variance() == 1.0);

        REQUIRE(f.ln_pdf(-0.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(0.0) == 0.0);
        REQUIRE(f.relative_ln_pdf(0.0) == 0.0);
        REQUIRE(f.ln_pdf(0.01) == -0.01);
        REQUIRE(f.relative_ln_pdf(0.01) == -0.01);
        REQUIRE(f.ln_pdf(1.0) == -1.0);
        REQUIRE(f.relative_ln_pdf(1.0) == -1.0);
        REQUIRE(f.ln_pdf(100.0) == -100.0);
        REQUIRE(f.relative_ln_pdf(100.0) == -100.0);

        REQUIRE(f.get_offset() == 0.0);
        REQUIRE(f.get_shape() == 1.0);
        REQUIRE(f.get_scale() == 1.0);

        RandomNumberGenerator rng = RandomNumberGenerator(1123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = f.draw(rng);
            mn = std::min(mn, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(f.get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(f.get_variance()).epsilon(0.01));
        REQUIRE(mn >= 0.0);
        REQUIRE(mn == Approx(0.0).epsilon(0.001));
    }
    SECTION("Testing constructor errors") {
        REQUIRE_THROWS_AS(GammaDistribution(0.0, 1.0), EcoevolityProbabilityDistributionError &);
        REQUIRE_THROWS_AS(GammaDistribution(1.0, 0.0), EcoevolityProbabilityDistributionError &);
        REQUIRE_THROWS_AS(GammaDistribution(-1.0, 1.0), EcoevolityProbabilityDistributionError &);
        REQUIRE_THROWS_AS(GammaDistribution(1.0, -1.0), EcoevolityProbabilityDistributionError &);
    }

    SECTION("Testing GammaDistribution(2, 4)") {
        GammaDistribution f = GammaDistribution(2.0, 4.0);
        REQUIRE(f.to_string() == "gamma(shape = 2, scale = 4)");
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_mean() == 8.0);
        REQUIRE(f.get_variance() == 32.0);
        REQUIRE(f.ln_pdf(0.0) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(1.0) == Approx(-3.0225887222397811));
        REQUIRE(f.ln_pdf(8.0) == Approx(-2.6931471805599454));
        REQUIRE(f.ln_pdf(100.0) == Approx(-23.167418536251692));

        REQUIRE(f.get_offset() == 0.0);
        REQUIRE(f.get_shape() == 2.0);
        REQUIRE(f.get_scale() == 4.0);

        RandomNumberGenerator rng = RandomNumberGenerator(1223);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = f.draw(rng);
            mn = std::min(mn, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(f.get_mean()).epsilon(0.01));
        REQUIRE(variance == Approx(f.get_variance()).epsilon(0.01));
        REQUIRE(mn >= 0.0);
        REQUIRE(mn == Approx(0.0).epsilon(0.01));
    }

    SECTION("Testing GammaDistribution(0.1, 100)") {
        GammaDistribution f = GammaDistribution(0.1, 100.0);
        REQUIRE(f.to_string() == "gamma(shape = 0.1, scale = 100)");
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_mean() == 10.0);
        REQUIRE(f.get_variance() == 1000.0);
        REQUIRE(f.ln_pdf(0.0) == std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(-0.1) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(1.0) == Approx(-2.7232296703330157));
        REQUIRE(f.ln_pdf(5.0) == Approx(-4.2117237915237062));
        REQUIRE(f.ln_pdf(10.0) == Approx(-4.8855562540276569));
        REQUIRE(f.ln_pdf(50.0) == Approx(-6.7340503752183469));

        REQUIRE(f.get_offset() == 0.0);
        REQUIRE(f.get_shape() == 0.1);
        REQUIRE(f.get_scale() == 100.0);

        RandomNumberGenerator rng = RandomNumberGenerator(123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = f.draw(rng);
            mn = std::min(mn, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(f.get_mean()).epsilon(0.01));
        REQUIRE(variance == Approx(f.get_variance()).epsilon(1.0));
        REQUIRE(mn >= 0.0);
        REQUIRE(mn == Approx(0.0).epsilon(0.001));
    }
}

TEST_CASE("Testing OffsetExponentialDistribution", "[OffsetExponentialDistribution]") {
    SECTION("Testing bare constructor") {
        OffsetExponentialDistribution f = OffsetExponentialDistribution();

        REQUIRE(f.to_string() == "exp(lambda = 1, offset = 0)");

        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_offset() == 0.0);
        REQUIRE(f.get_shape() == 1.0);
        REQUIRE(f.get_scale() == 1.0);
        REQUIRE(f.get_lambda() == 1.0);

        REQUIRE(f.get_mean() == 1.0);
        REQUIRE(f.get_variance() == 1.0);

        REQUIRE(f.ln_pdf(-0.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(0.0) == 0.0);
        REQUIRE(f.ln_pdf(0.01) == -0.01);
        REQUIRE(f.ln_pdf(1.0) == -1.0);
        REQUIRE(f.ln_pdf(100.0) == -100.0);

        RandomNumberGenerator rng = RandomNumberGenerator(123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = f.draw(rng);
            mn = std::min(mn, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(f.get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(f.get_variance()).epsilon(0.01));
        REQUIRE(mn >= 0.0);
        REQUIRE(mn == Approx(0.0).epsilon(0.001));
    }

    SECTION("Testing constructor errors") {
        REQUIRE_THROWS_AS(OffsetExponentialDistribution(0.0, 1.0), EcoevolityProbabilityDistributionError &);
        REQUIRE_THROWS_AS(OffsetExponentialDistribution(-0.01, 1.0), EcoevolityProbabilityDistributionError &);
    }

    SECTION("Testing OffsetExponentialDistribution(5, -5)") {
        OffsetExponentialDistribution f = OffsetExponentialDistribution(5.0, -5.0);

        REQUIRE(f.to_string() == "exp(lambda = 5, offset = -5)");

        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_min() == -5.0);
        REQUIRE(f.get_offset() == -5.0);
        REQUIRE(f.get_shape() == 1.0);
        REQUIRE(f.get_scale() == Approx(1.0/5.0));
        REQUIRE(f.get_lambda() == 5.0);

        REQUIRE(f.get_mean() == Approx(-4.8));
        REQUIRE(f.get_variance() == Approx(1.0/25.0));

        REQUIRE(f.ln_pdf(-5.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(-5.0) == Approx(1.6094379124341003));
        REQUIRE(f.ln_pdf(-4.99) == Approx(1.5594379124341002));
        REQUIRE(f.ln_pdf(-4.0) == Approx(-3.3905620875658995));
        REQUIRE(f.ln_pdf(95.0) == Approx(-498.39056208756591));

        RandomNumberGenerator rng = RandomNumberGenerator(123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = f.draw(rng);
            mn = std::min(mn, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(f.get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(f.get_variance()).epsilon(0.01));
        REQUIRE(mn >= -5.0);
        REQUIRE(mn == Approx(-5.0).epsilon(0.001));
    }
}

TEST_CASE("Testing ExponentialDistribution", "[ExponentialDistribution]") {
    SECTION("Testing bare constructor") {
        ExponentialDistribution f = ExponentialDistribution();

        REQUIRE(f.to_string() == "exp(lambda = 1)");

        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_offset() == 0.0);
        REQUIRE(f.get_shape() == 1.0);
        REQUIRE(f.get_scale() == 1.0);
        REQUIRE(f.get_lambda() == 1.0);

        REQUIRE(f.get_mean() == 1.0);
        REQUIRE(f.get_variance() == 1.0);

        REQUIRE(f.ln_pdf(-0.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(0.0) == 0.0);
        REQUIRE(f.ln_pdf(0.01) == -0.01);
        REQUIRE(f.ln_pdf(1.0) == -1.0);
        REQUIRE(f.ln_pdf(100.0) == -100.0);

        RandomNumberGenerator rng = RandomNumberGenerator(123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = f.draw(rng);
            mn = std::min(mn, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(f.get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(f.get_variance()).epsilon(0.01));
        REQUIRE(mn >= 0.0);
        REQUIRE(mn == Approx(0.0).epsilon(0.001));
    }

    SECTION("Testing constructor errors") {
        REQUIRE_THROWS_AS(ExponentialDistribution(0.0), EcoevolityProbabilityDistributionError &);
        REQUIRE_THROWS_AS(ExponentialDistribution(-0.01), EcoevolityProbabilityDistributionError &);
    }

    SECTION("Testing ExponentialDistribution(5)") {
        ExponentialDistribution f = ExponentialDistribution(5.0);

        REQUIRE(f.to_string() == "exp(lambda = 5)");

        REQUIRE(f.get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_offset() == 0.0);
        REQUIRE(f.get_shape() == 1.0);
        REQUIRE(f.get_scale() == Approx(1.0/5.0));
        REQUIRE(f.get_lambda() == 5.0);

        REQUIRE(f.get_mean() == Approx(0.2));
        REQUIRE(f.get_variance() == Approx(1.0/25.0));

        REQUIRE(f.ln_pdf(-0.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(0.0) == Approx(1.6094379124341003));
        REQUIRE(f.ln_pdf(0.01) == Approx(1.5594379124341002));
        REQUIRE(f.ln_pdf(1.0) == Approx(-3.3905620875658995));
        REQUIRE(f.ln_pdf(100.0) == Approx(-498.39056208756591));

        RandomNumberGenerator rng = RandomNumberGenerator(123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = f.draw(rng);
            mn = std::min(mn, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(f.get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(f.get_variance()).epsilon(0.01));
        REQUIRE(mn >= 0.0);
        REQUIRE(mn == Approx(0.0).epsilon(0.001));
    }

    SECTION("Testing ExponentialDistribution(5) with base pointer") {
        std::shared_ptr<ContinuousProbabilityDistribution> f = std::make_shared<ExponentialDistribution>(5.0);

        REQUIRE(f->to_string() == "exp(lambda = 5)");

        REQUIRE(f->get_max() == std::numeric_limits<double>::infinity());
        REQUIRE(f->get_min() == 0.0);

        REQUIRE(f->get_mean() == Approx(0.2));
        REQUIRE(f->get_variance() == Approx(1.0/25.0));

        REQUIRE(f->ln_pdf(-0.01) == -std::numeric_limits<double>::infinity());
        REQUIRE(f->ln_pdf(0.0) == Approx(1.6094379124341003));
        REQUIRE(f->ln_pdf(0.01) == Approx(1.5594379124341002));
        REQUIRE(f->ln_pdf(1.0) == Approx(-3.3905620875658995));
        REQUIRE(f->ln_pdf(100.0) == Approx(-498.39056208756591));

        RandomNumberGenerator rng = RandomNumberGenerator(123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = f->draw(rng);
            mn = std::min(mn, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(f->get_mean()).epsilon(0.001));
        REQUIRE(variance == Approx(f->get_variance()).epsilon(0.01));
        REQUIRE(mn >= 0.0);
        REQUIRE(mn == Approx(0.0).epsilon(0.001));
    }
}

TEST_CASE("Testing DirichletDistribution constructor errors",
        "[DirichletDistribution]") {

    SECTION("Testing constructor errors") {
        std::vector<double> parameters = {0.0, 2.0};
        REQUIRE_THROWS_AS(DirichletDistribution f = DirichletDistribution(parameters), EcoevolityProbabilityDistributionError &);
        parameters = {2.0, 0.0};
        REQUIRE_THROWS_AS(DirichletDistribution f = DirichletDistribution(parameters), EcoevolityProbabilityDistributionError &);
        parameters = {1.0, -1.0};
        REQUIRE_THROWS_AS(DirichletDistribution f = DirichletDistribution(parameters), EcoevolityProbabilityDistributionError &);
        parameters = {-1.0, 1.0};
        REQUIRE_THROWS_AS(DirichletDistribution f = DirichletDistribution(parameters), EcoevolityProbabilityDistributionError &);
    }
}

TEST_CASE("Testing DirichletDistribution(1, 1, 1) bare constructor",
        "[DirichletDistribution]") {

    SECTION("Testing bare constructor") {
        std::vector<double> parameters {1.0, 1.0, 1.0};
        double sum = 0.0;
        for (auto p: parameters) {
            sum += p;
        }
        double k = parameters.size();
        std::vector<double> expected_means(k, 0.0);
        std::vector<double> expected_variances(k, 0.0);
        for (unsigned int i = 0; i < k; ++i) {
            double p_i = parameters.at(i);
            expected_means.at(i) = p_i / sum;
            expected_variances.at(i) = (p_i * (sum - p_i)) / (sum * sum * (sum + 1));
        }

        DirichletDistribution f = DirichletDistribution();
        REQUIRE(f.to_string() == "dirichlet(1, 1, 1)");
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_max() == 1.0);
        REQUIRE(f.get_parameters() == parameters);
        std::vector<double> means = f.get_mean();
        std::vector<double> variances = f.get_variance();
        for (unsigned int i = 0; i < k; ++i) {
            REQUIRE(means.at(i) == Approx(expected_means.at(i)));
            REQUIRE(variances.at(i) == Approx(expected_variances.at(i)));
        }

        std::vector<double> values = {0.0, 0.5, 0.5};
        REQUIRE(f.ln_pdf(values) == -std::numeric_limits<double>::infinity());
        values = {1.0, 0.0, 0.0};
        REQUIRE(f.ln_pdf(values) == -std::numeric_limits<double>::infinity());

        // numbers from scipy.stats.dirichlet.logpdf
        values = {0.1, 0.3, 0.6};
        double expected1 = 0.69314718055994529;
        double ln_p1 = f.ln_pdf(values);
        double rln_p1 = f.relative_ln_pdf(values);
        REQUIRE(ln_p1 == Approx(expected1));
        REQUIRE(rln_p1 == Approx(0.0));

        values = {0.3, 0.4, 0.3};
        double expected2 = 0.69314718055994529;
        double ln_p2 = f.ln_pdf(values);
        double rln_p2 = f.relative_ln_pdf(values);
        REQUIRE(ln_p2 == Approx(expected2));
        REQUIRE(rln_p2 == Approx(0.0));


        std::vector< SampleSummarizer<double> > summaries(k);
        std::vector<double> x(k);
        unsigned int nsamples = 100000;
        RandomNumberGenerator rng(123);
        for (unsigned int i = 0; i < nsamples; ++i) {
            x = f.draw(rng);
            double s = 0.0;
            for (unsigned int j = 0; j < k; ++j) {
                summaries.at(j).add_sample(x.at(j));
                s += x.at(j);
            }
            if (i < 100) {
                REQUIRE(s == Approx(1.0));
            }
        }
        
        for (unsigned int i = 0; i < k; ++i) {
            REQUIRE(summaries.at(i).sample_size() == nsamples);
            REQUIRE(summaries.at(i).mean() == Approx(expected_means.at(i)).epsilon(0.001));
            REQUIRE(summaries.at(i).variance() == Approx(expected_variances.at(i)).epsilon(0.001));
            REQUIRE(summaries.at(i).min() > 0.0);
            REQUIRE(summaries.at(i).max() <= 1.0);
        }
    }
}

TEST_CASE("Testing DirichletDistribution(0.5, 3, 1, 5)",
        "[DirichletDistribution]") {

    SECTION("Testing p = 0.5, 3, 1, 5") {
        std::vector<double> parameters {0.5, 3.0, 1.0, 5.0};
        double sum = 0.0;
        for (auto p: parameters) {
            sum += p;
        }
        double k = parameters.size();
        std::vector<double> expected_means(k, 0.0);
        std::vector<double> expected_variances(k, 0.0);
        for (unsigned int i = 0; i < k; ++i) {
            double p_i = parameters.at(i);
            expected_means.at(i) = p_i / sum;
            expected_variances.at(i) = (p_i * (sum - p_i)) / (sum * sum * (sum + 1));
        }

        DirichletDistribution f = DirichletDistribution(parameters);
        REQUIRE(f.to_string() == "dirichlet(0.5, 3, 1, 5)");
        REQUIRE(f.get_min() == 0.0);
        REQUIRE(f.get_max() == 1.0);
        REQUIRE(f.get_parameters() == parameters);
        std::vector<double> means = f.get_mean();
        std::vector<double> variances = f.get_variance();
        for (unsigned int i = 0; i < k; ++i) {
            REQUIRE(means.at(i) == Approx(expected_means.at(i)));
            REQUIRE(variances.at(i) == Approx(expected_variances.at(i)));
        }

        std::vector<double> values = {0.0, 0.3, 0.3, 0.4};
        REQUIRE(f.ln_pdf(values) == -std::numeric_limits<double>::infinity());

        // numbers from scipy.stats.dirichlet.logpdf
        values = {0.25, 0.25, 0.25, 0.25};
        double expected1 = -0.37885151919472015;
        double ln_p1 = f.ln_pdf(values);
        double rln_p1 = f.relative_ln_pdf(values);
        REQUIRE(ln_p1 == Approx(expected1));

        values = {0.4, 0.1, 0.3, 0.2};
        double expected2 = -3.3390090028227366;
        double ln_p2 = f.ln_pdf(values);
        double rln_p2 = f.relative_ln_pdf(values);
        REQUIRE(ln_p2 == Approx(expected2));

        REQUIRE((rln_p1 - rln_p2) == Approx(ln_p1 - ln_p2));

        std::vector< SampleSummarizer<double> > summaries(k);
        std::vector<double> x(k);
        unsigned int nsamples = 100000;
        RandomNumberGenerator rng(123);
        for (unsigned int i = 0; i < nsamples; ++i) {
            x = f.draw(rng);
            double s = 0.0;
            for (unsigned int j = 0; j < k; ++j) {
                summaries.at(j).add_sample(x.at(j));
                s += x.at(j);
            }
            if (i < 100) {
                REQUIRE(s == Approx(1.0));
            }
        }
        
        for (unsigned int i = 0; i < k; ++i) {
            REQUIRE(summaries.at(i).sample_size() == nsamples);
            REQUIRE(summaries.at(i).mean() == Approx(expected_means.at(i)).epsilon(0.001));
            REQUIRE(summaries.at(i).variance() == Approx(expected_variances.at(i)).epsilon(0.001));
            REQUIRE(summaries.at(i).min() > 0.0);
            REQUIRE(summaries.at(i).max() <= 1.0);
        }
    }
}

TEST_CASE("Testing relative_ln_pdf of DirichletDistribution(1, 1, 1)",
        "[DirichletDistribution]") {

    SECTION("Testing relative_ln_pdf") {
        std::vector<double> alphas {1.0, 1.0, 1.0};
        DirichletDistribution f = DirichletDistribution(alphas);

        std::vector<double> v1 = {1.0/3.0, 1.0/3.0, 1.0/3.0};
        std::vector<double> v2 = {0.1, 0.2, 0.7};
        std::vector<double> v3 = {0.9, 0.05, 0.05};
        std::vector<double> v4 = {0.2, 0.6, 0.2};

        REQUIRE(f.ln_pdf(v1) != -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(v2) != -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(v3) != -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(v4) != -std::numeric_limits<double>::infinity());

        REQUIRE((f.ln_pdf(v1) - f.ln_pdf(v2)) == Approx(f.relative_ln_pdf(v1) - f.relative_ln_pdf(v2)));
        REQUIRE((f.ln_pdf(v1) - f.ln_pdf(v3)) == Approx(f.relative_ln_pdf(v1) - f.relative_ln_pdf(v3)));
        REQUIRE((f.ln_pdf(v1) - f.ln_pdf(v4)) == Approx(f.relative_ln_pdf(v1) - f.relative_ln_pdf(v4)));
        REQUIRE((f.ln_pdf(v2) - f.ln_pdf(v3)) == Approx(f.relative_ln_pdf(v2) - f.relative_ln_pdf(v3)));
        REQUIRE((f.ln_pdf(v2) - f.ln_pdf(v4)) == Approx(f.relative_ln_pdf(v2) - f.relative_ln_pdf(v4)));
        REQUIRE((f.ln_pdf(v3) - f.ln_pdf(v4)) == Approx(f.relative_ln_pdf(v3) - f.relative_ln_pdf(v4)));
    }
}

TEST_CASE("Testing relative_ln_pdf of DirichletDistribution(3, 1, 2)",
        "[DirichletDistribution]") {

    SECTION("Testing relative_ln_pdf") {
        std::vector<double> alphas {3.0, 1.0, 2.0};
        DirichletDistribution f = DirichletDistribution(alphas);

        std::vector<double> v1 = {1.0/3.0, 1.0/3.0, 1.0/3.0};
        std::vector<double> v2 = {0.1, 0.2, 0.7};
        std::vector<double> v3 = {0.9, 0.05, 0.05};
        std::vector<double> v4 = {0.2, 0.6, 0.2};

        REQUIRE(f.ln_pdf(v1) != -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(v2) != -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(v3) != -std::numeric_limits<double>::infinity());
        REQUIRE(f.ln_pdf(v4) != -std::numeric_limits<double>::infinity());

        REQUIRE((f.ln_pdf(v1) - f.ln_pdf(v2)) == Approx(f.relative_ln_pdf(v1) - f.relative_ln_pdf(v2)));
        REQUIRE((f.ln_pdf(v1) - f.ln_pdf(v3)) == Approx(f.relative_ln_pdf(v1) - f.relative_ln_pdf(v3)));
        REQUIRE((f.ln_pdf(v1) - f.ln_pdf(v4)) == Approx(f.relative_ln_pdf(v1) - f.relative_ln_pdf(v4)));
        REQUIRE((f.ln_pdf(v2) - f.ln_pdf(v3)) == Approx(f.relative_ln_pdf(v2) - f.relative_ln_pdf(v3)));
        REQUIRE((f.ln_pdf(v2) - f.ln_pdf(v4)) == Approx(f.relative_ln_pdf(v2) - f.relative_ln_pdf(v4)));
        REQUIRE((f.ln_pdf(v3) - f.ln_pdf(v4)) == Approx(f.relative_ln_pdf(v3) - f.relative_ln_pdf(v4)));
    }
}
