#include "catch.hpp"
#include "ecoevolity/stats_util.hpp"

TEST_CASE("Testing double SampleSummarizer", "[stats_util]") {

    SECTION("Testing empty double summarizer") {
        SampleSummarizer<double> ss;
        REQUIRE(ss.sample_size() == 0);
        REQUIRE_THROWS_AS(ss.min(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.max(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.mean(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.variance(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.population_variance(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.std_dev(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.std_error(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.skewness(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.kurtosis(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.excess_kurtosis(), EcoevolityError);
    }

    SECTION("Testing double summarizer with one sample 0") {
        SampleSummarizer<double> ss;
        ss.add_sample(0.0);
        REQUIRE(ss.sample_size() == 1);
        REQUIRE(ss.min() == 0.0);
        REQUIRE(ss.max() == 0.0);
        REQUIRE(ss.mean() == 0.0);
        REQUIRE(ss.variance() == std::numeric_limits<double>::infinity());
        REQUIRE(ss.population_variance() == 0.0);
        REQUIRE(ss.std_dev() == std::numeric_limits<double>::infinity());
        REQUIRE(ss.std_error() == std::numeric_limits<double>::infinity());
        REQUIRE(std::isnan(ss.skewness()));
        REQUIRE(std::isnan(ss.kurtosis()));
        REQUIRE(std::isnan(ss.excess_kurtosis()));
    }

    SECTION("Testing double summarizer with one sample -10") {
        SampleSummarizer<double> ss;
        ss.add_sample(-10.0);
        REQUIRE(ss.sample_size() == 1);
        REQUIRE(ss.min() == -10.0);
        REQUIRE(ss.max() == -10.0);
        REQUIRE(ss.mean() == -10.0);
        REQUIRE(ss.variance() == std::numeric_limits<double>::infinity());
        REQUIRE(ss.population_variance() == 0.0);
        REQUIRE(ss.std_dev() == std::numeric_limits<double>::infinity());
        REQUIRE(ss.std_error() == std::numeric_limits<double>::infinity());
        REQUIRE(std::isnan(ss.skewness()));
        REQUIRE(std::isnan(ss.kurtosis()));
        REQUIRE(std::isnan(ss.excess_kurtosis()));
    }

    SECTION("Testing double summarizer with one sample 10") {
        SampleSummarizer<double> ss;
        ss.add_sample(10.0);
        REQUIRE(ss.sample_size() == 1);
        REQUIRE(ss.min() == 10.0);
        REQUIRE(ss.max() == 10.0);
        REQUIRE(ss.mean() == 10.0);
        REQUIRE(ss.variance() == std::numeric_limits<double>::infinity());
        REQUIRE(ss.population_variance() == 0.0);
        REQUIRE(ss.std_dev() == std::numeric_limits<double>::infinity());
        REQUIRE(ss.std_error() == std::numeric_limits<double>::infinity());
        REQUIRE(std::isnan(ss.skewness()));
        REQUIRE(std::isnan(ss.kurtosis()));
        REQUIRE(std::isnan(ss.excess_kurtosis()));
    }

    SECTION("Testing double summarizer with 11 samples around 0") {
        SampleSummarizer<double> ss;
        ss.add_sample(0.0);
        ss.add_sample(1.0);
        ss.add_sample(-1.0);
        ss.add_sample(2.0);
        ss.add_sample(-2.0);
        ss.add_sample(3.0);
        ss.add_sample(-3.0);
        ss.add_sample(4.0);
        ss.add_sample(-4.0);
        ss.add_sample(5.0);
        ss.add_sample(-5.0);
        REQUIRE(ss.sample_size() == 11);
        REQUIRE(ss.min() == -5.0);
        REQUIRE(ss.max() == 5.0);
        REQUIRE(ss.mean() == 0.0);
        REQUIRE(ss.variance() == (110.0/10.0));
        REQUIRE(ss.population_variance() == (110.0/11.0));
        REQUIRE(ss.std_dev() == std::sqrt(110.0/10.0));
        REQUIRE(ss.std_error() == std::sqrt(110.0/10.0) / std::sqrt(11));
        REQUIRE(ss.skewness() == 0.0);
        REQUIRE(ss.kurtosis() == 1.78);
        REQUIRE(ss.excess_kurtosis() == -1.22);
    }

    SECTION("Testing double summarizer with 12 samples around 0 with skew") {
        SampleSummarizer<double> ss;
        ss.add_sample(0);
        ss.add_sample(1);
        ss.add_sample(-1);
        ss.add_sample(2);
        ss.add_sample(-2);
        ss.add_sample(3);
        ss.add_sample(-3);
        ss.add_sample(4);
        ss.add_sample(-4);
        ss.add_sample(5);
        ss.add_sample(-5);
        ss.add_sample(5);
        REQUIRE(ss.sample_size() == 12);
        REQUIRE(ss.min() == -5);
        REQUIRE(ss.max() == 5);
        REQUIRE(ss.mean() == 5.0/12.0);
        REQUIRE(ss.variance() == Approx(12.0833333));
        REQUIRE(ss.population_variance() == Approx(11.07639));
        REQUIRE(ss.std_dev() == Approx(3.476109));
        REQUIRE(ss.std_error() == Approx(1.003466));
        REQUIRE(ss.skewness() == Approx(-0.09497610248616362));
        REQUIRE(ss.kurtosis() == Approx(1.7077461896011248));
        REQUIRE(ss.excess_kurtosis() == Approx(-1.2922538103988752));
    }
}

TEST_CASE("Testing int SampleSummarizer", "[stats_util]") {

    SECTION("Testing empty int summarizer") {
        SampleSummarizer<int> ss;
        REQUIRE(ss.sample_size() == 0);
        REQUIRE_THROWS_AS(ss.min(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.max(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.mean(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.variance(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.population_variance(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.std_dev(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.std_error(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.skewness(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.kurtosis(), EcoevolityError);
        REQUIRE_THROWS_AS(ss.excess_kurtosis(), EcoevolityError);
    }

    SECTION("Testing int summarizer with one sample 0") {
        SampleSummarizer<int> ss;
        ss.add_sample(0);
        REQUIRE(ss.sample_size() == 1);
        REQUIRE(ss.min() == 0);
        REQUIRE(ss.max() == 0);
        REQUIRE(ss.mean() == 0.0);
        REQUIRE(ss.variance() == std::numeric_limits<double>::infinity());
        REQUIRE(ss.population_variance() == 0.0);
        REQUIRE(ss.std_dev() == std::numeric_limits<double>::infinity());
        REQUIRE(ss.std_error() == std::numeric_limits<double>::infinity());
        REQUIRE(std::isnan(ss.skewness()));
        REQUIRE(std::isnan(ss.kurtosis()));
        REQUIRE(std::isnan(ss.excess_kurtosis()));
    }

    SECTION("Testing int summarizer with one sample -10") {
        SampleSummarizer<int> ss;
        ss.add_sample(-10);
        REQUIRE(ss.sample_size() == 1);
        REQUIRE(ss.min() == -10);
        REQUIRE(ss.max() == -10);
        REQUIRE(ss.mean() == -10.0);
        REQUIRE(ss.variance() == std::numeric_limits<double>::infinity());
        REQUIRE(ss.population_variance() == 0.0);
        REQUIRE(ss.std_dev() == std::numeric_limits<double>::infinity());
        REQUIRE(ss.std_error() == std::numeric_limits<double>::infinity());
        REQUIRE(std::isnan(ss.skewness()));
        REQUIRE(std::isnan(ss.kurtosis()));
        REQUIRE(std::isnan(ss.excess_kurtosis()));
    }

    SECTION("Testing int summarizer with one sample 10") {
        SampleSummarizer<int> ss;
        ss.add_sample(10);
        REQUIRE(ss.sample_size() == 1);
        REQUIRE(ss.min() == 10);
        REQUIRE(ss.max() == 10);
        REQUIRE(ss.mean() == 10.0);
        REQUIRE(ss.variance() == std::numeric_limits<double>::infinity());
        REQUIRE(ss.population_variance() == 0.0);
        REQUIRE(ss.std_dev() == std::numeric_limits<double>::infinity());
        REQUIRE(ss.std_error() == std::numeric_limits<double>::infinity());
        REQUIRE(std::isnan(ss.skewness()));
        REQUIRE(std::isnan(ss.kurtosis()));
        REQUIRE(std::isnan(ss.excess_kurtosis()));
    }

    SECTION("Testing int summarizer with 11 samples around 0") {
        SampleSummarizer<int> ss;
        ss.add_sample(0);
        ss.add_sample(1);
        ss.add_sample(-1);
        ss.add_sample(2);
        ss.add_sample(-2);
        ss.add_sample(3);
        ss.add_sample(-3);
        ss.add_sample(4);
        ss.add_sample(-4);
        ss.add_sample(5);
        ss.add_sample(-5);
        REQUIRE(ss.sample_size() == 11);
        REQUIRE(ss.min() == -5);
        REQUIRE(ss.max() == 5);
        REQUIRE(ss.mean() == 0.0);
        REQUIRE(ss.variance() == (110.0/10.0));
        REQUIRE(ss.population_variance() == (110.0/11.0));
        REQUIRE(ss.std_dev() == std::sqrt(110.0/10.0));
        REQUIRE(ss.std_error() == std::sqrt(110.0/10.0) / std::sqrt(11));
        REQUIRE(ss.skewness() == 0.0);
        REQUIRE(ss.kurtosis() == 1.78);
        REQUIRE(ss.excess_kurtosis() == -1.22);
    }

    SECTION("Testing int summarizer with 12 samples around 0 with skew") {
        SampleSummarizer<int> ss;
        ss.add_sample(0);
        ss.add_sample(1);
        ss.add_sample(-1);
        ss.add_sample(2);
        ss.add_sample(-2);
        ss.add_sample(3);
        ss.add_sample(-3);
        ss.add_sample(4);
        ss.add_sample(-4);
        ss.add_sample(5);
        ss.add_sample(-5);
        ss.add_sample(5);
        REQUIRE(ss.sample_size() == 12);
        REQUIRE(ss.min() == -5);
        REQUIRE(ss.max() == 5);
        REQUIRE(ss.mean() == 5.0/12.0);
        REQUIRE(ss.variance() == Approx(12.0833333));
        REQUIRE(ss.population_variance() == Approx(11.07639));
        REQUIRE(ss.std_dev() == Approx(3.476109));
        REQUIRE(ss.std_error() == Approx(1.003466));
        REQUIRE(ss.skewness() == Approx(-0.09497610248616362));
        REQUIRE(ss.kurtosis() == Approx(1.7077461896011248));
        REQUIRE(ss.excess_kurtosis() == Approx(-1.2922538103988752));
    }
}
