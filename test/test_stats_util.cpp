#include "catch.hpp"
#include "ecoevolity/stats_util.hpp"
#include "ecoevolity/rng.hpp"

TEST_CASE("Testing double SampleSummarizer", "[stats_util]") {

    SECTION("Testing empty double summarizer") {
        SampleSummarizer<double> ss;
        REQUIRE(ss.sample_size() == 0);
        REQUIRE_THROWS_AS(ss.min(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.max(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.mean(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.variance(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.population_variance(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.std_dev(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.std_error(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.skewness(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.kurtosis(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.excess_kurtosis(), EcoevolityError &);
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
        REQUIRE_THROWS_AS(ss.min(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.max(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.mean(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.variance(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.population_variance(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.std_dev(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.std_error(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.skewness(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.kurtosis(), EcoevolityError &);
        REQUIRE_THROWS_AS(ss.excess_kurtosis(), EcoevolityError &);
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

TEST_CASE("Testing effective_sample_size", "[stats_util]") {

    SECTION("Testing effective_sample_size") {
        std::vector<double>x = {
                0.7977061294666541, 0.9150350307910423, 0.7209626707423714,
                0.5954848559944081, 0.18032194756853182, 0.210042410144069,
                0.3673333965443635, 0.8740467791825761, 0.6874289295702046,
                0.22144353794416716, 0.3233467553676893, 0.10398479380458114,
                0.5243615565040305, 0.5877894894599294, 0.42089823773318724,
                0.6266108731616019, 0.3343859686141625, 0.512551474670303,
                0.6446230257104236, 0.36282234951752024, 0.6228723575494212,
                0.7568718761184856, 0.3718316658814024, 0.6861537858829704,
                0.1257109245390987, 0.6412426639048084, 0.48211219814972295,
                0.593973829940721, 0.4036132973697879, 0.42477867300229544,
                0.31213513805943194, 0.7963559245316685, 0.6941826857155579,
                0.6805456463190873, 0.49143482763009017, 0.6290575158052324,
                0.08661756315679514, 0.7995156973771527, 0.27539069568104,
                0.3139293111140057, 0.32288336271807183, 0.2612070751385418,
                0.4545704301079062, 0.6359171147861155, 0.3737093467417866,
                0.9232642159501455, 0.8271543021690014, 0.34958286197540656,
                0.34815266170044323, 0.6056828909353177, 0.5011441473017468,
                0.8184372611091862, 0.06710536859043326, 0.019983484122365947,
                0.3176095570458911, 0.9800154385339, 0.5319803418547973,
                0.2523950819849151, 0.04169284733227552, 0.5240020836881362,
                0.040929832798068166, 0.5024077861662805, 0.7176655502585366,
                0.6306537858831496, 0.5774716670659389, 0.9104292864296849,
                0.35302437929192343, 0.8624334312505447, 0.6990861575487167,
                0.8394941343135478, 0.5795304077084198, 0.12535068024747653,
                0.7025132099214821, 0.177220279120623, 0.9070732428670005,
                0.7666417009611808, 0.08750652002252135, 0.9948532901833365,
                0.44265582277400917, 0.10322490371849158, 0.5094288068541217,
                0.13640416841602576, 0.20328541281100587, 0.7289261198868512,
                0.8040861608469766, 0.9670617517210303, 0.23243617749946088,
                0.25068739997092004, 0.2742590187495584, 0.307652725552081,
                0.8997811130977051, 0.35615376615317706, 0.0211059298791072,
                0.03263965076194353, 0.4416542975034954, 0.5586675733736068,
                0.21167935845287156, 0.47810451475326077, 0.7395889690656308,
                0.24135469373818985};
        // expected results calculated with R package mcmcse: Monte Carlo
        // Standard Errors for MCMC Version 1.2-1 created on 2016-03-24.
        double ess = effective_sample_size<double>(x, false);
        REQUIRE(ess == Approx(154.581627605617));

        ess = effective_sample_size<double>(x, true);
        REQUIRE(ess == Approx(100.0));
    }
}

TEST_CASE("Testing potential_scale_reduction_factor double", "[stats_util]") {
    SECTION("Testing potential_scale_reduction_factor<double>") {
        std::vector< std::vector<double> > chains {
            {1.1, 1.3, 1.2, 1.6, 1.5},
            {1.2, 1.7, 1.5, 1.9, 1.6}};
        double psrf = potential_scale_reduction_factor<double>(chains);
        // expectation calculated with commit aa83c8cc8584ba2d
        // of pymc.diagnostics.gelman_rubin
        // <https://github.com/pymc-devs/pymc/blob/master/pymc/diagnostics.py>
        double e_pymc = 1.2591483413222384;
        REQUIRE(psrf == Approx(e_pymc));
    }
}

TEST_CASE("Testing potential_scale_reduction_factor double, equal chains",
        "[stats_util]") {
    SECTION("Testing potential_scale_reduction_factor<double>, equal cains") {
        std::vector< std::vector<double> > chains {
            {1.1, 1.3, 1.2, 1.6, 1.5},
            {1.1, 1.3, 1.2, 1.6, 1.5}};
        double psrf = potential_scale_reduction_factor<double>(chains);
        // expectation calculated with commit aa83c8cc8584ba2d
        // of pymc.diagnostics.gelman_rubin
        // <https://github.com/pymc-devs/pymc/blob/master/pymc/diagnostics.py>
        double e_pymc = 0.89442719099991586;
        REQUIRE(psrf == Approx(e_pymc));
    }
}

TEST_CASE("Testing median with one double",
        "[stats_util]") {
    SECTION("Testing median with one double") {
        std::vector<double> samples {1.3};
        double m = get_median<double>(samples);
        REQUIRE(m == samples.at(0));
    }
}

TEST_CASE("Testing median with four doubles",
        "[stats_util]") {
    SECTION("Testing median with four doubles") {
        std::vector<double> samples {1.2, 0.1, 23.4, 1.3};
        double m = get_median<double>(samples);
        REQUIRE(m == Approx(1.25));
    }
}

TEST_CASE("Testing median with five doubles",
        "[stats_util]") {
    SECTION("Testing median with five doubles") {
        std::vector<double> samples {1.2, 0.1, 1.4, 23.1, 1.3};
        double m = get_median<double>(samples);
        REQUIRE(m == Approx(1.3));
    }
}

TEST_CASE("Testing median with one int",
        "[stats_util]") {
    SECTION("Testing median with one int") {
        std::vector<int> samples {13};
        double m = get_median<int>(samples);
        REQUIRE(m == samples.at(0));
    }
}

TEST_CASE("Testing median with four ints",
        "[stats_util]") {
    SECTION("Testing median with four doubles") {
        std::vector<int> samples {12, 1, 234, 13};
        double m = get_median<int>(samples);
        REQUIRE(m == Approx(12.5));
    }
}

TEST_CASE("Testing median with five ints",
        "[stats_util]") {
    SECTION("Testing median with five ints") {
        std::vector<int> samples {12, 1, 14, 231, 13};
        double m = get_median<int>(samples);
        REQUIRE(m == Approx(13));
    }
}

TEST_CASE("Testing HPD interval, quantile, percentile, and summary") {
    RandomNumberGenerator rng = RandomNumberGenerator(123);
    unsigned int nsamples = 1000000;
    std::vector<double> normal_samples(nsamples);
    std::vector<double> exponential_samples(nsamples);
    for (unsigned int i = 0; i < nsamples; ++i) {
        normal_samples.at(i) = rng.normal(0.0, 1.0);
        exponential_samples.at(i) = rng.gamma(1.0, 1.0);
    }

    SECTION("Testing standard normal HPD") {
        double eps = 0.002;
        std::pair<double, double> hpdi = get_hpd_interval<double>(normal_samples, 0.95);
        REQUIRE(hpdi.first == Approx(-1.96).epsilon(eps));
        REQUIRE(hpdi.second == Approx(1.96).epsilon(eps));
    }

    SECTION("Testing normal quantiles") {
        double eps = 0.002;
        std::pair<double, double> quants = quantiles_95<double>(normal_samples);
        double q025 = quantile<double>(normal_samples, 0.025);
        double q975 = quantile<double>(normal_samples, 0.975);
        REQUIRE(quants.first == Approx(q025));
        REQUIRE(quants.second == Approx(q975));
        REQUIRE(q025 == Approx(-1.96).epsilon(eps));
        REQUIRE(q975 == Approx(1.96).epsilon(eps));
    }

    SECTION("Testing exponential HPD") {
        double eps = 0.002;
        std::pair<double, double> hpdi = get_hpd_interval<double>(exponential_samples, 0.95);
        REQUIRE(hpdi.first == Approx(0.0).epsilon(eps));
        REQUIRE(hpdi.second == Approx(2.9957).epsilon(eps));
    }

    SECTION("Testing exponential quantiles") {
        double eps = 0.002;
        std::pair<double, double> quants = quantiles_95<double>(exponential_samples);
        double q025 = quantile<double>(exponential_samples, 0.025);
        double q975 = quantile<double>(exponential_samples, 0.975);
        REQUIRE(quants.first == Approx(q025));
        REQUIRE(quants.second == Approx(q975));
        REQUIRE(q025 == Approx(0.0253).epsilon(eps));
        REQUIRE(q975 == Approx(3.6889).epsilon(eps));
    }

    SECTION("Testing standard normal percentile") {
        double eps = 0.002;
        double r = percentile<double>(normal_samples, -1.96);
        REQUIRE(r == Approx(0.025).epsilon(eps));
        r = percentile<double>(normal_samples, 1.96);
        REQUIRE(r == Approx(0.975).epsilon(eps));
    }

    SECTION("Testing exponential percentile") {
        double eps = 0.002;
        double r = percentile<double>(exponential_samples, 0.0253);
        REQUIRE(r == Approx(0.025).epsilon(eps));
        r = percentile<double>(exponential_samples, 3.6889);
        REQUIRE(r == Approx(0.975).epsilon(eps));
    }

    SECTION("Testing standard normal SampleSummary") {
        double expected_min = std::numeric_limits<double>::max();
        double expected_max = -std::numeric_limits<double>::max();
        for (auto x : normal_samples) {
            if (x > expected_max) {
                expected_max = x;
            }
            if (x < expected_min) {
                expected_min = x;
            }
        }
        double eps = 0.002;
        SampleSummary<double> ss(normal_samples);
        REQUIRE(ss.sample_size() == normal_samples.size());
        REQUIRE(ss.min() == expected_min);
        REQUIRE(ss.max() == expected_max);
        REQUIRE(ss.mean() == Approx(0.0).epsilon(eps));
        REQUIRE(ss.median() == Approx(0.0).epsilon(eps));
        REQUIRE(ss.variance() == Approx(1.0).epsilon(eps));
        REQUIRE(ss.std_dev() == Approx(1.0).epsilon(eps));
        REQUIRE(ss.std_error() == Approx(1.0 / std::sqrt(normal_samples.size())).epsilon(eps));
        REQUIRE(ss.hpdi_95().first == Approx(-1.96).epsilon(eps));
        REQUIRE(ss.hpdi_95().second == Approx(1.96).epsilon(eps));
        REQUIRE(ss.qi_95().first == Approx(-1.96).epsilon(eps));
        REQUIRE(ss.qi_95().second == Approx(1.96).epsilon(eps));
    }
}
