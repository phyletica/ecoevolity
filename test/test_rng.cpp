#include "catch.hpp"
#include "ecoevolity/rng.hpp"
#include "ecoevolity/math_util.hpp"
#include "ecoevolity/stats_util.hpp"

#include <limits>
#include <map>
#include <set>

TEST_CASE("Testing seed constructor", "[RandomNumberGenerator]") {

    SECTION("Testing seed constructor") {
        RandomNumberGenerator rng = RandomNumberGenerator(123);
        REQUIRE(rng.get_seed() == 123);
    }
}

TEST_CASE("Testing bare constructor", "[RandomNumberGenerator]") {

    SECTION("Testing seed constructor") {
        RandomNumberGenerator rng = RandomNumberGenerator();
        REQUIRE(rng.get_seed() > 0);
    }
}

TEST_CASE("Testing uniform_real", "[RandomNumberGenerator]") {

    SECTION("Testing uniform_real()") {
        RandomNumberGenerator rng(123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = rng.uniform_real();
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(0.5).epsilon(0.001));
        REQUIRE(variance == Approx(1.0/12.0).epsilon(0.001));
        REQUIRE(mn >= 0.0);
        REQUIRE(mx < 1.0);
        REQUIRE(mn == Approx(0.0).epsilon(0.001));
        REQUIRE(mx == Approx(1.0).epsilon(0.001));
    }

    SECTION("Testing uniform_real(1.0, 3.0)") {
        RandomNumberGenerator rng(2123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = rng.uniform_real(1.0, 3.0);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            ++n;
            d = x - mean;
            d_n = d / n;
            mean += d_n;
            sum_devs += d * d_n * (n - 1);
        }
        double variance = sum_devs / (n - 1);
        
        REQUIRE(mean == Approx(2.0).epsilon(0.001));
        REQUIRE(variance == Approx(4.0/12.0).epsilon(0.001));
        REQUIRE(mn >= 1.0);
        REQUIRE(mx < 3.0);
        REQUIRE(mn == Approx(1.0).epsilon(0.001));
        REQUIRE(mx == Approx(3.0).epsilon(0.001));
    }
}

TEST_CASE("Testing uniform_int", "[RandomNumberGenerator]") {
    SECTION("Testing uniform_int") {
        RandomNumberGenerator rng(98729);
        unsigned int two_count = 0;
        unsigned int three_count = 0;
        for (unsigned int i = 0; i < 100000; ++i) {
            int x = rng.uniform_int(2,3);
            if (x == 2) {
                ++two_count;
            }
            else if (x == 3) {
                ++three_count;
            }
        }
        REQUIRE(two_count + three_count == 100000);
        REQUIRE(two_count / 100000.0 == Approx(0.5).epsilon(0.001));
    }
}

TEST_CASE("Testing uniform_positive_int", "[RandomNumberGenerator]") {
    SECTION("Testing uniform_positive_int with one arg") {
        RandomNumberGenerator rng(37229);
        unsigned int n0 = 0;
        unsigned int n1 = 0;
        unsigned int n2 = 0;
        for (unsigned int i = 0; i < 100000; ++i) {
            int x = rng.uniform_positive_int(2);
            if (x == 2) {
                ++n2;
            }
            else if (x == 1) {
                ++n1;
            }
            else if (x == 0) {
                ++n0;
            }
        }
        REQUIRE(n0 + n1 + n2 == 100000);
        REQUIRE(n0 / 100000.0 == Approx(1.0/3.0).epsilon(0.001));
    }

    SECTION("Testing uniform_positive_int with two args") {
        RandomNumberGenerator rng(37229);
        unsigned int n4 = 0;
        unsigned int n5 = 0;
        unsigned int n6 = 0;
        for (unsigned int i = 0; i < 100000; ++i) {
            int x = rng.uniform_positive_int(4,6);
            if (x == 6) {
                ++n6;
            }
            else if (x == 5) {
                ++n5;
            }
            else if (x == 4) {
                ++n4;
            }
        }
        REQUIRE(n4 + n5 + n6 == 100000);
        REQUIRE(n4 / 100000.0 == Approx(1.0/3.0).epsilon(0.001));
    }
}

TEST_CASE("Testing weighted_index", "[RandomNumberGenerator]") {

    SECTION("Testing weighted_index with even weights") {
        RandomNumberGenerator rng(11111);
        std::vector<double> weights = {0.25, 0.25, 0.25, 0.25};
        std::vector<unsigned int> counts(4, 0);

        for (unsigned int i = 0; i < 100000; ++i) {
            unsigned int x = rng.weighted_index(weights);
            ++counts.at(x);
        }

        REQUIRE(counts.at(0) / 100000.0 == Approx(0.25).epsilon(0.001));
    }

    SECTION("Testing weighted_index with uneven weights") {
        RandomNumberGenerator rng(11111);
        std::vector<double> weights = {0.2, 0.5, 0.3};
        std::vector<unsigned int> counts(3, 0);

        for (unsigned int i = 0; i < 100000; ++i) {
            unsigned int x = rng.weighted_index(weights);
            ++counts.at(x);
        }

        REQUIRE(counts.at(1) / 100000.0 == Approx(0.5).epsilon(0.001));
    }
}

TEST_CASE("Testing random_string", "[RandomNumberGenerator]") {
    SECTION("Testing expectation for first char in pool") {
        RandomNumberGenerator rng(11111);

        std::string s;
        unsigned int count = 0;
        unsigned int total = 0;
        unsigned int reps = 100000;
        unsigned int str_len = 10;
        for (unsigned int i = 0; i < reps; ++i) {
            s = rng.random_string(str_len);
            for (unsigned int j = 0; j < s.size(); ++j) {
                if (s.at(j) == '0') {
                    ++count;
                }
            }
            total += s.size();
        }
        REQUIRE(total == (reps * str_len));

        double p = count / (double)total;
        double e =  1.0 / 62.0;
        REQUIRE(p == Approx(e).epsilon(0.001));
    }

    SECTION("Testing expectation for last char in pool") {
        RandomNumberGenerator rng(11111);

        std::string s;
        unsigned int count = 0;
        unsigned int total = 0;
        unsigned int reps = 100000;
        unsigned int str_len = 10;
        for (unsigned int i = 0; i < reps; ++i) {
            s = rng.random_string(str_len);
            for (unsigned int j = 0; j < s.size(); ++j) {
                if (s.at(j) == 'Z') {
                    ++count;
                }
            }
            total += s.size();
        }
        REQUIRE(total == (reps * str_len));

        double p = count / (double)total;
        double e =  1.0 / 62.0;
        REQUIRE(p == Approx(e).epsilon(0.001));
    }
}

TEST_CASE("Testing random_subset_indices 1 from 2", "[RandomNumberGenerator]") {
    SECTION("Testing random sampling") {
        unsigned int number_of_elements = 2;
        unsigned int subset_size = 1;
        unsigned int nreps = 100000;
        RandomNumberGenerator rng(872349);
        std::vector<unsigned int> counts(number_of_elements, 0);

        std::vector<unsigned int> indices;
        for (unsigned int i = 0; i < nreps; ++i) {
            indices = rng.random_subset_indices(number_of_elements, subset_size);
            for (auto const idx: indices) {
                ++counts.at(idx);
            }
        }

        for (auto const cnt: counts) {
            REQUIRE(cnt / (double)(nreps * subset_size) == Approx(1.0/number_of_elements).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing random_subset_indices 1 from 10", "[RandomNumberGenerator]") {
    SECTION("Testing random sampling") {
        unsigned int number_of_elements = 10;
        unsigned int subset_size = 1;
        unsigned int nreps = 100000;
        RandomNumberGenerator rng(1234);
        std::vector<unsigned int> counts(number_of_elements, 0);

        std::vector<unsigned int> indices;
        for (unsigned int i = 0; i < nreps; ++i) {
            indices = rng.random_subset_indices(number_of_elements, subset_size);
            for (auto const idx: indices) {
                ++counts.at(idx);
            }
        }

        for (auto const cnt: counts) {
            REQUIRE(cnt / (double)(nreps * subset_size) == Approx(1.0/number_of_elements).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing random_subset_indices 2 from 10", "[RandomNumberGenerator]") {
    SECTION("Testing random sampling") {
        unsigned int number_of_elements = 10;
        unsigned int subset_size = 2;
        unsigned int nreps = 100000;
        RandomNumberGenerator rng(1234);
        std::vector<unsigned int> counts(number_of_elements, 0);

        std::vector<unsigned int> indices;
        for (unsigned int i = 0; i < nreps; ++i) {
            indices = rng.random_subset_indices(number_of_elements, subset_size);
            for (auto const idx: indices) {
                ++counts.at(idx);
            }
        }

        for (auto const cnt: counts) {
            REQUIRE(cnt / (double)(nreps * subset_size) == Approx(1.0/number_of_elements).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing random_subset_indices 3 from 10", "[RandomNumberGenerator]") {
    SECTION("Testing random sampling") {
        unsigned int number_of_elements = 10;
        unsigned int subset_size = 3;
        unsigned int nreps = 100000;
        RandomNumberGenerator rng(1234);
        std::vector<unsigned int> counts(number_of_elements, 0);

        std::vector<unsigned int> indices;
        for (unsigned int i = 0; i < nreps; ++i) {
            indices = rng.random_subset_indices(number_of_elements, subset_size);
            for (auto const idx: indices) {
                ++counts.at(idx);
            }
        }

        for (auto const cnt: counts) {
            REQUIRE(cnt / (double)(nreps * subset_size) == Approx(1.0/number_of_elements).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing random_subset_indices 4 from 10", "[RandomNumberGenerator]") {
    SECTION("Testing random sampling") {
        unsigned int number_of_elements = 10;
        unsigned int subset_size = 4;
        unsigned int nreps = 100000;
        RandomNumberGenerator rng(1234);
        std::vector<unsigned int> counts(number_of_elements, 0);

        std::vector<unsigned int> indices;
        for (unsigned int i = 0; i < nreps; ++i) {
            indices = rng.random_subset_indices(number_of_elements, subset_size);
            for (auto const idx: indices) {
                ++counts.at(idx);
            }
        }

        for (auto const cnt: counts) {
            REQUIRE(cnt / (double)(nreps * subset_size) == Approx(1.0/number_of_elements).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing random_subset_indices 5 from 10", "[RandomNumberGenerator]") {
    SECTION("Testing random sampling") {
        unsigned int number_of_elements = 10;
        unsigned int subset_size = 5;
        unsigned int nreps = 100000;
        RandomNumberGenerator rng(1234);
        std::vector<unsigned int> counts(number_of_elements, 0);

        std::vector<unsigned int> indices;
        for (unsigned int i = 0; i < nreps; ++i) {
            indices = rng.random_subset_indices(number_of_elements, subset_size);
            for (auto const idx: indices) {
                ++counts.at(idx);
            }
        }

        for (auto const cnt: counts) {
            REQUIRE(cnt / (double)(nreps * subset_size) == Approx(1.0/number_of_elements).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing random_subset_indices 6 from 10", "[RandomNumberGenerator]") {
    SECTION("Testing random sampling") {
        unsigned int number_of_elements = 10;
        unsigned int subset_size = 6;
        unsigned int nreps = 100000;
        RandomNumberGenerator rng(1234);
        std::vector<unsigned int> counts(number_of_elements, 0);

        std::vector<unsigned int> indices;
        for (unsigned int i = 0; i < nreps; ++i) {
            indices = rng.random_subset_indices(number_of_elements, subset_size);
            for (auto const idx: indices) {
                ++counts.at(idx);
            }
        }

        for (auto const cnt: counts) {
            REQUIRE(cnt / (double)(nreps * subset_size) == Approx(1.0/number_of_elements).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing random_subset_indices 7 from 10", "[RandomNumberGenerator]") {
    SECTION("Testing random sampling") {
        unsigned int number_of_elements = 10;
        unsigned int subset_size = 7;
        unsigned int nreps = 100000;
        RandomNumberGenerator rng(1234);
        std::vector<unsigned int> counts(number_of_elements, 0);

        std::vector<unsigned int> indices;
        for (unsigned int i = 0; i < nreps; ++i) {
            indices = rng.random_subset_indices(number_of_elements, subset_size);
            for (auto const idx: indices) {
                ++counts.at(idx);
            }
        }

        for (auto const cnt: counts) {
            REQUIRE(cnt / (double)(nreps * subset_size) == Approx(1.0/number_of_elements).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing random_subset_indices 8 from 10", "[RandomNumberGenerator]") {
    SECTION("Testing random sampling") {
        unsigned int number_of_elements = 10;
        unsigned int subset_size = 8;
        unsigned int nreps = 100000;
        RandomNumberGenerator rng(1234);
        std::vector<unsigned int> counts(number_of_elements, 0);

        std::vector<unsigned int> indices;
        for (unsigned int i = 0; i < nreps; ++i) {
            indices = rng.random_subset_indices(number_of_elements, subset_size);
            for (auto const idx: indices) {
                ++counts.at(idx);
            }
        }

        for (auto const cnt: counts) {
            REQUIRE(cnt / (double)(nreps * subset_size) == Approx(1.0/number_of_elements).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing random_subset_indices 9 from 10", "[RandomNumberGenerator]") {
    SECTION("Testing random sampling") {
        unsigned int number_of_elements = 10;
        unsigned int subset_size = 9;
        unsigned int nreps = 100000;
        RandomNumberGenerator rng(1234);
        std::vector<unsigned int> counts(number_of_elements, 0);

        std::vector<unsigned int> indices;
        for (unsigned int i = 0; i < nreps; ++i) {
            indices = rng.random_subset_indices(number_of_elements, subset_size);
            for (auto const idx: indices) {
                ++counts.at(idx);
            }
        }

        for (auto const cnt: counts) {
            REQUIRE(cnt / (double)(nreps * subset_size) == Approx(1.0/number_of_elements).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing beta(1,1)", "[RandomNumberGenerator]") {

    SECTION("Testing beta(1,1)") {
        double a = 1.0;
        double b = 1.0;
        double expected_mean = a / (a + b);
        double expected_variance = (a * b) / ((a + b) * (a + b) * (a + b + 1.0));

        RandomNumberGenerator rng(123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = rng.beta(a, b);
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
}

TEST_CASE("Testing beta(0.5,0.5)", "[RandomNumberGenerator]") {

    SECTION("Testing beta(1,1)") {
        double a = 0.5;
        double b = 0.5;
        double expected_mean = a / (a + b);
        double expected_variance = (a * b) / ((a + b) * (a + b) * (a + b + 1.0));

        RandomNumberGenerator rng(123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = rng.beta(a, b);
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
}

TEST_CASE("Testing beta(5,1)", "[RandomNumberGenerator]") {

    SECTION("Testing beta(5,1)") {
        double a = 5.0;
        double b = 1.0;
        double expected_mean = a / (a + b);
        double expected_variance = (a * b) / ((a + b) * (a + b) * (a + b + 1.0));

        RandomNumberGenerator rng(123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = rng.beta(a, b);
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
        // Unlikely to sample near zero
        // REQUIRE(mn == Approx(0.0).epsilon(0.001));
        REQUIRE(mx == Approx(1.0).epsilon(0.001));
    }
}

TEST_CASE("Testing beta(1,5)", "[RandomNumberGenerator]") {

    SECTION("Testing beta(1,5)") {
        double a = 1.0;
        double b = 5.0;
        double expected_mean = a / (a + b);
        double expected_variance = (a * b) / ((a + b) * (a + b) * (a + b + 1.0));

        RandomNumberGenerator rng(123);
        unsigned int n = 0;
        double mean = 0.0;
        double sum_devs = 0.0;
        double d;
        double d_n;
        double mn = std::numeric_limits<double>::max();
        double mx = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 100000; ++i) {
            double x = rng.beta(a, b);
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
        // Unlikely to sample near one
        // REQUIRE(mx == Approx(1.0).epsilon(0.001));
    }
}

TEST_CASE("Testing dirichet_process(6, 1.7)", "[RandomNumberGenerator]") {

    SECTION("Testing dirichlet_process(6, 1.7)") {
        unsigned int nsamples = 100000;
        unsigned int n = 6;
        double a = 1.7;

        std::vector<unsigned int> elements (6, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
                {4, 0},
                {5, 0},
                {6, 0}
        };

        RandomNumberGenerator rng(123);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.dirichlet_process(elements, a);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012345") == nevent_counts.at(6));

        for (std::string m: {"000000", "012345", "012344", "012340", "001234", "012314"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_dpp_log_prior_probability(m, a))).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing dirichet_process(3, 0.6)", "[RandomNumberGenerator]") {

    SECTION("Testing dirichlet_process(3, 0.6)") {
        unsigned int nsamples = 200000;
        unsigned int n = 3;
        double a = 0.6;

        std::vector<unsigned int> elements (3, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
        };

        RandomNumberGenerator rng(342254);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.dirichlet_process(elements, a);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_dpp_log_prior_probability(m, a))).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing pitman_yor_process(3, 0.6, 0.0)", "[RandomNumberGenerator]") {

    SECTION("Testing pitman_yor_process(3, 0.6, 0.0)") {
        unsigned int nsamples = 200000;
        unsigned int n = 3;
        double a = 0.6;
        double d = 0.0;

        std::vector<unsigned int> elements (3, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
        };

        RandomNumberGenerator rng(342254);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.pitman_yor_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_dpp_log_prior_probability(m, a))).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing weighted_discount_process(3, 0.6, 0.0)", "[RandomNumberGenerator]") {

    SECTION("Testing weighted_discount_process(3, 0.6, 0.0)") {
        unsigned int nsamples = 200000;
        unsigned int n = 3;
        double a = 0.6;
        double d = 0.0;

        std::vector<unsigned int> elements (3, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
        };

        RandomNumberGenerator rng(342254);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.weighted_discount_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_wdp_log_prior_probability(m, a, d))).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing dirichet_process(3, 4.2)", "[RandomNumberGenerator]") {

    SECTION("Testing dirichlet_process(3, 4.2)") {
        unsigned int nsamples = 200000;
        unsigned int n = 3;
        double a = 4.2;

        std::vector<unsigned int> elements (3, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
        };

        RandomNumberGenerator rng(42342);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.dirichlet_process(elements, a);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_dpp_log_prior_probability(m, a))).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing pitman_yor_process(3, 4.2, 0.0)", "[RandomNumberGenerator]") {

    SECTION("Testing pitman_yor_process(3, 4.2, 0.0)") {
        unsigned int nsamples = 200000;
        unsigned int n = 3;
        double a = 4.2;
        double d = 0.0;

        std::vector<unsigned int> elements (3, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
        };

        RandomNumberGenerator rng(42342);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.pitman_yor_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_dpp_log_prior_probability(m, a))).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing weighted_discount_process(3, 4.2, 0.0)", "[RandomNumberGenerator]") {

    SECTION("Testing weighted_discount_process(3, 4.2, 0.0)") {
        unsigned int nsamples = 200000;
        unsigned int n = 3;
        double a = 4.2;
        double d = 0.0;

        std::vector<unsigned int> elements (3, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
        };

        RandomNumberGenerator rng(42342);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.weighted_discount_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_wdp_log_prior_probability(m, a, d))).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing vector return version of dirichet_process(3, 0.6)",
        "[RandomNumberGenerator]") {

    SECTION("Testing dirichlet_process(3, 0.6)") {
        unsigned int nsamples = 200000;
        unsigned int n = 3;
        double a = 0.6;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
        };

        RandomNumberGenerator rng(342254);
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::pair< unsigned int, std::vector<unsigned int> > ret = rng.dirichlet_process(n, a);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: ret.second) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            unsigned int ncats = ret.first;
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_dpp_log_prior_probability(m, a))).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing vector return version of pitman_yor_process(3, 0.6, 0.0)",
        "[RandomNumberGenerator]") {

    SECTION("Testing pitman_yor_process(3, 0.6)") {
        unsigned int nsamples = 200000;
        unsigned int n = 3;
        double a = 0.6;
        double d = 0.0;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
        };

        RandomNumberGenerator rng(342254);
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::pair< unsigned int, std::vector<unsigned int> > ret = rng.pitman_yor_process(n, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: ret.second) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            unsigned int ncats = ret.first;
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_dpp_log_prior_probability(m, a))).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing vector return version of weighted_discount_process(3, 0.6, 0.0)",
        "[RandomNumberGenerator]") {

    SECTION("Testing weighted_discount_process(3, 0.6)") {
        unsigned int nsamples = 200000;
        unsigned int n = 3;
        double a = 0.6;
        double d = 0.0;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
        };

        RandomNumberGenerator rng(342254);
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::pair< unsigned int, std::vector<unsigned int> > ret = rng.weighted_discount_process(n, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: ret.second) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            unsigned int ncats = ret.first;
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_wdp_log_prior_probability(m, a, d))).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing random_subsets(1, 1)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_subsets(1, 1)") {
        unsigned int nsamples = 50;
        RandomNumberGenerator rng(123);
        std::vector< std::vector<unsigned int> > subsets;
        std::vector< std::vector<unsigned int> > expected_subsets;
        expected_subsets.push_back({0});
        for (unsigned int i = 0; i < nsamples; ++i) {
            subsets = rng.random_subsets(1, 1);
            REQUIRE(subsets == expected_subsets);
        }
    }
}

TEST_CASE("Testing random_subsets(2, 1)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_subsets(2, 1)") {
        unsigned int nsamples = 50;
        RandomNumberGenerator rng(123);
        std::vector< std::vector<unsigned int> > subsets;
        std::vector< std::vector<unsigned int> > expected_subsets;
        expected_subsets.push_back({0, 1});
        for (unsigned int i = 0; i < nsamples; ++i) {
            subsets = rng.random_subsets(2, 1);
            REQUIRE(subsets == expected_subsets);
        }
    }
}

TEST_CASE("Testing random_subsets(2, 2)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_subsets(2, 2)") {
        unsigned int nsamples = 50;
        RandomNumberGenerator rng(123);
        std::vector< std::vector<unsigned int> > subsets;
        std::vector< std::vector<unsigned int> > expected_subsets;
        expected_subsets.push_back({0});
        expected_subsets.push_back({1});
        for (unsigned int i = 0; i < nsamples; ++i) {
            subsets = rng.random_subsets(2, 2);
            REQUIRE(subsets == expected_subsets);
        }
    }
}

TEST_CASE("Testing random_subsets(3, 1)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_subsets(3, 1)") {
        unsigned int nsamples = 50;
        RandomNumberGenerator rng(123);
        std::vector< std::vector<unsigned int> > subsets;
        std::vector< std::vector<unsigned int> > expected_subsets;
        expected_subsets.push_back({0, 1, 2});
        for (unsigned int i = 0; i < nsamples; ++i) {
            subsets = rng.random_subsets(3, 1);
            REQUIRE(subsets == expected_subsets);
        }
    }
}

TEST_CASE("Testing random_subsets(3, 3)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_subsets(3, 3)") {
        unsigned int nsamples = 50;
        RandomNumberGenerator rng(123);
        std::vector< std::vector<unsigned int> > subsets;
        std::vector< std::vector<unsigned int> > expected_subsets;
        expected_subsets.push_back({0});
        expected_subsets.push_back({1});
        expected_subsets.push_back({2});
        for (unsigned int i = 0; i < nsamples; ++i) {
            subsets = rng.random_subsets(3, 3);
            REQUIRE(subsets == expected_subsets);
        }
    }
}

TEST_CASE("Testing random_subsets(3, 2)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_subsets(3, 2)") {
        unsigned int nsamples = 500000;
        RandomNumberGenerator rng(111);
        std::vector< std::vector<unsigned int> > subsets;
        std::vector< std::vector<unsigned int> > expected_subsets_1;
        std::vector< std::vector<unsigned int> > expected_subsets_2;
        std::vector< std::vector<unsigned int> > expected_subsets_3;
        expected_subsets_1.push_back({0, 1});
        expected_subsets_1.push_back({2});
        expected_subsets_2.push_back({0, 2});
        expected_subsets_2.push_back({1});
        expected_subsets_3.push_back({0});
        expected_subsets_3.push_back({1, 2});
        unsigned int count_1 = 0;
        unsigned int count_2 = 0;
        unsigned int count_3 = 0;
        unsigned int bad_count = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            subsets = rng.random_subsets(3, 2);
            if (subsets == expected_subsets_1) {
                ++count_1;
            }
            else if (subsets == expected_subsets_2) {
                ++count_2;
            }
            else if (subsets == expected_subsets_3) {
                ++count_3;
            }
            else {
                ++bad_count;
            }
        }
        REQUIRE(bad_count == 0);
        REQUIRE(count_1 + count_2 + count_3 == nsamples);
        REQUIRE(count_1 / (double)nsamples == Approx(1.0/3.0).epsilon(0.001));
        REQUIRE(count_2 / (double)nsamples == Approx(1.0/3.0).epsilon(0.001));
    }
}

TEST_CASE("Testing random_set_partition_with_k_subsets(1, 1)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition_with_k_subsets(1, 1)") {
        unsigned int nsamples = 50;
        RandomNumberGenerator rng(123);
        std::vector<unsigned int> subsets;
        std::vector<unsigned int> expected_subsets = {0};
        for (unsigned int i = 0; i < nsamples; ++i) {
            subsets = rng.random_set_partition_with_k_subsets(1, 1);
            REQUIRE(subsets == expected_subsets);
        }
    }
}

TEST_CASE("Testing random_set_partition_with_k_subsets(2, 1)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition_with_k_subsets(2, 1)") {
        unsigned int nsamples = 50;
        RandomNumberGenerator rng(123);
        std::vector<unsigned int> subsets;
        std::vector<unsigned int> expected_subsets = {0, 0};
        for (unsigned int i = 0; i < nsamples; ++i) {
            subsets = rng.random_set_partition_with_k_subsets(2, 1);
            REQUIRE(subsets == expected_subsets);
        }
    }
}

TEST_CASE("Testing random_set_partition_with_k_subsets(2, 2)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition_with_k_subsets(2, 2)") {
        unsigned int nsamples = 50;
        RandomNumberGenerator rng(123);
        std::vector<unsigned int> subsets;
        std::vector<unsigned int> expected_subsets = {0, 1};
        for (unsigned int i = 0; i < nsamples; ++i) {
            subsets = rng.random_set_partition_with_k_subsets(2, 2);
            REQUIRE(subsets == expected_subsets);
        }
    }
}

TEST_CASE("Testing random_set_partition_with_k_subsets(3, 1)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition_with_k_subsets(3, 1)") {
        unsigned int nsamples = 50;
        RandomNumberGenerator rng(123);
        std::vector<unsigned int> subsets;
        std::vector<unsigned int> expected_subsets = {0, 0, 0};
        for (unsigned int i = 0; i < nsamples; ++i) {
            subsets = rng.random_set_partition_with_k_subsets(3, 1);
            REQUIRE(subsets == expected_subsets);
        }
    }
}

TEST_CASE("Testing random_set_partition_with_k_subsets(3, 3)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition_with_k_subsets(3, 3)") {
        unsigned int nsamples = 50;
        RandomNumberGenerator rng(123);
        std::vector<unsigned int> subsets;
        std::vector<unsigned int> expected_subsets = {0, 1, 2};
        for (unsigned int i = 0; i < nsamples; ++i) {
            subsets = rng.random_set_partition_with_k_subsets(3, 3);
            REQUIRE(subsets == expected_subsets);
        }
    }
}


TEST_CASE("Testing random_set_partition_with_k_subsets(3, 2)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition_with_k_subsets(3, 2)") {
        unsigned int nsamples = 500000;
        RandomNumberGenerator rng(111);
        std::vector<unsigned int> subsets;
        std::vector<unsigned int> expected_subsets_1 = {0, 0, 1};
        std::vector<unsigned int> expected_subsets_2 = {0, 1, 1};
        std::vector<unsigned int> expected_subsets_3 = {0, 1, 0};
        unsigned int count_1 = 0;
        unsigned int count_2 = 0;
        unsigned int count_3 = 0;
        unsigned int bad_count = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            subsets = rng.random_set_partition_with_k_subsets(3, 2);
            if (subsets == expected_subsets_1) {
                ++count_1;
            }
            else if (subsets == expected_subsets_2) {
                ++count_2;
            }
            else if (subsets == expected_subsets_3) {
                ++count_3;
            }
            else {
                ++bad_count;
            }
        }
        REQUIRE(bad_count == 0);
        REQUIRE(count_1 + count_2 + count_3 == nsamples);
        REQUIRE(count_1 / (double)nsamples == Approx(1.0/3.0).epsilon(0.001));
        REQUIRE(count_2 / (double)nsamples == Approx(1.0/3.0).epsilon(0.001));
    }
}

TEST_CASE("Testing random_set_partition (1, 1.0)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition(1, 1.0)") {
        unsigned int nsamples = 100;
        unsigned int n = 1;
        double split_weight = 1.0;

        std::vector<unsigned int> expected_model(1, 0);
        unsigned int expected_nevents = 1;

        RandomNumberGenerator rng(654321);
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::pair< unsigned int, std::vector<unsigned int> > ret = rng.random_set_partition(n, split_weight);
            REQUIRE(ret.second == expected_model);
            REQUIRE(ret.first == expected_nevents);
        }
    }
}

TEST_CASE("Testing random_set_partition (1, 124.9)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition(1, 124.9)") {
        unsigned int nsamples = 100;
        unsigned int n = 1;
        double split_weight = 124.9;

        std::vector<unsigned int> expected_model(1, 0);
        unsigned int expected_nevents = 1;

        RandomNumberGenerator rng(654321);
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::pair< unsigned int, std::vector<unsigned int> > ret = rng.random_set_partition(n, split_weight);
            REQUIRE(ret.second == expected_model);
            REQUIRE(ret.first == expected_nevents);
        }
    }
}

TEST_CASE("Testing random_set_partition (2, 1.0)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition(2, 1.0)") {
        unsigned int nsamples = 200000;
        unsigned int n = 2;
        double split_weight = 1.0;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        RandomNumberGenerator rng(54321);
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::pair< unsigned int, std::vector<unsigned int> > ret = rng.random_set_partition(n, split_weight);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: ret.second) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            unsigned int nevents = ret.first;
            REQUIRE(max_index + 1 == nevents);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
            if (nevent_counts.count(nevents) < 1) {
                nevent_counts[nevents] = 1;
            }
            else {
                ++nevent_counts[nevents];
            }
        }

        REQUIRE(model_counts.at("00") == nevent_counts.at(1));
        REQUIRE(model_counts.at("01") == nevent_counts.at(2));
        unsigned int tally = 0;
        for (auto const & kv: model_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);
        tally = 0;
        for (auto const & kv: nevent_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: model_counts) {
            REQUIRE((kv.second / (double)nsamples) == Approx(1.0/model_counts.size()).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing random_set_partition (2, 0.5)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition(2, 0.5)") {
        unsigned int nsamples = 100000;
        unsigned int n = 2;
        double split_weight = 0.5;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        RandomNumberGenerator rng(4321);
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::pair< unsigned int, std::vector<unsigned int> > ret = rng.random_set_partition(n, split_weight);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: ret.second) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            unsigned int nevents = ret.first;
            REQUIRE(max_index + 1 == nevents);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
            if (nevent_counts.count(nevents) < 1) {
                nevent_counts[nevents] = 1;
            }
            else {
                ++nevent_counts[nevents];
            }
        }

        REQUIRE(model_counts.at("00") == nevent_counts.at(1));
        REQUIRE(model_counts.at("01") == nevent_counts.at(2));
        unsigned int tally = 0;
        for (auto const & kv: model_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);
        tally = 0;
        for (auto const & kv: nevent_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        REQUIRE(model_counts.at("00")/(double)nsamples == Approx(2.0/3.0).epsilon(0.001));
        REQUIRE(model_counts.at("01")/(double)nsamples == Approx(1.0/3.0).epsilon(0.001));
    }
}


TEST_CASE("Testing random_set_partition (3, 1.0)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition(3, 1.0)") {
        unsigned int nsamples = 750000;
        unsigned int n = 3;
        double split_weight = 1.0;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        RandomNumberGenerator rng(654321);
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::pair< unsigned int, std::vector<unsigned int> > ret = rng.random_set_partition(n, split_weight);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: ret.second) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            unsigned int nevents = ret.first;
            REQUIRE(max_index + 1 == nevents);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
            if (nevent_counts.count(nevents) < 1) {
                nevent_counts[nevents] = 1;
            }
            else {
                ++nevent_counts[nevents];
            }
        }

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));
        unsigned int tally = 0;
        for (auto const & kv: model_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);
        tally = 0;
        for (auto const & kv: nevent_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: model_counts) {
            REQUIRE((kv.second / (double)nsamples) == Approx(1.0/model_counts.size()).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing random_set_partition (3, 3.0)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition(3, 3.0)") {
        unsigned int nsamples = 500000;
        unsigned int n = 3;
        double split_weight = 3.0;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        RandomNumberGenerator rng(1234);
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::pair< unsigned int, std::vector<unsigned int> > ret = rng.random_set_partition(n, split_weight);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: ret.second) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            unsigned int nevents = ret.first;
            REQUIRE(max_index + 1 == nevents);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
            if (nevent_counts.count(nevents) < 1) {
                nevent_counts[nevents] = 1;
            }
            else {
                ++nevent_counts[nevents];
            }
        }

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));
        unsigned int tally = 0;
        for (auto const & kv: model_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);
        tally = 0;
        for (auto const & kv: nevent_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        REQUIRE(model_counts.at("000")/(double)nsamples == Approx(1.0/19.0).epsilon(0.001));
        REQUIRE(model_counts.at("001")/(double)nsamples == Approx(3.0/19.0).epsilon(0.001));
        REQUIRE(model_counts.at("010")/(double)nsamples == Approx(3.0/19.0).epsilon(0.001));
        REQUIRE(model_counts.at("011")/(double)nsamples == Approx(3.0/19.0).epsilon(0.001));
        REQUIRE(model_counts.at("012")/(double)nsamples == Approx(9.0/19.0).epsilon(0.001));
    }
}

TEST_CASE("Testing random_set_partition (3, 1/3.0)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition(3, 1/3.0)") {
        unsigned int nsamples = 500000;
        unsigned int n = 3;
        double split_weight = 1.0/3.0;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        RandomNumberGenerator rng(1234);
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::pair< unsigned int, std::vector<unsigned int> > ret = rng.random_set_partition(n, split_weight);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: ret.second) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            unsigned int nevents = ret.first;
            REQUIRE(max_index + 1 == nevents);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
            if (nevent_counts.count(nevents) < 1) {
                nevent_counts[nevents] = 1;
            }
            else {
                ++nevent_counts[nevents];
            }
        }

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));
        unsigned int tally = 0;
        for (auto const & kv: model_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);
        tally = 0;
        for (auto const & kv: nevent_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        REQUIRE(model_counts.at("012")/(double)nsamples == Approx(1.0/19.0).epsilon(0.001));
        REQUIRE(model_counts.at("001")/(double)nsamples == Approx(3.0/19.0).epsilon(0.001));
        REQUIRE(model_counts.at("010")/(double)nsamples == Approx(3.0/19.0).epsilon(0.001));
        REQUIRE(model_counts.at("011")/(double)nsamples == Approx(3.0/19.0).epsilon(0.001));
        REQUIRE(model_counts.at("000")/(double)nsamples == Approx(9.0/19.0).epsilon(0.001));
    }
}

TEST_CASE("Testing random_set_partition (4, 1.0)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition(4, 1.0)") {
        unsigned int nsamples = 1000000;
        unsigned int n = 4;
        double split_weight = 1.0;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        RandomNumberGenerator rng(333333);
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::pair< unsigned int, std::vector<unsigned int> > ret = rng.random_set_partition(n, split_weight);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: ret.second) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            unsigned int nevents = ret.first;
            REQUIRE(max_index + 1 == nevents);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
            if (nevent_counts.count(nevents) < 1) {
                nevent_counts[nevents] = 1;
            }
            else {
                ++nevent_counts[nevents];
            }
        }

        REQUIRE(model_counts.at("0000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("0123") == nevent_counts.at(4));
        unsigned int tally = 0;
        for (auto const & kv: model_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);
        tally = 0;
        for (auto const & kv: nevent_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: model_counts) {
            REQUIRE((kv.second / (double)nsamples) == Approx(1.0/model_counts.size()).epsilon(0.005));
        }
    }
}

TEST_CASE("Testing random_set_partition (4, 3.0)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition(4, 3.0)") {
        unsigned int nsamples = 1000000;
        unsigned int n = 4;
        double split_weight = 3.0;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        RandomNumberGenerator rng(333333);
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::pair< unsigned int, std::vector<unsigned int> > ret = rng.random_set_partition(n, split_weight);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: ret.second) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            unsigned int nevents = ret.first;
            REQUIRE(max_index + 1 == nevents);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
            if (nevent_counts.count(nevents) < 1) {
                nevent_counts[nevents] = 1;
            }
            else {
                ++nevent_counts[nevents];
            }
        }

        REQUIRE(model_counts.at("0000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("0123") == nevent_counts.at(4));
        unsigned int tally = 0;
        for (auto const & kv: model_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);
        tally = 0;
        for (auto const & kv: nevent_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }

        // Below are the expected results. 'k' is the number of categories,
        // 'S(n,k)' is the number of possible partitions of n elements into k
        // categories (Stirling number of second kind), 'w' is the relative split
        // weight for any one of those possible partitions, 'total_w' is the
        // overall weight of the k class (w * S(n,k)), and 'prob' is the
        // probability of any one possible partition in the k class.
        //
        // k    S(n,k)  w       total_w     prob
        // -------------------------------------
        // 1    1       1       1           1/103
        // 2    7       3       21          3/103
        // 3    6       9       54          9/103
        // 4    1       27      27          27/103
        //                      103

        REQUIRE((model_counts.at("0000") / (double)nsamples) == Approx(1.0/103.0).epsilon(0.001));
        REQUIRE((model_counts.at("0010") / (double)nsamples) == Approx(3.0/103.0).epsilon(0.001));
        REQUIRE((model_counts.at("0122") / (double)nsamples) == Approx(9.0/103.0).epsilon(0.001));
        REQUIRE((model_counts.at("0123") / (double)nsamples) == Approx(27.0/103.0).epsilon(0.001));
        REQUIRE((nevent_counts.at(2) / (double)nsamples) == Approx(21.0/103.0).epsilon(0.001));
        REQUIRE((nevent_counts.at(3) / (double)nsamples) == Approx(54.0/103.0).epsilon(0.001));
    }
}

TEST_CASE("Testing random_set_partition (4, 1/3.0)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition(4, 1/3.0)") {
        unsigned int nsamples = 1000000;
        unsigned int n = 4;
        double split_weight = 1.0/3.0;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        RandomNumberGenerator rng(333333);
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::pair< unsigned int, std::vector<unsigned int> > ret = rng.random_set_partition(n, split_weight);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: ret.second) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            unsigned int nevents = ret.first;
            REQUIRE(max_index + 1 == nevents);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
            if (nevent_counts.count(nevents) < 1) {
                nevent_counts[nevents] = 1;
            }
            else {
                ++nevent_counts[nevents];
            }
        }

        REQUIRE(model_counts.at("0000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("0123") == nevent_counts.at(4));
        unsigned int tally = 0;
        for (auto const & kv: model_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);
        tally = 0;
        for (auto const & kv: nevent_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }

        // Below are the expected results. 'k' is the number of categories,
        // 'S(n,k)' is the number of possible partitions of n elements into k
        // categories (Stirling number of second kind), 'w' is the relative split
        // weight for any one of those possible partitions, 'total_w' is the
        // overall weight of the k class (w * S(n,k)), and 'prob' is the
        // probability of any one possible partition in the k class.
        //
        // k    S(n,k)  w       total_w     prob
        // -------------------------------------
        // 1    1       27      27          27/109
        // 2    7       9       63          9/109
        // 3    6       3       18          3/109
        // 4    1       1       1           1/109
        //                      109

        REQUIRE((model_counts.at("0000") / (double)nsamples) == Approx(27.0/109.0).epsilon(0.001));
        REQUIRE((model_counts.at("0010") / (double)nsamples) == Approx(9.0/109.0).epsilon(0.001));
        REQUIRE((model_counts.at("0122") / (double)nsamples) == Approx(3.0/109.0).epsilon(0.001));
        REQUIRE((model_counts.at("0123") / (double)nsamples) == Approx(1.0/109.0).epsilon(0.001));
        REQUIRE((nevent_counts.at(2) / (double)nsamples) == Approx(63.0/109.0).epsilon(0.001));
        REQUIRE((nevent_counts.at(3) / (double)nsamples) == Approx(18.0/109.0).epsilon(0.001));
    }
}

TEST_CASE("Testing random_set_partition (5, 1.0)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition(5, 1.0)") {
        unsigned int nsamples = 1000000;
        unsigned int n = 5;
        double split_weight = 1.0;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        RandomNumberGenerator rng(333333);
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::pair< unsigned int, std::vector<unsigned int> > ret = rng.random_set_partition(n, split_weight);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: ret.second) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            unsigned int nevents = ret.first;
            REQUIRE(max_index + 1 == nevents);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
            if (nevent_counts.count(nevents) < 1) {
                nevent_counts[nevents] = 1;
            }
            else {
                ++nevent_counts[nevents];
            }
        }

        REQUIRE(model_counts.at("00000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("01234") == nevent_counts.at(5));
        unsigned int tally = 0;
        for (auto const & kv: model_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);
        tally = 0;
        for (auto const & kv: nevent_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: model_counts) {
            REQUIRE((kv.second / (double)nsamples) == Approx(1.0/model_counts.size()).epsilon(0.005));
        }
    }
}

TEST_CASE("Testing random_set_partition (6, 1.0)",
        "[RandomNumberGenerator]") {
    SECTION("Testing random_set_partition(6, 1.0)") {
        unsigned int nsamples = 1000000;
        unsigned int n = 6;
        double split_weight = 1.0;

        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts;

        RandomNumberGenerator rng(333333);
        for (unsigned int i = 0; i < nsamples; ++i) {
            std::pair< unsigned int, std::vector<unsigned int> > ret = rng.random_set_partition(n, split_weight);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: ret.second) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            unsigned int nevents = ret.first;
            REQUIRE(max_index + 1 == nevents);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
            if (nevent_counts.count(nevents) < 1) {
                nevent_counts[nevents] = 1;
            }
            else {
                ++nevent_counts[nevents];
            }
        }

        REQUIRE(model_counts.at("000000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012345") == nevent_counts.at(6));
        unsigned int tally = 0;
        for (auto const & kv: model_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);
        tally = 0;
        for (auto const & kv: nevent_counts) {
            tally += kv.second;
        }
        REQUIRE(tally == nsamples);

        for (auto const & kv: model_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: nevent_counts) {
            std::cout << kv.first << ": " << kv.second / (double)nsamples << "\n";
        }
        for (auto const & kv: model_counts) {
            REQUIRE((kv.second / (double)nsamples) == Approx(1.0/model_counts.size()).epsilon(0.005));
        }
    }
}

TEST_CASE("Testing dirichlet(1, 1)", "[RandomNumberGenerator]") {

    SECTION("Testing dirichlet({1, 1})") {
        std::vector<double> parameters {1.0, 1.0};
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

        std::vector< SampleSummarizer<double> > summaries(k);
        std::vector<double> x(k);
        unsigned int nsamples = 100000;
        RandomNumberGenerator rng(123);
        for (unsigned int i = 0; i < nsamples; ++i) {
            x = rng.dirichlet(parameters);
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

TEST_CASE("Testing dirichlet(3, 3)", "[RandomNumberGenerator]") {

    SECTION("Testing dirichlet({3, 3})") {
        std::vector<double> parameters {3.0, 3.0};
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

        std::vector< SampleSummarizer<double> > summaries(k);
        std::vector<double> x(k);
        unsigned int nsamples = 100000;
        RandomNumberGenerator rng(123);
        for (unsigned int i = 0; i < nsamples; ++i) {
            x = rng.dirichlet(parameters);
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

TEST_CASE("Testing dirichlet(0.3, 0.3)", "[RandomNumberGenerator]") {

    SECTION("Testing dirichlet({0.3, 0.3})") {
        std::vector<double> parameters {0.3, 0.3};
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

        std::vector< SampleSummarizer<double> > summaries(k);
        std::vector<double> x(k);
        unsigned int nsamples = 100000;
        RandomNumberGenerator rng(123);
        for (unsigned int i = 0; i < nsamples; ++i) {
            x = rng.dirichlet(parameters);
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

TEST_CASE("Testing dirichlet(6, 1)", "[RandomNumberGenerator]") {

    SECTION("Testing dirichlet({6, 1})") {
        std::vector<double> parameters {6.0, 1.0};
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

        std::vector< SampleSummarizer<double> > summaries(k);
        std::vector<double> x(k);
        unsigned int nsamples = 100000;
        RandomNumberGenerator rng(123);
        for (unsigned int i = 0; i < nsamples; ++i) {
            x = rng.dirichlet(parameters);
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

TEST_CASE("Testing dirichlet(1, 4)", "[RandomNumberGenerator]") {

    SECTION("Testing dirichlet({1, 4})") {
        std::vector<double> parameters {1.0, 4.0};
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

        std::vector< SampleSummarizer<double> > summaries(k);
        std::vector<double> x(k);
        unsigned int nsamples = 100000;
        RandomNumberGenerator rng(123);
        for (unsigned int i = 0; i < nsamples; ++i) {
            x = rng.dirichlet(parameters);
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

TEST_CASE("Testing dirichlet(1, 1, 1, 1, 1)", "[RandomNumberGenerator]") {

    SECTION("Testing dirichlet({1, 1, 1, 1, 1})") {
        std::vector<double> parameters {1.0, 1.0, 1.0, 1.0, 1.0};
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

        std::vector< SampleSummarizer<double> > summaries(k);
        std::vector<double> x(k);
        unsigned int nsamples = 100000;
        RandomNumberGenerator rng(123);
        for (unsigned int i = 0; i < nsamples; ++i) {
            x = rng.dirichlet(parameters);
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

TEST_CASE("Testing dirichlet(1, 0.5, 1, 10, 5)", "[RandomNumberGenerator]") {

    SECTION("Testing dirichlet({1, 0.5, 1, 10, 5})") {
        std::vector<double> parameters {1.0, 0.5, 1.0, 10.0, 5.0};
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

        std::vector< SampleSummarizer<double> > summaries(k);
        std::vector<double> x(k);
        unsigned int nsamples = 100000;
        RandomNumberGenerator rng(123);
        for (unsigned int i = 0; i < nsamples; ++i) {
            x = rng.dirichlet(parameters);
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

TEST_CASE("Testing pitman_yor_process(6, 1.7, 0.0)", "[RandomNumberGenerator]") {

    SECTION("Testing pitman_yor_process(6, 1.7, 0.0)") {
        unsigned int nsamples = 100000;
        unsigned int n = 6;
        double a = 1.7;
        double d = 0.0;

        std::vector<unsigned int> elements (6, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
                {4, 0},
                {5, 0},
                {6, 0}
        };

        RandomNumberGenerator rng(123);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.pitman_yor_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012345") == nevent_counts.at(6));

        for (std::string m: {"000000", "012345", "012344", "012340", "001234", "012314"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_dpp_log_prior_probability(m, a))).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing weighted_discount_process(6, 1.7, 0.0)", "[RandomNumberGenerator]") {

    SECTION("Testing weighted_discount_process(6, 1.7, 0.0)") {
        unsigned int nsamples = 100000;
        unsigned int n = 6;
        double a = 1.7;
        double d = 0.0;

        std::vector<unsigned int> elements (6, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
                {4, 0},
                {5, 0},
                {6, 0}
        };

        RandomNumberGenerator rng(123);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.weighted_discount_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012345") == nevent_counts.at(6));

        for (std::string m: {"000000", "012345", "012344", "012340", "001234", "012314"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_wdp_log_prior_probability(m, a, d))).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing pitman_yor_process(6, 1.7, 0.1)", "[RandomNumberGenerator]") {

    SECTION("Testing pitman_yor_process(6, 1.7, 0.1)") {
        unsigned int nsamples = 100000;
        unsigned int n = 6;
        double a = 1.7;
        double d = 0.1;

        std::vector<unsigned int> elements (6, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
                {4, 0},
                {5, 0},
                {6, 0}
        };

        RandomNumberGenerator rng(123);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.pitman_yor_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012345") == nevent_counts.at(6));

        for (std::string m: {"000000", "012345", "012344", "012340", "001234", "012314"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_pyp_log_prior_probability(m, a, d))).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing weighted_discount_process(6, 1.7, 0.1)", "[RandomNumberGenerator]") {

    SECTION("Testing weighted_discount_process(6, 1.7, 0.1)") {
        unsigned int nsamples = 100000;
        unsigned int n = 6;
        double a = 1.7;
        double d = 0.1;

        std::vector<unsigned int> elements (6, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
                {4, 0},
                {5, 0},
                {6, 0}
        };

        RandomNumberGenerator rng(123);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.weighted_discount_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012345") == nevent_counts.at(6));

        for (std::string m: {"000000", "012345", "012344", "012340", "001234", "012314"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_wdp_log_prior_probability(m, a, d))).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing pitman_yor_process(6, 1.7, 0.5)", "[RandomNumberGenerator]") {

    SECTION("Testing pitman_yor_process(6, 1.7, 0.5)") {
        unsigned int nsamples = 200000;
        unsigned int n = 6;
        double a = 1.7;
        double d = 0.5;

        std::vector<unsigned int> elements (6, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
                {4, 0},
                {5, 0},
                {6, 0}
        };

        RandomNumberGenerator rng(12345);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.pitman_yor_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012345") == nevent_counts.at(6));

        for (std::string m: {"000000", "012345", "012344", "012340", "001234", "012314"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_pyp_log_prior_probability(m, a, d))).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing weighted_discount_process(6, 1.7, 0.5)", "[RandomNumberGenerator]") {

    SECTION("Testing weighted_discount_process(6, 1.7, 0.5)") {
        unsigned int nsamples = 200000;
        unsigned int n = 6;
        double a = 1.7;
        double d = 0.5;

        std::vector<unsigned int> elements (6, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
                {4, 0},
                {5, 0},
                {6, 0}
        };

        RandomNumberGenerator rng(12345);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.weighted_discount_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012345") == nevent_counts.at(6));

        for (std::string m: {"000000", "012345", "012344", "012340", "001234", "012314"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_wdp_log_prior_probability(m, a, d))).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing pitman_yor_process(6, 1.7, 0.9)", "[RandomNumberGenerator]") {

    SECTION("Testing pitman_yor_process(6, 1.7, 0.9)") {
        unsigned int nsamples = 100000;
        unsigned int n = 6;
        double a = 1.7;
        double d = 0.9;

        std::vector<unsigned int> elements (6, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
                {4, 0},
                {5, 0},
                {6, 0}
        };

        RandomNumberGenerator rng(123);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.pitman_yor_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012345") == nevent_counts.at(6));

        for (std::string m: {"000000", "012345", "012344", "012340", "001234", "012314"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_pyp_log_prior_probability(m, a, d))).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing weighted_discount_process(6, 1.7, 0.9)", "[RandomNumberGenerator]") {

    SECTION("Testing weighted_discount_process(6, 1.7, 0.9)") {
        unsigned int nsamples = 100000;
        unsigned int n = 6;
        double a = 1.7;
        double d = 0.9;

        std::vector<unsigned int> elements (6, 0);
        std::map<std::string, int> model_counts;
        model_counts["000000"] = 0;
        model_counts["012345"] = 0;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
                {4, 0},
                {5, 0},
                {6, 0}
        };

        RandomNumberGenerator rng(123);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.weighted_discount_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012345") == nevent_counts.at(6));

        for (std::string m: {"000000", "012345", "012344", "012340", "001234", "012314"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_wdp_log_prior_probability(m, a, d))).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing pitman_yor_process(3, 0.6, 0.1)", "[RandomNumberGenerator]") {

    SECTION("Testing pitman_yor_process(3, 0.6, 0.1)") {
        unsigned int nsamples = 200000;
        unsigned int n = 3;
        double a = 0.6;
        double d = 0.1;

        std::vector<unsigned int> elements (3, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
        };

        RandomNumberGenerator rng(342254);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.pitman_yor_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_pyp_log_prior_probability(m, a, d))).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing weighted_discount_process(3, 0.6, 0.1)", "[RandomNumberGenerator]") {

    SECTION("Testing weighted_discount_process(3, 0.6, 0.1)") {
        unsigned int nsamples = 200000;
        unsigned int n = 3;
        double a = 0.6;
        double d = 0.1;

        std::vector<unsigned int> elements (3, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
        };

        RandomNumberGenerator rng(342254);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.weighted_discount_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_wdp_log_prior_probability(m, a, d))).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing pitman_yor_process(3, 0.6, 0.5)", "[RandomNumberGenerator]") {

    SECTION("Testing pitman_yor_process(3, 0.6, 0.5)") {
        unsigned int nsamples = 200000;
        unsigned int n = 3;
        double a = 0.6;
        double d = 0.5;

        std::vector<unsigned int> elements (3, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
        };

        RandomNumberGenerator rng(342254);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.pitman_yor_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_pyp_log_prior_probability(m, a, d))).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing weighted_discount_process(3, 0.6, 0.5)", "[RandomNumberGenerator]") {

    SECTION("Testing weighted_discount_process(3, 0.6, 0.5)") {
        unsigned int nsamples = 200000;
        unsigned int n = 3;
        double a = 0.6;
        double d = 0.5;

        std::vector<unsigned int> elements (3, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
        };

        RandomNumberGenerator rng(342254);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.weighted_discount_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_wdp_log_prior_probability(m, a, d))).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing pitman_yor_process(3, 0.6, 0.9)", "[RandomNumberGenerator]") {

    SECTION("Testing pitman_yor_process(3, 0.6, 0.9)") {
        unsigned int nsamples = 200000;
        unsigned int n = 3;
        double a = 0.6;
        double d = 0.9;

        std::vector<unsigned int> elements (3, 0);
        std::map<std::string, int> model_counts;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
        };

        RandomNumberGenerator rng(342254);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.pitman_yor_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_pyp_log_prior_probability(m, a, d))).epsilon(0.002));
        }
    }
}

TEST_CASE("Testing weighted_discount_process(3, 0.6, 0.9)", "[RandomNumberGenerator]") {

    SECTION("Testing weighted_discount_process(3, 0.6, 0.9)") {
        unsigned int nsamples = 200000;
        unsigned int n = 3;
        double a = 0.6;
        double d = 0.9;

        std::vector<unsigned int> elements (3, 0);
        std::map<std::string, int> model_counts;
        model_counts["000"] = 0;
        model_counts["012"] = 0;
        std::map<int, int> nevent_counts = {
                {1, 0},
                {2, 0},
                {3, 0},
        };

        RandomNumberGenerator rng(342254);
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int ncats = rng.weighted_discount_process(elements, a, d);
            std::ostringstream stream;
            unsigned int max_index = 0;
            for (auto e: elements) {
                stream << e;
                if (e > max_index) {
                    max_index = e;
                }
            }
            REQUIRE(max_index + 1 == ncats);
            ++nevent_counts.at(ncats);
            std::string model_str = stream.str();
            if (model_counts.count(model_str) < 1) {
                model_counts[model_str] = 1;
            }
            else {
                ++model_counts[model_str];
            }
        }

        int total = 0;
        for (auto const &kv: model_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);
        total = 0;
        for (auto const &kv: nevent_counts) {
            total += kv.second;
        }
        REQUIRE(total == nsamples);

        REQUIRE(model_counts.at("000") == nevent_counts.at(1));
        REQUIRE(model_counts.at("012") == nevent_counts.at(3));
        REQUIRE((model_counts.at("001") + model_counts.at("010") + model_counts.at("011")) == nevent_counts.at(2));

        for (std::string m: {"000", "001", "010", "011", "012"}) {
            REQUIRE((model_counts.at(m) / (double)nsamples) == Approx(std::exp(
                    get_wdp_log_prior_probability(m, a, d))).epsilon(0.002));
        }
    }
}
