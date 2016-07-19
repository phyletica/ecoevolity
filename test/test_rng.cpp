#include "catch.hpp"
#include "ecoevolity/rng.hpp"

#include <limits>

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
