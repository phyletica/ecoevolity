#include "catch.hpp"
#include "ecoevolity/rng.hpp"
#include "ecoevolity/math_util.hpp"
#include "ecoevolity/stats_util.hpp"

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=1, split_weight=1.0",
        "[math_util]") {
    SECTION("Testing n=1, weight=1.0") {
        unsigned int n = 1;
        double w = 1.0;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(1.0)));
    }
}

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=1, split_weight=12.9",
        "[math_util]") {
    SECTION("Testing n=1, weight=12.9") {
        unsigned int n = 1;
        double w = 12.9;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(1.0)));
    }
}

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=2, split_weight=1.0",
        "[math_util]") {
    SECTION("Testing n=2, weight=1.0") {
        unsigned int n = 2;
        double w = 1.0;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(0.5)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 2, w) == Approx(std::log(0.5)));
    }
}

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=2, split_weight=2.0",
        "[math_util]") {
    SECTION("Testing n=2, weight=2.0") {
        unsigned int n = 2;
        double w = 2.0;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(1.0/3.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 2, w) == Approx(std::log(2.0/3.0)));
    }
}

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=2, split_weight=2.5",
        "[math_util]") {
    SECTION("Testing n=2, weight=2.5") {
        unsigned int n = 2;
        double w = 2.5;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(2.0/7.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 2, w) == Approx(std::log(5.0/7.0)));
    }
}

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=2, split_weight=0.5",
        "[math_util]") {
    SECTION("Testing n=2, weight=0.5") {
        unsigned int n = 2;
        double w = 0.5;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(2.0/3.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 2, w) == Approx(std::log(1.0/3.0)));
    }
}

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=2, split_weight=1/2.5",
        "[math_util]") {
    SECTION("Testing n=2, weight=1/2.5") {
        unsigned int n = 2;
        double w = 1.0/2.5;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(5.0/7.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 2, w) == Approx(std::log(2.0/7.0)));
    }
}

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=3, split_weight=1.0",
        "[math_util]") {
    SECTION("Testing n=3, weight=1.0") {
        unsigned int n = 3;
        double w = 1.0;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(1.0/5.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 2, w) == Approx(std::log(1.0/5.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 3, w) == Approx(std::log(1.0/5.0)));
    }
}

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=3, split_weight=2.0",
        "[math_util]") {
    SECTION("Testing n=3, weight=2.0") {
        unsigned int n = 3;
        double w = 2.0;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(1.0/11.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 2, w) == Approx(std::log(2.0/11.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 3, w) == Approx(std::log(4.0/11.0)));
    }
}

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=3, split_weight=1/2.0",
        "[math_util]") {
    SECTION("Testing n=3, weight=1/2.0") {
        unsigned int n = 3;
        double w = 1.0/2.0;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(4.0/11.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 2, w) == Approx(std::log(2.0/11.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 3, w) == Approx(std::log(1.0/11.0)));
    }
}

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=3, split_weight=3.0",
        "[math_util]") {
    SECTION("Testing n=3, weight=3.0") {
        unsigned int n = 3;
        double w = 3.0;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(1.0/19.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 2, w) == Approx(std::log(3.0/19.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 3, w) == Approx(std::log(9.0/19.0)));
    }
}

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=3, split_weight=1/3.0",
        "[math_util]") {
    SECTION("Testing n=3, weight=1/3.0") {
        unsigned int n = 3;
        double w = 1.0/3.0;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(9.0/19.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 2, w) == Approx(std::log(3.0/19.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 3, w) == Approx(std::log(1.0/19.0)));
    }
}

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=4, split_weight=1.0",
        "[math_util]") {
    SECTION("Testing n=4, weight=1.0") {
        unsigned int n = 4;
        double w = 1.0;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(1.0/15.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 2, w) == Approx(std::log(1.0/15.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 3, w) == Approx(std::log(1.0/15.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 4, w) == Approx(std::log(1.0/15.0)));
    }
}

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=4, split_weight=2.0",
        "[math_util]") {
    // Below are the expected results. 'k' is the number of categories,
    // 'S(n,k)' is the number of possible partitions of n elements into k
    // categories (Stirling number of second kind), 'w' is the relative split
    // weight for any one of those possible partitions, 'total_w' is the
    // overall weight of the k class (w * S(n,k)), and 'prob' is the
    // probability of any one possible partition in the k class.
    //
    // k    S(n,k)  w       total_w     prob
    // -------------------------------------
    // 1    1       1       1           1/47
    // 2    7       2       14          2/47
    // 3    6       4       24          4/47
    // 4    1       8       8           8/47
    //                      47
    SECTION("Testing n=4, weight=2.0") {
        unsigned int n = 4;
        double w = 2.0;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(1.0/47.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 2, w) == Approx(std::log(2.0/47.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 3, w) == Approx(std::log(4.0/47.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 4, w) == Approx(std::log(8.0/47.0)));
    }
}

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=4, split_weight=1/2.0",
        "[math_util]") {
    // Below are the expected results. 'k' is the number of categories,
    // 'S(n,k)' is the number of possible partitions of n elements into k
    // categories (Stirling number of second kind), 'w' is the relative split
    // weight for any one of those possible partitions, 'total_w' is the
    // overall weight of the k class (w * S(n,k)), and 'prob' is the
    // probability of any one possible partition in the k class.
    //
    // k    S(n,k)  w       total_w     prob
    // -------------------------------------
    // 1    1       8       8           8/49
    // 2    7       4       28          4/49
    // 3    6       2       12          2/49
    // 4    1       1       1           1/49
    //                      49

    SECTION("Testing n=4, weight=1/2.0") {
        unsigned int n = 4;
        double w = 1.0/2.0;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(8.0/49.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 2, w) == Approx(std::log(4.0/49.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 3, w) == Approx(std::log(2.0/49.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 4, w) == Approx(std::log(1.0/49.0)));
    }
}

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=4, split_weight=3.0",
        "[math_util]") {
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

    SECTION("Testing n=4, weight=3.0") {
        unsigned int n = 4;
        double w = 3.0;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(1.0/103.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 2, w) == Approx(std::log(3.0/103.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 3, w) == Approx(std::log(9.0/103.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 4, w) == Approx(std::log(27.0/103.0)));
    }
}

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=4, split_weight=1/3.0",
        "[math_util]") {
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

    SECTION("Testing n=4, weight=1/3.0") {
        unsigned int n = 4;
        double w = 1.0/3.0;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(27.0/109.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 2, w) == Approx(std::log(9.0/109.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 3, w) == Approx(std::log(3.0/109.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 4, w) == Approx(std::log(1.0/109.0)));
    }
}

TEST_CASE("Testing get_uniform_model_log_prior_probability with n=5, split_weight=3.0",
        "[math_util]") {
    // Below are the expected results. 'k' is the number of categories,
    // 'S(n,k)' is the number of possible partitions of n elements into k
    // categories (Stirling number of second kind), 'w' is the relative split
    // weight for any one of those possible partitions, 'total_w' is the
    // overall weight of the k class (w * S(n,k)), and 'prob' is the
    // probability of any one possible partition in the k class.
    //
    // k    S(n,k)  w       total_w     prob
    // -------------------------------------
    // 1    1       1       1           1/622
    // 2    15      3       45          3/622
    // 3    25      9       225         9/622
    // 4    10      27      270         27/622
    // 5    1       81      81          81/622
    //                      622

    SECTION("Testing n=5, weight=3.0") {
        unsigned int n = 5;
        double w = 3.0;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(1.0/622.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 2, w) == Approx(std::log(3.0/622.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 3, w) == Approx(std::log(9.0/622.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 4, w) == Approx(std::log(27.0/622.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 5, w) == Approx(std::log(81.0/622.0)));
    }
}
TEST_CASE("Testing get_uniform_model_log_prior_probability with n=5, split_weight=1/3.0",
        "[math_util]") {
    // Below are the expected results. 'k' is the number of categories,
    // 'S(n,k)' is the number of possible partitions of n elements into k
    // categories (Stirling number of second kind), 'w' is the relative split
    // weight for any one of those possible partitions, 'total_w' is the
    // overall weight of the k class (w * S(n,k)), and 'prob' is the
    // probability of any one possible partition in the k class.
    //
    // k    S(n,k)  w       total_w     prob
    // -------------------------------------
    // 1    1       81      81          81/742
    // 2    15      27      405         27/742
    // 3    25      9       225         9/742
    // 4    10      3       30          3/742
    // 5    1       1       1           1/742
    //                      742

    SECTION("Testing n=5, weight=1/3.0") {
        unsigned int n = 5;
        double w = 1.0/3.0;
        REQUIRE(get_uniform_model_log_prior_probability(n, 1, w) == Approx(std::log(81.0/742.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 2, w) == Approx(std::log(27.0/742.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 3, w) == Approx(std::log(9.0/742.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 4, w) == Approx(std::log(3.0/742.0)));
        REQUIRE(get_uniform_model_log_prior_probability(n, 5, w) == Approx(std::log(1.0/742.0)));
    }
}

TEST_CASE("Testing normalize_log_likelihoods", "[math_util]") {

    SECTION("Testing round trip") {
        std::vector<double> v = {0.2, 0.2, 0.2, 0.2, 0.2};
        std::vector<double> v2(v.size());
        for (unsigned int i = 0; i < v.size(); ++i) {
            v2.at(i) = std::log(v.at(i));
        }
        normalize_log_likelihoods(v2);
        for (unsigned int i = 0; i < v.size(); ++i) {
            REQUIRE(v.at(i) == Approx(v2.at(i)));
        }
    }

    SECTION("Testing uneven round trip") {
        std::vector<double> v = {0.1, 0.2, 0.3, 0.4};
        std::vector<double> v2(v.size());
        for (unsigned int i = 0; i < v.size(); ++i) {
            v2.at(i) = std::log(v.at(i));
        }
        normalize_log_likelihoods(v2);
        for (unsigned int i = 0; i < v.size(); ++i) {
            REQUIRE(v.at(i) == Approx(v2.at(i)));
        }
    }

    SECTION("Testing uneven unnormalized round trip") {
        std::vector<double> v = {0.00001, 0.00002, 0.00003, 0.00004};
        std::vector<double> v2(v.size());
        for (unsigned int i = 0; i < v.size(); ++i) {
            v2.at(i) = std::log(v.at(i));
        }
        normalize_log_likelihoods(v2);
        for (unsigned int i = 0; i < v.size(); ++i) {
            REQUIRE(v.at(i) * 10000.0 == Approx(v2.at(i)));
        }
    }
}

TEST_CASE("Testing normalize_log_likelihoods of copy", "[math_util]") {
    SECTION("Testing copy normalization") {
        std::vector<double> v = {-23.45, -12.24, -13.46, -45.12};
        std::vector<double> v2(v);
        normalize_log_likelihoods(v2);
        REQUIRE(v.at(0) == -23.45);
        REQUIRE(v.at(1) == -12.24);
        REQUIRE(v.at(2) == -13.46);
        REQUIRE(v.at(3) == -45.12);
    }
}

TEST_CASE("Testing get_dpp_expected_number_of_categories 0.218 7282", "[math_util]") {
    double e = get_dpp_expected_number_of_categories(0.218, 7282);
    REQUIRE(e == Approx(3.0).epsilon(0.001));
}

TEST_CASE("Testing get_dpp_expected_number_of_categories 0.449 7282", "[math_util]") {
    double e = get_dpp_expected_number_of_categories(0.449, 7282);
    REQUIRE(e == Approx(5.0).epsilon(0.001));
}

TEST_CASE("Testing get_dpp_expected_number_of_categories 0.814 7282", "[math_util]") {
    double e = get_dpp_expected_number_of_categories(0.814, 7282);
    REQUIRE(e == Approx(8.0).epsilon(0.001));
}

TEST_CASE("Testing get_dpp_expected_number_of_categories 1.068 7282", "[math_util]") {
    double e = get_dpp_expected_number_of_categories(1.068, 7282);
    REQUIRE(e == Approx(10.0).epsilon(0.0001));
}

TEST_CASE("Testing get_dpp_expected_number_of_categories 2.9, 10", "[math_util]") {
    SECTION("Testing 2.9, 10") {
        unsigned int num_elements = 10;
        double a = 2.9;
        double e = get_dpp_expected_number_of_categories(a, num_elements);

        std::vector<unsigned int> elements (num_elements, 0);
        REQUIRE(elements.size() == num_elements);
        SampleSummarizer<double> ncats;
        RandomNumberGenerator rng(342254);

        unsigned int nsamples = 100000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int n = rng.dirichlet_process(elements, a);
            ncats.add_sample(float(n));
        }
        REQUIRE(ncats.mean() == Approx(e).epsilon(0.005));
    }
}

TEST_CASE("Testing get_dpp_concentration 3.0 7282", "[math_util]") {
    double e = get_dpp_concentration(3.0, 7282);
    REQUIRE(e == Approx(0.218).epsilon(0.001));
}
TEST_CASE("Testing get_dpp_concentration 5.0 7282", "[math_util]") {
    double e = get_dpp_concentration(5.0, 7282);
    REQUIRE(e == Approx(0.449).epsilon(0.001));
}
TEST_CASE("Testing get_dpp_concentration 8.0 7282", "[math_util]") {
    double e = get_dpp_concentration(8.0, 7282);
    REQUIRE(e == Approx(0.814).epsilon(0.0001));
}
TEST_CASE("Testing get_dpp_concentration 10.0 7282", "[math_util]") {
    double e = get_dpp_concentration(10.0, 7282);
    REQUIRE(e == Approx(1.068).epsilon(0.0001));
}

TEST_CASE("Testing get_dpp_gamma_scale", "[math_util]") {
    unsigned int n = 7282;
    double ncats = 10.0;
    double concentration = 1.068;
    std::vector<double> shapes {0.001, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0, 10000.0};
    for (auto shape : shapes) {
        double scale = get_dpp_gamma_scale(ncats, n, shape);
        REQUIRE(concentration == Approx(shape * scale).epsilon(0.001));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,0} alpha=1",
        "[math_util]") {
    double concentration = 1.0;
    double expected_prob = std::log(0.5);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "00";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_wdp_log_prior_probability with {0,0} alpha=1, discount=0.0",
        "[math_util]") {
    double concentration = 1.0;
    double discount = 0.0;
    double expected_prob = std::log(0.5);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0};
        REQUIRE(get_wdp_log_prior_probability<unsigned int>(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0};
        REQUIRE(get_wdp_log_prior_probability<int>(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0};
        REQUIRE(get_wdp_log_prior_probability<size_t>(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0'};
        REQUIRE(get_wdp_log_prior_probability<char>(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "00";
        REQUIRE(get_wdp_log_prior_probability(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_pyp_log_prior_probability with {0,0} alpha=1, discount=0.0",
        "[math_util]") {
    double concentration = 1.0;
    double discount = 0.0;
    double expected_prob = std::log(0.5);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0};
        REQUIRE(get_pyp_log_prior_probability<unsigned int>(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0};
        REQUIRE(get_pyp_log_prior_probability<int>(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0};
        REQUIRE(get_pyp_log_prior_probability<size_t>(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0'};
        REQUIRE(get_pyp_log_prior_probability<char>(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "00";
        REQUIRE(get_pyp_log_prior_probability(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,0} alpha=0.5",
        "[math_util]") {
    double concentration = 0.5;
    double expected_prob = std::log(2.0/3.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "00";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,1} alpha=0.5",
        "[math_util]") {
    double concentration = 0.5;
    double expected_prob = std::log(1.0/3.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "01";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_pyp_log_prior_probability with {0,1} alpha=0.5, discount=0.0",
        "[math_util]") {
    double concentration = 0.5;
    double discount = 0.0;
    double expected_prob = std::log(1.0/3.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1};
        REQUIRE(get_pyp_log_prior_probability<unsigned int>(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1};
        REQUIRE(get_pyp_log_prior_probability<int>(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1};
        REQUIRE(get_pyp_log_prior_probability<size_t>(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1'};
        REQUIRE(get_pyp_log_prior_probability<char>(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "01";
        REQUIRE(get_pyp_log_prior_probability(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_wdp_log_prior_probability with {0,1} alpha=0.5, discount=0.0",
        "[math_util]") {
    double concentration = 0.5;
    double discount = 0.0;
    double expected_prob = std::log(1.0/3.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1};
        REQUIRE(get_wdp_log_prior_probability<unsigned int>(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1};
        REQUIRE(get_wdp_log_prior_probability<int>(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1};
        REQUIRE(get_wdp_log_prior_probability<size_t>(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1'};
        REQUIRE(get_wdp_log_prior_probability<char>(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "01";
        REQUIRE(get_wdp_log_prior_probability(partition,
                    concentration,
                    discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,0} alpha=2.0",
        "[math_util]") {
    double concentration = 2.0;
    double expected_prob = std::log(1.0/3.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "00";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,1} alpha=2.0",
        "[math_util]") {
    double concentration = 2.0;
    double expected_prob = std::log(2.0/3.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "01";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,0,0} alpha=1.0",
        "[math_util]") {
    double concentration = 1.0;
    double expected_prob = std::log(2.0/6.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0, 0};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0, 0};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0, 0};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0', '0'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "000";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,0,1} alpha=1.0",
        "[math_util]") {
    double concentration = 1.0;
    double expected_prob = std::log(1.0/6.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0, 1};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0, 1};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0, 1};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0', '1'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "001";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,1,0} alpha=1.0",
        "[math_util]") {
    double concentration = 1.0;
    double expected_prob = std::log(1.0/6.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1, 0};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1, 0};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1, 0};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1', '0'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "010";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,1,1} alpha=1.0",
        "[math_util]") {
    double concentration = 1.0;
    double expected_prob = std::log(1.0/6.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1, 1};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1, 1};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1, 1};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1', '1'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "011";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,1,2} alpha=1.0",
        "[math_util]") {
    double concentration = 1.0;
    double expected_prob = std::log(1.0/6.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1, 2};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1, 2};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1, 2};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1', '2'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "012";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,0,0} alpha=3.0",
        "[math_util]") {
    double concentration = 3.0;
    double expected_prob = std::log(2.0/20.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0, 0};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0, 0};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0, 0};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0', '0'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "000";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,0,1} alpha=3.0",
        "[math_util]") {
    double concentration = 3.0;
    double expected_prob = std::log(3.0/20.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0, 1};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0, 1};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0, 1};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0', '1'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "001";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,1,0} alpha=3.0",
        "[math_util]") {
    double concentration = 3.0;
    double expected_prob = std::log(3.0/20.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1, 0};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1, 0};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1, 0};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1', '0'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "010";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,1,1} alpha=3.0",
        "[math_util]") {
    double concentration = 3.0;
    double expected_prob = std::log(3.0/20.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1, 1};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1, 1};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1, 1};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1', '1'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "011";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,1,2} alpha=3.0",
        "[math_util]") {
    double concentration = 3.0;
    double expected_prob = std::log(9.0/20.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1, 2};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1, 2};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1, 2};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1', '2'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "012";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,0,0,0,0} alpha=1.0",
        "[math_util]") {
    double concentration = 1.0;
    double expected_prob = std::log(24.0/120.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0, 0, 0, 0};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0, 0, 0, 0};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0, 0, 0, 0};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0', '0', '0', '0'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "00000";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,1,2,3,4} alpha=1.0",
        "[math_util]") {
    double concentration = 1.0;
    double expected_prob = std::log(1.0/120.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1, 2, 3, 4};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1, 2, 3, 4};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1, 2, 3, 4};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1', '2', '3', '4'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "01234";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability with {0,0,1,0,2} alpha=1.0",
        "[math_util]") {
    double concentration = 1.0;
    double expected_prob = std::log(2.0/120.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0, 1, 0, 2};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0, 1, 0, 2};
        REQUIRE(get_dpp_log_prior_probability<int>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0, 1, 0, 2};
        REQUIRE(get_dpp_log_prior_probability<size_t>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0', '1', '0', '2'};
        REQUIRE(get_dpp_log_prior_probability<char>(partition, concentration) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "00102";
        REQUIRE(get_dpp_log_prior_probability(partition, concentration) == Approx(expected_prob));
    }
}

TEST_CASE("Testing nonstandardized partitions for get_dpp_log_prior_probability",
        "[math_util]") {
    double concentration = 1.0;
    double expected_prob = std::log(2.0/120.0);
    SECTION("Testing standardized") {
    }
    SECTION("Testing non-standardized partitions") {
        std::vector<unsigned int> partition = {0, 0, 1, 0, 2};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) == Approx(expected_prob));
        std::vector<unsigned int> partition1 = {2, 2, 1, 2, 0};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) ==
                get_dpp_log_prior_probability<unsigned int>(partition1, concentration));
        std::vector<unsigned int> partition2 = {1, 1, 0, 1, 2};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) ==
                get_dpp_log_prior_probability<unsigned int>(partition2, concentration));
        std::vector<unsigned int> partition3 = {1, 1, 2, 1, 0};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) ==
                get_dpp_log_prior_probability<unsigned int>(partition3, concentration));
        std::vector<unsigned int> partition4 = {7, 7, 5, 7, 3};
        REQUIRE(get_dpp_log_prior_probability<unsigned int>(partition, concentration) ==
                get_dpp_log_prior_probability<unsigned int>(partition4, concentration));
    }
}

TEST_CASE("Testing stirling number functions", "[math_util]") {
    SECTION("Testing stirling2") {
        REQUIRE(stirling2(1, 1) == 1);
        REQUIRE(stirling2(2, 1) == 1);
        REQUIRE(stirling2(2, 2) == 1);
        REQUIRE(stirling2(3, 1) == 1);
        REQUIRE(stirling2(3, 2) == 3);
        REQUIRE(stirling2(3, 3) == 1);
        REQUIRE(stirling2(4, 1) == 1);
        REQUIRE(stirling2(4, 2) == 7);
        REQUIRE(stirling2(4, 3) == 6);
        REQUIRE(stirling2(4, 4) == 1);
        REQUIRE(stirling2(5, 1) == 1);
        REQUIRE(stirling2(5, 2) == 15);
        REQUIRE(stirling2(5, 3) == 25);
        REQUIRE(stirling2(5, 4) == 10);
        REQUIRE(stirling2(5, 5) == 1);
        REQUIRE(stirling2(9, 3) == 3025);
    }
    SECTION("Testing stirling2_float") {
        REQUIRE(stirling2_float(1, 1) == 1);
        REQUIRE(stirling2_float(2, 1) == 1);
        REQUIRE(stirling2_float(2, 2) == 1);
        REQUIRE(stirling2_float(3, 1) == 1);
        REQUIRE(stirling2_float(3, 2) == 3);
        REQUIRE(stirling2_float(3, 3) == 1);
        REQUIRE(stirling2_float(4, 1) == 1);
        REQUIRE(stirling2_float(4, 2) == 7);
        REQUIRE(stirling2_float(4, 3) == 6);
        REQUIRE(stirling2_float(4, 4) == 1);
        REQUIRE(stirling2_float(5, 1) == 1);
        REQUIRE(stirling2_float(5, 2) == 15);
        REQUIRE(stirling2_float(5, 3) == 25);
        REQUIRE(stirling2_float(5, 4) == 10);
        REQUIRE(stirling2_float(5, 5) == 1);
        REQUIRE(stirling2_float(9, 3) == 3025);
    }
}

TEST_CASE("Testing bell number functions", "[math_util]") {
    SECTION("Testing bell_number") {
        REQUIRE(bell_number(1) == 1);
        REQUIRE(bell_number(2) == 2);
        REQUIRE(bell_number(3) == 5);
        REQUIRE(bell_number(4) == 15);
        REQUIRE(bell_number(5) == 52);
        REQUIRE(bell_number(6) == 203);
        REQUIRE(bell_number(7) == 877);
        REQUIRE(bell_number(8) == 4140);
        REQUIRE(bell_number(9) == 21147);
        REQUIRE(bell_number(10) == 115975);
        REQUIRE(bell_number(15) == 1382958545);
    }
    SECTION("Testing bell_float") {
        REQUIRE(bell_float(1) == 1);
        REQUIRE(bell_float(2) == 2);
        REQUIRE(bell_float(3) == 5);
        REQUIRE(bell_float(4) == 15);
        REQUIRE(bell_float(5) == 52);
        REQUIRE(bell_float(6) == 203);
        REQUIRE(bell_float(7) == 877);
        REQUIRE(bell_float(8) == 4140);
        REQUIRE(bell_float(9) == 21147);
        REQUIRE(bell_float(10) == 115975);
        REQUIRE(bell_float(15) == 1382958545);
    }
}

TEST_CASE("Testing get_integer_partitions(2, 2)", "[math_util]") {
    SECTION("Testing get_integer_partitions") {
        std::vector< std::vector<unsigned int> > expected_partitions;
        expected_partitions.push_back({1, 1});
        REQUIRE(get_integer_partitions(2, 2) == expected_partitions);
    }
}

TEST_CASE("Testing get_integer_partitions(3, 2)", "[math_util]") {
    SECTION("Testing get_integer_partitions") {
        std::vector< std::vector<unsigned int> > expected_partitions;
        expected_partitions.push_back({2, 1});
        REQUIRE(get_integer_partitions(3, 2) == expected_partitions);
    }
}

TEST_CASE("Testing get_integer_partitions(3, 3)", "[math_util]") {
    SECTION("Testing get_integer_partitions") {
        std::vector< std::vector<unsigned int> > expected_partitions;
        expected_partitions.push_back({1, 1, 1});
        REQUIRE(get_integer_partitions(3, 3) == expected_partitions);
    }
}

TEST_CASE("Testing get_integer_partitions(4, 2)", "[math_util]") {
    SECTION("Testing get_integer_partitions") {
        std::vector< std::vector<unsigned int> > expected_partitions;
        expected_partitions.push_back({3, 1});
        expected_partitions.push_back({2, 2});
        REQUIRE(get_integer_partitions(4, 2) == expected_partitions);
    }
}

TEST_CASE("Testing get_integer_partitions(4, 3)", "[math_util]") {
    SECTION("Testing get_integer_partitions") {
        std::vector< std::vector<unsigned int> > expected_partitions;
        expected_partitions.push_back({2, 1, 1});
        REQUIRE(get_integer_partitions(4, 3) == expected_partitions);
    }
}

TEST_CASE("Testing get_integer_partitions(4, 4)", "[math_util]") {
    SECTION("Testing get_integer_partitions") {
        std::vector< std::vector<unsigned int> > expected_partitions;
        expected_partitions.push_back({1, 1, 1, 1});
        REQUIRE(get_integer_partitions(4, 4) == expected_partitions);
    }
}

TEST_CASE("Testing get_integer_partitions(5, 2)", "[math_util]") {
    SECTION("Testing get_integer_partitions") {
        std::vector< std::vector<unsigned int> > expected_partitions;
        expected_partitions.push_back({4, 1});
        expected_partitions.push_back({3, 2});
        REQUIRE(get_integer_partitions(5, 2) == expected_partitions);
    }
}

TEST_CASE("Testing get_integer_partitions(5, 3)", "[math_util]") {
    SECTION("Testing get_integer_partitions") {
        std::vector< std::vector<unsigned int> > expected_partitions;
        expected_partitions.push_back({3, 1, 1});
        expected_partitions.push_back({2, 2, 1});
        REQUIRE(get_integer_partitions(5, 3) == expected_partitions);
    }
}

TEST_CASE("Testing get_integer_partitions(5, 4)", "[math_util]") {
    SECTION("Testing get_integer_partitions") {
        std::vector< std::vector<unsigned int> > expected_partitions;
        expected_partitions.push_back({2, 1, 1, 1});
        REQUIRE(get_integer_partitions(5, 4) == expected_partitions);
    }
}

TEST_CASE("Testing get_integer_partitions(5, 5)", "[math_util]") {
    SECTION("Testing get_integer_partitions") {
        std::vector< std::vector<unsigned int> > expected_partitions;
        expected_partitions.push_back({1, 1, 1, 1, 1});
        REQUIRE(get_integer_partitions(5, 5) == expected_partitions);
    }
}

TEST_CASE("Testing get_integer_partitions(6, 2)", "[math_util]") {
    SECTION("Testing get_integer_partitions") {
        std::vector< std::vector<unsigned int> > expected_partitions;
        expected_partitions.push_back({5, 1});
        expected_partitions.push_back({4, 2});
        expected_partitions.push_back({3, 3});
        REQUIRE(get_integer_partitions(6, 2) == expected_partitions);
    }
}

TEST_CASE("Testing get_integer_partitions(6, 3)", "[math_util]") {
    SECTION("Testing get_integer_partitions") {
        std::vector< std::vector<unsigned int> > expected_partitions;
        expected_partitions.push_back({4, 1, 1});
        expected_partitions.push_back({3, 2, 1});
        expected_partitions.push_back({2, 2, 2});
        REQUIRE(get_integer_partitions(6, 3) == expected_partitions);
    }
}

TEST_CASE("Testing get_integer_partitions(6, 4)", "[math_util]") {
    SECTION("Testing get_integer_partitions") {
        std::vector< std::vector<unsigned int> > expected_partitions;
        expected_partitions.push_back({3, 1, 1, 1});
        expected_partitions.push_back({2, 2, 1, 1});
        REQUIRE(get_integer_partitions(6, 4) == expected_partitions);
    }
}

TEST_CASE("Testing get_integer_partitions(6, 5)", "[math_util]") {
    SECTION("Testing get_integer_partitions") {
        std::vector< std::vector<unsigned int> > expected_partitions;
        expected_partitions.push_back({2, 1, 1, 1, 1});
        REQUIRE(get_integer_partitions(6, 5) == expected_partitions);
    }
}

TEST_CASE("Testing get_integer_partitions(6, 6)", "[math_util]") {
    SECTION("Testing get_integer_partitions") {
        std::vector< std::vector<unsigned int> > expected_partitions;
        expected_partitions.push_back({1, 1, 1, 1, 1, 1});
        REQUIRE(get_integer_partitions(6, 6) == expected_partitions);
    }
}

TEST_CASE("Testing get_integer_partitions(11, 4)", "[math_util]") {
    SECTION("Testing get_integer_partitions") {
        std::vector< std::vector<unsigned int> > expected_partitions;
        expected_partitions.push_back({8, 1, 1, 1});
        expected_partitions.push_back({7, 2, 1, 1});
        expected_partitions.push_back({6, 3, 1, 1});
        expected_partitions.push_back({5, 4, 1, 1});
        expected_partitions.push_back({6, 2, 2, 1});
        expected_partitions.push_back({5, 3, 2, 1});
        expected_partitions.push_back({4, 4, 2, 1});
        expected_partitions.push_back({4, 3, 3, 1});
        expected_partitions.push_back({5, 2, 2, 2});
        expected_partitions.push_back({4, 3, 2, 2});
        expected_partitions.push_back({3, 3, 3, 2});
        REQUIRE(get_integer_partitions(11, 4) == expected_partitions);
    }
}

TEST_CASE("Testing log_factorial", "[math_util]") {
    SECTION("Testing log_factorial(0)") {
        REQUIRE(log_factorial(0) == 0.0);
    }
    SECTION("Testing log_factorial(1)") {
        REQUIRE(log_factorial(1) == 0.0);
    }
    SECTION("Testing log_factorial(2)") {
        REQUIRE(log_factorial(2) == Approx(std::log(2)));
    }
    SECTION("Testing log_factorial(3)") {
        REQUIRE(log_factorial(3) == Approx(std::log(6)));
    }
    SECTION("Testing log_factorial(4)") {
        REQUIRE(log_factorial(4) == Approx(std::log(24)));
    }
    SECTION("Testing log_factorial(5)") {
        REQUIRE(log_factorial(5) == Approx(std::log(120)));
    }
    SECTION("Testing log_factorial(6)") {
        REQUIRE(log_factorial(6) == Approx(std::log(720)));
    }
    SECTION("Testing log_factorial(12)") {
        REQUIRE(log_factorial(12) == Approx(std::log(479001600)));
    }
}

TEST_CASE("Testing log_multinomial_coefficient", "[math_util]") {
    SECTION("Testing log_factorial(1, 1)") {
        std::vector<unsigned int> x {1, 1};
        REQUIRE(log_multinomial_coefficient(x) == Approx(std::log(2.0)));
    }

    SECTION("Testing log_factorial(2, 1)") {
        std::vector<unsigned int> x {2, 1};
        REQUIRE(log_multinomial_coefficient(x) == Approx(std::log(3.0)));
    }

    SECTION("Testing log_factorial(1, 2)") {
        std::vector<unsigned int> x {1, 2};
        REQUIRE(log_multinomial_coefficient(x) == Approx(std::log(3.0)));
    }

    SECTION("Testing log_factorial(2, 2)") {
        std::vector<unsigned int> x {2, 2};
        REQUIRE(log_multinomial_coefficient(x) == Approx(std::log(6.0)));
    }

    SECTION("Testing log_factorial(3, 1)") {
        std::vector<unsigned int> x {3, 1};
        REQUIRE(log_multinomial_coefficient(x) == Approx(std::log(4.0)));
    }

    SECTION("Testing log_factorial(1, 3)") {
        std::vector<unsigned int> x {1, 3};
        REQUIRE(log_multinomial_coefficient(x) == Approx(std::log(4.0)));
    }

    SECTION("Testing log_factorial(3, 2)") {
        std::vector<unsigned int> x {3, 2};
        REQUIRE(log_multinomial_coefficient(x) == Approx(std::log(10.0)));
    }

    SECTION("Testing log_factorial(2, 3)") {
        std::vector<unsigned int> x {2, 3};
        REQUIRE(log_multinomial_coefficient(x) == Approx(std::log(10.0)));
    }

    SECTION("Testing log_factorial(3, 3)") {
        std::vector<unsigned int> x {3, 3};
        REQUIRE(log_multinomial_coefficient(x) == Approx(std::log(20.0)));
    }

    SECTION("Testing log_factorial(1, 1, 1)") {
        std::vector<unsigned int> x {1, 1, 1};
        REQUIRE(log_multinomial_coefficient(x) == Approx(std::log(6.0)));
    }

    SECTION("Testing log_factorial(4, 3, 2)") {
        std::vector<unsigned int> x {4, 3, 2};
        REQUIRE(log_multinomial_coefficient(x) == Approx(std::log(1260.0)));
    }

    SECTION("Testing log_factorial(2, 4, 3)") {
        std::vector<unsigned int> x {2, 4, 3};
        REQUIRE(log_multinomial_coefficient(x) == Approx(std::log(1260.0)));
    }
}


TEST_CASE("Testing get_pyp_log_prior_probability with {0,0,0} alpha=3.0, discount=0.0",
        "[math_util]") {
    double concentration = 3.0;
    double discount = 0.0;
    double expected_prob = std::log(2.0/20.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0, 0};
        REQUIRE(get_pyp_log_prior_probability<unsigned int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0, 0};
        REQUIRE(get_pyp_log_prior_probability<int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0, 0};
        REQUIRE(get_pyp_log_prior_probability<size_t>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0', '0'};
        REQUIRE(get_pyp_log_prior_probability<char>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "000";
        REQUIRE(get_pyp_log_prior_probability(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_wdp_log_prior_probability with {0,0,0} alpha=3.0, discount=0.0",
        "[math_util]") {
    double concentration = 3.0;
    double discount = 0.0;
    double expected_prob = std::log(2.0/20.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0, 0};
        REQUIRE(get_wdp_log_prior_probability<unsigned int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0, 0};
        REQUIRE(get_wdp_log_prior_probability<int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0, 0};
        REQUIRE(get_wdp_log_prior_probability<size_t>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0', '0'};
        REQUIRE(get_wdp_log_prior_probability<char>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "000";
        REQUIRE(get_wdp_log_prior_probability(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_pyp_log_prior_probability with {0,0,1} alpha=3.0, discount=0.0",
        "[math_util]") {
    double concentration = 3.0;
    double discount = 0.0;
    double expected_prob = std::log(3.0/20.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0, 1};
        REQUIRE(get_pyp_log_prior_probability<unsigned int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0, 1};
        REQUIRE(get_pyp_log_prior_probability<int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0, 1};
        REQUIRE(get_pyp_log_prior_probability<size_t>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0', '1'};
        REQUIRE(get_pyp_log_prior_probability<char>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "001";
        REQUIRE(get_pyp_log_prior_probability(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_wdp_log_prior_probability with {0,0,1} alpha=3.0, discount=0.0",
        "[math_util]") {
    double concentration = 3.0;
    double discount = 0.0;
    double expected_prob = std::log(3.0/20.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0, 1};
        REQUIRE(get_wdp_log_prior_probability<unsigned int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0, 1};
        REQUIRE(get_wdp_log_prior_probability<int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0, 1};
        REQUIRE(get_wdp_log_prior_probability<size_t>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0', '1'};
        REQUIRE(get_wdp_log_prior_probability<char>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "001";
        REQUIRE(get_wdp_log_prior_probability(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_pyp_log_prior_probability with {0,1,0} alpha=3.0, discount=0.0",
        "[math_util]") {
    double concentration = 3.0;
    double discount = 0.0;
    double expected_prob = std::log(3.0/20.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1, 0};
        REQUIRE(get_pyp_log_prior_probability<unsigned int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1, 0};
        REQUIRE(get_pyp_log_prior_probability<int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1, 0};
        REQUIRE(get_pyp_log_prior_probability<size_t>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1', '0'};
        REQUIRE(get_pyp_log_prior_probability<char>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "010";
        REQUIRE(get_pyp_log_prior_probability(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_wdp_log_prior_probability with {0,1,0} alpha=3.0, discount=0.0",
        "[math_util]") {
    double concentration = 3.0;
    double discount = 0.0;
    double expected_prob = std::log(3.0/20.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1, 0};
        REQUIRE(get_wdp_log_prior_probability<unsigned int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1, 0};
        REQUIRE(get_wdp_log_prior_probability<int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1, 0};
        REQUIRE(get_wdp_log_prior_probability<size_t>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1', '0'};
        REQUIRE(get_wdp_log_prior_probability<char>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "010";
        REQUIRE(get_wdp_log_prior_probability(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_pyp_log_prior_probability with {0,1,1} alpha=3.0, discount=0.0",
        "[math_util]") {
    double concentration = 3.0;
    double discount = 0.0;
    double expected_prob = std::log(3.0/20.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1, 1};
        REQUIRE(get_pyp_log_prior_probability<unsigned int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1, 1};
        REQUIRE(get_pyp_log_prior_probability<int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1, 1};
        REQUIRE(get_pyp_log_prior_probability<size_t>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1', '1'};
        REQUIRE(get_pyp_log_prior_probability<char>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "011";
        REQUIRE(get_pyp_log_prior_probability(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_wdp_log_prior_probability with {0,1,1} alpha=3.0, discount=0.0",
        "[math_util]") {
    double concentration = 3.0;
    double discount = 0.0;
    double expected_prob = std::log(3.0/20.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1, 1};
        REQUIRE(get_wdp_log_prior_probability<unsigned int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1, 1};
        REQUIRE(get_wdp_log_prior_probability<int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1, 1};
        REQUIRE(get_wdp_log_prior_probability<size_t>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1', '1'};
        REQUIRE(get_wdp_log_prior_probability<char>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "011";
        REQUIRE(get_wdp_log_prior_probability(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_pyp_log_prior_probability with {0,1,2} alpha=3.0, discount=0.0",
        "[math_util]") {
    double concentration = 3.0;
    double discount = 0.0;
    double expected_prob = std::log(9.0/20.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1, 2};
        REQUIRE(get_pyp_log_prior_probability<unsigned int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1, 2};
        REQUIRE(get_pyp_log_prior_probability<int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1, 2};
        REQUIRE(get_pyp_log_prior_probability<size_t>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1', '2'};
        REQUIRE(get_pyp_log_prior_probability<char>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "012";
        REQUIRE(get_pyp_log_prior_probability(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_wdp_log_prior_probability with {0,1,2} alpha=3.0, discount=0.0",
        "[math_util]") {
    double concentration = 3.0;
    double discount = 0.0;
    double expected_prob = std::log(9.0/20.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1, 2};
        REQUIRE(get_wdp_log_prior_probability<unsigned int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1, 2};
        REQUIRE(get_wdp_log_prior_probability<int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1, 2};
        REQUIRE(get_wdp_log_prior_probability<size_t>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1', '2'};
        REQUIRE(get_wdp_log_prior_probability<char>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "012";
        REQUIRE(get_wdp_log_prior_probability(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_pyp_log_prior_probability with {0,0,0,0,0} alpha=1.0, discount=0.0",
        "[math_util]") {
    double concentration = 1.0;
    double discount = 0.0;
    double expected_prob = std::log(24.0/120.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0, 0, 0, 0};
        REQUIRE(get_pyp_log_prior_probability<unsigned int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0, 0, 0, 0};
        REQUIRE(get_pyp_log_prior_probability<int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0, 0, 0, 0};
        REQUIRE(get_pyp_log_prior_probability<size_t>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0', '0', '0', '0'};
        REQUIRE(get_pyp_log_prior_probability<char>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "00000";
        REQUIRE(get_pyp_log_prior_probability(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_wdp_log_prior_probability with {0,0,0,0,0} alpha=1.0, discount=0.0",
        "[math_util]") {
    double concentration = 1.0;
    double discount = 0.0;
    double expected_prob = std::log(24.0/120.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0, 0, 0, 0};
        REQUIRE(get_wdp_log_prior_probability<unsigned int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0, 0, 0, 0};
        REQUIRE(get_wdp_log_prior_probability<int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0, 0, 0, 0};
        REQUIRE(get_wdp_log_prior_probability<size_t>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0', '0', '0', '0'};
        REQUIRE(get_wdp_log_prior_probability<char>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "00000";
        REQUIRE(get_wdp_log_prior_probability(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_pyp_log_prior_probability with {0,1,2,3,4} alpha=1.0, discount=0.0",
        "[math_util]") {
    double concentration = 1.0;
    double discount = 0.0;
    double expected_prob = std::log(1.0/120.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1, 2, 3, 4};
        REQUIRE(get_pyp_log_prior_probability<unsigned int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1, 2, 3, 4};
        REQUIRE(get_pyp_log_prior_probability<int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1, 2, 3, 4};
        REQUIRE(get_pyp_log_prior_probability<size_t>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1', '2', '3', '4'};
        REQUIRE(get_pyp_log_prior_probability<char>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "01234";
        REQUIRE(get_pyp_log_prior_probability(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_wdp_log_prior_probability with {0,1,2,3,4} alpha=1.0, discount=0.0",
        "[math_util]") {
    double concentration = 1.0;
    double discount = 0.0;
    double expected_prob = std::log(1.0/120.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 1, 2, 3, 4};
        REQUIRE(get_wdp_log_prior_probability<unsigned int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 1, 2, 3, 4};
        REQUIRE(get_wdp_log_prior_probability<int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 1, 2, 3, 4};
        REQUIRE(get_wdp_log_prior_probability<size_t>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '1', '2', '3', '4'};
        REQUIRE(get_wdp_log_prior_probability<char>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "01234";
        REQUIRE(get_wdp_log_prior_probability(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_pyp_log_prior_probability with {0,0,1,0,2} alpha=1.0, discount=0.0",
        "[math_util]") {
    double concentration = 1.0;
    double discount = 0.0;
    double expected_prob = std::log(2.0/120.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0, 1, 0, 2};
        REQUIRE(get_pyp_log_prior_probability<unsigned int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0, 1, 0, 2};
        REQUIRE(get_pyp_log_prior_probability<int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0, 1, 0, 2};
        REQUIRE(get_pyp_log_prior_probability<size_t>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0', '1', '0', '2'};
        REQUIRE(get_pyp_log_prior_probability<char>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "00102";
        REQUIRE(get_pyp_log_prior_probability(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_wdp_log_prior_probability with {0,0,1,0,2} alpha=1.0, discount=0.0",
        "[math_util]") {
    double concentration = 1.0;
    double discount = 0.0;
    double expected_prob = std::log(2.0/120.0);
    SECTION("Testing vector of unsigned int") {
        std::vector<unsigned int> partition = {0, 0, 1, 0, 2};
        REQUIRE(get_wdp_log_prior_probability<unsigned int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of int") {
        std::vector<int> partition = {0, 0, 1, 0, 2};
        REQUIRE(get_wdp_log_prior_probability<int>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of size_t") {
        std::vector<size_t> partition = {0, 0, 1, 0, 2};
        REQUIRE(get_wdp_log_prior_probability<size_t>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing vector of char") {
        std::vector<char> partition = {'0', '0', '1', '0', '2'};
        REQUIRE(get_wdp_log_prior_probability<char>(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
    SECTION("Testing string") {
        std::string partition = "00102";
        REQUIRE(get_wdp_log_prior_probability(partition,
                concentration,
                discount) == Approx(expected_prob));
    }
}

TEST_CASE("Testing get_pyp_expected_number_of_categories 0.218, 0.0, 7282", "[math_util]") {
    double e = get_pyp_expected_number_of_categories(0.218, 0.0, 7282);
    REQUIRE(e == Approx(3.0).epsilon(0.001));
}

TEST_CASE("Testing get_pyp_expected_number_of_categories 1.068, 0.0, 7282", "[math_util]") {
    double e = get_pyp_expected_number_of_categories(1.068, 0.0, 7282);
    REQUIRE(e == Approx(10.0).epsilon(0.0001));
}

TEST_CASE("Testing ln_pochhammer", "[math_util]") {
    SECTION("testing ln_pochhammer(1, 5)") {
        REQUIRE(ln_pochhammer(1, 5) == Approx(std::log(120.0)));
        REQUIRE(ln_pochhammer(2.9, 10) == Approx(std::log(203697425)));
    }
}

TEST_CASE("Testing get_pyp_expected_number_of_categories 2.9, 0.5, 10", "[math_util]") {
    SECTION("Testing 2.9, 0.5, 10") {
        unsigned int num_elements = 10;
        double a = 2.9;
        double d = 0.5;
        double e = get_pyp_expected_number_of_categories(a, d, num_elements);

        std::vector<unsigned int> elements (num_elements, 0);
        REQUIRE(elements.size() == num_elements);
        SampleSummarizer<double> ncats;
        RandomNumberGenerator rng(342254);

        unsigned int nsamples = 100000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int n = rng.pitman_yor_process(elements, a, d);
            ncats.add_sample(float(n));
        }
        REQUIRE(ncats.mean() == Approx(e).epsilon(0.005));
    }
}

TEST_CASE("Testing get_pyp_expected_number_of_categories 2.9, 0.1, 10", "[math_util]") {
    SECTION("Testing 2.9, 0.1, 10") {
        unsigned int num_elements = 10;
        double a = 2.9;
        double d = 0.1;
        double e = get_pyp_expected_number_of_categories(a, d, num_elements);

        std::vector<unsigned int> elements (num_elements, 0);
        REQUIRE(elements.size() == num_elements);
        SampleSummarizer<double> ncats;
        RandomNumberGenerator rng(342254);

        unsigned int nsamples = 100000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int n = rng.pitman_yor_process(elements, a, d);
            ncats.add_sample(float(n));
        }
        REQUIRE(ncats.mean() == Approx(e).epsilon(0.005));
    }
}

TEST_CASE("Testing get_pyp_expected_number_of_categories 2.9, 0.9, 10", "[math_util]") {
    SECTION("Testing 2.9, 0.9, 10") {
        unsigned int num_elements = 10;
        double a = 2.9;
        double d = 0.9;
        double e = get_pyp_expected_number_of_categories(a, d, num_elements);

        std::vector<unsigned int> elements (num_elements, 0);
        REQUIRE(elements.size() == num_elements);
        SampleSummarizer<double> ncats;
        RandomNumberGenerator rng(342254);

        unsigned int nsamples = 100000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            unsigned int n = rng.pitman_yor_process(elements, a, d);
            ncats.add_sample(float(n));
        }
        REQUIRE(ncats.mean() == Approx(e).epsilon(0.005));
    }
}

TEST_CASE("Testing get_pyp_concentration 3.0, 0.0,  7282", "[math_util]") {
    double e = get_pyp_concentration(3.0, 7282, 0.0);
    REQUIRE(e == Approx(0.218).epsilon(0.001));
}
TEST_CASE("Testing get_pyp_concentration 10.0, 0.0, 7282", "[math_util]") {
    double e = get_pyp_concentration(10.0, 7282, 0.0);
    REQUIRE(e == Approx(1.068).epsilon(0.0001));
}

TEST_CASE("Testing get_pyp_concentration 2.9, 0.1, 10", "[math_util]") {
    SECTION("Testing 2.9, 0.1, 10") {
        unsigned int num_elements = 10;
        double a = 2.9;
        double d = 0.1;
        double e = get_pyp_expected_number_of_categories(a, d, num_elements);

        double conc = get_pyp_concentration(e, num_elements, d);

        REQUIRE(conc == Approx(a).epsilon(0.001));
    }
}

TEST_CASE("Testing get_pyp_concentration 2.9, 0.5, 10", "[math_util]") {
    SECTION("Testing 2.9, 0.5, 10") {
        unsigned int num_elements = 10;
        double a = 2.9;
        double d = 0.5;
        double e = get_pyp_expected_number_of_categories(a, d, num_elements);

        double conc = get_pyp_concentration(e, num_elements, d);

        REQUIRE(conc == Approx(a).epsilon(0.001));
    }
}

TEST_CASE("Testing get_pyp_concentration 2.9, 0.9, 10", "[math_util]") {
    SECTION("Testing 2.9, 0.9, 10") {
        unsigned int num_elements = 10;
        double a = 2.9;
        double d = 0.9;
        double e = get_pyp_expected_number_of_categories(a, d, num_elements);

        double conc = get_pyp_concentration(e, num_elements, d);

        REQUIRE(conc == Approx(a).epsilon(0.001));
    }
}

TEST_CASE("Testing get_pyp_concentration_gamma_scale", "[math_util]") {
    SECTION("Testing 2.9, 0.5, 10") {
        unsigned int num_elements = 10;
        double a = 2.9;
        double d = 0.5;
        double e = get_pyp_expected_number_of_categories(a, d, num_elements);
        std::vector<double> shapes {0.001, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0, 10000.0};
        for (auto shape : shapes) {
            double scale = get_pyp_concentration_gamma_scale(e, num_elements, d, shape);
            REQUIRE(a == Approx(shape * scale).epsilon(0.001));
        }
    }
}

TEST_CASE("Testing get_dpp_log_prior_probability exchangeability with 6 elements sets of size 1, 2, 3",
        "[math_util]") {
    double concentration = 1.0;
    std::vector<unsigned int> partition = {0, 0, 0, 1, 1, 2};
    double log_p = get_dpp_log_prior_probability<unsigned int>(partition, concentration);
    std::vector< std::vector<unsigned int> > permutations = {
            {0, 1, 1, 2, 2, 2},
            {0, 0, 1, 2, 2, 2},
            {0, 1, 1, 1, 2, 2},
            {0, 0, 1, 1, 1, 2},
            {0, 1, 0, 1, 0, 2},
            {0, 1, 1, 0, 0, 2},
            {0, 1, 2, 0, 0, 1},
            {0, 1, 2, 1, 0, 0},
            {0, 1, 2, 0, 1, 0},
    };
    for (unsigned int i = 0; i < permutations.size(); ++i) {
        REQUIRE(log_p == get_dpp_log_prior_probability(permutations.at(i), concentration));
    }
}

TEST_CASE("Testing get_pyp_log_prior_probability exchangeability with 6 elements sets of size 1, 2, 3",
        "[math_util]") {
    double concentration = 1.0;
    double discount = 0.5;
    std::vector<unsigned int> partition = {0, 0, 0, 1, 1, 2};
    double log_p = get_pyp_log_prior_probability<unsigned int>(partition, concentration, discount);
    std::vector< std::vector<unsigned int> > permutations = {
            {0, 1, 1, 2, 2, 2},
            {0, 0, 1, 2, 2, 2},
            {0, 1, 1, 1, 2, 2},
            {0, 0, 1, 1, 1, 2},
            {0, 1, 0, 1, 0, 2},
            {0, 1, 1, 0, 0, 2},
            {0, 1, 2, 0, 0, 1},
            {0, 1, 2, 1, 0, 0},
            {0, 1, 2, 0, 1, 0},
    };
    for (unsigned int i = 0; i < permutations.size(); ++i) {
        REQUIRE(log_p == get_pyp_log_prior_probability(permutations.at(i), concentration, discount));
    }
}

// This test confirms that the elements of the weighted-discount process are
// not exchangeable, so the Gibbs sampler will not work.
// TEST_CASE("Testing get_wdp_log_prior_probability exchangeability with 6 elements sets of size 1, 2, 3",
//         "[math_util]") {
//     double concentration = 1.0;
//     double discount = 0.5;
//     std::vector<unsigned int> partition = {0, 0, 0, 1, 1, 2};
//     double log_p = get_wdp_log_prior_probability<unsigned int>(partition, concentration, discount);
//     std::vector< std::vector<unsigned int> > permutations = {
//             {0, 1, 1, 2, 2, 2},
//             {0, 0, 1, 2, 2, 2},
//             {0, 1, 1, 1, 2, 2},
//             {0, 0, 1, 1, 1, 2},
//             {0, 1, 0, 1, 0, 2},
//             {0, 1, 1, 0, 0, 2},
//             {0, 1, 2, 0, 0, 1},
//             {0, 1, 2, 1, 0, 0},
//             {0, 1, 2, 0, 1, 0},
//     };
//     for (unsigned int i = 0; i < permutations.size(); ++i) {
//         REQUIRE(log_p == get_wdp_log_prior_probability(permutations.at(i), concentration, discount));
//     }
// }
