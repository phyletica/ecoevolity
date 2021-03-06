/******************************************************************************
 * Copyright (C) 2015-2016 Jamie R. Oaks.
 *
 * This file is part of Ecoevolity.
 *
 * Ecoevolity is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Ecoevolity is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Ecoevolity.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

#ifndef ECOEVOLITY_RNG_HPP
#define ECOEVOLITY_RNG_HPP

#include <ctime>
#include <random>
#include <unordered_map>

#include "assert.hpp"
#include "math_util.hpp"

class RandomNumberGenerator {
    public:
        typedef typename std::mt19937 EngineType;
        typedef typename EngineType::result_type RandomSeedType;
        EngineType engine_;
        RandomSeedType seed_;
        std::uniform_real_distribution<double> uniform_real_rng_;
        std::uniform_int_distribution<long> uniform_int_rng_;
        std::uniform_int_distribution<unsigned long> uniform_positive_int_rng_;
        std::gamma_distribution<double> gamma_rng_;

        RandomNumberGenerator() {
            this->set_seed_from_time();
        }

        RandomNumberGenerator(RandomSeedType seed) {
            this->set_seed(seed);
        }

        ~RandomNumberGenerator() { }

        RandomSeedType get_seed() const {
            return this->seed_;
        }

        void set_seed(RandomSeedType seed) {
            this->seed_ = seed;
            this->engine_.seed(seed);
        }

        void set_seed_from_time() {
            this->set_seed(static_cast<RandomSeedType>(std::time(NULL)));
        }

        inline double uniform_real() {
            return this->uniform_real_rng_(
                    this->engine_,
                    typename decltype(this->uniform_real_rng_)::param_type(
                            0.0, 1.0));
        }

        inline double uniform_real(double a, double b) {
            return this->uniform_real_rng_(
                    this->engine_,
                    typename decltype(this->uniform_real_rng_)::param_type(
                            a, b));
        }

        inline long uniform_int(long a, long b) {
            return this->uniform_int_rng_(
                    this->engine_,
                    typename decltype(this->uniform_int_rng_)::param_type(
                            a, b));
        }

        inline long uniform_positive_int(unsigned long a, unsigned long b) {
            return this->uniform_positive_int_rng_(
                    this->engine_,
                    typename decltype(this->uniform_positive_int_rng_)::param_type(
                            a, b));
        }

        inline long uniform_positive_int(unsigned long b) {
            return this->uniform_positive_int_rng_(
                    this->engine_,
                    typename decltype(this->uniform_positive_int_rng_)::param_type(
                            0, b));
        }

        inline double gamma() {
            return this->gamma_rng_(
                    this->engine_,
                    typename decltype(this->gamma_rng_)::param_type(
                            1.0, 1.0));
        }
        inline double gamma(double shape, double scale) {
            return this->gamma_rng_(
                    this->engine_,
                    typename decltype(this->gamma_rng_)::param_type(
                            shape, scale));
        }

        inline double beta(double alpha, double beta) {
            double x = this->gamma(alpha, 1.0);
            double y = this->gamma(beta, 1.0);
            return (x / (x + y));
        }

        inline std::vector<double> dirichlet(const std::vector<double> & parameters) {
            ECOEVOLITY_ASSERT(parameters.size() > 1);
            std::vector<double> x;
            x.reserve(parameters.size());
            double sum = 0.0;
            for (unsigned int i = 0; i < parameters.size(); ++i) {
                double g = this->gamma(parameters.at(i), 1.0);
                x.push_back(g);
                sum += g;
            }
            for (unsigned int i = 0; i < parameters.size(); ++i) {
                x.at(i) = x.at(i) / sum;
            }
            return x;
        }

        inline unsigned int weighted_index(
                const std::vector<double>& probabilities) {
            double u = this->uniform_real();
            for (unsigned int i = 0; i < probabilities.size(); ++i) {
                u -= probabilities.at(i);
                if (u < 0.0) {
                    return i;
                }
            }
            ECOEVOLITY_ASSERT_APPROX_EQUAL(u, 0.0);
            return probabilities.size() - 1;
        }

        inline unsigned int weighted_index(
                const std::vector<long double>& probabilities) {
            long double u = this->uniform_real();
            for (unsigned int i = 0; i < probabilities.size(); ++i) {
                u -= probabilities.at(i);
                if (u < 0.0) {
                    return i;
                }
            }
            ECOEVOLITY_ASSERT_APPROX_EQUAL(u, 0.0);
            return probabilities.size() - 1;
        }

        inline std::vector<unsigned int> random_subset_indices(
                unsigned int number_of_elements,
                unsigned int subset_size) {
            ECOEVOLITY_ASSERT((number_of_elements > 0) &&
                              (subset_size > 0) &&
                              (number_of_elements >= subset_size));
            std::vector<unsigned int> indices;
            indices.reserve(subset_size);
            if (subset_size == number_of_elements) {
                for (unsigned int i = 0; i < number_of_elements; ++i) {
                    indices.push_back(i);
                }
                return indices;
            }

            unsigned int total = 0;
            double u;

            while (indices.size() < subset_size) {
                u = this->uniform_real();
                if (((number_of_elements - total) * u) >= (subset_size - indices.size())) {
                    ++total;
                }
                else {
                    indices.push_back(total);
                    ++total;
                }
            }
            return indices;
        }

        inline std::string random_string(
                unsigned int length,
                std::vector<char> char_pool = {
                    '0', '1', '2', '3', '4',
                    '5', '6', '7', '8', '9',
                    'a', 'b', 'c', 'd', 'e',
                    'f', 'g', 'h', 'i', 'j',
                    'k', 'l', 'm', 'n', 'o',
                    'p', 'q', 'r', 's', 't',
                    'u', 'v', 'w', 'x', 'y',
                    'z',
                    'A', 'B', 'C', 'D', 'E',
                    'F', 'G', 'H', 'I', 'J',
                    'K', 'L', 'M', 'N', 'O',
                    'P', 'Q', 'R', 'S', 'T',
                    'U', 'V', 'W', 'X', 'Y',
                    'Z'
                }) {
            ECOEVOLITY_ASSERT(char_pool.size() > 0);
            std::string s = "";
            for (unsigned int i = 0; i < length; ++i) {
                s += char_pool.at(this->uniform_int(0, char_pool.size() - 1));
            }
            return s;
        }

        /** 
         * A function for generating a random draw from a Dirichlet process.
         */
        inline unsigned int dirichlet_process(
                std::vector<unsigned int>& elements,
                double alpha) {
            ECOEVOLITY_ASSERT (alpha > 0.0);
            unsigned int n = elements.size();
            ECOEVOLITY_ASSERT(n > 0);
            double subset_prob;
            double new_subset_prob;
            double u;
            std::vector<unsigned int> subset_counts;
            subset_counts.reserve(n);
            subset_counts.push_back(1);
            elements.at(0) = 0;
            unsigned int num_subsets = 1;
            for (unsigned int i = 1; i < n; ++i) {
                new_subset_prob = (alpha / (alpha + (double)i));
                u = this->uniform_real();
                u -= new_subset_prob;
                if (u < 0.0) {
                    elements.at(i) = num_subsets;
                    subset_counts.push_back(1);
                    ++num_subsets;
                    continue;
                }
                for (unsigned int j = 0; j < num_subsets; ++j) {
                    subset_prob = ((double)subset_counts.at(j) / (alpha + (double)i));
                    u -= subset_prob;
                    if (u < 0.0) {
                        elements.at(i) = j;
                        ++subset_counts.at(j);
                        break;
                    }
                }
                if (u > 0.0) {
                    elements.at(i) = num_subsets - 1;
                    ++subset_counts.at(num_subsets - 1);
                }
            }
            return num_subsets;
        }

        inline std::pair< unsigned int, std::vector<unsigned int> > dirichlet_process(
                unsigned int number_of_elements,
                double alpha) {
            std::vector<unsigned int> elements (number_of_elements, 0);
            unsigned int num_subsets = this->dirichlet_process(elements, alpha);
            return std::make_pair(num_subsets, elements);
        }

        /**
         * A function for generating a random set partition.
         */
        inline unsigned int random_set_partition(
                std::vector<unsigned int>& elements,
                double split_weight = 1.0) {
            // This is very inefficent. It would be better to draw the sizes of
            // the 'ncats' categories weighted by how many possible partitions
            // there are for each (i.e., the sum of all possible set partitions
            // for each possible integer partition with 'ncats' categories.
            // NOTE: adding 'ncats' categories to random elements, and then
            // randomly assigning the remaining elements to categories (or
            // using shuffling) does NOT work. This produces a multinomial
            // distribution over partitions, which is not uniform.

            ECOEVOLITY_ASSERT (split_weight > 0.0);
            unsigned int n = elements.size();
            ECOEVOLITY_ASSERT(n > 0);
            std::vector<long double> ncat_probs;
            ncat_probs.reserve(n);
            long double denom = 0.0;
            long double p = 0.0;
            for (unsigned int k = 1; k <= n; ++k) {
                p = stirling2_base<long double>(n, k) * std::pow(split_weight, (k - 1));
                ncat_probs.push_back(p);
                denom += p;
            }
            for (unsigned int i = 0; i < n; ++i) {
                ncat_probs.at(i) = ncat_probs.at(i) / denom;
            }
            unsigned int ncats = this->weighted_index(ncat_probs) + 1;
            if (ncats == 1) {
                for (unsigned int i = 0; i < n; ++i) {
                    elements.at(i) = 0;
                }
                return ncats;
            }
            for (unsigned int i = 0; i < n; ++i) {
                elements.at(i) = i;
            }
            if (ncats == n) {
                for (unsigned int i = 0; i < n; ++i) {
                    elements.at(i) = i;
                }
                return ncats;
            }
            std::unordered_map<unsigned int, unsigned int> standardizing_map;
            standardizing_map.reserve(ncats);
            while(true) {
                standardizing_map.clear();
                unsigned int next_idx = 0;
                for (unsigned int i = 0; i < n; ++i) {
                    unsigned int raw_idx = this->uniform_int(0, ncats-1);
                    if (standardizing_map.count(raw_idx) == 0) {
                        standardizing_map[raw_idx] = next_idx;
                        ++next_idx;
                    }
                    elements.at(i) = standardizing_map.at(raw_idx);
                }
                if (next_idx == ncats) {
                    break;
                }
            }
            return ncats;
        }

        inline std::pair< unsigned int, std::vector<unsigned int> > random_set_partition(
                unsigned int number_of_elements,
                double split_weight = 1.0) {
            std::vector<unsigned int> elements (number_of_elements, 0);
            unsigned int num_subsets = this->random_set_partition(elements, split_weight);
            return std::make_pair(num_subsets, elements);
        }
};

#endif
