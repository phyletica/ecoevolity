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

        template <typename T>
        inline void shuffle(std::vector<T> & vec) {
            std::shuffle(std::begin(vec), std::end(vec), this->engine_);
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
         * A function for generating a random draw from a Pitman-Yor process.
         */
        inline unsigned int pitman_yor_process(
                std::vector<unsigned int>& elements,
                double alpha,
                double discount) {
            ECOEVOLITY_ASSERT ((discount >= 0.0) && (discount < 1.0));
            ECOEVOLITY_ASSERT (alpha > -discount);
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
                new_subset_prob = ((alpha + (discount * num_subsets)) /
                        (alpha + (double)i));
                u = this->uniform_real();
                u -= new_subset_prob;
                if (u < 0.0) {
                    elements.at(i) = num_subsets;
                    subset_counts.push_back(1);
                    ++num_subsets;
                    continue;
                }
                for (unsigned int j = 0; j < num_subsets; ++j) {
                    subset_prob = (((double)subset_counts.at(j) - discount) /
                            (alpha + (double)i));
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

        inline std::pair< unsigned int, std::vector<unsigned int> > pitman_yor_process(
                unsigned int number_of_elements,
                double alpha,
                double discount) {
            std::vector<unsigned int> elements (number_of_elements, 0);
            unsigned int num_subsets = this->pitman_yor_process(elements, alpha,
                    discount);
            return std::make_pair(num_subsets, elements);
        }

        /** 
         * A function for generating a random draw from a weighted-discount process.
         */
        inline unsigned int weighted_discount_process(
                std::vector<unsigned int>& elements,
                double alpha,
                double discount) {
            ECOEVOLITY_ASSERT ((discount >= 0.0) && (discount < 1.0));
            ECOEVOLITY_ASSERT (alpha > -discount);
            unsigned int n = elements.size();
            ECOEVOLITY_ASSERT(n > 0);
            double subset_prob;
            double new_subset_prob;
            double u;
            std::vector<unsigned int> subset_counts;
            subset_counts.reserve(n);
            subset_counts.push_back(1);
            elements.at(0) = 0;
            unsigned n_elements = 1;
            unsigned int num_subsets = 1;
            for (unsigned int i = 1; i < n; ++i) {
                new_subset_prob = ((alpha + (discount * n_elements)) /
                        (alpha + (double)i));
                u = this->uniform_real();
                u -= new_subset_prob;
                if (u < 0.0) {
                    elements.at(i) = num_subsets;
                    ++n_elements;
                    subset_counts.push_back(1);
                    ++num_subsets;
                    continue;
                }
                for (unsigned int j = 0; j < num_subsets; ++j) {
                    subset_prob = (((double)subset_counts.at(j) - (discount * (double)subset_counts.at(j))) /
                            (alpha + (double)i));
                    u -= subset_prob;
                    if (u < 0.0) {
                        elements.at(i) = j;
                        ++n_elements;
                        ++subset_counts.at(j);
                        break;
                    }
                }
                if (u > 0.0) {
                    elements.at(i) = num_subsets - 1;
                    ++n_elements;
                    ++subset_counts.at(num_subsets - 1);
                }
            }
            return num_subsets;
        }

        inline std::pair< unsigned int, std::vector<unsigned int> > weighted_discount_process(
                unsigned int number_of_elements,
                double alpha,
                double discount) {
            std::vector<unsigned int> elements (number_of_elements, 0);
            unsigned int num_subsets = this->weighted_discount_process(elements, alpha,
                    discount);
            return std::make_pair(num_subsets, elements);
        }

        inline void random_subsets_rev(
                std::vector< std::vector<unsigned int> > &subsets,
                unsigned int N,
                unsigned int k) {
            if (N == 0) {
                return;
            }
            else if (N == 1){
                subsets.push_back({0});
                return;
            }
            else if (
                    (stirling2_base<long double>(N-1, k-1) /
                     stirling2_base<long double>(N, k)) >
                    this->uniform_real()) {
                subsets.push_back({N-1});
                std::vector< std::vector<unsigned int> > remaining_subsets;
                remaining_subsets.reserve(k-1);
                this->random_subsets_rev(remaining_subsets, N-1, k-1);
                for (unsigned int subset_idx = 0;
                        subset_idx < remaining_subsets.size();
                        ++subset_idx) {
                    subsets.push_back(remaining_subsets.at(subset_idx));
                }

                return;
            }
            else {
                this->random_subsets_rev(subsets, N-1, k);
                subsets.at(this->uniform_positive_int(subsets.size() - 1)).push_back(N-1);
            }

        }

        inline void random_subsets(
                std::vector< std::vector<unsigned int> > &subsets,
                unsigned int number_of_elements,
                unsigned int number_of_subsets) {
            this->random_subsets_rev(subsets,
                    number_of_elements,
                    number_of_subsets);
            std::reverse(std::begin(subsets), std::end(subsets));
        }

        inline std::vector< std::vector<unsigned int> > random_subsets(
                unsigned int number_of_elements,
                unsigned int number_of_subsets) {
            std::vector< std::vector<unsigned int> > subsets;
            subsets.reserve(number_of_subsets);
            this->random_subsets(subsets,
                    number_of_elements,
                    number_of_subsets);
            return subsets;
        }

        inline void random_set_partition_with_k_subsets(
                std::vector<unsigned int>& elements,
                unsigned int number_of_subsets) {
            std::vector< std::vector<unsigned int> > subsets;
            subsets.reserve(number_of_subsets);
            this->random_subsets(subsets,
                    elements.size(),
                    number_of_subsets);
            for (unsigned int subset_idx = 0;
                    subset_idx < subsets.size();
                    ++subset_idx) {
                for (auto element_idx : subsets.at(subset_idx)) {
                    elements.at(element_idx) = subset_idx;
                }
            }
            return;
        }

        inline std::vector<unsigned int> random_set_partition_with_k_subsets(
                unsigned int number_of_elements,
                unsigned int number_of_subsets) {
            std::vector<unsigned int> elements(number_of_elements);
            this->random_set_partition_with_k_subsets(
                    elements,
                    number_of_subsets);
            return elements;
        }

        inline unsigned int random_number_of_subsets(
                unsigned int number_of_elements,
                double split_weight = 1.0) {
            ECOEVOLITY_ASSERT(split_weight > 0.0);
            ECOEVOLITY_ASSERT(number_of_elements > 0);
            std::vector<long double> ncat_probs;
            ncat_probs.reserve(number_of_elements);
            long double denom = 0.0;
            long double p = 0.0;
            for (unsigned int k = 1; k <= number_of_elements; ++k) {
                p = stirling2_base<long double>(number_of_elements, k) * std::pow(split_weight, (k - 1));
                ncat_probs.push_back(p);
                denom += p;
            }
            for (unsigned int i = 0; i < number_of_elements; ++i) {
                ncat_probs.at(i) = ncat_probs.at(i) / denom;
            }
            unsigned int ncats = this->weighted_index(ncat_probs) + 1;
            return ncats;
        }

        inline unsigned int restricted_random_number_of_subsets(
                unsigned int number_of_elements,
                const std::vector<unsigned int> & possible_numbers_of_subsets,
                double split_weight = 1.0) {
            ECOEVOLITY_ASSERT(split_weight > 0.0);
            ECOEVOLITY_ASSERT(number_of_elements > 0);
            ECOEVOLITY_ASSERT(possible_numbers_of_subsets.size() > 0);
            if (possible_numbers_of_subsets.size() == 1) {
                return possible_numbers_of_subsets.at(0);
            }
            std::vector<long double> ncat_probs;
            ncat_probs.reserve(possible_numbers_of_subsets.size());
            long double denom = 0.0;
            long double p = 0.0;
            for (auto k : possible_numbers_of_subsets) {
                ECOEVOLITY_ASSERT((k > 0) && (k <= number_of_elements));
                p = stirling2_base<long double>(number_of_elements, k) * std::pow(split_weight, (k - 1));
                ncat_probs.push_back(p);
                denom += p;
            }
            for (unsigned int i = 0; i < ncat_probs.size(); ++i) {
                ncat_probs.at(i) = ncat_probs.at(i) / denom;
            }
            unsigned int ncats = this->weighted_index(ncat_probs) + 1;
            return ncats;
        }

        /**
         * A function for generating a random set partition.
         */
        inline unsigned int random_set_partition(
                std::vector<unsigned int>& elements,
                double split_weight = 1.0) {
            unsigned int n = elements.size();
            unsigned int ncats = random_number_of_subsets(n, split_weight);
            if (ncats == 1) {
                for (unsigned int i = 0; i < n; ++i) {
                    elements.at(i) = 0;
                }
                return ncats;
            }
            if (ncats == n) {
                for (unsigned int i = 0; i < n; ++i) {
                    elements.at(i) = i;
                }
                return ncats;
            }
            // No longer need rejection hack to uniformly sample set partitions
            // std::unordered_map<unsigned int, unsigned int> standardizing_map;
            // standardizing_map.reserve(ncats);
            // while(true) {
            //     standardizing_map.clear();
            //     unsigned int next_idx = 0;
            //     for (unsigned int i = 0; i < n; ++i) {
            //         unsigned int raw_idx = this->uniform_int(0, ncats-1);
            //         if (standardizing_map.count(raw_idx) == 0) {
            //             standardizing_map[raw_idx] = next_idx;
            //             ++next_idx;
            //         }
            //         elements.at(i) = standardizing_map.at(raw_idx);
            //     }
            //     if (next_idx == ncats) {
            //         break;
            //     }
            // }
            this->random_set_partition_with_k_subsets(elements, ncats);
            return ncats;
        }

        inline std::pair< unsigned int, std::vector<unsigned int> > random_set_partition(
                unsigned int number_of_elements,
                double split_weight = 1.0) {
            std::vector<unsigned int> elements (number_of_elements, 0);
            unsigned int num_subsets = this->random_set_partition(elements, split_weight);
            return std::make_pair(num_subsets, elements);
        }


        /**
         * A function for generating a random set partition.
         */
        inline unsigned int random_set_partition_as_subsets(
                std::vector< std::vector<unsigned int> > &subsets,
                unsigned int number_of_elements,
                double split_weight = 1.0) {
            unsigned int n = number_of_elements;
            unsigned int ncats = random_number_of_subsets(n, split_weight);
            subsets.reserve(ncats);
            if (ncats == 1) {
                std::vector<unsigned int> sset(number_of_elements);
                for (unsigned int i = 0; i < n; ++i) {
                    sset.at(i) = i;
                }
                subsets.push_back(sset);
                return ncats;
            }
            if (ncats == n) {
                for (unsigned int i = 0; i < n; ++i) {
                    std::vector<unsigned int> sset = {i};
                    subsets.push_back(sset);
                }
                return ncats;
            }
            // No longer need rejection hack to uniformly sample set partitions
            // std::unordered_map<unsigned int, unsigned int> standardizing_map;
            // standardizing_map.reserve(ncats);
            // while(true) {
            //     standardizing_map.clear();
            //     unsigned int next_idx = 0;
            //     for (unsigned int i = 0; i < n; ++i) {
            //         unsigned int raw_idx = this->uniform_int(0, ncats-1);
            //         if (standardizing_map.count(raw_idx) == 0) {
            //             standardizing_map[raw_idx] = next_idx;
            //             ++next_idx;
            //         }
            //         elements.at(i) = standardizing_map.at(raw_idx);
            //     }
            //     if (next_idx == ncats) {
            //         break;
            //     }
            // }
            this->random_subsets(subsets, n, ncats);
            return ncats;
        }

        inline std::vector< std::vector<unsigned int> > random_set_partition_as_subsets(
                unsigned int number_of_elements,
                double split_weight = 1.0) {
            std::vector< std::vector<unsigned int> > subsets;
            unsigned int num_subsets = this->random_set_partition_as_subsets(
                    subsets,
                    number_of_elements,
                    split_weight);
            return subsets;
        }

        /**
         * A function for generating a random set partition conditional on a
         * restricted possible number of subsets.
         */
        inline unsigned int restricted_random_set_partition_as_subsets(
                std::vector< std::vector<unsigned int> > &subsets,
                unsigned int number_of_elements,
                const std::vector<unsigned int> & possible_number_of_subsets,
                double split_weight = 1.0) {
            unsigned int n = number_of_elements;
            unsigned int ncats = restricted_random_number_of_subsets(n,
                    possible_number_of_subsets, split_weight);
            subsets.reserve(ncats);
            if (ncats == 1) {
                std::vector<unsigned int> sset(number_of_elements);
                for (unsigned int i = 0; i < n; ++i) {
                    sset.at(i) = i;
                }
                subsets.push_back(sset);
                return ncats;
            }
            if (ncats == n) {
                for (unsigned int i = 0; i < n; ++i) {
                    std::vector<unsigned int> sset = {i};
                    subsets.push_back(sset);
                }
                return ncats;
            }
            this->random_subsets(subsets, n, ncats);
            return ncats;
        }

        inline std::vector< std::vector<unsigned int> > restricted_random_set_partition_as_subsets(
                unsigned int number_of_elements,
                const std::vector<unsigned int> & possible_number_of_subsets,
                double split_weight = 1.0) {
            std::vector< std::vector<unsigned int> > subsets;
            unsigned int num_subsets = this->restricted_random_set_partition_as_subsets(
                    subsets,
                    number_of_elements,
                    possible_number_of_subsets,
                    split_weight);
            return subsets;
        }
};

#endif
