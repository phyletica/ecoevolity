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

#include "assert.hpp"

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
};

#endif
