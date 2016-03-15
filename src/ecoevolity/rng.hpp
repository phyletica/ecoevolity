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

class RandomNumberGenerator {
    public:
        typedef typename std::mt19937 EngineType;
        typedef typename EngineType::result_type RandomSeedType;
        EngineType engine_;
        RandomSeedType seed_;
        std::uniform_real_distribution<double> uniform_real_rng_;
        std::uniform_int_distribution<long> uniform_int_rng_;
        std::uniform_int_distribution<unsigned long> uniform_positive_int_rng_;

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

        inline long uniform_positive_int(unsigned int a, unsigned long b) {
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
};

#endif
