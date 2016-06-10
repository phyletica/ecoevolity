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

#ifndef ECOEVOLITY_STATS_UTIL_HPP
#define ECOEVOLITY_STATS_UTIL_HPP

#include <vector>
#include <limits>

#include "assert.hpp"
#include "error.hpp"


template <typename T>
class SampleSummarizer {

    private:
        unsigned int n_ = 0;
        T min_ = std::numeric_limits<T>::max();
        T max_ = std::numeric_limits<T>::lowest();
        double mean_ = 0.0;
        double sum_devs_2_ = 0.0;
        double sum_devs_3_ = 0.0;
        double sum_devs_4_ = 0.0;

        void check_for_samples() const {
            if (this->n_ < 1) {
                throw EcoevolityError(
                    "Requesting a summary statistic from a SampleSummarizer "
                    "with no samples.");
            }
        }

    public:

        void add_sample(T x) {
            unsigned int n = this->n_ + 1;
            double d = ((double)x) - this->mean_;
            double d_n = d / n;
            double d_n2 = d_n * d_n;
            this->mean_ = this->mean_ + d_n;
            double first_term = d * d_n * this->n_;
            this->sum_devs_4_ += (first_term * d_n2 * ((n * n) - (3 * n) + 3)) +
                (6 * d_n2 * this->sum_devs_2_) - (4 * d_n * this->sum_devs_3_);
            this->sum_devs_3_ += (first_term * d_n * (n - 2)) -
                (3 * d_n * this->sum_devs_2_);
            this->sum_devs_2_ += first_term;
            this->n_ = n;
            if (x < this->min_) {
                this->min_ = x;
            }
            if (x > this->max_) {
                this->max_ = x;
            }
        }

        unsigned int sample_size() const {
            return this->n_;
        }

        T min() const {
            this->check_for_samples();
            return this->min_;
        }
        T max() const {
            this->check_for_samples();
            return this->max_;
        }
        T range() const {
            this->check_for_samples();
            return this->max_ - this->min_;
        }

        double mean() const {
            this->check_for_samples();
            return this->mean_;
        }

        double variance() const {
            this->check_for_samples();
            if (this->n_ == 1) {
                return std::numeric_limits<double>::infinity();
            }
            return (this->sum_devs_2_ / (this->n_ - 1));
        }

        double population_variance() const {
            this->check_for_samples();
            return (this->sum_devs_2_ / this->n_);
        }

        double std_dev() const {
            return std::sqrt(this->variance());
        }
        double std_error() const {
            return (this->std_dev() / std::sqrt(this->n_));
        }

        double skewness() const {
            this->check_for_samples();
            if (this->sum_devs_2_ == 0.0) {
                return std::numeric_limits<double>::quiet_NaN();
            }
            return ((this->sum_devs_3_ * std::sqrt(this->n_)) /
                    std::pow(this->sum_devs_2_, (3.0/2.0)));
        }
        double kurtosis() const {
            this->check_for_samples();
            if (this->sum_devs_2_ == 0.0) {
                return std::numeric_limits<double>::quiet_NaN();
            }
            return ((this->n_ * this->sum_devs_4_) /
                    std::pow(this->sum_devs_2_, 2.0));
        }
        double excess_kurtosis() const {
            return this->kurtosis() - 3.0;
        }

};

#endif
