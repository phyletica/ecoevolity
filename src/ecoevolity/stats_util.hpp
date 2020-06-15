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
#include <cmath>
#include <numeric>

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


/**
 * Calculate Monte Carlo standard error.
 * 
 * Adapted from 'mcse' function of the 'mcmcse' R package (Gnu GPL version 2).
 * See:
 * 
 * Flegal, J. M. and Jones, G. L. (2010) Batch means and spectral variance
 * estimators in Markov chain Monte Carlo. The Annals of Statistics,
 * 38:1034--1070.
 */
template <typename T>
inline std::pair<double, double> monte_carlo_standard_error(
        const std::vector<T> & samples) {
    unsigned int n = samples.size(); 
    unsigned int b = (unsigned int)std::floor(std::sqrt(n));
    unsigned int a = (unsigned int)std::floor((double)n / b);
    std::vector<double> y;
    for (unsigned int k = 1; k < a + 1; ++k) {
        unsigned int count = 0;
        double z_sum = 0.0;
        for (unsigned int i = ((k - 1) * b); i < (k * b); ++i) {
            z_sum += (double)samples.at(i);
            ++count;
        }
        y.push_back(z_sum / (double)count);
    }
    double s_sum = std::accumulate(samples.begin(), samples.end(), 0.0);
    double mu_hat = s_sum / (double)n;
    std::vector<double> sq_diffs(y.size());
    for (unsigned int i = 0; i < y.size(); ++i) {
        sq_diffs.at(i) = std::pow(y.at(i) - mu_hat, 2);
    }
    double sum_sq_diffs = std::accumulate(sq_diffs.begin(), sq_diffs.end(), 0.0);
    double var_hat = (double)b * sum_sq_diffs / ((double)a - 1.0);
    double se = std::sqrt(var_hat / (double)n);
    return std::make_pair(mu_hat, se);
}

/**
 * Estimate effective sample size of MCMC sample.
 * 
 * Adapted from 'ess' function of the 'mcmcse' R package (Gnu GPL version 2).
 * See:
 * 
 * Gong, Lei, and James M. Flegal. A practical sequential stopping rule for
 * high-dimensional MCMC and its application to spatial-temporal Bayesian
 * models. arXiv:1403.5536v1 [stat.CO].
 */
template <typename T>
inline double effective_sample_size(
        const std::vector<T> & samples,
        const bool limit_to_number_of_samples = true) {
    SampleSummarizer<T> ss;
    for (auto x : samples) {
        ss.add_sample(x);
    }
    std::pair<double, double> mu_sigma = monte_carlo_standard_error(samples);
    double sigma = mu_sigma.second;
    if ((ss.sample_size() == 0) or (sigma == 0.0)) {
        return 0.0;
    }
    double n = (double)ss.sample_size();
    double ess = ((n * ss.variance()) / (std::pow(sigma, 2) * n));
    if ((ess > n) and (limit_to_number_of_samples)) {
        return n;
    }
    return ess;
}

/**
 * Calculate the potential scale reduction factor.
 *
 * Returns the square root of Equation 1.1 in:
 *
 * Brooks, Stephen P. and Andrew Gelman. 1998. General Methods for Monitoring
 * Convergence of Iterative Simulations. Journal of Computational and
 * Graphical Statistics, Volume7, Number 4, Pages 434-455.
 */
template <typename T>
inline double potential_scale_reduction_factor(
        const std::vector< std::vector<T> > & chains) {
    unsigned int nchains = chains.size();
    ECOEVOLITY_ASSERT(nchains > 1);
    unsigned int nsamples = chains.at(0).size();
    std::vector<SampleSummarizer<T>> sample_summaries(nchains);
    for (unsigned int chain_idx = 0; chain_idx < chains.size(); ++chain_idx) {
        for (unsigned int sample_idx = 0;
                sample_idx < chains.at(chain_idx).size();
                ++sample_idx) {
            sample_summaries.at(chain_idx).add_sample(
                    chains.at(chain_idx).at(sample_idx));
        }
        ECOEVOLITY_ASSERT(sample_summaries.at(chain_idx).sample_size() == nsamples);
    }
    ECOEVOLITY_ASSERT(sample_summaries.size() == nchains);

    SampleSummarizer<double> summary_of_variances;
    SampleSummarizer<double> summary_of_means;
    for (auto ss : sample_summaries) {
        summary_of_variances.add_sample(ss.variance());
        summary_of_means.add_sample(ss.mean());
    }
    double within_chain_var = summary_of_variances.mean();
    double between_chain_var = summary_of_means.variance();
    double pooled_var_term1 = (1.0 - (1.0 / (double)nsamples)) * within_chain_var;
    double pooled_var = pooled_var_term1 + between_chain_var;
    double pooled_posterior_var = pooled_var + (between_chain_var / nchains);
    if (within_chain_var == 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    return std::sqrt(pooled_posterior_var / within_chain_var);
}


/**
 * Return the median value from a list of samples.
 */
template <typename T>
inline double get_median(const std::vector<T> & samples) {
    std::vector<T> s = samples;
    std::sort(s.begin(), s.end());
    unsigned int n = s.size();
    if (n < 1) {
        throw EcoevolityError("stats_util::median: empty samples");
    }
    if (n % 2 == 0) {
        return (s.at((n / 2) - 1) + s.at(n / 2)) / 2.0;
    }
    return (double)s.at((n - 1) / 2);
}

/**
 * Return tuple of highest posterior density (HPD).
 *
 * The interval is estimated via a sliding window to find the narrowest
 * interval that contains the specified proportion of samples.
 */
template <typename T>
inline std::pair<T, T> get_hpd_interval(
        const std::vector<T> & samples,
        const double interval_prob = 0.95) {
    if (samples.size() < 1){
        throw EcoevolityError("stats_util::get_hpd_interval: empty samples");
    }
    if (interval_prob <= 0.0) {
        throw EcoevolityError("stats_util::get_hpd_interval: prob interval must be > 0");
    }
    if (interval_prob > 1.0) {
        throw EcoevolityError("stats_util::get_hpd_interval: prob interval greater than 1");
    }
    std::vector<T> s = samples;
    std::sort(s.begin(), s.end());
    double tail_prob = 1 - interval_prob;
    unsigned int n = s.size();
    int num_exclude = (int)round((n * tail_prob) - 0.5);
    if (num_exclude < 1) {
        num_exclude = 0;
    }
    std::vector<T> widths;
    // sliding window to find possible interval widths
    for (int i = 0; i < num_exclude + 1; ++i) {
        T lower = s.at(i);
        T upper = s.at((n - 1) - num_exclude + i);
        widths.push_back(upper - lower);
    }
    typename std::vector<T>::iterator min_width_iter = std::min_element(
            std::begin(widths), std::end(widths));
    int min_index = std::distance(std::begin(widths), min_width_iter);
    int max_index = (n - 1) - num_exclude + min_index;
    return std::make_pair(s.at(min_index), s.at(max_index));
}

/**
 * Return quantile associated with probability `p`.
 *
 * Modified from code by Wai Yip Tung, licensed under PSF license and available
 * here:
 * http://code.activestate.com/recipes/511478-finding-the-percentile-of-the-values/
 */
template <typename T>
inline double quantile(
        const std::vector<T> & samples,
        const double p) {
    if (samples.size() < 1){
        throw EcoevolityError("stats_util::quantile: empty samples");
    }
    if (p < 0.0) {
        throw EcoevolityError("stats_util::quantile: p < 0");
    }
    if (p > 1.0) {
        throw EcoevolityError("stats_util::quantile: p > 1");
    }
    std::vector<T> s = samples;
    std::sort(s.begin(), s.end());
    double k = (s.size() - 1) * p;
    double f = std::floor(k);
    double c = std::ceil(k);
    if (f == c) {
        return (double)s.at((int)(k));
    }
    double d0 = s.at((int)f) * (c - k);
    double d1 = s.at((int)c) * (k - f);
    return d0 + d1;
}

/**
 * Return tuple of interval of 2.5% and 97.5% quantiles.
 */
template <typename T>
inline std::pair<double, double> quantiles_95(
        const std::vector<T> & samples) {
    return std::make_pair(quantile(samples, 0.025), quantile(samples, 0.975));
}

/**
* Return the percentile of value `v` (from 0.0 to 1.0) relative to a collection
* of samples.
*
* If the value is less than all the samples, the percentile is 0.0. If the
* value is greater than all the samples the percentile is 1.0. If the value is
* greater than or equal to 1 sample out of 10, the percentile is 0.1.
*
* The samples do not have to be sorted (they are sorted internally before the
* percentile is determined).
*/
template <typename T>
inline double percentile(
        const std::vector<T> & samples,
        const T v) {
    if (samples.size() < 1){
        throw EcoevolityError("stats_util::percentile: empty samples");
    }
    std::vector<T> s = samples;
    std::sort(s.begin(), s.end());
    unsigned int n = s.size();
    for (unsigned int i = 0; i < n; ++i){
        if (v < s.at(i)) {
            return i / (double)n;
        }
    }
    return 1.0;
}


template <typename T>
class SampleSummary {
    private:
        unsigned int n_ = 0;
        double mean_;
        double median_;
        double variance_;
        T min_;
        T max_;
        std::pair<T, T> hpdi_95_;
        std::pair<double, double> qi_95_;

        void check_for_samples() const {
            if (this->n_ < 1) {
                throw EcoevolityError(
                    "Requesting a summary statistic from a SampleSummary "
                    "with no samples.");
            }
        }

    public:
        SampleSummary(const std::vector<T> & samples) {
            SampleSummarizer<T> ss;
            for (auto x : samples) {
                ss.add_sample(x);
            }
            this->n_ = ss.sample_size();
            this->mean_ = ss.mean();
            this->median_ = get_median<T>(samples);
            this->variance_ = ss.variance();
            this->min_ = ss.min();
            this->max_ = ss.max();
            this->hpdi_95_ = get_hpd_interval<T>(samples, 0.95);
            this->qi_95_ = quantiles_95<T>(samples);
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

        double mean() const {
            this->check_for_samples();
            return this->mean_;
        }

        double median() const {
            this->check_for_samples();
            return this->median_;
        }

        double variance() const {
            this->check_for_samples();
            return this->variance_;
        }

        double std_dev() const {
            return std::sqrt(this->variance_);
        }

        double std_error() const {
            return (this->std_dev() / std::sqrt(this->n_));
        }

        const std::pair<T, T> & hpdi_95() const {
            return this->hpdi_95_;
        }

        const std::pair<double, double> & qi_95() const {
            return this->qi_95_;
        }
};

#endif
