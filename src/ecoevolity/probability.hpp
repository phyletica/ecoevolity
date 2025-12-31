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

#ifndef ECOEVOLITY_PROBABILITY_HPP
#define ECOEVOLITY_PROBABILITY_HPP

#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>

#include "error.hpp"
#include "assert.hpp"
#include "rng.hpp"


class ContinuousProbabilityDistribution {
    public:
        ContinuousProbabilityDistribution() { }
        virtual ~ContinuousProbabilityDistribution() { }
        ContinuousProbabilityDistribution& operator=(const ContinuousProbabilityDistribution& other) {
            return * this;
        }

        static double ln_gamma_function(double x) {
            return std::lgamma(x);
        }

        virtual double ln_pdf(double x) const = 0;
        virtual double relative_ln_pdf(double x) const = 0;
        virtual double get_mean() const = 0;
        virtual double get_variance() const = 0;

        virtual double get_min() const {
            return -std::numeric_limits<double>::infinity();
        }
        virtual double get_max() const {
            return std::numeric_limits<double>::infinity();
        }

        virtual std::string get_name() const = 0;
        virtual std::string to_string() const = 0;

        virtual double draw(RandomNumberGenerator & rng) const = 0;
        
        virtual bool is_within_support(double x) const = 0;
};

class ImproperUniformDistribution : public ContinuousProbabilityDistribution {
    public:
        ImproperUniformDistribution() { }
        ~ImproperUniformDistribution() { }
        ImproperUniformDistribution& operator=(const ImproperUniformDistribution& other) {
            return * this;
        }

        std::string get_name() const {
            return "uniform(-inf, +inf)";
        }
        std::string to_string() const {
            return this->get_name();
        }

        double relative_ln_pdf(double x) const {
            return 0.0;
        }
        double ln_pdf(double x) const {
            throw EcoevolityProbabilityDistributionError(
                    "The density is undefinded for ImproperUniformDistribution");
        }
        double get_mean() const {
            throw EcoevolityProbabilityDistributionError(
                    "The mean is undefinded for ImproperUniformDistribution");
        }
        double get_variance() const {
            throw EcoevolityProbabilityDistributionError(
                    "The variance is undefinded for ImproperUniformDistribution");
        }

        double draw(RandomNumberGenerator & rng) const {
            return rng.uniform_real(
                    std::numeric_limits<double>::lowest(),
                    std::numeric_limits<double>::max());
        }

        bool is_within_support(double x) const {
            return true;
        }
};

class ImproperPositiveUniformDistribution: public ImproperUniformDistribution {
    public:
        ImproperPositiveUniformDistribution() { }
        ~ImproperPositiveUniformDistribution() { }
        ImproperPositiveUniformDistribution& operator=(const ImproperPositiveUniformDistribution& other) {
            return * this;
        }

        std::string get_name() const {
            return "uniform(0, +inf)";
        }
        std::string to_string() const {
            return this->get_name();
        }

        double get_min() const {
            return 0.0;
        }

        bool is_within_support(double x) const {
            if (x < 0.0) {
                return false;
            }
            return true;
        }

        double relative_ln_pdf(double x) const {
            if (! this->is_within_support(x)) {
                return -std::numeric_limits<double>::infinity();
            }
            return 0.0;
        }

        double draw(RandomNumberGenerator & rng) const {
            return rng.uniform_real(
                    0.0,
                    std::numeric_limits<double>::max());
        }
};

class UniformDistribution : public ContinuousProbabilityDistribution {
    protected:
        double min_ = 0.0;
        double max_ = 1.0;
        double ln_density_ = 0.0;

    public:
        UniformDistribution() { }
        ~UniformDistribution() { }
        UniformDistribution(double a, double b) {
            if (b <= a) {
                throw EcoevolityProbabilityDistributionError(
                        "The upper limit must be greater than lower limit for UniformDistribution");
            }
            this->min_ = a;
            this->max_ = b;
            this->ln_density_ = -1.0 * std::log(b - a);
        }
        UniformDistribution& operator=(const UniformDistribution& other) {
            this->min_ = other.min_;
            this->max_ = other.max_;
            this->ln_density_ = other.ln_density_;
            return * this;
        }

        double relative_ln_pdf(double x) const {
            return this->ln_pdf(x);
        }

        bool is_within_support(double x) const {
            if ((x < this->min_) || (x > this->max_)) {
                return false;
            }
            return true;
        }

        double ln_pdf(double x) const {
            if (! this->is_within_support(x)) {
                return -std::numeric_limits<double>::infinity();
            }
            return this->ln_density_;
        }

        double get_mean() const {
            return ((this->min_ + this->max_) / 2.0);
        }

        double get_variance() const {
            return ((this->max_ - this->min_) * (this->max_ - this->min_) / 12.0);
        }

        double get_min() const {
            return this->min_;
        }
        double get_max() const {
            return this->max_;
        }

        std::string get_name() const {
            return "uniform";
        }
        std::string to_string() const {
            std::ostringstream ss;
            ss << this->get_name() << "(" << this->min_ << ", " << this->max_ << ")";
            return ss.str();
        }

        double draw(RandomNumberGenerator & rng) const {
            return rng.uniform_real(this->min_, this->max_);
        }
};

class BetaDistribution: public ContinuousProbabilityDistribution {
    protected:
        double min_ = 0.0;
        double max_ = 1.0;
        double alpha_ = 1.0;
        double beta_ = 1.0;

    public:
        BetaDistribution() { }
        ~BetaDistribution() { }
        BetaDistribution(double alpha, double beta) {
            if ((alpha <= 0.0) || (beta <= 0.0)) {
                throw EcoevolityProbabilityDistributionError(
                        "alpha and beta must be greater than zero for beta distribution");
            }
            this->alpha_ = alpha;
            this->beta_ = beta;
        }
        BetaDistribution& operator=(const BetaDistribution& other) {
            this->min_ = other.min_;
            this->max_ = other.max_;
            this->alpha_ = other.alpha_;
            this->beta_ = other.beta_;
            return * this;
        }

        bool is_within_support(double x) const {
            if ((x <= this->min_) || (x >= this->max_)) {
                return false;
            }
            return true;
        }

        static double get_ln_pdf(double x, double alpha, double beta) {
            if ((x <= 0.0) || (x >= 1.0)) {
                return -std::numeric_limits<double>::infinity();
            }
		    double lnp = ((alpha - 1.0) * std::log(x)) + ((beta - 1.0) * std::log(1.0 - x));

		    lnp += BetaDistribution::ln_gamma_function(alpha + beta);
		    lnp -= BetaDistribution::ln_gamma_function(alpha);
		    lnp -= BetaDistribution::ln_gamma_function(beta);

		    return lnp;
        }

        double ln_pdf(double x) const {
            return BetaDistribution::get_ln_pdf(x, this->alpha_, this->beta_);
        }

        static double get_scaled_ln_pdf(double x,
                double alpha,
                double beta,
                double scale_parameter) {
            double unscaled_x = x / scale_parameter;
            double lnp = BetaDistribution::get_ln_pdf(unscaled_x, alpha, beta);
            lnp -= std::log(scale_parameter);
		    return lnp;
        }

        double scaled_ln_pdf(double x, double scale_parameter) const {
            return BetaDistribution::get_scaled_ln_pdf(x,
                    this->alpha_,
                    this->beta_,
                    scale_parameter);
        }

        double relative_ln_pdf(double x) const {
            return this->ln_pdf(x);
        }

        double get_mean() const {
            return (this->alpha_ / (this->alpha_ + this->beta_));
        }
        double get_variance() const {
            const double ab = this->alpha_ + this->beta_;
            return (this->alpha_ * this->beta_) / (ab * ab * (ab + 1.0));
        }

        double get_min() const {
            return this->min_;
        }

        double get_max() const {
            return this->max_;
        }

        double get_alpha() const {
            return this->alpha_;
        }
        double get_beta() const {
            return this->beta_;
        }

        std::string get_name() const {
            return "beta";
        }
        std::string to_string() const {
            std::ostringstream ss;
            ss << this->get_name() << "(" << this->alpha_ << ", " << this->beta_ << ")";
            return ss.str();
        }

        static double get_draw(RandomNumberGenerator & rng,
                double alpha,
                double beta) {
            return rng.beta(alpha, beta);
        }
        double draw(RandomNumberGenerator & rng) const {
            return BetaDistribution::get_draw(rng,
                    this->alpha_,
                    this->beta_);
        }
};

class OffsetGammaDistribution : public ContinuousProbabilityDistribution {
    protected:
        double min_ = 0.0;
        double shape_ = 1.0;
        double scale_ = 1.0;
        double ln_constant_ = 0.0;

        void compute_ln_constant() {
            double ln_gamma = this->ln_gamma_function(this->shape_);
            this->ln_constant_ = -1.0*(this->shape_ * std::log(this->scale_) + ln_gamma);
        }

    public:
        OffsetGammaDistribution() { }
        ~OffsetGammaDistribution() { }
        OffsetGammaDistribution(double shape, double scale, double offset) {
            if ((shape <= 0.0) || (scale <= 0.0)) {
                throw EcoevolityProbabilityDistributionError(
                        "Shape and scale must be greater than zero for gamma distribution");
            }
            this->shape_ = shape;
            this->scale_ = scale;
            this->min_ = offset;
            this->compute_ln_constant();
        }
        OffsetGammaDistribution& operator=(const OffsetGammaDistribution& other) {
            this->min_ = other.min_;
            this->shape_ = other.shape_;
            this->scale_ = other.scale_;
            this->ln_constant_ = other.ln_constant_;
            return * this;
        }

        bool is_within_support(double x) const {
            if ((x < this->min_) || ((x == this->min_) && (this->shape_ > 1.0))) {
                return false;
            }
            return true;
        }

        double ln_pdf(double x) const {
            if (! this->is_within_support(x)) {
                return -std::numeric_limits<double>::infinity();
            }
            x -= this->min_;
            double term1 = (this->shape_ - 1.0);
            // corner case when x=0 and shape=1
            // the pdf should be 1.0, and so ln_pdf is 0.0
            // avoiding std::log(0) which will return nan.
            if (term1 != 0.0) {
                term1 *= std::log(x);
            }
            double term2 = -x / this->scale_;
            return (term1 + term2 + this->ln_constant_);
        }

        double relative_ln_pdf(double x) const {
            return this->ln_pdf(x);
        }

        double get_mean() const {
            return (this->shape_ * this->scale_) + this->min_;
        }
        double get_variance() const {
            return this->shape_ * this->scale_ * this->scale_;
        }

        double get_min() const {
            return this->min_;
        }
        double get_offset() const {
            return this->get_min();
        }

        double get_shape() const {
            return this->shape_;
        }
        double get_scale() const {
            return this->scale_;
        }

        std::string get_name() const {
            return "gamma";
        }
        std::string to_string() const {
            std::ostringstream ss;
            ss << this->get_name() << "(shape = " << this->shape_ << ", scale = " << this->scale_ << ", offset = " << this->min_ << ")";
            return ss.str();
        }

        double draw(RandomNumberGenerator & rng) const {
            return this->min_ + rng.gamma(this->shape_, this->scale_);
        }
};

class GammaDistribution : public OffsetGammaDistribution {
    public:
        GammaDistribution() : OffsetGammaDistribution() { }
        ~GammaDistribution() { }
        GammaDistribution(double shape, double scale) : OffsetGammaDistribution(shape, scale, 0.0) { }

        GammaDistribution& operator=(const GammaDistribution& other) {
            this->min_ = other.min_;
            this->shape_ = other.shape_;
            this->scale_ = other.scale_;
            this->ln_constant_ = other.ln_constant_;
            return * this;
        }

        std::string to_string() const {
            std::ostringstream ss;
            ss << this->get_name() << "(shape = " << this->shape_ << ", scale = " << this->scale_ << ")";
            return ss.str();
        }
};

class OffsetExponentialDistribution : public OffsetGammaDistribution {
    public:
        OffsetExponentialDistribution() : OffsetGammaDistribution() { }
        ~OffsetExponentialDistribution() { }
        OffsetExponentialDistribution(double lambda, double offset)
                : OffsetGammaDistribution(1.0, 1.0/lambda, offset) {
            if (lambda <= 0.0) {
                throw EcoevolityProbabilityDistributionError(
                        "lambda must be greater than 0 for exponential distribution");
            }
        }
        OffsetExponentialDistribution& operator=(const OffsetExponentialDistribution& other) {
            this->min_ = other.min_;
            this->shape_ = other.shape_;
            this->scale_ = other.scale_;
            this->ln_constant_ = other.ln_constant_;
            return * this;
        }

        double get_lambda() const {
            return 1.0 / this->scale_;
        }

        std::string get_name() const {
            return "exp";
        }

        std::string to_string() const {
            std::ostringstream ss;
            ss << this->get_name() << "(lambda = " << this->get_lambda() << ", offset = " << this->min_ << ")";
            return ss.str();
        }
};

class ExponentialDistribution: public OffsetExponentialDistribution {
    public:
        ExponentialDistribution() : OffsetExponentialDistribution() { }
        ~ExponentialDistribution() { }
        ExponentialDistribution(double lambda)
                : OffsetExponentialDistribution(lambda, 0.0) { }

        ExponentialDistribution& operator=(const ExponentialDistribution& other) {
            this->min_ = other.min_;
            this->shape_ = other.shape_;
            this->scale_ = other.scale_;
            this->ln_constant_ = other.ln_constant_;
            return * this;
        }

        std::string to_string() const {
            std::ostringstream ss;
            ss << this->get_name() << "(lambda = " << this->get_lambda() << ")";
            return ss.str();
        }
};

class DirichletDistribution {
    protected:
        double min_ = 0.0;
        double max_ = 1.0;
        std::vector<double> parameters_ = {1.0, 1.0, 1.0};

    public:
        DirichletDistribution() { }
        ~DirichletDistribution() { }
        DirichletDistribution(const std::vector<double> & parameters) {
            this->parameters_.clear();
            this->parameters_.reserve(parameters.size());
            for (auto p : parameters) {
                if (p <= 0.0) {
                    throw EcoevolityProbabilityDistributionError(
                            "all parameters must be greater than zero for dirichlet distribution");
                }
                this->parameters_.push_back(p);
            }
        }

        DirichletDistribution& operator=(const DirichletDistribution& other) {
            this->min_ = other.min_;
            this->max_ = other.max_;
            this->parameters_ = other.parameters_;
            return * this;
        }

        double ln_gamma_function(double x) const {
            return std::lgamma(x);
        }

        bool is_within_support(double x) const {
            if ((x <= this->min_) || (x >= this->max_)) {
                return false;
            }
            return true;
        }

        bool is_within_support(const std::vector<double> & values) const {
            for (const auto & x : values) {
                if (! this->is_within_support(x)) {
                    return false;
                }
            }
            return true;
        }

        double ln_pdf(const std::vector<double> & x) const {
            ECOEVOLITY_ASSERT(x.size() == this->parameters_.size());
            double r = 0.0;
            double sum_p = 0.0;
            double sum_x = 0.0;
            for (unsigned int i = 0; i < x.size(); ++i) {
                if (! this->is_within_support(x.at(i))) {
                    return -std::numeric_limits<double>::infinity();
                }
                double p_i = this->parameters_.at(i);
                sum_p += p_i;
                r += (p_i - 1.0) * std::log(x.at(i));
                r -= this->ln_gamma_function(p_i);
                sum_x += x.at(i);
            }
            ECOEVOLITY_ASSERT_APPROX_EQUAL(sum_x, 1.0);
            r += this->ln_gamma_function(sum_p);
            return r;
        }

        double relative_ln_pdf(const std::vector<double> & x) const {
            ECOEVOLITY_ASSERT(x.size() == this->parameters_.size());
            double r = 0.0;
            double sum_x = 0.0;
            for (unsigned int i = 0; i < x.size(); ++i) {
                if (! this->is_within_support(x.at(i))) {
                    return -std::numeric_limits<double>::infinity();
                }
                r += (this->parameters_.at(i) - 1.0) * std::log(x.at(i));
                sum_x += x.at(i);
            }
            ECOEVOLITY_ASSERT_APPROX_EQUAL(sum_x, 1.0);
            return r;
        }

        std::vector<double> get_mean() const {
            std::vector<double> means;
            means.reserve(this->parameters_.size());
            double p = this->get_parameter_sum();
            for (auto p_i : this->parameters_) {
                means.push_back(p_i / p);
            }
            return means;
        }

        double get_mean(unsigned int i) const {
            double p = this->get_parameter_sum();
            return this->parameters_.at(i) / p;
        }

        std::vector<double> get_variance() const {
            std::vector<double> variances;
            variances.reserve(this->parameters_.size());
            double p = this->get_parameter_sum();
            double denom = p * p * (p + 1.0);
            for (auto p_i : this->parameters_) {
                variances.push_back((p_i * (p - p_i)) / denom);
            }
            return variances;
        }

        double get_variance(unsigned int i) const {
            double p = this->get_parameter_sum();
            double denom = p * p * (p + 1.0);
            return (this->parameters_.at(i) * (p - this->parameters_.at(i))) / denom;
        }

        double get_parameter_sum() const{
            double s = 0.0;
            for (auto p : this->parameters_) {
                s += p;
            }
            return s;
        }

        double get_min() const {
            return this->min_;
        }

        double get_max() const {
            return this->max_;
        }

        const std::vector<double>& get_parameters() const {
            return this->parameters_;
        }

        std::string get_name() const {
            return "dirichlet";
        }
        std::string to_string() const {
            std::ostringstream ss;
            ss << this->get_name() << "(" << this->parameters_.at(0);
            for (unsigned int i = 1; i < this->parameters_.size(); ++i) {
                ss << ", " << this->parameters_.at(i);
            }
            ss << ")";
            return ss.str();
        }

        std::vector<double> draw(RandomNumberGenerator & rng) const {
            return rng.dirichlet(this->parameters_);
        }
};

#endif
