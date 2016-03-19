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


class ContinuousProbabilityDist {
    public:
        ContinuousProbabilityDist() { }
        virtual ~ContinuousProbabilityDist() { }

        double ln_gamma_function(double x) const {
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
};

class RandomValueGenerator {
    public:
        RandomNumberGenerator * rng_ = NULL;

        RandomValueGenerator() {
            this->rng_ = new RandomNumberGenerator();
        }
        RandomValueGenerator(long seed) {
            this->rng_ = new RandomNumberGenerator(seed);
        }
        RandomValueGenerator(RandomNumberGenerator * rng) {
            this->rng_ = rng;
        }
        virtual ~RandomValueGenerator() {
            delete this->rng_;
        }
        RandomValueGenerator& operator=(const RandomValueGenerator& other) {
            this->rng_ = other.rng_;
            return * this;
        }

        virtual double draw() const = 0;
};

        

class ImproperUniformDist : public ContinuousProbabilityDist {
    public:
        ImproperUniformDist() { }
        ~ImproperUniformDist() { }

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
            throw EcoevolityProbabilityDistError(
                    "The density is undefinded for ImproperUniformDist");
        }
        double get_mean() const {
            throw EcoevolityProbabilityDistError(
                    "The mean is undefinded for ImproperUniformDist");
        }
        double get_variance() const {
            throw EcoevolityProbabilityDistError(
                    "The variance is undefinded for ImproperUniformDist");
        }
};

class ImproperPositiveUniformDist: public ImproperUniformDist {
    public:
        ImproperPositiveUniformDist() { }
        ~ImproperPositiveUniformDist() { }

        std::string get_name() const {
            return "uniform(0, +inf)";
        }
        std::string to_string() const {
            return this->get_name();
        }

        double get_min() const {
            return 0.0;
        }

        double relative_ln_pdf(double x) const {
            if (x < 0.0) {
                return -std::numeric_limits<double>::infinity();
            }
            return 0.0;
        }
};

class UniformDist : public ContinuousProbabilityDist {
    protected:
        double min_ = 0.0;
        double max_ = 1.0;
        double ln_density_ = 0.0;

    public:
        UniformDist() { }
        ~UniformDist() { }
        UniformDist(double a, double b) {
            if (b <= a) {
                throw EcoevolityProbabilityDistError(
                        "The upper limit must be greater than lower limit for UniformDist");
            }
            this->min_ = a;
            this->max_ = b;
            this->ln_density_ = -1.0 * std::log(b - a);
        }
        UniformDist& operator=(const UniformDist& other) {
            this->min_ = other.min_;
            this->max_ = other.max_;
            this->ln_density_ = other.ln_density_;
            return * this;
        }

        double relative_ln_pdf(double x) const {
            return this->ln_pdf(x);
        }

        double ln_pdf(double x) const {
            if ((x < this->min_) || (x > this->max_)) {
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
};

class UniformDistribution : public UniformDist, public RandomValueGenerator {
    public:
        UniformDistribution() : UniformDist(), RandomValueGenerator() { }
        UniformDistribution(double a, double b) : UniformDist(a, b), RandomValueGenerator() { }
        UniformDistribution(double a, double b, long seed) : UniformDist(a, b), RandomValueGenerator(seed) { }
        UniformDistribution(long seed) : UniformDist(), RandomValueGenerator(seed) { }
        UniformDistribution(double a, double b, RandomNumberGenerator * rng) : UniformDist(a, b), RandomValueGenerator(rng) { }
        UniformDistribution(RandomNumberGenerator * rng) : UniformDist(), RandomValueGenerator(rng) { }
        UniformDistribution& operator=(const UniformDistribution& other) {
            this->rng_ = other.rng_;
            this->min_ = other.min_;
            this->max_ = other.max_;
            this->ln_density_ = other.ln_density_;
            return * this;
        }

        double draw() const {
            return this->rng_->uniform_real(this->min_, this->max_);
        }
};

class OffsetGammaDist : public ContinuousProbabilityDist {
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
        OffsetGammaDist() { }
        ~OffsetGammaDist() { }
        OffsetGammaDist(double shape, double scale, double offset) {
            if ((shape <= 0.0) || (scale <= 0.0)) {
                throw EcoevolityProbabilityDistError(
                        "Shape and scale must be greater than zero for gamma distribution");
            }
            this->shape_ = shape;
            this->scale_ = scale;
            this->min_ = offset;
            this->compute_ln_constant();
        }
        OffsetGammaDist& operator=(const OffsetGammaDist& other) {
            this->min_ = other.min_;
            this->shape_ = other.shape_;
            this->scale_ = other.scale_;
            this->ln_constant_ = other.ln_constant_;
            return * this;
        }

        double ln_pdf(double x) const {
            if ((x < this->min_) || ((x == this->min_) && (this->shape_ > 1.0))) {
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
};

class OffsetGammaDistribution : public OffsetGammaDist, public RandomValueGenerator {
    public:
        OffsetGammaDistribution() : OffsetGammaDist(), RandomValueGenerator() { }
        OffsetGammaDistribution(double shape, double scale, double offset)
                : OffsetGammaDist(shape, scale, offset),
                  RandomValueGenerator() { }
        OffsetGammaDistribution(double shape, double scale, double offset,
                long seed)
                : OffsetGammaDist(shape, scale, offset),
                  RandomValueGenerator(seed) { }
        OffsetGammaDistribution(long seed)
                : OffsetGammaDist(),
                  RandomValueGenerator(seed) { }
        OffsetGammaDistribution(double shape, double scale, double offset,
                RandomNumberGenerator * rng)
                : OffsetGammaDist(shape, scale, offset),
                  RandomValueGenerator(rng) { }
        OffsetGammaDistribution(RandomNumberGenerator * rng)
                : OffsetGammaDist(),
                  RandomValueGenerator(rng) { }
        OffsetGammaDistribution& operator=(const OffsetGammaDistribution& other) {
            this->rng_ = other.rng_;
            this->min_ = other.min_;
            this->shape_ = other.shape_;
            this->scale_ = other.scale_;
            this->ln_constant_ = other.ln_constant_;
            return * this;
        }

        double draw() const {
            return this->min_ + this->rng_->gamma(this->shape_, this->scale_);
        }
};

class GammaDist : public OffsetGammaDist {
    public:
        GammaDist() : OffsetGammaDist() { }
        ~GammaDist() { }
        GammaDist(double shape, double scale) : OffsetGammaDist(shape, scale, 0.0) { }

        GammaDist& operator=(const GammaDist& other) {
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

class GammaDistribution : public OffsetGammaDistribution {
    public:
        GammaDistribution() : OffsetGammaDistribution() { }
        GammaDistribution(double shape, double scale)
                : OffsetGammaDistribution(shape, scale, 0.0)
                { }
        GammaDistribution(double shape, double scale, long seed)
                : OffsetGammaDistribution(shape, scale, 0.0, seed)
                { }
        GammaDistribution(long seed)
                : OffsetGammaDistribution(seed)
                { }
        GammaDistribution(double shape, double scale,
                RandomNumberGenerator * rng)
                : OffsetGammaDistribution(shape, scale, 0.0, rng)
                { }
        GammaDistribution(RandomNumberGenerator * rng)
                : OffsetGammaDistribution(rng)
                { }
        GammaDistribution& operator=(const GammaDistribution& other) {
            this->rng_ = other.rng_;
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

class OffsetExponentialDist : public OffsetGammaDist {
    public:
        OffsetExponentialDist() : OffsetGammaDist() { }
        ~OffsetExponentialDist() { }
        OffsetExponentialDist(double lambda, double offset)
                : OffsetGammaDist(1.0, 1.0/lambda, offset) {
            if (lambda <= 0.0) {
                throw EcoevolityProbabilityDistError(
                        "lambda must be greater than 0 for exponential distribution");
            }
        }
        OffsetExponentialDist& operator=(const OffsetExponentialDist& other) {
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

class OffsetExponentialDistribution : public OffsetGammaDistribution {
    public:
        OffsetExponentialDistribution() : OffsetGammaDistribution() { }
        OffsetExponentialDistribution(double lambda, double offset)
                : OffsetGammaDistribution(1.0, 1.0/lambda, offset) {
            if (lambda <= 0.0) {
                throw EcoevolityProbabilityDistError(
                        "lambda must be greater than 0 for exponential distribution");
            }
        }
        OffsetExponentialDistribution(double lambda, double offset, long seed)
                : OffsetGammaDistribution(1.0, 1.0/lambda, offset, seed) {
            if (lambda <= 0.0) {
                throw EcoevolityProbabilityDistError(
                        "lambda must be greater than 0 for exponential distribution");
            }
        }
        OffsetExponentialDistribution(long seed)
                : OffsetGammaDistribution(seed)
                { }
        OffsetExponentialDistribution(double lambda, double offset,
                RandomNumberGenerator * rng)
                : OffsetGammaDistribution(1.0, 1.0/lambda, offset, rng) {
            if (lambda <= 0.0) {
                throw EcoevolityProbabilityDistError(
                        "lambda must be greater than 0 for exponential distribution");
            }
        }
        OffsetExponentialDistribution(RandomNumberGenerator * rng)
                : OffsetGammaDistribution(rng)
                { }
        OffsetExponentialDistribution& operator=(const OffsetExponentialDistribution& other) {
            this->rng_ = other.rng_;
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

class ExponentialDist: public OffsetExponentialDist {
    public:
        ExponentialDist() : OffsetExponentialDist() { }
        ~ExponentialDist() { }
        ExponentialDist(double lambda)
                : OffsetExponentialDist(lambda, 0.0) { }

        ExponentialDist& operator=(const ExponentialDist& other) {
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

class ExponentialDistribution : public OffsetExponentialDistribution {
    public:
        ExponentialDistribution() : OffsetExponentialDistribution() { }
        ExponentialDistribution(double lambda)
                : OffsetExponentialDistribution(lambda, 0.0)
                { }
        ExponentialDistribution(double lambda, long seed)
                : OffsetExponentialDistribution(lambda, 0.0, seed)
                { }
        ExponentialDistribution(long seed)
                : OffsetExponentialDistribution(seed)
                { }
        ExponentialDistribution(double lambda,
                RandomNumberGenerator * rng)
                : OffsetExponentialDistribution(lambda, 0.0, rng)
                { }
        ExponentialDistribution(RandomNumberGenerator * rng)
                : OffsetExponentialDistribution(rng)
                { }
        ExponentialDistribution& operator=(const ExponentialDistribution& other) {
            this->rng_ = other.rng_;
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

#endif

