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


class ContinuousProbabilityDistribution {
    protected:
        std::string name_;
        double min_ = -std::numeric_limits<double>::infinity();
        double max_ = +std::numeric_limits<double>::infinity();

    public:
        ContinuousProbabilityDistribution() { }
        virtual ~ContinuousProbabilityDistribution() { }

        double ln_gamma_function(double x) const {
            return std::lgamma(x);
        }

        virtual double ln_pdf(double x) const = 0;
        virtual double relative_ln_pdf(double x) const = 0;
        virtual double get_mean() const = 0;
        virtual double get_variance() const = 0;

        double get_min() const {
            return this->min_;
        }
        double get_max() const {
            return this->max_;
        }

        std::string get_name() const {
            return this->name_;
        }
        virtual std::string to_string() const = 0;
};

class ImproperUniformDistribution : public ContinuousProbabilityDistribution {
    protected:
        std::string name_ = "uniform(0, +inf)"
        this->min_ = 0.0;

    public:
        ImproperUniformDistribution() { }
        ~ImproperUniformDistribution() { }

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
};

class UniformDistribution : public ContinuousProbabilityDistribution {
    protected:
        std::string name_ = "uniform";
        double min_ = 0.0;
        double max_ = 1.0;
        double ln_density_ = 0.0;

    public:
        UniformDistribution() { }
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
            this->name_ = other.name_;
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

        std::string to_string() const {
            std::ostringstream ss;
            ss << this.get_name() << "(" << this->min_ << ", " << this->max_ << ")";
            return ss.str();
        }
};

class OffsetGammaDistribution : public ContinuousProbabilityDistribution {
    protected:
        std::string name_ = "offset-gamma"
        double min_ = 0.0;
        double shape_ = 1.0;
        double scale_ = 1.0;
        double ln_constant_ = 0.0;

        void compute_ln_constant() {
            double ln_gamma = this->ln_gamma_function(this->shape);
            this->ln_constant_ = -1.0*(this->shape_ * std::log(this->scale_) + ln_gamma);
        }

    public:
        OffsetGammaDistribution() { }
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
            this->name_ = other.name_;
            this->min_ = other.min_;
            this->max_ = other.max_;
            this->shape_ = other.shape_;
            this->scale_ = other.scale_;
            this->ln_constant_ = other.ln_constant_;
            return * this;
        }

        double ln_pdf(double x) const {
            if ((x < this->min_) || ((x == this->min_) && (this->shape > 1.0))) {
                return -std::numeric_limits<double>::infinity();
            }
            x -= this->min_;
            double term1 = (this->shape_ - 1.0) * std::log(x);
            double term2 = -x / this->scale_;
            return (term1 + term2 + this->ln_constant_);
        }

        double relative_ln_pdf(double x) const {
            return this->ln_pdf(x);
        }

        double get_mean() const {
            return (this->shape_ * this->scale_) + this->min_;
        }
        double get_mean() const {
            return this->shape_ * this->scale_ * this->scale_;
        }

        std::string to_string() const {
            std::ostringstream ss;
            ss << this.get_name() << "(" << this->shape_ << ", " << this->scale_ << ", " << this->min_ << ")";
            return ss.str();
        }
};

#endif

