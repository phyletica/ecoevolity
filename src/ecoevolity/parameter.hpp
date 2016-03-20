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

#ifndef ECOEVOLITY_PARAMETER_HPP
#define ECOEVOLITY_PARAMETER_HPP

#include "error.hpp"
#include "assert.hpp"
#include "rng.hpp"
#include "probability.hpp"


template<class VariableType>
class Variable {
    protected:
        VariableType value_;
        VariableType stored_value_;
        VariableType upper_;
        VariableType lower_;
        VariableType min_;
        VariableType max_;
        typedef Variable<VariableType> DerivedClass;

    public:
        // Constructors
        Variable() {
            this->upper_ = std::numeric_limits<VariableType>::max();
            this->lower_ = std::numeric_limits<VariableType>::lowest();
            this->max_ = this->upper_;
            if (std::numeric_limits<VariableType>::has_infinity) {
                this->max_ = std::numeric_limits<VariableType>::infinity();
            }
            this->min_ = 0;
            if (std::numeric_limits<VariableType>::is_signed) {
                this->min_ = -this->max_;
            }
        }
        Variable(const VariableType& value) : Variable() {
            this->set_value(value);
        }

        virtual ~Variable() { }

        DerivedClass& operator=(const DerivedClass& p) {
            this->value_ = p.value_;
            this->stored_value_ = p.stored_value_;
            this->upper_ = p.upper_;
            this->lower_ = p.lower_;
            this->min_ = p.min_;
            this->max_ = p.max_;
            return * this;
        }

        DerivedClass* clone() const {
            return new DerivedClass(static_cast<DerivedClass const &>(* this));
        }
        
        //Methods
        const VariableType& get_value() const { return this->value_; }
        const VariableType& get_stored_value() const { return this->stored_value_; }
        const VariableType& get_max() const { return this->max_; }
        const VariableType& get_min() const { return this->min_; }
        const VariableType& get_upper() const { return this->upper_; }
        const VariableType& get_lower() const { return this->lower_; }

        void update_value(const VariableType& value) {
            this->store();
            this->set_value(value);
        }
        void set_value(const VariableType& value) {
            if ((value < this->lower_) || (value > this->upper_)) {
                throw EcoevolityParameterValueError("value outside of parameter bounds");
            }
            this->value_ = value;
        }
        void set_upper(const VariableType& upper) {
            this->upper_ = upper;
        }
        void set_lower(const VariableType& lower) {
            this->lower_ = lower;
        }
        void set_bounds(const VariableType& lower, const VariableType& upper) {
            this->set_lower(lower);
            this->set_upper(upper);
        }

        void store() {
            this->stored_value_ = this->value_;
        }

        void restore() {
            this->value_ = this->stored_value_;
        }

        virtual std::string to_string() {
            return std::to_string(this->value_);
        }
};

class RealVariable: public Variable<double> {
    protected:
        typedef Variable<double> BaseClass;

    public:
        RealVariable() : BaseClass() { }
        RealVariable(double value) : BaseClass() {
            this->set_value(value);
        }
};

class RealParameter: public RealVariable {
    public:
        ContinuousProbabilityDistribution * prior = NULL;

        RealParameter()
                : RealVariable()
                { }
        RealParameter(ContinuousProbabilityDistribution * prior_ptr)
                : RealVariable()
        {
            this->prior = prior_ptr;
        }
        RealParameter(double value)
                : RealVariable(value)
                { }
        RealParameter(ContinuousProbabilityDistribution * prior_ptr, double value)
                : RealVariable(value)
        {
            this->prior = prior_ptr;
        }
        virtual ~RealParameter() {
            delete this->prior;
        }

        virtual void set_prior(ContinuousProbabilityDistribution * prior_ptr) {
            this->prior = prior_ptr;
        }
        virtual void check_prior() const {
            if (! this->prior) {
                throw EcoevolityNullPointerError("tried to use null prior of real parameter");
            }
        }

        virtual double draw_from_prior(RandomNumberGenerator & rng) {
            this->check_prior();
            return this->prior->draw(rng);
        }
        virtual void set_value_from_prior(RandomNumberGenerator & rng) {
            this->set_value(this->draw_from_prior(rng));
        }
        virtual void update_value_from_prior(RandomNumberGenerator & rng) {
            this->update_value(this->draw_from_prior(rng));
        }
        virtual double get_prior_mean() const {
            this->check_prior();
            return this->prior->get_mean();
        }
        virtual double get_prior_variance() const {
            this->check_prior();
            return this->prior->get_variance();
        }
        virtual double get_prior_min() const {
            this->check_prior();
            return this->prior->get_min();
        }
        virtual double get_prior_max() const {
            this->check_prior();
            return this->prior->get_max();
        }
        virtual std::string get_prior_name() const {
            this->check_prior();
            return this->prior->get_name();
        }
        virtual std::string get_prior_string() const {
            this->check_prior();
            return this->prior->to_string();
        }
        virtual double prior_ln_pdf() const {
            this->check_prior();
            return this->prior->ln_pdf(this->get_value());
        }
        virtual double prior_ln_pdf(double x) const {
            this->check_prior();
            return this->prior->ln_pdf(x);
        }
        virtual double relative_prior_ln_pdf() const {
            this->check_prior();
            return this->prior->relative_ln_pdf(this->get_value());
        }
        virtual double relative_prior_ln_pdf(double x) const {
            this->check_prior();
            return this->prior->relative_ln_pdf(x);
        }
};

class PositiveRealVariable: public RealVariable {
    public:
        PositiveRealVariable() : RealVariable() {
            this->set_lower(0.0);
        }
        PositiveRealVariable(double value) : RealVariable() {
            this->set_lower(0.0);
            this->set_value(value);
        }
};

class PositiveRealParameter: public RealParameter {
    public:
        PositiveRealParameter() : RealParameter()
        {
            this->set_lower(0.0);
        }
        PositiveRealParameter(ContinuousProbabilityDistribution * prior_ptr)
                : RealParameter(prior_ptr)
        {
            this->set_lower(0.0);
        }
        PositiveRealParameter(double value)
                : RealParameter()
        {
            this->set_lower(0.0);
            this->set_value(value);
        }
        PositiveRealParameter(ContinuousProbabilityDistribution * prior_ptr, double value)
                : RealParameter(prior_ptr)
        {
            this->set_lower(0.0);
            this->set_value(value);
        }
};

class IntVariable: public Variable<int> {
    protected:
        typedef Variable<int> BaseClass;

    public:
        IntVariable() : BaseClass() { }
        IntVariable(int value) : BaseClass() {
            this->set_value(value);
        }
};

class Probability: public RealVariable {
    public:
        Probability() : RealVariable() {
            this->set_bounds(0.0, 1.0);
        }
        Probability(double value) : RealVariable() {
            this->set_bounds(0.0, 1.0);
            this->set_value(value);
        }
};

class ProbabilityDensity: public PositiveRealVariable {
    public:
        ProbabilityDensity() : PositiveRealVariable() {}
        ProbabilityDensity(double value) : PositiveRealVariable(value) {}
};

class LogProbability: public RealVariable {
    public:
        LogProbability() : RealVariable() {
            this->set_upper(0.0);
        }
        LogProbability(double value) : RealVariable() {
            this->set_upper(0.0);
            this->set_value(value);
        }
};

class LogProbabilityDensity: public RealVariable {
    public:
        LogProbabilityDensity() : RealVariable() {}
        LogProbabilityDensity(double value) : RealVariable(value) {}
};

#endif
