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

#include <memory>

#include "error.hpp"
#include "assert.hpp"
#include "rng.hpp"
#include "probability.hpp"
#include "settings.hpp"


template<class VariableType>
class Variable {
    protected:
        VariableType value_;
        VariableType stored_value_;
        VariableType max_;
        VariableType min_;
        bool is_fixed_ = false;
        bool value_is_initialized_ = false;
        typedef Variable<VariableType> DerivedClass;

    public:
        // Constructors
        Variable() {
            this->max_ = std::numeric_limits<VariableType>::max();
            if (std::numeric_limits<VariableType>::has_infinity) {
                this->max_ = std::numeric_limits<VariableType>::infinity();
            }
            this->min_ = 0;
            if (std::numeric_limits<VariableType>::is_signed) {
                this->min_ = -this->max_;
            }
        }
        Variable(const VariableType& value, bool fix = false) : Variable() {
            this->set_value(value);
            this->is_fixed_ = fix;
        }

        virtual ~Variable() { }

        DerivedClass& operator=(const DerivedClass& p) {
            this->value_ = p.value_;
            this->stored_value_ = p.stored_value_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            return * this;
        }
        bool operator< (const DerivedClass & other) const {
            return this->get_value() < other.get_value();
        }
        bool operator> (const DerivedClass & other) const {
            return this->get_value() > other.get_value();
        }
        static bool sort_by_value(
                const std::shared_ptr<DerivedClass>& v1,
                const std::shared_ptr<DerivedClass>& v2) {
            return v1->get_value() < v2->get_value();
        }

        DerivedClass* clone() const {
            return new DerivedClass(static_cast<DerivedClass const &>(* this));
        }
        
        //Methods
        VariableType get_value() const { return this->value_; }
        VariableType get_stored_value() const { return this->stored_value_; }
        VariableType get_max() const { return this->max_; }
        VariableType get_min() const { return this->min_; }

        bool is_fixed() const { return this->is_fixed_; }
        void fix() {
            if (! this->value_is_initialized_) {
                throw EcoevolityParameterValueError(
                        "cannot fix parameter with uninitialized value");
            }
            this->is_fixed_ = true;
        }
        void estimate() {
            this->is_fixed_ = false;
        }

        void update_value(const VariableType& value) {
            if (this->is_fixed()) {
                return;
            }
            this->store();
            this->set_value(value);
        }
        void set_value(const VariableType& value) {
            if (this->is_fixed()) {
                return;
            }
            if ((value < this->min_) || (value > this->max_)) {
                throw EcoevolityParameterValueError("value outside of parameter bounds");
            }
            this->value_ = value;
            this->value_is_initialized_ = true;
        }
        void set_max(const VariableType& max) {
            this->max_ = max;
        }
        void set_min(const VariableType& min) {
            this->min_ = min;
        }
        void set_bounds(const VariableType& min, const VariableType& max) {
            this->set_min(min);
            this->set_max(max);
        }

        void store() {
            this->stored_value_ = this->value_;
        }

        void restore() {
            if (this->is_fixed()) {
                return;
            }
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
        RealVariable() : BaseClass() {
            this->value_ = std::numeric_limits<double>::quiet_NaN();
        }
        RealVariable(double value, bool fix = false) : BaseClass() {
            this->set_value(value);
            this->is_fixed_ = fix;
        }
        virtual ~RealVariable() { }
        RealVariable& operator=(const RealVariable& p) {
            this->value_ = p.value_;
            this->stored_value_ = p.stored_value_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            return * this;
        }
};

class RealParameter: public RealVariable {
    public:
        std::shared_ptr<ContinuousProbabilityDistribution> prior;

        RealParameter()
                : RealVariable()
                { }
        RealParameter(std::shared_ptr<ContinuousProbabilityDistribution> prior_ptr)
                : RealVariable()
        {
            this->prior = prior_ptr;
        }
        RealParameter(double value, bool fix = false)
                : RealVariable(value, fix)
                { }
        RealParameter(std::shared_ptr<ContinuousProbabilityDistribution> prior_ptr, double value, bool fix = false)
                : RealVariable(value, fix)
        {
            this->prior = prior_ptr;
        }
        virtual ~RealParameter() { }
        RealParameter& operator=(const RealParameter& p) {
            this->value_ = p.value_;
            this->stored_value_ = p.stored_value_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            this->prior = p.prior;
            return * this;
        }

        virtual std::shared_ptr<ContinuousProbabilityDistribution> get_prior() const {
            return this->prior;
        }
        virtual void set_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior_ptr) {
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
            if (this->is_fixed()) {
                return;
            }
            this->set_value(this->draw_from_prior(rng));
        }
        virtual void initialize_value(RandomNumberGenerator & rng) {
            if (std::isnan(this->value_)) {
                this->set_value_from_prior(rng);
            }
        }
        virtual void update_value_from_prior(RandomNumberGenerator & rng) {
            if (this->is_fixed()) {
                return;
            }
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
            if (this->is_fixed()) {
                return 0.0;
            }
            this->check_prior();
            return this->prior->ln_pdf(this->get_value());
        }
        virtual double prior_ln_pdf(double x) const {
            if (this->is_fixed()) {
                return 0.0;
            }
            this->check_prior();
            return this->prior->ln_pdf(x);
        }
        virtual double relative_prior_ln_pdf() const {
            if (this->is_fixed()) {
                return 0.0;
            }
            this->check_prior();
            return this->prior->relative_ln_pdf(this->get_value());
        }
        virtual double relative_prior_ln_pdf(double x) const {
            if (this->is_fixed()) {
                return 0.0;
            }
            this->check_prior();
            return this->prior->relative_ln_pdf(x);
        }
};

class PositiveRealVariable: public RealVariable {
    public:
        PositiveRealVariable() : RealVariable() {
            this->set_min(0.0);
        }
        PositiveRealVariable(double value, bool fix = false) : RealVariable() {
            this->set_min(0.0);
            this->set_value(value);
            this->is_fixed_ = fix;
        }
        virtual ~PositiveRealVariable() { }
        PositiveRealVariable& operator=(const PositiveRealVariable& p) {
            this->value_ = p.value_;
            this->stored_value_ = p.stored_value_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            return * this;
        }
};

class PositiveRealParameter: public RealParameter {
    public:
        PositiveRealParameter() : RealParameter()
        {
            this->set_min(0.0);
        }
        PositiveRealParameter(std::shared_ptr<ContinuousProbabilityDistribution> prior_ptr)
                : RealParameter(prior_ptr)
        {
            this->set_min(0.0);
        }
        PositiveRealParameter(double value, bool fix = false)
                : RealParameter()
        {
            this->set_min(0.0);
            this->set_value(value);
            this->is_fixed_ = fix;
        }
        PositiveRealParameter(std::shared_ptr<ContinuousProbabilityDistribution> prior_ptr, double value, bool fix = false)
                : RealParameter(prior_ptr)
        {
            this->set_min(0.0);
            this->set_value(value);
            this->is_fixed_ = fix;
        }
        PositiveRealParameter(
                const PositiveRealParameterSettings& settings,
                RandomNumberGenerator& rng) {
            this->set_min(0.0);
            this->set_value(settings.get_value());
            this->is_fixed_ = settings.is_fixed();
            if (! this->is_fixed()) {
                this->set_prior(settings.get_prior_settings().get_instance());
            }
            this->initialize_value(rng);
        }
        virtual ~PositiveRealParameter() { }
        PositiveRealParameter& operator=(const PositiveRealParameter& p) {
            this->value_ = p.value_;
            this->stored_value_ = p.stored_value_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            this->prior = p.prior;
            return * this;
        }
};

class DiscountParameter: public RealParameter {
    public:
        DiscountParameter() : RealParameter()
        {
            this->set_min(0.0);
            this->set_max(1.0);
        }
        DiscountParameter(std::shared_ptr<ContinuousProbabilityDistribution> prior_ptr)
                : RealParameter(prior_ptr)
        {
            this->set_min(0.0);
            this->set_max(1.0);
        }
        DiscountParameter(double value, bool fix = false)
                : RealParameter()
        {
            this->set_min(0.0);
            this->set_max(1.0);
            this->set_value(value);
            this->is_fixed_ = fix;
        }
        DiscountParameter(std::shared_ptr<ContinuousProbabilityDistribution> prior_ptr, double value, bool fix = false)
                : RealParameter(prior_ptr)
        {
            this->set_min(0.0);
            this->set_max(1.0);
            this->set_value(value);
            this->is_fixed_ = fix;
        }
        DiscountParameter(
                const PositiveRealParameterSettings& settings,
                RandomNumberGenerator& rng) {
            this->set_min(0.0);
            this->set_max(1.0);
            this->set_value(settings.get_value());
            this->is_fixed_ = settings.is_fixed();
            if (! this->is_fixed()) {
                this->set_prior(settings.get_prior_settings().get_instance());
            }
            this->initialize_value(rng);
        }
        virtual ~DiscountParameter() { }
        DiscountParameter& operator=(const DiscountParameter& p) {
            this->value_ = p.value_;
            this->stored_value_ = p.stored_value_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            this->prior = p.prior;
            return * this;
        }

        // override set_value method, because the discount parameter cannot be
        // equal to 1.0
        void set_value(const double& value) {
            if (this->is_fixed()) {
                return;
            }
            if ((value < this->min_) || (value >= this->max_)) {
                throw EcoevolityParameterValueError("value outside of parameter bounds");
            }
            this->value_ = value;
            this->value_is_initialized_ = true;
        }
};

class CoalescenceRateParameter: public PositiveRealParameter {
    public:
        CoalescenceRateParameter() : PositiveRealParameter() { }
        CoalescenceRateParameter(std::shared_ptr<ContinuousProbabilityDistribution> prior_ptr)
                : PositiveRealParameter(prior_ptr) { }
        CoalescenceRateParameter(double value, bool fix = false)
                : PositiveRealParameter(value, fix) { }
        CoalescenceRateParameter(std::shared_ptr<ContinuousProbabilityDistribution> prior_ptr,
                double value,
                bool fix = false)
                : PositiveRealParameter(prior_ptr, value, fix) { }
        virtual ~CoalescenceRateParameter() { }
        CoalescenceRateParameter& operator=(const CoalescenceRateParameter& p) {
            this->value_ = p.value_;
            this->stored_value_ = p.stored_value_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            this->prior = p.prior;
            return * this;
        }

        static double get_population_size_from_rate(double coalescence_rate) {
            if (coalescence_rate == 0.0) {
                return std::numeric_limits<double>::infinity();
            }
            return 2.0 / coalescence_rate;
        }
        static double get_rate_from_population_size(double population_size) {
            if (population_size == 0.0) {
                return std::numeric_limits<double>::infinity();
            }
            return 2.0 / population_size;
        }
        double get_population_size() const {
            return this->get_population_size_from_rate(this->get_value());
        }
        void set_population_size(double population_size) {
            this->set_value(this->get_rate_from_population_size(population_size));
        }

        virtual double draw_from_prior(RandomNumberGenerator & rng) {
            this->check_prior();
            return this->get_rate_from_population_size(this->prior->draw(rng));
        }
        virtual void set_value_from_prior(RandomNumberGenerator & rng) {
            if (this->is_fixed()) {
                return;
            }
            this->set_value(this->draw_from_prior(rng));
        }
        virtual void initialize_value(RandomNumberGenerator & rng) {
            if (std::isnan(this->value_)) {
                this->set_value_from_prior(rng);
            }
        }
        virtual void update_value_from_prior(RandomNumberGenerator & rng) {
            if (this->is_fixed()) {
                return;
            }
            this->update_value(this->draw_from_prior(rng));
        }

        double prior_ln_pdf() const {
            if (this->is_fixed()) {
                return 0.0;
            }
            this->check_prior();
            return this->prior->ln_pdf(this->get_population_size());
        }
        double prior_ln_pdf(double coalescence_rate) const {
            if (this->is_fixed()) {
                return 0.0;
            }
            this->check_prior();
            return this->prior->ln_pdf(
                    this->get_population_size_from_rate(coalescence_rate));
        }
        double relative_prior_ln_pdf() const {
            if (this->is_fixed()) {
                return 0.0;
            }
            this->check_prior();
            return this->prior->relative_ln_pdf(this->get_population_size());
        }
        double relative_prior_ln_pdf(double coalescence_rate) const {
            if (this->is_fixed()) {
                return 0.0;
            }
            this->check_prior();
            return this->prior->relative_ln_pdf(
                    this->get_population_size_from_rate(coalescence_rate));
        }
};

class IntVariable: public Variable<int> {
    protected:
        typedef Variable<int> BaseClass;

    public:
        IntVariable() : BaseClass() { }
        IntVariable(int value, bool fix = false) : BaseClass() {
            this->set_value(value);
            this->is_fixed_ = fix;
        }
        virtual ~IntVariable() { }
        IntVariable& operator=(const IntVariable& p) {
            this->value_ = p.value_;
            this->stored_value_ = p.stored_value_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            return * this;
        }
};

class Probability: public RealVariable {
    public:
        Probability() : RealVariable() {
            this->set_bounds(0.0, 1.0);
        }
        Probability(double value, bool fix = false) : RealVariable() {
            this->set_bounds(0.0, 1.0);
            this->set_value(value);
            this->is_fixed_ = fix;
        }
        virtual ~Probability() { }
        Probability& operator=(const Probability& p) {
            this->value_ = p.value_;
            this->stored_value_ = p.stored_value_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            return * this;
        }
};

class ProbabilityDensity: public PositiveRealVariable {
    public:
        ProbabilityDensity() : PositiveRealVariable() {}
        ProbabilityDensity(double value, bool fix = false) : PositiveRealVariable(value, fix) {}
        virtual ~ProbabilityDensity() { }
        ProbabilityDensity& operator=(const ProbabilityDensity& p) {
            this->value_ = p.value_;
            this->stored_value_ = p.stored_value_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            return * this;
        }
};

class LogProbability: public RealVariable {
    public:
        LogProbability() : RealVariable() {
            this->set_max(0.0);
        }
        LogProbability(double value, bool fix = false) : RealVariable() {
            this->set_max(0.0);
            this->set_value(value);
            this->is_fixed_ = fix;
        }
        virtual ~LogProbability() { }
        LogProbability& operator=(const LogProbability& p) {
            this->value_ = p.value_;
            this->stored_value_ = p.stored_value_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            return * this;
        }
};

class LogProbabilityDensity: public RealVariable {
    public:
        LogProbabilityDensity() : RealVariable() {}
        LogProbabilityDensity(double value, bool fix = false) : RealVariable(value, fix) {}
        virtual ~LogProbabilityDensity() { }
        LogProbabilityDensity& operator=(const LogProbabilityDensity& p) {
            this->value_ = p.value_;
            this->stored_value_ = p.stored_value_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            return * this;
        }
};

#endif
