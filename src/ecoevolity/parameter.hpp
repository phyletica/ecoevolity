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
#include "math_util.hpp"


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
            else {
                if (this->prior && (! this->prior->is_within_support(this->get_value()))) {
                    std::ostringstream ss;
                    ss << "Initial value "
                       << this->get_value()
                       << " is outside the support of prior "
                       << this->prior->to_string();
                    throw EcoevolityParameterValueError(ss.str());
                }
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

class ProportionParameter: public RealParameter {
    public:
        ProportionParameter() : RealParameter()
        {
            this->set_min(0.0);
            this->set_max(1.0);
        }
        ProportionParameter(std::shared_ptr<ContinuousProbabilityDistribution> prior_ptr)
                : RealParameter(prior_ptr)
        {
            this->set_min(0.0);
            this->set_max(1.0);
        }
        ProportionParameter(double value, bool fix = false)
                : RealParameter()
        {
            this->set_min(0.0);
            this->set_max(1.0);
            this->set_value(value);
            this->is_fixed_ = fix;
        }
        ProportionParameter(std::shared_ptr<ContinuousProbabilityDistribution> prior_ptr, double value, bool fix = false)
                : RealParameter(prior_ptr)
        {
            this->set_min(0.0);
            this->set_max(1.0);
            this->set_value(value);
            this->is_fixed_ = fix;
        }
        ProportionParameter(
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
        virtual ~ProportionParameter() { }
        ProportionParameter& operator=(const ProportionParameter& p) {
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
            else {
                if (this->prior && (! this->prior->is_within_support(this->get_value()))) {
                    std::ostringstream ss;
                    ss << "Initial value "
                       << this->get_value()
                       << " is outside the support of prior "
                       << this->prior->to_string();
                    throw EcoevolityParameterValueError(ss.str());
                }
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


class RealVariableVector {
    public:
        typedef double                      variable_t;
        typedef std::vector<variable_t>     vector_t;

    protected:
        vector_t                            values_;
        vector_t                            stored_values_;
        variable_t                          max_;
        variable_t                          min_;
        bool                                is_sum_constrained_ = false;
        variable_t                          sum_constraint_ = 0.0;
        variable_t                          sum_epsilon_ = 1e-6;
        bool                                is_fixed_ = false;
        bool                                values_are_initialized_ = false;

    public:
        // Constructors
        RealVariableVector() {
            this->max_ = std::numeric_limits<variable_t>::infinity();
            this->min_ = -this->max_;
        }
        RealVariableVector(const vector_t & values, bool fix = false) : RealVariableVector() {
            this->set_values(values);
            this->is_fixed_ = fix;
        }

        virtual ~RealVariableVector() { }

        RealVariableVector& operator=(const RealVariableVector& p) {
            this->values_ = p.values_;
            this->stored_values_ = p.stored_values_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            this->is_sum_constrained_ = p.is_sum_constrained_;
            this->sum_constraint_ = p.sum_constraint_;
            this->sum_epsilon_ = p.sum_epsilon_;
            return * this;
        }

        //Methods
        const vector_t & get_values() const { return this->values_; }
        const variable_t * get_raw_values() const { return this->values_.data(); }
        const vector_t & get_stored_values() const { return this->stored_values_; }
        variable_t get_max() const { return this->max_; }
        variable_t get_min() const { return this->min_; }
        variable_t get_sum_constraint() const { return this->sum_constraint_; }
        variable_t get_sum_epsilon() const { return this->sum_epsilon_; }
        bool is_sum_constrained() const {return is_sum_constrained_; }

        bool is_fixed() const { return this->is_fixed_; }
        void fix() {
            if (! this->values_are_initialized_) {
                throw EcoevolityParameterValueError(
                        "cannot fix parameter vector with uninitialized values");
            }
            this->is_fixed_ = true;
        }
        void estimate() {
            this->is_fixed_ = false;
        }

        void set_sum_constraint(variable_t c) {
            if ((c < this->min_) || (c > this->max_)) {
                throw EcoevolityParameterValueError(
                        "sum constraint is outside min/max");
            }
            this->sum_constraint_ = c;
            this->is_sum_constrained_ = true;
        }

        void set_sum_epsilon(variable_t epsilon) {
            if (epsilon < 0.0) {
                throw EcoevolityParameterValueError(
                        "sum epsilon cannot be negative");
            }
            this->sum_epsilon_ = epsilon;
        }

        void update_values(const vector_t & values) {
            if (this->is_fixed()) {
                return;
            }
            if (! this->values_are_initialized_) {
                throw EcoevolityParameterValueError(
                        "calling update_values but values have not been set yet");
            }
            this->store();
            this->set_values(values);
        }
        void set_values(const vector_t & values) {
            if (this->is_fixed()) {
                return;
            }
            for (const auto v : values) {
                if ((v < this->min_) || (v > this->max_)) {
                    throw EcoevolityParameterValueError("value outside of parameter bounds");
                }
            }
            if (! this->values_are_initialized_) {
                this->values_.resize(values.size(), 0.0);
            }
            else if (values.size() != this->values_.size()) {
                throw EcoevolityParameterValueError("value vector has incorrect number of values");
            }
            if (this->is_sum_constrained_) {
                variable_t sum = std::accumulate(values.begin(), values.end(), 0.0);
                double abs_diff = std::abs(sum - this->sum_constraint_);
                if (abs_diff > this->sum_epsilon_) {
                    throw EcoevolityParameterValueError("values do not sum to constraint");
                }
            }
            this->values_ = values;
            this->values_are_initialized_ = true;
        }
        void set_max(const variable_t & max) {
            this->max_ = max;
        }
        void set_min(const variable_t & min) {
            this->min_ = min;
        }
        void set_bounds(const variable_t & min, const variable_t & max) {
            this->set_min(min);
            this->set_max(max);
        }

        void store() {
            this->stored_values_ = this->values_;
        }

        void restore() {
            if (this->is_fixed()) {
                return;
            }
            this->values_ = this->stored_values_;
        }
};

template<class PriorType>
class RealParameterVector: public RealVariableVector {
    public:
        typedef RealVariableVector::variable_t  variable_t;
        typedef RealVariableVector::vector_t    vector_t;

        typedef RealParameterVector<PriorType>  DerivedClass;
        typedef RealVariableVector              BaseClass;

        std::shared_ptr<PriorType> prior;

        RealParameterVector()
                : BaseClass()
                { }
        RealParameterVector(std::shared_ptr<PriorType> prior_ptr)
                : BaseClass()
        {
            this->prior = prior_ptr;
        }
        RealParameterVector(BaseClass::vector_t values, bool fix = false)
                : BaseClass(values, fix)
                { }
        RealParameterVector(std::shared_ptr<PriorType> prior_ptr, BaseClass::vector_t values, bool fix = false)
                : BaseClass(values, fix)
        {
            this->prior = prior_ptr;
        }
        virtual ~RealParameterVector() { }
        DerivedClass& operator=(const DerivedClass& p) {
            this->values_ = p.values_;
            this->stored_values_ = p.stored_values_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            this->is_sum_constrained_ = p.is_sum_constrained_;
            this->sum_constraint_ = p.sum_constraint_;
            this->sum_epsilon_ = p.sum_epsilon_;
            this->prior = p.prior;
            return * this;
        }

        virtual std::shared_ptr<PriorType> get_prior() const {
            return this->prior;
        }
        virtual void set_prior(std::shared_ptr<PriorType> prior_ptr) {
            this->prior = prior_ptr;
        }
        virtual void check_prior() const {
            if (! this->prior) {
                throw EcoevolityNullPointerError("tried to use null prior of real parameter");
            }
        }

        virtual BaseClass::vector_t draw_from_prior(RandomNumberGenerator & rng) {
            this->check_prior();
            return this->prior->draw(rng);
        }
        virtual void set_values_from_prior(RandomNumberGenerator & rng) {
            if (this->is_fixed()) {
                return;
            }
            this->set_values(this->draw_from_prior(rng));
        }
        virtual void initialize_values(RandomNumberGenerator & rng) {
            if (! this->values_are_initialized_) {
                this->set_values_from_prior(rng);
            }
            else {
                if (this->prior && (! this->prior->is_within_support(this->get_values()))) {
                    throw EcoevolityParameterValueError("parameter vector value is outside support of prior");
                }
            }
        }
        virtual void update_values_from_prior(RandomNumberGenerator & rng) {
            if (this->is_fixed()) {
                return;
            }
            this->update_values(this->draw_from_prior(rng));
        }
        virtual std::vector<double> get_prior_mean() const {
            this->check_prior();
            return this->prior->get_mean();
        }
        virtual std::vector<double> get_prior_variance() const {
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
            return this->prior->ln_pdf(this->get_values());
        }
        virtual double prior_ln_pdf(std::vector<double> x) const {
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
            return this->prior->relative_ln_pdf(this->get_values());
        }
        virtual double relative_prior_ln_pdf(std::vector<double> x) const {
            if (this->is_fixed()) {
                return 0.0;
            }
            this->check_prior();
            return this->prior->relative_ln_pdf(x);
        }
};

class PositiveRealVariableVector: public RealVariableVector {
    public:
        typedef RealVariableVector::variable_t  variable_t;
        typedef RealVariableVector::vector_t    vector_t;

        PositiveRealVariableVector() : RealVariableVector() {
            this->set_min(0.0);
        }
        PositiveRealVariableVector(const vector_t & values, bool fix = false) : RealVariableVector() {
            this->set_min(0.0);
            this->set_values(values);
            this->is_fixed_ = fix;
        }
        virtual ~PositiveRealVariableVector() { }
        PositiveRealVariableVector& operator=(const PositiveRealVariableVector& p) {
            this->values_ = p.values_;
            this->stored_values_ = p.stored_values_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            this->is_sum_constrained_ = p.is_sum_constrained_;
            this->sum_constraint_ = p.sum_constraint_;
            this->sum_epsilon_ = p.sum_epsilon_;
            return * this;
        }
};

template<class PriorType>
class PositiveRealParameterVector: public RealParameterVector<PriorType> {
    public:
        typedef RealVariableVector::variable_t  variable_t;
        typedef RealVariableVector::vector_t    vector_t;

        typedef RealParameterVector<PriorType> BaseClass;
        typedef PositiveRealParameterVector<PriorType>  DerivedClass;

        PositiveRealParameterVector() : BaseClass()
        {
            this->set_min(0.0);
        }
        PositiveRealParameterVector(std::shared_ptr<PriorType> prior_ptr)
                : BaseClass(prior_ptr)
        {
            this->set_min(0.0);
        }
        PositiveRealParameterVector(const vector_t & values, bool fix = false)
                : BaseClass()
        {
            this->set_min(0.0);
            this->set_values(values);
            this->is_fixed_ = fix;
        }
        PositiveRealParameterVector(std::shared_ptr<PriorType> prior_ptr, const vector_t & values, bool fix = false)
                : BaseClass(prior_ptr)
        {
            this->set_min(0.0);
            this->set_values(values);
            this->is_fixed_ = fix;
        }
        PositiveRealParameterVector(
                const PositiveRealParameterSettings& settings,
                RandomNumberGenerator& rng) {
            this->set_min(0.0);
            this->set_values(settings.get_values());
            this->is_fixed_ = settings.is_fixed();
            if (! this->is_fixed()) {
                this->set_prior(settings.get_prior_settings().get_instance());
            }
            this->initialize_values(rng);
        }
        virtual ~PositiveRealParameterVector() { }
        DerivedClass& operator=(const DerivedClass& p) {
            this->values_ = p.values_;
            this->stored_values_ = p.stored_values_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            this->is_sum_constrained_ = p.is_sum_constrained_;
            this->sum_constraint_ = p.sum_constraint_;
            this->sum_epsilon_ = p.sum_epsilon_;
            this->prior = p.prior;
            return * this;
        }
};

class ProportionVariableVector: public PositiveRealVariableVector {
    public:
        typedef RealVariableVector::variable_t  variable_t;
        typedef RealVariableVector::vector_t    vector_t;

        ProportionVariableVector() : PositiveRealVariableVector() {
            this->set_min(0.0);
            this->set_max(1.0);
            this->set_sum_constraint(1.0);
        }
        ProportionVariableVector(const vector_t & values, bool fix = false) : PositiveRealVariableVector() {
            this->set_min(0.0);
            this->set_max(1.0);
            this->set_sum_constraint(1.0);
            this->set_values(values);
            this->is_fixed_ = fix;
        }
        virtual ~ProportionVariableVector() { }
        ProportionVariableVector& operator=(const ProportionVariableVector& p) {
            this->values_ = p.values_;
            this->stored_values_ = p.stored_values_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            this->is_sum_constrained_ = p.is_sum_constrained_;
            this->sum_constraint_ = p.sum_constraint_;
            this->sum_epsilon_ = p.sum_epsilon_;
            return * this;
        }

        void set_unnormalized_values(const vector_t & values) {
            for (const auto & v : values) {
                if (v < 0.0) {
                    throw EcoevolityParameterValueError("unnormalized values cannot be negative");
                }
            }
            vector_t new_values = get_sum_normalized_vector(values);
            this->set_values(new_values);
        }
};

template<class PriorType>
class BaseProportionParameterVector: public PositiveRealParameterVector<PriorType> {
    public:
        typedef RealVariableVector::variable_t  variable_t;
        typedef RealVariableVector::vector_t    vector_t;

        typedef PositiveRealParameterVector<PriorType> BaseClass;
        typedef BaseProportionParameterVector<PriorType> DerivedClass;

        BaseProportionParameterVector() : BaseClass()
        {
            this->set_min(0.0);
            this->set_max(1.0);
            this->set_sum_constraint(1.0);
        }
        BaseProportionParameterVector(std::shared_ptr<PriorType> prior_ptr)
                : BaseClass(prior_ptr)
        {
            this->set_min(0.0);
            this->set_max(1.0);
            this->set_sum_constraint(1.0);
        }
        BaseProportionParameterVector(const vector_t & values, bool fix = false)
                : BaseClass()
        {
            this->set_min(0.0);
            this->set_max(1.0);
            this->set_sum_constraint(1.0);
            this->set_values(values);
            this->is_fixed_ = fix;
        }
        BaseProportionParameterVector(std::shared_ptr<PriorType> prior_ptr, const vector_t & values, bool fix = false)
                : BaseClass(prior_ptr)
        {
            this->set_min(0.0);
            this->set_max(1.0);
            this->set_sum_constraint(1.0);
            this->set_values(values);
            this->is_fixed_ = fix;
        }
        virtual ~BaseProportionParameterVector() { }
        DerivedClass& operator=(const DerivedClass& p) {
            this->values_ = p.values_;
            this->stored_values_ = p.stored_values_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            this->is_sum_constrained_ = p.is_sum_constrained_;
            this->sum_constraint_ = p.sum_constraint_;
            this->sum_epsilon_ = p.sum_epsilon_;
            this->prior = p.prior;
            return * this;
        }

        void set_unnormalized_values(const vector_t & values) {
            for (const auto & v : values) {
                if (v < 0.0) {
                    throw EcoevolityParameterValueError("unnormalized values cannot be negative");
                }
            }
            vector_t new_values = get_sum_normalized_vector(values);
            this->set_values(new_values);
        }
};

class ProportionParameterVector: public BaseProportionParameterVector<DirichletDistribution> {
    public:

        typedef BaseProportionParameterVector<DirichletDistribution> BaseClass;

        ProportionParameterVector() : BaseClass() { }

        ProportionParameterVector(std::shared_ptr<DirichletDistribution> prior_ptr)
                : BaseClass(prior_ptr) { }
        ProportionParameterVector(const vector_t & values, bool fix = false)
                : BaseClass(values, fix) { }
        ProportionParameterVector(std::shared_ptr<DirichletDistribution> prior_ptr, const vector_t & values, bool fix = false)
                : BaseClass(prior_ptr, values, fix) { }
        ProportionParameterVector(
                const PositiveRealParameterSettings& settings,
                RandomNumberGenerator& rng) {
            this->set_min(0.0);
            this->set_max(1.0);
            this->set_sum_constraint(1.0);
            this->set_values(settings.get_values());
            this->is_fixed_ = settings.is_fixed();
            if (! this->is_fixed()) {
                this->set_prior(settings.get_prior_settings().get_dirichlet_distribution_instance());
            }
            this->initialize_values(rng);
        }
        virtual ~ProportionParameterVector() { }
        ProportionParameterVector& operator=(const ProportionParameterVector& p) {
            this->values_ = p.values_;
            this->stored_values_ = p.stored_values_;
            this->max_ = p.max_;
            this->min_ = p.min_;
            this->is_fixed_ = p.is_fixed_;
            this->is_sum_constrained_ = p.is_sum_constrained_;
            this->sum_constraint_ = p.sum_constraint_;
            this->sum_epsilon_ = p.sum_epsilon_;
            this->prior = p.prior;
            return * this;
        }
};

#endif
