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

template<class ParameterType>
class Parameter {
    protected:
        ParameterType value_;
        ParameterType stored_value_;
        ParameterType upper_;
        ParameterType lower_;
        ParameterType min_;
        ParameterType max_;
        typedef Parameter<ParameterType> DerivedClass;

    public:
        // Constructors
        Parameter() {
            this->upper_ = std::numeric_limits<ParameterType>::max();
            this->lower_ = std::numeric_limits<ParameterType>::lowest();
            this->max_ = this->upper_;
            if (std::numeric_limits<ParameterType>::has_infinity) {
                this->max_ = std::numeric_limits<ParameterType>::infinity();
            }
            this->min_ = 0;
            if (std::numeric_limits<ParameterType>::is_signed) {
                this->min_ = -this->max_;
            }
        }
        Parameter(const ParameterType& value) : Parameter() {
            this->set_value(value);
        }

        virtual ~Parameter() { }

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
        const ParameterType& get_value() const { return this->value_; }
        const ParameterType& get_stored_value() const { return this->stored_value_; }
        const ParameterType& get_max() const { return this->max_; }
        const ParameterType& get_min() const { return this->min_; }
        const ParameterType& get_upper() const { return this->upper_; }
        const ParameterType& get_lower() const { return this->lower_; }

        void update_value(const ParameterType& value) {
            this->store();
            this->set_value(value);
        }
        void set_value(const ParameterType& value) {
            if ((value < this->lower_) || (value > this->upper_)) {
                throw EcoevolityParameterValueError("value outside of parameter bounds");
            }
            this->value_ = value;
        }
        void set_upper(const ParameterType& upper) {
            this->upper_ = upper;
        }
        void set_lower(const ParameterType& lower) {
            this->lower_ = lower;
        }
        void set_bounds(const ParameterType& lower, const ParameterType& upper) {
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

class RealParameter: public Parameter<double> {
    protected:
        typedef Parameter<double> BaseClass;

    public:
        RealParameter() : BaseClass() { }
        RealParameter(double value) : BaseClass() {
            this->set_value(value);
        }
};

class PositiveRealParameter: public RealParameter {
    public:
        PositiveRealParameter() : RealParameter() {
            this->set_lower(0.0);
        }
        PositiveRealParameter(double value) : RealParameter() {
            this->set_lower(0.0);
            this->set_value(value);
        }
};

class IntParameter: public Parameter<int> {
    protected:
        typedef Parameter<int> BaseClass;

    public:
        IntParameter() : BaseClass() { }
        IntParameter(int value) : BaseClass() {
            this->set_value(value);
        }
};

class Probability: public RealParameter {
    public:
        Probability() : RealParameter() {
            this->set_bounds(0.0, 1.0);
        }
        Probability(double value) : RealParameter() {
            this->set_bounds(0.0, 1.0);
            this->set_value(value);
        }
};

class LogProbability: public RealParameter {
    public:
        LogProbability() : RealParameter() {
            this->set_upper(0.0);
        }
        LogProbability(double value) : RealParameter() {
            this->set_upper(0.0);
            this->set_value(value);
        }
};

#endif
