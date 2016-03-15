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

#ifndef ECOEVOLITY_OPERATOR_HPP
#define ECOEVOLITY_OPERATOR_HPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <limits>

#include "rng.hpp"
#include "assert.hpp"

class OperatorSchedule {
    protected:
        RandomNumberGenerator * rng_ = RandomNumberGenerator();
        std::vector<Operator> operators_;
        double total_weight_ = 0.0;
        std::vector<double> cumulative_probs_;
        unsigned int auto_optimize_delay_ = 10000;
        unsigned int auto_optimize_delay_count_ = 0;
        bool auto_optimize_ = true;

    public:

        void add_operator(Operator o) {
            this->operators_.push_back(o);
            o.set_operator_schedule(this);
            this->total_weight_ += o.get_weight();
            this->cumulative_probs_.push_back(0.0);
            ECOEVOLITY_ASSERT(this->operators_.size() == this->cumulative_probs_.size());
            this->cumulative_probs_.at(0) = this->operators_.at(0).get_weight() / this->total_weight_;
            for (unsigned int i = 1; i < this->operators_.size(); ++i) {
                this->cumulative_probs_.at(i) =
                        (this->operators_.at(i).get_weight() /
                        this->total_weight_) + 
                        this->cumulative_probs_.at(i - 1);
            }
        }

        Operator& draw_operator() {
            double u = this->rng_->uniform_real();
            for (unsigned int i = 0; i < this->cumulative_probs_.size(); ++i) {
                if (u <= this->cumulative_probs_.at(i)) {
                    return this->operators_.at(i);
                }
            return this->operators_.back();
        }

        double calc_delta(const Operator& operator, double log_alpha) {
            if ((this->get_auto_optimize_delay_count() < this->get_auto_optimize_delay()) ||
                    (! this->auto_optimize_)) {
                return 0.0;
            }
            double target = operator.get_target_acceptance_probability();
            double count = (operator.get_number_rejected_for_correction() +
                            operator.get_number_accepted_for_correction() +
                            1.0);
            double delta_p = ((1.0 / count) * (std::exp(std::min(log_alpha, 0)) - target));
            double mx = std::numeric_limits<double>::max();
            if ((delta_p > -mx) && (delta_p < mx)) {
                return delta_p;
            }
            return 0.0;
        }

        const double& get_total_weight() const {
            return this->total_weight_;
        }
        const unsigned int& get_auto_optimize_delay_count() const {
            return this->auto_optimize_delay_count_;
        }
        const unsigned int& get_auto_optimize_delay() const {
            return this->auto_optimize_delay_;
        }

        void write_operator_rates(std::ofstream out) {
            const Operator& op = this->operators_.at(0);
            out << op.header_string();
            out << op.to_string();
            for (unsigned int i = 1; i < this->operators_.size(); ++i) {
                out << this->operators_.at(0).to_string();
            }
        }
};

class Operator {
    public:
        Operator() { }
        virtual ~Operator() {
            delete this->operator_schedule_;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double proposal() = 0;

        virtual void optimize(double log_alpha) = 0;

        double get_target_acceptance_probability() {
            return 0.234;
        }

        void set_operator_schedule(OperatorSchedule * os) {
            this->operator_schedule_ = os;
        }

        virtual double get_coercable_parameter_value() {
            return std::numeric_limits<double>::quiet_NaN();
        }

        virtual void set_coercable_parameter_value(double value) { }

        const double& get_weight() const {
            return this->weight_;
        }

        void accept() {
            ++this->number_accepted_;
            if (this->operator_schedule_->get_auto_optimize_delay_count() >= this->operator_schedule_->get_auto_optimize_delay()) {
                ++this->number_accepted_for_correction_;
            }
        }

        void reject() {
            ++this->number_rejected_;
            if (this->operator_schedule_->get_auto_optimize_delay_count() >= this->operator_schedule_->get_auto_optimize_delay()) {
                ++this->number_rejected_for_correction_;
            }
        }

        const unsigned int& get_number_rejected() const {
            return this->number_rejected_;
        }
        const unsigned int& get_number_accepted() const {
            return this->number_accepted_;
        }
        const unsigned int& get_number_rejected_for_correction() const {
            return this->number_rejected_for_correction_;
        }
        const unsigned int& get_number_accepted_for_correction() const {
            return this->number_accepted_for_correction_;
        }

        const std::string& get_name() const {
            return this->name_;
        }

        std::string header_string() const {
            return "name\tnumber_accepted\tnumber_rejected\tweight\tweight_prob\ttuning_parameter\n";
        }
        std::string to_string() const {
            std::ostringstream ss;
            ss << this->get_name() << "\t" 
               << this->get_number_accepted() << "\t"
               << this->get_number_rejected() << "\t"
               << this->get_weight() << "\t";

            if (this->operator_schedule_->get_total_weight() > 0.0) {
                ss << this->get_weight() / this->operator_schedule_->get_total_weight() << "\t";
            }
            else {
                ss << "\t";
            }

            double tuning = this->get_coercable_parameter_value();
            if (std::isnan(tuning)) {
                ss << "\t";
            }
            else {
                ss << tuning << "\t";
            }
            ss << "\n";
            return ss.str();
        }

    protected:
        OperatorSchedule * operator_schedule_;
        double weight_ = 1.0;
        unsigned int number_rejected_ = 0;
        unsigned int number_accepted_ = 0;
        unsigned int number_rejected_for_correction = 0;
        unsigned int number_accepted_for_correction = 0;
        std::string name_ = "BaseOperator";

        double calc_delta(double log_alpha) {
            return operator_schedule_->calc_delta(this, log_alpha);
        }

};

#endif
