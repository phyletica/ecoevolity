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
#include <memory>

#include "rng.hpp"
#include "assert.hpp"

class OperatorSchedule {
    protected:
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

        Operator& draw_operator(RandomNumberGenerator& rng) {
            double u = this->rng.uniform_real();
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
        virtual ~Operator() { }
		enum TargetTypeEnum {
            ComparisonPopulationTree = 1,
            ComparisonPopulationTreeCollection = 2
        };

        virtual Operator::TargetTypeEnum get_target_type() const = 0;

        virtual void optimize(double log_alpha) = 0;

        double get_target_acceptance_probability() {
            return 0.234;
        }

        void set_operator_schedule(std::shared_ptr<OperatorSchedule> os) {
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
        std::shared_ptr<OperatorSchedule> operator_schedule_;
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

class ComparisonCollectionOperator : public Operator {
    public:
        ComparisonCollectionOperator() : Operator() { }
        virtual ~ComparisonCollectionOperator() { }

        Operator::TargetTypeEnum get_target_type() const {
            return Operator::TargetTypeEnum::ComparisonPopulationTreeCollection;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons) = 0;
};

class ComparisonTreeOperator : public Operator {
    public:
        ComparisonTreeOperator() : Operator() { }
        virtual ~ComparisonTreeOperator() { }

        Operator::TargetTypeEnum get_target_type() const {
            return Operator::TargetTypeEnum::ComparisonPopulationTree;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) = 0;
};

class MutationRateMover : public ComparisonTreeOperator {
    protected:
        double window_size_ = 0.1;

    public:
        MutationRateMover() : ComparisonTreeOperator() { }
        MutationRateMover(double window_size) : ComparisonTreeOperator() {
            this->set_window_size(window_size);
        }
        virtual ~MutationRateMover() { }

        void set_window_size(double window_size) {
            this->window_size_ = window_size;
        }
        double get_window_size() const {
            return this->window_size_;
        }
        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) {
            double red_freq = tree->get_v() / (tree->get_u() + tree->get_v());
            double x = this->window_size_ * ((2.0 * rng.uniform_real()) - 1.0);
            if ((red_freq + x >= 0.0) && (red_freq + x <= 1.0)) {
                red_freq += x;
                this->set_u(1.0 / (2.0 * red_freq));
                this->set_v(1.0 / (2.0 * (1.0 - red_freq)));
                return 0.0;
            }
            return -std::numeric_limits<double>::infinity();
        }
};

class ChildCoalescenceRateScaler : public ComparisonTreeOperator {
    protected:
        double scale_ = 0.5;

    public:
        CoalescenceRateScaler() : ComparisonTreeOperator() { }
        CoalescenceRateScaler(double scale) : ComparisonTreeOperator() {
            this->set_scale(scale);
        }
        virtual ~CoalescenceRateScaler() { }

        void set_scale(double scale) {
            this->scale_ = scale;
        }
        double get_scale() const {
            return this->scale_;
        }
        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) {
            int pop_idx = rng.uniform_int(0, tree->get_leaf_allele_count() - 1);
            double rate = tree->get_child_coalescence_rate(pop_idx);

            double multiplier = std::exp(this->scale_ * ((2.0 * rng.uniform_real()) - 1.0));

            tree->set_child_coalescence_rate(pop_idx, rate * multiplier);

            return std::log(multiplier);
        }

        void optimize(double log_alpha) {
            double delta = this->calc_delta(log_alpha);
            delta += std::log(this->scale_);
            this->scale_ = std::exp(delta);
        }

        double get_coercable_parameter_value() {
            return this->scale_;
        }

        void set_coercable_parameter_value(double value) {
            this->scale_ = value;
        }
};

class RootCoalescenceRateScaler : public ChildCoalescenceRateScaler {
        RootCoalescenceRateScaler() : ChildCoalescenceRateScaler() { }
        RootCoalescenceRateScaler(double scale) : ChildCoalescenceRateScaler(scale) { }
        virtual ~RootCoalescenceRateScaler() { }

        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) {
            double rate = tree->get_root_coalescence_rate();

            double multiplier = std::exp(this->scale_ * ((2.0 * rng.uniform_real()) - 1.0));

            tree->set_root_coalescence_rate(rate * multiplier);

            return std::log(multiplier);
        }
};


#endif
