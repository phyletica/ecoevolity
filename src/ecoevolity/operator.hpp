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

        const unsigned int& get_auto_optimize_delay_count() const {
            return this->auto_optimize_delay_count_;
        }
        const unsigned int& get_auto_optimize_delay() const {
            return this->auto_optimize_delay_;
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

        virtual double get_coercable_parameter_value() = 0;

        virtual void set_coercable_parameter_value(double value) = 0;

        const double& get_weight() const {
            return this->weight_;
        }
        void set_weight(double weight) {
            this->weight_ = weight;
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

    protected:
        OperatorSchedule * operator_schedule_;
        double weight_ = 1.0;
        unsigned int number_rejected_ = 0;
        unsigned int number_accepted_ = 0;
        unsigned int number_rejected_for_correction = 0;
        unsigned int number_accepted_for_correction = 0;

        double calc_delta(double log_alpha) {
            return operator_schedule_->calc_delta(this, log_alpha);
        }

};

#endif
