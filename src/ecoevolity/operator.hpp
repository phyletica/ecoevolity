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

        double calc_delta(const Operator& op, double log_alpha) {
            if ((this->get_auto_optimize_delay_count() < this->get_auto_optimize_delay()) ||
                    (! this->auto_optimize_)) {
                return 0.0;
            }
            double target = op.get_target_acceptance_probability();
            double count = (op.get_number_rejected_for_correction() +
                            op.get_number_accepted_for_correction() +
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


//////////////////////////////////////////////////////////////////////////////
// Operator base class
//////////////////////////////////////////////////////////////////////////////

class Operator {
    public:
        Operator() { }
        virtual ~Operator() { }
		enum OperatorTypeEnum {
            tree_operator = 1;
            time_operator = 2;
            model_operator = 3;
        };

        virtual Operator::OperatorTypeEnum get_type() const = 0;

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

        virtual std::string get_name() const {
            return "BaseOperator";
        }

        virtual std::string target_parameter() const = 0;

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

        double calc_delta(double log_alpha) {
            return operator_schedule_->calc_delta(this, log_alpha);
        }

};


//////////////////////////////////////////////////////////////////////////////
// Operator Type base class
//////////////////////////////////////////////////////////////////////////////

class ModelOperator : public Operator {
    public:
        ModelOperator() : Operator() { }
        virtual ~ModelOperator() { }

        Operator::OperatorTypeEnum get_type() const {
            return Operator::OperatorTypeEnum::model_operator;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons) = 0;

        std::string get_name() const {
            return "ModelOperator";
        }

        std::string target_parameter() const {;
            return "model";
        }
};

class ComparisonTreeOperator : public Operator {
    public:
        ComparisonTreeOperator() : Operator() { }
        virtual ~ComparisonTreeOperator() { }

        Operator::OperatorTypeEnum get_type() const {
            return Operator::OperatorTypeEnum::tree_operator;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) = 0;

        std::string get_name() const {
            return "ComparisonTreeOperator";
        }
};

class NodeHeightOperator : public Operator {
    public:
        NodeHeightOperator() : Operator() { }
        virtual ~NodeHeightOperator() { }

        Operator::OperatorTypeEnum get_type() const {
            return Operator::OperatorTypeEnum::time_operator;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& node_height) = 0;

        std::string get_name() const {
            return "NodeHeightOperator";
        }

        std::string target_parameter() const {;
            return "node height";
        }
};


//////////////////////////////////////////////////////////////////////////////
// Operator base classes for each parameter
//////////////////////////////////////////////////////////////////////////////

class ConcentrationOperator : public ModelOperator {
    public:
        ConcentrationOperator() : ModelOperator() { }
        virtual ~ConcentrationOperator() { }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons) = 0;

        std::string get_name() const {
            return "ConcentrationOperator";
        }

        std::string target_parameter() const {;
            return "concentration";
        }
};

class MutationRateOperator : public ComparisonTreeOperator {
    public:
        MutationRateOperator() : ComparisonTreeOperator() { }
        virtual ~MutationRateOperator() { }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& node_height) = 0;

        std::string get_name() const {
            return "MutationRateOperator";
        }

        std::string target_parameter() const {;
            return "mutation rate";
        }
};

class CoalescenceRateOperator : public ComparisonTreeOperator {
    public:
        CoalescenceRateOperator() : ComparisonTreeOperator() { }
        virtual ~CoalescenceRateOperator() { }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& node_height) = 0;

        std::string get_name() const {
            return "CoalescenceRateOperator";
        }

        std::string target_parameter() const {;
            return "coalescence rate";
        }
};

class NodeHeightMultiplierOperator : public ComparisonTreeOperator {
    public:
        NodeHeightMultiplierOperator() : ComparisonTreeOperator() { }
        virtual ~NodeHeightMultiplierOperator() { }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& node_height) = 0;

        std::string get_name() const {
            return "NodeHeightMultiplierOperator";
        }

        std::string target_parameter() const {;
            return "node height multiplier";
        }
};


//////////////////////////////////////////////////////////////////////////////
// Derived Operator classes
//////////////////////////////////////////////////////////////////////////////

class MutationRateMover : public MutationRateOperator {
    protected:
        double window_size_ = 0.1;

    public:
        MutationRateMover() : MutationRateOperator() { }
        MutationRateMover(double window_size) : MutationRateOperator() {
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

        std::string get_name() const {
            return "MutationRateMover";
        }
};

class ConcentrationScaler: public ConcentrationOperator {
    protected:
        double scale_ = 0.5;

    public:
        ConcentrationScaler() : ConcentrationOperator() { }
        ConcentrationScaler(double scale) : public ConcentrationOperator() {
            this->set_scale(scale);
        }
        virtual ~ConcentrationScaler() { }

        void set_scale(double scale) {
            this->scale_ = scale;
        }
        double get_scale() const {
            return this->scale_;
        }
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons) {
            double v = comparisons->get_concentration();
            double multiplier = std::exp(this->scale_ * ((2.0 * rng.uniform_real()) - 1.0));
            comparisons->set_concentration(v * multiplier);
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

        std::string get_name() const {
            return "ConcentrationScaler";
        }
};

class NodeHeightMultiplierScaler : public NodeHeightMultiplierOperator {
    protected:
        double scale_ = 0.3;

    public:
        NodeHeightMultiplierScaler() : NodeHeightMultiplierOperator() { }
        NodeHeightMultiplierScaler(double scale) : NodeHeightMultiplierOperator() {
            this->set_scale(scale);
        }
        virtual ~NodeHeightMultiplierScaler() { }

        void set_scale(double scale) {
            this->scale_ = scale;
        }
        double get_scale() const {
            return this->scale_;
        }
        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) {
            double v = tree->get_node_height_multiplier();
            double multiplier = std::exp(this->scale_ * ((2.0 * rng.uniform_real()) - 1.0));
            tree->set_node_height_multiplier(v * multiplier);
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

        std::string get_name() const {
            return "NodeHeightMultiplierScaler";
        }
};

class ChildCoalescenceRateScaler : public CoalescenceRateOperator {
    protected:
        double scale_ = 0.5;

    public:
        ChildCoalescenceRateScaler() : CoalescenceRateOperator() { }
        ChildCoalescenceRateScaler(double scale) : CoalescenceRateOperator() {
            this->set_scale(scale);
        }
        virtual ~ChildCoalescenceRateScaler() { }

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

        std::string get_name() const {
            return "ChildCoalescenceRateScaler";
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

        std::string get_name() const {
            return "RootCoalescenceRateScaler";
        }
};

class ComparisonHeightScaler : public NodeHeightOperator {
    protected:
        double scale_ = 0.5;

    public:
        ComparisonHeightScaler() : NodeHeightOperator() { }
        ComparisonHeightScaler(double scale) : NodeHeightOperator() {
            this->set_scale(scale);
        }
        virtual ~ComparisonHeightScaler() { }

        void set_scale(double scale) {
            this->scale_ = scale;
        }
        double get_scale() const {
            return this->scale_;
        }
        double propose(
                RandomNumberGenerator& rng,
                PositiveRealParameter& node_height) {
            double h = node_height->get_value();
            double multiplier = std::exp(this->scale_ * ((2.0 * rng.uniform_real()) - 1.0));
            node_height.set_value(h * multiplier);
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

        std::string get_name() const {
            return "ComparisonHeightScaler";
        }
};

class DirichletProcessGibbsSampler : public ModelOperator {
    public:
        DirichletProcessGibbsSampler() : ModelOperator() { }
        virtual ~DirichletProcessGibbsSampler() { }

        std::string get_name() const {
            return "DirichletProcessGibbsSampler";
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons) {

            const double ln_concentration_over_num_aux = std::log(
                    comparisons.get_concentration() /
                    comparisons.get_number_of_auxiliary_heights())

            for (unsigned int tree_idx = 0;
                    tree_idx < comparisons.get_number_of_comparisons();
                    ++tree_idx) {
                std::vector<unsigned int> other_height_indices = comparisons.get_other_height_indices(tree_idx);
                std::vector<double> ln_category_probs(other_height_indices.size() + comparisons.get_number_of_auxiliary_heights());

                for (auto height_idx : other_height_indices) {
                    unsigned int number_of_elements = comparisons.get_number_of_trees_mapped_to_height_index(height_idx);
                    if (height_idx == tree_idx) {
                        --number_of_elements;
                    }
                    comparisons.trees_.at(tree_idx).set_height(comparisons.node_heights_.at(height_idx).get_value());
                    double lnl = comparisons.trees_.at(tree_idx).compute_log_likelihood();
                    ln_category_probs.push_back(lnl + std::log(number_of_elements));
                }

                std::vector<double> auxiliary_heights(comparisons.get_number_of_auxiliary_heights());
                for (unsigned int i = 0; i < comparisons.get_number_of_auxiliary_heights(); ++i) {
                    double fresh_height = comparisons.node_height_prior.draw(rng);
                    auxiliary_heights.push_back(fresh_height);
                    comparisons.trees_.at(tree_idx).set_height(fresh_height);
                    double lnl = comparisons.trees_.at(tree_idx).compute_log_likelihood();
                    ln_category_probs.push_back(lnl + ln_concentration_over_num_aux);
                }

                normalize(ln_category_probs);



            }
        }
};


#endif
