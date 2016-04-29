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
#include "math_util.hpp"


class OperatorSchedule;

//////////////////////////////////////////////////////////////////////////////
// Operator base classes
//////////////////////////////////////////////////////////////////////////////

class Operator {
    public:
        Operator() { }
        Operator(double weight);
        virtual ~Operator() { }
		enum OperatorTypeEnum {
            tree_operator = 1,
            time_operator = 2,
            model_operator = 3
        };

        virtual Operator::OperatorTypeEnum get_type() const = 0;

        virtual void optimize(double log_alpha) = 0;

        double get_target_acceptance_probability() const;

        void set_operator_schedule(std::shared_ptr<OperatorSchedule> os);

        virtual double get_coercable_parameter_value();

        virtual void set_coercable_parameter_value(double value) { }

        double get_weight() const { return this->weight_; }

        void set_weight(double weight);

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const = 0;

        void accept();

        void reject();

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

        virtual std::string get_name() const;

        virtual std::string target_parameter() const = 0;

        std::string header_string() const;

        std::string to_string() const;

    protected:
        std::shared_ptr<OperatorSchedule> operator_schedule_;
        double weight_ = 1.0;
        unsigned int number_rejected_ = 0;
        unsigned int number_accepted_ = 0;
        unsigned int number_rejected_for_correction = 0;
        unsigned int number_accepted_for_correction = 0;

        double calc_delta(double log_alpha) const;
};

class ScaleOperator : public Operator {
    protected:
        double scale_ = 0.5;

    public:
        ScaleOperator() : Operator() { }
        ScaleOperator(double weight) : Operator(weight) { }
        ScaleOperator(double weight, double scale) : Operator(weight);
        virtual ~ScaleOperator() { }

        void set_scale(double scale);

        double get_scale() const {
            return this->scale_;
        }

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const;

        void optimize(double log_alpha);

        double get_coercable_parameter_value();

        void set_coercable_parameter_value(double value);

        std::string get_name() const;
};

class WindowOperator : public Operator {
    protected:
        double window_size_ = 0.1;

    public:
        WindowOperator() : Operator() { }
        WindowOperator(double weight) : Operator(weight) { }
        WindowOperator(double weight, double window_size) : Operator(weight) {
            this->set_window_size(window_size);
        }
        virtual ~WindowOperator() { }

        void set_window_size(double window_size) {
            ECOEVOLITY_ASSERT(window_size > 0.0);
            this->window_size_ = window_size;
        }
        double get_window_size() const {
            return this->window_size_;
        }

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const {
            double addend = (rng.uniform_real() * 2 * this->window_size_) - this->window_size;
            parameter_value += addend
            hastings_ratio = 0.0;
        }

        void optimize(double log_alpha) {
            double delta = this->calc_delta(log_alpha);
            delta += std::log(this->window_size_);
            this->set_window_size(std::exp(delta));
        }

        double get_coercable_parameter_value() {
            return this->window_size_;
        }

        void set_coercable_parameter_value(double value) {
            this->set_window_size(value);
        }

        std::string get_name() const {
            return "BaseWindowOperator";
        }
};


//////////////////////////////////////////////////////////////////////////////
// Operator Type base classes
//////////////////////////////////////////////////////////////////////////////

class ModelOperator : public Operator {
    public:
        ModelOperator() : Operator() { }
        ModelOperator(double weight) : Operator(weight) { }
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

class ComparisonTreeScaleOperator : public ScaleOperator {
    public:
        ComparisonTreeScaleOperator() : ScaleOperator() { }
        ComparisonTreeScaleOperator(double weight) : ScaleOperator(weight) { }
        ComparisonTreeScaleOperator(double weight, double scale) : ScaleOperator(weight, scale) { }
        virtual ~ComparisonTreeScaleOperator() { }

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
            return "ComparisonTreeScaleOperator";
        }
};

class ComparisonTreeWindowOperator : public WindowOperator {
    public:
        ComparisonTreeWindowOperator() : WindowOperator() { }
        ComparisonTreeWindowOperator(double weight) : WindowOperator(weight) { }
        ComparisonTreeWindowOperator(double weight, double window_size) : WindowOperator(weight, window_size) { }
        virtual ~ComparisonTreeWindowOperator() { }

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
            return "ComparisonTreeWindowOperator";
        }
};

class NodeHeightScaleOperator : public ScaleOperator {
    public:
        NodeHeightScaleOperator() : ScaleOperator() { }
        NodeHeightScaleOperator(double weight) : ScaleOperator(weight) { }
        NodeHeightScaleOperator(double weight, double scale) : ScaleOperator(weight, scale) { }
        virtual ~NodeHeightScaleOperator() { }

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
            return "NodeHeightScaleOperator";
        }

        std::string target_parameter() const {;
            return "node height";
        }
};

class NodeHeightWindowOperator : public WindowOperator {
    public:
        NodeHeightWindowOperator() : WindowOperator() { }
        NodeHeightWindowOperator(double weight) : WindowOperator(weight) { }
        NodeHeightWindowOperator(double weight, double window_size) : WindowOperator(weight, window_size) { }
        virtual ~NodeHeightWindowOperator() { }

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
            return "NodeHeightWindowOperator";
        }

        std::string target_parameter() const {;
            return "node height";
        }
};


//////////////////////////////////////////////////////////////////////////////
// Derived Operator classes
//////////////////////////////////////////////////////////////////////////////

class ConcentrationScaler: public ScaleOperator {
    public:
        ConcentrationScaler() : ScaleOperator() { }
        ConcentrationScaler(double weight) : ScaleOperator(weight) { }
        ConcentrationScaler(double weight, double scale) : ScaleOperator(weight, scale) { }
        virtual ~ConcentrationScaler() { }

        Operator::OperatorTypeEnum get_type() const {
            return Operator::OperatorTypeEnum::model_operator;
        }

        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons) {
            double v = comparisons->get_concentration();
            double hastings;
            this->update(rng, v, hastings);
            comparisons->set_concentration(v);
            return hastings;
        }

        std::string target_parameter() const {;
            return "concentration";
        }

        std::string get_name() const {
            return "ConcentrationScaler";
        }
};

class MutationRateMover : public ComparisonTreeWindowOperator {
    public:
        MutationRateMover() : ComparisonTreeWindowOperator() { }
        MutationRateMover(double weight) : ComparisonTreeWindowOperator(weight) { }
        MutationRateMover(double weight, double window_size) : ComparisonTreeWindowOperator(weight, window_size) { }
        virtual ~MutationRateMover() { }

        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) {
            double red_freq = tree->get_v();
            double hastings;
            this->update(rng, red_freq, hastings);
            if ((red_freq >= 0.0) && (red_freq <= 1.0)) {
                // v is also set here
                this->set_u(1.0 / (2.0 * red_freq));
                return hastings; 
            }
            return -std::numeric_limits<double>::infinity();
        }

        std::string target_parameter() const {;
            return "mutation rate";
        }

        std::string get_name() const {
            return "MutationRateMover";
        }
};

class ComparisonHeightMultiplierScaler : public ComparisonTreeScaleOperator {
    public:
        ComparisonHeightMultiplierScaler() : ComparisonTreeScaleOperator() { }
        ComparisonHeightMultiplierScaler(double weight) : ComparisonTreeScaleOperator(weight) { }
        ComparisonHeightMultiplierScaler(double weight, double scale) : ComparisonTreeScaleOperator(weight, scale) { }
        virtual ~ComparisonHeightMultiplierScaler() { }

        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) {
            double v = tree->get_node_height_multiplier();
            double hastings;
            this->update(rng, v, hastings);
            tree->set_node_height_multiplier(v);
            return hastings;
        }

        std::string target_parameter() const {;
            return "node height multiplier";
        }

        std::string get_name() const {
            return "ComparisonHeightMultiplierScaler";
        }
};

class ChildCoalescenceRateScaler : public ComparisonTreeScaleOperator {
    public:
        ChildCoalescenceRateScaler() : ComparisonTreeScaleOperator() { }
        ChildCoalescenceRateScaler(double weight) : ComparisonTreeScaleOperator(weight) { }
        ChildCoalescenceRateScaler(double weight, double scale) : ComparisonTreeScaleOperator(weight, scale) { }
        virtual ~ChildCoalescenceRateScaler() { }

        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) {
            int pop_idx = rng.uniform_int(0, tree->get_leaf_allele_count() - 1);
            double rate = tree->get_child_coalescence_rate(pop_idx);

            double hastings;
            this->update(rng, rate, hastings);

            tree->set_child_coalescence_rate(pop_idx, rate);
            return hastings;
        }

        std::string target_parameter() const {;
            return "coalescence rate";
        }

        std::string get_name() const {
            return "ChildCoalescenceRateScaler";
        }
};

class RootCoalescenceRateScaler : public ChildCoalescenceRateScaler {
        RootCoalescenceRateScaler() : ChildCoalescenceRateScaler() { }
        RootCoalescenceRateScaler(double weight) : ChildCoalescenceRateScaler(weight) { }
        RootCoalescenceRateScaler(double weight, double scale) : ChildCoalescenceRateScaler(weight, scale) { }
        virtual ~RootCoalescenceRateScaler() { }

        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) {
            double rate = tree->get_root_coalescence_rate();

            double hastings;
            this->update(rng, rate, hastings);

            tree->set_root_coalescence_rate(rate);

            return hastings;
        }

        std::string get_name() const {
            return "RootCoalescenceRateScaler";
        }
};

class ComparisonHeightScaler : public NodeHeightScaleOperator {
    public:
        ComparisonHeightScaler() : NodeHeightScaleOperator() { }
        ComparisonHeightScaler(double weight) : NodeHeightScaleOperator(weight) { }
        ComparisonHeightScaler(double weight, double scale) : NodeHeightScaleOperator(weight, scale) {
            this->set_scale(scale);
        }
        virtual ~ComparisonHeightScaler() { }

        double propose(
                RandomNumberGenerator& rng,
                PositiveRealParameter& node_height) {
            double h = node_height->get_value();
            double hastings;
            this->update(rng, h, hastings);
            node_height.set_value(h);
            return hastings;
        }

        std::string get_name() const {
            return "ComparisonHeightScaler";
        }
};

class DirichletProcessGibbsSampler : public ModelOperator {
    public:
        DirichletProcessGibbsSampler() : ModelOperator() { }
        DirichletProcessGibbsSampler(double weight) : ModelOperator(weight) { }
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
                    tree_idx < comparisons.get_number_of_trees();
                    ++tree_idx) {
                std::vector<unsigned int> other_height_indices = comparisons.get_other_height_indices(tree_idx);
                std::vector<double> ln_category_likelihoods;
                ln_category_likelihoods.reserve(other_height_indices.size() +
                        comparisons.get_number_of_auxiliary_heights());

                // store height associated with this tree
                comparisons.node_heights_.at(comparisons.get_height_index(tree_idx)).store();
                for (auto height_idx : other_height_indices) {
                    unsigned int number_of_elements = comparisons.get_number_of_trees_mapped_to_height(height_idx);
                    if (height_idx == comparisons.get_height_index(tree_idx)) {
                        --number_of_elements;
                        if (! comparisons.trees_.at(tree_idx).is_dirty()) {
                            double lnl = comparisons.trees_.at(tree_idx).get_log_likelihood_value();
                            ln_category_likelihoods.push_back(lnl + std::log(number_of_elements));
                            continue;
                        }
                    }
                    else {
                        comparisons.trees_.at(tree_idx).set_height(comparisons.node_heights_.at(height_idx).get_value());
                    }
                    double lnl = comparisons.trees_.at(tree_idx).compute_log_likelihood();
                    ln_category_likelihoods.push_back(lnl + std::log(number_of_elements));
                }

                std::vector<double> auxiliary_heights;
                auxiliary_heights.reserve(comparisons.get_number_of_auxiliary_heights());
                for (unsigned int i = 0; i < comparisons.get_number_of_auxiliary_heights(); ++i) {
                    double fresh_height = comparisons.node_height_prior_.draw(rng);
                    auxiliary_heights.push_back(fresh_height);
                    comparisons.trees_.at(tree_idx).set_height(fresh_height);
                    double lnl = comparisons.trees_.at(tree_idx).compute_log_likelihood();
                    ln_category_likelihoods.push_back(lnl + ln_concentration_over_num_aux);
                }

                // restore height associated with this tree
                comparisons.node_heights_.at(comparisons.get_height_index(tree_idx)).restore();

                std::vector<double> category_probs(ln_category_likelihoods);
                normalize_log_likelihoods(category_probs);
                unsigned int prob_index = rng.weighted_index(category_probs);
                if (prob_index < other_height_indices.size()) {
                    comparisons.remap_tree(
                            tree_idx,
                            other_height_indices.at(prob_idx),
                            ln_category_likelihoods.at(prob_idx));
                }
                else {
                    comparisons.map_tree_to_new_height(
                            tree_idx,
                            auxiliary_heights.at(new_height_index -
                                    other_height_indices.size()),
                            ln_category_likelihoods.at(prob_idx));
                }
            }
        }
};

class ReversibleJumpSampler : public ModelOperator {
    public:
        ReversibleJumpSampler() : ModelOperator() { }
        ReversibleJumpSampler(double weight) : ModelOperator(weight) { }
        virtual ~ReversibleJumpSampler() { }

        std::string get_name() const {
            return "ReversibleJumpSampler";
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons) {
            throw EcoevolityNotImplementedError("rjMCMC not implemented yet");
        }
};


//////////////////////////////////////////////////////////////////////////////
// Operator schedule 
//////////////////////////////////////////////////////////////////////////////

class OperatorSchedule {
    protected:
        std::vector< std::shared_ptr<Operator> > operators_;
        double total_weight_ = 0.0;
        std::vector<double> cumulative_probs_;
        unsigned int auto_optimize_delay_ = 10000;
        unsigned int auto_optimize_delay_count_ = 0;
        bool auto_optimize_ = true;

    public:
        OperatorSchedule() { }
        virtual ~OperatorSchedule() { }

        void add_operator(std::shared_ptr<Operator> o) {
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

        double get_total_weight() const {
            return this->total_weight_;
        }
        unsigned int get_auto_optimize_delay_count() const {
            return this->auto_optimize_delay_count_;
        }
        unsigned int get_auto_optimize_delay() const {
            return this->auto_optimize_delay_;
        }
        void set_auto_optimize_delay(unsigned int delay) {
            this->auto_optimize_delay_ = delay;
        }

        void write_operator_rates(std::ofstream out) {
            const Operator& op = this->operators_.at(0);
            out << op.header_string();
            out << op.to_string();
            for (unsigned int i = 1; i < this->operators_.size(); ++i) {
                out << this->operators_.at(i).to_string();
            }
        }

        bool auto_optimizing() const {
            return this->auto_optimize_;
        }
        void turn_on_auto_optimize() {
            this->auto_optimize_ = true;
        }
        void turn_off_auto_optimize() {
            this->auto_optimize_ = false;
        }
};

#endif
