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

#ifndef ECOEVOLITY_GENERAL_TREE_OPERATOR_HPP
#define ECOEVOLITY_GENERAL_TREE_OPERATOR_HPP

#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>
#include <memory>

#include "basetree.hpp"
#include "tree.hpp"
#include "rng.hpp"
#include "assert.hpp"
#include "math_util.hpp"

class BaseGeneralTreeOperatorTemplate {
    protected:
        double weight_ = 1.0;
        bool ignore_proposal_attempt_ = false;

    public:
		enum OperatorTypeEnum {
            topology_operator                   = 1,
            topology_model_operator             = 2,
            node_height_operator                = 3,
            root_height_operator                = 4,
            global_height_operator              = 5,
            node_height_prior_operator          = 6,
            derived_operator                    = 7,
        };
		enum OperatorScopeEnum {
            topology                            = 1,
            global                              = 2,
            node_height                         = 3,
            root_height                         = 4,
            branch                              = 5,
            hyper                               = 6,
        };

        BaseGeneralTreeOperatorTemplate() { }
        BaseGeneralTreeOperatorTemplate(double weight) {
            this->set_weight(weight);
        }

        virtual BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const = 0;
        virtual BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const = 0;

        double get_weight() const {
            return this->weight_;
        }

        void set_weight(double weight) {
            ECOEVOLITY_ASSERT(weight >= 0.0);
            this->weight_ = weight;
        }
        virtual void set_default_weight(unsigned int number_of_leaves) {
            if (this->get_type() == BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::topology_model_operator) {
                this->weight_ = (double)number_of_leaves;
            }
            else if (this->get_scope() == BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::topology) {
                this->weight_ = (double)number_of_leaves / 2.0;
            }
            else if (this->get_scope() == BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::node_height) {
                this->weight_ = (double)number_of_leaves / 2.0;
            }
            else if (this->get_scope() == BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::branch) {
                this->weight_ = (double)number_of_leaves / 2.0;
            }
            else {
                this->weight_ = 1.0;
            }
        }

        virtual std::string get_name() const = 0;

        virtual std::string target_parameter() const = 0;

        virtual void optimize(double log_alpha) = 0;

        virtual void accept() = 0;

        virtual void reject() = 0;

        virtual double get_coercable_parameter_value() const = 0;

        virtual void set_coercable_parameter_value(double value) = 0;

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const = 0;

        virtual double calc_delta(double log_alpha) const = 0;

        virtual double get_target_acceptance_probability() const = 0;

        virtual unsigned int get_number_rejected() const = 0;
        virtual unsigned int get_number_accepted() const = 0;
        virtual unsigned int get_number_rejected_for_correction() const = 0;
        virtual unsigned int get_number_accepted_for_correction() const = 0;

        virtual unsigned int get_number_of_attempts() const = 0;
        virtual unsigned int get_number_of_attempts_for_correction() const = 0;

        virtual unsigned int get_auto_optimize_delay() const = 0;
        virtual void set_auto_optimize_delay(unsigned int delay) = 0;
        virtual bool auto_optimizing() const = 0;
        virtual void turn_on_auto_optimize() = 0;
        virtual void turn_off_auto_optimize() = 0;

        virtual std::string header_string() const = 0;

        virtual std::string to_string() const = 0;
};

template<class TreeType>
class GeneralTreeOperatorTemplate : public BaseGeneralTreeOperatorTemplate {

    public:
        GeneralTreeOperatorTemplate() : BaseGeneralTreeOperatorTemplate() { }
        GeneralTreeOperatorTemplate(double weight) : BaseGeneralTreeOperatorTemplate(weight) { }

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate< TreeType > > > helper_ops;

        virtual void call_store_methods(
                TreeType * tree) const {
            tree->store_state();
        }
        virtual void call_restore_methods(
                TreeType * tree) const {
            tree->restore_state();
        }

        virtual void perform_move(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) = 0;

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) = 0;

        virtual bool is_operable(TreeType * tree) const { return true; }
        virtual void operate(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1,
                unsigned int number_of_moves = 1) = 0;

        virtual void operate_plus(RandomNumberGenerator& rng,
                TreeType * tree,
                std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > other_operators,
                unsigned int nthreads = 1,
                unsigned int number_of_moves = 1,
                unsigned int other_op_number_of_moves = 1) = 0;
        virtual void operate_with_helpers(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1,
                unsigned int number_of_moves = 1,
                unsigned int helper_op_number_of_moves = 1) = 0;
};


template<class TreeType, class OperatorType>
class GeneralTreeOperatorInterface : public GeneralTreeOperatorTemplate<TreeType> {

    public:
        OperatorType op_;

        GeneralTreeOperatorInterface() : GeneralTreeOperatorTemplate<TreeType>() { }
        GeneralTreeOperatorInterface(double weight) : GeneralTreeOperatorTemplate<TreeType>(weight) { }
        GeneralTreeOperatorInterface(
                double weight,
                double tuning_parameter
                ) : GeneralTreeOperatorTemplate<TreeType>(weight) {
            this->op_.set_coercable_parameter_value(tuning_parameter);
        }

        void perform_move(
                RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            this->call_store_methods(tree);

            // std::cout << "lnl before move: " << tree->get_log_likelihood_value() << "\n";
            // std::cout << "Calling propose on " << this->get_name() << "\n";
        
            double hastings_ratio = this->propose(rng, tree, nthreads);

            // Debug check for any negative branch lengths
            // std::string t = tree->to_parentheses(false);
            // if (t.find(":-") != std::string::npos) {
            //     std::cerr << t << "\n";
            //     throw EcoevolityError("neg branch length!");
            // }

            // If the Hasting's ratio is zero (i.e., the log is -inf), we can
            // reject before going any further and wasting computation on the
            // likelihood
            if (hastings_ratio == -std::numeric_limits<double>::infinity()) {
                if (! this->ignore_proposal_attempt_) {
                    this->reject();
                }
                this->call_restore_methods(tree);
                tree->make_clean();
                if (this->auto_optimizing() && (! this->ignore_proposal_attempt_)) {
                    this->optimize(hastings_ratio);
                }
                // std::cout << "bad move; ignored: " << this->ignore_proposal_attempt_ << "\n";
                this->ignore_proposal_attempt_ = false;
                return;
            }

            tree->compute_log_likelihood_and_prior(true);

            // std::cout << "lnl after comp: " << tree->get_log_likelihood_value() << "\n";
        
            double likelihood_ratio = 
                    tree->get_log_likelihood_value() -
                    tree->get_stored_log_likelihood_value();
            double prior_ratio = 
                    tree->get_log_prior_density_value() -
                    tree->get_stored_log_prior_density_value();
            double acceptance_probability =
                    likelihood_ratio + 
                    prior_ratio +
                    hastings_ratio;
            // std::cout << "ln(like ratio) = " << likelihood_ratio << "\n";
            // std::cout << "ln(prior ratio) = " << prior_ratio << "\n";
            // std::cout << "ln(hastings ratio) = " << hastings_ratio << "\n";
            // std::cout << "ln(p(accept)) = " << acceptance_probability << "\n";
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_probability)) {
                if (! this->ignore_proposal_attempt_) {
                    this->accept();
                }
                // std::cout << "ACCEPT!\n";
            }
            else {
                if (! this->ignore_proposal_attempt_) {
                    this->reject();
                }
                // std::cout << "REJECT!\n";
                this->call_restore_methods(tree);
            }
            tree->make_clean();
            if (this->auto_optimizing() && (! this->ignore_proposal_attempt_)) {
                this->optimize(acceptance_probability);
            }
            this->ignore_proposal_attempt_ = false;
            // std::cout << "lnl at end: " << tree->get_log_likelihood_value() << "\n";
        }

        void operate(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1,
                unsigned int number_of_moves = 1) {
            for (unsigned int i = 0; i < number_of_moves; ++i) {
                this->perform_move(rng, tree, nthreads);
            }
        }

        virtual void operate_plus(RandomNumberGenerator& rng,
                TreeType * tree,
                std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > other_operators,
                unsigned int nthreads = 1,
                unsigned int number_of_moves = 1,
                unsigned int other_op_number_of_moves = 1) {
            this->operate(rng, tree, nthreads, number_of_moves);
        }

        void operate_with_helpers(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1,
                unsigned int number_of_moves = 1,
                unsigned int helper_op_number_of_moves = 1) {
            for (unsigned int i = 0; i < number_of_moves; ++i) {
                // std::cout << "\ncalling: " << this->get_name() << "\n";
                this->perform_move(rng, tree, nthreads);
                for (auto helper_op : this->helper_ops) {
                    // std::cout << "\nhelper call: " << helper_op->get_name() << "\n";
                    helper_op->operate(rng, tree, nthreads, helper_op_number_of_moves);
                }
            }
        }

        virtual void optimize(double log_alpha) {
            this->op_.optimize(log_alpha);
        }

        void accept() {
            this->op_.accept();
        }

        void reject() {
            this->op_.reject();
        }

        virtual double get_coercable_parameter_value() const {
            return this->op_.get_coercable_parameter_value();
        }

        virtual void set_coercable_parameter_value(double value) {
            this->op_.set_coercable_parameter_value(value);
        }

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const {
            this->op_.update(rng, parameter_value, hastings_ratio);
        }

        double calc_delta(double log_alpha) const {
            return this->op_.calc_delta(log_alpha);
        }

        double get_target_acceptance_probability() const {
            return this->op_.get_target_acceptance_probability();
        }

        unsigned int get_number_rejected() const {
            return this->op_.get_number_rejected();
        }
        unsigned int get_number_accepted() const {
            return this->op_.get_number_accepted();
        }
        unsigned int get_number_rejected_for_correction() const {
            return this->op_.get_number_rejected_for_correction();
        }
        unsigned int get_number_accepted_for_correction() const {
            return this->op_.get_number_accepted_for_correction();
        }
        unsigned int get_number_of_attempts() const {
            return this->op_.get_number_of_attempts();
        }
        unsigned int get_number_of_attempts_for_correction() const {
            return this->op_.get_number_of_attempts_for_correction();
        }

        unsigned int get_auto_optimize_delay() const {
            return this->op_.get_auto_optimize_delay();
        }
        void set_auto_optimize_delay(unsigned int delay) {
            this->op_.set_auto_optimize_delay(delay);
        }
        bool auto_optimizing() const {
            return this->op_.auto_optimizing();
        }
        void turn_on_auto_optimize() {
            this->op_.turn_on_auto_optimize();
        }
        void turn_off_auto_optimize() {
            this->op_.turn_off_auto_optimize();
        }

        std::string header_string() const {
            return "name\tprop_accepted\tnumber_accepted\tnumber_rejected\tweight\ttuning_parameter\n";
        }

        virtual std::string to_string() const {
            unsigned int n_accepted = this->op_.get_number_accepted();
            unsigned int n_rejected = this->op_.get_number_rejected();
            double p_accepted = (double)n_accepted / (n_accepted + n_rejected);
            std::ostringstream ss;
            ss << this->get_name() << "\t"
               << p_accepted << "\t"
               << n_accepted << "\t"
               << n_rejected << "\t"
               << this->get_weight() << "\t";

            double tuning = this->op_.get_coercable_parameter_value();
            if (std::isnan(tuning)) {
                ss << "none\t";
            }
            else {
                ss << tuning << "\t";
            }
            ss << "\n";
            return ss.str();
        }
};


//////////////////////////////////////////////////////////////////////////////
// Operator base classes
//////////////////////////////////////////////////////////////////////////////

class Op {

    public:
        Op(const bool auto_optimize = false) :
            auto_optimize_(auto_optimize) { }
        Op(const unsigned int auto_optimize_delay,
                const bool auto_optimize = false) :
            auto_optimize_delay_(auto_optimize_delay),
            auto_optimize_(auto_optimize) { }

        virtual void optimize(double log_alpha) { }

        virtual double get_coercable_parameter_value() const {
            return std::numeric_limits<double>::quiet_NaN();
        }

        virtual void set_coercable_parameter_value(double value) { }

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const { }

        virtual double get_move_amount(RandomNumberGenerator& rng) const { return 0.0; }

        double calc_delta(double log_alpha) const {
            if ((this->get_number_of_attempts() <= this->get_auto_optimize_delay()) ||
                    (! this->auto_optimize_)) {
                return 0.0;
            }
            double target = this->get_target_acceptance_probability();
            double count = (this->get_number_rejected_for_correction() +
                            this->get_number_accepted_for_correction() +
                            1.0);
            double delta_p = ((1.0 / count) * (std::exp(std::min(log_alpha, 0.0)) - target));
            double mx = std::numeric_limits<double>::max();
            if ((delta_p > -mx) && (delta_p < mx)) {
                return delta_p;
            }
            return 0.0;
        }

        void accept() {
            ++this->number_accepted_;
            if (this->get_number_of_attempts() > this->get_auto_optimize_delay()) {
                ++this->number_accepted_for_correction_;
            }
        }

        void reject() {
            ++this->number_rejected_;
            if (this->get_number_of_attempts() > this->get_auto_optimize_delay()) {
                ++this->number_rejected_for_correction_;
            }
        }

        double get_target_acceptance_probability() const {
            // Some prelim tests confirm that an acceptance rate of 0.44 leads to
            // better mixing for the simple, univariate random variables that the
            // operators are updating.
            // return 0.234;
            return 0.44;
        }

        unsigned int get_number_rejected() const {
            return this->number_rejected_;
        }
        unsigned int get_number_accepted() const {
            return this->number_accepted_;
        }
        unsigned int get_number_rejected_for_correction() const {
            return this->number_rejected_for_correction_;
        }
        unsigned int get_number_accepted_for_correction() const {
            return this->number_accepted_for_correction_;
        }
        unsigned int get_number_of_attempts() const {
            return this->number_rejected_ + this->number_accepted_;
        }
        unsigned int get_number_of_attempts_for_correction() const {
            return this->number_rejected_for_correction_ + this->number_accepted_for_correction_;
        }

        unsigned int get_auto_optimize_delay() const {
            return this->auto_optimize_delay_;
        }
        void set_auto_optimize_delay(unsigned int delay) {
            this->auto_optimize_delay_ = delay;
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

    protected:
        unsigned int number_rejected_ = 0;
        unsigned int number_accepted_ = 0;
        unsigned int number_rejected_for_correction_ = 0;
        unsigned int number_accepted_for_correction_ = 0;
        unsigned int auto_optimize_delay_ = 1000;
        bool auto_optimize_ = true;
};


class BaseOptimizingOp : public Op {

    public:
        BaseOptimizingOp(const bool auto_optimize = true) :
            Op(auto_optimize) { }
        BaseOptimizingOp(const unsigned int auto_optimize_delay,
                const bool auto_optimize = true) :
            Op(auto_optimize_delay, auto_optimize) { }

        void optimize(double log_alpha) {
            double delta = this->calc_delta(log_alpha);
            if (delta == 0.0) {
                return;
            }
            delta += std::log(this->get_coercable_parameter_value());
            this->set_coercable_parameter_value(std::exp(delta));
        }
};


class ScaleOp : public BaseOptimizingOp {

    protected:
        double scale_ = 0.5;

    public:
        ScaleOp(const bool auto_optimize = true) : BaseOptimizingOp(auto_optimize) { }
        ScaleOp(const unsigned int auto_optimize_delay,
                const bool auto_optimize = true) :
                BaseOptimizingOp(auto_optimize_delay, auto_optimize) { }
        ScaleOp(double scale,
                const bool auto_optimize = true) :
                BaseOptimizingOp(auto_optimize) {
            this->set_scale(scale);
        }
        ScaleOp(double scale,
                const unsigned int auto_optimize_delay,
                const bool auto_optimize = true) :
                BaseOptimizingOp(auto_optimize_delay, auto_optimize) {
            this->set_scale(scale);
        }

        void set_scale(double scale) {
            ECOEVOLITY_ASSERT(scale > 0.0);
            this->scale_ = scale;
        }

        double get_scale() const {
            return this->scale_;
        }

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const {
            double multiplier = this->get_move_amount(rng);
            parameter_value *= multiplier;
            hastings_ratio = std::log(multiplier);
        }

        double get_move_amount(RandomNumberGenerator& rng) const {
            return std::exp(this->scale_ * ((2.0 * rng.uniform_real()) - 1.0));
        }

        double get_coercable_parameter_value() const {
            return this->scale_;
        }

        void set_coercable_parameter_value(double value) {
            this->set_scale(value);
        }
};


class WindowOp : public BaseOptimizingOp {

    protected:
        double window_size_ = 0.1;

    public:
        WindowOp(const bool auto_optimize = true) : BaseOptimizingOp(auto_optimize) { }
        WindowOp(const unsigned int auto_optimize_delay,
                const bool auto_optimize = true) :
                BaseOptimizingOp(auto_optimize_delay, auto_optimize) { }
        WindowOp(double window_size,
                const bool auto_optimize = true) :
                BaseOptimizingOp(auto_optimize) {
            this->set_window_size(window_size);
        }
        WindowOp(double window_size,
                const unsigned int auto_optimize_delay,
                const bool auto_optimize = true) :
                BaseOptimizingOp(auto_optimize_delay, auto_optimize) {
            this->set_window_size(window_size);
        }

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
            double addend = this->get_move_amount(rng);
            parameter_value += addend;
            hastings_ratio = 0.0;
        }

        double get_move_amount(RandomNumberGenerator& rng) const {
            return (rng.uniform_real() * 2 * this->window_size_) - this->window_size_;
        }

        double get_coercable_parameter_value() const {
            return this->window_size_;
        }

        void set_coercable_parameter_value(double value) {
            this->set_window_size(value);
        }
};


class DirichletOp : public BaseOptimizingOp {

    protected:
        double scale_ = 0.5;

    public:
        DirichletOp(const bool auto_optimize = true) : BaseOptimizingOp(auto_optimize) { }
        DirichletOp(const unsigned int auto_optimize_delay,
                const bool auto_optimize = true) :
                BaseOptimizingOp(auto_optimize_delay, auto_optimize) { }
        DirichletOp(double scale,
                const bool auto_optimize = true) :
                BaseOptimizingOp(auto_optimize) {
            this->set_scale(scale);
        }
        DirichletOp(double scale,
                const unsigned int auto_optimize_delay,
                const bool auto_optimize = true) :
                BaseOptimizingOp(auto_optimize_delay, auto_optimize) {
            this->set_scale(scale);
        }

        void set_scale(double scale) {
            ECOEVOLITY_ASSERT(scale > 0.0);
            this->scale_ = scale;
        }

        double get_scale() const {
            return this->scale_;
        }

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const {
            ECOEVOLITY_ASSERT(parameter_value <= 1.0);
            std::vector<double> values { parameter_value, 1.0 - parameter_value };
            this->update_vector(rng, values, hastings_ratio);
            parameter_value = values.at(0);
        }

        void update_vector(
                RandomNumberGenerator& rng,
                std::vector<double> & parameter_values,
                double& hastings_ratio) const {
            ECOEVOLITY_ASSERT(parameter_values.size() > 1);
            double sum = std::accumulate(parameter_values.begin(), parameter_values.end(), 0.0);
            ECOEVOLITY_ASSERT(almost_equal(sum, 1.0));
        
            std::vector<double> old_proportions = parameter_values;
            std::vector<double> forward_dir_parameters = old_proportions;
            for (unsigned int i = 0; i < forward_dir_parameters.size(); ++i) {
                forward_dir_parameters.at(i) = 1.0 + (forward_dir_parameters.at(i) * (1.0 / this->get_scale()));
            }
            DirichletDistribution dir_forward = DirichletDistribution(forward_dir_parameters);
            std::vector<double> new_proportions = dir_forward.draw(rng);
        
            for (unsigned int i = 0; i < new_proportions.size(); ++i) {
                parameter_values.at(i) = new_proportions.at(i);
            }
        
            std::vector<double> reverse_dir_parameters = new_proportions;
            for (unsigned int i = 0; i < reverse_dir_parameters.size(); ++i) {
                reverse_dir_parameters.at(i) = 1.0 + (reverse_dir_parameters.at(i) * (1.0 / this->get_scale()));
            }
            DirichletDistribution dir_reverse = DirichletDistribution(reverse_dir_parameters);
        
            hastings_ratio = (dir_reverse.ln_pdf(old_proportions) -
                    dir_forward.ln_pdf(new_proportions));
        }

        double get_move_amount(RandomNumberGenerator& rng) const {
            // We should never call this
            ECOEVOLITY_ASSERT(0 == 1);
            return 0.0;
        }

        double get_coercable_parameter_value() const {
            return this->scale_;
        }

        void set_coercable_parameter_value(double value) {
            this->set_scale(value);
        }
};


class BetaOp : public BaseOptimizingOp {

    protected:
        double scale_ = 0.5;

    public:
        BetaOp(const bool auto_optimize = true) : BaseOptimizingOp(auto_optimize) { }
        BetaOp(const unsigned int auto_optimize_delay,
                const bool auto_optimize = true) :
                BaseOptimizingOp(auto_optimize_delay, auto_optimize) { }
        BetaOp(double scale,
                const bool auto_optimize = true) :
                BaseOptimizingOp(auto_optimize) {
            this->set_scale(scale);
        }
        BetaOp(double scale,
                const unsigned int auto_optimize_delay,
                const bool auto_optimize = true) :
                BaseOptimizingOp(auto_optimize_delay, auto_optimize) {
            this->set_scale(scale);
        }

        void set_scale(double scale) {
            ECOEVOLITY_ASSERT(scale > 0.0);
            this->scale_ = scale;
        }

        double get_scale() const {
            return this->scale_;
        }

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const {
            ECOEVOLITY_ASSERT(parameter_value <= 1.0);

            // We are going to propose a random value from a beta distribution
            // centered on the current value (v). The mean of a beta is
            //
            //   v = alpha / (alpha + beta)
            //
            // The larger alpha and beta are, the smaller moves we will
            // propose, so we will set
            //
            //   alpha + beta = 1 / this->scale_
            // 
            // So that when the scale gets larger, we will propose larger moves.
            // The alpha and beta parameters of the beta distributions we will
            // draw from are then:
            //
            //   alpha = v * (alpha + beta)
            //   beta = (alpha + beta) - alpha
            double alpha_plus_beta = 1.0 / this->scale_;
            double alpha = parameter_value * alpha_plus_beta;
            double beta = alpha_plus_beta - alpha;
            double new_value = BetaDistribution::get_draw(rng, alpha, beta);
            double ln_prob_forward_move = BetaDistribution::get_ln_pdf(
                    new_value,
                    alpha,
                    beta);

            // Now we need prob of reverse move
            double rev_alpha = new_value * alpha_plus_beta;
            double rev_beta = alpha_plus_beta - rev_alpha;
            double ln_prob_reverse_move = BetaDistribution::get_ln_pdf(
                    parameter_value,
                    rev_alpha,
                    rev_beta);
            hastings_ratio = ln_prob_reverse_move - ln_prob_forward_move;
            parameter_value = new_value;

        }

        double get_move_amount(RandomNumberGenerator& rng) const {
            // We should never call this
            ECOEVOLITY_ASSERT(0 == 1);
            return 0.0;
        }

        double get_coercable_parameter_value() const {
            return this->scale_;
        }

        void set_coercable_parameter_value(double value) {
            this->set_scale(value);
        }
};


//////////////////////////////////////////////////////////////////////////////
// BaseTree operators 
//////////////////////////////////////////////////////////////////////////////

template<class TreeType>
class NodeHeightPriorAlphaScaler : public GeneralTreeOperatorInterface<TreeType, ScaleOp> {
    public:
        NodeHeightPriorAlphaScaler() : GeneralTreeOperatorInterface<TreeType, ScaleOp>() { }
        NodeHeightPriorAlphaScaler(double weight) : GeneralTreeOperatorInterface<TreeType, ScaleOp>(weight) { }
        NodeHeightPriorAlphaScaler(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<TreeType, ScaleOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "NodeHeightPriorAlphaScaler";
        }

        std::string target_parameter() const {
            return "alpha of node height beta distribution";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::node_height_prior_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::hyper;
        }

        bool is_operable(TreeType * tree) const {
            if (tree->alpha_of_node_height_beta_prior_is_fixed()) {
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            double new_alpha = tree->get_alpha_of_node_height_beta_prior();
            double ln_multiplier;
            this->update(rng, new_alpha, ln_multiplier);
            tree->set_alpha_of_node_height_beta_prior(new_alpha);
            return ln_multiplier;
        }
};


template<class TreeType>
class NodeHeightPriorAlphaMover : public GeneralTreeOperatorInterface<TreeType, WindowOp> {
    public:
        NodeHeightPriorAlphaMover() : GeneralTreeOperatorInterface<TreeType, WindowOp>() { }
        NodeHeightPriorAlphaMover(double weight) : GeneralTreeOperatorInterface<TreeType, WindowOp>(weight) { }
        NodeHeightPriorAlphaMover(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<TreeType, WindowOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "NodeHeightPriorAlphaMover";
        }

        std::string target_parameter() const {
            return "alpha of node height beta distribution";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::node_height_prior_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::hyper;
        }

        bool is_operable(TreeType * tree) const {
            if (tree->alpha_of_node_height_beta_prior_is_fixed()) {
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            double new_alpha = tree->get_alpha_of_node_height_beta_prior();
            double ln_hastings;
            this->update(rng, new_alpha, ln_hastings);
            if (new_alpha <= 0.0) {
                return -std::numeric_limits<double>::infinity();
            }
            tree->set_alpha_of_node_height_beta_prior(new_alpha);
            return ln_hastings;
        }
};


template<class TreeType>
class NodeHeightPriorBetaScaler : public GeneralTreeOperatorInterface<TreeType, ScaleOp> {
    public:
        NodeHeightPriorBetaScaler() : GeneralTreeOperatorInterface<TreeType, ScaleOp>() { }
        NodeHeightPriorBetaScaler(double weight) : GeneralTreeOperatorInterface<TreeType, ScaleOp>(weight) { }
        NodeHeightPriorBetaScaler(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<TreeType, ScaleOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "NodeHeightPriorBetaScaler";
        }

        std::string target_parameter() const {
            return "beta of node height beta distribution";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::node_height_prior_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::hyper;
        }

        bool is_operable(TreeType * tree) const {
            if (tree->beta_of_node_height_beta_prior_is_fixed()) {
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            double new_beta = tree->get_beta_of_node_height_beta_prior();
            double ln_multiplier;
            this->update(rng, new_beta, ln_multiplier);
            tree->set_beta_of_node_height_beta_prior(new_beta);
            return ln_multiplier;
        }
};


template<class TreeType>
class NodeHeightPriorBetaMover : public GeneralTreeOperatorInterface<TreeType, WindowOp> {
    public:
        NodeHeightPriorBetaMover() : GeneralTreeOperatorInterface<TreeType, WindowOp>() { }
        NodeHeightPriorBetaMover(double weight) : GeneralTreeOperatorInterface<TreeType, WindowOp>(weight) { }
        NodeHeightPriorBetaMover(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<TreeType, WindowOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "NodeHeightPriorBetaMover";
        }

        std::string target_parameter() const {
            return "beta of node height beta distribution";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::node_height_prior_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::hyper;
        }

        bool is_operable(TreeType * tree) const {
            if (tree->beta_of_node_height_beta_prior_is_fixed()) {
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            double new_beta = tree->get_beta_of_node_height_beta_prior();
            double ln_hastings;
            this->update(rng, new_beta, ln_hastings);
            if (new_beta <= 0.0) {
                return -std::numeric_limits<double>::infinity();
            }
            tree->set_beta_of_node_height_beta_prior(new_beta);
            return ln_hastings;
        }
};


template<class TreeType>
class TreeScaler : public GeneralTreeOperatorInterface<TreeType, ScaleOp> {

    public:
        TreeScaler() : GeneralTreeOperatorInterface<TreeType, ScaleOp>() { }
        TreeScaler(double weight) : GeneralTreeOperatorInterface<TreeType, ScaleOp>(weight) { }
        TreeScaler(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<TreeType, ScaleOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "TreeScaler";
        }

        std::string target_parameter() const {
            return "node heights";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::global_height_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::global;
        }

        bool is_operable(TreeType * tree) const {
            if (tree->root_height_is_fixed()) {
                return false;
            }
            return true;
        }
        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            unsigned int num_heights = tree->get_number_of_node_heights();
            double multiplier = this->op_.get_move_amount(rng);
            tree->scale_tree(multiplier);
            return std::log(multiplier) * num_heights;
        }
};


template<class TreeType>
class NodeHeightScaler : public GeneralTreeOperatorInterface<TreeType, ScaleOp> {
    public:
        NodeHeightScaler() : GeneralTreeOperatorInterface<TreeType, ScaleOp>() { }
        NodeHeightScaler(double weight) : GeneralTreeOperatorInterface<TreeType, ScaleOp>(weight) { }
        NodeHeightScaler(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<TreeType, ScaleOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "NodeHeightScaler";
        }

        std::string target_parameter() const {
            return "node heights";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::node_height_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::node_height;
        }

        bool is_operable(TreeType * tree) const {
            unsigned int num_heights = tree->get_number_of_node_heights();
            if (num_heights < 2) {
                // No non-root heights to operate on
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            unsigned int num_heights = tree->get_number_of_node_heights();
            unsigned int max_height_index = num_heights - 2;
            unsigned int height_index = rng.uniform_int(0,
                    max_height_index);
            double new_height = tree->get_height(height_index);
            double ln_multiplier;
            this->update(rng, new_height, ln_multiplier);
            if ((new_height < tree->get_height_of_oldest_child(height_index)) ||
                    (new_height > tree->get_height_of_youngest_parent(height_index))) {
                return -std::numeric_limits<double>::infinity();
            }
            tree->set_height(height_index, new_height);
            return ln_multiplier;
        }
};


template<class TreeType>
class RootHeightScaler : public GeneralTreeOperatorInterface<TreeType, ScaleOp> {
    public:
        RootHeightScaler() : GeneralTreeOperatorInterface<TreeType, ScaleOp>() { }
        RootHeightScaler(double weight) : GeneralTreeOperatorInterface<TreeType, ScaleOp>(weight) { }
        RootHeightScaler(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<TreeType, ScaleOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "RootHeightScaler";
        }

        std::string target_parameter() const {
            return "root height";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::root_height_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::root_height;
        }

        bool is_operable(TreeType * tree) const {
            if (tree->root_height_is_fixed()) {
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            unsigned int num_heights = tree->get_number_of_node_heights();
            unsigned int height_index = num_heights - 1;
            double height = tree->get_height(height_index);
            double ln_multiplier;
            this->update(rng, height, ln_multiplier);
            if (height < tree->get_height_of_oldest_child(height_index)) {
                return -std::numeric_limits<double>::infinity();
            }
            tree->set_height(height_index, height);
            return ln_multiplier;
        }
};


template<class TreeType>
class GlobalNodeHeightDirichletOperator : public GeneralTreeOperatorInterface<TreeType, DirichletOp> {

    public:
        GlobalNodeHeightDirichletOperator() : GeneralTreeOperatorInterface<TreeType, DirichletOp>() { }
        GlobalNodeHeightDirichletOperator(double weight) : GeneralTreeOperatorInterface<TreeType, DirichletOp>(weight) { }
        GlobalNodeHeightDirichletOperator(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<TreeType, DirichletOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "GlobalNodeHeightDirichletOperator";
        }

        std::string target_parameter() const {
            return "node heights";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::global_height_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::global;
        }

        bool is_operable(TreeType * tree) const {
            unsigned int num_heights = tree->get_number_of_node_heights();
            if (num_heights < 2) {
                // No non-root heights to operate on
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            unsigned int nheights = tree->get_number_of_node_heights();
            double root_height = tree->get_root_height();
            std::vector<double> rel_height_gaps(nheights);
            double last_rel_height = 0.0;
            for (unsigned int i = 0; i < nheights; ++i) {
                double rel_ht = tree->get_height(i) / root_height;
                rel_height_gaps.at(i) = rel_ht - last_rel_height;
                last_rel_height = rel_ht;
            }

            double ln_hastings;
            this->op_.update_vector(rng, rel_height_gaps, ln_hastings);

            last_rel_height = 0.0;
            // Update all heights except the root
            for (unsigned int i = 0; i < nheights - 1; ++i) {
                double rel_ht = rel_height_gaps.at(i) + last_rel_height;
                double abs_value = rel_ht * root_height;
                tree->node_heights_.at(i)->set_value(abs_value);
                last_rel_height = rel_ht;
            }
            tree->sort_node_heights();
            tree->make_dirty();
            return ln_hastings;
        }
};


template<class TreeType>
class NodeHeightDirichletOperator : public GeneralTreeOperatorInterface<TreeType, DirichletOp> {

    public:
        NodeHeightDirichletOperator() : GeneralTreeOperatorInterface<TreeType, DirichletOp>() { }
        NodeHeightDirichletOperator(double weight) : GeneralTreeOperatorInterface<TreeType, DirichletOp>(weight) { }
        NodeHeightDirichletOperator(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<TreeType, DirichletOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "NodeHeightDirichletOperator";
        }

        std::string target_parameter() const {
            return "node heights";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::node_height_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::node_height;
        }

        bool is_operable(TreeType * tree) const {
            unsigned int num_heights = tree->get_number_of_node_heights();
            if (num_heights < 2) {
                // No non-root heights to operate on
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            unsigned int num_heights = tree->get_number_of_node_heights();
            unsigned int max_height_index = num_heights - 2;
            unsigned int height_index = rng.uniform_positive_int(
                    max_height_index);
            double height = tree->get_height(height_index);
            double upper_height = tree->get_height_of_youngest_parent(height_index);
            double lower_height = tree->get_height_of_oldest_child(height_index);
            double shifted_upper_height = upper_height - lower_height;
            double shifted_height = height - lower_height;
            double rel_shifted_height = shifted_height / shifted_upper_height;

            double ln_hastings;
            this->op_.update(rng, rel_shifted_height, ln_hastings);

            double new_height = (rel_shifted_height * shifted_upper_height) + lower_height;

            tree->set_height(height_index, new_height);
            return ln_hastings;
        }
};


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// This operator is not working correctly
///////////////////////////////////////////////////////////////////////////////
// template<class TreeType>
// class NodeHeightBetaOperator : public GeneralTreeOperatorInterface<TreeType, BetaOp> {
// 
//     public:
//         NodeHeightBetaOperator() : GeneralTreeOperatorInterface<TreeType, BetaOp>() { }
//         NodeHeightBetaOperator(double weight) : GeneralTreeOperatorInterface<TreeType, BetaOp>(weight) { }
// 
//         std::string get_name() const {
//             return "NodeHeightBetaOperator";
//         }
// 
//         std::string target_parameter() const {
//             return "node heights";
//         }
// 
//         BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
//             return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::node_height_operator;
//         }
//         BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
//             return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::node_height;
//         }
// 
//         bool is_operable(TreeType * tree) const {
//             unsigned int num_heights = tree->get_number_of_node_heights();
//             if (num_heights < 2) {
//                 // No non-root heights to operate on
//                 return false;
//             }
//             return true;
//         }
// 
//         /**
//          * @brief   Propose a new state.
//          *
//          * @return  Log of Hastings Ratio.
//          */
//         double propose(RandomNumberGenerator& rng,
//                 TreeType * tree,
//                 unsigned int nthreads = 1) {
//             if (! this->is_operable(tree)) {
//                 this->ignore_proposal_attempt_ = true;
//                 return -std::numeric_limits<double>::infinity();
//             }
//             unsigned int num_heights = tree->get_number_of_node_heights();
//             unsigned int max_height_index = num_heights - 2;
//             unsigned int height_index = rng.uniform_positive_int(
//                     max_height_index);
//             double height = tree->get_height(height_index);
//             double upper_height = tree->get_height_of_youngest_parent(height_index);
//             double lower_height = tree->get_height_of_oldest_child(height_index);
//             double shifted_upper_height = upper_height - lower_height;
//             double shifted_height = height - lower_height;
//             double rel_shifted_height = shifted_height / shifted_upper_height;
// 
//             double ln_hastings;
//             this->op_.update(rng, rel_shifted_height, ln_hastings);
// 
//             double new_height = (rel_shifted_height * shifted_upper_height) + lower_height;
// 
//             tree->set_height(height_index, new_height);
//             return ln_hastings;
//         }
// };
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


template<class TreeType>
class NodeHeightMover : public GeneralTreeOperatorInterface<TreeType, WindowOp> {
    public:
        NodeHeightMover() : GeneralTreeOperatorInterface<TreeType, WindowOp>() { }
        NodeHeightMover(double weight) : GeneralTreeOperatorInterface<TreeType, WindowOp>(weight) { }
        NodeHeightMover(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<TreeType, WindowOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "NodeHeightMover";
        }

        std::string target_parameter() const {
            return "node heights";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::node_height_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::node_height;
        }

        bool is_operable(TreeType * tree) const {
            unsigned int num_heights = tree->get_number_of_node_heights();
            if (num_heights < 2) {
                // No non-root heights to operate on
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            unsigned int num_heights = tree->get_number_of_node_heights();
            unsigned int max_height_index = num_heights - 2;
            unsigned int height_index = rng.uniform_int(0,
                    max_height_index);
            double height = tree->get_height(height_index);
            double ln_hastings;
            this->update(rng, height, ln_hastings);
            if (height > tree->get_height_of_youngest_parent(height_index)) {
                return -std::numeric_limits<double>::infinity();
            }
            if (height < tree->get_height_of_oldest_child(height_index)) {
                return -std::numeric_limits<double>::infinity();
            }
            tree->set_height(height_index, height);
            return ln_hastings;
        }
};


template<class TreeType>
class NodeHeightSlideBumpScaler : public GeneralTreeOperatorInterface<TreeType, ScaleOp> {
    protected:
        virtual bool call_tree_method_(
                TreeType * tree,
                RandomNumberGenerator& rng,
                unsigned int height_index,
                double height) {
            return tree->slide_bump_height(rng,
                    height_index,
                    height);
        }

        bool operate_on_root_ = false;

    public:
        NodeHeightSlideBumpScaler() : GeneralTreeOperatorInterface<TreeType, ScaleOp>() { }
        NodeHeightSlideBumpScaler(double weight) : GeneralTreeOperatorInterface<TreeType, ScaleOp>(weight) { }
        NodeHeightSlideBumpScaler(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<TreeType, ScaleOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "NodeHeightSlideBumpScaler";
        }

        std::string target_parameter() const {
            return "node heights";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::node_height_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::node_height;
        }

        virtual void set_operate_on_root(bool operate_on_root) {
            this->operate_on_root_ = operate_on_root;
        }

        bool is_operable(TreeType * tree) const {
            unsigned int num_heights = tree->get_number_of_node_heights();
            if ((! this->operate_on_root_) && (num_heights < 2)) {
                // No non-root heights to operate on
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            unsigned int num_heights = tree->get_number_of_node_heights();
            unsigned int max_height_index = num_heights - 2;
            if (this->operate_on_root_) {
                max_height_index = num_heights - 1;
            }
            unsigned int height_index = rng.uniform_int(0,
                    max_height_index);
            double new_height = tree->get_height(height_index);
            double ln_multiplier;
            this->update(rng, new_height, ln_multiplier);
            if (new_height < 0) {
                return -std::numeric_limits<double>::infinity();
            }
            if (new_height > tree->get_root_height()) {
                if (tree->root_height_is_fixed() || (! this->operate_on_root_)) {
                    return -std::numeric_limits<double>::infinity();
                }
            }
            bool move_happened = this->call_tree_method_(
                    tree,
                    rng,
                    height_index,
                    new_height);
            if (! move_happened) {
                return -std::numeric_limits<double>::infinity();
            }
            return ln_multiplier;
        }
};


template<class TreeType>
class NodeHeightSlideBumpPermuteScaler : public NodeHeightSlideBumpScaler<TreeType> {
    protected:
        bool call_tree_method_(
                TreeType * tree,
                RandomNumberGenerator& rng,
                unsigned int height_index,
                double height) {
            return tree->slide_bump_permute_height(rng,
                    height_index,
                    height);
        }

    public:
        NodeHeightSlideBumpPermuteScaler() : NodeHeightSlideBumpScaler<TreeType>() { }
        NodeHeightSlideBumpPermuteScaler(double weight) : NodeHeightSlideBumpScaler<TreeType>(weight) { }
        NodeHeightSlideBumpPermuteScaler(double weight, double tuning_parameter)
            : NodeHeightSlideBumpScaler<TreeType>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "NodeHeightSlideBumpPermuteScaler";
        }
};


template<class TreeType>
class NodeHeightSlideBumpSwapScaler : public NodeHeightSlideBumpScaler<TreeType> {
    protected:
        bool call_tree_method_(
                TreeType * tree,
                RandomNumberGenerator& rng,
                unsigned int height_index,
                double height) {
            return tree->slide_bump_swap_height(rng,
                    height_index,
                    height);
        }

    public:
        NodeHeightSlideBumpSwapScaler() : NodeHeightSlideBumpScaler<TreeType>() { }
        NodeHeightSlideBumpSwapScaler(double weight) : NodeHeightSlideBumpScaler<TreeType>(weight) { }
        NodeHeightSlideBumpSwapScaler(double weight, double tuning_parameter)
            : NodeHeightSlideBumpScaler<TreeType>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "NodeHeightSlideBumpSwapScaler";
        }
};


template<class TreeType>
class NodeHeightSlideBumpSwapAllScaler : public NodeHeightSlideBumpScaler<TreeType> {
    protected:
        bool call_tree_method_(
                TreeType * tree,
                RandomNumberGenerator& rng,
                unsigned int height_index,
                double height) {
            return tree->slide_bump_swap_all_height(rng,
                    height_index,
                    height);
        }

    public:
        NodeHeightSlideBumpSwapAllScaler() : NodeHeightSlideBumpScaler<TreeType>() { }
        NodeHeightSlideBumpSwapAllScaler(double weight) : NodeHeightSlideBumpScaler<TreeType>(weight) { }
        NodeHeightSlideBumpSwapAllScaler(double weight, double tuning_parameter)
            : NodeHeightSlideBumpScaler<TreeType>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "NodeHeightSlideBumpSwapAllScaler";
        }
};


template<class TreeType>
class NodeHeightSlideBumpMover : public GeneralTreeOperatorInterface<TreeType, WindowOp> {
    protected:
        virtual bool call_tree_method_(
                TreeType * tree,
                RandomNumberGenerator& rng,
                unsigned int height_index,
                double height) {
            return tree->slide_bump_height(rng,
                    height_index,
                    height);
        }

        bool operate_on_root_ = false;

    public:
        NodeHeightSlideBumpMover() : GeneralTreeOperatorInterface<TreeType, WindowOp>() { }
        NodeHeightSlideBumpMover(double weight) : GeneralTreeOperatorInterface<TreeType, WindowOp>(weight) { }
        NodeHeightSlideBumpMover(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<TreeType, WindowOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "NodeHeightSlideBumpMover";
        }

        std::string target_parameter() const {
            return "node heights";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::node_height_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::node_height;
        }

        virtual void set_operate_on_root(bool operate_on_root) {
            this->operate_on_root_ = operate_on_root;
        }

        bool is_operable(TreeType * tree) const {
            unsigned int num_heights = tree->get_number_of_node_heights();
            if ((! this->operate_on_root_) && (num_heights < 2)) {
                // No non-root heights to operate on
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            unsigned int num_heights = tree->get_number_of_node_heights();
            unsigned int max_height_index = num_heights - 2;
            if (this->operate_on_root_) {
                max_height_index = num_heights - 1;
            }
            unsigned int height_index = rng.uniform_int(0,
                    max_height_index);
            double new_height = tree->get_height(height_index);
            double ln_hastings;
            this->update(rng, new_height, ln_hastings);
            if (new_height < 0) {
                return -std::numeric_limits<double>::infinity();
            }
            if (new_height > tree->get_root_height()) {
                if (tree->root_height_is_fixed() || (! this->operate_on_root_)) {
                    return -std::numeric_limits<double>::infinity();
                }
            }
            bool move_happened = this->call_tree_method_(
                    tree,
                    rng,
                    height_index,
                    new_height);
            if (! move_happened) {
                return -std::numeric_limits<double>::infinity();
            }
            return ln_hastings;
        }
};


template<class TreeType>
class NodeHeightSlideBumpPermuteMover : public NodeHeightSlideBumpMover<TreeType> {
    protected:
        bool call_tree_method_(
                TreeType * tree,
                RandomNumberGenerator& rng,
                unsigned int height_index,
                double height) {
            return tree->slide_bump_permute_height(rng,
                    height_index,
                    height);
        }

    public:
        NodeHeightSlideBumpPermuteMover() : NodeHeightSlideBumpMover<TreeType>() { }
        NodeHeightSlideBumpPermuteMover(double weight) : NodeHeightSlideBumpMover<TreeType>(weight) { }
        NodeHeightSlideBumpPermuteMover(double weight, double tuning_parameter)
            : NodeHeightSlideBumpMover<TreeType>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "NodeHeightSlideBumpPermuteMover";
        }
};


template<class TreeType>
class NodeHeightSlideBumpSwapMover : public NodeHeightSlideBumpMover<TreeType> {
    protected:
        bool call_tree_method_(
                TreeType * tree,
                RandomNumberGenerator& rng,
                unsigned int height_index,
                double height) {
            return tree->slide_bump_swap_height(rng,
                    height_index,
                    height);
        }

    public:
        NodeHeightSlideBumpSwapMover() : NodeHeightSlideBumpMover<TreeType>() { }
        NodeHeightSlideBumpSwapMover(double weight) : NodeHeightSlideBumpMover<TreeType>(weight) { }
        NodeHeightSlideBumpSwapMover(double weight, double tuning_parameter)
            : NodeHeightSlideBumpMover<TreeType>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "NodeHeightSlideBumpSwapMover";
        }
};


template<class TreeType>
class NodeHeightSlideBumpSwapAllMover : public NodeHeightSlideBumpMover<TreeType> {
    protected:
        bool call_tree_method_(
                TreeType * tree,
                RandomNumberGenerator& rng,
                unsigned int height_index,
                double height) {
            return tree->slide_bump_swap_all_height(rng,
                    height_index,
                    height);
        }

    public:
        NodeHeightSlideBumpSwapAllMover() : NodeHeightSlideBumpMover<TreeType>() { }
        NodeHeightSlideBumpSwapAllMover(double weight) : NodeHeightSlideBumpMover<TreeType>(weight) { }
        NodeHeightSlideBumpSwapAllMover(double weight, double tuning_parameter)
            : NodeHeightSlideBumpMover<TreeType>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "NodeHeightSlideBumpSwapAllMover";
        }
};


template<class TreeType>
class NeighborHeightNodePermute : public GeneralTreeOperatorInterface<TreeType, Op> {

    public:
        NeighborHeightNodePermute() : GeneralTreeOperatorInterface<TreeType, Op>() { }
        NeighborHeightNodePermute(double weight) : GeneralTreeOperatorInterface<TreeType, Op>(weight) { }

        std::string get_name() const {
            return "NeighborHeightNodePermute";
        }

        std::string target_parameter() const {
            return "topology";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::topology_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::topology;
        }

        bool is_operable(TreeType * tree) const {
            unsigned int num_node_heights = tree->get_number_of_node_heights();
            if (num_node_heights == 1) {
                // In comb state, so nothing to do; force rejection
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            unsigned int num_node_heights = tree->get_number_of_node_heights();
            unsigned int height_index = rng.uniform_int(0,
                    num_node_heights - 2);
            unsigned int parent_index = tree->get_index_of_youngest_parent(height_index);
            tree->collision_node_permute(rng,
                    parent_index,
                    height_index);
            return 0.0;
        }
};


template<class TreeType>
class NeighborHeightNodeSwap : public GeneralTreeOperatorInterface<TreeType, Op> {

    public:
        NeighborHeightNodeSwap() : GeneralTreeOperatorInterface<TreeType, Op>() { }
        NeighborHeightNodeSwap(double weight) : GeneralTreeOperatorInterface<TreeType, Op>(weight) { }

        std::string get_name() const {
            return "NeighborHeightNodeSwap";
        }

        std::string target_parameter() const {
            return "topology";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::topology_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::topology;
        }

        bool is_operable(TreeType * tree) const {
            unsigned int num_node_heights = tree->get_number_of_node_heights();
            if (num_node_heights == 1) {
                // In comb state, so nothing to do; force rejection
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            unsigned int num_node_heights = tree->get_number_of_node_heights();
            unsigned int height_index = rng.uniform_int(0,
                    num_node_heights - 2);
            unsigned int parent_index = tree->get_index_of_youngest_parent(height_index);
            tree->collision_node_swap(rng,
                    parent_index,
                    height_index);
            return 0.0;
        }
};

template<class TreeType>
class NeighborHeightNodeSwapAll : public NeighborHeightNodeSwap<TreeType> {

    public:
        NeighborHeightNodeSwapAll() : NeighborHeightNodeSwap<TreeType>() { }
        NeighborHeightNodeSwapAll(double weight) : NeighborHeightNodeSwap<TreeType>(weight) { }

        std::string get_name() const {
            return "NeighborHeightNodeSwapAll";
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            unsigned int num_node_heights = tree->get_number_of_node_heights();
            unsigned int height_index = rng.uniform_int(0,
                    num_node_heights - 2);
            unsigned int parent_index = tree->get_index_of_youngest_parent(height_index);
            tree->collision_node_swap_all(rng,
                    parent_index,
                    height_index);
            return 0.0;
        }
};


template<class TreeType>
class SplitLumpNodesRevJumpSampler : public GeneralTreeOperatorInterface<TreeType, Op> {
    protected:
        std::unordered_map<unsigned int, long double> stirling2_numbers;
        std::unordered_map<unsigned int, long double> bell_numbers;
        double beta_a_ = 1.0;
        double beta_b_ = 1.0;

    public:
        SplitLumpNodesRevJumpSampler() : GeneralTreeOperatorInterface<TreeType, Op>() { }
        SplitLumpNodesRevJumpSampler(double weight) : GeneralTreeOperatorInterface<TreeType, Op>(weight) { }

        std::string get_name() const {
            return "SplitLumpNodesRevJumpSampler";
        }

        std::string target_parameter() const {
            return "topology";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::topology_model_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::topology;
        }

        double get_stirling2(unsigned int n) {
            if (this->stirling2_numbers.count(n) < 1) {
                this->stirling2_numbers[n] = stirling2_float(n, 2);
            }
            return this->stirling2_numbers[n];
        }

        long double get_bell_number(unsigned int n) {
            if (this->bell_numbers.count(n) < 1) {
                this->bell_numbers[n] = bell_float(n);
            }
            return this->bell_numbers[n];
        }

        bool is_operable(TreeType * tree) const {
            return true;
        }

        void operate_plus(RandomNumberGenerator& rng,
                TreeType * tree,
                std::vector< std::shared_ptr< GeneralTreeOperatorTemplate< TreeType > > > other_operators,
                unsigned int nthreads = 1,
                unsigned int number_of_moves = 1,
                unsigned int other_op_number_of_moves = 1) {
            for (unsigned int i = 0; i < number_of_moves; ++i) {
                this->perform_move(rng, tree, nthreads);
                for (auto other_op : other_operators) {
                    other_op->operate(rng, tree, nthreads, other_op_number_of_moves);
                }
            }
        }

        double propose_split(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            const unsigned int num_heights = tree->get_number_of_node_heights();
            const bool in_comb_state_before = (num_heights == 1);
            std::vector<unsigned int> splittable_height_indices = tree->get_indices_of_splittable_heights();
            const unsigned int number_of_splittable_heights = splittable_height_indices.size();
            ECOEVOLITY_ASSERT(number_of_splittable_heights > 0);
            unsigned int splittable_vector_idx = rng.uniform_int(0, number_of_splittable_heights - 1);
            unsigned int split_height_idx = splittable_height_indices.at(splittable_vector_idx);
            double current_height = tree->get_height(split_height_idx);
            double height_lower_bound = -1.0;
            double proposed_height = -1.0;
            unsigned int number_of_mapped_nodes = 0;
            unsigned int number_of_nodes_in_split_subset;
            std::vector<unsigned int> moving_polytomy_sizes;
            bool mapped_nodes_include_polytomy;
            double ln_prob_of_drawing_new_node_states = tree->split_node_height_down(
                    rng,
                    this->beta_a_,
                    this->beta_b_,
                    split_height_idx,
                    proposed_height,
                    height_lower_bound,
                    number_of_mapped_nodes,
                    number_of_nodes_in_split_subset,
                    moving_polytomy_sizes,
                    mapped_nodes_include_polytomy,
                    false);
            ECOEVOLITY_ASSERT(number_of_mapped_nodes > 0);
            double height_window = current_height - height_lower_bound;
            ECOEVOLITY_ASSERT(height_window > 0.0);
            double beta_val = proposed_height - height_lower_bound;
            double ln_density_of_proposed_height = BetaDistribution::get_scaled_ln_pdf(
                    beta_val,
                    this->beta_a_,
                    this->beta_b_,
                    height_window);

            // ================================================================
            // IN CASE OF A SINGLETON POLYTOMY
            // ================================================================
            // The probability of forward split move (just proposed) is the
            // product of the probabilites of:
            //   1) choosing the splittable height to split
            //          = 1 / number of splittable heights
            //   2a) drawing the new height uniformly between the height we
            //      are splitting and the next younger height (or zero).
            //          = 1 / (height - younger neighbor height)
            //          = 1 / d
            //      NOTE: This can also be thought of as the Jacobian term
            //      for the reversible jump move.
            //   2b) more generally this can be the prob density of any
            //      proposal distribution.
            //          = pdens
            //      NOTE: This can also be thought of as the Jacobian term
            //      for the reversible jump move.
            //   3) randomly partitioning the children of the one polytomy;
            //          = 1 / (Bell(nchildren) - 2)
            //      The minus 2 is because we reject 2 possible
            //      partitions: (1) where every child is in its own subset,
            //      and (2) when all children end up together in one subst.
            //      The second partition would result in the polytomy node
            //      simply moving down to the new height (not getting split
            //      up). This is normally OK, when other nodes are
            //      remaining at the old height (this is guaranteed,
            //      because we always split nodes into two non-empty
            //      subsets), but in the singleton case that doesn't
            //      happen.
            //   4) drawing the values of parameters of any new nodes that were
            //      created (e.g., the effective population size of a
            //      PopulationNode). This is what is returned by
            //      'split_node_height_down'
            //          = z
            //
            // So the prob of the forward split move is
            // p(split move) = z * pdens * 1 / (number of splittable heights *
            //                       (Bell(nchildren) - 2))
            //
            // The probability of the reverse move is simply the
            // probability of randomly selecting the proposed (split)
            // height from among all node heights except the root.
            // p(reverse merge) = 1 / (number of heights before proposal + 1 - 1)
            //                  = 1 / (number of heights before the proposal)
            //
            // So, the Hasting ratio for the proposed split move is:
            // p(reverse merge) / p(proposed split) = 
            //       (number of splittable heights * (Bell(n) - 2))      1
            //     -------------------------------------------------- * ---
            //      pdens * (number of heights before the proposal)      z
            //
            // ================================================================
            // IN CASE OF SHARED BIFURCATING NODES (NO POLYTOMIES)
            // ================================================================
            // The probability of forward split move (just proposed) is the
            // product of the probabilites of:
            //   1) choosing the splittable height to split
            //          = 1 / number of splittable heights
            //   2) randomly splitting the subset out of the 'n' nodes
            //      mapped to the height
            //          = 1 / (2 * stirling2(n, 2))
            //          = 1 / (2^n - 2)
            //   3) drawing the new height uniformly between the height we
            //      are splitting and the next younger height (or zero).
            //          = 1 / (height - younger neighbor height)
            //          = 1 / d
            //          OR more generally
            //          = pdens
            //      NOTE: This can also be thought of as the Jacobian term
            //      for the reversible jump move.
            //
            // So the prob of the forward split move is
            // p(split move) =  pdens * 1 / (number of splittable heights *
            //                       (2^n - 2))
            //               =  pdens * 1 / (number of splittable heights * 
            //                       2 * stirling2(n, 2))
            //
            // The probability of the reverse move is simply the
            // probability of randomly selecting the proposed (split)
            // height from among all node heights except the root.
            // p(reverse merge) = 1 / (number of heights before proposal + 1 - 1)
            //                  = 1 / (number of heights before the proposal)
            //
            // So, the Hasting ratio for the proposed split move is:
            // p(reverse merge) / p(proposed split) = 
            //      (number of splittable heights * 2 * stirling2(n, 2))
            //     --------------------------------------------------------
            //         (number of heights before the proposal) * pdens
            //
            // ================================================================
            // IN CASE WHEN SUBSET OF SHARED NODES INCLUDE POLYTOMIES
            // ================================================================
            // The probability of forward split move (just proposed) is the
            // product of the probabilites of:
            //   1) choosing the splittable height to split
            //          = 1 / number of splittable heights
            //   2) randomly splitting the subset out of the 'n' nodes
            //      mapped to the height. Because there are polytomies,
            //      this has to be done differently. We need to allow
            //      all nodes ending up in the subset to be split down.
            //      We can do this by splitting the 'n' nodes into 2 subsets
            //      and ALLOWING an empty subset. There number of ways to do
            //      this is:
            //          stirling2(n, 2) + 1
            //      We randomly pick one of the 2 sets so there are now
            //          2 * (stirling2(n, 2) + 1)
            //      ways to do this. But, we won't allow the empty set
            //      to be sampled to move, so there are now
            //          (2 * (stirling2(n, 2) + 1)) - 1
            //          = (2 * stirling2(n, 2)) + 2 - 1
            //          = (2 * stirling2(n, 2)) + 1
            //      ways to choose our subset of nodes to be split/moved.
            //      We avoid the empty set via rejection during the move, so
            //      all the probability of all other ways remains uniform:
            //          = 1 / ((2 * stirling2(n, 2)) + 1)
            //   3) drawing the new height uniformly between the height we
            //      are splitting and the next younger height (or zero).
            //          = 1 / (height - younger neighbor height)
            //          = 1 / d
            //          OR more generally
            //          = pdens
            //      NOTE: This can also be thought of as the Jacobian term
            //      for the reversible jump move.
            //   4) randomly partitioning the children of each polytomy;
            //      each subset with 2 or more children is split off as a
            //      clade and mapped to the new, younger height. The
            //      number of ways to do this is:
            //          = Bell(nchildren) - 1
            //      where Bell(n) is the Bell number, and 1 is subtracted,
            //      because we reject the partition where every child is in
            //      its own subset.  We reject this, because then the node
            //      would not be split or moved to the new height. We have
            //      to take the product of this over all polytomy nodes
            //      that end up in the ``split'' subset of nodes. So, the
            //      number of ways across all polytomy nodes is:
            //          = \prod (Bell(nchildren) - 1)
            //      We sample these uniformly, so the prob is
            //          = 1/ \prod (Bell(nchildren) - 1)
            //   5) drawing the values of parameters of any new nodes that were
            //      created (e.g., the effective population size of a
            //      PopulationNode). This is what is returned by
            //      'split_node_height_down'
            //          = z
            //
            // So the prob of the forward split move is
            // p(split move) = z * pdens * 1 / (number of splittable heights *
            //                       ((2 * stirling2(n, 2)) + 1) *
            //                       \prod (Bell(nchildren) - 1))
            //
            // The probability of the reverse move is simply the
            // probability of randomly selecting the proposed (split)
            // height from among all node heights except the root.
            // p(reverse merge) = 1 / (number of heights before proposal + 1 - 1)
            //                  = 1 / (number of heights before the proposal)
            //
            // So, the Hasting ratio for the proposed split move is:
            // p(reverse merge) / p(proposed split) = (1/z) *
            //     (number of splittable heights * ((2 * stirling2(n, 2)) + 1) * \prod (Bell(nchild) - 1))
            //     ------------------------------------------------------------------------------------------
            //          (number of heights before the proposal) * pdens
            //
            // ================================================================
            // IN CASE WHEN ALL SHARED NODES IN MOVE SET AND INCLUDES
            // POLYTOMIES
            // ================================================================
            // The probability of forward split move (just proposed) is the
            // product of the probabilites of:
            //   1) choosing the splittable height to split
            //          = 1 / number of splittable heights
            //   2) randomly splitting the subset out of the 'n' nodes
            //      mapped to the height. Because there are polytomies,
            //      this has to be done differently. We need to allow
            //      all nodes ending up in the subset to be split down.
            //      We can do this by splitting the 'n' nodes into 2 subsets
            //      and ALLOWING an empty subset. There number of ways to do
            //      this is:
            //          stirling2(n, 2) + 1
            //      We randomly pick one of the 2 sets so there are now
            //          2 * (stirling2(n, 2) + 1)
            //      ways to do this. But, we won't allow the empty set
            //      to be sampled to move, so there are now
            //          (2 * (stirling2(n, 2) + 1)) - 1
            //          = (2 * stirling2(n, 2)) + 2 - 1
            //          = (2 * stirling2(n, 2)) + 1
            //      ways to choose our subset of nodes to be split/moved.
            //      We avoid the empty set via rejection during the move, so
            //      all the probability of all other ways remains uniform:
            //          = 1 / ((2 * stirling2(n, 2)) + 1)
            //   3) drawing the new height uniformly between the height we
            //      are splitting and the next younger height (or zero).
            //          = 1 / (height - younger neighbor height)
            //          = 1 / d
            //          OR more generally
            //          = pdens
            //      NOTE: This can also be thought of as the Jacobian term
            //      for the reversible jump move.
            //   4) randomly partitioning the children of each polytomy;
            //      each subset with 2 or more children is split off as a
            //      clade and mapped to the new, younger height. The
            //      number of ways to do this is:
            //          = Bell(nchildren) - 1
            //      where Bell(n) is the Bell number, and 1 is subtracted,
            //      because we reject the partition where every child is in
            //      its own subset.  We reject this, because then the node
            //      would not be split or moved to the new height. We have
            //      to take the product of this over all polytomy nodes
            //      that end up in the ``split'' subset of nodes. So, the
            //      number of ways across all polytomy nodes is:
            //          = \prod (Bell(nchildren) - 1)
            //      NOTE: we reject the situation (during the move) where all
            //      nodes get moved down in full (i.e., no nodes remain at old
            //      node height), because this doesn't result in a model jump
            //      and cannot be reversed by a merge move. In this special
            //      case, the number of ways to split up all the nodes is:
            //          = [\prod (Bell(nchildren) - 1)] - 1
            //   5) drawing the values of parameters of any new nodes that were
            //      created (e.g., the effective population size of a
            //      PopulationNode). This is what is returned by
            //      'split_node_height_down'
            //          = z
            //
            // So the prob of the forward split move is
            // p(split move) = z * pdens * 1 / (
            //                       number of splittable heights *
            //                       ((2 * stirling2(n, 2)) + 1) *
            //                       ((\prod (Bell(nchildren) - 1)) - 1)
            //                       )
            //
            // The probability of the reverse move is simply the
            // probability of randomly selecting the proposed (split)
            // height from among all node heights except the root.
            // p(reverse merge) = 1 / (number of heights before proposal + 1 - 1)
            //                  = 1 / (number of heights before the proposal)
            //
            // So, the Hasting ratio for the proposed split move is:
            // p(reverse merge) / p(proposed split) = (1/z)
            //     (number of splittable heights * ((2 * stirling2(n, 2)) + 1) * ((\prod (Bell(nchild) - 1)) - 1)
            //     ------------------------------------------------------------------------------------------
            //          (number of heights before the proposal) * pdens
            //

            double ln_hastings = 0.0;
            if (number_of_mapped_nodes == 1) {
                // We have the special case of only a single polytomy
                ECOEVOLITY_ASSERT(moving_polytomy_sizes.size() == 1);
                double ln_bell_num_minus_2 = std::log(
                        this->get_bell_number(moving_polytomy_sizes.at(0)) - 2.0);
                ln_hastings =
                        std::log(number_of_splittable_heights) +
                        ln_bell_num_minus_2;
                ln_hastings -= (std::log(num_heights) +
                        ln_prob_of_drawing_new_node_states +
                        ln_density_of_proposed_height);
            }
            else if (! mapped_nodes_include_polytomy) {
                // We have multiple bifurcating nodes mapped to height
                double ln_stirling2 = std::log(this->get_stirling2(number_of_mapped_nodes));
                ln_hastings =
                        std::log(number_of_splittable_heights) +
                        std::log(2.0) +
                        ln_stirling2;
                ln_hastings -= (std::log(num_heights) +
                        ln_prob_of_drawing_new_node_states +
                        ln_density_of_proposed_height);
            }
            else {
                // We have multiple nodes mapped to height and some are
                // polytomies
                bool hit_overflow = false;
                long double long_double_max = std::numeric_limits<long double>::max();
                long double stirl2_term = this->get_stirling2(number_of_mapped_nodes);
                // Check for multiplication overflow
                if (2.0 > (long_double_max / stirl2_term)) {
                    hit_overflow = true;
                } else {
                    stirl2_term *= 2.0;
                }
                // Check for addition overflow
                if (stirl2_term > (long_double_max - 1.0)) {
                    hit_overflow = true;
                } else {
                    stirl2_term += 1.0;
                }
                double ln_stirling2_term;
                if (hit_overflow) {
                    // If we hit overflow, we will ignore the tiny plus 1 to 
                    // (2 * Stirling)
                    ln_stirling2_term = std::log(2.0) + std::log(this->get_stirling2(number_of_mapped_nodes));
                } else {
                    ln_stirling2_term = std::log(stirl2_term);
                }

                double ln_bell_num_minus_1_sum = 0.0;
                for (auto polytomy_size : moving_polytomy_sizes) {
                    ln_bell_num_minus_1_sum += std::log(this->get_bell_number(polytomy_size) - 1.0);
                }
                double ln_bell_term = ln_bell_num_minus_1_sum;

                if (number_of_mapped_nodes == number_of_nodes_in_split_subset) {
                    ECOEVOLITY_ASSERT(moving_polytomy_sizes.size() > 0);
                    // If all nodes mapped to height ended up in the move set,
                    // we need to account for the case we reject where none of
                    // the polytomies get broken up (i.e., all node simply
                    // slide down and no parameter is added to model).
                    hit_overflow = false;
                    long double bell_num_minus_1_prod = 1.0;
                    for (auto polytomy_size : moving_polytomy_sizes) {
                        // Check for multiplication overflow
                        long double bell_minus_1 = this->get_bell_number(polytomy_size) - 1.0;
                        if (bell_minus_1 > (long_double_max / bell_num_minus_1_prod)) {
                            hit_overflow = true;
                            break;
                        }
                        bell_num_minus_1_prod *= bell_minus_1;
                    }
                    // If we overflowed, we will use ln_bell_num_minus_1_sum,
                    // because the minus 1 from the product will be tiny
                    // TODO: We could probably be less conservative and use
                    // ln_bell_num_minus_1_sum when the bell_num_minus_1_prod
                    // gets to a threshold much smaller than the numeric limit.
                    if (! hit_overflow) {
                        ln_bell_term = std::log(bell_num_minus_1_prod - 1.0);
                    }
                }
                ln_hastings =
                        std::log(number_of_splittable_heights) +
                        ln_stirling2_term +
                        ln_bell_term;
                ln_hastings -= (std::log(num_heights) +
                        ln_prob_of_drawing_new_node_states +
                        ln_density_of_proposed_height);
            }

            // For the Hastings ratio, we also have to consider the
            // probability of proposing a split move versus a merge move.
            // The ratio of the p(proposing reverse merge) / p(proposing
            // this split) is 1.0, exept for two corner cases:
            //  1. If we are in the comb state (nheights = 1) before
            //     the move, and NOT in the general state (i.e.,
            //     nheights = nleaves - 1) after the move, the ratio
            //     is:
            //
            //     p(propose rev merge) / p(propose split)
            //     = 0.5 / 1.0 = 0.5
            //
            //  2. If we are in general state (nheights = nleaves - 1)
            //     AFTER the move, and NOT in the comb state (nheigths
            //     = 1) BEFORE, the ratio is:
            //     
            //     p(propose rev merge) / p(propose split)
            //     = 1.0 / 0.5 = 2.0
            const bool in_general_state_after = (tree->get_number_of_node_heights() ==
                    (tree->get_leaf_node_count() - 1));

            if (in_comb_state_before && (! in_general_state_after)) {
                ln_hastings -= std::log(2.0);
            }
            else if (in_general_state_after && (! in_comb_state_before)) {
                ln_hastings += std::log(2.0);
            }
            ECOEVOLITY_ASSERT(! std::isnan(ln_hastings));
            return ln_hastings;
        }

        double propose_merge(RandomNumberGenerator& rng,
                TreeType * tree,
                const bool in_general_state_before,
                unsigned int nthreads = 1) {
            const unsigned int num_heights = tree->get_number_of_node_heights();
            ECOEVOLITY_ASSERT(in_general_state_before ==
                    (num_heights == tree->get_leaf_node_count() - 1));
            // MERGE MOVE
            //
            // For the merge move, we simply randomly select a height that
            // isn't the root and merge with the next older height.
            // So, the probability of the forward merge move is simply
            //   = 1 / number of non-root heights before merge
            //   = 1 / nheights before merge - 1
            //   = 1 / nheights after merge
            // 
            // The probability of the split move that would reverse
            // the merge is detailed above, and varies under some
            // conditions:
            // 
            // When only one node that is a polytomy is mapped to the
            // post-merge height:
            //   p(rev split move) =  z * pdens * 1 / (post-merge num of splittable heights *
            //                            (Bell(nchildren) - 2))
            //   so, HR =
            //                  z * nheights after merge * pdens
            //   ----------------------------------------------------------------
            //   (post-merge num of splittable heights * bell(nchildren) - 2)
            //
            // When only multiple bifurcating nodes are mapped to the
            // post-merge height:
            //   p(rev split move) =  pdens * 1 / (post-merge num of splittable heights *
            //                            2 * stirling2(n, 2))
            //   so, HR =
            //                   nheights after merge * pdens
            //   ----------------------------------------------------------------
            //   (post-merge num of splittable heights * 2 * Stirling2(n, 2))
            //
            // When multiple nodes are mapped to the post-merge height, and
            // include at least one polytomy:
            //   IF THE NUMBER OF NODES MAPPED TO OLD HEIGHT (PRE-MERGE) WAS
            //   EQUAL TO THE NUMBER THAT WAS MAPPED TO THE POST-MERGE HEIGHT
            //   (i.e., in the reverse split move, all nodes are included in
            //   the move subset):
            //     p(rev split move) = z * pdens * 1 / (
            //                           post-merge num of splittable heights *
            //                           ((2 * stirling2(n, 2)) + 1) *
            //                           ((\prod (Bell(nchildren) - 1)) - 1)
            //                           )
            //     So, HR =
            //               z * nheights after merge * pdens
            //     -------------------------------------------------
            //     (post-merge num of splittable heights *
            //         ((2 * stirling2(n, 2)) + 1) *
            //         ((\prod (Bell(nchildren) - 1)) - 1))
            //   ELSE:
            //     p(rev split move) = z * pdens * 1 / (
            //                           post-merge num of splittable heights *
            //                           ((2 * stirling2(n, 2)) + 1) *
            //                           \prod (Bell(nchildren) - 1)
            //                           )
            //     So, HR =
            //              z * nheights after merge * pdens
            //     -------------------------------------------------
            //     (post-merge num of splittable heights *
            //         ((2 * stirling2(n, 2)) + 1) *
            //         \prod (Bell(nchildren) - 1)
            //         )

            // std::cout << "Merging...\n";
            const unsigned int merge_height_idx = rng.uniform_int(0, num_heights - 2);
            const double original_height = tree->get_height(merge_height_idx);
            std::vector<unsigned int> sizes_of_mapped_polytomies_after_merge;
            unsigned int number_of_resulting_merged_nodes;
            double ln_prob_of_drawing_old_node_states = tree->merge_node_height_up(merge_height_idx,
                    sizes_of_mapped_polytomies_after_merge,
                    number_of_resulting_merged_nodes);
            const double older_height = tree->get_height(merge_height_idx);
            double younger_height = 0.0;
            if (merge_height_idx > 0) {
                younger_height = tree->get_height(merge_height_idx - 1);
            }
            const double height_window = older_height - younger_height;
            ECOEVOLITY_ASSERT(height_window > 0.0);
            double rev_beta_val = original_height - younger_height;
            double ln_density_of_rev_height = BetaDistribution::get_scaled_ln_pdf(
                    rev_beta_val,
                    this->beta_a_,
                    this->beta_b_,
                    height_window);
            const unsigned int post_num_splittable_heights = tree->get_number_of_splittable_heights();
            std::vector< typename TreeType::NodePtr > post_mapped_nodes = tree->get_mapped_nodes(merge_height_idx);
            const unsigned int post_num_mapped_nodes = post_mapped_nodes.size();
            unsigned int post_num_mapped_poly_nodes = 0; 
            for (auto node : post_mapped_nodes) {
                if (node->is_polytomy()) {
                    ++post_num_mapped_poly_nodes;
                }
            }
            double ln_hastings = 0.0;
            if (post_num_mapped_nodes == 1) {
                const unsigned int num_polytomy_children = post_mapped_nodes.at(0)->get_number_of_children();
                ECOEVOLITY_ASSERT(num_polytomy_children > 2);
                ln_hastings =
                        std::log(num_heights - 1) +
                        ln_prob_of_drawing_old_node_states +
                        ln_density_of_rev_height;
                ln_hastings -=
                        (std::log(post_num_splittable_heights) +
                        std::log(this->get_bell_number(num_polytomy_children) - 2.0));
            }
            else if (post_num_mapped_poly_nodes < 1) {
                // Only shared bifurcating nodes
                double ln_stirling2_num = std::log(this->get_stirling2(post_num_mapped_nodes));
                ln_hastings =
                        std::log(num_heights - 1) +
                        ln_prob_of_drawing_old_node_states +
                        ln_density_of_rev_height;
                ln_hastings -=
                        (std::log(post_num_splittable_heights) +
                        std::log(2.0) +
                        ln_stirling2_num);
            }
            // We have shared nodes that include at least on polytomy
            else {
                bool hit_overflow = false;
                long double long_double_max = std::numeric_limits<long double>::max();
                long double stirl2_term = this->get_stirling2(post_num_mapped_nodes);
                // Check for multiplication overflow
                if (2.0 > (long_double_max / stirl2_term)) {
                    hit_overflow = true;
                } else {
                    stirl2_term *= 2.0;
                }
                // Check for addition overflow
                if (stirl2_term > (long_double_max - 1.0)) {
                    hit_overflow = true;
                } else {
                    stirl2_term += 1.0;
                }
                double ln_stirling2_term;
                if (hit_overflow) {
                    // If we hit overflow, we will ignore the tiny plus 1 to 
                    // (2 * Stirling)
                    ln_stirling2_term = std::log(2.0) + std::log(this->get_stirling2(post_num_mapped_nodes));
                } else {
                    ln_stirling2_term = std::log(stirl2_term);
                }

                double ln_bell_num_minus_1_sum = 0.0;
                for (auto poly_size : sizes_of_mapped_polytomies_after_merge) {
                    ln_bell_num_minus_1_sum += std::log(
                            this->get_bell_number(poly_size) - 1.0);
                }
                double ln_bell_term = ln_bell_num_minus_1_sum;

                if (post_num_mapped_nodes == number_of_resulting_merged_nodes) {
                    ECOEVOLITY_ASSERT(sizes_of_mapped_polytomies_after_merge.size() > 0);
                    // If all nodes mapped to height need to end up in the move
                    // set for the reverse move, we need to account for the
                    // case we reject where none of the polytomies get broken
                    // up (i.e., all node simply slide down and no parameter is
                    // added to model).
                    hit_overflow = false;
                    long double bell_num_minus_1_prod = 1.0;
                    for (auto poly_size : sizes_of_mapped_polytomies_after_merge) {
                        // Check for multiplication overflow
                        long double bell_minus_1 = this->get_bell_number(poly_size) - 1.0;
                        if (bell_minus_1 > (long_double_max / bell_num_minus_1_prod)) {
                            hit_overflow = true;
                            break;
                        }
                        bell_num_minus_1_prod *= bell_minus_1;
                    }
                    // If we overflowed, we will use ln_bell_num_minus_1_sum,
                    // because the minus 1 from the product will be tiny
                    // TODO: We could probably be less conservative and use
                    // ln_bell_num_minus_1_sum when the bell_num_minus_1_prod
                    // gets to a threshold much smaller than the numeric limit.
                    if (! hit_overflow) {
                        ln_bell_term = std::log(bell_num_minus_1_prod - 1.0);
                    }
                }
                ln_hastings =
                        std::log(num_heights - 1) +
                        ln_prob_of_drawing_old_node_states +
                        ln_density_of_rev_height;
                ln_hastings -=
                        (std::log(post_num_splittable_heights) +
                        ln_stirling2_term +
                        ln_bell_term);
            }

            // For the hastings ratio we also have to include the ratio of the
            // probability of choosing to split for the reverse move over the
            // probability of choosing to merge. Normally, this cancels out,
            // because we randomly choose (50/50) to try a merge or split move.
            // However, there are two corner cases where this ratio is not
            // equal to 1:
            //  1) We are in the GENERAL state (nheights = nleaves - 1) BEFORE
            //     the move and NOT in the COMB state (nheights = 1) AFTER
            //     the move. Then this ratio is
            //      = 0.5 / 1.0 = 0.5
            //  2) We are in the COMB state AFTER the move and NOT in the
            //     GENERAL state BEFORE the move. Then the ratio is
            //      = 1.0 / 0.5 = 2.0

            const bool in_comb_state_after = (tree->get_number_of_node_heights() == 1);
            if (in_general_state_before && (! in_comb_state_after)) {
                ln_hastings -= std::log(2.0);
            }
            else if (in_comb_state_after && (! in_general_state_before)) {
                ln_hastings += std::log(2.0);
            }
            ECOEVOLITY_ASSERT(! std::isnan(ln_hastings));
            return ln_hastings;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                TreeType * tree,
                unsigned int nthreads = 1) {
            const unsigned int num_heights = tree->get_number_of_node_heights();
            const unsigned int num_leaves = tree->get_leaf_node_count();
            ECOEVOLITY_ASSERT(num_leaves > 2);
            const bool in_comb_state_before = (num_heights == 1);
            const bool in_general_state_before = (num_heights == (num_leaves - 1));
            const bool split_event = ((! in_general_state_before) &&
                    (in_comb_state_before || (rng.uniform_real() < 0.5)));
            if (split_event) {
                return this->propose_split(rng,
                        tree,
                        nthreads);
            }
            return this->propose_merge(rng,
                    tree,
                    in_general_state_before,
                    nthreads);
        }
};


//////////////////////////////////////////////////////////////////////////////
// BasePopulationTree operators 
//////////////////////////////////////////////////////////////////////////////

class GlobalPopSizeScaler : public GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp> {

    public:
        GlobalPopSizeScaler() : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>() { }
        GlobalPopSizeScaler(double weight) : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>(weight) { }
        GlobalPopSizeScaler(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "GlobalPopSizeScaler";
        }

        std::string target_parameter() const {
            return "population sizes";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::derived_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::global;
        }

        bool is_operable(BasePopulationTree * tree) const {
            if (tree->population_sizes_are_fixed()) {
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                BasePopulationTree * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            double multiplier = this->op_.get_move_amount(rng);
            unsigned int n_parameters_scaled;
            n_parameters_scaled = tree->scale_all_population_sizes(multiplier);
            return std::log(multiplier) * n_parameters_scaled;
        }
};

class PopSizeScaler : public GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp> {

    public:
        PopSizeScaler() : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>() { }
        PopSizeScaler(double weight) : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>(weight) { }
        PopSizeScaler(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "PopSizeScaler";
        }

        std::string target_parameter() const {
            return "population sizes";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::derived_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::branch;
        }

        bool is_operable(BasePopulationTree * tree) const {
            if (tree->population_sizes_are_fixed()) {
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                BasePopulationTree * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            unsigned int random_node_index = rng.uniform_positive_int(
                    tree->pre_ordered_nodes_.size() - 1);
            double pop_size = tree->pre_ordered_nodes_.at(
                    random_node_index)->get_population_size();

            double ln_multiplier;
            this->update(rng, pop_size, ln_multiplier);

            if (pop_size <= 0.0) {
                return -std::numeric_limits<double>::infinity();
            }

            tree->pre_ordered_nodes_.at(random_node_index)->set_population_size(
                    pop_size);

            return ln_multiplier;
        }
};

class MuRateScaler : public GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp> {

    public:
        MuRateScaler() : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>() { }
        MuRateScaler(double weight) : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>(weight) { }
        MuRateScaler(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "MuRateScaler";
        }

        std::string target_parameter() const {
            return "mutation rate";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::derived_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::global;
        }

        bool is_operable(BasePopulationTree * tree) const {
            if (tree->mutation_rate_is_fixed()) {
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                BasePopulationTree * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            double mutation_rate = tree->get_mutation_rate();
            double ln_multiplier;
            this->update(rng, mutation_rate, ln_multiplier);
            if (mutation_rate <= 0.0) {
                return -std::numeric_limits<double>::infinity();
            }
            tree->set_mutation_rate(mutation_rate);
            return ln_multiplier;
        }
};

class GlobalHeightSizeMixer : public GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp> {

    public:
        GlobalHeightSizeMixer() : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>() { }
        GlobalHeightSizeMixer(double weight) : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>(weight) { }
        GlobalHeightSizeMixer(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "GlobalHeightSizeMixer";
        }

        std::string target_parameter() const {
            return "node heights and population sizes";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::global_height_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::global;
        }

        bool is_operable(BasePopulationTree * tree) const {
            if (tree->root_height_is_fixed()) {
                return false;
            }
            if (tree->population_sizes_are_fixed()) {
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         *
         * The expected height of a gene tree node given the height of the
         * species tree node (t) is:
         *
         *   t*mu + ploidy*N*mu
         *
         * where mu is the mutation rate, and N is the effective size of the
         * ancestral population.
         *
         * In this proposal we want to mainatin this relationship:
         *
         *   t'*mu + ploidy*N'*mu = t*mu + ploidy*N*mu
         *
         * where (') represents the proposed (new) values. The mutation rates
         * cancel, giving us:
         *
         *   t' + ploidy*N' = t + ploidy*N
         *
         * Re-writing this in terms of the amount of change in t and N (denoted
         * as dt and dN, respectively), we have:
         *
         *   t + dt + ploidy(N + dN) = t + ploidy*N
         *
         * multiply out N and dN on the left:
         *
         *   t + dt + ploidy*N + ploidy*dN = t + ploidy*N
         *
         * subtract ploidy*N from both sides:
         *
         *   t + dt + ploidy*dN = t
         *
         * subtract t from both sides:
         *
         *   dt + ploidy*dN = 0
         *
         *   ploidy*dN = -dt
         *
         *   dN = -dt / ploidy
         *
         * So, we will scale all node heights, and then add (dt / ploidy) to
         * all ancestral pop sizes
         * 
         */
        double propose(RandomNumberGenerator& rng,
                BasePopulationTree * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            if (tree->population_sizes_are_constrained()) {
                // When pop sizes are constrained to be equal, it makes sense
                // to scale the shared pop size, and then change all the node
                // heights accordingly. We are (probably) less likely to
                // propose a negative node height value this way then we are to
                // get a negative pop size if we scaled the node heights.
                double old_size = tree->get_root_population_size();
                double multiplier = this->op_.get_move_amount(rng);
                double new_size = old_size * multiplier;
                tree->set_root_population_size(new_size);
                double size_diff = new_size - old_size;
                double height_change = -size_diff * tree->get_ploidy();
                for (auto ht_ptr : tree->get_node_height_pointers()) {
                    double new_height = ht_ptr->get_value() + height_change;
                    if (new_height < 0.0) {
                        return -std::numeric_limits<double>::infinity();
                    }
                    ht_ptr->set_value(new_height);
                }
                return std::log(multiplier);
            }
            // If we reach this point, the pop sizes are unconstrained.  In
            // this case, we will scale the smallest pop size in the tree
            // (ignoring leaf sizes), and then change all node heights and
            // other sizes according to expectations under the coalescent.
            // It's a bit convoluted, but by doing this, we are least likely
            // to propose negative values for any parameters (though it is
            // still possible for node heights). An alternative would be to
            // simply scale all the node heights and then adjust their mapped
            // nodes according to coalescent expectations, but this route seems
            // much more likely to propose negative pop size values (which
            // could be remedied by auto-tuning the magnitude of the proposals,
            // but scaling the smallest pop size should limit how much of this
            // is needed). The old code that scaled node heights is commented
            // out below
            double smallest_size = std::numeric_limits<double>::max();
            for (auto node : tree->pre_ordered_nodes_) {
                if (! node->is_leaf()) {
                    if (node->get_population_size() < smallest_size) {
                        smallest_size = node->get_population_size();
                    }
                }
            }
            double old_size = smallest_size;
            double multiplier = this->op_.get_move_amount(rng);
            double new_size = old_size * multiplier;
            double size_change = new_size - old_size;
            double height_change = -size_change * tree->get_ploidy();
            for (unsigned int i = 0;
                    i < tree->get_number_of_node_heights();
                    ++i) {
                double new_height = tree->get_height(i) + height_change;
                if (new_height < 0.0) {
                    return -std::numeric_limits<double>::infinity();
                }
                tree->node_heights_.at(i)->set_value(new_height);
            }
            for (auto node : tree->pre_ordered_nodes_) {
                if (! node->is_leaf()) {
                    double new_node_size = node->get_population_size() + size_change;
                    if (new_node_size <= 0.0) {
                        return -std::numeric_limits<double>::infinity();
                    }
                    node->set_population_size(new_node_size);
                }
            }
            return std::log(multiplier);

            // Old code that scaled all the node heights and then adjusted all
            // the mapped nodes according to coalescent expectations.  Scaling
            // the smallest pop size instead (as above), likely proposes less
            // negative values, and thus requires less autotuning to work well.
            // std::vector<double> old_heights = tree->get_node_heights();
            // unsigned int num_heights = old_heights.size();
            // double multiplier = this->op_.get_move_amount(rng);
            // tree->scale_tree(multiplier);
            // std::vector<double> new_heights = tree->get_node_heights();
            // // Change internal node pop sizes according to time changes
            // for (unsigned int i = 0; i < new_heights.size(); ++i) {
            //     double height_diff = new_heights.at(i) - old_heights.at(i);
            //     double pop_size_change = -height_diff / tree->get_ploidy();
            //     std::vector< std::shared_ptr<PopulationNode> > mapped_nodes = tree->get_mapped_nodes(i);
            //     for (auto node : mapped_nodes) {
            //         double pop_size = node->get_population_size();
            //         double new_size = pop_size + pop_size_change;
            //         if (new_size <= 0.0) {
            //             return -std::numeric_limits<double>::infinity();
            //         }
            //         node->set_population_size(pop_size + pop_size_change);
            //     }
            // }
            // return std::log(multiplier) * num_heights;
        }
};

class HeightSizeMixer : public GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp> {
    public:
        HeightSizeMixer() : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>() { }
        HeightSizeMixer(double weight) : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>(weight) { }
        HeightSizeMixer(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "HeightSizeMixer";
        }

        std::string target_parameter() const {
            return "node heights";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::node_height_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::node_height;
        }

        bool is_operable(BasePopulationTree * tree) const {
            if (tree->population_sizes_are_fixed() || tree->population_sizes_are_constrained()) {
                // It doesn't make sense to use this move if the pop sizes are
                // constrained or fixed
                return false;
            }
            unsigned int num_heights = tree->get_number_of_node_heights();
            if (num_heights < 2) {
                // No non-root heights to operate on
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         *
         * The expected height of a gene tree node given the height of the
         * species tree node (t) is:
         *
         *   t*mu + ploidy*N*mu
         *
         * where mu is the mutation rate, and N is the effective size of the
         * ancestral population.
         *
         * In this proposal we want to mainatin this relationship:
         *
         *   t'*mu + ploidy*N'*mu = t*mu + ploidy*N*mu
         *
         * where (') represents the proposed (new) values. The mutation rates
         * cancel, giving us:
         *
         *   t' + ploidy*N' = t + ploidy*N
         *
         * Re-writing this in terms of the amount of change in t and N (denoted
         * as dt and dN, respectively), we have:
         *
         *   t + dt + ploidy(N + dN) = t + ploidy*N
         *
         * multiply out N and dN on the left:
         *
         *   t + dt + ploidy*N + ploidy*dN = t + ploidy*N
         *
         * subtract ploidy*N from both sides:
         *
         *   t + dt + ploidy*dN = t
         *
         * subtract t from both sides:
         *
         *   dt + ploidy*dN = 0
         *
         *   ploidy*dN = -dt
         *
         *   dN = -dt / ploidy
         *
         * So, we will scale one node heights, and then add (-dt / ploidy) to
         * the pop sizes of all nodes that map to it.
         */
        double propose(RandomNumberGenerator& rng,
                BasePopulationTree * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            unsigned int num_heights = tree->get_number_of_node_heights();
            unsigned int max_height_index = num_heights - 2;
            unsigned int height_index = rng.uniform_positive_int(
                    max_height_index);

            // Pop sizes are unconstrained, so will pick a node height, find
            // the smallest pop size mapped to it and scale it. Then we will
            // update the node height and the pop sizes of all other nodes
            // mapped to it according to coalescent expectations.
            std::vector< std::shared_ptr<PopulationNode> > mapped_nodes = tree->get_mapped_nodes(
                    height_index);
            double smallest_size = std::numeric_limits<double>::max();
            for (auto node : mapped_nodes) {
                if (node->get_population_size() < smallest_size) {
                    smallest_size = node->get_population_size();
                }
            }
            double old_size = smallest_size;
            double multiplier = this->op_.get_move_amount(rng);
            double new_size = old_size * multiplier;
            double size_change = new_size - old_size;
            double height_change = -size_change * tree->get_ploidy();
            for (auto node : mapped_nodes) {
                double new_node_size = node->get_population_size() + size_change;
                if (new_node_size <= 0.0) {
                    return -std::numeric_limits<double>::infinity();
                }
                node->set_population_size(new_node_size);
            }
            double new_height = tree->get_height(height_index) + height_change;
            if (new_height < tree->get_height_of_oldest_child(height_index)) {
                return -std::numeric_limits<double>::infinity();
            }
            if (new_height > tree->get_height_of_youngest_parent(height_index)) {
                return -std::numeric_limits<double>::infinity();
            }
            tree->set_height(height_index, new_height);

            return std::log(multiplier);
        }
};

class RootHeightSizeMixer : public GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp> {
    public:
        RootHeightSizeMixer() : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>() { }
        RootHeightSizeMixer(double weight) : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>(weight) { }
        RootHeightSizeMixer(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "RootHeightSizeMixer";
        }

        std::string target_parameter() const {
            return "root height";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::root_height_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::root_height;
        }

        bool is_operable(BasePopulationTree * tree) const {
            if (tree->root_height_is_fixed()) {
                return false;
            }
            if (tree->population_sizes_are_fixed() || tree->population_sizes_are_constrained()) {
                // It doesn't make sense to use this move if the pop sizes are
                // constrained or fixed
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         *
         * The expected height of a gene tree node given the height of the
         * species tree node (t) is:
         *
         *   t*mu + ploidy*N*mu
         *
         * where mu is the mutation rate, and N is the effective size of the
         * ancestral population.
         *
         * In this proposal we want to mainatin this relationship:
         *
         *   t'*mu + ploidy*N'*mu = t*mu + ploidy*N*mu
         *
         * where (') represents the proposed (new) values. The mutation rates
         * cancel, giving us:
         *
         *   t' + ploidy*N' = t + ploidy*N
         *
         * Re-writing this in terms of the amount of change in t and N (denoted
         * as dt and dN, respectively), we have:
         *
         *   t + dt + ploidy(N + dN) = t + ploidy*N
         *
         * multiply out N and dN on the left:
         *
         *   t + dt + ploidy*N + ploidy*dN = t + ploidy*N
         *
         * subtract ploidy*N from both sides:
         *
         *   t + dt + ploidy*dN = t
         *
         * subtract t from both sides:
         *
         *   dt + ploidy*dN = 0
         *
         *   ploidy*dN = -dt
         *
         *   dN = -dt / ploidy
         *
         * So, we will scale one node heights, and then add (-dt / ploidy) to
         * the pop sizes of all nodes that map to it.
         */
        double propose(RandomNumberGenerator& rng,
                BasePopulationTree * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            double old_size = tree->get_root_population_size();
            double multiplier = this->op_.get_move_amount(rng);
            double new_size = old_size * multiplier;
            if (new_size <= 0.0) {
                return -std::numeric_limits<double>::infinity();
            }
            double size_change = new_size - old_size;
            double height_change = -size_change * tree->get_ploidy();
            double new_height = tree->get_root_height() + height_change;
            unsigned int root_index = tree->get_number_of_node_heights() - 1;
            if (new_height < tree->get_height_of_oldest_child(root_index)) {
                return -std::numeric_limits<double>::infinity();
            }
            tree->set_root_population_size(new_size);
            tree->set_root_height(new_height);
            return std::log(multiplier);
        }
};


class HeightSizeSlideBumpMixer : public GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp> {
    protected:
        virtual bool call_tree_method_(
                BasePopulationTree * tree,
                RandomNumberGenerator& rng,
                unsigned int height_index,
                double height) {
            return tree->slide_bump_height(rng,
                    height_index,
                    height);
        }

        bool operate_on_root_ = false;

    public:
        HeightSizeSlideBumpMixer() : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>() { }
        HeightSizeSlideBumpMixer(double weight) : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>(weight) { }
        HeightSizeSlideBumpMixer(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "HeightSizeSlideBumpMixer";
        }

        std::string target_parameter() const {
            return "node heights and population sizes";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::node_height_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::node_height;
        }

        virtual void set_operate_on_root(bool operate_on_root) {
            this->operate_on_root_ = operate_on_root;
        }

        bool is_operable(BasePopulationTree * tree) const {
            if (tree->population_sizes_are_fixed() || tree->population_sizes_are_constrained()) {
                // It doesn't make sense to use the move if the pop sizes are
                // constrained or fixed
                return false;
            }
            unsigned int num_heights = tree->get_number_of_node_heights();
            if ((! this->operate_on_root_) && (num_heights < 2)) {
                // No non-root heights to operate on
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                BasePopulationTree * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            unsigned int num_heights = tree->get_number_of_node_heights();
            unsigned int max_height_index = num_heights - 2;
            if (this->operate_on_root_) {
                max_height_index = num_heights - 1;
            }
            unsigned int height_index = rng.uniform_int(0,
                    max_height_index);
            double new_height = tree->get_height(height_index);
            double ln_multiplier;
            this->update(rng, new_height, ln_multiplier);
            if (new_height < 0.0) {
                return -std::numeric_limits<double>::infinity();
            }
            if (new_height > tree->get_root_height()) {
                if (tree->root_height_is_fixed() || (! this->operate_on_root_)) {
                    return -std::numeric_limits<double>::infinity();
                }
            }
            bool move_happened = this->call_tree_method_(
                    tree,
                    rng,
                    height_index,
                    new_height);
            if (! move_happened) {
                return -std::numeric_limits<double>::infinity();
            }
            // Change internal node pop sizes according to time changes
            for (unsigned int i = 0; i < max_height_index + 1; ++i) {
                double height_diff = tree->get_height(i) - tree->get_stored_height(i);
                if (almost_equal_abs(height_diff, 0.0, 1e-10)) {
                    continue;
                }
                double pop_size_change = -height_diff / tree->get_ploidy();
                std::vector< std::shared_ptr<PopulationNode> > mapped_nodes = tree->get_mapped_nodes(i);
                for (auto node : mapped_nodes) {
                    double new_pop_size = node->get_population_size() + pop_size_change;
                    if (new_pop_size <= 0.0) {
                        return -std::numeric_limits<double>::infinity();
                    }
                    node->set_population_size(new_pop_size);
                }
            }
            return ln_multiplier;
        }
};


class GlobalHeightSizeRateScaler : public GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp> {

    protected:
        bool scale_sizes_ = true;
        bool scale_rate_ = true;

    public:
        GlobalHeightSizeRateScaler() : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>() { }
        GlobalHeightSizeRateScaler(double weight) : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>(weight) { }
        GlobalHeightSizeRateScaler(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<BasePopulationTree, ScaleOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "GlobalHeightSizeRateScaler";
        }

        std::string target_parameter() const {
            return "node heights, population sizes, and mutation rate";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::global_height_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::global;
        }

        bool is_operable(BasePopulationTree * tree) const {
            if (tree->root_height_is_fixed()) {
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                BasePopulationTree * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            double multiplier = this->op_.get_move_amount(rng);
            tree->scale_tree(multiplier);
            int num_params_scaled = tree->get_number_of_node_heights();
            int num_params_inv_scaled = 0;
            if (this->scale_sizes_) {
                num_params_inv_scaled = tree->scale_all_population_sizes(1.0/multiplier);
            }
            if (this->scale_rate_ and (! tree->mutation_rate_is_fixed())) {
                tree->set_mutation_rate(
                        tree->get_mutation_rate() * (1.0/multiplier));
                ++num_params_inv_scaled;
            }
            return std::log(multiplier) * (num_params_scaled - num_params_inv_scaled);
        }
};


class GlobalHeightRateScaler : public GlobalHeightSizeRateScaler {

    public:
        GlobalHeightRateScaler() : GlobalHeightSizeRateScaler() {
            this->scale_sizes_ = false;
        }
        GlobalHeightRateScaler(double weight) : GlobalHeightSizeRateScaler(weight) {
            this->scale_sizes_ = false;
        }
        GlobalHeightRateScaler(double weight, double tuning_parameter)
            : GlobalHeightSizeRateScaler(weight, tuning_parameter) { }

        std::string get_name() const {
            return "GlobalHeightRateScaler";
        }

        std::string target_parameter() const {
            return "node heights and mutation rate";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::global_height_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::global;
        }
};


class GlobalHeightSizeScaler : public GlobalHeightSizeRateScaler {

    public:
        GlobalHeightSizeScaler() : GlobalHeightSizeRateScaler() {
            this->scale_rate_ = false;
        }
        GlobalHeightSizeScaler(double weight) : GlobalHeightSizeRateScaler(weight) {
            this->scale_rate_ = false;
        }
        GlobalHeightSizeScaler(double weight, double tuning_parameter)
            : GlobalHeightSizeRateScaler(weight, tuning_parameter) { }

        std::string get_name() const {
            return "GlobalHeightSizeScaler";
        }

        std::string target_parameter() const {
            return "node heights and population sizes";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::global_height_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::global;
        }
};


class StateFreqMover : public GeneralTreeOperatorInterface<BasePopulationTree, WindowOp> {
    public:
        StateFreqMover() : GeneralTreeOperatorInterface<BasePopulationTree, WindowOp>() { }
        StateFreqMover(double weight) : GeneralTreeOperatorInterface<BasePopulationTree, WindowOp>(weight) { }
        StateFreqMover(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<BasePopulationTree, WindowOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "StateFreqMover";
        }

        std::string target_parameter() const {
            return "freq 1";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::derived_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::global;
        }

        bool is_operable(BasePopulationTree * tree) const {
            if (tree->state_frequencies_are_fixed()) {
                return false;
            }
            if (tree->state_frequencies_are_constrained()) {
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                BasePopulationTree * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            double freq_1 = tree->get_freq_1();
            double ln_hastings;
            this->update(rng, freq_1, ln_hastings);
            if ((freq_1 <= 0.0) || (freq_1 >= 1.0)) {
                return -std::numeric_limits<double>::infinity();
            }
            tree->set_freq_1(freq_1);
            return ln_hastings; 
        }
};


class StateFreqDirichletOperator : public GeneralTreeOperatorInterface<BasePopulationTree, DirichletOp> {
    public:
        StateFreqDirichletOperator() : GeneralTreeOperatorInterface<BasePopulationTree, DirichletOp>() { }
        StateFreqDirichletOperator(double weight) : GeneralTreeOperatorInterface<BasePopulationTree, DirichletOp>(weight) { }
        StateFreqDirichletOperator(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<BasePopulationTree, DirichletOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "StateFreqDirichletOperator";
        }

        std::string target_parameter() const {
            return "freq 1";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::derived_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::global;
        }

        bool is_operable(BasePopulationTree * tree) const {
            if (tree->state_frequencies_are_fixed()) {
                return false;
            }
            if (tree->state_frequencies_are_constrained()) {
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                BasePopulationTree * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            double freq_1 = tree->get_freq_1();
            double ln_hastings;
            this->update(rng, freq_1, ln_hastings);
            if ((freq_1 <= 0.0) || (freq_1 >= 1.0)) {
                return -std::numeric_limits<double>::infinity();
            }
            tree->set_freq_1(freq_1);
            return ln_hastings; 
        }
};


class StateFreqBetaOperator : public GeneralTreeOperatorInterface<BasePopulationTree, BetaOp> {
    public:
        StateFreqBetaOperator() : GeneralTreeOperatorInterface<BasePopulationTree, BetaOp>() { }
        StateFreqBetaOperator(double weight) : GeneralTreeOperatorInterface<BasePopulationTree, BetaOp>(weight) { }
        StateFreqBetaOperator(double weight, double tuning_parameter)
            : GeneralTreeOperatorInterface<BasePopulationTree, BetaOp>(weight, tuning_parameter) { }

        std::string get_name() const {
            return "StateFreqBetaOperator";
        }

        std::string target_parameter() const {
            return "freq 1";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::derived_operator;
        }
        BaseGeneralTreeOperatorTemplate::OperatorScopeEnum get_scope() const {
            return BaseGeneralTreeOperatorTemplate::OperatorScopeEnum::global;
        }

        bool is_operable(BasePopulationTree * tree) const {
            if (tree->state_frequencies_are_fixed()) {
                return false;
            }
            if (tree->state_frequencies_are_constrained()) {
                return false;
            }
            return true;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                BasePopulationTree * tree,
                unsigned int nthreads = 1) {
            if (! this->is_operable(tree)) {
                this->ignore_proposal_attempt_ = true;
                return -std::numeric_limits<double>::infinity();
            }
            double freq_1 = tree->get_freq_1();
            double ln_hastings;
            this->update(rng, freq_1, ln_hastings);
            if ((freq_1 <= 0.0) || (freq_1 >= 1.0)) {
                return -std::numeric_limits<double>::infinity();
            }
            tree->set_freq_1(freq_1);
            return ln_hastings; 
        }
};

#endif
