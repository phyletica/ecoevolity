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

#include "tree.hpp"
#include "rng.hpp"
#include "assert.hpp"
#include "math_util.hpp"


template<class NodeType>
class GeneralTreeOperatorInterface {
    protected:
        double weight_ = 1.0;

    public:
		enum OperatorTypeEnum {
            node_height_operator = 1,
            topology_operator = 2,
            population_size_operator = 3,
            rj_operator = 4,
        };

        GeneralTreeOperatorInterface() { }
        GeneralTreeOperatorInterface(double weight) {
            this->set_weight(weight);
        }

        virtual GeneralTreeOperatorInterface::OperatorTypeEnum get_type() const = 0;

        double get_weight() const {
            return this->weight_;
        }

        void set_weight(double weight) {
            ECOEVOLITY_ASSERT(weight >= 0.0);
            this->weight_ = weight;
        }

        virtual void call_store_methods(
                BaseTree<NodeType> * tree) const {
            tree->store_state();
        }
        virtual void call_restore_methods(
                BaseTree<NodeType> * tree) const {
            tree->restore_state();
        }

        virtual std::string get_name() const = 0;

        virtual std::string target_parameter() const = 0;

        virtual void perform_move(RandomNumberGenerator& rng,
                BaseTree<NodeType> * tree,
                unsigned int nthreads) = 0;

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                BaseTree<NodeType> * tree,
                unsigned int nthreads) = 0;

        virtual void operate(RandomNumberGenerator& rng,
                BaseTree<NodeType> * tree,
                unsigned int nthreads = 1) = 0;

        void operate(RandomNumberGenerator& rng,
                BaseTree<NodeType> * tree,
                std::vector< std::shared_ptr< BaseGenTreeOperatorInterface<NodeType> > > node_height_operators,
                unsigned int nthreads = 1) = 0;

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

        virtual std::string header_string() const = 0;

        virtual std::string to_string() const = 0;
};


template<class NodeType>
template<class OperatorType>
class BaseGenTreeOperatorInterface : public GeneralTreeOperatorInterface<NodeType> {

    public:
        OperatorType op_;

        BaseGenTreeOperatorInterface() : GeneralTreeOperatorInterface<NodeType>() { }
        BaseGenTreeOperatorInterface(double weight) : GeneralTreeOperatorInterface<NodeType>(weight) { }

        void perform_move(
                RandomNumberGenerator& rng,
                BaseTree<NodeType> * tree,
                unsigned int nthreads) {
            this->call_store_methods(tree);
        
            double hastings_ratio = this->propose(rng, tree, nthreads);
            tree->compute_log_likelihood_and_prior(true);
        
            double likelihood_ratio = 
                    tree->get_log_likelihood() -
                    tree->get_stored_log_likelihood();
            // The prior ratio is included in the hastings ratio
            double acceptance_probability =
                    likelihood_ratio + 
                    hastings_ratio;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_probability)) {
                this->accept();
            }
            else {
                this->reject();
                this->call_restore_methods(tree);
            }
            tree->make_clean();
            this->optimize(acceptance_probability);
        }

        void operate(RandomNumberGenerator& rng,
                BaseTree<NodeType> * tree,
                unsigned int nthreads = 1) {
            this->perform_move(rng, tree, nthreads);
        }

        void operate(RandomNumberGenerator& rng,
                BaseTree<NodeType> * tree,
                std::vector< std::shared_ptr< BaseGenTreeOperatorInterface<NodeType> > > node_height_operators,
                unsigned int nthreads = 1) {
            this->operate(rng, tree, nthreads);
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

        std::string header_string() const {
            return "name\tnumber_accepted\tnumber_rejected\tweight\ttuning_parameter\n";
        }

        virtual std::string to_string() const {
            std::ostringstream ss;
            ss << this->get_name() << "\t" 
               << this->op_.get_number_accepted() << "\t"
               << this->op_.get_number_rejected() << "\t"
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

class Operator {

    public:
        Operator(const bool auto_optimize = false) :
            auto_optimize_(auto_optimize) { }
        Operator(const unsigned int auto_optimize_delay,
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
            if ((this->get_number_of_attempts() < this->get_auto_optimize_delay()) ||
                    (! this->auto_optimize_)) {
                ++this->auto_optimize_delay_count_;
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
            if (os.get_number_of_attempts() >= os.get_auto_optimize_delay()) {
                ++this->number_accepted_for_correction_;
            }
        }

        void reject() {
            ++this->number_rejected_;
            if (os.get_number_of_attempts() >= os.get_auto_optimize_delay()) {
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

        unsigned int get_auto_optimize_delay() const {
            return this->get_auto_optimize_delay_;
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


class ScaleOperator : public Operator {

    protected:
        double scale_ = 0.5;

    public:
        ScaleOperator(const bool auto_optimize = true) : Operator(auto_optimize) { }
        ScaleOperator(const unsigned int auto_optimize_delay,
                const bool auto_optimize = true) :
                Operator(auto_optimize_delay, auto_optimize) { }
        ScaleOperator(double scale,
                const bool auto_optimize = true) :
                Operator(auto_optimize) {
            this->set_scale(scale);
        }
        ScaleOperator(double scale,
                const unsigned int auto_optimize_delay,
                const bool auto_optimize = true) :
                Operator(auto_optimize_delay, auto_optimize) {
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

        void optimize(double log_alpha) {
            double delta = this->calc_delta(log_alpha);
            delta += std::log(this->scale_);
            this->set_scale(std::exp(delta));
        }

        double get_coercable_parameter_value() const {
            return this->scale_;
        }

        void set_coercable_parameter_value(double value) {
            this->set_scale(value);
        }
};


class WindowOperator : public Operator {

    protected:
        double window_size_ = 0.1;

    public:
        WindowOperator(const bool auto_optimize = true) : Operator(auto_optimize) { }
        WindowOperator(const unsigned int auto_optimize_delay,
                const bool auto_optimize = true) :
                Operator(auto_optimize_delay, auto_optimize) { }
        WindowOperator(double window_size,
                const bool auto_optimize = true) :
                Operator(auto_optimize) {
            this->set_window_size(window_size);
        }
        WindowOperator(double window_size,
                const unsigned int auto_optimize_delay,
                const bool auto_optimize = true) :
                Operator(auto_optimize_delay, auto_optimize) {
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

        void optimize(double log_alpha) {
            double delta = this->calc_delta(log_alpha);
            delta += std::log(this->window_size_);
            this->set_window_size(std::exp(delta));
        }

        double get_coercable_parameter_value() const {
            return this->window_size_;
        }

        void set_coercable_parameter_value(double value) {
            this->set_window_size(value);
        }
};


template<class NodeType>
class NodeHeightSlideBumpScaler : public BaseGenTreeOperatorInterface<NodeType><ScaleOperator> {

    public:
        NodeHeightSlideBumpScaler() : BaseGenTreeOperatorInterface<ScaleOperator>() { }
        NodeHeightSlideBumpScaler(double weight) : BaseGenTreeOperatorInterface<ScaleOperator>(weight) { }

        std::string get_name() const {
            return "NodeHeightSlideBumpScaler";
        }

        std::string target_parameter() const {
            return "node heights";
        }

        GeneralTreeOperatorInterface::OperatorTypeEnum get_type() const {
            return GeneralTreeOperatorInterface::OperatorTypeEnum::node_height_operator;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                BaseTree<NodeType> * tree,
                unsigned int nthreads) {
            unsigned int height_index = rng.uniform_int(0,
                    tree->get_number_of_node_heights() - 1);
            double height = tree->get_height(height_index);
            double ln_multiplier;
            this->update(rng, heigh, ln_multiplier);
        }
};


template<class NodeType>
class SplitLumpNodesRevJumpSampler : public BaseGenTreeOperatorInterface<NodeType><Operator> {

    public:
        SplitLumpNodesRevJumpSampler() : BaseGenTreeOperatorInterface<Operator>() { }
        SplitLumpNodesRevJumpSampler(double weight) : BaseGenTreeOperatorInterface<Operator>(weight) { }

        std::string get_name() const {
            return "SplitLumpNodesRevJumpSampler";
        }

        std::string target_parameter() const {
            return "topology";
        }

        GeneralTreeOperatorInterface::OperatorTypeEnum get_type() const {
            return GeneralTreeOperatorInterface::OperatorTypeEnum::rj_operator;
        }

        void operate(RandomNumberGenerator& rng,
                BaseTree<NodeType> * tree,
                std::vector< std::shared_ptr< BaseGenTreeOperatorInterface<NodeType> > > node_height_operators,
                unsigned int nthreads = 1) {
            this->perform_move(rng, tree, nthreads);
        
            // Perform sweep of univariate time moves
            for (std::shared_ptr< BaseGenTreeOperatorInterface<NodeType> > time_op : node_height_operators) {
                time_op->operate(rng, tree, nthreads);
            }
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                BaseTree<NodeType> * tree,
                unsigned int nthreads) {
            // TODO
        }
};

#endif
