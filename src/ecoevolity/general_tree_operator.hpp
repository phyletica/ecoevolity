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
#include "rng.hpp"
#include "assert.hpp"
#include "math_util.hpp"


class BaseGeneralTreeOperatorTemplate {
    protected:
        double weight_ = 1.0;

    public:
		enum OperatorTypeEnum {
            node_height_operator = 1,
            topology_operator = 2,
            population_size_operator = 3,
            rj_operator = 4,
        };

        BaseGeneralTreeOperatorTemplate() { }
        BaseGeneralTreeOperatorTemplate(double weight) {
            this->set_weight(weight);
        }

        virtual BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const = 0;

        double get_weight() const {
            return this->weight_;
        }

        void set_weight(double weight) {
            ECOEVOLITY_ASSERT(weight >= 0.0);
            this->weight_ = weight;
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

        virtual std::string header_string() const = 0;

        virtual std::string to_string() const = 0;
};


template<class NodeType>
class GeneralTreeOperatorTemplate : public BaseGeneralTreeOperatorTemplate {

    public:
        GeneralTreeOperatorTemplate() : BaseGeneralTreeOperatorTemplate() { }
        GeneralTreeOperatorTemplate(double weight) : BaseGeneralTreeOperatorTemplate(weight) { }

        virtual void call_store_methods(
                BaseTree<NodeType> * tree) const {
            tree->store_state();
        }
        virtual void call_restore_methods(
                BaseTree<NodeType> * tree) const {
            tree->restore_state();
        }

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

        // void operate(RandomNumberGenerator& rng,
        //         BaseTree<NodeType> * tree,
        //         std::vector< std::shared_ptr< GeneralTreeOperatorInterface<NodeType> > > node_height_operators,
        //         unsigned int nthreads = 1) = 0;

};


template<class NodeType, class OperatorType>
class GeneralTreeOperatorInterface : public GeneralTreeOperatorTemplate<NodeType> {

    public:
        OperatorType op_;

        GeneralTreeOperatorInterface() : GeneralTreeOperatorTemplate<NodeType>() { }
        GeneralTreeOperatorInterface(double weight) : GeneralTreeOperatorTemplate<NodeType>(weight) { }

        void perform_move(
                RandomNumberGenerator& rng,
                BaseTree<NodeType> * tree,
                unsigned int nthreads) {
            this->call_store_methods(tree);
        
            double hastings_ratio = this->propose(rng, tree, nthreads);
            tree->compute_log_likelihood_and_prior(true);
        
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
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_probability)) {
                this->accept();
            }
            else {
                this->reject();
                this->call_restore_methods(tree);
            }
            tree->make_clean();
            if (this->auto_optimizing()) {
                this->optimize(acceptance_probability);
            }
        }

        void operate(RandomNumberGenerator& rng,
                BaseTree<NodeType> * tree,
                unsigned int nthreads = 1) {
            this->perform_move(rng, tree, nthreads);
        }

        // virtual void operate(RandomNumberGenerator& rng,
        //         BaseTree<NodeType> * tree,
        //         std::vector< std::shared_ptr< GeneralTreeOperatorInterface<NodeType, OperatorType> > > node_height_operators,
        //         unsigned int nthreads = 1) {
        //     this->operate(rng, tree, nthreads);
        // }

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


template<class NodeType>
class NodeHeightSlideBumpScaler : public GeneralTreeOperatorInterface<NodeType, ScaleOp> {
    protected:
        virtual bool call_tree_method_(
                BaseTree<NodeType> * tree,
                RandomNumberGenerator& rng,
                unsigned int height_index,
                double height) {
            return tree->slide_bump_height(rng,
                    height_index,
                    height);
        }

    public:
        NodeHeightSlideBumpScaler() : GeneralTreeOperatorInterface<NodeType, ScaleOp>() { }
        NodeHeightSlideBumpScaler(double weight) : GeneralTreeOperatorInterface<NodeType, ScaleOp>(weight) { }

        std::string get_name() const {
            return "NodeHeightSlideBumpScaler";
        }

        std::string target_parameter() const {
            return "node heights";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::node_height_operator;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                BaseTree<NodeType> * tree,
                unsigned int nthreads) {
            unsigned int max_height_index = tree->get_number_of_node_heights() - 1;
            if (tree->root_height_is_fixed()) {
                --max_height_index;
            }
            unsigned int height_index = rng.uniform_int(0,
                    max_height_index);
            double height = tree->get_height(height_index);
            double ln_multiplier;
            this->update(rng, height, ln_multiplier);
            if (height < 0) {
                return -std::numeric_limits<double>::infinity();
            }
            if (tree->root_height_is_fixed() && (height > tree->get_root_height())) {
                return -std::numeric_limits<double>::infinity();
            }
            bool move_happened = this->call_tree_method_(
                    tree,
                    rng,
                    height_index,
                    height);
            if (! move_happened) {
                return -std::numeric_limits<double>::infinity();
            }
            return ln_multiplier;
        }
};


template<class NodeType>
class NodeHeightSlideBumpPermuteScaler : public NodeHeightSlideBumpScaler<NodeType> {
    protected:
        bool call_tree_method_(
                BaseTree<NodeType> * tree,
                RandomNumberGenerator& rng,
                unsigned int height_index,
                double height) {
            return tree->slide_bump_permute_height(rng,
                    height_index,
                    height);
        }

    public:
        NodeHeightSlideBumpPermuteScaler() : NodeHeightSlideBumpScaler<NodeType>() { }
        NodeHeightSlideBumpPermuteScaler(double weight) : NodeHeightSlideBumpScaler<NodeType>(weight) { }

        std::string get_name() const {
            return "NodeHeightSlideBumpPermuteScaler";
        }
};


template<class NodeType>
class NodeHeightSlideBumpSwapScaler : public NodeHeightSlideBumpScaler<NodeType> {
    protected:
        bool call_tree_method_(
                BaseTree<NodeType> * tree,
                RandomNumberGenerator& rng,
                unsigned int height_index,
                double height) {
            return tree->slide_bump_swap_height(rng,
                    height_index,
                    height);
        }

    public:
        NodeHeightSlideBumpSwapScaler() : NodeHeightSlideBumpScaler<NodeType>() { }
        NodeHeightSlideBumpSwapScaler(double weight) : NodeHeightSlideBumpScaler<NodeType>(weight) { }

        std::string get_name() const {
            return "NodeHeightSlideBumpSwapScaler";
        }
};


template<class NodeType>
class NodeHeightSlideBumpMover : public GeneralTreeOperatorInterface<NodeType, WindowOp> {
    protected:
        virtual bool call_tree_method_(
                BaseTree<NodeType> * tree,
                RandomNumberGenerator& rng,
                unsigned int height_index,
                double height) {
            return tree->slide_bump_height(rng,
                    height_index,
                    height);
        }

    public:
        NodeHeightSlideBumpMover() : GeneralTreeOperatorInterface<NodeType, WindowOp>() { }
        NodeHeightSlideBumpMover(double weight) : GeneralTreeOperatorInterface<NodeType, WindowOp>(weight) { }

        std::string get_name() const {
            return "NodeHeightSlideBumpMover";
        }

        std::string target_parameter() const {
            return "node heights";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::node_height_operator;
        }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                BaseTree<NodeType> * tree,
                unsigned int nthreads) {
            unsigned int max_height_index = tree->get_number_of_node_heights() - 1;
            if (tree->root_height_is_fixed()) {
                --max_height_index;
            }
            unsigned int height_index = rng.uniform_int(0,
                    max_height_index);
            double height = tree->get_height(height_index);
            double ln_hastings;
            this->update(rng, height, ln_hastings);
            if (height < 0) {
                return -std::numeric_limits<double>::infinity();
            }
            if (tree->root_height_is_fixed() && (height > tree->get_root_height())) {
                return -std::numeric_limits<double>::infinity();
            }
            bool move_happened = this->call_tree_method_(
                    tree,
                    rng,
                    height_index,
                    height);
            if (! move_happened) {
                return -std::numeric_limits<double>::infinity();
            }
            return ln_hastings;
        }
};


template<class NodeType>
class NodeHeightSlideBumpPermuteMover : public NodeHeightSlideBumpMover<NodeType> {
    protected:
        bool call_tree_method_(
                BaseTree<NodeType> * tree,
                RandomNumberGenerator& rng,
                unsigned int height_index,
                double height) {
            return tree->slide_bump_permute_height(rng,
                    height_index,
                    height);
        }

    public:
        NodeHeightSlideBumpPermuteMover() : NodeHeightSlideBumpMover<NodeType>() { }
        NodeHeightSlideBumpPermuteMover(double weight) : NodeHeightSlideBumpMover<NodeType>(weight) { }

        std::string get_name() const {
            return "NodeHeightSlideBumpPermuteMover";
        }
};


template<class NodeType>
class NodeHeightSlideBumpSwapMover : public NodeHeightSlideBumpMover<NodeType> {
    protected:
        bool call_tree_method_(
                BaseTree<NodeType> * tree,
                RandomNumberGenerator& rng,
                unsigned int height_index,
                double height) {
            return tree->slide_bump_swap_height(rng,
                    height_index,
                    height);
        }

    public:
        NodeHeightSlideBumpSwapMover() : NodeHeightSlideBumpMover<NodeType>() { }
        NodeHeightSlideBumpSwapMover(double weight) : NodeHeightSlideBumpMover<NodeType>(weight) { }

        std::string get_name() const {
            return "NodeHeightSlideBumpSwapMover";
        }
};



template<class NodeType>
class NeighborHeightNodePermute : public GeneralTreeOperatorInterface<NodeType, Op> {

    public:
        NeighborHeightNodePermute() : GeneralTreeOperatorInterface<NodeType, Op>() { }
        NeighborHeightNodePermute(double weight) : GeneralTreeOperatorInterface<NodeType, Op>(weight) { }

        std::string get_name() const {
            return "NeighborHeightNodePermute";
        }

        std::string target_parameter() const {
            return "topology";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::topology_operator;
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
                    tree->get_number_of_node_heights() - 2);
            tree->collision_node_permute(rng,
                    height_index + 1,
                    height_index);
            return 0.0;
        }
};


template<class NodeType>
class NeighborHeightNodeSwap : public GeneralTreeOperatorInterface<NodeType, Op> {

    public:
        NeighborHeightNodeSwap() : GeneralTreeOperatorInterface<NodeType, Op>() { }
        NeighborHeightNodeSwap(double weight) : GeneralTreeOperatorInterface<NodeType, Op>(weight) { }

        std::string get_name() const {
            return "NeighborHeightNodeSwap";
        }

        std::string target_parameter() const {
            return "topology";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::topology_operator;
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
                    tree->get_number_of_node_heights() - 2);
            tree->collision_node_swap(rng,
                    height_index + 1,
                    height_index);
            return 0.0;
        }
};


template<class NodeType>
class SplitLumpNodesRevJumpSampler : public GeneralTreeOperatorInterface<NodeType, Op> {
    protected:
        std::map<unsigned int, double> ln_stirling2_numbers;
        std::map<unsigned int, long double> bell_numbers;

    public:
        SplitLumpNodesRevJumpSampler() : GeneralTreeOperatorInterface<NodeType, Op>() { }
        SplitLumpNodesRevJumpSampler(double weight) : GeneralTreeOperatorInterface<NodeType, Op>(weight) { }

        std::string get_name() const {
            return "SplitLumpNodesRevJumpSampler";
        }

        std::string target_parameter() const {
            return "topology";
        }

        BaseGeneralTreeOperatorTemplate::OperatorTypeEnum get_type() const {
            return BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::rj_operator;
        }

        double get_ln_stirling2(unsigned int n) {
            if (this->ln_stirling2_numbers.count(n) < 1) {
                this->ln_stirling2_numbers.at(n) = std::log(
                        stirling2_float(n, 2));
            }
            return this->ln_stirling2_numbers.at(n);
        }

        long double get_bell_number(unsigned int n) {
            if (this->bell_numbers.count(n) < 1) {
                this->bell_numbers.at(n) = bell_float(n);
            }
            return this->bell_numbers.at(n);
        }

        // void operate(RandomNumberGenerator& rng,
        //         BaseTree<NodeType> * tree,
        //         unsigned int nthreads = 1) {
        //     this->perform_move(rng, tree, nthreads);
        // 
        //     // Perform sweep of univariate time moves
        //     for (std::shared_ptr< GeneralTreeOperatorInterface<NodeType> > time_op : node_height_operators) {
        //         time_op->operate(rng, tree, nthreads);
        //     }
        // }

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                BaseTree<NodeType> * tree,
                unsigned int nthreads) {
            const unsigned int num_heights = tree->get_number_of_node_heights();
            const unsigned int num_internal_nodes = tree->get_internal_node_count();
            const unsigned int num_leaves = tree->get_leaf_node_count();
            ECOEVOLITY_ASSERT(num_leaves > 2);
            const bool in_comb_state_before = (num_heights == 1);
            const bool in_general_state_before = (num_heights == (num_leaves - 1));
            const bool split_event = ((! in_general_state_before) &&
                    (in_comb_state_before || (rng.uniform_real() < 0.5)));
            if (split_event) {
                std::vector<unsigned int> splittable_height_indices = tree->get_indices_of_splittable_heights();
                const unsigned int number_of_splittable_heights = splittable_height_indices.size();
                ECOEVOLITY_ASSERT(number_of_splittable_heights > 0);
                unsigned int split_height_idx = rng.uniform_int(0, number_of_splittable_heights - 1);
                double current_height = tree->get_height(split_height_idx);
                double height_lower_bound = -1.0;
                unsigned int number_of_mapped_nodes = 0;
                std::vector<unsigned int> mapped_polytomy_sizes;
                tree->split_node_height_down(rng,
                        split_height_idx,
                        height_lower_bound,
                        number_of_mapped_nodes,
                        mapped_polytomy_sizes,
                        false);
                ECOEVOLITY_ASSERT(number_of_mapped_nodes > 0);
                double height_diff = current_height - height_lower_bound;
                ECOEVOLITY_ASSERT(height_diff > 0.0);

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
                //      NOTE: This can also be thought of as the Jacobian term
                //      for the reversible jump move.
                //   4) randomly partitioning the children of each polytomy;
                //      each subset with 2 or more children is split off as a
                //      clade and mapped to the new, younger height. The prob
                //      for each polytomy
                //          = \prod 1 / (Bell(nchildren) - 1)
                //      where Bell(n) is the Bell number, and 1 is subtracted,
                //      because we reject the partition where every child is in
                //      its own subset.  We reject this, because then the node
                //      would node be split or moved to the new height. We have
                //      to take the product of this over all polytomy nodes
                //      that end up in the ``split'' subset of nodes.
                //
                //      In the special case where there is only one polytomy
                //      mapped to the node we are splitting, this needs to be:
                //          = 1 / (Bell(nchildren) - 2)
                //      The minus 2 is because we now reject 2 possible
                //      partitions: (1) where every child is in its own subset,
                //      and (2) when all children end up together in one subst.
                //      The second partition would result in the polytomy node
                //      simply moving down to the new height (not getting split
                //      up). This is normally OK, when other nodes are
                //      remaining at the old height (this is guaranteed,
                //      because we always split nodes into two non-empty
                //      subsets), but in the singleton case that doesn't
                //      happen.
                //
                // So the prob of the forward split move is
                // p(split move) =  1 / (number of splittable heights *
                //                       (2^n - 2) * d *
                //                       \prod 1 / (Bell(nchildren) - X))
                //               =  1 / (number of splittable heights * 
                //                       2 * stirling2(n, 2) * d *
                //                       \prod 1 / (Bell(nchildren) - X))
                //
                // The probability of the reverse move is simply the
                // probability of randomly selecting the proposed (split)
                // height from among all node heights except the root.
                // p(reverse merge) = 1 / (number of heights before proposal + 1 - 1)
                //                  = 1 / (number of heights before the proposal)
                //
                // So, the Hasting ratio for the proposed split move is:
                // p(reverse merge) / p(proposed split) = 
                //     (number of splittable heights * 2 * stirling2(n, 2) * d * \prod (Bell(nchild) - 1))
                //     ----------------------------------------------------------------------------------
                //          (number of heights before the proposal)
                //
                // EXCEPT for the case of a singleton polytomy, which is:
                // p(reverse merge) / p(proposed split) = 
                //     (number of splittable heights * d * (Bell(n) - 2))
                //     --------------------------------------------------
                //         (number of heights before the proposal)

                double ln_hastings = 0.0;
                if (number_of_mapped_nodes == 1) {
                    // We have the special case of only a single polytomy
                    ECOEVOLITY_ASSERT(mapped_polytomy_sizes.size() == 1);
                    double ln_bell_num_minus_2 = std::log(
                            this->get_bell_number(mapped_polytomy_sizes.at(0)) - 2.0);
                    ln_hastings =
                            std::log(number_of_splittable_heights) +
                            std::log(height_diff) +
                            ln_bell_num_minus_2 -
                            std::log(num_heights);
                }
                else {
                    // We have multiple nodes mapped to height
                    double ln_stirling2 = this->get_ln_stirling2(number_of_mapped_nodes);
                    double ln_bell_num_minus_1_sum = 0.0;
                    for (auto polytomy_size : mapped_polytomy_sizes) {
                        ln_bell_num_minus_1_sum += std::log(this->get_bell_number(polytomy_size) - 1.0);
                    }
                    ln_hastings =
                            std::log(number_of_splittable_heights) +
                            std::log(2.0) +
                            ln_stirling2 +
                            std::log(height_diff) +
                            ln_bell_num_minus_1_sum -
                            std::log(num_heights);
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
                bool in_general_state_after = (tree->get_number_of_heights() ==
                        (tree->get_leaf_node_count() - 1));

                if (in_comb_state_before && (! in_general_state_after)) {
                    ln_hastings -= std::log(2.0);
                }
                else if (in_general_state_after && (! in_comb_state_before)) {
                    ln_hastings += std::log(2.0);
                }
                return ln_hastings;
            }
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
            // the merge is the product of the probabilities of:
            //  1) Randomly selecting the merged height from among the
            //     splittable heights that exist AFTER the merge move
            //       = 1 / number of splittable heights after merge
            //
            //  2) Randomly selecting the correct nodes to split down to a new
            //     height to get the pre-merge state. The number of ways you
            //     can split the $n$ nodes mapped to the merged height into two
            //     non-empy sets is equal to Stirling2(n, 2) (the Stirling
            //     number of the second kind). Then you have to select one
            //     of the two sets to move down to the new height, so
            //     the probability of this is
            //       = 1/2 * 1 / Stirling2(n, 2)
            //  3) The polytomy nodes selected to "split down" need to be
            //     resolved (or not) to their pre-merge state. For each
            //     node with k children the number of ways to split it is equal
            //     to the number of ways the k children can be partitioned into
            //     1--k-1 non-empty sets, or bell(k) - 1 (the Bell number).
            //     All subsets with more than 1 child is split into a new node
            //     that descends from the original polytomy node and is
            //     assigned the new younger height. The "minus 1" is 
            //     because we don't allow the partition where every child is
            //     assigned to its own subset, because this would result in
            //     nothing happening to the polytomy node. The probability
            //     of this is a product over each polytomy node in the
            //     "split subset"
            //      = \prod 1 / bell(k) - 1
            //  4) Randomly drawing the pre-merged younger height between
            //     the merged height and the next younger (post-merge) height
            //      = 1 / merged height - next younger height (or zero)
            //      = 1 / d
            // 
            // So, the overall probability of the reverse split move is
            //
            // = 1 / (post-merge num of splittable heights * 2 * Stirling2(n, 2) * d * \prod bell(k) - 1)
            //
            // However, this is not true in one corner case: when, after the
            // merge, there is only one node mapped to the merged height (which
            // is necessarily a polytomy). In this special case of a singleton
            // polytomy node mapped to the merged height, the probability of
            // the reverse move is:
            //
            // = 1 / (post-merge num of splittable heights * d * bell(k) - 2)
            //
            // The number of ways to split the mapped nodes into 2 non-empty
            // subsets is simply 1 (because there is only one node mapped to
            // the height. Also, we now subtract 2 from all the ways we can
            // partition the children of the polytomy into non-empty sets.
            // This is because, in addition to not allowing the partition where
            // all children are assigned to their own set, we also do not allow
            // the partition where they are all mapped to one subset (because
            // this would not break up the polytomy).
            //
            // The hastings ratio is then
            //
            // = prob of forward merge / prob of reverse split
            //
            // =
            // (post-merge num of splittable heights * 2 * Stirling2(n, 2) * d * \prod bell(k) - 1)
            // ------------------------------------------------------------------------------------
            //                 nheights after merge
            // 
            // EXCEPT if there is only one polytomy node mapped to the merged
            // height, in which case, the hastings ratio is:
            //
            // (post-merge num of splittable heights * d * \prod bell(k) - 2)
            // ---------------------------------------------------------------
            //                 nheights after merge

            unsigned int merge_height_idx = rng.uniform_int(0, num_heights - 2);
            const double younger_height = tree->get_height(merge_height_idx);
            tree->merge_node_height_up(merge_height_idx);
            const double older_height = tree->get_height(merge_height_idx);
            const double height_diff = older_height - younger_height;
            ECOEVOLITY_ASSERT(height_diff > 0.0);
            const unsigned int post_num_splittable_heights = tree->get_number_of_splittable_heights();
            std::vector< std::shared_ptr<NodeType> > post_mapped_nodes = tree->get_mapped_nodes(merge_height_idx);
            std::vector< std::shared_ptr<NodeType> > post_mapped_poly_nodes = tree->get_mapped_polytomy_nodes(merge_height_idx);
            const bool post_num_mapped_nodes = post_mapped_nodes.size();
            double ln_hastings = 0.0;
            if (post_num_mapped_nodes == 1) {
                const unsigned int num_polytomy_children = post_mapped_nodes.at(0)->get_number_of_children();
                ECOEVOLITY_ASSERT(num_polytomy_children > 2);
                ln_hastings =
                        std::log(post_num_splittable_heights) +
                        std::log(height_diff) +
                        std::log(this->get_bell_number(num_polytomy_children) - 2.0) -
                        std::log(num_heights - 1);

            }
            else {
                double ln_stirling2_num = this->get_ln_stirling2(post_num_mapped_nodes);
                double ln_bell_number_sum = 0.0;
                for (auto poly_node : post_mapped_poly_nodes) {
                    ln_bell_number_sum += std::log(
                            this->get_bell_number(poly_node->get_number_of_children()) - 1.0);
                }
                ln_hastings =
                        std::log(post_num_splittable_heights) +
                        std::log(2.0) +
                        ln_stirling2_num +
                        std::log(height_diff) +
                        ln_bell_number_sum -
                        std::log(num_heights - 1);
            }
            return 0.0;
            /* const unsigned int num_heights = tree->get_number_of_node_heights(); */
            /* const unsigned int num_internal_nodes = tree->get_internal_node_count(); */
            /* const unsigned int num_leaves = tree->get_leaf_node_count(); */
            /* const bool in_comb_state_before = (num_heights == 1); */
            /* const bool in_general_state_before = (num_heights == (num_leaves - 1)); */
        }
};

#endif
