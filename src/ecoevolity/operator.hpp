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
        WindowOperator(double weight, double window_size) : Operator(weight);
        virtual ~WindowOperator() { }

        void set_window_size(double window_size);
        double get_window_size() const;

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const;

        void optimize(double log_alpha);

        double get_coercable_parameter_value();

        void set_coercable_parameter_value(double value);

        std::string get_name() const;
};


//////////////////////////////////////////////////////////////////////////////
// Operator Type base classes
//////////////////////////////////////////////////////////////////////////////

class ModelOperator : public Operator {
    public:
        ModelOperator() : Operator() { }
        ModelOperator(double weight) : Operator(weight) { }
        virtual ~ModelOperator() { }

        Operator::OperatorTypeEnum get_type() const;

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons) const = 0;

        std::string get_name() const;

        std::string target_parameter() const;
};

class ComparisonTreeScaleOperator : public ScaleOperator {
    public:
        ComparisonTreeScaleOperator() : ScaleOperator() { }
        ComparisonTreeScaleOperator(double weight) : ScaleOperator(weight) { }
        ComparisonTreeScaleOperator(double weight, double scale) : ScaleOperator(weight, scale) { }
        virtual ~ComparisonTreeScaleOperator() { }

        Operator::OperatorTypeEnum get_type() const;

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const = 0;

        std::string get_name() const;
};

class ComparisonTreeWindowOperator : public WindowOperator {
    public:
        ComparisonTreeWindowOperator() : WindowOperator() { }
        ComparisonTreeWindowOperator(double weight) : WindowOperator(weight) { }
        ComparisonTreeWindowOperator(double weight, double window_size) : WindowOperator(weight, window_size) { }
        virtual ~ComparisonTreeWindowOperator() { }

        Operator::OperatorTypeEnum get_type() const;

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const = 0;

        std::string get_name() const;
};

class NodeHeightScaleOperator : public ScaleOperator {
    public:
        NodeHeightScaleOperator() : ScaleOperator() { }
        NodeHeightScaleOperator(double weight) : ScaleOperator(weight) { }
        NodeHeightScaleOperator(double weight, double scale) : ScaleOperator(weight, scale) { }
        virtual ~NodeHeightScaleOperator() { }

        Operator::OperatorTypeEnum get_type() const;

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& node_height) const = 0;

        std::string get_name() const;

        std::string target_parameter() const;
};

class NodeHeightWindowOperator : public WindowOperator {
    public:
        NodeHeightWindowOperator() : WindowOperator() { }
        NodeHeightWindowOperator(double weight) : WindowOperator(weight) { }
        NodeHeightWindowOperator(double weight, double window_size) : WindowOperator(weight, window_size) { }
        virtual ~NodeHeightWindowOperator() { }

        Operator::OperatorTypeEnum get_type() const;

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& node_height) const = 0;

        std::string get_name() const;

        std::string target_parameter() const;
};


//////////////////////////////////////////////////////////////////////////////
// Derived Operator classes
//////////////////////////////////////////////////////////////////////////////

class ConcentrationScaler : public ScaleOperator {
    public:
        ConcentrationScaler() : ScaleOperator() { }
        ConcentrationScaler(double weight) : ScaleOperator(weight) { }
        ConcentrationScaler(double weight, double scale) : ScaleOperator(weight, scale) { }
        virtual ~ConcentrationScaler() { }

        Operator::OperatorTypeEnum get_type() const;

        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons) const;

        std::string target_parameter() const;

        std::string get_name() const;
};

class MutationRateMover : public ComparisonTreeWindowOperator {
    public:
        MutationRateMover() : ComparisonTreeWindowOperator() { }
        MutationRateMover(double weight) : ComparisonTreeWindowOperator(weight) { }
        MutationRateMover(double weight, double window_size) : ComparisonTreeWindowOperator(weight, window_size) { }
        virtual ~MutationRateMover() { }

        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const;

        std::string target_parameter() const;

        std::string get_name() const;
};

class ComparisonHeightMultiplierScaler : public ComparisonTreeScaleOperator {
    public:
        ComparisonHeightMultiplierScaler() : ComparisonTreeScaleOperator() { }
        ComparisonHeightMultiplierScaler(double weight) : ComparisonTreeScaleOperator(weight) { }
        ComparisonHeightMultiplierScaler(double weight, double scale) : ComparisonTreeScaleOperator(weight, scale) { }
        virtual ~ComparisonHeightMultiplierScaler() { }

        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const;

        std::string target_parameter() const;

        std::string get_name() const;
};

class ChildCoalescenceRateScaler : public ComparisonTreeScaleOperator {
    public:
        ChildCoalescenceRateScaler() : ComparisonTreeScaleOperator() { }
        ChildCoalescenceRateScaler(double weight) : ComparisonTreeScaleOperator(weight) { }
        ChildCoalescenceRateScaler(double weight, double scale) : ComparisonTreeScaleOperator(weight, scale) { }
        virtual ~ChildCoalescenceRateScaler() { }

        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const;

        std::string target_parameter() const;

        std::string get_name() const;
};

class RootCoalescenceRateScaler : public ChildCoalescenceRateScaler {
    public:
        RootCoalescenceRateScaler() : ChildCoalescenceRateScaler() { }
        RootCoalescenceRateScaler(double weight) : ChildCoalescenceRateScaler(weight) { }
        RootCoalescenceRateScaler(double weight, double scale) : ChildCoalescenceRateScaler(weight, scale) { }
        virtual ~RootCoalescenceRateScaler() { }

        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const;

        std::string get_name() const;
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
                PositiveRealParameter& node_height) const;

        std::string get_name() const;
};

class DirichletProcessGibbsSampler : public ModelOperator {
    public:
        DirichletProcessGibbsSampler() : ModelOperator() { }
        DirichletProcessGibbsSampler(double weight) : ModelOperator(weight) { }
        virtual ~DirichletProcessGibbsSampler() { }

        std::string get_name() const;

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons) const;
};

class ReversibleJumpSampler : public ModelOperator {
    public:
        ReversibleJumpSampler() : ModelOperator() { }
        ReversibleJumpSampler(double weight) : ModelOperator(weight) { }
        virtual ~ReversibleJumpSampler() { }

        std::string get_name() const;

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons) const;
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

        void add_operator(std::shared_ptr<Operator> o);

        Operator& draw_operator(RandomNumberGenerator& rng);

        double calc_delta(const Operator& op, double log_alpha);

        double get_total_weight() const;
        unsigned int get_auto_optimize_delay_count() const;
        unsigned int get_auto_optimize_delay() const;
        void set_auto_optimize_delay(unsigned int delay);

        void write_operator_rates(std::ofstream out);

        bool auto_optimizing() const;
        void turn_on_auto_optimize();
        void turn_off_auto_optimize();
};

#endif
