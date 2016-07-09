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
#include <cmath>
#include <limits>
#include <memory>

#include "collection.hpp"
#include "rng.hpp"
#include "assert.hpp"
#include "math_util.hpp"
#include "operator_schedule.hpp"


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
            model_operator = 3,
            rj_operator = 4
        };

        virtual Operator::OperatorTypeEnum get_type() const = 0;

        virtual void optimize(OperatorSchedule& os, double log_alpha) = 0;

        double get_target_acceptance_probability() const;

        virtual double get_coercable_parameter_value() const;

        virtual void set_coercable_parameter_value(double value) { }

        double get_weight() const { return this->weight_; }

        void set_weight(double weight);

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const = 0;

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons) = 0;
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const = 0;
        virtual double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& node_height) const = 0;

        void accept(const OperatorSchedule& os);

        void reject(const OperatorSchedule& os);

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

        virtual std::string get_name() const;

        virtual std::string target_parameter() const = 0;

        std::string header_string() const;

        std::string to_string(const OperatorSchedule& os) const;

    protected:
        double weight_ = 1.0;
        unsigned int number_rejected_ = 0;
        unsigned int number_accepted_ = 0;
        unsigned int number_rejected_for_correction_ = 0;
        unsigned int number_accepted_for_correction_ = 0;

        double calc_delta(OperatorSchedule& os, double log_alpha) const;
};

class ScaleOperator : public Operator {
    protected:
        double scale_ = 0.5;

    public:
        ScaleOperator() : Operator() { }
        ScaleOperator(double weight) : Operator(weight) { }
        ScaleOperator(double weight, double scale);
        virtual ~ScaleOperator() { }

        void set_scale(double scale);

        double get_scale() const {
            return this->scale_;
        }

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const;

        void optimize(OperatorSchedule& os, double log_alpha);

        double get_coercable_parameter_value() const;

        void set_coercable_parameter_value(double value);

        std::string get_name() const;
};

class WindowOperator : public Operator {
    protected:
        double window_size_ = 0.1;

    public:
        WindowOperator() : Operator() { }
        WindowOperator(double weight) : Operator(weight) { }
        WindowOperator(double weight, double window_size);
        virtual ~WindowOperator() { }

        void set_window_size(double window_size);
        double get_window_size() const;

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const;

        void optimize(OperatorSchedule& os, double log_alpha);

        double get_coercable_parameter_value() const;

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
                ComparisonPopulationTreeCollection& comparisons) = 0;
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const {
            throw EcoevolityError("calling wrong propose signature");
        }
        double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& node_height) const {
            throw EcoevolityError("calling wrong propose signature");
        }

        std::string get_name() const;

        std::string target_parameter() const;

        void optimize(OperatorSchedule& os, double log_alpha) { }
        void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const { }
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
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons) {
            throw EcoevolityError("calling wrong propose signature");
        }
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const = 0;
        double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& node_height) const {
            throw EcoevolityError("calling wrong propose signature");
        }

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
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons) {
            throw EcoevolityError("calling wrong propose signature");
        }
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const = 0;
        double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& node_height) const {
            throw EcoevolityError("calling wrong propose signature");
        }

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
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons) {
            throw EcoevolityError("calling wrong propose signature");
        }
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const {
            throw EcoevolityError("calling wrong propose signature");
        }
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
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons) {
            throw EcoevolityError("calling wrong propose signature");
        }
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const {
            throw EcoevolityError("calling wrong propose signature");
        }
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
                ComparisonPopulationTreeCollection& comparisons);
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const {
            throw EcoevolityError("calling wrong propose signature");
        }
        double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& node_height) const {
            throw EcoevolityError("calling wrong propose signature");
        }

        std::string target_parameter() const;

        std::string get_name() const;
};


// TODO:
// NOTE: Failing to sample from the prior with no data; use at your own peril
// (or just use UScaler).
class UMover : public ComparisonTreeWindowOperator {

    using ComparisonTreeWindowOperator::propose;

    public:
        UMover() : ComparisonTreeWindowOperator() { }
        UMover(double weight) : ComparisonTreeWindowOperator(weight) { }
        UMover(double weight, double window_size) : ComparisonTreeWindowOperator(weight, window_size) { }
        virtual ~UMover() { }

        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const;

        std::string target_parameter() const;

        std::string get_name() const;
};

class UScaler : public ComparisonTreeScaleOperator {

    using ComparisonTreeScaleOperator::propose;

    public:
        UScaler() : ComparisonTreeScaleOperator() { }
        UScaler(double weight) : ComparisonTreeScaleOperator(weight) { }
        UScaler(double weight, double scale) : ComparisonTreeScaleOperator(weight, scale) { }
        virtual ~UScaler() { }

        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const;

        std::string target_parameter() const;

        std::string get_name() const;
};

class ComparisonMutationRateScaler : public ComparisonTreeScaleOperator {

    using ComparisonTreeScaleOperator::propose;

    public:
        ComparisonMutationRateScaler() : ComparisonTreeScaleOperator() { }
        ComparisonMutationRateScaler(double weight) : ComparisonTreeScaleOperator(weight) { }
        ComparisonMutationRateScaler(double weight, double scale) : ComparisonTreeScaleOperator(weight, scale) { }
        virtual ~ComparisonMutationRateScaler() { }

        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const;

        std::string target_parameter() const;

        std::string get_name() const;
};


class ChildPopulationSizeScaler : public ComparisonTreeScaleOperator {

    using ComparisonTreeScaleOperator::propose;

    public:
        ChildPopulationSizeScaler() : ComparisonTreeScaleOperator() { }
        ChildPopulationSizeScaler(double weight) : ComparisonTreeScaleOperator(weight) { }
        ChildPopulationSizeScaler(double weight, double scale) : ComparisonTreeScaleOperator(weight, scale) { }
        virtual ~ChildPopulationSizeScaler() { }

        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const;

        std::string target_parameter() const;

        std::string get_name() const;
};


class RootPopulationSizeScaler : public ComparisonTreeScaleOperator {

    using ComparisonTreeScaleOperator::propose;

    public:
        RootPopulationSizeScaler() : ComparisonTreeScaleOperator() { }
        RootPopulationSizeScaler(double weight) : ComparisonTreeScaleOperator(weight) { }
        RootPopulationSizeScaler(double weight, double scale) : ComparisonTreeScaleOperator(weight, scale) { }
        virtual ~RootPopulationSizeScaler() { }

        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const;

        std::string target_parameter() const;

        std::string get_name() const;
};

class ComparisonHeightScaler : public NodeHeightScaleOperator {

    using NodeHeightScaleOperator::propose;

    public:
        ComparisonHeightScaler() : NodeHeightScaleOperator() { }
        ComparisonHeightScaler(double weight) : NodeHeightScaleOperator(weight) { }
        ComparisonHeightScaler(double weight, double scale) : NodeHeightScaleOperator(weight, scale) { }
        virtual ~ComparisonHeightScaler() { }

        double propose(
                RandomNumberGenerator& rng,
                PositiveRealParameter& node_height) const;

        std::string get_name() const;
};

class ComparisonHeightMover : public NodeHeightWindowOperator {

    using NodeHeightWindowOperator::propose;

    public:
        ComparisonHeightMover() : NodeHeightWindowOperator() { }
        ComparisonHeightMover(double weight) : NodeHeightWindowOperator(weight) { }
        ComparisonHeightMover(double weight, double window_size) : NodeHeightWindowOperator(weight, window_size) { }
        virtual ~ComparisonHeightMover() { }

        double propose(
                RandomNumberGenerator& rng,
                PositiveRealParameter& node_height) const;

        std::string get_name() const;
};

class DirichletProcessGibbsSampler : public ModelOperator {

    using ModelOperator::propose;

    protected:
        unsigned int number_of_auxiliary_categories_ = 4;

    public:
        DirichletProcessGibbsSampler() : ModelOperator() { }
        DirichletProcessGibbsSampler(double weight) : ModelOperator(weight) { }
        DirichletProcessGibbsSampler(double weight,
                unsigned int number_of_auxiliary_categories);
        virtual ~DirichletProcessGibbsSampler() { }

        void set_number_of_auxiliary_categories(unsigned int n);
        unsigned int get_number_of_auxiliary_categories() const;

        std::string get_name() const;

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons);
};

class ReversibleJumpSampler : public ModelOperator {

    using ModelOperator::propose;

    protected:
        std::map<unsigned int, std::vector<double> > split_subset_size_probs_;
        std::map<unsigned int, double> ln_number_of_possible_splits_;
        void populate_split_subset_size_probabilities(
                unsigned int number_of_nodes_in_event);

    public:
        ReversibleJumpSampler() : ModelOperator() { }
        ReversibleJumpSampler(double weight) : ModelOperator(weight) { }
        virtual ~ReversibleJumpSampler() { }

        std::string get_name() const;

        Operator::OperatorTypeEnum get_type() const;

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons);
        virtual double propose_jump_to_prior(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons);
        virtual double propose_jump_to_gap(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons);

        const std::vector<double>& get_split_subset_size_probabilities(
                unsigned int number_of_nodes_in_event);

        void write_split_probabilities(std::ostream& out) const {
            for (auto const & set_size_probs: this->split_subset_size_probs_) {
                out << "Set size: " << set_size_probs.first << "\n";
                double ln_ways_to_split = this->ln_number_of_possible_splits_.at(set_size_probs.first);
                double ways_to_split = std::exp(ln_ways_to_split);
                out << "\tlog number of ways to split: " << ln_ways_to_split << "\n";
                out << "\tnumber of ways to split: " << ways_to_split << "\n";
                out << "\tprobability of split subset sizes:\n";
                unsigned int split_size = 1;
                for (auto const & split_size_prob: set_size_probs.second) {
                    out << "\t\t" << split_size << ": " << split_size_prob << "\n";
                    ++split_size;
                }
            }
        }

};

#endif
