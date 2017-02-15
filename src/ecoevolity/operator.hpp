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


class OperatorInterface {
    protected:
        double weight_ = 1.0;

    public:
        OperatorInterface() { }
        OperatorInterface(double weight);
		enum OperatorTypeEnum {
            tree_operator = 1,
            time_operator = 2,
            collection_operator = 3,
            rj_operator = 4,
        };

        virtual OperatorInterface::OperatorTypeEnum get_type() const = 0;

        double get_weight() const;

        void set_weight(double weight);

        virtual void call_store_methods(
                ComparisonPopulationTreeCollection& comparisons) const;
        virtual void call_restore_methods(
                ComparisonPopulationTreeCollection& comparisons) const;

        virtual std::string get_name() const = 0;

        virtual std::string target_parameter() const = 0;

        virtual void perform_collection_move(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1) = 0;

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& parameter) const = 0;
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const = 0;
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1) = 0;

        virtual void operate(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1) = 0;

        virtual void optimize(OperatorSchedule& os, double log_alpha) = 0;

        virtual void accept(const OperatorSchedule& os) = 0;

        virtual void reject(const OperatorSchedule& os) = 0;

        virtual double get_coercable_parameter_value() const = 0;

        virtual void set_coercable_parameter_value(double value) = 0;

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const = 0;

        virtual double calc_delta(OperatorSchedule& os, double log_alpha) const = 0;

        virtual double get_target_acceptance_probability() const = 0;

        virtual unsigned int get_number_rejected() const = 0;
        virtual unsigned int get_number_accepted() const = 0;
        virtual unsigned int get_number_rejected_for_correction() const = 0;
        virtual unsigned int get_number_accepted_for_correction() const = 0;

        virtual std::string header_string() const = 0;

        virtual std::string to_string(const OperatorSchedule& os) const = 0;
};

template<class DerivedOperatorType>
class BaseOperatorInterface : public OperatorInterface {

    public:
        DerivedOperatorType op_;

        BaseOperatorInterface() : OperatorInterface() { }
        BaseOperatorInterface(double weight) : OperatorInterface(weight) { }

        virtual void optimize(OperatorSchedule& os, double log_alpha);

        void accept(const OperatorSchedule& os);

        void reject(const OperatorSchedule& os);

        virtual double get_coercable_parameter_value() const;

        virtual void set_coercable_parameter_value(double value);

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const;

        double calc_delta(OperatorSchedule& os, double log_alpha) const;

        double get_target_acceptance_probability() const;

        unsigned int get_number_rejected() const;
        unsigned int get_number_accepted() const;
        unsigned int get_number_rejected_for_correction() const;
        unsigned int get_number_accepted_for_correction() const;

        std::string header_string() const;

        virtual std::string to_string(const OperatorSchedule& os) const;
};


template<class DerivedOperatorType>
class TimeOperatorInterface : public BaseOperatorInterface<DerivedOperatorType> {

    public:
        TimeOperatorInterface() : BaseOperatorInterface<DerivedOperatorType>() { }
        TimeOperatorInterface(double weight) : BaseOperatorInterface<DerivedOperatorType>(weight) { }

        virtual void perform_collection_move(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& parameter) const = 0;

        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const {
            throw EcoevolityError("calling wrong propose signature");
        }
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1) {
            throw EcoevolityError("calling wrong propose signature");
        }

        virtual OperatorInterface::OperatorTypeEnum get_type() const;

        virtual void operate(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1) = 0;
};

template<class DerivedOperatorType>
class TreeOperatorInterface : public BaseOperatorInterface<DerivedOperatorType> {

    public:
        TreeOperatorInterface() : BaseOperatorInterface<DerivedOperatorType>() { }
        TreeOperatorInterface(double weight) : BaseOperatorInterface<DerivedOperatorType>(weight) { }

        virtual void perform_collection_move(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const = 0;

        virtual double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& parameter) const {
            throw EcoevolityError("calling wrong propose signature");
        }
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1) {
            throw EcoevolityError("calling wrong propose signature");
        }

        virtual OperatorInterface::OperatorTypeEnum get_type() const;

        virtual void operate(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1) = 0;
};

template<class DerivedOperatorType>
class CollectionOperatorInterface : public BaseOperatorInterface<DerivedOperatorType> {

    public:
        CollectionOperatorInterface() : BaseOperatorInterface<DerivedOperatorType>() { }
        CollectionOperatorInterface(double weight) : BaseOperatorInterface<DerivedOperatorType>(weight) { }

        virtual void perform_collection_move(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1) = 0;

        virtual double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& parameter) const {
            throw EcoevolityError("calling wrong propose signature");
        }
        virtual double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const {
            throw EcoevolityError("calling wrong propose signature");
        }

        virtual OperatorInterface::OperatorTypeEnum get_type() const;

        virtual void operate(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1) = 0;
};


//////////////////////////////////////////////////////////////////////////////
// Operator base classes
//////////////////////////////////////////////////////////////////////////////

class Operator {

    public:
        Operator() { }
        virtual ~Operator() { }

        virtual void optimize(OperatorSchedule& os, double log_alpha) { }

        virtual double get_coercable_parameter_value() const;

        virtual void set_coercable_parameter_value(double value) { }

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const { }

        double calc_delta(OperatorSchedule& os, double log_alpha) const;

        void accept(const OperatorSchedule& os);

        void reject(const OperatorSchedule& os);

        double get_target_acceptance_probability() const;

        unsigned int get_number_rejected() const;
        unsigned int get_number_accepted() const;
        unsigned int get_number_rejected_for_correction() const;
        unsigned int get_number_accepted_for_correction() const;

    protected:
        unsigned int number_rejected_ = 0;
        unsigned int number_accepted_ = 0;
        unsigned int number_rejected_for_correction_ = 0;
        unsigned int number_accepted_for_correction_ = 0;
};


class ScaleOperator : public Operator {

    protected:
        double scale_ = 0.5;

    public:
        ScaleOperator() : Operator() { }
        ScaleOperator(double scale);
        virtual ~ScaleOperator() { }

        void set_scale(double scale);

        double get_scale() const;

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const;

        void optimize(OperatorSchedule& os, double log_alpha);

        double get_coercable_parameter_value() const;

        void set_coercable_parameter_value(double value);
};


class WindowOperator : public Operator {

    protected:
        double window_size_ = 0.1;

    public:
        WindowOperator() : Operator() { }
        WindowOperator(double window_size);
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
};


//////////////////////////////////////////////////////////////////////////////
// Derived Operator classes
//////////////////////////////////////////////////////////////////////////////

class ConcentrationScaler : public CollectionOperatorInterface<ScaleOperator> {

    public:
        ConcentrationScaler();
        ConcentrationScaler(double weight);
        ConcentrationScaler(double weight, double scale);
        // virtual ~ConcentrationScaler() { }

        void operate(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& parameter) const {
            throw EcoevolityError("calling wrong propose signature");
        }
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const {
            throw EcoevolityError("calling wrong propose signature");
        }

        std::string target_parameter() const;

        std::string get_name() const;
};


class FreqMover : public TreeOperatorInterface<WindowOperator> {

    public:
        FreqMover();
        FreqMover(double weight);
        FreqMover(double weight, double window_size);
        // virtual ~FreqMover() { }

        void operate(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);
        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const;

        double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& parameter) const {
            throw EcoevolityError("calling wrong propose signature");
        }
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1) {
            throw EcoevolityError("calling wrong propose signature");
        }

        std::string target_parameter() const;

        std::string get_name() const;
};

class ComparisonMutationRateScaler : public TreeOperatorInterface<ScaleOperator> {

    public:
        ComparisonMutationRateScaler();
        ComparisonMutationRateScaler(double weight);
        ComparisonMutationRateScaler(double weight, double scale);
        // virtual ~ComparisonMutationRateScaler() { }

        void operate(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);
        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const;

        double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& parameter) const {
            throw EcoevolityError("calling wrong propose signature");
        }
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1) {
            throw EcoevolityError("calling wrong propose signature");
        }

        std::string target_parameter() const;

        std::string get_name() const;
};


class ChildPopulationSizeScaler : public TreeOperatorInterface<ScaleOperator> {

    public:
        ChildPopulationSizeScaler();
        ChildPopulationSizeScaler(double weight);
        ChildPopulationSizeScaler(double weight, double scale);
        // virtual ~ChildPopulationSizeScaler() { }

        void operate(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);
        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const;

        double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& parameter) const {
            throw EcoevolityError("calling wrong propose signature");
        }
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1) {
            throw EcoevolityError("calling wrong propose signature");
        }

        std::string target_parameter() const;

        std::string get_name() const;
};


class RootPopulationSizeScaler : public TreeOperatorInterface<ScaleOperator> {

    public:
        RootPopulationSizeScaler();
        RootPopulationSizeScaler(double weight);
        RootPopulationSizeScaler(double weight, double scale);
        // virtual ~RootPopulationSizeScaler() { }

        void operate(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);
        double propose(
                RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const;

        double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& parameter) const {
            throw EcoevolityError("calling wrong propose signature");
        }
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1) {
            throw EcoevolityError("calling wrong propose signature");
        }

        std::string target_parameter() const;

        std::string get_name() const;
};

class ComparisonHeightScaler : public TimeOperatorInterface<ScaleOperator> {

    public:
        ComparisonHeightScaler();
        ComparisonHeightScaler(double weight);
        ComparisonHeightScaler(double weight, double scale);
        // virtual ~ComparisonHeightScaler() { }

        void operate(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);
        double propose(
                RandomNumberGenerator& rng,
                PositiveRealParameter& node_height) const;

        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const {
            throw EcoevolityError("calling wrong propose signature");
        }
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1) {
            throw EcoevolityError("calling wrong propose signature");
        }

        std::string get_name() const;
        std::string target_parameter() const;
};

class ComparisonHeightMover : public TimeOperatorInterface<WindowOperator> {

    public:
        ComparisonHeightMover();
        ComparisonHeightMover(double weight);
        ComparisonHeightMover(double weight, double window_size);
        // virtual ~ComparisonHeightMover() { }

        void operate(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);
        double propose(
                RandomNumberGenerator& rng,
                PositiveRealParameter& node_height) const;

        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const {
            throw EcoevolityError("calling wrong propose signature");
        }
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1) {
            throw EcoevolityError("calling wrong propose signature");
        }

        std::string get_name() const;
        std::string target_parameter() const;
};


class UnivariateCollectionScaler : public CollectionOperatorInterface<ScaleOperator> {

    public:
        UnivariateCollectionScaler();
        UnivariateCollectionScaler(double weight);
        UnivariateCollectionScaler(double weight, double scale);
        // virtual ~UnivariateCollectionScaler() { }

        void operate(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& parameter) const {
            throw EcoevolityError("calling wrong propose signature");
        }
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const {
            throw EcoevolityError("calling wrong propose signature");
        }

        std::string target_parameter() const;

        std::string get_name() const;
};

class UnivariateCompositeCollectionScaler : public CollectionOperatorInterface<Operator> {

    protected:
        ComparisonHeightScaler height_scaler_ = ComparisonHeightScaler(0.0, 0.5);
        RootPopulationSizeScaler root_size_scaler_ = RootPopulationSizeScaler(0.0, 0.5);
        ChildPopulationSizeScaler child_size_scaler_ = ChildPopulationSizeScaler(0.0, 0.5);

    public:
        UnivariateCompositeCollectionScaler() : CollectionOperatorInterface<Operator>() { }
        UnivariateCompositeCollectionScaler(double weight) : CollectionOperatorInterface<Operator>(weight) { }
        UnivariateCompositeCollectionScaler(double weight, double scale);
        virtual ~UnivariateCompositeCollectionScaler() { }

        void operate(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);
        void perform_collection_move(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& parameter) const {
            throw EcoevolityError("calling wrong propose signature");
        }
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const {
            throw EcoevolityError("calling wrong propose signature");
        }

        std::string target_parameter() const;

        std::string get_name() const;

        std::string to_string(const OperatorSchedule& os) const;

        virtual void optimize(OperatorSchedule& os, double log_alpha) { }

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const { };

};


class CollectionScaler : public CollectionOperatorInterface<ScaleOperator> {

    friend class ReversibleJumpSampler;

    protected:
        UnivariateCollectionScaler uni_collection_scaler_ = UnivariateCollectionScaler(0.0, 0.5);

    public:
        CollectionScaler();
        CollectionScaler(double weight);
        CollectionScaler(double weight, double scale);
        // virtual ~CollectionScaler() { }

        void operate(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& parameter) const {
            throw EcoevolityError("calling wrong propose signature");
        }
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const {
            throw EcoevolityError("calling wrong propose signature");
        }

        std::string target_parameter() const;

        std::string get_name() const;

        std::string to_string(const OperatorSchedule& os) const;
};


class CompositeCollectionScaler : public CollectionScaler {

    friend class ReversibleJumpSampler;

    using CollectionScaler::propose;

    protected:
        UnivariateCompositeCollectionScaler uni_composite_collection_scaler_ = UnivariateCompositeCollectionScaler(0.0);

    public:
        CompositeCollectionScaler() : CollectionScaler() { }
        CompositeCollectionScaler(double weight) : CollectionScaler(weight) { }
        CompositeCollectionScaler(double weight, double scale) : CollectionScaler(weight, scale) { }
        // virtual ~CompositeCollectionScaler() { }

        void operate(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);

        std::string get_name() const;

        std::string to_string(const OperatorSchedule& os) const;
};


class DirichletProcessGibbsSampler : public CollectionOperatorInterface<Operator> {

    protected:
        unsigned int number_of_auxiliary_categories_ = 4;

    public:
        DirichletProcessGibbsSampler() : CollectionOperatorInterface<Operator>() { }
        DirichletProcessGibbsSampler(double weight) : CollectionOperatorInterface<Operator>(weight) { }
        DirichletProcessGibbsSampler(double weight,
                unsigned int number_of_auxiliary_categories);
        virtual ~DirichletProcessGibbsSampler() { }

        void set_number_of_auxiliary_categories(unsigned int n);
        unsigned int get_number_of_auxiliary_categories() const;

        std::string get_name() const;

        void perform_collection_move(
                RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads);

        void operate(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);
        double propose_gibbs(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& parameter) const {
            throw EcoevolityError("calling wrong propose signature");
        }
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const {
            throw EcoevolityError("calling wrong propose signature");
        }

        std::string target_parameter() const;

        virtual void optimize(OperatorSchedule& os, double log_alpha) { }

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const { };
};

class ReversibleJumpSampler : public CollectionOperatorInterface<Operator> {

    protected:
        CollectionScaler collection_scaler_ = CollectionScaler(0.0, 0.5);
        // ComparisonHeightScaler collection_scaler_ = ComparisonHeightScaler(0.0, 0.5);
        std::map<unsigned int, std::vector<double> > split_subset_size_probs_;
        std::map<unsigned int, double> ln_number_of_possible_splits_;
        void populate_split_subset_size_probabilities(
                unsigned int number_of_nodes_in_event);

    public:
        ReversibleJumpSampler() : CollectionOperatorInterface<Operator>() { }
        ReversibleJumpSampler(double weight) : CollectionOperatorInterface<Operator>(weight) { }
        virtual ~ReversibleJumpSampler() { }

        std::string get_name() const;

        std::string to_string(const OperatorSchedule& os) const;

        virtual void call_store_methods(
                ComparisonPopulationTreeCollection& comparisons) const;
        virtual void call_restore_methods(
                ComparisonPopulationTreeCollection& comparisons) const;

        void operate(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);

        OperatorInterface::OperatorTypeEnum get_type() const;

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons,
                unsigned int nthreads = 1);
        virtual double propose_jump_to_prior(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons);
        virtual double propose_jump_to_gap(RandomNumberGenerator& rng,
                ComparisonPopulationTreeCollection& comparisons);

        double propose(RandomNumberGenerator& rng,
                PositiveRealParameter& parameter) const {
            throw EcoevolityError("calling wrong propose signature");
        }
        double propose(RandomNumberGenerator& rng,
                ComparisonPopulationTree& tree) const {
            throw EcoevolityError("calling wrong propose signature");
        }

        const std::vector<double>& get_split_subset_size_probabilities(
                unsigned int number_of_nodes_in_event);

        void write_split_probabilities(std::ostream& out) const;

        std::string target_parameter() const;

        virtual void optimize(OperatorSchedule& os, double log_alpha) { }

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const { };
};


///////////////////////////////////////////////////////////////////////////////
// This 'works' but is *very* sensitive to window size and will often fail to
// sample the target (prior distribution)
//
// class ReversibleJumpWindowOperator : public ReversibleJumpSampler {
// 
//     protected:
//         ComparisonHeightMover height_mover_ = ComparisonHeightMover(0.0, 0.5);
// 
//     public:
//         ReversibleJumpWindowOperator() : ReversibleJumpSampler() { }
//         ReversibleJumpWindowOperator(double weight) : ReversibleJumpSampler(weight) { }
//         ReversibleJumpWindowOperator(double weight, double window_size);
//         virtual ~ReversibleJumpWindowOperator() { }
// 
//         void operate(RandomNumberGenerator& rng,
//                 ComparisonPopulationTreeCollection& comparisons,
//                 unsigned int nthreads = 1);
// 
//         /**
//          * @brief   Propose a new state.
//          *
//          * @return  Log of Hastings Ratio.
//          */
//         double propose(RandomNumberGenerator& rng,
//                 ComparisonPopulationTreeCollection& comparisons,
//                 unsigned int nthreads = 1);
// 
//         double propose(RandomNumberGenerator& rng,
//                 PositiveRealParameter& parameter) const {
//             throw EcoevolityError("calling wrong propose signature");
//         }
//         double propose(RandomNumberGenerator& rng,
//                 ComparisonPopulationTree& tree) const {
//             throw EcoevolityError("calling wrong propose signature");
//         }
// 
//         void update_height(RandomNumberGenerator& rng,
//                 double& height,
//                 double& hastings,
//                 double window_size) const;
// 
//         void propose_height_moves(RandomNumberGenerator& rng,
//                 ComparisonPopulationTreeCollection& comparisons,
//                 unsigned int nthreads = 1);
// 
//         OperatorInterface::OperatorTypeEnum get_type() const;
// 
//         std::string get_name() const;
// 
//         std::string to_string(const OperatorSchedule& os) const;
// 
//         std::string target_parameter() const;
// 
//         virtual void optimize(OperatorSchedule& os, double log_alpha) { }
// 
//         virtual void update(
//                 RandomNumberGenerator& rng,
//                 double& parameter_value,
//                 double& hastings_ratio) const { };
// };
///////////////////////////////////////////////////////////////////////////////

#endif
