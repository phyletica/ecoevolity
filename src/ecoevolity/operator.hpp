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
            multivariate_time_operator = 5,
        };

        virtual OperatorInterface::OperatorTypeEnum get_type() const = 0;

        double get_weight() const;

        void set_weight(double weight);

        virtual void call_store_methods(
                BaseComparisonPopulationTreeCollection * comparisons) const;
        virtual void call_restore_methods(
                BaseComparisonPopulationTreeCollection * comparisons) const;

        virtual std::string get_name() const = 0;

        virtual std::string target_parameter() const = 0;

        virtual void perform_collection_move(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads) = 0;

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads) = 0;

        virtual void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
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

        virtual int get_tree_index() const { return -1; }
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

    protected:
        int tree_index_ = -1;

    public:
        TimeOperatorInterface() : BaseOperatorInterface<DerivedOperatorType>() { }
        TimeOperatorInterface(unsigned int tree_index) :
            BaseOperatorInterface<DerivedOperatorType>(),
            tree_index_(tree_index)
        { }
        TimeOperatorInterface(double weight) : BaseOperatorInterface<DerivedOperatorType>(weight) { }
        TimeOperatorInterface(
                unsigned int tree_index,
                double weight) :
            BaseOperatorInterface<DerivedOperatorType>(weight),
            tree_index_(tree_index)
        { }

        virtual void perform_collection_move(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        virtual void perform_global_collection_move(
                RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        virtual void perform_tree_specific_collection_move(
                RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int height_index) = 0;

        virtual OperatorInterface::OperatorTypeEnum get_type() const;

        virtual void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1) = 0;

        virtual int get_tree_index() const;
};

template<class DerivedOperatorType>
class TreeOperatorInterface : public BaseOperatorInterface<DerivedOperatorType> {

    protected:
        int tree_index_ = -1;

    public:
        TreeOperatorInterface() : BaseOperatorInterface<DerivedOperatorType>() { }
        TreeOperatorInterface(unsigned int tree_index) :
            BaseOperatorInterface<DerivedOperatorType>(),
            tree_index_(tree_index)
        { }
        TreeOperatorInterface(double weight) : BaseOperatorInterface<DerivedOperatorType>(weight) { }
        TreeOperatorInterface(
                unsigned int tree_index,
                double weight) :
            BaseOperatorInterface<DerivedOperatorType>(weight),
            tree_index_(tree_index)
        { }

        virtual void perform_collection_move(
                RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);
        
        void perform_global_collection_move(
                RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        void perform_tree_specific_collection_move(
                RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index) = 0;

        virtual OperatorInterface::OperatorTypeEnum get_type() const;

        virtual void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1) = 0;

        virtual int get_tree_index() const;
};


template<class DerivedOperatorType>
class CollectionOperatorInterface : public BaseOperatorInterface<DerivedOperatorType> {
    public:
        CollectionOperatorInterface() : BaseOperatorInterface<DerivedOperatorType>() { }
        CollectionOperatorInterface(double weight) : BaseOperatorInterface<DerivedOperatorType>(weight) { }

        virtual void perform_collection_move(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        virtual double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads) = 0;

        virtual OperatorInterface::OperatorTypeEnum get_type() const;

        virtual void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
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

        virtual double get_move_amount(RandomNumberGenerator& rng) const { return 0.0; }

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

        double get_move_amount(RandomNumberGenerator& rng) const;

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

        double get_move_amount(RandomNumberGenerator& rng) const;

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
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        std::string target_parameter() const;

        std::string get_name() const;
};


class DiscountScaler : public CollectionOperatorInterface<ScaleOperator> {

    public:
        DiscountScaler();
        DiscountScaler(double weight);
        DiscountScaler(double weight, double scale);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        std::string target_parameter() const;

        std::string get_name() const;
};


class DiscountMover : public CollectionOperatorInterface<WindowOperator> {

    public:
        DiscountMover();
        DiscountMover(double weight);
        DiscountMover(double weight, double window_size);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        std::string target_parameter() const;

        std::string get_name() const;
};


class FreqMover : public TreeOperatorInterface<WindowOperator> {

    public:
        FreqMover();
        FreqMover(unsigned int tree_index);
        FreqMover(double weight);
        FreqMover(unsigned int tree_index, double weight);
        FreqMover(double weight, double window_size);
        FreqMover(unsigned int tree_index, double weight, double window_size);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(
                RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index);

        std::string target_parameter() const;

        std::string get_name() const;
};

class MutationRateScaler : public TreeOperatorInterface<ScaleOperator> {

    public:
        MutationRateScaler();
        MutationRateScaler(unsigned int tree_index);
        MutationRateScaler(double weight);
        MutationRateScaler(unsigned int tree_index, double weight);
        MutationRateScaler(double weight, double scale);
        MutationRateScaler(unsigned int tree_index, double weight, double scale);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(
                RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index);

        std::string target_parameter() const;

        std::string get_name() const;
};

/**
 * Propose new population size multipliers by sampling from a Dirichlet target
 * distribution.
 *
 * @note    This move is adapted from the 'DirichletMove' class in:
 *              Phycas
 *              <http://www.phycas.org/>
 *              <https://github.com/plewis/phycas>
 *              License:    Gnu GPL Version 2
 *              Authors:    Paul Lewis, Mark Holder, and David Swofford
 */
class RelativePopulationSizeMixer : public TreeOperatorInterface<ScaleOperator> {

    public:
        RelativePopulationSizeMixer();
        RelativePopulationSizeMixer(unsigned int tree_index);
        RelativePopulationSizeMixer(double weight);
        RelativePopulationSizeMixer(unsigned int tree_index, double weight);
        RelativePopulationSizeMixer(double weight, double scale);
        RelativePopulationSizeMixer(unsigned int tree_index, double weight, double scale);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(
                RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index);

        std::string target_parameter() const;

        std::string get_name() const;
};

class MeanPopulationSizeScaler : public TreeOperatorInterface<ScaleOperator> {

    public:
        MeanPopulationSizeScaler();
        MeanPopulationSizeScaler(unsigned int tree_index);
        MeanPopulationSizeScaler(double weight);
        MeanPopulationSizeScaler(unsigned int tree_index, double weight);
        MeanPopulationSizeScaler(double weight, double scale);
        MeanPopulationSizeScaler(unsigned int tree_index, double weight, double scale);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(
                RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index);

        std::string target_parameter() const;

        std::string get_name() const;
};

class LeafPopulationSizeScaler : public TreeOperatorInterface<ScaleOperator> {

    public:
        LeafPopulationSizeScaler();
        LeafPopulationSizeScaler(unsigned int tree_index);
        LeafPopulationSizeScaler(double weight);
        LeafPopulationSizeScaler(unsigned int tree_index, double weight);
        LeafPopulationSizeScaler(double weight, double scale);
        LeafPopulationSizeScaler(unsigned int tree_index, double weight, double scale);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(
                RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index);

        std::string target_parameter() const;

        std::string get_name() const;
};


class RootPopulationSizeScaler : public TreeOperatorInterface<ScaleOperator> {

    public:
        RootPopulationSizeScaler();
        RootPopulationSizeScaler(unsigned int tree_index);
        RootPopulationSizeScaler(double weight);
        RootPopulationSizeScaler(unsigned int tree_index, double weight);
        RootPopulationSizeScaler(double weight, double scale);
        RootPopulationSizeScaler(unsigned int tree_index, double weight, double scale);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(
                RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index);

        std::string target_parameter() const;

        std::string get_name() const;
};


class RootRelativePopulationSizeMover : public TreeOperatorInterface<WindowOperator> {

    public:
        RootRelativePopulationSizeMover();
        RootRelativePopulationSizeMover(unsigned int tree_index);
        RootRelativePopulationSizeMover(double weight);
        RootRelativePopulationSizeMover(unsigned int tree_index, double weight);
        RootRelativePopulationSizeMover(double weight, double scale);
        RootRelativePopulationSizeMover(unsigned int tree_index, double weight, double scale);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(
                RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index);

        std::string target_parameter() const;

        std::string get_name() const;
};

class LeafRelativePopulationSizeMover : public TreeOperatorInterface<WindowOperator> {

    public:
        LeafRelativePopulationSizeMover();
        LeafRelativePopulationSizeMover(unsigned int tree_index);
        LeafRelativePopulationSizeMover(double weight);
        LeafRelativePopulationSizeMover(unsigned int tree_index, double weight);
        LeafRelativePopulationSizeMover(double weight, double scale);
        LeafRelativePopulationSizeMover(unsigned int tree_index, double weight, double scale);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(
                RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index);

        std::string target_parameter() const;

        std::string get_name() const;
};


class EventTimeScaler : public TimeOperatorInterface<ScaleOperator> {

    public:
        EventTimeScaler();
        EventTimeScaler(unsigned int tree_index);
        EventTimeScaler(double weight);
        EventTimeScaler(unsigned int tree_index, double weight);
        EventTimeScaler(double weight, double scale);
        EventTimeScaler(unsigned int tree_index, double weight, double scale);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int index);

        double propose_by_height(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int height_index);

        double propose_by_tree(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index);

        std::string get_name() const;
        std::string target_parameter() const;
};

class EventTimeMover : public TimeOperatorInterface<WindowOperator> {

    public:
        EventTimeMover();
        EventTimeMover(unsigned int tree_index);
        EventTimeMover(double weight);
        EventTimeMover(unsigned int tree_index, double weight);
        EventTimeMover(double weight, double window_size);
        EventTimeMover(unsigned int tree_index, double weight, double window_size);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int index);

        double propose_by_height(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int height_index);

        double propose_by_tree(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index);

        std::string get_name() const;
        std::string target_parameter() const;
};


class UnivariateTimeSizeScaler : public CollectionOperatorInterface<ScaleOperator> {

    public:
        UnivariateTimeSizeScaler();
        UnivariateTimeSizeScaler(double weight);
        UnivariateTimeSizeScaler(double weight, double scale);
        // virtual ~UnivariateTimeSizeScaler() { }

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        std::string target_parameter() const;

        std::string get_name() const;
};

class UnivariateTimeSizeRateScaler : public CollectionOperatorInterface<ScaleOperator> {

    public:
        UnivariateTimeSizeRateScaler();
        UnivariateTimeSizeRateScaler(double weight);
        UnivariateTimeSizeRateScaler(double weight, double scale);
        // virtual ~UnivariateTimeSizeRateScaler() { }

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);


        std::string target_parameter() const;

        std::string get_name() const;
};

class UnivariateTimeMeanSizeRateScaler : public CollectionOperatorInterface<ScaleOperator> {

    public:
        UnivariateTimeMeanSizeRateScaler();
        UnivariateTimeMeanSizeRateScaler(double weight);
        UnivariateTimeMeanSizeRateScaler(double weight, double scale);
        // virtual ~UnivariateTimeSizeRateScaler() { }

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);


        std::string target_parameter() const;

        std::string get_name() const;
};

class UnivariateCompositeTimeSizeScaler : public CollectionOperatorInterface<Operator> {

    protected:
        EventTimeScaler height_scaler_ = EventTimeScaler(0.0, 0.5);
        RootPopulationSizeScaler root_size_scaler_ = RootPopulationSizeScaler(0.0, 0.5);
        LeafPopulationSizeScaler child_size_scaler_ = LeafPopulationSizeScaler(0.0, 0.5);

    public:
        UnivariateCompositeTimeSizeScaler() : CollectionOperatorInterface<Operator>() { }
        UnivariateCompositeTimeSizeScaler(double weight) : CollectionOperatorInterface<Operator>(weight) { }
        UnivariateCompositeTimeSizeScaler(double weight, double scale);
        virtual ~UnivariateCompositeTimeSizeScaler() { }

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);
        void perform_collection_move(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);
        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        void scale_heights(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);
        void scale_root_population_sizes(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);
        void scale_child_population_sizes(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        std::string target_parameter() const;

        std::string get_name() const;

        std::string to_string(const OperatorSchedule& os) const;

        virtual void optimize(OperatorSchedule& os, double log_alpha) { }

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const { };

};

class UnivariateCompositeTimeSizeRateScaler : public CollectionOperatorInterface<Operator> {

    protected:
        EventTimeScaler height_scaler_ = EventTimeScaler(0.0, 0.5);
        RootPopulationSizeScaler root_size_scaler_ = RootPopulationSizeScaler(0.0, 0.5);
        LeafPopulationSizeScaler child_size_scaler_ = LeafPopulationSizeScaler(0.0, 0.5);
        MutationRateScaler mutation_rate_scaler_ = MutationRateScaler(0.0, 0.5);

    public:
        UnivariateCompositeTimeSizeRateScaler() : CollectionOperatorInterface<Operator>() { }
        UnivariateCompositeTimeSizeRateScaler(double weight) : CollectionOperatorInterface<Operator>(weight) { }
        UnivariateCompositeTimeSizeRateScaler(double weight, double scale);
        virtual ~UnivariateCompositeTimeSizeRateScaler() { }

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);
        void perform_collection_move(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);
        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        void scale_heights(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);
        void scale_root_population_sizes(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);
        void scale_child_population_sizes(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);
        void scale_mutation_rates(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        std::string target_parameter() const;

        std::string get_name() const;

        std::string to_string(const OperatorSchedule& os) const;

        virtual void optimize(OperatorSchedule& os, double log_alpha) { }

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const { };

};

class UnivariateCompositeTimeMeanSizeRateScaler : public CollectionOperatorInterface<Operator> {

    protected:
        EventTimeScaler height_scaler_ = EventTimeScaler(0.0, 0.5);
        MeanPopulationSizeScaler size_scaler_ = MeanPopulationSizeScaler(0.0, 0.5);
        RelativePopulationSizeMixer size_multiplier_mixer_ = RelativePopulationSizeMixer(0.0, 0.05);
        MutationRateScaler mutation_rate_scaler_ = MutationRateScaler(0.0, 0.5);

    public:
        UnivariateCompositeTimeMeanSizeRateScaler() : CollectionOperatorInterface<Operator>() { }
        UnivariateCompositeTimeMeanSizeRateScaler(double weight) : CollectionOperatorInterface<Operator>(weight) { }
        UnivariateCompositeTimeMeanSizeRateScaler(double weight, double scale);
        virtual ~UnivariateCompositeTimeMeanSizeRateScaler() { }

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);
        void perform_collection_move(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);
        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        void scale_heights(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);
        void scale_population_sizes(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);
        void mix_population_size_multipliers(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);
        void scale_mutation_rates(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        std::string target_parameter() const;

        std::string get_name() const;

        std::string to_string(const OperatorSchedule& os) const;

        virtual void optimize(OperatorSchedule& os, double log_alpha) { }

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const { };
};


class TimeSizeMixer : public TimeOperatorInterface<ScaleOperator> {

    public:
        TimeSizeMixer();
        TimeSizeMixer(unsigned int tree_index);
        TimeSizeMixer(double weight);
        TimeSizeMixer(unsigned int tree_index, double weight);
        TimeSizeMixer(double weight, double scale);
        TimeSizeMixer(unsigned int tree_index, double weight, double scale);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int index);

        double propose_by_height(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int height_index);

        double propose_by_tree(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index);

        std::string get_name() const;

        std::string target_parameter() const;

        std::string to_string(const OperatorSchedule& os) const;

        OperatorInterface::OperatorTypeEnum get_type() const {
            return OperatorInterface::OperatorTypeEnum::multivariate_time_operator;
        }
};


class TimeRootSizeMixer : public TimeOperatorInterface<ScaleOperator> {

    public:
        TimeRootSizeMixer();
        TimeRootSizeMixer(unsigned int tree_index);
        TimeRootSizeMixer(double weight);
        TimeRootSizeMixer(unsigned int tree_index, double weight);
        TimeRootSizeMixer(double weight, double scale);
        TimeRootSizeMixer(unsigned int tree_index, double weight, double scale);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int index);

        double propose_by_height(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int height_index);

        double propose_by_tree(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index);

        std::string get_name() const;

        std::string target_parameter() const;

        std::string to_string(const OperatorSchedule& os) const;

        OperatorInterface::OperatorTypeEnum get_type() const {
            return OperatorInterface::OperatorTypeEnum::multivariate_time_operator;
        }
};


class TimeSizeScaler : public TimeOperatorInterface<ScaleOperator> {

    public:
        TimeSizeScaler();
        TimeSizeScaler(unsigned int tree_index);
        TimeSizeScaler(double weight);
        TimeSizeScaler(unsigned int tree_index, double weight);
        TimeSizeScaler(double weight, double scale);
        TimeSizeScaler(unsigned int tree_index, double weight, double scale);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int index);

        double propose_by_height(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int height_index);

        double propose_by_tree(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index);

        std::string get_name() const;

        std::string target_parameter() const;

        std::string to_string(const OperatorSchedule& os) const;

        OperatorInterface::OperatorTypeEnum get_type() const {
            return OperatorInterface::OperatorTypeEnum::multivariate_time_operator;
        }
};


class TimeSizeRateMixer : public TimeOperatorInterface<ScaleOperator> {

    public:
        TimeSizeRateMixer();
        TimeSizeRateMixer(unsigned int tree_index);
        TimeSizeRateMixer(double weight);
        TimeSizeRateMixer(unsigned int tree_index, double weight);
        TimeSizeRateMixer(double weight, double scale);
        TimeSizeRateMixer(unsigned int tree_index, double weight, double scale);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int index);

        double propose_by_height(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int height_index);

        double propose_by_tree(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index);

        std::string get_name() const;

        std::string target_parameter() const;

        std::string to_string(const OperatorSchedule& os) const;

        OperatorInterface::OperatorTypeEnum get_type() const {
            return OperatorInterface::OperatorTypeEnum::multivariate_time_operator;
        }
};

class TimeMeanSizeRateMixer : public TimeOperatorInterface<ScaleOperator> {

    public:
        TimeMeanSizeRateMixer();
        TimeMeanSizeRateMixer(unsigned int tree_index);
        TimeMeanSizeRateMixer(double weight);
        TimeMeanSizeRateMixer(unsigned int tree_index, double weight);
        TimeMeanSizeRateMixer(double weight, double scale);
        TimeMeanSizeRateMixer(unsigned int tree_index, double weight, double scale);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int index);

        double propose_by_height(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int height_index);

        double propose_by_tree(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index);

        std::string get_name() const;

        std::string target_parameter() const;

        std::string to_string(const OperatorSchedule& os) const;

        OperatorInterface::OperatorTypeEnum get_type() const {
            return OperatorInterface::OperatorTypeEnum::multivariate_time_operator;
        }
};


class TimeSizeRateScaler : public TimeOperatorInterface<ScaleOperator> {

    public:
        TimeSizeRateScaler();
        TimeSizeRateScaler(unsigned int tree_index);
        TimeSizeRateScaler(double weight);
        TimeSizeRateScaler(unsigned int tree_index, double weight);
        TimeSizeRateScaler(double weight, double scale);
        TimeSizeRateScaler(unsigned int tree_index, double weight, double scale);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int index);

        double propose_by_height(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int height_index);

        double propose_by_tree(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index);

        std::string get_name() const;

        std::string target_parameter() const;

        std::string to_string(const OperatorSchedule& os) const;

        OperatorInterface::OperatorTypeEnum get_type() const {
            return OperatorInterface::OperatorTypeEnum::multivariate_time_operator;
        }
};


class TimeMeanSizeRateScaler : public TimeOperatorInterface<ScaleOperator> {

    public:
        TimeMeanSizeRateScaler();
        TimeMeanSizeRateScaler(unsigned int tree_index);
        TimeMeanSizeRateScaler(double weight);
        TimeMeanSizeRateScaler(unsigned int tree_index, double weight);
        TimeMeanSizeRateScaler(double weight, double scale);
        TimeMeanSizeRateScaler(unsigned int tree_index, double weight, double scale);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int index);

        double propose_by_height(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int height_index);

        double propose_by_tree(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int tree_index);

        std::string get_name() const;

        std::string target_parameter() const;

        std::string to_string(const OperatorSchedule& os) const;

        OperatorInterface::OperatorTypeEnum get_type() const {
            return OperatorInterface::OperatorTypeEnum::multivariate_time_operator;
        }
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

        std::string to_string(const OperatorSchedule& os) const;

        void perform_collection_move(
                RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        double propose_gibbs(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        std::string target_parameter() const;

        virtual void optimize(OperatorSchedule& os, double log_alpha) { }

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const { };
};


class PitmanYorProcessGibbsSampler : public CollectionOperatorInterface<Operator> {

    protected:
        unsigned int number_of_auxiliary_categories_ = 4;

    public:
        PitmanYorProcessGibbsSampler() : CollectionOperatorInterface<Operator>() { }
        PitmanYorProcessGibbsSampler(double weight) : CollectionOperatorInterface<Operator>(weight) { }
        PitmanYorProcessGibbsSampler(double weight,
                unsigned int number_of_auxiliary_categories);
        virtual ~PitmanYorProcessGibbsSampler() { }

        void set_number_of_auxiliary_categories(unsigned int n);
        unsigned int get_number_of_auxiliary_categories() const;

        std::string get_name() const;

        std::string to_string(const OperatorSchedule& os) const;

        void perform_collection_move(
                RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        double propose_gibbs(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);

        std::string target_parameter() const;

        virtual void optimize(OperatorSchedule& os, double log_alpha) { }

        virtual void update(
                RandomNumberGenerator& rng,
                double& parameter_value,
                double& hastings_ratio) const { };
};


class ReversibleJumpSampler : public CollectionOperatorInterface<Operator> {

    protected:
        EventTimeScaler time_scaler_ = EventTimeScaler(0.0, 0.5);
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
                BaseComparisonPopulationTreeCollection * comparisons) const;
        virtual void call_restore_methods(
                BaseComparisonPopulationTreeCollection * comparisons) const;

        void perform_collection_move(
                RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        void operate(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads = 1);

        OperatorInterface::OperatorTypeEnum get_type() const;

        /**
         * @brief   Propose a new state.
         *
         * @return  Log of Hastings Ratio.
         */
        double propose(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons,
                unsigned int nthreads);
        virtual double propose_jump_to_prior(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons);
        virtual double propose_jump_to_gap(RandomNumberGenerator& rng,
                BaseComparisonPopulationTreeCollection * comparisons);

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
//         EventTimeMover height_mover_ = EventTimeMover(0.0, 0.5);
// 
//     public:
//         ReversibleJumpWindowOperator() : ReversibleJumpSampler() { }
//         ReversibleJumpWindowOperator(double weight) : ReversibleJumpSampler(weight) { }
//         ReversibleJumpWindowOperator(double weight, double window_size);
//         virtual ~ReversibleJumpWindowOperator() { }
// 
//         void operate(RandomNumberGenerator& rng,
//                 BaseComparisonPopulationTreeCollection * comparisons,
//                 unsigned int nthreads = 1);
// 
//         /**
//          * @brief   Propose a new state.
//          *
//          * @return  Log of Hastings Ratio.
//          */
//         double propose(RandomNumberGenerator& rng,
//                 BaseComparisonPopulationTreeCollection * comparisons,
//                 unsigned int nthreads);
// 
//         void update_height(RandomNumberGenerator& rng,
//                 double& height,
//                 double& hastings,
//                 double window_size) const;
// 
//         void propose_height_moves(RandomNumberGenerator& rng,
//                 BaseComparisonPopulationTreeCollection * comparisons,
//                 unsigned int nthreads);
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
