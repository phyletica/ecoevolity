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

#include "operator.hpp"


//////////////////////////////////////////////////////////////////////////////
// OperatorInterface methods
//////////////////////////////////////////////////////////////////////////////

OperatorInterface::OperatorInterface(double weight) {
    this->set_weight(weight);
}

void OperatorInterface::set_weight(double weight) {
    ECOEVOLITY_ASSERT(weight >= 0.0);
    this->weight_ = weight;
}

double OperatorInterface::get_weight() const {
    return this->weight_;
}

void OperatorInterface::call_store_methods(
        BaseComparisonPopulationTreeCollection * comparisons) const {
    comparisons->store_state();
}

void OperatorInterface::call_restore_methods(
        BaseComparisonPopulationTreeCollection * comparisons) const {
    comparisons->restore_state();
}


//////////////////////////////////////////////////////////////////////////////
// BaseOperatorInterface methods
//////////////////////////////////////////////////////////////////////////////

template<class DerivedOperatorType>
void BaseOperatorInterface<DerivedOperatorType>::optimize(
        OperatorSchedule& os,
        double log_alpha) {
    this->op_.optimize(os, log_alpha);
}

template<class DerivedOperatorType>
void BaseOperatorInterface<DerivedOperatorType>::accept(
        const OperatorSchedule& os) {
    this->op_.accept(os);
}

template<class DerivedOperatorType>
void BaseOperatorInterface<DerivedOperatorType>::reject(
        const OperatorSchedule& os) {
    this->op_.reject(os);
}

template<class DerivedOperatorType>
std::string BaseOperatorInterface<DerivedOperatorType>::header_string() const {
    return "name\tnumber_accepted\tnumber_rejected\tweight\tweight_prob\ttuning_parameter\n";
}

template<class DerivedOperatorType>
std::string BaseOperatorInterface<DerivedOperatorType>::to_string(
        const OperatorSchedule& os) const {
    std::ostringstream ss;
    ss << this->get_name() << "\t" 
       << this->op_.get_number_accepted() << "\t"
       << this->op_.get_number_rejected() << "\t"
       << this->get_weight() << "\t";

    if (os.get_total_weight() > 0.0) {
        ss << this->get_weight() / os.get_total_weight() << "\t";
    }
    else {
        ss << "nan\t";
    }

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

template<class DerivedOperatorType>
double BaseOperatorInterface<DerivedOperatorType>::get_coercable_parameter_value() const {
    return this->op_.get_coercable_parameter_value();
}

template<class DerivedOperatorType>
void BaseOperatorInterface<DerivedOperatorType>::set_coercable_parameter_value(
        double value) {
    this->op_.set_coercable_parameter_value(value);
}

template<class DerivedOperatorType>
void BaseOperatorInterface<DerivedOperatorType>::update(
        RandomNumberGenerator& rng,
        double& parameter_value,
        double& hastings_ratio) const {
    this->op_.update(rng, parameter_value, hastings_ratio);
}

template<class DerivedOperatorType>
double BaseOperatorInterface<DerivedOperatorType>::calc_delta(
        OperatorSchedule& os,
        double log_alpha) const {
    return this->op_.calc_delta(os, log_alpha);
}

template<class DerivedOperatorType>
double BaseOperatorInterface<DerivedOperatorType>::get_target_acceptance_probability() const {
    return this->op_.get_target_acceptance_probability();
}

template<class DerivedOperatorType>
unsigned int BaseOperatorInterface<DerivedOperatorType>::get_number_rejected() const {
    return this->op_.get_number_rejected();
}

template<class DerivedOperatorType>
unsigned int BaseOperatorInterface<DerivedOperatorType>::get_number_accepted() const {
    return this->op_.get_number_accepted();
}

template<class DerivedOperatorType>
unsigned int BaseOperatorInterface<DerivedOperatorType>::get_number_rejected_for_correction() const {
    return this->op_.get_number_rejected_for_correction();
}

template<class DerivedOperatorType>
unsigned int BaseOperatorInterface<DerivedOperatorType>::get_number_accepted_for_correction() const {
    return this->op_.get_number_accepted_for_correction();
}


//////////////////////////////////////////////////////////////////////////////
// Operator methods
//////////////////////////////////////////////////////////////////////////////

double Operator::get_coercable_parameter_value() const {
    return std::numeric_limits<double>::quiet_NaN();
}

double Operator::calc_delta(OperatorSchedule& os,
        double log_alpha) const {
    return os.calc_delta(*this, log_alpha);
}

void Operator::accept(const OperatorSchedule& os) {
    ++this->number_accepted_;
    if (os.get_auto_optimize_delay_count() >= os.get_auto_optimize_delay()) {
        ++this->number_accepted_for_correction_;
    }
}

void Operator::reject(const OperatorSchedule& os) {
    ++this->number_rejected_;
    if (os.get_auto_optimize_delay_count() >= os.get_auto_optimize_delay()) {
        ++this->number_rejected_for_correction_;
    }
}

double Operator::get_target_acceptance_probability() const {
    // Some prelim tests confirm that an acceptance rate of 0.44 leads to
    // better mixing for the simple, univariate random variables that the
    // operators are updating.
    // return 0.234;
    return 0.44;
}

unsigned int Operator::get_number_rejected() const {
    return this->number_rejected_;
}
unsigned int Operator::get_number_accepted() const {
    return this->number_accepted_;
}
unsigned int Operator::get_number_rejected_for_correction() const {
    return this->number_rejected_for_correction_;
}
unsigned int Operator::get_number_accepted_for_correction() const {
    return this->number_accepted_for_correction_;
}


//////////////////////////////////////////////////////////////////////////////
// ScaleOperator methods
//////////////////////////////////////////////////////////////////////////////

ScaleOperator::ScaleOperator(double scale) : Operator() {
    this->set_scale(scale);
}

void ScaleOperator::set_scale(double scale) {
    ECOEVOLITY_ASSERT(scale > 0.0);
    this->scale_ = scale;
}

double ScaleOperator::get_scale() const {
    return this->scale_;
}

double ScaleOperator::get_move_amount(RandomNumberGenerator& rng) const {
    return std::exp(this->scale_ * ((2.0 * rng.uniform_real()) - 1.0));
}

void ScaleOperator::update(
        RandomNumberGenerator& rng,
        double& parameter_value,
        double& hastings_ratio) const {
    double multiplier = this->get_move_amount(rng);
    parameter_value *= multiplier;
    hastings_ratio = std::log(multiplier);
}

void ScaleOperator::optimize(OperatorSchedule& os, double log_alpha) {
    double delta = this->calc_delta(os, log_alpha);
    delta += std::log(this->scale_);
    this->set_scale(std::exp(delta));
}

double ScaleOperator::get_coercable_parameter_value() const {
    return this->scale_;
}

void ScaleOperator::set_coercable_parameter_value(double value) {
    this->set_scale(value);
}


//////////////////////////////////////////////////////////////////////////////
// WindowOperator methods
//////////////////////////////////////////////////////////////////////////////

WindowOperator::WindowOperator(double window_size) : Operator() {
    this->set_window_size(window_size);
}
void WindowOperator::set_window_size(double window_size) {
    ECOEVOLITY_ASSERT(window_size > 0.0);
    this->window_size_ = window_size;
}
double WindowOperator::get_window_size() const {
    return this->window_size_;
}

double WindowOperator::get_move_amount(RandomNumberGenerator& rng) const {
    return (rng.uniform_real() * 2 * this->window_size_) - this->window_size_;
}

void WindowOperator::update(
        RandomNumberGenerator& rng,
        double& parameter_value,
        double& hastings_ratio) const {
    double addend = this->get_move_amount(rng);
    parameter_value += addend;
    hastings_ratio = 0.0;
}

void WindowOperator::optimize(OperatorSchedule& os, double log_alpha) {
    double delta = this->calc_delta(os, log_alpha);
    delta += std::log(this->window_size_);
    this->set_window_size(std::exp(delta));
}

double WindowOperator::get_coercable_parameter_value() const {
    return this->window_size_;
}

void WindowOperator::set_coercable_parameter_value(double value) {
    this->set_window_size(value);
}


//////////////////////////////////////////////////////////////////////////////
// TimeOperatorInterface methods
//////////////////////////////////////////////////////////////////////////////

template<class DerivedOperatorType>
OperatorInterface::OperatorTypeEnum TimeOperatorInterface<DerivedOperatorType>::get_type() const {
    return OperatorInterface::OperatorTypeEnum::time_operator;
}

template<class DerivedOperatorType>
void TimeOperatorInterface<DerivedOperatorType>::perform_collection_move(
        RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    if (this->tree_index_ < 0) {
        this->perform_global_collection_move(rng, comparisons, nthreads);
        return;
    }
    this->perform_tree_specific_collection_move(rng, comparisons, nthreads);
}

template<class DerivedOperatorType>
void TimeOperatorInterface<DerivedOperatorType>::perform_global_collection_move(
        RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->call_store_methods(comparisons);

    std::vector<double> hastings_ratios;
    hastings_ratios.reserve(comparisons->get_number_of_events());
    for (unsigned int height_idx = 0; height_idx < comparisons->get_number_of_events(); ++height_idx) {
        hastings_ratios.push_back(this->propose(rng, comparisons, height_idx));
    }
    comparisons->make_trees_dirty();
    comparisons->compute_tree_partials();
    for (unsigned int height_idx = 0; height_idx < comparisons->get_number_of_events(); ++height_idx) {
        double old_lnl = 0.0;
        double new_lnl = 0.0;
        double old_prior_ln_pdf = 0.0;
        double new_prior_ln_pdf = 0.0;
        for (unsigned int tree_idx = 0; tree_idx < comparisons->get_number_of_trees(); ++tree_idx) {
            if (comparisons->get_height_index(tree_idx) == height_idx) {
                std::shared_ptr<PopulationTree> tree = comparisons->get_tree(tree_idx);
                old_lnl += tree->get_stored_log_likelihood_value();
                new_lnl += tree->get_log_likelihood_value();
                old_prior_ln_pdf += tree->get_stored_log_prior_density_value();
                new_prior_ln_pdf += tree->get_log_prior_density_value();
            }
        }
        old_prior_ln_pdf += comparisons->get_height_parameter(height_idx)->relative_prior_ln_pdf(
                comparisons->get_height_parameter(height_idx)->get_stored_value());
        new_prior_ln_pdf += comparisons->get_height_parameter(height_idx)->relative_prior_ln_pdf();
        double likelihood_ratio = new_lnl - old_lnl;
        double prior_ratio = new_prior_ln_pdf - old_prior_ln_pdf;
        double hastings_ratio = hastings_ratios.at(height_idx);
        double acceptance_probability =
                likelihood_ratio + 
                prior_ratio +
                hastings_ratio;
        double u = rng.uniform_real();
        if (u < std::exp(acceptance_probability)) {
            this->accept(comparisons->get_operator_schedule());
        }
        else {
            this->reject(comparisons->get_operator_schedule());
            comparisons->get_height_parameter(height_idx)->restore();
            for (unsigned int tree_idx = 0; tree_idx < comparisons->get_number_of_trees(); ++tree_idx) {
                if (comparisons->get_height_index(tree_idx) == height_idx) {
                    comparisons->get_tree(tree_idx)->restore_state();
                }
            }
        }
        this->optimize(comparisons->get_operator_schedule(), acceptance_probability);
    }
    comparisons->make_trees_clean();
    comparisons->compute_log_likelihood_and_prior(false);
}

template<class DerivedOperatorType>
void TimeOperatorInterface<DerivedOperatorType>::perform_tree_specific_collection_move(
        RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->call_store_methods(comparisons);

    double hastings_ratio = this->propose(rng, comparisons, this->tree_index_);
    comparisons->compute_log_likelihood_and_prior(true);

    double likelihood_ratio = 
        comparisons->get_log_likelihood() -
        comparisons->get_stored_log_likelihood();
    double prior_ratio = 
        comparisons->get_log_prior_density() -
        comparisons->get_stored_log_prior_density();
    double acceptance_probability =
            likelihood_ratio + 
            prior_ratio +
            hastings_ratio;
    double u = rng.uniform_real();
    if (u < std::exp(acceptance_probability)) {
        this->accept(comparisons->get_operator_schedule());
    }
    else {
        this->reject(comparisons->get_operator_schedule());
        this->call_restore_methods(comparisons);
    }
    comparisons->make_trees_clean();
    this->optimize(comparisons->get_operator_schedule(), acceptance_probability);
}

template<class DerivedOperatorType>
int TimeOperatorInterface<DerivedOperatorType>::get_tree_index() const {
    return this->tree_index_;
}


//////////////////////////////////////////////////////////////////////////////
// TreeOperatorInterface methods
//////////////////////////////////////////////////////////////////////////////

template<class DerivedOperatorType>
OperatorInterface::OperatorTypeEnum TreeOperatorInterface<DerivedOperatorType>::get_type() const {
    return OperatorInterface::OperatorTypeEnum::tree_operator;
}

template<class DerivedOperatorType>
void TreeOperatorInterface<DerivedOperatorType>::perform_collection_move(
        RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    if (this->tree_index_ < 0) {
        this->perform_global_collection_move(rng, comparisons, nthreads);
        return;
    }
    this->perform_tree_specific_collection_move(rng, comparisons, nthreads);
}

template<class DerivedOperatorType>
void TreeOperatorInterface<DerivedOperatorType>::perform_global_collection_move(
        RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->call_store_methods(comparisons);

    std::vector<double> hastings_ratios;
    hastings_ratios.reserve(comparisons->get_number_of_trees());
    for (unsigned int tree_idx = 0; tree_idx < comparisons->get_number_of_trees(); ++tree_idx) {
        hastings_ratios.push_back(this->propose(rng, comparisons, tree_idx));
    }
    comparisons->compute_tree_partials();
    for (unsigned int tree_idx = 0; tree_idx < comparisons->get_number_of_trees(); ++tree_idx) {
        std::shared_ptr<PopulationTree> tree = comparisons->get_tree(tree_idx);
        // Check to see if we updated a fixed parameter. If so, do
        // nothing and continue to next tree (to avoid counting toward
        // operator acceptance ratio).
        if (((this->target_parameter() == "population size") ||
                (this->target_parameter() == "population size multipliers")) &&
                        (tree->population_sizes_are_fixed())) {
            ECOEVOLITY_ASSERT(! tree->is_dirty());
            continue;
        }
        if (tree->using_population_size_multipliers()) {
            if ((this->target_parameter() == "population size") &&
                    (tree->mean_population_size_is_fixed())) {
                ECOEVOLITY_ASSERT(! tree->is_dirty());
                continue;
            }
            if ((this->target_parameter() == "population size multipliers") &&
                    (tree->population_size_multipliers_are_fixed())) {
                ECOEVOLITY_ASSERT(! tree->is_dirty());
                continue;
            }
        }
        if ((this->target_parameter() == "freq 1") &&
                (tree->state_frequencies_are_fixed())) {
            ECOEVOLITY_ASSERT(! tree->is_dirty());
            continue;
        }
        if ((this->target_parameter() == "mutation rate") &&
                (tree->mutation_rate_is_fixed())) {
            ECOEVOLITY_ASSERT(! tree->is_dirty());
            continue;
        }
        double likelihood_ratio =
                tree->get_log_likelihood_value() -
                tree->get_stored_log_likelihood_value();
        double prior_ratio =
                tree->get_log_prior_density_value() -
                tree->get_stored_log_prior_density_value();
        double hastings_ratio = hastings_ratios.at(tree_idx);
        double acceptance_probability =
                likelihood_ratio + 
                prior_ratio +
                hastings_ratio;
        double u = rng.uniform_real();
        if (u < std::exp(acceptance_probability)) {
            this->accept(comparisons->get_operator_schedule());
        }
        else {
            this->reject(comparisons->get_operator_schedule());
            tree->restore_state();
        }
        tree->make_clean();
        this->optimize(comparisons->get_operator_schedule(), acceptance_probability);
    }
    // Make sure the likelihood and prior of comparisons is up to date
    comparisons->compute_log_likelihood_and_prior(false);
}

template<class DerivedOperatorType>
void TreeOperatorInterface<DerivedOperatorType>::perform_tree_specific_collection_move(
        RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    std::shared_ptr<PopulationTree> tree = comparisons->get_tree(this->tree_index_);
    tree->store_state();

    double hastings_ratio = this->propose(rng, comparisons, this->tree_index_);

    tree->compute_log_likelihood_and_prior(nthreads);

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
        this->accept(comparisons->get_operator_schedule());
    }
    else {
        this->reject(comparisons->get_operator_schedule());
        tree->restore_state();
    }
    tree->make_clean();
    this->optimize(comparisons->get_operator_schedule(), acceptance_probability);
    // Make sure the likelihood and prior of comparisons is up to date
    comparisons->compute_log_likelihood_and_prior(false);
}

template<class DerivedOperatorType>
int TreeOperatorInterface<DerivedOperatorType>::get_tree_index() const {
    return this->tree_index_;
}


//////////////////////////////////////////////////////////////////////////////
// CollectionOperatorInterface methods
//////////////////////////////////////////////////////////////////////////////

template<class DerivedOperatorType>
OperatorInterface::OperatorTypeEnum CollectionOperatorInterface<DerivedOperatorType>::get_type() const {
    return OperatorInterface::OperatorTypeEnum::collection_operator;
}

template<class DerivedOperatorType>
void CollectionOperatorInterface<DerivedOperatorType>::perform_collection_move(
        RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->call_store_methods(comparisons);

    double hastings_ratio = this->propose(rng, comparisons, nthreads);
    comparisons->compute_log_likelihood_and_prior(true);

    double likelihood_ratio = 
        comparisons->get_log_likelihood() -
        comparisons->get_stored_log_likelihood();
    double prior_ratio = 
        comparisons->get_log_prior_density() -
        comparisons->get_stored_log_prior_density();
    double acceptance_probability =
            likelihood_ratio + 
            prior_ratio +
            hastings_ratio;
    double u = rng.uniform_real();
    if (u < std::exp(acceptance_probability)) {
        this->accept(comparisons->get_operator_schedule());
    }
    else {
        this->reject(comparisons->get_operator_schedule());
        this->call_restore_methods(comparisons);
    }
    comparisons->make_trees_clean();
    this->optimize(comparisons->get_operator_schedule(), acceptance_probability);
}


//////////////////////////////////////////////////////////////////////////////
// UnivariateTimeSizeScaler methods
//////////////////////////////////////////////////////////////////////////////

UnivariateTimeSizeScaler::UnivariateTimeSizeScaler(
        ) : CollectionOperatorInterface<ScaleOperator>() {
    this->op_ = ScaleOperator();
}

UnivariateTimeSizeScaler::UnivariateTimeSizeScaler(
        double weight) : CollectionOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator();
}

UnivariateTimeSizeScaler::UnivariateTimeSizeScaler(
        double weight,
        double scale) : CollectionOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator(scale);
}

void UnivariateTimeSizeScaler::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);
}

double UnivariateTimeSizeScaler::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    double hastings_ratio = 0.0;
    double size;
    double height;
    double hastings;
    for (unsigned int tree_idx = 0;
            tree_idx < comparisons->get_number_of_trees();
            ++tree_idx) {
        std::shared_ptr<PopulationTree> tree = comparisons->get_tree(tree_idx);
        if (! tree->population_sizes_are_fixed()) {
            size = tree->get_root_population_size();
            this->update(rng, size, hastings);
            hastings_ratio += hastings;
            tree->set_root_population_size(size);
            if (! tree->population_sizes_are_constrained()) {
                for (unsigned int i = 0; i < tree->get_leaf_node_count(); ++i) {
                    size = tree->get_child_population_size(i);
                    this->update(rng, size, hastings);
                    hastings_ratio += hastings;
                    tree->set_child_population_size(i, size);
                }
            }
        }
    }
    for (unsigned int height_idx = 0;
            height_idx < comparisons->get_number_of_events();
            ++height_idx) {
        height = comparisons->get_height(height_idx);
        this->update(rng, height, hastings);
        hastings_ratio += hastings;
        comparisons->set_height(height_idx, height);
    }
    comparisons->make_trees_dirty();
    return hastings_ratio;
}

std::string UnivariateTimeSizeScaler::target_parameter() const {
    return "node heights and population sizes";
}

std::string UnivariateTimeSizeScaler::get_name() const {
    return "UnivariateTimeSizeScaler";
}


//////////////////////////////////////////////////////////////////////////////
// UnivariateTimeSizeRateScaler methods
//////////////////////////////////////////////////////////////////////////////

UnivariateTimeSizeRateScaler::UnivariateTimeSizeRateScaler(
        ) : CollectionOperatorInterface<ScaleOperator>() {
    this->op_ = ScaleOperator();
}

UnivariateTimeSizeRateScaler::UnivariateTimeSizeRateScaler(
        double weight) : CollectionOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator();
}

UnivariateTimeSizeRateScaler::UnivariateTimeSizeRateScaler(
        double weight,
        double scale) : CollectionOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator(scale);
}

void UnivariateTimeSizeRateScaler::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);
}

double UnivariateTimeSizeRateScaler::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    double hastings_ratio = 0.0;
    double size;
    double height;
    double hastings;
    double rate;
    for (unsigned int tree_idx = 0;
            tree_idx < comparisons->get_number_of_trees();
            ++tree_idx) {
        std::shared_ptr<PopulationTree> tree = comparisons->get_tree(tree_idx);
        if (! tree->population_sizes_are_fixed()) {
            size = tree->get_root_population_size();
            this->update(rng, size, hastings);
            hastings_ratio += hastings;
            tree->set_root_population_size(size);
            if (! tree->population_sizes_are_constrained()) {
                for (unsigned int i = 0; i < tree->get_leaf_node_count(); ++i) {
                    size = tree->get_child_population_size(i);
                    this->update(rng, size, hastings);
                    hastings_ratio += hastings;
                    tree->set_child_population_size(i, size);
                }
            }
        }
        if (! tree->mutation_rate_is_fixed()) {
            rate = tree->get_mutation_rate();
            this->update(rng, rate, hastings);
            hastings_ratio += hastings;
            tree->set_mutation_rate(rate);
        }
    }
    for (unsigned int height_idx = 0;
            height_idx < comparisons->get_number_of_events();
            ++height_idx) {
        height = comparisons->get_height(height_idx);
        this->update(rng, height, hastings);
        hastings_ratio += hastings;
        comparisons->set_height(height_idx, height);
    }
    comparisons->make_trees_dirty();
    return hastings_ratio;
}

std::string UnivariateTimeSizeRateScaler::target_parameter() const {
    return "node heights, population sizes, and mutation rates";
}

std::string UnivariateTimeSizeRateScaler::get_name() const {
    return "UnivariateTimeSizeRateScaler";
}


//////////////////////////////////////////////////////////////////////////////
// UnivariateTimeMeanSizeRateScaler methods
//////////////////////////////////////////////////////////////////////////////

UnivariateTimeMeanSizeRateScaler::UnivariateTimeMeanSizeRateScaler(
        ) : CollectionOperatorInterface<ScaleOperator>() {
    this->op_ = ScaleOperator();
}

UnivariateTimeMeanSizeRateScaler::UnivariateTimeMeanSizeRateScaler(
        double weight) : CollectionOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator();
}

UnivariateTimeMeanSizeRateScaler::UnivariateTimeMeanSizeRateScaler(
        double weight,
        double scale) : CollectionOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator(scale);
}

void UnivariateTimeMeanSizeRateScaler::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);
}

double UnivariateTimeMeanSizeRateScaler::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    double hastings_ratio = 0.0;
    double size;
    double height;
    double hastings;
    double rate;
    for (unsigned int tree_idx = 0;
            tree_idx < comparisons->get_number_of_trees();
            ++tree_idx) {
        std::shared_ptr<PopulationTree> tree = comparisons->get_tree(tree_idx);
        if (! tree->mean_population_size_is_fixed()) {
            size = tree->get_mean_population_size();
            this->update(rng, size, hastings);
            hastings_ratio += hastings;
            tree->set_mean_population_size(size);
        }
        if (! tree->mutation_rate_is_fixed()) {
            rate = tree->get_mutation_rate();
            this->update(rng, rate, hastings);
            hastings_ratio += hastings;
            tree->set_mutation_rate(rate);
        }
    }
    for (unsigned int height_idx = 0;
            height_idx < comparisons->get_number_of_events();
            ++height_idx) {
        height = comparisons->get_height(height_idx);
        this->update(rng, height, hastings);
        hastings_ratio += hastings;
        comparisons->set_height(height_idx, height);
    }
    comparisons->make_trees_dirty();
    return hastings_ratio;
}

std::string UnivariateTimeMeanSizeRateScaler::target_parameter() const {
    return "node heights, population sizes, and mutation rates";
}

std::string UnivariateTimeMeanSizeRateScaler::get_name() const {
    return "UnivariateTimeMeanSizeRateScaler";
}


//////////////////////////////////////////////////////////////////////////////
// UnivariateCompositeTimeSizeScaler methods
//////////////////////////////////////////////////////////////////////////////

UnivariateCompositeTimeSizeScaler::UnivariateCompositeTimeSizeScaler(
        double weight,
        double scale) : CollectionOperatorInterface<Operator>(weight) {
    this->height_scaler_.op_.set_scale(scale);
    this->root_size_scaler_.op_.set_scale(scale);
    this->child_size_scaler_.op_.set_scale(scale);
}

void UnivariateCompositeTimeSizeScaler::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);
}

void UnivariateCompositeTimeSizeScaler::perform_collection_move(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->height_scaler_.operate(rng, comparisons, nthreads);
    this->root_size_scaler_.operate(rng, comparisons, nthreads);
    this->child_size_scaler_.operate(rng, comparisons, nthreads);
}

double UnivariateCompositeTimeSizeScaler::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->height_scaler_.operate(rng, comparisons, nthreads);
    this->root_size_scaler_.operate(rng, comparisons, nthreads);
    this->child_size_scaler_.operate(rng, comparisons, nthreads);
    return std::numeric_limits<double>::infinity();
}

void UnivariateCompositeTimeSizeScaler::scale_heights(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->height_scaler_.operate(rng, comparisons, nthreads);
}

void UnivariateCompositeTimeSizeScaler::scale_root_population_sizes(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->root_size_scaler_.operate(rng, comparisons, nthreads);
}

void UnivariateCompositeTimeSizeScaler::scale_child_population_sizes(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->child_size_scaler_.operate(rng, comparisons, nthreads);
}

std::string UnivariateCompositeTimeSizeScaler::target_parameter() const {
    return "node heights and population sizes";
}

std::string UnivariateCompositeTimeSizeScaler::get_name() const {
    return "UnivariateCompositeTimeSizeScaler";
}

std::string UnivariateCompositeTimeSizeScaler::to_string(const OperatorSchedule& os) const {
    std::ostringstream ss;
    ss << this->get_name() << "\t" 
       << this->get_number_accepted() << "\t"
       << this->get_number_rejected() << "\t"
       << this->get_weight() << "\t";

    if (os.get_total_weight() > 0.0) {
        ss << this->get_weight() / os.get_total_weight() << "\t";
    }
    else {
        ss << "nan\t";
    }

    double tuning = this->get_coercable_parameter_value();
    if (std::isnan(tuning)) {
        ss << "none\t";
    }
    else {
        ss << tuning << "\t";
    }
    ss << "\n";
    ss << this->height_scaler_.to_string(os);
    ss << this->root_size_scaler_.to_string(os);
    ss << this->child_size_scaler_.to_string(os);
    return ss.str();
}


//////////////////////////////////////////////////////////////////////////////
// UnivariateCompositeTimeSizeRateScaler methods
//////////////////////////////////////////////////////////////////////////////

UnivariateCompositeTimeSizeRateScaler::UnivariateCompositeTimeSizeRateScaler(
        double weight,
        double scale) : CollectionOperatorInterface<Operator>(weight) {
    this->height_scaler_.op_.set_scale(scale);
    this->root_size_scaler_.op_.set_scale(scale);
    this->child_size_scaler_.op_.set_scale(scale);
    this->mutation_rate_scaler_.op_.set_scale(scale);
}

void UnivariateCompositeTimeSizeRateScaler::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);
}

void UnivariateCompositeTimeSizeRateScaler::perform_collection_move(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->height_scaler_.operate(rng, comparisons, nthreads);
    this->root_size_scaler_.operate(rng, comparisons, nthreads);
    this->child_size_scaler_.operate(rng, comparisons, nthreads);
    this->mutation_rate_scaler_.operate(rng, comparisons, nthreads);
}

double UnivariateCompositeTimeSizeRateScaler::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->height_scaler_.operate(rng, comparisons, nthreads);
    this->root_size_scaler_.operate(rng, comparisons, nthreads);
    this->child_size_scaler_.operate(rng, comparisons, nthreads);
    this->mutation_rate_scaler_.operate(rng, comparisons, nthreads);
    return std::numeric_limits<double>::infinity();
}

void UnivariateCompositeTimeSizeRateScaler::scale_heights(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->height_scaler_.operate(rng, comparisons, nthreads);
}

void UnivariateCompositeTimeSizeRateScaler::scale_root_population_sizes(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->root_size_scaler_.operate(rng, comparisons, nthreads);
}

void UnivariateCompositeTimeSizeRateScaler::scale_child_population_sizes(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->child_size_scaler_.operate(rng, comparisons, nthreads);
}

void UnivariateCompositeTimeSizeRateScaler::scale_mutation_rates(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->mutation_rate_scaler_.operate(rng, comparisons, nthreads);
}

std::string UnivariateCompositeTimeSizeRateScaler::target_parameter() const {
    return "node heights, population sizes, and mutation rates";
}

std::string UnivariateCompositeTimeSizeRateScaler::get_name() const {
    return "UnivariateCompositeTimeSizeRateScaler";
}

std::string UnivariateCompositeTimeSizeRateScaler::to_string(const OperatorSchedule& os) const {
    std::ostringstream ss;
    ss << this->get_name() << "\t" 
       << this->get_number_accepted() << "\t"
       << this->get_number_rejected() << "\t"
       << this->get_weight() << "\t";

    if (os.get_total_weight() > 0.0) {
        ss << this->get_weight() / os.get_total_weight() << "\t";
    }
    else {
        ss << "nan\t";
    }

    double tuning = this->get_coercable_parameter_value();
    if (std::isnan(tuning)) {
        ss << "none\t";
    }
    else {
        ss << tuning << "\t";
    }
    ss << "\n";
    ss << this->height_scaler_.to_string(os);
    ss << this->root_size_scaler_.to_string(os);
    ss << this->child_size_scaler_.to_string(os);
    ss << this->mutation_rate_scaler_.to_string(os);
    return ss.str();
}


//////////////////////////////////////////////////////////////////////////////
// UnivariateCompositeTimeMeanSizeRateScaler methods
//////////////////////////////////////////////////////////////////////////////

UnivariateCompositeTimeMeanSizeRateScaler::UnivariateCompositeTimeMeanSizeRateScaler(
        double weight,
        double scale) : CollectionOperatorInterface<Operator>(weight) {
    this->height_scaler_.op_.set_scale(scale);
    this->size_scaler_.op_.set_scale(scale);
    this->size_multiplier_mixer_.op_.set_scale(scale);
    this->mutation_rate_scaler_.op_.set_scale(scale);
}

void UnivariateCompositeTimeMeanSizeRateScaler::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);
}

void UnivariateCompositeTimeMeanSizeRateScaler::perform_collection_move(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->height_scaler_.operate(rng, comparisons, nthreads);
    this->size_scaler_.operate(rng, comparisons, nthreads);
    this->size_multiplier_mixer_.operate(rng, comparisons, nthreads);
    this->mutation_rate_scaler_.operate(rng, comparisons, nthreads);
}

double UnivariateCompositeTimeMeanSizeRateScaler::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->height_scaler_.operate(rng, comparisons, nthreads);
    this->size_scaler_.operate(rng, comparisons, nthreads);
    this->size_multiplier_mixer_.operate(rng, comparisons, nthreads);
    this->mutation_rate_scaler_.operate(rng, comparisons, nthreads);
    return std::numeric_limits<double>::infinity();
}

void UnivariateCompositeTimeMeanSizeRateScaler::scale_heights(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->height_scaler_.operate(rng, comparisons, nthreads);
}

void UnivariateCompositeTimeMeanSizeRateScaler::scale_population_sizes(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->size_scaler_.operate(rng, comparisons, nthreads);
}

void UnivariateCompositeTimeMeanSizeRateScaler::mix_population_size_multipliers(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->size_multiplier_mixer_.operate(rng, comparisons, nthreads);
}

void UnivariateCompositeTimeMeanSizeRateScaler::scale_mutation_rates(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->mutation_rate_scaler_.operate(rng, comparisons, nthreads);
}

std::string UnivariateCompositeTimeMeanSizeRateScaler::target_parameter() const {
    return "node heights, population sizes, and mutation rates";
}

std::string UnivariateCompositeTimeMeanSizeRateScaler::get_name() const {
    return "UnivariateCompositeTimeMeanSizeRateScaler";
}

std::string UnivariateCompositeTimeMeanSizeRateScaler::to_string(const OperatorSchedule& os) const {
    std::ostringstream ss;
    ss << this->get_name() << "\t" 
       << this->get_number_accepted() << "\t"
       << this->get_number_rejected() << "\t"
       << this->get_weight() << "\t";

    if (os.get_total_weight() > 0.0) {
        ss << this->get_weight() / os.get_total_weight() << "\t";
    }
    else {
        ss << "nan\t";
    }

    double tuning = this->get_coercable_parameter_value();
    if (std::isnan(tuning)) {
        ss << "none\t";
    }
    else {
        ss << tuning << "\t";
    }
    ss << "\n";
    ss << this->height_scaler_.to_string(os);
    ss << this->size_scaler_.to_string(os);
    ss << this->size_multiplier_mixer_.to_string(os);
    ss << this->mutation_rate_scaler_.to_string(os);
    return ss.str();
}


//////////////////////////////////////////////////////////////////////////////
// ConcentrationScaler methods
//////////////////////////////////////////////////////////////////////////////

ConcentrationScaler::ConcentrationScaler() : CollectionOperatorInterface<ScaleOperator>() {
    this->op_ = ScaleOperator();
}

ConcentrationScaler::ConcentrationScaler(
        double weight) : CollectionOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator();
}

ConcentrationScaler::ConcentrationScaler(
        double weight,
        double scale) : CollectionOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator(scale);
}

void ConcentrationScaler::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);
}

double ConcentrationScaler::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    double v = comparisons->get_concentration();
    double hastings;
    this->update(rng, v, hastings);
    comparisons->set_concentration(v);
    return hastings;
}

std::string ConcentrationScaler::target_parameter() const {
    return "concentration";
}

std::string ConcentrationScaler::get_name() const {
    return "ConcentrationScaler";
}


//////////////////////////////////////////////////////////////////////////////
// FreqMover methods
//////////////////////////////////////////////////////////////////////////////

FreqMover::FreqMover() : TreeOperatorInterface<WindowOperator>() {
    this->op_ = WindowOperator();
}

FreqMover::FreqMover(
        unsigned int tree_index) : TreeOperatorInterface<WindowOperator>(tree_index) {
    this->op_ = WindowOperator();
}

FreqMover::FreqMover(
        double weight) : TreeOperatorInterface<WindowOperator>(weight) {
    this->op_ = WindowOperator();
}

FreqMover::FreqMover(
        unsigned int tree_index,
        double weight) : TreeOperatorInterface<WindowOperator>(tree_index, weight) {
    this->op_ = WindowOperator();
}

FreqMover::FreqMover(
        double weight,
        double window_size) : TreeOperatorInterface<WindowOperator>(weight) {
    this->op_ = WindowOperator(window_size);
}

FreqMover::FreqMover(
        unsigned int tree_index,
        double weight,
        double window_size) : TreeOperatorInterface<WindowOperator>(tree_index, weight) {
    this->op_ = WindowOperator(window_size);
}

void FreqMover::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);
}

double FreqMover::propose(
        RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int tree_index) {
    double freq_1 = comparisons->get_tree(tree_index)->get_freq_1();
    double hastings;
    this->update(rng, freq_1, hastings);
    if ((freq_1 <= 0.0) || (freq_1 > 1.0)) {
        return -std::numeric_limits<double>::infinity();
    }
    comparisons->get_tree(tree_index)->set_freq_1(freq_1);
    return hastings; 
}

std::string FreqMover::target_parameter() const {
    return "freq 1";
}

std::string FreqMover::get_name() const {
    std::string name = "FreqMover";
    if (this->tree_index_ > -1) {
        name += std::to_string(this->tree_index_);
    }
    return name;
}


//////////////////////////////////////////////////////////////////////////////
// MutationRateScaler methods
//////////////////////////////////////////////////////////////////////////////

MutationRateScaler::MutationRateScaler(
        ) : TreeOperatorInterface<ScaleOperator>() {
    this->op_ = ScaleOperator();
}

MutationRateScaler::MutationRateScaler(
        unsigned int tree_index) : TreeOperatorInterface<ScaleOperator>(tree_index) {
    this->op_ = ScaleOperator();
}

MutationRateScaler::MutationRateScaler(
        double weight) : TreeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator();
}

MutationRateScaler::MutationRateScaler(
        unsigned int tree_index,
        double weight) : TreeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator();
}

MutationRateScaler::MutationRateScaler(
        double weight,
        double scale) : TreeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator(scale);
}

MutationRateScaler::MutationRateScaler(
        unsigned int tree_index,
        double weight,
        double scale) : TreeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator(scale);
}

void MutationRateScaler::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);
}

double MutationRateScaler::propose(
        RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int tree_index) {
    double v = comparisons->get_tree(tree_index)->get_mutation_rate();
    double hastings;
    this->update(rng, v, hastings);
    comparisons->get_tree(tree_index)->set_mutation_rate(v);
    return hastings;
}

std::string MutationRateScaler::target_parameter() const {
    return "mutation rate";
}

std::string MutationRateScaler::get_name() const {
    std::string name = "MutationRateScaler";
    if (this->tree_index_ > -1) {
        name += std::to_string(this->tree_index_);
    }
    return name;
}


//////////////////////////////////////////////////////////////////////////////
// RelativePopulationSizeMixer methods
//////////////////////////////////////////////////////////////////////////////

RelativePopulationSizeMixer::RelativePopulationSizeMixer(
        ) : TreeOperatorInterface<ScaleOperator>() {
    this->op_ = ScaleOperator();
}

RelativePopulationSizeMixer::RelativePopulationSizeMixer(
        unsigned int tree_index) : TreeOperatorInterface<ScaleOperator>(tree_index) {
    this->op_ = ScaleOperator();
}

RelativePopulationSizeMixer::RelativePopulationSizeMixer(
        double weight) : TreeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator();
}

RelativePopulationSizeMixer::RelativePopulationSizeMixer(
        unsigned int tree_index,
        double weight) : TreeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator();
}

RelativePopulationSizeMixer::RelativePopulationSizeMixer(
        double weight,
        double scale) : TreeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator(scale);
}

RelativePopulationSizeMixer::RelativePopulationSizeMixer(
        unsigned int tree_index,
        double weight,
        double scale) : TreeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator(scale);
}

void RelativePopulationSizeMixer::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);
}

double RelativePopulationSizeMixer::propose(
        RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int tree_index) {
    std::vector<double> old_pop_proportions = comparisons->get_tree(tree_index)->get_population_sizes_as_proportions();
    ECOEVOLITY_ASSERT(old_pop_proportions.size() > 1);
    std::vector<double> forward_dir_parameters = old_pop_proportions;
    for (unsigned int i = 0; i < forward_dir_parameters.size(); ++i) {
        forward_dir_parameters.at(i) = 1.0 + (forward_dir_parameters.at(i) * (1.0 / this->op_.get_scale()));
    }
    DirichletDistribution dir_forward = DirichletDistribution(forward_dir_parameters);
    std::vector<double> new_pop_proportions = dir_forward.draw(rng);

    comparisons->get_tree(tree_index)->set_population_sizes_as_proportions(new_pop_proportions);

    std::vector<double> reverse_dir_parameters = new_pop_proportions;
    for (unsigned int i = 0; i < reverse_dir_parameters.size(); ++i) {
        reverse_dir_parameters.at(i) = 1.0 + (reverse_dir_parameters.at(i) * (1.0 / this->op_.get_scale()));
    }
    DirichletDistribution dir_reverse = DirichletDistribution(reverse_dir_parameters);

    double hastings = (dir_reverse.ln_pdf(old_pop_proportions) -
            dir_forward.ln_pdf(new_pop_proportions));

    return hastings;
}

std::string RelativePopulationSizeMixer::target_parameter() const {
    return "population size multipliers";
}

std::string RelativePopulationSizeMixer::get_name() const {
    std::string name = "RelativePopulationSizeMixer";
    if (this->tree_index_ > -1) {
        name += std::to_string(this->tree_index_);
    }
    return name;
}


//////////////////////////////////////////////////////////////////////////////
// MeanPopulationSizeScaler methods
//////////////////////////////////////////////////////////////////////////////

MeanPopulationSizeScaler::MeanPopulationSizeScaler(
        ) : TreeOperatorInterface<ScaleOperator>() {
    this->op_ = ScaleOperator();
}

MeanPopulationSizeScaler::MeanPopulationSizeScaler(
        unsigned int tree_index) : TreeOperatorInterface<ScaleOperator>(tree_index) {
    this->op_ = ScaleOperator();
}

MeanPopulationSizeScaler::MeanPopulationSizeScaler(
        double weight) : TreeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator();
}

MeanPopulationSizeScaler::MeanPopulationSizeScaler(
        unsigned int tree_index,
        double weight) : TreeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator();
}

MeanPopulationSizeScaler::MeanPopulationSizeScaler(
        double weight,
        double scale) : TreeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator(scale);
}

MeanPopulationSizeScaler::MeanPopulationSizeScaler(
        unsigned int tree_index,
        double weight,
        double scale) : TreeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator(scale);
}

void MeanPopulationSizeScaler::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);
}

double MeanPopulationSizeScaler::propose(
        RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int tree_index) {
    double size = comparisons->get_tree(tree_index)->get_mean_population_size();

    double hastings;
    this->update(rng, size, hastings);

    // avoid zero division to get coalescence rate
    if (size <= 0.0) {
        return -std::numeric_limits<double>::infinity();
    }

    comparisons->get_tree(tree_index)->set_mean_population_size(size);
    return hastings;
}

std::string MeanPopulationSizeScaler::target_parameter() const {
    return "population size";
}

std::string MeanPopulationSizeScaler::get_name() const {
    std::string name = "MeanPopulationSizeScaler";
    if (this->tree_index_ > -1) {
        name += std::to_string(this->tree_index_);
    }
    return name;
}


//////////////////////////////////////////////////////////////////////////////
// LeafPopulationSizeScaler methods
//////////////////////////////////////////////////////////////////////////////

LeafPopulationSizeScaler::LeafPopulationSizeScaler(
        ) : TreeOperatorInterface<ScaleOperator>() {
    this->op_ = ScaleOperator();
}

LeafPopulationSizeScaler::LeafPopulationSizeScaler(
        unsigned int tree_index) : TreeOperatorInterface<ScaleOperator>(tree_index) {
    this->op_ = ScaleOperator();
}

LeafPopulationSizeScaler::LeafPopulationSizeScaler(
        double weight) : TreeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator();
}

LeafPopulationSizeScaler::LeafPopulationSizeScaler(
        unsigned int tree_index,
        double weight) : TreeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator();
}

LeafPopulationSizeScaler::LeafPopulationSizeScaler(
        double weight,
        double scale) : TreeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator(scale);
}

LeafPopulationSizeScaler::LeafPopulationSizeScaler(
        unsigned int tree_index,
        double weight,
        double scale) : TreeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator(scale);
}

void LeafPopulationSizeScaler::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);
}

double LeafPopulationSizeScaler::propose(
        RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int tree_index) {
    int pop_idx = rng.uniform_int(0, comparisons->get_tree(tree_index)->get_leaf_node_count() - 1);
    double size = comparisons->get_tree(tree_index)->get_child_population_size(pop_idx);

    double hastings;
    this->update(rng, size, hastings);

    // avoid zero division to get coalescence rate
    if (size <= 0.0) {
        return -std::numeric_limits<double>::infinity();
    }

    comparisons->get_tree(tree_index)->set_child_population_size(pop_idx, size);
    return hastings;
}

std::string LeafPopulationSizeScaler::target_parameter() const {
    return "population size";
}

std::string LeafPopulationSizeScaler::get_name() const {
    std::string name = "LeafPopulationSizeScaler";
    if (this->tree_index_ > -1) {
        name += std::to_string(this->tree_index_);
    }
    return name;
}


//////////////////////////////////////////////////////////////////////////////
// RootPopulationSizeScaler methods
//////////////////////////////////////////////////////////////////////////////

RootPopulationSizeScaler::RootPopulationSizeScaler(
        ) : TreeOperatorInterface<ScaleOperator>() {
    this->op_ = ScaleOperator();
}

RootPopulationSizeScaler::RootPopulationSizeScaler(
        unsigned int tree_index) : TreeOperatorInterface<ScaleOperator>(tree_index) {
    this->op_ = ScaleOperator();
}

RootPopulationSizeScaler::RootPopulationSizeScaler(
        double weight) : TreeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator();
}

RootPopulationSizeScaler::RootPopulationSizeScaler(
        unsigned int tree_index,
        double weight) : TreeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator();
}

RootPopulationSizeScaler::RootPopulationSizeScaler(
        double weight,
        double scale) : TreeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator(scale);
}

RootPopulationSizeScaler::RootPopulationSizeScaler(
        unsigned int tree_index,
        double weight,
        double scale) : TreeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator(scale);
}

void RootPopulationSizeScaler::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);
}

double RootPopulationSizeScaler::propose(
        RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int tree_index) {
    double size = comparisons->get_tree(tree_index)->get_root_population_size();

    double hastings;
    this->update(rng, size, hastings);

    // avoid zero division to get coalescence rate
    if (size <= 0.0) {
        return -std::numeric_limits<double>::infinity();
    }

    comparisons->get_tree(tree_index)->set_root_population_size(size);

    return hastings;
}

std::string RootPopulationSizeScaler::target_parameter() const {
    return "population size";
}

std::string RootPopulationSizeScaler::get_name() const {
    std::string name = "RootPopulationSizeScaler";
    if (this->tree_index_ > -1) {
        name += std::to_string(this->tree_index_);
    }
    return name;
}


//////////////////////////////////////////////////////////////////////////////
// RootRelativePopulationSizeMover methods
//////////////////////////////////////////////////////////////////////////////

RootRelativePopulationSizeMover::RootRelativePopulationSizeMover(
        ) : TreeOperatorInterface<WindowOperator>() {
    this->op_ = WindowOperator();
}

RootRelativePopulationSizeMover::RootRelativePopulationSizeMover(
        unsigned int tree_index) : TreeOperatorInterface<WindowOperator>(tree_index) {
    this->op_ = WindowOperator();
}

RootRelativePopulationSizeMover::RootRelativePopulationSizeMover(
        double weight) : TreeOperatorInterface<WindowOperator>(weight) {
    this->op_ = WindowOperator();
}

RootRelativePopulationSizeMover::RootRelativePopulationSizeMover(
        unsigned int tree_index,
        double weight) : TreeOperatorInterface<WindowOperator>(tree_index, weight) {
    this->op_ = WindowOperator();
}

RootRelativePopulationSizeMover::RootRelativePopulationSizeMover(
        double weight,
        double scale) : TreeOperatorInterface<WindowOperator>(weight) {
    this->op_ = WindowOperator(scale);
}

RootRelativePopulationSizeMover::RootRelativePopulationSizeMover(
        unsigned int tree_index,
        double weight,
        double scale) : TreeOperatorInterface<WindowOperator>(tree_index, weight) {
    this->op_ = WindowOperator(scale);
}

void RootRelativePopulationSizeMover::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);
}

double RootRelativePopulationSizeMover::propose(
        RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int tree_index) {

    std::vector<double> sizes = comparisons->get_tree(tree_index)->get_population_sizes();
    ECOEVOLITY_ASSERT(sizes.size() > 1);

    double addend = this->op_.get_move_amount(rng);

    unsigned int num_nonroot_branches = sizes.size() - 1;
    sizes.at(num_nonroot_branches) += addend;
    if (sizes.at(num_nonroot_branches) <= 0.0) {
        return -std::numeric_limits<double>::infinity();
    }
    for (unsigned int i = 0; i < num_nonroot_branches; ++i) {
        sizes.at(i) -= (addend / (double)num_nonroot_branches);
        if (sizes.at(i) <= 0.0) {
            return -std::numeric_limits<double>::infinity();
        }
    }

    comparisons->get_tree(tree_index)->set_population_sizes(sizes);

    return 0.0;
}

std::string RootRelativePopulationSizeMover::target_parameter() const {
    return "population size multipliers";
}

std::string RootRelativePopulationSizeMover::get_name() const {
    std::string name = "RootRelativePopulationSizeMover";
    if (this->tree_index_ > -1) {
        name += std::to_string(this->tree_index_);
    }
    return name;
}


//////////////////////////////////////////////////////////////////////////////
// LeafRelativePopulationSizeMover methods
//////////////////////////////////////////////////////////////////////////////

LeafRelativePopulationSizeMover::LeafRelativePopulationSizeMover(
        ) : TreeOperatorInterface<WindowOperator>() {
    this->op_ = WindowOperator();
}

LeafRelativePopulationSizeMover::LeafRelativePopulationSizeMover(
        unsigned int tree_index) : TreeOperatorInterface<WindowOperator>(tree_index) {
    this->op_ = WindowOperator();
}

LeafRelativePopulationSizeMover::LeafRelativePopulationSizeMover(
        double weight) : TreeOperatorInterface<WindowOperator>(weight) {
    this->op_ = WindowOperator();
}

LeafRelativePopulationSizeMover::LeafRelativePopulationSizeMover(
        unsigned int tree_index,
        double weight) : TreeOperatorInterface<WindowOperator>(tree_index, weight) {
    this->op_ = WindowOperator();
}

LeafRelativePopulationSizeMover::LeafRelativePopulationSizeMover(
        double weight,
        double scale) : TreeOperatorInterface<WindowOperator>(weight) {
    this->op_ = WindowOperator(scale);
}

LeafRelativePopulationSizeMover::LeafRelativePopulationSizeMover(
        unsigned int tree_index,
        double weight,
        double scale) : TreeOperatorInterface<WindowOperator>(tree_index, weight) {
    this->op_ = WindowOperator(scale);
}

void LeafRelativePopulationSizeMover::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);
}

double LeafRelativePopulationSizeMover::propose(
        RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int tree_index) {

    std::vector<double> sizes = comparisons->get_tree(tree_index)->get_population_sizes();
    ECOEVOLITY_ASSERT(sizes.size() > 1);

    double addend = this->op_.get_move_amount(rng);

    unsigned int num_nonroot_branches = sizes.size() - 1;
    unsigned int pop_idx = rng.uniform_int(0, num_nonroot_branches - 1);
    for (unsigned int i = 0; i < sizes.size(); ++i) {
        if (i == pop_idx) {
            sizes.at(i) += addend;
            if (sizes.at(i) <= 0.0) {
                return -std::numeric_limits<double>::infinity();
            }
        }
        else {
            sizes.at(i) -= (addend / (double)num_nonroot_branches);
            if (sizes.at(i) <= 0.0) {
                return -std::numeric_limits<double>::infinity();
            }
        }
    }

    comparisons->get_tree(tree_index)->set_population_sizes(sizes);

    return 0.0;
}

std::string LeafRelativePopulationSizeMover::target_parameter() const {
    return "population size multipliers";
}

std::string LeafRelativePopulationSizeMover::get_name() const {
    std::string name = "LeafRelativePopulationSizeMover";
    if (this->tree_index_ > -1) {
        name += std::to_string(this->tree_index_);
    }
    return name;
}


//////////////////////////////////////////////////////////////////////////////
// EventTimeScaler methods
//////////////////////////////////////////////////////////////////////////////

EventTimeScaler::EventTimeScaler(
        ) : TimeOperatorInterface<ScaleOperator>() {
    this->op_ = ScaleOperator();
}

EventTimeScaler::EventTimeScaler(
        unsigned int tree_index) : TimeOperatorInterface<ScaleOperator>(tree_index) {
    this->op_ = ScaleOperator();
}

EventTimeScaler::EventTimeScaler(
        double weight) : TimeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator();
}

EventTimeScaler::EventTimeScaler(
        unsigned int tree_index,
        double weight) : TimeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator();
}

EventTimeScaler::EventTimeScaler(
        double weight,
        double scale) : TimeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator(scale);
}

EventTimeScaler::EventTimeScaler(
        unsigned int tree_index,
        double weight,
        double scale) : TimeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator(scale);
}

void EventTimeScaler::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);
}

double EventTimeScaler::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int index) {
    if (this->tree_index_ < 0) {
        return this->propose_by_height(rng, comparisons, index);
    }
    return this->propose_by_tree(rng, comparisons, index);
}

double EventTimeScaler::propose_by_height(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int height_index) {
    double h = comparisons->get_height(height_index);
    double hastings;
    this->update(rng, h, hastings);
    comparisons->set_height(height_index, h);
    return hastings;
}

double EventTimeScaler::propose_by_tree(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int tree_index) {
    double h = comparisons->get_height_of_tree(tree_index);
    double hastings;
    this->update(rng, h, hastings);
    comparisons->set_height_of_tree(tree_index, h);
    return hastings;
}

std::string EventTimeScaler::get_name() const {
    std::string name = "EventTimeScaler";
    if (this->tree_index_ > -1) {
        name += std::to_string(this->tree_index_);
    }
    return name;
}

std::string EventTimeScaler::target_parameter() const {
    return "node height";
}


//////////////////////////////////////////////////////////////////////////////
// EventTimeMover methods
//////////////////////////////////////////////////////////////////////////////

EventTimeMover::EventTimeMover(
        ) : TimeOperatorInterface<WindowOperator>() {
    this->op_ = WindowOperator();
}

EventTimeMover::EventTimeMover(
        unsigned int tree_index) : TimeOperatorInterface<WindowOperator>(tree_index) {
    this->op_ = WindowOperator();
}

EventTimeMover::EventTimeMover(
        double weight) : TimeOperatorInterface<WindowOperator>(weight) {
    this->op_ = WindowOperator();
}

EventTimeMover::EventTimeMover(
        unsigned int tree_index,
        double weight) : TimeOperatorInterface<WindowOperator>(tree_index, weight) {
    this->op_ = WindowOperator();
}

EventTimeMover::EventTimeMover(
        double weight,
        double window_size) : TimeOperatorInterface<WindowOperator>(weight) {
    this->op_ = WindowOperator(window_size);
}

EventTimeMover::EventTimeMover(
        unsigned int tree_index,
        double weight,
        double window_size) : TimeOperatorInterface<WindowOperator>(tree_index, weight) {
    this->op_ = WindowOperator(window_size);
}

void EventTimeMover::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);
}

double EventTimeMover::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int index) {
    if (this->tree_index_ < 0) {
        return this->propose_by_height(rng, comparisons, index);
    }
    return this->propose_by_tree(rng, comparisons, index);
}

double EventTimeMover::propose_by_height(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int height_index) {
    double h = comparisons->get_height(height_index);
    double hastings;
    this->update(rng, h, hastings);
    if (h < 0.0) {
        return -std::numeric_limits<double>::infinity();
    }
    comparisons->set_height(height_index, h);
    return hastings;
}

double EventTimeMover::propose_by_tree(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int tree_index) {
    double h = comparisons->get_height_of_tree(tree_index);
    double hastings;
    this->update(rng, h, hastings);
    if (h < 0.0) {
        return -std::numeric_limits<double>::infinity();
    }
    comparisons->set_height_of_tree(tree_index, h);
    return hastings;
}

std::string EventTimeMover::get_name() const {
    std::string name = "EventTimeMover";
    if (this->tree_index_ > -1) {
        name += std::to_string(this->tree_index_);
    }
    return name;
}

std::string EventTimeMover::target_parameter() const {
    return "node height";
}


//////////////////////////////////////////////////////////////////////////////
// TimeSizeMixer methods
//////////////////////////////////////////////////////////////////////////////

TimeSizeMixer::TimeSizeMixer(
        ) : TimeOperatorInterface<ScaleOperator>() {
    this->op_ = ScaleOperator();
}

TimeSizeMixer::TimeSizeMixer(
        unsigned int tree_index) : TimeOperatorInterface<ScaleOperator>(tree_index) {
    this->op_ = ScaleOperator();
}

TimeSizeMixer::TimeSizeMixer(
        double weight) : TimeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator();
}

TimeSizeMixer::TimeSizeMixer(
        unsigned int tree_index,
        double weight) : TimeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator();
}

TimeSizeMixer::TimeSizeMixer(
        double weight,
        double scale) : TimeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator(scale);
}

TimeSizeMixer::TimeSizeMixer(
        unsigned int tree_index,
        double weight,
        double scale) : TimeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator(scale);
}

void TimeSizeMixer::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);

    // Perform sweep of univariate moves
    if (this->tree_index_ < 0) {
        for (std::shared_ptr<OperatorInterface> time_op : comparisons->get_time_operators()) {
            time_op->operate(rng, comparisons, nthreads);
        }
        for (std::shared_ptr<OperatorInterface> tree_op : comparisons->get_tree_operators()) {
            tree_op->operate(rng, comparisons, nthreads);
        }
        return;
    }
    std::vector< std::shared_ptr<OperatorInterface> > time_operators = comparisons->get_time_operators(this->tree_index_);
    if (time_operators.size() < 1) {
        time_operators = comparisons->get_time_operators();
    }
    for (std::shared_ptr<OperatorInterface> time_op : time_operators) {
        time_op->operate(rng, comparisons, nthreads);
    }
    for (std::shared_ptr<OperatorInterface> tree_op : comparisons->get_tree_operators(this->tree_index_)) {
        tree_op->operate(rng, comparisons, nthreads);
    }
}

double TimeSizeMixer::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int index) {
    if (this->tree_index_ < 0) {
        return this->propose_by_height(rng, comparisons, index);
    }
    return this->propose_by_tree(rng, comparisons, index);
}

double TimeSizeMixer::propose_by_height(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int height_index) {
    double multiplier = this->op_.get_move_amount(rng);

    double new_height = comparisons->get_height(height_index) * multiplier;
    comparisons->set_height(height_index, new_height);

    int ndimensions = 1; // for height scaled above
    int number_of_free_parameters_scaled = 0;
    int number_of_free_parameters_inverse_scaled = 0;
    std::vector<unsigned int> tree_indices = comparisons->get_indices_of_mapped_trees(height_index);
    for (auto tree_idx : tree_indices) {
        std::shared_ptr<PopulationTree> tree = comparisons->get_tree(tree_idx);
        if (! tree->population_sizes_are_fixed()) {
            if (tree->population_sizes_are_constrained()) {
                number_of_free_parameters_scaled += tree->scale_root_population_size(multiplier);
            }
            else {
                number_of_free_parameters_inverse_scaled += tree->scale_root_population_size(1.0/multiplier);
                unsigned int nleaves = tree->get_leaf_node_count();
                for (unsigned int i = 0; i < nleaves; ++i) { 
                    tree->set_child_population_size(i, 
                            tree->get_child_population_size(i) * multiplier);
                }
                number_of_free_parameters_scaled += nleaves;
            }
        }

    }
    ndimensions += (number_of_free_parameters_scaled - number_of_free_parameters_inverse_scaled);
    return std::log(multiplier) * ndimensions;
}

double TimeSizeMixer::propose_by_tree(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int tree_index) {
    double multiplier = this->op_.get_move_amount(rng);

    double new_height = comparisons->get_height_of_tree(tree_index) * multiplier;

    comparisons->set_height_of_tree(tree_index, new_height);

    std::shared_ptr<PopulationTree> tree = comparisons->get_tree(this->tree_index_);

    int ndimensions = 1; // for height scaled above
    int number_of_free_parameters_scaled = 0;
    int number_of_free_parameters_inverse_scaled = 0;
    if (! tree->population_sizes_are_fixed()) {
        if (tree->population_sizes_are_constrained()) {
            number_of_free_parameters_scaled += tree->scale_root_population_size(multiplier);
        }
        else {
            number_of_free_parameters_inverse_scaled += tree->scale_root_population_size(1.0/multiplier);
            unsigned int nleaves = tree->get_leaf_node_count();
            for (unsigned int i = 0; i < nleaves; ++i) { 
                tree->set_child_population_size(i, 
                        tree->get_child_population_size(i) * multiplier);
            }
            number_of_free_parameters_scaled += nleaves;
        }
    }
    ndimensions += (number_of_free_parameters_scaled - number_of_free_parameters_inverse_scaled);
    return std::log(multiplier) * ndimensions;
}

std::string TimeSizeMixer::get_name() const {
    std::string name = "TimeSizeMixer";
    if (this->tree_index_ > -1) {
        name += std::to_string(this->tree_index_);
    }
    return name;
}

std::string TimeSizeMixer::target_parameter() const {
    return "node heights and population sizes";
}

std::string TimeSizeMixer::to_string(const OperatorSchedule& os) const {
    std::ostringstream ss;
    ss << this->get_name() << "\t" 
       << this->get_number_accepted() << "\t"
       << this->get_number_rejected() << "\t"
       << this->get_weight() << "\t";

    if (os.get_total_weight() > 0.0) {
        ss << this->get_weight() / os.get_total_weight() << "\t";
    }
    else {
        ss << "nan\t";
    }

    double tuning = this->get_coercable_parameter_value();
    if (std::isnan(tuning)) {
        ss << "none\t";
    }
    else {
        ss << tuning << "\t";
    }
    ss << "\n";
    return ss.str();
}


//////////////////////////////////////////////////////////////////////////////
// TimeRootSizeMixer methods
//////////////////////////////////////////////////////////////////////////////

TimeRootSizeMixer::TimeRootSizeMixer(
        unsigned int tree_index) : TimeOperatorInterface<ScaleOperator>(tree_index) {
    this->op_ = ScaleOperator();
}

TimeRootSizeMixer::TimeRootSizeMixer(
        unsigned int tree_index,
        double weight) : TimeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator();
}

TimeRootSizeMixer::TimeRootSizeMixer(
        unsigned int tree_index,
        double weight,
        double scale) : TimeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator(scale);
}

void TimeRootSizeMixer::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);

    // Perform sweep of univariate moves
    if (this->tree_index_ < 0) {
        for (std::shared_ptr<OperatorInterface> time_op : comparisons->get_time_operators()) {
            time_op->operate(rng, comparisons, nthreads);
        }
        for (std::shared_ptr<OperatorInterface> tree_op : comparisons->get_tree_operators()) {
            tree_op->operate(rng, comparisons, nthreads);
        }
        return;
    }
    std::vector< std::shared_ptr<OperatorInterface> > time_operators = comparisons->get_time_operators(this->tree_index_);
    if (time_operators.size() < 1) {
        time_operators = comparisons->get_time_operators();
    }
    for (std::shared_ptr<OperatorInterface> time_op : time_operators) {
        time_op->operate(rng, comparisons, nthreads);
    }
    for (std::shared_ptr<OperatorInterface> tree_op : comparisons->get_tree_operators(this->tree_index_)) {
        tree_op->operate(rng, comparisons, nthreads);
    }
}

double TimeRootSizeMixer::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int index) {
    if (this->tree_index_ < 0) {
        return this->propose_by_height(rng, comparisons, index);
    }
    return this->propose_by_tree(rng, comparisons, index);
}

double TimeRootSizeMixer::propose_by_height(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int height_index) {
    return -std::numeric_limits<double>::infinity();
}

double TimeRootSizeMixer::propose_by_tree(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int tree_index) {
    std::shared_ptr<PopulationTree> tree = comparisons->get_tree(tree_index);
    if (tree->root_population_size_is_fixed()) {
        return -std::numeric_limits<double>::infinity();
    }

    double multiplier = this->op_.get_move_amount(rng);

    double old_size = tree->get_root_population_size();
    double new_size = old_size * multiplier;
    double delta_Nmu = old_size - new_size;
    double new_height = comparisons->get_height_of_tree(tree_index) + (2.0 * delta_Nmu);

    if (new_height < 0.0) {
        return -std::numeric_limits<double>::infinity();
    }

    tree->set_root_population_size(new_size);
    comparisons->set_height_of_tree(tree_index, new_height);

    return std::log(multiplier);
}

std::string TimeRootSizeMixer::get_name() const {
    std::string name = "TimeRootSizeMixer";
    if (this->tree_index_ > -1) {
        name += std::to_string(this->tree_index_);
    }
    return name;
}

std::string TimeRootSizeMixer::target_parameter() const {
    return "node heights and population sizes";
}

std::string TimeRootSizeMixer::to_string(const OperatorSchedule& os) const {
    std::ostringstream ss;
    ss << this->get_name() << "\t" 
       << this->get_number_accepted() << "\t"
       << this->get_number_rejected() << "\t"
       << this->get_weight() << "\t";

    if (os.get_total_weight() > 0.0) {
        ss << this->get_weight() / os.get_total_weight() << "\t";
    }
    else {
        ss << "nan\t";
    }

    double tuning = this->get_coercable_parameter_value();
    if (std::isnan(tuning)) {
        ss << "none\t";
    }
    else {
        ss << tuning << "\t";
    }
    ss << "\n";
    return ss.str();
}


//////////////////////////////////////////////////////////////////////////////
// TimeSizeScaler methods
//////////////////////////////////////////////////////////////////////////////

TimeSizeScaler::TimeSizeScaler(
        ) : TimeOperatorInterface<ScaleOperator>() {
    this->op_ = ScaleOperator();
}

TimeSizeScaler::TimeSizeScaler(
        unsigned int tree_index) : TimeOperatorInterface<ScaleOperator>(tree_index) {
    this->op_ = ScaleOperator();
}

TimeSizeScaler::TimeSizeScaler(
        double weight) : TimeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator();
}

TimeSizeScaler::TimeSizeScaler(
        unsigned int tree_index,
        double weight) : TimeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator();
}

TimeSizeScaler::TimeSizeScaler(
        double weight,
        double scale) : TimeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator(scale);
}

TimeSizeScaler::TimeSizeScaler(
        unsigned int tree_index,
        double weight,
        double scale) : TimeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator(scale);
}

void TimeSizeScaler::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);

    // Perform sweep of univariate moves
    if (this->tree_index_ < 0) {
        for (std::shared_ptr<OperatorInterface> time_op : comparisons->get_time_operators()) {
            time_op->operate(rng, comparisons, nthreads);
        }
        for (std::shared_ptr<OperatorInterface> tree_op : comparisons->get_tree_operators()) {
            tree_op->operate(rng, comparisons, nthreads);
        }
        return;
    }
    std::vector< std::shared_ptr<OperatorInterface> > time_operators = comparisons->get_time_operators(this->tree_index_);
    if (time_operators.size() < 1) {
        time_operators = comparisons->get_time_operators();
    }
    for (std::shared_ptr<OperatorInterface> time_op : time_operators) {
        time_op->operate(rng, comparisons, nthreads);
    }
    for (std::shared_ptr<OperatorInterface> tree_op : comparisons->get_tree_operators(this->tree_index_)) {
        tree_op->operate(rng, comparisons, nthreads);
    }
}

double TimeSizeScaler::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int index) {
    if (this->tree_index_ < 0) {
        return this->propose_by_height(rng, comparisons, index);
    }
    return this->propose_by_tree(rng, comparisons, index);
}

double TimeSizeScaler::propose_by_height(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int height_index) {
    double multiplier = this->op_.get_move_amount(rng);
    unsigned int number_of_free_parameters_scaled = 0;

    double new_height = comparisons->get_height(height_index) * multiplier;
    comparisons->set_height(height_index, new_height);
    ++number_of_free_parameters_scaled;
    std::vector<unsigned int> tree_indices = comparisons->get_indices_of_mapped_trees(height_index);

    for (auto tree_idx : tree_indices) {
        unsigned int nparameters_scaled = comparisons->get_tree(tree_idx)->scale_all_population_sizes(multiplier);
        number_of_free_parameters_scaled += nparameters_scaled;
    }

    return std::log(multiplier) * number_of_free_parameters_scaled;
}

double TimeSizeScaler::propose_by_tree(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int tree_index) {
    double multiplier = this->op_.get_move_amount(rng);
    unsigned int number_of_free_parameters_scaled = 0;

    double new_height = comparisons->get_height_of_tree(tree_index) * multiplier;
    comparisons->set_height_of_tree(tree_index, new_height);
    ++number_of_free_parameters_scaled;

    unsigned int nparameters_scaled = comparisons->get_tree(tree_index)->scale_all_population_sizes(multiplier);
    number_of_free_parameters_scaled += nparameters_scaled;

    return std::log(multiplier) * number_of_free_parameters_scaled;
}

std::string TimeSizeScaler::get_name() const {
    std::string name = "TimeSizeScaler";
    if (this->tree_index_ > -1) {
        name += std::to_string(this->tree_index_);
    }
    return name;
}

std::string TimeSizeScaler::target_parameter() const {
    return "node heights and population sizes";
}

std::string TimeSizeScaler::to_string(const OperatorSchedule& os) const {
    std::ostringstream ss;
    ss << this->get_name() << "\t" 
       << this->get_number_accepted() << "\t"
       << this->get_number_rejected() << "\t"
       << this->get_weight() << "\t";

    if (os.get_total_weight() > 0.0) {
        ss << this->get_weight() / os.get_total_weight() << "\t";
    }
    else {
        ss << "nan\t";
    }

    double tuning = this->get_coercable_parameter_value();
    if (std::isnan(tuning)) {
        ss << "none\t";
    }
    else {
        ss << tuning << "\t";
    }
    ss << "\n";
    return ss.str();
}


//////////////////////////////////////////////////////////////////////////////
// TimeSizeRateMixer methods
//////////////////////////////////////////////////////////////////////////////

TimeSizeRateMixer::TimeSizeRateMixer(
        ) : TimeOperatorInterface<ScaleOperator>() {
    this->op_ = ScaleOperator();
}

TimeSizeRateMixer::TimeSizeRateMixer(
        unsigned int tree_index) : TimeOperatorInterface<ScaleOperator>(tree_index) {
    this->op_ = ScaleOperator();
}

TimeSizeRateMixer::TimeSizeRateMixer(
        double weight) : TimeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator();
}

TimeSizeRateMixer::TimeSizeRateMixer(
        unsigned int tree_index,
        double weight) : TimeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator();
}

TimeSizeRateMixer::TimeSizeRateMixer(
        double weight,
        double scale) : TimeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator(scale);
}

TimeSizeRateMixer::TimeSizeRateMixer(
        unsigned int tree_index,
        double weight,
        double scale) : TimeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator(scale);
}

void TimeSizeRateMixer::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);

    // Perform sweep of univariate moves
    if (this->tree_index_ < 0) {
        for (std::shared_ptr<OperatorInterface> time_op : comparisons->get_time_operators()) {
            time_op->operate(rng, comparisons, nthreads);
        }
        for (std::shared_ptr<OperatorInterface> tree_op : comparisons->get_tree_operators()) {
            tree_op->operate(rng, comparisons, nthreads);
        }
        return;
    }
    std::vector< std::shared_ptr<OperatorInterface> > time_operators = comparisons->get_time_operators(this->tree_index_);
    if (time_operators.size() < 1) {
        time_operators = comparisons->get_time_operators();
    }
    for (std::shared_ptr<OperatorInterface> time_op : time_operators) {
        time_op->operate(rng, comparisons, nthreads);
    }
    for (std::shared_ptr<OperatorInterface> tree_op : comparisons->get_tree_operators(this->tree_index_)) {
        tree_op->operate(rng, comparisons, nthreads);
    }
}

double TimeSizeRateMixer::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int index) {
    if (this->tree_index_ < 0) {
        return this->propose_by_height(rng, comparisons, index);
    }
    return this->propose_by_tree(rng, comparisons, index);
}

double TimeSizeRateMixer::propose_by_height(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int height_index) {
    double multiplier = this->op_.get_move_amount(rng);

    double new_height = comparisons->get_height(height_index) * multiplier;
    comparisons->set_height(height_index, new_height);

    int ndimensions = 1; // for height scaled above
    int number_of_free_parameters_scaled = 0;
    int number_of_free_parameters_inverse_scaled = 0;
    std::vector<unsigned int> tree_indices = comparisons->get_indices_of_mapped_trees(height_index);
    for (auto tree_idx : tree_indices) {
        std::shared_ptr<PopulationTree> tree = comparisons->get_tree(tree_idx);
        if (! tree->population_sizes_are_fixed()) {
            if (tree->population_sizes_are_constrained()) {
                number_of_free_parameters_scaled += tree->scale_root_population_size(multiplier);
            }
            else {
                number_of_free_parameters_inverse_scaled += tree->scale_root_population_size(1.0/multiplier);
                unsigned int nleaves = tree->get_leaf_node_count();
                for (unsigned int i = 0; i < nleaves; ++i) { 
                    if (! tree->get_child_population_size_parameter(i)->is_fixed()) {
                        tree->set_child_population_size(i, 
                                tree->get_child_population_size(i) * multiplier);
                        ++number_of_free_parameters_scaled;
                    }
                }
            }
        }
        if (! tree->mutation_rate_is_fixed()) {
            tree->set_mutation_rate(tree->get_mutation_rate() * (1.0/multiplier));
            ++number_of_free_parameters_inverse_scaled;
        }
    }
    ndimensions += (number_of_free_parameters_scaled - number_of_free_parameters_inverse_scaled);
    return std::log(multiplier) * ndimensions;
}

double TimeSizeRateMixer::propose_by_tree(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int tree_index) {
    double multiplier = this->op_.get_move_amount(rng);

    double new_height = comparisons->get_height_of_tree(tree_index) * multiplier;
    comparisons->set_height_of_tree(tree_index, new_height);

    int ndimensions = 1; // for height scaled above
    int number_of_free_parameters_scaled = 0;
    int number_of_free_parameters_inverse_scaled = 0;
    std::shared_ptr<PopulationTree> tree = comparisons->get_tree(tree_index);
    if (! tree->population_sizes_are_fixed()) {
        if (tree->population_sizes_are_constrained()) {
            number_of_free_parameters_scaled += tree->scale_root_population_size(multiplier);
        }
        else {
            number_of_free_parameters_inverse_scaled += tree->scale_root_population_size(1.0/multiplier);
            unsigned int nleaves = tree->get_leaf_node_count();
            for (unsigned int i = 0; i < nleaves; ++i) { 
                if (! tree->get_child_population_size_parameter(i)->is_fixed()) {
                    tree->set_child_population_size(i, 
                            tree->get_child_population_size(i) * multiplier);
                    ++number_of_free_parameters_scaled;
                }
            }
        }
    }
    if (! tree->mutation_rate_is_fixed()) {
        tree->set_mutation_rate(tree->get_mutation_rate() * (1.0/multiplier));
        ++number_of_free_parameters_inverse_scaled;
    }
    ndimensions += (number_of_free_parameters_scaled - number_of_free_parameters_inverse_scaled);
    return std::log(multiplier) * ndimensions;
}

std::string TimeSizeRateMixer::get_name() const {
    std::string name = "TimeSizeRateMixer";
    if (this->tree_index_ > -1) {
        name += std::to_string(this->tree_index_);
    }
    return name;
}

std::string TimeSizeRateMixer::target_parameter() const {
    return "node heights, population sizes, and mutation rates";
}

std::string TimeSizeRateMixer::to_string(const OperatorSchedule& os) const {
    std::ostringstream ss;
    ss << this->get_name() << "\t" 
       << this->get_number_accepted() << "\t"
       << this->get_number_rejected() << "\t"
       << this->get_weight() << "\t";

    if (os.get_total_weight() > 0.0) {
        ss << this->get_weight() / os.get_total_weight() << "\t";
    }
    else {
        ss << "nan\t";
    }

    double tuning = this->get_coercable_parameter_value();
    if (std::isnan(tuning)) {
        ss << "none\t";
    }
    else {
        ss << tuning << "\t";
    }
    ss << "\n";
    return ss.str();
}


//////////////////////////////////////////////////////////////////////////////
// TimeMeanSizeRateMixer methods
//////////////////////////////////////////////////////////////////////////////

TimeMeanSizeRateMixer::TimeMeanSizeRateMixer(
        ) : TimeOperatorInterface<ScaleOperator>() {
    this->op_ = ScaleOperator();
}

TimeMeanSizeRateMixer::TimeMeanSizeRateMixer(
        unsigned int tree_index) : TimeOperatorInterface<ScaleOperator>(tree_index) {
    this->op_ = ScaleOperator();
}

TimeMeanSizeRateMixer::TimeMeanSizeRateMixer(
        double weight) : TimeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator();
}

TimeMeanSizeRateMixer::TimeMeanSizeRateMixer(
        unsigned int tree_index,
        double weight) : TimeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator();
}

TimeMeanSizeRateMixer::TimeMeanSizeRateMixer(
        double weight,
        double scale) : TimeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator(scale);
}

TimeMeanSizeRateMixer::TimeMeanSizeRateMixer(
        unsigned int tree_index,
        double weight,
        double scale) : TimeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator(scale);
}

void TimeMeanSizeRateMixer::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);

    // Perform sweep of univariate moves
    if (this->tree_index_ < 0) {
        for (std::shared_ptr<OperatorInterface> time_op : comparisons->get_time_operators()) {
            time_op->operate(rng, comparisons, nthreads);
        }
        for (std::shared_ptr<OperatorInterface> tree_op : comparisons->get_tree_operators()) {
            tree_op->operate(rng, comparisons, nthreads);
        }
        return;
    }
    std::vector< std::shared_ptr<OperatorInterface> > time_operators = comparisons->get_time_operators(this->tree_index_);
    if (time_operators.size() < 1) {
        time_operators = comparisons->get_time_operators();
    }
    for (std::shared_ptr<OperatorInterface> time_op : time_operators) {
        time_op->operate(rng, comparisons, nthreads);
    }
    for (std::shared_ptr<OperatorInterface> tree_op : comparisons->get_tree_operators(this->tree_index_)) {
        tree_op->operate(rng, comparisons, nthreads);
    }
}

double TimeMeanSizeRateMixer::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int index) {
    if (this->tree_index_ < 0) {
        return this->propose_by_height(rng, comparisons, index);
    }
    return this->propose_by_tree(rng, comparisons, index);
}

double TimeMeanSizeRateMixer::propose_by_height(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int height_index) {
    double multiplier = this->op_.get_move_amount(rng);

    double new_height = comparisons->get_height(height_index) * multiplier;
    comparisons->set_height(height_index, new_height);

    int ndimensions = 1; // for height scaled above
    std::vector<unsigned int> tree_indices = comparisons->get_indices_of_mapped_trees(height_index);
    for (auto tree_idx : tree_indices) {
        int number_of_free_parameters_scaled = 0;
        int number_of_free_parameters_inverse_scaled = 0;
        std::shared_ptr<PopulationTree> tree = comparisons->get_tree(tree_idx);
        if (! (tree->mean_population_size_is_fixed() && tree->population_size_multipliers_are_fixed())) {
            if (tree->population_size_multipliers_are_fixed()) {
                tree->scale_all_population_sizes(multiplier);
                ++number_of_free_parameters_scaled;
            }
            else if (tree->mean_population_size_is_fixed()) {
                double mean_size = tree->get_mean_population_size();
                bool decrease_root = false;
                if (multiplier > 1.0) {
                    decrease_root = true;
                }
                // This is a hack to keep this move symmetric on average
                int magnitude_toggle = rng.uniform_int(0, 1);
                double delta = std::abs(mean_size - (mean_size * multiplier));
                if (magnitude_toggle == 0) {
                    delta = std::abs(mean_size - (mean_size * (1.0/multiplier)));
                }
                if (decrease_root) {
                    delta *= -1.0;
                }
                std::vector<double> sizes = tree->get_population_sizes();
                ECOEVOLITY_ASSERT(sizes.size() > 1);
                unsigned int num_nonroot_branches = sizes.size() - 1;
                sizes.at(num_nonroot_branches) += delta;
                if (sizes.at(num_nonroot_branches) <= 0.0) {
                    return -std::numeric_limits<double>::infinity();
                }
                /* int leaf_idx = rng.uniform_int(0, tree->get_leaf_node_count() - 1); */
                /* sizes.at(leaf_idx) -= delta; */
                /* if (sizes.at(leaf_idx) < 0.0) { */
                /*     return -std::numeric_limits<double>::infinity(); */
                /* } */
                for (unsigned int i = 0; i < num_nonroot_branches; ++i) {
                    sizes.at(i) -= (delta/num_nonroot_branches);
                    if (sizes.at(i) <= 0.0) {
                        return -std::numeric_limits<double>::infinity();
                    }
                }
                // ECOEVOLITY_ASSERT_APPROX_EQUAL(tree->get_mean_population_size(), mean_size);
                tree->set_population_sizes(sizes);
            }
            else {
                tree->scale_all_population_sizes(1.0/multiplier);
                ++number_of_free_parameters_inverse_scaled;

                tree->scale_root_population_size(1.0/multiplier);
                ++number_of_free_parameters_inverse_scaled;
                unsigned int nleaves = tree->get_leaf_node_count();
                for (unsigned int i = 0; i < nleaves; ++i) { 
                    tree->set_child_population_size(i, 
                            tree->get_child_population_size(i) * multiplier);
                    ++number_of_free_parameters_scaled;
                }
            }
        }
        if (! tree->mutation_rate_is_fixed()) {
            tree->set_mutation_rate(tree->get_mutation_rate() * (1.0/multiplier));
            ++number_of_free_parameters_inverse_scaled;
        }
        if ((number_of_free_parameters_inverse_scaled < 1) ||
                ((number_of_free_parameters_scaled - number_of_free_parameters_inverse_scaled) < 2)) {
            ndimensions += (number_of_free_parameters_scaled - number_of_free_parameters_inverse_scaled);
        }
        else {
            ndimensions += (number_of_free_parameters_scaled - number_of_free_parameters_inverse_scaled - 2);
        }
    }
    return std::log(multiplier) * ndimensions;
}

double TimeMeanSizeRateMixer::propose_by_tree(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int tree_index) {
    double multiplier = this->op_.get_move_amount(rng);

    double new_height = comparisons->get_height_of_tree(tree_index) * multiplier;
    comparisons->set_height_of_tree(tree_index, new_height);

    int ndimensions = 1; // for height scaled above
    int number_of_free_parameters_scaled = 0;
    int number_of_free_parameters_inverse_scaled = 0;
    std::shared_ptr<PopulationTree> tree = comparisons->get_tree(tree_index);
    if (! (tree->mean_population_size_is_fixed() && tree->population_size_multipliers_are_fixed())) {
        if (tree->population_size_multipliers_are_fixed()) {
            tree->scale_all_population_sizes(multiplier);
            ++number_of_free_parameters_scaled;
        }
        else if (tree->mean_population_size_is_fixed()) {
            double mean_size = tree->get_mean_population_size();
            bool decrease_root = false;
            if (multiplier > 1.0) {
                decrease_root = true;
            }
            // This is a hack to keep this move symmetric on average
            int magnitude_toggle = rng.uniform_int(0, 1);
            double delta = std::abs(mean_size - (mean_size * multiplier));
            if (magnitude_toggle == 0) {
                delta = std::abs(mean_size - (mean_size * (1.0/multiplier)));
            }
            if (decrease_root) {
                delta *= -1.0;
            }
            std::vector<double> sizes = tree->get_population_sizes();
            ECOEVOLITY_ASSERT(sizes.size() > 1);
            unsigned int num_nonroot_branches = sizes.size() - 1;
            sizes.at(num_nonroot_branches) += delta;
            if (sizes.at(num_nonroot_branches) <= 0.0) {
                return -std::numeric_limits<double>::infinity();
            }
            /* int leaf_idx = rng.uniform_int(0, tree->get_leaf_node_count() - 1); */
            /* sizes.at(leaf_idx) -= delta; */
            /* if (sizes.at(leaf_idx) < 0.0) { */
            /*     return -std::numeric_limits<double>::infinity(); */
            /* } */
            for (unsigned int i = 0; i < num_nonroot_branches; ++i) {
                sizes.at(i) -= (delta/num_nonroot_branches);
                if (sizes.at(i) <= 0.0) {
                    return -std::numeric_limits<double>::infinity();
                }
            }
            // ECOEVOLITY_ASSERT_APPROX_EQUAL(tree->get_mean_population_size(), mean_size);
            tree->set_population_sizes(sizes);
        }
        else {
            tree->scale_all_population_sizes(1.0/multiplier);
            ++number_of_free_parameters_inverse_scaled;

            tree->scale_root_population_size(1.0/multiplier);
            ++number_of_free_parameters_inverse_scaled;
            unsigned int nleaves = tree->get_leaf_node_count();
            for (unsigned int i = 0; i < nleaves; ++i) { 
                tree->set_child_population_size(i, 
                        tree->get_child_population_size(i) * multiplier);
                ++number_of_free_parameters_scaled;
            }
        }
    }
    if (! tree->mutation_rate_is_fixed()) {
        tree->set_mutation_rate(tree->get_mutation_rate() * (1.0/multiplier));
        ++number_of_free_parameters_inverse_scaled;
    }
    if ((number_of_free_parameters_inverse_scaled < 1) ||
            ((number_of_free_parameters_scaled - number_of_free_parameters_inverse_scaled) < 2)) {
        ndimensions += (number_of_free_parameters_scaled - number_of_free_parameters_inverse_scaled);
    }
    else {
        ndimensions += (number_of_free_parameters_scaled - number_of_free_parameters_inverse_scaled - 2);
    }
    return std::log(multiplier) * ndimensions;
}

std::string TimeMeanSizeRateMixer::get_name() const {
    std::string name = "TimeMeanSizeRateMixer";
    if (this->tree_index_ > -1) {
        name += std::to_string(this->tree_index_);
    }
    return name;
}

std::string TimeMeanSizeRateMixer::target_parameter() const {
    return "node heights, population sizes, and mutation rates";
}

std::string TimeMeanSizeRateMixer::to_string(const OperatorSchedule& os) const {
    std::ostringstream ss;
    ss << this->get_name() << "\t" 
       << this->get_number_accepted() << "\t"
       << this->get_number_rejected() << "\t"
       << this->get_weight() << "\t";

    if (os.get_total_weight() > 0.0) {
        ss << this->get_weight() / os.get_total_weight() << "\t";
    }
    else {
        ss << "nan\t";
    }

    double tuning = this->get_coercable_parameter_value();
    if (std::isnan(tuning)) {
        ss << "none\t";
    }
    else {
        ss << tuning << "\t";
    }
    ss << "\n";
    return ss.str();
}


//////////////////////////////////////////////////////////////////////////////
// TimeSizeRateScaler methods
//////////////////////////////////////////////////////////////////////////////

TimeSizeRateScaler::TimeSizeRateScaler(
        ) : TimeOperatorInterface<ScaleOperator>() {
    this->op_ = ScaleOperator();
}

TimeSizeRateScaler::TimeSizeRateScaler(
        unsigned int tree_index) : TimeOperatorInterface<ScaleOperator>(tree_index) {
    this->op_ = ScaleOperator();
}

TimeSizeRateScaler::TimeSizeRateScaler(
        double weight) : TimeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator();
}

TimeSizeRateScaler::TimeSizeRateScaler(
        unsigned int tree_index,
        double weight) : TimeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator();
}

TimeSizeRateScaler::TimeSizeRateScaler(
        double weight,
        double scale) : TimeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator(scale);
}

TimeSizeRateScaler::TimeSizeRateScaler(
        unsigned int tree_index,
        double weight,
        double scale) : TimeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator(scale);
}


void TimeSizeRateScaler::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);

    // Perform sweep of univariate moves
    if (this->tree_index_ < 0) {
        for (std::shared_ptr<OperatorInterface> time_op : comparisons->get_time_operators()) {
            time_op->operate(rng, comparisons, nthreads);
        }
        for (std::shared_ptr<OperatorInterface> tree_op : comparisons->get_tree_operators()) {
            tree_op->operate(rng, comparisons, nthreads);
        }
        return;
    }
    std::vector< std::shared_ptr<OperatorInterface> > time_operators = comparisons->get_time_operators(this->tree_index_);
    if (time_operators.size() < 1) {
        time_operators = comparisons->get_time_operators();
    }
    for (std::shared_ptr<OperatorInterface> time_op : time_operators) {
        time_op->operate(rng, comparisons, nthreads);
    }
    for (std::shared_ptr<OperatorInterface> tree_op : comparisons->get_tree_operators(this->tree_index_)) {
        tree_op->operate(rng, comparisons, nthreads);
    }
}

double TimeSizeRateScaler::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int index) {
    if (this->tree_index_ < 0) {
        return this->propose_by_height(rng, comparisons, index);
    }
    return this->propose_by_tree(rng, comparisons, index);
}

double TimeSizeRateScaler::propose_by_height(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int height_index) {
    double multiplier = this->op_.get_move_amount(rng);

    double new_height = comparisons->get_height(height_index) * multiplier;
    comparisons->set_height(height_index, new_height);

    std::vector<unsigned int> tree_indices = comparisons->get_indices_of_mapped_trees(height_index);
    int ndimensions = 1; // for height scaled above
    for (auto tree_idx : tree_indices) {
        std::shared_ptr<PopulationTree> tree = comparisons->get_tree(tree_idx);
        unsigned int nparameters_scaled = tree->scale_all_population_sizes(multiplier);
        ndimensions += nparameters_scaled;
        if (! tree->mutation_rate_is_fixed()) {
            tree->set_mutation_rate(tree->get_mutation_rate() * (1.0/multiplier));
            --ndimensions;
        }
    }

    return std::log(multiplier) * ndimensions;
}

double TimeSizeRateScaler::propose_by_tree(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int tree_index) {
    double multiplier = this->op_.get_move_amount(rng);

    double new_height = comparisons->get_height_of_tree(tree_index) * multiplier;
    comparisons->set_height_of_tree(tree_index, new_height);

    int ndimensions = 1; // for height scaled above
    std::shared_ptr<PopulationTree> tree = comparisons->get_tree(tree_index);
    unsigned int nparameters_scaled = tree->scale_all_population_sizes(multiplier);
    ndimensions += nparameters_scaled;
    if (! tree->mutation_rate_is_fixed()) {
        tree->set_mutation_rate(tree->get_mutation_rate() * (1.0/multiplier));
        --ndimensions;
    }

    return std::log(multiplier) * ndimensions;
}

std::string TimeSizeRateScaler::get_name() const {
    std::string name = "TimeSizeRateScaler";
    if (this->tree_index_ > -1) {
        name += std::to_string(this->tree_index_);
    }
    return name;
}

std::string TimeSizeRateScaler::target_parameter() const {
    return "node heights, population sizes, and mutation rates";
}

std::string TimeSizeRateScaler::to_string(const OperatorSchedule& os) const {
    std::ostringstream ss;
    ss << this->get_name() << "\t" 
       << this->get_number_accepted() << "\t"
       << this->get_number_rejected() << "\t"
       << this->get_weight() << "\t";

    if (os.get_total_weight() > 0.0) {
        ss << this->get_weight() / os.get_total_weight() << "\t";
    }
    else {
        ss << "nan\t";
    }

    double tuning = this->get_coercable_parameter_value();
    if (std::isnan(tuning)) {
        ss << "none\t";
    }
    else {
        ss << tuning << "\t";
    }
    ss << "\n";
    return ss.str();
}


//////////////////////////////////////////////////////////////////////////////
// TimeMeanSizeRateScaler methods
//////////////////////////////////////////////////////////////////////////////

TimeMeanSizeRateScaler::TimeMeanSizeRateScaler(
        ) : TimeOperatorInterface<ScaleOperator>() {
    this->op_ = ScaleOperator();
}

TimeMeanSizeRateScaler::TimeMeanSizeRateScaler(
        unsigned int tree_index) : TimeOperatorInterface<ScaleOperator>(tree_index) {
    this->op_ = ScaleOperator();
}

TimeMeanSizeRateScaler::TimeMeanSizeRateScaler(
        double weight) : TimeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator();
}

TimeMeanSizeRateScaler::TimeMeanSizeRateScaler(
        unsigned int tree_index,
        double weight) : TimeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator();
}

TimeMeanSizeRateScaler::TimeMeanSizeRateScaler(
        double weight,
        double scale) : TimeOperatorInterface<ScaleOperator>(weight) {
    this->op_ = ScaleOperator(scale);
}

TimeMeanSizeRateScaler::TimeMeanSizeRateScaler(
        unsigned int tree_index,
        double weight,
        double scale) : TimeOperatorInterface<ScaleOperator>(tree_index, weight) {
    this->op_ = ScaleOperator(scale);
}

void TimeMeanSizeRateScaler::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);

    // Perform sweep of univariate moves
    if (this->tree_index_ < 0) {
        for (std::shared_ptr<OperatorInterface> time_op : comparisons->get_time_operators()) {
            time_op->operate(rng, comparisons, nthreads);
        }
        for (std::shared_ptr<OperatorInterface> tree_op : comparisons->get_tree_operators()) {
            tree_op->operate(rng, comparisons, nthreads);
        }
        return;
    }
    std::vector< std::shared_ptr<OperatorInterface> > time_operators = comparisons->get_time_operators(this->tree_index_);
    if (time_operators.size() < 1) {
        time_operators = comparisons->get_time_operators();
    }
    for (std::shared_ptr<OperatorInterface> time_op : time_operators) {
        time_op->operate(rng, comparisons, nthreads);
    }
    for (std::shared_ptr<OperatorInterface> tree_op : comparisons->get_tree_operators(this->tree_index_)) {
        tree_op->operate(rng, comparisons, nthreads);
    }
}

double TimeMeanSizeRateScaler::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int index) {
    if (this->tree_index_ < 0) {
        return this->propose_by_height(rng, comparisons, index);
    }
    return this->propose_by_tree(rng, comparisons, index);
}

double TimeMeanSizeRateScaler::propose_by_height(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int height_index) {
    double multiplier = this->op_.get_move_amount(rng);

    double new_height = comparisons->get_height(height_index) * multiplier;
    comparisons->set_height(height_index, new_height);

    std::vector<unsigned int> tree_indices = comparisons->get_indices_of_mapped_trees(height_index);
    int ndimensions = 1; // for height scaled above
    for (auto tree_idx : tree_indices) {
        std::shared_ptr<PopulationTree> tree = comparisons->get_tree(tree_idx);
        if (! tree->mean_population_size_is_fixed()) {
            tree->set_mean_population_size(tree->get_mean_population_size() * multiplier);
            ++ndimensions;
        }
        if (! tree->mutation_rate_is_fixed()) {
            tree->set_mutation_rate(tree->get_mutation_rate() * (1.0/multiplier));
            --ndimensions;
        }
    }

    return std::log(multiplier) * ndimensions;
}

double TimeMeanSizeRateScaler::propose_by_tree(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int tree_index) {
    double multiplier = this->op_.get_move_amount(rng);

    double new_height = comparisons->get_height_of_tree(tree_index) * multiplier;
    comparisons->set_height_of_tree(tree_index, new_height);

    int ndimensions = 1; // for height scaled above
    std::shared_ptr<PopulationTree> tree = comparisons->get_tree(tree_index);
    if (! tree->mean_population_size_is_fixed()) {
        tree->set_mean_population_size(tree->get_mean_population_size() * multiplier);
        ++ndimensions;
    }
    if (! tree->mutation_rate_is_fixed()) {
        tree->set_mutation_rate(tree->get_mutation_rate() * (1.0/multiplier));
        --ndimensions;
    }

    return std::log(multiplier) * ndimensions;
}

std::string TimeMeanSizeRateScaler::get_name() const {
    std::string name = "TimeMeanSizeRateScaler";
    if (this->tree_index_ > -1) {
        name += std::to_string(this->tree_index_);
    }
    return name;
}

std::string TimeMeanSizeRateScaler::target_parameter() const {
    return "node heights, population sizes, and mutation rates";
}

std::string TimeMeanSizeRateScaler::to_string(const OperatorSchedule& os) const {
    std::ostringstream ss;
    ss << this->get_name() << "\t" 
       << this->get_number_accepted() << "\t"
       << this->get_number_rejected() << "\t"
       << this->get_weight() << "\t";

    if (os.get_total_weight() > 0.0) {
        ss << this->get_weight() / os.get_total_weight() << "\t";
    }
    else {
        ss << "nan\t";
    }

    double tuning = this->get_coercable_parameter_value();
    if (std::isnan(tuning)) {
        ss << "none\t";
    }
    else {
        ss << tuning << "\t";
    }
    ss << "\n";
    return ss.str();
}


//////////////////////////////////////////////////////////////////////////////
// DirichletProcessGibbsSampler methods
//////////////////////////////////////////////////////////////////////////////

DirichletProcessGibbsSampler::DirichletProcessGibbsSampler(
        double weight,
        unsigned int number_of_auxiliary_categories) : CollectionOperatorInterface<Operator>(weight) {
    this->set_number_of_auxiliary_categories(number_of_auxiliary_categories);
}

void DirichletProcessGibbsSampler::set_number_of_auxiliary_categories(
        unsigned int n) {
    this->number_of_auxiliary_categories_ = n;
}
unsigned int DirichletProcessGibbsSampler::get_number_of_auxiliary_categories() const {
    return this->number_of_auxiliary_categories_;
}

std::string DirichletProcessGibbsSampler::get_name() const {
    return "DirichletProcessGibbsSampler";
}

std::string DirichletProcessGibbsSampler::target_parameter() const {
    return "model";
}

std::string DirichletProcessGibbsSampler::to_string(const OperatorSchedule& os) const {
    std::ostringstream ss;
    ss << this->get_name() << "\t" 
       << this->get_number_accepted() << "\t"
       << this->get_number_rejected() << "\t"
       << this->get_weight() << "\t";

    if (os.get_total_weight() > 0.0) {
        ss << this->get_weight() / os.get_total_weight() << "\t";
    }
    else {
        ss << "nan\t";
    }

    double tuning = this->get_coercable_parameter_value();
    if (std::isnan(tuning)) {
        ss << "none\t";
    }
    else {
        ss << tuning << "\t";
    }
    ss << "\n";
    return ss.str();
}

void DirichletProcessGibbsSampler::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    this->perform_collection_move(rng, comparisons, nthreads);

    // Perform sweep of univariate time moves
    for (std::shared_ptr<OperatorInterface> time_op : comparisons->get_time_operators()) {
        time_op->operate(rng, comparisons, nthreads);
    }
}

void DirichletProcessGibbsSampler::perform_collection_move(
        RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    // this->call_store_methods(comparisons);

    // comparisons->make_trees_dirty();
    // comparisons->compute_log_likelihood_and_prior(true);
    this->propose(rng, comparisons, nthreads);
    // comparisons->make_trees_dirty();
    // comparisons->compute_log_likelihood_and_prior(true);

    // Likelihoods are clean, but update priors
    comparisons->compute_log_likelihood_and_prior(false);

    this->accept(comparisons->get_operator_schedule());
    comparisons->make_trees_clean();
}

double DirichletProcessGibbsSampler::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {

    const double ln_concentration_over_num_aux = std::log(
            comparisons->get_concentration() /
            this->get_number_of_auxiliary_categories());

    for (unsigned int tree_idx = 0;
            tree_idx < comparisons->get_number_of_trees();
            ++tree_idx) {
        std::shared_ptr<PopulationTree> tree = comparisons->get_tree(tree_idx);
        const unsigned int original_height_index = comparisons->get_height_index(tree_idx);
        const double original_height = tree->get_root_height();
        const double original_likelihood = tree->get_log_likelihood_value();

        std::vector<unsigned int> other_height_indices = comparisons->get_other_height_indices(tree_idx);
        std::vector<double> ln_category_likelihoods;
        std::vector<double> ln_tree_likelihoods;
        unsigned int number_of_aux_categories = this->get_number_of_auxiliary_categories();
        bool tree_in_singleton_category = false;
        if (other_height_indices.size() < comparisons->get_number_of_events()) {
            // Tree is in a singleton category. Need to consider singleton
            // category as one of the auxillary categories.
            tree_in_singleton_category = true;
            other_height_indices.push_back(original_height_index);
            --number_of_aux_categories;
        }
        ln_category_likelihoods.reserve(other_height_indices.size() +
                number_of_aux_categories);
        ln_tree_likelihoods.reserve(other_height_indices.size() +
                number_of_aux_categories);

        for (auto height_idx : other_height_indices) {
            unsigned int number_of_elements = comparisons->get_number_of_trees_mapped_to_height(height_idx);
            if (height_idx == original_height_index) {
                if (tree_in_singleton_category) {
                    // Considering singleton category as one of the auxillary categories
                    ln_category_likelihoods.push_back(original_likelihood + ln_concentration_over_num_aux);
                    ln_tree_likelihoods.push_back(original_likelihood);
                }
                else {
                    --number_of_elements;
                    ln_category_likelihoods.push_back(original_likelihood + std::log(number_of_elements));
                    ln_tree_likelihoods.push_back(original_likelihood);
                }
                continue;
            }
            tree->set_root_height(comparisons->get_height(height_idx));
            double lnl = tree->compute_log_likelihood(nthreads);
            ln_category_likelihoods.push_back(lnl + std::log(number_of_elements));
            ln_tree_likelihoods.push_back(lnl);
        }

        std::vector<double> auxiliary_heights;
        auxiliary_heights.reserve(number_of_aux_categories);
        for (unsigned int i = 0; i < number_of_aux_categories; ++i) {
            double fresh_height = comparisons->get_draw_from_node_height_prior(rng);
            auxiliary_heights.push_back(fresh_height);
            tree->set_root_height(fresh_height);
            double lnl = tree->compute_log_likelihood(nthreads);
            ln_category_likelihoods.push_back(lnl + ln_concentration_over_num_aux);
            ln_tree_likelihoods.push_back(lnl);
        }

        // restore height associated with this tree
        tree->set_root_height(original_height);

        std::vector<double> category_probs(ln_category_likelihoods);
        normalize_log_likelihoods(category_probs);
        unsigned int prob_index = rng.weighted_index(category_probs);
        if (prob_index < other_height_indices.size()) {
            comparisons->remap_tree(
                    tree_idx,
                    other_height_indices.at(prob_index),
                    ln_tree_likelihoods.at(prob_index));
        }
        else {
            comparisons->map_tree_to_new_height(
                    tree_idx,
                    auxiliary_heights.at(prob_index -
                            other_height_indices.size()),
                    ln_tree_likelihoods.at(prob_index));
        }
    }
    // Always accept, so returning inf
    return std::numeric_limits<double>::infinity();
}


//////////////////////////////////////////////////////////////////////////////
// ReversibleJumpSampler methods
//////////////////////////////////////////////////////////////////////////////

std::string ReversibleJumpSampler::get_name() const {
    return "ReversibleJumpSampler";
}

std::string ReversibleJumpSampler::target_parameter() const {
    return "model";
}

void ReversibleJumpSampler::call_store_methods(
        BaseComparisonPopulationTreeCollection * comparisons) const {
    comparisons->store_state();
    comparisons->store_model_state();
}

void ReversibleJumpSampler::call_restore_methods(
        BaseComparisonPopulationTreeCollection * comparisons) const {
    comparisons->restore_model_state();
}

OperatorInterface::OperatorTypeEnum ReversibleJumpSampler::get_type() const {
    return OperatorInterface::OperatorTypeEnum::rj_operator;
}

void ReversibleJumpSampler::write_split_probabilities(std::ostream& out) const {
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

std::string ReversibleJumpSampler::to_string(const OperatorSchedule& os) const {
    std::ostringstream ss;
    ss << this->get_name() << "\t" 
       << this->get_number_accepted() << "\t"
       << this->get_number_rejected() << "\t"
       << this->get_weight() << "\t";

    if (os.get_total_weight() > 0.0) {
        ss << this->get_weight() / os.get_total_weight() << "\t";
    }
    else {
        ss << "nan\t";
    }

    double tuning = this->get_coercable_parameter_value();
    if (std::isnan(tuning)) {
        ss << "none\t";
    }
    else {
        ss << tuning << "\t";
    }
    ss << "\n";
    ss << this->time_scaler_.to_string(os);
    return ss.str();
}

void ReversibleJumpSampler::operate(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    for (unsigned int i = 0; i < comparisons->get_number_of_trees(); ++i) {
        this->perform_collection_move(rng, comparisons, nthreads);

        // Perform sweep of univariate time moves
        this->time_scaler_.operate(rng, comparisons, nthreads);
        // for (std::shared_ptr<OperatorInterface> time_op : comparisons->get_time_operators()) {
        //     time_op->operate(rng, comparisons, nthreads);
        // }
    }
}

double ReversibleJumpSampler::propose(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons,
        unsigned int nthreads) {
    return this->propose_jump_to_gap(rng, comparisons);
}

double ReversibleJumpSampler::propose_jump_to_prior(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons) {
    throw EcoevolityNotImplementedError(
            "The reversible 'jump to prior' move is currently not "
            "implemented");
    // TODO:
    // This currently only works for an exponential distribution on node
    // heights (the jacobian needs to be solved for a gamma distribution), and
    // it seems to have some bad corner cases (e.g., with only two comparisons,
    // this gets accepted every time).
    // It needs work/debugging. For now, we will just use the 'jump to gap'
    // move, which might perform better anyway, especially if the prior on
    // heights is misspecified.
    const unsigned int nnodes = comparisons->get_number_of_trees();
    const unsigned int nevents = comparisons->get_number_of_events();
    const bool in_general_state_before = (nnodes == nevents);
    const bool in_shared_state_before = (nevents == 1);
    const bool split_event = ((! in_general_state_before) &&
            (in_shared_state_before || (rng.uniform_real() < 0.5)));
    double mean_height = comparisons->get_node_height_prior_mean();
    if (split_event) {
        std::vector<unsigned int> shared_indices =
                comparisons->get_shared_event_indices();
        unsigned int i = rng.uniform_int(0, shared_indices.size() - 1);
        unsigned int event_index = shared_indices.at(i);
        double new_height = comparisons->get_draw_from_node_height_prior(rng);

        std::vector<unsigned int> tree_indices = comparisons->get_indices_of_mapped_trees(event_index);
        unsigned int num_mapped_nodes = tree_indices.size();
        const std::vector<double>& split_size_probs = 
                this->get_split_subset_size_probabilities(num_mapped_nodes);

        unsigned int subset_size = rng.weighted_index(split_size_probs) + 1;
        ECOEVOLITY_ASSERT((subset_size > 0) && (subset_size < num_mapped_nodes));

        std::vector<unsigned int> random_indices = rng.random_subset_indices(
                num_mapped_nodes,
                subset_size);
        std::vector<unsigned int> subset_indices;
        subset_indices.reserve(subset_size);
        for (auto const random_idx: random_indices) {
            subset_indices.push_back(tree_indices.at(random_idx));
        }
        comparisons->map_trees_to_new_height(subset_indices, new_height);

        // TODO: check this
        double ln_jacobian = std::log(mean_height) + (new_height / mean_height);

        // The probability of forward split move (just proposed) is the product
        // of the probabilites of
        //   1) choosing the shared event to split
        //          = 1 / number of shared events
        //   2) randomly splitting the subset out of the 'n' nodes sharing the event
        //          = 1 / (2 * stirling2(n, 2))
        //          = 1 / (2^n - 2)
        // So the prob of the forward split move is
        // p(split move) =  1 / (number of shared events * (2^n - 2))
        //                  1 / (number of shared events * 2 * stirling2(n, 2))
        //
        // The probability of the reverse move is the product of the
        // probability of
        //   1) randomly selecting the proposed (split) event from among all events.
        //          = 1 / (number of events before proposal + 1)
        //          = 1 / (number of events after the proposal)
        //   2) randomly selecting the correct target to which to merge
        //          = 1 / (number of events before the proposal)
        //          = 1 / (number of events after the proposal - 1)
        //
        // So, the Hasting ratio for the proposed split move is:
        // p(reverse merge) / p(proposed split) = 
        //                 (number of shared events * 2 * stirling2(n, 2))
        //    ----------------------------------------------------------------------------
        //    (number of events after the proposal * number of events before the proposal)
        double ln_hastings =
                this->ln_number_of_possible_splits_.at(num_mapped_nodes) +
                std::log(shared_indices.size()) -
                (std::log(nevents + 1) + std::log(nevents));

        const bool in_general_state_after = (comparisons->get_number_of_trees() ==
                comparisons->get_number_of_events());
        // Account for probability of choosing to split
        // This is 1.0 (or zero on log scale) except for these two corner
        // cases:
        if (in_shared_state_before && (! in_general_state_after)) {
            ln_hastings -= std::log(2.0);
        }
        else if (in_general_state_after && (! in_shared_state_before)) {
            ln_hastings += std::log(2.0);
        }
        return ln_hastings + ln_jacobian;
    }
    // Merge move
    std::vector<unsigned int> height_indices =  rng.random_subset_indices(nevents, 2);
    unsigned int move_height_index = height_indices.at(0);
    unsigned int target_height_index = height_indices.at(1);
    if (rng.uniform_real() < 0.5) {
        move_height_index = height_indices.at(1);
        target_height_index = height_indices.at(0);
    }
    double removed_height = comparisons->get_height(move_height_index);
    unsigned int new_merged_event_index = comparisons->merge_height(move_height_index, target_height_index);
    unsigned int nnodes_in_merged_event = comparisons->get_number_of_trees_mapped_to_height(new_merged_event_index);
    // Don't need the returned probability vector, but need to make sure we
    // update the stored number of ways to make the reverse split of this
    // number of nodes.
    this->get_split_subset_size_probabilities(nnodes_in_merged_event);

    // TODO: check this
    double ln_jacobian = -removed_height / mean_height - std::log(mean_height);

    // The probability of the forward merge move is the product of the
    // probability of
    //   1) randomly selecting the height to move from among all events
    //          = 1 / (number of events before proposal)
    //   2) randomly selecting the target to merge with
    //          = 1 / (number of events before proposal - 1)
    // p(forward merge) = 1 / (number of events before proposal * (number of events before proposal - 1))
    //
    // The probability of reverse split move is the product of the probabilites of
    //   1) choosing the merged shared event to split
    //          = 1 / number of shared events after the proposal
    //   2) randomly splitting the subset out of the 'n' nodes in the merged
    //   event (the number of nodes in the event after the proposal).
    //          = 1 / (2 * stirling2(n, 2))
    //          = 1 / (2^n - 2)
    // So the prob of the reverse split move is
    // p(reverse split move) =  1 / (number of shared events after the proposal * (2^n - 2))
    //                  1 / (number of shared events after the proposal * 2 * stirling2(n, 2))
    //
    // So, the Hasting ratio for the proposed split move is:
    // p(reverse split) / p(proposed merge) = 
    //     (number of events after the proposal) * (number of events before proposal - 1)
    //     ------------------------------------------------------------------------------
    //          (number of shared events after the proposal * 2 * stirling2(n, 2))
    ECOEVOLITY_ASSERT(comparisons->get_number_of_events() + 1 == nevents);
    unsigned int nshared_after = comparisons->get_shared_event_indices().size();
    double ln_hastings = std::log(nevents) + std::log(nevents - 1);
    ln_hastings -= (
            std::log(nshared_after) +
            this->ln_number_of_possible_splits_.at(nnodes_in_merged_event));

    const bool in_shared_state_after = (comparisons->get_number_of_events() == 1);
    // Account for probability of choosing to merge
    // This is 1.0 (or zero on log scale) except for these two corner
    // cases:
    if (in_general_state_before && (! in_shared_state_after)) {
        ln_hastings -= std::log(2.0);
    }
    else if (in_shared_state_after && (! in_general_state_before)) {
        ln_hastings += std::log(2.0);
    }
    return ln_hastings + ln_jacobian;
}

double ReversibleJumpSampler::propose_jump_to_gap(RandomNumberGenerator& rng,
        BaseComparisonPopulationTreeCollection * comparisons) {
    const unsigned int nnodes = comparisons->get_number_of_trees();
    const unsigned int nevents = comparisons->get_number_of_events();
    const bool in_general_state_before = (nnodes == nevents);
    const bool in_shared_state_before = (nevents == 1);
    const bool split_event = ((! in_general_state_before) &&
            (in_shared_state_before || (rng.uniform_real() < 0.5)));
    if (split_event) {
        std::vector<unsigned int> shared_indices =
                comparisons->get_shared_event_indices();
        unsigned int i = rng.uniform_int(0, shared_indices.size() - 1);
        unsigned int event_index = shared_indices.at(i);
        double event_height = comparisons->get_height(event_index);
        double lower_bound = comparisons->get_nearest_smaller_height(event_index);
        double new_height = rng.uniform_real(lower_bound, event_height);

        std::vector<unsigned int> tree_indices = comparisons->get_indices_of_mapped_trees(event_index);
        unsigned int num_mapped_nodes = tree_indices.size();
        const std::vector<double>& split_size_probs = 
                this->get_split_subset_size_probabilities(num_mapped_nodes);

        unsigned int subset_size = rng.weighted_index(split_size_probs) + 1;
        ECOEVOLITY_ASSERT((subset_size > 0) && (subset_size < num_mapped_nodes));

        std::vector<unsigned int> random_indices = rng.random_subset_indices(
                num_mapped_nodes,
                subset_size);
        std::vector<unsigned int> subset_indices;
        subset_indices.reserve(subset_size);
        for (auto const random_idx: random_indices) {
            subset_indices.push_back(tree_indices.at(random_idx));
        }
        comparisons->map_trees_to_new_height(subset_indices, new_height);

        // TODO: check this
        double ln_jacobian = 0.0;

        double ln_model_prior_ratio = std::log(comparisons->get_concentration());

        // The probability of forward split move (just proposed) is the product of the probabilites of
        //   1) choosing the shared event to split
        //          = 1 / number of shared events
        //   2) randomly splitting the subset out of the 'n' nodes sharing the event
        //          = 1 / (2 * stirling2(n, 2))
        //          = 1 / (2^n - 2)
        //   3) drawing the new height uniformly between the height of the
        //      shared event and the next younger height (or zero).
        //          = 1 / (height - younger neighbor height)
        //          = 1 / d
        // So the prob of the forward split move is
        // p(split move) =  1 / (number of shared events * (2^n - 2) * d)
        //                  1 / (number of shared events * 2 * stirling2(n, 2) * d)
        //
        // The probability of the reverse move is simply the probability of
        // randomly selecting the proposed (split) event from among all but the
        // oldest event.
        // p(reverse merge) = 1 / (number of events before proposal + 1 - 1)
        //                  = 1 / (number of events before the proposal)
        //
        // So, the Hasting ratio for the proposed split move is:
        // p(reverse merge) / p(proposed split) = 
        //     (number of shared events * 2 * stirling2(n, 2) * d)
        //     ---------------------------------------------------
        //          (number of events before the proposal)
        double ln_hastings =
                this->ln_number_of_possible_splits_.at(num_mapped_nodes) +
                std::log(shared_indices.size()) +
                std::log(event_height - lower_bound) -
                std::log(nevents);

        const bool in_general_state_after = (comparisons->get_number_of_trees() ==
                comparisons->get_number_of_events());
        // Account for probability of choosing to split
        // This is 1.0 (or zero on log scale) except for these two corner
        // cases:
        if (in_shared_state_before && (! in_general_state_after)) {
            ln_hastings -= std::log(2.0);
        }
        else if (in_general_state_after && (! in_shared_state_before)) {
            ln_hastings += std::log(2.0);
        }
        return ln_model_prior_ratio + ln_hastings + ln_jacobian;
    }
    // Merge move
    std::vector<unsigned int> candidate_indices =
            comparisons->get_height_indices_sans_largest();
    unsigned int i = rng.uniform_int(0, candidate_indices.size() - 1);
    unsigned int height_index = candidate_indices.at(i);
    unsigned int target_height_index = comparisons->get_nearest_larger_height_index(height_index);
    unsigned int new_merged_event_index = comparisons->merge_height(height_index, target_height_index);
    unsigned int nnodes_in_merged_event = comparisons->get_number_of_trees_mapped_to_height(new_merged_event_index);
    // Don't need the returned probability vector, but need to make sure we
    // update the stored number of ways to make the reverse split of this
    // number of nodes.
    this->get_split_subset_size_probabilities(nnodes_in_merged_event);

    // TODO: check this
    double ln_jacobian = 0.0;

    double ln_model_prior_ratio = std::log(1.0 / comparisons->get_concentration());

    // The probability of the forward merge move is simply the probability of
    // randomly selecting the height to merge from among all but the
    // oldest height.
    // p(reverse merge) = 1 / (number of events before proposal - 1)
    //                  = 1 / (number of events after the proposal)
    //
    // The probability of reverse split move is the product of the probabilites of
    //   1) choosing the merged shared event to split
    //          = 1 / number of shared events after the proposal
    //   2) randomly splitting the subset out of the 'n' nodes in the merged
    //   event (the number of nodes in the event after the proposal).
    //          = 1 / (2 * stirling2(n, 2))
    //          = 1 / (2^n - 2)
    //   3) drawing the original height uniformly between the height of the
    //   merged event and the next younger height (or zero).
    //          = 1 / (merged height - younger neighbor height)
    //          = 1 / d
    // So the prob of the reverse split move is
    // p(split move) =  1 / (number of shared events after the proposal * (2^n - 2) * d)
    //                  1 / (number of shared events after the proposal * 2 * stirling2(n, 2) * d)
    //
    // So, the Hasting ratio for the proposed split move is:
    // p(reverse merge) / p(proposed split) = 
    //                  (number of events after the proposal)
    //     ----------------------------------------------------------------------
    //     (number of shared events after the proposal * 2 * stirling2(n, 2) * d)
    ECOEVOLITY_ASSERT(comparisons->get_number_of_events() + 1 == nevents);
    unsigned int nshared_after = comparisons->get_shared_event_indices().size();
    double new_merged_event_height = comparisons->get_height(new_merged_event_index);
    double lower_bound = comparisons->get_nearest_smaller_height(new_merged_event_index);
    double ln_hastings = std::log(comparisons->get_number_of_events());
    ln_hastings -= (
            std::log(nshared_after) +
            this->ln_number_of_possible_splits_.at(nnodes_in_merged_event) +
            std::log(new_merged_event_height - lower_bound));

    const bool in_shared_state_after = (comparisons->get_number_of_events() == 1);
    // Account for probability of choosing to merge
    // This is 1.0 (or zero on log scale) except for these two corner
    // cases:
    if (in_general_state_before && (! in_shared_state_after)) {
        ln_hastings -= std::log(2.0);
    }
    else if (in_shared_state_after && (! in_general_state_before)) {
        ln_hastings += std::log(2.0);
    }
    return ln_model_prior_ratio + ln_hastings + ln_jacobian;
}

const std::vector<double>& ReversibleJumpSampler::get_split_subset_size_probabilities(
        unsigned int number_of_nodes_in_event) {
    unsigned int n = number_of_nodes_in_event;
    if (this->split_subset_size_probs_.count(n) > 0) {
        return this->split_subset_size_probs_.at(n);
    }
    double ln_n_factorial = std::lgamma(n + 1);
    double ln_twice_stirling_num = std::log(
            stirling2_float(number_of_nodes_in_event, 2)) + std::log(2.0);
    this->ln_number_of_possible_splits_[n] = ln_twice_stirling_num;
    this->split_subset_size_probs_[n];
    this->split_subset_size_probs_.at(n).reserve(n - 1);
    double total = 0.0;
    for (unsigned int k = 1; k <= n - 1; ++k) {
        double ln_n_choose_k = ln_n_factorial - std::lgamma(k + 1) - std::lgamma(n - k + 1);
        double prob_k = std::exp(ln_n_choose_k - ln_twice_stirling_num);
        this->split_subset_size_probs_.at(n).push_back(prob_k);
        total += prob_k;
    }
    ECOEVOLITY_ASSERT_APPROX_EQUAL(total, 1.0);
    return this->split_subset_size_probs_.at(n);
}


//////////////////////////////////////////////////////////////////////////////
// ReversibleJumpWindowOperator methods
//////////////////////////////////////////////////////////////////////////////
// 
// ReversibleJumpWindowOperator::ReversibleJumpWindowOperator(
//         double weight,
//         double window_size) : ReversibleJumpSampler(weight) {
//     this->height_mover_.op_.set_window_size(window_size);
// }
// 
// std::string ReversibleJumpWindowOperator::get_name() const {
//     return "ReversibleJumpWindowOperator";
// }
// 
// std::string ReversibleJumpWindowOperator::target_parameter() const {
//     return "model";
// }
// 
// OperatorInterface::OperatorTypeEnum ReversibleJumpWindowOperator::get_type() const {
//     return OperatorInterface::OperatorTypeEnum::rj_operator;
// }
// 
// std::string ReversibleJumpWindowOperator::to_string(const OperatorSchedule& os) const {
//     std::ostringstream ss;
//     ss << this->get_name() << "\t" 
//        << this->get_number_accepted() << "\t"
//        << this->get_number_rejected() << "\t"
//        << this->get_weight() << "\t";
// 
//     if (os.get_total_weight() > 0.0) {
//         ss << this->get_weight() / os.get_total_weight() << "\t";
//     }
//     else {
//         ss << "nan\t";
//     }
// 
//     double tuning = this->height_mover_.get_coercable_parameter_value();
//     if (std::isnan(tuning)) {
//         ss << "none\t";
//     }
//     else {
//         ss << tuning << "\t";
//     }
//     ss << "\n";
//     ss << this->height_mover_.to_string(os);
//     return ss.str();
// }
// 
// void ReversibleJumpWindowOperator::operate(RandomNumberGenerator& rng,
//         BaseComparisonPopulationTreeCollection * comparisons,
//         unsigned int nthreads) {
//     this->propose_height_moves(rng, comparisons, nthreads);
//     this->perform_collection_move(rng, comparisons, nthreads);
// }
// 
// void ReversibleJumpWindowOperator::propose_height_moves(RandomNumberGenerator& rng,
//         BaseComparisonPopulationTreeCollection * comparisons,
//         unsigned int nthreads) {
//     this->height_mover_.operate(rng, comparisons, nthreads);
// }
// 
// void ReversibleJumpWindowOperator::update_height(RandomNumberGenerator& rng,
//         double& height,
//         double& hastings,
//         double window_size) const {
//     double addend = (rng.uniform_real() * 2 * window_size) - window_size;
//     height += addend;
//     hastings = 0.0;
// }
// 
// double ReversibleJumpWindowOperator::propose(RandomNumberGenerator& rng,
//         BaseComparisonPopulationTreeCollection * comparisons,
//         unsigned int nthreads) {
//     const unsigned int nnodes = comparisons->get_number_of_trees();
//     const unsigned int nevents = comparisons->get_number_of_events();
//     const bool in_general_state_before = (nnodes == nevents);
//     const bool in_shared_state_before = (nevents == 1);
//     const bool split_event = ((! in_general_state_before) &&
//             (in_shared_state_before || (rng.uniform_real() < 0.5)));
//     double window_size = (0.5 * this->height_mover_.op_.get_window_size());
//     // double window_size = this->height_mover_.op_.get_window_size();
//     if (split_event) {
//         std::vector<unsigned int> shared_indices =
//                 comparisons->get_shared_event_indices();
//         unsigned int i = rng.uniform_int(0, shared_indices.size() - 1);
//         unsigned int event_index = shared_indices.at(i);
//         double old_height = comparisons->get_height(event_index);
//         double new_height = old_height;
//         double new_height_ln_hastings = 0.0;
//         this->update_height(rng, new_height, new_height_ln_hastings, window_size);
//         // this->height_mover_.update(rng, new_height, new_height_ln_hastings);
//         if (new_height <= 0.0) {
//             return -std::numeric_limits<double>::infinity();
//         }
// 
//         double lower_bound = 0.0;
//         double upper_bound = 0.0;
//         unsigned int next_proximal_height_index = event_index;
//         if (old_height > new_height) {
//             next_proximal_height_index = comparisons->get_nearest_larger_height_index(event_index, true);
//             lower_bound = old_height;
//             upper_bound = comparisons->get_height(next_proximal_height_index);
//             if ((next_proximal_height_index == event_index) ||
//                ((upper_bound - new_height) > window_size)) {
//                 upper_bound = new_height + window_size;
//             }
//         }
//         else {
//             next_proximal_height_index = comparisons->get_nearest_smaller_height_index(event_index, true);
//             upper_bound = old_height;
//             lower_bound = comparisons->get_height(next_proximal_height_index);
//             if ((next_proximal_height_index == event_index) ||
//                ((new_height - lower_bound) > window_size)) {
//                 lower_bound = new_height - window_size;
//                 // if (lower_bound < 0.0) {
//                 //     lower_bound = 0.0;
//                 // }
//             }
//         }
//         double merge_gap = upper_bound - lower_bound;
// 
//         // std::cout << "Splitting:\n";
//         // std::cout << "old height: " << old_height << " " << comparisons->get_height(event_index) << "\n";
//         // std::cout << "new height: " << new_height << "\n";
//         // std::cout << "next proximal height: " << comparisons->get_height(next_proximal_height_index) << "\n";
//         // std::cout << "next_proximal_height_index: " << next_proximal_height_index << "\n";
//         // std::cout << "window size: " << window_size << "\n";
//         // std::cout << "upper bond: " << upper_bound << "\n";
//         // std::cout << "lower bond: " << lower_bound << "\n";
//         // std::cout << "merge gap: " << merge_gap << "\n";
//         // std::cout << "\n";
// 
//         ECOEVOLITY_ASSERT(merge_gap > 0.0);
// 
//         std::vector<unsigned int> tree_indices = comparisons->get_indices_of_mapped_trees(event_index);
//         unsigned int num_mapped_nodes = tree_indices.size();
//         const std::vector<double>& split_size_probs = 
//                 this->get_split_subset_size_probabilities(num_mapped_nodes);
// 
//         unsigned int subset_size = rng.weighted_index(split_size_probs) + 1;
//         ECOEVOLITY_ASSERT((subset_size > 0) && (subset_size < num_mapped_nodes));
// 
//         std::vector<unsigned int> random_indices = rng.random_subset_indices(
//                 num_mapped_nodes,
//                 subset_size);
//         std::vector<unsigned int> subset_indices;
//         subset_indices.reserve(subset_size);
//         for (auto const random_idx: random_indices) {
//             subset_indices.push_back(tree_indices.at(random_idx));
//         }
//         comparisons->map_trees_to_new_height(subset_indices, new_height);
// 
// 
//         // TODO: check this
//         double ln_jacobian = 0.0;
// 
//         // The probability of forward split move (just proposed) is the product of the probabilites of
//         //   1) choosing the shared event to split
//         //          = 1 / number of shared events
//         //   2) randomly splitting the subset out of the 'n' nodes sharing the event
//         //          = 1 / (2 * stirling2(n, 2))
//         //          = 1 / (2^n - 2)
//         //   3) drawing the new height
//         //          = 1 / 2w
//         //      where w = window size
//         // So the prob of the forward split move is
//         // p(split move) = 1 / (number of shared events * (2^n - 2) * 2w)
//         //               = 1 / (number of shared events * 2 * stirling2(n, 2) * 2w)
//         //
//         // The probability of the reverse move is simply the probability of
//         // randomly selecting the proposed (split) event from among all events.
//         // Times the probability of drawing a height that will lead to a merge
//         // p(selecting proposed event) = 1 / (number of events before proposal + 1)
//         // p(drawing merge height) = width of window that leads to merge / 2w
//         //                           = merge gap / 2w 
//         // p(reverse merge) = merge gap / (2w * (number of events before proposal + 1))
//         //
//         // So, the Hasting ratio for the proposed split move is:
//         // p(reverse merge) / p(proposed split) = 
//         //     (number of shared events * 2 * stirling2(n, 2) * 2w * merge gap)
//         //     -----------------------------------------------------------------
//         //               2w * (number of events before the proposal + 1)
//         //  =
//         //     (number of shared events * 2 * stirling2(n, 2) * merge gap)
//         //     -----------------------------------------------------------------
//         //               (number of events before the proposal + 1)
// 
//         double ln_hastings =
//                 this->ln_number_of_possible_splits_.at(num_mapped_nodes) +
//                 std::log(shared_indices.size()) +
//                 std::log(merge_gap) -
//                 std::log(nevents + 1);
// 
//         // Account for hastings of proposed height
//         ln_hastings += new_height_ln_hastings;
// 
//         const bool in_general_state_after = (comparisons->get_number_of_trees() ==
//                 comparisons->get_number_of_events());
//         // Account for probability of choosing to split
//         // This is 1.0 (or zero on log scale) except for these two corner
//         // cases:
//         if (in_shared_state_before && (! in_general_state_after)) {
//             ln_hastings -= std::log(2.0);
//         }
//         else if (in_general_state_after && (! in_shared_state_before)) {
//             ln_hastings += std::log(2.0);
//         }
// 
//         // double merge_overhang = window_size - std::abs(new_height - old_height);
//         // std::cout << "splitting:\n";
//         // std::cout << "hastings: " << std::exp(ln_hastings) << "\n";
//         // std::cout << "gap: " << merge_gap << "\n";
//         // std::cout << "overhang: " << merge_overhang << "\n";
//         // std::cout << "\n";
// 
//         return ln_hastings + ln_jacobian;
//     }
//     // Merge move
//     unsigned int height_index = rng.uniform_int(0, comparisons->get_number_of_events() - 1);
//     double old_height = comparisons->get_height(height_index);
//     double new_height = old_height;
//     double new_height_ln_hastings = 0.0;
//     this->update_height(rng, new_height, new_height_ln_hastings, window_size);
//     // this->height_mover_.update(rng, new_height, new_height_ln_hastings);
//     if (new_height <= 0.0) {
//         return -std::numeric_limits<double>::infinity();
//     }
// 
//     unsigned int target_height_index = comparisons->get_distal_height_index_within_move(height_index, (new_height - old_height));
//     if (target_height_index == height_index) {
//         comparisons->get_height_parameter(height_index)->set_value(new_height);
//         return new_height_ln_hastings;
//     }
// 
//     double lower_bound = 0.0;
//     double upper_bound = 0.0;
//     unsigned int next_proximal_height_index = target_height_index;
//     if (old_height < new_height) {
//         next_proximal_height_index = comparisons->get_nearest_larger_height_index(target_height_index, true);
//         lower_bound = comparisons->get_height(target_height_index);
//         upper_bound = comparisons->get_height(next_proximal_height_index);
//         if ((next_proximal_height_index == target_height_index) ||
//            ((upper_bound - old_height) > window_size)) {
//             upper_bound = old_height + window_size;
//         }
//     }
//     else {
//         next_proximal_height_index = comparisons->get_nearest_smaller_height_index(target_height_index, true);
//         upper_bound = comparisons->get_height(target_height_index);
//         lower_bound = comparisons->get_height(next_proximal_height_index);
//         if ((next_proximal_height_index == target_height_index) ||
//            ((old_height - lower_bound) > window_size)) {
//             lower_bound = old_height - window_size;
//             // if (lower_bound < 0.0) {
//             //     lower_bound = 0.0;
//             // }
//         }
//     }
//     double merge_gap = upper_bound - lower_bound;
// 
//     // std::cout << "Merging:\n";
//     // std::cout << "old height: " << old_height << " " << comparisons->get_height(height_index) << "\n";
//     // std::cout << "new height: " << new_height << "\n";
//     // std::cout << "targe height: " << comparisons->get_height(target_height_index) << "\n";
//     // std::cout << "next height: " << comparisons->get_height(next_proximal_height_index) << "\n";
//     // std::cout << "height index: " << height_index << "\n";
//     // std::cout << "target_height index: " << target_height_index << "\n";
//     // std::cout << "next_proximal_height_index: " << next_proximal_height_index << "\n";
//     // std::cout << "window size: " << window_size << "\n";
//     // std::cout << "upper bond: " << upper_bound << "\n";
//     // std::cout << "lower bond: " << lower_bound << "\n";
//     // std::cout << "merge gap: " << merge_gap << "\n";
//     // std::cout << "\n";
// 
//     ECOEVOLITY_ASSERT(merge_gap > 0.0);
// 
//     unsigned int new_merged_event_index = comparisons->merge_height(height_index, target_height_index);
//     unsigned int nnodes_in_merged_event = comparisons->get_number_of_trees_mapped_to_height(new_merged_event_index);
//     // Don't need the returned probability vector, but need to make sure we
//     // update the stored number of ways to make the reverse split of this
//     // number of nodes.
//     this->get_split_subset_size_probabilities(nnodes_in_merged_event);
// 
//     // TODO: check this
//     double ln_jacobian = 0.0;
// 
//     // The probability of the forward merge move is simply the probability of
//     // randomly selecting the height to merge from among all heights times the
//     // probability of drawing a height that will lead to the merge
//     // p(selecting height) = 1 / (number of events before proposal)
//     // p(drawing a merge height) = width of window that leads to merge / 2w
//     //                           = merge gap / 2w 
//     // where w = window size
//     //
//     // p(merge) = merge gap / (2w * number of events before proposal)
//     //
//     // The probability of reverse split move is the product of the probabilites of
//     //   1) choosing the merged shared event to split
//     //          = 1 / number of shared events after the proposal
//     //   2) randomly splitting the subset out of the 'n' nodes in the merged
//     //      event (the number of nodes in the event after the proposal).
//     //          = 1 / (2 * stirling2(n, 2))
//     //          = 1 / (2^n - 2)
//     //   3) randomly drawing the new time = 1 / 2w
//     // So the prob of the reverse split move is
//     // p(split move) =  1 / (number of shared events after the proposal * (2^n - 2) * 2w)
//     //                  1 / (number of shared events after the proposal * 2 * stirling2(n, 2) * 2w)
//     //
//     // So, the Hasting ratio for the proposed split move is:
//     // p(reverse split) / p(merge) = 
//     //                  (number of events before the proposal * 2w)
//     //     ------------------------------------------------------------------------------------
//     //     (number of shared events after the proposal * 2 * stirling2(n, 2) * 2w * merge gap)
//     // =
//     //                  (number of events before the proposal)
//     //     ------------------------------------------------------------------------------------
//     //     (number of shared events after the proposal * 2 * stirling2(n, 2) * merge gap)
//     ECOEVOLITY_ASSERT(comparisons->get_number_of_events() + 1 == nevents);
//     unsigned int nshared_after = comparisons->get_shared_event_indices().size();
//     double ln_hastings = std::log(nevents);
//     ln_hastings -= (
//             std::log(nshared_after) +
//             this->ln_number_of_possible_splits_.at(nnodes_in_merged_event) +
//             std::log(merge_gap));
// 
//     // Account for hastings of proposed height
//     ln_hastings += new_height_ln_hastings;
// 
//     const bool in_shared_state_after = (comparisons->get_number_of_events() == 1);
//     // Account for probability of choosing to merge
//     // This is 1.0 (or zero on log scale) except for these two corner
//     // cases:
//     if (in_general_state_before && (! in_shared_state_after)) {
//         ln_hastings -= std::log(2.0);
//     }
//     else if (in_shared_state_after && (! in_general_state_before)) {
//         ln_hastings += std::log(2.0);
//     }
// 
//     // double merge_overhang = window_size - std::abs(old_height - comparisons->get_height(new_merged_event_index));
//     // std::cout << "merging:\n";
//     // std::cout << "hastings: " << std::exp(ln_hastings) << "\n";
//     // std::cout << "1/gap: " << 1.0/merge_gap << "\n";
//     // std::cout << "1/overhang: " << 1.0/merge_overhang << "\n";
//     // std::cout << "\n";
// 
//     return ln_hastings + ln_jacobian;
// }
