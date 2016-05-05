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
// Operator methods
//////////////////////////////////////////////////////////////////////////////

Operator::Operator(double weight) {
    this->set_weight(weight);
}

double Operator::get_target_acceptance_probability() const {
    return 0.234;
}

double Operator::get_coercable_parameter_value() const {
    return std::numeric_limits<double>::quiet_NaN();
}

void Operator::set_weight(double weight) {
    ECOEVOLITY_ASSERT(weight > 0.0);
    this->weight_ = weight;
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

std::string Operator::get_name() const {
    return "BaseOperator";
}

std::string Operator::header_string() const {
    return "name\tnumber_accepted\tnumber_rejected\tweight\tweight_prob\ttuning_parameter\n";
}

std::string Operator::to_string(const OperatorSchedule& os) const {
    std::ostringstream ss;
    ss << this->get_name() << "\t" 
       << this->get_number_accepted() << "\t"
       << this->get_number_rejected() << "\t"
       << this->get_weight() << "\t";

    if (os.get_total_weight() > 0.0) {
        ss << this->get_weight() / os.get_total_weight() << "\t";
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

double Operator::calc_delta(OperatorSchedule& os,
        double log_alpha) const {
    return os.calc_delta(shared_from_this(), log_alpha);
}


//////////////////////////////////////////////////////////////////////////////
// ScaleOperator methods
//////////////////////////////////////////////////////////////////////////////

ScaleOperator::ScaleOperator(double weight, double scale) : Operator(weight) {
    this->set_scale(scale);
}

void ScaleOperator::set_scale(double scale) {
    ECOEVOLITY_ASSERT(scale > 0.0);
    this->scale_ = scale;
}

void ScaleOperator::update(
        RandomNumberGenerator& rng,
        double& parameter_value,
        double& hastings_ratio) const {
    double multiplier = std::exp(this->scale_ * ((2.0 * rng.uniform_real()) - 1.0));
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

std::string ScaleOperator::get_name() const {
    return "BaseScaleOperator";
}


//////////////////////////////////////////////////////////////////////////////
// WindowOperator methods
//////////////////////////////////////////////////////////////////////////////

WindowOperator::WindowOperator(double weight, double window_size) : Operator(weight) {
    this->set_window_size(window_size);
}
void WindowOperator::set_window_size(double window_size) {
    ECOEVOLITY_ASSERT(window_size > 0.0);
    this->window_size_ = window_size;
}
double WindowOperator::get_window_size() const {
    return this->window_size_;
}

void WindowOperator::update(
        RandomNumberGenerator& rng,
        double& parameter_value,
        double& hastings_ratio) const {
    double addend = (rng.uniform_real() * 2 * this->window_size_) - this->window_size_;
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

std::string WindowOperator::get_name() const {
    return "BaseWindowOperator";
}


//////////////////////////////////////////////////////////////////////////////
// ModelOperator methods
//////////////////////////////////////////////////////////////////////////////

Operator::OperatorTypeEnum ModelOperator::get_type() const {
    return Operator::OperatorTypeEnum::model_operator;
}

std::string ModelOperator::get_name() const {
    return "ModelOperator";
}

std::string ModelOperator::target_parameter() const {
    return "model";
}


//////////////////////////////////////////////////////////////////////////////
// ComparisonTreeScaleOperator methods
//////////////////////////////////////////////////////////////////////////////

Operator::OperatorTypeEnum ComparisonTreeScaleOperator::get_type() const {
    return Operator::OperatorTypeEnum::tree_operator;
}

std::string ComparisonTreeScaleOperator::get_name() const {
    return "ComparisonTreeScaleOperator";
}

//////////////////////////////////////////////////////////////////////////////
// ComparisonTreeWindowOperator methods
//////////////////////////////////////////////////////////////////////////////

Operator::OperatorTypeEnum ComparisonTreeWindowOperator::get_type() const {
    return Operator::OperatorTypeEnum::tree_operator;
}

std::string ComparisonTreeWindowOperator::get_name() const {
    return "ComparisonTreeWindowOperator";
}


//////////////////////////////////////////////////////////////////////////////
// NodeHeightScaleOperator methods
//////////////////////////////////////////////////////////////////////////////

Operator::OperatorTypeEnum NodeHeightScaleOperator::get_type() const {
    return Operator::OperatorTypeEnum::time_operator;
}

std::string NodeHeightScaleOperator::get_name() const {
    return "NodeHeightScaleOperator";
}

std::string NodeHeightScaleOperator::target_parameter() const {
    return "node height";
}


//////////////////////////////////////////////////////////////////////////////
// NodeHeightWindowOperator methods
//////////////////////////////////////////////////////////////////////////////

Operator::OperatorTypeEnum NodeHeightWindowOperator::get_type() const {
    return Operator::OperatorTypeEnum::time_operator;
}

std::string NodeHeightWindowOperator::get_name() const {
    return "NodeHeightWindowOperator";
}

std::string NodeHeightWindowOperator::target_parameter() const {
    return "node height";
}


//////////////////////////////////////////////////////////////////////////////
// ConcentrationScaler methods
//////////////////////////////////////////////////////////////////////////////

Operator::OperatorTypeEnum ConcentrationScaler::get_type() const {
    return Operator::OperatorTypeEnum::model_operator;
}

double ConcentrationScaler::propose(RandomNumberGenerator& rng,
        ComparisonPopulationTreeCollection& comparisons) const {
    double v = comparisons.get_concentration();
    double hastings;
    this->update(rng, v, hastings);
    comparisons.set_concentration(v);
    return hastings;
}

std::string ConcentrationScaler::target_parameter() const {
    return "concentration";
}

std::string ConcentrationScaler::get_name() const {
    return "ConcentrationScaler";
}


//////////////////////////////////////////////////////////////////////////////
// MutationRateMover methods
//////////////////////////////////////////////////////////////////////////////

double MutationRateMover::propose(
        RandomNumberGenerator& rng,
        ComparisonPopulationTree& tree) const {
    double red_freq = tree.get_v();
    double hastings;
    this->update(rng, red_freq, hastings);
    if ((red_freq >= 0.0) && (red_freq <= 1.0)) {
        // v is also set here
        tree.set_u(1.0 / (2.0 * red_freq));
        return hastings; 
    }
    return -std::numeric_limits<double>::infinity();
}

std::string MutationRateMover::target_parameter() const {
    return "mutation rate";
}

std::string MutationRateMover::get_name() const {
    return "MutationRateMover";
}


//////////////////////////////////////////////////////////////////////////////
// ComparisonHeightMultiplierScaler methods
//////////////////////////////////////////////////////////////////////////////

double ComparisonHeightMultiplierScaler::propose(
        RandomNumberGenerator& rng,
        ComparisonPopulationTree& tree) const {
    double v = tree.get_node_height_multiplier();
    double hastings;
    this->update(rng, v, hastings);
    tree.set_node_height_multiplier(v);
    return hastings;
}

std::string ComparisonHeightMultiplierScaler::target_parameter() const {
    return "node height multiplier";
}

std::string ComparisonHeightMultiplierScaler::get_name() const {
    return "ComparisonHeightMultiplierScaler";
}


//////////////////////////////////////////////////////////////////////////////
// ChildCoalescenceRateScaler methods
//////////////////////////////////////////////////////////////////////////////

double ChildCoalescenceRateScaler::propose(
        RandomNumberGenerator& rng,
        ComparisonPopulationTree& tree) const {
    int pop_idx = rng.uniform_int(0, tree.get_leaf_node_count() - 1);
    double rate = tree.get_child_coalescence_rate(pop_idx);

    double hastings;
    this->update(rng, rate, hastings);

    tree.set_child_coalescence_rate(pop_idx, rate);
    return hastings;
}

std::string ChildCoalescenceRateScaler::target_parameter() const {
    return "coalescence rate";
}

std::string ChildCoalescenceRateScaler::get_name() const {
    return "ChildCoalescenceRateScaler";
}


//////////////////////////////////////////////////////////////////////////////
// RootCoalescenceRateScaler methods
//////////////////////////////////////////////////////////////////////////////

double RootCoalescenceRateScaler::propose(
        RandomNumberGenerator& rng,
        ComparisonPopulationTree& tree) const {
    double rate = tree.get_root_coalescence_rate();

    double hastings;
    this->update(rng, rate, hastings);

    tree.set_root_coalescence_rate(rate);

    return hastings;
}

std::string RootCoalescenceRateScaler::target_parameter() const {
    return "coalescence rate";
}

std::string RootCoalescenceRateScaler::get_name() const {
    return "RootCoalescenceRateScaler";
}


//////////////////////////////////////////////////////////////////////////////
// ComparisonHeightScaler methods
//////////////////////////////////////////////////////////////////////////////

double ComparisonHeightScaler::propose(
        RandomNumberGenerator& rng,
        PositiveRealParameter& node_height) const {
    double h = node_height.get_value();
    double hastings;
    this->update(rng, h, hastings);
    node_height.set_value(h);
    return hastings;
}

std::string ComparisonHeightScaler::get_name() const {
    return "ComparisonHeightScaler";
}


//////////////////////////////////////////////////////////////////////////////
// DirichletProcessGibbsSampler methods
//////////////////////////////////////////////////////////////////////////////

DirichletProcessGibbsSampler::DirichletProcessGibbsSampler(
        double weight,
        unsigned int number_of_auxiliary_categories) : ModelOperator(weight) {
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

double DirichletProcessGibbsSampler::propose(RandomNumberGenerator& rng,
        ComparisonPopulationTreeCollection& comparisons) const {

    const double ln_concentration_over_num_aux = std::log(
            comparisons.get_concentration() /
            this->get_number_of_auxiliary_categories());

    for (unsigned int tree_idx = 0;
            tree_idx < comparisons.get_number_of_trees();
            ++tree_idx) {
        std::vector<unsigned int> other_height_indices = comparisons.get_other_height_indices(tree_idx);
        std::vector<double> ln_category_likelihoods;
        ln_category_likelihoods.reserve(other_height_indices.size() +
                this->get_number_of_auxiliary_categories());

        // store height associated with this tree
        comparisons.node_heights_.at(comparisons.get_height_index(tree_idx))->store();
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
                comparisons.trees_.at(tree_idx).set_height(comparisons.node_heights_.at(height_idx)->get_value());
            }
            double lnl = comparisons.trees_.at(tree_idx).compute_log_likelihood();
            ln_category_likelihoods.push_back(lnl + std::log(number_of_elements));
        }

        std::vector<double> auxiliary_heights;
        auxiliary_heights.reserve(this->get_number_of_auxiliary_categories());
        for (unsigned int i = 0; i < this->get_number_of_auxiliary_categories(); ++i) {
            double fresh_height = comparisons.node_height_prior_->draw(rng);
            auxiliary_heights.push_back(fresh_height);
            comparisons.trees_.at(tree_idx).set_height(fresh_height);
            double lnl = comparisons.trees_.at(tree_idx).compute_log_likelihood();
            ln_category_likelihoods.push_back(lnl + ln_concentration_over_num_aux);
        }

        // restore height associated with this tree
        comparisons.node_heights_.at(comparisons.get_height_index(tree_idx))->restore();

        std::vector<double> category_probs(ln_category_likelihoods);
        normalize_log_likelihoods(category_probs);
        unsigned int prob_index = rng.weighted_index(category_probs);
        if (prob_index < other_height_indices.size()) {
            comparisons.remap_tree(
                    tree_idx,
                    other_height_indices.at(prob_index),
                    ln_category_likelihoods.at(prob_index));
        }
        else {
            comparisons.map_tree_to_new_height(
                    tree_idx,
                    auxiliary_heights.at(prob_index -
                            other_height_indices.size()),
                    ln_category_likelihoods.at(prob_index));
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

double ReversibleJumpSampler::propose(RandomNumberGenerator& rng,
        ComparisonPopulationTreeCollection& comparisons) const {
    throw EcoevolityNotImplementedError("rjMCMC not implemented yet");
}
