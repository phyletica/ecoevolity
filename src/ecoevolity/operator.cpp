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
    // Some prelim tests confirm that an acceptance rate of 0.44 leads to
    // better mixing for the simple, univariate random variables that the
    // operators are updating.
    // return 0.234;
    return 0.44;
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

double Operator::calc_delta(OperatorSchedule& os,
        double log_alpha) const {
    return os.calc_delta(*this, log_alpha);
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
// CollectionScaler methods
//////////////////////////////////////////////////////////////////////////////

Operator::OperatorTypeEnum CollectionScaler::get_type() const {
    return Operator::OperatorTypeEnum::model_operator;
}

double CollectionScaler::propose(RandomNumberGenerator& rng,
        ComparisonPopulationTreeCollection& comparisons,
        unsigned int nthreads) {
    double multiplier = std::exp(this->scale_ * ((2.0 * rng.uniform_real()) - 1.0));
    unsigned int number_of_free_parameters_scaled = 0;
    for (unsigned int tree_idx = 0;
            tree_idx < comparisons.trees_.size();
            ++tree_idx) {
        number_of_free_parameters_scaled += comparisons.trees_.at(tree_idx).scale_population_sizes(multiplier);
    }
    for (unsigned int height_idx = 0;
            height_idx < comparisons.node_heights_.size();
            ++height_idx) {
        comparisons.node_heights_.at(height_idx)->set_value(
                comparisons.node_heights_.at(height_idx)->get_value() * multiplier);
        ++number_of_free_parameters_scaled;
    }
    comparisons.make_trees_dirty();

    return std::log(multiplier) * number_of_free_parameters_scaled;
}

// The following element-wise proposal worked in the sense that it sampled from
// the expected prior when there were no data. However, it did not improve the
// mixing of node heights and pop sizes when there are no constant sits
//
// double CollectionScaler::propose(RandomNumberGenerator& rng,
//         ComparisonPopulationTreeCollection& comparisons,
//         unsigned int nthreads) {
//     double ln_hastings = 0.0;
//     double size;
//     double height;
//     double hastings;
//     for (unsigned int tree_idx = 0;
//             tree_idx < comparisons.trees_.size();
//             ++tree_idx) {
//         size = comparisons.trees_.at(tree_idx).get_root_population_size();
//         this->update(rng, size, hastings);
//         ln_hastings += hastings;
//         comparisons.trees_.at(tree_idx).set_root_population_size(size);
//         if (! comparisons.trees_.at(tree_idx).population_sizes_are_constrained()) {
//             for (unsigned int i = 0; i < comparisons.trees_.at(tree_idx).get_leaf_node_count(); ++i) {
//                 size = comparisons.trees_.at(tree_idx).get_child_population_size(i);
//                 this->update(rng, size, hastings);
//                 ln_hastings += hastings;
//                 comparisons.trees_.at(tree_idx).set_child_population_size(i, size);
//             }
//         }
//     }
//     for (unsigned int height_idx = 0;
//             height_idx < comparisons.node_heights_.size();
//             ++height_idx) {
//         height = comparisons.node_heights_.at(height_idx)->get_value();
//         this->update(rng, height, hastings);
//         ln_hastings += hastings;
//         comparisons.node_heights_.at(height_idx)->set_value(height);
//     }
//     comparisons.make_trees_dirty();
//
//     return ln_hastings;
// }

std::string CollectionScaler::target_parameter() const {
    return "node heights and population sizes";
}

std::string CollectionScaler::get_name() const {
    return "CollectionScaler";
}


//////////////////////////////////////////////////////////////////////////////
// ConcentrationScaler methods
//////////////////////////////////////////////////////////////////////////////

Operator::OperatorTypeEnum ConcentrationScaler::get_type() const {
    return Operator::OperatorTypeEnum::model_operator;
}

double ConcentrationScaler::propose(RandomNumberGenerator& rng,
        ComparisonPopulationTreeCollection& comparisons,
        unsigned int nthreads) {
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
// FreqMover methods
//////////////////////////////////////////////////////////////////////////////

double FreqMover::propose(
        RandomNumberGenerator& rng,
        ComparisonPopulationTree& tree) const {
    double freq_1 = tree.get_freq_1();
    double hastings;
    this->update(rng, freq_1, hastings);
    if ((freq_1 <= 0.0) || (freq_1 > 1.0)) {
        return -std::numeric_limits<double>::infinity();
    }
    tree.set_freq_1(freq_1);
    return hastings; 
}

std::string FreqMover::target_parameter() const {
    return "freq 1";
}

std::string FreqMover::get_name() const {
    return "FreqMover";
}



//////////////////////////////////////////////////////////////////////////////
// UScaler methods
//////////////////////////////////////////////////////////////////////////////

// u/v rates no longer being proposed
// double UScaler::propose(
//         RandomNumberGenerator& rng,
//         ComparisonPopulationTree& tree) const {
//     double val = tree.get_u();
//     double hastings;
//     this->update(rng, val, hastings);
//     if (val < 0.5) {
//         return -std::numeric_limits<double>::infinity();
//     }
//     // v is also set here
//     tree.set_u(val);
//     return hastings;
// }
// 
// std::string UScaler::target_parameter() const {
//     return "u rate";
// }
// 
// std::string UScaler::get_name() const {
//     return "UScaler";
// }


//////////////////////////////////////////////////////////////////////////////
// ComparisonMutationRateScaler methods
//////////////////////////////////////////////////////////////////////////////

double ComparisonMutationRateScaler::propose(
        RandomNumberGenerator& rng,
        ComparisonPopulationTree& tree) const {
    double v = tree.get_mutation_rate();
    double hastings;
    this->update(rng, v, hastings);
    tree.set_mutation_rate(v);
    return hastings;
}

std::string ComparisonMutationRateScaler::target_parameter() const {
    return "mutation rate";
}

std::string ComparisonMutationRateScaler::get_name() const {
    return "ComparisonMutationRateScaler";
}


//////////////////////////////////////////////////////////////////////////////
// ChildPopulationSizeScaler methods
//////////////////////////////////////////////////////////////////////////////

double ChildPopulationSizeScaler::propose(
        RandomNumberGenerator& rng,
        ComparisonPopulationTree& tree) const {
    int pop_idx = rng.uniform_int(0, tree.get_leaf_node_count() - 1);
    double size = tree.get_child_population_size(pop_idx);

    double hastings;
    this->update(rng, size, hastings);

    // avoid zero division to get coalescence rate
    if (size <= 0.0) {
        return -std::numeric_limits<double>::infinity();
    }

    tree.set_child_population_size(pop_idx, size);
    return hastings;
}

std::string ChildPopulationSizeScaler::target_parameter() const {
    return "population size";
}

std::string ChildPopulationSizeScaler::get_name() const {
    return "ChildPopulationSizeScaler";
}


//////////////////////////////////////////////////////////////////////////////
// RootPopulationSizeScaler methods
//////////////////////////////////////////////////////////////////////////////

double RootPopulationSizeScaler::propose(
        RandomNumberGenerator& rng,
        ComparisonPopulationTree& tree) const {
    double size = tree.get_root_population_size();

    double hastings;
    this->update(rng, size, hastings);

    // avoid zero division to get coalescence rate
    if (size <= 0.0) {
        return -std::numeric_limits<double>::infinity();
    }

    tree.set_root_population_size(size);

    return hastings;
}

std::string RootPopulationSizeScaler::target_parameter() const {
    return "population size";
}

std::string RootPopulationSizeScaler::get_name() const {
    return "RootPopulationSizeScaler";
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
// ComparisonHeightMover methods
//////////////////////////////////////////////////////////////////////////////

double ComparisonHeightMover::propose(
        RandomNumberGenerator& rng,
        PositiveRealParameter& node_height) const {
    double h = node_height.get_value();
    double hastings;
    this->update(rng, h, hastings);
    if (h < 0.0) {
        return -std::numeric_limits<double>::infinity();
    }
    node_height.set_value(h);
    return hastings;
}

std::string ComparisonHeightMover::get_name() const {
    return "ComparisonHeightMover";
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
        ComparisonPopulationTreeCollection& comparisons,
        unsigned int nthreads) {

    const double ln_concentration_over_num_aux = std::log(
            comparisons.get_concentration() /
            this->get_number_of_auxiliary_categories());

    for (unsigned int tree_idx = 0;
            tree_idx < comparisons.get_number_of_trees();
            ++tree_idx) {
        std::vector<unsigned int> other_height_indices = comparisons.get_other_height_indices(tree_idx);
        std::vector<double> ln_category_likelihoods;
        std::vector<double> ln_tree_likelihoods;
        ln_category_likelihoods.reserve(other_height_indices.size() +
                this->get_number_of_auxiliary_categories());
        ln_tree_likelihoods.reserve(other_height_indices.size() +
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
                    ln_tree_likelihoods.push_back(lnl);
                    continue;
                }
            }
            else {
                comparisons.trees_.at(tree_idx).set_height(comparisons.node_heights_.at(height_idx)->get_value());
            }
            double lnl = comparisons.trees_.at(tree_idx).compute_log_likelihood(nthreads);
            ln_category_likelihoods.push_back(lnl + std::log(number_of_elements));
            ln_tree_likelihoods.push_back(lnl);
        }

        std::vector<double> auxiliary_heights;
        auxiliary_heights.reserve(this->get_number_of_auxiliary_categories());
        for (unsigned int i = 0; i < this->get_number_of_auxiliary_categories(); ++i) {
            double fresh_height = comparisons.node_height_prior_->draw(rng);
            auxiliary_heights.push_back(fresh_height);
            comparisons.trees_.at(tree_idx).set_height(fresh_height);
            double lnl = comparisons.trees_.at(tree_idx).compute_log_likelihood(nthreads);
            ln_category_likelihoods.push_back(lnl + ln_concentration_over_num_aux);
            ln_tree_likelihoods.push_back(lnl);
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
                    ln_tree_likelihoods.at(prob_index));
        }
        else {
            comparisons.map_tree_to_new_height(
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

Operator::OperatorTypeEnum ReversibleJumpSampler::get_type() const {
    return Operator::OperatorTypeEnum::rj_operator;
}

double ReversibleJumpSampler::propose(RandomNumberGenerator& rng,
        ComparisonPopulationTreeCollection& comparisons,
        unsigned int nthreads) {
    return this->propose_jump_to_gap(rng, comparisons);
}

double ReversibleJumpSampler::propose_jump_to_prior(RandomNumberGenerator& rng,
        ComparisonPopulationTreeCollection& comparisons) {
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
    const unsigned int nnodes = comparisons.get_number_of_trees();
    const unsigned int nevents = comparisons.get_number_of_events();
    const bool in_general_state_before = (nnodes == nevents);
    const bool in_shared_state_before = (nevents == 1);
    const bool split_event = ((! in_general_state_before) &&
            (in_shared_state_before || (rng.uniform_real() < 0.5)));
    double mean_height = comparisons.node_height_prior_->get_mean();
    if (split_event) {
        std::vector<unsigned int> shared_indices =
                comparisons.get_shared_event_indices();
        unsigned int i = rng.uniform_int(0, shared_indices.size() - 1);
        unsigned int event_index = shared_indices.at(i);
        double new_height = comparisons.node_height_prior_->draw(rng);

        std::vector<unsigned int> tree_indices = comparisons.get_indices_of_mapped_trees(event_index);
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
        comparisons.map_trees_to_new_height(subset_indices, new_height);

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

        const bool in_general_state_after = (comparisons.get_number_of_trees() ==
                comparisons.get_number_of_events());
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
    double removed_height = comparisons.get_height(move_height_index);
    unsigned int new_merged_event_index = comparisons.merge_height(move_height_index, target_height_index);
    unsigned int nnodes_in_merged_event = comparisons.get_number_of_trees_mapped_to_height(new_merged_event_index);
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
    ECOEVOLITY_ASSERT(comparisons.get_number_of_events() + 1 == nevents);
    unsigned int nshared_after = comparisons.get_shared_event_indices().size();
    double ln_hastings = std::log(nevents) + std::log(nevents - 1);
    ln_hastings -= (
            std::log(nshared_after) +
            this->ln_number_of_possible_splits_.at(nnodes_in_merged_event));

    const bool in_shared_state_after = (comparisons.get_number_of_events() == 1);
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
        ComparisonPopulationTreeCollection& comparisons) {
    const unsigned int nnodes = comparisons.get_number_of_trees();
    const unsigned int nevents = comparisons.get_number_of_events();
    const bool in_general_state_before = (nnodes == nevents);
    const bool in_shared_state_before = (nevents == 1);
    const bool split_event = ((! in_general_state_before) &&
            (in_shared_state_before || (rng.uniform_real() < 0.5)));
    if (split_event) {
        std::vector<unsigned int> shared_indices =
                comparisons.get_shared_event_indices();
        unsigned int i = rng.uniform_int(0, shared_indices.size() - 1);
        unsigned int event_index = shared_indices.at(i);
        double event_height = comparisons.get_height(event_index);
        double lower_bound = comparisons.get_nearest_smaller_height(event_index);
        double new_height = rng.uniform_real(lower_bound, event_height);

        std::vector<unsigned int> tree_indices = comparisons.get_indices_of_mapped_trees(event_index);
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
        comparisons.map_trees_to_new_height(subset_indices, new_height);

        // TODO: check this
        double ln_jacobian = 0.0;

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

        const bool in_general_state_after = (comparisons.get_number_of_trees() ==
                comparisons.get_number_of_events());
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
    std::vector<unsigned int> candidate_indices =
            comparisons.get_height_indices_sans_largest();
    unsigned int i = rng.uniform_int(0, candidate_indices.size() - 1);
    unsigned int height_index = candidate_indices.at(i);
    unsigned int target_height_index = comparisons.get_nearest_larger_height_index(height_index);
    unsigned int new_merged_event_index = comparisons.merge_height(height_index, target_height_index);
    unsigned int nnodes_in_merged_event = comparisons.get_number_of_trees_mapped_to_height(new_merged_event_index);
    // Don't need the returned probability vector, but need to make sure we
    // update the stored number of ways to make the reverse split of this
    // number of nodes.
    this->get_split_subset_size_probabilities(nnodes_in_merged_event);

    // TODO: check this
    double ln_jacobian = 0.0;

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
    ECOEVOLITY_ASSERT(comparisons.get_number_of_events() + 1 == nevents);
    unsigned int nshared_after = comparisons.get_shared_event_indices().size();
    double new_merged_event_height = comparisons.get_height(new_merged_event_index);
    double lower_bound = comparisons.get_nearest_smaller_height(new_merged_event_index);
    double ln_hastings = std::log(comparisons.get_number_of_events());
    ln_hastings -= (
            std::log(nshared_after) +
            this->ln_number_of_possible_splits_.at(nnodes_in_merged_event) +
            std::log(new_merged_event_height - lower_bound));

    const bool in_shared_state_after = (comparisons.get_number_of_events() == 1);
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
