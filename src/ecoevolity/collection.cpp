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

#include "collection.hpp"
#include "operator.hpp"

ComparisonPopulationTreeCollection::ComparisonPopulationTreeCollection(
        const CollectionSettings & settings,
        RandomNumberGenerator & rng
        ) {
    this->state_log_path_ = settings.get_state_log_path();
    this->operator_log_path_ = settings.get_operator_log_path();
    this->node_height_prior_ = settings.get_time_prior_settings().get_instance();
    this->concentration_ = std::make_shared<PositiveRealParameter>(
            settings.get_concentration_settings(),
            rng);
    this->operator_schedule_ = OperatorSchedule(
            settings.get_operator_schedule_settings(),
            settings.using_dpp());
    this->init_trees(settings.get_comparison_settings(), rng);
}

void ComparisonPopulationTreeCollection::init_trees(
    const std::vector<ComparisonSettings> & comparison_settings,
    RandomNumberGenerator & rng) {
    double fresh_height;
    for (unsigned int tree_idx = 0;
            tree_idx < comparison_settings.size();
            ++tree_idx) {
        fresh_height = this->node_height_prior_->draw(rng);
        std::shared_ptr<PositiveRealParameter> new_height_parameter = std::make_shared<PositiveRealParameter>(this->node_height_prior_, fresh_height);
        ComparisonPopulationTree new_tree = ComparisonPopulationTree(
                comparison_settings.at(tree_idx),
                rng);
        new_tree.set_node_height_prior(this->node_height_prior_);
        new_tree.set_height_parameter(new_height_parameter);
        this->node_heights_.push_back(new_height_parameter);
        this->trees_.push_back(new_tree);
        this->node_height_indices_.push_back(tree_idx);
        ECOEVOLITY_ASSERT(
                this->trees_.at(tree_idx).get_height_parameter() ==
                this->node_heights_.at(this->node_height_indices_.at(tree_idx))
                );
    }
    ECOEVOLITY_ASSERT(this->trees_.size() == this->node_heights_.size());
    ECOEVOLITY_ASSERT(this->trees_.size() == this->node_height_indices_.size());
}

void ComparisonPopulationTreeCollection::store_state() {
    this->log_likelihood_.store();
    this->log_prior_density_.store();
    for (unsigned int i = 0; i < this->trees_.size(); ++i) {
        this->trees_.at(i).store_state();
    }
}
void ComparisonPopulationTreeCollection::restore_state() {
    this->log_likelihood_.restore();
    this->log_prior_density_.restore();
    for (unsigned int i = 0; i < this->trees_.size(); ++i) {
        this->trees_.at(i).restore_state();
    }
}
void ComparisonPopulationTreeCollection::compute_log_likelihood_and_prior(bool compute_partials) {
    double lnl = 0.0;
    double lnp = 0.0;
    if (compute_partials) {
        if ((this->use_multithreading_) && (this->trees_.size() > 1)) {
            this->compute_tree_partials_threaded();
        }
        else {
            this->compute_tree_partials();
        }
    }
    for (unsigned int i = 0; i < this->trees_.size(); ++i) {
        lnl += this->trees_.at(i).get_log_likelihood_value();
        lnp += this->trees_.at(i).get_log_prior_density_value();
    }
    for (unsigned int h = 0; h < this->node_heights_.size(); ++h) {
        lnp += this->node_heights_.at(h)->relative_prior_ln_pdf();
    }

    this->log_likelihood_.set_value(lnl);
    this->log_prior_density_.set_value(lnp);
}

void ComparisonPopulationTreeCollection::compute_tree_partials() {
    for (unsigned int i = 0; i < this->trees_.size(); ++i) {
        this->trees_.at(i).compute_log_likelihood_and_prior();
    }
}

void ComparisonPopulationTreeCollection::compute_tree_partials_threaded() {
    this->compute_tree_partials();
    // TODO: get multithreading working
    // std::vector<std::thread> threads;
    // threads.reserve(this->trees_.size());

    // for (unsigned int i = 0; i < this->trees_.size(); ++i) {
    //     threads.push_back(std::thread(
    //                 &ComparisonPopulationTree::compute_log_likelihood_and_prior,
    //                 this->trees_.at(i)));
    // }

    // for (unsigned int i = 0; i < threads.size(); ++i) {
    //     threads.at(i).join()
    // }
}

void ComparisonPopulationTreeCollection::make_trees_clean() {
    for (unsigned int i = 0; i < this->trees_.size(); ++i) {
        this->trees_.at(i).make_clean();
    }
}
void ComparisonPopulationTreeCollection::make_trees_dirty() {
    for (unsigned int i = 0; i < this->trees_.size(); ++i) {
        this->trees_.at(i).make_dirty();
    }
}

std::vector<unsigned int> ComparisonPopulationTreeCollection::get_standardized_height_indices() const {
    std::vector<unsigned int> standardized_indices;
    standardized_indices.reserve(this->node_height_indices_.size());
    std::map<unsigned int, unsigned int> standardizing_map;
    unsigned int next_idx = 0;
    for (auto const raw_idx: this->node_height_indices_) {
        if (standardizing_map.count(raw_idx) == 0) {
            standardizing_map[raw_idx] = next_idx;
            ++next_idx;
        }
        standardized_indices.push_back(standardizing_map.at(raw_idx));
    }
    return standardized_indices;
}

unsigned int ComparisonPopulationTreeCollection::get_number_of_trees_mapped_to_height(
        unsigned int height_index) const {
    unsigned int count = 0;
    for (auto h_index_iter : this->node_height_indices_) {
        if (h_index_iter == height_index) {
            ++count;
        }
    }
    return count;
}

std::vector<unsigned int> ComparisonPopulationTreeCollection::get_other_height_indices(
        unsigned int tree_index) const {
    std::vector<unsigned int> others;
    others.reserve(this->node_heights_.size());
    if (this->get_number_of_partners(tree_index) > 0) {
        for (unsigned int i = 0; i < this->node_heights_.size(); ++i) {
            others.push_back(i);
        }
    }
    else {
        unsigned int removed_height_index = this->get_height_index(tree_index);
        for (unsigned int i = 0; i < this->node_heights_.size(); ++i) {
            if (i != removed_height_index) {
                others.push_back(i);
            }
        }
    }
    return others;
}

unsigned int ComparisonPopulationTreeCollection::get_number_of_partners(
        unsigned int tree_index) const {
    unsigned int h_index = this->node_height_indices_.at(tree_index);
    unsigned int count = this->get_number_of_trees_mapped_to_height(h_index);
    // should always find at least itself
    ECOEVOLITY_ASSERT(count > 0);
    return count - 1;
}

void ComparisonPopulationTreeCollection::remap_tree(
        unsigned int tree_index,
        unsigned int height_index,
        double log_likelihood) {
    if (! (this->get_height_index(tree_index) == height_index)) {
        unsigned int num_partners = this->get_number_of_partners(tree_index);
        unsigned int old_height_index = this->get_height_index(tree_index);
        this->node_height_indices_.at(tree_index) = height_index;
        this->trees_.at(tree_index).set_height_parameter(
                this->get_height_parameter(height_index));
        if (num_partners == 0) {
            this->remove_height(old_height_index);
        }
    }
    this->trees_.at(tree_index).log_likelihood_.set_value(log_likelihood);
    this->trees_.at(tree_index).make_clean();
}

void ComparisonPopulationTreeCollection::map_tree_to_new_height(
        unsigned int tree_index,
        double height,
        double log_likelihood) {
    unsigned int num_partners = this->get_number_of_partners(tree_index);
    unsigned int old_height_index = this->get_height_index(tree_index);
    std::vector<unsigned int> mapped_tree_indices = {tree_index};
    this->add_height(height, mapped_tree_indices);
    if (num_partners == 0) {
        this->remove_height(old_height_index);
    }
    this->trees_.at(tree_index).log_likelihood_.set_value(log_likelihood);
    this->trees_.at(tree_index).make_clean();
}

void ComparisonPopulationTreeCollection::remove_height(
        unsigned int height_index) {
    ECOEVOLITY_ASSERT(this->get_number_of_trees_mapped_to_height(height_index) == 0);
    ECOEVOLITY_ASSERT(this->node_heights_.at(height_index).use_count() < 2);
    this->node_heights_.erase(this->node_heights_.begin() + height_index);
    for (unsigned int i = 0; i < this->node_height_indices_.size(); ++i) {
        ECOEVOLITY_ASSERT(this->node_height_indices_.at(i) != height_index);
        if (this->node_height_indices_.at(i) > height_index) {
            --this->node_height_indices_.at(i);
        }
    }
}

void ComparisonPopulationTreeCollection::add_height(
        double height,
        const std::vector<unsigned int>& mapped_tree_indices) {
    this->node_heights_.push_back(std::make_shared<PositiveRealParameter>(this->node_height_prior_, height));
    for (auto tree_idx : mapped_tree_indices) {
        this->node_height_indices_.at(tree_idx) = this->node_heights_.size() - 1;
        this->trees_.at(tree_idx).set_height_parameter(this->node_heights_.back());
    }
}

void ComparisonPopulationTreeCollection::write_state_log_header(
        std::ostream& out,
        bool short_summary) const {
    out << "generation" << this->logging_delimiter_
        << "ln_likelihood" << this->logging_delimiter_
        << "number_of_events";
    if (this->using_dpp()) {
        out << this->logging_delimiter_ << "concentration";
    }
    if (short_summary) {
        out << std::endl;
        return;
    }
    for (unsigned int tree_idx = 0;
            tree_idx < this->trees_.size();
            ++tree_idx) {
        out << this->logging_delimiter_;
        this->trees_.at(tree_idx).write_state_log_header(out,
                true,
                this->logging_delimiter_);
    }
    out << std::endl;
}

void ComparisonPopulationTreeCollection::log_state(std::ostream& out,
        unsigned int generation_index,
        bool short_summary) const {
    out << generation_index << this->logging_delimiter_
        << this->log_likelihood_.get_value() << this->logging_delimiter_
        << this->get_number_of_events();
    if (this->using_dpp()) {
        out << this->logging_delimiter_ << this->get_concentration();
    }
    if (short_summary) {
        out << std::endl;
        return;
    }
    std::vector<unsigned int> standardized_height_indices =
            this->get_standardized_height_indices();
    for (unsigned int tree_idx = 0;
            tree_idx < this->trees_.size();
            ++tree_idx) {
        out << this->logging_delimiter_;
        this->trees_.at(tree_idx).log_state(out,
                standardized_height_indices.at(tree_idx),
                this->logging_delimiter_);
    }
    out << std::endl;
}

void ComparisonPopulationTreeCollection::update_log_paths(
        unsigned int max_number_of_attempts) {
    if (! path::exists(this->get_state_log_path())) {
        return;
    }
    unsigned int ntries = 0;
    while (true) {
        this->increment_log_paths();
        ++ntries;
        if (! path::exists(this->get_state_log_path())) {
            break;
        }
        if (ntries > max_number_of_attempts) {
            throw EcoevolityError("Could not generate unique output files");
        }
    }
}

void ComparisonPopulationTreeCollection::increment_log_paths() {
    std::vector<std::string> path_elements;
    std::pair<std::string, std::string> prefix_ext;
    std::string new_suffix;
    int run_number;

    prefix_ext = path::splitext(this->get_state_log_path());
    path_elements = string_util::split(prefix_ext.first, '-');
    run_number = std::stoi(path_elements.back());
    path_elements.pop_back();
    ++run_number;
    
    new_suffix = "-" + std::to_string(run_number) + prefix_ext.second;

    this->set_state_log_path(
            string_util::join(path_elements, "-") + new_suffix);

    prefix_ext = path::splitext(this->get_operator_log_path());
    path_elements = string_util::split(prefix_ext.first, '-');
    path_elements.pop_back();
    this->set_operator_log_path(
            string_util::join(path_elements, "-") + new_suffix);
}

void ComparisonPopulationTreeCollection::mcmc(
        RandomNumberGenerator& rng,
        unsigned int chain_length,
        unsigned int sample_frequency) {

    std::ofstream state_log_stream;
    std::ofstream operator_log_stream;
    this->update_log_paths();
    if (path::exists(this->get_state_log_path())) {
        std::ostringstream message;
        message << "ERROR: The parameter log file \'"
                << this->get_state_log_path()
                << "\' already exists!\n";
        throw EcoevolityError(message.str());
    }
    if (path::exists(this->get_operator_log_path())) {
        std::ostringstream message;
        message << "ERROR: The operator log file \'"
                << this->get_operator_log_path()
                << "\' already exists!\n";
        throw EcoevolityError(message.str());
    }
    state_log_stream.open(this->get_state_log_path());
    operator_log_stream.open(this->get_operator_log_path());
    
    if (! state_log_stream.is_open()) {
        std::ostringstream message;
        message << "ERROR: Could not open parameter log file \'"
                << this->get_state_log_path()
                << "\'\n";
        throw EcoevolityError(message.str());
    }
    if (! operator_log_stream.is_open()) {
        std::ostringstream message;
        message << "ERROR: Could not open operator log file \'"
                << this->get_operator_log_path()
                << "\'\n";
        throw EcoevolityError(message.str());
    }

    std::cout << "State log path: " << this->get_state_log_path() << std::endl;
    std::cout << "Operator log path: " << this->get_operator_log_path() << std::endl;
    
    state_log_stream.precision(this->get_logging_precision());
    operator_log_stream.precision(this->get_logging_precision());

    this->write_state_log_header(state_log_stream);
    this->write_state_log_header(std::cout, true);

    this->make_trees_dirty();
    this->compute_log_likelihood_and_prior(true);
    this->log_state(state_log_stream, 0);
    this->log_state(std::cout, 0, true);

    for (unsigned int gen = 0; gen < chain_length; ++gen) {
        this->store_state();
        Operator& op = this->operator_schedule_.draw_operator(rng);
        if (op.get_type() == Operator::OperatorTypeEnum::tree_operator) {
            std::vector<double> hastings_ratios;
            hastings_ratios.reserve(this->trees_.size());
            for (unsigned int tree_idx = 0; tree_idx < this->trees_.size(); ++tree_idx) {
                hastings_ratios.push_back(op.propose(rng, this->trees_.at(tree_idx)));
            }
            this->compute_tree_partials();
            for (unsigned int tree_idx = 0; tree_idx < this->trees_.size(); ++tree_idx) {
                // Check to see if we updated a fixed parameter. If so, do
                // nothing and continue to next tree (to avoid counting toward
                // operator acceptance ratio).
                if ((op.target_parameter() == "coalescence rate") &&
                        (this->trees_.at(tree_idx).coalescence_rates_are_fixed())) {
                    ECOEVOLITY_ASSERT(! this->trees_.at(tree_idx).is_dirty());
                    continue;
                }
                if ((op.target_parameter() == "mutation rate") &&
                        (this->trees_.at(tree_idx).mutation_rates_are_fixed())) {
                    ECOEVOLITY_ASSERT(! this->trees_.at(tree_idx).is_dirty());
                    continue;
                }
                if ((op.target_parameter() == "node height multiplier") &&
                        (this->trees_.at(tree_idx).node_height_multiplier_is_fixed())) {
                    ECOEVOLITY_ASSERT(! this->trees_.at(tree_idx).is_dirty());
                    continue;
                }
                double likelihood_ratio =
                        this->trees_.at(tree_idx).get_log_likelihood_value() -
                        this->trees_.at(tree_idx).get_stored_log_likelihood_value();
                double prior_ratio =
                        this->trees_.at(tree_idx).get_log_prior_density_value() -
                        this->trees_.at(tree_idx).get_stored_log_prior_density_value();
                double hastings_ratio = hastings_ratios.at(tree_idx);
                double acceptance_probability =
                        likelihood_ratio + 
                        prior_ratio +
                        hastings_ratio;
                double u = rng.uniform_real();
                if (u < std::exp(acceptance_probability)) {
                    op.accept(this->operator_schedule_);
                }
                else {
                    op.reject(this->operator_schedule_);
                    this->trees_.at(tree_idx).restore_state();
                }
                this->trees_.at(tree_idx).make_clean();
                op.optimize(this->operator_schedule_, acceptance_probability);
            }
            this->compute_log_likelihood_and_prior(false);
        }
        else if (op.get_type() == Operator::OperatorTypeEnum::time_operator) {
            std::vector<double> hastings_ratios;
            hastings_ratios.reserve(this->node_heights_.size());
            for (unsigned int height_idx = 0; height_idx < this->node_heights_.size(); ++height_idx) {
                hastings_ratios.push_back(op.propose(rng, *(this->node_heights_.at(height_idx))));
            }
            this->make_trees_dirty();
            this->compute_tree_partials();
            for (unsigned int height_idx = 0; height_idx < this->node_heights_.size(); ++height_idx) {
                double old_lnl = 0.0;
                double new_lnl = 0.0;
                for (unsigned int tree_idx = 0; tree_idx < this->node_height_indices_.size(); ++tree_idx) {
                    if (this->node_height_indices_.at(tree_idx) == height_idx) {
                        old_lnl += this->trees_.at(tree_idx).get_stored_log_likelihood_value();
                        new_lnl += this->trees_.at(tree_idx).get_log_likelihood_value();
                    }
                }
                double likelihood_ratio = new_lnl - old_lnl;
                double prior_ratio =
                        this->node_heights_.at(height_idx)->relative_prior_ln_pdf() -
                        this->node_heights_.at(height_idx)->relative_prior_ln_pdf(
                                this->node_heights_.at(height_idx)->get_stored_value());
                double hastings_ratio = hastings_ratios.at(height_idx);
                double acceptance_probability =
                        likelihood_ratio + 
                        prior_ratio +
                        hastings_ratio;
                double u = rng.uniform_real();
                if (u < std::exp(acceptance_probability)) {
                    op.accept(this->operator_schedule_);
                }
                else {
                    op.reject(this->operator_schedule_);
                    this->node_heights_.at(height_idx)->restore();
                    for (unsigned int tree_idx = 0; tree_idx < this->node_height_indices_.size(); ++tree_idx) {
                        if (this->node_height_indices_.at(tree_idx) == height_idx) {
                            this->trees_.at(tree_idx).restore_likelihood();
                            // this->trees_.at(tree_idx).restore_prior_density();
                        }
                    }
                }
                op.optimize(this->operator_schedule_, acceptance_probability);
            }
            this->make_trees_clean();
            this->compute_log_likelihood_and_prior(false);
        }
        else if (op.get_type() == Operator::OperatorTypeEnum::model_operator) {
            double hastings_ratio = op.propose(rng, *this);
            this->compute_log_likelihood_and_prior(true);
            double likelihood_ratio = 
                this->log_likelihood_.get_value() -
                this->log_likelihood_.get_stored_value();
            double prior_ratio = 
                this->log_prior_density_.get_value() -
                this->log_prior_density_.get_stored_value();
            double acceptance_probability =
                    likelihood_ratio + 
                    prior_ratio +
                    hastings_ratio;
            double u = rng.uniform_real();
            if (u < std::exp(acceptance_probability)) {
                op.accept(this->operator_schedule_);
            }
            else {
                op.reject(this->operator_schedule_);
                this->restore_state();
            }
            this->make_trees_clean();
            op.optimize(this->operator_schedule_, acceptance_probability);
        }
        else {
            state_log_stream.close();
            operator_log_stream.close();
            throw EcoevolityError("unexpected operator");
        }

        if ((gen + 1) % sample_frequency == 0) {
            this->operator_schedule_.write_operator_rates(operator_log_stream);
            this->log_state(state_log_stream, gen + 1);
            this->log_state(std::cout, gen + 1, true);
        }

        // Check if the chain has gone astray
        ECOEVOLITY_DEBUG(
            if ((gen + 1) % (sample_frequency * 10) == 0) {
                double chain_ln_likelihood = this->log_likelihood_.get_value();
                this->make_trees_dirty();
                this->compute_log_likelihood_and_prior(true);
                double expected_ln_likelihood = this->log_likelihood_.get_value();
                if (fabs(chain_ln_likelihood - expected_ln_likelihood) > 1e-6) {
                    std::ostringstream message;
                    message << "ERROR: The likelihood at generation "
                            << gen
                            << " is " << chain_ln_likelihood
                            << "; expected " << expected_ln_likelihood
                            << "; last operator " << op.get_name();
                    state_log_stream.close();
                    operator_log_stream.close();
                    throw EcoevolityError(message.str());
                }
            }
        )
    }
    state_log_stream.close();
    operator_log_stream.close();
}
