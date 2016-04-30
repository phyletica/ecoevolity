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
    for (auto tree_iter : this->trees_) {
        tree_iter.store_state();
    }
}
void ComparisonPopulationTreeCollection::restore_state() {
    this->log_likelihood_.restore();
    this->log_prior_density_.restore();
    for (auto tree_iter : this->trees_) {
        tree_iter.restore_state();
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
    for (auto tree_iter : this->trees_) {
        lnl += tree_iter.get_log_likelihood_value();
        lnp += tree_iter.get_log_prior_density_value();
    }
    for (auto height_iter : this->node_heights_) {
        lnp += height_iter->relative_prior_ln_pdf();
    }

    this->log_likelihood_.set_value(lnl);
    this->log_prior_density_.set_value(lnp);
}

void ComparisonPopulationTreeCollection::compute_tree_partials() {
    for (auto tree_iter : this->trees_) {
        tree_iter.compute_log_likelihood_and_prior();
    }
}

void ComparisonPopulationTreeCollection::compute_tree_partials_threaded() {
    this->compute_tree_partials();
    // TODO: get multithreading working
    // std::vector<std::thread> threads;
    // threads.reserve(this->trees_.size());

    // for (auto tree_iter : this->trees_) {
    //     threads.push_back(std::thread(
    //                 &ComparisonPopulationTree::compute_log_likelihood_and_prior,
    //                 tree_iter));
    // }

    // for (auto &t : threads) {
    //     t.join();
    // }
}

void ComparisonPopulationTreeCollection::make_trees_clean() {
    for (auto tree_iter : this->trees_) {
        tree_iter.make_clean();
    }
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
    // TODO:
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
    std::shared_ptr<PositiveRealParameter> new_height = std::make_shared<PositiveRealParameter>(this->node_height_prior_, height);
    this->node_heights_.push_back(new_height);
    for (auto tree_idx : mapped_tree_indices) {
        this->node_height_indices_.at(tree_idx) = this->node_heights_.size() - 1;
        this->trees_.at(tree_idx).set_height_parameter(new_height);
    }
}

void ComparisonPopulationTreeCollection::mcmc(
        RandomNumberGenerator& rng,
        unsigned int chain_length) {
    // TODO:  initialize state
    std::shared_ptr<Operator> op;
    for (unsigned int gen = 0; gen < chain_length; ++gen) {
        this->store_state();
        op = this->operator_schedule_.draw_operator(rng);
        if (op->get_type() == Operator::OperatorTypeEnum::tree_operator) {
            std::vector<double> hastings_ratios;
            hastings_ratios.reserve(this->trees_.size());
            for (unsigned int tree_idx = 0; tree_idx < this->trees_.size(); ++tree_idx) {
                hastings_ratios.push_back(op->propose(rng, this->trees_.at(tree_idx)));
            }
            this->compute_tree_partials();
            for (unsigned int tree_idx = 0; tree_idx < this->trees_.size(); ++tree_idx) {
                // Check to see if we updated a fixed parameter. If so, do
                // nothing and continue to next tree (to avoid counting toward
                // operator acceptance ratio).
                if ((op->target_parameter() == "coalescence rate") &&
                        (this->trees_.at(tree_idx).coalescence_rates_are_fixed())) {
                    ECOEVOLITY_ASSERT(! this->trees_.at(tree_idx).is_dirty());
                    continue;
                }
                if ((op->target_parameter() == "mutation rate") &&
                        (this->trees_.at(tree_idx).mutation_rates_are_fixed())) {
                    ECOEVOLITY_ASSERT(! this->trees_.at(tree_idx).is_dirty());
                    continue;
                }
                if ((op->target_parameter() == "node height multiplier") &&
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
                if (u < acceptance_probability) {
                    op->accept();
                }
                else {
                    op->reject();
                    this->trees_.at(tree_idx).restore_state();
                }
                this->trees_.at(tree_idx).make_clean();
            }
            this->compute_log_likelihood_and_prior(false);
        }
        else if (op->get_type() == Operator::OperatorTypeEnum::time_operator) {
            std::vector<double> hastings_ratios;
            hastings_ratios.reserve(this->node_heights_.size());
            for (unsigned int height_idx = 0; height_idx < this->node_heights_.size(); ++height_idx) {
                hastings_ratios.push_back(op->propose(rng, *(this->node_heights_.at(height_idx))));
            }
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
                if (u < acceptance_probability) {
                    op->accept();
                }
                else {
                    op->reject();
                    this->node_heights_.at(height_idx)->restore();
                    for (unsigned int tree_idx = 0; tree_idx < this->node_height_indices_.size(); ++tree_idx) {
                        if (this->node_height_indices_.at(tree_idx) == height_idx) {
                            this->trees_.at(tree_idx).restore_likelihood();
                            // this->trees_.at(tree_idx).restore_prior_density();
                        }
                    }
                }
            }
            this->make_trees_clean();
            this->compute_log_likelihood_and_prior(false);
        }
        else if (op->get_type() == Operator::OperatorTypeEnum::model_operator) {
            // TODO: for DPP nothing to accept/reject
            double hastings_ratio = op->propose(rng, *this);
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
            if (u < acceptance_probability) {
                op->accept();
            }
            else {
                op->reject();
                this->restore_state();
            }
            this->make_trees_clean();
        }
        else {
            throw EcoevolityError("unexpected operator");
        }
        // TODO: check if its a logging generation
    }
}
