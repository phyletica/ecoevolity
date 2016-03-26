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

double ComparisonPopulationTreeCollection::compute_log_prior_density_of_node_heights() {
}

void ComparisonPopulationTreeCollection::store_state() {
    this->log_likelihood_.store();
    this->log_prior_density_.store();
    for (auto tree_iter : this->trees_) {
        tree_iter->store_state();
    }
}
void ComparisonPopulationTreeCollection::restore_state() {
    this->log_likelihood_.restore();
    this->log_prior_density_.restore();
    for (auto tree_iter : this->trees_) {
        tree_iter->restore_state();
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
        lnl += tree_iter->get_log_likelihood_value();
        lnp += tree_iter->get_log_prior_density_value();
    }
    for (auto height_iter : this->node_heights_) {
        lnp += this->node_heights_->relative_prior_ln_pdf();
    }

    this->log_likelihood_.set_value(lnl);
    this->log_prior_density_.set_value(lnp);
}

void ComparisonPopulationTreeCollection::compute_tree_partials() {
    for (auto tree_iter : this->trees_) {
        tree_iter->compute_log_likelihood_and_prior();
    }
}

void ComparisonPopulationTreeCollection::compute_tree_partials_threaded() {
    std::vector<std::thread> threads(this->trees_.size());

    for (auto tree_iter : this->trees_) {
        threads.push_back(std::thread(tree_iter->compute_log_likelihood_and_prior));
    }

    for (auto &t : threads) {
        t.join();
    }
}

void ComparisonPopulationTreeCollection::make_trees_clean() {
    for (auto tree_iter : this->trees_) {
        tree_iter->make_clean();
    }
}

void ComparisonPopulationTreeCollection::mcmc() {
    // TODO:  initialize state
    for (unsigned int gen = 0; gen < ngens; ++gen) {
        this->store_state();
        Operator& op = this->operator_schedule_.draw_operator(this->rng_);
        if (op.get_target_type() == Operator::TargetTypeEnum::ComparisonPopulationTree) {
            std::vector<double> hastings_ratios(this->trees_.size());
            for (unsigned int tree_idx = 0; tree_idx < this->trees_.size(); ++ tree_idx) {
                hastings_ratios.push_back(op.propose(this->rng_, this->trees_.at(tree_idx)));
            }
            this->compute_tree_partials();
            for (unsigned int tree_idx = 0; tree_idx < this->trees_.size(); ++ tree_idx) {
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
                double u = this->rng_.uniform_real();
                if (u < acceptance_probability) {
                    op.accept();
                }
                else {
                    op.reject();
                    this->trees_.at(tree_idx).restore_state();
                }
                this->trees_.at(tree_idx).make_clean();
            }
            this->compute_log_likelihood_and_prior(false);
        }
        else if (op.get_target_type() == Operator::TargetTypeEnum::ComparisonPopulationTreeCollection) {
            double hastings_ratio = op.propose(this->rng_, this);
            this->compute_log_likelihood_and_prior(true);
            double likelihood_ratio = 
                this->log_likelihood_.get_value() -
                this->log_likelihood_.get_store_value();
            double prior_ratio = 
                this->log_prior_density_.get_value() -
                this->log_prior_density_.get_store_value();
            double acceptance_probability =
                    likelihood_ratio + 
                    prior_ratio +
                    hastings_ratio;
            double u = this->rng_.uniform_real();
            if (u < acceptance_probability) {
                op.accept();
            }
            else {
                op.reject();
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

