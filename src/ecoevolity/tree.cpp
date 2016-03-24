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

#include "tree.hpp"

PopulationTree::PopulationTree(
        const std::string path, 
        const char population_name_delimiter,
        const bool population_name_is_prefix,
        const bool genotypes_are_diploid,
        const bool markers_are_dominant,
        const bool validate) {
    this->init(path,
               population_name_delimiter,
               population_name_is_prefix,
               genotypes_are_diploid,
               markers_are_dominant,
               validate);
}

void PopulationTree::init(
        const std::string path, 
        const char population_name_delimiter,
        const bool population_name_is_prefix,
        const bool genotypes_are_diploid,
        const bool markers_are_dominant,
        const bool validate) {
    this->data_.init(
            path,
            population_name_delimiter,
            population_name_is_prefix,
            genotypes_are_diploid,
            markers_are_dominant,
            validate);
    if (this->data_.get_number_of_populations() < 1) {
        throw EcoevolityError("PopulationTree(); no populations were found");
    }
    unsigned int number_of_missing_patterns_removed = this->data_.remove_missing_population_patterns();
    if (this->correct_for_constant_patterns_) {
        unsigned int number_of_constant_patterns_removed = this->data_.remove_constant_patterns();
    }

    this->init_tree();

    this->root_->resize_all();

    this->pattern_likelihoods_.assign(this->data_.get_number_of_patterns(), 0.0);
}

void PopulationTree::init_tree() {
    if (this->data_.get_number_of_populations() < 3) {
        this->root_ = new PopulationNode(0.0);
        for (unsigned int pop_idx = 0;
                pop_idx < this->data_.get_number_of_populations();
                ++pop_idx) {
            PopulationNode * tip = new PopulationNode(
                    this->data_.get_population_label(pop_idx),
                    0.0,
                    this->data_.get_max_allele_count(pop_idx));
            tip->fix_node_height();
            this->root_->add_child(tip);
        }
        return;
    }
    PopulationNode * ancestor = new PopulationNode(0.0);
    ancestor->add_child(new PopulationNode(this->data_.get_population_label(0)));
    ancestor->add_child(new PopulationNode(this->data_.get_population_label(1)));
    for (unsigned int pop_idx = 2;
            pop_idx < this->data_.get_number_of_populations();
            ++pop_idx) {
        PopulationNode * next_ancestor = new PopulationNode(0.0);
        next_ancestor->add_child(ancestor);
        PopulationNode * tip = new PopulationNode(
                this->data_.get_population_label(pop_idx),
                0.0,
                this->data_.get_max_allele_count(pop_idx));
        tip->fix_node_height();
        next_ancestor->add_child(tip);
        ancestor = next_ancestor;
    }
    this->root_ = ancestor;
}

void PopulationTree::compute_leaf_partials(
        int pattern_index,
        PopulationNode * node) {
    unsigned int pop_idx = this->data_.get_population_index(node->get_label());
    unsigned int allele_count = 0;
    unsigned int red_allele_count = 0;
    if (pattern_index == -1) {
        allele_count = this->data_.get_max_allele_count(pop_idx);
        red_allele_count = 0;
    }
    else if (pattern_index == -2) {
        allele_count = this->data_.get_max_allele_count(pop_idx);
        red_allele_count = allele_count;
    }
    else if (pattern_index > -1) {
        allele_count = this->data_.get_allele_count(pattern_index, pop_idx);
        red_allele_count = this->data_.get_red_allele_count(pattern_index, pop_idx);
    }
    else {
        throw EcoevolityError("PopulationTree::compute_leaf_partials(): Unexpected negative pattern index");
    }
    if (this->data_.markers_are_dominant()) {
        unsigned int n = allele_count;
        unsigned int n_reds = red_allele_count;
        allele_count = allele_count * 2;
        if (red_allele_count > 0) {
            BiallelicPatternProbabilityMatrix m(allele_count);
            double p_r_k_n = 1.0;
            for (unsigned int r = 1; r <= n_reds; ++r) {
                p_r_k_n = (p_r_k_n * 2.0 * (n - r + 1.0)) / ((2.0 * n) - r + 1.0);
            }
            for (unsigned int k = n_reds; k <= (2 * n_reds); ++k) {
                if (k > n_reds) {
                    p_r_k_n = (p_r_k_n * ((2.0 * n_reds) - k + 1) * k) /
                              (2.0 * ( k - n_reds) * ((2.0 * n) - k + 1.0));
                }
                m.set_pattern_probability(allele_count, k, p_r_k_n);
            }
            node->copy_bottom_pattern_probs(m);
            return;
        }
        BiallelicPatternProbabilityMatrix m(allele_count, n_reds);
        node->copy_bottom_pattern_probs(m);
        return;
    }
    BiallelicPatternProbabilityMatrix m(allele_count, red_allele_count);
    node->copy_bottom_pattern_probs(m);
    return;
}

void PopulationTree::compute_top_of_branch_partials(
        PopulationNode * node) {
    if (node->get_allele_count() == 0) {
        node->copy_top_pattern_probs(node->get_bottom_pattern_probs());
        return;
    }

    BiallelicPatternProbabilityMatrix m = matrix_exponentiator.expQTtx(
            node->get_allele_count(),
            this->u_->get_value(),
            this->v_->get_value(),
            node->get_coalescence_rate(),
            node->get_length() * this->node_height_multiplier_->get_value(),
            node->get_bottom_pattern_probs());
    node->copy_top_pattern_probs(m);
}

void PopulationTree::compute_internal_partials(
        PopulationNode * node) {
    if (node->get_number_of_children() == 1) {
        node->copy_bottom_pattern_probs(node->get_child(0)->get_top_pattern_probs());
        return;
    }
    if (node->get_child(0)->get_allele_count() == 0) {
        node->copy_bottom_pattern_probs(node->get_child(1)->get_top_pattern_probs());
        return;
    }
    if (node->get_child(1)->get_allele_count() == 0) {
        node->copy_bottom_pattern_probs(node->get_child(0)->get_top_pattern_probs());
        return;
    }
    unsigned int allele_count_child1 = node->get_child(0)->get_allele_count();
    unsigned int allele_count_child2 = node->get_child(1)->get_allele_count();

    std::vector<double> pattern_probs_child1 = node->get_child(0)->get_top_pattern_probs().get_pattern_prob_matrix();
    std::vector<double> pattern_probs_child2 = node->get_child(1)->get_top_pattern_probs().get_pattern_prob_matrix();

    for (unsigned int n = 1; n <= allele_count_child1; ++n) {
        double b_nr = 1.0;
        for (unsigned int r = 0; r <= n; ++r) {
            pattern_probs_child1.at(((n*(n+1))/2)-1+r) *= b_nr;
            b_nr *= ((double)n - r)/(r+1);
        }
    }

    for (unsigned int n = 1; n<= allele_count_child2; ++n) {
        double b_nr = 1.0;
        for (unsigned int r = 0; r <= n; ++r) {
            pattern_probs_child2.at(((n*(n+1))/2)-1+r) *= b_nr;
            b_nr *= ((double)n - r)/(r+1);
        }
    }

    unsigned int allele_count = allele_count_child1 + allele_count_child2;
    std::vector<double> pattern_probs; 
    pattern_probs.assign(
            ((((allele_count + 1) * (allele_count + 2))/2) - 1),
            0);
    for (unsigned int n1 = 1; n1 <= allele_count_child1; ++n1) {
        for (unsigned int r1 = 0; r1 <= n1; ++r1) {
            double f11 = pattern_probs_child1.at(n1*(n1+1)/2-1+r1);
            for (unsigned int n2 = 1; n2 <= allele_count_child2; ++n2) {
                for (unsigned int r2 = 0; r2 <= n2; ++r2) {
                    pattern_probs.at((n1+n2)*(n1+n2+1)/2-1+(r1+r2)) += f11 * pattern_probs_child2.at(n2*(n2+1)/2-1+r2);
                }
            }
        }
    }

    for (unsigned int n = 1; n <= allele_count; ++n) {
        double b_nr = 1.0;
        for (unsigned int r = 0; r <= n; ++r) {
            double f_nr = pattern_probs.at(n*(n+1)/2-1+r);
            f_nr /= b_nr;
            // TODO: better way to fix this?
            f_nr = std::max(f_nr, 0.0);
            pattern_probs.at(n*(n+1)/2-1+r) = f_nr;
            b_nr *= ((double)n - r)/(r+1);

        }
    }
    BiallelicPatternProbabilityMatrix m(allele_count, pattern_probs);
    node->copy_bottom_pattern_probs(m);
}

void PopulationTree::compute_pattern_partials(
        int pattern_index,
        PopulationNode * node) {
    if (node->is_leaf()) {
        this->compute_leaf_partials(pattern_index, node);
    }
    else if (node->get_number_of_children() == 1) {
        compute_pattern_partials(pattern_index, node->get_child(0));
        compute_top_of_branch_partials(node->get_child(0));
        compute_internal_partials(node);
    }
    else if (node->get_number_of_children() == 2) {
        compute_pattern_partials(pattern_index, node->get_child(0));
        compute_pattern_partials(pattern_index, node->get_child(1));
        compute_top_of_branch_partials(node->get_child(0));
        compute_top_of_branch_partials(node->get_child(1));
        compute_internal_partials(node);
    }
    else {
        throw EcoevolityError(
            "PopulationTree::compute_pattern_probability(); unexpected number of children");
    }
}

std::vector< std::vector<double> > PopulationTree::compute_root_probabilities() {
    unsigned int N = this->root_->get_allele_count();
    std::vector< std::vector<double> > x (N + 1); 
    QMatrix q = QMatrix(
            N,
            this->u_->get_value(),
            this->v_->get_value(),
            this->root_->get_coalescence_rate());
    std::vector<double> xcol = q.find_orthogonal_vector();

    // ECOEVOLITY_DEBUG(
    //     std::cerr << "xcol = [";
    //     for (unsigned int i = 0; i < xcol.size(); ++i) {
    //         std::cerr << xcol.at(i) << " ";
    //     }
    //     std::cerr << "]" << std::endl;
    // )

    unsigned int index = 1;
    for (unsigned int n = 1; n <= N; ++n) {
        x.at(n).resize(n + 1, 0.0);
        double row_sum = 0.0;
        for (unsigned int r = 0; r <= n; ++r) {
            double xcol_index = std::max(xcol.at(index), 0.0);
            row_sum += xcol_index;
            x.at(n).at(r) = xcol_index;
            index++;
        }
        for (unsigned int r = 0; r <= n; ++r) {
            x.at(n).at(r) = x.at(n).at(r) / row_sum;
        }
    }
    return x;
}

double PopulationTree::compute_root_likelihood() {
    unsigned int N = this->root_->get_allele_count();
    std::vector< std::vector<double> > conditionals = this->compute_root_probabilities();

    // ECOEVOLITY_DEBUG(
    //     for (unsigned int n = 1; n <= N; ++n) {
    //         for (unsigned int r = 0; r <= n; ++r) {
    //             std::cerr << "root height: " << this->root_->get_height() << std::endl;
    //             std::cerr << "conditional[" << n << ", " << r << "] = " << conditionals.at(n).at(r) << std::endl;
    //             std::cerr << "bottom_pattern_probs[" << n << ", " << r << "] = " << this->root_->get_bottom_pattern_probability(n, r) << std::endl;
    //         }
    //     }
    // )

    double sum = 0.0;
    for (unsigned int n = 1; n <= N; ++n) {
        for (unsigned int r = 0; r <= n; ++r) {
            double term = conditionals.at(n).at(r) * this->root_->get_bottom_pattern_probability(n, r);
            sum += term;
            if (sum < 0.0) {
                throw EcoevolityError("PopulationTree::compute_root_likelihood(): Numerical error");
            }
        }
    }

    // ECOEVOLITY_DEBUG(
    //     std::cerr << "root likelihood: " << sum << std::endl;
    // )
    return sum;
}

double PopulationTree::compute_pattern_likelihood(int pattern_index) {
    this->compute_pattern_partials(pattern_index, this->root_);
    return this->compute_root_likelihood();
}

void PopulationTree::compute_pattern_likelihoods() {
    for (unsigned int pattern_idx = 0;
            pattern_idx < this->data_.get_number_of_patterns();
            ++pattern_idx) {
        this->pattern_likelihoods_.at(pattern_idx) = this->compute_pattern_likelihood(pattern_idx);
    }
    double all_green_pattern_likelihood = this->compute_pattern_likelihood(-1);
    double all_red_pattern_likelihood = all_green_pattern_likelihood;
    if (! this->mutation_rates_are_constrained()) {
        all_red_pattern_likelihood = this->compute_pattern_likelihood(-2);
    }
    this->all_green_pattern_likelihood_.set_value(all_green_pattern_likelihood);
    this->all_red_pattern_likelihood_.set_value(all_red_pattern_likelihood);
}

void PopulationTree::calculate_likelihood_correction() {
    double log_likelihood_correction = 0.0;
    for (unsigned int pattern_idx = 0;
            pattern_idx < this->data_.get_number_of_patterns();
            ++pattern_idx) {
        for (unsigned int pop_idx = 0;
                pop_idx < this->data_.get_number_of_populations();
                ++pop_idx) {
            log_likelihood_correction -= (
                    this->calculate_log_binomial(
                        this->data_.get_red_allele_count(pattern_idx, pop_idx),
                        this->data_.get_allele_count(pattern_idx, pop_idx)) *
                    this->data_.get_pattern_weight(pattern_idx)
                    );
        }
    }
    this->log_likelihood_correction_.set_value(log_likelihood_correction);
    this->likelihood_correction_was_calculated_ = true;
    // ECOEVOLITY_DEBUG(
    //     std::cerr << "Log likelihood correction: " << this->log_likelihood_correction_ << std::endl;
    // )
}

double PopulationTree::get_likelihood_correction(bool force) {
    if ((! this->likelihood_correction_was_calculated_) || (force)) {
        this->calculate_likelihood_correction();
    }
    return this->log_likelihood_correction_.get_value();
}

double PopulationTree::calculate_log_binomial(
        unsigned int red_allele_count,
        unsigned int allele_count) const {
    double f = 0.0;
    for (unsigned int i = red_allele_count + 1; i <= allele_count; ++i) {
        f += std::log(i) - std::log(allele_count - i + 1);
    }
    return f;
}

bool PopulationTree::constant_site_counts_were_provided() {
    if ((this->number_of_constant_green_sites_ > -1) && (this->number_of_constant_red_sites_ > -1)) {
        return true;
    }
    return false;
}

double PopulationTree::compute_log_likelihood() {
    double log_likelihood = 0.0;
    this->compute_pattern_likelihoods();
    for (unsigned int pattern_idx = 0;
            pattern_idx < this->data_.get_number_of_patterns();
            ++pattern_idx) {
        double pattern_likelihood = this->pattern_likelihoods_.at(pattern_idx);
        double weight = (double) this->data_.get_pattern_weight(pattern_idx);
        if (pattern_likelihood ==  0.0) {
            log_likelihood = -10e100;
            break;
        }
        log_likelihood += weight * std::log(pattern_likelihood);
    }

    if (this->correct_for_constant_patterns_) {
        if (this->constant_site_counts_were_provided()) {
            double constant_log_likelihood =
                    ((double)this->number_of_constant_green_sites_ * std::log(this->all_green_pattern_likelihood_.get_value())) +
                    ((double)this->number_of_constant_red_sites_ * std::log(this->all_red_pattern_likelihood_.get_value()));
            log_likelihood += constant_log_likelihood;
        }
        else if (this->use_removed_constant_site_counts_){
            double constant_log_likelihood =
                    ((double)this->data_.get_number_of_constant_green_sites_removed() *
                    std::log(this->all_green_pattern_likelihood_.get_value())) +
                    ((double)this->data_.get_number_of_constant_red_sites_removed() *
                    std::log(this->all_red_pattern_likelihood_.get_value()));
            log_likelihood += constant_log_likelihood;
        }
        else {
            log_likelihood -= ((double)this->data_.get_number_of_sites() * 
                    std::log(1.0 - this->all_green_pattern_likelihood_.get_value() -
                            this->all_red_pattern_likelihood_.get_value()));
        }
    }


    if (this->correct_for_full_likelihood_) {
        log_likelihood += this->get_likelihood_correction();
    }

    // ECOEVOLITY_DEBUG(
    //     std::cerr << "PopulationTree::compute_log_likelihood(): " << log_likelihood << std::endl;
    // )

    this->log_likelihood_.set_value(log_likelihood);
    ++this->number_of_likelihood_calculations_;
    return log_likelihood;
}

double PopulationTree::get_log_likelihood_value() const {
    return this->log_likelihood_.get_value();
}
double PopulationTree::get_stored_log_likelihood_value() const {
    return this->log_likelihood_.get_stored_value();
}

void PopulationTree::fold_patterns() {
    if (! this->mutation_rates_are_constrained()) {
        std::cerr << 
            "WARNING: Site patterns are being folded when foward/backward\n" <<
            "         mutation rates are not constrained." << std::endl;
    }
    this->data_.fold_patterns();
    this->make_dirty();
}

void PopulationTree::set_root_height(double height) {
    this->root_->set_height(height);
}
void PopulationTree::update_root_height(double height) {
    this->root_->update_height(height);
}
const double& PopulationTree::get_root_height() const {
    return this->root_->get_height();
}
void PopulationTree::store_root_height() {
    this->root_->store_height();
}
void PopulationTree::restore_root_height() {
    this->root_->restore_height();
}
void PopulationTree::set_root_height_parameter(PositiveRealParameter * h) {
    this->root_->set_height_parameter(h);
}
PositiveRealParameter * PopulationTree::get_root_height_parameter() const {
    return this->root_->get_height_parameter();
}

void PopulationTree::set_u(double u) {
    this->u_->set_value(u);
    this->make_dirty();
}
void PopulationTree::update_u(double u) {
    this->u_->update_value(u);
    this->make_dirty();
}
void PopulationTree::set_v(double v) {
    this->v_->set_value(v);
    this->make_dirty();
}
void PopulationTree::update_v(double v) {
    this->v_->update_value(v);
    this->make_dirty();
}
const double& PopulationTree::get_u() const {
    return this->u_->get_value();
}
const double& PopulationTree::get_v() const {
    return this->v_->get_value();
}
void PopulationTree::store_u() {
    this->u_->store();
}
void PopulationTree::store_v() {
    this->v_->store();
}
void PopulationTree::restore_u() {
    this->u_->restore();
    this->make_dirty();
}
void PopulationTree::restore_v() {
    this->v_->restore();
    this->make_dirty();
}

void PopulationTree::set_u_parameter(PositiveRealParameter * u) {
    this->u_ = u;
    this->make_dirty();
}
void PopulationTree::set_v_parameter(PositiveRealParameter * v) {
    this->v_ = v;
    this->make_dirty();
}
PositiveRealParameter * PopulationTree::get_u_parameter() const {
    return this->u_;
}
PositiveRealParameter * PopulationTree::get_v_parameter() const {
    return this->v_;
}

void PopulationTree::set_root_coalescence_rate(double rate) {
    this->root_->set_coalescence_rate(rate);
}
void PopulationTree::set_coalescence_rate(double rate) {
    this->root_->set_all_coalescence_rates(rate);
}
double PopulationTree::get_root_coalescence_rate() const {
    return this->root_->get_coalescence_rate();
}

void PopulationTree::store_state() {
    this->store_likelihood();
    this->store_prior_density();
    this->store_parameters();
}
void PopulationTree::store_likelihood() {
    this->log_likelihood_.store();
    this->all_green_pattern_likelihood_.store();
    this->all_red_pattern_likelihood_.store();
}
void PopulationTree::store_prior_density() {
    this->log_prior_density_.store();
}
void PopulationTree::store_parameters() {
    this->store_u();
    this->store_v();
    this->store_all_coalescence_rates();
    this->store_all_heights();
}
void PopulationTree::store_all_coalescence_rates() {
    this->root_->store_all_coalescence_rates();
}
void PopulationTree::store_all_heights() {
    this->root_->store_all_heights();
}

void PopulationTree::restore_state() {
    this->restore_likelihood();
    this->restore_prior_density();
    this->restore_parameters();
}
void PopulationTree::restore_likelihood() {
    this->log_likelihood_.restore();
    this->all_green_pattern_likelihood_.restore();
    this->all_red_pattern_likelihood_.restore();
}
void PopulationTree::restore_prior_density() {
    this->log_prior_density_.restore();
}
void PopulationTree::restore_parameters() {
    this->restore_u();
    this->restore_v();
    this->restore_all_coalescence_rates();
    this->restore_all_heights();
}
void PopulationTree::restore_all_coalescence_rates() {
    this->root_->restore_all_coalescence_rates();
}
void PopulationTree::restore_all_heights() {
    this->root_->restore_all_heights();
}

void PopulationTree::set_node_height_prior(ContinuousProbabilityDistribution * prior) {
    this->node_height_prior_ = prior;
    this->root_->set_all_node_height_priors(prior);
}

void PopulationTree::set_population_size_prior(ContinuousProbabilityDistribution * prior) {
    delete this->population_size_prior_;
    this->population_size_prior_ = prior;
    this->root_->set_all_population_size_priors(prior);
}

void PopulationTree::set_u_prior(ContinuousProbabilityDistribution * prior) {
    delete this->u_prior_;
    this->u_prior_ = prior;
    this->u_->set_prior(prior);
    this->make_dirty();
}
void PopulationTree::set_v_prior(ContinuousProbabilityDistribution * prior) {
    delete this->v_prior_;
    this->v_prior_ = prior;
    this->v_->set_prior(prior);
    this->make_dirty();
}

double PopulationTree::compute_log_prior_density() {
    double d = 0.0;
    d += this->compute_log_prior_density_of_mutation_rates();
    d += this->compute_log_prior_density_of_node_heights();
    d += this->compute_log_prior_density_of_coalescence_rates();
    this->log_prior_density_.set_value(d);
    return d;
}
double PopulationTree::compute_log_prior_density_of_mutation_rates() const {
    double d = this->u_->relative_prior_ln_pdf();
    if (! this->mutation_rates_are_constrained()) {
        d += this->v_->relative_prior_ln_pdf();
    }
    return d;
}
double PopulationTree::compute_log_prior_density_of_node_heights() const {
    return this->root_->calculate_ln_relative_node_height_prior_density();
}
double PopulationTree::compute_log_prior_density_of_coalescence_rates() const {
    return this->root_->calculate_ln_relative_coalescence_rate_prior_density();
}

double PopulationTree::get_log_prior_density_value() const {
    return this->log_prior_density_.get_value();
}
double PopulationTree::get_stored_log_prior_density_value() const {
    return this->log_prior_density_.get_stored_value();
}

bool PopulationTree::is_dirty() const {
    if (this->is_dirty_) {
        return true;
    }
    return this->root_->clade_has_dirt();
}

void PopulationTree::make_dirty() {
    this->is_dirty_ = true;
}
void PopulationTree::make_clean() {
    this->is_dirty_ = false;
    this->root_->make_all_clean();
}


ComparisonPopulationTree::ComparisonPopulationTree(
        const std::string path, 
        const char population_name_delimiter,
        const bool population_name_is_prefix,
        const bool genotypes_are_diploid,
        const bool markers_are_dominant,
        const bool validate) {
    this->init(path,
               population_name_delimiter,
               population_name_is_prefix,
               genotypes_are_diploid,
               markers_are_dominant,
               validate);
    if (this->data_.get_number_of_populations() > 2) {
        throw EcoevolityError("ComparisonPopulationTree(); does not support more than 2 populations");
    }
    this->root_->set_label("root-" + this->root_->get_child(0)->get_label());
}

void ComparisonPopulationTree::set_child_coalescence_rate(
        unsigned int child_index,
        double rate) {
    this->root_->get_child(child_index)->set_coalescence_rate(rate);
}
void ComparisonPopulationTree::update_child_coalescence_rate(
        unsigned int child_index,
        double rate) {
    this->root_->get_child(child_index)->update_coalescence_rate(rate);
}
const double& ComparisonPopulationTree::get_child_coalescence_rate(
        unsigned int child_index) const {
    return this->root_->get_child(child_index)->get_coalescence_rate();
}
void ComparisonPopulationTree::store_child_coalescence_rate(
        unsigned int child_index) {
    this->root_->get_child(child_index)->store_coalescence_rate();
}
void ComparisonPopulationTree::restore_child_coalescence_rate(
        unsigned int child_index) {
    this->root_->get_child(child_index)->restore_coalescence_rate();
}
void ComparisonPopulationTree::set_child_coalescence_rate_parameter(
        unsigned int child_index,
        CoalescenceRateParameter * r) {
    this->root_->get_child(child_index)->set_coalescence_rate_parameter(r);
}
CoalescenceRateParameter * ComparisonPopulationTree::get_child_coalescence_rate_parameter(
        unsigned int child_index) const {
    return this->root_->get_child(child_index)->get_coalescence_rate_parameter();
}

// Node eight sharing needs to be dealt with in next level up in
// class hierarchy
double ComparisonPopulationTree::compute_log_prior_density() {
    double d = 0.0;
    d += this->compute_log_prior_density_of_mutation_rates();
    // d += this->compute_log_prior_density_of_node_heights();
    d += this->compute_log_prior_density_of_coalescence_rates();
    this->log_prior_density_.set_value(d);
    return d;
}
