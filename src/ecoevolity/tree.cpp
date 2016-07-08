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
        std::string path, 
        char population_name_delimiter,
        bool population_name_is_prefix,
        bool genotypes_are_diploid,
        bool markers_are_dominant,
        bool constant_sites_removed,
        bool validate,
        bool strict_on_constant_sites,
        bool strict_on_missing_sites) {
    this->init(path,
               population_name_delimiter,
               population_name_is_prefix,
               genotypes_are_diploid,
               markers_are_dominant,
               constant_sites_removed,
               validate,
               strict_on_constant_sites,
               strict_on_missing_sites);
}

void PopulationTree::init(
        std::string path, 
        char population_name_delimiter,
        bool population_name_is_prefix,
        bool genotypes_are_diploid,
        bool markers_are_dominant,
        bool constant_sites_removed,
        bool validate,
        bool strict_on_constant_sites,
        bool strict_on_missing_sites) {
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
    if (number_of_missing_patterns_removed > 0) {
        if (strict_on_missing_sites) {
            std::ostringstream message;
            message << "\n#######################################################################\n"
                    <<   "###############################  ERROR  ###############################\n"
                    << this->data_.get_number_of_missing_sites_removed()
                    << " sites from the alignment in:\n    \'"
                    << path << "\'\n"
                    << "have no data for at least one population.\n"
                    << "#######################################################################\n";
            throw EcoevolityMissingDataError(message.str(), path);
        }
        else {
            std::ostringstream message;
            message << "\n#######################################################################\n"
                    <<   "##############################  WARNING  ##############################\n"
                    << this->data_.get_number_of_missing_sites_removed()
                    << " sites will be ignored from the alignment in:\n    \'"
                    << path << "\'\n"
                    << "due to at least one population with no data.\n"
                    << "#######################################################################\n";
            std::cerr << message.str() << std::endl;
        }
    }
    this->constant_sites_removed_ = constant_sites_removed;
    if (this->constant_sites_removed_) {
        // Have to make sure there are no missing sites
        unsigned int number_of_constant_patterns_removed = this->data_.remove_constant_patterns();
        if (number_of_constant_patterns_removed > 0) {
            if (strict_on_constant_sites) {
                std::ostringstream message;
                message << "\n#######################################################################\n"
                        <<   "###############################  ERROR  ###############################\n"
                        << this->data_.get_number_of_constant_sites_removed()
                        << " constant sites were found in the alignment in:\n"
                        << "    \'" << path << "\'\n"
                        << "but you indicated that such sites were already removed with option:\n"
                        << "    constant_sites_removed = true\n"
                        << "If you intended to remove them, please do so and re-run the analysis.\n"
                        << "If you intended for constant sites to be used in the likelihood\n"
                        << "calculations, you should set \'constant_sites_removed\' to false for\n"
                        << "this alignment and re-run the analysis.\n"
                        << "#######################################################################\n";
                throw EcoevolityConstantSitesError(message.str(), path);
            }
            else {
                std::ostringstream message;
                message << "\n#######################################################################\n"
                        <<   "##############################  WARNING  ##############################\n"
                        << this->data_.get_number_of_constant_sites_removed()
                        << " constant sites were found in the alignment in:\n"
                        << "    \'" << path << "\'\n"
                        << "but you indicated that such sites were already removed with option:\n"
                        << "    constant_sites_removed = true\n"
                        << "These sites have been removed, so if you intended to remove them, but\n"
                        << "missed them, all is well. However, if you intended for the constant\n"
                        << "sites to be used in the likelihood calculations, you should set\n"
                        << "\'constant_sites_removed\' to false for this alignment and re-run this\n"
                        << "analysis.\n"
                        << "#######################################################################\n";
                std::cerr << message.str() << std::endl;
            }
        }
    }

    this->init_tree();

    this->root_->resize_all();

    this->pattern_likelihoods_.assign(this->data_.get_number_of_patterns(), 0.0);
}

void PopulationTree::init_tree() {
    if (this->data_.get_number_of_populations() < 3) {
        this->root_ = std::make_shared<PopulationNode>(0.0);
        for (unsigned int pop_idx = 0;
                pop_idx < this->data_.get_number_of_populations();
                ++pop_idx) {
            std::shared_ptr<PopulationNode> tip = std::make_shared<PopulationNode>(
                    this->data_.get_population_label(pop_idx),
                    0.0,
                    this->data_.get_max_allele_count(pop_idx));
            tip->fix_node_height();
            this->root_->add_child(tip);
        }
        return;
    }
    std::shared_ptr<PopulationNode> ancestor = std::make_shared<PopulationNode>(0.0);
    ancestor->add_child(std::make_shared<PopulationNode>(this->data_.get_population_label(0)));
    ancestor->add_child(std::make_shared<PopulationNode>(this->data_.get_population_label(1)));
    for (unsigned int pop_idx = 2;
            pop_idx < this->data_.get_number_of_populations();
            ++pop_idx) {
        std::shared_ptr<PopulationNode> next_ancestor = std::make_shared<PopulationNode>(0.0);
        next_ancestor->add_child(ancestor);
        std::shared_ptr<PopulationNode> tip = std::make_shared<PopulationNode>(
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
        PopulationNode& node) {
    unsigned int pop_idx = this->data_.get_population_index(node.get_label());
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
            node.copy_bottom_pattern_probs(m);
            return;
        }
        BiallelicPatternProbabilityMatrix m(allele_count, n_reds);
        node.copy_bottom_pattern_probs(m);
        return;
    }
    BiallelicPatternProbabilityMatrix m(allele_count, red_allele_count);
    node.copy_bottom_pattern_probs(m);
    return;
}

void PopulationTree::compute_top_of_branch_partials(
        PopulationNode& node) {
    if (node.get_allele_count() == 0) {
        node.copy_top_pattern_probs(node.get_bottom_pattern_probs());
        return;
    }

    BiallelicPatternProbabilityMatrix m = matrix_exponentiator.expQTtx(
            node.get_allele_count(),
            this->u_->get_value(),
            this->get_v(),
            node.get_population_size() * this->get_rate_multiplier(),
            node.get_length() * this->get_rate_multiplier(),
            node.get_bottom_pattern_probs());
    node.copy_top_pattern_probs(m);
}

void PopulationTree::compute_internal_partials(
        PopulationNode& node) {
    if (node.get_number_of_children() == 1) {
        node.copy_bottom_pattern_probs(node.get_child(0)->get_top_pattern_probs());
        return;
    }
    if (node.get_child(0)->get_allele_count() == 0) {
        node.copy_bottom_pattern_probs(node.get_child(1)->get_top_pattern_probs());
        return;
    }
    if (node.get_child(1)->get_allele_count() == 0) {
        node.copy_bottom_pattern_probs(node.get_child(0)->get_top_pattern_probs());
        return;
    }
    unsigned int allele_count_child1 = node.get_child(0)->get_allele_count();
    unsigned int allele_count_child2 = node.get_child(1)->get_allele_count();

    std::vector<double> pattern_probs_child1 = node.get_child(0)->get_top_pattern_probs().get_pattern_prob_matrix();
    std::vector<double> pattern_probs_child2 = node.get_child(1)->get_top_pattern_probs().get_pattern_prob_matrix();

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
    node.copy_bottom_pattern_probs(m);
}

void PopulationTree::compute_pattern_partials(
        int pattern_index,
        PopulationNode& node) {
    if (node.is_leaf()) {
        this->compute_leaf_partials(pattern_index, node);
    }
    else if (node.get_number_of_children() == 1) {
        compute_pattern_partials(pattern_index, *node.get_child(0));
        compute_top_of_branch_partials(*node.get_child(0));
        compute_internal_partials(node);
    }
    else if (node.get_number_of_children() == 2) {
        compute_pattern_partials(pattern_index, *node.get_child(0));
        compute_pattern_partials(pattern_index, *node.get_child(1));
        compute_top_of_branch_partials(*node.get_child(0));
        compute_top_of_branch_partials(*node.get_child(1));
        compute_internal_partials(node);
    }
    else {
        std::ostringstream message;
        message << "PopulationTree::compute_pattern_probability(); "
                << "unexpected number of children: "
                << node.get_number_of_children();
        throw EcoevolityError(message.str());
    }
}

std::vector< std::vector<double> > PopulationTree::compute_root_probabilities() {
    unsigned int N = this->root_->get_allele_count();
    std::vector< std::vector<double> > x (N + 1); 
    QMatrix q = QMatrix(
            N,
            this->u_->get_value(),
            this->get_v(),
            this->root_->get_population_size() * this->get_rate_multiplier());
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
    this->compute_pattern_partials(pattern_index, *this->root_);
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
    if (! this->u_v_rates_are_constrained()) {
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
    if ((red_allele_count == 0) || (red_allele_count == allele_count)) {
        return 0.0;
    }
    double f = 0.0;
    for (unsigned int i = red_allele_count + 1; i <= allele_count; ++i) {
        f += std::log(i) - std::log(allele_count - i + 1);
    }
    return f;
}

bool PopulationTree::constant_site_counts_were_provided() {
    if ((this->provided_number_of_constant_green_sites_ > -1) && (this->provided_number_of_constant_red_sites_ > -1)) {
        return true;
    }
    return false;
}

double PopulationTree::compute_log_likelihood() {
    if (this->ignore_data_) {
        this->log_likelihood_.set_value(0.0);
        return 0.0;
    }
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

    if (this->constant_sites_removed_) {
        if (this->constant_site_counts_were_provided()) {
            double constant_log_likelihood =
                    ((double)this->provided_number_of_constant_green_sites_ * std::log(this->all_green_pattern_likelihood_.get_value())) +
                    ((double)this->provided_number_of_constant_red_sites_ * std::log(this->all_red_pattern_likelihood_.get_value()));
            log_likelihood += constant_log_likelihood;
        }
        //////////////////////////////////////////////////////////////////////
        // No reason to use removed site counts. Simply leave constant sites in
        // and calc likelihood without correction. This is better, because it
        // doesn't treat all constant site patterns equally (i.e., it accounts
        // for constant patterns with missing data).
        // else if (this->use_removed_constant_site_counts_){
        //     double constant_log_likelihood =
        //             ((double)this->data_.get_number_of_constant_green_sites_removed() *
        //             std::log(this->all_green_pattern_likelihood_.get_value())) +
        //             ((double)this->data_.get_number_of_constant_red_sites_removed() *
        //             std::log(this->all_red_pattern_likelihood_.get_value()));
        //     log_likelihood += constant_log_likelihood;
        // }
        //////////////////////////////////////////////////////////////////////
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
    if (! this->u_v_rates_are_constrained()) {
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
double PopulationTree::get_root_height() const {
    return this->root_->get_height();
}
void PopulationTree::store_root_height() {
    this->root_->store_height();
}
void PopulationTree::restore_root_height() {
    this->root_->restore_height();
}
void PopulationTree::set_root_height_parameter(std::shared_ptr<PositiveRealParameter> h) {
    ECOEVOLITY_ASSERT(h->prior == this->node_height_prior_);
    this->root_->set_height_parameter(h);
}
std::shared_ptr<PositiveRealParameter> PopulationTree::get_root_height_parameter() const {
    return this->root_->get_height_parameter();
}

void PopulationTree::set_u(double u) {
    if (this->u_v_rates_are_fixed()) {
        return;
    }
    ECOEVOLITY_ASSERT(u >= 0.5);
    this->u_->set_value(u);
    ECOEVOLITY_DEBUG(
        double v = this->get_v();
        ECOEVOLITY_ASSERT_APPROX_EQUAL(2*u*v/(u+v), 1.0);
    )
    this->make_dirty();
}
void PopulationTree::update_u(double u) {
    if (this->u_v_rates_are_fixed()) {
        return;
    }
    ECOEVOLITY_ASSERT(u >= 0.5);
    this->u_->update_value(u);
    ECOEVOLITY_DEBUG(
        double v = this->get_v();
        ECOEVOLITY_ASSERT_APPROX_EQUAL(2*u*v/(u+v), 1.0);
    )
    this->make_dirty();
}
double PopulationTree::get_u() const {
    return this->u_->get_value();
}
double PopulationTree::get_v() const {
    double u = this->get_u();
    return u / ((2.0 * u) - 1.0);
}
void PopulationTree::store_u() {
    this->u_->store();
}
void PopulationTree::restore_u() {
    this->u_->restore();
    this->make_dirty();
}

void PopulationTree::set_rate_multiplier(double m) {
    if (this->rate_multiplier_is_fixed()) {
        return;
    }
    this->rate_multiplier_->set_value(m);
    this->make_dirty();
}
void PopulationTree::update_rate_multiplier(double m) {
    if (this->rate_multiplier_is_fixed()) {
        return;
    }
    this->rate_multiplier_->update_value(m);
    this->make_dirty();
}
double PopulationTree::get_rate_multiplier() const {
    return this->rate_multiplier_->get_value();
}
void PopulationTree::store_rate_multiplier() {
    this->rate_multiplier_->store();
}
void PopulationTree::restore_rate_multiplier() {
    this->rate_multiplier_->restore();
    this->make_dirty();
}

std::shared_ptr<PositiveRealParameter> PopulationTree::get_u_parameter() const {
    return this->u_;
}

void PopulationTree::set_rate_multiplier_parameter(std::shared_ptr<PositiveRealParameter> h) {
    this->rate_multiplier_ = h;
    this->make_dirty();
}
std::shared_ptr<PositiveRealParameter> PopulationTree::get_rate_multiplier_parameter() const {
    return this->rate_multiplier_;
}

void PopulationTree::set_root_population_size(double size) {
    if (this->population_sizes_are_fixed()) {
        return;
    }
    this->root_->set_population_size(size);
}
void PopulationTree::set_population_size(double size) {
    if (this->population_sizes_are_fixed()) {
        return;
    }
    this->root_->set_all_population_sizes(size);
}

double PopulationTree::get_root_population_size() const {
    return this->root_->get_population_size();
}
std::shared_ptr<PositiveRealParameter> PopulationTree::get_root_population_size_parameter() const {
    return this->root_->get_population_size_parameter();
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
    this->store_rate_multiplier();
    this->store_all_population_sizes();
    this->store_all_heights();
}
void PopulationTree::store_all_population_sizes() {
    this->root_->store_all_population_sizes();
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
    this->restore_rate_multiplier();
    this->restore_all_population_sizes();
    this->restore_all_heights();
}
void PopulationTree::restore_all_population_sizes() {
    this->root_->restore_all_population_sizes();
}
void PopulationTree::restore_all_heights() {
    this->root_->restore_all_heights();
}

void PopulationTree::set_node_height_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
    this->node_height_prior_ = prior;
    this->root_->set_all_node_height_priors(prior);
}

void PopulationTree::set_population_size_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
    this->population_size_prior_ = prior;
    this->root_->set_all_population_size_priors(prior);
}

void PopulationTree::set_u_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
    this->u_->set_prior(prior);
    this->make_dirty();
}
void PopulationTree::set_rate_multiplier_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
    this->rate_multiplier_->set_prior(prior);
    this->make_dirty();
}

double PopulationTree::compute_log_prior_density() {
    double d = 0.0;
    d += this->compute_log_prior_density_of_u_v_rates();
    d += this->compute_log_prior_density_of_rate_multiplier();
    d += this->compute_log_prior_density_of_node_heights();
    d += this->compute_log_prior_density_of_population_sizes();
    this->log_prior_density_.set_value(d);
    return d;
}
double PopulationTree::compute_log_prior_density_of_u_v_rates() const {
    return this->u_->relative_prior_ln_pdf();
}
double PopulationTree::compute_log_prior_density_of_rate_multiplier() const {
    return this->rate_multiplier_->relative_prior_ln_pdf();
}
double PopulationTree::compute_log_prior_density_of_node_heights() const {
    return this->root_->calculate_ln_relative_node_height_prior_density();
}
double PopulationTree::compute_log_prior_density_of_population_sizes() const {
    return this->root_->calculate_ln_relative_population_size_prior_density();
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

void PopulationTree::provide_number_of_constant_sites(
                unsigned int number_all_red,
                unsigned int number_all_green) {
    if (! this->constant_sites_removed_) {
        throw EcoevolityError(
                "Trying to provide number of constant sites, but they haven't been removed");
    }
    this->provided_number_of_constant_red_sites_ = number_all_red;
    this->provided_number_of_constant_green_sites_ = number_all_green;
}


ComparisonPopulationTree::ComparisonPopulationTree(
        std::string path, 
        char population_name_delimiter,
        bool population_name_is_prefix,
        bool genotypes_are_diploid,
        bool markers_are_dominant,
        bool constant_sites_removed,
        bool validate,
        bool strict_on_constant_sites,
        bool strict_on_missing_sites) {
    this->init(path,
               population_name_delimiter,
               population_name_is_prefix,
               genotypes_are_diploid,
               markers_are_dominant,
               constant_sites_removed,
               validate,
               strict_on_constant_sites,
               strict_on_missing_sites);
    if (this->data_.get_number_of_populations() > 2) {
        throw EcoevolityError("ComparisonPopulationTree(); does not support more than 2 populations");
    }
    this->root_->set_label("root-" + this->root_->get_child(0)->get_label());
}
ComparisonPopulationTree::ComparisonPopulationTree(
        const ComparisonSettings& settings,
        RandomNumberGenerator& rng,
        bool strict_on_constant_sites,
        bool strict_on_missing_sites) {
    this->init(settings.get_path(),
               settings.get_population_name_delimiter(),
               settings.population_name_is_prefix(),
               settings.genotypes_are_diploid(),
               settings.markers_are_dominant(),
               settings.constant_sites_removed(),
               true, // validate
               strict_on_constant_sites,
               strict_on_missing_sites);
    if (settings.constrain_u_v_rates()) {
        this->constrain_u_v_rates();
        this->fold_patterns();
    }
    this->set_population_size_prior(
            settings.get_population_size_settings().get_prior_settings().get_instance());
    if (settings.constrain_population_sizes()) {
        this->constrain_population_sizes();
    }
    PositiveRealParameter p = PositiveRealParameter(
            settings.get_population_size_settings(),
            rng);
    this->set_population_size(p.get_value());
    if (settings.get_population_size_settings().is_fixed()) {
        this->fix_population_sizes();
    }
    
    this->set_u_prior(settings.get_u_settings().get_prior_settings().get_instance());
    if (settings.constrain_u_v_rates()) {
        this->constrain_u_v_rates();
    }
    else {
        PositiveRealParameter u = PositiveRealParameter(
                settings.get_u_settings(),
                rng);
        this->set_u(u.get_value());
        if (u.is_fixed()) {
            this->fix_u_v_rates();
        }
    }
    this->set_rate_multiplier_parameter(
            std::make_shared<PositiveRealParameter>(
                    settings.get_rate_multiplier_settings(),
                    rng));
    if (
        (this->data_.get_number_of_populations() == 1) &&
        (
            (this->population_sizes_are_fixed()) ||
            (this->population_sizes_are_constrained())
        )
    ) {
        std::ostringstream message;
        message << "\n#######################################################################\n"
                <<   "###############################  ERROR  ###############################\n"
                << "The alignment in:\n    \'"
                << this->data_.get_path() << "\'\n"
                << "contains only a single population, but you have fixed and/or "
                << "constrained the population sizes for this comparison. The timing of "
                << "population expansion/contraction cannot be estimated if the ancestral "
                << "and descendant population sizes for this comparison are either "
                << "fixed or constrained to be equal. Please update your configuration "
                << "file to estimate the unconstrained population sizes for this "
                << "comparison and re-run the analysis.\n"
                << "#######################################################################\n";
        throw EcoevolityComparisonSettingError(message.str(), this->data_.get_path());
    }
}

void ComparisonPopulationTree::set_child_population_size(
        unsigned int child_index,
        double size) {
    if (this->population_sizes_are_fixed()) {
        return;
    }
    this->root_->get_child(child_index)->set_population_size(size);
}
void ComparisonPopulationTree::update_child_population_size(
        unsigned int child_index,
        double size) {
    if (this->population_sizes_are_fixed()) {
        return;
    }
    this->root_->get_child(child_index)->update_population_size(size);
}
double ComparisonPopulationTree::get_child_population_size(
        unsigned int child_index) const {
    return this->root_->get_child(child_index)->get_population_size();
}
void ComparisonPopulationTree::store_child_population_size(
        unsigned int child_index) {
    this->root_->get_child(child_index)->store_population_size();
}
void ComparisonPopulationTree::restore_child_population_size(
        unsigned int child_index) {
    this->root_->get_child(child_index)->restore_population_size();
}
std::shared_ptr<PositiveRealParameter> ComparisonPopulationTree::get_child_population_size_parameter(
        unsigned int child_index) const {
    return this->root_->get_child(child_index)->get_population_size_parameter();
}

// Node height sharing needs to be dealt with in next level up in
// class hierarchy (ComparisonPopulationTreeCollection)
double ComparisonPopulationTree::compute_log_prior_density() {
    double d = 0.0;
    d += this->compute_log_prior_density_of_u_v_rates();
    d += this->compute_log_prior_density_of_rate_multiplier();
    // d += this->compute_log_prior_density_of_node_heights();
    d += this->compute_log_prior_density_of_population_sizes();
    this->log_prior_density_.set_value(d);
    return d;
}

// Node height (re)storing is managed by ComparisonPopulationTree.
void ComparisonPopulationTree::store_parameters() {
    this->store_u();
    this->store_rate_multiplier();
    this->store_all_population_sizes();
    // this->store_all_heights();
}
void ComparisonPopulationTree::restore_parameters() {
    this->restore_u();
    this->restore_rate_multiplier();
    this->restore_all_population_sizes();
    // this->restore_all_heights();
}

void ComparisonPopulationTree::write_state_log_header(
        std::ostream& out,
        bool include_event_index,
        const std::string& delimiter) const {
    std::string suffix = "_" + this->root_->get_child(0)->get_label();
    if (include_event_index) {
        out << "root_height_index" << suffix << delimiter;
    }
    out << "ln_likelihood" << suffix << delimiter
        << "ln_prior" << suffix << delimiter
        << "root_height" << suffix << delimiter
        << "rate_multiplier" << suffix << delimiter
        << "u" << suffix << delimiter
        << "v" << suffix << delimiter
        << "pop_size" << suffix << delimiter;
    if (this->root_->get_number_of_children() > 1) {
        out << "pop_size" << "_" << this->root_->get_child(1)->get_label() << delimiter;
    }
    out << "pop_size_root" << suffix;
}

void ComparisonPopulationTree::log_state(
        std::ostream& out,
        const std::string& delimiter) const {
    out << this->log_likelihood_.get_value() << delimiter
        << this->log_prior_density_.get_value() << delimiter
        << this->get_height() << delimiter
        << this->get_rate_multiplier() << delimiter
        << this->get_u() << delimiter
        << this->get_v() << delimiter
        << this->get_child_population_size(0) << delimiter;
    if (this->root_->get_number_of_children() > 1) {
        out << this->get_child_population_size(1) << delimiter;
    }
    out << this->get_root_population_size();
}
void ComparisonPopulationTree::log_state(
        std::ostream& out,
        unsigned int event_index,
        const std::string& delimiter) const {
    out << event_index << delimiter;
    this->log_state(out, delimiter);
}

std::string ComparisonPopulationTree::get_state_header_string(
        const std::string& delimiter) const {
    std::ostringstream ss;
    this->write_state_log_header(ss, false, delimiter);
    return ss.str();
}

std::string ComparisonPopulationTree::get_state_string(
        const std::string& delimiter,
        unsigned int precision) const {
    std::ostringstream ss;
    ss.precision(precision);
    this->log_state(ss, delimiter);
    return ss.str();
}

// TODO: This is a hack. The more general solution would be a recursive method
// of PopulationTree
std::shared_ptr<GeneTreeSimNode> ComparisonPopulationTree::simulate_gene_tree(
        const unsigned int pattern_index,
        RandomNumberGenerator& rng) const {

    std::vector< std::shared_ptr<GeneTreeSimNode> > left_lineages;
    std::vector< std::shared_ptr<GeneTreeSimNode> > right_lineages;
    std::vector< std::shared_ptr<GeneTreeSimNode> > root_lineages;
    std::vector<std::string> tip_labels;
    unsigned int allele_count;
    tip_labels = this->data_.get_sequence_labels(
            this->data_.get_population_index(
                    this->root_->get_child(0)->get_label()));
    allele_count = this->data_.get_allele_count(
            pattern_index,
            this->data_.get_population_index(
                    this->root_->get_child(0)->get_label()));
    if (this->data_.markers_are_dominant()) {
        allele_count *= 2;
    }
    left_lineages.reserve(allele_count);
    for (unsigned int tip_idx = 0; tip_idx < allele_count; ++tip_idx) {
        std::shared_ptr<GeneTreeSimNode> tip = std::make_shared<GeneTreeSimNode>(
                    tip_labels.at(tip_idx),
                    0.0);
            tip->fix_node_height();
            left_lineages.push_back(tip);
    }

    double top_of_branch_height = this->get_height() * this->get_rate_multiplier();
    double current_height = 0.0;
    double last_left_coal_height = this->coalesce_in_branch(
            left_lineages,
            this->get_child_population_size(0) * this->get_rate_multiplier(),
            rng,
            current_height,
            top_of_branch_height
            );

    if (this->root_->get_number_of_children() > 1) {
        tip_labels = this->data_.get_sequence_labels(
                this->data_.get_population_index(
                        this->root_->get_child(1)->get_label()));
        allele_count = this->data_.get_allele_count(
                pattern_index,
                this->data_.get_population_index(
                        this->root_->get_child(1)->get_label()));
        if (this->data_.markers_are_dominant()) {
            allele_count *= 2;
        }
        right_lineages.reserve(allele_count);
        for (unsigned int tip_idx = 0; tip_idx < allele_count; ++tip_idx) {
            std::shared_ptr<GeneTreeSimNode> tip = std::make_shared<GeneTreeSimNode>(
                        tip_labels.at(tip_idx),
                        0.0);
                tip->fix_node_height();
                right_lineages.push_back(tip);
        }

        double last_right_coal_height = this->coalesce_in_branch(
                right_lineages,
                this->get_child_population_size(1) * this->get_rate_multiplier(),
                rng,
                current_height,
                top_of_branch_height
                );
    }

    for (unsigned int i = 0; i < left_lineages.size(); ++i) {
        root_lineages.push_back(left_lineages.at(i));
    }
    left_lineages.clear();
    for (unsigned int i = 0; i < right_lineages.size(); ++i) {
        root_lineages.push_back(right_lineages.at(i));
    }
    right_lineages.clear();
    ECOEVOLITY_ASSERT(root_lineages.size() > 0);
    if (root_lineages.size() == 1) {
        return root_lineages.at(0);
    }
    double last_root_coal_height = this->coalesce_in_branch(
            root_lineages,
            this->get_root_population_size() * this->get_rate_multiplier(),
            rng,
            top_of_branch_height,
            std::numeric_limits<double>::infinity()
            );
    ECOEVOLITY_ASSERT(root_lineages.size() == 1);
    return root_lineages.at(0);
}

double ComparisonPopulationTree::coalesce_in_branch(
        std::vector< std::shared_ptr<GeneTreeSimNode> >& lineages,
        double population_size,
        RandomNumberGenerator& rng,
        double bottom_of_branch_height,
        double top_of_branch_height
        ) {
    ECOEVOLITY_ASSERT(lineages.size() > 0);
    ECOEVOLITY_ASSERT(bottom_of_branch_height < top_of_branch_height);
    double current_height = bottom_of_branch_height;
    unsigned int k = lineages.size();
    while (true) {
        if (k == 1) {
            ECOEVOLITY_ASSERT(lineages.size() == k);
            break;
        }
        double scale = population_size / (((double)k) * (k - 1.0));
        double wait = rng.gamma(1.0, scale);
        
        if ((current_height + wait) >= top_of_branch_height) {
            break;
        }
        current_height += wait;
        std::shared_ptr<GeneTreeSimNode> mrca = std::make_shared<GeneTreeSimNode>(
                current_height);
        for (int i = 0; i < 2; ++i) {
            int idx = rng.uniform_int(0, lineages.size() - 1);
            mrca->add_child(lineages.at(idx));
            lineages.erase(lineages.begin() + idx);
        }
        lineages.push_back(mrca);
        --k;
        ECOEVOLITY_ASSERT(lineages.size() == k);
    }
    return current_height;
}

BiallelicData ComparisonPopulationTree::simulate_biallelic_data_set(
        RandomNumberGenerator& rng,
        bool validate) const {
    BiallelicData sim_data = this->data_.get_empty_copy();
    const bool filtering_constant_sites = this->constant_sites_removed_;
    std::unordered_map<std::string, unsigned int> seq_label_to_pop_index_map;
    for (unsigned int pop_idx = 0;
            pop_idx < this->data_.get_number_of_populations();
            ++pop_idx) {
        for (auto seq_label: this->data_.get_sequence_labels(pop_idx)) {
            ECOEVOLITY_ASSERT(
                this->data_.get_population_index_from_seq_label(seq_label) ==
                pop_idx);
            seq_label_to_pop_index_map[seq_label] = pop_idx;
        }
    }
    // Looping over patterns to make sure simulated dataset has exact same
    // sample configuration (i.e., the same pattern of missing data) as the
    // member dataset.
    for (unsigned int pattern_idx = 0;
            pattern_idx < this->data_.get_number_of_patterns();
            ++pattern_idx) {
        for (unsigned int i = 0;
                i < this->data_.get_pattern_weight(pattern_idx);
                ++i) {
            bool site_added = false;
            while (! site_added) {
                auto pattern_tree = this->simulate_biallelic_site(
                        pattern_idx,
                        seq_label_to_pop_index_map,
                        rng);
                auto pattern = pattern_tree.first;
                std::vector<unsigned int> red_allele_counts = pattern.first;
                std::vector<unsigned int> allele_counts = pattern.second;
                std::shared_ptr<GeneTreeSimNode> gtree = pattern_tree.second;
                site_added = sim_data.add_site(red_allele_counts,
                        allele_counts,
                        filtering_constant_sites);
            }
        }
    }
    sim_data.update_pattern_booleans();
    if (validate) {
        sim_data.validate();
    }
    // What about gene trees? Write to stream/path (don't want to store/return
    // giant vector)
    return sim_data;
}

std::pair<
        std::pair<std::vector<unsigned int>, std::vector<unsigned int> >,
        std::shared_ptr<GeneTreeSimNode> >
ComparisonPopulationTree::simulate_biallelic_site(
        const unsigned int pattern_idx,
        std::unordered_map<std::string, unsigned int> seq_label_to_pop_index_map,
        RandomNumberGenerator& rng) const {
    double freq_0 = this->get_u() / (this->get_u() + this->get_v());

    std::shared_ptr<GeneTreeSimNode> gene_tree = this->simulate_gene_tree(pattern_idx, rng);
    gene_tree->compute_binary_transition_probabilities(this->get_u(), this->get_v());
    gene_tree->simulate_binary_character(freq_0, rng);

    const std::vector<unsigned int>& expected_allele_counts = this->data_.get_allele_counts(pattern_idx);
    std::vector<unsigned int> allele_counts(expected_allele_counts.size(), 0);
    std::vector<unsigned int> red_allele_counts(expected_allele_counts.size(), 0);
    if (this->data_.markers_are_dominant()) {
        std::vector<int> last_allele(expected_allele_counts.size(), -1);
        gene_tree->get_allele_counts(
                seq_label_to_pop_index_map,
                allele_counts,
                red_allele_counts,
                last_allele);
    }
    else {
        gene_tree->get_allele_counts(
                seq_label_to_pop_index_map,
                allele_counts,
                red_allele_counts);
    }
    ECOEVOLITY_ASSERT(allele_counts == expected_allele_counts);
    std::pair<std::vector<unsigned int>, std::vector<unsigned int> > pattern = 
            std::make_pair(red_allele_counts, allele_counts);
    return std::make_pair(pattern, gene_tree);
}
