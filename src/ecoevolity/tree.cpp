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


BasePopulationTree::BasePopulationTree(
        std::string path, 
        char population_name_delimiter,
        bool population_name_is_prefix,
        bool genotypes_are_diploid,
        bool markers_are_dominant,
        bool constant_sites_removed,
        bool validate,
        bool strict_on_constant_sites,
        bool strict_on_missing_sites,
        bool strict_on_triallelic_sites,
        double ploidy,
        bool store_seq_loci_info
        ) : BaseTree<PopulationNode>() {
    this->init(path,
               population_name_delimiter,
               population_name_is_prefix,
               genotypes_are_diploid,
               markers_are_dominant,
               constant_sites_removed,
               validate,
               strict_on_constant_sites,
               strict_on_missing_sites,
               strict_on_triallelic_sites,
               ploidy,
               store_seq_loci_info);
}

BasePopulationTree::BasePopulationTree(
        std::shared_ptr<PopulationNode> root,
        unsigned int number_of_loci,
        unsigned int length_of_loci,
        bool validate_data) : BaseTree<PopulationNode>() {
    const std::vector< std::shared_ptr<PopulationNode> >& leaves = root->get_leaves();
    std::vector<std::string> pop_labels;
    std::vector<unsigned int> haploid_sample_sizes;
    std::vector<unsigned int> leaf_pop_indices;
    std::vector<unsigned int> internal_pop_indices;
    std::vector<unsigned int> expected_leaf_pop_indices;
    std::vector<unsigned int> expected_internal_pop_indices;
    pop_labels.reserve(leaves.size());
    haploid_sample_sizes.reserve(leaves.size());
    leaf_pop_indices.reserve(leaves.size());
    expected_leaf_pop_indices.reserve(leaves.size());
    for (unsigned int i = 0; i < leaves.size(); ++i) {
        pop_labels.push_back(leaves.at(i)->get_label());
        haploid_sample_sizes.push_back(leaves.at(i)->get_allele_count());
        expected_leaf_pop_indices.push_back(i);
    }
    root->get_node_indices(internal_pop_indices, leaf_pop_indices);
    for (unsigned int i = leaves.size();
            i < (leaves.size() + internal_pop_indices.size());
            ++i) {
        expected_internal_pop_indices.push_back(i);
    }
    ECOEVOLITY_ASSERT(std::is_permutation(
                leaf_pop_indices.begin(), leaf_pop_indices.end(),
                expected_leaf_pop_indices.begin()));
    ECOEVOLITY_ASSERT(std::is_permutation(
                internal_pop_indices.begin(), internal_pop_indices.end(),
                expected_internal_pop_indices.begin()));
    BiallelicData bd(pop_labels,
            haploid_sample_sizes,
            number_of_loci,
            length_of_loci,
            validate_data);
    this->data_ = bd;
    this->set_root(root);
    this->constant_sites_removed_ = false;
    this->root_->resize_all();
    std::shared_ptr<ContinuousProbabilityDistribution> root_height_prior = std::make_shared<ExponentialDistribution>(100.0);
    this->root_->set_all_node_height_priors(root_height_prior);
    this->root_->set_all_population_size_priors(this->population_size_prior_);
    this->update_unique_allele_counts();
    this->is_dirty_ = true;
    this->number_of_likelihood_calculations_ = 0;
    this->likelihood_correction_was_calculated_ = false;
    if (validate_data) {
        this->data_.validate();
    }
    this->update_node_heights();
}

void BasePopulationTree::init(
        std::string path, 
        char population_name_delimiter,
        bool population_name_is_prefix,
        bool genotypes_are_diploid,
        bool markers_are_dominant,
        bool constant_sites_removed,
        bool validate,
        bool strict_on_constant_sites,
        bool strict_on_missing_sites,
        bool strict_on_triallelic_sites,
        double ploidy,
        bool store_seq_loci_info) {
    this->set_ploidy(ploidy);
    try {
        // First, try to parse data as YAML formatted
        if (store_seq_loci_info) {
            // If we are expecting charset info, we need to abort YAML
            // processing and try nexus parsing.
            throw EcoevolityYamlDataError(
                    "Charsets option specified, abort YAML processing");
        }
        this->data_.init_from_yaml_path(
                path,
                validate);
    }
    catch (...) {
        // If YAML parsing failed, try parsing nexus data
        this->data_ = BiallelicData();
        if (genotypes_are_diploid && (ploidy != 2.0)) {
            throw EcoevolityBiallelicDataError(
                    "Genotypes cannot be diploid if ploidy is not 2",
                    path);
        }
        this->data_.init(
                path,
                population_name_delimiter,
                population_name_is_prefix,
                genotypes_are_diploid,
                markers_are_dominant,
                validate,
                store_seq_loci_info);
    }
    this->constant_sites_removed_ = constant_sites_removed;

    this->process_and_vet_initialized_data(
            strict_on_constant_sites,
            strict_on_missing_sites,
            strict_on_triallelic_sites);

    this->init_tree();

    this->root_->resize_all();

    // Store unique allele counts and weights for calculating the likelihood
    // correction term for constant sites.
    // At this point, no data manipulation should happen that would require
    // these to be updated (only pattern folding, which will not change the
    // unique allele counts or weights).
    this->update_unique_allele_counts();
    this->update_node_heights();
}

void BasePopulationTree::process_and_vet_initialized_data(
        bool strict_on_constant_sites,
        bool strict_on_missing_sites,
        bool strict_on_triallelic_sites) {
    if (this->data_.get_number_of_populations() < 1) {
        throw EcoevolityError("BasePopulationTree(); no populations were found");
    }
    unsigned int number_of_missing_patterns_removed = this->data_.remove_missing_population_patterns();
    if (this->data_.has_recoded_triallelic_sites()) {
        if (strict_on_triallelic_sites) {
            std::ostringstream message;
            message << "\n#######################################################################\n"
                    <<   "###############################  ERROR  ###############################\n"
                    << this->data_.get_number_of_triallelic_sites_recoded()
                    << " sites from the alignment in:\n    \'"
                    << this->data_.get_path() << "\'\n"
                    << "have more than two character states.\n"
                    << "#######################################################################\n";
            throw EcoevolityTriallelicDataError(message.str(), this->data_.get_path());
        }
        else {
            std::ostringstream message;
            message << "\n#######################################################################\n"
                    <<   "##############################  WARNING  ##############################\n"
                    << this->data_.get_number_of_triallelic_sites_recoded()
                    << " sites had more than two nucleotide states from the alignment in:\n    \'"
                    << this->data_.get_path() << "\'.\n"
                    << "These sites have been recoded as biallelic, by treating the first\n"
                    << "nucleotide as 0 and all others as 1. If you would prefer to ignore\n"
                    << "these sites, please remove all sites with more than two nucleotide\n"
                    << "states from your DNA alignments and re-run the analysis.\n"
                    << "#######################################################################\n";
            std::cerr << message.str() << std::endl;
        }
    }
    if (number_of_missing_patterns_removed > 0) {
        if (strict_on_missing_sites) {
            std::ostringstream message;
            message << "\n#######################################################################\n"
                    <<   "###############################  ERROR  ###############################\n"
                    << this->data_.get_number_of_missing_sites_removed()
                    << " sites from the alignment in:\n    \'"
                    << this->data_.get_path() << "\'\n"
                    << "have no data for at least one population.\n"
                    << "#######################################################################\n";
            throw EcoevolityMissingDataError(message.str(), this->data_.get_path());
        }
        else {
            std::ostringstream message;
            message << "\n#######################################################################\n"
                    <<   "##############################  WARNING  ##############################\n"
                    << this->data_.get_number_of_missing_sites_removed()
                    << " sites will be ignored from the alignment in:\n    \'"
                    << this->data_.get_path() << "\'\n"
                    << "due to at least one population with no data.\n"
                    << "#######################################################################\n";
            std::cerr << message.str() << std::endl;
        }
    }
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
                        << "    \'" << this->data_.get_path() << "\'\n"
                        << "but you indicated that such sites were already removed with option:\n"
                        << "    constant_sites_removed = true\n"
                        << "If you intended to remove them, please do so and re-run the analysis.\n"
                        << "If you intended for constant sites to be used in the likelihood\n"
                        << "calculations, you should set \'constant_sites_removed\' to false for\n"
                        << "this alignment and re-run the analysis.\n"
                        << "#######################################################################\n";
                throw EcoevolityConstantSitesError(message.str(), this->data_.get_path());
            }
            else {
                std::ostringstream message;
                message << "\n#######################################################################\n"
                        <<   "##############################  WARNING  ##############################\n"
                        << this->data_.get_number_of_constant_sites_removed()
                        << " constant sites were found in the alignment in:\n"
                        << "    \'" << this->data_.get_path() << "\'\n"
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
}

void BasePopulationTree::init_tree() {
    if (this->data_.get_number_of_populations() < 3) {
        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(this->data_.get_number_of_populations(), 0.0);
        for (unsigned int pop_idx = 0;
                pop_idx < this->data_.get_number_of_populations();
                ++pop_idx) {
            std::shared_ptr<PopulationNode> tip = std::make_shared<PopulationNode>(
                    pop_idx,
                    this->data_.get_population_label(pop_idx),
                    0.0,
                    this->data_.get_max_allele_count(pop_idx));
            tip->fix_node_height();
            root->add_child(tip);
        }
        this->set_root(root);
        return;
    }
    unsigned int next_index = this->data_.get_number_of_populations();
    std::shared_ptr<PopulationNode> ancestor = std::make_shared<PopulationNode>(next_index, 0.0);
    ++next_index;
    ancestor->add_child(std::make_shared<PopulationNode>(
                0,
                this->data_.get_population_label(0),
                0.0,
                this->data_.get_max_allele_count(0)
                ));
    ancestor->add_child(std::make_shared<PopulationNode>(
                1,
                this->data_.get_population_label(1),
                0.0,
                this->data_.get_max_allele_count(1)
                ));
    for (unsigned int pop_idx = 2;
            pop_idx < this->data_.get_number_of_populations();
            ++pop_idx) {
        std::shared_ptr<PopulationNode> next_ancestor = std::make_shared<PopulationNode>(next_index, 0.0);
        ++next_index;
        next_ancestor->add_child(ancestor);
        std::shared_ptr<PopulationNode> tip = std::make_shared<PopulationNode>(
                pop_idx,
                this->data_.get_population_label(pop_idx),
                0.0,
                this->data_.get_max_allele_count(pop_idx));
        tip->fix_node_height();
        next_ancestor->add_child(tip);
        ancestor = next_ancestor;
    }
    this->set_root(ancestor);
}

void BasePopulationTree::set_data(const BiallelicData & data, bool constant_sites_removed) {
    const std::vector< std::shared_ptr<PopulationNode> >& leaves = this->root_->get_leaves();
    ECOEVOLITY_ASSERT(this->data_.get_number_of_populations() == leaves.size());
    for (unsigned int i = 0; i < leaves.size(); ++i) {
        if (leaves.at(i)->get_allele_count() != data.get_max_allele_count(leaves.at(i)->get_index())) {
            leaves.at(i)->resize(data.get_max_allele_count(leaves.at(i)->get_index()));
        }
        if (leaves.at(i)->get_label() != data.get_population_label(leaves.at(i)->get_index())) {
            leaves.at(i)->set_label(data.get_population_label(leaves.at(i)->get_index()));
        }
    }
    this->root_->resize_all();
    this->data_ = data;
    this->constant_sites_removed_ = constant_sites_removed;
    this->update_unique_allele_counts();
    this->is_dirty_ = true;
    this->number_of_likelihood_calculations_ = 0;
    this->likelihood_correction_was_calculated_ = false;
}

void BasePopulationTree::update_unique_allele_counts() {
    this->unique_allele_counts_.clear();
    this->unique_allele_count_weights_.clear();
    std::map<std::vector<unsigned int>, unsigned int> unique_allele_count_map = this->data_.get_unique_allele_counts();
    for (auto const & kv: unique_allele_count_map) {
        this->unique_allele_counts_.push_back(kv.first);
        this->unique_allele_count_weights_.push_back(kv.second);
    }
}

void BasePopulationTree::calculate_likelihood_correction() {
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

double BasePopulationTree::get_likelihood_correction(bool force) {
    if ((! this->likelihood_correction_was_calculated_) || (force)) {
        this->calculate_likelihood_correction();
    }
    return this->log_likelihood_correction_.get_value();
}

double BasePopulationTree::calculate_log_binomial(
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

// bool BasePopulationTree::constant_site_counts_were_provided() {
//     if ((this->provided_number_of_constant_green_sites_ > -1) && (this->provided_number_of_constant_red_sites_ > -1)) {
//         return true;
//     }
//     return false;
// }

void BasePopulationTree::fold_patterns() {
    if (! this->state_frequencies_are_constrained()) {
        std::cerr << "WARNING: Site patterns are being folded when u/v rates are not constrained." << std::endl;
    }
    this->data_.fold_patterns();
    this->make_dirty();
}

void BasePopulationTree::set_freq_1(double p) {
    if (this->state_frequencies_are_fixed()) {
        return;
    }
    this->freq_1_->set_value(p);
    this->make_dirty();
}
double BasePopulationTree::get_freq_1() const {
    return this->freq_1_->get_value();
}
double BasePopulationTree::get_freq_0() const {
    return 1.0 - this->get_freq_1();
}
double BasePopulationTree::get_u() const {
    ECOEVOLITY_ASSERT(this->get_freq_1() > 0.0);
    return 1.0 / (2.0 * this->get_freq_1());
}
double BasePopulationTree::get_v() const {
    // double u = this->get_u();
    // return u / ((2.0 * u) - 1.0);
    ECOEVOLITY_ASSERT(this->get_freq_0() > 0.0);
    return 1.0 / (2.0 * this->get_freq_0());
}
void BasePopulationTree::store_freq_1() {
    this->freq_1_->store();
}
void BasePopulationTree::restore_freq_1() {
    this->freq_1_->restore();
    this->make_dirty();
}

void BasePopulationTree::set_mutation_rate(double m) {
    if (this->mutation_rate_is_fixed()) {
        return;
    }
    this->mutation_rate_->set_value(m);
    this->make_dirty();
}
double BasePopulationTree::get_mutation_rate() const {
    return this->mutation_rate_->get_value();
}
void BasePopulationTree::store_mutation_rate() {
    this->mutation_rate_->store();
}
void BasePopulationTree::restore_mutation_rate() {
    this->mutation_rate_->restore();
    this->make_dirty();
}

std::shared_ptr<PositiveRealParameter> BasePopulationTree::get_freq_1_parameter() const {
    return this->freq_1_;
}

void BasePopulationTree::set_mutation_rate_parameter(std::shared_ptr<PositiveRealParameter> h) {
    this->mutation_rate_ = h;
    this->make_dirty();
}
std::shared_ptr<PositiveRealParameter> BasePopulationTree::get_mutation_rate_parameter() const {
    return this->mutation_rate_;
}

void BasePopulationTree::set_root_population_size(double size) {
    if (this->root_->population_size_is_fixed()) {
        return;
    }
    this->root_->set_population_size(size);
}
void BasePopulationTree::set_all_population_sizes(double size) {
    if (this->population_sizes_are_fixed()) {
        return;
    }
    this->root_->set_all_population_sizes(size);
}
unsigned int BasePopulationTree::scale_all_population_sizes(double scale) {
    if (this->population_sizes_are_fixed()) {
        return 0;
    }
    unsigned int number_of_free_parameters_scaled = this->root_->scale_all_population_sizes(scale);
    return number_of_free_parameters_scaled;
}
unsigned int BasePopulationTree::scale_root_population_size(double scale) {
    if (this->root_->population_size_is_fixed()) {
        return 0;
    }
    this->root_->set_population_size(this->root_->get_population_size() * scale);
    return 1;
}

double BasePopulationTree::get_root_population_size() const {
    return this->root_->get_population_size();
}
std::shared_ptr<PositiveRealParameter> BasePopulationTree::get_root_population_size_parameter() const {
    return this->root_->get_population_size_parameter();
}

void BasePopulationTree::set_population_sizes(
        std::shared_ptr<PopulationNode> node,
        const std::vector<double> & sizes) {
    node->set_population_size(sizes.at(node->get_index()));
    for (unsigned int i = 0; i < node->get_number_of_children(); ++i) {
        this->set_population_sizes(node->get_child(i), sizes);
    }
}

void BasePopulationTree::get_population_sizes(
        std::shared_ptr<PopulationNode> node,
        std::vector<double> & sizes) const {
    sizes.at(node->get_index()) = node->get_population_size();
    for (unsigned int i = 0; i < node->get_number_of_children(); ++i) {
        this->get_population_sizes(node->get_child(i), sizes);
    }
}

std::vector<double> BasePopulationTree::get_population_sizes() const {
    std::vector<double> sizes(this->get_node_count(), 0.0);
    this->get_population_sizes(this->root_, sizes);
    return sizes;
}

std::vector< std::shared_ptr<PositiveRealParameter> > BasePopulationTree::get_pointers_to_population_sizes() const {
    return this->root_->get_all_population_size_parameters();
}

void BasePopulationTree::set_population_sizes(
        const std::vector<double> & sizes) {
    ECOEVOLITY_ASSERT(sizes.size() == this->get_node_count());
    this->set_population_sizes(this->root_, sizes);
}

std::vector<double> BasePopulationTree::get_population_sizes_as_proportions() const {
    unsigned int num_nodes = this->get_node_count();
    std::vector<double> sizes(num_nodes, 0.0);
    this->get_population_sizes(this->root_, sizes);
    double sum_size = std::accumulate(sizes.begin(), sizes.end(), 0.0);
    for (unsigned int i = 0; i < sizes.size(); ++i) {
        sizes.at(i) /= sum_size;
    }
    return sizes;
}

std::vector<double> BasePopulationTree::get_population_sizes_as_multipliers() const {
    unsigned int num_nodes = this->get_node_count();
    std::vector<double> sizes(num_nodes, 0.0);
    this->get_population_sizes(this->root_, sizes);
    double sum_size = std::accumulate(sizes.begin(), sizes.end(), 0.0);
    double mean_size = sum_size / (double)num_nodes;
    for (unsigned int i = 0; i < sizes.size(); ++i) {
        sizes.at(i) /= mean_size;
    }
    return sizes;
}

void BasePopulationTree::set_population_sizes_as_proportions(const std::vector<double> & proportions) {
    unsigned int nnodes = this->get_node_count();
    ECOEVOLITY_ASSERT(proportions.size() == nnodes);
    ECOEVOLITY_ASSERT_APPROX_EQUAL(
            std::accumulate(proportions.begin(), proportions.end(), 0.0),
            1.0);
    std::vector<double> sizes(nnodes, 0.0);
    this->get_population_sizes(this->root_, sizes);
    double sum_size = std::accumulate(sizes.begin(), sizes.end(), 0.0);
    std::vector<double> new_sizes = proportions;
    for (unsigned int i = 0; i < new_sizes.size(); ++i) {
        new_sizes.at(i) *= sum_size;
    }
    this->set_population_sizes(new_sizes);
}

double BasePopulationTree::get_mean_population_size() const {
    unsigned int num_nodes = this->get_node_count();
    std::vector<double> sizes(num_nodes, 0.0);
    this->get_population_sizes(this->root_, sizes);
    double sum_size = std::accumulate(sizes.begin(), sizes.end(), 0.0);
    return sum_size / (double)num_nodes;
}

double BasePopulationTree::get_leaf_mean_population_size() const {
    unsigned int num_nodes = this->get_node_count();
    std::vector<double> sizes(num_nodes, 0.0);
    this->get_population_sizes(this->root_, sizes);
    unsigned int nleaves = this->get_leaf_node_count();
    double sum_size = 0.0;
    for (unsigned int i = 0; i < nleaves; ++i) {
        sum_size += sizes.at(i);
    }
    return sum_size / (double)nleaves;
}

void BasePopulationTree::set_mean_population_size(double size) {
    double scale = size / this->get_mean_population_size();
    this->scale_all_population_sizes(scale);
}


void BasePopulationTree::store_derived_class_parameters() {
    // Storing of pop sizes is taken care of by PopulationNode
    this->store_freq_1();
    this->store_mutation_rate();
}
void BasePopulationTree::store_all_population_sizes() {
    this->root_->store_all_population_sizes();
}

void BasePopulationTree::restore_derived_class_parameters() {
    // Restoring of pop sizes is taken care of by PopulationNode
    this->restore_freq_1();
    this->restore_mutation_rate();
}
void BasePopulationTree::restore_all_population_sizes() {
    this->root_->restore_all_population_sizes();
}

void BasePopulationTree::set_population_size_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
    this->population_size_prior_ = prior;
    this->root_->set_all_population_size_priors(prior);
}

void BasePopulationTree::set_freq_1_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
    this->freq_1_->set_prior(prior);
    this->make_dirty();
}
void BasePopulationTree::set_mutation_rate_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
    this->mutation_rate_->set_prior(prior);
    this->make_dirty();
}

double BasePopulationTree::get_derived_class_component_of_log_prior_density() const {
    double d = 0.0;
    d += this->compute_log_prior_density_of_state_frequencies();
    d += this->compute_log_prior_density_of_mutation_rate();
    d += this->compute_log_prior_density_of_population_sizes();
    return d;
}

double BasePopulationTree::compute_log_prior_density_of_state_frequencies() const {
    return this->freq_1_->relative_prior_ln_pdf();
}
double BasePopulationTree::compute_log_prior_density_of_mutation_rate() const {
    return this->mutation_rate_->relative_prior_ln_pdf();
}
double BasePopulationTree::compute_log_prior_density_of_population_sizes() const {
    return this->root_->calculate_ln_relative_population_size_prior_density();
}

// void BasePopulationTree::provide_number_of_constant_sites(
//                 unsigned int number_all_red,
//                 unsigned int number_all_green) {
//     if (! this->constant_sites_removed_) {
//         throw EcoevolityError(
//                 "Trying to provide number of constant sites, but they haven't been removed");
//     }
//     this->provided_number_of_constant_red_sites_ = number_all_red;
//     this->provided_number_of_constant_green_sites_ = number_all_green;
// }

void BasePopulationTree::simulate_gene_tree(
        const std::shared_ptr<PopulationNode> node,
        std::unordered_map<unsigned int, std::vector< std::shared_ptr<GeneTreeSimNode> > > & branch_lineages,
        const unsigned int pattern_index,
        RandomNumberGenerator & rng,
        const bool use_max_allele_counts) const {

    // std::cout << "Top of simulate_gene_tree\n";
    // std::cout << "node index: " << node->get_index() << "\n";
    // std::cout << "node label: " << node->get_label() << "\n";
    std::vector< std::shared_ptr<GeneTreeSimNode> > lineages;
    if (node->has_children()) {
        // Handle internal branch: must get uncoalesced lineages from children
        // std::cout << "Has children!\n";
        for (unsigned int i = 0;
                i < node->get_number_of_children();
                ++i) {
            std::shared_ptr<PopulationNode> child = node->get_child(i);
            unsigned int child_idx = child->get_index();
            // std::cout << "child has index: " << "\n";
            // std::cout << "child has nlineages: " << branch_lineages.count(child_idx) << "\n";
            // ECOEVOLITY_ASSERT(child_idx >= 0);
            if (branch_lineages.count(child_idx) < 1) {
                this->simulate_gene_tree(child,
                        branch_lineages,
                        pattern_index,
                        rng,
                        use_max_allele_counts);
            }
            for (unsigned int j = 0; j < branch_lineages.at(child_idx).size(); ++j) {
                lineages.push_back(branch_lineages.at(child_idx).at(j));
            }
            branch_lineages.at(child_idx).clear();
        }
        if (node->get_number_of_parents() < 1) {
            // At root node: coalesce until 1 gene lineage and return
            // std::cout << "At root!\n";
            if (lineages.size() < 2) {
                branch_lineages[node->get_index()] = lineages;
                return;
            }
            double node_height = this->get_node_height_in_subs_per_site(*node);
            this->coalesce_in_branch(
                    lineages,
                    this->get_node_theta(*node),
                    rng,
                    node_height,
                    std::numeric_limits<double>::infinity(),
                    node->get_index()
                    );
            branch_lineages[node->get_index()] = lineages;
            return;
        }
        // Internal branch that is not the root: Coalesce from bottom to top of
        // branch
        // std::cout << "Internal node that's not root!\n";
        double node_length = this->get_node_length_in_subs_per_site(*node);
        double node_height = this->get_node_height_in_subs_per_site(*node);
        this->coalesce_in_branch(
                lineages,
                this->get_node_theta(*node),
                rng,
                node_height,
                node_height + node_length,
                node->get_index()
                );
        branch_lineages[node->get_index()] = lineages;
    } else {
        // Handle terminal branch: create genealogy tips and coalesce to top of
        // branch
        // std::cout << "Terminal node!\n";
        unsigned int allele_count;
        if (use_max_allele_counts) {
            allele_count = this->data_.get_max_allele_count(
                    node->get_index());
        }
        else {
            allele_count = this->data_.get_allele_count(
                    pattern_index,
                    node->get_index());
        }
        if (this->data_.markers_are_dominant()) {
            allele_count *= 2;
        }
        lineages.reserve(allele_count);
        for (unsigned int tip_idx = 0; tip_idx < allele_count; ++tip_idx) {
            std::shared_ptr<GeneTreeSimNode> tip = std::make_shared<GeneTreeSimNode>(
                        node->get_index(),
                        0.0);
                tip->fix_node_height();
                lineages.push_back(tip);
        }

        double node_length = this->get_node_length_in_subs_per_site(*node);
        double node_height = this->get_node_height_in_subs_per_site(*node);
        this->coalesce_in_branch(
                lineages,
                this->get_node_theta(*node),
                rng,
                node_height,
                node_height + node_length,
                node->get_index()
                );
        branch_lineages[node->get_index()] = lineages;
    }
}

std::shared_ptr<GeneTreeSimNode> BasePopulationTree::simulate_gene_tree(
        const unsigned int pattern_index,
        RandomNumberGenerator& rng,
        const bool use_max_allele_counts) const {
    std::unordered_map<unsigned int, std::vector< std::shared_ptr<GeneTreeSimNode> > > branch_lineages;
    branch_lineages.reserve(this->get_node_count());
    this->simulate_gene_tree(
            this->root_,
            branch_lineages,
            pattern_index,
            rng,
            use_max_allele_counts);
    ECOEVOLITY_ASSERT(branch_lineages.at(this->root_->get_index()).size() == 1);
    return branch_lineages.at(this->root_->get_index()).at(0);
}

double BasePopulationTree::coalesce_in_branch(
        std::vector< std::shared_ptr<GeneTreeSimNode> >& lineages,
        double population_size,
        RandomNumberGenerator& rng,
        double bottom_of_branch_height,
        double top_of_branch_height,
        unsigned int branch_index
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
                branch_index, current_height);
        for (int i = 0; i < 2; ++i) {
            int idx = rng.uniform_int(0, lineages.size() - 1);
            mrca->add_child(lineages.at(idx));
            lineages.erase(lineages.begin() + idx);
        }
        // std::cout << "coalescence!\n";
        lineages.push_back(mrca);
        --k;
        ECOEVOLITY_ASSERT(lineages.size() == k);
    }
    return current_height;
}

bool BasePopulationTree::sample_pattern(
        RandomNumberGenerator& rng,
        const float singleton_sample_probability,
        const std::vector<unsigned int>& red_allele_counts,
        const std::vector<unsigned int>& allele_counts) const {
    ECOEVOLITY_ASSERT(red_allele_counts.size() == allele_counts.size());
    if (singleton_sample_probability >= 1.0) {
        return true;
    }
    unsigned int number_of_reds = 0;
    unsigned int number_of_alleles = 0;
    for (unsigned int count_idx = 0;
            count_idx < allele_counts.size();
            ++count_idx) {
        number_of_reds += red_allele_counts.at(count_idx);
        number_of_alleles += allele_counts.at(count_idx);
    }
    if (number_of_alleles < 4) {
        throw EcoevolityBiallelicDataError(
                "Singleton filtering requires 4 or more sampled alleles for each site",
                this->data_.get_path());
    }
    if ((number_of_reds == 1) || (number_of_reds == (number_of_alleles - 1))) {
        double u = rng.uniform_real();
        if (u > singleton_sample_probability) {
            return false;
        }
    }
    return true;
}

BiallelicData BasePopulationTree::simulate_biallelic_data_set(
        RandomNumberGenerator& rng,
        float singleton_sample_probability,
        bool validate) const {
    BiallelicData sim_data = this->data_.get_empty_copy();
    const bool filtering_constant_sites = this->constant_sites_removed_;
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
                        rng,
                        false);
                auto pattern = pattern_tree.first;
                std::vector<unsigned int> red_allele_counts = pattern.first;
                std::vector<unsigned int> allele_counts = pattern.second;
                if (singleton_sample_probability < 1.0) {
                    bool sample_pattern = this->sample_pattern(
                            rng,
                            singleton_sample_probability,
                            red_allele_counts,
                            allele_counts);
                    if (! sample_pattern) {
                        continue;
                    }
                }
                std::shared_ptr<GeneTreeSimNode> gtree = pattern_tree.second;
                site_added = sim_data.add_site(red_allele_counts,
                        allele_counts,
                        filtering_constant_sites,
                        false);
            }
        }
    }
    sim_data.update_max_allele_counts();
    sim_data.update_pattern_booleans();
    if (validate) {
        sim_data.validate();
    }
    // What about gene trees? Write to stream/path (don't want to store/return
    // giant vector)
    return sim_data;
}

BiallelicData BasePopulationTree::simulate_linked_biallelic_data_set(
        RandomNumberGenerator& rng,
        float singleton_sample_probability,
        bool max_one_variable_site_per_locus,
        bool validate) const {
    ECOEVOLITY_ASSERT(this->data_.has_seq_loci_info());
    BiallelicData sim_data = this->data_.get_empty_copy();
    bool filtering_constant_sites = this->constant_sites_removed_;
    if (filtering_constant_sites) {
        // It doesn't make sense to simulate contiguous loci using a template
        // dataset that had constant sites removed
        throw EcoevolityBiallelicDataError(
                "Cannot simulate contiguous loci from charsets when constant sites were removed",
                this->data_.get_path());
    }
    sim_data.start_storing_seq_loci_info();
    if (max_one_variable_site_per_locus) {
        filtering_constant_sites = true;
        sim_data.stop_storing_seq_loci_info();
    }
    const std::vector<unsigned int> & locus_end_indices = this->data_.get_locus_end_indices();
    unsigned int site_idx = 0;
    for (unsigned int locus_idx = 0; locus_idx < locus_end_indices.size(); ++locus_idx) {
        std::shared_ptr<GeneTreeSimNode> gene_tree;
        gene_tree = this->simulate_gene_tree(0, rng, true);
        std::pair<std::vector<unsigned int>, std::vector<unsigned int> > pattern;
        while (site_idx <= locus_end_indices.at(locus_idx)) {
            bool site_added = false;
            pattern = this->simulate_biallelic_site_sans_missing(
                    gene_tree,
                    this->data_.get_allele_counts(this->data_.get_pattern_index_for_site(site_idx)),
                    rng);
            std::vector<unsigned int> red_allele_counts = pattern.first;
            std::vector<unsigned int> allele_counts = pattern.second;
            if (singleton_sample_probability < 1.0) {
                bool sample_pattern = this->sample_pattern(
                        rng,
                        singleton_sample_probability,
                        red_allele_counts,
                        allele_counts);
                if (! sample_pattern) {
                    continue;
                }
            }
            bool end_of_locus = (site_idx == locus_end_indices.at(locus_idx));
            site_added = sim_data.add_site(red_allele_counts,
                    allele_counts,
                    filtering_constant_sites,
                    end_of_locus);
            if (site_added && max_one_variable_site_per_locus) {
                site_idx = locus_end_indices.at(locus_idx);
            }
            ++site_idx;
        }
    }
    ECOEVOLITY_ASSERT(site_idx == this->data_.get_number_of_sites());
    if (! max_one_variable_site_per_locus) {
        // Debugging output
        // std::cout << "\n";
        // std::cout << "Length of end indices: "
        //         << sim_data.get_locus_end_indices().size()
        //         << " (" << this->data_.get_locus_end_indices().size()
        //         << ")\n";
        // std::cout << "End indices:\n";
        // for (auto end_idx : sim_data.get_locus_end_indices()) {
        //     std::cout << " " << end_idx;
        // }
        // std::cout << "\n";
        // std::cout << "Expected end indices:\n";
        // for (auto end_idx : this->data_.get_locus_end_indices()) {
        //     std::cout << " " << end_idx;
        // }
        // std::cout << "\n";
        ECOEVOLITY_ASSERT(this->data_.get_locus_end_indices() == sim_data.get_locus_end_indices());
    }
    sim_data.update_max_allele_counts();
    sim_data.update_pattern_booleans();
    if (validate) {
        sim_data.validate();
    }
    return sim_data;
}

std::pair<BiallelicData, unsigned int>
BasePopulationTree::simulate_complete_biallelic_data_set(
        RandomNumberGenerator& rng,
        unsigned int locus_size,
        float singleton_sample_probability,
        bool validate) const {
    ECOEVOLITY_ASSERT(locus_size > 0);
    BiallelicData sim_data = this->data_.get_empty_copy();
    const bool filtering_constant_sites = this->constant_sites_removed_;
    if (filtering_constant_sites) {
        // It doesn't make sense to simulate contiguous loci using a template
        // dataset that had constant sites removed
        throw EcoevolityBiallelicDataError(
                "Cannot simulate contiguous loci when constant sites were removed",
                this->data_.get_path());
    }
    sim_data.start_storing_seq_loci_info();
    if (locus_size < 2) {
        sim_data.stop_storing_seq_loci_info();
    }
    std::shared_ptr<GeneTreeSimNode> current_gene_tree;
    auto pattern_tree = this->simulate_biallelic_site(0, rng, true);
    auto pattern = pattern_tree.first;
    current_gene_tree = pattern_tree.second;
    unsigned int number_of_loci = 1;
    unsigned int locus_site_count = 0;
    for (unsigned int site_idx = 0;
            site_idx < this->data_.get_number_of_sites();
            ++site_idx) {
        bool site_added = false;
        while (! site_added) {
            if (locus_site_count < locus_size) {
                pattern = this->simulate_biallelic_site(
                        current_gene_tree,
                        rng);
            }
            else {
                pattern_tree = this->simulate_biallelic_site(0, rng, true);
                pattern = pattern_tree.first;
                current_gene_tree = pattern_tree.second;
                locus_site_count = 0;
                ++number_of_loci;
            }
            std::vector<unsigned int> red_allele_counts = pattern.first;
            std::vector<unsigned int> allele_counts = pattern.second;
            if (singleton_sample_probability < 1.0) {
                bool sample_pattern = this->sample_pattern(
                        rng,
                        singleton_sample_probability,
                        red_allele_counts,
                        allele_counts);
                if (! sample_pattern) {
                    continue;
                }
            }
            bool end_of_locus = (
                    (
                        (locus_site_count == (locus_size - 1)) ||
                        (site_idx == (this->data_.get_number_of_sites() - 1))
                    ) &&
                    (locus_size > 1)
                    );
            site_added = sim_data.add_site(red_allele_counts,
                    allele_counts,
                    filtering_constant_sites,
                    end_of_locus);
        }
        ++locus_site_count;
    }
    sim_data.update_max_allele_counts();
    sim_data.update_pattern_booleans();
    if (validate) {
        sim_data.validate();
    }
    return std::make_pair(sim_data, number_of_loci);
}

std::pair<BiallelicData, unsigned int>
BasePopulationTree::simulate_data_set_max_one_variable_site_per_locus(
        RandomNumberGenerator& rng,
        unsigned int locus_size,
        float singleton_sample_probability,
        bool validate) const {
    ECOEVOLITY_ASSERT(locus_size > 0);
    BiallelicData sim_data = this->data_.get_empty_copy();
    if (this->constant_sites_removed_) {
        throw EcoevolityError(
                "To simulated one variable site per locus, the \"template\" "
                "dataset must include constant sites");
    }
    const bool filtering_constant_sites = true;
    sim_data.stop_storing_seq_loci_info();
    std::shared_ptr<GeneTreeSimNode> current_gene_tree;
    auto pattern_tree = this->simulate_biallelic_site(0, rng, true);
    auto pattern = pattern_tree.first;
    current_gene_tree = pattern_tree.second;
    unsigned int number_of_loci = 1;
    unsigned int locus_site_count = 0;
    bool site_added = false;
    for (unsigned int site_idx = 0;
            site_idx < this->data_.get_number_of_sites();
            ++site_idx) {
        if (locus_site_count < locus_size) {
            pattern = this->simulate_biallelic_site(
                    current_gene_tree,
                    rng);
        }
        else {
            pattern_tree = this->simulate_biallelic_site(0, rng, true);
            pattern = pattern_tree.first;
            current_gene_tree = pattern_tree.second;
            locus_site_count = 0;
            ++number_of_loci;
        }
        std::vector<unsigned int> red_allele_counts = pattern.first;
        std::vector<unsigned int> allele_counts = pattern.second;
        if (singleton_sample_probability < 1.0) {
            bool sample_pattern = this->sample_pattern(
                    rng,
                    singleton_sample_probability,
                    red_allele_counts,
                    allele_counts);
            if (! sample_pattern) {
                site_idx -= 1;
                continue;
            }
        }
        site_added = sim_data.add_site(red_allele_counts,
                allele_counts,
                filtering_constant_sites,
                false);
        if (site_added) {
            // We have a variable site for this locus, so need to skip over
            // the remaining sites of this locus.
            // Subtracting 1 for the loop increment.
            int remaining_locus_sites = locus_size - locus_site_count;
            site_idx += (remaining_locus_sites - 1);
            locus_site_count += remaining_locus_sites;
            continue;
        }
        ++locus_site_count;
    }
    sim_data.update_max_allele_counts();
    sim_data.update_pattern_booleans();
    if (validate) {
        sim_data.validate();
    }
    return std::make_pair(sim_data, number_of_loci);
}

std::pair<
        std::pair<std::vector<unsigned int>, std::vector<unsigned int> >,
        std::shared_ptr<GeneTreeSimNode> >
BasePopulationTree::simulate_biallelic_site(
        const unsigned int pattern_idx,
        RandomNumberGenerator& rng,
        const bool use_max_allele_counts) const {

    std::shared_ptr<GeneTreeSimNode> gene_tree = this->simulate_gene_tree(
            pattern_idx,
            rng,
            use_max_allele_counts);
    gene_tree->compute_binary_transition_probabilities(this->get_u(), this->get_v());
    gene_tree->simulate_binary_character(this->get_freq_0(), rng);

    std::vector<unsigned int> expected_allele_counts;
    if (use_max_allele_counts) {
        expected_allele_counts = this->data_.get_max_allele_counts();
    }
    else {
        expected_allele_counts = this->data_.get_allele_counts(pattern_idx);
    }
    std::vector<unsigned int> allele_counts(expected_allele_counts.size(), 0);
    std::vector<unsigned int> red_allele_counts(expected_allele_counts.size(), 0);
    if (this->data_.markers_are_dominant()) {
        std::vector<int> last_allele(expected_allele_counts.size(), -1);
        gene_tree->get_allele_counts(
                allele_counts,
                red_allele_counts,
                last_allele);
    }
    else {
        gene_tree->get_allele_counts(
                allele_counts,
                red_allele_counts);
    }

    ECOEVOLITY_ASSERT(allele_counts == expected_allele_counts);
    std::pair<std::vector<unsigned int>, std::vector<unsigned int> > pattern = 
            std::make_pair(red_allele_counts, allele_counts);
    return std::make_pair(pattern, gene_tree);
}

std::pair<std::vector<unsigned int>, std::vector<unsigned int> >
BasePopulationTree::simulate_biallelic_site(
        std::shared_ptr<GeneTreeSimNode> gene_tree,
        RandomNumberGenerator& rng) const {
    gene_tree->compute_binary_transition_probabilities(this->get_u(), this->get_v());
    gene_tree->simulate_binary_character(this->get_freq_0(), rng);

    const std::vector<unsigned int>& expected_allele_counts = this->data_.get_max_allele_counts();
    std::vector<unsigned int> allele_counts(expected_allele_counts.size(), 0);
    std::vector<unsigned int> red_allele_counts(expected_allele_counts.size(), 0);
    if (this->data_.markers_are_dominant()) {
        std::vector<int> last_allele(expected_allele_counts.size(), -1);
        gene_tree->get_allele_counts(
                allele_counts,
                red_allele_counts,
                last_allele);
    }
    else {
        gene_tree->get_allele_counts(
                allele_counts,
                red_allele_counts);
    }

    ECOEVOLITY_ASSERT(allele_counts == expected_allele_counts);
    std::pair<std::vector<unsigned int>, std::vector<unsigned int> > pattern = 
            std::make_pair(red_allele_counts, allele_counts);
    return pattern;
}

std::pair<std::vector<unsigned int>, std::vector<unsigned int> >
BasePopulationTree::simulate_biallelic_site_sans_missing(
        std::shared_ptr<GeneTreeSimNode> gene_tree,
        const std::vector<unsigned int> & site_allele_counts,
        RandomNumberGenerator& rng) const {
    auto pattern = this->simulate_biallelic_site(
            gene_tree,
            rng);
    std::vector<unsigned int> red_allele_counts = pattern.first;
    std::vector<unsigned int> allele_counts = pattern.second;
    ECOEVOLITY_ASSERT(site_allele_counts.size() == allele_counts.size());
    std::vector<unsigned int> sampled_allele_counts (allele_counts);
    std::vector<unsigned int> sampled_red_allele_counts (red_allele_counts);
    for (unsigned int pop_idx = 0; pop_idx < allele_counts.size(); ++pop_idx) {
        if (allele_counts.at(pop_idx) == site_allele_counts.at(pop_idx)) {
            // No gene copies to prune
            continue;
        }
        double prob_pick_red = (double)sampled_red_allele_counts.at(pop_idx) / (double)sampled_allele_counts.at(pop_idx);
        for (unsigned int missing_idx = 0;
                missing_idx < (allele_counts.at(pop_idx) - site_allele_counts.at(pop_idx));
                ++missing_idx) {
            double u = rng.uniform_real();
            if (u < prob_pick_red) {
                --sampled_red_allele_counts.at(pop_idx);
            }
            --sampled_allele_counts.at(pop_idx);
            prob_pick_red = (double)sampled_red_allele_counts.at(pop_idx) / (double)sampled_allele_counts.at(pop_idx);
        }
    }
    ECOEVOLITY_ASSERT(site_allele_counts == sampled_allele_counts);
    std::pair<std::vector<unsigned int>, std::vector<unsigned int> > sampled_pattern =
            std::make_pair(sampled_red_allele_counts, sampled_allele_counts);
    return sampled_pattern;
}

double BasePopulationTree::compute_log_likelihood(
        unsigned int nthreads) {
    if (this->ignoring_data()) {
        this->log_likelihood_.set_value(0.0);
        return 0.0;
    }
    double constant_pattern_lnl_correction = 0.0;
    double log_likelihood = get_log_likelihood(
            this->get_mutable_root(),
            this->data_.get_red_allele_count_matrix(),
            this->data_.get_allele_count_matrix(),
            this->data_.get_pattern_weights(),
            this->unique_allele_counts_,
            this->unique_allele_count_weights_,
            this->get_u(),
            this->get_v(),
            this->get_mutation_rate(),
            this->get_ploidy(),
            this->data_.markers_are_dominant(),
            this->state_frequencies_are_constrained(),
            this->constant_sites_removed(),
            constant_pattern_lnl_correction,
            nthreads);

    if (this->constant_sites_removed()) {
        //////////////////////////////////////////////////////////////////////
        // No reason to use removed site counts. Simply leave constant sites in
        // and calc likelihood without correction. This is better, because it
        // doesn't treat all constant site patterns equally (i.e., it accounts
        // for constant patterns with missing data).
        // if (this->constant_site_counts_were_provided()) {
        //     double constant_log_likelihood =
        //             ((double)this->get_provided_number_of_constant_green_sites() * std::log(all_green_pattern_prob)) +
        //             ((double)this->get_provided_number_of_constant_red_sites() * std::log(all_red_pattern_prob));
        //     log_likelihood += constant_log_likelihood;
        // }
        // else if (this->use_removed_constant_site_counts_){
        //     double constant_log_likelihood =
        //             ((double)this->data_.get_number_of_constant_green_sites_removed() *
        //             std::log(this->all_green_pattern_likelihood_.get_value())) +
        //             ((double)this->data_.get_number_of_constant_red_sites_removed() *
        //             std::log(this->all_red_pattern_likelihood_.get_value()));
        //     log_likelihood += constant_log_likelihood;
        // }
        //////////////////////////////////////////////////////////////////////
        // else {
        if (constant_pattern_lnl_correction == -std::numeric_limits<double>::infinity()) {
            // Rather than throw an error return -inf and let the MCMC machinery reject
            // TODO: Is there a better way to handle this? Technically, the log likelihood
            // would be inf or NAN (not -inf)
            log_likelihood = -std::numeric_limits<double>::infinity();
            double root_height = this->get_node_height_in_subs_per_site(this->get_root());
            std::vector<double> pop_sizes = this->get_population_sizes();
            std::ostringstream message;
            message << "\n"
                    << "\n#######################################################################\n"
                    <<   "##############################  WARNING  ##############################\n"
                    <<   "The probability of a variable character is zero for the current state\n"
                    <<   "of the population-tree model for the data in:\n    \'"
                    <<   this->data_.get_path() << "\'.\n"
                    <<   "Correcting the likelihood for missing constant characters would thus\n"
                    <<   "result an infinite likelihood for any character pattern.\n"
                    <<   "This is likely due to the event time and population sizes being very\n"
                    <<   "small. The current height of the root node in expected subsitutions\n"
                    <<   "per site is:\n    "
                    <<   root_height * this->get_mutation_rate() << "\n"
                    <<   "The current population sizes (scaled by the mutation rate) are:\n    "
                    <<   pop_sizes.at(0) * this->get_mutation_rate();
            for (unsigned int i = 1; i < pop_sizes.size(); ++i) {
                message << " " << pop_sizes.at(i) * this->get_mutation_rate();
            }
            message << "\n"
                    << "This state will be rejected by the Metropolis-Hastings algorithm,\n"
                    << "however, the MCMC exploring such parameter space could indicate a\n"
                    << "larger problem, such as a prior specified in the wrong units\n"
                    << "#######################################################################\n\n";
            std::cerr << message.str();
            // throw EcoevolityError(message.str());
        }
        else {
            log_likelihood -= constant_pattern_lnl_correction;
        }
        // }
    }

    log_likelihood += this->get_likelihood_correction();

    // ECOEVOLITY_DEBUG(
    //     std::cerr << "compute_log_likelihood(): " << log_likelihood << std::endl;
    // )
    if (std::isnan(log_likelihood)) {
        throw EcoevolityError("BasePopulationTree::compute_log_likelihood resulted in a NAN likelihood");
    }

    this->log_likelihood_.set_value(log_likelihood);
    return log_likelihood;
}




void PopulationTree::set_population_size_multiplier_prior(std::shared_ptr<DirichletDistribution> prior) {
    this->population_size_multiplier_prior_ = prior;
    this->make_dirty();
}

void PopulationTree::set_population_sizes(
        const std::vector<double> & sizes) {
    if (this->population_size_multipliers_are_fixed()) {
        return;
    }
    BasePopulationTree::set_population_sizes(sizes);
}

void PopulationTree::set_population_sizes_as_proportions(const std::vector<double> & proportions) {
    if (this->population_size_multipliers_are_fixed()) {
        return;
    }
    BasePopulationTree::set_population_sizes_as_proportions(proportions);
}

void PopulationTree::set_mean_population_size(double size) {
    if (this->mean_population_size_is_fixed()) {
        return;
    }
    BasePopulationTree::set_mean_population_size(size);
}

void PopulationTree::set_root_height_parameter(std::shared_ptr<PositiveRealParameter> h) {
    ECOEVOLITY_ASSERT(h->prior == this->root_node_height_prior_);
    this->root_->set_height_parameter(h);
}
std::shared_ptr<PositiveRealParameter> PopulationTree::get_root_height_parameter() const {
    return this->root_->get_height_parameter();
}

void PopulationTree::store_root_height() {
    this->root_->store_height();
}
void PopulationTree::restore_root_height() {
    this->root_->restore_height();
}

void PopulationTree::store_parameters() {
    this->store_freq_1();
    this->store_mutation_rate();
    this->store_all_population_sizes();
    this->store_all_heights();
}

void PopulationTree::restore_parameters() {
    this->restore_freq_1();
    this->restore_mutation_rate();
    this->restore_all_population_sizes();
    this->restore_all_heights();
}

double PopulationTree::compute_log_prior_density() {
    double d = 0.0;
    d += this->compute_log_prior_density_of_state_frequencies();
    d += this->compute_log_prior_density_of_mutation_rate();
    d += this->compute_log_prior_density_of_node_heights();
    d += this->compute_log_prior_density_of_population_sizes();
    this->log_prior_density_.set_value(d);
    return d;
}

double PopulationTree::compute_log_prior_density_of_node_heights() const {
    return this->root_->calculate_ln_relative_node_height_prior_density();
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
        bool strict_on_missing_sites,
        bool strict_on_triallelic_sites,
        double ploidy,
        bool store_seq_loci_info) {
    this->comparison_init(path,
               population_name_delimiter,
               population_name_is_prefix,
               genotypes_are_diploid,
               markers_are_dominant,
               constant_sites_removed,
               validate,
               strict_on_constant_sites,
               strict_on_missing_sites,
               strict_on_triallelic_sites,
               ploidy,
               store_seq_loci_info);
}
void ComparisonPopulationTree::comparison_init(
        std::string path, 
        char population_name_delimiter,
        bool population_name_is_prefix,
        bool genotypes_are_diploid,
        bool markers_are_dominant,
        bool constant_sites_removed,
        bool validate,
        bool strict_on_constant_sites,
        bool strict_on_missing_sites,
        bool strict_on_triallelic_sites,
        double ploidy,
        bool store_seq_loci_info) {
    this->init(path,
               population_name_delimiter,
               population_name_is_prefix,
               genotypes_are_diploid,
               markers_are_dominant,
               constant_sites_removed,
               validate,
               strict_on_constant_sites,
               strict_on_missing_sites,
               strict_on_triallelic_sites,
               ploidy,
               store_seq_loci_info);
    if (this->data_.get_number_of_populations() > 2) {
        throw EcoevolityComparisonSettingError(
                "ComparisonPopulationTree() does not support more than 2 populations",
                path);
    }
    this->root_->set_label("root-" + this->root_->get_child(0)->get_label());
}
ComparisonPopulationTree::ComparisonPopulationTree(
        const ComparisonSettings& settings,
        RandomNumberGenerator& rng,
        bool strict_on_constant_sites,
        bool strict_on_missing_sites, 
        bool strict_on_triallelic_sites,
        bool store_seq_loci_info) {
    this->comparison_init(settings.get_path(),
               settings.get_population_name_delimiter(),
               settings.population_name_is_prefix(),
               settings.genotypes_are_diploid(),
               settings.markers_are_dominant(),
               settings.constant_sites_removed(),
               true, // validate
               strict_on_constant_sites,
               strict_on_missing_sites,
               strict_on_triallelic_sites,
               settings.get_ploidy(),
               store_seq_loci_info);
    if (settings.constrain_state_frequencies()) {
        this->constrain_state_frequencies();
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
    this->set_all_population_sizes(p.get_value());
    if (settings.get_population_size_settings().is_fixed()) {
        this->fix_population_sizes();
    }
    
    this->set_freq_1_prior(settings.get_freq_1_settings().get_prior_settings().get_instance());
    if (settings.constrain_state_frequencies()) {
        this->constrain_state_frequencies();
    }
    else {
        PositiveRealParameter freq_1 = PositiveRealParameter(
                settings.get_freq_1_settings(),
                rng);
        this->set_freq_1(freq_1.get_value());
        if (freq_1.is_fixed()) {
            this->fix_state_frequencies();
        }
    }
    this->set_mutation_rate_parameter(
            std::make_shared<PositiveRealParameter>(
                    settings.get_mutation_rate_settings(),
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
                << "contains only a single population, but you have fixed and/or\n"
                << "constrained the population sizes for this comparison. The timing of\n"
                << "population expansion/contraction cannot be estimated if the ancestral\n"
                << "and descendant population sizes for this comparison are either\n"
                << "fixed or constrained to be equal. Please update your configuration\n"
                << "file to estimate the unconstrained population sizes for this\n"
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
double ComparisonPopulationTree::get_child_population_size(
        unsigned int child_index) const {
    return this->root_->get_child(child_index)->get_population_size();
}
std::shared_ptr<PositiveRealParameter> ComparisonPopulationTree::get_child_population_size_parameter(
        unsigned int child_index) const {
    return this->root_->get_child(child_index)->get_population_size_parameter();
}

// Node height sharing needs to be dealt with in next level up in
// class hierarchy (ComparisonPopulationTreeCollection)
double ComparisonPopulationTree::compute_log_prior_density() {
    double d = 0.0;
    d += this->compute_log_prior_density_of_state_frequencies();
    d += this->compute_log_prior_density_of_mutation_rate();
    // d += this->compute_log_prior_density_of_node_heights();
    d += this->compute_log_prior_density_of_population_sizes();
    this->log_prior_density_.set_value(d);
    return d;
}

// Node height (re)storing is managed by ComparisonPopulationTree.
void ComparisonPopulationTree::store_parameters() {
    this->store_freq_1();
    this->store_mutation_rate();
    this->store_all_population_sizes();
    // this->store_all_heights();
}
void ComparisonPopulationTree::restore_parameters() {
    this->restore_freq_1();
    this->restore_mutation_rate();
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
        << "mutation_rate" << suffix << delimiter
        << "freq_1" << suffix << delimiter
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
        << this->get_root_height() << delimiter
        << this->get_mutation_rate() << delimiter
        << this->get_freq_1() << delimiter
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

void ComparisonPopulationTree::draw_from_prior(RandomNumberGenerator& rng) {
    if ((! this->state_frequencies_are_fixed()) && (! this->state_frequencies_are_constrained())) {
        this->freq_1_->set_value_from_prior(rng);
    }
    if (! this->mutation_rate_is_fixed()) {
        this->mutation_rate_->set_value_from_prior(rng);
    }
    if (! this->population_sizes_are_fixed()) {
        this->set_root_population_size(this->population_size_prior_->draw(rng));
        if (! this->population_sizes_are_constrained()) {
            this->set_child_population_size(0, this->population_size_prior_->draw(rng));
            if (this->root_->get_number_of_children() > 1) {
                this->set_child_population_size(1, this->population_size_prior_->draw(rng));
            }
        }
    }
}



DirichletPopulationTree::DirichletPopulationTree(
        std::string path, 
        char population_name_delimiter,
        bool population_name_is_prefix,
        bool genotypes_are_diploid,
        bool markers_are_dominant,
        bool constant_sites_removed,
        bool validate,
        bool strict_on_constant_sites,
        bool strict_on_missing_sites,
        bool strict_on_triallelic_sites,
        double ploidy,
        bool store_seq_loci_info) {
    this->init(path,
               population_name_delimiter,
               population_name_is_prefix,
               genotypes_are_diploid,
               markers_are_dominant,
               constant_sites_removed,
               validate,
               strict_on_constant_sites,
               strict_on_missing_sites,
               strict_on_triallelic_sites,
               ploidy,
               store_seq_loci_info);
    std::vector<double> dirichlet_parameters (this->get_node_count(), 1.0);
    this->set_population_size_multiplier_prior(std::make_shared<DirichletDistribution>(dirichlet_parameters));
}

double DirichletPopulationTree::compute_log_prior_density_of_population_sizes() const {
    if (this->population_sizes_are_fixed()) {
        return 0.0;
    }
    double ln_p = 0.0;
    if (! this->mean_population_size_is_fixed()) {
        ln_p += this->population_size_prior_->relative_ln_pdf(
                this->get_mean_population_size());
    }
    if (! this->population_size_multipliers_are_fixed()) {
        ln_p += this->population_size_multiplier_prior_->relative_ln_pdf(
                this->get_population_sizes_as_proportions());
    }
    return ln_p;
}



ComparisonDirichletPopulationTree::ComparisonDirichletPopulationTree(
        std::string path, 
        char population_name_delimiter,
        bool population_name_is_prefix,
        bool genotypes_are_diploid,
        bool markers_are_dominant,
        bool constant_sites_removed,
        bool validate,
        bool strict_on_constant_sites,
        bool strict_on_missing_sites,
        bool strict_on_triallelic_sites,
        double ploidy,
        bool store_seq_loci_info) {
    this->comparison_init(path,
               population_name_delimiter,
               population_name_is_prefix,
               genotypes_are_diploid,
               markers_are_dominant,
               constant_sites_removed,
               validate,
               strict_on_constant_sites,
               strict_on_missing_sites,
               strict_on_triallelic_sites,
               ploidy,
               store_seq_loci_info);
    std::vector<double> dirichlet_parameters (this->get_node_count(), 1.0);
    this->set_population_size_multiplier_prior(std::make_shared<DirichletDistribution>(dirichlet_parameters));
}

void ComparisonDirichletPopulationTree::comparison_init(
        std::string path, 
        char population_name_delimiter,
        bool population_name_is_prefix,
        bool genotypes_are_diploid,
        bool markers_are_dominant,
        bool constant_sites_removed,
        bool validate,
        bool strict_on_constant_sites,
        bool strict_on_missing_sites,
        bool strict_on_triallelic_sites,
        double ploidy,
        bool store_seq_loci_info) {
    this->init(path,
               population_name_delimiter,
               population_name_is_prefix,
               genotypes_are_diploid,
               markers_are_dominant,
               constant_sites_removed,
               validate,
               strict_on_constant_sites,
               strict_on_missing_sites,
               strict_on_triallelic_sites,
               ploidy,
               store_seq_loci_info);
    if (this->data_.get_number_of_populations() > 2) {
        throw EcoevolityComparisonSettingError(
                "ComparisonDirichletPopulationTree() does not support more than 2 populations",
                path);
    }
    this->root_->set_label("root-" + this->root_->get_child(0)->get_label());
}

ComparisonDirichletPopulationTree::ComparisonDirichletPopulationTree(
        const DirichletComparisonSettings& settings,
        RandomNumberGenerator& rng,
        bool strict_on_constant_sites,
        bool strict_on_missing_sites,
        bool strict_on_triallelic_sites,
        bool store_seq_loci_info) {
    this->comparison_init(settings.get_path(),
               settings.get_population_name_delimiter(),
               settings.population_name_is_prefix(),
               settings.genotypes_are_diploid(),
               settings.markers_are_dominant(),
               settings.constant_sites_removed(),
               true, // validate
               strict_on_constant_sites,
               strict_on_missing_sites,
               strict_on_triallelic_sites,
               settings.get_ploidy(),
               store_seq_loci_info);
    if (settings.constrain_state_frequencies()) {
        this->constrain_state_frequencies();
        this->fold_patterns();
    }
    this->set_population_size_prior(
            settings.get_population_size_settings().get_prior_settings().get_instance());
    PositiveRealParameter p = PositiveRealParameter(
            settings.get_population_size_settings(),
            rng);
    this->set_mean_population_size(p.get_value());
    if (settings.get_population_size_settings().is_fixed()) {
        this->fix_mean_population_size();
    }

    this->set_population_size_multiplier_prior(
            settings.get_population_size_multiplier_settings().get_prior_settings().get_dirichlet_distribution_instance());
    std::vector<double> multiplier_values = settings.get_population_size_multiplier_settings().get_values();
    if (multiplier_values.size() < 1) {
        multiplier_values = this->population_size_multiplier_prior_->draw(rng);
        this->set_population_sizes_as_proportions(multiplier_values);
    }
    else {
        for (unsigned int i = 0; i < multiplier_values.size(); ++i) {
            multiplier_values.at(i) /= multiplier_values.size();
        }
        this->set_population_sizes_as_proportions(multiplier_values);
        if (settings.get_population_size_multiplier_settings().is_fixed()) {
            this->fix_population_size_multipliers();
        }
    }
    
    this->set_freq_1_prior(settings.get_freq_1_settings().get_prior_settings().get_instance());
    if (settings.constrain_state_frequencies()) {
        this->constrain_state_frequencies();
    }
    else {
        PositiveRealParameter freq_1 = PositiveRealParameter(
                settings.get_freq_1_settings(),
                rng);
        this->set_freq_1(freq_1.get_value());
        if (freq_1.is_fixed()) {
            this->fix_state_frequencies();
        }
    }
    this->set_mutation_rate_parameter(
            std::make_shared<PositiveRealParameter>(
                    settings.get_mutation_rate_settings(),
                    rng));
    if ((this->data_.get_number_of_populations() == 1) &&
            (this->population_size_multipliers_are_fixed()) &&
            (this->get_root_population_size() == this->get_child_population_size(0)))
    {
        std::ostringstream message;
        message << "\n#######################################################################\n"
                <<   "###############################  ERROR  ###############################\n"
                << "The alignment in:\n    \'"
                << this->data_.get_path() << "\'\n"
                << "contains only a single population, but you have fixed the population\n"
                << "sizes for this comparison to be equal. The timing of population\n"
                << "expansion/contraction cannot be estimated if the ancestral and\n"
                << "descendant population sizes for this comparison are fixed to be equal.\n"
                << "Please update your configuration for this comparison accordingly and\n"
                << "re-run the analysis.\n"
                << "#######################################################################\n";
        throw EcoevolityComparisonSettingError(message.str(), this->data_.get_path());
    }
}

void ComparisonDirichletPopulationTree::write_state_log_header(
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
        << "mutation_rate" << suffix << delimiter
        << "freq_1" << suffix << delimiter
        << "pop_size" << suffix << delimiter
        << "pop_size_multiplier" << suffix << delimiter;
    if (this->root_->get_number_of_children() > 1) {
        out << "pop_size_multiplier" << "_" << this->root_->get_child(1)->get_label() << delimiter;
    }
    out << "pop_size_multiplier_root" << suffix;
}

void ComparisonDirichletPopulationTree::log_state(
        std::ostream& out,
        const std::string& delimiter) const {
    std::vector<double> multipliers = this->get_population_sizes_as_multipliers();
    out << this->log_likelihood_.get_value() << delimiter
        << this->log_prior_density_.get_value() << delimiter
        << this->get_root_height() << delimiter
        << this->get_mutation_rate() << delimiter
        << this->get_freq_1() << delimiter
        << this->get_mean_population_size() << delimiter
        << multipliers.at(0) << delimiter;
    unsigned int root_index = 1;
    if (this->root_->get_number_of_children() > 1) {
        out << multipliers.at(1) << delimiter;
        ++root_index;
    }
    out << multipliers.at(root_index);
}

void ComparisonDirichletPopulationTree::log_state(
        std::ostream& out,
        unsigned int event_index,
        const std::string& delimiter) const {
    out << event_index << delimiter;
    this->log_state(out, delimiter);
}

void ComparisonDirichletPopulationTree::draw_from_prior(RandomNumberGenerator& rng) {
    if ((! this->state_frequencies_are_fixed()) && (! this->state_frequencies_are_constrained())) {
        this->freq_1_->set_value_from_prior(rng);
    }
    if (! this->mutation_rate_is_fixed()) {
        this->mutation_rate_->set_value_from_prior(rng);
    }
    if (! this->mean_population_size_is_fixed()) {
        this->set_mean_population_size(this->population_size_prior_->draw(rng));
    }
    if (! this->population_size_multipliers_are_fixed()) {
        this->set_population_sizes_as_proportions(this->population_size_multiplier_prior_->draw(rng));
    }
}

double ComparisonDirichletPopulationTree::compute_log_prior_density_of_population_sizes() const {
    if (this->population_sizes_are_fixed()) {
        return 0.0;
    }
    double ln_p = 0.0;
    if (! this->mean_population_size_is_fixed()) {
        ln_p += this->population_size_prior_->relative_ln_pdf(
                this->get_mean_population_size());
    }
    if (! this->population_size_multipliers_are_fixed()) {
        ln_p += this->population_size_multiplier_prior_->relative_ln_pdf(
                this->get_population_sizes_as_proportions());
    }
    return ln_p;
}



//////////////////////////////////////////////////////////////////////////////
// RelativeRootPopulationTree methods
//////////////////////////////////////////////////////////////////////////////

void RelativeRootPopulationTree::update_root_population_size() {
    this->root_->set_population_size(
            this->relative_root_population_size_->get_value() *
            this->get_leaf_mean_population_size());
    this->make_dirty();
}

void RelativeRootPopulationTree::update_relative_root_population_size() {
    this->relative_root_population_size_->set_value(
            this->root_->get_population_size() /
            this->get_leaf_mean_population_size());
    this->make_dirty();
}

void RelativeRootPopulationTree::set_root_population_size(double size) {
    if (this->population_sizes_are_constrained()) {
        PopulationTree::set_root_population_size(size);
        return;
    }
    if (this->relative_root_population_size_is_fixed()) {
        return;
    }
    ECOEVOLITY_ASSERT(! this->root_->population_size_is_fixed());
    PopulationTree::set_root_population_size(size);
    this->update_relative_root_population_size();
}

double RelativeRootPopulationTree::get_root_population_size() const {
    if (this->population_sizes_are_constrained()) {
        return PopulationTree::get_root_population_size();
    }
    double size = this->root_->get_population_size();
    ECOEVOLITY_ASSERT_APPROX_EQUAL(size, 
            (this->relative_root_population_size_->get_value() *
             this->get_leaf_mean_population_size()));
    return size;
}

void RelativeRootPopulationTree::set_all_population_sizes(
        double size) {
    PopulationTree::set_all_population_sizes(size);
    if (this->population_sizes_are_fixed() || this->population_sizes_are_constrained()) {
        return;
    }
    if (this->relative_root_population_size_is_fixed()) {
        this->update_root_population_size();
        return;
    }
    this->update_relative_root_population_size();
}

unsigned int RelativeRootPopulationTree::scale_all_population_sizes(
        double scale) {
    unsigned int ret = PopulationTree::scale_all_population_sizes(scale);
    if (this->population_sizes_are_fixed() || this->population_sizes_are_constrained()) {
        return ret;
    }
    if (! this->leaf_population_sizes_are_fixed_) {
        --ret;
    }
    if (this->relative_root_population_size_is_fixed()) {
        this->update_root_population_size();
        return ret;
    }
    this->update_relative_root_population_size();
    return ret;
}

unsigned int RelativeRootPopulationTree::scale_root_population_size(
        double scale) {
    if (this->population_sizes_are_constrained()) {
        return PopulationTree::scale_root_population_size(scale);
    }
    if (this->relative_root_population_size_is_fixed()) {
        return 0;
    }
    ECOEVOLITY_ASSERT(! this->root_->population_size_is_fixed());
    unsigned int ret = PopulationTree::scale_root_population_size(scale);
    this->update_relative_root_population_size();
    return ret;
}

void RelativeRootPopulationTree::set_mean_population_size(double size) {
    PopulationTree::set_mean_population_size(size);
    if (this->population_sizes_are_fixed() || this->population_sizes_are_constrained()) {
        return;
    }
    if (this->relative_root_population_size_is_fixed()) {
        this->update_root_population_size();
        return;
    }
    this->update_relative_root_population_size();
}

void RelativeRootPopulationTree::set_population_sizes_as_proportions(
        const std::vector<double> & proportions) {
    PopulationTree::set_population_sizes_as_proportions(proportions);
    if (this->population_sizes_are_fixed() || this->population_sizes_are_constrained()) {
        return;
    }
    if (this->relative_root_population_size_is_fixed()) {
        this->update_root_population_size();
        return;
    }
    this->update_relative_root_population_size();
}

void RelativeRootPopulationTree::set_population_sizes(const std::vector<double> & sizes) {
    PopulationTree::set_population_sizes(sizes);
    if (this->population_sizes_are_fixed() || this->population_sizes_are_constrained()) {
        return;
    }
    if (this->relative_root_population_size_is_fixed()) {
        this->update_root_population_size();
        return;
    }
    this->update_relative_root_population_size();
}

void RelativeRootPopulationTree::set_relative_root_population_size_prior(
        std::shared_ptr<ContinuousProbabilityDistribution> prior) {
    this->relative_root_population_size_->prior = prior;
    this->make_dirty();
}

void RelativeRootPopulationTree::store_all_population_sizes() {
    this->root_->store_all_population_sizes();
    this->relative_root_population_size_->store();
}

void RelativeRootPopulationTree::restore_all_population_sizes() {
    this->root_->restore_all_population_sizes();
    this->relative_root_population_size_->restore();
}

double RelativeRootPopulationTree::compute_log_prior_density_of_population_sizes() const {
    if (this->population_sizes_are_fixed() || this->population_sizes_are_constrained()) {
        return PopulationTree::compute_log_prior_density_of_population_sizes();
    }
    double ln_p = 0.0;
    if (! this->leaf_population_sizes_are_fixed_) {
        // Can't call relative_ln_pdf when pop sizes were fixed, because root
        // size will not have a prior
        ln_p = this->root_->calculate_ln_relative_population_size_prior_density();
        ln_p -= this->root_->get_population_size_relative_prior_ln_pdf();
    }
    ln_p += this->relative_root_population_size_->relative_prior_ln_pdf();
    return ln_p;
}


//////////////////////////////////////////////////////////////////////////////
// ComparisonRelativeRootPopulationTree methods
//////////////////////////////////////////////////////////////////////////////

ComparisonRelativeRootPopulationTree::ComparisonRelativeRootPopulationTree(
        std::string path, 
        char population_name_delimiter,
        bool population_name_is_prefix,
        bool genotypes_are_diploid,
        bool markers_are_dominant,
        bool constant_sites_removed,
        bool validate,
        bool strict_on_constant_sites,
        bool strict_on_missing_sites,
        bool strict_on_triallelic_sites,
        double ploidy,
        bool store_seq_loci_info) {
    this->comparison_init(path,
               population_name_delimiter,
               population_name_is_prefix,
               genotypes_are_diploid,
               markers_are_dominant,
               constant_sites_removed,
               validate,
               strict_on_constant_sites,
               strict_on_missing_sites,
               strict_on_triallelic_sites,
               ploidy,
               store_seq_loci_info);
}
void ComparisonRelativeRootPopulationTree::comparison_init(
        std::string path, 
        char population_name_delimiter,
        bool population_name_is_prefix,
        bool genotypes_are_diploid,
        bool markers_are_dominant,
        bool constant_sites_removed,
        bool validate,
        bool strict_on_constant_sites,
        bool strict_on_missing_sites,
        bool strict_on_triallelic_sites,
        double ploidy,
        bool store_seq_loci_info) {
    this->init(path,
               population_name_delimiter,
               population_name_is_prefix,
               genotypes_are_diploid,
               markers_are_dominant,
               constant_sites_removed,
               validate,
               strict_on_constant_sites,
               strict_on_missing_sites,
               strict_on_triallelic_sites,
               ploidy,
               store_seq_loci_info);
    if (this->data_.get_number_of_populations() > 2) {
        throw EcoevolityComparisonSettingError(
                "ComparisonRelativeRootPopulationTree() does not support more than 2 populations",
                path);
    }
    this->root_->set_label("root-" + this->root_->get_child(0)->get_label());
}

ComparisonRelativeRootPopulationTree::ComparisonRelativeRootPopulationTree(
        const RelativeRootComparisonSettings& settings,
        RandomNumberGenerator& rng,
        bool strict_on_constant_sites,
        bool strict_on_missing_sites,
        bool strict_on_triallelic_sites,
        bool store_seq_loci_info) {
    this->comparison_init(settings.get_path(),
               settings.get_population_name_delimiter(),
               settings.population_name_is_prefix(),
               settings.genotypes_are_diploid(),
               settings.markers_are_dominant(),
               settings.constant_sites_removed(),
               true, // validate
               strict_on_constant_sites,
               strict_on_missing_sites,
               strict_on_triallelic_sites,
               settings.get_ploidy(),
               store_seq_loci_info);
    if (settings.constrain_state_frequencies()) {
        this->constrain_state_frequencies();
        this->fold_patterns();
    }

    PositiveRealParameter p = PositiveRealParameter(
            settings.get_population_size_settings(),
            rng);
    this->set_all_population_sizes(p.get_value());

    this->relative_root_population_size_ = std::make_shared<PositiveRealParameter>(
            settings.get_relative_root_population_size_settings(),
            rng);
    this->update_root_population_size();

    this->set_population_size_prior(
            settings.get_population_size_settings().get_prior_settings().get_instance());

    if (settings.constrain_population_sizes()) {
        this->constrain_population_sizes();
    }

    if (settings.get_population_size_settings().is_fixed()) {
        this->fix_population_sizes();
    }
    
    this->set_freq_1_prior(settings.get_freq_1_settings().get_prior_settings().get_instance());
    if (settings.constrain_state_frequencies()) {
        this->constrain_state_frequencies();
    }
    else {
        PositiveRealParameter freq_1 = PositiveRealParameter(
                settings.get_freq_1_settings(),
                rng);
        this->set_freq_1(freq_1.get_value());
        if (freq_1.is_fixed()) {
            this->fix_state_frequencies();
        }
    }
    this->set_mutation_rate_parameter(
            std::make_shared<PositiveRealParameter>(
                    settings.get_mutation_rate_settings(),
                    rng));
    if (
        (this->data_.get_number_of_populations() == 1) &&
        (
            (this->population_sizes_are_constrained()) ||
            (
                settings.get_population_size_settings().is_fixed() &&
                this->relative_root_population_size_is_fixed() &&
                (this->get_root_population_size() == this->get_child_population_size(0))
            )
        )
    ) {
        std::ostringstream message;
        message << "\n#######################################################################\n"
                <<   "###############################  ERROR  ###############################\n"
                << "The alignment in:\n    \'"
                << this->data_.get_path() << "\'\n"
                << "contains only a single population, but you have fixed and/or\n"
                << "constrained the population sizes for this comparison. The timing of\n"
                << "population expansion/contraction cannot be estimated if the ancestral\n"
                << "and descendant population sizes for this comparison are either\n"
                << "fixed or constrained to be equal. Please update your configuration\n"
                << "file to estimate the unconstrained population sizes for this\n"
                << "comparison and re-run the analysis.\n"
                << "#######################################################################\n";
        throw EcoevolityComparisonSettingError(message.str(), this->data_.get_path());
    }
}

void ComparisonRelativeRootPopulationTree::set_child_population_size(
        unsigned int child_index,
        double size) {
    if (this->root_->get_child(child_index)->population_size_is_fixed()) {
        return;
    }
    this->root_->get_child(child_index)->set_population_size(size);
    ECOEVOLITY_ASSERT(! this->root_->population_size_is_fixed());
    this->update_root_population_size();
}

double ComparisonRelativeRootPopulationTree::get_child_population_size(
        unsigned int child_index) const {
    return this->root_->get_child(child_index)->get_population_size();
}

std::shared_ptr<PositiveRealParameter> ComparisonRelativeRootPopulationTree::get_child_population_size_parameter(
        unsigned int child_index) const {
    return this->root_->get_child(child_index)->get_population_size_parameter();
}

// Node height sharing needs to be dealt with in next level up in
// class hierarchy (ComparisonRelativeRootPopulationTreeCollection)
double ComparisonRelativeRootPopulationTree::compute_log_prior_density() {
    double d = 0.0;
    d += this->compute_log_prior_density_of_state_frequencies();
    d += this->compute_log_prior_density_of_mutation_rate();
    // d += this->compute_log_prior_density_of_node_heights();
    d += this->compute_log_prior_density_of_population_sizes();
    this->log_prior_density_.set_value(d);
    return d;
}

// Node height (re)storing is managed by ComparisonRelativeRootPopulationTree.
void ComparisonRelativeRootPopulationTree::store_parameters() {
    this->store_freq_1();
    this->store_mutation_rate();
    this->store_all_population_sizes();
    // this->store_all_heights();
}
void ComparisonRelativeRootPopulationTree::restore_parameters() {
    this->restore_freq_1();
    this->restore_mutation_rate();
    this->restore_all_population_sizes();
    // this->restore_all_heights();
}

void ComparisonRelativeRootPopulationTree::write_state_log_header(
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
        << "mutation_rate" << suffix << delimiter
        << "freq_1" << suffix << delimiter
        << "pop_size" << suffix << delimiter;
    if (this->root_->get_number_of_children() > 1) {
        out << "pop_size" << "_" << this->root_->get_child(1)->get_label() << delimiter;
    }
    out << "pop_size_root" << suffix;
}

void ComparisonRelativeRootPopulationTree::log_state(
        std::ostream& out,
        const std::string& delimiter) const {
    out << this->log_likelihood_.get_value() << delimiter
        << this->log_prior_density_.get_value() << delimiter
        << this->get_root_height() << delimiter
        << this->get_mutation_rate() << delimiter
        << this->get_freq_1() << delimiter
        << this->get_child_population_size(0) << delimiter;
    if (this->root_->get_number_of_children() > 1) {
        out << this->get_child_population_size(1) << delimiter;
    }
    out << this->get_root_population_size();
}
void ComparisonRelativeRootPopulationTree::log_state(
        std::ostream& out,
        unsigned int event_index,
        const std::string& delimiter) const {
    out << event_index << delimiter;
    this->log_state(out, delimiter);
}

void ComparisonRelativeRootPopulationTree::draw_from_prior(RandomNumberGenerator& rng) {
    if ((! this->state_frequencies_are_fixed()) && (! this->state_frequencies_are_constrained())) {
        this->freq_1_->set_value_from_prior(rng);
    }
    if (! this->mutation_rate_is_fixed()) {
        this->mutation_rate_->set_value_from_prior(rng);
    }
    if (! this->population_sizes_are_fixed()) {
        if (this->population_sizes_are_constrained()) {
            this->set_root_population_size(this->population_size_prior_->draw(rng));
        }
        else {
            if (! this->root_->get_child(0)->population_size_is_fixed()) {
                this->set_child_population_size(0, this->population_size_prior_->draw(rng));
            }
            if (this->root_->get_number_of_children() > 1) {
                if (! this->root_->get_child(1)->population_size_is_fixed()) {
                    this->set_child_population_size(1, this->population_size_prior_->draw(rng));
                }
            }
            if (! this->relative_root_population_size_is_fixed()) {
                this->relative_root_population_size_->set_value_from_prior(rng);
                this->update_root_population_size();
            }
        }
    }
}
