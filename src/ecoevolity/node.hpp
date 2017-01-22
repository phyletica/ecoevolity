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

#ifndef ECOEVOLITY_NODE_HPP
#define ECOEVOLITY_NODE_HPP

#include <memory>

#include "basenode.hpp"
#include "parameter.hpp"

/**
 * Base class for a node of a phylogenetic tree.
 *
 * @note    Many of the class' methods modified from:
 *              BasicTNode class of Bio++ Library
 *              <http://biopp.univ-montp2.fr/wiki/index.php/Main_Page>
 *              License:    CeCILL <http://www.cecill.info>
 *              Author:     Sylvain Gaillard
 *              Copyright:  CNRS, (January 12, 2011)
 */
class Node: public BaseNode<Node>{
    private:
        typedef BaseNode<Node> BaseClass;

    public:
        Node() { }
        Node(std::string label) : BaseClass(label) { }
        Node(double height) : BaseClass(height) { }
        Node(std::string label, double height) : BaseClass(label, height) { }
};

class GeneTreeSimNode : public BaseNode<GeneTreeSimNode> {
    private:
        typedef BaseNode<GeneTreeSimNode> BaseClass;
        int character_state_ = -1;
        unsigned int population_index_ = 0;
        double p[2][2];

    public:
        GeneTreeSimNode() { }
        GeneTreeSimNode(std::string label) : BaseClass(label) { }
        GeneTreeSimNode(double height) : BaseClass(height) { }
        GeneTreeSimNode(unsigned int population_index) :
            BaseClass(),
            population_index_(population_index)
            { }
        GeneTreeSimNode(std::string label, double height) :
            BaseClass(label, height)
            { }
        GeneTreeSimNode(unsigned int population_index, double height) :
            BaseClass(height),
            population_index_(population_index)
            { }
        GeneTreeSimNode(unsigned int population_index, std::string label) :
            BaseClass(label),
            population_index_(population_index)
            { }
        GeneTreeSimNode(unsigned int population_index, std::string label,
                double height) :
            BaseClass(label, height),
            population_index_(population_index)
            { }

        int get_character_state() const {
            return this->character_state_;
        }
        void set_character_state(int state) {
            this->character_state_ = state;
        }
        unsigned int get_population_index() const {
            return this->population_index_;
        }
        void set_population_index(unsigned int index) {
            this->population_index_ = index;
        }

        void compute_branch_binary_transition_probabilities(
                const double u,
                const double v) {
            double freq_0 = u / (u + v);
            double freq_1 = 1.0 - freq_0;
            double t = this->get_length();
            double r = 1.0 - std::exp(-(u + v) * t);
            this->p[1][0] = freq_0 * r;
            this->p[1][1] = 1.0 - (freq_0 * r);
            this->p[0][1] = freq_1 * r;
            this->p[0][0] = 1.0 - (freq_1 * r);
        }

        void compute_binary_transition_probabilities(
                const double u,
                const double v) {
            this->compute_branch_binary_transition_probabilities(u, v);
            for (auto child_iter: this->children_) {
                child_iter->compute_binary_transition_probabilities(u, v);
            }
        }

        void draw_binary_character(
                const double freq_0,
                RandomNumberGenerator& rng) {
            double u = rng.uniform_real();
            if (this->has_parent()) {
                int i = this->get_parent()->get_character_state();
                ECOEVOLITY_ASSERT((i == 0) || (i == 1));
                if (u < this->p[i][0]) {
                    this->set_character_state(0);
                }
                else {
                    this->set_character_state(1);
                }
                return;
            }
            if (u < freq_0) {
                this->set_character_state(0);
            }
            else {
                this->set_character_state(1);
            }
        }

        void simulate_binary_character(
                const double freq_0,
                RandomNumberGenerator& rng) {
            this->draw_binary_character(freq_0, rng);
            for (auto child_iter: this->children_) {
                child_iter->simulate_binary_character(freq_0, rng);
            }

        }

        void get_allele_counts(
                std::vector<unsigned int>& allele_counts,
                std::vector<unsigned int>& red_allele_counts) const {
            if (this->is_leaf()) {
                unsigned int pop_idx = this->get_population_index();
                ++allele_counts.at(pop_idx);
                if (this->get_character_state() == 1) {
                    ++red_allele_counts.at(pop_idx);
                }
            }
            for (auto child_iter: this->children_) {
                child_iter->get_allele_counts(
                        allele_counts,
                        red_allele_counts);
            }
        }

        void get_allele_counts(
                std::vector<unsigned int>& allele_counts,
                std::vector<unsigned int>& red_allele_counts,
                std::vector<int>& last_allele) const {
            if (this->is_leaf()) {
                unsigned int pop_idx = this->get_population_index();
                if (last_allele.at(pop_idx) < 0) {
                    last_allele.at(pop_idx) = this->get_character_state();
                    return;
                }
                if ((last_allele.at(pop_idx) == 1) || (this->get_character_state() == 1)) {
                    ++red_allele_counts.at(pop_idx);
                }
                ++allele_counts.at(pop_idx);
                last_allele.at(pop_idx) = -1;
                return;
            }
            for (auto child_iter: this->children_) {
                child_iter->get_allele_counts(
                        allele_counts,
                        red_allele_counts,
                        last_allele);
            }
        }
};

class PopulationNode: public BaseNode<PopulationNode>{
    private:
        typedef BaseNode<PopulationNode> BaseClass;
        int population_index_ = -1;
        BiallelicPatternProbabilityMatrix bottom_pattern_probs_;
        BiallelicPatternProbabilityMatrix top_pattern_probs_;
        std::shared_ptr<PositiveRealParameter> population_size_ = std::make_shared<PositiveRealParameter>(0.001);

        void add_ln_relative_population_size_prior_density(
                double& density,
                std::vector< std::shared_ptr<PositiveRealParameter> >& parameters) const {
            bool parameter_found = false;
            for (auto parameter_iter : parameters) {
                if (parameter_iter == this->population_size_) {
                    parameter_found = true;
                    break;
                }
            }
            if (! parameter_found) {
                density += this->population_size_->relative_prior_ln_pdf();
                parameters.push_back(this->population_size_);
            }
            for (unsigned int i = 0; i < this->children_.size(); ++i) {
                this->children_.at(i)->add_ln_relative_population_size_prior_density(density, parameters);
            }
        }

        void scale_all_population_size_parameters(
                double scale,
                std::vector< std::shared_ptr<PositiveRealParameter> >& parameters,
                unsigned int & number_of_free_parameters_scaled) {
            bool parameter_found = false;
            for (auto parameter_iter : parameters) {
                if (parameter_iter == this->population_size_) {
                    parameter_found = true;
                    break;
                }
            }
            if (! parameter_found) {
                this->set_population_size(this->get_population_size() * scale);
                parameters.push_back(this->population_size_);
                ++number_of_free_parameters_scaled;
            }
            for (unsigned int i = 0; i < this->children_.size(); ++i) {
                this->children_.at(i)->scale_all_population_size_parameters(
                        scale,
                        parameters,
                        number_of_free_parameters_scaled);
            }
        }

    public:
        PopulationNode() { }
        PopulationNode(std::string label) : BaseClass(label) { }
        PopulationNode(double height) : BaseClass(height) { }
        PopulationNode(std::string label, double height) :
            BaseClass(label, height)
            { }
        PopulationNode(unsigned int allele_count) : BaseClass() {
            this->bottom_pattern_probs_.resize(allele_count);
            this->top_pattern_probs_.resize(allele_count);
        }
        PopulationNode(std::string label, unsigned int allele_count) :
            BaseClass(label)
        {
            this->bottom_pattern_probs_.resize(allele_count);
            this->top_pattern_probs_.resize(allele_count);
        }
        PopulationNode(double height, unsigned int allele_count) :
            BaseClass(height)
        {
            this->bottom_pattern_probs_.resize(allele_count);
            this->top_pattern_probs_.resize(allele_count);
        }
        PopulationNode(std::string label,
                double height,
                unsigned int allele_count) :
            BaseClass(label, height)
        {
            this->bottom_pattern_probs_.resize(allele_count);
            this->top_pattern_probs_.resize(allele_count);
        }
        PopulationNode(
                unsigned int population_index,
                std::string label,
                double height,
                unsigned int allele_count) :
            BaseClass(label, height)
        {
            this->bottom_pattern_probs_.resize(allele_count);
            this->top_pattern_probs_.resize(allele_count);
            this->population_index_ = population_index;
        }
        PopulationNode(const PopulationNode& node) :
            BaseClass(node.label_, node.height_)
        {
            this->population_size_ = node.population_size_;
            this->bottom_pattern_probs_ = node.bottom_pattern_probs_;
            this->top_pattern_probs_ = node.top_pattern_probs_;
            this->population_index_ = node.population_index_;
        }

        // overload copy operator
        // PopulationNode& operator=(const PopulationNode& node) {
        //     this->children_ = node.children_;
        //     this->parent_ = node.parent_;
        //     this->height_->set_value(node.height_->get_value());
        //     this->label_ = node.label_;
        //     this->is_dirty_ = node.is_dirty_;
        //     this->bottom_pattern_probs_ = node.bottom_pattern_probs_;
        //     this->top_pattern_probs_ = node.top_pattern_probs_;
        //     return * this;
        // }

        int get_population_index() const {
            return this->population_index_;
        }

        // methods for accessing/changing pattern probabilities
        unsigned int get_allele_count() const {
            return this->bottom_pattern_probs_.get_allele_count();
        }

        unsigned int get_leaf_allele_count() const {
            if (this->is_leaf()) {
                return this->get_allele_count();
            }
            unsigned int n = 0;
            for (auto child_iter: this->children_) {
                n += child_iter->get_leaf_allele_count();
            }
            return n;
        }

        void resize(unsigned int allele_count) {
            this->bottom_pattern_probs_.resize(allele_count);
            this->top_pattern_probs_.resize(allele_count);
        }
        void reset(unsigned int allele_count) {
            this->bottom_pattern_probs_.reset(allele_count);
            this->top_pattern_probs_.reset(allele_count);
        }
        void resize_all() {
            this->resize(this->get_leaf_allele_count());
            for (auto child_iter: this->children_) {
                child_iter->resize_all();
            }
        }

        std::shared_ptr<PopulationNode> get_clone() const {
            return std::make_shared<PopulationNode>(*this);
        }

        std::shared_ptr<PopulationNode> get_clade_clone() const {
            if (this->is_leaf()) {
                return this->get_clone();
            }
            std::shared_ptr<PopulationNode> n = this->get_clone();
            for (auto child_iter: this->children_) {
                n->add_child(child_iter->get_clade_clone());
            }
            return n;
        }

        const BiallelicPatternProbabilityMatrix& get_bottom_pattern_probs() const{
            return this->bottom_pattern_probs_;
        }
        const BiallelicPatternProbabilityMatrix& get_top_pattern_probs() const{
            return this->top_pattern_probs_;
        }
        BiallelicPatternProbabilityMatrix* clone_bottom_pattern_probs() const{
            return this->bottom_pattern_probs_.clone();
        }
        BiallelicPatternProbabilityMatrix* clone_top_pattern_probs() const{
            return this->top_pattern_probs_.clone();
        }

        void copy_bottom_pattern_probs(const BiallelicPatternProbabilityMatrix& m) {
            // No check here; bottom probs will be updated first and can differ
            // in size until the top is also updated.
            // if (m.get_allele_count() != this->top_pattern_probs_.get_allele_count()) {
            //     throw EcoevolityError(
            //             "PopulationNode:copy_bottom_pattern_probs(); allele counts must be the same between top and bottom of branch");
            // }
            this->bottom_pattern_probs_.copy(m);
        }
        void copy_top_pattern_probs(const BiallelicPatternProbabilityMatrix& m) {
            if (m.get_allele_count() != this->bottom_pattern_probs_.get_allele_count()) {
                throw EcoevolityError(
                        "PopulationNode:copy_top_pattern_probs(); allele counts must be the same between top and bottom of branch");
            }
            this->top_pattern_probs_.copy(m);
        }
        void copy_pattern_probs(
                const BiallelicPatternProbabilityMatrix& bottom_probs,
                const BiallelicPatternProbabilityMatrix& top_probs) {
            if (bottom_probs.get_allele_count() != top_probs.get_allele_count()) {
                throw EcoevolityError(
                        "PopulationNode:copy_pattern_probs(); allele counts must be the same between top and bottom of branch");
            }
            this->bottom_pattern_probs_.copy(bottom_probs);
            this->top_pattern_probs_.copy(top_probs);
        }

        double get_bottom_pattern_probability(
                unsigned int allele_count,
                unsigned int red_allele_count) const {
            return this->bottom_pattern_probs_.get_pattern_probability(
                    allele_count,
                    red_allele_count);
        }
        double get_top_pattern_probability(
                unsigned int allele_count,
                unsigned int red_allele_count) const {
            return this->top_pattern_probs_.get_pattern_probability(
                    allele_count,
                    red_allele_count);
        }
        void set_bottom_pattern_probability(
                unsigned int allele_count,
                unsigned int red_allele_count,
                double probability) {
            this->bottom_pattern_probs_.set_pattern_probability(
                    allele_count,
                    red_allele_count,
                    probability);
        }
        void set_top_pattern_probability(
                unsigned int allele_count,
                unsigned int red_allele_count,
                double probability) {
            this->top_pattern_probs_.set_pattern_probability(
                    allele_count,
                    red_allele_count,
                    probability);
        }

        double get_population_size() const {
            return this->population_size_->get_value();
        }
        std::shared_ptr<PositiveRealParameter> get_population_size_parameter() const {
            return this->population_size_;
        }
        void set_population_size_parameter(std::shared_ptr<PositiveRealParameter> size) {
            this->population_size_ = size;
            this->make_all_dirty();
        }
        void set_all_population_size_parameters(std::shared_ptr<PositiveRealParameter> size) {
            this->population_size_ = size;
            this->make_dirty();
            for (auto child_iter: this->children_) {
                child_iter->set_all_population_size_parameters(size);
            }
        }
        void set_all_population_size_parameters() {
            std::shared_ptr<PositiveRealParameter> size = this->population_size_;
            this->make_dirty();
            for (auto child_iter: this->children_) {
                child_iter->set_all_population_size_parameters(size);
            }
        }

        void set_population_size(double size) {
            this->population_size_->set_value(size);
            this->make_all_dirty();
        }
        void set_all_population_sizes(double size) {
            this->population_size_->set_value(size);
            this->make_dirty();
            for (auto child_iter: this->children_) {
                child_iter->set_all_population_sizes(size);
            }
        }
        void update_population_size(double size) {
            this->population_size_->update_value(size);
            this->make_all_dirty();
        }
        void update_all_population_sizes(double size) {
            this->population_size_->update_value(size);
            this->make_dirty();
            for (auto child_iter: this->children_) {
                child_iter->update_all_population_sizes(size);
            }
        }

        void store_population_size() {
            this->population_size_->store();
        }
        void restore_population_size() {
            this->population_size_->restore();
            this->make_all_dirty();
        }
        void store_all_population_sizes() {
            this->population_size_->store();
            for (auto child_iter: this->children_) {
                child_iter->store_all_population_sizes();
            }
        }
        void restore_all_population_sizes() {
            this->population_size_->restore();
            this->make_dirty();
            for (auto child_iter: this->children_) {
                child_iter->restore_all_population_sizes();
            }
        }

        void set_population_size_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            this->population_size_->set_prior(prior);
            this->make_all_dirty();
        }
        void set_all_population_size_priors(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            this->population_size_->set_prior(prior);
            this->make_dirty();
            for (auto child_iter: this->children_) {
                child_iter->set_all_population_size_priors(prior);
            }
        }

        void fix_population_size() {
            this->population_size_->fix();
        }
        void fix_all_population_sizes() {
            this->population_size_->fix();
            for (auto child_iter: this->children_) {
                child_iter->fix_all_population_sizes();
            }
        }
        void estimate_population_size() {
            this->population_size_->estimate();
        }
        void estimate_all_population_sizes() {
            this->population_size_->estimate();
            for (auto child_iter: this->children_) {
                child_iter->estimate_all_population_sizes();
            }
        }

        double calculate_ln_relative_population_size_prior_density() const {
            double d = 0.0;
            std::vector< std::shared_ptr<PositiveRealParameter> > parameters;
            parameters.reserve(this->get_node_count());
            this->add_ln_relative_population_size_prior_density(d, parameters);
            return d;
        }

        unsigned int scale_all_population_size_parameters(double scale) {
            std::vector< std::shared_ptr<PositiveRealParameter> > parameters;
            parameters.reserve(this->get_node_count());
            unsigned int number_of_free_parameters_scaled = 0;
            this->scale_all_population_size_parameters(
                    scale,
                    parameters,
                    number_of_free_parameters_scaled);
            return number_of_free_parameters_scaled;
        }

        bool population_size_is_fixed() const {
            return this->population_size_->is_fixed();
        }

        bool all_population_sizes_are_fixed() const {
            if (! this->population_size_is_fixed()) {
                return false;
            }
            for (auto child_iter: this->children_) {
                if (! child_iter->all_population_sizes_are_fixed()) {
                    return false;
                }
            }
            return true;
        }

};

#endif
