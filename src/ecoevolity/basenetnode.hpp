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

#ifndef ECOEVOLITY_BASENETNODE_HPP
#define ECOEVOLITY_BASENETNODE_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include <stdexcept>
#include <queue>

#include "matrix.hpp"
#include "parameter.hpp"
#include "probability.hpp"
#include "debug.hpp"
#include "assert.hpp"
#include "error.hpp"
#include "split.hpp"
#include "rng.hpp"

/**
 * Base class for a node of a phylogenetic network.
 *
 * @note    Many of the class' methods modified from:
 *              BasicTNode class of Bio++ Library
 *              <http://biopp.univ-montp2.fr/wiki/index.php/Main_Page>
 *              License:    CeCILL <http://www.cecill.info>
 *              Author:     Sylvain Gaillard
 *              Copyright:  CNRS, (January 12, 2011)
 */
template<class DerivedNodeT>
class BaseNetNode : public std::enable_shared_from_this<DerivedNodeT> {
    protected:
        std::vector< std::shared_ptr<DerivedNodeT> > children_;
        std::vector<std::weak_ptr<DerivedNodeT> > parents_;
        std::string label_ = "";
        int index_ = -1;
        std::shared_ptr<PositiveRealParameter> height_ = std::make_shared<PositiveRealParameter>(0.0);
        std::shared_ptr<PositiveRealParameter> stored_height_ = std::make_shared<PositiveRealParameter>(0.0);
        bool is_dirty_ = true;

        void add_ln_relative_node_height_prior_density(
                double& density,
                std::vector< std::shared_ptr<PositiveRealParameter> >& parameters,
                std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) const {
            this->visit(visited_nodes);
            bool parameter_found = false;
            for (auto parameter_iter : parameters) {
                if (parameter_iter == this->height_) {
                    parameter_found = true;
                    break;
                }
            }
            if (! parameter_found) {
                density += this->height_->relative_prior_ln_pdf();
                parameters.push_back(this->height_);
            }
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->add_ln_relative_node_height_prior_density(density, parameters, visited_nodes);
                }
            }
        }

    public:
        Split split_;

        // Constructors
        BaseNetNode() { }
        BaseNetNode(int index) {
            this->index_ = index;
        }
        BaseNetNode(std::string label) {
            this->label_ = label;
        }
        BaseNetNode(int index, std::string label) {
            this->index_ = index;
            this->label_ = label;
        }
        BaseNetNode(double height) {
            this->height_->set_value(height);
        }
        BaseNetNode(int index, double height) {
            this->index_ = index;
            this->height_->set_value(height);
        }
        BaseNetNode(std::string label, double height) {
            this->label_ = label;
            this->height_->set_value(height);
        }
        BaseNetNode(int index, std::string label, double height) {
            this->index_ = index;
            this->label_ = label;
            this->height_->set_value(height);
        }
        BaseNetNode(std::shared_ptr<PositiveRealParameter> height) {
            this->height_ = height;
        }
        BaseNetNode(int index, std::shared_ptr<PositiveRealParameter> height) {
            this->index_ = index;
            this->height_ = height;
        }
        BaseNetNode(std::string label, std::shared_ptr<PositiveRealParameter> height) {
            this->label_ = label;
            this->height_ = height;
        }
        BaseNetNode(int index, std::string label, std::shared_ptr<PositiveRealParameter> height) {
            this->index_ = index;
            this->label_ = label;
            this->height_ = height;
        }

        // Destructor
        // virtual ~BaseNetNode() { }

        /* DerivedNodeT* clone() const { */
        /*     return new DerivedNodeT(static_cast<DerivedNodeT const &>(* this)); */
        /* } */

        bool operator< (const DerivedNodeT & other) const {
            return this->get_height() < other.get_height();
        }

        std::shared_ptr<DerivedNodeT> get_copy() const {
            std::shared_ptr<DerivedNodeT> c = std::make_shared<DerivedNodeT>();
            c->label_ = this->label_;
            c->index_ = this->index_;
            c->split_ = this->split_;
            // Keep height parameter copies "shallow"
            c->height_ = this->height_;
            c->stored_height_ = this->stored_height_;
            c->is_dirty_ = this->is_dirty_;
            this->copy_node_type_specific_members(c);
            // Copy from this node to leaves; do not copy parent
            // if (this->has_parent()) {
            //     c->add_parent(this->get_parent()->get_copy());
            // }
            if (this->has_children()) {
                for (auto child : this->children_) {
                    c->add_child(child->get_copy());
                }
            }
            return c;
        }
        virtual void copy_node_type_specific_members(std::shared_ptr<DerivedNodeT> other) const { }

        std::shared_ptr<DerivedNodeT> get_deep_copy() const {
            std::shared_ptr<DerivedNodeT> c = std::make_shared<DerivedNodeT>();
            c->label_ = this->label_;
            c->index_ = this->index_;
            c->split_ = this->split_;
            c->height_ = std::make_shared<PositiveRealParameter>(*this->height_);
            c->stored_height_ = std::make_shared<PositiveRealParameter>(*this->stored_height_);
            c->is_dirty_ = this->is_dirty_;
            this->deep_copy_node_type_specific_members(c);
            // Copy from this node to leaves; do not copy parent
            if (this->has_children()) {
                for (auto child : this->children_) {
                    c->add_child(child->get_deep_copy());
                }
            }
            return c;
        }
        virtual void deep_copy_node_type_specific_members(std::shared_ptr<DerivedNodeT> other) const { }

        //Methods

        void visit(
                std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) const {
            if (this->has_multiple_parents()) {
                ECOEVOLITY_ASSERT(visited_nodes.count(this->shared_from_this()) < 1);
                visited_nodes.insert(this->shared_from_this());
            }
        }

        virtual void get_parameter_map(
                std::map<std::string, double> & parameter_map) const {
            parameter_map["height"] = this->get_height();
            parameter_map["length"] = this->get_length();
        }

        /**
         * Method to populate non-height related data (e.g., pop size) from
         * node comments in the form of a map.  Nothing to do for BaseNetNode, but
         * this is intended for derived classes to override
         */
        virtual void extract_data_from_node_comments(
                const std::map<std::string, std::string> & comment_map) { }

        unsigned int degree() const {
            unsigned int d = this->children_.size() + this->parents_.size();
            return d;
        }

        void resize_splits(unsigned int number_of_leaves) {
            this->split_.resize(number_of_leaves);
            for (auto child_iter: this->children_) {
                child_iter->resize_splits(number_of_leaves);
            }
        }

        unsigned int get_split_size() const {
            return this->split_.size();
        }

        // bool has_parent() const { return this->parent_.expired() ? false : true; }
        bool has_parent() const { return !this->parents_.empty(); }
        bool has_multiple_parents() const {
            return (this->parents_.size() > 1);
        }
        // bool is_root() const { return this->parent_.expired() ? true : false; }
        bool is_root() const { return this->parents_.empty(); }

        // unsigned int get_number_of_parents() const { return this->parent_.expired() ? 0 : 1; }
        unsigned int get_number_of_parents() const { return this->parents_.size(); }

        const std::shared_ptr<DerivedNodeT>& get_parent(const unsigned int index) const {
            if (index >= this->parents_.size()) {
                throw std::out_of_range("BaseNetNode::get_parent() index out of range");
            }
            return this->parents_.at(index).lock();
        }

        std::shared_ptr<DerivedNodeT> get_parent(const unsigned int index) {
            if (index >= this->parents_.size()) {
                throw std::out_of_range("BaseNetNode::get_parent() index out of range");
            }
            return this->parents_.at(index).lock();
        }
        const std::shared_ptr<DerivedNodeT>& get_parent() const {
            return this->get_parent(0);
        }
        std::shared_ptr<DerivedNodeT> get_parent() {
            return this->get_parent(0);
        }

        bool is_parent(const std::shared_ptr<DerivedNodeT>& node) const {
            for (unsigned int i = 0; i < this->parents_.size(); ++i) {
                if (this->parents_.at(i).lock() == node) {
                    return true;
                }
            }
            return false;
        }

        bool is_ancestor(const std::shared_ptr<DerivedNodeT>& node) const {
            if (this->is_root()) {
                return false;
            }
            if (this->is_parent(node)) {
                return true;
            }
            for (unsigned int i = 0; i < this->parents_.size(); ++i) {
                if (this->parents_.at(i).lock()->is_ancestor(node)) {
                    return true;
                }
            }
            return false;
        }

        void add_parent(std::shared_ptr<DerivedNodeT> node) {
            if (!node) {
                throw EcoevolityNullPointerError("BaseNetNode::add_parent(), empty node given");
            }
            if (! this->is_parent(node)) {
                this->parents_.push_back(node);
            }
            if (! node->is_child(this->shared_from_this())) {
                node->add_child(this->shared_from_this());
            }
            node->make_dirty();
        }

        std::shared_ptr<DerivedNodeT> remove_parent() {
            ECOEVOLITY_ASSERT(this->parents_.size() < 2);
            if (this->has_parent()) {
                return this->remove_parent(0);
            }
            return nullptr;
        }
        unsigned int remove_parent(std::shared_ptr<DerivedNodeT> node) {
            if (!node) {
                throw EcoevolityNullPointerError("BaseNetNode::remove_parent(), empty node given");
            }
            for (unsigned int i = 0; i < this->parents_.size(); ++i) {
                if (this->parents_.at(i).lock() == node) {
                    this->parents_.at(i).reset();
                    this->parents_.erase(this->parents_.begin() + i);
                    node->remove_child(this->shared_from_this());
                    node->make_dirty();
                    this->make_dirty();
                    return i;
                }
            }
            throw EcoevolityError("BaseNetNode::remove_parent(node); Node is not parent!");
        }
        std::shared_ptr<DerivedNodeT> remove_parent(const unsigned int index) {
            if (index >= this->parents_.size()) {
                throw std::out_of_range("BaseNetNode::remove_parent() index out of range");
            }
            std::shared_ptr<DerivedNodeT> p = this->parents_.at(index).lock();
            this->parents_.at(index).reset();
            this->parents_.erase(this->parents_.begin() + index);
            p->remove_child(this->shared_from_this());
            p->make_dirty();
            this->make_dirty();
            return p;
        }

        unsigned int collapse() {
            if (! this->has_parent()) {
                throw EcoevolityError("BaseNetNode::collapse(), node has no parent");
            }
            if (this->parents_.size() > 1) {
                throw EcoevolityError("BaseNetNode::collapse(), node has multiple parents");
            }
            std::shared_ptr<DerivedNodeT> grand_parent = this->get_parent();
            while (this->has_children()) {
                std::shared_ptr<DerivedNodeT> child = this->get_child(0);
                child->remove_parent(this->shared_from_this());
                child->add_parent(grand_parent);
            }
            while (this->has_parent()) {
                this->remove_parent(0);
            }
            grand_parent->make_dirty();
            return grand_parent->get_number_of_children();
        }

        std::shared_ptr<DerivedNodeT> split_children_from_polytomy(
                RandomNumberGenerator & rng,
                std::vector< std::shared_ptr<DerivedNodeT> >& children_to_split,
                std::shared_ptr<PositiveRealParameter> new_height_parameter,
                unsigned int number_of_leaves_in_tree = 0) {
            ECOEVOLITY_ASSERT(children_to_split.size() > 1);
            ECOEVOLITY_ASSERT(this->is_polytomy());
            std::shared_ptr<DerivedNodeT> new_node = std::make_shared<DerivedNodeT>(new_height_parameter);
            new_node->make_dirty();
            this->make_dirty();
            if (number_of_leaves_in_tree > 0) {
                new_node->split_.resize(number_of_leaves_in_tree);
            }
            for (auto child : children_to_split) {
                child->remove_parent(this->shared_from_this());
                child->add_parent(new_node);
            }
            // If all children were assigned to the new node, we need to remove
            // the old polytomy node that now has no children
            if (! this->has_children()) {
                while (this->has_parent()) {
                    this->get_parent()->add_child(new_node);
                    this->remove_parent(0);
                }
            }
            else {
                this->add_child(new_node);
            }
            new_node->finish_initializing_inserted_internal_node(rng);
            return new_node;
        }

        virtual void finish_initializing_inserted_internal_node(
                RandomNumberGenerator & rng) { }

        virtual double get_ln_prob_of_drawing_state() { return 0.0; }

        bool has_children() const { return !this->children_.empty(); }
        bool is_leaf() const { return this->children_.empty(); }

        unsigned int get_number_of_children() const { return this->children_.size(); }
        const std::shared_ptr<DerivedNodeT>& get_child(unsigned int index) const {
            if (index >= this->children_.size()) {
                throw std::out_of_range("BaseNetNode::get_child() index out of range");
            }
            return this->children_.at(index);
        }

        std::shared_ptr<DerivedNodeT> get_child(unsigned int index) {
            if (index >= this->children_.size()) {
                throw std::out_of_range("BaseNetNode::get_child() index out of range");
            }
            return this->children_.at(index);
        }

        std::shared_ptr<DerivedNodeT> get_oldest_child() {
            if (this->is_leaf()) {
                throw EcoevolityError("called BaseNetNode::get_oldest_child() on a leaf");
            }
            std::shared_ptr<DerivedNodeT> oldest_child = this->children_.at(0);
            for (unsigned int i = 1; i < this->children_.size(); ++i) {
                if (this->children_.at(i)->get_height() > oldest_child->get_height()) {
                    oldest_child = this->children_.at(i);
                }
            }
            return oldest_child;
        }

        std::vector< std::shared_ptr<DerivedNodeT> > get_children(
                std::vector<unsigned int> indices) {
            std::vector< std::shared_ptr<DerivedNodeT> > children_ptrs(indices.size());
            for (unsigned int i = 0; i < indices.size(); ++i) {
                children_ptrs.at(i) = this->get_child(indices.at(i));
            }
            return children_ptrs;
        }

        std::vector< std::shared_ptr<DerivedNodeT> > get_leaf_children() {
            std::vector< std::shared_ptr<DerivedNodeT> > children_ptrs;
            for (auto child : this->children_) {
                if (child->is_leaf()) {
                    children_ptrs.push_back(child);
                }
            }
            return children_ptrs;
        }

        std::vector< std::shared_ptr<DerivedNodeT> > get_internal_children() {
            std::vector< std::shared_ptr<DerivedNodeT> > children_ptrs;
            for (auto child : this->children_) {
                if (! child->is_leaf()) {
                    children_ptrs.push_back(child);
                }
            }
            return children_ptrs;
        }

        const std::vector< std::shared_ptr<DerivedNodeT> >& get_all_children() const {
            return this->children_;
        }

        bool is_child(const std::shared_ptr<DerivedNodeT>& node) const {
            for (auto child_iter: this->children_) {
                if (child_iter == node) {
                    return true;
                }
            }
            return false;
        }
        bool is_child(const std::string& label) const {
            for (auto child_iter: this->children_) {
                if (child_iter->get_label() == label) {
                    return true;
                }
            }
            return false;
        }

        void add_child(std::shared_ptr<DerivedNodeT> node) {
            if (!node) {
                throw EcoevolityNullPointerError("BaseNetNode::add_child(), empty node given");
            }
            if (! this->is_child(node)) {
                this->children_.push_back(node);
            }
            if (! node->is_parent(this->shared_from_this())) {
                node->add_parent(this->shared_from_this());
            }
            this->make_dirty();
        }

        void remove_child(std::shared_ptr<DerivedNodeT> node) {
            if (!node) {
                throw EcoevolityNullPointerError("BaseNetNode::remove_child(), empty node given");
            }
            for (unsigned int i = 0; i < this->children_.size(); ++i) {
                if (this->children_.at(i) == node) {
                    this->children_.at(i).reset();
                    this->children_.erase(this->children_.begin() + i);
                    node->remove_parent(this->shared_from_this());
                }
            }
            this->make_dirty();
        }

        std::shared_ptr<DerivedNodeT> remove_child(unsigned int index) {
            if (index >= this->children_.size()) {
                throw std::out_of_range("BaseNetNode::remove_child() index out of range");
            }
            std::shared_ptr<DerivedNodeT> c = this->children_.at(index);
            this->children_.at(index).reset();
            this->children_.erase(this->children_.begin() + index);
            c->remove_parent(this->shared_from_this());
            this->make_dirty();
            return c;
        }

        bool is_polytomy() const {
            return (this->get_number_of_children() > 2);
        }

        double get_height() const {
            return this->height_->get_value();
        }

        void set_height(double height) {
            this->height_->set_value(height);
            this->make_all_dirty();
        }

        void set_height_parameter(std::shared_ptr<PositiveRealParameter> height_parameter) {
            this->height_ = height_parameter;
            this->make_all_dirty();
        }

        void scale(const double multiplier, std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) {
            this->set_height(this->get_height() * multiplier);
            this->visit(visited_nodes);
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->_scale(multiplier, visited_nodes);
                }
            }
        }
        void scale(const double multiplier) {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->scale(multiplier, visited_nodes);
        }

        std::shared_ptr<PositiveRealParameter> get_height_parameter() const {
            return this->height_;
        }

        void update_height(double height) {
            this->height_->update_value(height);
            this->make_all_dirty();
        }

        void get_node(const std::string& label,
                std::shared_ptr<DerivedNodeT>& node_ptr,
                std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) {
            this->visit(visited_nodes);
            if (node_ptr != nullptr) {
                return;
            }
            if (this->is_label(label)) {
                node_ptr = this->shared_from_this();
                return;
            }
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->get_node(label, node_ptr, visited_nodes);
                }
            }
        }
        std::shared_ptr<DerivedNodeT> get_node(const std::string& label) {
            std::shared_ptr<DerivedNodeT> node_ptr = nullptr;
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->get_node(label, node_ptr, visited_nodes);
            return node_ptr;
        }

        void store_height() {
            this->height_->store();
        }
        void store_all_heights(std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) {
            this->height_->store();
            this->visit(visited_nodes);
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->store_all_heights(visited_nodes);
                }
            }
        }
        void store_all_heights() {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->store_all_heights(visited_nodes);
        }
        virtual void store_all_parameter_values() {
            this->store_all_heights();
        }

        void restore_height() {
            this->height_->restore();
            this->make_all_dirty();
        }
        void restore_all_heights(std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) {
            this->height_->restore();
            this->make_dirty();
            this->visit(visited_nodes);
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->restore_all_heights(visited_nodes);
                }
            }
        }
        void restore_all_heights() {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->restore_all_heights(visited_nodes);
        }
        virtual void restore_all_parameter_values() {
            this->restore_all_heights();
        }

        void store_height_pointer() {
            this->stored_height_ = this->height_;
        }
        void store_all_height_pointers(std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) {
            this->stored_height_ = this->height_;
            this->visit(visited_nodes);
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->store_all_height_pointers(visited_nodes);
                }
            }
        }
        void store_all_height_pointers() {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->store_all_height_pointers(visited_nodes);
        }
        virtual void store_all_parameter_pointers() {
            this->store_all_height_pointers();
        }

        void restore_height_pointer() {
            this->height_ = this->stored_height_;
            this->make_all_dirty();
        }
        void restore_all_height_pointers(std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) {
            this->height_ = this->stored_height_;
            this->visit(visited_nodes);
            this->make_dirty();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->restore_all_height_pointers(visited_nodes);
                }
            }
        }
        void restore_all_height_pointers() {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->restore_all_height_pointers(visited_nodes);
        }
        virtual void restore_all_parameter_pointers() {
            this->restore_all_height_pointers();
        }

        double get_length(const unsigned int index) const {
            return this->parents_.at(index).lock()->get_height() - this->get_height();
        }
        double get_length(const std::shared_ptr<DerivedNodeT>& node) const {
            if (! this->is_parent(node)) {
                throw EcoevolityError("BaseNetNode::get_length(), passed node is not a parent");
            }
            return node->get_height() - this->get_height();
        }
        double get_length() const {
            if (! this->has_parent()) {
                return 0.0;
            }
            return this->get_length(0);
        }

        void sum_clade_length(double & length,
                std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) const  {
            length += this->get_length();
            this->visit(visited_nodes);
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->sum_clade_length(length, visited_nodes);
                }
            }
        }
        void sum_clade_length(double & length) const  {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->sum_clade_length(length, visited_nodes);
        }

        double get_clade_length() const {
            double l = 0.0;
            this->sum_clade_length(l);
            return l;
        }

        bool node_height_is_valid() const {
            double height = this->get_height();
            if (this->has_parent()) {
                for (unsigned int i = 0; i < this->parents_.size(); ++i) {
                    if (this->parents_.at(i).lock()->get_height() < height) {
                        return false;
                    }
                }
            }
            if (this->has_children()) {
                for (unsigned int i = 0; i < this->get_number_of_children(); ++i) {
                    if (this->get_child(i)->get_height() > height) {
                        return false;
                    }
                }
            }
            return true;
        }

        bool node_heights_are_valid(
                std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) const  {
            this->visit(visited_nodes);
            if (! this->node_height_is_valid()) {
                return false;
            }
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    if (! child_iter->node_heights_are_valid(visited_nodes)) {
                        return false;
                    }
                }
            }
            return true;
        }
        bool node_heights_are_valid() const {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            return this->node_heights_are_valid(visited_nodes);
        }

        virtual int get_index() const {
            ECOEVOLITY_ASSERT(this->index_ > -1);
            return this->index_;
        }
        void set_index(int i) {
            this->index_ = i;
        }
        const std::string& get_label() const {
            return this->label_;
        }
        void set_label(std::string label) {
            this->label_ = label;
        }
        bool is_label(const std::string& label) {
            return (this->label_ == label);
        }

        bool is_dirty() const {
            return this->is_dirty_;
        }

        bool clade_has_dirt(std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) const {
            this->visit(visited_nodes);
            if (this->is_dirty()) {
                return true;
            }
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    if (child_iter->clade_has_dirt(visited_nodes)) {
                        return true;
                    }
                }
            }
            return false;
        }
        bool clade_has_dirt() const {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            return this->clade_has_dirt(visited_nodes);
        }

        void make_dirty() {
            this->is_dirty_ = true;
        }
        void make_all_dirty(std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) {
            this->is_dirty_ = true;
            this->visit(visited_nodes);
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->make_all_dirty(visited_nodes);
                }
            }
        }
        void make_all_dirty() {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->make_all_dirty(visited_nodes);
        }
        void make_clean() {
            this->is_dirty_ = false;
        }
        void make_all_clean(std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) {
            this->is_dirty_ = false;
            this->visit(visited_nodes);
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->make_all_clean(visited_nodes);
                }
            }
        }
        void make_all_clean() {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->make_all_clean(visited_nodes);
        }

        unsigned int get_node_count(std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) const {
            unsigned int n = 1;
            this->visit(visited_nodes);
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    n += child_iter->get_node_count(visited_nodes);
                }
            }
            return n;
        }
        unsigned int get_node_count() const {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            return this->get_node_count(visited_nodes);
        }
        unsigned int get_leaf_node_count(std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) const {
            this->visit(visited_nodes);
            if (this->is_leaf()) {
                return 1;
            }
            unsigned int n = 0;
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    n += child_iter->get_leaf_node_count(visited_nodes);
                }
            }
            return n;
        }
        unsigned int get_leaf_node_count() const {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            return this->get_leaf_node_count(visited_nodes);
        }
        unsigned int get_internal_node_count(std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) const {
            this->visit(visited_nodes);
            if (this->is_leaf()) {
                return 0;
            }
            unsigned int n = 1;
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    n += child_iter->get_internal_node_count(visited_nodes);
                }
            }
            return n;
        }
        unsigned int get_internal_node_count() const {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            return this->get_internal_node_count(visited_nodes);
        }
        unsigned int get_polytomy_node_count(std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) const {
            this->visit(visited_nodes);
            if (this->is_polytomy()) {
                return 1;
            }
            unsigned int n = 0;
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    n += child_iter->get_polytomy_node_count(visited_nodes);
                }
            }
            return n;
        }
        unsigned int get_polytomy_node_count() const {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            return this->get_polytomy_node_count(visited_nodes);
        }
        unsigned int get_mapped_node_count(
                std::shared_ptr<PositiveRealParameter> node_height_pointer,
                std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) const {
            this->visit(visited_nodes);
            if (this->get_height_parameter() == node_height_pointer) {
                return 1;
            }
            unsigned int n = 0;
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    n += child_iter->get_mapped_node_count(node_height_pointer, visited_nodes);
                }
            }
            return n;
        }
        unsigned int get_mapped_node_count(
                std::shared_ptr<PositiveRealParameter> node_height_pointer) const {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            return this->get_mapped_node_count(node_height_pointer, visited_nodes);
        }
        unsigned int get_mapped_polytomy_node_count(
                std::shared_ptr<PositiveRealParameter> node_height_pointer,
                std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) const {
            this->visit(visited_nodes);
            if ((this->get_height_parameter() == node_height_pointer) &&
                    (this->is_polytomy())) {
                return 1;
            }
            unsigned int n = 0;
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    n += child_iter->get_mapped_polytomy_node_count(node_height_pointer, visited_nodes);
                }
            }
            return n;
        }
        unsigned int get_mapped_polytomy_node_count(
                std::shared_ptr<PositiveRealParameter> node_height_pointer) const {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            return this->get_mapped_polytomy_node_count(node_height_pointer, visited_nodes);
        }

        void get_nodes(std::vector< std::shared_ptr<DerivedNodeT> >& nodes,
                std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) {
            this->visit(visited_nodes);
            nodes.push_back(this->shared_from_this());
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->get_nodes(nodes, visited_nodes);
                }
            }
        }
        void get_nodes(std::vector< std::shared_ptr<DerivedNodeT> >& nodes) {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->get_nodes(nodes, visited_nodes);
        }
        std::vector< std::shared_ptr<DerivedNodeT> > get_nodes() {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            std::vector< std::shared_ptr<DerivedNodeT> > nodes;
            this->get_nodes(nodes, visited_nodes);
            return nodes;
        }

        void pre_order(std::vector< std::shared_ptr<DerivedNodeT> >& nodes) {
            nodes.clear();
            this->get_nodes(nodes);
        }

        void level_order(std::vector< std::shared_ptr<DerivedNodeT> >& nodes) {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            nodes.clear();
            std::queue< std::shared_ptr<DerivedNodeT> > q;
            std::shared_ptr<DerivedNodeT> nd = this->shared_from_this();
            q.push(nd);
            while(! q.empty()) {
                nd = q.front(); q.pop();
                nodes.push_back(nd);
                this->visit(visited_nodes);
                if (! nd->is_leaf()) {
                    for (auto child : nd->children_) {
                        if (visited_nodes.count(child) < 1) {
                            q.push(child);
                        }
                    }
                }
            }
        }

        void get_internal_nodes(std::vector< std::shared_ptr<DerivedNodeT> >& internal_nodes,
                std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) {
            this->visit(visited_nodes);
            if (! this->is_leaf()) {
                internal_nodes.push_back(this->shared_from_this());
            }
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->get_internal_nodes(internal_nodes, visited_nodes);
                }
            }
        }
        void get_internal_nodes(std::vector< std::shared_ptr<DerivedNodeT> >& internal_nodes) {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->get_internal_nodes(internal_nodes, visited_nodes);
        }
        std::vector< std::shared_ptr<DerivedNodeT> > get_internal_nodes() {
            std::vector< std::shared_ptr<DerivedNodeT> > internal_nodes;
            internal_nodes.reserve(this->get_internal_node_count());
            this->get_internal_nodes(internal_nodes);
            return internal_nodes;
        }

        void get_polytomy_nodes(std::vector< std::shared_ptr<DerivedNodeT> >& polytomy_nodes,
                std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) {
            this->visit(visited_nodes);
            if (this->is_polytomy()) {
                polytomy_nodes.push_back(this->shared_from_this());
            }
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->get_polytomy_nodes(polytomy_nodes, visited_nodes);
                }
            }
        }
        void get_polytomy_nodes(std::vector< std::shared_ptr<DerivedNodeT> >& polytomy_nodes) {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->get_polytomy_nodes(polytomy_nodes, visited_nodes);
        }
        std::vector< std::shared_ptr<DerivedNodeT> > get_polytomy_nodes() {
            std::vector< std::shared_ptr<DerivedNodeT> > nodes;
            this->get_polytomy_nodes(nodes);
            return nodes;
        }

        void get_mapped_nodes(std::vector< std::shared_ptr<DerivedNodeT> >& mapped_nodes,
                std::shared_ptr<PositiveRealParameter> node_height_pointer,
                std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) {
            this->visit(visited_nodes);
            if (this->get_height_parameter() == node_height_pointer) {
                mapped_nodes.push_back(this->shared_from_this());
            }
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->get_mapped_nodes(mapped_nodes, node_height_pointer, visited_nodes);
                }
            }
        }
        void get_mapped_nodes(std::vector< std::shared_ptr<DerivedNodeT> >& mapped_nodes,
                std::shared_ptr<PositiveRealParameter> node_height_pointer) {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->get_mapped_nodes(mapped_nodes, node_height_pointer, visited_nodes);
        }
        std::vector< std::shared_ptr<DerivedNodeT> > get_mapped_nodes(
                std::shared_ptr<PositiveRealParameter> node_height_pointer) {
            std::vector< std::shared_ptr<DerivedNodeT> > mapped_nodes;
            this->get_mapped_nodes(mapped_nodes, node_height_pointer);
            return mapped_nodes;
        }

        void get_mapped_polytomy_nodes(std::vector< std::shared_ptr<DerivedNodeT> >& mapped_nodes,
                std::shared_ptr<PositiveRealParameter> node_height_pointer,
                std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) {
            this->visit(visited_nodes);
            if ((this->get_height_parameter() == node_height_pointer) &&
                    (this->is_polytomy())) {
                mapped_nodes.push_back(this->shared_from_this());
            }
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->get_mapped_polytomy_nodes(mapped_nodes, node_height_pointer, visited_nodes);
                }
            }
        }
        void get_mapped_polytomy_nodes(std::vector< std::shared_ptr<DerivedNodeT> >& mapped_nodes,
                std::shared_ptr<PositiveRealParameter> node_height_pointer) {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->get_mapped_polytomy_nodes(mapped_nodes, node_height_pointer, visited_nodes);
        }
        std::vector< std::shared_ptr<DerivedNodeT> > get_mapped_polytomy_nodes(
                std::shared_ptr<PositiveRealParameter> node_height_pointer) {
            std::vector< std::shared_ptr<DerivedNodeT> > mapped_nodes;
            this->get_mapped_polytomy_nodes(mapped_nodes, node_height_pointer);
            return mapped_nodes;
        }

        void get_leaves(std::vector< std::shared_ptr<DerivedNodeT> >& leaves,
                std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) {
            this->visit(visited_nodes);
            if (this->is_leaf()) {
                leaves.push_back(this->shared_from_this());
            }
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->get_leaves(leaves, visited_nodes);
                }
            }
        }
        void get_leaves(std::vector< std::shared_ptr<DerivedNodeT> >& leaves) {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->get_leaves(leaves, visited_nodes);
        }
        std::vector< std::shared_ptr<DerivedNodeT> > get_leaves() {
            std::vector< std::shared_ptr<DerivedNodeT> > leaves;
            leaves.reserve(this->get_leaf_node_count());
            this->get_leaves(leaves);
            return leaves;
        }

        void get_leaf_labels(std::vector<std::string>& leaf_labels,
                std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) const {
            this->visit(visited_nodes);
            if (this->is_leaf()) {
                leaf_labels.push_back(this->get_label());
            }
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->get_leaf_labels(leaf_labels, visited_nodes);
                }
            }
        }
        void get_leaf_labels(std::vector<std::string>& leaf_labels) const {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->get_leaf_labels(leaf_labels, visited_nodes);
        }
        std::vector<std::string> get_leaf_labels() const {
            std::vector<std::string> leaf_labels;
            leaf_labels.reserve(this->get_leaf_node_count());
            this->get_leaf_labels(leaf_labels);
            return leaf_labels;
        }

        void get_leaf_label_set(std::set<std::string>& leaf_labels,
                std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) const {
            this->visit(visited_nodes);
            if (this->is_leaf()) {
                leaf_labels.insert(this->get_label());
            }
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->get_leaf_label_set(leaf_labels, visited_nodes);
                }
            }
        }
        void get_leaf_label_set(std::set<std::string>& leaf_labels) const {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->get_leaf_label_set(leaf_labels, visited_nodes);
        }
        std::set<std::string> get_leaf_label_set() const {
            std::set<std::string> leaf_labels;
            this->get_leaf_label_set(leaf_labels);
            return leaf_labels;
        }

        std::shared_ptr<ContinuousProbabilityDistribution> get_node_height_prior() const {
            return this->height_->get_prior();
        }
        void set_node_height_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            this->height_->set_prior(prior);
            this->make_all_dirty();
        }
        void set_all_node_height_priors(std::shared_ptr<ContinuousProbabilityDistribution> prior,
                std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) {
            this->visit(visited_nodes);
            this->height_->set_prior(prior);
            this->make_dirty();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->set_all_node_height_priors(prior, visited_nodes);
                }
            }
        }
        void set_all_node_height_priors(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->set_all_node_height_priors(prior, visited_nodes);
        }

        void fix_node_height() {
            this->height_->fix();
        }
        void fix_all_node_heights(std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) {
            this->visit(visited_nodes);
            this->height_->fix();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->fix_all_node_heights(visited_nodes);
                }
            }
        }
        void fix_all_node_heights() {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->fix_all_node_heights(visited_nodes);
        }

        void estimate_node_height() {
            this->height_->estimate();
        }
        void estimate_all_node_heights(std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) {
            this->visit(visited_nodes);
            this->height_->estimate();
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    child_iter->estimate_all_node_heights(visited_nodes);
                }
            }
        }
        void estimate_all_node_heights() {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->estimate_all_node_heights(visited_nodes);
        }

        bool node_height_is_fixed() const {
            return this->height_->is_fixed();
        }

        bool all_node_heights_are_fixed(std::set< std::shared_ptr<const DerivedNodeT> > & visited_nodes) const {
            this->visit(visited_nodes);
            if (! this->node_height_is_fixed()) {
                return false;
            }
            for (auto child_iter: this->children_) {
                if (visited_nodes.count(child_iter) < 1) {
                    if (! child_iter->all_node_heights_are_fixed(visited_nodes)) {
                        return false;
                    }
                }
            }
            return true;
        }
        bool all_node_heights_are_fixed() const {
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            return this->all_node_heights_are_fixed(visited_nodes);
        }

        double get_height_relative_prior_ln_pdf() const {
            return this->height_->relative_prior_ln_pdf();
        }

        double calculate_ln_relative_node_height_prior_density() const {
            double d = 0.0;
            std::vector< std::shared_ptr<PositiveRealParameter> > parameters;
            parameters.reserve(this->get_node_count());
            std::set< std::shared_ptr<const DerivedNodeT> > visited_nodes;
            this->add_ln_relative_node_height_prior_density(d, parameters, visited_nodes);
            return d;
        }

        // TODO:: This needs to be updated to handle reticulating nodes
        std::string to_parentheses(unsigned int precision = 12,
                const bool label_internal_nodes = false) const {
            std::ostringstream s;
            s.precision(precision);
            if (this->is_leaf()) {
                s << this->get_label();
            }
            else {
                unsigned int child_idx = 0;
                s << "(";
                for (auto child_iter: this->children_) {
                    if (child_idx > 0) {
                        s << ",";
                    }
                    s << child_iter->to_parentheses(precision,
                            label_internal_nodes);
                    ++child_idx;
                }
                s << ")";
                if (label_internal_nodes) {
                    s << this->get_label();
                }
            }
            s << ":" << this->get_length();
            return s.str();
        }

        std::string get_comment_data_string(
                unsigned int precision = 12) const {
            std::ostringstream s;
            s.precision(precision);
            std::string additional_comment_str = this->get_additional_comment_data_string(precision);
            s << "height="
              << this->get_height();
            if (additional_comment_str.length() > 0) {
                s << ","
                  << additional_comment_str;
            }
            return s.str();
        }

        // Descendant classes can override this method to populate additional
        // node data (in addition to height. E.g., population size).
        virtual std::string get_additional_comment_data_string(
                unsigned int precision = 12) const { return ""; }
};

#endif
