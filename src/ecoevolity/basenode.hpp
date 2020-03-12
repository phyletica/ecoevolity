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

#ifndef ECOEVOLITY_BASENODE_HPP
#define ECOEVOLITY_BASENODE_HPP

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
template<class DerivedNodeT>
class BaseNode : public std::enable_shared_from_this<DerivedNodeT> {
    protected:
        std::vector< std::shared_ptr<DerivedNodeT> > children_;
        std::weak_ptr<DerivedNodeT> parent_;
        std::string label_ = "";
        int index_ = -1;
        std::shared_ptr<PositiveRealParameter> height_ = std::make_shared<PositiveRealParameter>(0.0);
        std::shared_ptr<PositiveRealParameter> stored_height_ = std::make_shared<PositiveRealParameter>(0.0);
        bool is_dirty_ = true;

        void add_ln_relative_node_height_prior_density(
                double& density,
                std::vector< std::shared_ptr<PositiveRealParameter> >& parameters) const {
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
                child_iter->add_ln_relative_node_height_prior_density(density, parameters);
            }
        }

    public:
        Split split_;

        // Constructors
        BaseNode() { }
        BaseNode(int index) {
            this->index_ = index;
        }
        BaseNode(std::string label) {
            this->label_ = label;
        }
        BaseNode(int index, std::string label) {
            this->index_ = index;
            this->label_ = label;
        }
        BaseNode(double height) {
            this->height_->set_value(height);
        }
        BaseNode(int index, double height) {
            this->index_ = index;
            this->height_->set_value(height);
        }
        BaseNode(std::string label, double height) {
            this->label_ = label;
            this->height_->set_value(height);
        }
        BaseNode(int index, std::string label, double height) {
            this->index_ = index;
            this->label_ = label;
            this->height_->set_value(height);
        }
        BaseNode(std::shared_ptr<PositiveRealParameter> height) {
            this->height_ = height;
        }
        BaseNode(int index, std::shared_ptr<PositiveRealParameter> height) {
            this->index_ = index;
            this->height_ = height;
        }
        BaseNode(std::string label, std::shared_ptr<PositiveRealParameter> height) {
            this->label_ = label;
            this->height_ = height;
        }
        BaseNode(int index, std::string label, std::shared_ptr<PositiveRealParameter> height) {
            this->index_ = index;
            this->label_ = label;
            this->height_ = height;
        }

        // Destructor
        // virtual ~BaseNode() { }

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

        /**
         * Method to populate non-height related data (e.g., pop size) from
         * node comments in the form of a map.  Nothing to do for BaseNode, but
         * this is intended for derived classes to override
         */
        virtual void extract_data_from_node_comments(
                std::map<std::string, std::string> comment_map) { }

        unsigned int degree() const {
            unsigned int d = children_.size();
            if (this->has_parent()) {
                d += 1;
            }
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

        bool has_parent() const { return this->parent_.expired() ? false : true; }
        bool is_root() const { return this->parent_.expired() ? true : false; }

        unsigned int get_number_of_parents() const { return this->parent_.expired() ? 0 : 1; }

        const std::shared_ptr<DerivedNodeT>& get_parent() const {
            return this->parent_.lock();
        }
        std::shared_ptr<DerivedNodeT> get_parent() {
            return this->parent_.lock();
        }

        bool is_parent(const std::shared_ptr<DerivedNodeT>& node) const {
            if (this->parent_.lock() == node) {
                return true;
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
            return (this->parent_.lock()->is_ancestor(node));
        }

        void add_parent(std::shared_ptr<DerivedNodeT> node) {
            if (!node) {
                throw EcoevolityNullPointerError("BaseNode::add_parent(), empty node given");
            }
            if (! this->parent_.expired()) {
                throw EcoevolityError("BaseNode::add_parent(), this node already has a parent");
            }
            this->parent_ = node;
            if (! node->is_child(this->shared_from_this())) {
                node->add_child(this->shared_from_this());
            }
            node->make_dirty();
        }

        std::shared_ptr<DerivedNodeT> remove_parent() {
            if (this->has_parent()) {
                std::shared_ptr<DerivedNodeT> p = this->parent_.lock();
                this->parent_.reset();
                p->remove_child(this->shared_from_this());
                p->make_dirty();
                return p;
            }
            return nullptr;
        }

        unsigned int collapse() {
            if (! this->has_parent()) {
                throw EcoevolityError("BaseNode::collapse(), node has no parent");
            }
            std::shared_ptr<DerivedNodeT> grand_parent = this->get_parent();
            while (this->has_children()) {
                std::shared_ptr<DerivedNodeT> child = this->get_child(0);
                child->remove_parent();
                child->add_parent(grand_parent);
            }
            this->remove_parent();
            grand_parent->make_dirty();
            return grand_parent->get_number_of_children();
        }

        void split_children_from_polytomy(
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
                child->remove_parent();
                child->add_parent(new_node);
            }
            // If all children were assigned to the new node, we need to remove
            // the old polytomy node that now has no children
            if (! this->has_children()) {
                this->get_parent()->add_child(new_node);
                this->remove_parent();
            }
            else {
                this->add_child(new_node);
            }
            new_node->finish_initializing_inserted_internal_node();
        }

        virtual void finish_initializing_inserted_internal_node() { }

        bool has_children() const { return !this->children_.empty(); }
        bool is_leaf() const { return this->children_.empty(); }

        unsigned int get_number_of_children() const { return this->children_.size(); }
        const std::shared_ptr<DerivedNodeT>& get_child(unsigned int index) const {
            if (index >= this->children_.size()) {
                throw std::out_of_range("BaseNode::get_child() index out of range");
            }
            return this->children_.at(index);
        }

        std::shared_ptr<DerivedNodeT> get_child(unsigned int index) {
            if (index >= this->children_.size()) {
                throw std::out_of_range("BaseNode::get_child() index out of range");
            }
            return this->children_.at(index);
        }

        std::shared_ptr<DerivedNodeT> get_oldest_child() {
            if (this->is_leaf()) {
                throw EcoevolityError("called BaseNode::get_oldest_child() on a leaf");
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
                throw EcoevolityNullPointerError("BaseNode::add_child(), empty node given");
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
                throw EcoevolityNullPointerError("BaseNode::remove_child(), empty node given");
            }
            for (unsigned int i = 0; i < this->children_.size(); ++i) {
                if (this->children_.at(i) == node) {
                    this->children_.at(i).reset();
                    this->children_.erase(this->children_.begin() + i);
                    node->remove_parent();
                }
            }
            this->make_dirty();
        }

        std::shared_ptr<DerivedNodeT> remove_child(unsigned int index) {
            if (index >= this->children_.size()) {
                throw std::out_of_range("BaseNode::remove_child() index out of range");
            }
            std::shared_ptr<DerivedNodeT> c = this->children_.at(index);
            this->children_.at(index).reset();
            this->children_.erase(this->children_.begin() + index);
            c->remove_parent();
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

        void scale(double multiplier) {
            this->set_height(this->get_height() * multiplier);
            for (auto child_iter: this->children_) {
                child_iter->scale(multiplier);
            }
        }

        std::shared_ptr<PositiveRealParameter> get_height_parameter() const {
            return this->height_;
        }

        void update_height(double height) {
            this->height_->update_value(height);
            this->make_all_dirty();
        }

        void get_node(const std::string& label,
                std::shared_ptr<DerivedNodeT>& node_ptr) {
            if (node_ptr != nullptr) {
                return;
            }
            if (this->is_label(label)) {
                node_ptr = this->shared_from_this();
                return;
            }
            for (auto child_iter: this->children_) {
                if (child_iter->is_label(label)) {
                    node_ptr = child_iter;
                    return;
                }
                child_iter->get_node(label, node_ptr);
            }
        }
        std::shared_ptr<DerivedNodeT> get_node(const std::string& label) {
            std::shared_ptr<DerivedNodeT> node_ptr = nullptr;
            this->get_node(label, node_ptr);
            return node_ptr;
        }

        void store_height() {
            this->height_->store();
        }
        void store_all_heights() {
            this->height_->store();
            for (auto child_iter: this->children_) {
                child_iter->store_all_heights();
            }
        }
        virtual void store_all_parameter_values() {
            this->store_all_heights();
        }

        void restore_height() {
            this->height_->restore();
            this->make_all_dirty();
        }
        void restore_all_heights() {
            this->height_->restore();
            this->make_dirty();
            for (auto child_iter: this->children_) {
                child_iter->restore_all_heights();
            }
        }
        virtual void restore_all_parameter_values() {
            this->restore_all_heights();
        }

        void store_height_pointer() {
            this->stored_height_ = this->height_;
        }
        void store_all_height_pointers() {
            this->stored_height_ = this->height_;
            for (auto child_iter: this->children_) {
                child_iter->store_all_height_pointers();
            }
        }
        virtual void store_all_parameter_pointers() {
            this->store_all_height_pointers();
        }

        void restore_height_pointer() {
            this->height_ = this->stored_height_;
            this->make_all_dirty();
        }
        void restore_all_height_pointers() {
            this->height_ = this->stored_height_;
            this->make_dirty();
            for (auto child_iter: this->children_) {
                child_iter->restore_all_height_pointers();
            }
        }
        virtual void restore_all_parameter_pointers() {
            this->restore_all_height_pointers();
        }

        double get_length() const {
            if (this->has_parent()) {
                return this->parent_.lock()->get_height() - this->get_height();
            }
            return 0.0;
        }

        void sum_clade_length(double & length) const  {
            length += this->get_length();
            for (auto child_iter: this->children_) {
                child_iter->sum_clade_length(length);
            }
        }

        double get_clade_length() const {
            double l = 0.0;
            this->sum_clade_length(l);
            return l;
        }

        bool node_height_is_valid() const {
            double height = this->get_height();
            if (this->has_parent()) {
                if (this->parent_.lock()->get_height() < height) {
                    return false;
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

        bool node_heights_are_valid() const {
            if (! this->node_height_is_valid()) {
                return false;
            }
            for (auto child_iter: this->children_) {
                if (! child_iter->node_heights_are_valid()) {
                    return false;
                }
            }
            return true;
        }

        virtual int get_index() const {
            ECOEVOLITY_ASSERT(this->index_ > -1);
            return this->index_;
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

        bool clade_has_dirt() const {
            if (this->is_dirty()) {
                return true;
            }
            for (auto child_iter: this->children_) {
                if (child_iter->clade_has_dirt()) {
                    return true;
                }
            }
            return false;
        }

        void make_dirty() {
            this->is_dirty_ = true;
        }
        void make_all_dirty() {
            this->is_dirty_ = true;
            for (auto child_iter: this->children_) {
                child_iter->make_all_dirty();
            }
        }
        void make_clean() {
            this->is_dirty_ = false;
        }
        void make_all_clean() {
            this->is_dirty_ = false;
            for (auto child_iter: this->children_) {
                child_iter->make_all_clean();
            }
        }

        unsigned int get_node_count() const {
            unsigned int n = 1;
            for (auto child_iter: this->children_) {
                n += child_iter->get_node_count();
            }
            return n;
        }
        unsigned int get_leaf_node_count() const {
            if (this->is_leaf()) {
                return 1;
            }
            unsigned int n = 0;
            for (auto child_iter: this->children_) {
                n += child_iter->get_leaf_node_count();
            }
            return n;
        }
        unsigned int get_internal_node_count() const {
            if (this->is_leaf()) {
                return 0;
            }
            unsigned int n = 1;
            for (auto child_iter: this->children_) {
                n += child_iter->get_internal_node_count();
            }
            return n;
        }
        unsigned int get_polytomy_node_count() const {
            if (this->is_polytomy()) {
                return 1;
            }
            unsigned int n = 0;
            for (auto child_iter: this->children_) {
                n += child_iter->get_polytomy_node_count();
            }
            return n;
        }
        unsigned int get_mapped_node_count(
                std::shared_ptr<PositiveRealParameter> node_height_pointer) const {
            if (this->get_height_parameter() == node_height_pointer) {
                return 1;
            }
            unsigned int n = 0;
            for (auto child_iter: this->children_) {
                n += child_iter->get_mapped_node_count(node_height_pointer);
            }
            return n;
        }
        unsigned int get_mapped_polytomy_node_count(
                std::shared_ptr<PositiveRealParameter> node_height_pointer) const {
            if ((this->get_height_parameter() == node_height_pointer) &&
                    (this->is_polytomy())) {
                return 1;
            }
            unsigned int n = 0;
            for (auto child_iter: this->children_) {
                n += child_iter->get_mapped_polytomy_node_count(node_height_pointer);
            }
            return n;
        }

        void get_nodes(std::vector< std::shared_ptr<DerivedNodeT> >& nodes) {
            nodes.push_back(this->shared_from_this());
            for (auto child_iter: this->children_) {
                child_iter->get_nodes(nodes);
            }
        }
        std::vector< std::shared_ptr<DerivedNodeT> > get_nodes() {
            std::vector< std::shared_ptr<DerivedNodeT> > nodes;
            this->get_nodes(nodes);
            return nodes;
        }

        void pre_order(std::vector< std::shared_ptr<DerivedNodeT> >& nodes) {
            nodes.clear();
            this->get_nodes(nodes);
        }

        void level_order(std::vector< std::shared_ptr<DerivedNodeT> >& nodes) {
            nodes.clear();
            std::queue< std::shared_ptr<DerivedNodeT> > q;
            std::shared_ptr<DerivedNodeT> nd = this->shared_from_this();
            q.push(nd);
            while(! q.empty()) {
                nd = q.front(); q.pop();
                nodes.push_back(nd);
                if (! nd->is_leaf()) {
                    for (auto child : nd->children_) {
                        q.push(child);
                    }
                }
            }
        }

        void get_internal_nodes(std::vector< std::shared_ptr<DerivedNodeT> >& internal_nodes) {
            if (! this->is_leaf()) {
                internal_nodes.push_back(this->shared_from_this());
            }
            for (auto child_iter: this->children_) {
                child_iter->get_internal_nodes(internal_nodes);
            }
        }
        std::vector< std::shared_ptr<DerivedNodeT> > get_internal_nodes() {
            std::vector< std::shared_ptr<DerivedNodeT> > internal_nodes;
            internal_nodes.reserve(this->get_internal_node_count());
            this->get_internal_nodes(internal_nodes);
            return internal_nodes;
        }

        void get_polytomy_nodes(std::vector< std::shared_ptr<DerivedNodeT> >& polytomy_nodes) {
            if (this->is_polytomy()) {
                polytomy_nodes.push_back(this->shared_from_this());
            }
            for (auto child_iter: this->children_) {
                child_iter->get_polytomy_nodes(polytomy_nodes);
            }
        }
        std::vector< std::shared_ptr<DerivedNodeT> > get_polytomy_nodes() {
            std::vector< std::shared_ptr<DerivedNodeT> > nodes;
            this->get_polytomy_nodes(nodes);
            return nodes;
        }

        void get_mapped_nodes(std::vector< std::shared_ptr<DerivedNodeT> >& mapped_nodes,
                std::shared_ptr<PositiveRealParameter> node_height_pointer) {
            if (this->get_height_parameter() == node_height_pointer) {
                mapped_nodes.push_back(this->shared_from_this());
            }
            for (auto child_iter: this->children_) {
                child_iter->get_mapped_nodes(mapped_nodes, node_height_pointer);
            }
        }
        std::vector< std::shared_ptr<DerivedNodeT> > get_mapped_nodes(
                std::shared_ptr<PositiveRealParameter> node_height_pointer) {
            std::vector< std::shared_ptr<DerivedNodeT> > mapped_nodes;
            this->get_mapped_nodes(mapped_nodes, node_height_pointer);
            return mapped_nodes;
        }

        void get_mapped_polytomy_nodes(std::vector< std::shared_ptr<DerivedNodeT> >& mapped_nodes,
                std::shared_ptr<PositiveRealParameter> node_height_pointer) {
            if ((this->get_height_parameter() == node_height_pointer) &&
                    (this->is_polytomy())) {
                mapped_nodes.push_back(this->shared_from_this());
            }
            for (auto child_iter: this->children_) {
                child_iter->get_mapped_polytomy_nodes(mapped_nodes, node_height_pointer);
            }
        }
        std::vector< std::shared_ptr<DerivedNodeT> > get_mapped_polytomy_nodes(
                std::shared_ptr<PositiveRealParameter> node_height_pointer) {
            std::vector< std::shared_ptr<DerivedNodeT> > mapped_nodes;
            this->get_mapped_polytomy_nodes(mapped_nodes, node_height_pointer);
            return mapped_nodes;
        }

        void get_leaves(std::vector< std::shared_ptr<DerivedNodeT> >& leaves) {
            if (this->is_leaf()) {
                leaves.push_back(this->shared_from_this());
            }
            for (auto child_iter: this->children_) {
                child_iter->get_leaves(leaves);
            }
        }
        std::vector< std::shared_ptr<DerivedNodeT> > get_leaves() {
            std::vector< std::shared_ptr<DerivedNodeT> > leaves;
            leaves.reserve(this->get_leaf_node_count());
            this->get_leaves(leaves);
            return leaves;
        }

        void get_leaf_labels(std::vector<std::string>& leaf_labels) const {
            if (this->is_leaf()) {
                leaf_labels.push_back(this->get_label());
            }
            for (auto child_iter: this->children_) {
                child_iter->get_leaf_labels(leaf_labels);
            }
        }
        std::vector<std::string> get_leaf_labels() const {
            std::vector<std::string> leaf_labels;
            leaf_labels.reserve(this->get_leaf_node_count());
            this->get_leaf_labels(leaf_labels);
            return leaf_labels;
        }

        std::shared_ptr<ContinuousProbabilityDistribution> get_node_height_prior() const {
            return this->height_->get_prior();
        }
        void set_node_height_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            this->height_->set_prior(prior);
            this->make_all_dirty();
        }
        void set_all_node_height_priors(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            this->height_->set_prior(prior);
            this->make_dirty();
            for (auto child_iter: this->children_) {
                child_iter->set_all_node_height_priors(prior);
            }
        }

        void fix_node_height() {
            this->height_->fix();
        }
        void fix_all_node_heights() {
            this->height_->fix();
            for (auto child_iter: this->children_) {
                child_iter->fix_all_node_heights();
            }
        }
        void estimate_node_height() {
            this->height_->estimate();
        }
        void estimate_all_node_heights() {
            this->height_->estimate();
            for (auto child_iter: this->children_) {
                child_iter->estimate_all_node_heights();
            }
        }

        bool node_height_is_fixed() const {
            return this->height_->is_fixed();
        }

        bool all_node_heights_are_fixed() const {
            if (! this->node_height_is_fixed()) {
                return false;
            }
            for (auto child_iter: this->children_) {
                if (! child_iter->all_node_heights_are_fixed()) {
                    return false;
                }
            }
            return true;
        }

        double get_height_relative_prior_ln_pdf() const {
            return this->height_->relative_prior_ln_pdf();
        }

        double calculate_ln_relative_node_height_prior_density() const {
            double d = 0.0;
            std::vector< std::shared_ptr<PositiveRealParameter> > parameters;
            parameters.reserve(this->get_node_count());
            this->add_ln_relative_node_height_prior_density(d, parameters);
            return d;
        }

        std::string to_parentheses(unsigned int precision = 12) const {
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
                    s << child_iter->to_parentheses(precision);
                    ++child_idx;
                }
                s << ")";
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
