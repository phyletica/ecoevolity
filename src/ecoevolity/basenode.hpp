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

#include "matrix.hpp"
#include "parameter.hpp"
#include "debug.hpp"
#include "assert.hpp"
#include "error.hpp"

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
class BaseNode {
    protected:
        std::vector< DerivedNodeT * > children_;
        DerivedNodeT * parent_ = 0;
        std::string label_ = "";
        PositiveRealParameter * height_ = new PositiveRealParameter(0.0);
        bool is_dirty_ = true;

    public:
        // Constructors
        BaseNode() { }
        BaseNode(const DerivedNodeT& node) {
            this->children_ = node.children_;
            this->parent_ = node.parent_;
            this->height_->set_value(node.height_->get_value());
            this->label_ = node.label_;
            this->is_dirty_ = node.is_dirty_;
        }
        BaseNode(std::string label) {
            this->label_ = label;
        }
        BaseNode(double height) {
            this->height_->set_value(height);
        }
        BaseNode(std::string label, double height) {
            this->label_ = label;
            this->height_->set_value(height);
        }

        // Destructor
        virtual ~BaseNode() {
            delete this->height_;
            if (this->parent_) {
                this->parent_->remove_child(static_cast<DerivedNodeT *>(this));
            }
            for (auto child_iter: this->children_) {
                child_iter->remove_parent();
            }
        }

        DerivedNodeT& operator=(const DerivedNodeT& node) {
            this->children_ = node.children_;
            this->parent_ = node.parent_;
            this->height_->set_value(node.height_->get_value());
            this->label_ = node.label_;
            this->is_dirty_ = node.is_dirty_;
            return * this;
        }

        DerivedNodeT* clone() const {
            return new DerivedNodeT(static_cast<DerivedNodeT const &>(* this));
        }
        
        //Methods
        unsigned int degree() const {
            unsigned int d = children_.size();
            if (this->has_parent()) {
                d += 1;
            }
            return d;
        }

        bool has_parent() const { return parent_ ? true : false; }
        bool is_root() const { return parent_ ? false : true; }

        unsigned int get_number_of_parents() const { return parent_ ? 1 : 0; }

        const DerivedNodeT* get_parent() const {
            return this->parent_;
        }
        DerivedNodeT* get_parent() {
            return this->parent_;
        }

        bool is_parent(const DerivedNodeT* node) const {
            if (this->parent_ ==  node) {
                return true;
            }
            return false;
        }

        void add_parent(DerivedNodeT* node) {
            if (!node) {
                throw EcoevolityNullPointerError("BaseNode::add_parent(), empty node given");
            }
            if (this->parent_) {
                throw EcoevolityError("BaseNode::add_parent(), this node already has a parent");
            }
            this->parent_ = node;
            if (! node->is_child(static_cast<DerivedNodeT *>(this))) {
                node->add_child(static_cast<DerivedNodeT *>(this));
            }
        }

        DerivedNodeT* remove_parent() {
            if (this->has_parent()) {
                DerivedNodeT* p = this->parent_;
                this->parent_ = 0;
                p->remove_child(static_cast<DerivedNodeT *>(this));
                return p;
            }
            return 0;
        }

        bool has_children() const { return !this->children_.empty(); }
        bool is_leaf() const { return this->children_.empty(); }

        unsigned int get_number_of_children() const { return this->children_.size(); }
        const DerivedNodeT* get_child(unsigned int index) const {
            if (index >= this->children_.size()) {
                throw std::out_of_range("BaseNode::get_child() index out of range");
            }
            return this->children_.at(index);
        }

        DerivedNodeT* get_child(unsigned int index) {
            if (index >= this->children_.size()) {
                throw std::out_of_range("BaseNode::get_child() index out of range");
            }
            return this->children_.at(index);
        }

        bool is_child(const DerivedNodeT* node) const {
            for (auto child_iter: this->children_) {
                if (child_iter == node) {
                    return true;
                }
            }
            return false;
        }

        void add_child(DerivedNodeT* node) {
            if (!node) {
                throw EcoevolityNullPointerError("BaseNode::add_child(), empty node given");
            }
            if (! this->is_child(node)) {
                this->children_.push_back(node);
            }
            if (! node->is_parent(static_cast<DerivedNodeT *>(this))) {
                node->add_parent(static_cast<DerivedNodeT *>(this));
            }
        }

        void remove_child(DerivedNodeT* node) {
            if (!node) {
                throw EcoevolityNullPointerError("BaseNode::remove_child(), empty node given");
            }
            for (unsigned int i = 0; i < this->children_.size(); ++i) {
                if (this->children_.at(i) == node) {
                    this->children_.erase(this->children_.begin() + i);
                    node->remove_parent();
                }
            }
        }

        DerivedNodeT* remove_child(unsigned int index) {
            if (index >= this->children_.size()) {
                throw std::out_of_range("BaseNode::remove_child() index out of range");
            }
            DerivedNodeT* c = this->children_.at(index);
            this->children_.erase(this->children_.begin() + index);
            c->remove_parent();
            return c;
        }

        const double& get_height() const {
            return this->height_->get_value();
        }

        void set_height(double height) {
            this->height_->set_value(height);
            this->make_dirty();
            for (auto child_iter: this->children_) {
                child_iter->make_dirty();
            }
        }

        void set_height_parameter(PositiveRealParameter * height_parameter) {
            this->height_ = height_parameter;
        }

        PositiveRealParameter * get_height_parameter() {
            return this->height_;
        }

        void update_height(double height) {
            this->height_->update_value(height);
            this->make_dirty();
            for (auto child_iter: this->children_) {
                child_iter->make_dirty();
            }
        }

        void store_height() {
            this->height_->store();
        }

        void restore_height() {
            this->height_->restore();
        }

        double get_length() const {
            if (this->has_parent()) {
                return this->get_parent()->get_height() - this->get_height();
            }
            return 0.0;
        }

        const std::string& get_label() const {
            return this->label_;
        }
        void set_label(std::string label) {
            this->label_ = label;
        }

        const bool& is_dirty() const {
            return this->is_dirty_;
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
};

#endif
