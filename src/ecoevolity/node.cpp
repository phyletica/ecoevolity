/******************************************************************************
 * Copyright (C) 2016 Jamie R. Oaks.
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

#include "node.hpp"

Node::Node(const Node& node) {
    this->children_ = node.children_;
    this->parent_ = node.parent_;
    this->height_ = node.height_;
    this->label_ = node.label_;
    this->is_dirty_ = node.is_dirty_;
}

Node::Node(std::string label) {
    this->label_ = label;
}

Node::~Node() {
    if (this->parent_) {
        this->parent_->remove_child(this);
    }
    for (auto child_iter: this->children_) {
        child_iter->remove_parent();
    }
}

Node::Node(double height) {
    this->height_ = height;
}

Node::Node(std::string label, double height) {
    this->label_ = label;
    this->height_ = height;
}

Node& Node::operator=(const Node& node) {
    this->children_ = node.children_;
    this->parent_ = node.parent_;
    this->height_ = node.height_;
    this->label_ = node.label_;
    this->is_dirty_ = node.is_dirty_;
    return * this;
}

unsigned int Node::degree() const {
    unsigned int d = children_.size();
    if (this->has_parent()) {
        d += 1;
    }
    return d;
}

const Node* Node::get_parent() const {
    return this->parent_;
}
Node* Node::get_parent() {
    return this->parent_;
}

bool Node::is_parent(const Node* node) const {
    if (this->parent_ ==  node) {
        return true;
    }
    return false;
}

void Node::add_parent(Node* node) {
    if (!node) {
        throw EcoevolityNullPointerError("Node::add_parent(), empty node given");
    }
    if (this->parent_) {
        throw EcoevolityError("Node::add_parent(), this node already has a parent");
    }
    this->parent_ = node;
    if (! node->is_child(this)) {
        node->add_child(this);
    }
}

Node* Node::remove_parent() {
    if (this->has_parent()) {
        Node* p = this->parent_;
        this->parent_ = 0;
        p->remove_child(this);
        return p;
    }
    return 0;
}

const Node* Node::get_child(unsigned int index) const {
    if (index >= this->children_.size()) {
        throw std::out_of_range("Node::get_child() index out of range");
    }
    return this->children_.at(index);
}

Node* Node::get_child(unsigned int index) {
    if (index >= this->children_.size()) {
        throw std::out_of_range("Node::get_child() index out of range");
    }
    return this->children_.at(index);
}

bool Node::is_child(const Node* node) const {
    for (auto child_iter: this->children_) {
        if (child_iter == node) {
            return true;
        }
    }
    return false;
}

void Node::add_child(Node* node) {
    if (!node) {
        throw EcoevolityNullPointerError("Node::add_child(), empty node given");
    }
    if (! this->is_child(node)) {
        this->children_.push_back(node);
    }
    if (! node->is_parent(this)) {
        node->add_parent(this);
    }
}

void Node::remove_child(Node* node) {
    if (!node) {
        throw EcoevolityNullPointerError("Node::remove_child(), empty node given");
    }
    for (unsigned int i = 0; i < this->children_.size(); ++i) {
        if (this->children_.at(i) == node) {
            this->children_.erase(this->children_.begin() + i);
            node->remove_parent();
        }
    }
}

Node* Node::remove_child(unsigned int index) {
    if (index >= this->children_.size()) {
        throw std::out_of_range("Node::remove_child() index out of range");
    }
    Node* c = this->children_.at(index);
    this->children_.erase(this->children_.begin() + index);
    c->remove_parent();
    return c;
}

const bool& Node::is_dirty() const {
    return this->is_dirty_;
}

void Node::make_dirty() {
    this->is_dirty_ = true;
}

void Node::make_all_dirty() {
    this->is_dirty_ = true;
    for (auto child_iter: this->children_) {
        child_iter->make_all_dirty();
    }
}

void Node::make_clean() {
    this->is_dirty_ = false;
}

void Node::make_all_clean() {
    this->is_dirty_ = false;
    for (auto child_iter: this->children_) {
        child_iter->make_all_clean();
    }
}

const double& Node::get_height() const {
    return this->height_;
}

void Node::set_height(double height) {
    this->height_ = height;
    this->make_dirty();
    for (auto child_iter: this->children_) {
        child_iter->make_dirty();
    }
}

double Node::get_length() const {
    if (this->has_parent()) {
        return this->get_parent()->get_height() - this->get_height();
    }
    return 0.0;
}

const std::string& Node::get_label() const {
    return this->label_;
}

void Node::set_label(std::string label) {
    this->label_ = label;
}

unsigned int Node::get_node_count() const {
    unsigned int n = 1;
    for (auto child_iter: this->children_) {
        n += child_iter->get_node_count();
    }
    return n;
}

unsigned int Node::get_leaf_node_count() const {
    if (this->is_leaf()) {
        return 1;
    }
    unsigned int n = 0;
    for (auto child_iter: this->children_) {
        n += child_iter->get_leaf_node_count();
    }
    return n;
}

unsigned int Node::get_internal_node_count() const {
    if (this->is_leaf()) {
        return 0;
    }
    unsigned int n = 1;
    for (auto child_iter: this->children_) {
        n += child_iter->get_internal_node_count();
    }
    return n;
}
