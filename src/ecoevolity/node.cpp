#include "node.hpp"

Node::~Node() {
    if (this->parent_) {
        this->parent_->remove_child(this);
    }
    for (auto child_iter: this->children_) {
        child_iter->remove_parent();
    }
}

Node::Node(const Node& node) {
    this->children_ = node.children_;
    this->parent_ = node.parent_;
}

Node& Node::operator=(const Node& node) {
    this->children_ = node.children_;
    this->parent_ = node.parent_;
    return * this;
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
        throw EvoevolityNullPointerError("Node::add_parent(), empty node given");
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
    if (index < 0 || index >= this->children_.size()) {
        throw std::out_of_range("Node::get_child() index out of range");
    }
    return this->children_.at(index);
}

Node* Node::get_child(unsigned int index) {
    if (index < 0 || index >= this->children_.size()) {
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
        throw EvoevolityNullPointerError("Node::add_child(), empty node given");
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
        throw EvoevolityNullPointerError("Node::remove_child(), empty node given");
    }
    for (unsigned int i = 0; i < this->children_.size(); ++i) {
        if (this->children_.at(i) == node) {
            this->children_.erase(this->children_.begin() + i);
            node->remove_parent();
        }
    }
}

Node* Node::remove_child(unsigned int index) {
    if (index < 0 || index >= this->children_.size()) {
        throw std::out_of_range("Node::remove_child() index out of range");
    }
    Node* c = this->children_.at(index);
    this->children_.erase(this->children_.begin() + index);
    c->remove_parent();
    return c;
}
