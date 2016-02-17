#ifndef ECOEVOLITY_NODE_HPP
#define ECOEVOLITY_NODE_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "matrix.hpp"
#include "debug.hpp"
#include "assert.hpp"

class Node {
    public:
        // Constructor
        Node();
        Node(std::string label);
        // Destructor
        // ~Node();
        
        //Methods
        const double& get_height() const;
        void set_height(double height);

        const Node& get_parent() const;
        void set_parent(Node & parent);

        const std::vector<Node>& get_children() const;
        const Node& get_child(unsigned int child_index) const;
        void set_child(unsigned int child_index, Node & child);
        void add_child(Node & child);
        void remove_child(unsigned int child_index);
        const Node& get_left_child() const {
            if (this->children_.empty()) {
                return NULL;
            }
            return this->children_.at(0);
        }
        const Node& get_right_child() const {
            if (this->children_.size() < 2) {
                return NULL;
            }
            return this->children_.at(1);
        }
        unsigned int get_number_of_children();

        bool is_root() const {
            return (this->parent_ == NULL);
        }

        bool is_leaf() const {
            return (this->children_.empty());
        }

        const bool& is_dirty() const;

        void make_dirty();
        void make_all_dirty();
        void make_clean();

        double get_length() const;

        unsigned int get_node_count() const;
        unsigned int get_leaf_node_count() const;
        unsigned int get_internal_node_count() const;
        // from NodeData
        unsigned int get_number_of_alleles() const;


    private:
        std::string label_;
        double height_;
        std::vector<Node> children_;
        Node parent_ = NULL;
        BiallelicPatternProbabilityMatrix pattern_probs_bottom_;
        BiallelicPatternProbabilityMatrix pattern_probs_top_;

        bool is_dirty_ = false;

        //Methods
};

#endif