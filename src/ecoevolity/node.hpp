#ifndef ECOEVOLITY_NODE_HPP
#define ECOEVOLITY_NODE_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "matrix.hpp"
#include "debug.hpp"
#include "assert.hpp"
#include "error.hpp"

class Node {
    private:
        std::vector< Node * > children_;
        Node * parent_ = 0;
        std::string label_ = "";
        double height_ = 0.0;
        bool is_dirty_ = true;

    public:
        // Constructors
        Node() { };
        Node(const Node& node);
        Node(std::string label);
        Node(double height);
        Node(std::string label, double height);

        // Destructor
        virtual ~Node();

        Node& operator=(const Node& node);

        Node* clone() const {
            return new Node(* this);
        }
        
        //Methods
        unsigned int degree() const { return children_.size() + parent_ ? 1 : 0; }

        bool has_parent() const { return parent_ ? true : false; }
        bool is_root() const { return parent_ ? false : true; }

        unsigned int get_number_of_parents() const { return parent_ ? 1 : 0; }

        const Node* get_parent() const;
        Node* get_parent();

        virtual bool is_parent(const Node* node) const;

        virtual void add_parent(Node* node);

        virtual Node* remove_parent();

        bool has_children() const { return !this->children_.empty(); }
        bool is_leaf() const { return this->children_.empty(); }

        unsigned int get_number_of_children() const { return this->children_.size(); }
        const Node* get_child(unsigned int index) const;
        Node* get_child(unsigned int index);

        virtual bool is_child(const Node* node) const;

        virtual void add_child(Node* node);

        virtual void remove_child(Node* node);

        virtual Node* remove_child(unsigned int index);

        const double& get_height() const;
        void set_height(double height);
        double get_length() const;

        const std::string& get_label() const;
        void set_label(std::string label);

        const bool& is_dirty() const;

        void make_dirty();
        void make_all_dirty();
        void make_clean();
        void make_all_clean();

        unsigned int get_node_count() const;
        unsigned int get_leaf_node_count() const;
        unsigned int get_internal_node_count() const;
};


// class PopulationNode: public Node {
//     private:
//         BiallelicPatternProbabilityMatrix pattern_probs_bottom_;
//         BiallelicPatternProbabilityMatrix pattern_probs_top_;
// 
//     public:
//         // Node(unsigned int allele_count);
//         // Node(std::string label, unsigned int allele_count);
// 
//         // from NodeData
//         unsigned int get_allele_count() const;
//         void resize(unsigned int allele_count);
//         void reset(unsigned int allele_count);
// }

#endif
