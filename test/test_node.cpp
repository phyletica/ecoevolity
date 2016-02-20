#include "catch.hpp"
#include "ecoevolity/node.hpp"

TEST_CASE("Testing constructors of Node", "[Node]") {

    SECTION("Testing bare constructor") {
        Node n = Node();
        REQUIRE(n.get_height() == 0.0);
        REQUIRE(n.get_length() == 0.0);
        REQUIRE(n.get_label() == "");
        REQUIRE(n.is_dirty());

        REQUIRE(n.degree() == 0);
        REQUIRE(n.has_parent() == false);
        REQUIRE(n.get_number_of_parents() == 0);
        REQUIRE(n.has_children() == false);
        REQUIRE(n.get_number_of_children() == 0);

        REQUIRE(n.is_leaf() == true);
        REQUIRE(n.is_root() == true);
        REQUIRE(n.get_node_count() == 1);
        REQUIRE(n.get_leaf_node_count() == 1);
        REQUIRE(n.get_internal_node_count() == 0);

        Node * p = n.get_parent();
        REQUIRE(p == NULL);

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);

        delete p;
    }

    SECTION("Testing label constructor") {
        Node n = Node("leaf1");
        REQUIRE(n.get_height() == 0.0);
        REQUIRE(n.get_length() == 0.0);
        REQUIRE(n.get_label() == "leaf1");
        REQUIRE(n.is_dirty());

        REQUIRE(n.degree() == 0);
        REQUIRE(n.has_parent() == false);
        REQUIRE(n.get_number_of_parents() == 0);
        REQUIRE(n.has_children() == false);
        REQUIRE(n.get_number_of_children() == 0);

        REQUIRE(n.is_leaf() == true);
        REQUIRE(n.is_root() == true);
        REQUIRE(n.get_node_count() == 1);
        REQUIRE(n.get_leaf_node_count() == 1);
        REQUIRE(n.get_internal_node_count() == 0);

        Node * p = n.get_parent();
        REQUIRE(p == NULL);

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);

        delete p;
    }

    SECTION("Testing height constructor") {
        Node n = Node(0.03);
        REQUIRE(n.get_height() == 0.03);
        REQUIRE(n.get_length() == 0.0);
        REQUIRE(n.get_label() == "");
        REQUIRE(n.is_dirty());

        REQUIRE(n.degree() == 0);
        REQUIRE(n.has_parent() == false);
        REQUIRE(n.get_number_of_parents() == 0);
        REQUIRE(n.has_children() == false);
        REQUIRE(n.get_number_of_children() == 0);

        REQUIRE(n.is_leaf() == true);
        REQUIRE(n.is_root() == true);
        REQUIRE(n.get_node_count() == 1);
        REQUIRE(n.get_leaf_node_count() == 1);
        REQUIRE(n.get_internal_node_count() == 0);

        Node * p = n.get_parent();
        REQUIRE(p == NULL);

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);

        delete p;
    }

    SECTION("Testing label and height constructor") {
        Node n = Node("leaf1", 0.02);
        REQUIRE(n.get_height() == 0.02);
        REQUIRE(n.get_length() == 0.0);
        REQUIRE(n.get_label() == "leaf1");
        REQUIRE(n.is_dirty());

        REQUIRE(n.degree() == 0);
        REQUIRE(n.has_parent() == false);
        REQUIRE(n.get_number_of_parents() == 0);
        REQUIRE(n.has_children() == false);
        REQUIRE(n.get_number_of_children() == 0);

        REQUIRE(n.is_leaf() == true);
        REQUIRE(n.is_root() == true);
        REQUIRE(n.get_node_count() == 1);
        REQUIRE(n.get_leaf_node_count() == 1);
        REQUIRE(n.get_internal_node_count() == 0);

        Node * p = n.get_parent();
        REQUIRE(p == NULL);

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);

        delete p;
    }

    SECTION("Testing node ref constructor") {
        Node n2 = Node("leaf1", 0.02);
        Node n = Node(n2);
        REQUIRE(n.get_height() == 0.02);
        REQUIRE(n.get_length() == 0.0);
        REQUIRE(n.get_label() == "leaf1");
        REQUIRE(n.is_dirty());

        REQUIRE(n.degree() == 0);
        REQUIRE(n.has_parent() == false);
        REQUIRE(n.get_number_of_parents() == 0);
        REQUIRE(n.has_children() == false);
        REQUIRE(n.get_number_of_children() == 0);

        REQUIRE(n.is_leaf() == true);
        REQUIRE(n.is_root() == true);
        REQUIRE(n.get_node_count() == 1);
        REQUIRE(n.get_leaf_node_count() == 1);
        REQUIRE(n.get_internal_node_count() == 0);

        Node * p = n.get_parent();
        REQUIRE(p == NULL);

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
        
        n.set_label("leaf2");
        REQUIRE(n.get_label() == "leaf2");
        REQUIRE(n2.get_label() == "leaf1");

        n2.set_height(0.06);
        REQUIRE(n.get_height() == 0.02);
        REQUIRE(n2.get_height() == 0.06);
    }

}

TEST_CASE("Testing copy operator of Node", "[Node]") {

    SECTION("Testing copy operator") {
        Node n2 = Node("leaf1", 0.02);
        Node n = n2;
        REQUIRE(n.get_height() == 0.02);
        REQUIRE(n.get_length() == 0.0);
        REQUIRE(n.get_label() == "leaf1");
        REQUIRE(n.is_dirty());

        REQUIRE(n.degree() == 0);
        REQUIRE(n.has_parent() == false);
        REQUIRE(n.get_number_of_parents() == 0);
        REQUIRE(n.has_children() == false);
        REQUIRE(n.get_number_of_children() == 0);

        REQUIRE(n.is_leaf() == true);
        REQUIRE(n.is_root() == true);
        REQUIRE(n.get_node_count() == 1);
        REQUIRE(n.get_leaf_node_count() == 1);
        REQUIRE(n.get_internal_node_count() == 0);

        Node * p = n.get_parent();
        REQUIRE(p == NULL);

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
        
        n.set_label("leaf2");
        REQUIRE(n.get_label() == "leaf2");
        REQUIRE(n2.get_label() == "leaf1");

        n2.set_height(0.06);
        REQUIRE(n.get_height() == 0.02);
        REQUIRE(n2.get_height() == 0.06);
    }

}

TEST_CASE("Testing clone method of Node", "[Node]") {

    SECTION("Testing clone") {
        Node * n2 = new Node("leaf1", 0.02);
        Node * n = n2->clone();
        REQUIRE(n->get_height() == 0.02);
        REQUIRE(n->get_length() == 0.0);
        REQUIRE(n->get_label() == "leaf1");
        REQUIRE(n->is_dirty());

        REQUIRE(n->degree() == 0);
        REQUIRE(n->has_parent() == false);
        REQUIRE(n->get_number_of_parents() == 0);
        REQUIRE(n->has_children() == false);
        REQUIRE(n->get_number_of_children() == 0);

        REQUIRE(n->is_leaf() == true);
        REQUIRE(n->is_root() == true);
        REQUIRE(n->get_node_count() == 1);
        REQUIRE(n->get_leaf_node_count() == 1);
        REQUIRE(n->get_internal_node_count() == 0);

        Node * p = n->get_parent();
        REQUIRE(p == NULL);

        REQUIRE_THROWS_AS(n->get_child(0), std::out_of_range);
        
        n->set_label("leaf2");
        REQUIRE(n->get_label() == "leaf2");
        REQUIRE(n2->get_label() == "leaf1");

        n2->set_height(0.06);
        REQUIRE(n->get_height() == 0.02);
        REQUIRE(n2->get_height() == 0.06);

        delete n;
        delete n2;
    }

}

TEST_CASE("Testing parent methods of Node", "[Node]") {

    SECTION("Testing parent methods") {
        Node c = Node();
        Node p = Node();

        REQUIRE(c.degree() == 0);
        REQUIRE(p.degree() == 0);
        REQUIRE(c.has_parent() == false);
        REQUIRE(c.get_number_of_parents() == 0);

        REQUIRE(c.get_parent() == NULL);
        REQUIRE_THROWS_AS(p.get_child(0), std::out_of_range);
        REQUIRE(c.is_parent(&p) == false);
        REQUIRE(p.is_child(&c) == false);

        c.add_parent(&p);
        REQUIRE(c.has_parent() == true);
        REQUIRE(c.is_parent(&p) == true);
        REQUIRE(p.is_child(&c) == true);
        REQUIRE(c.get_number_of_parents() == 1);
        REQUIRE(c.degree() == 1);
        REQUIRE(p.degree() == 1);

        Node * r = c.remove_parent();
        REQUIRE(r != NULL);
        REQUIRE(r == &p);
        REQUIRE(c.has_parent() == false);
        REQUIRE(c.is_parent(&p) == false);
        REQUIRE(p.is_child(&c) == false);
        REQUIRE(c.get_number_of_parents() == 0);
        REQUIRE(c.degree() == 0);
        REQUIRE(p.degree() == 0);
    }
}

TEST_CASE("Testing child methods of Node", "[Node]") {

    SECTION("Testing child methods") {
        Node c = Node();
        Node p = Node();

        REQUIRE(c.degree() == 0);
        REQUIRE(p.degree() == 0);
        REQUIRE(p.has_children() == false);
        REQUIRE(p.get_number_of_children() == 0);

        REQUIRE(c.get_parent() == NULL);
        REQUIRE_THROWS_AS(p.get_child(0), std::out_of_range);
        REQUIRE(c.is_parent(&p) == false);
        REQUIRE(p.is_child(&c) == false);

        p.add_child(&c);
        REQUIRE(p.has_children() == true);
        REQUIRE(c.is_parent(&p) == true);
        REQUIRE(p.is_child(&c) == true);
        REQUIRE(p.get_number_of_children() == 1);
        REQUIRE(c.degree() == 1);
        REQUIRE(p.degree() == 1);

        unsigned int i = 1;
        REQUIRE_THROWS_AS(p.remove_child(i), std::out_of_range);
        i = 0;
        Node * r = p.remove_child(i);
        REQUIRE(r == &c);
        REQUIRE(p.has_children() == false);
        REQUIRE(c.is_parent(&p) == false);
        REQUIRE(p.is_child(&c) == false);
        REQUIRE(p.get_number_of_children() == 0);
        REQUIRE(c.degree() == 0);
        REQUIRE(p.degree() == 0);

        p.add_child(&c);
        REQUIRE(p.has_children() == true);
        REQUIRE(c.is_parent(&p) == true);
        REQUIRE(p.is_child(&c) == true);
        REQUIRE(p.get_number_of_children() == 1);
        REQUIRE(c.degree() == 1);
        REQUIRE(p.degree() == 1);

        p.remove_child(&c);
        REQUIRE(p.has_children() == false);
        REQUIRE(c.is_parent(&p) == false);
        REQUIRE(p.is_child(&c) == false);
        REQUIRE(p.get_number_of_children() == 0);
        REQUIRE(c.degree() == 0);
        REQUIRE(p.degree() == 0);
        
    }
}
