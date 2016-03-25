#include <typeinfo>

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
        REQUIRE(typeid(&n).hash_code() == typeid(p).hash_code());

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
        REQUIRE(typeid(&n).hash_code() == typeid(p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);

        delete p;
    }

    SECTION("Testing height constructor") {
        Node n = Node(0.03);
        REQUIRE(n.get_height() == Approx(0.03));
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
        REQUIRE(typeid(&n).hash_code() == typeid(p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);

        delete p;
    }

    SECTION("Testing label and height constructor") {
        Node n = Node("leaf1", 0.02);
        REQUIRE(n.get_height() == Approx(0.02));
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
        REQUIRE(typeid(&n).hash_code() == typeid(p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);

        delete p;
    }

    SECTION("Testing node ref constructor") {
        Node n2 = Node("leaf1", 0.02);
        Node n = Node(n2);
        REQUIRE(typeid(n).hash_code() == typeid(n2).hash_code());
        REQUIRE(n.get_height() == Approx(0.02));
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
        REQUIRE(typeid(&n).hash_code() == typeid(p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
        
        n.set_label("leaf2");
        REQUIRE(n.get_label() == "leaf2");
        REQUIRE(n2.get_label() == "leaf1");

        n2.set_height(0.06);
        REQUIRE(n.get_height() == Approx(0.02));
        REQUIRE(n2.get_height() == Approx(0.06));

        delete p;
    }

}

TEST_CASE("Testing copy operator of Node", "[Node]") {

    SECTION("Testing copy operator") {
        Node n2 = Node("leaf1", 0.02);
        Node n = n2;
        REQUIRE(typeid(n).hash_code() == typeid(n2).hash_code());
        REQUIRE(n.get_height() == Approx(0.02));
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
        REQUIRE(typeid(&n).hash_code() == typeid(p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
        
        n.set_label("leaf2");
        REQUIRE(n.get_label() == "leaf2");
        REQUIRE(n2.get_label() == "leaf1");

        n2.set_height(0.06);
        REQUIRE(n.get_height() == Approx(0.02));
        REQUIRE(n2.get_height() == Approx(0.06));

        delete p;
    }

}

TEST_CASE("Testing clone method of Node", "[Node]") {

    SECTION("Testing clone") {
        Node * n2 = new Node("leaf1", 0.02);
        Node * n = n2->clone();
        REQUIRE(typeid(n).hash_code() == typeid(n2).hash_code());
        REQUIRE(n->get_height() == Approx(0.02));
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
        REQUIRE(typeid(n).hash_code() == typeid(p).hash_code());

        REQUIRE_THROWS_AS(n->get_child(0), std::out_of_range);
        
        n->set_label("leaf2");
        REQUIRE(n->get_label() == "leaf2");
        REQUIRE(n2->get_label() == "leaf1");

        n2->set_height(0.06);
        REQUIRE(n->get_height() == Approx(0.02));
        REQUIRE(n2->get_height() == Approx(0.06));

        delete p;
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
        REQUIRE(typeid(r).hash_code() == typeid(&c).hash_code());
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
        REQUIRE(typeid(c).hash_code() == typeid(p).hash_code());

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
        REQUIRE(typeid(&c).hash_code() == typeid(r).hash_code());
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

TEST_CASE("Test simple tree building with Node", "[Node]") {

    SECTION("Testing tree building") {
        Node root = Node("root", 1.0);
        Node root_child1 = Node("root child 1", 0.8);
        Node root_child2 = Node("root child 2", 0.3);
        Node leaf1 = Node("leaf 1");
        Node leaf2 = Node("leaf 2");
        Node leaf3 = Node("leaf 3");
        Node leaf4 = Node("leaf 4");
        Node leaf5 = Node("leaf 5");

        root.add_child(&root_child1);
        root.add_child(&root_child2);

        leaf1.add_parent(&root_child1);
        leaf2.add_parent(&root_child1);
        leaf3.add_parent(&root_child1);

        leaf4.add_parent(&root_child2);
        leaf5.add_parent(&root_child2);

        REQUIRE(root.degree() == 2);
        REQUIRE(root_child1.degree() == 4);
        REQUIRE(root_child1.get_number_of_children() == 3);
        REQUIRE(root_child2.degree() == 3);
        REQUIRE(leaf1.degree() == 1);
        REQUIRE(leaf2.degree() == 1);
        REQUIRE(leaf3.degree() == 1);
        REQUIRE(leaf4.degree() == 1);
        REQUIRE(leaf5.degree() == 1);

        REQUIRE(root.is_root());
        REQUIRE(! root.is_leaf());
        REQUIRE(root.has_children());
        REQUIRE(! root.has_parent());
        REQUIRE(root.is_dirty());
        REQUIRE(root.get_node_count() == 8);
        REQUIRE(root.get_leaf_node_count() == 5);
        REQUIRE(root.get_internal_node_count() == 3);
        REQUIRE(root.get_number_of_children() == 2);
        REQUIRE(root.get_number_of_parents() == 0);
        REQUIRE(root.get_height() == Approx(1.0));
        REQUIRE(root.get_length() == Approx(0.0));
        REQUIRE(root.is_child(&root_child1));
        REQUIRE(root.is_child(&root_child2));
        REQUIRE(! root.is_child(&leaf1));

        REQUIRE(! root_child1.is_root());
        REQUIRE(! root_child1.is_leaf());
        REQUIRE(root_child1.has_parent());
        REQUIRE(root_child1.has_children());
        REQUIRE(root_child1.is_dirty());
        REQUIRE(root_child1.get_node_count() == 4);
        REQUIRE(root_child1.get_leaf_node_count() == 3);
        REQUIRE(root_child1.get_internal_node_count() == 1);
        REQUIRE(root_child1.get_number_of_children() == 3);
        REQUIRE(root_child1.get_number_of_parents() == 1);
        REQUIRE(root_child1.get_height() == Approx(0.8));
        REQUIRE(root_child1.get_length() == Approx(0.2));
        REQUIRE(root_child1.is_child(&leaf1));
        REQUIRE(root_child1.is_child(&leaf2));
        REQUIRE(root_child1.is_child(&leaf3));
        REQUIRE(! root_child1.is_child(&leaf4));
        REQUIRE(! root_child1.is_child(&leaf5));
        REQUIRE(root_child1.is_parent(&root));
        REQUIRE(! root_child1.is_parent(&root_child2));

        REQUIRE(! root_child2.is_root());
        REQUIRE(! root_child2.is_leaf());
        REQUIRE(root_child2.has_parent());
        REQUIRE(root_child2.has_children());
        REQUIRE(root_child2.is_dirty());
        REQUIRE(root_child2.get_node_count() == 3);
        REQUIRE(root_child2.get_leaf_node_count() == 2);
        REQUIRE(root_child2.get_internal_node_count() == 1);
        REQUIRE(root_child2.get_number_of_children() == 2);
        REQUIRE(root_child2.get_number_of_parents() == 1);
        REQUIRE(root_child2.get_height() == Approx(0.3));
        REQUIRE(root_child2.get_length() == Approx(0.7));
        REQUIRE(! root_child2.is_child(&leaf1));
        REQUIRE(! root_child2.is_child(&leaf2));
        REQUIRE(! root_child2.is_child(&leaf3));
        REQUIRE(root_child2.is_child(&leaf4));
        REQUIRE(root_child2.is_child(&leaf5));
        REQUIRE(root_child2.is_parent(&root));
        REQUIRE(! root_child2.is_parent(&root_child1));

        REQUIRE(! leaf1.is_root());
        REQUIRE(leaf1.is_leaf());
        REQUIRE(leaf1.has_parent());
        REQUIRE(! leaf1.has_children());
        REQUIRE(leaf1.is_dirty());
        REQUIRE(leaf1.get_node_count() == 1);
        REQUIRE(leaf1.get_leaf_node_count() == 1);
        REQUIRE(leaf1.get_internal_node_count() == 0);
        REQUIRE(leaf1.get_number_of_children() == 0);
        REQUIRE(leaf1.get_number_of_parents() == 1);
        REQUIRE(leaf1.get_height() == Approx(0.0));
        REQUIRE(leaf1.get_length() == Approx(0.8));
        REQUIRE(leaf1.is_parent(&root_child1));

        REQUIRE(! leaf2.is_root());
        REQUIRE(leaf2.is_leaf());
        REQUIRE(leaf2.has_parent());
        REQUIRE(! leaf2.has_children());
        REQUIRE(leaf2.is_dirty());
        REQUIRE(leaf2.get_node_count() == 1);
        REQUIRE(leaf2.get_leaf_node_count() == 1);
        REQUIRE(leaf2.get_internal_node_count() == 0);
        REQUIRE(leaf2.get_number_of_children() == 0);
        REQUIRE(leaf2.get_number_of_parents() == 1);
        REQUIRE(leaf2.get_height() == Approx(0.0));
        REQUIRE(leaf2.get_length() == Approx(0.8));
        REQUIRE(leaf2.is_parent(&root_child1));

        REQUIRE(! leaf3.is_root());
        REQUIRE(leaf3.is_leaf());
        REQUIRE(leaf3.has_parent());
        REQUIRE(! leaf3.has_children());
        REQUIRE(leaf3.is_dirty());
        REQUIRE(leaf3.get_node_count() == 1);
        REQUIRE(leaf3.get_leaf_node_count() == 1);
        REQUIRE(leaf3.get_internal_node_count() == 0);
        REQUIRE(leaf3.get_number_of_children() == 0);
        REQUIRE(leaf3.get_number_of_parents() == 1);
        REQUIRE(leaf3.get_height() == Approx(0.0));
        REQUIRE(leaf3.get_length() == Approx(0.8));
        REQUIRE(leaf3.is_parent(&root_child1));

        REQUIRE(! leaf4.is_root());
        REQUIRE(leaf4.is_leaf());
        REQUIRE(leaf4.has_parent());
        REQUIRE(! leaf4.has_children());
        REQUIRE(leaf4.is_dirty());
        REQUIRE(leaf4.get_node_count() == 1);
        REQUIRE(leaf4.get_leaf_node_count() == 1);
        REQUIRE(leaf4.get_internal_node_count() == 0);
        REQUIRE(leaf4.get_number_of_children() == 0);
        REQUIRE(leaf4.get_number_of_parents() == 1);
        REQUIRE(leaf4.get_height() == Approx(0.0));
        REQUIRE(leaf4.get_length() == Approx(0.3));
        REQUIRE(leaf4.is_parent(&root_child2));

        REQUIRE(! leaf5.is_root());
        REQUIRE(leaf5.is_leaf());
        REQUIRE(leaf5.has_parent());
        REQUIRE(! leaf5.has_children());
        REQUIRE(leaf5.is_dirty());
        REQUIRE(leaf5.get_node_count() == 1);
        REQUIRE(leaf5.get_leaf_node_count() == 1);
        REQUIRE(leaf5.get_internal_node_count() == 0);
        REQUIRE(leaf5.get_number_of_children() == 0);
        REQUIRE(leaf5.get_number_of_parents() == 1);
        REQUIRE(leaf5.get_height() == Approx(0.0));
        REQUIRE(leaf5.get_length() == Approx(0.3));
        REQUIRE(leaf5.is_parent(&root_child2));
    }
}

TEST_CASE("Test simple tree cleaning Node", "[Node]") {

    SECTION("Testing tree cleaning") {
        Node root = Node("root", 1.0);
        Node root_child1 = Node("root child 1", 0.8);
        Node root_child2 = Node("root child 2", 0.3);
        Node leaf1 = Node("leaf 1");
        Node leaf2 = Node("leaf 2");
        Node leaf3 = Node("leaf 3");
        Node leaf4 = Node("leaf 4");
        Node leaf5 = Node("leaf 5");

        root.add_child(&root_child1);
        root.add_child(&root_child2);

        leaf1.add_parent(&root_child1);
        leaf2.add_parent(&root_child1);
        leaf3.add_parent(&root_child1);

        leaf4.add_parent(&root_child2);
        leaf5.add_parent(&root_child2);

        REQUIRE(root.is_dirty());
        REQUIRE(root_child1.is_dirty());
        REQUIRE(root_child2.is_dirty());
        REQUIRE(leaf1.is_dirty());
        REQUIRE(leaf2.is_dirty());
        REQUIRE(leaf3.is_dirty());
        REQUIRE(leaf4.is_dirty());
        REQUIRE(leaf5.is_dirty());
        REQUIRE(root.clade_has_dirt());
        REQUIRE(root_child1.clade_has_dirt());
        REQUIRE(root_child2.clade_has_dirt());
        REQUIRE(leaf1.clade_has_dirt());
        REQUIRE(leaf2.clade_has_dirt());
        REQUIRE(leaf3.clade_has_dirt());
        REQUIRE(leaf4.clade_has_dirt());
        REQUIRE(leaf5.clade_has_dirt());

        root.make_all_clean();

        REQUIRE(! root.is_dirty());
        REQUIRE(! root_child1.is_dirty());
        REQUIRE(! root_child2.is_dirty());
        REQUIRE(! leaf1.is_dirty());
        REQUIRE(! leaf2.is_dirty());
        REQUIRE(! leaf3.is_dirty());
        REQUIRE(! leaf4.is_dirty());
        REQUIRE(! leaf5.is_dirty());
        REQUIRE(! root.clade_has_dirt());
        REQUIRE(! root_child1.clade_has_dirt());
        REQUIRE(! root_child2.clade_has_dirt());
        REQUIRE(! leaf1.clade_has_dirt());
        REQUIRE(! leaf2.clade_has_dirt());
        REQUIRE(! leaf3.clade_has_dirt());
        REQUIRE(! leaf4.clade_has_dirt());
        REQUIRE(! leaf5.clade_has_dirt());

        root_child1.make_all_dirty();

        REQUIRE(! root.is_dirty());
        REQUIRE(root_child1.is_dirty());
        REQUIRE(! root_child2.is_dirty());
        REQUIRE(leaf1.is_dirty());
        REQUIRE(leaf2.is_dirty());
        REQUIRE(leaf3.is_dirty());
        REQUIRE(! leaf4.is_dirty());
        REQUIRE(! leaf5.is_dirty());
        REQUIRE(root.clade_has_dirt());
        REQUIRE(root_child1.clade_has_dirt());
        REQUIRE(! root_child2.clade_has_dirt());
        REQUIRE(leaf1.clade_has_dirt());
        REQUIRE(leaf2.clade_has_dirt());
        REQUIRE(leaf3.clade_has_dirt());
        REQUIRE(! leaf4.clade_has_dirt());
        REQUIRE(! leaf5.clade_has_dirt());

        root.make_all_clean();

        REQUIRE(! root.is_dirty());
        REQUIRE(! root_child1.is_dirty());
        REQUIRE(! root_child2.is_dirty());
        REQUIRE(! leaf1.is_dirty());
        REQUIRE(! leaf2.is_dirty());
        REQUIRE(! leaf3.is_dirty());
        REQUIRE(! leaf4.is_dirty());
        REQUIRE(! leaf5.is_dirty());
        REQUIRE(! root.clade_has_dirt());
        REQUIRE(! root_child1.clade_has_dirt());
        REQUIRE(! root_child2.clade_has_dirt());
        REQUIRE(! leaf1.clade_has_dirt());
        REQUIRE(! leaf2.clade_has_dirt());
        REQUIRE(! leaf3.clade_has_dirt());
        REQUIRE(! leaf4.clade_has_dirt());
        REQUIRE(! leaf5.clade_has_dirt());

        root_child2.make_all_dirty();

        REQUIRE(! root.is_dirty());
        REQUIRE(! root_child1.is_dirty());
        REQUIRE(root_child2.is_dirty());
        REQUIRE(! leaf1.is_dirty());
        REQUIRE(! leaf2.is_dirty());
        REQUIRE(! leaf3.is_dirty());
        REQUIRE(leaf4.is_dirty());
        REQUIRE(leaf5.is_dirty());
        REQUIRE(root.clade_has_dirt());
        REQUIRE(! root_child1.clade_has_dirt());
        REQUIRE(root_child2.clade_has_dirt());
        REQUIRE(! leaf1.clade_has_dirt());
        REQUIRE(! leaf2.clade_has_dirt());
        REQUIRE(! leaf3.clade_has_dirt());
        REQUIRE(leaf4.clade_has_dirt());
        REQUIRE(leaf5.clade_has_dirt());

        root_child1.make_all_clean();

        REQUIRE(! root.is_dirty());
        REQUIRE(! root_child1.is_dirty());
        REQUIRE(root_child2.is_dirty());
        REQUIRE(! leaf1.is_dirty());
        REQUIRE(! leaf2.is_dirty());
        REQUIRE(! leaf3.is_dirty());
        REQUIRE(leaf4.is_dirty());
        REQUIRE(leaf5.is_dirty());
        REQUIRE(root.clade_has_dirt());
        REQUIRE(! root_child1.clade_has_dirt());
        REQUIRE(root_child2.clade_has_dirt());
        REQUIRE(! leaf1.clade_has_dirt());
        REQUIRE(! leaf2.clade_has_dirt());
        REQUIRE(! leaf3.clade_has_dirt());
        REQUIRE(leaf4.clade_has_dirt());
        REQUIRE(leaf5.clade_has_dirt());

        root_child2.make_all_clean();

        REQUIRE(! root.is_dirty());
        REQUIRE(! root_child1.is_dirty());
        REQUIRE(! root_child2.is_dirty());
        REQUIRE(! leaf1.is_dirty());
        REQUIRE(! leaf2.is_dirty());
        REQUIRE(! leaf3.is_dirty());
        REQUIRE(! leaf4.is_dirty());
        REQUIRE(! leaf5.is_dirty());
        REQUIRE(! root.clade_has_dirt());
        REQUIRE(! root_child1.clade_has_dirt());
        REQUIRE(! root_child2.clade_has_dirt());
        REQUIRE(! leaf1.clade_has_dirt());
        REQUIRE(! leaf2.clade_has_dirt());
        REQUIRE(! leaf3.clade_has_dirt());
        REQUIRE(! leaf4.clade_has_dirt());
        REQUIRE(! leaf5.clade_has_dirt());

        leaf1.make_all_dirty();

        REQUIRE(! root.is_dirty());
        REQUIRE(! root_child1.is_dirty());
        REQUIRE(! root_child2.is_dirty());
        REQUIRE(leaf1.is_dirty());
        REQUIRE(! leaf2.is_dirty());
        REQUIRE(! leaf3.is_dirty());
        REQUIRE(! leaf4.is_dirty());
        REQUIRE(! leaf5.is_dirty());
        REQUIRE(root.clade_has_dirt());
        REQUIRE(root_child1.clade_has_dirt());
        REQUIRE(! root_child2.clade_has_dirt());
        REQUIRE(leaf1.clade_has_dirt());
        REQUIRE(! leaf2.clade_has_dirt());
        REQUIRE(! leaf3.clade_has_dirt());
        REQUIRE(! leaf4.clade_has_dirt());
        REQUIRE(! leaf5.clade_has_dirt());

        root.make_all_dirty();

        REQUIRE(root.is_dirty());
        REQUIRE(root_child1.is_dirty());
        REQUIRE(root_child2.is_dirty());
        REQUIRE(leaf1.is_dirty());
        REQUIRE(leaf2.is_dirty());
        REQUIRE(leaf3.is_dirty());
        REQUIRE(leaf4.is_dirty());
        REQUIRE(leaf5.is_dirty());
        REQUIRE(root.clade_has_dirt());
        REQUIRE(root_child1.clade_has_dirt());
        REQUIRE(root_child2.clade_has_dirt());
        REQUIRE(leaf1.clade_has_dirt());
        REQUIRE(leaf2.clade_has_dirt());
        REQUIRE(leaf3.clade_has_dirt());
        REQUIRE(leaf4.clade_has_dirt());
        REQUIRE(leaf5.clade_has_dirt());

        root.make_all_clean();

        REQUIRE(! root.is_dirty());
        REQUIRE(! root_child1.is_dirty());
        REQUIRE(! root_child2.is_dirty());
        REQUIRE(! leaf1.is_dirty());
        REQUIRE(! leaf2.is_dirty());
        REQUIRE(! leaf3.is_dirty());
        REQUIRE(! leaf4.is_dirty());
        REQUIRE(! leaf5.is_dirty());
        REQUIRE(! root.clade_has_dirt());
        REQUIRE(! root_child1.clade_has_dirt());
        REQUIRE(! root_child2.clade_has_dirt());
        REQUIRE(! leaf1.clade_has_dirt());
        REQUIRE(! leaf2.clade_has_dirt());
        REQUIRE(! leaf3.clade_has_dirt());
        REQUIRE(! leaf4.clade_has_dirt());
        REQUIRE(! leaf5.clade_has_dirt());

        root_child2.make_dirty();

        REQUIRE(! root.is_dirty());
        REQUIRE(! root_child1.is_dirty());
        REQUIRE(root_child2.is_dirty());
        REQUIRE(! leaf1.is_dirty());
        REQUIRE(! leaf2.is_dirty());
        REQUIRE(! leaf3.is_dirty());
        REQUIRE(! leaf4.is_dirty());
        REQUIRE(! leaf5.is_dirty());
        REQUIRE(root.clade_has_dirt());
        REQUIRE(! root_child1.clade_has_dirt());
        REQUIRE(root_child2.clade_has_dirt());
        REQUIRE(! leaf1.clade_has_dirt());
        REQUIRE(! leaf2.clade_has_dirt());
        REQUIRE(! leaf3.clade_has_dirt());
        REQUIRE(! leaf4.clade_has_dirt());
        REQUIRE(! leaf5.clade_has_dirt());
    }
}

TEST_CASE("Testing bare constructor of PopulationNode", "[PopulationNode]") {

    SECTION("Testing bare constructor") {
        PopulationNode n = PopulationNode();
        REQUIRE(n.get_height() == 0.0);
        REQUIRE(n.get_length() == 0.0);
        REQUIRE(n.get_label() == "");
        REQUIRE(n.get_allele_count() == 0);
        REQUIRE(n.is_dirty());

        REQUIRE_THROWS_AS(n.get_bottom_pattern_probability(1, 0), std::out_of_range);
        REQUIRE_THROWS_AS(n.get_top_pattern_probability(1, 0), std::out_of_range);

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

        PopulationNode * p = n.get_parent();
        REQUIRE(p == NULL);
        REQUIRE(typeid(&n).hash_code() == typeid(p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);

        delete p;
    }
}

TEST_CASE("Testing label constructor of PopulationNode", "[PopulationNode]") {

    SECTION("Testing label constructor") {
        PopulationNode n = PopulationNode("leaf1");
        REQUIRE(n.get_height() == 0.0);
        REQUIRE(n.get_length() == 0.0);
        REQUIRE(n.get_label() == "leaf1");
        REQUIRE(n.get_allele_count() == 0);
        REQUIRE(n.is_dirty());

        REQUIRE_THROWS_AS(n.get_bottom_pattern_probability(1, 0), std::out_of_range);
        REQUIRE_THROWS_AS(n.get_top_pattern_probability(1, 0), std::out_of_range);

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

        PopulationNode * p = n.get_parent();
        REQUIRE(p == NULL);
        REQUIRE(typeid(&n).hash_code() == typeid(p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);

        delete p;
    }
}

TEST_CASE("Testing height constructor of PopulationNode", "[PopulationNode]") {
    SECTION("Testing height constructor") {
        PopulationNode n = PopulationNode(0.03);
        REQUIRE(n.get_height() == Approx(0.03));
        REQUIRE(n.get_length() == 0.0);
        REQUIRE(n.get_label() == "");
        REQUIRE(n.get_allele_count() == 0);
        REQUIRE(n.is_dirty());

        REQUIRE_THROWS_AS(n.get_bottom_pattern_probability(1, 0), std::out_of_range);
        REQUIRE_THROWS_AS(n.get_top_pattern_probability(1, 0), std::out_of_range);

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

        PopulationNode * p = n.get_parent();
        REQUIRE(p == NULL);
        REQUIRE(typeid(&n).hash_code() == typeid(p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);

        delete p;
    }
}

TEST_CASE("Testing label and height constructor of PopulationNode", "[PopulationNode]") {
    SECTION("Testing label and height constructor") {
        PopulationNode n = PopulationNode("leaf1", 0.02);
        REQUIRE(n.get_height() == Approx(0.02));
        REQUIRE(n.get_length() == 0.0);
        REQUIRE(n.get_label() == "leaf1");
        REQUIRE(n.get_allele_count() == 0);
        REQUIRE(n.is_dirty());

        REQUIRE_THROWS_AS(n.get_bottom_pattern_probability(1, 0), std::out_of_range);
        REQUIRE_THROWS_AS(n.get_top_pattern_probability(1, 0), std::out_of_range);

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

        PopulationNode * p = n.get_parent();
        REQUIRE(p == NULL);
        REQUIRE(typeid(&n).hash_code() == typeid(p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);

        delete p;
    }
}

TEST_CASE("Testing allele count constructor of PopulationNode", "[PopulationNode]") {
    SECTION("Testing allele count constructor") {
        PopulationNode n = PopulationNode((unsigned int)3);
        REQUIRE(n.get_height() == Approx(0.0));
        REQUIRE(n.get_length() == 0.0);
        REQUIRE(n.get_label() == "");
        REQUIRE(n.get_allele_count() == 3);
        REQUIRE(n.is_dirty());

        REQUIRE(n.get_bottom_pattern_probability(2, 1) == Approx(0.0));
        REQUIRE_THROWS_AS(n.get_bottom_pattern_probability(4, 0), std::out_of_range);
        REQUIRE(n.get_top_pattern_probability(2, 1) == Approx(0.0));
        REQUIRE_THROWS_AS(n.get_top_pattern_probability(4, 0), std::out_of_range);

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

        PopulationNode * p = n.get_parent();
        REQUIRE(p == NULL);
        REQUIRE(typeid(&n).hash_code() == typeid(p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);

        delete p;
    }
}

TEST_CASE("Testing height and allele count constructor of PopulationNode", "[PopulationNode]") {
    SECTION("Testing height and allele count constructor") {
        PopulationNode n = PopulationNode(0.03, (unsigned int)3);
        REQUIRE(n.get_height() == Approx(0.03));
        REQUIRE(n.get_length() == 0.0);
        REQUIRE(n.get_label() == "");
        REQUIRE(n.get_allele_count() == 3);
        REQUIRE(n.is_dirty());

        REQUIRE(n.get_bottom_pattern_probability(2, 1) == Approx(0.0));
        REQUIRE_THROWS_AS(n.get_bottom_pattern_probability(4, 0), std::out_of_range);
        REQUIRE(n.get_top_pattern_probability(2, 1) == Approx(0.0));
        REQUIRE_THROWS_AS(n.get_top_pattern_probability(4, 0), std::out_of_range);

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

        PopulationNode * p = n.get_parent();
        REQUIRE(p == NULL);
        REQUIRE(typeid(&n).hash_code() == typeid(p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);

        delete p;
    }
}

TEST_CASE("Testing label and allele count constructor of PopulationNode", "[PopulationNode]") {
    SECTION("Testing label and allele count constructor") {
        PopulationNode n = PopulationNode("leaf1", (unsigned int)2);
        REQUIRE(n.get_height() == Approx(0.0));
        REQUIRE(n.get_length() == 0.0);
        REQUIRE(n.get_label() == "leaf1");
        REQUIRE(n.get_allele_count() == 2);
        REQUIRE(n.is_dirty());

        REQUIRE(n.get_bottom_pattern_probability(2, 1) == Approx(0.0));
        REQUIRE_THROWS_AS(n.get_bottom_pattern_probability(3, 0), std::out_of_range);
        REQUIRE(n.get_top_pattern_probability(2, 1) == Approx(0.0));
        REQUIRE_THROWS_AS(n.get_top_pattern_probability(3, 0), std::out_of_range);

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

        PopulationNode * p = n.get_parent();
        REQUIRE(p == NULL);
        REQUIRE(typeid(&n).hash_code() == typeid(p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);

        delete p;
    }
}

TEST_CASE("Testing label, height, and allele count constructor of PopulationNode", "[PopulationNode]") {
    SECTION("Testing label, height, and allele count constructor") {
        PopulationNode n = PopulationNode("leaf1", 0.03, (unsigned int)3);
        REQUIRE(n.get_height() == Approx(0.03));
        REQUIRE(n.get_length() == 0.0);
        REQUIRE(n.get_label() == "leaf1");
        REQUIRE(n.get_allele_count() == 3);
        REQUIRE(n.is_dirty());

        REQUIRE(n.get_bottom_pattern_probability(2, 1) == Approx(0.0));
        REQUIRE_THROWS_AS(n.get_bottom_pattern_probability(4, 0), std::out_of_range);
        REQUIRE(n.get_top_pattern_probability(2, 1) == Approx(0.0));
        REQUIRE_THROWS_AS(n.get_top_pattern_probability(4, 0), std::out_of_range);

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

        PopulationNode * p = n.get_parent();
        REQUIRE(p == NULL);
        REQUIRE(typeid(&n).hash_code() == typeid(p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);

        delete p;
    }
}


TEST_CASE("Testing node ref constructor of PopulationNode", "[PopulationNode]") {
    SECTION("Testing node ref constructor") {
        PopulationNode n2 = PopulationNode("leaf1", 0.03, (unsigned int)3);
        n2.set_bottom_pattern_probability(2, 1, 1.0);
        PopulationNode n = PopulationNode(n2);
        REQUIRE(typeid(n).hash_code() == typeid(n2).hash_code());
        REQUIRE(n.get_height() == Approx(0.03));
        REQUIRE(n.get_length() == 0.0);
        REQUIRE(n.get_label() == "leaf1");
        REQUIRE(n.get_allele_count() == 3);
        REQUIRE(n.get_bottom_pattern_probability(2, 1) == Approx(1.0));
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

        PopulationNode * p = n.get_parent();
        REQUIRE(p == NULL);
        REQUIRE(typeid(&n).hash_code() == typeid(p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
        
        n.set_label("leaf2");
        REQUIRE(n.get_label() == "leaf2");
        REQUIRE(n2.get_label() == "leaf1");

        n2.set_height(0.06);
        n2.reset(2);
        n2.set_bottom_pattern_probability(1, 1, 1.0);
        REQUIRE(n.get_height() == Approx(0.03));
        REQUIRE(n2.get_height() == Approx(0.06));
        REQUIRE(n.get_allele_count() == 3);
        REQUIRE(n.get_bottom_pattern_probability(2, 1) == Approx(1.0));
        REQUIRE(n.get_bottom_pattern_probability(1, 1) == Approx(0.0));
        REQUIRE(n2.get_allele_count() == 2);
        REQUIRE(n2.get_bottom_pattern_probability(2, 1) == Approx(0.0));
        REQUIRE(n2.get_bottom_pattern_probability(1, 1) == Approx(1.0));

        std::vector<double> e_bottom_n = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        std::vector<double> e_top_n = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        std::vector<double> e_bottom_n2 = {0.0, 1.0, 0.0, 0.0, 0.0};
        std::vector<double> e_top_n2 = {0.0, 0.0, 0.0, 0.0, 0.0};
        REQUIRE(n.get_bottom_pattern_probs().get_pattern_prob_matrix() == e_bottom_n);
        REQUIRE(n.get_top_pattern_probs().get_pattern_prob_matrix() == e_top_n);
        REQUIRE(n2.get_bottom_pattern_probs().get_pattern_prob_matrix() == e_bottom_n2);
        REQUIRE(n2.get_top_pattern_probs().get_pattern_prob_matrix() == e_top_n2);

        delete p;
    }
}

TEST_CASE("Testing copy operator of PopulationNode", "[PopulationNode]") {
    SECTION("Testing copy operator") {
        PopulationNode n2 = PopulationNode("leaf1", 0.03, (unsigned int)3);
        n2.set_bottom_pattern_probability(2, 1, 1.0);
        PopulationNode n = n2;
        REQUIRE(typeid(n).hash_code() == typeid(n2).hash_code());
        REQUIRE(n.get_height() == Approx(0.03));
        REQUIRE(n.get_length() == 0.0);
        REQUIRE(n.get_label() == "leaf1");
        REQUIRE(n.get_allele_count() == 3);
        REQUIRE(n.get_bottom_pattern_probability(2, 1) == Approx(1.0));
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

        PopulationNode * p = n.get_parent();
        REQUIRE(p == NULL);
        REQUIRE(typeid(&n).hash_code() == typeid(p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
        
        n.set_label("leaf2");
        REQUIRE(n.get_label() == "leaf2");
        REQUIRE(n2.get_label() == "leaf1");

        n2.set_height(0.06);
        n2.reset(2);
        n2.set_bottom_pattern_probability(1, 1, 1.0);
        REQUIRE(n.get_height() == Approx(0.03));
        REQUIRE(n2.get_height() == Approx(0.06));
        REQUIRE(n.get_allele_count() == 3);
        REQUIRE(n.get_bottom_pattern_probability(2, 1) == Approx(1.0));
        REQUIRE(n.get_bottom_pattern_probability(1, 1) == Approx(0.0));
        REQUIRE(n2.get_allele_count() == 2);
        REQUIRE(n2.get_bottom_pattern_probability(2, 1) == Approx(0.0));
        REQUIRE(n2.get_bottom_pattern_probability(1, 1) == Approx(1.0));

        std::vector<double> e_bottom_n = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        std::vector<double> e_top_n = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        std::vector<double> e_bottom_n2 = {0.0, 1.0, 0.0, 0.0, 0.0};
        std::vector<double> e_top_n2 = {0.0, 0.0, 0.0, 0.0, 0.0};
        REQUIRE(n.get_bottom_pattern_probs().get_pattern_prob_matrix() == e_bottom_n);
        REQUIRE(n.get_top_pattern_probs().get_pattern_prob_matrix() == e_top_n);
        REQUIRE(n2.get_bottom_pattern_probs().get_pattern_prob_matrix() == e_bottom_n2);
        REQUIRE(n2.get_top_pattern_probs().get_pattern_prob_matrix() == e_top_n2);

        delete p;
    }
}

TEST_CASE("Testing clone method of PopulationNode", "[PopulationNode]") {

    SECTION("Testing clone") {
        PopulationNode * n2 = new PopulationNode("leaf1", 0.03, (unsigned int)3);
        n2->set_bottom_pattern_probability(2, 1, 1.0);
        PopulationNode * n = n2->clone();
        REQUIRE(typeid(n).hash_code() == typeid(n2).hash_code());
        REQUIRE(n->get_height() == Approx(0.03));
        REQUIRE(n->get_length() == 0.0);
        REQUIRE(n->get_label() == "leaf1");
        REQUIRE(n->get_allele_count() == 3);
        REQUIRE(n->get_bottom_pattern_probability(2, 1) == Approx(1.0));
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

        PopulationNode * p = n->get_parent();
        REQUIRE(p == NULL);
        REQUIRE(typeid(n).hash_code() == typeid(p).hash_code());

        REQUIRE_THROWS_AS(n->get_child(0), std::out_of_range);
        
        n->set_label("leaf2");
        REQUIRE(n->get_label() == "leaf2");
        REQUIRE(n2->get_label() == "leaf1");

        n2->set_height(0.06);
        n2->reset(2);
        n2->set_bottom_pattern_probability(1, 1, 1.0);
        REQUIRE(n->get_height() == Approx(0.03));
        REQUIRE(n2->get_height() == Approx(0.06));
        REQUIRE(n->get_allele_count() == 3);
        REQUIRE(n->get_bottom_pattern_probability(2, 1) == Approx(1.0));
        REQUIRE(n->get_bottom_pattern_probability(1, 1) == Approx(0.0));
        REQUIRE(n2->get_allele_count() == 2);
        REQUIRE(n2->get_bottom_pattern_probability(2, 1) == Approx(0.0));
        REQUIRE(n2->get_bottom_pattern_probability(1, 1) == Approx(1.0));

        std::vector<double> e_bottom_n = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        std::vector<double> e_top_n = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        std::vector<double> e_bottom_n2 = {0.0, 1.0, 0.0, 0.0, 0.0};
        std::vector<double> e_top_n2 = {0.0, 0.0, 0.0, 0.0, 0.0};
        REQUIRE(n->get_bottom_pattern_probs().get_pattern_prob_matrix() == e_bottom_n);
        REQUIRE(n->get_top_pattern_probs().get_pattern_prob_matrix() == e_top_n);
        REQUIRE(n2->get_bottom_pattern_probs().get_pattern_prob_matrix() == e_bottom_n2);
        REQUIRE(n2->get_top_pattern_probs().get_pattern_prob_matrix() == e_top_n2);

        delete p;
        delete n;
        delete n2;
    }
}

TEST_CASE("Testing parent methods of PopulationNode", "[PopulationNode]") {

    SECTION("Testing parent methods") {
        PopulationNode c = PopulationNode();
        PopulationNode p = PopulationNode();

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

        PopulationNode * r = c.remove_parent();
        REQUIRE(r != NULL);
        REQUIRE(r == &p);
        REQUIRE(typeid(r).hash_code() == typeid(&c).hash_code());
        REQUIRE(c.has_parent() == false);
        REQUIRE(c.is_parent(&p) == false);
        REQUIRE(p.is_child(&c) == false);
        REQUIRE(c.get_number_of_parents() == 0);
        REQUIRE(c.degree() == 0);
        REQUIRE(p.degree() == 0);
    }
}

TEST_CASE("Testing child methods of PopulationNode", "[PopulationNode]") {

    SECTION("Testing child methods") {
        PopulationNode c = PopulationNode();
        PopulationNode p = PopulationNode();
        REQUIRE(typeid(c).hash_code() == typeid(p).hash_code());

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
        PopulationNode * r = p.remove_child(i);
        REQUIRE(r == &c);
        REQUIRE(typeid(&c).hash_code() == typeid(r).hash_code());
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

TEST_CASE("Test simple tree building with PopulationNode", "[PopulationNode]") {

    SECTION("Testing tree building") {
        PopulationNode root = PopulationNode("root", 1.0);
        PopulationNode root_child1 = PopulationNode("root child 1", 0.8);
        PopulationNode root_child2 = PopulationNode("root child 2", 0.3);
        PopulationNode leaf1 = PopulationNode("leaf 1", (unsigned int)3);
        PopulationNode leaf2 = PopulationNode("leaf 2", (unsigned int)3);
        PopulationNode leaf3 = PopulationNode("leaf 3", (unsigned int)3);
        PopulationNode leaf4 = PopulationNode("leaf 4", (unsigned int)3);
        PopulationNode leaf5 = PopulationNode("leaf 5", (unsigned int)3);

        root.add_child(&root_child1);
        root.add_child(&root_child2);

        leaf1.add_parent(&root_child1);
        leaf2.add_parent(&root_child1);
        leaf3.add_parent(&root_child1);

        leaf4.add_parent(&root_child2);
        leaf5.add_parent(&root_child2);

        REQUIRE(root.degree() == 2);
        REQUIRE(root_child1.degree() == 4);
        REQUIRE(root_child1.get_number_of_children() == 3);
        REQUIRE(root_child2.degree() == 3);
        REQUIRE(leaf1.degree() == 1);
        REQUIRE(leaf2.degree() == 1);
        REQUIRE(leaf3.degree() == 1);
        REQUIRE(leaf4.degree() == 1);
        REQUIRE(leaf5.degree() == 1);

        REQUIRE(root.is_root());
        REQUIRE(! root.is_leaf());
        REQUIRE(root.has_children());
        REQUIRE(! root.has_parent());
        REQUIRE(root.is_dirty());
        REQUIRE(root.get_node_count() == 8);
        REQUIRE(root.get_leaf_node_count() == 5);
        REQUIRE(root.get_internal_node_count() == 3);
        REQUIRE(root.get_number_of_children() == 2);
        REQUIRE(root.get_number_of_parents() == 0);
        REQUIRE(root.get_height() == Approx(1.0));
        REQUIRE(root.get_length() == Approx(0.0));
        REQUIRE(root.is_child(&root_child1));
        REQUIRE(root.is_child(&root_child2));
        REQUIRE(! root.is_child(&leaf1));

        REQUIRE(! root_child1.is_root());
        REQUIRE(! root_child1.is_leaf());
        REQUIRE(root_child1.has_parent());
        REQUIRE(root_child1.has_children());
        REQUIRE(root_child1.is_dirty());
        REQUIRE(root_child1.get_node_count() == 4);
        REQUIRE(root_child1.get_leaf_node_count() == 3);
        REQUIRE(root_child1.get_internal_node_count() == 1);
        REQUIRE(root_child1.get_number_of_children() == 3);
        REQUIRE(root_child1.get_number_of_parents() == 1);
        REQUIRE(root_child1.get_height() == Approx(0.8));
        REQUIRE(root_child1.get_length() == Approx(0.2));
        REQUIRE(root_child1.is_child(&leaf1));
        REQUIRE(root_child1.is_child(&leaf2));
        REQUIRE(root_child1.is_child(&leaf3));
        REQUIRE(! root_child1.is_child(&leaf4));
        REQUIRE(! root_child1.is_child(&leaf5));
        REQUIRE(root_child1.is_parent(&root));
        REQUIRE(! root_child1.is_parent(&root_child2));

        REQUIRE(! root_child2.is_root());
        REQUIRE(! root_child2.is_leaf());
        REQUIRE(root_child2.has_parent());
        REQUIRE(root_child2.has_children());
        REQUIRE(root_child2.is_dirty());
        REQUIRE(root_child2.get_node_count() == 3);
        REQUIRE(root_child2.get_leaf_node_count() == 2);
        REQUIRE(root_child2.get_internal_node_count() == 1);
        REQUIRE(root_child2.get_number_of_children() == 2);
        REQUIRE(root_child2.get_number_of_parents() == 1);
        REQUIRE(root_child2.get_height() == Approx(0.3));
        REQUIRE(root_child2.get_length() == Approx(0.7));
        REQUIRE(! root_child2.is_child(&leaf1));
        REQUIRE(! root_child2.is_child(&leaf2));
        REQUIRE(! root_child2.is_child(&leaf3));
        REQUIRE(root_child2.is_child(&leaf4));
        REQUIRE(root_child2.is_child(&leaf5));
        REQUIRE(root_child2.is_parent(&root));
        REQUIRE(! root_child2.is_parent(&root_child1));

        REQUIRE(! leaf1.is_root());
        REQUIRE(leaf1.is_leaf());
        REQUIRE(leaf1.has_parent());
        REQUIRE(! leaf1.has_children());
        REQUIRE(leaf1.is_dirty());
        REQUIRE(leaf1.get_node_count() == 1);
        REQUIRE(leaf1.get_leaf_node_count() == 1);
        REQUIRE(leaf1.get_internal_node_count() == 0);
        REQUIRE(leaf1.get_number_of_children() == 0);
        REQUIRE(leaf1.get_number_of_parents() == 1);
        REQUIRE(leaf1.get_height() == Approx(0.0));
        REQUIRE(leaf1.get_length() == Approx(0.8));
        REQUIRE(leaf1.is_parent(&root_child1));

        REQUIRE(! leaf2.is_root());
        REQUIRE(leaf2.is_leaf());
        REQUIRE(leaf2.has_parent());
        REQUIRE(! leaf2.has_children());
        REQUIRE(leaf2.is_dirty());
        REQUIRE(leaf2.get_node_count() == 1);
        REQUIRE(leaf2.get_leaf_node_count() == 1);
        REQUIRE(leaf2.get_internal_node_count() == 0);
        REQUIRE(leaf2.get_number_of_children() == 0);
        REQUIRE(leaf2.get_number_of_parents() == 1);
        REQUIRE(leaf2.get_height() == Approx(0.0));
        REQUIRE(leaf2.get_length() == Approx(0.8));
        REQUIRE(leaf2.is_parent(&root_child1));

        REQUIRE(! leaf3.is_root());
        REQUIRE(leaf3.is_leaf());
        REQUIRE(leaf3.has_parent());
        REQUIRE(! leaf3.has_children());
        REQUIRE(leaf3.is_dirty());
        REQUIRE(leaf3.get_node_count() == 1);
        REQUIRE(leaf3.get_leaf_node_count() == 1);
        REQUIRE(leaf3.get_internal_node_count() == 0);
        REQUIRE(leaf3.get_number_of_children() == 0);
        REQUIRE(leaf3.get_number_of_parents() == 1);
        REQUIRE(leaf3.get_height() == Approx(0.0));
        REQUIRE(leaf3.get_length() == Approx(0.8));
        REQUIRE(leaf3.is_parent(&root_child1));

        REQUIRE(! leaf4.is_root());
        REQUIRE(leaf4.is_leaf());
        REQUIRE(leaf4.has_parent());
        REQUIRE(! leaf4.has_children());
        REQUIRE(leaf4.is_dirty());
        REQUIRE(leaf4.get_node_count() == 1);
        REQUIRE(leaf4.get_leaf_node_count() == 1);
        REQUIRE(leaf4.get_internal_node_count() == 0);
        REQUIRE(leaf4.get_number_of_children() == 0);
        REQUIRE(leaf4.get_number_of_parents() == 1);
        REQUIRE(leaf4.get_height() == Approx(0.0));
        REQUIRE(leaf4.get_length() == Approx(0.3));
        REQUIRE(leaf4.is_parent(&root_child2));

        REQUIRE(! leaf5.is_root());
        REQUIRE(leaf5.is_leaf());
        REQUIRE(leaf5.has_parent());
        REQUIRE(! leaf5.has_children());
        REQUIRE(leaf5.is_dirty());
        REQUIRE(leaf5.get_node_count() == 1);
        REQUIRE(leaf5.get_leaf_node_count() == 1);
        REQUIRE(leaf5.get_internal_node_count() == 0);
        REQUIRE(leaf5.get_number_of_children() == 0);
        REQUIRE(leaf5.get_number_of_parents() == 1);
        REQUIRE(leaf5.get_height() == Approx(0.0));
        REQUIRE(leaf5.get_length() == Approx(0.3));
        REQUIRE(leaf5.is_parent(&root_child2));

        REQUIRE(leaf1.get_allele_count() == 3);
        REQUIRE(leaf2.get_allele_count() == 3);
        REQUIRE(leaf3.get_allele_count() == 3);
        REQUIRE(leaf4.get_allele_count() == 3);
        REQUIRE(leaf5.get_allele_count() == 3);
        REQUIRE(root_child1.get_allele_count() == 0);
        REQUIRE(root_child2.get_allele_count() == 0);
        REQUIRE(root.get_allele_count() == 0);
        
        root.resize_all();
        REQUIRE(leaf1.get_allele_count() == 3);
        REQUIRE(leaf2.get_allele_count() == 3);
        REQUIRE(leaf3.get_allele_count() == 3);
        REQUIRE(leaf4.get_allele_count() == 3);
        REQUIRE(leaf5.get_allele_count() == 3);
        REQUIRE(root_child1.get_allele_count() == 9);
        REQUIRE(root_child2.get_allele_count() == 6);
        REQUIRE(root.get_allele_count() == 15);
    }
}

TEST_CASE("Test node height prior", "[Node]") {

    SECTION("Testing prior") {
        Node root = Node("root", 1.0);
        Node root_child1 = Node("root child 1", 0.8);
        Node root_child2 = Node("root child 2", 0.3);
        Node leaf1 = Node("leaf 1", 0.0);
        Node leaf2 = Node("leaf 2", 0.0);
        Node leaf3 = Node("leaf 3", 0.0);
        Node leaf4 = Node("leaf 4", 0.0);
        Node leaf5 = Node("leaf 5", 0.0);

        root.add_child(&root_child1);
        root.add_child(&root_child2);

        leaf1.add_parent(&root_child1);
        leaf2.add_parent(&root_child1);
        leaf3.add_parent(&root_child1);

        leaf4.add_parent(&root_child2);
        leaf5.add_parent(&root_child2);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(1.0, 1.0);
        root.set_all_node_height_priors(prior);
        REQUIRE(root.calculate_ln_relative_node_height_prior_density() == Approx(-2.1));

        leaf1.fix_node_height();
        leaf2.fix_node_height();
        leaf3.fix_node_height();
        leaf4.fix_node_height();
        leaf5.fix_node_height();

        REQUIRE(root.calculate_ln_relative_node_height_prior_density() == Approx(-2.1));


        std::shared_ptr<GammaDistribution> prior2 = std::make_shared<GammaDistribution>(1.0, 0.01);
        root.set_all_node_height_priors(prior2);
        REQUIRE(root.calculate_ln_relative_node_height_prior_density() == Approx(-196.18448944203573));

        leaf1.estimate_node_height();
        leaf2.estimate_node_height();
        leaf3.estimate_node_height();
        leaf4.estimate_node_height();
        leaf5.estimate_node_height();

        REQUIRE(root.calculate_ln_relative_node_height_prior_density() == Approx(-173.15863851209528));
    }
}

TEST_CASE("Test node height and coalescence rate priors", "[PopulationNode]") {

    SECTION("Testing prior") {
        PopulationNode root = PopulationNode("root", 1.0);
        PopulationNode root_child1 = PopulationNode("root child 1", 0.8);
        PopulationNode root_child2 = PopulationNode("root child 2", 0.3);
        PopulationNode leaf1 = PopulationNode("leaf 1", 0.0);
        PopulationNode leaf2 = PopulationNode("leaf 2", 0.0);
        PopulationNode leaf3 = PopulationNode("leaf 3", 0.0);
        PopulationNode leaf4 = PopulationNode("leaf 4", 0.0);
        PopulationNode leaf5 = PopulationNode("leaf 5", 0.0);

        root.add_child(&root_child1);
        root.add_child(&root_child2);

        leaf1.add_parent(&root_child1);
        leaf2.add_parent(&root_child1);
        leaf3.add_parent(&root_child1);

        leaf4.add_parent(&root_child2);
        leaf5.add_parent(&root_child2);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(1.0, 1.0);
        root.set_all_node_height_priors(prior);
        root.set_all_population_size_priors(prior);
        REQUIRE(root.calculate_ln_relative_node_height_prior_density() == Approx(-2.1));
        REQUIRE(root.calculate_ln_relative_coalescence_rate_prior_density() == Approx(-0.2*8));

        leaf1.fix_node_height();
        leaf2.fix_node_height();
        leaf3.fix_node_height();
        leaf4.fix_node_height();
        leaf5.fix_node_height();

        REQUIRE(root.calculate_ln_relative_node_height_prior_density() == Approx(-2.1));
        REQUIRE(root.calculate_ln_relative_coalescence_rate_prior_density() == Approx(-0.2*8));


        std::shared_ptr<GammaDistribution> prior2 = std::make_shared<GammaDistribution>(1.0, 0.01);
        root.set_all_node_height_priors(prior2);
        root.set_all_population_size_priors(prior2);
        REQUIRE(root.calculate_ln_relative_node_height_prior_density() == Approx(-196.18448944203573));
        REQUIRE(root.calculate_ln_relative_coalescence_rate_prior_density() == Approx(-15.39482981401191*8));

        leaf1.estimate_node_height();
        leaf2.estimate_node_height();
        leaf3.estimate_node_height();
        leaf4.estimate_node_height();
        leaf5.estimate_node_height();

        REQUIRE(root.calculate_ln_relative_node_height_prior_density() == Approx(-173.15863851209528));
        REQUIRE(root.calculate_ln_relative_coalescence_rate_prior_density() == Approx(-15.39482981401191*8));

        root.set_all_coalescence_rates(100.0);

        REQUIRE(root.calculate_ln_relative_node_height_prior_density() == Approx(-173.15863851209528));
        REQUIRE(root.calculate_ln_relative_coalescence_rate_prior_density() == Approx(2.6051701859880909*8));

        root.set_all_coalescence_rate_parameters();
        root_child1.set_height_parameter(root_child2.get_height_parameter());

        REQUIRE(root.calculate_ln_relative_node_height_prior_density() == Approx(-97.76380869808338));
        REQUIRE(root.calculate_ln_relative_coalescence_rate_prior_density() == Approx(2.6051701859880909));

        leaf1.fix_node_height();
        leaf2.fix_node_height();
        leaf3.fix_node_height();
        leaf4.fix_node_height();
        leaf5.fix_node_height();

        REQUIRE(root.calculate_ln_relative_node_height_prior_density() == Approx(-120.78965962802383));
        REQUIRE(root.calculate_ln_relative_coalescence_rate_prior_density() == Approx(2.6051701859880909));
    }
}
