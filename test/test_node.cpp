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

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);

        delete p;
    }

    SECTION("Testing node ref constructor") {
        Node n2 = Node("leaf1", 0.02);
        Node n = Node(n2);
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

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
        
        n.set_label("leaf2");
        REQUIRE(n.get_label() == "leaf2");
        REQUIRE(n2.get_label() == "leaf1");

        n2.set_height(0.06);
        REQUIRE(n.get_height() == Approx(0.02));
        REQUIRE(n2.get_height() == Approx(0.06));
    }

}

TEST_CASE("Testing copy operator of Node", "[Node]") {

    SECTION("Testing copy operator") {
        Node n2 = Node("leaf1", 0.02);
        Node n = n2;
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

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
        
        n.set_label("leaf2");
        REQUIRE(n.get_label() == "leaf2");
        REQUIRE(n2.get_label() == "leaf1");

        n2.set_height(0.06);
        REQUIRE(n.get_height() == Approx(0.02));
        REQUIRE(n2.get_height() == Approx(0.06));
    }

}

TEST_CASE("Testing clone method of Node", "[Node]") {

    SECTION("Testing clone") {
        Node * n2 = new Node("leaf1", 0.02);
        Node * n = n2->clone();
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

        REQUIRE_THROWS_AS(n->get_child(0), std::out_of_range);
        
        n->set_label("leaf2");
        REQUIRE(n->get_label() == "leaf2");
        REQUIRE(n2->get_label() == "leaf1");

        n2->set_height(0.06);
        REQUIRE(n->get_height() == Approx(0.02));
        REQUIRE(n2->get_height() == Approx(0.06));

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

        root.make_all_clean();

        REQUIRE(! root.is_dirty());
        REQUIRE(! root_child1.is_dirty());
        REQUIRE(! root_child2.is_dirty());
        REQUIRE(! leaf1.is_dirty());
        REQUIRE(! leaf2.is_dirty());
        REQUIRE(! leaf3.is_dirty());
        REQUIRE(! leaf4.is_dirty());
        REQUIRE(! leaf5.is_dirty());

        root_child1.make_all_dirty();

        REQUIRE(! root.is_dirty());
        REQUIRE(root_child1.is_dirty());
        REQUIRE(! root_child2.is_dirty());
        REQUIRE(leaf1.is_dirty());
        REQUIRE(leaf2.is_dirty());
        REQUIRE(leaf3.is_dirty());
        REQUIRE(! leaf4.is_dirty());
        REQUIRE(! leaf5.is_dirty());

        root.make_all_clean();

        REQUIRE(! root.is_dirty());
        REQUIRE(! root_child1.is_dirty());
        REQUIRE(! root_child2.is_dirty());
        REQUIRE(! leaf1.is_dirty());
        REQUIRE(! leaf2.is_dirty());
        REQUIRE(! leaf3.is_dirty());
        REQUIRE(! leaf4.is_dirty());
        REQUIRE(! leaf5.is_dirty());

        root_child2.make_all_dirty();

        REQUIRE(! root.is_dirty());
        REQUIRE(! root_child1.is_dirty());
        REQUIRE(root_child2.is_dirty());
        REQUIRE(! leaf1.is_dirty());
        REQUIRE(! leaf2.is_dirty());
        REQUIRE(! leaf3.is_dirty());
        REQUIRE(leaf4.is_dirty());
        REQUIRE(leaf5.is_dirty());

        root_child1.make_all_clean();

        REQUIRE(! root.is_dirty());
        REQUIRE(! root_child1.is_dirty());
        REQUIRE(root_child2.is_dirty());
        REQUIRE(! leaf1.is_dirty());
        REQUIRE(! leaf2.is_dirty());
        REQUIRE(! leaf3.is_dirty());
        REQUIRE(leaf4.is_dirty());
        REQUIRE(leaf5.is_dirty());

        root_child2.make_all_clean();

        REQUIRE(! root.is_dirty());
        REQUIRE(! root_child1.is_dirty());
        REQUIRE(! root_child2.is_dirty());
        REQUIRE(! leaf1.is_dirty());
        REQUIRE(! leaf2.is_dirty());
        REQUIRE(! leaf3.is_dirty());
        REQUIRE(! leaf4.is_dirty());
        REQUIRE(! leaf5.is_dirty());

        leaf1.make_all_dirty();

        REQUIRE(! root.is_dirty());
        REQUIRE(! root_child1.is_dirty());
        REQUIRE(! root_child2.is_dirty());
        REQUIRE(leaf1.is_dirty());
        REQUIRE(! leaf2.is_dirty());
        REQUIRE(! leaf3.is_dirty());
        REQUIRE(! leaf4.is_dirty());
        REQUIRE(! leaf5.is_dirty());

        root.make_all_dirty();

        REQUIRE(root.is_dirty());
        REQUIRE(root_child1.is_dirty());
        REQUIRE(root_child2.is_dirty());
        REQUIRE(leaf1.is_dirty());
        REQUIRE(leaf2.is_dirty());
        REQUIRE(leaf3.is_dirty());
        REQUIRE(leaf4.is_dirty());
        REQUIRE(leaf5.is_dirty());
    }
}
