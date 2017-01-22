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

        std::shared_ptr<Node> p = n.get_parent();
        REQUIRE(p == nullptr);
        REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
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

        std::shared_ptr<Node> p = n.get_parent();
        REQUIRE(p == nullptr);
        REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
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

        std::shared_ptr<Node> p = n.get_parent();
        REQUIRE(p == nullptr);
        REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
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

        std::shared_ptr<Node> p = n.get_parent();
        REQUIRE(p == nullptr);
        REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
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

        std::shared_ptr<Node> p = n.get_parent();
        REQUIRE(p == nullptr);
        REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
        
        n.set_label("leaf2");
        REQUIRE(n.get_label() == "leaf2");
        REQUIRE(n2.get_label() == "leaf1");

        n2.set_height(0.06);
        REQUIRE(n.get_height() == Approx(0.06));
        REQUIRE(n2.get_height() == Approx(0.06));
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

        std::shared_ptr<Node> p = n.get_parent();
        REQUIRE(p == nullptr);
        REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
        
        n.set_label("leaf2");
        REQUIRE(n.get_label() == "leaf2");
        REQUIRE(n2.get_label() == "leaf1");

        n2.set_height(0.06);
        REQUIRE(n.get_height() == Approx(0.06));
        REQUIRE(n2.get_height() == Approx(0.06));
    }

}

/* TEST_CASE("Testing clone method of Node", "[Node]") { */

/*     SECTION("Testing clone") { */
/*         Node * n2 = new Node("leaf1", 0.02); */
/*         Node * n = n2->clone(); */
/*         REQUIRE(typeid(n).hash_code() == typeid(n2).hash_code()); */
/*         REQUIRE(n->get_height() == Approx(0.02)); */
/*         REQUIRE(n->get_length() == 0.0); */
/*         REQUIRE(n->get_label() == "leaf1"); */
/*         REQUIRE(n->is_dirty()); */

/*         REQUIRE(n->degree() == 0); */
/*         REQUIRE(n->has_parent() == false); */
/*         REQUIRE(n->get_number_of_parents() == 0); */
/*         REQUIRE(n->has_children() == false); */
/*         REQUIRE(n->get_number_of_children() == 0); */

/*         REQUIRE(n->is_leaf() == true); */
/*         REQUIRE(n->is_root() == true); */
/*         REQUIRE(n->get_node_count() == 1); */
/*         REQUIRE(n->get_leaf_node_count() == 1); */
/*         REQUIRE(n->get_internal_node_count() == 0); */

/*         std::shared_ptr<Node> p = n->get_parent(); */
/*         REQUIRE(p == nullptr); */
/*         REQUIRE(typeid(n).hash_code() == typeid(p).hash_code()); */

/*         REQUIRE_THROWS_AS(n->get_child(0), std::out_of_range); */
        
/*         n->set_label("leaf2"); */
/*         REQUIRE(n->get_label() == "leaf2"); */
/*         REQUIRE(n2->get_label() == "leaf1"); */

/*         n2->set_height(0.06); */
/*         REQUIRE(n->get_height() == Approx(0.02)); */
/*         REQUIRE(n2->get_height() == Approx(0.06)); */
/*         delete n; */
/*         delete n2; */
/*     } */

/* } */

TEST_CASE("Testing parent methods of Node", "[Node]") {

    SECTION("Testing parent methods") {
        std::shared_ptr<Node> c = std::make_shared<Node>();
        std::shared_ptr<Node> p = std::make_shared<Node>();

        REQUIRE(c->degree() == 0);
        REQUIRE(p->degree() == 0);
        REQUIRE(c->has_parent() == false);
        REQUIRE(c->get_number_of_parents() == 0);

        REQUIRE(c->get_parent() == nullptr);
        REQUIRE_THROWS_AS(p->get_child(0), std::out_of_range);
        REQUIRE(c->is_parent(p) == false);
        REQUIRE(p->is_child(c) == false);

        c->add_parent(p);
        REQUIRE(c->has_parent() == true);
        REQUIRE(c->is_parent(p) == true);
        REQUIRE(p->is_child(c) == true);
        REQUIRE(c->get_number_of_parents() == 1);
        REQUIRE(c->degree() == 1);
        REQUIRE(p->degree() == 1);

        std::shared_ptr<Node> r = c->remove_parent();
        REQUIRE(r != nullptr);
        REQUIRE(r == p);
        REQUIRE(typeid(r).hash_code() == typeid(c).hash_code());
        REQUIRE(c->has_parent() == false);
        REQUIRE(c->is_parent(p) == false);
        REQUIRE(p->is_child(c) == false);
        REQUIRE(c->get_number_of_parents() == 0);
        REQUIRE(c->degree() == 0);
        REQUIRE(p->degree() == 0);
    }
}

TEST_CASE("Testing child methods of Node", "[Node]") {

    SECTION("Testing child methods") {
        std::shared_ptr<Node> c = std::make_shared<Node>();
        std::shared_ptr<Node> p = std::make_shared<Node>();
        REQUIRE(typeid(c).hash_code() == typeid(p).hash_code());

        REQUIRE(c->degree() == 0);
        REQUIRE(p->degree() == 0);
        REQUIRE(p->has_children() == false);
        REQUIRE(p->get_number_of_children() == 0);

        REQUIRE(c->get_parent() == nullptr);
        REQUIRE_THROWS_AS(p->get_child(0), std::out_of_range);
        REQUIRE(c->is_parent(p) == false);
        REQUIRE(p->is_child(c) == false);

        p->add_child(c);
        REQUIRE(p->has_children() == true);
        REQUIRE(c->is_parent(p) == true);
        REQUIRE(p->is_child(c) == true);
        REQUIRE(p->get_number_of_children() == 1);
        REQUIRE(c->degree() == 1);
        REQUIRE(p->degree() == 1);

        unsigned int i = 1;
        REQUIRE_THROWS_AS(p->remove_child(i), std::out_of_range);
        i = 0;
        std::shared_ptr<Node> r = p->remove_child(i);
        REQUIRE(r == c);
        REQUIRE(typeid(c).hash_code() == typeid(r).hash_code());
        REQUIRE(typeid(*c).hash_code() == typeid(*r).hash_code());
        REQUIRE(p->has_children() == false);
        REQUIRE(c->is_parent(p) == false);
        REQUIRE(p->is_child(c) == false);
        REQUIRE(p->get_number_of_children() == 0);
        REQUIRE(c->degree() == 0);
        REQUIRE(p->degree() == 0);

        p->add_child(c);
        REQUIRE(p->has_children() == true);
        REQUIRE(c->is_parent(p) == true);
        REQUIRE(p->is_child(c) == true);
        REQUIRE(p->get_number_of_children() == 1);
        REQUIRE(c->degree() == 1);
        REQUIRE(p->degree() == 1);

        p->remove_child(c);
        REQUIRE(p->has_children() == false);
        REQUIRE(c->is_parent(p) == false);
        REQUIRE(p->is_child(c) == false);
        REQUIRE(p->get_number_of_children() == 0);
        REQUIRE(c->degree() == 0);
        REQUIRE(p->degree() == 0);
        
    }
}

TEST_CASE("Test simple tree building with Node", "[Node]") {

    SECTION("Testing tree building") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> root_child1 = std::make_shared<Node>("root child 1", 0.8);
        std::shared_ptr<Node> root_child2 = std::make_shared<Node>("root child 2", 0.3);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf 1");
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf 2");
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf 3");
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf 4");
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf 5");

        root->add_child(root_child1);
        root->add_child(root_child2);

        leaf1->add_parent(root_child1);
        leaf2->add_parent(root_child1);
        leaf3->add_parent(root_child1);

        leaf4->add_parent(root_child2);
        leaf5->add_parent(root_child2);

        REQUIRE(root->get_clade_length() == Approx(3.9));

        REQUIRE(root->degree() == 2);
        REQUIRE(root_child1->degree() == 4);
        REQUIRE(root_child1->get_number_of_children() == 3);
        REQUIRE(root_child2->degree() == 3);
        REQUIRE(leaf1->degree() == 1);
        REQUIRE(leaf2->degree() == 1);
        REQUIRE(leaf3->degree() == 1);
        REQUIRE(leaf4->degree() == 1);
        REQUIRE(leaf5->degree() == 1);

        REQUIRE(root->is_root());
        REQUIRE(! root->is_leaf());
        REQUIRE(root->has_children());
        REQUIRE(! root->has_parent());
        REQUIRE(root->is_dirty());
        REQUIRE(root->get_node_count() == 8);
        REQUIRE(root->get_leaf_node_count() == 5);
        REQUIRE(root->get_internal_node_count() == 3);
        REQUIRE(root->get_number_of_children() == 2);
        REQUIRE(root->get_number_of_parents() == 0);
        REQUIRE(root->get_height() == Approx(1.0));
        REQUIRE(root->get_length() == Approx(0.0));
        REQUIRE(root->is_child(root_child1));
        REQUIRE(root->is_child(root_child2));
        REQUIRE(! root->is_child(leaf1));

        REQUIRE(! root_child1->is_root());
        REQUIRE(! root_child1->is_leaf());
        REQUIRE(root_child1->has_parent());
        REQUIRE(root_child1->has_children());
        REQUIRE(root_child1->is_dirty());
        REQUIRE(root_child1->get_node_count() == 4);
        REQUIRE(root_child1->get_leaf_node_count() == 3);
        REQUIRE(root_child1->get_internal_node_count() == 1);
        REQUIRE(root_child1->get_number_of_children() == 3);
        REQUIRE(root_child1->get_number_of_parents() == 1);
        REQUIRE(root_child1->get_height() == Approx(0.8));
        REQUIRE(root_child1->get_length() == Approx(0.2));
        REQUIRE(root_child1->is_child(leaf1));
        REQUIRE(root_child1->is_child(leaf2));
        REQUIRE(root_child1->is_child(leaf3));
        REQUIRE(! root_child1->is_child(leaf4));
        REQUIRE(! root_child1->is_child(leaf5));
        REQUIRE(root_child1->is_parent(root));
        REQUIRE(! root_child1->is_parent(root_child2));

        REQUIRE(! root_child2->is_root());
        REQUIRE(! root_child2->is_leaf());
        REQUIRE(root_child2->has_parent());
        REQUIRE(root_child2->has_children());
        REQUIRE(root_child2->is_dirty());
        REQUIRE(root_child2->get_node_count() == 3);
        REQUIRE(root_child2->get_leaf_node_count() == 2);
        REQUIRE(root_child2->get_internal_node_count() == 1);
        REQUIRE(root_child2->get_number_of_children() == 2);
        REQUIRE(root_child2->get_number_of_parents() == 1);
        REQUIRE(root_child2->get_height() == Approx(0.3));
        REQUIRE(root_child2->get_length() == Approx(0.7));
        REQUIRE(! root_child2->is_child(leaf1));
        REQUIRE(! root_child2->is_child(leaf2));
        REQUIRE(! root_child2->is_child(leaf3));
        REQUIRE(root_child2->is_child(leaf4));
        REQUIRE(root_child2->is_child(leaf5));
        REQUIRE(root_child2->is_parent(root));
        REQUIRE(! root_child2->is_parent(root_child1));

        REQUIRE(! leaf1->is_root());
        REQUIRE(leaf1->is_leaf());
        REQUIRE(leaf1->has_parent());
        REQUIRE(! leaf1->has_children());
        REQUIRE(leaf1->is_dirty());
        REQUIRE(leaf1->get_node_count() == 1);
        REQUIRE(leaf1->get_leaf_node_count() == 1);
        REQUIRE(leaf1->get_internal_node_count() == 0);
        REQUIRE(leaf1->get_number_of_children() == 0);
        REQUIRE(leaf1->get_number_of_parents() == 1);
        REQUIRE(leaf1->get_height() == Approx(0.0));
        REQUIRE(leaf1->get_length() == Approx(0.8));
        REQUIRE(leaf1->is_parent(root_child1));

        REQUIRE(! leaf2->is_root());
        REQUIRE(leaf2->is_leaf());
        REQUIRE(leaf2->has_parent());
        REQUIRE(! leaf2->has_children());
        REQUIRE(leaf2->is_dirty());
        REQUIRE(leaf2->get_node_count() == 1);
        REQUIRE(leaf2->get_leaf_node_count() == 1);
        REQUIRE(leaf2->get_internal_node_count() == 0);
        REQUIRE(leaf2->get_number_of_children() == 0);
        REQUIRE(leaf2->get_number_of_parents() == 1);
        REQUIRE(leaf2->get_height() == Approx(0.0));
        REQUIRE(leaf2->get_length() == Approx(0.8));
        REQUIRE(leaf2->is_parent(root_child1));

        REQUIRE(! leaf3->is_root());
        REQUIRE(leaf3->is_leaf());
        REQUIRE(leaf3->has_parent());
        REQUIRE(! leaf3->has_children());
        REQUIRE(leaf3->is_dirty());
        REQUIRE(leaf3->get_node_count() == 1);
        REQUIRE(leaf3->get_leaf_node_count() == 1);
        REQUIRE(leaf3->get_internal_node_count() == 0);
        REQUIRE(leaf3->get_number_of_children() == 0);
        REQUIRE(leaf3->get_number_of_parents() == 1);
        REQUIRE(leaf3->get_height() == Approx(0.0));
        REQUIRE(leaf3->get_length() == Approx(0.8));
        REQUIRE(leaf3->is_parent(root_child1));

        REQUIRE(! leaf4->is_root());
        REQUIRE(leaf4->is_leaf());
        REQUIRE(leaf4->has_parent());
        REQUIRE(! leaf4->has_children());
        REQUIRE(leaf4->is_dirty());
        REQUIRE(leaf4->get_node_count() == 1);
        REQUIRE(leaf4->get_leaf_node_count() == 1);
        REQUIRE(leaf4->get_internal_node_count() == 0);
        REQUIRE(leaf4->get_number_of_children() == 0);
        REQUIRE(leaf4->get_number_of_parents() == 1);
        REQUIRE(leaf4->get_height() == Approx(0.0));
        REQUIRE(leaf4->get_length() == Approx(0.3));
        REQUIRE(leaf4->is_parent(root_child2));

        REQUIRE(! leaf5->is_root());
        REQUIRE(leaf5->is_leaf());
        REQUIRE(leaf5->has_parent());
        REQUIRE(! leaf5->has_children());
        REQUIRE(leaf5->is_dirty());
        REQUIRE(leaf5->get_node_count() == 1);
        REQUIRE(leaf5->get_leaf_node_count() == 1);
        REQUIRE(leaf5->get_internal_node_count() == 0);
        REQUIRE(leaf5->get_number_of_children() == 0);
        REQUIRE(leaf5->get_number_of_parents() == 1);
        REQUIRE(leaf5->get_height() == Approx(0.0));
        REQUIRE(leaf5->get_length() == Approx(0.3));
        REQUIRE(leaf5->is_parent(root_child2));
    }
}

TEST_CASE("Test simple tree cleaning Node", "[Node]") {

    SECTION("Testing tree cleaning") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> root_child1 = std::make_shared<Node>("root child 1", 0.8);
        std::shared_ptr<Node> root_child2 = std::make_shared<Node>("root child 2", 0.3);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf 1");
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf 2");
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf 3");
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf 4");
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf 5");

        root->add_child(root_child1);
        root->add_child(root_child2);

        leaf1->add_parent(root_child1);
        leaf2->add_parent(root_child1);
        leaf3->add_parent(root_child1);

        leaf4->add_parent(root_child2);
        leaf5->add_parent(root_child2);

        REQUIRE(root->is_dirty());
        REQUIRE(root_child1->is_dirty());
        REQUIRE(root_child2->is_dirty());
        REQUIRE(leaf1->is_dirty());
        REQUIRE(leaf2->is_dirty());
        REQUIRE(leaf3->is_dirty());
        REQUIRE(leaf4->is_dirty());
        REQUIRE(leaf5->is_dirty());
        REQUIRE(root->clade_has_dirt());
        REQUIRE(root_child1->clade_has_dirt());
        REQUIRE(root_child2->clade_has_dirt());
        REQUIRE(leaf1->clade_has_dirt());
        REQUIRE(leaf2->clade_has_dirt());
        REQUIRE(leaf3->clade_has_dirt());
        REQUIRE(leaf4->clade_has_dirt());
        REQUIRE(leaf5->clade_has_dirt());

        root->make_all_clean();

        REQUIRE(! root->is_dirty());
        REQUIRE(! root_child1->is_dirty());
        REQUIRE(! root_child2->is_dirty());
        REQUIRE(! leaf1->is_dirty());
        REQUIRE(! leaf2->is_dirty());
        REQUIRE(! leaf3->is_dirty());
        REQUIRE(! leaf4->is_dirty());
        REQUIRE(! leaf5->is_dirty());
        REQUIRE(! root->clade_has_dirt());
        REQUIRE(! root_child1->clade_has_dirt());
        REQUIRE(! root_child2->clade_has_dirt());
        REQUIRE(! leaf1->clade_has_dirt());
        REQUIRE(! leaf2->clade_has_dirt());
        REQUIRE(! leaf3->clade_has_dirt());
        REQUIRE(! leaf4->clade_has_dirt());
        REQUIRE(! leaf5->clade_has_dirt());

        root_child1->make_all_dirty();

        REQUIRE(! root->is_dirty());
        REQUIRE(root_child1->is_dirty());
        REQUIRE(! root_child2->is_dirty());
        REQUIRE(leaf1->is_dirty());
        REQUIRE(leaf2->is_dirty());
        REQUIRE(leaf3->is_dirty());
        REQUIRE(! leaf4->is_dirty());
        REQUIRE(! leaf5->is_dirty());
        REQUIRE(root->clade_has_dirt());
        REQUIRE(root_child1->clade_has_dirt());
        REQUIRE(! root_child2->clade_has_dirt());
        REQUIRE(leaf1->clade_has_dirt());
        REQUIRE(leaf2->clade_has_dirt());
        REQUIRE(leaf3->clade_has_dirt());
        REQUIRE(! leaf4->clade_has_dirt());
        REQUIRE(! leaf5->clade_has_dirt());

        root->make_all_clean();

        REQUIRE(! root->is_dirty());
        REQUIRE(! root_child1->is_dirty());
        REQUIRE(! root_child2->is_dirty());
        REQUIRE(! leaf1->is_dirty());
        REQUIRE(! leaf2->is_dirty());
        REQUIRE(! leaf3->is_dirty());
        REQUIRE(! leaf4->is_dirty());
        REQUIRE(! leaf5->is_dirty());
        REQUIRE(! root->clade_has_dirt());
        REQUIRE(! root_child1->clade_has_dirt());
        REQUIRE(! root_child2->clade_has_dirt());
        REQUIRE(! leaf1->clade_has_dirt());
        REQUIRE(! leaf2->clade_has_dirt());
        REQUIRE(! leaf3->clade_has_dirt());
        REQUIRE(! leaf4->clade_has_dirt());
        REQUIRE(! leaf5->clade_has_dirt());

        root_child2->make_all_dirty();

        REQUIRE(! root->is_dirty());
        REQUIRE(! root_child1->is_dirty());
        REQUIRE(root_child2->is_dirty());
        REQUIRE(! leaf1->is_dirty());
        REQUIRE(! leaf2->is_dirty());
        REQUIRE(! leaf3->is_dirty());
        REQUIRE(leaf4->is_dirty());
        REQUIRE(leaf5->is_dirty());
        REQUIRE(root->clade_has_dirt());
        REQUIRE(! root_child1->clade_has_dirt());
        REQUIRE(root_child2->clade_has_dirt());
        REQUIRE(! leaf1->clade_has_dirt());
        REQUIRE(! leaf2->clade_has_dirt());
        REQUIRE(! leaf3->clade_has_dirt());
        REQUIRE(leaf4->clade_has_dirt());
        REQUIRE(leaf5->clade_has_dirt());

        root_child1->make_all_clean();

        REQUIRE(! root->is_dirty());
        REQUIRE(! root_child1->is_dirty());
        REQUIRE(root_child2->is_dirty());
        REQUIRE(! leaf1->is_dirty());
        REQUIRE(! leaf2->is_dirty());
        REQUIRE(! leaf3->is_dirty());
        REQUIRE(leaf4->is_dirty());
        REQUIRE(leaf5->is_dirty());
        REQUIRE(root->clade_has_dirt());
        REQUIRE(! root_child1->clade_has_dirt());
        REQUIRE(root_child2->clade_has_dirt());
        REQUIRE(! leaf1->clade_has_dirt());
        REQUIRE(! leaf2->clade_has_dirt());
        REQUIRE(! leaf3->clade_has_dirt());
        REQUIRE(leaf4->clade_has_dirt());
        REQUIRE(leaf5->clade_has_dirt());

        root_child2->make_all_clean();

        REQUIRE(! root->is_dirty());
        REQUIRE(! root_child1->is_dirty());
        REQUIRE(! root_child2->is_dirty());
        REQUIRE(! leaf1->is_dirty());
        REQUIRE(! leaf2->is_dirty());
        REQUIRE(! leaf3->is_dirty());
        REQUIRE(! leaf4->is_dirty());
        REQUIRE(! leaf5->is_dirty());
        REQUIRE(! root->clade_has_dirt());
        REQUIRE(! root_child1->clade_has_dirt());
        REQUIRE(! root_child2->clade_has_dirt());
        REQUIRE(! leaf1->clade_has_dirt());
        REQUIRE(! leaf2->clade_has_dirt());
        REQUIRE(! leaf3->clade_has_dirt());
        REQUIRE(! leaf4->clade_has_dirt());
        REQUIRE(! leaf5->clade_has_dirt());

        leaf1->make_all_dirty();

        REQUIRE(! root->is_dirty());
        REQUIRE(! root_child1->is_dirty());
        REQUIRE(! root_child2->is_dirty());
        REQUIRE(leaf1->is_dirty());
        REQUIRE(! leaf2->is_dirty());
        REQUIRE(! leaf3->is_dirty());
        REQUIRE(! leaf4->is_dirty());
        REQUIRE(! leaf5->is_dirty());
        REQUIRE(root->clade_has_dirt());
        REQUIRE(root_child1->clade_has_dirt());
        REQUIRE(! root_child2->clade_has_dirt());
        REQUIRE(leaf1->clade_has_dirt());
        REQUIRE(! leaf2->clade_has_dirt());
        REQUIRE(! leaf3->clade_has_dirt());
        REQUIRE(! leaf4->clade_has_dirt());
        REQUIRE(! leaf5->clade_has_dirt());

        root->make_all_dirty();

        REQUIRE(root->is_dirty());
        REQUIRE(root_child1->is_dirty());
        REQUIRE(root_child2->is_dirty());
        REQUIRE(leaf1->is_dirty());
        REQUIRE(leaf2->is_dirty());
        REQUIRE(leaf3->is_dirty());
        REQUIRE(leaf4->is_dirty());
        REQUIRE(leaf5->is_dirty());
        REQUIRE(root->clade_has_dirt());
        REQUIRE(root_child1->clade_has_dirt());
        REQUIRE(root_child2->clade_has_dirt());
        REQUIRE(leaf1->clade_has_dirt());
        REQUIRE(leaf2->clade_has_dirt());
        REQUIRE(leaf3->clade_has_dirt());
        REQUIRE(leaf4->clade_has_dirt());
        REQUIRE(leaf5->clade_has_dirt());

        root->make_all_clean();

        REQUIRE(! root->is_dirty());
        REQUIRE(! root_child1->is_dirty());
        REQUIRE(! root_child2->is_dirty());
        REQUIRE(! leaf1->is_dirty());
        REQUIRE(! leaf2->is_dirty());
        REQUIRE(! leaf3->is_dirty());
        REQUIRE(! leaf4->is_dirty());
        REQUIRE(! leaf5->is_dirty());
        REQUIRE(! root->clade_has_dirt());
        REQUIRE(! root_child1->clade_has_dirt());
        REQUIRE(! root_child2->clade_has_dirt());
        REQUIRE(! leaf1->clade_has_dirt());
        REQUIRE(! leaf2->clade_has_dirt());
        REQUIRE(! leaf3->clade_has_dirt());
        REQUIRE(! leaf4->clade_has_dirt());
        REQUIRE(! leaf5->clade_has_dirt());

        root_child2->make_dirty();

        REQUIRE(! root->is_dirty());
        REQUIRE(! root_child1->is_dirty());
        REQUIRE(root_child2->is_dirty());
        REQUIRE(! leaf1->is_dirty());
        REQUIRE(! leaf2->is_dirty());
        REQUIRE(! leaf3->is_dirty());
        REQUIRE(! leaf4->is_dirty());
        REQUIRE(! leaf5->is_dirty());
        REQUIRE(root->clade_has_dirt());
        REQUIRE(! root_child1->clade_has_dirt());
        REQUIRE(root_child2->clade_has_dirt());
        REQUIRE(! leaf1->clade_has_dirt());
        REQUIRE(! leaf2->clade_has_dirt());
        REQUIRE(! leaf3->clade_has_dirt());
        REQUIRE(! leaf4->clade_has_dirt());
        REQUIRE(! leaf5->clade_has_dirt());
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

        std::shared_ptr<PopulationNode> p = n.get_parent();
        REQUIRE(p == nullptr);
        REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
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

        std::shared_ptr<PopulationNode> p = n.get_parent();
        REQUIRE(p == nullptr);
        REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
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

        std::shared_ptr<PopulationNode> p = n.get_parent();
        REQUIRE(p == nullptr);
        REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
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

        std::shared_ptr<PopulationNode> p = n.get_parent();
        REQUIRE(p == nullptr);
        REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
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

        std::shared_ptr<PopulationNode> p = n.get_parent();
        REQUIRE(p == nullptr);
        REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
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

        std::shared_ptr<PopulationNode> p = n.get_parent();
        REQUIRE(p == nullptr);
        REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
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

        std::shared_ptr<PopulationNode> p = n.get_parent();
        REQUIRE(p == nullptr);
        REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
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

        std::shared_ptr<PopulationNode> p = n.get_parent();
        REQUIRE(p == nullptr);
        REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
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

        std::shared_ptr<PopulationNode> p = n.get_parent();
        REQUIRE(p == nullptr);
        REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
        
        n.set_label("leaf2");
        REQUIRE(n.get_label() == "leaf2");
        REQUIRE(n2.get_label() == "leaf1");

        n2.set_height(0.06);
        n2.reset(2);
        n2.set_bottom_pattern_probability(1, 1, 1.0);
        REQUIRE(n.get_height() == Approx(0.06));
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

        std::shared_ptr<PopulationNode> p = n.get_parent();
        REQUIRE(p == nullptr);
        REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code());

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range);
        
        n.set_label("leaf2");
        REQUIRE(n.get_label() == "leaf2");
        REQUIRE(n2.get_label() == "leaf1");

        n2.set_height(0.06);
        n2.reset(2);
        n2.set_bottom_pattern_probability(1, 1, 1.0);
        REQUIRE(n.get_height() == Approx(0.06));
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
    }
}

/* TEST_CASE("Testing clone method of PopulationNode", "[PopulationNode]") { */

/*     SECTION("Testing clone") { */
/*         PopulationNode * n2 = new PopulationNode("leaf1", 0.03, (unsigned int)3); */
/*         n2->set_bottom_pattern_probability(2, 1, 1.0); */
/*         PopulationNode * n = n2->clone(); */
/*         REQUIRE(typeid(n).hash_code() == typeid(n2).hash_code()); */
/*         REQUIRE(n->get_height() == Approx(0.03)); */
/*         REQUIRE(n->get_length() == 0.0); */
/*         REQUIRE(n->get_label() == "leaf1"); */
/*         REQUIRE(n->get_allele_count() == 3); */
/*         REQUIRE(n->get_bottom_pattern_probability(2, 1) == Approx(1.0)); */
/*         REQUIRE(n->is_dirty()); */

/*         REQUIRE(n->degree() == 0); */
/*         REQUIRE(n->has_parent() == false); */
/*         REQUIRE(n->get_number_of_parents() == 0); */
/*         REQUIRE(n->has_children() == false); */
/*         REQUIRE(n->get_number_of_children() == 0); */

/*         REQUIRE(n->is_leaf() == true); */
/*         REQUIRE(n->is_root() == true); */
/*         REQUIRE(n->get_node_count() == 1); */
/*         REQUIRE(n->get_leaf_node_count() == 1); */
/*         REQUIRE(n->get_internal_node_count() == 0); */

/*         std::shared_ptr<PopulationNode> p = n.get_parent(); */
/*         REQUIRE(p == nullptr); */
/*         REQUIRE(typeid(n).hash_code() == typeid(p).hash_code()); */

/*         REQUIRE_THROWS_AS(n->get_child(0), std::out_of_range); */
        
/*         n->set_label("leaf2"); */
/*         REQUIRE(n->get_label() == "leaf2"); */
/*         REQUIRE(n2->get_label() == "leaf1"); */

/*         n2->set_height(0.06); */
/*         n2->reset(2); */
/*         n2->set_bottom_pattern_probability(1, 1, 1.0); */
/*         REQUIRE(n->get_height() == Approx(0.03)); */
/*         REQUIRE(n2->get_height() == Approx(0.06)); */
/*         REQUIRE(n->get_allele_count() == 3); */
/*         REQUIRE(n->get_bottom_pattern_probability(2, 1) == Approx(1.0)); */
/*         REQUIRE(n->get_bottom_pattern_probability(1, 1) == Approx(0.0)); */
/*         REQUIRE(n2->get_allele_count() == 2); */
/*         REQUIRE(n2->get_bottom_pattern_probability(2, 1) == Approx(0.0)); */
/*         REQUIRE(n2->get_bottom_pattern_probability(1, 1) == Approx(1.0)); */

/*         std::vector<double> e_bottom_n = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0}; */
/*         std::vector<double> e_top_n = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; */
/*         std::vector<double> e_bottom_n2 = {0.0, 1.0, 0.0, 0.0, 0.0}; */
/*         std::vector<double> e_top_n2 = {0.0, 0.0, 0.0, 0.0, 0.0}; */
/*         REQUIRE(n->get_bottom_pattern_probs().get_pattern_prob_matrix() == e_bottom_n); */
/*         REQUIRE(n->get_top_pattern_probs().get_pattern_prob_matrix() == e_top_n); */
/*         REQUIRE(n2->get_bottom_pattern_probs().get_pattern_prob_matrix() == e_bottom_n2); */
/*         REQUIRE(n2->get_top_pattern_probs().get_pattern_prob_matrix() == e_top_n2); */

/*         delete n; */
/*         delete n2; */
/*     } */
/* } */

TEST_CASE("Testing parent methods of PopulationNode", "[PopulationNode]") {

    SECTION("Testing parent methods") {
        std::shared_ptr<PopulationNode> c = std::make_shared<PopulationNode>();
        std::shared_ptr<PopulationNode> p = std::make_shared<PopulationNode>();

        REQUIRE(c->degree() == 0);
        REQUIRE(p->degree() == 0);
        REQUIRE(c->has_parent() == false);
        REQUIRE(c->get_number_of_parents() == 0);

        REQUIRE(c->get_parent() == nullptr);
        REQUIRE_THROWS_AS(p->get_child(0), std::out_of_range);
        REQUIRE(c->is_parent(p) == false);
        REQUIRE(p->is_child(c) == false);

        c->add_parent(p);
        REQUIRE(c->has_parent() == true);
        REQUIRE(c->is_parent(p) == true);
        REQUIRE(p->is_child(c) == true);
        REQUIRE(c->get_number_of_parents() == 1);
        REQUIRE(c->degree() == 1);
        REQUIRE(p->degree() == 1);

        std::shared_ptr<PopulationNode> r = c->remove_parent();
        REQUIRE(r != nullptr);
        REQUIRE(r == p);
        REQUIRE(typeid(r).hash_code() == typeid(c).hash_code());
        REQUIRE(c->has_parent() == false);
        REQUIRE(c->is_parent(p) == false);
        REQUIRE(p->is_child(c) == false);
        REQUIRE(c->get_number_of_parents() == 0);
        REQUIRE(c->degree() == 0);
        REQUIRE(p->degree() == 0);
    }
}

TEST_CASE("Testing child methods of PopulationNode", "[PopulationNode]") {

    SECTION("Testing child methods") {
        std::shared_ptr<PopulationNode> c = std::make_shared<PopulationNode>();
        std::shared_ptr<PopulationNode> p = std::make_shared<PopulationNode>();
        REQUIRE(typeid(c).hash_code() == typeid(p).hash_code());

        REQUIRE(c->degree() == 0);
        REQUIRE(p->degree() == 0);
        REQUIRE(p->has_children() == false);
        REQUIRE(p->get_number_of_children() == 0);

        REQUIRE(c->get_parent() == nullptr);
        REQUIRE_THROWS_AS(p->get_child(0), std::out_of_range);
        REQUIRE(c->is_parent(p) == false);
        REQUIRE(p->is_child(c) == false);

        p->add_child(c);
        REQUIRE(p->has_children() == true);
        REQUIRE(c->is_parent(p) == true);
        REQUIRE(p->is_child(c) == true);
        REQUIRE(p->get_number_of_children() == 1);
        REQUIRE(c->degree() == 1);
        REQUIRE(p->degree() == 1);

        unsigned int i = 1;
        REQUIRE_THROWS_AS(p->remove_child(i), std::out_of_range);
        i = 0;
        std::shared_ptr<PopulationNode> r = p->remove_child(i);
        REQUIRE(r == c);
        REQUIRE(typeid(c).hash_code() == typeid(r).hash_code());
        REQUIRE(p->has_children() == false);
        REQUIRE(c->is_parent(p) == false);
        REQUIRE(p->is_child(c) == false);
        REQUIRE(p->get_number_of_children() == 0);
        REQUIRE(c->degree() == 0);
        REQUIRE(p->degree() == 0);

        p->add_child(c);
        REQUIRE(p->has_children() == true);
        REQUIRE(c->is_parent(p) == true);
        REQUIRE(p->is_child(c) == true);
        REQUIRE(p->get_number_of_children() == 1);
        REQUIRE(c->degree() == 1);
        REQUIRE(p->degree() == 1);

        p->remove_child(c);
        REQUIRE(p->has_children() == false);
        REQUIRE(c->is_parent(p) == false);
        REQUIRE(p->is_child(c) == false);
        REQUIRE(p->get_number_of_children() == 0);
        REQUIRE(c->degree() == 0);
        REQUIRE(p->degree() == 0);
        
    }
}

TEST_CASE("Test simple tree building with PopulationNode", "[PopulationNode]") {

    SECTION("Testing tree building") {
        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>("root", 1.0);
        std::shared_ptr<PopulationNode> root_child1 = std::make_shared<PopulationNode>("root child 1", 0.8);
        std::shared_ptr<PopulationNode> root_child2 = std::make_shared<PopulationNode>("root child 2", 0.3);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>("leaf 1", (unsigned int)3);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>("leaf 2", (unsigned int)3);
        std::shared_ptr<PopulationNode> leaf3 = std::make_shared<PopulationNode>("leaf 3", (unsigned int)3);
        std::shared_ptr<PopulationNode> leaf4 = std::make_shared<PopulationNode>("leaf 4", (unsigned int)3);
        std::shared_ptr<PopulationNode> leaf5 = std::make_shared<PopulationNode>("leaf 5", (unsigned int)3);

        root->add_child(root_child1);
        root->add_child(root_child2);

        leaf1->add_parent(root_child1);
        leaf2->add_parent(root_child1);
        leaf3->add_parent(root_child1);

        leaf4->add_parent(root_child2);
        leaf5->add_parent(root_child2);

        REQUIRE(root->degree() == 2);
        REQUIRE(root_child1->degree() == 4);
        REQUIRE(root_child1->get_number_of_children() == 3);
        REQUIRE(root_child2->degree() == 3);
        REQUIRE(leaf1->degree() == 1);
        REQUIRE(leaf2->degree() == 1);
        REQUIRE(leaf3->degree() == 1);
        REQUIRE(leaf4->degree() == 1);
        REQUIRE(leaf5->degree() == 1);

        REQUIRE(root->is_root());
        REQUIRE(! root->is_leaf());
        REQUIRE(root->has_children());
        REQUIRE(! root->has_parent());
        REQUIRE(root->is_dirty());
        REQUIRE(root->get_node_count() == 8);
        REQUIRE(root->get_leaf_node_count() == 5);
        REQUIRE(root->get_internal_node_count() == 3);
        REQUIRE(root->get_number_of_children() == 2);
        REQUIRE(root->get_number_of_parents() == 0);
        REQUIRE(root->get_height() == Approx(1.0));
        REQUIRE(root->get_length() == Approx(0.0));
        REQUIRE(root->is_child(root_child1));
        REQUIRE(root->is_child(root_child2));
        REQUIRE(! root->is_child(leaf1));

        REQUIRE(! root_child1->is_root());
        REQUIRE(! root_child1->is_leaf());
        REQUIRE(root_child1->has_parent());
        REQUIRE(root_child1->has_children());
        REQUIRE(root_child1->is_dirty());
        REQUIRE(root_child1->get_node_count() == 4);
        REQUIRE(root_child1->get_leaf_node_count() == 3);
        REQUIRE(root_child1->get_internal_node_count() == 1);
        REQUIRE(root_child1->get_number_of_children() == 3);
        REQUIRE(root_child1->get_number_of_parents() == 1);
        REQUIRE(root_child1->get_height() == Approx(0.8));
        REQUIRE(root_child1->get_length() == Approx(0.2));
        REQUIRE(root_child1->is_child(leaf1));
        REQUIRE(root_child1->is_child(leaf2));
        REQUIRE(root_child1->is_child(leaf3));
        REQUIRE(! root_child1->is_child(leaf4));
        REQUIRE(! root_child1->is_child(leaf5));
        REQUIRE(root_child1->is_parent(root));
        REQUIRE(! root_child1->is_parent(root_child2));

        REQUIRE(! root_child2->is_root());
        REQUIRE(! root_child2->is_leaf());
        REQUIRE(root_child2->has_parent());
        REQUIRE(root_child2->has_children());
        REQUIRE(root_child2->is_dirty());
        REQUIRE(root_child2->get_node_count() == 3);
        REQUIRE(root_child2->get_leaf_node_count() == 2);
        REQUIRE(root_child2->get_internal_node_count() == 1);
        REQUIRE(root_child2->get_number_of_children() == 2);
        REQUIRE(root_child2->get_number_of_parents() == 1);
        REQUIRE(root_child2->get_height() == Approx(0.3));
        REQUIRE(root_child2->get_length() == Approx(0.7));
        REQUIRE(! root_child2->is_child(leaf1));
        REQUIRE(! root_child2->is_child(leaf2));
        REQUIRE(! root_child2->is_child(leaf3));
        REQUIRE(root_child2->is_child(leaf4));
        REQUIRE(root_child2->is_child(leaf5));
        REQUIRE(root_child2->is_parent(root));
        REQUIRE(! root_child2->is_parent(root_child1));

        REQUIRE(! leaf1->is_root());
        REQUIRE(leaf1->is_leaf());
        REQUIRE(leaf1->has_parent());
        REQUIRE(! leaf1->has_children());
        REQUIRE(leaf1->is_dirty());
        REQUIRE(leaf1->get_node_count() == 1);
        REQUIRE(leaf1->get_leaf_node_count() == 1);
        REQUIRE(leaf1->get_internal_node_count() == 0);
        REQUIRE(leaf1->get_number_of_children() == 0);
        REQUIRE(leaf1->get_number_of_parents() == 1);
        REQUIRE(leaf1->get_height() == Approx(0.0));
        REQUIRE(leaf1->get_length() == Approx(0.8));
        REQUIRE(leaf1->is_parent(root_child1));

        REQUIRE(! leaf2->is_root());
        REQUIRE(leaf2->is_leaf());
        REQUIRE(leaf2->has_parent());
        REQUIRE(! leaf2->has_children());
        REQUIRE(leaf2->is_dirty());
        REQUIRE(leaf2->get_node_count() == 1);
        REQUIRE(leaf2->get_leaf_node_count() == 1);
        REQUIRE(leaf2->get_internal_node_count() == 0);
        REQUIRE(leaf2->get_number_of_children() == 0);
        REQUIRE(leaf2->get_number_of_parents() == 1);
        REQUIRE(leaf2->get_height() == Approx(0.0));
        REQUIRE(leaf2->get_length() == Approx(0.8));
        REQUIRE(leaf2->is_parent(root_child1));

        REQUIRE(! leaf3->is_root());
        REQUIRE(leaf3->is_leaf());
        REQUIRE(leaf3->has_parent());
        REQUIRE(! leaf3->has_children());
        REQUIRE(leaf3->is_dirty());
        REQUIRE(leaf3->get_node_count() == 1);
        REQUIRE(leaf3->get_leaf_node_count() == 1);
        REQUIRE(leaf3->get_internal_node_count() == 0);
        REQUIRE(leaf3->get_number_of_children() == 0);
        REQUIRE(leaf3->get_number_of_parents() == 1);
        REQUIRE(leaf3->get_height() == Approx(0.0));
        REQUIRE(leaf3->get_length() == Approx(0.8));
        REQUIRE(leaf3->is_parent(root_child1));

        REQUIRE(! leaf4->is_root());
        REQUIRE(leaf4->is_leaf());
        REQUIRE(leaf4->has_parent());
        REQUIRE(! leaf4->has_children());
        REQUIRE(leaf4->is_dirty());
        REQUIRE(leaf4->get_node_count() == 1);
        REQUIRE(leaf4->get_leaf_node_count() == 1);
        REQUIRE(leaf4->get_internal_node_count() == 0);
        REQUIRE(leaf4->get_number_of_children() == 0);
        REQUIRE(leaf4->get_number_of_parents() == 1);
        REQUIRE(leaf4->get_height() == Approx(0.0));
        REQUIRE(leaf4->get_length() == Approx(0.3));
        REQUIRE(leaf4->is_parent(root_child2));

        REQUIRE(! leaf5->is_root());
        REQUIRE(leaf5->is_leaf());
        REQUIRE(leaf5->has_parent());
        REQUIRE(! leaf5->has_children());
        REQUIRE(leaf5->is_dirty());
        REQUIRE(leaf5->get_node_count() == 1);
        REQUIRE(leaf5->get_leaf_node_count() == 1);
        REQUIRE(leaf5->get_internal_node_count() == 0);
        REQUIRE(leaf5->get_number_of_children() == 0);
        REQUIRE(leaf5->get_number_of_parents() == 1);
        REQUIRE(leaf5->get_height() == Approx(0.0));
        REQUIRE(leaf5->get_length() == Approx(0.3));
        REQUIRE(leaf5->is_parent(root_child2));

        REQUIRE(leaf1->get_allele_count() == 3);
        REQUIRE(leaf2->get_allele_count() == 3);
        REQUIRE(leaf3->get_allele_count() == 3);
        REQUIRE(leaf4->get_allele_count() == 3);
        REQUIRE(leaf5->get_allele_count() == 3);
        REQUIRE(root_child1->get_allele_count() == 0);
        REQUIRE(root_child2->get_allele_count() == 0);
        REQUIRE(root->get_allele_count() == 0);
        
        root->resize_all();
        REQUIRE(leaf1->get_allele_count() == 3);
        REQUIRE(leaf2->get_allele_count() == 3);
        REQUIRE(leaf3->get_allele_count() == 3);
        REQUIRE(leaf4->get_allele_count() == 3);
        REQUIRE(leaf5->get_allele_count() == 3);
        REQUIRE(root_child1->get_allele_count() == 9);
        REQUIRE(root_child2->get_allele_count() == 6);
        REQUIRE(root->get_allele_count() == 15);
    }
}

TEST_CASE("Test node height prior", "[Node]") {

    SECTION("Testing prior") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> root_child1 = std::make_shared<Node>("root child 1", 0.8);
        std::shared_ptr<Node> root_child2 = std::make_shared<Node>("root child 2", 0.3);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf 1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf 2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf 3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf 4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf 5", 0.0);

        root->add_child(root_child1);
        root->add_child(root_child2);

        leaf1->add_parent(root_child1);
        leaf2->add_parent(root_child1);
        leaf3->add_parent(root_child1);

        leaf4->add_parent(root_child2);
        leaf5->add_parent(root_child2);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(1.0, 1.0);
        root->set_all_node_height_priors(prior);
        REQUIRE(root->calculate_ln_relative_node_height_prior_density() == Approx(-2.1));

        leaf1->fix_node_height();
        leaf2->fix_node_height();
        leaf3->fix_node_height();
        leaf4->fix_node_height();
        leaf5->fix_node_height();

        REQUIRE(root->calculate_ln_relative_node_height_prior_density() == Approx(-2.1));


        std::shared_ptr<GammaDistribution> prior2 = std::make_shared<GammaDistribution>(1.0, 0.01);
        root->set_all_node_height_priors(prior2);
        REQUIRE(root->calculate_ln_relative_node_height_prior_density() == Approx(-196.18448944203573));

        leaf1->estimate_node_height();
        leaf2->estimate_node_height();
        leaf3->estimate_node_height();
        leaf4->estimate_node_height();
        leaf5->estimate_node_height();

        REQUIRE(root->calculate_ln_relative_node_height_prior_density() == Approx(-173.15863851209528));
    }
}

TEST_CASE("Test node height and population size priors", "[PopulationNode]") {

    SECTION("Testing prior") {
        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>("root", 1.0);
        std::shared_ptr<PopulationNode> root_child1 = std::make_shared<PopulationNode>("root child 1", 0.8);
        std::shared_ptr<PopulationNode> root_child2 = std::make_shared<PopulationNode>("root child 2", 0.3);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>("leaf 1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>("leaf 2", 0.0);
        std::shared_ptr<PopulationNode> leaf3 = std::make_shared<PopulationNode>("leaf 3", 0.0);
        std::shared_ptr<PopulationNode> leaf4 = std::make_shared<PopulationNode>("leaf 4", 0.0);
        std::shared_ptr<PopulationNode> leaf5 = std::make_shared<PopulationNode>("leaf 5", 0.0);

        root->add_child(root_child1);
        root->add_child(root_child2);

        leaf1->add_parent(root_child1);
        leaf2->add_parent(root_child1);
        leaf3->add_parent(root_child1);

        leaf4->add_parent(root_child2);
        leaf5->add_parent(root_child2);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(1.0, 1.0);
        root->set_all_node_height_priors(prior);
        root->set_all_population_size_priors(prior);
        REQUIRE(root->calculate_ln_relative_node_height_prior_density() == Approx(-2.1));
        REQUIRE(root->calculate_ln_relative_population_size_prior_density() == Approx(-0.001*8));

        leaf1->fix_node_height();
        leaf2->fix_node_height();
        leaf3->fix_node_height();
        leaf4->fix_node_height();
        leaf5->fix_node_height();

        REQUIRE(root->calculate_ln_relative_node_height_prior_density() == Approx(-2.1));
        REQUIRE(root->calculate_ln_relative_population_size_prior_density() == Approx(-0.001*8));


        std::shared_ptr<GammaDistribution> prior2 = std::make_shared<GammaDistribution>(1.0, 0.01);
        root->set_all_node_height_priors(prior2);
        root->set_all_population_size_priors(prior2);
        REQUIRE(root->calculate_ln_relative_node_height_prior_density() == Approx(-196.18448944203573));
        REQUIRE(root->calculate_ln_relative_population_size_prior_density() == Approx(4.5051701859880913*8));

        leaf1->estimate_node_height();
        leaf2->estimate_node_height();
        leaf3->estimate_node_height();
        leaf4->estimate_node_height();
        leaf5->estimate_node_height();

        REQUIRE(root->calculate_ln_relative_node_height_prior_density() == Approx(-173.15863851209528));
        REQUIRE(root->calculate_ln_relative_population_size_prior_density() == Approx(4.5051701859880913*8));

        root->set_all_population_sizes(0.02);

        REQUIRE(root->calculate_ln_relative_node_height_prior_density() == Approx(-173.15863851209528));
        REQUIRE(root->calculate_ln_relative_population_size_prior_density() == Approx(2.6051701859880909*8));

        root->set_all_population_size_parameters();
        root_child1->set_height_parameter(root_child2->get_height_parameter());

        REQUIRE(root->calculate_ln_relative_node_height_prior_density() == Approx(-97.76380869808338));
        REQUIRE(root->calculate_ln_relative_population_size_prior_density() == Approx(2.6051701859880909));

        leaf1->fix_node_height();
        leaf2->fix_node_height();
        leaf3->fix_node_height();
        leaf4->fix_node_height();
        leaf5->fix_node_height();

        REQUIRE(root->calculate_ln_relative_node_height_prior_density() == Approx(-120.78965962802383));
        REQUIRE(root->calculate_ln_relative_population_size_prior_density() == Approx(2.6051701859880909));
    }
}

TEST_CASE("Test node height fixing", "[Node]") {

    SECTION("Testing fixing of node height parameters") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> root_child1 = std::make_shared<Node>("root child 1", 0.8);
        std::shared_ptr<Node> root_child2 = std::make_shared<Node>("root child 2", 0.3);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf 1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf 2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf 3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf 4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf 5", 0.0);

        root->add_child(root_child1);
        root->add_child(root_child2);

        leaf1->add_parent(root_child1);
        leaf2->add_parent(root_child1);
        leaf3->add_parent(root_child1);

        leaf4->add_parent(root_child2);
        leaf5->add_parent(root_child2);

        REQUIRE(root->node_height_is_fixed() == false);
        REQUIRE(root->all_node_heights_are_fixed() == false);

        root->fix_all_node_heights();

        REQUIRE(root->node_height_is_fixed() == true);
        REQUIRE(root->all_node_heights_are_fixed() == true);

        root->estimate_all_node_heights();

        REQUIRE(root->node_height_is_fixed() == false);
        REQUIRE(root->all_node_heights_are_fixed() == false);

        root_child1->fix_node_height();

        REQUIRE(root_child1->node_height_is_fixed() == true);
        REQUIRE(root_child1->all_node_heights_are_fixed() == false);
        REQUIRE(leaf1->node_height_is_fixed() == false);
        REQUIRE(leaf2->node_height_is_fixed() == false);
        REQUIRE(leaf3->node_height_is_fixed() == false);

        root_child1->fix_all_node_heights();

        REQUIRE(root_child1->node_height_is_fixed() == true);
        REQUIRE(root_child1->all_node_heights_are_fixed() == true);
        REQUIRE(leaf1->node_height_is_fixed() == true);
        REQUIRE(leaf2->node_height_is_fixed() == true);
        REQUIRE(leaf3->node_height_is_fixed() == true);

        leaf2->estimate_node_height();
        REQUIRE(root_child1->node_height_is_fixed() == true);
        REQUIRE(root_child1->all_node_heights_are_fixed() == false);
        REQUIRE(leaf1->node_height_is_fixed() == true);
        REQUIRE(leaf2->node_height_is_fixed() == false);
        REQUIRE(leaf3->node_height_is_fixed() == true);
    }
}

TEST_CASE("Test population size fixing", "[PopulationNode]") {

    SECTION("Testing fixing of pop size parameters") {
        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>("root", 1.0);
        std::shared_ptr<PopulationNode> root_child1 = std::make_shared<PopulationNode>("root child 1", 0.8);
        std::shared_ptr<PopulationNode> root_child2 = std::make_shared<PopulationNode>("root child 2", 0.3);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>("leaf 1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>("leaf 2", 0.0);
        std::shared_ptr<PopulationNode> leaf3 = std::make_shared<PopulationNode>("leaf 3", 0.0);
        std::shared_ptr<PopulationNode> leaf4 = std::make_shared<PopulationNode>("leaf 4", 0.0);
        std::shared_ptr<PopulationNode> leaf5 = std::make_shared<PopulationNode>("leaf 5", 0.0);

        root->add_child(root_child1);
        root->add_child(root_child2);

        leaf1->add_parent(root_child1);
        leaf2->add_parent(root_child1);
        leaf3->add_parent(root_child1);

        leaf4->add_parent(root_child2);
        leaf5->add_parent(root_child2);

        REQUIRE(root->population_size_is_fixed() == false);
        REQUIRE(root->all_population_sizes_are_fixed() == false);

        root->fix_all_population_sizes();

        REQUIRE(root->population_size_is_fixed() == true);
        REQUIRE(root->all_population_sizes_are_fixed() == true);

        root->estimate_all_population_sizes();

        REQUIRE(root->population_size_is_fixed() == false);
        REQUIRE(root->all_population_sizes_are_fixed() == false);

        root_child1->fix_population_size();

        REQUIRE(root_child1->population_size_is_fixed() == true);
        REQUIRE(root_child1->all_population_sizes_are_fixed() == false);
        REQUIRE(leaf1->population_size_is_fixed() == false);
        REQUIRE(leaf2->population_size_is_fixed() == false);
        REQUIRE(leaf3->population_size_is_fixed() == false);

        root_child1->fix_all_population_sizes();

        REQUIRE(root_child1->population_size_is_fixed() == true);
        REQUIRE(root_child1->all_population_sizes_are_fixed() == true);
        REQUIRE(leaf1->population_size_is_fixed() == true);
        REQUIRE(leaf2->population_size_is_fixed() == true);
        REQUIRE(leaf3->population_size_is_fixed() == true);

        leaf2->estimate_population_size();
        REQUIRE(root_child1->population_size_is_fixed() == true);
        REQUIRE(root_child1->all_population_sizes_are_fixed() == false);
        REQUIRE(leaf1->population_size_is_fixed() == true);
        REQUIRE(leaf2->population_size_is_fixed() == false);
        REQUIRE(leaf3->population_size_is_fixed() == true);
    }
}

TEST_CASE("Test to_parentheses with Node", "[Node]") {

    SECTION("Testing newick string") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> root_child1 = std::make_shared<Node>("root child 1", 0.8);
        std::shared_ptr<Node> a = std::make_shared<Node>("A");
        std::shared_ptr<Node> b = std::make_shared<Node>("B");
        std::shared_ptr<Node> c = std::make_shared<Node>("C");

        root->add_child(root_child1);
        root->add_child(c);

        a->add_parent(root_child1);
        b->add_parent(root_child1);
        std::string p = root->to_parentheses();
        REQUIRE(p == "((A:0.8,B:0.8):0.2,C:1):0");

        std::shared_ptr<Node> d = std::make_shared<Node>("D");
        d->add_parent(root_child1);

        p = root->to_parentheses();
        REQUIRE(p == "((A:0.8,B:0.8,D:0.8):0.2,C:1):0");

        std::shared_ptr<Node> e = std::make_shared<Node>("E");
        std::shared_ptr<Node> f = std::make_shared<Node>("F");
        std::shared_ptr<Node> root_child2 = std::make_shared<Node>("root child 1", 0.3);
        root_child2->add_child(e);
        root_child2->add_child(f);
        root->add_child(root_child2);

        p = root->to_parentheses();
        REQUIRE(p == "((A:0.8,B:0.8,D:0.8):0.2,C:1,(E:0.3,F:0.3):0.7):0");
    }

    SECTION("Testing singleton") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);

        std::string p = root->to_parentheses();
        REQUIRE(p == "root:0");

        std::shared_ptr<Node> a = std::make_shared<Node>("A");

        root->add_child(a);

        p = root->to_parentheses();
        REQUIRE(p == "(A:1):0");
    }
}

TEST_CASE("Test clade cloning with PopulationNode", "[PopulationNode]") {

    SECTION("Testing clade cloning") {
        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>("root", 1.0);
        std::shared_ptr<PopulationNode> root_child1 = std::make_shared<PopulationNode>("root child 1", 0.8);
        std::shared_ptr<PopulationNode> root_child2 = std::make_shared<PopulationNode>("root child 2", 0.3);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>("leaf 1", (unsigned int)3);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>("leaf 2", (unsigned int)3);
        std::shared_ptr<PopulationNode> leaf3 = std::make_shared<PopulationNode>("leaf 3", (unsigned int)3);
        std::shared_ptr<PopulationNode> leaf4 = std::make_shared<PopulationNode>("leaf 4", (unsigned int)3);
        std::shared_ptr<PopulationNode> leaf5 = std::make_shared<PopulationNode>("leaf 5", (unsigned int)3);

        root->add_child(root_child1);
        root->add_child(root_child2);

        leaf1->add_parent(root_child1);
        leaf2->add_parent(root_child1);
        leaf3->add_parent(root_child1);

        leaf4->add_parent(root_child2);
        leaf5->add_parent(root_child2);

        std::shared_ptr<PopulationNode> clone = root->get_clade_clone();

        REQUIRE(root->degree() == 2);
        REQUIRE(root->is_root());
        REQUIRE(! root->is_leaf());
        REQUIRE(root->has_children());
        REQUIRE(! root->has_parent());
        REQUIRE(root->get_number_of_children() == 2);

        REQUIRE(clone->degree() == 2);
        REQUIRE(clone->is_root());
        REQUIRE(! clone->is_leaf());
        REQUIRE(clone->has_children());
        REQUIRE(! clone->has_parent());
        REQUIRE(clone->get_number_of_children() == 2);

        REQUIRE(root != clone);
        REQUIRE(root->get_height_parameter() == clone->get_height_parameter());
        REQUIRE(root->get_population_size_parameter() == clone->get_population_size_parameter());

        std::shared_ptr<PopulationNode> clone_child1 = clone->get_child(0);
        std::shared_ptr<PopulationNode> clone_child2 = clone->get_child(1);

        REQUIRE(root_child1->degree() == 4);
        REQUIRE(root_child1->get_number_of_children() == 3);
        REQUIRE(! root_child1->is_root());
        REQUIRE(! root_child1->is_leaf());
        REQUIRE(root_child1->has_children());
        REQUIRE(root_child1->has_parent());

        REQUIRE(clone_child1->degree() == 4);
        REQUIRE(clone_child1->get_number_of_children() == 3);
        REQUIRE(! clone_child1->is_root());
        REQUIRE(! clone_child1->is_leaf());
        REQUIRE(clone_child1->has_children());
        REQUIRE(clone_child1->has_parent());

        REQUIRE(root_child1 != clone_child1);
        REQUIRE(root_child1->get_height_parameter() == clone_child1->get_height_parameter());
        REQUIRE(root_child1->get_population_size_parameter() == clone_child1->get_population_size_parameter());

        REQUIRE(root_child2->degree() == 3);
        REQUIRE(root_child2->get_number_of_children() == 2);
        REQUIRE(! root_child2->is_root());
        REQUIRE(! root_child2->is_leaf());
        REQUIRE(root_child2->has_children());
        REQUIRE(root_child2->has_parent());

        REQUIRE(clone_child2->degree() == 3);
        REQUIRE(clone_child2->get_number_of_children() == 2);
        REQUIRE(! clone_child2->is_root());
        REQUIRE(! clone_child2->is_leaf());
        REQUIRE(clone_child2->has_children());
        REQUIRE(clone_child2->has_parent());

        REQUIRE(root_child2 != clone_child2);
        REQUIRE(root_child2->get_height_parameter() == clone_child2->get_height_parameter());
        REQUIRE(root_child2->get_population_size_parameter() == clone_child2->get_population_size_parameter());

        std::shared_ptr<PopulationNode> leaf1_clone = clone_child1->get_child(0);
        std::shared_ptr<PopulationNode> leaf2_clone = clone_child1->get_child(1);
        std::shared_ptr<PopulationNode> leaf3_clone = clone_child1->get_child(2);
        std::shared_ptr<PopulationNode> leaf4_clone = clone_child2->get_child(0);
        std::shared_ptr<PopulationNode> leaf5_clone = clone_child2->get_child(1);

        REQUIRE(leaf1->degree() == 1);
        REQUIRE(leaf1->get_number_of_children() == 0);
        REQUIRE(! leaf1->is_root());
        REQUIRE(leaf1->is_leaf());
        REQUIRE(! leaf1->has_children());
        REQUIRE(leaf1->has_parent());

        REQUIRE(leaf1_clone->degree() == 1);
        REQUIRE(leaf1_clone->get_number_of_children() == 0);
        REQUIRE(! leaf1_clone->is_root());
        REQUIRE(leaf1_clone->is_leaf());
        REQUIRE(! leaf1_clone->has_children());
        REQUIRE(leaf1_clone->has_parent());

        REQUIRE(leaf1 != leaf1_clone);
        REQUIRE(leaf1->get_height_parameter() == leaf1_clone->get_height_parameter());
        REQUIRE(leaf1->get_population_size_parameter() == leaf1_clone->get_population_size_parameter());

        REQUIRE(leaf2->degree() == 1);
        REQUIRE(leaf2->get_number_of_children() == 0);
        REQUIRE(! leaf2->is_root());
        REQUIRE(leaf2->is_leaf());
        REQUIRE(! leaf2->has_children());
        REQUIRE(leaf2->has_parent());

        REQUIRE(leaf2_clone->degree() == 1);
        REQUIRE(leaf2_clone->get_number_of_children() == 0);
        REQUIRE(! leaf2_clone->is_root());
        REQUIRE(leaf2_clone->is_leaf());
        REQUIRE(! leaf2_clone->has_children());
        REQUIRE(leaf2_clone->has_parent());

        REQUIRE(leaf2 != leaf2_clone);
        REQUIRE(leaf2->get_height_parameter() == leaf2_clone->get_height_parameter());
        REQUIRE(leaf2->get_population_size_parameter() == leaf2_clone->get_population_size_parameter());

        REQUIRE(leaf3->degree() == 1);
        REQUIRE(leaf3->get_number_of_children() == 0);
        REQUIRE(! leaf3->is_root());
        REQUIRE(leaf3->is_leaf());
        REQUIRE(! leaf3->has_children());
        REQUIRE(leaf3->has_parent());

        REQUIRE(leaf3_clone->degree() == 1);
        REQUIRE(leaf3_clone->get_number_of_children() == 0);
        REQUIRE(! leaf3_clone->is_root());
        REQUIRE(leaf3_clone->is_leaf());
        REQUIRE(! leaf3_clone->has_children());
        REQUIRE(leaf3_clone->has_parent());

        REQUIRE(leaf3 != leaf3_clone);
        REQUIRE(leaf3->get_height_parameter() == leaf3_clone->get_height_parameter());
        REQUIRE(leaf3->get_population_size_parameter() == leaf3_clone->get_population_size_parameter());

        REQUIRE(leaf4->degree() == 1);
        REQUIRE(leaf4->get_number_of_children() == 0);
        REQUIRE(! leaf4->is_root());
        REQUIRE(leaf4->is_leaf());
        REQUIRE(! leaf4->has_children());
        REQUIRE(leaf4->has_parent());

        REQUIRE(leaf4_clone->degree() == 1);
        REQUIRE(leaf4_clone->get_number_of_children() == 0);
        REQUIRE(! leaf4_clone->is_root());
        REQUIRE(leaf4_clone->is_leaf());
        REQUIRE(! leaf4_clone->has_children());
        REQUIRE(leaf4_clone->has_parent());

        REQUIRE(leaf4 != leaf4_clone);
        REQUIRE(leaf4->get_height_parameter() == leaf4_clone->get_height_parameter());
        REQUIRE(leaf4->get_population_size_parameter() == leaf4_clone->get_population_size_parameter());

        REQUIRE(leaf5->degree() == 1);
        REQUIRE(leaf5->get_number_of_children() == 0);
        REQUIRE(! leaf5->is_root());
        REQUIRE(leaf5->is_leaf());
        REQUIRE(! leaf5->has_children());
        REQUIRE(leaf5->has_parent());

        REQUIRE(leaf5_clone->degree() == 1);
        REQUIRE(leaf5_clone->get_number_of_children() == 0);
        REQUIRE(! leaf5_clone->is_root());
        REQUIRE(leaf5_clone->is_leaf());
        REQUIRE(! leaf5_clone->has_children());
        REQUIRE(leaf5_clone->has_parent());

        REQUIRE(leaf5 != leaf5_clone);
        REQUIRE(leaf5->get_height_parameter() == leaf5_clone->get_height_parameter());
        REQUIRE(leaf5->get_population_size_parameter() == leaf5_clone->get_population_size_parameter());


        REQUIRE(leaf1->get_allele_count() == 3);
        REQUIRE(leaf2->get_allele_count() == 3);
        REQUIRE(leaf3->get_allele_count() == 3);
        REQUIRE(leaf4->get_allele_count() == 3);
        REQUIRE(leaf5->get_allele_count() == 3);
        REQUIRE(root_child1->get_allele_count() == 0);
        REQUIRE(root_child2->get_allele_count() == 0);
        REQUIRE(root->get_allele_count() == 0);

        REQUIRE(leaf1_clone->get_allele_count() == 3);
        REQUIRE(leaf2_clone->get_allele_count() == 3);
        REQUIRE(leaf3_clone->get_allele_count() == 3);
        REQUIRE(leaf4_clone->get_allele_count() == 3);
        REQUIRE(leaf5_clone->get_allele_count() == 3);
        REQUIRE(clone_child1->get_allele_count() == 0);
        REQUIRE(clone_child2->get_allele_count() == 0);
        REQUIRE(clone->get_allele_count() == 0);
        
        root->resize_all();

        REQUIRE(leaf1->get_allele_count() == 3);
        REQUIRE(leaf2->get_allele_count() == 3);
        REQUIRE(leaf3->get_allele_count() == 3);
        REQUIRE(leaf4->get_allele_count() == 3);
        REQUIRE(leaf5->get_allele_count() == 3);
        REQUIRE(root_child1->get_allele_count() == 9);
        REQUIRE(root_child2->get_allele_count() == 6);
        REQUIRE(root->get_allele_count() == 15);

        REQUIRE(leaf1_clone->get_allele_count() == 3);
        REQUIRE(leaf2_clone->get_allele_count() == 3);
        REQUIRE(leaf3_clone->get_allele_count() == 3);
        REQUIRE(leaf4_clone->get_allele_count() == 3);
        REQUIRE(leaf5_clone->get_allele_count() == 3);
        REQUIRE(clone_child1->get_allele_count() == 0);
        REQUIRE(clone_child2->get_allele_count() == 0);
        REQUIRE(clone->get_allele_count() == 0);

        clone->resize_all();

        REQUIRE(leaf1_clone->get_allele_count() == 3);
        REQUIRE(leaf2_clone->get_allele_count() == 3);
        REQUIRE(leaf3_clone->get_allele_count() == 3);
        REQUIRE(leaf4_clone->get_allele_count() == 3);
        REQUIRE(leaf5_clone->get_allele_count() == 3);
        REQUIRE(clone_child1->get_allele_count() == 9);
        REQUIRE(clone_child2->get_allele_count() == 6);
        REQUIRE(clone->get_allele_count() == 15);
    }
}

TEST_CASE("Test population size setting and scaling", "[PopulationNode]") {

    SECTION("Testing setting and scaling population size") {
        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>("root", 1.0);
        std::shared_ptr<PopulationNode> root_child1 = std::make_shared<PopulationNode>("root child 1", 0.8);
        std::shared_ptr<PopulationNode> root_child2 = std::make_shared<PopulationNode>("root child 2", 0.3);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>("leaf 1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>("leaf 2", 0.0);
        std::shared_ptr<PopulationNode> leaf3 = std::make_shared<PopulationNode>("leaf 3", 0.0);
        std::shared_ptr<PopulationNode> leaf4 = std::make_shared<PopulationNode>("leaf 4", 0.0);
        std::shared_ptr<PopulationNode> leaf5 = std::make_shared<PopulationNode>("leaf 5", 0.0);

        root->add_child(root_child1);
        root->add_child(root_child2);

        leaf1->add_parent(root_child1);
        leaf2->add_parent(root_child1);
        leaf3->add_parent(root_child1);

        leaf4->add_parent(root_child2);
        leaf5->add_parent(root_child2);

        std::shared_ptr<ContinuousProbabilityDistribution> prior = std::make_shared<GammaDistribution>(1.0, 1.0);
        root->set_all_node_height_priors(prior);
        root->set_all_population_size_priors(prior);

        root->set_all_population_sizes(2.0);

        REQUIRE(root->get_population_size() == 2.0);
        REQUIRE(root_child1->get_population_size() == 2.0);
        REQUIRE(root_child2->get_population_size() == 2.0);
        REQUIRE(leaf1->get_population_size() == 2.0);
        REQUIRE(leaf2->get_population_size() == 2.0);
        REQUIRE(leaf3->get_population_size() == 2.0);
        REQUIRE(leaf4->get_population_size() == 2.0);
        REQUIRE(leaf5->get_population_size() == 2.0);

        unsigned int num_pop_sizes_scaled = root->scale_all_population_sizes(0.5);

        REQUIRE(num_pop_sizes_scaled == 8);

        REQUIRE(root->get_population_size() == 1.0);
        REQUIRE(root_child1->get_population_size() == 1.0);
        REQUIRE(root_child2->get_population_size() == 1.0);
        REQUIRE(leaf1->get_population_size() == 1.0);
        REQUIRE(leaf2->get_population_size() == 1.0);
        REQUIRE(leaf3->get_population_size() == 1.0);
        REQUIRE(leaf4->get_population_size() == 1.0);
        REQUIRE(leaf5->get_population_size() == 1.0);

        std::shared_ptr<PositiveRealParameter> size = std::make_shared<PositiveRealParameter>(prior, 4.0, false);

        root->set_all_population_size_parameters(size);

        REQUIRE(root->get_population_size() == 4.0);
        REQUIRE(root_child1->get_population_size() == 4.0);
        REQUIRE(root_child2->get_population_size() == 4.0);
        REQUIRE(leaf1->get_population_size() == 4.0);
        REQUIRE(leaf2->get_population_size() == 4.0);
        REQUIRE(leaf3->get_population_size() == 4.0);
        REQUIRE(leaf4->get_population_size() == 4.0);
        REQUIRE(leaf5->get_population_size() == 4.0);

        num_pop_sizes_scaled = root->scale_all_population_sizes(2.0);

        REQUIRE(num_pop_sizes_scaled == 1);

        REQUIRE(root->get_population_size() == 8.0);
        REQUIRE(root_child1->get_population_size() == 8.0);
        REQUIRE(root_child2->get_population_size() == 8.0);
        REQUIRE(leaf1->get_population_size() == 8.0);
        REQUIRE(leaf2->get_population_size() == 8.0);
        REQUIRE(leaf3->get_population_size() == 8.0);
        REQUIRE(leaf4->get_population_size() == 8.0);
        REQUIRE(leaf5->get_population_size() == 8.0);

        root->fix_all_population_sizes();

        num_pop_sizes_scaled = root->scale_all_population_sizes(0.5);

        REQUIRE(num_pop_sizes_scaled == 0);

        REQUIRE(root->get_population_size() == 8.0);
        REQUIRE(root_child1->get_population_size() == 8.0);
        REQUIRE(root_child2->get_population_size() == 8.0);
        REQUIRE(leaf1->get_population_size() == 8.0);
        REQUIRE(leaf2->get_population_size() == 8.0);
        REQUIRE(leaf3->get_population_size() == 8.0);
        REQUIRE(leaf4->get_population_size() == 8.0);
        REQUIRE(leaf5->get_population_size() == 8.0);
    }
}
