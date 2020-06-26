#include <typeinfo>

#include "catch.hpp"
#include "ecoevolity/node.hpp"
#include "ecoevolity/rng.hpp"

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
        REQUIRE(n.is_polytomy() == false);
        REQUIRE(n.get_node_count() == 1);
        REQUIRE(n.get_leaf_node_count() == 1);
        REQUIRE(n.get_internal_node_count() == 0);

        std::shared_ptr<Node> p = n.get_parent();
        REQUIRE(p == nullptr);
        /* REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code()); */

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range &);
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
        /* REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code()); */

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range &);
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
        /* REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code()); */

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range &);
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
        /* REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code()); */

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range &);
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
        /* REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code()); */

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range &);
        
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
        /* REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code()); */

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range &);
        
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

/*         REQUIRE_THROWS_AS(n->get_child(0), std::out_of_range &); */
        
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
        REQUIRE(p->is_polytomy() == false);
        REQUIRE(c->is_polytomy() == false);

        REQUIRE(c->get_parent() == nullptr);
        REQUIRE_THROWS_AS(p->get_child(0), std::out_of_range &);
        REQUIRE(c->is_parent(p) == false);
        REQUIRE(p->is_child(c) == false);

        c->add_parent(p);
        REQUIRE(c->has_parent() == true);
        REQUIRE(c->is_parent(p) == true);
        REQUIRE(p->is_child(c) == true);
        REQUIRE(c->get_number_of_parents() == 1);
        REQUIRE(c->degree() == 1);
        REQUIRE(p->degree() == 1);
        REQUIRE(p->is_polytomy() == false);
        REQUIRE(c->is_polytomy() == false);

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
        REQUIRE(p->is_polytomy() == false);
        REQUIRE(c->is_polytomy() == false);
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
        REQUIRE_THROWS_AS(p->get_child(0), std::out_of_range &);
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
        REQUIRE_THROWS_AS(p->remove_child(i), std::out_of_range &);
        i = 0;
        std::shared_ptr<Node> r = p->remove_child(i);
        REQUIRE(r == c);
        REQUIRE(typeid(c).hash_code() == typeid(r).hash_code());
        /* REQUIRE(typeid(*c).hash_code() == typeid(*r).hash_code()); */
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
        REQUIRE(root->is_polytomy() == false);

        leaf1->add_parent(root_child1);
        REQUIRE(root_child1->is_polytomy() == false);
        leaf2->add_parent(root_child1);
        REQUIRE(root_child1->is_polytomy() == false);
        leaf3->add_parent(root_child1);
        REQUIRE(root_child1->is_polytomy() == true);

        leaf4->add_parent(root_child2);
        leaf5->add_parent(root_child2);
        REQUIRE(root_child2->is_polytomy() == false);

        REQUIRE(root->get_polytomy_node_count() == 1);

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

        REQUIRE_THROWS_AS(n.get_bottom_pattern_probability(1, 0), std::out_of_range &);
        REQUIRE_THROWS_AS(n.get_top_pattern_probability(1, 0), std::out_of_range &);

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
        /* REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code()); */

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range &);
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

        REQUIRE_THROWS_AS(n.get_bottom_pattern_probability(1, 0), std::out_of_range &);
        REQUIRE_THROWS_AS(n.get_top_pattern_probability(1, 0), std::out_of_range &);

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
        /* REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code()); */

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range &);
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

        REQUIRE_THROWS_AS(n.get_bottom_pattern_probability(1, 0), std::out_of_range &);
        REQUIRE_THROWS_AS(n.get_top_pattern_probability(1, 0), std::out_of_range &);

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
        /* REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code()); */

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range &);
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

        REQUIRE_THROWS_AS(n.get_bottom_pattern_probability(1, 0), std::out_of_range &);
        REQUIRE_THROWS_AS(n.get_top_pattern_probability(1, 0), std::out_of_range &);

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
        /* REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code()); */

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range &);
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
        REQUIRE_THROWS_AS(n.get_bottom_pattern_probability(4, 0), std::out_of_range &);
        REQUIRE(n.get_top_pattern_probability(2, 1) == Approx(0.0));
        REQUIRE_THROWS_AS(n.get_top_pattern_probability(4, 0), std::out_of_range &);

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
        /* REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code()); */

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range &);
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
        REQUIRE_THROWS_AS(n.get_bottom_pattern_probability(4, 0), std::out_of_range &);
        REQUIRE(n.get_top_pattern_probability(2, 1) == Approx(0.0));
        REQUIRE_THROWS_AS(n.get_top_pattern_probability(4, 0), std::out_of_range &);

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
        /* REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code()); */

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range &);
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
        REQUIRE_THROWS_AS(n.get_bottom_pattern_probability(3, 0), std::out_of_range &);
        REQUIRE(n.get_top_pattern_probability(2, 1) == Approx(0.0));
        REQUIRE_THROWS_AS(n.get_top_pattern_probability(3, 0), std::out_of_range &);

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
        /* REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code()); */

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range &);
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
        REQUIRE_THROWS_AS(n.get_bottom_pattern_probability(4, 0), std::out_of_range &);
        REQUIRE(n.get_top_pattern_probability(2, 1) == Approx(0.0));
        REQUIRE_THROWS_AS(n.get_top_pattern_probability(4, 0), std::out_of_range &);

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
        /* REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code()); */

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range &);
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
        /* REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code()); */

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range &);
        
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
        /* REQUIRE(typeid(n).hash_code() == typeid(*p).hash_code()); */

        REQUIRE_THROWS_AS(n.get_child(0), std::out_of_range &);
        
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

/*         REQUIRE_THROWS_AS(n->get_child(0), std::out_of_range &); */
        
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
        REQUIRE_THROWS_AS(p->get_child(0), std::out_of_range &);
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
        REQUIRE_THROWS_AS(p->get_child(0), std::out_of_range &);
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
        REQUIRE_THROWS_AS(p->remove_child(i), std::out_of_range &);
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

TEST_CASE("Test get_leaves from one leaf", "[Node]") {

    SECTION("Testing get_leaves") {
        std::shared_ptr<Node> n = std::make_shared<Node>("root", 1.0);
        std::vector< std::shared_ptr<Node> > leaves = n->get_leaves();
        REQUIRE(leaves.size() == 1);
        REQUIRE(n == leaves.at(0));
    }
}

TEST_CASE("Test get_leaves from root with one leaf", "[Node]") {

    SECTION("Testing get_leaves") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf 1");
        root->add_child(leaf1);
        std::vector< std::shared_ptr<Node> > leaves_from_root = root->get_leaves();
        std::vector< std::shared_ptr<Node> > leaves_from_leaf = leaf1->get_leaves();
        REQUIRE(leaves_from_root.size() == 1);
        REQUIRE(leaves_from_leaf.size() == 1);
        REQUIRE(leaf1 == leaves_from_root.at(0));
        REQUIRE(leaf1 == leaves_from_leaf.at(0));
        REQUIRE(leaves_from_root == leaves_from_leaf);
    }
}

TEST_CASE("Test get_leaves from root with multiple leaves", "[Node]") {

    SECTION("Testing get_leaves") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::vector< std::shared_ptr<Node> > l;
        l.push_back(std::make_shared<Node>("leaf 1"));
        l.push_back(std::make_shared<Node>("leaf 2"));
        l.push_back(std::make_shared<Node>("leaf 3"));
        l.push_back(std::make_shared<Node>("leaf 4"));
        root->add_child(l.at(0));
        std::vector< std::shared_ptr<Node> > leaves = root->get_leaves();
        REQUIRE(leaves.size() == 1);
        REQUIRE(l.at(0) == leaves.at(0));

        root->add_child(l.at(1));
        leaves = root->get_leaves();
        REQUIRE(leaves.size() == 2);
        REQUIRE(l.at(0) == leaves.at(0));
        REQUIRE(l.at(1) == leaves.at(1));

        root->add_child(l.at(2));
        leaves = root->get_leaves();
        REQUIRE(leaves.size() == 3);
        REQUIRE(l.at(0) == leaves.at(0));
        REQUIRE(l.at(1) == leaves.at(1));
        REQUIRE(l.at(2) == leaves.at(2));

        root->add_child(l.at(3));
        leaves = root->get_leaves();
        REQUIRE(leaves.size() == 4);
        REQUIRE(l.at(0) == leaves.at(0));
        REQUIRE(l.at(1) == leaves.at(1));
        REQUIRE(l.at(2) == leaves.at(2));
        REQUIRE(l.at(3) == leaves.at(3));

        REQUIRE(l == leaves);
    }
}

TEST_CASE("Test get_leaves from root with multiple internal daughters", "[Node]") {

    SECTION("Testing get_leaves") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal 1", 1.0);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal 2", 1.0);
        std::vector< std::shared_ptr<Node> > l;
        l.push_back(std::make_shared<Node>("leaf 1"));
        l.push_back(std::make_shared<Node>("leaf 2"));
        l.push_back(std::make_shared<Node>("leaf 3"));
        l.push_back(std::make_shared<Node>("leaf 4"));
        l.push_back(std::make_shared<Node>("leaf 5"));

        root->add_child(internal1);

        internal1->add_child(l.at(0));
        internal1->add_child(l.at(1));
        internal1->add_child(l.at(2));

        root->add_child(internal2);

        internal2->add_child(l.at(3));
        internal2->add_child(l.at(4));

        std::vector< std::shared_ptr<Node> > leaves = root->get_leaves();
        REQUIRE(leaves.size() == 5);
        REQUIRE(l == leaves);
    }
}

TEST_CASE("Test get_leaf_labels from root with multiple internal daughters", "[Node]") {

    SECTION("Testing get_leaf_labels") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal 1", 1.0);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal 2", 1.0);
        std::vector< std::shared_ptr<Node> > l;
        l.push_back(std::make_shared<Node>("leaf 1"));
        l.push_back(std::make_shared<Node>("leaf 2"));
        l.push_back(std::make_shared<Node>("leaf 3"));
        l.push_back(std::make_shared<Node>("leaf 4"));
        l.push_back(std::make_shared<Node>("leaf 5"));

        root->add_child(internal1);

        internal1->add_child(l.at(0));
        internal1->add_child(l.at(1));
        internal1->add_child(l.at(2));

        root->add_child(internal2);

        internal2->add_child(l.at(3));
        internal2->add_child(l.at(4));

        std::vector<std::string> expected_leaf_labels = {
                "leaf 1",
                "leaf 2",
                "leaf 3",
                "leaf 4",
                "leaf 5"};
        std::vector<std::string> leaf_labels = root->get_leaf_labels();
        REQUIRE(leaf_labels.size() == 5);
        REQUIRE(leaf_labels == expected_leaf_labels);
    }
}


TEST_CASE("Test collapse Node to root", "[Node]") {
    SECTION("Testing collapse to root") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> root_child1 = std::make_shared<Node>("root child 1", 0.8);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf 1");
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf 2");
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf 3");
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf 4");

        root->add_child(root_child1);
        root->add_child(leaf4);
        REQUIRE(root->is_polytomy() == false);

        leaf1->add_parent(root_child1);
        REQUIRE(root_child1->is_polytomy() == false);
        leaf2->add_parent(root_child1);
        REQUIRE(root_child1->is_polytomy() == false);
        leaf3->add_parent(root_child1);
        REQUIRE(root_child1->is_polytomy() == true);

        REQUIRE(root->get_node_count() == 6);
        REQUIRE(root->get_internal_node_count() == 2);
        REQUIRE(root->get_polytomy_node_count() == 1);
        REQUIRE(root->get_leaf_node_count() == 4);
        std::vector< std::shared_ptr<Node> > expected_nodes = {root_child1};
        REQUIRE(root->get_polytomy_nodes() == expected_nodes);

        unsigned int poly_size = root_child1->collapse();
        REQUIRE(poly_size == 4);
        REQUIRE(root->is_polytomy() == true);
        REQUIRE(root->degree() == 4);
        REQUIRE(root->get_number_of_children() == 4);
        REQUIRE(root_child1->degree() == 0);
        REQUIRE(root_child1->has_parent() == false);
        REQUIRE(root_child1->has_children() == false);
        REQUIRE(leaf1->get_parent() == root);
        REQUIRE(leaf2->get_parent() == root);
        REQUIRE(leaf3->get_parent() == root);
        REQUIRE(leaf4->get_parent() == root);

        REQUIRE(root->get_internal_node_count() == 1);
        REQUIRE(root->get_polytomy_node_count() == 1);
        REQUIRE(root->get_leaf_node_count() == 4);
        REQUIRE(root->get_node_count() == 5);
        expected_nodes = {root};
        REQUIRE(root->get_polytomy_nodes() == expected_nodes);
    }
}

TEST_CASE("Test collapse Node", "[Node]") {

    SECTION("Testing collapse") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> root_child1 = std::make_shared<Node>("root child 1", 0.8);
        std::shared_ptr<Node> child2 = std::make_shared<Node>("child 2", 0.3);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf 1");
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf 2");
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf 3");
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf 4");

        root->add_child(root_child1);
        root->add_child(leaf4);
        REQUIRE(root->is_polytomy() == false);

        root_child1->add_child(leaf1);
        root_child1->add_child(child2);
        REQUIRE(root_child1->is_polytomy() == false);

        child2->add_child(leaf2);
        child2->add_child(leaf3);
        REQUIRE(child2->is_polytomy() == false);

        REQUIRE(root->get_node_count() == 7);
        REQUIRE(root->get_internal_node_count() == 3);
        REQUIRE(root->get_polytomy_node_count() == 0);
        REQUIRE(root->get_leaf_node_count() == 4);
        REQUIRE(root->get_polytomy_nodes().size() == 0);

        unsigned int poly_size = child2->collapse();
        REQUIRE(poly_size == 3);
        REQUIRE(child2->has_parent() == false);
        REQUIRE(child2->has_children() == false);

        REQUIRE(root_child1->is_polytomy() == true);
        REQUIRE(root_child1->degree() == 4);
        REQUIRE(root_child1->get_number_of_children() == 3);

        REQUIRE(leaf1->get_parent() == root_child1);
        REQUIRE(leaf2->get_parent() == root_child1);
        REQUIRE(leaf3->get_parent() == root_child1);

        REQUIRE(root->get_node_count() == 6);
        REQUIRE(root->get_internal_node_count() == 2);
        REQUIRE(root->get_polytomy_node_count() == 1);
        REQUIRE(root->get_leaf_node_count() == 4);
        std::vector< std::shared_ptr<Node> > expected_nodes = {root_child1};
        REQUIRE(root->get_polytomy_nodes() == expected_nodes);
    }
}

TEST_CASE("Test Node getters", "[Node]") {

    SECTION("Testing node getter methods") {
        std::shared_ptr<PositiveRealParameter> root_ht = std::make_shared<PositiveRealParameter>();
        root_ht->set_value(1.0);
        std::shared_ptr<PositiveRealParameter> int_ht = std::make_shared<PositiveRealParameter>();
        int_ht->set_value(0.5);
        std::shared_ptr<Node> root = std::make_shared<Node>("root");
        root->set_height_parameter(root_ht);
        std::shared_ptr<Node> root_child1 = std::make_shared<Node>("root child 1");
        root_child1->set_height_parameter(int_ht);
        std::shared_ptr<Node> root_child2 = std::make_shared<Node>("root child 2");
        root_child2->set_height_parameter(int_ht);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf 1");
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf 2");
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf 3");
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf 4");
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf 5");

        root->add_child(root_child1);
        root->add_child(root_child2);
        root_child1->add_child(leaf1);
        root_child1->add_child(leaf2);
        root_child2->add_child(leaf3);
        root_child2->add_child(leaf4);

        std::vector< std::shared_ptr<Node> > expected_mapped_nodes = {root_child1, root_child2};
        std::vector< std::shared_ptr<Node> > mapped_nodes = root->get_mapped_nodes(int_ht);
        REQUIRE(mapped_nodes.size() == expected_mapped_nodes.size());
        REQUIRE(std::is_permutation(
                    mapped_nodes.begin(), mapped_nodes.end(),
                    expected_mapped_nodes.begin()));
        REQUIRE(root->get_mapped_node_count(int_ht) == 2);

        expected_mapped_nodes = {root};
        mapped_nodes = root->get_mapped_nodes(root_ht);
        REQUIRE(mapped_nodes == expected_mapped_nodes);
        REQUIRE(root->get_mapped_node_count(root_ht) == 1);

        REQUIRE(root->get_mapped_polytomy_node_count(int_ht) == 0);
        REQUIRE(root->get_mapped_polytomy_node_count(root_ht) == 0);
        REQUIRE(root->get_mapped_polytomy_nodes(int_ht).size() == 0);
        REQUIRE(root->get_mapped_polytomy_nodes(root_ht).size() == 0);
        REQUIRE(root->get_polytomy_nodes().size() == 0);


        root_child2->add_child(leaf5);

        expected_mapped_nodes = {root_child1, root_child2};
        mapped_nodes = root->get_mapped_nodes(int_ht);
        REQUIRE(mapped_nodes.size() == expected_mapped_nodes.size());
        REQUIRE(std::is_permutation(
                    mapped_nodes.begin(), mapped_nodes.end(),
                    expected_mapped_nodes.begin()));
        REQUIRE(root->get_mapped_node_count(int_ht) == 2);

        expected_mapped_nodes = {root};
        mapped_nodes = root->get_mapped_nodes(root_ht);
        REQUIRE(mapped_nodes == expected_mapped_nodes);
        REQUIRE(root->get_mapped_node_count(root_ht) == 1);

        REQUIRE(root->get_mapped_polytomy_node_count(root_ht) == 0);
        REQUIRE(root->get_mapped_polytomy_nodes(root_ht).size() == 0);
        REQUIRE(root->get_mapped_polytomy_node_count(int_ht) == 1);
        REQUIRE(root->get_mapped_polytomy_nodes(int_ht).size() == 1);
        expected_mapped_nodes = {root_child2};
        mapped_nodes = root->get_mapped_polytomy_nodes(int_ht);
        REQUIRE(mapped_nodes == expected_mapped_nodes);
        REQUIRE(root->get_polytomy_nodes() == expected_mapped_nodes);


        unsigned int poly_size = root_child1->collapse();
        REQUIRE(poly_size == 3);

        expected_mapped_nodes = {root_child2};
        mapped_nodes = root->get_mapped_nodes(int_ht);
        REQUIRE(mapped_nodes == expected_mapped_nodes);
        REQUIRE(root->get_mapped_node_count(int_ht) == 1);

        expected_mapped_nodes = {root};
        mapped_nodes = root->get_mapped_nodes(root_ht);
        REQUIRE(mapped_nodes == expected_mapped_nodes);
        REQUIRE(root->get_mapped_node_count(root_ht) == 1);

        REQUIRE(root->get_mapped_polytomy_node_count(root_ht) == 1);
        REQUIRE(root->get_mapped_polytomy_nodes(root_ht).size() == 1);
        expected_mapped_nodes = {root};
        mapped_nodes = root->get_mapped_polytomy_nodes(root_ht);
        REQUIRE(mapped_nodes == expected_mapped_nodes);

        REQUIRE(root->get_mapped_polytomy_node_count(int_ht) == 1);
        REQUIRE(root->get_mapped_polytomy_nodes(int_ht).size() == 1);
        expected_mapped_nodes = {root_child2};
        mapped_nodes = root->get_mapped_polytomy_nodes(int_ht);
        REQUIRE(mapped_nodes == expected_mapped_nodes);

        expected_mapped_nodes = {root, root_child2};
        mapped_nodes = root->get_polytomy_nodes();
        REQUIRE(mapped_nodes.size() == expected_mapped_nodes.size());
        REQUIRE(std::is_permutation(
                    mapped_nodes.begin(), mapped_nodes.end(),
                    expected_mapped_nodes.begin()));
    }
}

TEST_CASE("Test split root polytomy", "[Node]") {
    SECTION("Testing split root polytomy") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf 1");
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf 2");
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf 3");
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf 4");

        root->add_child(leaf1);
        root->add_child(leaf2);
        root->add_child(leaf3);
        root->add_child(leaf4);

        REQUIRE(root->get_node_count() == 5);
        REQUIRE(root->get_leaf_node_count() == 4);
        REQUIRE(root->get_internal_node_count() == 1);
        REQUIRE(root->degree() == 4);
        
        std::vector< std::shared_ptr<Node> > children_to_split = {leaf1, leaf2};
        std::shared_ptr<PositiveRealParameter> new_height = std::make_shared<PositiveRealParameter>(0.5);

        RandomNumberGenerator rng(1234);
        root->split_children_from_polytomy(rng, children_to_split, new_height);

        REQUIRE(root->get_node_count() == 6);
        REQUIRE(root->get_leaf_node_count() == 4);
        REQUIRE(root->get_internal_node_count() == 2);
        REQUIRE(root->degree() == 3);

        std::shared_ptr<Node> new_node;
        for (unsigned int i = 0; i < root->get_number_of_children(); ++i) {
            if (! root->get_child(i)->is_leaf()) {
                new_node = root->get_child(i);
            }
        }

        REQUIRE(new_node->degree() == 3);
        REQUIRE(new_node->get_number_of_children() == 2);
        REQUIRE(new_node->is_parent(root));
        REQUIRE(new_node->is_child(leaf1));
        REQUIRE(new_node->is_child(leaf2));
        REQUIRE(new_node->get_height() == 0.5);
        REQUIRE(root->get_height() == 1.0);
    }
}

TEST_CASE("Test split nonroot polytomy", "[Node]") {
    SECTION("Testing split nonroot polytomy") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 2.0);
        std::shared_ptr<Node> root_child = std::make_shared<Node>("root child", 1.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf 1");
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf 2");
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf 3");
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf 4");
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf 5");

        root_child->add_child(leaf1);
        root_child->add_child(leaf2);
        root_child->add_child(leaf3);
        root_child->add_child(leaf4);

        root->add_child(root_child);
        root->add_child(leaf5);

        REQUIRE(root->get_node_count() == 7);
        REQUIRE(root->get_leaf_node_count() == 5);
        REQUIRE(root->get_internal_node_count() == 2);
        REQUIRE(root_child->degree() == 5);
        
        std::vector< std::shared_ptr<Node> > children_to_split = {leaf1, leaf2};
        std::shared_ptr<PositiveRealParameter> new_height = std::make_shared<PositiveRealParameter>(0.5);

        RandomNumberGenerator rng(123456);
        root_child->split_children_from_polytomy(rng, children_to_split, new_height);

        REQUIRE(root->get_node_count() == 8);
        REQUIRE(root->get_leaf_node_count() == 5);
        REQUIRE(root->get_internal_node_count() == 3);
        REQUIRE(root->degree() == 2);
        REQUIRE(root_child->degree() == 4);

        std::shared_ptr<Node> new_node;
        for (unsigned int i = 0; i < root_child->get_number_of_children(); ++i) {
            if (! root_child->get_child(i)->is_leaf()) {
                new_node = root_child->get_child(i);
            }
        }

        REQUIRE(root_child->degree() == 4);
        REQUIRE(new_node->degree() == 3);
        REQUIRE(new_node->get_number_of_children() == 2);
        REQUIRE(new_node->is_parent(root_child));
        REQUIRE(new_node->is_child(leaf1));
        REQUIRE(new_node->is_child(leaf2));
        REQUIRE(new_node->get_height() == 0.5);
        REQUIRE(root->get_height() == 2.0);
        REQUIRE(root_child->get_height() == 1.0);
    }
}


TEST_CASE("Test node_height_is_valid", "[Node]") {
    SECTION("Testing node_height_is_valid") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> root_child = std::make_shared<Node>("root child", 1.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf 1");
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf 2");
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf 3");

        REQUIRE(root->node_height_is_valid());
        REQUIRE(root_child->node_height_is_valid());
        REQUIRE(leaf1->node_height_is_valid());
        REQUIRE(leaf2->node_height_is_valid());
        REQUIRE(leaf3->node_height_is_valid());

        root_child->add_child(leaf1);
        root_child->add_child(leaf2);

        REQUIRE(root->node_height_is_valid());
        REQUIRE(root_child->node_height_is_valid());
        REQUIRE(leaf1->node_height_is_valid());
        REQUIRE(leaf2->node_height_is_valid());
        REQUIRE(leaf3->node_height_is_valid());

        root->add_child(root_child);
        root->add_child(leaf3);

        REQUIRE(root->node_height_is_valid());
        REQUIRE(root_child->node_height_is_valid());
        REQUIRE(leaf1->node_height_is_valid());
        REQUIRE(leaf2->node_height_is_valid());
        REQUIRE(leaf3->node_height_is_valid());

        leaf3->set_height(1.1);

        REQUIRE(! root->node_height_is_valid());
        REQUIRE(root_child->node_height_is_valid());
        REQUIRE(leaf1->node_height_is_valid());
        REQUIRE(leaf2->node_height_is_valid());
        REQUIRE(! leaf3->node_height_is_valid());

        leaf3->set_height(0.0);

        REQUIRE(root->node_height_is_valid());
        REQUIRE(root_child->node_height_is_valid());
        REQUIRE(leaf1->node_height_is_valid());
        REQUIRE(leaf2->node_height_is_valid());
        REQUIRE(leaf3->node_height_is_valid());

        leaf2->set_height(1.1);

        REQUIRE(root->node_height_is_valid());
        REQUIRE(! root_child->node_height_is_valid());
        REQUIRE(leaf1->node_height_is_valid());
        REQUIRE(! leaf2->node_height_is_valid());
        REQUIRE(leaf3->node_height_is_valid());

        leaf2->set_height(0.0);

        REQUIRE(root->node_height_is_valid());
        REQUIRE(root_child->node_height_is_valid());
        REQUIRE(leaf1->node_height_is_valid());
        REQUIRE(leaf2->node_height_is_valid());
        REQUIRE(leaf3->node_height_is_valid());

        root_child->set_height(1.1);

        REQUIRE(! root->node_height_is_valid());
        REQUIRE(! root_child->node_height_is_valid());
        REQUIRE(leaf1->node_height_is_valid());
        REQUIRE(leaf2->node_height_is_valid());
        REQUIRE(leaf3->node_height_is_valid());

        root_child->set_height(0.5);

        REQUIRE(root->node_height_is_valid());
        REQUIRE(root_child->node_height_is_valid());
        REQUIRE(leaf1->node_height_is_valid());
        REQUIRE(leaf2->node_height_is_valid());
        REQUIRE(leaf3->node_height_is_valid());
    }
}

TEST_CASE("Test node_heights_are_valid", "[Node]") {
    SECTION("Testing node_heights_are_valid") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> root_child = std::make_shared<Node>("root child", 1.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf 1");
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf 2");
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf 3");

        REQUIRE(root->node_heights_are_valid());
        REQUIRE(root_child->node_heights_are_valid());
        REQUIRE(leaf1->node_heights_are_valid());
        REQUIRE(leaf2->node_heights_are_valid());
        REQUIRE(leaf3->node_heights_are_valid());

        root_child->add_child(leaf1);
        root_child->add_child(leaf2);

        REQUIRE(root->node_heights_are_valid());
        REQUIRE(root_child->node_heights_are_valid());
        REQUIRE(leaf1->node_heights_are_valid());
        REQUIRE(leaf2->node_heights_are_valid());
        REQUIRE(leaf3->node_heights_are_valid());

        root->add_child(root_child);
        root->add_child(leaf3);

        REQUIRE(root->node_heights_are_valid());
        REQUIRE(root_child->node_heights_are_valid());
        REQUIRE(leaf1->node_heights_are_valid());
        REQUIRE(leaf2->node_heights_are_valid());
        REQUIRE(leaf3->node_heights_are_valid());

        leaf3->set_height(1.1);

        REQUIRE(! root->node_heights_are_valid());
        REQUIRE(root_child->node_heights_are_valid());
        REQUIRE(leaf1->node_heights_are_valid());
        REQUIRE(leaf2->node_heights_are_valid());
        REQUIRE(! leaf3->node_heights_are_valid());

        leaf3->set_height(0.0);

        REQUIRE(root->node_heights_are_valid());
        REQUIRE(root_child->node_heights_are_valid());
        REQUIRE(leaf1->node_heights_are_valid());
        REQUIRE(leaf2->node_heights_are_valid());
        REQUIRE(leaf3->node_heights_are_valid());

        leaf2->set_height(1.1);

        REQUIRE(! root->node_heights_are_valid());
        REQUIRE(! root_child->node_heights_are_valid());
        REQUIRE(leaf1->node_heights_are_valid());
        REQUIRE(! leaf2->node_heights_are_valid());
        REQUIRE(leaf3->node_heights_are_valid());

        leaf2->set_height(0.0);

        REQUIRE(root->node_heights_are_valid());
        REQUIRE(root_child->node_heights_are_valid());
        REQUIRE(leaf1->node_heights_are_valid());
        REQUIRE(leaf2->node_heights_are_valid());
        REQUIRE(leaf3->node_heights_are_valid());

        root_child->set_height(1.1);

        REQUIRE(! root->node_heights_are_valid());
        REQUIRE(! root_child->node_heights_are_valid());
        REQUIRE(leaf1->node_heights_are_valid());
        REQUIRE(leaf2->node_heights_are_valid());
        REQUIRE(leaf3->node_heights_are_valid());

        root_child->set_height(0.5);

        REQUIRE(root->node_heights_are_valid());
        REQUIRE(root_child->node_heights_are_valid());
        REQUIRE(leaf1->node_heights_are_valid());
        REQUIRE(leaf2->node_heights_are_valid());
        REQUIRE(leaf3->node_heights_are_valid());
    }
}

TEST_CASE("Test get_copy", "[Node]") {
    SECTION("Testing get_copy") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 1.0);
        std::shared_ptr<Node> root_child = std::make_shared<Node>("root child", 0.5);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf 1");
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf 2");
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf 3");

        root_child->add_child(leaf1);
        root_child->add_child(leaf2);

        root->add_child(root_child);
        root->add_child(leaf3);

        std::shared_ptr<Node> root_copy = root->get_copy();
        REQUIRE(root != root_copy);
        REQUIRE(root->get_label() == root_copy->get_label());
        REQUIRE(root->get_height() == root_copy->get_height());
        // node height pointers are shallow copied
        REQUIRE(root->get_height_parameter() == root_copy->get_height_parameter());

        std::shared_ptr<Node> root_child_copy = root_copy->get_child(0);
        REQUIRE(root_child->get_label() == root_child_copy->get_label());
        REQUIRE(root_child->get_height() == root_child_copy->get_height());
        REQUIRE(root_child->get_height_parameter() == root_child_copy->get_height_parameter());
        REQUIRE(root_child != root_child_copy);

        std::shared_ptr<Node> leaf1_copy = root_child_copy->get_child(0);
        REQUIRE(leaf1->get_label() == leaf1_copy->get_label());
        REQUIRE(leaf1->get_height() == leaf1_copy->get_height());
        REQUIRE(leaf1->get_height_parameter() == leaf1_copy->get_height_parameter());
        REQUIRE(leaf1 != leaf1_copy);

        std::shared_ptr<Node> leaf2_copy = root_child_copy->get_child(1);
        REQUIRE(leaf2->get_label() == leaf2_copy->get_label());
        REQUIRE(leaf2->get_height() == leaf2_copy->get_height());
        REQUIRE(leaf2->get_height_parameter() == leaf2_copy->get_height_parameter());
        REQUIRE(leaf2 != leaf2_copy);

        std::shared_ptr<Node> leaf3_copy = root_copy->get_child(1);
        REQUIRE(leaf3->get_label() == leaf3_copy->get_label());
        REQUIRE(leaf3->get_height() == leaf3_copy->get_height());
        REQUIRE(leaf3->get_height_parameter() == leaf3_copy->get_height_parameter());
        REQUIRE(leaf3 != leaf3_copy);
    }
}

TEST_CASE("Testing get_node(label)", "[Node]") {

    SECTION("Testing get_node(label)") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.1);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.07);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.05);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        internal1->add_child(leaf2);
        internal1->add_child(leaf3);
        root->add_child(internal0);
        root->add_child(internal1);
        root->add_child(leaf4);

        std::shared_ptr<Node> returned_node;

        returned_node = root->get_node("does not exist");
        REQUIRE(returned_node == nullptr);

        returned_node = leaf2->get_node("leaf2");
        REQUIRE(returned_node == leaf2);
        returned_node = internal1->get_node("leaf2");
        REQUIRE(returned_node == leaf2);

        returned_node = root->get_node("leaf0");
        REQUIRE(returned_node == leaf0);
        returned_node = root->get_node("internal0");
        REQUIRE(returned_node == internal0);
        returned_node = root->get_node("internal1");
        REQUIRE(returned_node == internal1);
        returned_node = root->get_node("leaf2");
        REQUIRE(returned_node == leaf2);
        returned_node = root->get_node("leaf4");
        REQUIRE(returned_node == leaf4);
        returned_node = root->get_node("root");
        REQUIRE(returned_node == root);
    }
}

TEST_CASE("Testing Node::get_oldest_child with leaf", "[Node]") {
    SECTION("Testing Node::get_oldest_child with leaf") {
        Node n = Node(0.03);
        REQUIRE_THROWS_AS(n.get_oldest_child(), EcoevolityError &);
    }
}

TEST_CASE("Testing Node::get_oldest_child with one child", "[Node]") {
    SECTION("Testing Node::get_oldest_child with one child") {
        std::shared_ptr<Node> n = std::make_shared<Node>("node", 0.3);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        n->add_child(leaf1);
        REQUIRE(n->get_oldest_child() == leaf1);
    }
}

TEST_CASE("Testing Node::get_oldest_child with multiple children", "[Node]") {
    SECTION("Testing Node::get_oldest_child with multiple children") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> n1 = std::make_shared<Node>("node1", 0.3);
        std::shared_ptr<Node> n2 = std::make_shared<Node>("node2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.2);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.1);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
        std::shared_ptr<Node> leaf6 = std::make_shared<Node>("leaf6", 0.0);
        std::shared_ptr<Node> leaf7 = std::make_shared<Node>("leaf7", 0.0);
        n2->add_child(leaf1);
        n2->add_child(leaf2);
        internal1->add_child(leaf3);
        internal1->add_child(leaf4);
        internal2->add_child(leaf5);
        internal2->add_child(leaf6);
        n1->add_child(internal1);
        n1->add_child(internal2);
        n1->add_child(leaf7);
        root->add_child(n1);
        root->add_child(n2);
        REQUIRE(n1->get_oldest_child() == internal1);
        REQUIRE(root->get_oldest_child() == n2);
    }
}

TEST_CASE("Testing Node::is_ancestor", "[Node]") {
    SECTION("Testing Node::is_ancestor") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> n1 = std::make_shared<Node>("node1", 0.3);
        std::shared_ptr<Node> n2 = std::make_shared<Node>("node2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.2);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.1);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
        std::shared_ptr<Node> leaf6 = std::make_shared<Node>("leaf6", 0.0);
        std::shared_ptr<Node> leaf7 = std::make_shared<Node>("leaf7", 0.0);
        n2->add_child(leaf1);
        n2->add_child(leaf2);
        internal1->add_child(leaf3);
        internal1->add_child(leaf4);
        internal2->add_child(leaf5);
        internal2->add_child(leaf6);
        n1->add_child(internal1);
        n1->add_child(internal2);
        n1->add_child(leaf7);
        root->add_child(n1);
        root->add_child(n2);
        REQUIRE(leaf1->is_ancestor(n2));
        REQUIRE(leaf1->is_ancestor(root));
        REQUIRE(! leaf1->is_ancestor(internal1));
        REQUIRE(! leaf1->is_ancestor(internal2));
        REQUIRE(! leaf1->is_ancestor(n1));
        REQUIRE(! leaf3->is_ancestor(n2));
        REQUIRE(leaf3->is_ancestor(internal1));
        REQUIRE(leaf3->is_ancestor(n1));
        REQUIRE(leaf3->is_ancestor(root));
        REQUIRE(! n1->is_ancestor(n2));
        REQUIRE(! n2->is_ancestor(n1));
        REQUIRE(! n1->is_ancestor(internal1));
        REQUIRE(internal1->is_ancestor(n1));
        REQUIRE(! n1->is_ancestor(internal2));
        REQUIRE(internal2->is_ancestor(n1));
        REQUIRE(! n2->is_ancestor(internal1));
        REQUIRE(! internal1->is_ancestor(n2));
        REQUIRE(! n2->is_ancestor(internal2));
        REQUIRE(! internal2->is_ancestor(n2));
        REQUIRE(n2->is_ancestor(root));
        REQUIRE(n1->is_ancestor(root));
        REQUIRE(internal1->is_ancestor(root));
        REQUIRE(internal2->is_ancestor(root));
        REQUIRE(! internal1->is_ancestor(leaf3));
    }
}

TEST_CASE("Testing Node::pre_order", "[xNode]") {
    SECTION("Testing Node::pre_order") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> n1 = std::make_shared<Node>("node1", 0.3);
        std::shared_ptr<Node> n2 = std::make_shared<Node>("node2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.2);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.1);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
        std::shared_ptr<Node> leaf6 = std::make_shared<Node>("leaf6", 0.0);
        std::shared_ptr<Node> leaf7 = std::make_shared<Node>("leaf7", 0.0);
        n2->add_child(leaf1);
        n2->add_child(leaf2);
        internal1->add_child(leaf3);
        internal1->add_child(leaf4);
        internal2->add_child(leaf5);
        internal2->add_child(leaf6);
        n1->add_child(internal1);
        n1->add_child(internal2);
        n1->add_child(leaf7);
        root->add_child(n1);
        root->add_child(n2);

        std::vector< std::shared_ptr<Node> > expected_nodes = {
            root,
            n1,
            internal1,
            leaf3,
            leaf4,
            internal2,
            leaf5,
            leaf6,
            leaf7,
            n2,
            leaf1,
            leaf2
        };
        std::vector< std::shared_ptr<Node> > slightly_wrong_nodes = {
            root,
            n1,
            internal1,
            leaf3,
            leaf4,
            internal2,
            leaf5,
            leaf6,
            leaf7,
            n2,
            leaf2,
            leaf1
        };
        std::vector< std::shared_ptr<Node> > nodes;
        root->pre_order(nodes);
        REQUIRE(nodes == expected_nodes);
        REQUIRE(nodes != slightly_wrong_nodes);
    }
}

TEST_CASE("Testing Node::level_order", "[xNode]") {
    SECTION("Testing Node::level_order") {
        std::shared_ptr<Node> root = std::make_shared<Node>("root", 0.5);
        std::shared_ptr<Node> n1 = std::make_shared<Node>("node1", 0.3);
        std::shared_ptr<Node> n2 = std::make_shared<Node>("node2", 0.4);
        std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.2);
        std::shared_ptr<Node> internal2 = std::make_shared<Node>("internal2", 0.1);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>("leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>("leaf5", 0.0);
        std::shared_ptr<Node> leaf6 = std::make_shared<Node>("leaf6", 0.0);
        std::shared_ptr<Node> leaf7 = std::make_shared<Node>("leaf7", 0.0);
        n2->add_child(leaf1);
        n2->add_child(leaf2);
        internal1->add_child(leaf3);
        internal1->add_child(leaf4);
        internal2->add_child(leaf5);
        internal2->add_child(leaf6);
        n1->add_child(internal1);
        n1->add_child(internal2);
        n1->add_child(leaf7);
        root->add_child(n1);
        root->add_child(n2);

        std::vector< std::shared_ptr<Node> > expected_nodes = {
            root,
            n1,
            n2,
            internal1,
            internal2,
            leaf7,
            leaf1,
            leaf2,
            leaf3,
            leaf4,
            leaf5,
            leaf6
        };
        std::vector< std::shared_ptr<Node> > slightly_wrong_nodes = {
            root,
            n1,
            n2,
            internal1,
            internal2,
            leaf7,
            leaf1,
            leaf2,
            leaf3,
            leaf4,
            leaf6,
            leaf5
        };
        std::vector< std::shared_ptr<Node> > nodes;
        root->level_order(nodes);
        REQUIRE(nodes == expected_nodes);
        REQUIRE(nodes != slightly_wrong_nodes);
    }
}
