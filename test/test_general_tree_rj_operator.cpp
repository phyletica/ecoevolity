#include "catch.hpp"
#include "ecoevolity/general_tree_operator.hpp"

#include <limits>
#include <memory>

#include "ecoevolity/tree.hpp"
#include "ecoevolity/node.hpp"
#include "ecoevolity/probability.hpp"
#include "ecoevolity/parameter.hpp"
#include "ecoevolity/stats_util.hpp"
#include "ecoevolity/rng.hpp"


inline void write_r_script(const std::vector<unsigned int> & counts,
        const std::string & path) {
    std::pair<std::string, std::string> prefix_ext = path::splitext(path::basename(path));
    std::ofstream os;
    os.open(path);
    os << "#! /usr/bin/env Rscript\n\n"
       << "plot_binom_on_hist <- function(counts, p) {\n"
       << "    k = seq(from = min(counts), to = max(counts), by = 1)\n"
       << "    n = sum(counts)\n"
       << "    binom_probs = dbinom(k, n, p)\n"
       << "    hist(counts, freq = F)\n"
       << "    lines(k, binom_probs, type = 'l')\n"
       << "}\n\n";

    os << "counts = c(";
    for (unsigned int i = 0; i < counts.size(); ++i) {
        if (i == 0) {
            os << counts.at(i);
        }
        else {
            if ((i + 1) % 1000 == 0) {
                os << ",\n        " << counts.at(i);
            }
            else {
                os << ", " << counts.at(i);
            }
        }
    }
    os << ")\n"
       << "number_of_topologies = length(counts)\n"
       << "binomial_prob = 1.0 / number_of_topologies\n"
       << "topology = seq(from = 1, to = number_of_topologies, by = 1)\n"
       << "d = data.frame(topology = topology, count = counts)\n\n";
    std::string plot_path = prefix_ext.first + "-topo-vs-count.pdf";
    os << "pdf(\"" << plot_path << "\")\n"
       << "plot(x = d$topology, y = d$count, ylim = c(0, max(counts)))\n"
       << "dev.off()\n\n";
    plot_path = prefix_ext.first + "-count-distribution.pdf";
    os << "pdf(\"" << plot_path << "\")\n"
       << "plot_binom_on_hist(counts = counts, p = binomial_prob)\n"
       << "dev.off()\n\n"
       << "chisq.test(counts)\n";
    os.close();
}

inline void write_r_script(
        const std::map< std::set< std::set<Split> >, unsigned int> & split_counts,
        const std::string & path) {
    std::vector<unsigned int> counts;
    counts.reserve(split_counts.size());
    for (auto splitset_count : split_counts) {
        counts.push_back(splitset_count.second);
    }
    write_r_script(counts, path);
}


TEST_CASE("Testing SplitLumpNodesRevJumpSampler::merge with 3 leaves",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing merge with 3 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(30);

        // MOVE FROM: ((A:0.1,B:0.1):0.2,C:0.3)
        // MOVE TO:   (A:0.3,B:0.3,C:0.3)
        //
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to merge) * pr(choose height)
        //                  =  1 * 1 = 1
        // pr(reverse move) = pr(choose to split) * pr(choose height) * pr(choosing nodes to move down) * pr(new height)
        //                  = 1 * 1 * 1/3 * 1/0.3
        //                  = 1/(3*0.3)
        // HR = 1/(3*0.3) / 1 = 1/(3*0.3)
        
        unsigned int nsamples = 50;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.3;
            std::shared_ptr<Node> root = std::make_shared<Node>("root", root_ht);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.1);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

            internal0->add_child(leaf0);
            internal0->add_child(leaf1);
            root->add_child(internal0);
            root->add_child(leaf2);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_number_of_node_heights() == 2);
            REQUIRE(tree.get_number_of_splittable_heights() == 0);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / pow(root_ht, 1))).epsilon(1e-8));

            double ln_hastings = op.propose(rng,
                    &tree);
            double exp_ln_hastings = std::log(1.0 / (3.0 * 0.3));
            REQUIRE(ln_hastings == Approx(exp_ln_hastings).epsilon(1e-8));
            REQUIRE(tree.get_number_of_node_heights() == 1);
            REQUIRE(tree.get_number_of_splittable_heights() == 1);
            REQUIRE(tree.get_height(0) == 0.3);
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / pow(root_ht, 0))).epsilon(1e-8));
        }
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler::split with 3 leaves",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing split with 3 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(30);

        // MOVE FROM: (A:0.3,B:0.3,C:0.3)
        // MOVE TO:   ((A:<0.3,B:<0.3):<0.3,C:0.3)
        //            OR
        // MOVE TO:   ((A:<0.3,C:<0.3):<0.3,B:0.3)
        //            OR
        // MOVE TO:   ((B:<0.3,C:<0.3):<0.3,C:0.3)
        //
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to split) * pr(choose height) * pr(choosing nodes to move down) * pr(new height)
        //                  = 1 * 1 * 1/3 * 1/0.3
        //                  = 1/(3*0.3)
        // pr(reverse move) = pr(choose to merge) * pr(choose height)
        //                  =  1 * 1 = 1
        // HR = 1 / 1/(3*0.3) = (3*0.3)
        
        unsigned int nsamples = 10000;
        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_12 = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.3;
            std::shared_ptr<Node> root = std::make_shared<Node>("root", root_ht);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

            root->add_child(leaf0);
            root->add_child(leaf1);
            root->add_child(leaf2);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_number_of_node_heights() == 1);
            REQUIRE(tree.get_number_of_splittable_heights() == 1);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / pow(root_ht, 0))).epsilon(1e-8));

            double ln_hastings = op.propose(rng,
                    &tree);
            double exp_ln_hastings = std::log(3.0 * 0.3);
            REQUIRE(ln_hastings == Approx(exp_ln_hastings).epsilon(1e-8));
            REQUIRE(tree.get_number_of_node_heights() == 2);
            REQUIRE(tree.get_number_of_splittable_heights() == 0);
            REQUIRE(tree.get_height(1) == 0.3);
            REQUIRE(tree.get_height(0) < 0.3);
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / pow(root_ht, 1))).epsilon(1e-8));
            if (tree.get_root_ptr()->is_child("leaf2")) {
                ++count_01;
            }
            if (tree.get_root_ptr()->is_child("leaf1")) {
                ++count_02;
            }
            if (tree.get_root_ptr()->is_child("leaf0")) {
                ++count_12;
            }
        }
        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        double eps = 0.01;
        std::cout << "Freq of ((0,1),2): " << count_01 / (double)nsamples << "\n";
        std::cout << "Freq of ((0,2),1): " << count_02 / (double)nsamples << "\n";
        std::cout << "Freq of ((1,2),0): " << count_12 / (double)nsamples << "\n";
        REQUIRE(count_01 / (double)nsamples == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(count_02 / (double)nsamples == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(count_12 / (double)nsamples == Approx(1.0/3.0).epsilon(eps));
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler with 3 leaves and fixed root",
        "[xSplitLumpNodesRevJumpSampler]") {

    SECTION("Testing 3 leaves with fixed root") {
        RandomNumberGenerator rng = RandomNumberGenerator(18);

        double root_ht = 0.2;
        std::shared_ptr<Node> root = std::make_shared<Node>("root", root_ht);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.1);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        SplitLumpNodesRevJumpSampler<Node> op;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_height_summary;

        unsigned int count_012 = 0;
        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_12 = 0;

        unsigned int niterations = 2000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            // std::cout << tree.to_parentheses() << "\n";
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                if (tree.get_root_ptr()->get_number_of_children() == 3) {
                    ++count_012;
                }
                else {
                    if (tree.get_root_ptr()->is_child("leaf0")) {
                        ++count_12;
                    }
                    if (tree.get_root_ptr()->is_child("leaf1")) {
                        ++count_02;
                    }
                    if (tree.get_root_ptr()->is_child("leaf2")) {
                        ++count_01;
                    }
                    internal_height_summary.add_sample(tree.get_height(0));
                }
                REQUIRE(tree.get_root_height() == root_ht);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        REQUIRE((count_01 + count_02 + count_12 + count_012) == nsamples);
        REQUIRE(internal_height_summary.sample_size() == (count_01 + count_02 + count_12));

        double freq_012 = count_012 / (double)nsamples;
        double freq_01 = count_01 / (double)nsamples;
        double freq_02 = count_02 / (double)nsamples;
        double freq_12 = count_12 / (double)nsamples;
        std::cout << "Freq of (0,1,2): " << freq_012 << "\n";
        std::cout << "Freq of ((0,1),2): " << freq_01 << "\n";
        std::cout << "Freq of ((0,2),1): " << freq_02 << "\n";
        std::cout << "Freq of ((1,2),0): " << freq_12 << "\n";

        std::vector<unsigned int> counts {count_012, count_01, count_02, count_12};
        write_r_script(counts, "../3-leaf-general-tree-test.r");

        double eps = 0.001;

        REQUIRE(freq_012 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_01 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_02 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_12 == Approx(0.25).epsilon(eps));
        
        UniformDistribution prior(0.0, root_ht);

        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler with 3 leaves, fixed root and operate_plus",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing 3 leaves, fixed root and operate_plus") {
        RandomNumberGenerator rng = RandomNumberGenerator(19191919);

        double root_ht = 0.2;
        std::shared_ptr<Node> root = std::make_shared<Node>("root", root_ht);
        std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.1);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);

        internal0->add_child(leaf0);
        internal0->add_child(leaf1);
        root->add_child(internal0);
        root->add_child(leaf2);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        SplitLumpNodesRevJumpSampler<Node> op;
        std::shared_ptr< NodeHeightScaler<Node> > node_height_op = std::make_shared<NodeHeightScaler<Node> >();
        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate< BaseTree<Node> > > > other_ops;
        other_ops.push_back(node_height_op);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_height_summary;

        unsigned int count_012 = 0;
        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_12 = 0;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            // std::cout << tree.to_parentheses() << "\n";
            op.operate_plus(rng, &tree, other_ops, 1, 2);
            if ((i + 1) % sample_freq == 0) {
                if (tree.get_root_ptr()->get_number_of_children() == 3) {
                    ++count_012;
                }
                else {
                    if (tree.get_root_ptr()->is_child("leaf0")) {
                        ++count_12;
                    }
                    if (tree.get_root_ptr()->is_child("leaf1")) {
                        ++count_02;
                    }
                    if (tree.get_root_ptr()->is_child("leaf2")) {
                        ++count_01;
                    }
                    internal_height_summary.add_sample(tree.get_height(0));
                }
                REQUIRE(tree.get_root_height() == root_ht);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << node_height_op->to_string();

        REQUIRE((count_01 + count_02 + count_12 + count_012) == nsamples);
        REQUIRE(internal_height_summary.sample_size() == (count_01 + count_02 + count_12));

        double freq_012 = count_012 / (double)nsamples;
        double freq_01 = count_01 / (double)nsamples;
        double freq_02 = count_02 / (double)nsamples;
        double freq_12 = count_12 / (double)nsamples;
        std::cout << "Freq of (0,1,2): " << freq_012 << "\n";
        std::cout << "Freq of ((0,1),2): " << freq_01 << "\n";
        std::cout << "Freq of ((0,2),1): " << freq_02 << "\n";
        std::cout << "Freq of ((1,2),0): " << freq_12 << "\n";

        double eps = 0.001;

        REQUIRE(freq_012 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_01 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_02 == Approx(0.25).epsilon(eps));
        REQUIRE(freq_12 == Approx(0.25).epsilon(eps));
        
        UniformDistribution prior(0.0, root_ht);

        REQUIRE(internal_height_summary.mean() == Approx(prior.get_mean()).epsilon(eps));
        REQUIRE(internal_height_summary.variance() == Approx(prior.get_variance()).epsilon(eps));
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler with 4 leaves and fixed root",
        "[xxSplitLumpNodesRevJumpSampler]") {

    SECTION("Testing 4 leaves with fixed root") {
        RandomNumberGenerator rng = RandomNumberGenerator(20);

        double root_ht = 0.5;
        std::shared_ptr<Node> root = std::make_shared<Node>(4, "root", root_ht);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);

        root->add_child(leaf0);
        root->add_child(leaf1);
        root->add_child(leaf2);
        root->add_child(leaf3);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        SplitLumpNodesRevJumpSampler<Node> op;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_height_summary;

        unsigned int count_0123 = 0;
        unsigned int count_0_ = 0;
        unsigned int count_1_ = 0;
        unsigned int count_2_ = 0;
        unsigned int count_3_ = 0;
        unsigned int count_01_ = 0;
        unsigned int count_02_ = 0;
        unsigned int count_03_ = 0;
        unsigned int count_12_ = 0;
        unsigned int count_13_ = 0;
        unsigned int count_23_ = 0;
        unsigned int count_01_23 = 0;
        unsigned int count_02_13 = 0;
        unsigned int count_03_12 = 0;
        unsigned int count_gen_01_23 = 0;
        unsigned int count_gen_02_13 = 0;
        unsigned int count_gen_03_12 = 0;
        unsigned int count_gen_01_2_3 = 0;
        unsigned int count_gen_01_3_2 = 0;
        unsigned int count_gen_02_1_3 = 0;
        unsigned int count_gen_02_3_1 = 0;
        unsigned int count_gen_03_1_2 = 0;
        unsigned int count_gen_03_2_1 = 0;
        unsigned int count_gen_12_0_3 = 0;
        unsigned int count_gen_12_3_0 = 0;
        unsigned int count_gen_13_0_2 = 0;
        unsigned int count_gen_13_2_0 = 0;
        unsigned int count_gen_23_0_1 = 0;
        unsigned int count_gen_23_1_0 = 0;
        unsigned int count_3_heights = 0;
        unsigned int count_2_heights = 0;
        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int niterations = 5000000;
        unsigned int sample_freq = 20;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {

                std::set< std::set<Split> > splits = tree.get_splits(false);
                if (split_counts.count(splits) > 0) {
                    ++split_counts[splits];
                }
                else {
                    split_counts[splits] = 1;
                }

                if (tree.get_number_of_node_heights() == 1) {
                    ++count_0123;
                    REQUIRE(tree.get_root_ptr()->get_number_of_children() == 4);
                }
                else if (tree.get_number_of_node_heights() == 2) {
                    ++count_2_heights;
                    if (tree.get_root_ptr()->get_number_of_children() == 2) {
                        if (tree.get_root_ptr()->get_child(0)->is_leaf() ||
                                tree.get_root_ptr()->get_child(1)->is_leaf()) {
                            if(tree.get_root_ptr()->is_child("leaf0")) {
                                ++count_0_;
                            }
                            else if(tree.get_root_ptr()->is_child("leaf1")) {
                                ++count_1_;
                            }
                            else if(tree.get_root_ptr()->is_child("leaf2")) {
                                ++count_2_;
                            }
                            else if(tree.get_root_ptr()->is_child("leaf3")) {
                                ++count_3_;
                            }
                            else {
                                REQUIRE(0 == 1);
                            }
                        }
                        else {
                            if (
                                    (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(0)->is_child("leaf1"))
                                    ||
                                    (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                               ) {
                                ++count_01_23;
                            }
                            else if (
                                    (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(0)->is_child("leaf2"))
                                    ||
                                    (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                               ) {
                                ++count_02_13;
                            }
                            else if (
                                    (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(0)->is_child("leaf3"))
                                    ||
                                    (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                               ) {
                                ++count_03_12;
                            }
                            else {
                                REQUIRE(0 == 1);
                            }
                        }
                    }
                    else if (tree.get_root_ptr()->get_number_of_children() == 3) {
                        if (
                                tree.get_root_ptr()->is_child("leaf0") &&
                                tree.get_root_ptr()->is_child("leaf1")
                            )
                        {
                            ++count_23_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf0") &&
                                tree.get_root_ptr()->is_child("leaf2")
                            )
                        {
                            ++count_13_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf0") &&
                                tree.get_root_ptr()->is_child("leaf3")
                            )
                        {
                            ++count_12_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf1") &&
                                tree.get_root_ptr()->is_child("leaf2")
                            )
                        {
                            ++count_03_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf1") &&
                                tree.get_root_ptr()->is_child("leaf3")
                            )
                        {
                            ++count_02_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf2") &&
                                tree.get_root_ptr()->is_child("leaf3")
                            )
                        {
                            ++count_01_;

                        }
                        else {
                            REQUIRE(0 == 1);
                        }
                    }
                    else {
                        REQUIRE(0 == 1);
                    }
                }
                else if (tree.get_number_of_node_heights() == 3) {
                    ++count_3_heights;
                    if ((! tree.get_root_ptr()->get_child(0)->is_leaf()) &&
                        (! tree.get_root_ptr()->get_child(1)->is_leaf())) {
                        if (
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(0)->is_child("leaf1"))
                                ||
                                (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                           ) {
                            ++count_gen_01_23;
                        }
                        else if (
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(0)->is_child("leaf2"))
                                ||
                                (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                           ) {
                            ++count_gen_02_13;
                        }
                        else if (
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(0)->is_child("leaf3"))
                                ||
                                (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                           ) {
                            ++count_gen_03_12;
                        }
                        else {
                            REQUIRE(0 == 1);
                        }
                    }
                    else {
                        // general ladderized topology
                        if (tree.get_root_ptr()->is_child("leaf3") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf2") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                            )
                        {
                            ++count_gen_01_2_3;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf2") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf3") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                            )
                        {
                            ++count_gen_01_3_2;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf3") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf1") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                            )
                        {
                            ++count_gen_02_1_3;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf1") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf3") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                            )
                        {
                            ++count_gen_02_3_1;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf2") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf1") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                            )
                        {
                            ++count_gen_03_1_2;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf1") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf2") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                            )
                        {
                            ++count_gen_03_2_1;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf3") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf0"))
                            )
                        {
                            ++count_gen_12_0_3;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf0") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf3") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                            )
                        {
                            ++count_gen_12_3_0;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf2") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf0"))
                            )
                        {
                            ++count_gen_13_0_2;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf0") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf2") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                            )
                        {
                            ++count_gen_13_2_0;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf1") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf0"))
                            )
                        {
                            ++count_gen_23_0_1;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf0") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf1") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                            )
                        {
                            ++count_gen_23_1_0;
                        }
                        else {
                            REQUIRE(0 == 1);
                        }
                    }
                }
                else {
                    REQUIRE(0 == 1);
                }
                REQUIRE(tree.get_root_height() == root_ht);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        REQUIRE((count_0123 + count_2_heights + count_3_heights) == nsamples);
        REQUIRE((count_01_23 +
                count_02_13 +
                count_03_12 +
                count_0_ +
                count_1_ +
                count_2_ +
                count_3_ +
                count_01_ +
                count_02_ +
                count_03_ +
                count_12_ +
                count_13_ +
                count_23_) == count_2_heights);
        REQUIRE((count_gen_01_23 +
                count_gen_02_13 +
                count_gen_03_12 +
                count_gen_01_2_3 +
                count_gen_01_3_2 +
                count_gen_02_1_3 +
                count_gen_02_3_1 +
                count_gen_03_1_2 +
                count_gen_03_2_1 +
                count_gen_12_0_3 +
                count_gen_12_3_0 +
                count_gen_13_0_2 +
                count_gen_13_2_0 +
                count_gen_23_0_1 +
                count_gen_23_1_0) == count_3_heights);

        double freq_2_heights = count_2_heights / (double)nsamples;
        double freq_3_heights = count_3_heights / (double)nsamples;
        double freq_0123 = count_0123 / (double)nsamples;
        double freq_01_23 = count_01_23 / (double)nsamples;
        double freq_02_13 = count_02_13 / (double)nsamples;
        double freq_03_12 = count_03_12 / (double)nsamples;
        double freq_gen_01_23 = count_gen_01_23 / (double)nsamples;
        double freq_gen_02_13 = count_gen_02_13 / (double)nsamples;
        double freq_gen_03_12 = count_gen_03_12 / (double)nsamples;
        double freq_0_ = count_0_ / (double)nsamples;
        double freq_1_ = count_1_ / (double)nsamples;
        double freq_2_ = count_2_ / (double)nsamples;
        double freq_3_ = count_3_ / (double)nsamples;
        double freq_01_ = count_01_ / (double)nsamples;
        double freq_02_ = count_02_ / (double)nsamples;
        double freq_03_ = count_03_ / (double)nsamples;
        double freq_12_ = count_12_ / (double)nsamples;
        double freq_13_ = count_13_ / (double)nsamples;
        double freq_23_ = count_23_ / (double)nsamples;
        double freq_gen_01_2_3 = count_gen_01_2_3 / (double)nsamples;
        double freq_gen_01_3_2 = count_gen_01_3_2 / (double)nsamples;
        double freq_gen_02_1_3 = count_gen_02_1_3 / (double)nsamples;
        double freq_gen_02_3_1 = count_gen_02_3_1 / (double)nsamples;
        double freq_gen_03_1_2 = count_gen_03_1_2 / (double)nsamples;
        double freq_gen_03_2_1 = count_gen_03_2_1 / (double)nsamples;
        double freq_gen_12_0_3 = count_gen_12_0_3 / (double)nsamples;
        double freq_gen_12_3_0 = count_gen_12_3_0 / (double)nsamples;
        double freq_gen_13_0_2 = count_gen_13_0_2 / (double)nsamples;
        double freq_gen_13_2_0 = count_gen_13_2_0 / (double)nsamples;
        double freq_gen_23_0_1 = count_gen_23_0_1 / (double)nsamples;
        double freq_gen_23_1_0 = count_gen_23_1_0 / (double)nsamples;

        double exp_freq = 1.0/29.0;

        std::cout << "Freq of (0,1,2,3): " << freq_0123 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((0,1),2,3): " << freq_01_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((0,2),1,3): " << freq_02_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((0,3),1,2): " << freq_03_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((1,2),0,3): " << freq_12_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((1,3),0,2): " << freq_13_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((2,3),0,1): " << freq_23_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of (0,(1,2,3)): " << freq_0_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of (1,(0,2,3)): " << freq_1_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of (2,(1,0,3)): " << freq_2_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of (3,(1,2,0)): " << freq_3_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of shared ((0,1),(2,3)): " << freq_01_23 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of shared ((0,2),(1,3)): " << freq_02_13 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of shared ((0,3),(1,2)): " << freq_03_12 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen ((0,1),(2,3)): " << freq_gen_01_23 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen ((0,2),(1,3)): " << freq_gen_02_13 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen ((0,3),(1,2)): " << freq_gen_03_12 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,1),2),3): " << freq_gen_01_2_3 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,1),3),2): " << freq_gen_01_3_2 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,2),1),3): " << freq_gen_02_1_3 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,2),3),1): " << freq_gen_02_3_1 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,3),1),2): " << freq_gen_03_1_2 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,3),2),1): " << freq_gen_03_2_1 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((1,2),0),3): " << freq_gen_12_0_3 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((1,2),3),0): " << freq_gen_12_3_0 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((1,3),0),2): " << freq_gen_13_0_2 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((1,3),2),0): " << freq_gen_13_2_0 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((2,3),0),1): " << freq_gen_23_0_1 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((2,3),1),0): " << freq_gen_23_1_0 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of 2 heights: " << freq_2_heights << " (expected " << 13 * exp_freq << ")\n";
        std::cout << "Freq of 3 heights: " << freq_3_heights << " (expected " << 15 * exp_freq << ")\n";

        write_r_script(split_counts, "../4-leaf-general-tree-test.r");

        double eps = 0.001;

        REQUIRE(freq_2_heights == Approx(13 * exp_freq).epsilon(eps));
        REQUIRE(freq_3_heights == Approx(15 * exp_freq).epsilon(eps));

        REQUIRE(freq_0123 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_01_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_02_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_03_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_12_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_13_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_23_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_0_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_1_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_2_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_3_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_01_23 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_02_13 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_03_12 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_01_23 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_02_13 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_03_12 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_01_2_3 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_01_3_2 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_02_1_3 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_02_3_1 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_03_1_2 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_03_2_1 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_12_0_3 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_12_3_0 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_13_0_2 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_13_2_0 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_23_0_1 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_23_1_0 == Approx(exp_freq).epsilon(eps));
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler with 4 leaves, fixed root, and operate_plus",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing 4 leaves with fixed root and operate_plus") {
        RandomNumberGenerator rng = RandomNumberGenerator(21);

        double root_ht = 0.5;
        std::shared_ptr<Node> root = std::make_shared<Node>("root", root_ht);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);

        root->add_child(leaf0);
        root->add_child(leaf1);
        root->add_child(leaf2);
        root->add_child(leaf3);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        SplitLumpNodesRevJumpSampler<Node> op;
        std::shared_ptr< NodeHeightScaler<Node> > node_height_op = std::make_shared<NodeHeightScaler<Node> >();
        /* std::shared_ptr< NeighborHeightNodeSwap<Node> > node_swap_op = std::make_shared<NeighborHeightNodeSwap<Node> >(); */
        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate< BaseTree<Node> > > > other_ops;
        other_ops.push_back(node_height_op);
        /* other_ops.push_back(node_swap_op); */

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        SampleSummarizer<double> internal_height_summary;

        unsigned int count_0123 = 0;
        unsigned int count_0_ = 0;
        unsigned int count_1_ = 0;
        unsigned int count_2_ = 0;
        unsigned int count_3_ = 0;
        unsigned int count_01_ = 0;
        unsigned int count_02_ = 0;
        unsigned int count_03_ = 0;
        unsigned int count_12_ = 0;
        unsigned int count_13_ = 0;
        unsigned int count_23_ = 0;
        unsigned int count_01_23 = 0;
        unsigned int count_02_13 = 0;
        unsigned int count_03_12 = 0;
        unsigned int count_gen_01_23 = 0;
        unsigned int count_gen_02_13 = 0;
        unsigned int count_gen_03_12 = 0;
        unsigned int count_gen_01_2_3 = 0;
        unsigned int count_gen_01_3_2 = 0;
        unsigned int count_gen_02_1_3 = 0;
        unsigned int count_gen_02_3_1 = 0;
        unsigned int count_gen_03_1_2 = 0;
        unsigned int count_gen_03_2_1 = 0;
        unsigned int count_gen_12_0_3 = 0;
        unsigned int count_gen_12_3_0 = 0;
        unsigned int count_gen_13_0_2 = 0;
        unsigned int count_gen_13_2_0 = 0;
        unsigned int count_gen_23_0_1 = 0;
        unsigned int count_gen_23_1_0 = 0;
        unsigned int count_3_heights = 0;
        unsigned int count_2_heights = 0;

        unsigned int niterations = 1000000;
        unsigned int sample_freq = 10;
        unsigned int nsamples = niterations / sample_freq;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate_plus(rng, &tree, other_ops, 1, 2, 2);
            if ((i + 1) % sample_freq == 0) {
                /* std::cout << "prior: " << tree.get_log_prior_density_value() << "\n"; */
                if (tree.get_number_of_node_heights() == 1) {
                    ++count_0123;
                    REQUIRE(tree.get_root_ptr()->get_number_of_children() == 4);
                }
                else if (tree.get_number_of_node_heights() == 2) {
                    ++count_2_heights;
                    if (tree.get_root_ptr()->get_number_of_children() == 2) {
                        if (tree.get_root_ptr()->get_child(0)->is_leaf() ||
                                tree.get_root_ptr()->get_child(1)->is_leaf()) {
                            if(tree.get_root_ptr()->is_child("leaf0")) {
                                ++count_0_;
                            }
                            else if(tree.get_root_ptr()->is_child("leaf1")) {
                                ++count_1_;
                            }
                            else if(tree.get_root_ptr()->is_child("leaf2")) {
                                ++count_2_;
                            }
                            else if(tree.get_root_ptr()->is_child("leaf3")) {
                                ++count_3_;
                            }
                            else {
                                REQUIRE(0 == 1);
                            }
                        }
                        else {
                            if (
                                    (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(0)->is_child("leaf1"))
                                    ||
                                    (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                               ) {
                                ++count_01_23;
                            }
                            else if (
                                    (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(0)->is_child("leaf2"))
                                    ||
                                    (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                               ) {
                                ++count_02_13;
                            }
                            else if (
                                    (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(0)->is_child("leaf3"))
                                    ||
                                    (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                    tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                               ) {
                                ++count_03_12;
                            }
                            else {
                                REQUIRE(0 == 1);
                            }
                        }
                    }
                    else if (tree.get_root_ptr()->get_number_of_children() == 3) {
                        if (
                                tree.get_root_ptr()->is_child("leaf0") &&
                                tree.get_root_ptr()->is_child("leaf1")
                            )
                        {
                            ++count_23_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf0") &&
                                tree.get_root_ptr()->is_child("leaf2")
                            )
                        {
                            ++count_13_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf0") &&
                                tree.get_root_ptr()->is_child("leaf3")
                            )
                        {
                            ++count_12_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf1") &&
                                tree.get_root_ptr()->is_child("leaf2")
                            )
                        {
                            ++count_03_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf1") &&
                                tree.get_root_ptr()->is_child("leaf3")
                            )
                        {
                            ++count_02_;

                        }
                        else if (
                                tree.get_root_ptr()->is_child("leaf2") &&
                                tree.get_root_ptr()->is_child("leaf3")
                            )
                        {
                            ++count_01_;

                        }
                        else {
                            REQUIRE(0 == 1);
                        }
                    }
                    else {
                        REQUIRE(0 == 1);
                    }
                }
                else if (tree.get_number_of_node_heights() == 3) {
                    ++count_3_heights;
                    if ((! tree.get_root_ptr()->get_child(0)->is_leaf()) &&
                        (! tree.get_root_ptr()->get_child(1)->is_leaf())) {
                        if (
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(0)->is_child("leaf1"))
                                ||
                                (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                           ) {
                            ++count_gen_01_23;
                        }
                        else if (
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(0)->is_child("leaf2"))
                                ||
                                (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                           ) {
                            ++count_gen_02_13;
                        }
                        else if (
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(0)->is_child("leaf3"))
                                ||
                                (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                                tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                           ) {
                            ++count_gen_03_12;
                        }
                        else {
                            REQUIRE(0 == 1);
                        }
                    }
                    else {
                        // general ladderized topology
                        if (tree.get_root_ptr()->is_child("leaf3") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf2") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                            )
                        {
                            ++count_gen_01_2_3;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf2") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf3") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                            )
                        {
                            ++count_gen_01_3_2;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf3") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf1") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                            )
                        {
                            ++count_gen_02_1_3;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf1") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf3") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                            )
                        {
                            ++count_gen_02_3_1;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf2") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf1") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                            )
                        {
                            ++count_gen_03_1_2;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf1") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf2") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                            )
                        {
                            ++count_gen_03_2_1;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf3") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf0"))
                            )
                        {
                            ++count_gen_12_0_3;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf0") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf3") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                            )
                        {
                            ++count_gen_12_3_0;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf2") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf0"))
                            )
                        {
                            ++count_gen_13_0_2;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf0") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf2") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                            )
                        {
                            ++count_gen_13_2_0;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf1") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf0") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf0"))
                            )
                        {
                            ++count_gen_23_0_1;
                        }
                        else if (tree.get_root_ptr()->is_child("leaf0") && 
                                (tree.get_root_ptr()->get_child(0)->is_child("leaf1") ||
                                 tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                            )
                        {
                            ++count_gen_23_1_0;
                        }
                        else {
                            REQUIRE(0 == 1);
                        }
                    }
                }
                else {
                    REQUIRE(0 == 1);
                }
                REQUIRE(tree.get_root_height() == root_ht);
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << node_height_op->to_string();
        /* std::cout << node_swap_op->to_string(); */


        REQUIRE((count_0123 + count_2_heights + count_3_heights) == nsamples);
        REQUIRE((count_01_23 +
                count_02_13 +
                count_03_12 +
                count_0_ +
                count_1_ +
                count_2_ +
                count_3_ +
                count_01_ +
                count_02_ +
                count_03_ +
                count_12_ +
                count_13_ +
                count_23_) == count_2_heights);
        REQUIRE((count_gen_01_23 +
                count_gen_02_13 +
                count_gen_03_12 +
                count_gen_01_2_3 +
                count_gen_01_3_2 +
                count_gen_02_1_3 +
                count_gen_02_3_1 +
                count_gen_03_1_2 +
                count_gen_03_2_1 +
                count_gen_12_0_3 +
                count_gen_12_3_0 +
                count_gen_13_0_2 +
                count_gen_13_2_0 +
                count_gen_23_0_1 +
                count_gen_23_1_0) == count_3_heights);

        double freq_2_heights = count_2_heights / (double)nsamples;
        double freq_3_heights = count_3_heights / (double)nsamples;
        double freq_0123 = count_0123 / (double)nsamples;
        double freq_01_23 = count_01_23 / (double)nsamples;
        double freq_02_13 = count_02_13 / (double)nsamples;
        double freq_03_12 = count_03_12 / (double)nsamples;
        double freq_gen_01_23 = count_gen_01_23 / (double)nsamples;
        double freq_gen_02_13 = count_gen_02_13 / (double)nsamples;
        double freq_gen_03_12 = count_gen_03_12 / (double)nsamples;
        double freq_0_ = count_0_ / (double)nsamples;
        double freq_1_ = count_1_ / (double)nsamples;
        double freq_2_ = count_2_ / (double)nsamples;
        double freq_3_ = count_3_ / (double)nsamples;
        double freq_01_ = count_01_ / (double)nsamples;
        double freq_02_ = count_02_ / (double)nsamples;
        double freq_03_ = count_03_ / (double)nsamples;
        double freq_12_ = count_12_ / (double)nsamples;
        double freq_13_ = count_13_ / (double)nsamples;
        double freq_23_ = count_23_ / (double)nsamples;
        double freq_gen_01_2_3 = count_gen_01_2_3 / (double)nsamples;
        double freq_gen_01_3_2 = count_gen_01_3_2 / (double)nsamples;
        double freq_gen_02_1_3 = count_gen_02_1_3 / (double)nsamples;
        double freq_gen_02_3_1 = count_gen_02_3_1 / (double)nsamples;
        double freq_gen_03_1_2 = count_gen_03_1_2 / (double)nsamples;
        double freq_gen_03_2_1 = count_gen_03_2_1 / (double)nsamples;
        double freq_gen_12_0_3 = count_gen_12_0_3 / (double)nsamples;
        double freq_gen_12_3_0 = count_gen_12_3_0 / (double)nsamples;
        double freq_gen_13_0_2 = count_gen_13_0_2 / (double)nsamples;
        double freq_gen_13_2_0 = count_gen_13_2_0 / (double)nsamples;
        double freq_gen_23_0_1 = count_gen_23_0_1 / (double)nsamples;
        double freq_gen_23_1_0 = count_gen_23_1_0 / (double)nsamples;

        double exp_freq = 1.0/29.0;

        std::cout << "Freq of (0,1,2,3): " << freq_0123 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((0,1),2,3): " << freq_01_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((0,2),1,3): " << freq_02_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((0,3),1,2): " << freq_03_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((1,2),0,3): " << freq_12_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((1,3),0,2): " << freq_13_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of ((2,3),0,1): " << freq_23_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of (0,(1,2,3)): " << freq_0_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of (1,(0,2,3)): " << freq_1_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of (2,(1,0,3)): " << freq_2_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of (3,(1,2,0)): " << freq_3_ << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of shared ((0,1),(2,3)): " << freq_01_23 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of shared ((0,2),(1,3)): " << freq_02_13 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of shared ((0,3),(1,2)): " << freq_03_12 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen ((0,1),(2,3)): " << freq_gen_01_23 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen ((0,2),(1,3)): " << freq_gen_02_13 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen ((0,3),(1,2)): " << freq_gen_03_12 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,1),2),3): " << freq_gen_01_2_3 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,1),3),2): " << freq_gen_01_3_2 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,2),1),3): " << freq_gen_02_1_3 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,2),3),1): " << freq_gen_02_3_1 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,3),1),2): " << freq_gen_03_1_2 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((0,3),2),1): " << freq_gen_03_2_1 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((1,2),0),3): " << freq_gen_12_0_3 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((1,2),3),0): " << freq_gen_12_3_0 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((1,3),0),2): " << freq_gen_13_0_2 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((1,3),2),0): " << freq_gen_13_2_0 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((2,3),0),1): " << freq_gen_23_0_1 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of gen (((2,3),1),0): " << freq_gen_23_1_0 << " (expected " << exp_freq << ")\n";
        std::cout << "Freq of 2 heights: " << freq_2_heights << " (expected " << 13 * exp_freq << ")\n";
        std::cout << "Freq of 3 heights: " << freq_3_heights << " (expected " << 15 * exp_freq << ")\n";

        double eps = 0.001;

        REQUIRE(freq_2_heights == Approx(13 * exp_freq).epsilon(eps));
        REQUIRE(freq_3_heights == Approx(15 * exp_freq).epsilon(eps));

        REQUIRE(freq_0123 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_01_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_02_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_03_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_12_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_13_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_23_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_0_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_1_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_2_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_3_ == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_01_23 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_02_13 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_03_12 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_01_23 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_02_13 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_03_12 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_01_2_3 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_01_3_2 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_02_1_3 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_02_3_1 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_03_1_2 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_03_2_1 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_12_0_3 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_12_3_0 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_13_0_2 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_13_2_0 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_23_0_1 == Approx(exp_freq).epsilon(eps));
        REQUIRE(freq_gen_23_1_0 == Approx(exp_freq).epsilon(eps));
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler::merge with ladderized tree with 4 leaves",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing merge with ladderized tree with 4") {
        RandomNumberGenerator rng = RandomNumberGenerator(22);

        // MOVE FROM: (((A:0.1,B:0.1):0.1,C:0.2):0.2,D:0.4)
        // MOVE TO:   ((A:0.2,B:0.2,C:0.2):0.2,D:0.4)
        //            OR
        //            ((A:0.1,B:0.1):0.3,C:0.4,D:0.4)
        //
        // HR of move to ((A:0.2,B:0.2,C:0.2):0.2,D:0.4)
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to merge) * pr(choose height)
        //                  =  1 * 1/2 = 1/2
        // pr(reverse move) = pr(choose to split) * pr(choose height) * pr(choosing nodes to move down) * pr(new height)
        //                  = 1/2 * 1 * 1/3 * 1/0.2
        //                  = 1/(6*0.2)
        // HR = 1/(6*0.2) / 1/2 = 2/(6*0.2) = 1 / (3*0.2)
        //
        // HR of move to ((A:0.1,B:0.1):0.3,C:0.4,D:0.4)
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to merge) * pr(choose height)
        //                  =  1 * 1/2 = 1/2
        // pr(reverse move) = pr(choose to split) * pr(choose height) * pr(choosing nodes to move down) * pr(new height)
        //                  = 1/2 * 1 * 1/3 * 1/0.3
        //                  = 1/(6*0.3)
        // HR = 1/(6*0.3) / 1/2 = 2/(6*0.3) = 1 / (3*0.3)

        unsigned int nsamples = 10000;
        unsigned int root_poly_count = 0;
        unsigned int internal_poly_count = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.4;
            std::shared_ptr<Node> root = std::make_shared<Node>("root", root_ht);
            std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.2);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.1);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);

            internal0->add_child(leaf0);
            internal0->add_child(leaf1);
            internal1->add_child(leaf2);
            internal1->add_child(internal0);
            root->add_child(leaf3);
            root->add_child(internal1);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            (1.0 / 0.4) *
                            (1.0 / 0.2)
                            )).epsilon(1e-8));

            double ln_hastings = op.propose(rng,
                    &tree);
            REQUIRE(tree.get_number_of_node_heights() == 2);
            REQUIRE(tree.get_number_of_splittable_heights() == 1);
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / root_ht)).epsilon(1e-8));
            // Pr(forward move) = pr(merging) * pr(picking node) = 1 * 1/2 = 1/2
            // Pr(reverse move) = pr(splitting) * pr(picking node) * pr(paritioning nodes) * pr(new height)
            // = 1/2 * 1 * 1/3 * 1/height diff = 1/6d
            // HR = 1/6d / 1/2 = 2 / 6d = 1 / 3d
            if (tree.get_root_ptr()->get_number_of_children() == 2) {
                double exp_ln_hastings = std::log(1.0 / (3.0 * 0.2));
                REQUIRE(ln_hastings == Approx(exp_ln_hastings).epsilon(1e-8));
                REQUIRE(tree.get_root_ptr()->is_child("leaf3"));
                REQUIRE(! tree.get_root_ptr()->is_child("leaf2"));
                ++internal_poly_count;
            }
            else if (tree.get_root_ptr()->get_number_of_children() == 3) {
                double exp_ln_hastings = std::log(1.0 / (3.0 * 0.3));
                REQUIRE(ln_hastings == Approx(exp_ln_hastings).epsilon(1e-8));
                REQUIRE(tree.get_root_ptr()->is_child("leaf3"));
                REQUIRE(tree.get_root_ptr()->is_child("leaf2"));
                ++root_poly_count;
            }
            else {
                REQUIRE(0 == 1);
            }
        }
        REQUIRE(root_poly_count + internal_poly_count == nsamples);
        REQUIRE(root_poly_count / (double)nsamples == Approx(0.5).epsilon(0.01));
        REQUIRE(internal_poly_count / (double)nsamples == Approx(0.5).epsilon(0.01));
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler::split to ladderized tree with 4 leaves",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing split to ladderized tree with 4") {
        RandomNumberGenerator rng = RandomNumberGenerator(23);

        // MOVE FROM: ((A:0.2,B:0.2,C:0.2):0.1,D:0.3)
        // MOVE TO:   (((A:<0.2,B:<0.2):?,C:0.2):0.1,D:0.3)
        //            OR
        //            (((A:<0.2,C:<0.2):?,B:0.2):0.1,D:0.3)
        //            OR
        //            (((B:<0.2,C:<0.2):?,A:0.2):0.1,D:0.3)
        //
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to split) * pr(choose height) * pr(choosing nodes to move down) * pr(new height)
        //                  = 1/2 * 1 * 1/3 * 1/(0.2)
        //                  = 1/(6*0.2)
        // pr(reverse move) = pr(choose to merge) * pr(choose height)
        //                  =  1 * 1/2 = 1/2
        // HR = 1/2 / 1/(6*0.2) = (6*0.2) / 2 = (3*0.2)

        unsigned int nsamples = 10000;
        unsigned int count_01 = 0;
        unsigned int count_02 = 0;
        unsigned int count_12 = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.3;
            std::shared_ptr<Node> root = std::make_shared<Node>("root", root_ht);
            std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.2);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);

            internal1->add_child(leaf0);
            internal1->add_child(leaf1);
            internal1->add_child(leaf2);
            root->add_child(leaf3);
            root->add_child(internal1);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / root_ht)).epsilon(1e-8));

            double ln_hastings = op.propose_split(rng,
                    &tree);
            double exp_ln_hastings = std::log(3.0 * 0.2);
            REQUIRE(ln_hastings == Approx(exp_ln_hastings).epsilon(1e-8));
            REQUIRE(tree.get_number_of_node_heights() == 3);
            REQUIRE(tree.get_number_of_splittable_heights() == 0);
            REQUIRE(tree.get_root_ptr()->is_child("leaf3"));
            REQUIRE(tree.get_root_ptr()->get_number_of_children() == 2);
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            (1.0 / root_ht) *
                            (1.0 / 0.2)
                            )).epsilon(1e-8));
            if (tree.get_root_ptr()->get_child(0)->is_child("leaf2") ||
                tree.get_root_ptr()->get_child(1)->is_child("leaf2")) {
                ++count_01;
            }
            else if (tree.get_root_ptr()->get_child(0)->is_child("leaf1") ||
                     tree.get_root_ptr()->get_child(1)->is_child("leaf1")) {
                ++count_02;
            }
            else if (tree.get_root_ptr()->get_child(0)->is_child("leaf0") ||
                     tree.get_root_ptr()->get_child(1)->is_child("leaf0")) {
                ++count_12;
            }
            else {
                REQUIRE(0 == 1);
            }
        }
        double eps = 0.01;
        REQUIRE(count_01 + count_02 + count_12 == nsamples);
        REQUIRE(count_01 / (double)nsamples == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(count_02 / (double)nsamples == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(count_12 / (double)nsamples == Approx(1.0/3.0).epsilon(eps));
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler::split from root poly to general tree with 4 leaves",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing split from root poly to general tree with 4 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(24);

        // MOVE FROM: ((A:0.2,B:0.2):0.1,C:0.3,D:0.3)
        // MOVE TO:   ((A:0.2,B:0.2):0.1,(C:0.2<x<0.3, D:0.2<x<0.3):<0.1)
        //            OR
        //            (((A:0.2,B:0.2):<0.1,C:<0.3),D:0.3)
        //            OR
        //            (((A:0.2,B:0.2):<0.1,D:<0.3),C:0.3)
        //
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to split) * pr(choose height) * pr(choosing nodes to move down) * pr(new height)
        //                  = 1/2 * 1 * 1/3 * 1/(0.3-0.2)
        //                  = 1/(6*0.1)
        // pr(reverse move) = pr(choose to merge) * pr(choose height)
        //                  =  1 * 1/2 = 1/2
        // HR = 1/2 / 1/(6*0.1) = (6*0.1) / 2 = (3*0.1)

        unsigned int nsamples = 10000;
        unsigned int ladder_count_2 = 0;
        unsigned int ladder_count_3 = 0;
        unsigned int balanced_count = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.3;
            std::shared_ptr<Node> root = std::make_shared<Node>("root", root_ht);
            std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.2);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);

            internal1->add_child(leaf0);
            internal1->add_child(leaf1);
            root->add_child(leaf2);
            root->add_child(leaf3);
            root->add_child(internal1);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / root_ht)).epsilon(1e-8));

            double ln_hastings = op.propose_split(rng,
                    &tree);
            double exp_ln_hastings = std::log(3.0 * 0.1);
            REQUIRE(ln_hastings == Approx(exp_ln_hastings).epsilon(1e-8));
            REQUIRE(tree.get_number_of_node_heights() == 3);
            REQUIRE(tree.get_number_of_splittable_heights() == 0);
            REQUIRE(tree.get_root_ptr()->get_number_of_children() == 2);
            REQUIRE(tree.get_height(0) == 0.2);
            REQUIRE(tree.get_height(2) == 0.3);
            REQUIRE(
                    ((tree.get_height(1) < 0.3) && (tree.get_height(1) > 0.2))
                    );
            tree.compute_log_likelihood_and_prior(true);
            if (tree.get_root_ptr()->get_child(0)->is_leaf() ||
                tree.get_root_ptr()->get_child(1)->is_leaf()) {
                if (tree.get_root_ptr()->is_child("leaf2")) {
                    ++ladder_count_2;
                }
                else if (tree.get_root_ptr()->is_child("leaf3")) {
                    ++ladder_count_3;
                }
                else {
                    REQUIRE(0 == 1);
                }
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                                (1.0 / root_ht) *
                                (1.0 / tree.get_height(1))
                                )).epsilon(1e-8));
            }
            else if ((! tree.get_root_ptr()->get_child(0)->is_leaf()) &&
                (! tree.get_root_ptr()->get_child(1)->is_leaf())) {
                ++balanced_count;
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                                (1.0 / root_ht) *
                                (1.0 / root_ht)
                                )).epsilon(1e-8));
            }
            else {
                REQUIRE(0 == 1);
            }
        }
        std::cout << "Freq of (((0,1),2),3): " << ladder_count_3 / (double)nsamples << "\n";
        std::cout << "Freq of (((0,1),3),2): " << ladder_count_2 / (double)nsamples << "\n";
        std::cout << "Freq of ((0,1),(2,3)): " << balanced_count / (double)nsamples << "\n";
        double eps = 0.01;
        REQUIRE(ladder_count_2 + ladder_count_3 + balanced_count == nsamples);
        REQUIRE(ladder_count_2 / (double)nsamples == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(ladder_count_3 / (double)nsamples == Approx(1.0/3.0).epsilon(eps));
        REQUIRE(balanced_count / (double)nsamples == Approx(1.0/3.0).epsilon(eps));
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler::split from comb with 4 leaves",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing split from comb with 4 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(25);

        // MOVE FROM: (A:0.3,B:0.3,C:0.3,D:0.3)
        // MOVE TO:   13 possible trees: 
        //            (A:0.3,B:0.3,(C:<0.3,D:<0.3):<0.3)
        //            (A:0.3,C:0.3,(B:<0.3,D:<0.3):<0.3)
        //            (C:0.3,B:0.3,(A:<0.3,D:<0.3):<0.3)
        //            (A:0.3,D:0.3,(B:<0.3,C:<0.3):<0.3)
        //            (B:0.3,D:0.3,(A:<0.3,C:<0.3):<0.3)
        //            (C:0.3,D:0.3,(A:<0.3,B:<0.3):<0.3)
        //            (A:0.3,(B:<0.3,C:<0.3,D:<0.3):<0.3)
        //            (B:0.3,(A:<0.3,C:<0.3,D:<0.3):<0.3)
        //            (C:0.3,(B:<0.3,A:<0.3,D:<0.3):<0.3)
        //            (D:0.3,(B:<0.3,C:<0.3,A:<0.3):<0.3)
        //            ((A:<0.3,B:<0.3):<0.3,(C:<0.3,D:<0.3):<0.3)
        //            ((A:<0.3,C:<0.3):<0.3,(B:<0.3,D:<0.3):<0.3)
        //            ((A:<0.3,D:<0.3):<0.3,(C:<0.3,B:<0.3):<0.3)
        //
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to split) * pr(choose height) * pr(choosing nodes to move down) * pr(new height)
        //                  = 1 * 1 * 1/13 * 1/0.3
        // pr(reverse move) = pr(choose to merge) * pr(choose height)
        //                  =  1/2 * 1 = 1/2
        // HR = 1/2 / 1 / (13*0.3) = (13*0.3) / 2

        unsigned int nsamples = 100000;
        unsigned int balanced_01_23_count = 0;
        unsigned int balanced_02_13_count = 0;
        unsigned int balanced_03_12_count = 0;
        unsigned int root_poly_01_count = 0;
        unsigned int root_poly_02_count = 0;
        unsigned int root_poly_03_count = 0;
        unsigned int root_poly_12_count = 0;
        unsigned int root_poly_13_count = 0;
        unsigned int root_poly_23_count = 0;
        unsigned int internal_poly_0_count = 0;
        unsigned int internal_poly_1_count = 0;
        unsigned int internal_poly_2_count = 0;
        unsigned int internal_poly_3_count = 0;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.3;
            std::shared_ptr<Node> root = std::make_shared<Node>("root", root_ht);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);

            root->add_child(leaf0);
            root->add_child(leaf1);
            root->add_child(leaf2);
            root->add_child(leaf3);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / pow(root_ht, 0))).epsilon(1e-8));

            // Calling propose, because split should be only possibility
            double ln_hastings = op.propose(rng,
                    &tree);
            double exp_ln_hastings = std::log((13.0 * 0.3) / 2.0);
            REQUIRE(ln_hastings == Approx(exp_ln_hastings).epsilon(1e-8));
            REQUIRE(tree.get_number_of_node_heights() == 2);
            REQUIRE(tree.get_number_of_splittable_heights() == 1);
            REQUIRE(tree.get_height(1) == root_ht);
            REQUIRE(tree.get_height(0) < root_ht);
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / pow(root_ht, 1))).epsilon(1e-8));
            if (tree.get_root_ptr()->get_number_of_children() == 3) {
                if (tree.get_root_ptr()->is_child("leaf2") &&
                    tree.get_root_ptr()->is_child("leaf3")) {
                    ++root_poly_01_count;
                }
                else if (tree.get_root_ptr()->is_child("leaf1") &&
                    tree.get_root_ptr()->is_child("leaf3")) {
                    ++root_poly_02_count;
                }
                else if (tree.get_root_ptr()->is_child("leaf1") &&
                    tree.get_root_ptr()->is_child("leaf2")) {
                    ++root_poly_03_count;
                }
                else if (tree.get_root_ptr()->is_child("leaf0") &&
                    tree.get_root_ptr()->is_child("leaf3")) {
                    ++root_poly_12_count;
                }
                else if (tree.get_root_ptr()->is_child("leaf0") &&
                    tree.get_root_ptr()->is_child("leaf2")) {
                    ++root_poly_13_count;
                }
                else if (tree.get_root_ptr()->is_child("leaf0") &&
                    tree.get_root_ptr()->is_child("leaf1")) {
                    ++root_poly_23_count;
                }
                else {
                    REQUIRE(0 == 1);
                }
            }
            else if (tree.get_root_ptr()->get_number_of_children() == 2) {
                if (tree.get_root_ptr()->get_child(0)->is_leaf() ||
                    tree.get_root_ptr()->get_child(1)->is_leaf()) {
                    if (tree.get_root_ptr()->is_child("leaf0")) {
                        ++internal_poly_0_count;
                    }
                    else if (tree.get_root_ptr()->is_child("leaf1")) {
                        ++internal_poly_1_count;
                    }
                    else if (tree.get_root_ptr()->is_child("leaf2")) {
                        ++internal_poly_2_count;
                    }
                    else if (tree.get_root_ptr()->is_child("leaf3")) {
                        ++internal_poly_3_count;
                    }
                    else {
                        REQUIRE(0 == 1);
                    }
                }
                else if ((! tree.get_root_ptr()->get_child(0)->is_leaf()) &&
                    (! tree.get_root_ptr()->get_child(1)->is_leaf())) {
                    if (
                            (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                             tree.get_root_ptr()->get_child(0)->is_child("leaf1")) ||
                            (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                             tree.get_root_ptr()->get_child(1)->is_child("leaf1"))
                        ) {
                        ++balanced_01_23_count;
                    }
                    else if (
                            (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                             tree.get_root_ptr()->get_child(0)->is_child("leaf2")) ||
                            (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                             tree.get_root_ptr()->get_child(1)->is_child("leaf2"))
                        ) {
                        ++balanced_02_13_count;
                    }
                    else if (
                            (tree.get_root_ptr()->get_child(0)->is_child("leaf0") &&
                             tree.get_root_ptr()->get_child(0)->is_child("leaf3")) ||
                            (tree.get_root_ptr()->get_child(1)->is_child("leaf0") &&
                             tree.get_root_ptr()->get_child(1)->is_child("leaf3"))
                        ) {
                        ++balanced_03_12_count;
                    }
                    else {
                        REQUIRE(0 == 1);
                    }
                }
                else {
                    REQUIRE(0 == 1);
                }
            }
            else {
                REQUIRE(0 == 1);
            }
        }

        double eps = 0.003;
        REQUIRE(balanced_01_23_count +
                balanced_02_13_count +
                balanced_03_12_count +
                internal_poly_0_count +
                internal_poly_1_count +
                internal_poly_2_count +
                internal_poly_3_count +
                root_poly_01_count +
                root_poly_02_count +
                root_poly_03_count +
                root_poly_12_count +
                root_poly_13_count +
                root_poly_23_count == nsamples);
        // There are 3 balanced-shared node trees
        //           6 trees with trichotomy at root
        //           4 trees with trichotomy that is a child of the root
        std::cout << "Freq of ((0,1),(2,3)): " << balanced_01_23_count / (double)nsamples << "\n";
        std::cout << "Freq of ((0,2),(1,3)): " << balanced_02_13_count / (double)nsamples << "\n";
        std::cout << "Freq of ((0,3),(1,2)): " << balanced_03_12_count / (double)nsamples << "\n";
        std::cout << "Freq of ((0,1,2),3): "   << internal_poly_3_count / (double)nsamples << "\n";
        std::cout << "Freq of ((0,1,3),2): "   << internal_poly_2_count / (double)nsamples << "\n";
        std::cout << "Freq of ((0,3,2),1): "   << internal_poly_1_count / (double)nsamples << "\n";
        std::cout << "Freq of ((3,1,2),0): "   << internal_poly_0_count / (double)nsamples << "\n";
        std::cout << "Freq of ((0,1),2,3): "   << root_poly_01_count / (double)nsamples << "\n";
        std::cout << "Freq of ((0,2),1,3): "   << root_poly_02_count / (double)nsamples << "\n";
        std::cout << "Freq of ((0,3),2,1): "   << root_poly_03_count / (double)nsamples << "\n";
        std::cout << "Freq of ((1,2),0,3): "   << root_poly_12_count / (double)nsamples << "\n";
        std::cout << "Freq of ((1,3),0,2): "   << root_poly_13_count / (double)nsamples << "\n";
        std::cout << "Freq of ((2,3),0,1): "   << root_poly_23_count / (double)nsamples << "\n";
        REQUIRE(balanced_01_23_count / (double)nsamples == Approx(1.0/13.0).epsilon(eps));
        REQUIRE(balanced_02_13_count / (double)nsamples == Approx(1.0/13.0).epsilon(eps));
        REQUIRE(balanced_03_12_count / (double)nsamples == Approx(1.0/13.0).epsilon(eps));
        REQUIRE(internal_poly_0_count / (double)nsamples == Approx(1.0/13.0).epsilon(eps));
        REQUIRE(internal_poly_1_count / (double)nsamples == Approx(1.0/13.0).epsilon(eps));
        REQUIRE(internal_poly_2_count / (double)nsamples == Approx(1.0/13.0).epsilon(eps));
        REQUIRE(internal_poly_3_count / (double)nsamples == Approx(1.0/13.0).epsilon(eps));
        REQUIRE(root_poly_01_count / (double)nsamples == Approx(1.0/13.0).epsilon(eps));
        REQUIRE(root_poly_02_count / (double)nsamples == Approx(1.0/13.0).epsilon(eps));
        REQUIRE(root_poly_03_count / (double)nsamples == Approx(1.0/13.0).epsilon(eps));
        REQUIRE(root_poly_12_count / (double)nsamples == Approx(1.0/13.0).epsilon(eps));
        REQUIRE(root_poly_13_count / (double)nsamples == Approx(1.0/13.0).epsilon(eps));
        REQUIRE(root_poly_23_count / (double)nsamples == Approx(1.0/13.0).epsilon(eps));
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler::merge from balanced to comb with 4 leaves",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing merge from balanced to comb with 4 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(26);

        // MOVE FROM: ((A:0.2,B:0.2):0.1,(C:0.2, D:0.2):0.1)
        // MOVE TO:   (A:0.3,B:0.3,C:0.3,D:0.3)
        //
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to merge) * pr(choose height)
        //                  =  1/2 * 1 = 1/2
        // pr(reverse move) = pr(choose to split) * pr(choose height) * pr(choosing nodes to move down) * pr(new height)
        //                  = 1 * 1 * 1/13 * 1/0.3
        //                  there are 13 ways to split up the root polytomy
        // HR = 1/(13*0.3) / 1/2 = 2/(13*0.3)

        unsigned int nsamples = 20;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.3;
            std::shared_ptr<Node> root = std::make_shared<Node>("root", root_ht);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.2);
            std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.2);
            internal0->set_height_parameter(internal1->get_height_parameter());
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);

            internal0->add_child(leaf0);
            internal0->add_child(leaf1);
            internal1->add_child(leaf2);
            internal1->add_child(leaf3);
            root->add_child(internal0);
            root->add_child(internal1);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_number_of_node_heights() == 2);
            REQUIRE(tree.get_number_of_splittable_heights() == 1);
            REQUIRE(tree.get_height(0) == 0.2);
            REQUIRE(tree.get_height(1) == 0.3);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / pow(root_ht, 1))).epsilon(1e-8));

            double ln_hastings = op.propose_merge(rng,
                    &tree,
                    false);
            double exp_ln_hastings = std::log(2.0 / (13.0 * 0.3));
            REQUIRE(ln_hastings == Approx(exp_ln_hastings).epsilon(1e-8));
            REQUIRE(tree.get_number_of_node_heights() == 1);
            REQUIRE(tree.get_number_of_splittable_heights() == 1);
            REQUIRE(tree.get_root_ptr()->get_number_of_children() == 4);
            REQUIRE(tree.get_height(0) == 0.3);
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / pow(root_ht, 0))).epsilon(1e-8));
        }
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler::merge from internal poly to comb with 4 leaves",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing merge from internal poly to comb with 4 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(27);

        // MOVE FROM: ((A:0.2,B:0.2,C:0.2):0.1,D:0.3)
        // MOVE TO:   (A:0.3,B:0.3,C:0.3,D:0.3)
        //
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to merge) * pr(choose height)
        //                  =  1/3 * 1 = 1/2
        // pr(reverse move) = pr(choose to split) * pr(choose height) * pr(choosing nodes to move down) * pr(new height)
        //                  = 1 * 1 * 1/13 * 1/0.3
        //                  there are 13 ways to split up the root polytomy
        // HR = 1/(13*0.3) / 1/2 = 2/(13*0.3)

        unsigned int nsamples = 20;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.3;
            std::shared_ptr<Node> root = std::make_shared<Node>("root", root_ht);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.2);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);

            internal0->add_child(leaf0);
            internal0->add_child(leaf1);
            internal0->add_child(leaf2);
            root->add_child(leaf3);
            root->add_child(internal0);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_number_of_node_heights() == 2);
            REQUIRE(tree.get_number_of_splittable_heights() == 1);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / pow(root_ht, 1))).epsilon(1e-8));

            double ln_hastings = op.propose_merge(rng,
                    &tree,
                    false);
            double exp_ln_hastings = std::log(2.0 / (13.0 * 0.3));
            REQUIRE(ln_hastings == Approx(exp_ln_hastings).epsilon(1e-8));
            REQUIRE(tree.get_number_of_node_heights() == 1);
            REQUIRE(tree.get_number_of_splittable_heights() == 1);
            REQUIRE(tree.get_root_ptr()->get_number_of_children() == 4);
            REQUIRE(tree.get_height(0) == 0.3);
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / pow(root_ht, 0))).epsilon(1e-8));
        }
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler::merge from root poly to comb with 4 leaves",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing merge from root poly to comb with 4 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(28);

        // MOVE FROM: ((A:0.2,B:0.2):0.1,C:0.3,D:0.3)
        // MOVE TO:   (A:0.3,B:0.3,C:0.3,D:0.3)
        //
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to merge) * pr(choose height)
        //                  =  1/2 * 1 = 1/2
        // pr(reverse move) = pr(choose to split) * pr(choose height) * pr(choosing nodes to move down) * pr(new height)
        //                  = 1 * 1 * 1/13 * 1/0.3
        //                  There are 13 possible ways to break up the root polytomy
        //                  = 1/(13*0.3)
        // HR = 1/(13*0.3) / 1/2 = 2/(13*0.3)

        unsigned int nsamples = 20;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.3;
            std::shared_ptr<Node> root = std::make_shared<Node>("root", root_ht);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.2);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);

            internal0->add_child(leaf0);
            internal0->add_child(leaf1);
            root->add_child(leaf2);
            root->add_child(leaf3);
            root->add_child(internal0);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_number_of_node_heights() == 2);
            REQUIRE(tree.get_number_of_splittable_heights() == 1);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / pow(root_ht, 1))).epsilon(1e-8));

            double ln_hastings = op.propose_merge(rng,
                    &tree,
                    false);
            double exp_ln_hastings = std::log(2.0 / (13.0 * 0.3));
            REQUIRE(ln_hastings == Approx(exp_ln_hastings).epsilon(1e-8));
            REQUIRE(tree.get_number_of_node_heights() == 1);
            REQUIRE(tree.get_number_of_splittable_heights() == 1);
            REQUIRE(tree.get_root_ptr()->get_number_of_children() == 4);
            REQUIRE(tree.get_height(0) == 0.3);
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / pow(root_ht, 0))).epsilon(1e-8));
        }
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler::split from shared to balanced with 4 leaves",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing split from shared to balanced with 4 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(29);

        // MOVE FROM: ((A:0.2,B:0.2):0.1,(C:0.2, D:0.2):0.1)
        // MOVE TO:   ((A:<0.2,B:<0.2):>0.1,(C:0.2, D:0.2):0.1)
        //            OR
        //            ((A:0.2,B:0.2):0.1,(C:<0.2, D:<0.2):>0.1)
        //
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to split) * pr(choose height) * pr(choosing nodes to move down) * pr(new height)
        //                  = 1/2 * 1 * 1/2 * 1/0.2
        // pr(reverse move) = pr(choose to merge) * pr(choose height)
        //                  =  1 * 1/2 = 1/2
        // HR = 1/2 / 1/(4*0.2) = 2*0.2

        unsigned int split0_count = 0;
        unsigned int split1_count = 0;
        unsigned int nsamples = 10000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.3;
            std::shared_ptr<Node> root = std::make_shared<Node>("root", root_ht);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.2);
            std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.2);
            internal0->set_height_parameter(internal1->get_height_parameter());
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);

            internal0->add_child(leaf0);
            internal0->add_child(leaf1);
            internal1->add_child(leaf2);
            internal1->add_child(leaf3);
            root->add_child(internal0);
            root->add_child(internal1);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_number_of_node_heights() == 2);
            REQUIRE(tree.get_number_of_splittable_heights() == 1);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / pow(root_ht, 1))).epsilon(1e-8));

            double ln_hastings = op.propose_split(rng,
                    &tree);
            double exp_ln_hastings = std::log(2.0 * 0.2);
            REQUIRE(ln_hastings == Approx(exp_ln_hastings).epsilon(1e-8));
            REQUIRE(tree.get_number_of_node_heights() == 3);
            REQUIRE(tree.get_number_of_splittable_heights() == 0);
            REQUIRE(tree.get_root_ptr()->get_number_of_children() == 2);
            REQUIRE(tree.get_height(2) == 0.3);
            REQUIRE(tree.get_height(1) == 0.2);
            REQUIRE(tree.get_height(0) < 0.2);
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / pow(root_ht, 2))).epsilon(1e-8));
            if (tree.get_root_ptr()->get_child(0)->get_height() < 0.2) {
                ++split0_count;
            }
            else if (tree.get_root_ptr()->get_child(1)->get_height() < 0.2) {
                ++split1_count;
            }
            else {
                REQUIRE(0 == 1);
            }
        }
        REQUIRE(split0_count + split1_count == nsamples);
        double eps = 0.01;
        REQUIRE(split0_count / (double)nsamples == Approx(0.5).epsilon(eps));
        REQUIRE(split1_count / (double)nsamples == Approx(0.5).epsilon(eps));
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler::merge from balanced general with 4 leaves",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing merge from balanced general with 4 leaves") {
        RandomNumberGenerator rng = RandomNumberGenerator(30);

        // MOVE FROM: ((A:0.2,B:0.2):0.2,(C:0.1, D:0.1):0.3)
        // MOVE TO:   ((A:0.2,B:0.2):0.2,(C:0.2, D:0.2):0.2)
        //            OR
        //            (A:0.4,B:0.4,(C:0.1, D:0.1):0.3)
        //
        // HR of move to ((A:0.2,B:0.2):0.2,(C:0.2, D:0.2):0.2)
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to merge) * pr(choose height)
        //                  =  1 * 1/2 = 1/2
        // pr(reverse move) = pr(choose to split) * pr(choose height) * pr(choosing nodes to move down) * pr(new height)
        //                  = 1/2 * 1 * 1/2 * 1/0.2
        //                  = 1/(4*0.2)
        // HR = 1/(4*0.2) / 1/2 = 2/(4*0.2) = 1 / (2*0.2)
        //
        // HR of move to ((A:0.4,B:0.4,(C:0.1, D:0.1):0.3)
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to merge) * pr(choose height)
        //                  =  1 * 1/2 = 1/2
        // pr(reverse move) = pr(choose to split) * pr(choose height) * pr(choosing nodes to move down) * pr(new height)
        //                  = 1/2 * 1 * 1/3 * 1/0.3
        //                  = 1/(6*0.3)
        // HR = 1/(6*0.3) / 1/2 = 2/(6*0.3) = 1 / (3*0.3)
        unsigned int balanced_count = 0;
        unsigned int root_poly_count = 0;
        unsigned int nsamples = 10000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.4;
            std::shared_ptr<Node> root = std::make_shared<Node>("root", root_ht);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>("internal0", 0.2);
            std::shared_ptr<Node> internal1 = std::make_shared<Node>("internal1", 0.1);
            std::shared_ptr<Node> leaf0 = std::make_shared<Node>("leaf0", 0.0);
            std::shared_ptr<Node> leaf1 = std::make_shared<Node>("leaf1", 0.0);
            std::shared_ptr<Node> leaf2 = std::make_shared<Node>("leaf2", 0.0);
            std::shared_ptr<Node> leaf3 = std::make_shared<Node>("leaf3", 0.0);

            internal0->add_child(leaf0);
            internal0->add_child(leaf1);
            internal1->add_child(leaf2);
            internal1->add_child(leaf3);
            root->add_child(internal0);
            root->add_child(internal1);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_number_of_node_heights() == 3);
            REQUIRE(tree.get_number_of_splittable_heights() == 0);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / pow(root_ht, 2))).epsilon(1e-8));

            double ln_hastings = op.propose(rng,
                    &tree);
            if (tree.get_root_ptr()->get_number_of_children() == 3) {
                double exp_ln_hastings = std::log(1.0 / (3.0 * 0.3));
                REQUIRE(ln_hastings == Approx(exp_ln_hastings).epsilon(1e-8));
                REQUIRE(tree.get_number_of_node_heights() == 2);
                REQUIRE(tree.get_number_of_splittable_heights() == 1);
                REQUIRE(tree.get_height(0) == 0.1);
                REQUIRE(tree.get_height(1) == 0.4);
                ++root_poly_count;
            }
            else if (tree.get_root_ptr()->get_number_of_children() == 2) {
                // Pr(forward) = Pr(merge) * Pr(pick node) = 1 * 1/2 = 1/2
                // Pr(reverse) = Pr(split) * Pr(pick height) * Pr(parition nodes) * Pr(new height)
                //             = 1/2 * 1 * 1/2 * 1/d = 1/4d
                // HR = 1/4d / 1/2 = 2/4d = 1/2d
                double exp_ln_hastings = std::log(1.0 / (2.0 * 0.2));
                REQUIRE(ln_hastings == Approx(exp_ln_hastings).epsilon(1e-8));
                REQUIRE(tree.get_number_of_node_heights() == 2);
                REQUIRE(tree.get_number_of_splittable_heights() == 1);
                REQUIRE(tree.get_height(0) == 0.2);
                REQUIRE(tree.get_height(1) == 0.4);
                ++balanced_count;
            }
            else {
                REQUIRE(0 == 1);
            }
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / pow(root_ht, 1))).epsilon(1e-8));
        }
        REQUIRE(root_poly_count + balanced_count == nsamples);
        double eps = 0.01;
        REQUIRE(root_poly_count / (double)nsamples == Approx(0.5).epsilon(eps));
        REQUIRE(balanced_count / (double)nsamples == Approx(0.5).epsilon(eps));
    }
}



// Expectations for tree with 5 leaves:
// ----------------------------------------------------------------------------
// # of unlabeled                            # of trees (labelings)      
//  topologies
// ----------------------------------------------------------------------------
// (a,b,c,d,e)                                                     = 1
// (a,(b,c,d,e))           = 5 choose 4                            = 5
// (a,b,(c,d,e))           = 5 choose 3                            = 10 
// (a,b,c,(d,e))           = 5 choose 2                            = 10
// !((a,b,c),(d,e))        = 5 choose 3                            = 10
// (a,(b,(c,d,e)))         = (5 choose 3) * 2!                     = 20 
// !((a,b),(c,d),e)        = ((5 choose 2) * (3 choose 2)) / 2     = 15
// !(((a,b),(c,d)),e       = ((5 choose 2) * (3 choose 2)) / 2     = 15
// !(((a,b),c),(d,e))      = (5 choose 2) * (3 choose 2)           = 30
// (a,(b,(c,(d,e))))       = (5 choose 2) * 3!                     = 60
// (a,b,(c,(d,e)))         = (5 choose 2) * (3 choose 2)           = 30
// (a,(b,c,(d,e)))         = (5 choose 2) * (3 choose 2)           = 30
// ----------------------------------------------------------------------------
// TOTAL                                                           = 236
//
// 236 matches Felsenstein 1978, but we need to account for shared node
// heights. The topologies above prefixed with '!' are topologies that have
// potentially shared node heights. For each shared node configuration of these
// topologies, we have to add that many additional trees, which we do below.
//
// ----------------------------------------------------------------------------
// # of unlabeled                            # of trees (labelings)      
//  topologies
// ----------------------------------------------------------------------------
// ((a,b,c)*,(d,e)*)       = 5 choose 3                            = 10
// ((a,b)*,(c,d)*,e)       = ((5 choose 2) * (3 choose 2)) / 2     = 15
// (((a,b)*,(c,d)*),e      = ((5 choose 2) * (3 choose 2)) / 2     = 15
// (((a,b),c)*,(d,e)*)     = (5 choose 2) * (3 choose 2)           = 30
// (((a,b)*,c),(d,e)*)     = (5 choose 2) * (3 choose 2)           = 30
// ----------------------------------------------------------------------------
// GRAND TOTAL # OF TREE MODELS                                    = 336
//
// The asterisks in the topologies above indicated shared node heights.
TEST_CASE("Testing SplitLumpNodesRevJumpSampler with 5 leaves and fixed root",
        "[xSplitLumpNodesRevJumpSampler]") {

    SECTION("Testing 5 leaves with fixed root") {
        RandomNumberGenerator rng = RandomNumberGenerator(2193647912);

        double root_ht = 0.5;
        std::shared_ptr<Node> root = std::make_shared<Node>(5, "root", root_ht);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);

        root->add_child(leaf0);
        root->add_child(leaf1);
        root->add_child(leaf2);
        root->add_child(leaf3);
        root->add_child(leaf4);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        SplitLumpNodesRevJumpSampler<Node> op;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int count_nheights_1 = 0;
        unsigned int count_nheights_2 = 0;
        unsigned int count_nheights_3 = 0;
        unsigned int count_nheights_4 = 0;

        unsigned int niterations = 50000000;
        unsigned int sample_freq = 50;
        unsigned int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 10000;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                if (tree.get_number_of_node_heights() == 1) {
                    ++count_nheights_1;
                }
                else if (tree.get_number_of_node_heights() == 2) {
                    ++count_nheights_2;
                }
                else if (tree.get_number_of_node_heights() == 3) {
                    ++count_nheights_3;
                }
                else if (tree.get_number_of_node_heights() == 4) {
                    ++count_nheights_4;
                }
                std::set< std::set<Split> > splits = tree.get_splits(false);
                if (split_counts.count(splits) > 0) {
                    ++split_counts[splits];
                }
                else {
                    split_counts[splits] = 1;
                }
                ++sample_count;
                if (sample_count % report_freq == 0) {
                    std::cout << "Sampled " << sample_count << " of " << nsamples << std::endl;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        REQUIRE((count_nheights_1 + count_nheights_2 + count_nheights_3 + count_nheights_4) == nsamples);

        double freq_nheights_1 = count_nheights_1 / (double)nsamples;
        double freq_nheights_2 = count_nheights_2 / (double)nsamples;
        double freq_nheights_3 = count_nheights_3 / (double)nsamples;
        double freq_nheights_4 = count_nheights_4 / (double)nsamples;

        double exp_freq = 1.0/336.0;
        double exp_count = nsamples/336.0;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.1;
        unsigned int total_trees_sampled = 0;
        std::map< std::set< std::set<Split> >, double> split_freqs;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            split_freqs[s_c.first] = s_c.second / (double)nsamples;
            std::cout << "Tree:\n";
            for (auto splitset : s_c.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            std::cout << "  nsamples: " << s_c.second << "\n";
            std::cout << "  prop error: " << prop_error << "\n";
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }

        double quantile_chi_sq_335_10 = 368.6;
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        std::cout << "Chi-square(335) 0.9 quantile: " << quantile_chi_sq_335_10 << "\n";

        std::cout << "BAD SPLITS (proportional error > " << prop_error_threshold << ")\n";
        for (auto s_e : bad_splits) {
            std::cout << "\nTree:\n";
            for (auto splitset : s_e.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            std::cout << "  prop error: " << s_e.second << "\n";
        }

        write_r_script(split_counts, "../5-leaf-general-tree-test.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        REQUIRE(split_counts.size() == 336);

        double eps = 0.001;

        for (auto s_f : split_freqs) {
            REQUIRE(s_f.second == Approx(exp_freq).epsilon(eps));
        }

        REQUIRE(chi_sq_test_statistic < quantile_chi_sq_335_10);
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler::merge with 3-2 shared tree with 5 leaves",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing merge 3-2 shared tree") {
        RandomNumberGenerator rng = RandomNumberGenerator(22);

        // MOVE FROM: ((A:0.1,B:0.1)*:0.1,(C:0.1,D:0.1,E:0.1)*:0.1)
        // MOVE TO:   (A:0.2,B:0.2,C:0.2,D:0.2,E:0.2)
        //
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to merge) * pr(choose height)
        //                  =  1/2 * 1 = 1/2
        // pr(reverse move) = pr(choose to split) * pr(choose height) * pr(partitioning branches to move down) * pr(new height)
        //                  = 1 * 1 * (1/(52-2)) * 1/0.2
        //                  = 1/(50*0.2)
        // HR = 1/(50*0.2) / 1/2 = 2/(50*0.2)
        
        std::string expected_tree_str = "(A:0.2,B:0.2,C:0.2,D:0.2,E:0.2)[&height_index=0,height=0.2];";
        BaseTree<Node> expected_tree(expected_tree_str);

        unsigned int nsamples = 20;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.2;
            std::shared_ptr<Node> root = std::make_shared<Node>(7, "root", root_ht);
            std::shared_ptr<Node> internal1 = std::make_shared<Node>(6, "internal1", 0.1);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>(5, "internal0", 0.1);
            internal0->set_height_parameter(internal1->get_height_parameter());
            std::shared_ptr<Node> leaf_a = std::make_shared<Node>(0, "A", 0.0);
            std::shared_ptr<Node> leaf_b = std::make_shared<Node>(1, "B", 0.0);
            std::shared_ptr<Node> leaf_c = std::make_shared<Node>(2, "C", 0.0);
            std::shared_ptr<Node> leaf_d = std::make_shared<Node>(3, "D", 0.0);
            std::shared_ptr<Node> leaf_e = std::make_shared<Node>(4, "E", 0.0);

            internal0->add_child(leaf_a);
            internal0->add_child(leaf_b);
            internal1->add_child(leaf_c);
            internal1->add_child(leaf_d);
            internal1->add_child(leaf_e);
            root->add_child(internal0);
            root->add_child(internal1);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / 0.2)).epsilon(1e-8));

            double ln_hastings = op.propose_merge(rng,
                    &tree,
                    false);
            REQUIRE(tree.get_number_of_node_heights() == 1);
            REQUIRE(tree.get_number_of_splittable_heights() == 1);
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0)).epsilon(1e-8));

            double exp_ln_hastings = std::log(2.0 / (50.0 * 0.2));
            REQUIRE(ln_hastings == Approx(exp_ln_hastings).epsilon(1e-8));
            REQUIRE(tree.get_splits_by_height_index() == expected_tree.get_splits_by_height_index());
        }
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler::split with 3-2 shared tree with 5 leaves",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing merge 3-2 shared tree") {
        RandomNumberGenerator rng = RandomNumberGenerator(765492347);

        // MOVE FROM: ((A:0.1,B:0.1)*:0.1,(C:0.1,D:0.1,E:0.1)*:0.1)
        // MOVE TO:   ((A:0.05,B:0.05):0.15,(C:0.1,D:0.1,E:0.1):0.1)
        //            OR
        //            ((A:0.1,B:0.1):0.1,(C:0.05,D:0.05,E:0.05):0.15)
        //            OR
        //            ((A:0.1,B:0.1)*:0.1,(C:0.1,(D:0.05,E:0.05):0.05)*:0.1)
        //            ((A:0.1,B:0.1)*:0.1,(D:0.1,(C:0.05,E:0.05):0.05)*:0.1)
        //            ((A:0.1,B:0.1)*:0.1,(E:0.1,(C:0.05,D:0.05):0.05)*:0.1)
        //            OR
        //            ((A:0.05,B:0.05)*:0.15,(C:0.1,(D:0.05,E:0.05)*:0.05):0.1)
        //            ((A:0.05,B:0.05)*:0.15,(D:0.1,(C:0.05,E:0.05)*:0.05):0.1)
        //            ((A:0.05,B:0.05)*:0.15,(E:0.1,(C:0.05,D:0.05)*:0.05):0.1)
        //
        // MOVE TO:   ((A:0.05,B:0.05):0.15,(C:0.1,D:0.1,E:0.1):0.1)
        // HR = pr(reverse move) / pr(forward move)
        // pr(reverse move) = pr(choose to merge) * pr(choose height)
        //                  =  1/2 * 1/2 = 1/4
        // pr(forward move) = pr(choose to split) * pr(choose height) * pr(choose to move down) * pr(new height)
        //                  = 1/2 * 1 * 1/3 (can choose one node, the other, or both to move down) * 1/0.1
        //                  = 1/0.6
        // HR = 1/4 / 1/0.6 = 0.6/4 = 0.3/2 = 0.15
        //
        //
        // MOVE TO:   ((A:0.1,B:0.1):0.1,(C:0.05,D:0.05,E:0.05):0.15)
        //            OR
        //            ((A:0.1,B:0.1)*:0.1,(C:0.1,(D:0.05,E:0.05):0.05)*:0.1)
        //            ((A:0.1,B:0.1)*:0.1,(D:0.1,(C:0.05,E:0.05):0.05)*:0.1)
        //            ((A:0.1,B:0.1)*:0.1,(E:0.1,(C:0.05,D:0.05):0.05)*:0.1)
        // HR = pr(reverse move) / pr(forward move)
        // pr(reverse move) = pr(choose to merge) * pr(choose height)
        //                  =  1/2 * 1/2 = 1/4
        // pr(forward move) = pr(choose to split) * pr(choose height) * pr(choose to move down) * pr(partition poly) * pr(new height)
        //                  = 1/2 * 1 * 1/3 * 1/4 * 1/0.1
        //                  = 1/(24*0.1) = 1/2.4
        // HR = 1/4 / 1/2.4 = 2.4/4 = 0.6
        //
        //
        // MOVE TO:   ((A:0.05,B:0.05)*:0.15,(C:0.1,(D:0.05,E:0.05)*:0.05):0.1)
        //            ((A:0.05,B:0.05)*:0.15,(D:0.1,(C:0.05,E:0.05)*:0.05):0.1)
        //            ((A:0.05,B:0.05)*:0.15,(E:0.1,(C:0.05,D:0.05)*:0.05):0.1)
        // HR = pr(reverse move) / pr(forward move)
        // pr(reverse move) = pr(choose to merge) * pr(choose height)
        //                  =  1/2 * 1/2 = 1/4
        // pr(forward move) = pr(choose to split) * pr(choose height) * pr(choose to move down) * pr(partition poly) * pr(new height)
        //                  = 1/2 * 1 * 1/3 * 1/3 * 1/0.1
        //                  = 1/(18*0.1) = 1/1.8
        // HR = 1/4 / 1/1.8 = 1.8/4 = 0.45
        
        BaseTree<Node> tree_cde("((A:0.05,B:0.05)[&height_index=0,height=0.05]:0.15,(C:0.1,D:0.1,E:0.1)[&height_index=1,height=0.1]:0.1)[&height_index=2,height=0.2];");
        BaseTree<Node> tree_ab("((A:0.1,B:0.1)[&height_index=1,height=0.1]:0.1,(C:0.05,D:0.05,E:0.05)[&height_index=0,height=0.05]:0.15)[&height_index=2,height=0.2];");

        BaseTree<Node> tree_abc("((A:0.1,B:0.1)[&height_index=1,height=0.1]:0.1,(C:0.1,(D:0.05,E:0.05)[&height_index=0,height=0.05]:0.05)[&height_index=1,height=0.01]:0.1)[&height_index=2,height=0.2];");
        BaseTree<Node> tree_abd("((A:0.1,B:0.1)[&height_index=1,height=0.1]:0.1,(D:0.1,(C:0.05,E:0.05)[&height_index=0,height=0.05]:0.05)[&height_index=1,height=0.01]:0.1)[&height_index=2,height=0.2];");
        BaseTree<Node> tree_abe("((A:0.1,B:0.1)[&height_index=1,height=0.1]:0.1,(E:0.1,(C:0.05,D:0.05)[&height_index=0,height=0.05]:0.05)[&height_index=1,height=0.01]:0.1)[&height_index=2,height=0.2];");

        BaseTree<Node> tree_c("((A:0.05,B:0.05)[&height_index=0,height=0.05]:0.15,(C:0.1,(D:0.05,E:0.05)[&height_index=0,height=0.05]:0.05)[&height_index=1,height=0.1]:0.1)[&height_index=2,height=0.2];");
        BaseTree<Node> tree_d("((A:0.05,B:0.05)[&height_index=0,height=0.05]:0.15,(D:0.1,(C:0.05,E:0.05)[&height_index=0,height=0.05]:0.05)[&height_index=1,height=0.1]:0.1)[&height_index=2,height=0.2];");
        BaseTree<Node> tree_e("((A:0.05,B:0.05)[&height_index=0,height=0.05]:0.15,(E:0.1,(C:0.05,D:0.05)[&height_index=0,height=0.05]:0.05)[&height_index=1,height=0.1]:0.1)[&height_index=2,height=0.2];");

        std::set<std::map< unsigned int, std::set<Split> > > choose_ab;
        choose_ab.insert(tree_cde.get_splits_by_height_index());
        double choose_ab_hr = std::log(0.15);

        std::set<std::map< unsigned int, std::set<Split> > > choose_cde;
        choose_cde.insert(tree_ab.get_splits_by_height_index());
        choose_cde.insert(tree_abc.get_splits_by_height_index());
        choose_cde.insert(tree_abd.get_splits_by_height_index());
        choose_cde.insert(tree_abe.get_splits_by_height_index());
        double choose_cde_hr = std::log(0.6);

        std::set<std::map< unsigned int, std::set<Split> > > choose_both;
        choose_both.insert(tree_c.get_splits_by_height_index());
        choose_both.insert(tree_d.get_splits_by_height_index());
        choose_both.insert(tree_e.get_splits_by_height_index());
        double choose_both_hr = std::log(0.45);

        std::map<std::map< unsigned int, std::set<Split> >, unsigned int> counts;
        for (auto splits : choose_ab) {
            counts[splits] = 0;
        }
        for (auto splits : choose_cde) {
            counts[splits] = 0;
        }
        for (auto splits : choose_both) {
            counts[splits] = 0;
        }

        std::map< unsigned int, std::set<Split> > splits;
        unsigned int nsamples = 50000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.2;
            std::shared_ptr<Node> root = std::make_shared<Node>(7, "root", root_ht);
            std::shared_ptr<Node> internal1 = std::make_shared<Node>(6, "internal1", 0.1);
            std::shared_ptr<Node> internal0 = std::make_shared<Node>(5, "internal0", 0.1);
            internal0->set_height_parameter(internal1->get_height_parameter());
            std::shared_ptr<Node> leaf_a = std::make_shared<Node>(0, "A", 0.0);
            std::shared_ptr<Node> leaf_b = std::make_shared<Node>(1, "B", 0.0);
            std::shared_ptr<Node> leaf_c = std::make_shared<Node>(2, "C", 0.0);
            std::shared_ptr<Node> leaf_d = std::make_shared<Node>(3, "D", 0.0);
            std::shared_ptr<Node> leaf_e = std::make_shared<Node>(4, "E", 0.0);

            internal0->add_child(leaf_a);
            internal0->add_child(leaf_b);
            internal1->add_child(leaf_c);
            internal1->add_child(leaf_d);
            internal1->add_child(leaf_e);
            root->add_child(internal0);
            root->add_child(internal1);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(1.0 / 0.2)).epsilon(1e-8));

            double ln_hastings = op.propose_split(rng,
                    &tree);
            REQUIRE(tree.get_number_of_node_heights() == 3);
            REQUIRE(tree.get_number_of_splittable_heights() == 1);
            tree.compute_log_likelihood_and_prior(true);

            splits = tree.get_splits_by_height_index();
            if (choose_ab.count(splits) > 0) {
                /* std::cout << "Proposed move AB\n"; */
                REQUIRE(ln_hastings == Approx(choose_ab_hr).epsilon(1e-8));
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            (1.0 / 0.2) * (1.0 / 0.2)
                        )).epsilon(1e-8));
                /* std::cout << "Passed move AB\n"; */
            }
            else if (choose_cde.count(splits) > 0) {
                /* std::cout << "Proposed move CDE\n"; */
                REQUIRE(ln_hastings == Approx(choose_cde_hr).epsilon(1e-8));
                if (tree.get_mapped_nodes(1).size() > 1) {
                    REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                                (1.0 / 0.2) * (1.0 / 0.1)
                            )).epsilon(1e-8));
                }
                else {
                    REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                                (1.0 / 0.2) * (1.0 / 0.2)
                            )).epsilon(1e-8));
                }
                /* std::cout << "Passed move CDE\n"; */
            }
            else if (choose_both.count(splits) > 0) {
                /* std::cout << "Proposed move BOTH\n"; */
                REQUIRE(ln_hastings == Approx(choose_both_hr).epsilon(1e-8));
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            (1.0 / 0.2) * (1.0 / 0.1)
                        )).epsilon(1e-8));
                /* std::cout << "Passed move BOTH\n"; */
            }
            else {
                std::cout << "Unexpected tree\n";
                REQUIRE(0 == 1);
            }
            ++counts[splits];
        }

        for (auto splits_count : counts) {
            if (choose_ab.count(splits_count.first) > 0) {
                std::cout << "Expected freq: " << 1.0/3.0 << "\n";
                std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
            }
            else if (choose_cde.count(splits_count.first) > 0) {
                std::cout << "Expected freq: " << 1.0/12.0 << "\n";
                std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
            }
            else if (choose_both.count(splits_count.first) > 0) {
                std::cout << "Expected freq: " << 1.0/9.0 << "\n";
                std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
            }
            else {
                std::cout << "Unexpected tree\n";
                REQUIRE(0 == 1);
            }
        }

        double eps = 0.005;
        for (auto splits_count : counts) {
            if (choose_ab.count(splits_count.first) > 0) {
                REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/3.0).epsilon(eps));
            }
            else if (choose_cde.count(splits_count.first) > 0) {
                REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/12.0).epsilon(eps));
            }
            else if (choose_both.count(splits_count.first) > 0) {
                REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/9.0).epsilon(eps));
            }
            else {
                std::cout << "Unexpected tree\n";
                REQUIRE(0 == 1);
            }
        }
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler::propose from 2-1-1 shared tree with 5 leaves",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing merge 2-1-1 shared tree") {
        RandomNumberGenerator rng = RandomNumberGenerator(25179271593);

        // MOVE FROM: ((A:0.05,B:0.05)*:0.15,(C:0.1,(D:0.05,E:0.05)*:0.05):0.1)
        // MOVE TO:   ((A:0.1,B:0.1)*:0.1,(C:0.1,D:0.1,E:0.1)*:0.1)
        //            OR
        //            ((A:0.05,B:0.05)*:0.15,(D:0.05,E:0.05)*:0.15,C:0.2)
        //            OR
        //            ((A:0.02,B:0.02):0.18,(C:0.1,(D:0.05,E:0.05):0.05):0.1)
        //            ((A:0.05,B:0.05):0.15,(C:0.1,(D:0.02,E:0.02):0.08):0.1)
        //
        // MOVE TO:   ((A:0.1,B:0.1)*:0.1,(C:0.1,D:0.1,E:0.1)*:0.1)
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to merge) * pr(choose height)
        //                  =  1/2 * 1/2 = 1/4
        // pr(reverse move) = pr(choose to split) * pr(choose height) * pr(choose to move down) * pr(partition poly) * pr(new height)
        //                  = 1/2 * 1 * 1/3 * 1/3 * 1/0.1
        //                  = 1/(18*0.1) = 1/1.8
        // HR = 1/1.8 / 1/4 = 4/1.8
        //
        // MOVE TO:   ((A:0.05,B:0.05)*:0.15,(D:0.05,E:0.05)*:0.15,C:0.2)
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to merge) * pr(choose height)
        //                  =  1/2 * 1/2 = 1/4
        // pr(reverse move) = pr(choose to split) * pr(choose height) * pr(partition poly) * pr(new height)
        //                  = 1/2 * 1/2 * 1/3 * 1/(0.2-0.05)
        //                  = 1/2 * 1/2 * 1/3 * 1/(0.15)
        //                  = 1/(12*0.15) = 1/1.8
        // HR = 1/1.8 / 1/4 = 4/1.8
        //
        // MOVE TO:   ((A:0.02,B:0.02):0.18,(C:0.1,(D:0.05,E:0.05):0.05):0.1)
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to split) * pr(choose height) * pr(choose move node) * pr(new height)
        //                  =  1/2 * 1 * 1/2 * 1/0.05 = 1/4 * 20 = 5
        // pr(reverse move) = pr(choose to merge) * pr(choose height)
        //                  = 1 * 1/3 = 1/6
        // HR = 1/3 / 5 = 1/3 * 1/5 = 1/15
        
        BaseTree<Node> tree_ab_cde("((A:0.1,B:0.1)[&height_index=0,height=0.1]:0.1,(C:0.1,D:0.1,E:0.1)[&height_index=0,height=0.1]:0.1)[&height_index=1,height=0.2];");
        BaseTree<Node> tree_ab_de_c("((A:0.05,B:0.05)[&height_index=0,height=0.05]:0.15,C:0.2,(D:0.05,E:0.05)[&height_index=0,height=0.05]:0.15)[&height_index=1,height=0.2];");
        BaseTree<Node> tree_gen_ab("((A:0.02,B:0.02)[&height_index=0,height=0.02]:0.18,(C:0.1,(D:0.05,E:0.05)[&height_index=1,height=0.05]:0.05)[&height_index=2,height=0.1]:0.1)[&height_index=3,height=0.2];");
        BaseTree<Node> tree_gen_de("((A:0.05,B:0.05)[&height_index=1,height=0.05]:0.15,(C:0.1,(D:0.02,E:0.02)[&height_index=0,height=0.02]:0.08)[&height_index=2,height=0.1]:0.1)[&height_index=3,height=0.2];");

        std::set<std::map< unsigned int, std::set<Split> > > merge_0;
        merge_0.insert(tree_ab_cde.get_splits_by_height_index());
        double merge_0_hr = std::log(4.0 / 1.8);

        std::set<std::map< unsigned int, std::set<Split> > > merge_1;
        merge_1.insert(tree_ab_de_c.get_splits_by_height_index());
        double merge_1_hr = std::log(4.0 / 1.8);

        std::set<std::map< unsigned int, std::set<Split> > > split_ab;
        split_ab.insert(tree_gen_ab.get_splits_by_height_index());
        double split_hr = std::log(1.0 / 15.0);

        std::set<std::map< unsigned int, std::set<Split> > > split_de;
        split_ab.insert(tree_gen_de.get_splits_by_height_index());

        std::map<std::map< unsigned int, std::set<Split> >, unsigned int> counts;
        for (auto splits : merge_0) {
            counts[splits] = 0;
        }
        for (auto splits : merge_1) {
            counts[splits] = 0;
        }
        for (auto splits : split_ab) {
            counts[splits] = 0;
        }
        for (auto splits : split_de) {
            counts[splits] = 0;
        }

        std::map< unsigned int, std::set<Split> > splits;
        unsigned int nsamples = 50000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.2;
            std::shared_ptr<Node> root = std::make_shared<Node>(8, "root", root_ht);
            std::shared_ptr<Node> internal_ab = std::make_shared<Node>(7, "internal_ab", 0.05);
            std::shared_ptr<Node> internal_de = std::make_shared<Node>(6, "internal_de", 0.05);
            std::shared_ptr<Node> internal_cde = std::make_shared<Node>(5, "internal_cde", 0.1);
            internal_ab->set_height_parameter(internal_de->get_height_parameter());
            std::shared_ptr<Node> leaf_a = std::make_shared<Node>(0, "A", 0.0);
            std::shared_ptr<Node> leaf_b = std::make_shared<Node>(1, "B", 0.0);
            std::shared_ptr<Node> leaf_c = std::make_shared<Node>(2, "C", 0.0);
            std::shared_ptr<Node> leaf_d = std::make_shared<Node>(3, "D", 0.0);
            std::shared_ptr<Node> leaf_e = std::make_shared<Node>(4, "E", 0.0);

            internal_ab->add_child(leaf_a);
            internal_ab->add_child(leaf_b);
            internal_de->add_child(leaf_d);
            internal_de->add_child(leaf_e);
            internal_cde->add_child(leaf_c);
            internal_cde->add_child(internal_de);
            root->add_child(internal_ab);
            root->add_child(internal_cde);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            (1.0 / 0.2) * (1.0 / 0.1)
                            )).epsilon(1e-8));
            REQUIRE(tree.get_number_of_node_heights() == 3);

            double ln_hastings = op.propose(rng,
                    &tree);
            tree.compute_log_likelihood_and_prior(true);

            splits = tree.get_splits_by_height_index();
            if (merge_0.count(splits) > 0) {
                /* std::cout << "Proposed merge 0\n"; */
                REQUIRE(tree.get_number_of_node_heights() == 2);
                REQUIRE(tree.get_number_of_splittable_heights() == 1);
                REQUIRE(ln_hastings == Approx(merge_0_hr).epsilon(1e-8));
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            1.0 / 0.2
                        )).epsilon(1e-8));
                /* std::cout << "Passed merge 0\n"; */
            }
            else if (merge_1.count(splits) > 0) {
                /* std::cout << "Proposed merge 1\n"; */
                REQUIRE(tree.get_number_of_node_heights() == 2);
                REQUIRE(tree.get_number_of_splittable_heights() == 2);
                REQUIRE(ln_hastings == Approx(merge_1_hr).epsilon(1e-8));
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            1.0 / 0.2
                        )).epsilon(1e-8));
                /* std::cout << "Passed merge 1\n"; */
            }
            else if (split_ab.count(splits) > 0) {
                /* std::cout << "Proposed split AB\n"; */
                REQUIRE(tree.get_number_of_node_heights() == 4);
                REQUIRE(tree.get_number_of_splittable_heights() == 0);
                REQUIRE(ln_hastings == Approx(split_hr).epsilon(1e-8));
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            (1.0 / 0.2) * (1.0 / 0.2) * (1.0 / 0.1) 
                        )).epsilon(1e-8));
                /* std::cout << "Passed split AB\n"; */
            }
            else if (split_de.count(splits) > 0) {
                /* std::cout << "Proposed split DE\n"; */
                REQUIRE(tree.get_number_of_node_heights() == 4);
                REQUIRE(tree.get_number_of_splittable_heights() == 0);
                REQUIRE(ln_hastings == Approx(split_hr).epsilon(1e-8));
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            (1.0 / 0.2) * (1.0 / 0.2) * (1.0 / 0.1) 
                        )).epsilon(1e-8));
                /* std::cout << "Passed merge DE\n"; */
            }
            else {
                std::cout << "Unexpected tree proposed\n";
                for (auto height_splits : splits) {
                    std::cout << "  " << height_splits.first << ": ";
                    unsigned int s_count = 0;
                    for (auto split : height_splits.second) {
                        if (s_count > 0) {
                            std::cout << "   ";
                        }
                        std::cout << "  " << split.as_string() << "\n";
                        ++s_count;
                    }
                }
                REQUIRE(0 == 1);
            }
            ++counts[splits];
        }

        for (auto splits_count : counts) {
            if (merge_0.count(splits_count.first) > 0) {
                std::cout << "Expected freq: " << 1.0/4.0 << "\n";
                std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
            }
            else if (merge_1.count(splits_count.first) > 0) {
                std::cout << "Expected freq: " << 1.0/4.0 << "\n";
                std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
            }
            else if (split_ab.count(splits_count.first) > 0) {
                std::cout << "Expected freq: " << 1.0/4.0 << "\n";
                std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
            }
            else if (split_de.count(splits_count.first) > 0) {
                std::cout << "Expected freq: " << 1.0/4.0 << "\n";
                std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
            }
            else {
                std::cout << "Unexpected tree sampled " << splits_count.second << "times:\n";
                for (auto height_splits : splits_count.first) {
                    std::cout << "  " << height_splits.first << ": ";
                    unsigned int s_count = 0;
                    for (auto split : height_splits.second) {
                        if (s_count > 0) {
                            std::cout << "   ";
                        }
                        std::cout << "  " << split.as_string() << "\n";
                        ++s_count;
                    }
                }
                REQUIRE(0 == 1);
            }
        }

        double eps = 0.005;
        for (auto splits_count : counts) {
            if (merge_0.count(splits_count.first) > 0) {
                REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/4.0).epsilon(eps));
            }
            else if (merge_1.count(splits_count.first) > 0) {
                REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/4.0).epsilon(eps));
            }
            else if (split_ab.count(splits_count.first) > 0) {
                REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/4.0).epsilon(eps));
            }
            else if (split_de.count(splits_count.first) > 0) {
                REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/4.0).epsilon(eps));
            }
            else {
                std::cout << "Unexpected tree\n";
                REQUIRE(0 == 1);
            }
        }
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler::propose from 1-2-1 shared tree with 5 leaves",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing merge 2-1-1 shared tree") {
        RandomNumberGenerator rng = RandomNumberGenerator(25179271593);

        // MOVE FROM: ((A:0.1,B:0.1)*:0.1,(C:0.1,(D:0.05,E:0.05):0.05)*:0.1)
        // MERGE TO:  ((A:0.1,B:0.1)*:0.1,(C:0.1,D:0.1,E:0.1)*:0.1)
        //            OR
        //            (A:0.2,B:0.2,C:0.2,(D:0.05,E:0.05):0.15)
        //            OR
        // SPLIT TO:  ((A:0.08,B:0.08):0.12,(C:0.1,(D:0.05,E:0.05):0.05)*:0.1)
        //            OR
        //            ((A:0.1,B:0.1):0.1,(C:0.08,(D:0.05,E:0.05):0.03)*:0.12)
        //
        // MERGE TO:  ((A:0.1,B:0.1)*:0.1,(C:0.1,D:0.1,E:0.1)*:0.1)
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to merge) * pr(choose height)
        //                  =  1/2 * 1/2 = 1/4
        // pr(reverse move) = pr(choose to split) * pr(choose height) * pr(choose to move down) * pr(partition poly) * pr(new height)
        //                  = 1/2 * 1 * 1/3 * 1/4 * 1/0.1
        //                  = 1/(24*0.1) = 1/2.4
        // HR = 1/2.4 / 1/4 = 4/2.4
        //
        // MERGE TO:  (A:0.2,B:0.2,C:0.2,(D:0.05,E:0.05):0.15)
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to merge) * pr(choose height)
        //                  =  1/2 * 1/2 = 1/4
        // pr(reverse move) = pr(choose to split) * pr(choose height) * pr(choose to move down) * pr(partition poly) * pr(new height)
        //                  = 1/2 * 1 * 1 * 1/13 * 1/(0.2-0.05)
        //                  = 1/2 * 1 * 1 * 1/13 * 1/0.15
        //                  = 1/(26*0.15) = 1/3.9
        // HR = 1/3.9 / 1/4 = 4/3.9
        //
        // SPLIT TO:  ((A:0.08,B:0.08):0.12,(C:0.1,(D:0.05,E:0.05):0.05)*:0.1)
        // HR = pr(reverse move) / pr(forward move)
        // pr(reverse move) = pr(choose to merge) * pr(choose height)
        //                  =  1 * 1/3 = 1/3
        // pr(forward move) = pr(choose to split) * pr(choose height) * pr(choose to move down) * pr(new height)
        //                  = 1/2 * 1 * 1/2 * 1/(0.1-0.05)
        //                  = 1/2 * 1 * 1/2 * 1/(0.05)
        //                  = 1/(4*0.05) = 1/0.2 = 5
        // HR = 1/3 / 5 = 1/15
        //
        // SPLIT TO:  ((A:0.1,B:0.1):0.1,(C:0.08,(D:0.05,E:0.05):0.03)*:0.12)
        // HR = pr(reverse move) / pr(forward move)
        // pr(reverse move) = pr(choose to merge) * pr(choose height)
        //                  =  1 * 1/3 = 1/3
        // pr(forward move) = pr(choose to split) * pr(choose height) * pr(choose to move down) * pr(new height)
        //                  = 1/2 * 1 * 1/2 * 1/(0.1-0.05)
        //                  = 1/2 * 1 * 1/2 * 1/(0.05)
        //                  = 1/(4*0.05) = 1/0.2 = 5
        // HR = 1/3 / 5 = 1/15
        
        BaseTree<Node> tree_merge_0("((A:0.1,B:0.1)[&height_index=0,height=0.1]:0.1,(C:0.1,D:0.1,E:0.1)[&height_index=0,height=0.1]:0.1)[&height_index=1,height=0.2];");
        BaseTree<Node> tree_merge_1("(A:0.2,B:0.2,C:0.2,(D:0.05,E:0.05)[&height_index=0,height=0.05]:0.15)[&height_index=1,height=0.2];");
        BaseTree<Node> tree_split_ab("((A:0.08,B:0.08)[&height_index=1,height=0.08]:0.12,(C:0.1,(D:0.05,E:0.05)[&height_index=0,height=0.05]:0.05)[&height_index=2,height=0.1]:0.1)[&height_index=3,height=0.2];");
        BaseTree<Node> tree_split_cde("((A:0.1,B:0.1)[&height_index=2,height=0.1]:0.1,(C:0.08,(D:0.05,E:0.05)[&height_index=0,height=0.05]:0.03)[&height_index=1,height=0.08]:0.12)[&height_index=3,height=0.2];");

        std::set<std::map< unsigned int, std::set<Split> > > merge_0;
        merge_0.insert(tree_merge_0.get_splits_by_height_index());
        double merge_0_hr = std::log(4.0 / 2.4);

        std::set<std::map< unsigned int, std::set<Split> > > merge_1;
        merge_1.insert(tree_merge_1.get_splits_by_height_index());
        double merge_1_hr = std::log(4.0 / 3.9);

        std::set<std::map< unsigned int, std::set<Split> > > split_ab;
        split_ab.insert(tree_split_ab.get_splits_by_height_index());
        double split_hr = std::log(1.0 / 15.0);

        std::set<std::map< unsigned int, std::set<Split> > > split_cde;
        split_cde.insert(tree_split_cde.get_splits_by_height_index());

        std::map<std::map< unsigned int, std::set<Split> >, unsigned int> counts;
        for (auto splits : merge_0) {
            counts[splits] = 0;
        }
        for (auto splits : merge_1) {
            counts[splits] = 0;
        }
        for (auto splits : split_ab) {
            counts[splits] = 0;
        }
        for (auto splits : split_cde) {
            counts[splits] = 0;
        }

        std::map< unsigned int, std::set<Split> > splits;
        unsigned int nsamples = 50000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.2;
            std::shared_ptr<Node> root = std::make_shared<Node>(8, "root", root_ht);
            std::shared_ptr<Node> internal_ab = std::make_shared<Node>(7, "internal_ab", 0.1);
            std::shared_ptr<Node> internal_de = std::make_shared<Node>(6, "internal_de", 0.05);
            std::shared_ptr<Node> internal_cde = std::make_shared<Node>(5, "internal_cde", 0.1);
            internal_ab->set_height_parameter(internal_cde->get_height_parameter());
            std::shared_ptr<Node> leaf_a = std::make_shared<Node>(0, "A", 0.0);
            std::shared_ptr<Node> leaf_b = std::make_shared<Node>(1, "B", 0.0);
            std::shared_ptr<Node> leaf_c = std::make_shared<Node>(2, "C", 0.0);
            std::shared_ptr<Node> leaf_d = std::make_shared<Node>(3, "D", 0.0);
            std::shared_ptr<Node> leaf_e = std::make_shared<Node>(4, "E", 0.0);

            internal_ab->add_child(leaf_a);
            internal_ab->add_child(leaf_b);
            internal_de->add_child(leaf_d);
            internal_de->add_child(leaf_e);
            internal_cde->add_child(leaf_c);
            internal_cde->add_child(internal_de);
            root->add_child(internal_ab);
            root->add_child(internal_cde);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            (1.0 / 0.2) * (1.0 / 0.1)
                            )).epsilon(1e-8));
            REQUIRE(tree.get_number_of_node_heights() == 3);

            double ln_hastings = op.propose(rng,
                    &tree);
            tree.compute_log_likelihood_and_prior(true);

            splits = tree.get_splits_by_height_index();
            if (merge_0.count(splits) > 0) {
                /* std::cout << "Proposed merge 0\n"; */
                REQUIRE(tree.get_number_of_node_heights() == 2);
                REQUIRE(tree.get_number_of_splittable_heights() == 1);
                REQUIRE(ln_hastings == Approx(merge_0_hr).epsilon(1e-8));
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            1.0 / 0.2
                        )).epsilon(1e-8));
                /* std::cout << "Passed merge 0\n"; */
            }
            else if (merge_1.count(splits) > 0) {
                /* std::cout << "Proposed merge 1\n"; */
                REQUIRE(tree.get_number_of_node_heights() == 2);
                REQUIRE(tree.get_number_of_splittable_heights() == 1);
                REQUIRE(ln_hastings == Approx(merge_1_hr).epsilon(1e-8));
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            1.0 / 0.2
                        )).epsilon(1e-8));
                /* std::cout << "Passed merge 1\n"; */
            }
            else if (split_ab.count(splits) > 0) {
                /* std::cout << "Proposed split AB\n"; */
                REQUIRE(tree.get_number_of_node_heights() == 4);
                REQUIRE(tree.get_number_of_splittable_heights() == 0);
                REQUIRE(ln_hastings == Approx(split_hr).epsilon(1e-8));
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            (1.0 / 0.2) * (1.0 / 0.2) * (1.0 / 0.1) 
                        )).epsilon(1e-8));
                /* std::cout << "Passed split AB\n"; */
            }
            else if (split_cde.count(splits) > 0) {
                /* std::cout << "Proposed split DE\n"; */
                REQUIRE(tree.get_number_of_node_heights() == 4);
                REQUIRE(tree.get_number_of_splittable_heights() == 0);
                REQUIRE(ln_hastings == Approx(split_hr).epsilon(1e-8));
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            (1.0 / 0.2) * (1.0 / 0.2) * (1.0 / tree.get_height(1))
                        )).epsilon(1e-8));
                /* std::cout << "Passed merge DE\n"; */
            }
            else {
                std::cout << "Unexpected tree proposed\n";
                for (auto height_splits : splits) {
                    std::cout << "  " << height_splits.first << ": ";
                    unsigned int s_count = 0;
                    for (auto split : height_splits.second) {
                        if (s_count > 0) {
                            std::cout << "   ";
                        }
                        std::cout << "  " << split.as_string() << "\n";
                        ++s_count;
                    }
                }
                REQUIRE(0 == 1);
            }
            ++counts[splits];
        }

        for (auto splits_count : counts) {
            if (merge_0.count(splits_count.first) > 0) {
                std::cout << "Expected freq: " << 1.0/4.0 << "\n";
                std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
            }
            else if (merge_1.count(splits_count.first) > 0) {
                std::cout << "Expected freq: " << 1.0/4.0 << "\n";
                std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
            }
            else if (split_ab.count(splits_count.first) > 0) {
                std::cout << "Expected freq: " << 1.0/4.0 << "\n";
                std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
            }
            else if (split_cde.count(splits_count.first) > 0) {
                std::cout << "Expected freq: " << 1.0/4.0 << "\n";
                std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
            }
            else {
                std::cout << "Unexpected tree sampled " << splits_count.second << "times:\n";
                for (auto height_splits : splits_count.first) {
                    std::cout << "  " << height_splits.first << ": ";
                    unsigned int s_count = 0;
                    for (auto split : height_splits.second) {
                        if (s_count > 0) {
                            std::cout << "   ";
                        }
                        std::cout << "  " << split.as_string() << "\n";
                        ++s_count;
                    }
                }
                REQUIRE(0 == 1);
            }
        }

        double eps = 0.005;
        for (auto splits_count : counts) {
            if (merge_0.count(splits_count.first) > 0) {
                REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/4.0).epsilon(eps));
            }
            else if (merge_1.count(splits_count.first) > 0) {
                REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/4.0).epsilon(eps));
            }
            else if (split_ab.count(splits_count.first) > 0) {
                REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/4.0).epsilon(eps));
            }
            else if (split_cde.count(splits_count.first) > 0) {
                REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/4.0).epsilon(eps));
            }
            else {
                std::cout << "Unexpected tree\n";
                REQUIRE(0 == 1);
            }
        }
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler::propose from ((A,B),(C,D,E)*) tree",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing merge ((A,B),(C,D,E)*) tree") {
        RandomNumberGenerator rng = RandomNumberGenerator(14875215);

        // MOVE FROM: ((A:0.1,B:0.1):0.1,(C:0.05,D:0.05,E:0.05):0.15)
        // MERGE TO:  ((A:0.1,B:0.1)*:0.1,(C:0.1,D:0.1,E:0.1)*:0.1)
        //            OR
        //            (A:0.2,B:0.2,(C:0.05,D:0.05,E:0.05):0.15)
        //            OR
        // SPLIT TO:  ((A:0.1,B:0.1):0.1,(C:0.05,(D:0.02,E:0.02)):0.15)
        //            OR
        //            ((A:0.1,B:0.1):0.1,(D:0.05,(C:0.02,E:0.02)):0.15)
        //            OR
        //            ((A:0.1,B:0.1):0.1,(E:0.05,(C:0.02,D:0.02)):0.15)
        //
        // MERGE TO:  ((A:0.1,B:0.1)*:0.1,(C:0.1,D:0.1,E:0.1)*:0.1)
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to merge) * pr(choose height)
        //                  =  1/2 * 1/2 = 1/4
        // pr(reverse move) = pr(choose to split) * pr(choose height) * pr(choose to move down) * pr(partition poly) * pr(new height)
        //                  = 1/2 * 1 * 1/3 * 1/4 * 1/0.1
        //                  = 1/(24*0.1) = 1/2.4
        // HR = 1/2.4 / 1/4 = 4/2.4
        //
        // MERGE TO:  (A:0.2,B:0.2,(C:0.05,D:0.05,E:0.05):0.15)
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to merge) * pr(choose height)
        //                  =  1/2 * 1/2 = 1/4
        // pr(reverse move) = pr(choose to split) * pr(choose height) * pr(choose to move down) * pr(partition poly) * pr(new height)
        //                  = 1/2 * 1/2 * 1 * 1/3 * 1/(0.2-0.05)
        //                  = 1/2 * 1/2 * 1 * 1/3 * 1/0.15
        //                  = 1/(12*0.15) = 1/1.8
        // HR = 1/1.8 / 1/4 = 4/1.8
        //
        // SPLIT TO:  ((A:0.1,B:0.1):0.1,(C:0.05,(D:0.02,E:0.02)):0.15) [OR OTHER 2]
        // HR = pr(reverse move) / pr(forward move)
        // pr(reverse move) = pr(choose to merge) * pr(choose height)
        //                  =  1 * 1/3 = 1/3
        // pr(forward move) = pr(choose to split) * pr(choose height) * pr(choose to move down) * pr(partition poly) * pr(new height)
        //                  = 1/2 * 1 * 1 * 1/3 * 1/(0.05)
        //                  = 1/(6*0.05) = 1/0.3
        // HR = 1/3 / 1/0.3 = 0.3/3 = 0.1
        
        BaseTree<Node> tree_merge_0("((A:0.1,B:0.1)[&height_index=0,height=0.1]:0.1,(C:0.1,D:0.1,E:0.1)[&height_index=0,height=0.1]:0.1)[&height_index=1,height=0.2];");
        BaseTree<Node> tree_merge_1("(A:0.2,B:0.2,(C:0.05,D:0.05,E:0.05)[&height_index=0,height=0.05]:0.15)[&height_index=1,height=0.2];");
        BaseTree<Node> tree_split_de("((A:0.1,B:0.1)[&height_index=2,height=0.1]:0.1,(C:0.05,(D:0.02,E:0.02)[&height_index=0,height=0.02]:0.03)[&height_index=1,height=0.05]:0.15)[&height_index=3,height=0.2];");
        BaseTree<Node> tree_split_ce("((A:0.1,B:0.1)[&height_index=2,height=0.1]:0.1,(D:0.05,(C:0.02,E:0.02)[&height_index=0,height=0.02]:0.03)[&height_index=1,height=0.05]:0.15)[&height_index=3,height=0.2];");
        BaseTree<Node> tree_split_cd("((A:0.1,B:0.1)[&height_index=2,height=0.1]:0.1,(E:0.05,(C:0.02,D:0.02)[&height_index=0,height=0.02]:0.03)[&height_index=1,height=0.05]:0.15)[&height_index=3,height=0.2];");

        std::set<std::map< unsigned int, std::set<Split> > > merge_0;
        merge_0.insert(tree_merge_0.get_splits_by_height_index());
        double merge_0_hr = std::log(4.0 / 2.4);

        std::set<std::map< unsigned int, std::set<Split> > > merge_1;
        merge_1.insert(tree_merge_1.get_splits_by_height_index());
        double merge_1_hr = std::log(4.0 / 1.8);

        std::set<std::map< unsigned int, std::set<Split> > > split_poly;
        split_poly.insert(tree_split_de.get_splits_by_height_index());
        split_poly.insert(tree_split_ce.get_splits_by_height_index());
        split_poly.insert(tree_split_cd.get_splits_by_height_index());
        double split_hr = std::log(0.1);

        std::map<std::map< unsigned int, std::set<Split> >, unsigned int> counts;
        for (auto splits : merge_0) {
            counts[splits] = 0;
        }
        for (auto splits : merge_1) {
            counts[splits] = 0;
        }
        for (auto splits : split_poly) {
            counts[splits] = 0;
        }

        std::map< unsigned int, std::set<Split> > splits;
        unsigned int nsamples = 50000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.2;
            std::shared_ptr<Node> root = std::make_shared<Node>(7, "root", root_ht);
            std::shared_ptr<Node> internal_ab = std::make_shared<Node>(6, "internal_ab", 0.1);
            std::shared_ptr<Node> internal_cde = std::make_shared<Node>(5, "internal_cde", 0.05);
            std::shared_ptr<Node> leaf_a = std::make_shared<Node>(0, "A", 0.0);
            std::shared_ptr<Node> leaf_b = std::make_shared<Node>(1, "B", 0.0);
            std::shared_ptr<Node> leaf_c = std::make_shared<Node>(2, "C", 0.0);
            std::shared_ptr<Node> leaf_d = std::make_shared<Node>(3, "D", 0.0);
            std::shared_ptr<Node> leaf_e = std::make_shared<Node>(4, "E", 0.0);

            internal_ab->add_child(leaf_a);
            internal_ab->add_child(leaf_b);
            internal_cde->add_child(leaf_d);
            internal_cde->add_child(leaf_e);
            internal_cde->add_child(leaf_c);
            root->add_child(internal_ab);
            root->add_child(internal_cde);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            (1.0 / 0.2) * (1.0 / 0.2)
                            )).epsilon(1e-8));
            REQUIRE(tree.get_number_of_node_heights() == 3);

            double ln_hastings = op.propose(rng,
                    &tree);
            tree.compute_log_likelihood_and_prior(true);

            splits = tree.get_splits_by_height_index();
            if (merge_0.count(splits) > 0) {
                /* std::cout << "Proposed merge 0\n"; */
                REQUIRE(tree.get_number_of_node_heights() == 2);
                REQUIRE(tree.get_number_of_splittable_heights() == 1);
                REQUIRE(ln_hastings == Approx(merge_0_hr).epsilon(1e-8));
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            1.0 / 0.2
                        )).epsilon(1e-8));
                /* std::cout << "Passed merge 0\n"; */
            }
            else if (merge_1.count(splits) > 0) {
                /* std::cout << "Proposed merge 1\n"; */
                REQUIRE(tree.get_number_of_node_heights() == 2);
                REQUIRE(tree.get_number_of_splittable_heights() == 2);
                REQUIRE(ln_hastings == Approx(merge_1_hr).epsilon(1e-8));
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            1.0 / 0.2
                        )).epsilon(1e-8));
                /* std::cout << "Passed merge 1\n"; */
            }
            else if (split_poly.count(splits) > 0) {
                /* std::cout << "Proposed split AB\n"; */
                REQUIRE(tree.get_number_of_node_heights() == 4);
                REQUIRE(tree.get_number_of_splittable_heights() == 0);
                REQUIRE(ln_hastings == Approx(split_hr).epsilon(1e-8));
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            (1.0 / 0.2) * (1.0 / 0.2) * (1.0 / 0.05) 
                        )).epsilon(1e-8));
                /* std::cout << "Passed split AB\n"; */
            }
            else {
                std::cout << "Unexpected tree proposed\n";
                for (auto height_splits : splits) {
                    std::cout << "  " << height_splits.first << ": ";
                    unsigned int s_count = 0;
                    for (auto split : height_splits.second) {
                        if (s_count > 0) {
                            std::cout << "   ";
                        }
                        std::cout << "  " << split.as_string() << "\n";
                        ++s_count;
                    }
                }
                REQUIRE(0 == 1);
            }
            ++counts[splits];
        }

        for (auto splits_count : counts) {
            if (merge_0.count(splits_count.first) > 0) {
                std::cout << "Expected freq: " << 1.0/4.0 << "\n";
                std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
            }
            else if (merge_1.count(splits_count.first) > 0) {
                std::cout << "Expected freq: " << 1.0/4.0 << "\n";
                std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
            }
            else if (split_poly.count(splits_count.first) > 0) {
                std::cout << "Expected freq: " << 1.0/6.0 << "\n";
                std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
            }
            else {
                std::cout << "Unexpected tree sampled " << splits_count.second << "times:\n";
                for (auto height_splits : splits_count.first) {
                    std::cout << "  " << height_splits.first << ": ";
                    unsigned int s_count = 0;
                    for (auto split : height_splits.second) {
                        if (s_count > 0) {
                            std::cout << "   ";
                        }
                        std::cout << "  " << split.as_string() << "\n";
                        ++s_count;
                    }
                }
                REQUIRE(0 == 1);
            }
        }

        double eps = 0.005;
        for (auto splits_count : counts) {
            if (merge_0.count(splits_count.first) > 0) {
                REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/4.0).epsilon(eps));
            }
            else if (merge_1.count(splits_count.first) > 0) {
                REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/4.0).epsilon(eps));
            }
            else if (split_poly.count(splits_count.first) > 0) {
                REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/6.0).epsilon(eps));
            }
            else {
                std::cout << "Unexpected tree\n";
                REQUIRE(0 == 1);
            }
        }
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler::propose from ((A,B)*,(C,D,E)) tree",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing merge ((A,B)*,(C,D,E)) tree") {
        RandomNumberGenerator rng = RandomNumberGenerator(5421523145);

        // MOVE FROM: ((A:0.05,B:0.05):0.15,(C:0.1,D:0.1,E:0.1):0.1)
        // MERGE TO:  ((A:0.1,B:0.1)*:0.1,(C:0.1,D:0.1,E:0.1)*:0.1)
        //            OR
        //            ((A:0.05,B:0.05):0.15,C:0.2,D:0.2,E:0.2)
        //            OR
        // SPLIT TO:  ((A:0.05,B:0.05):0.15,(C:0.1,(D:0.08,E:0.08):0.02):0.1)
        //            OR
        //            ((A:0.05,B:0.05):0.15,(D:0.1,(C:0.08,E:0.08):0.02):0.1)
        //            OR
        //            ((A:0.05,B:0.05):0.15,(E:0.1,(D:0.08,C:0.08):0.02):0.1)
        //
        // MERGE TO:  ((A:0.1,B:0.1)*:0.1,(C:0.1,D:0.1,E:0.1)*:0.1)
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to merge) * pr(choose height)
        //                  =  1/2 * 1/2 = 1/4
        // pr(reverse move) = pr(choose to split) * pr(choose height) * pr(choose to move down) * pr(partition poly) * pr(new height)
        //                  = 1/2 * 1 * 1/3 * 1/0.1
        //                  = 1/(6*0.1) = 1/0.6
        // HR = 1/0.6 / 1/4 = 4/0.6 = 1/0.15
        //
        // MERGE TO:  ((A:0.05,B:0.05):0.15,C:0.2,D:0.2,E:0.2)
        // HR = pr(reverse move) / pr(forward move)
        // pr(forward move) = pr(choose to merge) * pr(choose height)
        //                  =  1/2 * 1/2 = 1/4
        // pr(reverse move) = pr(choose to split) * pr(choose height) * pr(choose to move down) * pr(partition poly) * pr(new height)
        //                  = 1/2 * 1 * 1 * 1/13 * 1/(0.2-0.05)
        //                  = 1/2 * 1 * 1 * 1/13 * 1/0.15
        //                  = 1/(26*0.15) = 1/3.9
        // HR = 1/3.9 / 1/4 = 4/3.9
        //
        // SPLIT TO:  ((A:0.05,B:0.05):0.15,(C:0.1,(D:0.08,E:0.08):0.02):0.1) [OR OTHER 2]
        // HR = pr(reverse move) / pr(forward move)
        // pr(reverse move) = pr(choose to merge) * pr(choose height)
        //                  =  1 * 1/3 = 1/3
        // pr(forward move) = pr(choose to split) * pr(choose height) * pr(choose to move down) * pr(partition poly) * pr(new height)
        //                  = 1/2 * 1 * 1 * 1/3 * 1/(0.1-0.05)
        //                  = 1/(6*0.05) = 1/0.3
        // HR = 1/3 / 1/0.3 = 0.3/3 = 0.1
        
        BaseTree<Node> tree_merge_0("((A:0.1,B:0.1)[&height_index=0,height=0.1]:0.1,(C:0.1,D:0.1,E:0.1)[&height_index=0,height=0.1]:0.1)[&height_index=1,height=0.2];");
        BaseTree<Node> tree_merge_1("((A:0.05,B:0.05)[&height_index=0,height=0.05]:0.15,C:0.2,D:0.2,E:0.2)[&height_index=1,height=0.2];");
        BaseTree<Node> tree_split_de("((A:0.05,B:0.05)[&height_index=0,height=0.05]:0.15,(C:0.1,(D:0.08,E:0.08)[&height_index=1,height=0.08]:0.02)[&height_index=2,height=0.1]:0.1)[&height_index=3,height=0.2];");
        BaseTree<Node> tree_split_ce("((A:0.05,B:0.05)[&height_index=0,height=0.05]:0.15,(D:0.1,(C:0.08,E:0.08)[&height_index=1,height=0.08]:0.02)[&height_index=2,height=0.1]:0.1)[&height_index=3,height=0.2];");
        BaseTree<Node> tree_split_cd("((A:0.05,B:0.05)[&height_index=0,height=0.05]:0.15,(E:0.1,(D:0.08,C:0.08)[&height_index=1,height=0.08]:0.02)[&height_index=2,height=0.1]:0.1)[&height_index=3,height=0.2];");

        std::set<std::map< unsigned int, std::set<Split> > > merge_0;
        merge_0.insert(tree_merge_0.get_splits_by_height_index());
        double merge_0_hr = std::log(1.0 / 0.15);

        std::set<std::map< unsigned int, std::set<Split> > > merge_1;
        merge_1.insert(tree_merge_1.get_splits_by_height_index());
        double merge_1_hr = std::log(4.0 / 3.9);

        std::set<std::map< unsigned int, std::set<Split> > > split_poly;
        split_poly.insert(tree_split_de.get_splits_by_height_index());
        split_poly.insert(tree_split_ce.get_splits_by_height_index());
        split_poly.insert(tree_split_cd.get_splits_by_height_index());
        double split_hr = std::log(0.1);

        std::map<std::map< unsigned int, std::set<Split> >, unsigned int> counts;
        for (auto splits : merge_0) {
            counts[splits] = 0;
        }
        for (auto splits : merge_1) {
            counts[splits] = 0;
        }
        for (auto splits : split_poly) {
            counts[splits] = 0;
        }

        std::map< unsigned int, std::set<Split> > splits;
        unsigned int nsamples = 50000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.2;
            std::shared_ptr<Node> root = std::make_shared<Node>(7, "root", root_ht);
            std::shared_ptr<Node> internal_ab = std::make_shared<Node>(6, "internal_ab", 0.05);
            std::shared_ptr<Node> internal_cde = std::make_shared<Node>(5, "internal_cde", 0.1);
            std::shared_ptr<Node> leaf_a = std::make_shared<Node>(0, "A", 0.0);
            std::shared_ptr<Node> leaf_b = std::make_shared<Node>(1, "B", 0.0);
            std::shared_ptr<Node> leaf_c = std::make_shared<Node>(2, "C", 0.0);
            std::shared_ptr<Node> leaf_d = std::make_shared<Node>(3, "D", 0.0);
            std::shared_ptr<Node> leaf_e = std::make_shared<Node>(4, "E", 0.0);

            internal_ab->add_child(leaf_a);
            internal_ab->add_child(leaf_b);
            internal_cde->add_child(leaf_d);
            internal_cde->add_child(leaf_e);
            internal_cde->add_child(leaf_c);
            root->add_child(internal_ab);
            root->add_child(internal_cde);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            (1.0 / 0.2) * (1.0 / 0.2)
                            )).epsilon(1e-8));
            REQUIRE(tree.get_number_of_node_heights() == 3);

            double ln_hastings = op.propose(rng,
                    &tree);
            tree.compute_log_likelihood_and_prior(true);

            splits = tree.get_splits_by_height_index();
            if (merge_0.count(splits) > 0) {
                /* std::cout << "Proposed merge 0\n"; */
                REQUIRE(tree.get_number_of_node_heights() == 2);
                REQUIRE(tree.get_number_of_splittable_heights() == 1);
                REQUIRE(ln_hastings == Approx(merge_0_hr).epsilon(1e-8));
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            1.0 / 0.2
                        )).epsilon(1e-8));
                /* std::cout << "Passed merge 0\n"; */
            }
            else if (merge_1.count(splits) > 0) {
                /* std::cout << "Proposed merge 1\n"; */
                REQUIRE(tree.get_number_of_node_heights() == 2);
                REQUIRE(tree.get_number_of_splittable_heights() == 1);
                REQUIRE(ln_hastings == Approx(merge_1_hr).epsilon(1e-8));
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            1.0 / 0.2
                        )).epsilon(1e-8));
                /* std::cout << "Passed merge 1\n"; */
            }
            else if (split_poly.count(splits) > 0) {
                /* std::cout << "Proposed split AB\n"; */
                REQUIRE(tree.get_number_of_node_heights() == 4);
                REQUIRE(tree.get_number_of_splittable_heights() == 0);
                REQUIRE(ln_hastings == Approx(split_hr).epsilon(1e-8));
                REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            (1.0 / 0.2) * (1.0 / 0.2) * (1.0 / 0.1) 
                        )).epsilon(1e-8));
                /* std::cout << "Passed split AB\n"; */
            }
            else {
                std::cout << "Unexpected tree proposed\n";
                for (auto height_splits : splits) {
                    std::cout << "  " << height_splits.first << ": ";
                    unsigned int s_count = 0;
                    for (auto split : height_splits.second) {
                        if (s_count > 0) {
                            std::cout << "   ";
                        }
                        std::cout << "  " << split.as_string() << "\n";
                        ++s_count;
                    }
                }
                REQUIRE(0 == 1);
            }
            ++counts[splits];
        }

        for (auto splits_count : counts) {
            if (merge_0.count(splits_count.first) > 0) {
                std::cout << "Expected freq: " << 1.0/4.0 << "\n";
                std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
            }
            else if (merge_1.count(splits_count.first) > 0) {
                std::cout << "Expected freq: " << 1.0/4.0 << "\n";
                std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
            }
            else if (split_poly.count(splits_count.first) > 0) {
                std::cout << "Expected freq: " << 1.0/6.0 << "\n";
                std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
            }
            else {
                std::cout << "Unexpected tree sampled " << splits_count.second << "times:\n";
                for (auto height_splits : splits_count.first) {
                    std::cout << "  " << height_splits.first << ": ";
                    unsigned int s_count = 0;
                    for (auto split : height_splits.second) {
                        if (s_count > 0) {
                            std::cout << "   ";
                        }
                        std::cout << "  " << split.as_string() << "\n";
                        ++s_count;
                    }
                }
                REQUIRE(0 == 1);
            }
        }

        double eps = 0.005;
        for (auto splits_count : counts) {
            if (merge_0.count(splits_count.first) > 0) {
                REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/4.0).epsilon(eps));
            }
            else if (merge_1.count(splits_count.first) > 0) {
                REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/4.0).epsilon(eps));
            }
            else if (split_poly.count(splits_count.first) > 0) {
                REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/6.0).epsilon(eps));
            }
            else {
                std::cout << "Unexpected tree\n";
                REQUIRE(0 == 1);
            }
        }
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler::propose from (A,B,C,D,E) tree",
        "[SplitLumpNodesRevJumpSampler]") {

    SECTION("Testing merge (A,B,C,D,E) tree") {
        RandomNumberGenerator rng = RandomNumberGenerator(458262135);

        // MOVE FROM: (A:0.1,B:0.1,C:0.1,D:0.1,E:0.1)
        // SPLIT TO:  One of 50 different ways to break up polytomy
        // HR = pr(reverse move) / pr(forward move)
        // pr(reverse move) = pr(choose to merge) * pr(choose height)
        //                  =  1/2 * 1 = 1/2
        // pr(forward move) = pr(choose to split) * pr(choose height) * pr(choose to move down) * pr(partition poly) * pr(new height)
        //                  = 1 * 1 * 1 * 1/50 * 1/0.1
        //                  = 1/5
        // HR = 1/2 / 1/5 = 5/2 = 2.5
        //
        std::map<std::map< unsigned int, std::set<Split> >, unsigned int> counts;

        std::map< unsigned int, std::set<Split> > splits;
        unsigned int nsamples = 50000;
        for (unsigned int i = 0; i < nsamples; ++i) {
            double root_ht = 0.1;
            std::shared_ptr<Node> root = std::make_shared<Node>(5, "root", root_ht);
            std::shared_ptr<Node> leaf_a = std::make_shared<Node>(0, "A", 0.0);
            std::shared_ptr<Node> leaf_b = std::make_shared<Node>(1, "B", 0.0);
            std::shared_ptr<Node> leaf_c = std::make_shared<Node>(2, "C", 0.0);
            std::shared_ptr<Node> leaf_d = std::make_shared<Node>(3, "D", 0.0);
            std::shared_ptr<Node> leaf_e = std::make_shared<Node>(4, "E", 0.0);

            root->add_child(leaf_a);
            root->add_child(leaf_b);
            root->add_child(leaf_d);
            root->add_child(leaf_e);
            root->add_child(leaf_c);

            BaseTree<Node> tree(root);

            tree.ignore_data();
            tree.fix_root_height();

            SplitLumpNodesRevJumpSampler<Node> op;

            // Initialize prior probs
            tree.compute_log_likelihood_and_prior(true);
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                            1.0
                            )).epsilon(1e-8));
            REQUIRE(tree.get_number_of_node_heights() == 1);
            double expected_ln_hastings = std::log(2.5);

            double ln_hastings = op.propose(rng,
                    &tree);
            tree.compute_log_likelihood_and_prior(true);

            REQUIRE(tree.get_number_of_node_heights() == 2);
            REQUIRE(tree.get_number_of_splittable_heights() > 0);
            REQUIRE(tree.get_number_of_splittable_heights() < 3);
            REQUIRE(ln_hastings == Approx(expected_ln_hastings).epsilon(1e-8));
            REQUIRE(tree.get_log_prior_density_value() == Approx(std::log(
                        1.0 / 0.1
                    )).epsilon(1e-8));

            splits = tree.get_splits_by_height_index();
            if (counts.count(splits) > 0) {
                ++counts[splits];
            }
            else {
                counts[splits] = 1;
            }
        }

        for (auto splits_count : counts) {
            std::cout << "Expected freq: " << 1.0/50.0 << "\n";
            std::cout << "Sample freq: " << splits_count.second / (double)nsamples << "\n";
        }

        REQUIRE(counts.size() == 50);

        double eps = 0.005;
        for (auto splits_count : counts) {
            REQUIRE((splits_count.second / (double)nsamples) == Approx(1.0/50.0).epsilon(eps));
        }
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler with 6 leaves and fixed root",
        "[xSplitLumpNodesRevJumpSampler]") {

    SECTION("Testing 6 leaves with fixed root") {
        /* RandomNumberGenerator rng = RandomNumberGenerator(25478465); */
        RandomNumberGenerator rng = RandomNumberGenerator(7892471234);

        double root_ht = 0.5;
        std::shared_ptr<Node> root = std::make_shared<Node>(6, "root", root_ht);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>(5, "leaf5", 0.0);

        root->add_child(leaf0);
        root->add_child(leaf1);
        root->add_child(leaf2);
        root->add_child(leaf3);
        root->add_child(leaf4);
        root->add_child(leaf5);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        SplitLumpNodesRevJumpSampler<Node> op;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int niterations = 1000000000;
        unsigned int sample_freq = 100;
        unsigned int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 100000;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                std::set< std::set<Split> > splits = tree.get_splits(false);
                if (split_counts.count(splits) > 0) {
                    ++split_counts[splits];
                }
                else {
                    split_counts[splits] = 1;
                }
                ++sample_count;
                if (sample_count % report_freq == 0) {
                    std::cout << "Sampled " << sample_count << " of " << nsamples << std::endl;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        // TODO: Figure this out
        unsigned int num_tree_models = split_counts.size();
        double exp_freq = 1.0/num_tree_models;
        double exp_count = nsamples/num_tree_models;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.2;
        unsigned int total_trees_sampled = 0;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            std::cout << "Tree:\n";
            for (auto splitset : s_c.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            std::cout << "  nsamples: " << s_c.second << "\n";
            std::cout << "  prop error: " << prop_error << "\n";
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }


        double quantile_chi_sq_5627_10 = 5763.4;
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        std::cout << "Chi-square(5627) 0.9 quantile: " << quantile_chi_sq_5627_10 << "\n";

        std::cout << "BAD SPLITS (proportional error > " << prop_error_threshold << ")\n";
        for (auto s_e : bad_splits) {
            std::cout << "\nTree:\n";
            for (auto splitset : s_e.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            std::cout << "  prop error: " << s_e.second << "\n";
        }

        write_r_script(split_counts, "../6-leaf-general-tree-test.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        // REQUIRE(split_counts.size() == ???);

        REQUIRE(chi_sq_test_statistic < quantile_chi_sq_5627_10);
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler with 7 leaves and fixed root",
        "[xSplitLumpNodesRevJumpSampler]") {

    SECTION("Testing 7 leaves with fixed root") {
        RandomNumberGenerator rng = RandomNumberGenerator(2347243665);

        double root_ht = 0.5;
        std::shared_ptr<Node> root = std::make_shared<Node>(7, "root", root_ht);
        std::shared_ptr<Node> leaf0 = std::make_shared<Node>(0, "leaf0", 0.0);
        std::shared_ptr<Node> leaf1 = std::make_shared<Node>(1, "leaf1", 0.0);
        std::shared_ptr<Node> leaf2 = std::make_shared<Node>(2, "leaf2", 0.0);
        std::shared_ptr<Node> leaf3 = std::make_shared<Node>(3, "leaf3", 0.0);
        std::shared_ptr<Node> leaf4 = std::make_shared<Node>(4, "leaf4", 0.0);
        std::shared_ptr<Node> leaf5 = std::make_shared<Node>(5, "leaf5", 0.0);
        std::shared_ptr<Node> leaf6 = std::make_shared<Node>(6, "leaf6", 0.0);

        root->add_child(leaf0);
        root->add_child(leaf1);
        root->add_child(leaf2);
        root->add_child(leaf3);
        root->add_child(leaf4);
        root->add_child(leaf5);
        root->add_child(leaf6);

        BaseTree<Node> tree(root);

        tree.ignore_data();
        tree.fix_root_height();

        SplitLumpNodesRevJumpSampler<Node> op;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned long long int niterations = 5000000000;
        unsigned int sample_freq = 100;
        unsigned long long int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 100000;
        for (unsigned long long int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                std::set< std::set<Split> > splits = tree.get_splits(false);
                if (split_counts.count(splits) > 0) {
                    ++split_counts[splits];
                }
                else {
                    split_counts[splits] = 1;
                }
                ++sample_count;
                if (sample_count % report_freq == 0) {
                    std::cout << "Sampled " << sample_count << " of " << nsamples << std::endl;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        // niterations is likely beyond limit of unsigned int
        // REQUIRE(op.get_number_of_attempts() == niterations);

        // TODO: Figure this out
        unsigned int num_tree_models = split_counts.size();
        double exp_freq = 1.0/num_tree_models;
        double exp_count = nsamples/num_tree_models;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.2;
        unsigned long long int total_trees_sampled = 0;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            std::cout << "Tree:\n";
            for (auto splitset : s_c.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            std::cout << "  nsamples: " << s_c.second << "\n";
            std::cout << "  prop error: " << prop_error << "\n";
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }

        /* double quantile_chi_sq_5627_10 = 5763.4; */
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        /* std::cout << "Chi-square(5627) 0.9 quantile: " << quantile_chi_sq_5627_10 << "\n"; */

        std::cout << "BAD SPLITS (proportional error > " << prop_error_threshold << ")\n";
        for (auto s_e : bad_splits) {
            std::cout << "\nTree:\n";
            for (auto splitset : s_e.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            std::cout << "  prop error: " << s_e.second << "\n";
        }

        write_r_script(split_counts, "../7-leaf-general-tree-test.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        // REQUIRE(split_counts.size() == ???);

        /* REQUIRE(chi_sq_test_statistic < quantile_chi_sq_5627_10); */
    }
}


// 5 leaf test with BasePopulationTree
// Expectations for tree with 5 leaves:
// ----------------------------------------------------------------------------
// # of unlabeled                            # of trees (labelings)      
//  topologies
// ----------------------------------------------------------------------------
// (a,b,c,d,e)                                                     = 1
// (a,(b,c,d,e))           = 5 choose 4                            = 5
// (a,b,(c,d,e))           = 5 choose 3                            = 10 
// (a,b,c,(d,e))           = 5 choose 2                            = 10
// !((a,b,c),(d,e))        = 5 choose 3                            = 10
// (a,(b,(c,d,e)))         = (5 choose 3) * 2!                     = 20 
// !((a,b),(c,d),e)        = ((5 choose 2) * (3 choose 2)) / 2     = 15
// !(((a,b),(c,d)),e       = ((5 choose 2) * (3 choose 2)) / 2     = 15
// !(((a,b),c),(d,e))      = (5 choose 2) * (3 choose 2)           = 30
// (a,(b,(c,(d,e))))       = (5 choose 2) * 3!                     = 60
// (a,b,(c,(d,e)))         = (5 choose 2) * (3 choose 2)           = 30
// (a,(b,c,(d,e)))         = (5 choose 2) * (3 choose 2)           = 30
// ----------------------------------------------------------------------------
// TOTAL                                                           = 236
//
// 236 matches Felsenstein 1978, but we need to account for shared node
// heights. The topologies above prefixed with '!' are topologies that have
// potentially shared node heights. For each shared node configuration of these
// topologies, we have to add that many additional trees, which we do below.
//
// ----------------------------------------------------------------------------
// # of unlabeled                            # of trees (labelings)      
//  topologies
// ----------------------------------------------------------------------------
// ((a,b,c)*,(d,e)*)       = 5 choose 3                            = 10
// ((a,b)*,(c,d)*,e)       = ((5 choose 2) * (3 choose 2)) / 2     = 15
// (((a,b)*,(c,d)*),e      = ((5 choose 2) * (3 choose 2)) / 2     = 15
// (((a,b),c)*,(d,e)*)     = (5 choose 2) * (3 choose 2)           = 30
// (((a,b)*,c),(d,e)*)     = (5 choose 2) * (3 choose 2)           = 30
// ----------------------------------------------------------------------------
// GRAND TOTAL # OF TREE MODELS                                    = 336
//
// The asterisks in the topologies above indicated shared node heights.
TEST_CASE("Testing SplitLumpNodesRevJumpSampler with BasePopulationTree, 5 leaves, fixed sizes, and fixed root",
        "[BasePopulationTree]") {

    SECTION("Testing 5 leaves with BasePopulationTree, fixed sizes") {
        RandomNumberGenerator rng = RandomNumberGenerator(4157514221);

        double pop_size_shape = 10.0;
        double pop_size_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        double root_ht = 0.5;
        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(5, "root", root_ht);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);
        std::shared_ptr<PopulationNode> leaf3 = std::make_shared<PopulationNode>(3, "leaf3", 0.0);
        std::shared_ptr<PopulationNode> leaf4 = std::make_shared<PopulationNode>(4, "leaf4", 0.0);

        root->add_child(leaf0);
        root->add_child(leaf1);
        root->add_child(leaf2);
        root->add_child(leaf3);
        root->add_child(leaf4);

        BasePopulationTree tree(root);

        tree.set_population_size_prior(pop_size_prior);
        tree.ignore_data();
        tree.fix_root_height();
        tree.constrain_population_sizes();
        tree.fix_population_sizes();

        SplitLumpNodesRevJumpSampler<PopulationNode> op;

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int count_nheights_1 = 0;
        unsigned int count_nheights_2 = 0;
        unsigned int count_nheights_3 = 0;
        unsigned int count_nheights_4 = 0;

        unsigned int niterations = 10000000;
        unsigned int sample_freq = 50;
        unsigned int nsamples = niterations / sample_freq;

        unsigned int sample_count = 0;
        unsigned int report_freq = 10000;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1);
            if ((i + 1) % sample_freq == 0) {
                if (tree.get_number_of_node_heights() == 1) {
                    ++count_nheights_1;
                }
                else if (tree.get_number_of_node_heights() == 2) {
                    ++count_nheights_2;
                }
                else if (tree.get_number_of_node_heights() == 3) {
                    ++count_nheights_3;
                }
                else if (tree.get_number_of_node_heights() == 4) {
                    ++count_nheights_4;
                }
                std::set< std::set<Split> > splits = tree.get_splits(false);
                if (split_counts.count(splits) > 0) {
                    ++split_counts[splits];
                }
                else {
                    split_counts[splits] = 1;
                }
                ++sample_count;
                if (sample_count % report_freq == 0) {
                    std::cout << "Sampled " << sample_count << " of " << nsamples << std::endl;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        REQUIRE((count_nheights_1 + count_nheights_2 + count_nheights_3 + count_nheights_4) == nsamples);

        double freq_nheights_1 = count_nheights_1 / (double)nsamples;
        double freq_nheights_2 = count_nheights_2 / (double)nsamples;
        double freq_nheights_3 = count_nheights_3 / (double)nsamples;
        double freq_nheights_4 = count_nheights_4 / (double)nsamples;

        double exp_freq = 1.0/336.0;
        double exp_count = nsamples/336.0;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.2;
        unsigned int total_trees_sampled = 0;
        std::map< std::set< std::set<Split> >, double> split_freqs;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            split_freqs[s_c.first] = s_c.second / (double)nsamples;
            std::cout << "Tree:\n";
            for (auto splitset : s_c.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            std::cout << "  nsamples: " << s_c.second << "\n";
            std::cout << "  prop error: " << prop_error << "\n";
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }

        double quantile_chi_sq_335_10 = 368.6;
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        std::cout << "Chi-square(335) 0.9 quantile: " << quantile_chi_sq_335_10 << "\n";

        std::cout << "BAD SPLITS (proportional error > " << prop_error_threshold << ")\n";
        for (auto s_e : bad_splits) {
            std::cout << "\nTree:\n";
            for (auto splitset : s_e.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            std::cout << "  prop error: " << s_e.second << "\n";
        }

        write_r_script(split_counts, "../5-leaf-general-tree-test-fixed-pop-sizes.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        REQUIRE(split_counts.size() == 336);

        double eps = 0.001;

        for (auto s_f : split_freqs) {
            REQUIRE(s_f.second == Approx(exp_freq).epsilon(eps));
        }

        REQUIRE(chi_sq_test_statistic < quantile_chi_sq_335_10);
    }
}

TEST_CASE("Testing SplitLumpNodesRevJumpSampler with BasePopulationTree, 5 leaves, constrained sizes, and fixed root",
        "[BasePopulationTree]") {

    SECTION("Testing 5 leaves with BasePopulationTree, constrained sizes") {
        RandomNumberGenerator rng = RandomNumberGenerator(4164626423);

        double pop_size_shape = 10.0;
        double pop_size_scale = 0.05;
        std::shared_ptr<ContinuousProbabilityDistribution> pop_size_prior = std::make_shared<GammaDistribution>(
                pop_size_shape,
                pop_size_scale);

        double root_ht = 0.5;
        std::shared_ptr<PopulationNode> root = std::make_shared<PopulationNode>(5, "root", root_ht);
        std::shared_ptr<PopulationNode> leaf0 = std::make_shared<PopulationNode>(0, "leaf0", 0.0);
        std::shared_ptr<PopulationNode> leaf1 = std::make_shared<PopulationNode>(1, "leaf1", 0.0);
        std::shared_ptr<PopulationNode> leaf2 = std::make_shared<PopulationNode>(2, "leaf2", 0.0);
        std::shared_ptr<PopulationNode> leaf3 = std::make_shared<PopulationNode>(3, "leaf3", 0.0);
        std::shared_ptr<PopulationNode> leaf4 = std::make_shared<PopulationNode>(4, "leaf4", 0.0);

        root->add_child(leaf0);
        root->add_child(leaf1);
        root->add_child(leaf2);
        root->add_child(leaf3);
        root->add_child(leaf4);

        BasePopulationTree tree(root);

        tree.set_population_size_prior(pop_size_prior);
        tree.ignore_data();
        tree.fix_root_height();
        tree.constrain_population_sizes();

        SplitLumpNodesRevJumpSampler<PopulationNode> op;

        GlobalPopSizeScaler op2;
        op2.turn_on_auto_optimize();
        op2.set_auto_optimize_delay(100);

        // Initialize prior probs
        tree.compute_log_likelihood_and_prior(true);

        std::map< std::set< std::set<Split> >, unsigned int> split_counts;

        unsigned int count_nheights_1 = 0;
        unsigned int count_nheights_2 = 0;
        unsigned int count_nheights_3 = 0;
        unsigned int count_nheights_4 = 0;

        unsigned int niterations = 10000000;
        unsigned int sample_freq = 50;
        unsigned int nsamples = niterations / sample_freq;

        SampleSummarizer<double> pop_size_summary;

        std::vector< std::shared_ptr<PositiveRealParameter> > pop_sizes;

        unsigned int sample_count = 0;
        unsigned int report_freq = 10000;
        for (unsigned int i = 0; i < niterations; ++i) {
            op.operate(rng, &tree, 1, 1);
            op2.operate(rng, &tree, 1, 1);
            if ((i + 1) % sample_freq == 0) {
                pop_sizes = tree.get_pointers_to_population_sizes();
                REQUIRE(pop_sizes.size() == 1);
                pop_size_summary.add_sample(pop_sizes.at(0)->get_value());

                if (tree.get_number_of_node_heights() == 1) {
                    ++count_nheights_1;
                }
                else if (tree.get_number_of_node_heights() == 2) {
                    ++count_nheights_2;
                }
                else if (tree.get_number_of_node_heights() == 3) {
                    ++count_nheights_3;
                }
                else if (tree.get_number_of_node_heights() == 4) {
                    ++count_nheights_4;
                }
                std::set< std::set<Split> > splits = tree.get_splits(false);
                if (split_counts.count(splits) > 0) {
                    ++split_counts[splits];
                }
                else {
                    split_counts[splits] = 1;
                }
                ++sample_count;
                if (sample_count % report_freq == 0) {
                    std::cout << "Sampled " << sample_count << " of " << nsamples << std::endl;
                }
            }
        }
        std::cout << op.header_string();
        std::cout << op.to_string();
        std::cout << op2.header_string();
        std::cout << op2.to_string();

        REQUIRE(op.get_number_of_attempts() == niterations);

        REQUIRE(op2.get_number_of_attempts() == niterations);
        REQUIRE(op2.get_number_of_attempts_for_correction() == (niterations - 100));

        REQUIRE(pop_size_summary.sample_size() == nsamples);

        REQUIRE((count_nheights_1 + count_nheights_2 + count_nheights_3 + count_nheights_4) == nsamples);

        double freq_nheights_1 = count_nheights_1 / (double)nsamples;
        double freq_nheights_2 = count_nheights_2 / (double)nsamples;
        double freq_nheights_3 = count_nheights_3 / (double)nsamples;
        double freq_nheights_4 = count_nheights_4 / (double)nsamples;

        double exp_freq = 1.0/336.0;
        double exp_count = nsamples/336.0;
        std::map< std::set< std::set<Split> >, double> bad_splits;

        double prop_error_threshold = 0.2;
        unsigned int total_trees_sampled = 0;
        std::map< std::set< std::set<Split> >, double> split_freqs;
        double chi_sq_test_statistic = 0.0;
        std::cout << "Total tree topologies sampled: " << split_counts.size() << "\n";
        for (auto s_c : split_counts) {
            total_trees_sampled += s_c.second;
            split_freqs[s_c.first] = s_c.second / (double)nsamples;
            std::cout << "Tree:\n";
            for (auto splitset : s_c.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            double prop_error = ((double)s_c.second - exp_count) / exp_count;
            std::cout << "  nsamples: " << s_c.second << "\n";
            std::cout << "  prop error: " << prop_error << "\n";
            if (fabs(prop_error) > prop_error_threshold) {
                bad_splits[s_c.first] = prop_error;
            }
            double count_diff = s_c.second - exp_count;
            chi_sq_test_statistic += (count_diff * count_diff) / exp_count;
        }

        double quantile_chi_sq_335_10 = 368.6;
        std::cout << "Chi-square test statistic: " << chi_sq_test_statistic << "\n";
        std::cout << "Chi-square(335) 0.9 quantile: " << quantile_chi_sq_335_10 << "\n";

        std::cout << "BAD SPLITS (proportional error > " << prop_error_threshold << ")\n";
        for (auto s_e : bad_splits) {
            std::cout << "\nTree:\n";
            for (auto splitset : s_e.first) {
                unsigned int s_count = 0;
                for (auto split : splitset) {
                    if (s_count > 0) {
                        // Indent shared splits
                        std::cout << "  ";
                    }
                    std::cout << "  " << split.as_string() << "\n";
                    ++s_count;
                }
            }
            std::cout << "  prop error: " << s_e.second << "\n";
        }

        write_r_script(split_counts, "../5-leaf-general-tree-test-constrained-pop-sizes.r");

        REQUIRE(total_trees_sampled == nsamples);

        // We should sample every possible tree
        REQUIRE(split_counts.size() == 336);

        double eps = 0.001;
        
        REQUIRE(pop_size_summary.mean() == Approx(pop_size_prior->get_mean()).epsilon(eps));
        REQUIRE(pop_size_summary.variance() == Approx(pop_size_prior->get_variance()).epsilon(eps));

        for (auto s_f : split_freqs) {
            REQUIRE(s_f.second == Approx(exp_freq).epsilon(eps));
        }

        REQUIRE(chi_sq_test_statistic < quantile_chi_sq_335_10);
    }
}
