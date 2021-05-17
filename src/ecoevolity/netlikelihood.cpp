/******************************************************************************
 * Copyright (C) 2015-2016 Jamie R. Oaks.
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

#include "netlikelihood.hpp"

namespace netlikelihood {

const static MatrixExponentiator matrix_exponentiator;

void compute_leaf_partials(
        std::shared_ptr<PopulationNetNode> node,
        const unsigned int red_allele_count,
        const unsigned int allele_count,
        const bool markers_are_dominant
        ) {
    if (allele_count == 0) {
        BiallelicPatternProbabilityMatrix m(0, 0);
        node->copy_bottom_pattern_probs(m);
        return;
    }
    if (markers_are_dominant) {
        unsigned int n = allele_count;
        unsigned int n_reds = red_allele_count;
        const unsigned int a_count = allele_count * 2;
        if (red_allele_count > 0) {
            BiallelicPatternProbabilityMatrix m(a_count);
            double p_r_k_n = 1.0;
            for (unsigned int r = 1; r <= n_reds; ++r) {
                p_r_k_n = (p_r_k_n * 2.0 * (n - r + 1.0)) / ((2.0 * n) - r + 1.0);
            }
            for (unsigned int k = n_reds; k <= (2 * n_reds); ++k) {
                if (k > n_reds) {
                    p_r_k_n = (p_r_k_n * ((2.0 * n_reds) - k + 1) * k) /
                              (2.0 * ( k - n_reds) * ((2.0 * n) - k + 1.0));
                }
                m.set_pattern_probability(a_count, k, p_r_k_n);
            }
            node->copy_bottom_pattern_probs(m);
            return;
        }
        BiallelicPatternProbabilityMatrix m(a_count, n_reds);
        node->copy_bottom_pattern_probs(m);
        return;
    }
    BiallelicPatternProbabilityMatrix m(allele_count, red_allele_count);
    node->copy_bottom_pattern_probs(m);
    return;
}

void compute_top_of_branch_partials(
        std::shared_ptr<PopulationNetNode> node,
        const unsigned int branch_index,
        const double u,
        const double v,
        const double mutation_rate, 
        const double ploidy
        ) {
    if (node->get_allele_count() == 0) {
        node->copy_top_pattern_probs(branch_index, node->get_bottom_pattern_probs(branch_index));
        ECOEVOLITY_DEBUG(
            std::cerr << "compute_top_of_branch_partials: Node " << node->get_label() << "has no alleles\n";
            std::cerr << "top " << branch_index << ":\n";
            std::cerr << node->get_top_pattern_probs(branch_index).to_string() << "\n";
        )
        return;
    }

    // Nested for loops here over discrete lengths and pop sizes (percentiles
    // of discretized distributions). Average the probs inside the
    // internal for loop, and sum averaged probs in outer for loop
    double theta = 2 * ploidy * node->get_population_size(branch_index) * mutation_rate;
    double length = node->get_length(branch_index) * mutation_rate;
    BiallelicPatternProbabilityMatrix m = matrix_exponentiator.expQTtx(
            node->get_allele_count(),
            u,
            v,
            theta,
            length,
            node->get_bottom_pattern_probs(branch_index));
    // Get average probs after outer for loop finishes, and copy
    // the final average probs to the top patterns
    m.set_pattern_probability(0, 0, node->get_bottom_pattern_probability(branch_index, 0, 0));
    node->copy_top_pattern_probs(branch_index, m);
    ECOEVOLITY_DEBUG(
        std::cerr << "compute_top_of_branch_partials: Node " << node->get_label() << "\n";
        std::cerr << "top " << branch_index << ":\n";
        std::cerr << node->get_top_pattern_probs(branch_index).to_string() << "\n";
    )
}

/**
 * @brief   Take conditional probabilities at the top of a daugther branch with
 *          two parents and propagate them to the conditional probabilites at
 *          the bottom of each parent branch.
 *
 * For each possible allele pattern at the top of the daughter branch, we need
 * to consider all the ways the alleles can be split between the parents and
 * calc the prob of each. For example, let's say the probability that an allele
 * came from the left / right parent is 4/12 / 8/12, respectively, and we have
 * the following conditional probs at the top of the daugher branch:
 *
 *   Daughter prob  Ways to split                       Prob (over 1728)
 *   0,0 = 1/12
 *                  0,0 | 0,0 = 1/12                    144
 *   1,0 = 1/12
 *                  1,0 | 0,0 = 1/12 * 4/12             48
 *                  0,0 | 1,0 = 1/12 * 8/12             96
 *   0,1 = 1/12
 *                  0,1 | 0,0 = 1/12 * 4/12             48
 *                  0,0 | 0,1 = 1/12 * 8/12             96
 *   2,0 = 2/12
 *                  2,0 | 0,0 = 2/12 * 4/12 * 4/12      32
 *                  0,0 | 2,0 = 2/12 * 8/12 * 8/12      128
 *                  1,0 | 1,0 = 2(2/12 * 4/12 * 8/12)   128
 *                              x 2 above because 2
 *                              ways for red alleles
 *                              to end up in each
 *                              parent
 *   1,1 = 3/12
 *                  1,1 | 0,0 = 3/12 * 4/12 * 4/12      48
 *                  0,0 | 1,1 = 3/12 * 8/12 * 8/12      192
 *                  1,0 | 0,1 = 3/12 * 4/12 * 8/12      96
 *                  0,1 | 1,0 = 3/12 * 4/12 * 8/12      96
 *   0,2 = 4/12
 *                  0,2 | 0,0 = 4/12 * 4/12 * 4/12      64
 *                  0,0 | 0,2 = 4/12 * 8/12 * 8/12      256
 *                  0,1 | 0,1 = 2(4/12 * 4/12 * 8/12)   256
 *                              x 2 above because 2
 *                              ways for green alleles
 *                              to end up in each
 *                              parent
 *                                                      Total = 1728/1728
 *
 * From above, we can calculate the conditional probability of all possible
 * allele patterns at the bottom of each parent branch (which is the goal of
 * this function)
 *
 *   Left parent bottom probs (over 1728)
 *   0,0 = 144+96+96+128+192+256 = 912
 *   1,0 = 48+128+96             = 272
 *   0,1 = 48+96+256             = 400
 *   2,0 = 32                    = 32
 *   1,1 = 48                    = 48
 *   0,2 = 64                    = 64
 *   total                       = 1728 / 1728
 *
 *   Right parent bottom probs (over 1728)
 *   0,0 = 144+48+48+32+48+64    = 384
 *   1,0 = 96+128+96             = 320
 *   0,1 = 96+96+256             = 448
 *   2,0 = 128                   = 128
 *   1,1 = 192                   = 192
 *   0,2 = 256                   = 256
 *   total                       = 1728 / 1728
 *
 * And if the probability that an allele came from the left / right parent is
 * 0/12 / 12/12, we have:
 *
 *   Daughter prob  Ways to split                       Prob (over 1728)
 *   0,0 = 1/12
 *                  0,0 | 0,0 = 1/12                    144
 *   1,0 = 1/12
 *                  1,0 | 0,0 = 1/12 * 0                0
 *                  0,0 | 1,0 = 1/12 * 1                144
 *   0,1 = 1/12
 *                  0,1 | 0,0 = 1/12 * 0                0
 *                  0,0 | 0,1 = 1/12 * 1                144
 *   2,0 = 2/12
 *                  2,0 | 0,0 = 2/12 * 0 * 0            0
 *                  0,0 | 2,0 = 2/12 * 1 * 1            288
 *                  1,0 | 1,0 = 2(2/12 * 0 * 1)         0
 *                              x 2 above because 2
 *                              ways for red alleles
 *                              to end up in each
 *                              parent
 *   1,1 = 3/12
 *                  1,1 | 0,0 = 3/12 * 0 * 0            0
 *                  0,0 | 1,1 = 3/12 * 1 * 1            432
 *                  1,0 | 0,1 = 3/12 * 0 * 1            0
 *                  0,1 | 1,0 = 3/12 * 0 * 1            0
 *   0,2 = 4/12
 *                  0,2 | 0,0 = 4/12 * 0 * 0            0
 *                  0,0 | 0,2 = 4/12 * 1 * 1            576
 *                  0,1 | 0,1 = 2(4/12 * 0 * 1)         0
 *                              x 2 above because 2
 *                              ways for green alleles
 *                              to end up in each
 *                              parent
 *                                                      Total = 1728/1728
 *
 * From above, we can calculate the conditional probability of all possible
 * allele patterns at the bottom of each parent branch (which is the goal of
 * this function)
 *
 *   Left parent bottom probs (over 1728)
 *   0,0 = 144+144+144+288+432+576 = 1728
 *   1,0 =                         = 0
 *   0,1 =                         = 0
 *   2,0 =                         = 0
 *   1,1 =                         = 0
 *   0,2 =                         = 0
 *   total                         = 1728 / 1728
 *
 *   Right parent bottom probs (over 1728)
 *   0,0 = 144                   = 144
 *   1,0 = 144                   = 144
 *   0,1 = 144                   = 144
 *   2,0 = 288                   = 288
 *   1,1 = 432                   = 432
 *   0,2 = 576                   = 576
 *   total                       = 1728 / 1728
 */
void split_top_of_branch_partials(
        const unsigned int max_num_alleles,
        const BiallelicPatternProbabilityMatrix & top_child_partials,
        const double prob_to_parent1,
        const double prob_to_parent2,
        BiallelicPatternProbabilityMatrix & bottom_partials_1,
        BiallelicPatternProbabilityMatrix & bottom_partials_2) {
    // TODO: Do we need to rescale probabilities by hypergeometric
    // probabilities before and after splitting like we do before and after
    // merging in merge_top_of_branch_partials?

    bottom_partials_1.reset(max_num_alleles);
    bottom_partials_2.reset(max_num_alleles);
    if (max_num_alleles < 1) {
        return;
    }

    unsigned int n_alleles, n_red, n_green, n_r_p1, n_r_p2, n_g_p1, n_g_p2;
    double p, nr_choose_nrp1, ng_choose_ngp1;
    for (n_alleles = 0; n_alleles <= max_num_alleles; ++n_alleles) {
        for (n_red = 0; n_red <= n_alleles; ++n_red) {
            n_green = n_alleles - n_red;
            for (n_r_p1 = 0; n_r_p1 <= n_red; ++n_r_p1) {
                n_r_p2 = n_red - n_r_p1;
                for (n_g_p1 = 0; n_g_p1 <= n_green; ++n_g_p1) {
                    n_g_p2 = n_green - n_g_p1;
                    p = (top_child_partials.get_pattern_probability(n_alleles, n_red)
                            * std::pow(prob_to_parent1, (n_r_p1 + n_g_p1))
                            * std::pow(prob_to_parent2, (n_r_p2 + n_g_p2)));
                    if ((n_r_p1 > 0) && (n_r_p2 > 0)) {
                        // If red alleles went to both parents we need to
                        // account for all the ways one parent can choose its
                        // red alleles from the child. We only have to worry
                        // about one parent (and it doesn't matter which one),
                        // because the other always gets what is left over.
                        // This is:
                        //   "number of red alleles"
                        //   choose
                        //   "number of reds chosen by parent 1"
                        nr_choose_nrp1 = std::exp(ln_n_choose_k(n_red, n_r_p1));
                        p *= nr_choose_nrp1;
                    }
                    if ((n_g_p1 > 0) && (n_g_p2 > 0)) {
                        // If green alleles went to both parents we need to
                        // account for the combinatorics like we did for the
                        // red alleles above.
                        ng_choose_ngp1 = std::exp(ln_n_choose_k(n_green, n_g_p1));
                        p *= ng_choose_ngp1;
                    }
                    // std::cout << n_g_p1 << "," << n_r_p1 << " | " << n_g_p2 << "," << n_r_p2 << " " << p << "\n";
                    bottom_partials_1.set_pattern_probability(
                            (n_r_p1 + n_g_p1), n_r_p1,
                            (bottom_partials_1.get_pattern_probability(
                                    (n_r_p1 + n_g_p1), n_r_p1) + p));
                    bottom_partials_2.set_pattern_probability(
                            (n_r_p2 + n_g_p2), n_r_p2,
                            (bottom_partials_2.get_pattern_probability(
                                    (n_r_p2 + n_g_p2), n_r_p2) + p));
                }
            }
        }
    }
}

void merge_top_of_branch_partials(
        const unsigned int allele_count_child1,
        const unsigned int allele_count_child2,
        const double prob_no_alleles_child1,
        const double prob_no_alleles_child2,
        std::vector<double> & top_partials_child1,
        std::vector<double> & top_partials_child2,
        unsigned int & merged_allele_count,
        std::vector<double> & merged_pattern_probs,
        double & merged_prob_no_alleles,
        const bool do_hypergeom_scaling) {

    if (do_hypergeom_scaling) {
        for (unsigned int n = 1; n <= allele_count_child1; ++n) {
            double b_nr = 1.0;
            for (unsigned int r = 0; r <= n; ++r) {
                top_partials_child1.at(((n*(n+1))/2)-1+r) *= b_nr;
                b_nr *= ((double)n - r)/(r+1);
            }
        }

        for (unsigned int n = 1; n <= allele_count_child2; ++n) {
            double b_nr = 1.0;
            for (unsigned int r = 0; r <= n; ++r) {
                top_partials_child2.at(((n*(n+1))/2)-1+r) *= b_nr;
                b_nr *= ((double)n - r)/(r+1);
            }
        }
    }

    unsigned int allele_count = allele_count_child1 + allele_count_child2;
    merged_pattern_probs.assign(
            ((((allele_count + 1) * (allele_count + 2))/2) - 1),
            0);
    merged_prob_no_alleles = 0.0;
    unsigned int n1, r1, n2, r2;
    double top_prob_child1, top_prob_child2;
    for (n1 = 0; n1 <= allele_count_child1; ++n1) {
        for (r1 = 0; r1 <= n1; ++r1) {
            if (n1 == 0) {
                top_prob_child1 = prob_no_alleles_child1;
            } else {
                top_prob_child1 = top_partials_child1.at(n1*(n1+1)/2-1+r1);
            }
            for (n2 = 0; n2 <= allele_count_child2; ++n2) {
                for (r2 = 0; r2 <= n2; ++r2) {
                    if (n2 == 0) {
                        top_prob_child2 = prob_no_alleles_child2;
                    } else {
                        top_prob_child2 = top_partials_child2.at(n2*(n2+1)/2-1+r2);
                    }
                    if ((n1 + n2) < 1) {
                        merged_prob_no_alleles += top_prob_child1 * top_prob_child2;
                    }
                    else {
                        merged_pattern_probs.at((n1+n2)*(n1+n2+1)/2-1+(r1+r2)) += top_prob_child1 * top_prob_child2;
                    }
                }
            }
        }
    }

    if (do_hypergeom_scaling) {
        for (unsigned int n = 1; n <= allele_count; ++n) {
            double b_nr = 1.0;
            for (unsigned int r = 0; r <= n; ++r) {
                double f_nr = merged_pattern_probs.at(n*(n+1)/2-1+r);
                f_nr /= b_nr;
                // TODO: better way to fix this?
                f_nr = std::max(f_nr, 0.0);
                merged_pattern_probs.at(n*(n+1)/2-1+r) = f_nr;
                b_nr *= ((double)n - r)/(r+1);

            }
        }
    }
    merged_allele_count = allele_count;
}

void compute_internal_partials(
        std::shared_ptr<PopulationNetNode> node,
        std::set< std::shared_ptr<const PopulationNetNode> > & retic_nodes_visited) {

    if (retic_nodes_visited.count(node) > 0) {
        ECOEVOLITY_DEBUG(
            std::cerr << "compute_internal_partials: Retic node " << node->get_label() << " already visited\n";
        )
        return;
    }
    node->visit(retic_nodes_visited);

    unsigned int number_of_children_with_alleles = 0;
    std::vector<unsigned int> indices_of_children_with_alleles;
    for (unsigned int child_idx = 0; child_idx < node->get_number_of_children(); ++child_idx) {
        if (node->get_child(child_idx)->get_allele_count() > 0) {
            ++number_of_children_with_alleles;
            indices_of_children_with_alleles.push_back(child_idx);
        }
    }
    if (number_of_children_with_alleles < 1) {
        // std::ostringstream message;
        // message << "compute_internal_partials; "
        //         << "no children have alleles!";
        // throw EcoevolityError(message.str());
        BiallelicPatternProbabilityMatrix m(0, 0);
        ECOEVOLITY_DEBUG(
            std::cerr << "compute_internal_partials: Node " << node->get_label() << " has not children with alleles:\n";
        )
        for (unsigned int i = 0; i < node->get_number_of_parents(); ++i) {
            node->copy_bottom_pattern_probs(i, m);
            ECOEVOLITY_DEBUG(
                std::cerr << "bottom " << i << ":\n";
                std::cerr << node->get_bottom_pattern_probs(i).to_string() << "\n";
            )
        }
        return;
    }

    /* std::cout << "Processing " << node->get_label() << "\n"; */
    if ((node->get_number_of_children() == 1) && (node->has_multiple_parents())) {
        // Handle reticulation
        ECOEVOLITY_ASSERT(node->get_number_of_parents() == 2);
        BiallelicPatternProbabilityMatrix bottom_partials_1;
        BiallelicPatternProbabilityMatrix bottom_partials_2;
        split_top_of_branch_partials(
                node->get_child(0)->get_allele_count(),
                node->get_child(0)->get_top_pattern_probs(node),
                node->get_inheritance_proportion(0),
                node->get_inheritance_proportion(1),
                bottom_partials_1,
                bottom_partials_2);
        node->copy_bottom_pattern_probs(0, bottom_partials_1);
        node->copy_bottom_pattern_probs(1, bottom_partials_2);
        /* retic_nodes_visited.insert(node); */
        ECOEVOLITY_DEBUG(
            std::cerr << "compute_internal_partials: Splitting retic node " << node->get_label() << "\n";
            std::cerr << "bottom 0:\n";
            std::cerr << node->get_bottom_pattern_probs(0).to_string() << "\n";
            std::cerr << "bottom 1:\n";
            std::cerr << node->get_bottom_pattern_probs(1).to_string() << "\n";
        )
        return;
    }

    if (number_of_children_with_alleles == 1) {
        node->copy_bottom_pattern_probs(node->get_child(indices_of_children_with_alleles.at(0))->get_top_pattern_probs(node));
        ECOEVOLITY_DEBUG(
            std::cerr << "compute_internal_partials: Handling in-1-out-1 node " << node->get_label() << "\n";
            std::cerr << "bottom:\n";
            std::cerr << node->get_bottom_pattern_probs().to_string() << "\n";
        )
        return;
    }

    unsigned int merged_allele_count = node->get_child(indices_of_children_with_alleles.at(0))->get_allele_count();

    std::vector<double> merged_pattern_probs = node->get_child(
            indices_of_children_with_alleles.at(0))->get_top_pattern_probs(node).get_pattern_prob_matrix();
    double merged_prob_no_alleles = node->get_child(
            indices_of_children_with_alleles.at(0))->get_top_pattern_probability(
                    node, 0, 0);

    for (unsigned int i = 1; i < node->get_number_of_children(); ++i) {
        std::vector<double> pattern_probs_child1 = merged_pattern_probs;
        double prob_no_alleles_child1 = merged_prob_no_alleles;
        unsigned int allele_count_child1 = merged_allele_count;
        unsigned int allele_count_child2 = node->get_child(indices_of_children_with_alleles.at(i))->get_allele_count();
        double prob_no_alleles_child2 = node->get_child(indices_of_children_with_alleles.at(i))->get_top_pattern_probability(node,0,0);
        std::vector<double> pattern_probs_child2 = node->get_child(
                indices_of_children_with_alleles.at(i))->get_top_pattern_probs(node
                ).get_pattern_prob_matrix();
        merge_top_of_branch_partials(
                allele_count_child1,
                allele_count_child2,
                prob_no_alleles_child1,
                prob_no_alleles_child2,
                pattern_probs_child1,
                pattern_probs_child2,
                merged_allele_count,
                merged_pattern_probs,
                merged_prob_no_alleles,
                true);
    }
    BiallelicPatternProbabilityMatrix m(merged_allele_count, merged_pattern_probs);
    m.set_pattern_probability(0, 0, merged_prob_no_alleles);
    node->copy_bottom_pattern_probs(m);
    ECOEVOLITY_DEBUG(
        std::cerr << "compute_internal_partials: Merging node " << node->get_label() << "\n";
        std::cerr << "bottom:\n";
        std::cerr << node->get_bottom_pattern_probs().to_string() << "\n";
    )
}

void compute_pattern_partials(
        std::shared_ptr<PopulationNetNode> node,
        const std::vector<unsigned int>& red_allele_counts,
        const std::vector<unsigned int>& allele_counts,
        const double u,
        const double v,
        const double mutation_rate,
        const double ploidy,
        const bool markers_are_dominant,
        std::set< std::shared_ptr<const PopulationNetNode> > & retic_nodes_visited
        ) {
    ECOEVOLITY_DEBUG(
        std::cerr << "*** " << node->get_label() << " ***\n";
    )
    ECOEVOLITY_ASSERT(red_allele_counts.size() == allele_counts.size());
    if (node->is_leaf()) {
        compute_leaf_partials(
                node,
                red_allele_counts.at(node->get_index()),
                allele_counts.at(node->get_index()),
                markers_are_dominant);
        ECOEVOLITY_DEBUG(
            std::cerr << "Bottom of leaf " << node->get_label() << ":\n";
            std::cerr << node->get_bottom_pattern_probs().to_string() << "\n";
        )
        return;
    }
    else if (node->get_number_of_children() == 1) {
        if (retic_nodes_visited.count(node) > 0) {
            ECOEVOLITY_DEBUG(
                std::cerr << "compute_pattern_partials: already visited retic node " << node->get_label() << "... Skipping!\n";
            )
            return;
        }
        compute_pattern_partials(node->get_child(0),
                red_allele_counts,
                allele_counts,
                u,
                v,
                mutation_rate,
                ploidy,
                markers_are_dominant,
                retic_nodes_visited);
        compute_top_of_branch_partials(
                node->get_child(0),
                node->get_child(0)->get_parent_index(node),
                u, v, mutation_rate, ploidy);
        compute_internal_partials(node, retic_nodes_visited);
        return;
    }
    else if (node->get_number_of_children() > 1) {
        for (unsigned int child_idx = 0; child_idx < node->get_number_of_children(); ++child_idx) {
            compute_pattern_partials(node->get_child(child_idx),
                    red_allele_counts,
                    allele_counts,
                    u,
                    v,
                    mutation_rate,
                    ploidy,
                    markers_are_dominant,
                    retic_nodes_visited);
            compute_top_of_branch_partials(
                    node->get_child(child_idx),
                    node->get_child(child_idx)->get_parent_index(node),
                    u, v, mutation_rate, ploidy);
        }
        compute_internal_partials(node, retic_nodes_visited);
        return;
    }
    std::ostringstream message;
    message << "compute_pattern_partials(); "
            << "unexpected number of children: "
            << node->get_number_of_children();
    throw EcoevolityError(message.str());
}

std::vector< std::vector<double> > compute_root_probabilities(
        std::shared_ptr<PopulationNetNode> root,
        const double u,
        const double v,
        const double mutation_rate,
        const double ploidy
        ) {
    unsigned int N = root->get_allele_count();
    std::vector< std::vector<double> > x (N + 1); 
    QMatrix q = QMatrix(
            N,
            u,
            v,
            2 * ploidy * root->get_population_size() * mutation_rate);
    std::vector<double> xcol = q.find_orthogonal_vector();

    // ECOEVOLITY_DEBUG(
    //     std::cerr << "xcol = [";
    //     for (unsigned int i = 0; i < xcol.size(); ++i) {
    //         std::cerr << xcol.at(i) << " ";
    //     }
    //     std::cerr << "]" << std::endl;
    // )

    unsigned int index = 1;
    for (unsigned int n = 1; n <= N; ++n) {
        x.at(n).resize(n + 1, 0.0);
        double row_sum = 0.0;
        for (unsigned int r = 0; r <= n; ++r) {
            double xcol_index = std::max(xcol.at(index), 0.0);
            row_sum += xcol_index;
            x.at(n).at(r) = xcol_index;
            index++;
        }
        for (unsigned int r = 0; r <= n; ++r) {
            x.at(n).at(r) = x.at(n).at(r) / row_sum;
        }
    }
    return x;
}

double compute_root_likelihood(
        std::shared_ptr<PopulationNetNode> root,
        const double u,
        const double v,
        const double mutation_rate,
        const double ploidy
        ) {
    unsigned int N = root->get_allele_count();
    std::vector< std::vector<double> > conditionals = compute_root_probabilities(root, u, v, mutation_rate, ploidy);

    // ECOEVOLITY_DEBUG(
    //     for (unsigned int n = 1; n <= N; ++n) {
    //         for (unsigned int r = 0; r <= n; ++r) {
    //             std::cerr << "root height: " << root->get_height() << std::endl;
    //             std::cerr << "conditional[" << n << ", " << r << "] = " << conditionals.at(n).at(r) << std::endl;
    //             std::cerr << "bottom_pattern_probs[" << n << ", " << r << "] = " << root->get_bottom_pattern_probability(n, r) << std::endl;
    //         }
    //     }
    // )

    // NOTE about analytically integrating over pop sizes using a discretized
    // distibution: To integrate over pop size of the root branch, before this
    // point, we need to sum over `condtionals`. The bottom pattern probs
    // (`root.get_bottom_pattern_probability(n, r)`) should already be
    // integrated over pop sizes up to the bottom of the root branch, so only
    // the `conditionals` need to be integrated here.
    double sum = 0.0;
    for (unsigned int n = 1; n <= N; ++n) {
        for (unsigned int r = 0; r <= n; ++r) {
            double term = conditionals.at(n).at(r) * root->get_bottom_pattern_probability(n, r);
            sum += term;
            // if (sum < 0.0) {
            //     std::ostringstream message;
            //     message << "compute_root_likelihood(): Numerical error; "
            //             << "root likelihood is negative: "
            //             << sum;
            //     throw EcoevolityError(message.str());
            // }
        }
    }

    // ECOEVOLITY_DEBUG(
    //     std::cerr << "root likelihood: " << sum << std::endl;
    // )
    if (sum < 0.0) {
        return 0.0;
    }
    // TODO: There's probably a better way to deal with NANs before this point
    // This is likely the result of underflow when probabilities are getting
    // extremely small.
    // When conditionals are on the order of 1e-305, the bottom pattern probs
    // can be NAN or -NAN.
    if (std::isnan(sum)) {
        return 0.0;
    }
    return sum;
}

double compute_pattern_likelihood(
        std::shared_ptr<PopulationNetNode> root,
        const std::vector<unsigned int>& red_allele_counts,
        const std::vector<unsigned int>& allele_counts,
        const double u,
        const double v,
        const double mutation_rate,
        const double ploidy,
        const bool markers_are_dominant
        ) {
    std::set< std::shared_ptr<const PopulationNetNode> > retic_nodes_visited;
    compute_pattern_partials(root,
            red_allele_counts,
            allele_counts,
            u,
            v,
            mutation_rate,
            ploidy,
            markers_are_dominant,
            retic_nodes_visited);
    return compute_root_likelihood(root, u, v, mutation_rate, ploidy);
}

void compute_constant_pattern_log_likelihood_correction(
        std::shared_ptr<PopulationNetNode> root,
        const std::vector< std::vector<unsigned int> > & unique_allele_counts,
        const std::vector<unsigned int> & unique_allele_count_weights,
        const double u,
        const double v,
        const double mutation_rate,
        const double ploidy,
        const bool markers_are_dominant,
        const bool state_frequencies_are_constrained,
        double& constant_log_likelihood_correction 
        ) {
    double lnl_correction = 0.0;
    for (unsigned int pattern_idx = 0;
            pattern_idx < unique_allele_count_weights.size();
            ++pattern_idx) {
        std::vector<unsigned int> red_allele_counts (
                unique_allele_counts.at(pattern_idx).size(),
                0); 
        double all_green_likelihood = compute_pattern_likelihood(root,
                red_allele_counts,
                unique_allele_counts.at(pattern_idx),
                u,
                v,
                mutation_rate,
                ploidy,
                markers_are_dominant);
        double all_red_likelihood = all_green_likelihood;
        if (! state_frequencies_are_constrained) {
            all_red_likelihood = compute_pattern_likelihood(root,
                    unique_allele_counts.at(pattern_idx),
                    unique_allele_counts.at(pattern_idx),
                    u,
                    v,
                    mutation_rate,
                    ploidy,
                    markers_are_dominant);
        }
        double variable_likelihood = 1.0 - all_green_likelihood - all_red_likelihood;
        if (variable_likelihood <= 0.0) {
            lnl_correction = -std::numeric_limits<double>::infinity();
            break;
        }
        lnl_correction += (unique_allele_count_weights.at(pattern_idx) *
                std::log(variable_likelihood));
    }
    constant_log_likelihood_correction = lnl_correction;
}

double get_log_likelihood_for_pattern_range(
        std::shared_ptr<PopulationNetNode> root,
        const std::vector< std::vector<unsigned int> >& red_allele_count_matrix,
        const std::vector< std::vector<unsigned int> >& allele_count_matrix,
        const std::vector<unsigned int>& pattern_weights,
        const unsigned int start_index,
        const unsigned int stop_index,
        const double u,
        const double v,
        const double mutation_rate,
        const double ploidy,
        const bool markers_are_dominant
        ) {
    ECOEVOLITY_ASSERT((red_allele_count_matrix.size() == allele_count_matrix.size()) &&
                      (red_allele_count_matrix.size() == pattern_weights.size()));
    double log_likelihood = 0.0;
    for (unsigned int pattern_idx = start_index;
            pattern_idx < stop_index;
            ++pattern_idx) {
        double pattern_likelihood = compute_pattern_likelihood(root,
                red_allele_count_matrix.at(pattern_idx),
                allele_count_matrix.at(pattern_idx),
                u,
                v,
                mutation_rate,
                ploidy,
                markers_are_dominant);
        if (pattern_likelihood <= 0.0) {
            return -std::numeric_limits<double>::infinity();
        }
        double weight = (double) pattern_weights.at(pattern_idx);
        log_likelihood += weight * std::log(pattern_likelihood);
    }
    return log_likelihood;
}

double get_log_likelihood(
        std::shared_ptr<PopulationNetNode> root,
        const std::vector< std::vector<unsigned int> >& red_allele_count_matrix,
        const std::vector< std::vector<unsigned int> >& allele_count_matrix,
        const std::vector<unsigned int>& pattern_weights,
        const std::vector< std::vector<unsigned int> > & unique_allele_counts,
        const std::vector<unsigned int> & unique_allele_count_weights,
        const double u,
        const double v,
        const double mutation_rate,
        const double ploidy,
        const bool markers_are_dominant,
        const bool state_frequencies_are_constrained,
        const bool constant_sites_removed,
        double& constant_log_likelihood_correction,
        unsigned int nthreads
        ) {
#ifdef BUILD_WITH_THREADS
    if (nthreads < 2) {
#endif
        if (constant_sites_removed) {
            compute_constant_pattern_log_likelihood_correction(
                    root,
                    unique_allele_counts,
                    unique_allele_count_weights,
                    u,
                    v,
                    mutation_rate,
                    ploidy,
                    markers_are_dominant,
                    state_frequencies_are_constrained,
                    constant_log_likelihood_correction);
        }
        return get_log_likelihood_for_pattern_range(
                root,
                red_allele_count_matrix,
                allele_count_matrix,
                pattern_weights,
                0,
                pattern_weights.size(),
                u,
                v,
                mutation_rate,
                ploidy,
                markers_are_dominant);
#ifdef BUILD_WITH_THREADS
    }
    double log_likelihood = 0.0;
    const unsigned int npatterns = pattern_weights.size();
    if (npatterns < nthreads) {
        nthreads = npatterns;
    }
    const unsigned int batch_size = npatterns / nthreads;
    unsigned int start_idx = 0;
    std::vector< std::future<double> > threads(nthreads - 1);
    std::vector< std::shared_ptr<PopulationNetNode> > root_clones(nthreads - 1);

    // Launch nthreads - 1 threads
    for (unsigned int i = 0; i < (nthreads - 1); ++i) {
        root_clones.at(i) = root->get_clade_clone();
        threads.at(i) = std::async(
                std::launch::async,
                get_log_likelihood_for_pattern_range,
                root_clones.at(i),
                std::cref(red_allele_count_matrix),
                std::cref(allele_count_matrix),
                std::cref(pattern_weights),
                start_idx,
                start_idx + batch_size,
                u,
                v,
                mutation_rate,
                ploidy,
                markers_are_dominant
                );
        start_idx += batch_size;
    }

    // Use the main thread as the last thread
    if (constant_sites_removed) {
        compute_constant_pattern_log_likelihood_correction(
                root,
                unique_allele_counts,
                unique_allele_count_weights,
                u,
                v,
                mutation_rate,
                ploidy,
                markers_are_dominant,
                state_frequencies_are_constrained,
                constant_log_likelihood_correction);
    }
    log_likelihood += get_log_likelihood_for_pattern_range(
            root,
            red_allele_count_matrix,
            allele_count_matrix,
            pattern_weights,
            start_idx,
            pattern_weights.size(),
            u,
            v,
            mutation_rate,
            ploidy,
            markers_are_dominant);

    // Join the launched threads
    for (auto &t : threads) {
        log_likelihood += t.get();
    }
    return log_likelihood;
#endif
}

}
