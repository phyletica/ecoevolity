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

#ifndef ECOEVOLITY_TREE_HPP
#define ECOEVOLITY_TREE_HPP

#include <algorithm>

#include "data.hpp"
#include "node.hpp"
#include "matrix.hpp"
#include "parameter.hpp"
#include "error.hpp"
#include "assert.hpp"

class PopulationTree {
    protected:
        BiallelicData data_;
        PopulationNode * root_;
        MatrixExponentiator matrix_exponentiator;
        PositiveRealParameter * u_ = new PositiveRealParameter(1.0);
        PositiveRealParameter * v_ = new PositiveRealParameter(1.0);
        std::vector<double> pattern_likelihoods_;
        LogProbability log_likelihood_ = LogProbability(0.0);
        LogProbability log_likelihood_correction_ = LogProbability(0.0);
        bool likelihood_correction_was_calculated_ = false;
        Probability all_green_pattern_likelihood_ = Probability(0.0);
        Probability all_red_pattern_likelihood_ = Probability(0.0);
        bool correct_for_full_likelihood_ = true;
        bool correct_for_constant_patterns_ = true;
        int number_of_constant_red_sites_ = -1;
        int number_of_constant_green_sites_ = -1;
        bool use_removed_constant_site_counts_ = false;

        // methods
        void init_tree();
        bool constant_site_counts_were_provided();
        void calculate_likelihood_correction();

        double calculate_log_binomial(
                unsigned int red_allele_count,
                unsigned int allele_count) const;

        double compute_pattern_likelihood(int pattern_index);
        void compute_pattern_likelihoods();

        std::vector< std::vector<double> > compute_root_probabilities();
        double compute_root_likelihood();

        void compute_pattern_partials(
                int pattern_index,
                PopulationNode * node);
        void compute_internal_partials(PopulationNode * node);
        void compute_top_of_branch_partials(PopulationNode * node);
        void compute_leaf_partials(
                int pattern_index,
                PopulationNode * node);
        

    public:
        PopulationTree() { }
        PopulationTree(
                const std::string path, 
                const char population_name_delimiter = '_',
                const bool population_name_is_prefix = true,
                const bool genotypes_are_diploid = true,
                const bool markers_are_dominant = false,
                const bool validate = true);
        ~PopulationTree () { }

        void init(
                const std::string path, 
                const char population_name_delimiter = '_',
                const bool population_name_is_prefix = true,
                const bool genotypes_are_diploid = true,
                const bool markers_are_dominant = false,
                const bool validate = true);

        void fold_patterns();

        PopulationNode * get_root() const {return this->root_;}
        void set_root_height(double height);
        void update_root_height(double height);
        const double& get_root_height() const;
        void store_root_height();
        void restore_root_height();
        void set_root_height_parameter(PositiveRealParameter * h);
        PositiveRealParameter * get_root_height_parameter() const;

        void set_u(double u);
        void set_v(double v);
        void update_u(double u);
        void update_v(double v);
        const double& get_u() const;
        const double& get_v() const;
        void store_u();
        void store_v();
        void restore_u();
        void restore_v();
        void set_u_parameter(PositiveRealParameter * u);
        void set_v_parameter(PositiveRealParameter * v);
        PositiveRealParameter * get_u_parameter() const;
        PositiveRealParameter * get_v_parameter() const;

        void set_root_coalescence_rate(double rate);
        void set_coalescence_rate(double rate);

        double get_likelihood_correction(bool force = false);

        double compute_log_likelihood();

        void store_state();
        void store_likelihood();
        void store_parameters();
        void store_all_coalescence_rates();
        virtual void store_all_heights();
        void restore_state();
        void restore_likelihood();
        void restore_parameters();
        void restore_all_coalescence_rates();
        virtual void restore_all_heights();
};

class ComparisonPopulationTree: public PopulationTree {
    public:
        ComparisonPopulationTree() { }
        ComparisonPopulationTree(
                const std::string path, 
                const char population_name_delimiter = '_',
                const bool population_name_is_prefix = true,
                const bool genotypes_are_diploid = true,
                const bool markers_are_dominant = false,
                const bool validate = true);

        void set_child_coalescence_rate(unsigned int child_index, double rate);
        void update_child_coalescence_rate(unsigned int child_index, double rate);
        const double& get_child_coalescence_rate(unsigned int child_index);
        void store_child_coalescence_rate(unsigned int child_index);
        void restore_child_coalescence_rate(unsigned int child_index);
        void set_child_coalescence_rate_parameter(unsigned int child_index,
                PositiveRealParameter * r);
        PositiveRealParameter * get_child_coalescence_rate_parameter(
                unsigned int child_index) const;

        void set_height(double height) {this->set_root_height(height);}
        void update_height(double height) {this->update_root_height(height);}
        const double& get_height() const {return this->get_root_height();}
        void store_height() {this->store_root_height();}
        void restore_height() {this->restore_root_height();}
        void set_height_parameter(PositiveRealParameter * h) {
            this->set_root_height_parameter(h);
        }
        PositiveRealParameter * get_height_parameter() const {
            return this->get_root_height_parameter();
        }

        void store_all_heights() {
            this->store_height();
        }
        void restore_all_heights() {
            this->restore_height();
        }
};

#endif
