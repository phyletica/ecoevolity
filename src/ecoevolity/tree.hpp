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

#pragma once

#include "basetree.hpp"
#include "simple_nuc_data.hpp"
#include "node.hpp"
#include "parameter.hpp"
#include "probability.hpp"
#include "error.hpp"
#include "debug.hpp"
#include "assert.hpp"
#include "general_tree_settings.hpp"


template<class NodeType>
class SeqTree : public BaseTree<NodeType> {
    protected:
        ecoevolity::NucData::SharedPtr          data_;
        ecoevolity::Likelihood::SharedPtr       likelihood_;
        Model::SharedPtr                        model_;

        // methods

        void init();
        void init_data();
        void init_model();
        void init_likelihood();

        void set_tree_to_upgma();

    public:
        SeqTree() { }
        SeqTree(
                const NucTreeAnalysisSettings& settings,
                RandomNumberGenerator& rng
                );

        SeqTree(
                std::shared_ptr<NodeType> root
            ) : BaseTree<NodeType>(root) { }

        SeqTree(
                std::istream & tree_stream,
                const std::string & ncl_file_format
            ) : BaseTree<NodeType>(tree_stream, ncl_file_format) { }

        SeqTree(
                const NxsFullTreeDescription & tree_description,
                NxsTaxaBlock * taxa_block,
                const double ultrametricity_tolerance,
                const double multiplier = -1.0
            ) : BaseTree<NodeType>(tree_description, taxa_block,
                        ultrametricity_tolerance, multiplier) { }

        double get_ln_prob_of_drawing_node_state(
                std::shared_ptr<NodeType> node) const;

        void check_if_leaf_labels_match_data() const;

        bool initialized() const {return (bool)this->root_;}

        // const ecoevolity::NucData& get_data() const {
        //     return this->data_;
        // }

        double compute_log_likelihood(const unsigned int nthreads = 1);

        // double get_derived_class_component_of_log_prior_density() const;

        // double compute_log_prior_density_of_state_frequencies() const;

        // void store_derived_class_parameters();
        // void restore_derived_class_parameters();

        // void draw_from_prior(RandomNumberGenerator& rng);

        // void write_state_log_header(std::ostream& out,
        //         const std::string& delimiter = "\t",
        //         const bool short_summary = false) const;
        // void log_state(std::ostream& out,
        //         const unsigned int generation_index,
        //         const std::string& delimiter = "\t",
        //         const bool short_summary = false) const;
};
