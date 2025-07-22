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
#include "pll_likelihood.hpp"

namespace ecoevolity {

    template<class NodeType>
    class SeqTree : public BaseTree<NodeType> {
        protected:
            NucData::SharedPtr                          data_;
            std::shared_ptr< Likelihood<NodeType> >     likelihood_;
            Model::SharedPtr                            model_;
    
            // methods
    
            void init(const std::string & path);
            void init_data(const std::string & path);
            void init_model(const NucTreeAnalysisSettings& settings);
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
    
            // const NucData& get_data() const {
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
    
    
    template<class NodeType>
    SeqTree<NodeType>::SeqTree(
            const NucTreeAnalysisSettings& settings,
            RandomNumberGenerator& rng) {
        this->init(settings.data_settings.get_path());
        this->establish_tree_and_node_heights(settings.tree_model_settings, rng);
        try {
            this->check_if_leaf_labels_match_data();
        }
        catch (const EcoevolityError & e) {
            this->data_->replace_label_underscores_with_spaces();
            this->check_if_leaf_labels_match_data();
        }
    
        this->init_model(settings);
        this->init_likelihood();
    }
    
    template<class NodeType>
    void SeqTree<NodeType>::init(const std::string & path) {
        this->init_data(path);
        this->init_tree_from_leaf_labels(this->data_->get_seq_labels());
    }
    
    template<class NodeType>
    void SeqTree<NodeType>::init_data(const std::string & path) {
        this->data_ = std::make_shared<NucData>();
        this->data_->init_from_phylip_path(path);
    }
    
    template<class NodeType>
    void SeqTree<NodeType>::init_model(const NucTreeAnalysisSettings& settings) {
        this->model_ = std::make_shared<Model>();
    
        QMatrixNucleotide::SharedPtr q = std::make_shared<QMatrixNucleotide>();
        q->set_state_freqs(settings.state_freq_settings.get_values());
        q->set_exchangeabilities(settings.rate_matrix_settings.get_values());
    
        this->model_->init(q,
                settings.asrv_num_cats,
                settings.asrv_one_over_shape_settings.get_value());
    }
    
    template<class NodeType>
    void SeqTree<NodeType>::init_likelihood() {
        this->likelihood_ = std::make_shared< Likelihood<NodeType> >();
        this->likelihood_->init(
                this->data_,
                this->model_,
                this->root_);
    }

    template<class NodeType>
    void SeqTree<NodeType>::set_tree_to_upgma() {
        throw EcoevolityNotImplementedError(
                "Tree::set_tree_to_upgma is not implemented");
    }
    
    template<class NodeType>
    double SeqTree<NodeType>::get_ln_prob_of_drawing_node_state(
                    std::shared_ptr<NodeType> node) const {
        return node->get_ln_prob_of_drawing_state();
    }

    template<class NodeType>
    void SeqTree<NodeType>::check_if_leaf_labels_match_data() const {
        // Make sure labels on tree match the data
        std::vector<std::string> leaf_labels = this->root_->get_leaf_labels();
        std::vector<std::string> seq_labels = this->data_->get_seq_labels();
        if ((seq_labels.size() != leaf_labels.size()) ||
                (! std::is_permutation(seq_labels.begin(), seq_labels.end(), leaf_labels.begin()))) {
            std::ostringstream message;
            message << "Tip labels in provided tree do not match labels in data file";
            throw EcoevolityError(message.str());
        }
    }

    template<class NodeType>
    double SeqTree<NodeType>::compute_log_likelihood(
            const unsigned int nthreads) {
        if (this->ignoring_data()) {
            this->log_likelihood_.set_value(0.0);
            return 0.0;
        }
        double log_like = this->likelihood_->calc_log_likelihood();
        this->log_likelihood_.set_value(log_like);
        return log_like;
    }

}
