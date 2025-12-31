#pragma once

#include "pll.h"

#include "error.hpp"
#include "simple_nuc_data.hpp"
#include "model.hpp"
#include "node.hpp"

namespace ecoevolity {

    template<class NodeType>
    class Likelihood {
        // TODO:: Need to update to accommodate partitioned data and model
        public:
                                                    Likelihood();
                                                    ~Likelihood();

            void                                    clear();
            void                                    init(
                                                        NucData::SharedPtr data,
                                                        Model::SharedPtr model,
                                                        std::shared_ptr<NodeType> root_node);
            
            unsigned                                calc_num_edges_in_fully_resolved_tree() const;
            unsigned                                calc_num_internals_in_fully_resolved_tree() const;
            unsigned                                calc_num_nodes_in_fully_resolved_tree() const;
            unsigned                                get_num_seqs() const;

            double                                  calc_log_likelihood();

        private:

            void                                    update_state_freqs();
            void                                    update_exchangeabilities();
            void                                    update_asrv_rates();
            void                                    update_proportion_invariable_sites();
            void                                    set_tip_states();
            void                                    update_prob_matrices();
            void                                    update_operations();
            void                                    update_partials();

            std::vector<pll_operation_t>            _operations;
            pll_partition_t *                       _pll_partition;
            NucData::SharedPtr                      _data;
            Model::SharedPtr                        _model;
            std::shared_ptr<NodeType>               _root;
    }; 

    template<class NodeType>
    inline Likelihood<NodeType>::Likelihood() {
        //std::cout << "Constructing a Likelihood" << std::endl;
        this->clear();
    }

    template<class NodeType>
    inline Likelihood<NodeType>::~Likelihood() {
        //std::cout << "Destroying a Likelihood" << std::endl;
        if (_pll_partition) {
            pll_partition_destroy(_pll_partition);
        }
        _pll_partition = nullptr;
    }
    
    template<class NodeType>
    inline unsigned Likelihood<NodeType>::calc_num_edges_in_fully_resolved_tree() const {
        ECOEVOLITY_ASSERT(this->get_num_seqs() > 0);
        return ((2 * this->get_num_seqs()) - 2);
    }
    
    template<class NodeType>
    inline unsigned Likelihood<NodeType>::calc_num_internals_in_fully_resolved_tree() const {
        ECOEVOLITY_ASSERT(this->get_num_seqs() > 0);
        return (this->get_num_seqs() - 1);
    }

    template<class NodeType>
    inline unsigned Likelihood<NodeType>::calc_num_nodes_in_fully_resolved_tree() const {
        ECOEVOLITY_ASSERT(this->get_num_seqs() > 0);
        return ((2 * this->get_num_seqs()) - 1);
    }

    template<class NodeType>
    inline void Likelihood<NodeType>::clear() {
        this->_data                       = std::make_shared<NucData>();
        this->_model                      = std::make_shared<Model>();
        this->_root                       = std::make_shared<NodeType>();
        
        this->_operations.clear();
        if (_pll_partition) {
            pll_partition_destroy(_pll_partition);
        }
        _pll_partition = nullptr;
    }
    
    template<class NodeType>
    inline unsigned Likelihood<NodeType>::get_num_seqs() const {
        return this->_data->get_num_seqs();
    }

    template<class NodeType>
    inline void Likelihood<NodeType>::init(NucData::SharedPtr data, Model::SharedPtr model, std::shared_ptr<NodeType> root_node) {
        ECOEVOLITY_ASSERT(data->get_num_seqs() > 0);
        ECOEVOLITY_ASSERT(data->get_num_sites() > 0);
        this->_data = data;
        this->_model = model;
        this->_root = root_node;
        if (this->_data->get_num_seqs() != root_node->get_leaf_node_count()) {
            throw EcoevolityError("Likelihood<NodeType>::init: number of seqs and leaves must match");
        }

        unsigned arch;
        arch = PLL_ATTRIB_ARCH_CPU | PLL_ATTRIB_ARCH_SSE | PLL_ATTRIB_ARCH_AVX;
        _pll_partition = pll_partition_create(
                this->_data->get_num_seqs(),                        // Number of tip seqs
                this->calc_num_internals_in_fully_resolved_tree(),  // Extra number of partials (in addition to tips)
                this->_data->get_num_states(),                      // Number of states in data
                this->_data->get_num_sites(),                       // Length of alignment
                1,                                                  // Number of substitution models
                this->calc_num_edges_in_fully_resolved_tree(),      // Number of probability matrices to allocate
                this->_model->get_asrv()->get_num_categ(),          // Number of discrete rate categories
                this->calc_num_internals_in_fully_resolved_tree(),  // Number of scale buffers to allocate
                PLL_ATTRIB_ARCH_CPU);
                // arch);                                              // What hardware accel to us
        this->update_state_freqs();
        this->update_exchangeabilities();
        this->update_asrv_rates();
        this->set_tip_states();
        if (this->_model->get_asrv()->is_invariable_sites_model()) {
            // Update partition with info about which sites are invariant
            pll_update_invariant_sites(this->_pll_partition);
        }
    }

    template<class NodeType>
    inline void Likelihood<NodeType>::update_state_freqs() {
        pll_set_frequencies(
                _pll_partition,
                0,
                this->_model->get_q_matrix()->get_state_freqs());
    }

    template<class NodeType>
    inline void Likelihood<NodeType>::update_exchangeabilities() {
        pll_set_subst_params(
                _pll_partition,
                0,
                this->_model->get_q_matrix()->get_exchangeabilities());
    }

    template<class NodeType>
    inline void Likelihood<NodeType>::update_asrv_rates() {
        std::vector<double> asrv_rates(this->_model->get_asrv()->get_num_categ());
        pll_compute_gamma_cats(
                this->_model->get_asrv()->get_shape(),
                this->_model->get_asrv()->get_num_categ(),
                asrv_rates.data(),
                PLL_GAMMA_RATES_MEAN);
        pll_set_category_rates(_pll_partition, asrv_rates.data());
    }

    template<class NodeType>
    inline void Likelihood<NodeType>::update_proportion_invariable_sites() {
        pll_update_invariant_sites_proportion(
                _pll_partition,
                0,
                this->_model->get_asrv()->get_proportion_invariable_sites());
    }

    template<class NodeType>
    inline void Likelihood<NodeType>::set_tip_states() {
        for (auto n : this->_root->get_leaves()) {
            pll_set_tip_states(
                    _pll_partition,
                    n->get_index(),
                    pll_map_nt,
                    this->_data->get_seq(n->get_label()).c_str()
            );
        }
    }

    template<class NodeType>
    inline void Likelihood<NodeType>::update_prob_matrices() {
        unsigned max_n_edges = this->calc_num_edges_in_fully_resolved_tree();

        std::vector< std::shared_ptr<NodeType> > nodes;
        // nodes.reserve(max_n_edges + 1);
        this->_root->level_order(nodes);
        ECOEVOLITY_ASSERT(nodes.size() <= (max_n_edges + 1));

        std::vector<double> branch_lengths;
        branch_lengths.reserve(max_n_edges);
        std::vector<unsigned> matrix_indices;
        matrix_indices.reserve(max_n_edges);

        // subtract 1 for the root node, which does not have a trans prob
        // matrix
        unsigned zero_len_matrix_idx = nodes.size() - 1;
        unsigned n_children;
        std::shared_ptr<NodeType> n;
        for (int i = (nodes.size() - 1); i >= 1; --i) { // stop at 1 to skip root node
            n = nodes.at(i);
            // n_children = n->get_number_of_children();
            // for (unsigned j = 2; j < n_children; ++j) {
            //     branch_lengths.push_back(0.0);
            //     matrix_indices.push_back(zero_len_matrix_idx);
            // }
            // if (i == 0) {
            //     // This is the root node
            //     // We need to process it above in case it's a polytomy, but we
            //     // don't want to add it to the branch_lengths and
            //     // matrix_indices below
            //     assert(n->is_root());
            //     break;
            // }
            branch_lengths.push_back(n->get_length());
            matrix_indices.push_back(n->get_index());
        }
        // Only update the zero-length matrix once, because this matrix can be
        // used by all operations below that involve calculating the
        // probability of evolution along zero-length branches
        if (nodes.size() < (max_n_edges + 1)) {
            branch_lengths.push_back(0.0);
            matrix_indices.push_back(zero_len_matrix_idx);
        }
        ECOEVOLITY_ASSERT(branch_lengths.size() <= max_n_edges);
        ECOEVOLITY_ASSERT(matrix_indices.size() == branch_lengths.size());

        std::vector<unsigned> q_matrix_indices_of_asrv_cats(this->_model->get_asrv()->get_num_categ(), 0);

        pll_update_prob_matrices(
                _pll_partition,
                q_matrix_indices_of_asrv_cats.data(),
                matrix_indices.data(),
                branch_lengths.data(),
                branch_lengths.size());
    }

    template<class NodeType>
    inline void Likelihood<NodeType>::update_operations() {
        unsigned max_n_edges = this->calc_num_edges_in_fully_resolved_tree();
        unsigned max_n_internals = this->calc_num_internals_in_fully_resolved_tree();

        std::vector< std::shared_ptr<NodeType> > nodes;
        nodes.reserve(max_n_edges + 1);
        this->_root->level_order(nodes);
        ECOEVOLITY_ASSERT(nodes.size() <= (max_n_edges + 1));

        // Minus 1 for root node, which lacks a trans prob matrix
        unsigned zero_len_matrix_idx = nodes.size() - 1;

        this->_operations.clear();
        this->_operations.resize(max_n_internals);

        unsigned num_leaves = this->get_num_seqs();
        std::shared_ptr<NodeType> n;
        std::shared_ptr<NodeType> child1;
        std::shared_ptr<NodeType> child2;
        unsigned n_children;
        unsigned next_dummy_internal_idx = nodes.size();
        unsigned op_idx = 0;
        for (int i = (nodes.size() - 1); i >= 0; --i) {
            if (nodes.at(i)->is_leaf()) {
                continue;
            }
            n = nodes.at(i);

            unsigned node_idx = n->get_index();
            ECOEVOLITY_ASSERT(node_idx >= num_leaves);
            unsigned scaler_idx = node_idx - num_leaves;
            n_children = n->get_number_of_children();

            child1 = n->get_child(0);
            unsigned child1_idx = child1->get_index();
            unsigned child1_scaler_idx = PLL_SCALE_BUFFER_NONE;
            if (! child1->is_leaf()) {
                ECOEVOLITY_ASSERT(child1_idx >= num_leaves);
                child1_scaler_idx = child1_idx - num_leaves;
            }

            child2 = n->get_child(1);
            unsigned child2_idx = child2->get_index();
            unsigned child2_scaler_idx = PLL_SCALE_BUFFER_NONE;
            if (! child2->is_leaf()) {
                ECOEVOLITY_ASSERT(child2_idx >= num_leaves);
                child2_scaler_idx = child2_idx - num_leaves;
            }

            _operations.at(op_idx).parent_clv_index     = node_idx;
            _operations.at(op_idx).parent_scaler_index  = scaler_idx;
            _operations.at(op_idx).child1_clv_index     = child1_idx;
            _operations.at(op_idx).child2_clv_index     = child2_idx;
            _operations.at(op_idx).child1_matrix_index  = child1_idx;
            _operations.at(op_idx).child2_matrix_index  = child2_idx;
            _operations.at(op_idx).child1_scaler_index  = child1_scaler_idx;
            _operations.at(op_idx).child2_scaler_index  = child2_scaler_idx;
            ++op_idx;

            if (n_children > 2) {
                unsigned current_child_clade_idx = node_idx;
                std::shared_ptr<NodeType> next_child;
                for (unsigned j = 2; j < n_children; ++j) {
                    next_child = n->get_child(j);
                    unsigned next_child_idx = next_child->get_index();
                    unsigned next_child_scaler_idx = PLL_SCALE_BUFFER_NONE;
                    if (! next_child->is_leaf()) {
                        ECOEVOLITY_ASSERT(next_child_idx >= num_leaves);
                        next_child_scaler_idx = next_child_idx - num_leaves;
                    }
                    _operations.at(op_idx).parent_clv_index     = next_dummy_internal_idx;
                    _operations.at(op_idx).parent_scaler_index  = next_dummy_internal_idx - num_leaves;
                    _operations.at(op_idx).child1_clv_index     = current_child_clade_idx;
                    _operations.at(op_idx).child2_clv_index     = next_child_idx;
                    _operations.at(op_idx).child1_matrix_index  = zero_len_matrix_idx;
                    _operations.at(op_idx).child2_matrix_index  = next_child_idx;
                    _operations.at(op_idx).child1_scaler_index  = current_child_clade_idx - num_leaves;
                    _operations.at(op_idx).child2_scaler_index  = next_child_scaler_idx;
                    ++op_idx;

                    current_child_clade_idx = next_dummy_internal_idx;
                    ++next_dummy_internal_idx;
                }
            }
        }
        ECOEVOLITY_ASSERT(op_idx == max_n_internals);
        ECOEVOLITY_ASSERT(next_dummy_internal_idx == max_n_edges + 1);
    }

    template<class NodeType>
    inline void Likelihood<NodeType>::update_partials() {
        pll_update_partials(
                _pll_partition,
                _operations.data(),
                _operations.size());
    }

    template<class NodeType>
    inline double Likelihood<NodeType>::calc_log_likelihood() {
        this->update_state_freqs();
        this->update_exchangeabilities();
        this->update_asrv_rates();
        if (this->_model->get_asrv()->is_invariable_sites_model()) {
            this->update_proportion_invariable_sites();
        }
        this->update_prob_matrices();
        this->update_operations();
        this->update_partials();

        std::vector<unsigned> q_matrix_indices_of_asrv_cats(this->_model->get_asrv()->get_num_categ(), 0);
        unsigned n_nodes = this->calc_num_nodes_in_fully_resolved_tree();
        unsigned n_leaves = this->get_num_seqs();
        unsigned root_idx = n_nodes - 1;
        unsigned root_scaler_idx = root_idx - n_leaves;
        double logl = pll_compute_root_loglikelihood(
                _pll_partition,
                root_idx,                          // Index of CLV used to integrate likelihood
                root_scaler_idx,                   // Index of scaler for the node to integrate 
                q_matrix_indices_of_asrv_cats.data(), // Indices of q matrix to use for each ASRV category
                NULL);                             // Array to store per-site log-like values
        return logl;
    }
}
