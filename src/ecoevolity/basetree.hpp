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

#ifndef ECOEVOLITY_BASETREE_HPP
#define ECOEVOLITY_BASETREE_HPP

#include <ncl/nxsmultiformat.h>

#include "split.hpp"
#include "parameter.hpp"
#include "probability.hpp"
#include "error.hpp"
#include "assert.hpp"

class PopSizeScaler;
class GlobalHeightSizeMixer;
template<class TreeType> class GlobalNodeHeightDirichletOperator;

template<class NodeType>
class BaseTree {
        friend class PopSizeScaler;
        friend class GlobalHeightSizeMixer;
        template<class TreeType>
        friend class GlobalNodeHeightDirichletOperator;

    protected:
        std::shared_ptr<NodeType> root_;
        std::shared_ptr<NodeType> stored_root_;
        std::vector< std::shared_ptr<PositiveRealParameter> > node_heights_;
        LogProbabilityDensity log_likelihood_ = LogProbabilityDensity(0.0);
        LogProbabilityDensity log_prior_density_ = LogProbabilityDensity(0.0);
        bool is_dirty_ = true;

        // Alpha and beta parameters of the beta-distributed prior on
        // (non-root) internal node heights
        std::shared_ptr<PositiveRealParameter> alpha_of_node_height_beta_prior_ = std::make_shared<PositiveRealParameter>(
                std::make_shared<GammaDistribution>(10.0, 0.1), // prior
                1.0,   // value of parameter
                true); // is fixed?
        std::shared_ptr<PositiveRealParameter> beta_of_node_height_beta_prior_ = std::make_shared<PositiveRealParameter>(
                std::make_shared<GammaDistribution>(10.0, 0.1), // prior
                1.0,   // value of parameter
                true); // is fixed?

        bool ignore_data_ = false;
        unsigned int number_of_likelihood_calculations_ = 0;

        std::vector< std::shared_ptr<NodeType> > pre_ordered_nodes_;
        std::vector< std::shared_ptr<NodeType> > level_ordered_nodes_;

        // The case where a singleton polytomy is mapped to the height we are
        // splitting; we have to handle this a little differently.
        // Because there is only one node mapped to the height, we have to make
        // sure we break the polytomy into 2 sets of nodes, because if we
        // don't, we are not adding a new parameter to the model (we are simply
        // sliding the node heightdown). So, we can't assign all the children
        // to one subset.
        // This CAN happen when there are 2 or more nodes mapped to the height,
        // because we will always break them up into 2 subsets (which means we
        // will always maintain the old node height AND end up with a new
        // younger height).
        double split_singleton_polytomy(
                RandomNumberGenerator & rng,
                std::shared_ptr<NodeType> polytomy_node,
                std::shared_ptr<PositiveRealParameter> new_height_parameter,
                const bool refresh_node_heights = false,
                const bool refresh_node_ordering = true) {
            unsigned int n_children = polytomy_node->get_number_of_children();
            std::vector< std::vector<unsigned int> > child_subsets;
            // Need to avoid the partitions where all children are assigned to
            // their own subset or they are all assigned to one subset. In
            // these cases, there would be no split.
            while ((child_subsets.size() == 0) ||
                    (child_subsets.size() == 1) ||
                    (child_subsets.size() == n_children)) {
                child_subsets = rng.random_set_partition_as_subsets(
                        n_children);
            }
            // Need to get the pointers to the children to split (can't
            // work with child indices, because these will change as
            // they are split off from polytomy
            std::vector< std::vector< std::shared_ptr<NodeType> > > child_node_subsets;
            for (auto child_subset : child_subsets) {
                // Any singleton children remain at polytomy node
                if (child_subset.size() < 2) {
                    continue;
                }
                // Any groups of children are split from polytomy node
                child_node_subsets.push_back(polytomy_node->get_children(child_subset));
            }
            double ln_prob_new_vals = 0.0;
            for (unsigned int subset_index = 0;
                    subset_index < child_node_subsets.size();
                    ++subset_index) {
                std::shared_ptr<NodeType> new_node = polytomy_node->split_children_from_polytomy(
                        rng,
                        child_node_subsets.at(subset_index),
                        new_height_parameter,
                        this->get_leaf_node_count());
                ln_prob_new_vals += this->get_ln_prob_of_drawing_node_state(new_node);
            }
            if (refresh_node_heights) {
                this->update_node_heights();
            }
            else {
                this->node_heights_.push_back(new_height_parameter);
                this->sort_node_heights();
            }
            if (refresh_node_ordering) {
                this->refresh_ordered_nodes();
            }
            return ln_prob_new_vals;
        }

        bool slide_bump_height_(
                RandomNumberGenerator & rng,
                const unsigned int height_index,
                const double new_height,
                unsigned int collisions_swap_nodes = 0) {
            ECOEVOLITY_ASSERT(new_height >= 0.0);
            if (this->root_height_is_fixed() && (new_height > this->get_root_height())) {
                return false;
            }
            bool getting_older = true;
            if (this->get_height(height_index) > new_height) {
                getting_older = false;
            }
            unsigned int current_idx = height_index;
            unsigned int next_height_idx;
            if (getting_older) {
                while (true) {
                    if (current_idx == (this->get_number_of_node_heights() - 1)) {
                        this->node_heights_.at(current_idx)->set_value(new_height);
                        this->sort_node_heights();
                        this->make_dirty();
                        return true;
                    }
                    next_height_idx = this->get_index_of_youngest_parent(current_idx);
                    if (this->get_height(next_height_idx) > new_height) {
                        this->node_heights_.at(current_idx)->set_value(new_height);
                        this->sort_node_heights();
                        this->make_dirty();
                        return true;
                    }
                    this->node_heights_.at(current_idx)->set_value(
                            this->node_heights_.at(next_height_idx)->get_value()
                            );
                    if (collisions_swap_nodes == 1) {
                        this->slide_collision_node_permute(rng, next_height_idx, current_idx);
                    }
                    else if (collisions_swap_nodes == 2) {
                        this->slide_collision_node_swap(rng, next_height_idx, current_idx);
                    }
                    else if (collisions_swap_nodes == 3) {
                        this->slide_collision_node_swap_all(rng, next_height_idx, current_idx);
                    }
                    current_idx = next_height_idx;
                }
            }
            // Younger nodes to bump down
            while (true) {
                if (current_idx == 0) {
                    this->node_heights_.at(current_idx)->set_value(new_height);
                    this->sort_node_heights();
                    this->make_dirty();
                    return true;
                }
                std::shared_ptr<NodeType> oldest_child = this->get_oldest_child(current_idx);
                if (oldest_child->is_leaf()) {
                    this->node_heights_.at(current_idx)->set_value(new_height);
                    this->sort_node_heights();
                    this->make_dirty();
                    return true;
                }
                next_height_idx = this->get_node_height_index(
                        oldest_child->get_height_parameter());
                if (this->get_height(next_height_idx) < new_height) {
                    this->node_heights_.at(current_idx)->set_value(new_height);
                    this->sort_node_heights();
                    this->make_dirty();
                    return true;
                }
                this->node_heights_.at(current_idx)->set_value(
                        this->node_heights_.at(next_height_idx)->get_value()
                        );
                if (collisions_swap_nodes == 1) {
                    this->slide_collision_node_permute(rng, current_idx, next_height_idx);
                }
                else if (collisions_swap_nodes == 2) {
                    this->slide_collision_node_swap(rng, current_idx, next_height_idx);
                }
                else if (collisions_swap_nodes == 3) {
                    this->slide_collision_node_swap_all(rng, current_idx, next_height_idx);
                }
                current_idx = next_height_idx;
            }
        }

        void build_from_path_(const std::string & path,
                const std::string & ncl_file_format) {
            std::ifstream in_stream;
            in_stream.open(path);
            if (! in_stream.is_open()) {
                throw EcoevolityParsingError(
                        "Could not open tree file",
                        path);
            }
            try {
                this->build_from_stream_(in_stream, ncl_file_format);
            }
            catch(...) {
                std::cerr << "ERROR: Problem parsing tree file path: "
                        << path << "\n";
                throw;
            }
        }

        void build_from_stream_(std::istream & tree_stream,
                const std::string & ncl_file_format = "relaxedphyliptree",
                double ultrametricity_tolerance = 1e-6) {
            MultiFormatReader nexus_reader(-1, NxsReader::WARNINGS_TO_STDERR);
            try {
                nexus_reader.ReadStream(tree_stream, ncl_file_format.c_str());
            }
            catch(...) {
                nexus_reader.DeleteBlocksFromFactories();
                throw;
            }
            unsigned int num_taxa_blocks = nexus_reader.GetNumTaxaBlocks();
            ECOEVOLITY_ASSERT(num_taxa_blocks == 1);
            NxsTaxaBlock * taxa_block = nexus_reader.GetTaxaBlock(0);

            unsigned int num_tree_blocks = nexus_reader.GetNumTreesBlocks(taxa_block);
            ECOEVOLITY_ASSERT(num_tree_blocks == 1);

            NxsTreesBlock * tree_block = nexus_reader.GetTreesBlock(taxa_block, 0);
            unsigned int num_trees = tree_block->GetNumTrees();
            ECOEVOLITY_ASSERT(num_trees > 0);

            const NxsFullTreeDescription & tree_description = tree_block->GetFullTreeDescription(0);
            // std::cout << "Tree is processed: " << tree_description.IsProcessed() << "\n";
            // std::cout << "Tree is rooted: " << tree_description.IsRooted() << "\n";
            this->build_from_ncl_tree_description_(tree_description,
                    taxa_block,
                    ultrametricity_tolerance);

            nexus_reader.DeleteBlocksFromFactories();
        }

        void build_from_ncl_tree_description_(
                const NxsFullTreeDescription & tree_description,
                NxsTaxaBlock * taxa_block,
                double ultrametricity_tolerance = 1e-6) {
            // Tree must be processed to create the NxsSimpleTree
            if (! tree_description.IsProcessed()) {
                throw EcoevolityError("Input tree was not processed by NCL");
            }

            int default_int_edge_length = 0;
            double default_double_edge_length = 0.0;
            NxsSimpleTree simple_tree(tree_description,
                    default_int_edge_length,
                    default_double_edge_length);
            bool tree_is_ultrametric = this->parsed_tree_is_ultrametric_(
                    simple_tree,
                    ultrametricity_tolerance);
            if (! tree_is_ultrametric) {
                throw EcoevolityError("Input tree not ultrametric");
            }
            const std::vector<NxsSimpleNode *> & leaves = simple_tree.GetLeavesRef();
            unsigned int num_leaves = leaves.size();
            int num_leaves_int = leaves.size();
            // std::cout << "Number of leaves: " << num_leaves << "\n";
            std::vector<std::string> leaf_labels;
            leaf_labels.reserve(num_leaves);
            for (auto leaf_node : leaves) {
                NxsString label = taxa_block->GetTaxonLabel(leaf_node->GetTaxonIndex());
                leaf_labels.push_back(label);
            }
            // Sort labels to ensure we always give the same label the same
            // leaf node index
            std::sort(leaf_labels.begin(), leaf_labels.end());
            // Create map of labels to node indices
            std::unordered_map<std::string, int> leaf_label_to_index_map;
            leaf_label_to_index_map.reserve(num_leaves);
            for (int i = 0; i < num_leaves_int; ++i) {
                ECOEVOLITY_ASSERT(leaf_label_to_index_map.count(leaf_labels.at(i)) == 0);
                leaf_label_to_index_map[leaf_labels.at(i)] = i;
            }
            const NxsSimpleNode * simple_root = simple_tree.GetRootConst();
            NxsSimpleEdge root_edge = simple_root->GetEdgeToParent();
            std::map<std::string, std::string> root_info;
            for (auto nxs_comment : root_edge.GetUnprocessedComments()) {
                // std::cout << "Root comment: " << nxs_comment.GetText() << "\n";
                std::string comment = string_util::strip(nxs_comment.GetText());
                if (string_util::startswith(comment, "&")) {
                    std::string root_info_str = comment.substr(1);
                    string_util::parse_map(
                            root_info_str,
                            root_info,
                            ',',
                            '=');
                }
            }
            // Preferably, we want to get the node heights and node height
            // indices from the metadata in the comments.
            // We can check the root if this information exists.
            // If so, we shoud throw an error if any other nodes lacks these
            // data.
            // If the root does not have these data, we could still parse the
            // tree and calculate heights from the edge lengths; without height
            // indices, the resulting tree would have no shared node heights.
            bool using_height_comments = false;
            if (root_info.count("height_index") > 0) {
                using_height_comments = true;
            }
            // For storing heights so they can be shared; only needed if using
            // comments
            std::map<unsigned int, std::shared_ptr<PositiveRealParameter> > indices_to_heights;
            int next_internal_index = num_leaves;

            // Make the root node
            std::shared_ptr<NodeType> root = this->create_internal_node_(
                    simple_root,
                    indices_to_heights,
                    using_height_comments,
                    next_internal_index);
            ++next_internal_index;

            // Recursively add nodes down the tree toward the leaves
            this->add_child_nodes_(simple_root,
                    taxa_block,
                    root,
                    next_internal_index,
                    indices_to_heights,
                    leaf_label_to_index_map,
                    using_height_comments);

            this->set_root(root);
        }

        bool parsed_tree_is_ultrametric_(const NxsSimpleTree & simple_tree,
                double proportional_tolerance = 1e-6) {
            std::vector< std::vector<double> > pairwise_dists = simple_tree.GetDblPathDistances(false);
            std::vector< std::vector<double> > dists_to_mrca = simple_tree.GetDblPathDistances(true);
            unsigned int num_leaves = dists_to_mrca.size();
            ECOEVOLITY_ASSERT(dists_to_mrca.at(0).size() == num_leaves);

            // First find the maximum distance between two tips
            double max_pairwise_dist = -1.0;
            for (unsigned int i = 0; i < num_leaves - 1; ++i) {
                for (unsigned int j = i + 1; j < num_leaves; ++j) {
                    if (pairwise_dists.at(i).at(j) > max_pairwise_dist) {
                        max_pairwise_dist = pairwise_dists.at(i).at(j);
                    }
                }
            }
            // std::cout << "Max pairwise distance: " << max_pairwise_dist << "\n";

            // Get tolerance proportional to the maximum root height
            double max_root_height = max_pairwise_dist / 2.0;
            double abs_tol = max_root_height * proportional_tolerance;

            // Check if distance from leaf i to the MRCA of i and j is (almost)
            // equal to the distance from leaf j to the MRCA of i and j. If
            // this is true for all pairs of tips, then the tree is
            // ultrametric.
            for (unsigned int i = 0; i < num_leaves - 1; ++i) {
                for (unsigned int j = i + 1; j < num_leaves; ++j) {
                    double height_diff = fabs(dists_to_mrca.at(i).at(j) - dists_to_mrca.at(j).at(i));
                    if (height_diff > abs_tol) {
                        return false;
                    }
                }
            }
            return true;
        }

        void add_child_nodes_(const NxsSimpleNode * parent_ncl_node,
                const NxsTaxaBlock * taxa_block,
                std::shared_ptr<NodeType> parent_node,
                int & next_internal_index,
                std::map<unsigned int, std::shared_ptr<PositiveRealParameter> > & indices_to_heights,
                std::unordered_map<std::string, int> & leaf_label_to_index_map,
                const bool using_height_comments) {
            for (auto child_ncl_node : parent_ncl_node->GetChildren()) {
                if (! child_ncl_node->GetFirstChild()) {
                   // No children; this is a leaf
                    std::shared_ptr<NodeType> leaf = this->create_leaf_node_(
                            child_ncl_node,
                            taxa_block,
                            leaf_label_to_index_map);
                    parent_node->add_child(leaf);
                }
                else {
                    std::shared_ptr<NodeType> node = this->create_internal_node_(child_ncl_node,
                            indices_to_heights,
                            using_height_comments,
                            next_internal_index);
                    ++next_internal_index;
                    parent_node->add_child(node);
                    // Need to continue to recurse down the tree, toward the
                    // tips
                    this->add_child_nodes_(child_ncl_node, taxa_block,
                            node,
                            next_internal_index,
                            indices_to_heights,
                            leaf_label_to_index_map,
                            using_height_comments);
                }
            }
        }

        void parse_node_comments_(
                const NxsSimpleNode * ncl_node,
                std::map<std::string, std::string> & comment_map) {
            NxsSimpleEdge ncl_edge = ncl_node->GetEdgeToParent();
            for (auto nxs_comment : ncl_edge.GetUnprocessedComments()) {
                std::string raw_comment = string_util::strip(nxs_comment.GetText());
                if (string_util::startswith(raw_comment, "&")) {
                    std::string comment = raw_comment.substr(1);
                    string_util::parse_map(
                            comment,
                            comment_map,
                            ',',
                            '=');
                }
            }
        }

        std::shared_ptr<NodeType> create_leaf_node_(
                const NxsSimpleNode * ncl_node,
                const NxsTaxaBlock * taxa_block,
                std::unordered_map<std::string, int> & leaf_label_to_index_map) {
            ECOEVOLITY_ASSERT(! ncl_node->GetFirstChild());
            std::map<std::string, std::string> comment_map;
            this->parse_node_comments_(ncl_node, comment_map);
            NxsString leaf_label = taxa_block->GetTaxonLabel(ncl_node->GetTaxonIndex());
            // std::cout << leaf_label << "\n";
            std::shared_ptr<NodeType> leaf = std::make_shared<NodeType>(
                    leaf_label_to_index_map[leaf_label],
                    leaf_label,
                    0.0);
            leaf->fix_node_height();
            leaf->extract_data_from_node_comments(comment_map);
            return leaf;
        }

        std::shared_ptr<NodeType> create_internal_node_(
                const NxsSimpleNode * ncl_node,
                std::map<unsigned int, std::shared_ptr<PositiveRealParameter> > & indices_to_heights,
                const bool using_height_comments,
                const unsigned int internal_index) {
            ECOEVOLITY_ASSERT(ncl_node->GetFirstChild());
            std::map<std::string, std::string> comment_map;
            this->parse_node_comments_(ncl_node, comment_map);
            double height;
            unsigned int height_index;
            if (using_height_comments) {
                std::stringstream h_converter(comment_map["height"]);
                if (! (h_converter >> height)) {
                    throw EcoevolityError("could not convert node height \'" +
                            h_converter.str() + "\'");
                }
                std::stringstream i_converter(comment_map["height_index"]);
                if (! (i_converter >> height_index)) {
                    throw EcoevolityError("could not convert node height index \'" +
                            i_converter.str() + "\'");
                }
            }
            else {
                // Need to get height from edge lengths
                height = this->get_simple_node_height_(ncl_node);
            }
            std::shared_ptr<NodeType> node = std::make_shared<NodeType>(
                    internal_index,
                    height);
            if (using_height_comments) {
                if (indices_to_heights.count(height_index) > 0) {
                    node->set_height_parameter(indices_to_heights[height_index]);
                }
                else {
                    indices_to_heights[height_index] = node->get_height_parameter();
                }
            }
            node->extract_data_from_node_comments(comment_map);
            return node;
        }
        
        double get_simple_node_height_(const NxsSimpleNode * node) {
            double height = 0.0;
            NxsSimpleNode * child = node->GetFirstChild();
            while (child) {
                height += child->GetEdgeToParent().GetDblEdgeLen();
                child = child->GetFirstChild();
            }
            return height;
        }

        void add_edge_length(const NxsSimpleNode * node, double & height) {
            double l = node->GetEdgeToParent().GetDblEdgeLen();
            height += l;
        }


    public:
        BaseTree() { }
        BaseTree(std::shared_ptr<NodeType> root) {
            this->set_root(root);
        }
        BaseTree(std::istream & tree_stream,
                const std::string & ncl_file_format) : BaseTree() {
            this->build_from_stream_(tree_stream, ncl_file_format);
        }
        BaseTree(const std::string & path,
                const std::string & ncl_file_format) : BaseTree() {
            this->build_from_path_(path, ncl_file_format);
        }
        BaseTree(const std::string & newick_tree_string) : BaseTree() {
            std::istringstream tree_stream(newick_tree_string);
            try {
                this->build_from_stream_(tree_stream, "relaxedphyliptree");
            }
            catch(...) {
                std::cerr << "ERROR: Problem parsing newick tree string:\n"
                        << newick_tree_string << "\n";
                throw;
            }
        }

        typedef std::shared_ptr<NodeType> NodePtr;

        template<class TreeType>
        static void get_trees(
                std::istream & tree_stream,
                const std::string & ncl_file_format,
                std::vector<TreeType> & trees,
                unsigned int skip = 0,
                double ultrametricity_tolerance = 1e-6
                ) {
            MultiFormatReader nexus_reader(-1, NxsReader::WARNINGS_TO_STDERR);
            try {
                nexus_reader.ReadStream(tree_stream, ncl_file_format.c_str());
            }
            catch(...) {
                nexus_reader.DeleteBlocksFromFactories();
                throw;
            }
            unsigned int num_taxa_blocks = nexus_reader.GetNumTaxaBlocks();
            ECOEVOLITY_ASSERT(num_taxa_blocks == 1);
            NxsTaxaBlock * taxa_block = nexus_reader.GetTaxaBlock(0);

            unsigned int num_tree_blocks = nexus_reader.GetNumTreesBlocks(taxa_block);
            ECOEVOLITY_ASSERT(num_tree_blocks == 1);

            NxsTreesBlock * tree_block = nexus_reader.GetTreesBlock(taxa_block, 0);
            unsigned int num_trees = tree_block->GetNumTrees();
            ECOEVOLITY_ASSERT(num_trees > 0);

            for (unsigned int i = skip; i < num_trees; ++i) {
                const NxsFullTreeDescription & tree_description = tree_block->GetFullTreeDescription(i);
                TreeType t;
                t.build_from_ncl_tree_description_(tree_description,
                        taxa_block,
                        ultrametricity_tolerance);
                trees.push_back(t);
            }
            nexus_reader.DeleteBlocksFromFactories();
        }
        template<class TreeType>
        static std::vector<TreeType> get_trees(
                std::istream & tree_stream,
                const std::string & ncl_file_format,
                unsigned int skip = 0,
                double ultrametricity_tolerance = 1e-6
                ) {
            std::vector<TreeType> trees;
            BaseTree::get_trees(tree_stream,
                    ncl_file_format,
                    trees,
                    skip,
                    ultrametricity_tolerance);
            return trees;
        }
        template<class TreeType>
        static void get_trees(
                const std::string & path,
                const std::string & ncl_file_format,
                std::vector<TreeType> & trees,
                unsigned int skip = 0,
                double ultrametricity_tolerance = 1e-6
                ) {
            std::ifstream in_stream;
            in_stream.open(path);
            if (! in_stream.is_open()) {
                throw EcoevolityParsingError(
                        "Could not open tree file",
                        path);
            }
            try {
                BaseTree::get_trees(in_stream,
                        ncl_file_format,
                        trees,
                        skip,
                        ultrametricity_tolerance);
            }
            catch(...) {
                std::cerr << "ERROR: Problem parsing tree file path: "
                        << path << "\n";
                throw;
            }
        }
        template<class TreeType>
        static std::vector<TreeType> get_trees(
                const std::string & path,
                const std::string & ncl_file_format,
                unsigned int skip = 0,
                double ultrametricity_tolerance = 1e-6
                ) {
            std::vector<TreeType> trees;
            BaseTree::get_trees(path,
                    ncl_file_format,
                    trees,
                    skip,
                    ultrametricity_tolerance);
            return trees;
        }
        template<class TreeType>
        static void get_trees(
                const std::vector<std::string> & paths,
                const std::string & ncl_file_format,
                std::vector<TreeType> & trees,
                unsigned int skip = 0,
                double ultrametricity_tolerance = 1e-6
                ) {
            for (auto path : paths) {
                BaseTree::get_trees(path,
                        ncl_file_format,
                        trees,
                        skip,
                        ultrametricity_tolerance);
            }
        }
        template<class TreeType>
        static std::vector<TreeType> get_trees(
                const std::vector<std::string> & paths,
                const std::string & ncl_file_format,
                unsigned int skip = 0,
                double ultrametricity_tolerance = 1e-6
                ) {
            std::vector<TreeType> trees;
            BaseTree::get_trees(paths,
                    ncl_file_format,
                    trees,
                    skip,
                    ultrametricity_tolerance);
            return trees;
        }

        virtual double get_ln_prob_of_drawing_node_state(
                std::shared_ptr<NodeType>) const {
            return 0.0;
        }
        // Used to signify that the likelihood needs recalculating.
        // Node methods are expected to make nodes dirty. However, whenever the
        // node height parameters are updated, this needs to be called, because
        // nodes are bypassed (and thus don't end up dirty)
        void make_dirty() {
            this->is_dirty_ = true;
        }
        bool is_dirty() const {
            if (this->is_dirty_) {
                return true;
            }
            return this->root_->clade_has_dirt();
        }
        void make_clean() {
            this->is_dirty_ = false;
            this->root_->make_all_clean();
        }

        void refresh_pre_ordered_nodes() {
            this->root_->pre_order(this->pre_ordered_nodes_);
        }
        void refresh_level_ordered_nodes() {
            this->root_->level_order(this->level_ordered_nodes_);
        }
        void refresh_ordered_nodes() {
            this->refresh_pre_ordered_nodes();
            this->refresh_level_ordered_nodes();
        }

        void update_node_heights() {
            this->node_heights_.clear();
            std::vector< std::shared_ptr<NodeType> > internal_nodes = this->root_->get_internal_nodes();
            for (unsigned int i = 0; i < internal_nodes.size(); ++i) {
                bool exists = false;
                // Check if we already have a pointer to the same place in memory
                for (unsigned int j = 0; j < this->node_heights_.size(); ++j) {
                    if (internal_nodes.at(i)->get_height_parameter() == this->node_heights_.at(j)) {
                        exists = true;
                    }
                }
                if (! exists) {
                    this->node_heights_.push_back(internal_nodes.at(i)->get_height_parameter());
                }
            }
            this->sort_node_heights();
        }

        void sort_node_heights() {
            std::sort(this->node_heights_.begin(), this->node_heights_.end(), PositiveRealParameter::sort_by_value);
        }

        bool root_has_parent() const {
            if (this->root_->has_parent()) {
                return true;
            }
            return false;
        }
        bool root_has_children() const {
            if (this->root_->get_number_of_children() < 1) {
                return false;
            }
            return true;
        }
        bool root_is_valid() const {
            if (this->root_has_parent()) {
                return false;
            }
            if (! this->root_has_children()) {
                return false;
            }
            return true;
        }
        bool node_heights_are_valid() const {
            return this->root_->node_heights_are_valid();
        }
        bool tree_is_valid() const {
            if (! this->root_is_valid()) {
                return false;
            }
            if (! this->node_heights_are_valid()) {
                return false;
            }
            return true;
        }

        void vet_root() const {
            if (this->root_has_parent()) {
                throw EcoevolityError("Root has a parent");
            }
            if (! this->root_has_children()) {
                throw EcoevolityError("Root has no children");
            }
        }
        void vet_tree() const {
            this->vet_root();
            if (! this->node_heights_are_valid()) {
                throw EcoevolityError("Node ages are not valid");
            }
        }

        bool root_height_is_fixed() const {
            return this->root_->node_height_is_fixed();
        }
        void fix_root_height() {
            this->root_->fix_node_height();
        }
        void estimate_root_height() {
            this->root_->estimate_node_height();
        }

        unsigned int get_node_height_index(const std::shared_ptr<PositiveRealParameter> height) const {
            for (unsigned int i = 0; i < this->node_heights_.size(); ++i) {
                if (this->node_heights_.at(i) == height) {
                    return i;
                }
            }
            throw EcoevolityError("Node height does not exist");
        }

        std::vector< std::shared_ptr<NodeType> > get_mapped_nodes(const unsigned int height_index) const {
            return this->root_->get_mapped_nodes(this->node_heights_.at(height_index));
        }

        std::vector< std::shared_ptr<NodeType> > get_mapped_polytomy_nodes(const unsigned int height_index) const {
            return this->root_->get_mapped_polytomy_nodes(this->node_heights_.at(height_index));
        }

        std::vector< std::shared_ptr<NodeType> > get_polytomy_nodes() {
            return this->root_->get_polytomy_nodes();
        }

        double merge_node_height_up(const unsigned int height_index,
                std::vector<unsigned int> & sizes_of_mapped_polytomies_after_merge,
                unsigned int & number_of_resulting_merged_nodes,
                const bool refresh_node_heights = false,
                const bool refresh_node_ordering = true) {
            // Make sure we aren't dealing with the root node
            ECOEVOLITY_ASSERT(height_index < (this->get_number_of_node_heights() - 1));

            sizes_of_mapped_polytomies_after_merge.clear();
            number_of_resulting_merged_nodes = 0;
            double ln_prob_of_drawing_state_of_removed_nodes = 0.0;

            std::set<std::shared_ptr<NodeType> > nodes_of_polytomies_created;
            std::shared_ptr<PositiveRealParameter> new_height = this->node_heights_.at(height_index + 1);
            std::vector< std::shared_ptr<NodeType> > mapped_nodes = this->get_mapped_nodes(height_index);
            for (unsigned int i = 0; i < mapped_nodes.size(); ++i) {
                // If the parent of the node we are moving up is assigned to
                // the next larger node height, we need to add the child to a
                // polytomy
                if (mapped_nodes.at(i)->get_parent()->get_height_parameter() == new_height) {
                    nodes_of_polytomies_created.insert(mapped_nodes.at(i)->get_parent());
                    ln_prob_of_drawing_state_of_removed_nodes += this->get_ln_prob_of_drawing_node_state(
                            mapped_nodes.at(i));
                    mapped_nodes.at(i)->collapse();
                }
                else {
                    mapped_nodes.at(i)->set_height_parameter(new_height);
                    ++number_of_resulting_merged_nodes;
                    unsigned int n_children = mapped_nodes.at(i)->get_number_of_children();
                    ECOEVOLITY_ASSERT(n_children > 1);
                    if (n_children > 2) {
                        // We moved an existing polytomy up; we need to keep
                        // track of this polytomies size for calculating the
                        // probability of the reverse split move.
                        sizes_of_mapped_polytomies_after_merge.push_back(n_children);
                    }
                }
            }
            for (auto poly_node : nodes_of_polytomies_created) {
                sizes_of_mapped_polytomies_after_merge.push_back(poly_node->get_number_of_children());
                ++number_of_resulting_merged_nodes;
            }
            if (refresh_node_heights) {
                this->update_node_heights();
            }
            else {
                this->node_heights_.erase(this->node_heights_.begin() + height_index);
            }
            if (refresh_node_ordering) {
                this->refresh_ordered_nodes();
            }
            return ln_prob_of_drawing_state_of_removed_nodes;
        }

        double split_node_height_down(
                RandomNumberGenerator & rng,
                const unsigned int height_index,
                double & height_lower_bound,
                unsigned int & number_of_mapped_nodes,
                unsigned int & number_of_nodes_in_split_subset,
                std::vector<unsigned int> & moving_polytomy_sizes,
                bool & mapped_nodes_include_polytomy,
                const bool refresh_node_heights = false,
                const bool refresh_node_ordering = true) {
            moving_polytomy_sizes.clear();
            std::vector< std::shared_ptr<NodeType> > mapped_nodes = this->get_mapped_nodes(height_index);
            number_of_mapped_nodes = mapped_nodes.size();
            if (mapped_nodes.size() < 1) {
                number_of_nodes_in_split_subset = 0;
                mapped_nodes_include_polytomy = false;
                return 0.0;
            }
            double max_height = this->node_heights_.at(height_index)->get_value();
            double min_height = 0.0;
            if (height_index > 0) {
                min_height = this->node_heights_.at(height_index - 1)->get_value();
            }
            height_lower_bound = min_height;
            double new_height = rng.uniform_real(min_height, max_height);
            std::shared_ptr<PositiveRealParameter> new_height_parameter = std::make_shared<PositiveRealParameter>(new_height);
            if (mapped_nodes.size() == 1) {
                // If we only have a single polytomy, we need to handle the
                // splitting differently
                ECOEVOLITY_ASSERT(mapped_nodes.at(0)->is_polytomy());
                moving_polytomy_sizes.push_back(mapped_nodes.at(0)->get_number_of_children());
                number_of_nodes_in_split_subset = 1;
                mapped_nodes_include_polytomy = false;
                return this->split_singleton_polytomy(rng,
                        mapped_nodes.at(0),
                        new_height_parameter,
                        refresh_node_heights,
                        refresh_node_ordering);
            }

            mapped_nodes_include_polytomy = false;
            for (auto n : mapped_nodes) {
                if (n->is_polytomy()) {
                    mapped_nodes_include_polytomy = true;
                    break;
                }
            }

            if (! mapped_nodes_include_polytomy) {
                // We only have shared bifurcating nodes
                std::vector< std::vector<unsigned int> > subsets = rng.random_subsets(
                        mapped_nodes.size(), 2);

                unsigned int move_subset_index = rng.uniform_positive_int(0, 1);

                number_of_nodes_in_split_subset = subsets.at(move_subset_index).size();
                for (auto node_index : subsets.at(move_subset_index)) {
                    mapped_nodes.at(node_index)->set_height_parameter(new_height_parameter);
                }
                if (refresh_node_heights) {
                    this->update_node_heights();
                }
                else {
                    this->node_heights_.push_back(new_height_parameter);
                    this->sort_node_heights();
                }
                if (refresh_node_ordering) {
                    this->refresh_ordered_nodes();
                }
                return 0.0;
            }

            // We have at least one polytomy node included in the shared nodes.
            // So, we need to allow all nodes to end up in the move set.
            bool all_nodes_in_move_set = true;
            bool all_poly_nodes_will_move_in_full = true;
            std::vector<unsigned int> move_subset;
            std::vector<unsigned int> allowed_numbers_of_subsets {1, 2};
            std::map<unsigned int, std::vector< std::vector<unsigned int> > > node_index_to_child_subsets;
            while (all_nodes_in_move_set && all_poly_nodes_will_move_in_full) {
                // This while loop ensures that we do not end up simply sliding
                // all the nodes down to the new height. This would not add a
                // parameter to the model and cannot be reversed by a merge
                // move. If this happens we simply reject and continue, so that
                // all other possible ways to split will be sampled uniformly.
                while (move_subset.size() < 1) {
                    // This while loop ensures that we do not sample an empty
                    // set to move. If we do we continue until we don't. By
                    // using rejection, all other possible move sets will be
                    // sample uniformly.
                    all_nodes_in_move_set = false;
                    std::vector< std::vector<unsigned int> > subsets = rng.restricted_random_set_partition_as_subsets(
                            mapped_nodes.size(),
                            allowed_numbers_of_subsets);
                    unsigned int move_subset_index = rng.uniform_positive_int(0, 1);
                    if (subsets.size() == 1) {
                        all_nodes_in_move_set = true;
                        if (move_subset_index == 1) {
                            // We have sampled an empty set to move, which we avoid by
                            // rejecting and continuing the loop until we get a
                            // non-empty move set
                            continue;
                        }
                    }
                    move_subset = subsets.at(move_subset_index);
                }

                for (auto node_index : move_subset) {
                    // If node is a polytomy, we need to deal with possibility of
                    // splitting it up
                    if (mapped_nodes.at(node_index)->is_polytomy()) {
                        unsigned int n_children = mapped_nodes.at(node_index)->get_number_of_children();
                        std::vector< std::vector<unsigned int> > child_subsets;
                        // Need to avoid the partition where all children are
                        // assigned to their own subset. This would result in the
                        // no nodes being moved to the new height (i.e., the
                        // polytomy remains as is), which complicates Hasting's
                        // ratio, because then there are multiple ways polytomy
                        // nodes are not included in the split (it can either not
                        // be selected in the move pool, or selected, but all the
                        // children are in singleton subsets and thus not split
                        // off).
                        while ((child_subsets.size() == 0) || (child_subsets.size() == n_children)) {
                            child_subsets = rng.random_set_partition_as_subsets(
                                    n_children);
                        }
                        if (child_subsets.size() > 1) {
                            // We have at least one polytomy node that will leave a
                            // node at the old node height (we will have a
                            // model-jump move)
                            all_poly_nodes_will_move_in_full = false;
                        }
                        node_index_to_child_subsets[node_index] = child_subsets;
                    }
                }
            }

            double ln_prob_new_vals = 0.0;
            number_of_nodes_in_split_subset = move_subset.size();
            for (auto node_index : move_subset) {
                // If node is a polytomy, we need to deal with splitting it up
                if (mapped_nodes.at(node_index)->is_polytomy()) {
                    unsigned int n_children = mapped_nodes.at(node_index)->get_number_of_children();
                    moving_polytomy_sizes.push_back(n_children);
                    std::vector< std::vector<unsigned int> > child_subsets = node_index_to_child_subsets[node_index];
                    std::vector< std::vector< std::shared_ptr<NodeType> > > child_node_subsets;
                    // Need to get the pointers to the children to split (can't
                    // work with child indices, because these will change as
                    // they are split off from polytomy
                    for (auto child_subset : child_subsets) {
                        // Any singleton children remain at polytomy node
                        if (child_subset.size() < 2) {
                            continue;
                        }
                        // Any groups of children are split from polytomy node
                        child_node_subsets.push_back(mapped_nodes.at(node_index)->get_children(child_subset));
                    }
                    for (unsigned int subset_index = 0;
                            subset_index < child_node_subsets.size();
                            ++subset_index) {
                        std::shared_ptr<NodeType> new_node = mapped_nodes.at(
                                node_index)->split_children_from_polytomy(
                                    rng,
                                    child_node_subsets.at(subset_index),
                                    new_height_parameter,
                                    this->get_leaf_node_count());
                        ln_prob_new_vals += this->get_ln_prob_of_drawing_node_state(new_node);
                    }
                }
                // Node is not a polytomy, so we simply assign it to the new height
                else {
                    mapped_nodes.at(node_index)->set_height_parameter(new_height_parameter);
                }
            }
            if (refresh_node_heights) {
                this->update_node_heights();
            }
            else {
                this->node_heights_.push_back(new_height_parameter);
                this->sort_node_heights();
            }
            if (refresh_node_ordering) {
                this->refresh_ordered_nodes();
            }
            return ln_prob_new_vals;
        }

        std::vector<unsigned int> get_indices_of_splittable_heights() const {
            std::vector<unsigned int> splittable_heights;
            for (unsigned int i = 0; i < this->node_heights_.size(); ++i) {
                if (this->height_is_splittable(i)) {
                    splittable_heights.push_back(i);
                }
            }
            return splittable_heights;
        }

        unsigned int get_number_of_splittable_heights() const {
            unsigned int num_splittable_heights = 0;
            for (unsigned int i = 0; i < this->node_heights_.size(); ++i) {
                if (this->height_is_splittable(i)) {
                    ++num_splittable_heights;
                }
            }
            return num_splittable_heights;
        }

        bool height_is_splittable(const unsigned int height_index) const {
            unsigned int mapped_node_count = this->root_->get_mapped_node_count(
                    this->node_heights_.at(height_index));
            if (mapped_node_count > 1) {
                return true;
            }
            unsigned int mapped_polytomy_node_count = this->root_->get_mapped_polytomy_node_count(
                    this->node_heights_.at(height_index));
            if (mapped_polytomy_node_count > 0) {
                return true;
            }
            return false;
        }

        unsigned int get_mapped_node_count(const unsigned int height_index) const {
            return this->root_->get_mapped_node_count(this->node_heights_.at(height_index));
        }
        unsigned int get_mapped_polytomy_node_count(const unsigned int height_index) const {
            return this->root_->get_mapped_polytomy_node_count(this->node_heights_.at(height_index));
        }

        std::vector< std::shared_ptr<NodeType> > get_collision_parents(
                const unsigned int older_height_index,
                const unsigned int younger_height_index) {
            std::vector< std::shared_ptr<NodeType> > collision_parents;
            std::vector< std::shared_ptr<NodeType> > older_nodes = this->get_mapped_nodes(older_height_index);
            std::vector< std::shared_ptr<NodeType> > younger_nodes = this->get_mapped_nodes(younger_height_index);
            for (auto older_nd : older_nodes) {
                for (auto younger_nd : younger_nodes) {
                    if (older_nd->is_child(younger_nd)) {
                        collision_parents.push_back(older_nd);
                        break;
                    }
                }
            }
            return collision_parents;
        }

        void slide_collision_node_permute(
                RandomNumberGenerator & rng,
                const unsigned int older_height_index,
                const unsigned int younger_height_index,
                const bool refresh_node_ordering = true) {
            std::vector< std::shared_ptr<NodeType> > collision_parents = this->get_collision_parents(
                    older_height_index,
                    younger_height_index);
            if (collision_parents.size() < 1) {
                std::ostringstream message;
                message << "ERROR: slide_collision_node_permute: "
                        << "Height " << older_height_index
                        << " has no parents of " << younger_height_index;
                throw EcoevolityError(message.str());
            }
            bool change_happened = false;
            while (! change_happened) {
                std::map< std::shared_ptr<NodeType>, std::shared_ptr<NodeType> > parent_swap_node_map;
                for (auto parent_nd : collision_parents) {
                    std::vector< std::shared_ptr<NodeType> > swap_node_pool;
                    std::vector< std::shared_ptr<NodeType> > child_colliding_nodes;
                    std::vector< std::shared_ptr<NodeType> > child_non_colliding_nodes;
                    for (auto child_nd : parent_nd->get_all_children()) {
                        // Each colliding child contributes one random grandchild to
                        // the swap pool
                        if ((! child_nd->is_leaf()) &&
                                (this->get_node_height_index(child_nd->get_height_parameter()) ==
                                 younger_height_index)) {
                            ECOEVOLITY_ASSERT(child_nd->has_children());
                            child_colliding_nodes.push_back(child_nd);
                            unsigned int random_child_index = rng.uniform_positive_int(child_nd->get_number_of_children() - 1);
                            swap_node_pool.push_back(child_nd->get_child(random_child_index));
                            parent_swap_node_map[child_nd] = child_nd->get_child(random_child_index);
                        }
                        else {
                            child_non_colliding_nodes.push_back(child_nd);
                        }
                    }
                    // If parent node in collision has non-colliding descendants,
                    // it contributes one of these descendants (randomly) to the
                    // swap pool
                    if (child_non_colliding_nodes.size() > 0) {
                        unsigned int random_child_index = rng.uniform_positive_int(child_non_colliding_nodes.size() - 1);
                        swap_node_pool.push_back(child_non_colliding_nodes.at(random_child_index));
                        parent_swap_node_map[parent_nd] = child_non_colliding_nodes.at(random_child_index);
                    }

                    // Randomize the swap pool
                    std::shuffle(std::begin(swap_node_pool), std::end(swap_node_pool), rng.engine_);

                    std::shared_ptr<NodeType> random_swap_node;
                    // If the parent contributed a non-colliding child to the swap
                    // pool, randomly assign it a new child from the swap pool
                    if (child_non_colliding_nodes.size() > 0) {
                        std::shared_ptr<NodeType> random_swap_node = swap_node_pool.back();
                        swap_node_pool.pop_back();
                        random_swap_node->remove_parent();
                        random_swap_node->add_parent(parent_nd);
                    }
                    for (auto child_colliding_nd : child_colliding_nodes) {
                        random_swap_node = swap_node_pool.back();
                        swap_node_pool.pop_back();
                        random_swap_node->remove_parent();
                        random_swap_node->add_parent(child_colliding_nd);
                    }
                    ECOEVOLITY_ASSERT(swap_node_pool.size() == 0);
                }
                // The most probable outcome is for the swap nodes to be
                // reassigned back to their original parents (because this is
                // always possible for every set of nodes in the swap pool).
                // To avoid returning with no change having occurred, we'll
                // use rejection sampling until we get a change.
                //
                // The simplest example is the case of:
                //   ((A,B),C)
                // C will always be added to the swap pool along with A (prob
                // 1/2) or B (prob 1/2). So the probability of all three
                // possible outcomes are:
                //   ((A,C),B)  p(B added to pool) X p(root draws B) = 1/2 * 1/2 = 1/4
                //   ((B,C),A)  p(A added to pool) X p(root draws A) = 1/2 * 1/2 = 1/4
                //   ((A,B),C)  p(A OR B added to pool) X p(root draws C) = 1 * 1/2 = 1/2
                //
                // NOTE: this simplest example is also the maximum probability
                // of ending up in the unchanged state. Since in the worst case
                // scenario has a prob of no change of only 1/2, rejection
                // sampling is a pretty cheap way of handling this.
                //
                // ALSO NOTE: in this simplest example this is equivalent to a
                // NNI move
                for (auto parent_child : parent_swap_node_map) {
                    if (! parent_child.first->is_child(parent_child.second)) {
                        change_happened = true;
                        break;
                    }
                }
            }
            if (refresh_node_ordering) {
                this->refresh_ordered_nodes();
            }
            this->root_->resize_splits(this->get_leaf_node_count());
        }

        void collision_node_permute(
                RandomNumberGenerator & rng,
                const unsigned int older_height_index,
                const unsigned int younger_height_index,
                const bool refresh_node_ordering = true) {
            std::vector< std::shared_ptr<NodeType> > collision_parents = this->get_collision_parents(
                    older_height_index,
                    younger_height_index);
            if (collision_parents.size() < 1) {
                std::ostringstream message;
                message << "ERROR: collision_node_permute: "
                        << "Height " << older_height_index
                        << " has no parents of " << younger_height_index;
                throw EcoevolityError(message.str());
            }
            bool change_happened = false;
            while (! change_happened) {
                std::map< std::shared_ptr<NodeType>, std::shared_ptr<NodeType> > parent_swap_node_map;
                for (auto parent_nd : collision_parents) {
                    std::vector< std::shared_ptr<NodeType> > swap_node_pool;
                    std::vector< std::shared_ptr<NodeType> > child_colliding_nodes;
                    std::vector< std::shared_ptr<NodeType> > child_non_colliding_nodes;
                    std::shared_ptr<NodeType> non_colliding_node_parent;
                    for (auto child_nd : parent_nd->get_all_children()) {
                        // Each colliding child contributes one random grandchild to
                        // the swap pool
                        if ((! child_nd->is_leaf()) &&
                                (this->get_node_height_index(child_nd->get_height_parameter()) ==
                                 younger_height_index)) {
                            ECOEVOLITY_ASSERT(child_nd->has_children());
                            child_colliding_nodes.push_back(child_nd);
                            unsigned int random_child_index = rng.uniform_positive_int(child_nd->get_number_of_children() - 1);
                            swap_node_pool.push_back(child_nd->get_child(random_child_index));
                            parent_swap_node_map[child_nd] = child_nd->get_child(random_child_index);
                        }
                        else {
                            child_non_colliding_nodes.push_back(child_nd);
                        }
                    }
                    // If parent node in collision has non-colliding descendants,
                    // it contributes one of these descendants (randomly) to the
                    // swap pool
                    if (child_non_colliding_nodes.size() > 0) {
                        unsigned int random_child_index = rng.uniform_positive_int(child_non_colliding_nodes.size() - 1);
                        std::shared_ptr<NodeType> candidate_swap_node = child_non_colliding_nodes.at(random_child_index);
                        // Unlike for the slide_collision_node_permute, we have to make
                        // sure the node we add to the swap pool is younger than the
                        // younger height. We keep randomly grabbing children down the
                        // tree until this is true
                        while (candidate_swap_node->get_height() > this->get_height(younger_height_index)) {
                            random_child_index = rng.uniform_positive_int(candidate_swap_node->get_number_of_children() - 1);
                            candidate_swap_node = candidate_swap_node->get_child(random_child_index);
                        }
                        swap_node_pool.push_back(candidate_swap_node);
                        non_colliding_node_parent = candidate_swap_node->get_parent();
                        parent_swap_node_map[non_colliding_node_parent] = candidate_swap_node;
                    }

                    // Randomize the swap pool
                    std::shuffle(std::begin(swap_node_pool), std::end(swap_node_pool), rng.engine_);

                    std::shared_ptr<NodeType> random_swap_node;
                    // If the parent contributed a non-colliding child to the swap
                    // pool, randomly assign it a new child from the swap pool
                    if (child_non_colliding_nodes.size() > 0) {
                        std::shared_ptr<NodeType> random_swap_node = swap_node_pool.back();
                        swap_node_pool.pop_back();
                        random_swap_node->remove_parent();
                        random_swap_node->add_parent(non_colliding_node_parent);
                    }
                    for (auto child_colliding_nd : child_colliding_nodes) {
                        random_swap_node = swap_node_pool.back();
                        swap_node_pool.pop_back();
                        random_swap_node->remove_parent();
                        random_swap_node->add_parent(child_colliding_nd);
                    }
                    ECOEVOLITY_ASSERT(swap_node_pool.size() == 0);
                }
                // The most probable outcome is for the swap nodes to be
                // reassigned back to their original parents (because this is
                // always possible for every set of nodes in the swap pool).
                // To avoid returning with no change having occurred, we'll
                // use rejection sampling until we get a change.
                //
                // The simplest example is the case of:
                //   ((A,B),C)
                // C will always be added to the swap pool along with A (prob
                // 1/2) or B (prob 1/2). So the probability of all three
                // possible outcomes are:
                //   ((A,C),B)  p(B added to pool) X p(root draws B) = 1/2 * 1/2 = 1/4
                //   ((B,C),A)  p(A added to pool) X p(root draws A) = 1/2 * 1/2 = 1/4
                //   ((A,B),C)  p(A OR B added to pool) X p(root draws C) = 1 * 1/2 = 1/2
                //
                // NOTE: this simplest example is also the maximum probability
                // of ending up in the unchanged state. Since in the worst case
                // scenario has a prob of no change of only 1/2, rejection
                // sampling is a pretty cheap way of handling this.
                //
                // ALSO NOTE: in this simplest example this is equivalent to a
                // NNI move
                for (auto parent_child : parent_swap_node_map) {
                    if (! parent_child.first->is_child(parent_child.second)) {
                        change_happened = true;
                        break;
                    }
                }
            }
            if (refresh_node_ordering) {
                this->refresh_ordered_nodes();
            }
            this->root_->resize_splits(this->get_leaf_node_count());
        }

        void slide_collision_node_swap_all(
                RandomNumberGenerator & rng,
                const unsigned int older_height_index,
                const unsigned int younger_height_index,
                const bool refresh_node_ordering =  true) {
            std::vector< std::shared_ptr<NodeType> > collision_parents = this->get_collision_parents(
                    older_height_index,
                    younger_height_index);
            if (collision_parents.size() < 1) {
                std::ostringstream message;
                message << "ERROR: slide_collision_node_swap_all: "
                        << "Height " << older_height_index
                        << " has no parents of " << younger_height_index;
                throw EcoevolityError(message.str());
            }
            for (auto parent_nd : collision_parents) {
                std::vector< std::shared_ptr<NodeType> > swap_node_pool;
                std::vector< std::shared_ptr<NodeType> > child_colliding_nodes;
                std::vector< std::shared_ptr<NodeType> > child_non_colliding_nodes;
                for (auto child_nd : parent_nd->get_all_children()) {
                    // Each colliding child contributes one random grandchild to
                    // the swap pool
                    if ((! child_nd->is_leaf()) &&
                            (this->get_node_height_index(child_nd->get_height_parameter()) ==
                             younger_height_index)) {
                        ECOEVOLITY_ASSERT(child_nd->has_children());
                        child_colliding_nodes.push_back(child_nd);
                        unsigned int random_child_index = rng.uniform_positive_int(child_nd->get_number_of_children() - 1);
                        swap_node_pool.push_back(child_nd->get_child(random_child_index));
                    }
                    else {
                        child_non_colliding_nodes.push_back(child_nd);
                    }
                }
                // If parent node in collision has non-colliding descendants,
                // it contributes one of these descendants (randomly) to the
                // swap pool
                if (child_non_colliding_nodes.size() > 0) {
                    unsigned int random_child_index = rng.uniform_positive_int(child_non_colliding_nodes.size() - 1);
                    swap_node_pool.push_back(child_non_colliding_nodes.at(random_child_index));
                }

                // Randomly pick two nodes to swap
                std::vector<unsigned int> indices_of_nodes_to_swap = rng.random_subset_indices(
                        swap_node_pool.size(),
                        2);
                ECOEVOLITY_ASSERT(indices_of_nodes_to_swap.size() == 2);
                std::shared_ptr<NodeType> swap_node1 = swap_node_pool.at(indices_of_nodes_to_swap.at(0));
                std::shared_ptr<NodeType> swap_node2 = swap_node_pool.at(indices_of_nodes_to_swap.at(1));
                std::shared_ptr<NodeType> parent1 = swap_node1->get_parent();
                std::shared_ptr<NodeType> parent2 = swap_node2->get_parent();
                swap_node1->remove_parent();
                swap_node2->remove_parent();
                swap_node1->add_parent(parent2);
                swap_node2->add_parent(parent1);
            }
            if (refresh_node_ordering) {
                this->refresh_ordered_nodes();
            }
            this->root_->resize_splits(this->get_leaf_node_count());
        }

        void collision_node_swap_all(
                RandomNumberGenerator & rng,
                const unsigned int older_height_index,
                const unsigned int younger_height_index,
                const bool refresh_node_ordering =  true) {
            std::vector< std::shared_ptr<NodeType> > collision_parents = this->get_collision_parents(
                    older_height_index,
                    younger_height_index);
            if (collision_parents.size() < 1) {
                std::ostringstream message;
                message << "ERROR: collision_node_swap_all: "
                        << "Height " << older_height_index
                        << " has no parents of " << younger_height_index;
                throw EcoevolityError(message.str());
            }
            for (auto parent_nd : collision_parents) {
                std::vector< std::shared_ptr<NodeType> > swap_node_pool;
                std::vector< std::shared_ptr<NodeType> > child_colliding_nodes;
                std::vector< std::shared_ptr<NodeType> > child_non_colliding_nodes;
                for (auto child_nd : parent_nd->get_all_children()) {
                    // Each colliding child contributes one random grandchild to
                    // the swap pool
                    if ((! child_nd->is_leaf()) &&
                            (this->get_node_height_index(child_nd->get_height_parameter()) ==
                             younger_height_index)) {
                        ECOEVOLITY_ASSERT(child_nd->has_children());
                        child_colliding_nodes.push_back(child_nd);
                        unsigned int random_child_index = rng.uniform_positive_int(child_nd->get_number_of_children() - 1);
                        swap_node_pool.push_back(child_nd->get_child(random_child_index));
                    }
                    else {
                        child_non_colliding_nodes.push_back(child_nd);
                    }
                }
                // If parent node in collision has non-colliding descendants,
                // it contributes one of these descendants (randomly) to the
                // swap pool
                if (child_non_colliding_nodes.size() > 0) {
                    unsigned int random_child_index = rng.uniform_positive_int(child_non_colliding_nodes.size() - 1);
                    std::shared_ptr<NodeType> candidate_swap_node = child_non_colliding_nodes.at(random_child_index);
                    // Unlike for the slide_collision_node_swap, we have to make
                    // sure the node we add to the swap pool is younger than the
                    // younger height. We keep randomly grabbing children down the
                    // tree until this is true
                    while (candidate_swap_node->get_height() > this->get_height(younger_height_index)) {
                        random_child_index = rng.uniform_positive_int(candidate_swap_node->get_number_of_children() - 1);
                        candidate_swap_node = candidate_swap_node->get_child(random_child_index);
                    }
                    swap_node_pool.push_back(candidate_swap_node);
                }

                // Randomly pick two nodes to swap
                std::vector<unsigned int> indices_of_nodes_to_swap = rng.random_subset_indices(
                        swap_node_pool.size(),
                        2);
                ECOEVOLITY_ASSERT(indices_of_nodes_to_swap.size() == 2);
                std::shared_ptr<NodeType> swap_node1 = swap_node_pool.at(indices_of_nodes_to_swap.at(0));
                std::shared_ptr<NodeType> swap_node2 = swap_node_pool.at(indices_of_nodes_to_swap.at(1));
                std::shared_ptr<NodeType> parent1 = swap_node1->get_parent();
                std::shared_ptr<NodeType> parent2 = swap_node2->get_parent();
                swap_node1->remove_parent();
                swap_node2->remove_parent();
                swap_node1->add_parent(parent2);
                swap_node2->add_parent(parent1);
            }
            if (refresh_node_ordering) {
                this->refresh_ordered_nodes();
            }
            this->root_->resize_splits(this->get_leaf_node_count());
        }

        void slide_collision_node_swap(
                RandomNumberGenerator & rng,
                const unsigned int older_height_index,
                const unsigned int younger_height_index,
                const bool refresh_node_ordering =  true) {
            std::vector< std::shared_ptr<NodeType> > collision_parents = this->get_collision_parents(
                    older_height_index,
                    younger_height_index);
            if (collision_parents.size() < 1) {
                std::ostringstream message;
                message << "ERROR: slide_collision_node_swap: "
                        << "Height " << older_height_index
                        << " has no parents of " << younger_height_index;
                throw EcoevolityError(message.str());
            }
            unsigned int random_parent_index = rng.uniform_positive_int(collision_parents.size() - 1);
            std::shared_ptr<NodeType> parent_nd = collision_parents.at(random_parent_index);
            std::vector< std::shared_ptr<NodeType> > swap_node_pool;
            std::vector< std::shared_ptr<NodeType> > child_colliding_nodes;
            std::vector< std::shared_ptr<NodeType> > child_non_colliding_nodes;
            for (auto child_nd : parent_nd->get_all_children()) {
                // Each colliding child contributes one random grandchild to
                // the swap pool
                if ((! child_nd->is_leaf()) &&
                        (this->get_node_height_index(child_nd->get_height_parameter()) ==
                         younger_height_index)) {
                    ECOEVOLITY_ASSERT(child_nd->has_children());
                    child_colliding_nodes.push_back(child_nd);
                    unsigned int random_child_index = rng.uniform_positive_int(child_nd->get_number_of_children() - 1);
                    swap_node_pool.push_back(child_nd->get_child(random_child_index));
                }
                else {
                    child_non_colliding_nodes.push_back(child_nd);
                }
            }
            // If parent node in collision has non-colliding descendants,
            // it contributes one of these descendants (randomly) to the
            // swap pool
            if (child_non_colliding_nodes.size() > 0) {
                unsigned int random_child_index = rng.uniform_positive_int(child_non_colliding_nodes.size() - 1);
                swap_node_pool.push_back(child_non_colliding_nodes.at(random_child_index));
            }

            // Randomly pick two nodes to swap
            std::vector<unsigned int> indices_of_nodes_to_swap = rng.random_subset_indices(
                    swap_node_pool.size(),
                    2);
            ECOEVOLITY_ASSERT(indices_of_nodes_to_swap.size() == 2);
            std::shared_ptr<NodeType> swap_node1 = swap_node_pool.at(indices_of_nodes_to_swap.at(0));
            std::shared_ptr<NodeType> swap_node2 = swap_node_pool.at(indices_of_nodes_to_swap.at(1));
            std::shared_ptr<NodeType> parent1 = swap_node1->get_parent();
            std::shared_ptr<NodeType> parent2 = swap_node2->get_parent();
            swap_node1->remove_parent();
            swap_node2->remove_parent();
            swap_node1->add_parent(parent2);
            swap_node2->add_parent(parent1);
            if (refresh_node_ordering) {
                this->refresh_ordered_nodes();
            }
            this->root_->resize_splits(this->get_leaf_node_count());
        }

        void collision_node_swap(
                RandomNumberGenerator & rng,
                const unsigned int older_height_index,
                const unsigned int younger_height_index,
                const bool refresh_node_ordering =  true) {
            std::vector< std::shared_ptr<NodeType> > collision_parents = this->get_collision_parents(
                    older_height_index,
                    younger_height_index);
            if (collision_parents.size() < 1) {
                std::ostringstream message;
                message << "ERROR: collision_node_swap: "
                        << "Height " << older_height_index
                        << " has no parents of " << younger_height_index;
                throw EcoevolityError(message.str());
            }
            unsigned int random_parent_index = rng.uniform_positive_int(collision_parents.size() - 1);
            std::shared_ptr<NodeType> parent_nd = collision_parents.at(random_parent_index);
            std::vector< std::shared_ptr<NodeType> > swap_node_pool;
            std::vector< std::shared_ptr<NodeType> > child_colliding_nodes;
            std::vector< std::shared_ptr<NodeType> > child_non_colliding_nodes;
            for (auto child_nd : parent_nd->get_all_children()) {
                // Each colliding child contributes one random grandchild to
                // the swap pool
                if ((! child_nd->is_leaf()) &&
                        (this->get_node_height_index(child_nd->get_height_parameter()) ==
                         younger_height_index)) {
                    ECOEVOLITY_ASSERT(child_nd->has_children());
                    child_colliding_nodes.push_back(child_nd);
                    unsigned int random_child_index = rng.uniform_positive_int(child_nd->get_number_of_children() - 1);
                    swap_node_pool.push_back(child_nd->get_child(random_child_index));
                }
                else {
                    child_non_colliding_nodes.push_back(child_nd);
                }
            }
            // If parent node in collision has non-colliding descendants,
            // it contributes one of these descendants (randomly) to the
            // swap pool
            if (child_non_colliding_nodes.size() > 0) {
                unsigned int random_child_index = rng.uniform_positive_int(child_non_colliding_nodes.size() - 1);
                std::shared_ptr<NodeType> candidate_swap_node = child_non_colliding_nodes.at(random_child_index);
                // Unlike for the slide_collision_node_swap, we have to make
                // sure the node we add to the swap pool is younger than the
                // younger height. We keep randomly grabbing children down the
                // tree until this is true
                while (candidate_swap_node->get_height() > this->get_height(younger_height_index)) {
                    random_child_index = rng.uniform_positive_int(candidate_swap_node->get_number_of_children() - 1);
                    candidate_swap_node = candidate_swap_node->get_child(random_child_index);
                }
                swap_node_pool.push_back(candidate_swap_node);
            }

            // Randomly pick two nodes to swap
            std::vector<unsigned int> indices_of_nodes_to_swap = rng.random_subset_indices(
                    swap_node_pool.size(),
                    2);
            ECOEVOLITY_ASSERT(indices_of_nodes_to_swap.size() == 2);
            std::shared_ptr<NodeType> swap_node1 = swap_node_pool.at(indices_of_nodes_to_swap.at(0));
            std::shared_ptr<NodeType> swap_node2 = swap_node_pool.at(indices_of_nodes_to_swap.at(1));
            std::shared_ptr<NodeType> parent1 = swap_node1->get_parent();
            std::shared_ptr<NodeType> parent2 = swap_node2->get_parent();
            swap_node1->remove_parent();
            swap_node2->remove_parent();
            swap_node1->add_parent(parent2);
            swap_node2->add_parent(parent1);
            if (refresh_node_ordering) {
                this->refresh_ordered_nodes();
            }
            this->root_->resize_splits(this->get_leaf_node_count());
        }

        std::vector<unsigned int> get_intervening_height_indices(
                const unsigned int height_index,
                const double value) {
            std::vector<unsigned int> indices;
            // value is larger than this height
            if (this->get_height(height_index) < value) {
                if (height_index == (this->get_number_of_node_heights() - 1)) {
                    return indices;
                }
                for (unsigned int i = (height_index + 1);
                        i < this->get_number_of_node_heights();
                        ++i) {
                    if (this->get_height(i) > value) {
                        break;
                    }
                    indices.push_back(i);
                }
                return indices;
            }
            // value is less than this height
            if (height_index == 0) {
                return indices;
            }
            // decrementing an unsigned int to zero is a bit tricky; this stops
            // before we hit zero in the evaluation, but decrements before
            // jumping into loop body
            for (unsigned int i = height_index;
                    i-- > 0;
                ) {
                if (this->get_height(i) < value) {
                    break;
                }
                indices.push_back(i);
            }
            return indices;
        }

        std::vector<unsigned int> get_indices_of_intervening_nodes(
                const unsigned int height_index,
                const double value) {
            std::vector<unsigned int> indices;
            // value is larger than this height
            if (this->get_height(height_index) < value) {
                if (height_index == (this->get_number_of_node_heights() - 1)) {
                    return indices;
                }
                for (unsigned int i = (height_index + 1);
                        i < this->get_number_of_node_heights();
                        ++i) {
                    if (this->get_height(i) > value) {
                        break;
                    }
                    if (this->maps_ancestors_of(height_index, i)) {
                        indices.push_back(i);
                        continue;
                    }
                    for (auto idx : indices) {
                        if (this->maps_ancestors_of(idx, i)) {
                            indices.push_back(i);
                            break;
                        }
                    }
                }
                return indices;
            }
            // value is less than this height
            if (height_index == 0) {
                return indices;
            }
            // decrementing an unsigned int to zero is a bit tricky; this stops
            // before we hit zero in the evaluation, but decrements before
            // jumping into loop body
            for (unsigned int i = height_index;
                    i-- > 0;
                ) {
                if (this->get_height(i) < value) {
                    break;
                }
                if (this->maps_ancestors_of(i, height_index)) {
                    indices.push_back(i);
                    continue;
                }
                for (auto idx : indices) {
                    if (this->maps_ancestors_of(i, idx)) {
                        indices.push_back(i);
                        break;
                    }
                }
            }
            return indices;
        }

        unsigned int get_nearest_height_index(const double value) {
            if (this->get_number_of_node_heights() < 2) {
                return 0;
            }
            if (this->get_height(0) > value) {
                return 0;
            }
            if (this->get_height(this->get_number_of_node_heights() - 1) < value) {
                return this->get_number_of_node_heights() - 1;
            }
            double diff = fabs(this->get_height(0) - value);
            double current_diff = diff;
            for (unsigned int i = 1; i < this->node_heights_.size(); ++i) {
                current_diff = fabs(this->get_height(i) - value);
                if (current_diff > diff) {
                    return i - 1;
                }
                diff = current_diff;
            }
            return this->get_number_of_node_heights() - 1;
        }

        bool slide_bump_height(
                RandomNumberGenerator & rng,
                const unsigned int height_index,
                const double new_height) {
            return this->slide_bump_height_(rng, height_index, new_height, 0);
        }
        bool slide_bump_permute_height(
                RandomNumberGenerator & rng,
                const unsigned int height_index,
                const double new_height) {
            return this->slide_bump_height_(rng, height_index, new_height, 1);
        }
        bool slide_bump_swap_height(
                RandomNumberGenerator & rng,
                const unsigned int height_index,
                const double new_height) {
            return this->slide_bump_height_(rng, height_index, new_height, 2);
        }
        bool slide_bump_swap_all_height(
                RandomNumberGenerator & rng,
                const unsigned int height_index,
                const double new_height) {
            return this->slide_bump_height_(rng, height_index, new_height, 3);
        }

        virtual void set_root(std::shared_ptr<NodeType> root) {
            this->root_ = root;
            this->vet_tree();
            this->update_node_heights();
            this->refresh_ordered_nodes();
            this->root_->resize_splits(this->get_leaf_node_count());
            this->root_->make_all_dirty();
            this->make_dirty();
        }

        const NodeType& get_root() const {return *this->root_;}
        NodeType& get_mutable_root() const {return *this->root_;}

        std::shared_ptr<NodeType> get_root_ptr() const {return this->root_;}

        std::shared_ptr<NodeType> get_node(std::string label) const {
            std::shared_ptr<NodeType> n = this->root_->get_node(label);
            if (! n) {
                throw EcoevolityError(
                        "Node with label " + label + " not in tree");
            }
            return this->root_->get_node(label);
        }

        virtual void set_root_node_height_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            this->root_->set_node_height_prior(prior);
        }
        virtual std::shared_ptr<ContinuousProbabilityDistribution> get_root_node_height_prior() const {
            return this->root_->get_node_height_prior();
        }

        double get_height(const unsigned int height_index) const {
            return this->node_heights_.at(height_index)->get_value();
        }
        double get_stored_height(const unsigned int height_index) const {
            return this->node_heights_.at(height_index)->get_stored_value();
        }

        void set_height(const unsigned int height_index, double height) {
            double youngest_parent = std::numeric_limits<double>::infinity();
            if (height_index < (this->node_heights_.size() - 1)) {
                youngest_parent = this->get_height_of_youngest_parent(height_index);
            }
            double oldest_child = this->get_height_of_oldest_child(height_index);
            ECOEVOLITY_ASSERT(height < youngest_parent);
            ECOEVOLITY_ASSERT(height > oldest_child);
            this->node_heights_.at(height_index)->set_value(height);
            this->make_dirty();
            this->sort_node_heights();
        }

        double get_relative_height(const unsigned int height_index) const {
            return this->get_height(height_index) / this->get_root_height();
        }
        double get_height_relative_to_youngest_parent(
                const unsigned int height_index) const {
            return this->get_height(height_index) / this->get_height_of_youngest_parent(height_index);
        }

        std::shared_ptr<PositiveRealParameter> get_height_parameter(const unsigned int height_index) const {
            return this->node_heights_.at(height_index);
        }

        std::shared_ptr<NodeType> get_youngest_parent(const unsigned int height_index) const {
            if (height_index == (this->node_heights_.size() - 1)) {
                throw EcoevolityError("called get_height_index_of_youngest_parent with root index");
            }
            std::vector< std::shared_ptr<NodeType> > mapped_nodes = this->get_mapped_nodes(height_index);
            std::shared_ptr<NodeType> youngest_parent = mapped_nodes.at(0)->get_parent();
            for (unsigned int i = 1; i < mapped_nodes.size(); ++i) {
                if (mapped_nodes.at(i)->get_parent()->get_height() < youngest_parent->get_height()) {
                    youngest_parent = mapped_nodes.at(i)->get_parent();
                }
            }
            return youngest_parent;
        }

        double get_height_of_youngest_parent(const unsigned int height_index) const {
            return this->get_youngest_parent(height_index)->get_height();
        }

        double get_relative_height_of_youngest_parent(const unsigned int height_index) const {
            return this->get_height_of_youngest_parent(height_index) / this->get_root_height();
        }

        unsigned int get_index_of_youngest_parent(const unsigned int height_index) const {
            return this->get_node_height_index(this->get_youngest_parent(height_index)->get_height_parameter());
        }

        bool maps_ancestors_of(const unsigned int younger_index,
                const unsigned int older_index) {
            ECOEVOLITY_ASSERT(younger_index < older_index);
            std::vector< std::shared_ptr<NodeType> > older_nodes = this->get_mapped_nodes(older_index);
            std::vector< std::shared_ptr<NodeType> > younger_nodes = this->get_mapped_nodes(younger_index);
            for (auto young_node: younger_nodes) {
                for (auto old_node : older_nodes) {
                    if (young_node->is_ancestor(old_node)) {
                        return true;
                    }
                }
            }
            return false;
        }

        std::shared_ptr<NodeType> get_oldest_child(const unsigned int height_index) const {
            std::vector< std::shared_ptr<NodeType> > mapped_nodes = this->get_mapped_nodes(height_index);
            std::shared_ptr<NodeType> oldest_child = mapped_nodes.at(0)->get_oldest_child();
            for (unsigned int i = 1; i < mapped_nodes.size(); ++i) {
                if (mapped_nodes.at(i)->get_oldest_child()->get_height() > oldest_child->get_height()) {
                    oldest_child = mapped_nodes.at(i)->get_oldest_child();
                }
            }
            return oldest_child;
        }

        double get_height_of_oldest_child(const unsigned int height_index) const {
            return this->get_oldest_child(height_index)->get_height();
        }

        double get_relative_height_of_oldest_child(const unsigned int height_index) const {
            return this->get_height_of_oldest_child(height_index) / this->get_root_height();
        }

        void set_root_height(double height) {
            this->root_->set_height(height);
        }
        double get_root_height() const {
            return this->root_->get_height();
        }

        void scale_tree(double multiplier) {
            ECOEVOLITY_ASSERT(multiplier >= 0.0);
            ECOEVOLITY_ASSERT(! this->root_height_is_fixed());
            for (unsigned int i = 0; i < this->node_heights_.size(); ++i) {
                this->node_heights_.at(i)->set_value(
                        multiplier * this->node_heights_.at(i)->get_value());
            }
            this->make_dirty();
        }

        unsigned int get_degree_of_root() const {
            return this->root_->degree();
        }

        unsigned int get_leaf_node_count() const {
            return this->root_->get_leaf_node_count();
        }
        unsigned int get_internal_node_count() const {
            return this->root_->get_internal_node_count();
        }
        unsigned int get_node_count() const {
            return this->root_->get_node_count();
        }

        unsigned int get_number_of_node_heights() const {
            return this->node_heights_.size();
        }

        const std::vector< std::shared_ptr<PositiveRealParameter> >& get_node_height_pointers() const {
            return this->node_heights_;
        }

        std::vector<double> get_node_heights() const {
            std::vector<double> heights (this->node_heights_.size());
            for (unsigned int i = 0; i < this->node_heights_.size(); ++i) {
                heights.at(i) = this->node_heights_.at(i)->get_value();
            }
            return heights;
        }

        void ignore_data() {
            this->ignore_data_ = true;
        }
        void use_data() {
            this->ignore_data_ = false;
        }
        bool ignoring_data() const {
            return this->ignore_data_;
        }

        unsigned int get_number_of_likelihood_calculations() {
            return this->number_of_likelihood_calculations_;
        }

        virtual void compute_log_likelihood_and_prior(unsigned int nthreads = 1) {
            if (this->is_dirty()) {
                this->compute_log_likelihood(nthreads);
            }
            this->compute_log_prior_density();
            this->make_clean();
            return;
        }

        virtual double compute_log_likelihood(unsigned int nthreads = 1) {
            ++this->number_of_likelihood_calculations_;
            this->log_likelihood_.set_value(0.0);
            return 0.0;
        }

        void set_log_likelihood_value(double value) {
            this->log_likelihood_.set_value(value);
        }
        double get_log_likelihood_value() const {
            return this->log_likelihood_.get_value();
        }
        double get_stored_log_likelihood_value() const {
            return this->log_likelihood_.get_stored_value();
        }

        virtual double compute_log_prior_density() {
            double d = 0.0;
            d += this->compute_log_prior_density_of_parameters_of_beta_prior_on_node_heights();
            d += this->compute_log_prior_density_of_node_heights();
            d += this->compute_relative_log_prior_density_of_topology();
            d += this->get_derived_class_component_of_log_prior_density();
            this->log_prior_density_.set_value(d);
            return d;
        }
        virtual double get_derived_class_component_of_log_prior_density() const {
            return 0.0;
        }

        virtual double compute_log_prior_density_of_parameters_of_beta_prior_on_node_heights() {
            double d = 0.0;
            // prior_ln_pdf() method checks if parameter is fixed (and returns
            // 0 if so)
            d += this->alpha_of_node_height_beta_prior_->prior_ln_pdf();
            d += this->beta_of_node_height_beta_prior_->prior_ln_pdf();
            return d;
        }

        virtual double compute_log_prior_density_of_node_heights() const {
            // double root_height = this->root_->get_height();
            double d = 0.0;
            d += this->root_->get_height_relative_prior_ln_pdf();
            // I was using a Uniform(0, root_height) on all unique node
            // heights, but this caused a problem. Basically, heights mapped ot
            // nodes that descend from a non-root internal node were higher
            // order statistics that could not sample from their support
            // interval.
            // A Dirichlet distribution on unique heights doesn't fix this
            // either, because this causes problems with non-nested nodes
            // (e.g., a Dirichlet is not a good distribution for a balanced,
            // bifurcating, 4-tipped tree with independent node heights).
            // Kishino, Thorne, and Bruno (2001) had a solution for this, but
            // it's not easy to adapt to cases of shared node heights across
            // the tree.
            // The solution below is a bit of a hack, but works.
            // We assume that each unique node height is uniformly distributed
            // between zero and the height of the youngest parent of a node
            // mapped to the height.
            // This makes both nested and balanced tree shapes coherent.
            // This distribution is funky, but if we used scaled betas rather
            // than uniform, and put a hyper prior on the alpha shape parameter
            // of the betas, it would make it flexible.
            //
            // Old uniform(0, root_height) code:
            //     // prior prob density of non-root internal nodes = (1 / root_height)^n,
            //     // where n = the number of non-root internal nodes, so on
            //     // log scale = n * (log(1)-log(root_height)) = n * -log(root_height)
            //     double internal_node_height_prior_density = -std::log(root_height);
            //     d += internal_node_height_prior_density * (this->get_number_of_node_heights() - 1);
            // The conditional uniform solution (uniform(0, youngest parent)):
            for (unsigned int i = 0; i < this->get_number_of_node_heights() - 1; ++i) {
                ///////////////////////////////////////////////////////////////
                // Prior on the absolute ages of non-root internal nodes
                double youngest_parent_height = this->get_height_of_youngest_parent(i);
                // d -= std::log(youngest_parent_height);
                d += BetaDistribution::get_scaled_ln_pdf(this->get_height(i),
                        this->alpha_of_node_height_beta_prior_->get_value(),
                        this->beta_of_node_height_beta_prior_->get_value(),
                        youngest_parent_height);
                ///////////////////////////////////////////////////////////////
                // Prior on the relative ages of non-root internal nodes
                // double youngest_parent_rel_height = this->get_relative_height_of_youngest_parent(i);
                // d -= std::log(youngest_parent_rel_height);
            }
            return d;
        }

        virtual double compute_relative_log_prior_density_of_topology() const {
            return 0.0;
        }

        double get_log_prior_density_value() const {
            return this->log_prior_density_.get_value();
        }
        double get_stored_log_prior_density_value() const {
            return this->log_prior_density_.get_stored_value();
        }

        void store_state() {
            this->store_likelihood();
            this->store_prior_density();
            this->store_parameters();
        }
        void store_likelihood() {
            this->log_likelihood_.store();
        }
        void store_prior_density() {
            this->log_prior_density_.store();
        }
        virtual void store_parameters() {
            this->store_node_parameter_values();
            this->store_node_parameter_pointers();
            this->store_topology();
            this->alpha_of_node_height_beta_prior_->store();
            this->beta_of_node_height_beta_prior_->store();
            // Derived classes can override this
            this->store_derived_class_parameters();
        }
        virtual void store_all_heights() {
            this->root_->store_all_heights();
        }
        virtual void store_all_height_pointers() {
            this->root_->store_all_height_pointers();
        }
        virtual void store_node_parameter_values() {
            this->root_->store_all_parameter_values();
        }
        virtual void store_node_parameter_pointers() {
            this->root_->store_all_parameter_pointers();
        }
        virtual void store_topology() {
            this->stored_root_ = this->root_->get_copy();
        }
        // Derived classes can override this
        virtual void store_derived_class_parameters() { }
        void restore_state() {
            this->restore_likelihood();
            this->restore_prior_density();
            this->restore_parameters();
        }
        void restore_likelihood() {
            this->log_likelihood_.restore();
        }
        void restore_prior_density() {
            this->log_prior_density_.restore();
        }
        virtual void restore_parameters() {
            this->restore_topology();
            this->restore_node_parameter_pointers();
            this->restore_node_parameter_values();
            this->update_node_heights();
            this->alpha_of_node_height_beta_prior_->restore();
            this->beta_of_node_height_beta_prior_->restore();
            // Derived classes can override this
            this->restore_derived_class_parameters();
        }
        virtual void restore_all_heights() {
            this->root_->restore_all_heights();
        }
        virtual void restore_all_height_pointers() {
            this->root_->restore_all_height_pointers();
        }
        virtual void restore_node_parameter_values() {
            this->root_->restore_all_parameter_values();
        }
        virtual void restore_node_parameter_pointers() {
            this->root_->restore_all_parameter_pointers();
        }
        virtual void restore_topology() {
            this->root_ = this->stored_root_;
            this->refresh_ordered_nodes();
        }
        // Derived classes can override this
        virtual void restore_derived_class_parameters() { }

        std::string to_parentheses(
                const bool include_comments = true,
                const unsigned int precision = 12) const {
            if (! include_comments) {
                return this->root_->to_parentheses(precision);
            }
            return this->to_parentheses(this->root_, precision);
        }

        std::string to_parentheses(
                std::shared_ptr<NodeType> node,
                const unsigned int precision = 12) const {
            std::ostringstream s;
            s.precision(precision);
            if (node->is_leaf()) {
                s << node->get_label();
            }
            else {
                unsigned int child_idx = 0;
                s << "(";
                for (auto child_iter: node->get_all_children()) {
                    if (child_idx > 0) {
                        s << ",";
                    }
                    s << this->to_parentheses(child_iter, precision);
                    ++child_idx;
                }
                s << ")";
            }
            // annotate node with comment string
            s << this->get_comment_data_string(node)
              << ":" << node->get_length();
            return s.str();
        }

        std::string get_comment_data_string(
                std::shared_ptr<NodeType> node,
                const unsigned int precision = 12) const {
            std::ostringstream s;
            s.precision(precision);
            s << "[&";
            if (! node->is_leaf()) {
                unsigned int height_idx = this->get_node_height_index(node->get_height_parameter());
                s << "height_index="
                  << height_idx
                  << ",";
            }
            s << node->get_comment_data_string(precision);
            s << "]";
            return s.str();
        }

        virtual void write_state_log_header(std::ostream& out,
                const std::string& delimiter = "\t",
                const bool short_summary = false) const {
            throw EcoevolityError("write_state_log_header called from base BaseTree class");
        }
        virtual void log_state(std::ostream& out,
                const unsigned int generation_index,
                const std::string& delimiter = "\t",
                const bool short_summary = false) const {
            throw EcoevolityError("log_state called from base BaseTree class");
        }
        virtual void log_nexus_tree(std::ostream& out,
                const unsigned int generation_index,
                const bool include_comments = true,
                const unsigned int precision = 12) const {
            out << "    TREE gen" << generation_index << " = [&R] "
                << this->to_parentheses(include_comments, precision)
                << ";" << std::endl;
        }
        void write_nexus_taxa_block(std::ostream& out) const {
            std::vector<std::string> leaf_labels = this->root_->get_leaf_labels();
            out << "BEGIN TAXA;\n"
                << "    DIMENSIONS NTAX=" << leaf_labels.size() << ";\n"
                << "    TAXLABELS\n";
            for (auto label : leaf_labels) {
                out << "        " << label << "\n";
            }
            out << "    ;\n"
                << "END;" << std::endl;
        }

        virtual void draw_from_prior(RandomNumberGenerator& rng) {
            this->alpha_of_node_height_beta_prior_->set_value_from_prior(rng);
            this->beta_of_node_height_beta_prior_->set_value_from_prior(rng);
            unsigned int num_heights = this->get_number_of_node_heights();
            if (! this->root_height_is_fixed()) {
                this->node_heights_.at(num_heights - 1)->set_value(
                        this->get_root_node_height_prior()->draw(rng));
            }
            for (int i = (num_heights - 2); i >= 0; --i) {
                double youngest_parent_height = this->get_height_of_youngest_parent(i);
                double beta_draw = rng.beta(
                        this->alpha_of_node_height_beta_prior_->get_value(),
                        this->beta_of_node_height_beta_prior_->get_value());
                double height = youngest_parent_height * beta_draw;
                this->node_heights_.at(i)->set_value(height);
            }
            this->sort_node_heights();
            this->refresh_ordered_nodes();
        }

        void store_splits_by_height_index(
                std::map< unsigned int, std::set<Split> > & split_map,
                bool resize_splits = false) const {
            if (resize_splits) {
                this->root_->resize_splits(this->get_leaf_node_count());
            }
            for (auto node = this->pre_ordered_nodes_.rbegin();
                    node != this->pre_ordered_nodes_.rend();
                    ++node) {
                if (! (*node)->is_leaf()) { 
                    // add this internal node's split to split set
                    unsigned int height_idx = this->get_node_height_index((*node)->get_height_parameter());
                    split_map[height_idx].insert((*node)->split_);
                }
                else {
                    // Set bit for this leaf node's index
                    (*node)->split_.set_leaf_bit((*node)->get_index());
                }
                if ((*node)->has_parent()) {
                    (*node)->get_parent()->split_.add_split((*node)->split_);
                }
            }
        }

        std::map< unsigned int, std::set<Split> > get_splits_by_height_index(
                bool resize_splits = false) const {
            std::map< unsigned int, std::set<Split> > split_map;
            this->store_splits_by_height_index(split_map, resize_splits);
            return split_map;
        }

        void store_splits(
                std::set< std::set<Split> > & split_set,
                bool resize_splits = false) const {
            std::map< unsigned int, std::set<Split> > split_map = this->get_splits_by_height_index(resize_splits);
            for (auto item : split_map) {
                split_set.insert(item.second);
            }
        }

        std::set< std::set<Split> > get_splits(
                bool resize_splits = false) const {
            std::set< std::set<Split> > split_set;
            this->store_splits(split_set, resize_splits);
            return split_set;
        }

        void store_splits_heights_parameters(
                std::set< std::set<Split> > & split_set,
                std::map<std::set<Split>, double> heights,
                std::map<Split, std::map<std::string, double> > parameters,
                bool resize_splits = false) const {
            std::map< unsigned int, std::set<Split> > split_map;
            if (resize_splits) {
                this->root_->resize_splits(this->get_leaf_node_count());
            }
            for (auto node = this->pre_ordered_nodes_.rbegin();
                    node != this->pre_ordered_nodes_.rend();
                    ++node) {
                if (! (*node)->is_leaf()) {
                    // add this internal node's split to split set
                    unsigned int height_idx = this->get_node_height_index((*node)->get_height_parameter());
                    split_map[height_idx].insert((*node)->split_);
                    std::map<std::string, double> parameter_map;
                    (*node)->get_parameter_map(parameter_map);
                    if (parameters.size() > 0) {
                        parameters[(*node)->split_] = parameter_map;
                    }
                }
                else {
                    // Set bit for this leaf node's index
                    (*node)->split_.set_leaf_bit((*node)->get_index());
                }
                if ((*node)->has_parent()) {
                    (*node)->get_parent()->split_.add_split((*node)->split_);
                }
            }
            for (auto item : split_map) {
                split_set.insert(item.second);
                heights[item.second] = this->get_height(item.first);
            }
        }

        double get_alpha_of_node_height_beta_prior() const {
            return this->alpha_of_node_height_beta_prior_->get_value();
        }
        void set_alpha_of_node_height_beta_prior(double value) {
            this->alpha_of_node_height_beta_prior_->set_value(value);
        }
        bool alpha_of_node_height_beta_prior_is_fixed() const {
            return this->alpha_of_node_height_beta_prior_->is_fixed();
        }
        void fix_alpha_of_node_height_beta_prior() {
            this->alpha_of_node_height_beta_prior_->fix();
        }
        void estimate_alpha_of_node_height_beta_prior() {
            this->alpha_of_node_height_beta_prior_->estimate();
        }
        void set_prior_on_alpha_of_node_height_beta_prior(
                std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            this->alpha_of_node_height_beta_prior_->set_prior(prior);
        }
        std::shared_ptr<ContinuousProbabilityDistribution> get_prior_on_alpha_of_node_height_beta_prior() const {
            return this->alpha_of_node_height_beta_prior_->get_prior();
        }

        double get_beta_of_node_height_beta_prior() const {
            return this->beta_of_node_height_beta_prior_->get_value();
        }
        void set_beta_of_node_height_beta_prior(double value) {
            this->beta_of_node_height_beta_prior_->set_value(value);
        }
        bool beta_of_node_height_beta_prior_is_fixed() const {
            return this->beta_of_node_height_beta_prior_->is_fixed();
        }
        void fix_beta_of_node_height_beta_prior() {
            this->beta_of_node_height_beta_prior_->fix();
        }
        void estimate_beta_of_node_height_beta_prior() {
            this->beta_of_node_height_beta_prior_->estimate();
        }
        void set_prior_on_beta_of_node_height_beta_prior(
                std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            this->beta_of_node_height_beta_prior_->set_prior(prior);
        }
        std::shared_ptr<ContinuousProbabilityDistribution> get_prior_on_beta_of_node_height_beta_prior() const {
            return this->beta_of_node_height_beta_prior_->get_prior();
        }
};

#endif
