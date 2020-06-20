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

#ifndef ECOEVOLITY_TREESUM_HPP
#define ECOEVOLITY_TREESUM_HPP

#include <iostream>
#include <sstream>
#include <cmath>

#include "assert.hpp"
#include "stats_util.hpp"
#include "basetree.hpp"
#include "treecomp.hpp"


namespace treesum {

class BaseSamples {
    protected:
        unsigned int n_ = 0;
        std::vector<unsigned int> source_indices_;
        std::vector<unsigned int> tree_indices_;

        void tally_sample_(
                unsigned int tree_index,
                unsigned int source_index) {
            this->tree_indices_.push_back(tree_index);
            this->source_indices_.push_back(source_index);
            ++this->n_;
        }

    public:
        bool operator< (const BaseSamples & other) const {
            return this->n_ < other.n_;
        }
        bool operator> (const BaseSamples & other) const {
            return this->n_ > other.n_;
        }
        static bool sort_by_n(
                const std::shared_ptr<BaseSamples> & s1,
                const std::shared_ptr<BaseSamples> & s2) {
            return s1->n_ < s2->n_;
        }
        static bool reverse_sort_by_n(
                const std::shared_ptr<BaseSamples> & s1,
                const std::shared_ptr<BaseSamples> & s2) {
            return s1->n_ > s2->n_;
        }

        unsigned int get_sample_size() const {
            return this->n_;
        }

        const std::vector<unsigned int> & get_source_indices() const {
            return this->source_indices_;
        }

        const std::vector<unsigned int> & get_tree_indices() const {
            return this->tree_indices_;
        }
};

class TopologySamples : public BaseSamples {
    protected:
        std::set< std::set<Split> > split_set_;
        std::map<std::set<Split>, std::vector<double> > heights_;

    public:
        void add_sample(
                const std::map<std::set<Split>, double> & height_map,
                unsigned int tree_index,
                unsigned int source_index) { 
            std::set< std::set<Split> > s_set;
            for (auto splits_height : height_map) {
                s_set.insert(splits_height.first);
                this->heights_[splits_height.first].push_back(
                        splits_height.second);
            }
            if (this->n_ == 0) {
                this->split_set_ = s_set;
            }
            else {
                ECOEVOLITY_ASSERT(s_set == this->split_set_);
            }
            this->tally_sample_(tree_index, source_index);
        }
};

class HeightSamples : public BaseSamples {
    protected:
        std::set<Split> split_set_;
        std::vector<double> heights_;

    public:
        void add_sample(
                const std::set<Split> & set_of_splits,
                double height,
                unsigned int tree_index,
                unsigned int source_index) {
            if (this->n_ == 0) {
                this->split_set_ = set_of_splits;
            }
            else {
                ECOEVOLITY_ASSERT(set_of_splits == this->split_set_);
            }
            this->heights_.push_back(height);
            this->tally_sample_(tree_index, source_index);
        }
};

class SplitSamples : public BaseSamples {
    protected:
        Split split_;
        std::map<std::string, std::vector<double> > parameters_;

    public:
        void add_sample(
                const Split & s,
                const std::map<std::string, double> & p,
                unsigned int tree_index,
                unsigned int source_index) {
            if (this->n_ == 0) {
                this->split_ = s;
            }
            else {
                ECOEVOLITY_ASSERT(s == this->split_);
            }
            for (auto pname_value : p) {
                this->parameters_[pname_value.first].push_back(pname_value.second);
            }
            this->tally_sample_(tree_index, source_index);
        }
};

template<class NodeType>
class TreeSample {
    public:
        typedef BaseTree<NodeType> tree_type;

    protected:
        std::vector< std::shared_ptr<TopologySamples> > topologies_;
        std::vector< std::shared_ptr<HeightSamples> > heights_;
        std::vector< std::shared_ptr<SplitSamples> > splits_;
        std::map< std::set< std::set<Split> >, std::shared_ptr<TopologySamples> > topologies_map_;
        std::map< std::set<Split>,             std::shared_ptr<HeightSamples>   > heights_map_;
        std::map< Split,                       std::shared_ptr<SplitSamples>    > splits_map_;
        std::vector<double> tree_lengths_;
        std::vector<std::string> source_paths_;
        std::vector<unsigned int> source_sample_sizes_;
        Split root_split_;
        std::vector<Split> leaf_splits_;
        std::vector<std::string> leaf_labels_;
        unsigned int sample_size_ = 0;
        bool target_tree_provided_ = false;
        tree_type target_tree_;
        std::vector<double> target_euclidean_distances_;

        void reverse_sort_samples_by_freq_() {
            std::sort(this->topologies_.begin(), this->topologies_.end(),
                    BaseSamples::reverse_sort_by_n);
            std::sort(this->heights_.begin(), this->heights_.end(),
                    BaseSamples::reverse_sort_by_n);
            std::sort(this->splits_.begin(), this->splits_.end(),
                    BaseSamples::reverse_sort_by_n);
        }

        void update_splits_and_labels_(const tree_type & tree) {
            unsigned int nleaves = tree.get_leaf_node_count();
            ECOEVOLITY_ASSERT(nleaves > 0);
            tree.get_leaf_labels(this->leaf_labels_);
            ECOEVOLITY_ASSERT(this->leaf_labels_.size() == nleaves);
            // Leaf labels are sorted for every tree is parsed by
            // BaseTree::get_trees.  Sorting here to make sure order is
            // consistent with sampled trees (parsed by BaseTree)
            std::sort(std::begin(this->leaf_labels_), std::end(this->leaf_labels_));
            std::set<std::string> leaf_label_set(
                    std::begin(this->leaf_labels_),
                    std::end(this->leaf_labels_));
            ECOEVOLITY_ASSERT(leaf_label_set.size() == nleaves);
            this->root_split_.resize(nleaves);
            for (unsigned int i = 0; i < nleaves; ++i) {
                this->root_split_.set_leaf_bit(i);
                Split leaf_split;
                leaf_split.resize(nleaves);
                leaf_split.set_leaf_bit(i);
                // Because leaf_labels_ are sorted, the leaf_splits_ will be in
                // order to match leaf_labels_
                this->leaf_splits_.push_back(leaf_split);
            }
        }

    public:

        TreeSample() { }
        TreeSample(
                const std::vector<std::string> & paths,
                const std::string & ncl_file_format,
                unsigned int skip = 0,
                double ultrametricity_tolerance = 1e-6) {
            for (auto path : paths) {
                this->add_trees(path, ncl_file_format, skip,
                        ultrametricity_tolerance);
            }
        }
        TreeSample(
                const std::string & target_tree_path,
                const std::vector<std::string> & paths,
                const std::string & ncl_file_format,
                unsigned int skip = 0,
                double ultrametricity_tolerance = 1e-6) {
            this->set_target_tree(target_tree_path, ncl_file_format);
            for (auto path : paths) {
                this->add_trees(path, ncl_file_format, skip,
                        ultrametricity_tolerance);
            }
        }

        void add_trees(
                const std::string & path,
                const std::string & ncl_file_format,
                unsigned int skip = 0,
                double ultrametricity_tolerance = 1e-6) {
            this->source_paths_.push_back(path);
            std::vector<tree_type> trees;
            get_trees<tree_type>(
                    path,
                    ncl_file_format,
                    trees,
                    skip,
                    ultrametricity_tolerance);
            this->add_trees(trees);
        }

        void add_trees(
                std::istream & tree_stream,
                const std::string & ncl_file_format,
                unsigned int skip = 0,
                double ultrametricity_tolerance = 1e-6) {
            std::vector<tree_type> trees;
            get_trees<tree_type>(
                    tree_stream,
                    ncl_file_format,
                    trees,
                    skip,
                    ultrametricity_tolerance);
            this->add_trees(trees);
        }

        void add_trees(
                const std::vector<tree_type> & trees) {
            this->source_sample_sizes_.push_back(0);
            std::set< std::set<Split> > split_set;
            std::map<std::set<Split>, double> heights;
            std::map<Split, std::map<std::string, double> > parameters;
            std::map<std::string, double> parameter_map;
            std::vector< std::shared_ptr<NodeType> > leaves;
            unsigned int source_index = this->source_sample_sizes_.size();
            for (auto tree : trees) {
                if (this->sample_size_ == 0) {
                    this->update_splits_and_labels_(tree);
                }
                unsigned int tree_index = this->sample_size_;
                split_set.clear();
                heights.clear();
                parameters.clear();
                tree.store_splits_heights_parameters(
                        split_set,
                        heights,
                        parameters,
                        false);
                if (this->topologies_map_.count(split_set) > 0) {
                    this->topologies_map_[split_set]->add_sample(
                            heights,
                            tree_index,
                            source_index);
                }
                else {
                    std::shared_ptr<TopologySamples> ts;
                    ts->add_sample(heights, tree_index, source_index);
                    this->topologies_.push_back(ts);
                    this->topologies_map_[split_set] = ts;
                }
                for (auto splits_height : heights) {
                    if (this->heights_map_.count(splits_height.first) > 0) {
                        this->heights[splits_height.first]->add_sample(
                                splits_height.first,
                                splits_height.second,
                                tree_index,
                                source_index);
                    }
                    else {
                        std::shared_ptr<HeightSamples> hs;
                        hs->add_sample(
                                splits_height.first,
                                splits_height.second,
                                tree_index,
                                source_index);
                        this->heights_.push_back(hs);
                        this->heights_map_[splits_height.first] = hs;
                    }
                }
                for (auto split_pmap : parameters) {
                    if (this->splits_map_.count(split_pmap.first) > 0) {
                        this->splits[split_pmap.first]->add_sample(
                                split_pmap.first,
                                split_pmap.second,
                                tree_index,
                                source_index);
                    }
                    else {
                        std::shared_ptr<SplitSamples> ss;
                        ss->add_sample(
                                split_pmap.first,
                                split_pmap.second,
                                tree_index,
                                source_index);
                        this->splits_.push_back(ss);
                        this->splits_map_[split_pmap.first] = ss;
                    }
                }
                if (this->target_tree_provided_) {
                    this->target_euclidean_distances_.push_back(
                            treecomp::euclidean_distance<tree_type>(
                                this->target_tree_,
                                tree,
                                false)
                            );
                }
                this->tree_lengths_.push_back(tree.get_tree_length());
                ++this->sample_size_;
                ++this->source_sample_sizes_.back();
            }
            unsigned int source_total = 0;
            for (unsigned int n : this->source_sample_sizes_) {
                source_total += n;
            }
            ECOEVOLITY_ASSERT(source_total == this->sample_size_);
            this->reverse_sort_samples_by_freq_();
        }

        void set_target_tree(
                std::istream & tree_stream,
                const std::string & ncl_file_format) {
            this->target_tree_ = tree_type(tree_stream, ncl_file_format);
            this->target_tree_provided_ = true;
        }

        void set_target_tree(
                const std::string & tree_path,
                const std::string & ncl_file_format) {
            this->target_tree_ = tree_type(tree_path, ncl_file_format);
            this->target_tree_provided_ = true;
        }

        void set_target_tree(
                const std::string & newick_tree_string) {
            this->target_tree_ = tree_type(newick_tree_string);
            this->target_tree_provided_ = true;
        }

        unsigned int get_number_of_sources() const {
            return this->source_sample_sizes_.size();
        }

        unsigned int get_number_of_leaves() const {
            return this->leaf_splits_.size();
        }

        double get_average_std_dev_of_split_freqs(
                double credible_set_cutoff = 0.9) const {
            if (this->get_number_of_sources() < 2) {
                throw EcoevolityError("Calculating the ASDSF requires multiple chains");
            }
            double cumulative_freq = 0.0;
            SampleSummarizer<double> std_devs_of_split_freqs;
            for (auto ss : this->splits_) {
                if (cumulative_freq > credible_set_cutoff) {
                    break;
                }
                cumulative_freq += ss->get_sample_size() / (double)this->sample_size_;
                std::vector<unsigned int> split_counts(this->get_number_of_sources(), 0);
                SampleSummarizer<double> split_freqs;
                for (auto source_idx : ss->get_source_indices()) {
                    ++split_counts.at(source_idx);
                }
                unsigned int total = 0;
                for (unsigned int source_idx = 0;
                        source_idx < this->get_number_of_sources();
                        ++ source_idx) {
                    split_freqs.add_sample(split_counts.at(source_idx) / (double)this->sample_size_);
                    total += split_counts.at(source_idx);
                }
                ECOEVOLITY_ASSERT(total == this->sample_size_);
                std_devs_of_split_freqs.add_sample(split_freqs.std_dev());
            }
            return std_devs_of_split_freqs.mean();
        }
};

} // treesum

#endif
