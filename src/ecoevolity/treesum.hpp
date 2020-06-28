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
        std::map< std::set<Split>, std::vector<double> > heights_;

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

        const std::set< std::set<Split> > & get_split_set() const {
            return this->split_set_;
        }

        const std::map< std::set<Split>, std::vector<double> > & get_height_map() const {
            return this->heights_;
        }

        const std::vector<double> & get_heights(const std::set<Split> & split_set) const {
            return this->heights_.at(split_set);
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

        const std::set<Split> & get_split_set() const {
            return this->split_set_;
        }

        const std::vector<double> & get_heights() const {
            return this->heights_;
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

        const Split & get_split() const {
            return this->split_;
        }
        std::string get_split_as_string(
                const char unset_char = '0',
                const char set_char = '1') const {
            return this->split_.as_string(unset_char, set_char);
        }

        const std::map<std::string, std::vector<double> > & get_parameter_map() const {
            return this->parameters_;
        }

        const std::vector<double> & get_values(
                const std::string & parameter_name) const {
            return this->parameters_.at(parameter_name);
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
        std::vector< std::shared_ptr<SplitSamples> > non_trivial_splits_;
        std::map< std::set< std::set<Split> >, std::shared_ptr<TopologySamples> > topologies_map_;
        std::map< std::set<Split>,             std::shared_ptr<HeightSamples>   > heights_map_;
        std::map< Split,                       std::shared_ptr<SplitSamples>    > splits_map_;
        std::vector<double> tree_lengths_;
        std::vector<std::string> source_paths_;
        std::vector<unsigned int> source_sample_sizes_;
        std::set<Split> trivial_splits_;
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
            std::sort(this->non_trivial_splits_.begin(), this->non_trivial_splits_.end(),
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
                this->trivial_splits_.insert(leaf_split);
            }
            this->trivial_splits_.insert(this->root_split_);
        }

        void check_leaf_labels_(const tree_type & tree) {
            std::vector<std::string> l = tree.get_leaf_labels();
            if (
                    (l.size() != this->leaf_labels_.size()) ||
                    (! std::is_permutation(std::begin(l), std::end(l), std::begin(this->leaf_labels_)))
               ){
                throw EcoevolityError("Tip labels in trees do not match");
            }
        }

        void check_target_tree_() {
            if (this->sample_size_ < 1) {
                this->update_splits_and_labels_(this->target_tree_);
            }
            else {
                this->check_leaf_labels_(this->target_tree_);
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
            std::set< std::set<Split> > split_set;
            std::map<std::set<Split>, double> heights;
            std::map<Split, std::map<std::string, double> > parameters;
            std::map<std::string, double> parameter_map;
            std::vector< std::shared_ptr<NodeType> > leaves;
            unsigned int source_index = this->source_sample_sizes_.size();
            this->source_sample_sizes_.push_back(0);
            for (auto tree : trees) {
                if ((this->sample_size_ < 1) && (this->leaf_labels_.size() < 1)) {
                    this->update_splits_and_labels_(tree);
                }
                /* else { */
                /*     this->check_leaf_labels_(tree); */
                /* } */
                this->check_leaf_labels_(tree);

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
                    std::shared_ptr<TopologySamples> ts = std::make_shared<TopologySamples>();
                    ts->add_sample(heights, tree_index, source_index);
                    this->topologies_.push_back(ts);
                    this->topologies_map_[split_set] = ts;
                }
                for (auto splits_height : heights) {
                    if (this->heights_map_.count(splits_height.first) > 0) {
                        this->heights_map_[splits_height.first]->add_sample(
                                splits_height.first,
                                splits_height.second,
                                tree_index,
                                source_index);
                    }
                    else {
                        std::shared_ptr<HeightSamples> hs = std::make_shared<HeightSamples>();
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
                        this->splits_map_[split_pmap.first]->add_sample(
                                split_pmap.first,
                                split_pmap.second,
                                tree_index,
                                source_index);
                    }
                    else {
                        std::shared_ptr<SplitSamples> ss = std::make_shared<SplitSamples>();
                        ss->add_sample(
                                split_pmap.first,
                                split_pmap.second,
                                tree_index,
                                source_index);
                        this->splits_.push_back(ss);
                        this->splits_map_[split_pmap.first] = ss;
                        if (this->trivial_splits_.count(ss->get_split()) < 1) {
                            this->non_trivial_splits_.push_back(ss);
                        }
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
            this->check_target_tree_();
        }

        void set_target_tree(
                const std::string & tree_path,
                const std::string & ncl_file_format) {
            this->target_tree_ = tree_type(tree_path, ncl_file_format);
            this->target_tree_provided_ = true;
            this->check_target_tree_();
        }

        void set_target_tree(
                const std::string & newick_tree_string) {
            this->target_tree_ = tree_type(newick_tree_string);
            this->target_tree_provided_ = true;
            this->check_target_tree_();
        }

        unsigned int get_number_of_sources() const {
            return this->source_sample_sizes_.size();
        }

        unsigned int get_number_of_leaves() const {
            return this->leaf_splits_.size();
        }

        unsigned int get_sample_size() const {
            return this->sample_size_;
        }

        const std::vector<unsigned int> & get_source_sample_sizes() const {
            return this->source_sample_sizes_;
        }

        unsigned int get_source_sample_size(unsigned int source_index) const {
            return this->source_sample_sizes_.at(source_index);
        }

        const std::vector< std::shared_ptr<TopologySamples> > & get_topologies() const {
            return this->topologies_;
        }
        const std::vector< std::shared_ptr<HeightSamples> > & get_heights() const {
            return this->heights_;
        }
        const std::vector< std::shared_ptr<SplitSamples> > & get_splits() const {
            return this->splits_;
        }
        const std::vector< std::shared_ptr<SplitSamples> > & get_non_trivial_splits() const {
            return this->non_trivial_splits_;
        }

        std::shared_ptr<TopologySamples> get_topology(
                const std::set< std::set<Split> > & topology
                ) const {
            if (this->topologies_map_.count(topology) < 1) {
                return nullptr;
            }
            return this->topologies_map_.at(topology);
        }
        std::shared_ptr<HeightSamples> get_height(
                const std::set<Split> & split_set 
                ) const {
            if (this->heights_map_.count(split_set) < 1) {
                return nullptr;
            }
            return this->heights_map_.at(split_set);
        }
        std::shared_ptr<SplitSamples> get_split(
                const Split & split
                ) const {
            if (this->splits_map_.count(split) < 1) {
                return nullptr;
            }
            return this->splits_map_.at(split);
        }

        unsigned int get_topology_count(
                const std::set< std::set<Split> > & topology
                ) const {
            if (this->topologies_map_.count(topology) < 1) {
                return 0;
            }
            return this->topologies_map_.at(topology)->get_sample_size();
        }
        double get_topology_frequency(
                const std::set< std::set<Split> > & topology
                ) const {
            return this->get_topology_count(topology) / (double)this->sample_size_;
        }

        unsigned int get_height_count(
                const std::set<Split> & split_set
                ) const {
            if (this->heights_map_.count(split_set) < 1) {
                return 0;
            }
            return this->heights_map_.at(split_set)->get_sample_size();
        }
        double get_height_frequency(
                const std::set<Split> & split_set
                ) const {
            return this->get_height_count(split_set) / (double)this->sample_size_;
        }

        unsigned int get_split_count(
                const Split & split
                ) const {
            if (this->splits_map_.count(split) < 1) {
                return 0;
            }
            return this->splits_map_.at(split)->get_sample_size();
        }
        double get_split_frequency(
                const Split & split
                ) const {
            return this->get_split_count(split) / (double)this->sample_size_;
        }

        double get_average_std_dev_of_split_freqs(
                double min_frequency = 0.1) const {
            if (this->get_number_of_sources() < 2) {
                throw EcoevolityError("Calculating the ASDSF requires multiple chains");
            }
            SampleSummarizer<double> std_devs_of_split_freqs;
            for (auto ss : this->non_trivial_splits_) {
                if (this->get_split_frequency(ss->get_split()) < min_frequency) {
                    break;
                }
                std::vector<unsigned int> split_counts(this->get_number_of_sources(), 0);
                SampleSummarizer<double> split_freqs;
                for (auto source_idx : ss->get_source_indices()) {
                    ++split_counts.at(source_idx);
                }
                for (unsigned int source_idx = 0;
                        source_idx < this->get_number_of_sources();
                        ++source_idx) {
                    split_freqs.add_sample(
                            split_counts.at(source_idx) /
                            (double)this->get_source_sample_size(source_idx));
                }
                std_devs_of_split_freqs.add_sample(split_freqs.std_dev());
            }
            return std_devs_of_split_freqs.mean();
        }
};

} // treesum

#endif
