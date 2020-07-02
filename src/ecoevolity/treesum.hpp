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

        template <typename T>
        void _write_summary_of_values(
                const std::vector<T> & values,
                std::ostream & out,
                const std::string & parameter_name = "",
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            out.precision(precision);
            SampleSummary<T> summary(values);
            std::string p_name = parameter_name;
            if (parameter_name != "") {
                p_name += "_";
            }
            out << margin << p_name << "mean: " << summary.mean() << "\n"
                << margin << p_name << "median: " << summary.median() << "\n"
                << margin << p_name << "std_dev: " << summary.std_dev() << "\n"
                << margin << p_name << "range: ["
                                    << summary.min() << ", "
                                    << summary.max() << "]\n"
                << margin << p_name << "eti_95: ["
                                    << summary.qi_95().first << ", "
                                    << summary.qi_95().second << "]\n"
                << margin << p_name << "hpdi_95: ["
                                    << summary.hpdi_95().first << ", "
                                    << summary.hpdi_95().second << "]"
                << std::endl;
        }

        template <typename T>
        void _write_node_annotations(
                const std::vector<T> & values,
                std::ostream & out,
                const std::string & parameter_name,
                const unsigned int precision = 12) const {
            out.precision(precision);
            std::string p_name = parameter_name;
            if (parameter_name != "") {
                p_name += "_";
            }
            if (values.size() < 1) {
                double nan = std::numeric_limits<double>::quiet_NaN();
                out << p_name << "mean=" << nan << ","
                    << p_name << "median=" << nan << ","
                    << p_name << "std_dev=" << nan << ","
                    << p_name << "range={"
                              << nan << ","
                              << nan << "},"
                    << p_name << "eti_95={"
                              << nan << ","
                              << nan << "}\n"
                    << p_name << "hpdi_95={"
                              << nan << ","
                              << nan << "}";
                return;
            }
            SampleSummary<T> summary(values);
            out << p_name << "mean=" << summary.mean() << ","
                << p_name << "median=" << summary.median() << ","
                << p_name << "std_dev=" << summary.std_dev() << ","
                << p_name << "range={"
                          << summary.min() << ","
                          << summary.max() << "},"
                << p_name << "eti_95={"
                          << summary.qi_95().first << ","
                          << summary.qi_95().second << "}\n"
                << p_name << "hpdi_95={"
                          << summary.hpdi_95().first << ","
                          << summary.hpdi_95().second << "}";
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

        static bool sort_height_keys(
                const std::set<Split> & split_set1,
                const std::set<Split> & split_set2) {
            bool s1_in_s2 = false;
            bool s2_in s1 = false;
            for (auto s1 : split_set1) {
                for (auto s2 : split_set2) {
                    if (s1.is_subset_of(s2)) {
                        s1_in_s2 = true;
                    }
                    if (s2.is_subset_of(s1)) {
                        s2_in_s1 = true;
                    }
                }
            }
            ECOEVOLITY_ASSERT(! (s1_in_s2 && s2_in_s1));
            if (s1_in_s2) {
                return true;
            }
            if (s2_in_s1) {
                return false;
            }
            // They are not nested, so we need to use heights
            if (this->heights_map_.count(split_set1) < 1) {
                return true;
            }
            if (this->heights_map_.count(split_set2) < 1) {
                return false;
            }
            SampleSummarizer<double> s1_sum();
            SampleSummarizer<double> s2_sum();
            for (auto h : this->get_height(split_set1)->get_heights()) {
                s1_sum.add_sample(h);
            }
            for (auto h : this->get_height(split_set2)->get_heights()) {
                s2_sum.add_sample(h);
            }
            return s1_sum.mean() < s2_sum.mean();
        }

        static bool reverse_sort_height_keys(
                const std::set<Split> & split_set1,
                const std::set<Split> & split_set2) {
            bool s1_in_s2 = false;
            bool s2_in s1 = false;
            for (auto s1 : split_set1) {
                for (auto s2 : split_set2) {
                    if (s1.is_subset_of(s2)) {
                        s1_in_s2 = true;
                    }
                    if (s2.is_subset_of(s1)) {
                        s2_in_s1 = true;
                    }
                }
            }
            ECOEVOLITY_ASSERT(! (s1_in_s2 && s2_in_s1));
            if (s1_in_s2) {
                return false;
            }
            if (s2_in_s1) {
                return true;
            }
            // They are not nested, so we need to use heights
            if (this->heights_map_.count(split_set1) < 1) {
                return false;
            }
            if (this->heights_map_.count(split_set2) < 1) {
                return true;
            }
            SampleSummarizer<double> s1_sum();
            SampleSummarizer<double> s2_sum();
            for (auto h : this->get_height(split_set1)->get_heights()) {
                s1_sum.add_sample(h);
            }
            for (auto h : this->get_height(split_set2)->get_heights()) {
                s2_sum.add_sample(h);
            }
            return s1_sum.mean() > s2_sum.mean();
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

        void write_summary_of_splits(std::ostream & out,
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            std::string indent = string_util::get_indent(1);
            out.precision(precision);

            std::string item_margin = margin + indent + indent;

            out << margin << "splits:\n";
            out << margin << indent << "root_split:\n";
            this->write_summary_of_split(this->root_split_,
                    out,
                    false,
                    item_margin,
                    precision);
            item_margin += "  ";
            out << margin << indent << "leaf_splits:\n";
            for (auto s : this->leaf_splits_) {
                out << margin << indent << indent << "-\n";
                this->write_summary_of_split(s,
                        out,
                        true,
                        item_margin,
                        precision);
            }
            out << margin << indent << "nontrivial_splits:\n";
            for (auto split_samples : this->get_non_trivial_splits()) {
                out << margin << indent << indent << "-\n";
                this->write_summary_of_split(split_samples->get_split(),
                        out,
                        true,
                        item_margin,
                        precision);
            }
        }

        void write_summary_of_split(const Split & split,
                std::ostream & out,
                const bool include_leaf_indices = true,
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            out.precision(precision);
            if (include_leaf_indices) {
                std::vector<unsigned int> leaf_indices = split.get_leaf_indices();
                out << margin << "leaf_indices: [" << leaf_indices.at(0);
                for (unsigned int i = 1; i < leaf_indices.size(); ++i) {
                    out << ", " << leaf_indices.at(i);
                }
                out << "]\n";
            }
            out << margin << "count: " << this->get_split_count(split) << "\n"
                << margin << "frequency: " << this->get_split_frequency(split) << "\n";
            for (auto param_vals : this->get_split(split)->get_parameter_map()) {
                this->_write_summary_of_values<double>(
                        param_vals.second,
                        out,
                        param_vals.first,
                        margin,
                        precision);
            }
        }

        void write_summary_of_heights(std::ostream & out,
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            std::string indent = string_util::get_indent(1);
            out.precision(precision);

            std::set<Split> root_split_set;
            root_split_set.insert(this->root_split_);
            out << margin << "heights:\n";
            for (auto height_samples : this->get_heights()) {
                if (height_samples->get_split_set() == root_split_set) {
                    continue;
                }
                out << margin << indent << "-\n";
                this->write_summary_of_height(height_samples->get_split_set(),
                        out,
                        margin + indent + "  ",
                        precision);
            }
        }

        void write_summary_of_height(const std::set<Split> & split_set,
                std::ostream & out,
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            std::string indent = string_util::get_indent(1);
            out.precision(precision);

            out << margin << "number_of_nodes: " << split_set.size() << "\n"
                << margin << "splits:\n";
            for (auto split : split_set) {
                std::vector<unsigned int> leaf_indices = split.get_leaf_indices();
                out << margin << indent << "- leaf_indices: [" << leaf_indices.at(0);
                for (unsigned int i = 1; i < leaf_indices.size(); ++i) {
                    out << ", " << leaf_indices.at(i);
                }
                out << "]\n";
            }
            out << margin << "count: " << this->get_height_count(split_set) << "\n"
                << margin << "frequency: " << this->get_height_frequency(split_set) << "\n";
            this->_write_summary_of_values<double>(
                    this->get_height(split_set)->get_heights(),
                    out,
                    "",
                    margin,
                    precision);
        }

        void write_summary_of_topologies(std::ostream & out,
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            std::string indent = string_util::get_indent(1);
            out.precision(precision);

            double cumulative_freq = 0.0;
            std::string topo_margin = margin + "  ";
            out << margin << "topologies:\n";
            for (auto height_samples : this->get_heights()) {
                out << margin << indent << "-\n";
                cumulative_freq += this->write_summary_of_topology(
                        height_samples->get_split_set(),
                        out,
                        cumulative_freq,
                        topo_margin,
                        precision);
            }
        }

        double write_summary_of_topology(const std::set< std::set<Split> > & split_set,
                std::ostream & out,
                const double cumulative_freq = 0.0,
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            std::string indent = string_util::get_indent(1);
            out.precision(precision);
            out << margin << "number_of_heights: " << split_set.size() << "\n"
                << margin << "heights:\n";
            std::string h_margin = margin + indent + "  ";
            for (auto s_set : split_set) {
                out << margin << indent << "- number_of_nodes: " << s_set.size() << "\n";
                out << h_margin << "splits:\n";
                for (auto split : s_set) {
                    std::vector<unsigned int> leaf_indices = split.get_leaf_indices();
                    out << h_margin << indent << "- leaf_indices: [" << leaf_indices.at(0);
                    for (unsigned int i = 1; i < leaf_indices.size(); ++i) {
                        out << ", " << leaf_indices.at(i);
                    }
                    out << "]\n";
                }
                this->_write_summary_of_values<double>(
                        this->get_topology(split_set)->get_heights(s_set),
                        out,
                        "",
                        h_margin,
                        precision);
            }
            double freq = this->get_topology_frequency(split_set);
            out << margin << "count: " << this->get_topology_count(split_set) << "\n"
                << margin << "frequency: " << freq << "\n"
                << margin << "cumulative_frequency: " << freq + cumulative_freq
                << std::endl;
            return freq + cumulative_freq;
        }

        std::string to_parentheses(const std::set< std::set<Split> > & split_set,
                const bool use_median_heights = false,
                const unsigned int precision = 12) const {
            std::vector<std::string> parameter_keys;
            const std::map<std::string, std::vector<double> > & split_parameter_map;
            split_parameter_map = this->get_split(this->root_split_)->get_parameter_map();
            std::ostringstream root_label;
            root_label << "[&height_index=" << split_set.size() - 1 << ",";
            bool first_pass = true
            for (auto p_values : split_parameter_map) {
                if (! first_pass) {
                    root_label << ","
                }
                first_pass = false;
                this->_write_node_annotations<double>(p_values.second,
                        root_label,
                        p_values.first,
                        precision);
                parameter_keys.push_back(p_values.first);
            }
            SampleSummary<double> root_height_summary(split_parameter_map["height"]);
            double root_height = root_height_summary.mean();
            if (use_median_heights) {
                root_height = root_height_summary.median();
            }
            std::shared_ptr<Node> root_node = std::make_shared<Node>(
                    root_label.str(), root_height);

            std::vector< std::shared_ptr<Node> > leaf_nodes;
            for (unsigned int i = 0; i < this->leaf_splits_.size(); ++i) {
                std::ostringstream label;
                label << this->leaf_labels_.at(i) << "[&"
                unsigned int leaf_index = this->leaf_splits_.at(i).get_leaf_indices().at(0);
                ECOEVOLITY_ASSERT(leaf_index == i);
                split_parameter_map = this->get_split(this->leaf_splits_.at(i))->get_parameter_map();
                SampleSummary<double> leaf_height_summary(split_parameter_map["height"]);
                double leaf_height = leaf_height_summary.mean();
                if (use_median_heights) {
                    leaf_height = leaf_height_summary.median();
                }
                first_pass = true
                for (auto p_values : split_parameter_map) {
                    if (! first_pass) {
                        label << ","
                    }
                    first_pass = false;
                    this->_write_node_annotations<double>(p_values.second,
                            label,
                            p_values.first,
                            precision);
                }
                std::shared_ptr leaf_nd = std::make_shared<Node>(
                            leaf_index, label.str(), height);
                leaf_nodes.push_back(leaf_nd);
                root_node->add_child(leaf_nd);
            }

            std::vector< std::set<Split> > height_keys(std::begin(split_set),
                    std::end(split_set));
            std::sort(std::begin(height_keys), std::end(height_keys),
                    TreeSample::reverse_sort_height_keys);
            ECOEVOLITY_ASSERT((height_keys.at(0).size() == 1)
                    && (height_keys.at(0).count(this->root_split_) > 0));
            int height_idx = height_keys.size() - 2;
            for (unsigned int key_idx = 1; key_idx < height_keys.size(); ++key_idx) {
                std::shared_ptr<HeightSamples> height_sample = this->get_height(height_keys.at(key_idx));
                std::vector<double> node_heights;
                double node_height = std::numeric_limits<double>::quiet_NaN();
                if (height_sample) {
                    node_heights = height_sample->get_heights();
                    SampleSummary<double> node_height_summary(node_heights);
                    node_height = node_height_summary.mean();
                    if (use_median_heights) {
                        node_height = node_height_summary.median();
                    }
                }
                for (auto split : height_keys.at(key_idx)) {
                    std::ostringstream node_label;
                    node_label << "[&height_index=" << height_idx << ",";
                    this->_write_node_annotations<double>(node_heights,
                            node_label,
                            "height",
                            precision);
                    split_parameter_map.clear();
                    for (auto param_key : parameter_keys) {
                        split_parameter_map[param_key] = {};
                    }
                    std::shared_ptr<SplitSamples> split_sample = this->get_split(split);
                    if (split_sample) {
                        split_parameter_map = split_sample->get_parameter_map();
                    }
                    for (auto p_values : split_parameter_map) {
                        if (p_values.first == "height") {
                            continue;
                        }
                        node_label << ","
                        first_pass = false;
                        this->_write_node_annotations<double>(p_values.second,
                                node_label,
                                p_values.first,
                                precision);
                    }
                    std::shared_ptr<Node> new_node = std::make_shared<Node>(
                            node_label.str(), node_height);
                    std::vector<unsigned int> leaf_indices = split.get_leaf_indices();
                    std::shared_ptr<Node> grand_parent_node = leaf_nodes.at(leaf_indices.at(0))->get_parent();
                    for (auto leaf_idx : leaf_indices) {
                        ECOEVOLITY_ASSERT(leaf_nodes.at(leaf_idx)->get_parent() == grand_parent_node);
                        leaf_nodes.at(leaf_idx)->remove_parent();
                        leaf_nodes.at(leaf_idx)->add_parent(new_node);
                        grand_parent_node->add_child(new_node);
                    }
                }
                --height_idx;
            }
            return root_node->to_parentheses(precision, true);
        }
};

} // treesum

#endif
