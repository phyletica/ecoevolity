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
#include "node.hpp"
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

class NumberOfHeightsSamples : public BaseSamples {
    protected:
        unsigned int number_of_heights_;

    public:
        void add_sample(
                const unsigned int number_of_heights,
                unsigned int tree_index,
                unsigned int source_index) { 
            if (this->n_ == 0) {
                this->number_of_heights_ = number_of_heights;
            }
            else {
                ECOEVOLITY_ASSERT(number_of_heights == this->number_of_heights_);
            }
            this->tally_sample_(tree_index, source_index);
        }

        unsigned int get_number_of_heights() const {
            return this->number_of_heights_;
        }
};

class TopologySamples : public BaseSamples {
    protected:
        std::set< std::set<Split> > split_set_;
        std::map< std::set<Split>, std::vector<double> > heights_;
        std::map< Split, std::set<Split> > split_to_node_map_;
        std::map< std::set<Split>, Split > node_to_split_map_;
        std::map< Split, std::set<Split> > split_to_height_split_set_map_;

    public:
        void add_sample(
                const std::map<std::set<Split>, double> & height_map,
                const std::map< Split, std::set<Split> > & node_map,
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
                this->split_to_node_map_ = node_map;
                for (auto split_node : node_map) {
                    this->node_to_split_map_[split_node.second] = split_node.first;
                }
                for (auto splt_set : s_set) {
                    for (auto splt : splt_set) {
                        this->split_to_height_split_set_map_[splt] = splt_set;
                    }
                }
            }
            else {
                ECOEVOLITY_ASSERT(s_set == this->split_set_);
                ECOEVOLITY_ASSERT(node_map == this->split_to_node_map_);
            }
            this->tally_sample_(tree_index, source_index);
        }

        const std::set< std::set<Split> > & get_split_set() const {
            return this->split_set_;
        }

        unsigned int get_number_of_heights() const {
            return this->split_set_.size();
        }
        unsigned int get_number_of_non_leaf_splits() const {
            return this->split_to_node_map_.size();
        }

        const std::set< Split > & get_height_split_set(
                const Split & split) const {
            return this->split_to_height_split_set_map_.at(split);
        }

        const std::map< std::set<Split>, std::vector<double> > & get_height_map() const {
            return this->heights_;
        }

        const std::vector<double> & get_heights(const std::set<Split> & split_set) const {
            return this->heights_.at(split_set);
        }

        const std::map< Split, std::set<Split> > & get_split_to_node_map() const {
            return this->split_to_node_map_;
        }

        const std::map< std::set<Split>, Split > & get_node_to_split_map() const {
            return this->node_to_split_map_;
        }

        bool has_split(const Split & split) const {
            return (this->split_to_node_map_.count(split) > 0);
        }
        bool has_node(const std::set<Split> & split_set) const {
            return (this->node_to_split_map_.count(split_set) > 0);
        }

        const std::set<Split> & get_node_split_set(const Split & split) const {
            ECOEVOLITY_ASSERT(this->has_split(split));
            return this->split_to_node_map_.at(split);
        }
        const Split & get_split(const std::set<Split> & node_split_set) const {
            ECOEVOLITY_ASSERT(this->has_node(node_split_set));
            return this->node_to_split_map_.at(node_split_set);
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
            ECOEVOLITY_ASSERT(set_of_splits.size() > 0);
            this->set_split_set(set_of_splits);
            this->heights_.push_back(height);
            this->tally_sample_(tree_index, source_index);
        }

        const std::set<Split> & get_split_set() const {
            return this->split_set_;
        }

        void set_split_set(const std::set<Split> & set_of_splits) {
            if ((this->n_ > 0) || (this->split_set_.size() > 0)) {
                ECOEVOLITY_ASSERT(set_of_splits == this->split_set_);
            }
            else {
                if (set_of_splits.size() > 1) {
                    // Make sure splits are not nested
                    ECOEVOLITY_ASSERT(Split::can_be_siblings(set_of_splits));
                }
                this->split_set_ = set_of_splits;
            }
        }

        const std::vector<double> & get_heights() const {
            return this->heights_;
        }
};

class NodeHeightSamples : public BaseSamples {
    protected:
        std::set< std::set<Split> > node_set_;
        std::vector<double> heights_;

    public:
        void add_sample(
                const std::set< std::set<Split> > & set_of_nodes,
                double height,
                unsigned int tree_index,
                unsigned int source_index) {
            ECOEVOLITY_ASSERT(set_of_nodes.size() > 0);
            this->set_node_set(set_of_nodes);
            this->heights_.push_back(height);
            this->tally_sample_(tree_index, source_index);
        }

        const std::set< std::set<Split> > & get_node_set() const {
            return this->node_set_;
        }

        void set_node_set(const std::set< std::set<Split> > & set_of_nodes) {
            if ((this->n_ > 0) || (this->node_set_.size() > 0)) {
                ECOEVOLITY_ASSERT(set_of_nodes == this->node_set_);
            }
            else {
                this->node_set_ = set_of_nodes;
            }
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


class NodeSamples : public BaseSamples {
    protected:
        std::set<Split> split_set_;
        std::map<std::string, std::vector<double> > parameters_;

    public:
        void add_sample(
                const std::set<Split> & set_of_splits,
                const std::map<std::string, double> & p,
                unsigned int tree_index,
                unsigned int source_index) {
            ECOEVOLITY_ASSERT(set_of_splits.size() > 0);
            this->set_split_set(set_of_splits);
            for (auto pname_value : p) {
                this->parameters_[pname_value.first].push_back(pname_value.second);
            }
            this->tally_sample_(tree_index, source_index);
        }

        const std::set<Split> & get_split_set() const {
            return this->split_set_;
        }

        Split get_split() const {
            Split ret_split;
            bool split_resized = false;
            for (auto split : this->split_set_) {
                if (! split_resized) {
                    ret_split.resize(split.size());
                    split_resized = true;
                }
                ret_split.add_split(split);
            }
            return ret_split;
        }

        void set_split_set(const std::set<Split> & set_of_splits) {
            if ((this->n_ > 0) || (this->split_set_.size() > 0)) {
                ECOEVOLITY_ASSERT(set_of_splits == this->split_set_);
            }
            else {
                // Make sure splits can be siblings
                ECOEVOLITY_ASSERT(Split::can_be_siblings(set_of_splits));
                this->split_set_ = set_of_splits;
            }
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
        std::vector< std::shared_ptr<NodeHeightSamples> > node_heights_;
        std::vector< std::shared_ptr<SplitSamples> > splits_;
        std::vector< std::shared_ptr<NodeSamples> > nodes_;
        std::vector< std::shared_ptr<SplitSamples> > non_trivial_splits_;
        std::vector< std::shared_ptr<NumberOfHeightsSamples> > num_heights_;
        std::map< std::set< std::set<Split> >, std::shared_ptr<TopologySamples> > topologies_map_;
        std::map< std::set<Split>,             std::shared_ptr<HeightSamples>   > heights_map_;
        std::map< std::set< std::set<Split> >, std::shared_ptr<NodeHeightSamples> > node_heights_map_;
        std::map< Split,                       std::shared_ptr<SplitSamples>    > splits_map_;
        std::map< std::set<Split>,             std::shared_ptr<NodeSamples>    > nodes_map_;
        std::map< unsigned int,                std::shared_ptr<NumberOfHeightsSamples> > num_heights_map_;
        std::vector<double> tree_lengths_;
        std::vector<std::string> source_paths_;
        std::vector<unsigned int> source_sample_sizes_;
        std::vector<unsigned int> source_num_skipped_;
        std::set<Split> trivial_splits_;
        Split root_split_;
        std::vector<Split> leaf_splits_;
        std::vector<std::string> leaf_labels_;
        unsigned int sample_size_ = 0;
        bool target_tree_provided_ = false;
        tree_type target_tree_;
        TopologySamples target_topology_sample_;
        std::set< std::set<Split> > target_topology_;
        std::map<std::set<Split>, double> target_heights_;
        std::map<std::set< std::set<Split> >, double> target_node_heights_;
        std::map<Split, std::map<std::string, double> > target_split_parameters_;
        std::map<std::set<Split>, std::map<std::string, double> > target_node_parameters_;
        std::vector<double> target_euclidean_distances_;
        std::set<std::string> constrained_node_parameters_;

        void reverse_sort_samples_by_freq_() {
            std::sort(this->topologies_.begin(), this->topologies_.end(),
                    BaseSamples::reverse_sort_by_n);
            std::sort(this->heights_.begin(), this->heights_.end(),
                    BaseSamples::reverse_sort_by_n);
            std::sort(this->node_heights_.begin(), this->node_heights_.end(),
                    BaseSamples::reverse_sort_by_n);
            std::sort(this->splits_.begin(), this->splits_.end(),
                    BaseSamples::reverse_sort_by_n);
            std::sort(this->non_trivial_splits_.begin(), this->non_trivial_splits_.end(),
                    BaseSamples::reverse_sort_by_n);
            std::sort(this->nodes_.begin(), this->nodes_.end(),
                    BaseSamples::reverse_sort_by_n);
            std::sort(this->num_heights_.begin(), this->num_heights_.end(),
                    BaseSamples::reverse_sort_by_n);
        }

        void update_constrained_node_parameters_() {
            this->constrained_node_parameters_.clear();
            const std::map< std::string, std::vector<double> > & root_parameters =
                this->get_split(this->root_split_)->get_parameter_map();
            const std::map< std::string, std::vector<double> > & leaf_parameters =
                this->get_split(this->leaf_splits_.at(0))->get_parameter_map();
            for (auto pname_v : root_parameters) {
                ECOEVOLITY_ASSERT(root_parameters.at(pname_v.first).size() ==
                    leaf_parameters.at(pname_v.first).size());
                // If the root and first leaf share the same sampled values, we
                // are going to assume this parameter is constrained across the
                // tree
                if (root_parameters.at(pname_v.first) == leaf_parameters.at(pname_v.first)) {
                    this->constrained_node_parameters_.insert(pname_v.first);
                }
            }
        }

        void update_splits_and_labels_(const tree_type & tree) {
            unsigned int nleaves = tree.get_leaf_node_count();
            ECOEVOLITY_ASSERT(nleaves > 0);
            tree.get_leaf_labels(this->leaf_labels_);
            ECOEVOLITY_ASSERT(this->leaf_labels_.size() == nleaves);
            // Leaf labels are sorted for every tree parsed by BaseTree.
            // Sorting here to make sure order is consistent with sampled trees
            // (parsed by BaseTree)
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

        void _update_target_tree() {
            this->target_topology_.clear();
            this->target_heights_.clear();
            this->target_node_heights_.clear();
            this->target_split_parameters_.clear();
            this->target_node_parameters_.clear();
            std::map< Split, std::set<Split> > node_map;
            this->target_tree_.store_splits_heights_parameters(
                    this->target_topology_,
                    this->target_heights_,
                    this->target_node_heights_,
                    node_map,
                    this->target_split_parameters_,
                    this->target_node_parameters_,
                    false);
            this->target_tree_provided_ = true;
            this->check_target_tree_();
            this->target_topology_sample_ = TopologySamples();
            this->target_topology_sample_.add_sample(
                    this->target_heights_,
                    node_map,
                    0, 0);
        }

        template <typename T>
        void _write_summary_of_values(
                const std::vector<T> & values,
                std::ostream & out,
                const std::string & parameter_name = "",
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            out.precision(precision);
            std::string p_name = parameter_name;
            if (parameter_name != "") {
                p_name += "_";
            }
            if (values.size() < 1) {
                double nan = std::numeric_limits<double>::quiet_NaN();
                out << margin << p_name << "n: " << 0 << "\n"
                    << margin << p_name << "ess: " << 0 << "\n"
                    << margin << p_name << "mean: " << nan << "\n"
                    << margin << p_name << "median: " << nan << "\n"
                    << margin << p_name << "std_dev: " << nan << "\n"
                    << margin << p_name << "range: ["
                                        << nan << ", "
                                        << nan << "]\n"
                    << margin << p_name << "eti_95: ["
                                        << nan << ", "
                                        << nan << "]\n"
                    << margin << p_name << "hpdi_95: ["
                                        << nan << ", "
                                        << nan << "]"
                    << std::endl;
                return;
            }
            SampleSummary<T> summary(values);
            double ess = effective_sample_size<T>(values, true);
            out << margin << p_name << "n: " << summary.sample_size() << "\n"
                << margin << p_name << "ess: " << ess << "\n"
                << margin << p_name << "mean: " << summary.mean() << "\n"
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
        void _write_summary_of_values(
                const std::vector< std::vector<T> > & values,
                std::ostream & out,
                const std::string & parameter_name = "",
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            out.precision(precision);
            std::string p_name = parameter_name;
            if (parameter_name != "") {
                p_name += "_";
            }
            unsigned int nvals = 0;
            for (auto vec : values) {
                nvals += vec.size();
            }
            std::vector<T> s;
            s.reserve(nvals);
            for (auto vec : values) {
                for (auto val : vec) {
                    s.push_back(val);
                }
            }
            this->_write_summary_of_values(s, out, parameter_name, margin, precision);
            double psrf = potential_scale_reduction_factor<T>(values);
            out << margin << p_name << "psrf: " << psrf << std::endl;
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
                out << p_name << "n=" << 0 << ","
                    << p_name << "ess=" << 0 << ","
                    << p_name << "mean=" << nan << ","
                    << p_name << "median=" << nan << ","
                    << p_name << "std_dev=" << nan << ","
                    << p_name << "range={"
                              << nan << ","
                              << nan << "},"
                    << p_name << "eti_95={"
                              << nan << ","
                              << nan << "},"
                    << p_name << "hpdi_95={"
                              << nan << ","
                              << nan << "}";
                return;
            }
            SampleSummary<T> summary(values);
            double ess = effective_sample_size<T>(values, true);
            out << p_name << "n=" << summary.sample_size() << ","
                << p_name << "ess=" << ess << ","
                << p_name << "mean=" << summary.mean() << ","
                << p_name << "median=" << summary.median() << ","
                << p_name << "std_dev=" << summary.std_dev() << ","
                << p_name << "range={"
                          << summary.min() << ","
                          << summary.max() << "},"
                << p_name << "eti_95={"
                          << summary.qi_95().first << ","
                          << summary.qi_95().second << "},"
                << p_name << "hpdi_95={"
                          << summary.hpdi_95().first << ","
                          << summary.hpdi_95().second << "}";
        }

        template <typename T>
        std::vector< std::vector<T> > _get_values_by_source(
                const std::vector<T> & values
                ) const {
            ECOEVOLITY_ASSERT(values.size() == this->get_sample_size());
            std::vector< std::vector<T> > source_values(
                    this->source_sample_sizes_.size());
            for (unsigned int i = 0; i < this->source_sample_sizes_.size(); ++i) {
                source_values.at(i).reserve(this->source_sample_sizes_.at(i));
            }

            const std::vector<unsigned int> & source_indices = this->get_source_indices();
            ECOEVOLITY_ASSERT(source_indices.size() == this->get_sample_size());
            for (unsigned int i = 0; i < source_indices.size(); ++i) {
                source_values.at(source_indices.at(i)).push_back(values.at(i));
            }
            return source_values;
        }

        void _add_tree(
                const tree_type & tree,
                const unsigned int tree_index,
                const unsigned int source_index) {
            if ((this->sample_size_ < 1) && (this->leaf_labels_.size() < 1)) {
                this->update_splits_and_labels_(tree);
            }
            this->check_leaf_labels_(tree);

            std::set< std::set<Split> > split_set;
            std::map<std::set<Split>, double> heights;
            std::map<std::set< std::set<Split> >, double> node_heights;
            std::map<Split, std::map<std::string, double> > split_parameters;
            std::map<std::set<Split>, std::map<std::string, double> > node_parameters;
            std::map<Split, std::set< Split > > node_map;
            tree.store_splits_heights_parameters(
                    split_set,
                    heights,
                    node_heights,
                    node_map,
                    split_parameters,
                    node_parameters,
                    false);
            unsigned int nheights = heights.size();
            if (this->num_heights_map_.count(nheights) > 0) {
                this->num_heights_map_[nheights]->add_sample(
                        nheights,
                        tree_index,
                        source_index);
            }
            else {
                std::shared_ptr<NumberOfHeightsSamples> nhs = std::make_shared<NumberOfHeightsSamples>();
                nhs->add_sample(nheights, tree_index, source_index);
                this->num_heights_.push_back(nhs);
                this->num_heights_map_[nheights] = nhs;
            }
            if (this->topologies_map_.count(split_set) > 0) {
                this->topologies_map_[split_set]->add_sample(
                        heights,
                        node_map,
                        tree_index,
                        source_index);
            }
            else {
                std::shared_ptr<TopologySamples> ts = std::make_shared<TopologySamples>();
                ts->add_sample(heights, node_map, tree_index, source_index);
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
            for (auto node_height : node_heights) {
                if (this->node_heights_map_.count(node_height.first) > 0) {
                    this->node_heights_map_[node_height.first]->add_sample(
                            node_height.first,
                            node_height.second,
                            tree_index,
                            source_index);
                }
                else {
                    std::shared_ptr<NodeHeightSamples> nhs = std::make_shared<NodeHeightSamples>();
                    nhs->add_sample(
                            node_height.first,
                            node_height.second,
                            tree_index,
                            source_index);
                    this->node_heights_.push_back(nhs);
                    this->node_heights_map_[node_height.first] = nhs;
                }
            }
            for (auto split_pmap : split_parameters) {
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
            for (auto split_set_pmap : node_parameters) {
                if (this->nodes_map_.count(split_set_pmap.first) > 0) {
                    this->nodes_map_[split_set_pmap.first]->add_sample(
                            split_set_pmap.first,
                            split_set_pmap.second,
                            tree_index,
                            source_index);
                }
                else {
                    std::shared_ptr<NodeSamples> ns = std::make_shared<NodeSamples>();
                    ns->add_sample(
                            split_set_pmap.first,
                            split_set_pmap.second,
                            tree_index,
                            source_index);
                    this->nodes_.push_back(ns);
                    this->nodes_map_[split_set_pmap.first] = ns;
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

        void _annotate_nodes(
                const TopologySamples & topo_sample,
                std::shared_ptr<NodeType> root_node,
                const unsigned int precision = 12) const {
            // Let's get sorted node heights so that we can use them for
            // looking up height indices
            std::vector< std::shared_ptr<PositiveRealParameter> > node_heights;
            std::vector< std::shared_ptr<NodeType> > internal_nodes = root_node->get_internal_nodes();
            for (unsigned int i = 0; i < internal_nodes.size(); ++i) {
                bool exists = false;
                // Check if we already have a pointer to the same place in memory
                for (unsigned int j = 0; j < node_heights.size(); ++j) {
                    if (internal_nodes.at(i)->get_height_parameter() == node_heights.at(j)) {
                        exists = true;
                    }
                }
                if (! exists) {
                    node_heights.push_back(internal_nodes.at(i)->get_height_parameter());
                }
            }
            std::sort(node_heights.begin(), node_heights.end(), PositiveRealParameter::sort_by_value);
            std::map< std::shared_ptr<PositiveRealParameter>, unsigned int > heights_to_index_map;
            for (unsigned int i = 0; i < node_heights.size(); ++i) {
                heights_to_index_map[node_heights.at(i)] = i;
            }

            // Collect values of parameters that are constrained across all
            // branches of the tree; we can do this just using the root split
            std::map<std::string, std::vector<double> > constrained_params;
            for (auto param : this->constrained_node_parameters_) {
                constrained_params[param] = this->get_split(this->root_split_)->get_values(param);
            }
            // Get parameter map using root split; we need the keys of
            // this map below for each internal node
            std::map<std::string, std::vector<double> > split_parameter_map;
            split_parameter_map = this->get_split(this->root_split_)->get_parameter_map();

            // Do a reverse level-order traversal of tree in order to annotate
            // the internal node labels.
            root_node->resize_splits(this->get_number_of_leaves());
            std::vector< std::shared_ptr<NodeType> > level_ordered_nodes;
            root_node->level_order(level_ordered_nodes);

            for (auto node = level_ordered_nodes.rbegin();
                    node != level_ordered_nodes.rend();
                    ++node) {
                if (! (*node)->is_leaf()) {
                    Split split = (*node)->split_;
                    std::set<Split> height_split_set = topo_sample.get_height_split_set(split);
                    unsigned int height_idx = heights_to_index_map.at((*node)->get_height_parameter());

                    // Get splits that descend from this node and use them to
                    // get summaries of node samples
                    std::set<Split> node_split_set;
                    for (auto child : (*node)->get_all_children()) {
                        node_split_set.insert(child->split_);
                    }

                    std::shared_ptr<HeightSamples> height_sample = this->get_height(height_split_set);
                    std::vector<double> index_heights;
                    if (height_sample) {
                        index_heights = height_sample->get_heights();
                        SampleSummary<double> index_height_summary(index_heights);
                    }

                    std::ostringstream node_label;
                    node_label.precision(precision);
                    node_label << "[&height_index=" << height_idx
                        << ",index_freq=" << this->get_height_frequency(height_split_set)
                        << ",split_freq=" << this->get_split_frequency(split)
                        << ",";
                    this->_write_node_annotations<double>(index_heights,
                            node_label,
                            "index_height",
                            precision);
                    for (auto key_params : split_parameter_map) {
                        key_params.second.clear();
                    }
                    std::shared_ptr<SplitSamples> split_sample = this->get_split(split);
                    if (split_sample) {
                        split_parameter_map = split_sample->get_parameter_map();
                    }
                    for (auto p_values : split_parameter_map) {
                        node_label << ",";
                        if (p_values.first == "height") {
                            this->_write_node_annotations<double>(p_values.second,
                                    node_label,
                                    "split_height",
                                    precision);
                        }
                        else if (constrained_params.count(p_values.first) > 0) {
                            this->_write_node_annotations<double>(
                                    constrained_params.at(p_values.first),
                                    node_label,
                                    p_values.first,
                                    precision);
                        }
                        else {
                            this->_write_node_annotations<double>(p_values.second,
                                    node_label,
                                    p_values.first,
                                    precision);
                        }
                    }
                    node_label << ",node_freq="
                               << this->get_node_frequency(node_split_set)
                               << "]";
                    (*node)->set_label(node_label.str());
                }
                else {
                    // Set bit for this leaf node's index
                    (*node)->split_.set_leaf_bit((*node)->get_index());

                    Split split = (*node)->split_;
                    std::ostringstream label;
                    label.precision(precision);
                    unsigned int leaf_index = split.get_leaf_indices().at(0);
                    ECOEVOLITY_ASSERT(split == this->leaf_splits_.at(leaf_index));
                    label << this->leaf_labels_.at(leaf_index) << "[&";
                    split_parameter_map = this->get_split(split)->get_parameter_map();
                    SampleSummary<double> leaf_height_summary(split_parameter_map.at("height"));
                    bool first_pass = true;
                    for (auto p_values : split_parameter_map) {
                        if (! first_pass) {
                            label << ",";
                        }
                        first_pass = false;
                        if (p_values.first == "height") {
                            this->_write_node_annotations<double>(p_values.second,
                                    label,
                                    "index_height",
                                    precision);
                            label << ",";
                            this->_write_node_annotations<double>(p_values.second,
                                    label,
                                    "split_height",
                                    precision);
                        }
                        else {
                            this->_write_node_annotations<double>(p_values.second,
                                    label,
                                    p_values.first,
                                    precision);
                        }
                    }
                    label << "]";
                    (*node)->set_label(label.str());
                }
                if ((*node)->has_parent()) {
                    (*node)->get_parent()->split_.add_split((*node)->split_);
                }
            }
        }

        void _process_node_splits(
                const Split & split,
                const TopologySamples & topo_sample,
                std::shared_ptr<NodeType> root_node,
                const std::vector< std::shared_ptr<NodeType> > & leaf_nodes,
                std::map< std::set<Split>, std::shared_ptr<PositiveRealParameter> > & split_set_to_height_parameter_map,
                const bool use_median_heights = false) const {
            for (auto s : topo_sample.get_node_split_set(split)) {
                if (s.get_leaf_node_count() < 2) {
                    // We have a leaf; do not process
                    continue;
                }
                std::shared_ptr<HeightSamples> height_sample = this->get_height(
                        topo_sample.get_height_split_set(s));
                std::vector<double> index_heights;
                double index_height = std::numeric_limits<double>::quiet_NaN();
                if (height_sample) {
                    index_heights = height_sample->get_heights();
                    SampleSummary<double> index_height_summary(index_heights);
                    index_height = index_height_summary.mean();
                    if (use_median_heights) {
                        index_height = index_height_summary.median();
                    }
                }
                std::shared_ptr<NodeType> new_node = std::make_shared<NodeType>(index_height);
                if (split_set_to_height_parameter_map.count(topo_sample.get_height_split_set(s)) < 1) {
                    split_set_to_height_parameter_map[topo_sample.get_height_split_set(s)] = new_node->get_height_parameter();
                }
                else {
                    new_node->set_height_parameter(split_set_to_height_parameter_map.at(topo_sample.get_height_split_set(s)));
                }
                std::vector<unsigned int> leaf_indices = s.get_leaf_indices();
                std::shared_ptr<NodeType> grand_parent_node = leaf_nodes.at(leaf_indices.at(0))->get_parent();
                for (auto leaf_idx : leaf_indices) {
                    ECOEVOLITY_ASSERT(leaf_nodes.at(leaf_idx)->get_parent() == grand_parent_node);
                    leaf_nodes.at(leaf_idx)->remove_parent();
                    leaf_nodes.at(leaf_idx)->add_parent(new_node);
                }
                grand_parent_node->add_child(new_node);

                this->_process_node_splits(s,
                        topo_sample,
                        root_node,
                        leaf_nodes,
                        split_set_to_height_parameter_map,
                        use_median_heights);
            }
        }

        void _write_to_nexus(
                std::vector< std::shared_ptr<NodeType> > trees,
                std::ostream & out,
                const unsigned int precision = 12) const {
            out.precision(precision);
            out<< "#NEXUS\n\n"
                << "BEGIN TAXA;\n"
                << "    DIMENSIONS NTAX=" << this->leaf_labels_.size() << ";\n"
                << "    TAXLABELS\n";
            for (auto label : this->leaf_labels_) {
                out << "        " << label << "\n";
            }
            out << "    ;\n"
                << "END;\n\n"
                << "BEGIN TREES;\n";
            unsigned int tree_idx = 0;
            for (auto t : trees) {
                out << "    TREE tree" << ++tree_idx << " = [&R] "
                    << t->to_parentheses(precision, true)
                    << ";\n";
            }
            out << "END;" << std::endl;
        }


    public:

        TreeSample() { }
        TreeSample(
                const std::vector<std::string> & paths,
                const std::string & ncl_file_format,
                const unsigned int skip = 0,
                const double ultrametricity_tolerance = 1e-6,
                const double multiplier = -1.0) {
            for (auto path : paths) {
                this->add_trees(path, ncl_file_format, skip,
                        ultrametricity_tolerance,
                        multiplier);
            }
        }
        TreeSample(
                const std::string & target_tree_path,
                const std::vector<std::string> & paths,
                const std::string & target_ncl_file_format,
                const std::string & ncl_file_format,
                const unsigned int skip = 0,
                const double ultrametricity_tolerance = 1e-6,
                const double multiplier = -1.0) {
            this->set_target_tree(target_tree_path, target_ncl_file_format);
            for (auto path : paths) {
                this->add_trees(path, ncl_file_format, skip,
                        ultrametricity_tolerance,
                        multiplier);
            }
        }

        void add_trees(
                const std::string & path,
                const std::string & ncl_file_format,
                const unsigned int skip = 0,
                const double ultrametricity_tolerance = 1e-6,
                const double multiplier = -1.0) {
            this->source_paths_.push_back(path);
            std::ifstream in_stream;
            in_stream.open(path);
            if (! in_stream.is_open()) {
                throw EcoevolityParsingError(
                        "Could not open tree file",
                        path);
            }
            try {
                this->add_trees(in_stream,
                        ncl_file_format,
                        skip,
                        ultrametricity_tolerance,
                        multiplier);
            }
            catch(...) {
                std::cerr << "ERROR: Problem parsing tree file path: "
                        << path << "\n";
                throw;
            }
        }

        void add_trees(
                std::istream & tree_stream,
                const std::string & ncl_file_format,
                const unsigned int skip = 0,
                const double ultrametricity_tolerance = 1e-6,
                const double multiplier = -1.0) {
            this->source_num_skipped_.push_back(skip);
            unsigned int source_index = this->source_sample_sizes_.size();
            this->source_sample_sizes_.push_back(0);

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
                tree_type t(tree_description,
                        taxa_block,
                        ultrametricity_tolerance,
                        multiplier);
                this->_add_tree(t, i, source_index);
            }
            nexus_reader.DeleteBlocksFromFactories();

            unsigned int source_total = 0;
            for (unsigned int n : this->source_sample_sizes_) {
                source_total += n;
            }
            ECOEVOLITY_ASSERT(source_total == this->sample_size_);
            this->reverse_sort_samples_by_freq_();
            this->update_constrained_node_parameters_();
        }

        void set_target_tree(
                std::istream & tree_stream,
                const std::string & ncl_file_format) {
            this->target_tree_ = tree_type(tree_stream, ncl_file_format);
            this->_update_target_tree();
        }

        void set_target_tree(
                const std::string & tree_path,
                const std::string & ncl_file_format) {
            this->target_tree_ = tree_type(tree_path, ncl_file_format);
            this->_update_target_tree();
        }

        void set_target_tree(
                const std::string & newick_tree_string) {
            this->target_tree_ = tree_type(newick_tree_string);
            this->_update_target_tree();
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

        unsigned int get_source_burnin(unsigned int source_index) const {
            return this->source_num_skipped_.at(source_index);
        }

        bool equal_source_sample_sizes() const {
            unsigned int n = this->source_sample_sizes_.at(0);
            for (unsigned int i = 1; i < this->source_sample_sizes_.size(); ++i) {
                if (this->source_sample_sizes_.at(i) != n) {
                    return false;
                }
            }
            return true;
        }

        bool has_multiple_sources_with_equal_n() const {
            return ((this->get_number_of_sources() > 1) && this->equal_source_sample_sizes());
        }

        const std::vector<unsigned int> & get_source_indices() const {
            std::shared_ptr<SplitSamples> root_sample = this->get_split(this->root_split_);
            ECOEVOLITY_ASSERT(root_sample->get_sample_size() == this->get_sample_size());
            return root_sample->get_source_indices();
        }

        const std::vector<double> & get_tree_lengths() const {
            return this->tree_lengths_;
        }
        std::vector< std::vector<double> > get_tree_lengths_by_source() const {
            return this->_get_values_by_source(this->tree_lengths_);
        }

        const std::vector<double> & get_root_parameter_values(
                const std::string & parameter_name) const {
            return this->get_split(this->root_split_)->get_values(parameter_name);
        }
        std::vector< std::vector<double> > get_root_parameter_values_by_source(
                const std::string & parameter_name) const {
            return this->_get_values_by_source(this->get_root_parameter_values(parameter_name));
        }

        const std::vector< std::shared_ptr<TopologySamples> > & get_topologies() const {
            return this->topologies_;
        }
        const std::vector< std::shared_ptr<HeightSamples> > & get_heights() const {
            return this->heights_;
        }
        const std::vector< std::shared_ptr<NodeHeightSamples> > & get_node_heights() const {
            return this->node_heights_;
        }
        const std::vector< std::shared_ptr<SplitSamples> > & get_splits() const {
            return this->splits_;
        }
        const std::vector< std::shared_ptr<SplitSamples> > & get_non_trivial_splits() const {
            return this->non_trivial_splits_;
        }
        const std::vector< std::shared_ptr<NodeSamples> > & get_nodes() const {
            return this->nodes_;
        }
        const std::vector< std::shared_ptr<NumberOfHeightsSamples> > & get_all_numbers_of_heights() const {
            return this->num_heights_;
        }
        std::vector< std::vector<unsigned int> > get_all_numbers_of_heights_by_source() const {
            std::vector< std::vector<unsigned int> > source_values(
                    this->source_sample_sizes_.size());
            for (unsigned int i = 0; i < this->source_sample_sizes_.size(); ++i) {
                source_values.at(i).reserve(this->source_sample_sizes_.at(i));
            }

            for (auto nhs : this->get_all_numbers_of_heights()) {
                for (auto src_idx : nhs->get_source_indices()) {
                    source_values.at(src_idx).push_back(nhs->get_number_of_heights());
                }
            }
            for (unsigned int i = 0; i < this->source_sample_sizes_.size(); ++i) {
                ECOEVOLITY_ASSERT(source_values.at(i).size() == this->source_sample_sizes_.at(i));
            }
            return source_values;
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
        std::shared_ptr<NodeSamples> get_node(
                const std::set<Split> & split_set
                ) const {
            if (this->nodes_map_.count(split_set) < 1) {
                return nullptr;
            }
            return this->nodes_map_.at(split_set);
        }
        std::vector< std::shared_ptr<NodeSamples> > get_nodes_of_split(
                const Split & split
                ) const {
            std::vector< std::shared_ptr<NodeSamples> > node_samples;
            for (auto ns : this->nodes_) {
                if (ns->get_split() == split) {
                    node_samples.push_back(ns);
                }
            }
            return node_samples;
        }
        std::set<Split> get_target_node_split_set_of_split(
                const Split & split
                ) const {
            for (auto splitset_params : this->target_node_parameters_) {
                if (split.is_parent_of(splitset_params.first)) {
                    return splitset_params.first;
                }
            }
            throw EcoevolityError("Split not in target tree");
        }
        std::shared_ptr<NumberOfHeightsSamples> get_number_of_heights(
                const unsigned int number_of_heights
                ) const {
            if (this->num_heights_map_.count(number_of_heights) < 1) {
                return nullptr;
            }
            return this->num_heights_map_.at(number_of_heights);
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

        unsigned int get_node_height_count(
                const std::set< std::set<Split> > & node_set
                ) const {
            if (this->node_heights_map_.count(node_set) < 1) {
                return 0;
            }
            return this->node_heights_map_.at(node_set)->get_sample_size();
        }
        double get_node_height_frequency(
                const std::set< std::set<Split> > & node_set
                ) const {
            return this->get_node_height_count(node_set) / (double)this->sample_size_;
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

        unsigned int get_node_count(
                const std::set<Split> & split_set
                ) const {
            if (this->nodes_map_.count(split_set) < 1) {
                return 0;
            }
            return this->nodes_map_.at(split_set)->get_sample_size();
        }
        double get_node_frequency(
                const std::set<Split> & split_set
                ) const {
            return this->get_node_count(split_set) / (double)this->sample_size_;
        }

        unsigned int get_number_of_heights_count(
                const unsigned int number_of_heights
                ) const {
            if (this->num_heights_map_.count(number_of_heights) < 1) {
                return 0;
            }
            return this->num_heights_map_.at(number_of_heights)->get_sample_size();
        }
        double get_number_of_heights_frequency(
                const unsigned int number_of_heights
                ) const {
            return this->get_number_of_heights_count(number_of_heights) / (double)this->sample_size_;
        }

        double get_topology_credibility_level(
                const std::set< std::set<Split> > & topology
                ) const {
            unsigned int topo_count = this->get_topology_count(topology);
            if (topo_count < 1) {
                return 0.0;
            }
            double total_prob = 0.0;
            for (auto topo_sample : this->get_topologies()) {
                if (topo_count == topo_sample->get_sample_size()) {
                    return 1.0 - total_prob;
                }
                total_prob += this->get_topology_frequency(
                        topo_sample->get_split_set());
            }
            ECOEVOLITY_ASSERT(fabs(total_prob - 1.0) < 1e-6);
            return 0.0;
        }

        double get_number_of_heights_credibility_level(
                const unsigned int number_of_heights
                ) const {
            unsigned int count = this->get_number_of_heights_count(number_of_heights);
            if (count < 1) {
                return 0.0;
            }
            double total_prob = 0.0;
            for (auto nh_sample : this->get_all_numbers_of_heights()) {
                if (count == nh_sample->get_sample_size()) {
                    return 1.0 - total_prob;
                }
                total_prob += this->get_number_of_heights_frequency(
                        nh_sample->get_number_of_heights());
            }
            ECOEVOLITY_ASSERT(fabs(total_prob - 1.0) < 1e-6);
            return 0.0;
        }

        SampleSummarizer<double> get_summary_of_split_freq_std_devs(
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
            return std_devs_of_split_freqs;
        }

        double get_average_std_dev_of_split_freqs(
                double min_frequency = 0.1) const {
            if (this->get_number_of_sources() < 2) {
                throw EcoevolityError("Calculating the ASDSF requires multiple chains");
            }
            SampleSummarizer<double> std_devs_of_split_freqs = this->get_summary_of_split_freq_std_devs(
                    min_frequency);
            return std_devs_of_split_freqs.mean();
        }

        void write_summary_of_splits(std::ostream & out,
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            std::string indent = string_util::get_indent(1);
            out.precision(precision);

            std::string item_margin = margin + indent + indent;

            out << margin << "splits:\n";
            out << margin << indent << "root:\n";
            this->write_summary_of_trivial_split(this->root_split_,
                    out,
                    false,
                    true,
                    item_margin,
                    precision);
            item_margin += "  ";
            out << margin << indent << "leaves:\n";
            for (auto s : this->leaf_splits_) {
                out << margin << indent << indent << "-\n";
                this->write_summary_of_trivial_split(s,
                        out,
                        true,
                        false,
                        item_margin,
                        precision);
            }
            out << margin << indent << "nontrivial_splits:\n";
            for (auto split_samples : this->get_non_trivial_splits()) {
                out << margin << indent << indent << "-\n";
                this->write_summary_of_nontrivial_split(split_samples->get_split(),
                        out,
                        true,
                        true,
                        item_margin,
                        precision);
            }
        }

        void write_summary_of_trivial_split(const Split & split,
                std::ostream & out,
                const bool include_leaf_indices = true,
                const bool include_nodes = false,
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            if (! this->has_multiple_sources_with_equal_n()) {
                this->write_summary_of_nontrivial_split(split,
                        out,
                        include_leaf_indices,
                        include_nodes,
                        margin,
                        precision);
                return;
            }
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
                std::vector< std::vector<double> > vals_by_source = this->_get_values_by_source(
                        param_vals.second);
                this->_write_summary_of_values<double>(
                        vals_by_source,
                        out,
                        param_vals.first,
                        margin,
                        precision);
            }
            if ((split.get_leaf_node_count() > 2) && include_nodes) {
                out << margin << "nodes:\n";
                std::string indent = string_util::get_indent(1);
                std::string item_margin = margin + indent + "  ";
                unsigned int cumulative_count = 0;
                for (auto node_samples : this->get_nodes_of_split(split)) {
                    out << margin << indent << "-\n";
                    cumulative_count += this->write_summary_of_node(
                            node_samples->get_split_set(),
                            out,
                            item_margin,
                            precision);
                }
                ECOEVOLITY_ASSERT(cumulative_count == this->get_split_count(split));
            }
        }

        void write_summary_of_nontrivial_split(const Split & split,
                std::ostream & out,
                const bool include_leaf_indices = true,
                const bool include_nodes = false,
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
            if (this->splits_map_.count(split) < 1) {
                // This split was not sampled (this can happen for target
                // tree), so we can't call get_split on this split, and we need
                // to feed _write_summary_of_values empty vectors of values
                std::vector<double> empty_values;
                for (auto param_vals : this->get_split(this->leaf_splits_.at(0))->get_parameter_map()) {
                    this->_write_summary_of_values<double>(
                            empty_values,
                            out,
                            param_vals.first,
                            margin,
                            precision);
                }
                return;
            }
            for (auto param_vals : this->get_split(split)->get_parameter_map()) {
                this->_write_summary_of_values<double>(
                        param_vals.second,
                        out,
                        param_vals.first,
                        margin,
                        precision);
            }
            if ((split.get_leaf_node_count() > 2) && include_nodes) {
                out << margin << "nodes:\n";
                std::string indent = string_util::get_indent(1);
                std::string item_margin = margin + indent + "  ";
                unsigned int cumulative_count = 0;
                for (auto node_samples : this->get_nodes_of_split(split)) {
                    out << margin << indent << "-\n";
                    cumulative_count += this->write_summary_of_node(
                            node_samples->get_split_set(),
                            out,
                            item_margin,
                            precision);
                }
                ECOEVOLITY_ASSERT(cumulative_count == this->get_split_count(split));
            }
        }

        unsigned int write_summary_of_node(const std::set<Split> & split_set,
                std::ostream & out,
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            std::string indent = string_util::get_indent(1);
            out.precision(precision);
            out << margin << "number_of_descendants: " << split_set.size() << "\n"
                << margin << "descendant_splits:\n";
            for (auto split : split_set) {
                std::vector<unsigned int> leaf_indices = split.get_leaf_indices();
                out << margin << indent << "- [" << leaf_indices.at(0);
                for (unsigned int i = 1; i < leaf_indices.size(); ++i) {
                    out << ", " << leaf_indices.at(i);
                }
                out << "]\n";
            }
            double count = this->get_node_count(split_set);
            out << margin << "count: " << count << "\n"
                << margin << "frequency: " << this->get_node_frequency(split_set) << "\n";
            if (this->nodes_map_.count(split_set) < 1) {
                // This node was not sampled (this can happen for target tree),
                // so we can't call get_split_set on this node, and we need to
                // feed _write_summary_of_values empty vectors of values
                std::vector<double> empty_values;
                for (auto param_vals : this->nodes_.at(0)->get_parameter_map()) {
                    this->_write_summary_of_values<double>(
                            empty_values,
                            out,
                            param_vals.first,
                            margin,
                            precision);
                }
                return count;
            }
            for (auto param_vals : this->get_node(split_set)->get_parameter_map()) {
                this->_write_summary_of_values<double>(
                        param_vals.second,
                        out,
                        param_vals.first,
                        margin,
                        precision);
            }
            return count;
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
            if (this->heights_map_.count(split_set) < 1) {
                // This height was not sampled (this can happen for target
                // tree), so we can't call get_height, and we need to feed
                // _write_summary_of_values an empty vector of values
                std::vector<double> empty_heights;
                this->_write_summary_of_values<double>(
                        empty_heights,
                        out,
                        "",
                        margin,
                        precision);
                return;
            }
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
            std::string topo_margin = margin + indent + "  ";
            out << margin << "topologies:\n";
            for (auto topo_samples : this->get_topologies()) {
                out << margin << indent << "-\n";
                cumulative_freq += this->write_summary_of_topology(
                        topo_samples->get_split_set(),
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
            double freq = this->get_topology_frequency(split_set);
            out << margin << "count: " << this->get_topology_count(split_set) << "\n"
                << margin << "frequency: " << freq << "\n"
                << margin << "cumulative_frequency: " << freq + cumulative_freq << "\n"
                << margin << "number_of_heights: " << split_set.size() << "\n"
                << margin << "heights:\n";
            std::string h_margin = margin + indent + "  ";
            std::string n_margin = h_margin + indent + "  ";
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
                    if (leaf_indices.size() > 2) {
                        out << n_margin << "node:\n";
                        out << n_margin << indent << "descendant_splits:\n";
                        for (auto node_split : this->get_topology(split_set)->get_node_split_set(split)) {
                            leaf_indices = node_split.get_leaf_indices();
                            out << n_margin << indent << indent << "- [" << leaf_indices.at(0);
                            for (unsigned int i = 1; i < leaf_indices.size(); ++i) {
                                out << ", " << leaf_indices.at(i);
                            }
                            out << "]\n";
                        }
                    }
                }
                this->_write_summary_of_values<double>(
                        this->get_topology(split_set)->get_heights(s_set),
                        out,
                        "",
                        h_margin,
                        precision);
            }
            return freq;
        }

        void write_summary_of_all_numbers_of_heights(std::ostream & out,
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            std::string indent = string_util::get_indent(1);
            out.precision(precision);

            out << margin << "number_of_heights_summary:\n";
            std::vector< std::vector<unsigned int> > n_heights = this->get_all_numbers_of_heights_by_source();
            std::vector<unsigned int> n_hts;
            n_hts.reserve(this->get_sample_size());
            for (auto n_vec : n_heights) {
                for (auto num_ht : n_vec) {
                    n_hts.push_back(num_ht);
                }
            }
            ECOEVOLITY_ASSERT(n_hts.size() == this->get_sample_size());
            this->_write_summary_of_values<unsigned int>(
                    n_hts,
                    out,
                    "",
                    margin + indent,
                    precision);

            out << margin << "numbers_of_heights:\n";
            double cumulative_freq = 0.0;
            std::string nh_margin = margin + indent + "  ";
            for (auto nh_samples : this->get_all_numbers_of_heights()) {
                out << margin << indent << "-\n";
                cumulative_freq += this->write_summary_of_number_of_heights(
                        nh_samples->get_number_of_heights(),
                        out,
                        cumulative_freq,
                        nh_margin,
                        precision);
            }
        }

        double write_summary_of_number_of_heights(
                const unsigned int number_of_heights,
                std::ostream & out,
                const double cumulative_freq = 0.0,
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            std::string indent = string_util::get_indent(1);
            out.precision(precision);
            double freq = this->get_number_of_heights_frequency(number_of_heights);
            out << margin << "number_of_heights: " << number_of_heights << "\n"
                << margin << "count: " << this->get_number_of_heights_count(number_of_heights) << "\n"
                << margin << "frequency: " << freq << "\n"
                << margin << "cumulative_frequency: " << freq + cumulative_freq << std::endl;
            return freq;
        }

        void write_summary_of_tree_lengths(std::ostream & out,
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            std::string indent = string_util::get_indent(1);
            out.precision(precision);

            double cumulative_freq = 0.0;
            out << margin << "tree_length:\n";
            if (this->has_multiple_sources_with_equal_n()) {
                this->_write_summary_of_values<double>(
                        this->get_tree_lengths_by_source(),
                        out,
                        "",
                        indent,
                        precision);
            }
            else {
                this->_write_summary_of_values<double>(
                        this->get_tree_lengths(),
                        out,
                        "",
                        indent,
                        precision);
            }
        }

        void write_summary_of_source_data(std::ostream & out,
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            std::string indent = string_util::get_indent(1);
            out.precision(precision);
            out << margin << "summary_of_tree_sources:\n"
                << margin << indent << "total_number_of_trees_sampled: " << this->get_sample_size() << "\n"
                << margin << indent << "sources:\n";
            unsigned int n_sources = this->get_number_of_sources();
            ECOEVOLITY_ASSERT(n_sources = this->source_sample_sizes_.size());
            ECOEVOLITY_ASSERT(n_sources = this->source_num_skipped_.size());
            bool source_paths_provided = false;
            if (this->source_paths_.size() > 0) {
                source_paths_provided = true;
                ECOEVOLITY_ASSERT(n_sources = this->source_paths_.size());
            }
            std::string src_indent = margin + indent + indent + "  ";
            for (unsigned int i = 0; i < this->get_number_of_sources(); ++i) {
                out << margin << indent << indent << "-\n";
                if (source_paths_provided) {
                    out << src_indent << "path: " << this->source_paths_.at(i) << "\n";
                }
                out << src_indent << "number_of_trees_skipped: "
                    << this->source_num_skipped_.at(i) << "\n"
                    << src_indent << "number_of_trees_sampled: "
                    << this->source_sample_sizes_.at(i) << "\n";
            }
        }

        bool is_a_map_node(const std::set< Split > & split_set) const {
            Split split;
            split.resize(this->get_number_of_leaves());
            for (auto s : split_set) {
                split.add_split(s);
            }
            std::vector< std::shared_ptr<NodeSamples> > node_samples = this->get_nodes_of_split(split);
            if (node_samples.empty()) {
                return false;
            }
            unsigned int map_count = node_samples.at(0)->get_sample_size();
            for (auto ns : node_samples) {
                if (ns->get_sample_size() < map_count) {
                    return false;
                }
                if (ns->get_split_set() == split_set) {
                    return true;
                }
            }
            return false;
        }

        bool is_a_map_tree(const std::set< std::set<Split> > & split_set) const {
            unsigned int map_count = this->topologies_.at(0)->get_sample_size();
            for (auto topo_samples : this->topologies_) {
                if (topo_samples->get_sample_size() < map_count) {
                    return false;
                }
                if (topo_samples->get_split_set() == split_set) {
                    return true;
                }
            }
            return false;
        }

        void write_summary_of_target_tree(
                std::ostream & out,
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            if (! this->target_tree_provided_) {
                return;
            }
            out.precision(precision);
            out << std::boolalpha;
            const std::set< std::set<Split> > & split_set = this->target_topology_;
            double target_tree_length = this->target_tree_.get_tree_length();
            double target_tree_length_percentile = percentile(
                    this->get_tree_lengths(),
                    target_tree_length);
            bool target_is_a_map_tree = this->is_a_map_tree(this->target_topology_);
            std::shared_ptr<NodeType> root_node = this->target_tree_.get_root_ptr();
            this->_annotate_nodes(this->target_topology_sample_,
                    root_node,
                    precision);
            std::string indent = string_util::get_indent(1);
            out.precision(precision);
            out << margin << "count: " << this->get_topology_count(split_set) << "\n"
                << margin << "frequency: " << this->get_topology_frequency(split_set) << "\n"
                << margin << "credibility_level: " << this->get_topology_credibility_level(split_set) << "\n"
                << margin << "is_a_map_topology: " << target_is_a_map_tree << "\n"
                << margin << "number_of_heights: " << split_set.size() << "\n"
                << margin << "number_of_heights_count: " << this->get_number_of_heights_count(split_set.size()) << "\n"
                << margin << "number_of_heights_frequency: " << this->get_number_of_heights_frequency(split_set.size()) << "\n"
                << margin << "number_of_heights_credibility_level: " << this->get_number_of_heights_credibility_level(split_set.size()) << "\n"
                << margin << "tree_length: " << target_tree_length << "\n"
                << margin << "tree_length_percentile: " << target_tree_length_percentile << "\n"
                << margin << "newick: " << root_node->to_parentheses(precision, true) << "\n";

            std::string item_margin = margin + indent + indent;

            out << margin << "splits:\n";
            out << margin << indent << "root:\n";
            for (auto param_val : this->target_split_parameters_.at(this->root_split_)) {
                double perc = percentile(
                        this->get_split(this->root_split_)->get_values(param_val.first),
                        param_val.second);
                out << item_margin << param_val.first << ": " << param_val.second << "\n"
                    << item_margin << param_val.first << "_percentile: "
                    << perc << "\n";
            }
            std::set<Split> node_split_set = this->get_target_node_split_set_of_split(
                    this->root_split_);
            out << item_margin << "node:\n";
            this->write_summary_of_node(
                    node_split_set,
                    out,
                    item_margin + indent,
                    precision);
            out << item_margin << indent << "is_a_map_node_given_split: "
                << this->is_a_map_node(node_split_set) << "\n";

            item_margin += "  ";
            out << margin << indent << "leaves:\n";
            for (auto s : this->leaf_splits_) {
                out << margin << indent << indent << "-\n";
                std::vector<unsigned int> leaf_indices = s.get_leaf_indices();
                ECOEVOLITY_ASSERT(leaf_indices.size() == 1);
                out << item_margin << "leaf_indices: [" << leaf_indices.at(0) << "]\n";
                for (auto param_val : this->target_split_parameters_.at(s)) {
                    double perc = percentile(
                            this->get_split(s)->get_values(param_val.first),
                            param_val.second);
                    out << item_margin << param_val.first << ": " << param_val.second << "\n"
                        << item_margin << param_val.first << "_percentile: "
                        << perc << "\n";
                }
            }
            out << margin << indent << "nontrivial_splits:\n";
            for (auto s_set : split_set) {
                for (auto split : s_set) {
                    if (split == this->root_split_) {
                        continue;
                    }
                    out << margin << indent << indent << "-\n";
                    this->write_summary_of_nontrivial_split(split,
                            out,
                            true,
                            false,
                            item_margin,
                            precision);
                    for (auto param_val : this->target_split_parameters_.at(split)) {
                        out << item_margin << param_val.first << ": "
                            << param_val.second << "\n";
                    }

                    if (split.get_leaf_node_count() > 2) {
                        node_split_set = this->get_target_node_split_set_of_split(split);
                        out << item_margin << "node:\n";
                        std::string node_margin = item_margin + indent;
                        this->write_summary_of_node(
                                node_split_set,
                                out,
                                node_margin,
                                precision);
                        out << node_margin << "is_a_map_node: "
                            << this->is_a_map_node(node_split_set) << "\n";
                    }
                }
            }
            std::set<Split> root_split_set;
            root_split_set.insert(this->root_split_);
            std::string h_margin = margin + indent + "  ";
            out << margin << "heights:\n";
            for (auto s_set : split_set) {
                if (s_set == root_split_set) {
                    continue;
                }
                out << margin << indent << "-\n";
                this->write_summary_of_height(s_set, out, h_margin, precision);
                out << h_margin << "height: " << this->target_heights_.at(s_set) << "\n";
            }
            if (this->target_euclidean_distances_.size() > 0) {
                out << margin << "euclidean_distance:\n";
                this->_write_summary_of_values<double>(
                        this->target_euclidean_distances_,
                        out,
                        "",
                        margin + indent,
                        precision);
            }
        }

        void write_summary_of_map_trees(std::ostream & out,
                const bool use_median_heights = false,
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            out.precision(precision);
            unsigned int map_count = this->topologies_.at(0)->get_sample_size();
            for (auto topo_samples : this->topologies_) {
                if (topo_samples->get_sample_size() < map_count) {
                    break;
                }
                out << margin << "- count: " << topo_samples->get_sample_size() << "\n"
                    << margin << "  frequency: " << (topo_samples->get_sample_size() /
                            (double)this->get_sample_size()) << "\n"
                    << margin << "  newick: " << this->to_parentheses(
                            *topo_samples,
                            use_median_heights,
                            precision) << "\n";
            }
        }

        std::shared_ptr<NodeType> get_tree(
                const TopologySamples & topo_sample,
                const bool use_median_heights = false,
                const unsigned int precision = 12) const {
            std::map<std::string, std::vector<double> > split_parameter_map;
            split_parameter_map = this->get_split(this->root_split_)->get_parameter_map();
            SampleSummary<double> root_height_summary(split_parameter_map.at("height"));
            double root_height = root_height_summary.mean();
            if (use_median_heights) {
                root_height = root_height_summary.median();
            }
            std::shared_ptr<NodeType> root_node = std::make_shared<NodeType>(
                    root_height);

            std::vector< std::shared_ptr<NodeType> > leaf_nodes;
            for (unsigned int i = 0; i < this->leaf_splits_.size(); ++i) {
                unsigned int leaf_index = this->leaf_splits_.at(i).get_leaf_indices().at(0);
                ECOEVOLITY_ASSERT(leaf_index == i);
                std::ostringstream label;
                label << this->leaf_labels_.at(i);
                split_parameter_map = this->get_split(this->leaf_splits_.at(i))->get_parameter_map();
                SampleSummary<double> leaf_height_summary(split_parameter_map.at("height"));
                double leaf_height = leaf_height_summary.mean();
                if (use_median_heights) {
                    leaf_height = leaf_height_summary.median();
                }
                std::shared_ptr<NodeType> leaf_nd = std::make_shared<NodeType>(
                            leaf_index, label.str(), leaf_height);
                leaf_nodes.push_back(leaf_nd);
                root_node->add_child(leaf_nd);
            }
            std::map< std::set<Split>, std::shared_ptr<PositiveRealParameter> > split_set_to_height_parameter_map;
            this->_process_node_splits(this->root_split_,
                    topo_sample,
                    root_node,
                    leaf_nodes,
                    split_set_to_height_parameter_map,
                    use_median_heights);
            ECOEVOLITY_ASSERT(root_node->get_leaf_node_count() == this->get_number_of_leaves());
            root_node->resize_splits(this->get_number_of_leaves());
            ECOEVOLITY_ASSERT(root_node->get_internal_node_count() == topo_sample.get_number_of_non_leaf_splits());

            this->_annotate_nodes(topo_sample,
                    root_node,
                    precision);
            return root_node;
        }

        std::shared_ptr<NodeType> get_tree(
                const std::set< std::set<Split> > & split_set,
                const bool use_median_heights = false,
                const unsigned int precision = 12) const {
            if (this->get_topology_count(split_set) < 1) {
                throw EcoevolityError("TreeSample::get_tree called with unsampled topology");
            }
            return this->get_tree(
                    this->get_topology(split_set),
                    use_median_heights,
                    precision);
        }

        std::string to_parentheses(
                const TopologySamples & topo_sample,
                const bool use_median_heights = false,
                const unsigned int precision = 12) const {

            std::shared_ptr<NodeType> root_node = this->get_tree(topo_sample,
                    use_median_heights,
                    precision);
            return root_node->to_parentheses(precision, true);
        }

        // std::string to_parentheses(const std::set< std::set<Split> > & split_set,
        //         const bool use_median_heights = false,
        //         const unsigned int precision = 12) const {
        //     if (this->get_topology_count(split_set) < 1) {
        //         throw EcoevolityError("TreeSample::to_parentheses called with unsampled topology");
        //     }
        //     return this->to_parentheses(
        //             this->get_topology_count(split_set),
        //             use_median_heights,
        //             precision);
        // }

        void write_summary_of_leaf_labels(std::ostream & out,
                const std::string & margin = "") const {
            std::string indent = string_util::get_indent(1);
            out << "leaf_label_map:\n";
            for (unsigned int i = 0; i < this->leaf_splits_.size(); ++i) {
                std::vector<unsigned int> leaf_indices = this->leaf_splits_.at(i).get_leaf_indices();
                ECOEVOLITY_ASSERT(leaf_indices.size() == 1);
                out << indent << leaf_indices.at(0) << ": " << this->leaf_labels_.at(i) << "\n";
            }
        }

        void write_summary(std::ostream & out,
                const bool use_median_heights = false,
                const double min_freq_for_asdsf = 0.1,
                const std::string & margin = "",
                const unsigned int precision = 12) const {
            std::string indent = string_util::get_indent(1);
            out.precision(precision);
            out << "---\n";
            this->write_summary_of_leaf_labels(out, margin);
            this->write_summary_of_source_data(out, margin, precision);
            if (this->get_number_of_sources() > 1) {
                SampleSummarizer<double> sdsf_summary =
                    this->get_summary_of_split_freq_std_devs(min_freq_for_asdsf);
                out << "summary_of_split_freq_std_deviations:\n"
                    << indent << "min_frequency: " << min_freq_for_asdsf << "\n"
                    << indent << "n: " << sdsf_summary.sample_size() << "\n"
                    << indent << "mean: " << sdsf_summary.mean() << "\n"
                    << indent << "max: " << sdsf_summary.max() << "\n";
            }
            this->write_summary_of_tree_lengths(out, margin, precision);
            out << "summary_of_map_topologies:\n";
            this->write_summary_of_map_trees(out,
                    use_median_heights,
                    indent,
                    precision);
            if (this->target_tree_provided_) {
                out << "summary_of_target_tree:\n";
                this->write_summary_of_target_tree(out,
                    indent,
                    precision);
            }
            this->write_summary_of_topologies(out, "", precision);
            this->write_summary_of_heights(out, "", precision);
            this->write_summary_of_splits(out, "", precision);
            this->write_summary_of_all_numbers_of_heights(out, "", precision);
        }

        void write_to_nexus(
                std::vector< std::set< std::set<Split> > >::const_iterator split_sets_start,
                std::vector< std::set< std::set<Split> > >::const_iterator split_sets_end,
                std::ostream & out,
                const bool use_median_heights = false,
                const unsigned int precision = 12) const {
            std::vector< std::shared_ptr<NodeType> > trees;
            while (split_sets_start != split_sets_end) {
                std::shared_ptr<NodeType> t = this->get_tree(
                        *(split_sets_start++),
                        use_median_heights,
                        precision);
                trees.push_back(t);
            }
            this->_write_to_nexus(
                    trees,
                    out,
                    precision);
        }

        void write_target_tree_to_nexus(
                std::ostream & out,
                const unsigned int precision = 12) const {
            std::shared_ptr<NodeType> root_node = this->target_tree_.get_root_ptr();
            this->_annotate_nodes(this->target_topology_sample_,
                    root_node,
                    precision);
            std::vector< std::shared_ptr<NodeType> > target_tree_vec = {root_node};
            this->_write_to_nexus(
                    target_tree_vec,
                    out,
                    precision);
        }

        void write_map_trees_to_nexus(
                std::ostream & out,
                const bool use_median_heights = false,
                const unsigned int precision = 12) const {
            std::vector< std::shared_ptr<NodeType> > trees;
            unsigned int map_count = this->topologies_.at(0)->get_sample_size();
            for (auto topo_samples : this->topologies_) {
                if (topo_samples->get_sample_size() < map_count) {
                    break;
                }
                std::shared_ptr<NodeType> t = this->get_tree(
                        *topo_samples,
                        use_median_heights,
                        precision);
                trees.push_back(t);
            }
            this->_write_to_nexus(
                    trees,
                    out,
                    precision);
        }

        void write_summary_of_merged_target_heights(std::ostream & out,
                const std::string & margin = "",
                const unsigned int precision = 12) {
            if (! this->target_tree_provided_) {
                return;
            }
            unsigned int nheights = this->target_tree_.get_number_of_node_heights();
            if (nheights < 2) {
                // No heights to merge
                return;
            }
            std::string indent = string_util::get_indent(1);
            out.precision(precision);

            std::vector<unsigned int> sizes_of_mapped_polytomies;
            unsigned int number_of_nodes_merged_up;
            std::map< unsigned int, std::set< std::set<Split> > > height_index_to_node_set;

            std::string item_margin = margin + indent + "  ";
            out << margin << "merged_target_heights:\n";

            for (unsigned int ht_idx = 0; ht_idx < (nheights - 1); ++ht_idx) {
                this->target_tree_.store_state();

                double younger_height = this->target_tree_.get_height(ht_idx);
                double older_height = this->target_tree_.get_height(ht_idx + 1);
                unsigned int n_nodes_mapped = this->target_tree_.get_mapped_node_count(ht_idx);
                unsigned int older_n_nodes_mapped = this->target_tree_.get_mapped_node_count(ht_idx + 1);

                this->target_tree_.merge_node_height_up(ht_idx,
                        sizes_of_mapped_polytomies,
                        number_of_nodes_merged_up,
                        true);
                ECOEVOLITY_ASSERT(this->target_tree_.get_height(ht_idx) == older_height);
                ECOEVOLITY_ASSERT(this->target_tree_.get_number_of_node_heights() == (nheights - 1));
                unsigned int post_merge_num_nodes = this->target_tree_.get_mapped_node_count(ht_idx);
                height_index_to_node_set.clear();
                this->target_tree_.store_nodes_by_height_index(
                        height_index_to_node_set,
                        false); // Don't update split size (i.e., # of leaves)
                ECOEVOLITY_ASSERT(height_index_to_node_set.at(ht_idx).size() == post_merge_num_nodes);
                double height_freq = this->get_node_height_frequency(height_index_to_node_set.at(ht_idx));
                double height_count = this->get_node_height_count(height_index_to_node_set.at(ht_idx));

                out << margin << indent << "-\n";
                out << item_margin << "younger_height_index: " << ht_idx << "\n";
                out << item_margin << "younger_height: " << younger_height << "\n";
                out << item_margin << "older_height: " << older_height << "\n";
                out << item_margin << "younger_height_number_of_nodes: " << n_nodes_mapped << "\n";
                out << item_margin << "older_height_number_of_nodes: " << older_n_nodes_mapped << "\n";
                out << item_margin << "merged_height_number_of_nodes: " << post_merge_num_nodes << "\n";
                out << item_margin << "count: " << height_count << "\n";
                out << item_margin << "frequency: " << height_freq << "\n";

                this->target_tree_.restore_state();
                ECOEVOLITY_ASSERT(this->target_tree_.get_number_of_node_heights() == nheights);
            }
        }
};

} // treesum

#endif
