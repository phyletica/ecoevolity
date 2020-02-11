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

#ifndef ECOEVOLITY_DATA_HPP
#define ECOEVOLITY_DATA_HPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <unordered_map>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <ncl/nxsmultiformat.h>

#include "yaml_util.hpp"

#include "string_util.hpp"
#include "debug.hpp"
#include "assert.hpp"
#include "error.hpp"

/**
 * Class for storing biallelic site patterns.
 *
 */
class BiallelicData {
    public:
        // Constructor
        BiallelicData() { }
        BiallelicData(std::string path,
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool validate = true,
                bool store_seq_loci_info = false);
        void init(const std::string path,
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool validate = true,
                bool store_seq_loci_info = false);
        void init_from_yaml_stream(
                std::istream& stream,
                const std::string& path,
                bool validate = true);
        void init_from_yaml_file(
                const std::string& path,
                bool validate = true);

        // Destructor
        // ~BiallelicData();
        
        //Methods
        BiallelicData get_empty_copy() const;
        const std::vector<unsigned int>& get_red_allele_counts(unsigned int pattern_index) const;
        const std::vector<unsigned int>& get_allele_counts(unsigned int pattern_index) const;

        unsigned int get_red_allele_count(unsigned int pattern_index,
                unsigned int population_index) const;
        unsigned int get_allele_count(unsigned int pattern_index,
                unsigned int population_index) const;

        unsigned int get_pattern_weight(unsigned int pattern_index) const;
        unsigned int get_max_allele_count(unsigned int population_index) const;
        const std::vector<unsigned int>& get_max_allele_counts() const;
        unsigned int get_population_index(std::string population_label) const;
        unsigned int get_population_index_from_seq_label(std::string seq_label) const;
        const std::string& get_population_label(unsigned int population_index) const;
        const std::vector<std::string>& get_population_labels() const {
            return this->population_labels_;
        }
        const std::vector<std::string>& get_sequence_labels(unsigned int population_index) const;
        const std::string& get_path() const;

        unsigned int get_number_of_patterns() const;
        unsigned int get_number_of_sites() const;
        unsigned int get_number_of_variable_sites() const;
        unsigned int get_number_of_populations() const;

        bool markers_are_dominant() const;
        bool genotypes_are_diploid() const;

        bool has_constant_patterns() const;
        bool has_missing_population_patterns() const;
        bool has_mirrored_patterns() const;
        bool patterns_are_folded() const;
        bool has_recoded_triallelic_sites() const;

        static bool pattern_is_constant(
                const std::vector<unsigned int>& red_allele_counts,
                const std::vector<unsigned int>& allele_counts);

        void get_pattern_index(
                bool& was_found,
                unsigned int& pattern_index,
                const std::vector<unsigned int>& red_allele_counts,
                const std::vector<unsigned int>& allele_counts) const;

        const std::vector< std::vector<unsigned int> > get_mirrored_pattern(unsigned int pattern_index) const;

        unsigned int remove_constant_patterns(const bool validate = true);
        unsigned int remove_missing_population_patterns(const bool validate = true);
        unsigned int fold_patterns(const bool validate = true);

        unsigned int get_number_of_constant_sites_removed() const;
        unsigned int get_number_of_constant_red_sites_removed() const;
        unsigned int get_number_of_constant_green_sites_removed() const;
        unsigned int get_number_of_missing_sites_removed() const;
        unsigned int get_number_of_triallelic_sites_recoded() const;

        double get_proportion_of_red_alleles() const;
        double get_proportion_1() const {
            return this->get_proportion_of_red_alleles();
        }
        void get_empirical_u_v_rates(double& u, double& v) const;

        void validate() const;

        bool add_site(
                const std::vector<unsigned int>& red_allele_counts,
                const std::vector<unsigned int>& allele_counts,
                bool filtering_constant_patterns = true,
                bool end_of_locus = false);

        void update_pattern_booleans();
        void update_max_allele_counts();

        std::vector< std::vector<std::string> > get_alignment() const;
        void write_nexus(
                std::ostream& out,
                char population_name_delimiter) const;
        void write_charsets(
                std::ostream& out) const;
        void write_alignment(
                std::ostream& out,
                char population_name_delimiter) const;
        void write_yaml(std::ostream& out) const;

        void write_summary(
                std::ostream& out,
                unsigned int indent_level = 0) const;

        const std::vector< std::vector<unsigned int> >& get_red_allele_count_matrix() const {
            return this->red_allele_counts_;
        }
        const std::vector< std::vector<unsigned int> >& get_allele_count_matrix() const {
            return this->allele_counts_;
        }
        const std::vector<unsigned int>& get_pattern_weights() const {
            return this->pattern_weights_;
        }
        const std::vector<unsigned int>& get_contiguous_pattern_indices() const {
            return this->contiguous_pattern_indices_;
        }
        const std::vector<unsigned int>& get_locus_end_indices() const {
            return this->locus_end_indices_;
        }
        bool has_seq_loci_info() const {
            return this->storing_seq_loci_info_;
        }
        void start_storing_seq_loci_info() {
            if (! this->appendable_) {
                throw EcoevolityBiallelicDataError(
                        "Cannot toggle storing of loci info for a parsed dataset",
                        this->path_);
            }
            this->contiguous_pattern_indices_.clear();
            this->locus_end_indices_.clear();
            this->storing_seq_loci_info_ = true;
        }
        void stop_storing_seq_loci_info() {
            if (! this->appendable_) {
                throw EcoevolityBiallelicDataError(
                        "Cannot toggle storing of loci info for a parsed dataset",
                        this->path_);
            }
            this->contiguous_pattern_indices_.clear();
            this->locus_end_indices_.clear();
            this->storing_seq_loci_info_ = false;
        }
        unsigned int get_pattern_index_for_site(unsigned int site_index) const {
            return this->contiguous_pattern_indices_.at(site_index);
        }

        /**
         * Get the unique allele counts.
         *
         * Normally, the BiallelicData class works with allele patterns and
         * their weights (the number of sites that have the pattern), which
         * include the red allele counts and total allele counts from each
         * population. Here what we want is only the unique total allele counts
         * and their weights (so, ignoring the red allele counts).
         *
         * Returns a map of unique allele count vectors to
         * the weights of those allele counts.
         *
         */
        std::map< std::vector<unsigned int>, unsigned int >
        get_unique_allele_counts() const;


    private:
        unsigned int number_of_constant_red_sites_removed_ = 0;
        unsigned int number_of_constant_green_sites_removed_ = 0;
        unsigned int number_of_missing_sites_removed_ = 0;
        unsigned int number_of_triallelic_sites_recoded_ = 0;
        bool markers_are_dominant_ = true;
        bool genotypes_are_diploid_ = true;
        bool has_missing_population_patterns_ = false;
        bool has_constant_patterns_ = false;
        bool has_mirrored_patterns_ = true;
        bool patterns_are_folded_ = false;
        bool appendable_ = false;
        std::string path_ = "";
        std::vector< std::vector<unsigned int> > red_allele_counts_;
        std::vector< std::vector<unsigned int> > allele_counts_;
        std::vector<unsigned int> pattern_weights_;
        std::vector<unsigned int> max_allele_counts_;
        std::vector<std::string> population_labels_;
        std::vector< std::vector<std::string> > sequence_labels_;
        std::unordered_map<std::string, std::string> seq_label_to_pop_label_map_;
        std::unordered_map<std::string, unsigned int> pop_label_to_index_map_;

        bool storing_seq_loci_info_ = false;
        std::vector<unsigned int> locus_end_indices_;
        std::vector<unsigned int> contiguous_pattern_indices_;

        //Methods
        void remove_pattern(unsigned int pattern_index);
        void update_has_missing_population_patterns();
        void update_has_constant_patterns();
        void update_has_mirrored_patterns();
        void update_patterns_are_folded();
        void remove_first_constant_pattern(
                bool& was_removed,
                unsigned int& removed_index);
        void remove_first_missing_population_pattern(
                bool& was_removed,
                unsigned int& removed_index);
        void fold_first_mirrored_pattern(
                bool& was_folded,
                unsigned int& folded_index);
        
        void parse_yaml_data(std::istream& yaml_stream, bool validate = true);
        void parse_yaml_top_level(const YAML::Node& top_level_node);
        void parse_yaml_marker_dominance(const YAML::Node& node);
        void parse_yaml_population_labels(const YAML::Node& node);
        void parse_yaml_pattern_weights(const YAML::Node& node);
        void parse_yaml_all_allele_count_patterns(const YAML::Node& node);
        void parse_yaml_allele_count_pattern(const YAML::Node& node,
                std::vector<unsigned int> & red_allele_counts,
                std::vector<unsigned int> & allele_counts);
        void parse_yaml_allele_count(const YAML::Node& node,
                std::vector<unsigned int> & red_allele_counts,
                std::vector<unsigned int> & allele_counts);
};

#endif
