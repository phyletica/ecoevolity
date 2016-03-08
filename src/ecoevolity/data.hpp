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
#include <assert.h>
#include <unordered_map>
#include <map>
#include <string>
#include <vector>
#include <ncl/nxsmultiformat.h>

#include "util.hpp"
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
        BiallelicData(const std::string path,
                const char population_name_delimiter = '_',
                const bool population_name_is_prefix = true,
                const bool genotypes_are_diploid = true,
                const bool markers_are_dominant = false,
                const bool validate = true);
        // Destructor
        // ~BiallelicData();
        
        //Methods
        const std::vector<unsigned int>& get_red_allele_counts(unsigned int pattern_index) const;
        const std::vector<unsigned int>& get_allele_counts(unsigned int pattern_index) const;

        const unsigned int& get_red_allele_count(unsigned int pattern_index,
                unsigned int population_index) const;
        const unsigned int& get_allele_count(unsigned int pattern_index,
                unsigned int population_index) const;

        const unsigned int& get_pattern_weight(unsigned int pattern_index) const;
        const unsigned int& get_max_allele_count(unsigned int population_index) const;
        const std::vector<unsigned int>& get_max_allele_counts() const;
        const unsigned int& get_population_index(std::string population_label) const;
        const unsigned int& get_population_index_from_seq_label(std::string seq_label) const;
        const std::string& get_population_label(unsigned int population_index) const;
        const std::vector<std::string>& get_sequence_labels(unsigned int population_index) const;
        const std::string& get_path() const;

        unsigned int get_number_of_patterns() const;
        unsigned int get_number_of_populations() const;

        const bool& markers_are_dominant() const;
        const bool& genotypes_are_diploid() const;

        const bool& has_constant_patterns() const;
        const bool& has_missing_population_patterns() const;
        const bool& has_mirrored_patterns() const;
        const bool& patterns_are_folded() const;

        int get_pattern_index(
                const std::vector<unsigned int> red_allele_counts,
                const std::vector<unsigned int> allele_counts) const;

        const std::vector< std::vector<unsigned int> > get_mirrored_pattern(unsigned int pattern_index) const;

        unsigned int remove_constant_patterns(const bool validate = true);
        unsigned int remove_missing_population_patterns(const bool validate = true);
        unsigned int fold_patterns(const bool validate = true);

        unsigned int get_number_of_constant_sites_removed() const;
        unsigned int get_number_of_missing_sites_removed() const;

        void validate() const;

    private:
        unsigned int number_of_constant_sites_removed_ = 0;
        unsigned int number_of_missing_sites_removed_ = 0;
        bool markers_are_dominant_ = true;
        bool genotypes_are_diploid_ = true;
        bool has_missing_population_patterns_ = false;
        bool has_constant_patterns_ = false;
        bool has_mirrored_patterns_ = true;
        bool patterns_are_folded_ = false;
        std::string path_;
        std::vector< std::vector<unsigned int> > red_allele_counts_;
        std::vector< std::vector<unsigned int> > allele_counts_;
        std::vector<unsigned int> pattern_weights_;
        std::vector<unsigned int> max_allele_counts_;
        std::vector<std::string> population_labels_;
        std::vector< std::vector<std::string> > sequence_labels_;
        std::unordered_map<std::string, std::string> seq_label_to_pop_label_map_;
        std::unordered_map<std::string, unsigned int> pop_label_to_index_map_;

        //Methods
        void remove_pattern(unsigned int pattern_index);
        void update_has_missing_population_patterns();
        void update_has_constant_patterns();
        void update_has_mirrored_patterns();
        void update_patterns_are_folded();
        void update_pattern_booleans();
        int remove_first_constant_pattern();
        int remove_first_missing_population_pattern();
        int fold_first_mirrored_pattern();
        void update_max_allele_counts();
};

#endif

