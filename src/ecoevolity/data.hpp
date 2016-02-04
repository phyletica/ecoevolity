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
                const bool markers_are_dominant = false);
        // Destructor
        // ~BiallelicData();
        
        //Methods
        std::vector<unsigned int> get_red_allele_counts(unsigned int pattern_index) const;
        std::vector<unsigned int> get_allele_counts(unsigned int pattern_index) const;
        unsigned int get_pattern_weight(unsigned int pattern_index) const;
        unsigned int get_number_of_patterns() const;
        unsigned int get_number_of_populations() const;
        unsigned int get_population_index(std::string population_label) const;
        unsigned int get_population_index_from_seq_label(std::string seq_label) const;
        int get_pattern_index(
                const std::vector<unsigned int> red_allele_counts,
                const std::vector<unsigned int> allele_counts) const;


        void remove_constant_patterns();
        void remove_patterns_with_missing_taxa();

    private:
        bool markers_are_dominant_ = true;
        bool genotypes_are_diploid_ = true;
        std::vector< std::vector<unsigned int> > red_allele_counts_;
        std::vector< std::vector<unsigned int> > allele_counts_;
        std::vector<unsigned int> pattern_weights_;
        std::vector<std::string> population_labels_;
        std::vector< std::vector<std::string> > sequence_labels_;
        std::unordered_map<std::string, std::string> seq_label_to_pop_label_map_;
        std::unordered_map<std::string, unsigned int> pop_label_to_index_map_;
};

#endif

