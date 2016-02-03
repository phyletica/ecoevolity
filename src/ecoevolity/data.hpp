#ifndef ECOEVOLITY_DATA_HPP
#define ECOEVOLITY_DATA_HPP

#include <iostream>
#include <assert.h>
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
        std::vector<unsigned int> get_number_of_red_alleles(unsigned int pattern_index);
        std::vector<unsigned int> get_number_of_alleles(unsigned int pattern_index);
        unsigned int get_pattern_weight(unsigned int pattern_index);
        unsigned int get_number_of_patterns();
        unsigned int get_number_of_populations();

        void remove_constant_patterns();
        void remove_patterns_with_missing_taxa();

    private:
        bool markers_are_dominant_ = true;
        bool genotypes_are_diploid_ = true;
        std::vector< std::vector<unsigned int> > number_of_red_alleles_;
        std::vector< std::vector<unsigned int> > number_of_alleles_;
        std::vector<unsigned int> pattern_weights_;
        std::vector<std::string> population_labels_;
        std::vector< std::vector<std::string> > sequence_labels_;
};

#endif

