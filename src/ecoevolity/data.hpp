#ifndef DATA_HPP
#define DATA_HPP

#include <iostream>
#include <ncl/nxsmultiformat.h>

/**
 * Class for storing biallelic site patterns.
 *
 */
class BiallelicData {
    public:
        // Constructor
        BiallelicData(const std::string path,
                const char population_name_delimiter = '_',
                const bool population_name_is_prefix = true);
        // Destructor
        // ~BiallelicData();
        
        //Methods
        std::vector<unsigned int> get_num_red_alleles(unsigned int pattern_index);
        std::vector<unsigned int> get_num_alleles(unsigned int pattern_index);
        unsigned int get_pattern_weight(unsigned int pattern_index);
        unsigned int get_number_of_patterns();
        unsigned int get_number_of_taxa();

        void remove_constant_patterns();
        void remove_patterns_with_missing_taxa();

    private:
        std::vector< std::vector<unsigned int> > num_red_alleles_;
        std::vector< std::vector<unsigned int> > num_alleles_;
        std::vector<unsigned int> pattern_weights_;
        std::vector<std::string> taxon_names_;
};

#endif

