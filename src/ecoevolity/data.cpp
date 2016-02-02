#include "data.hpp"

BiallelicData::BiallelicData(
        const std::string path, 
        const char population_name_delimiter,
        const bool population_name_is_prefix) {
    char pop_name_delimiter = population_name_delimiter;
    if (population_name_delimiter == '_') {
        pop_name_delimiter = ' ';
    }
    std::cout << "population name delimiter: " << population_name_delimiter << std::endl;
    std::cout << "pop name delimiter: " << pop_name_delimiter << std::endl;
    std::cout << "pop name is prefix: " << population_name_is_prefix << std::endl;

    MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);
    nexusReader.ReadFilepath(path.c_str(), MultiFormatReader::NEXUS_FORMAT);
    
    unsigned int numTaxaBlocks = nexusReader.GetNumTaxaBlocks();
    std::cout << numTaxaBlocks << " TAXA block(s) read.\n";
    for (unsigned int i = 0; i < numTaxaBlocks; ++i) {
        NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(i);
        std::string taxaBlockTitle = taxaBlock->GetTitle();
        unsigned int num_taxa = taxaBlock->GetNTax();
        std::cout << "Taxa block index " << i << " has the Title \"" << taxaBlockTitle << "\"\n";
        std::cout << "Taxa block index " << i << " has the " << num_taxa << " taxa\n";

        std::cout << "Taxa block index " << i << " has the labesl:\n";
        for (unsigned int j = 0; j < num_taxa; ++j) {
            NxsString taxon_label = taxaBlock->GetTaxonLabel(j);
            std::cout << taxon_label << std::endl;
            std::vector<std::string> taxon_label_elements = split(taxon_label, pop_name_delimiter);
            assert (! taxon_label_elements.empty());
            std::string taxon_name = taxon_label_elements.front();
            if (! population_name_is_prefix) {
                taxon_name = taxon_label_elements.back();
            }
            std::cout << taxon_name << std::endl;
        }
    
        const unsigned int nCharBlocks = nexusReader.GetNumCharactersBlocks(taxaBlock);
        std::cout  <<  nCharBlocks << " CHARACTERS/DATA block(s) refer to this TAXA block\n";
        for (unsigned int j = 0; j < nCharBlocks; ++j) {
            const NxsCharactersBlock * charBlock = nexusReader.GetCharactersBlock(taxaBlock, j);
            std::string charBlockTitle = charBlock->GetTitle();
            //unsigned int dtype = charBlock->GetDataType();
            NxsCharactersBlock::DataTypesEnum dtype;
            dtype = charBlock->GetDataType();
            bool data_is_standard;
            data_is_standard = (dtype == NxsCharactersBlock::DataTypesEnum::standard); 
            std::cout << std::boolalpha; // write booleans as true/false
            std::cout << "Char block index " << j << " has the Title \"" << charBlockTitle << "\"" << std::endl;
            std::cout << "Char block index " << j << " has data type: \"" << (int)dtype << "\"" << std::endl;
            std::cout << "Char block index " << j << " data type \"" << (int)dtype << "\" == standard: " << data_is_standard << std::endl;
        }
    }
}


std::vector<unsigned int> BiallelicData::get_number_of_red_alleles(unsigned int pattern_index) {
    std::vector<unsigned int> v (2, 0);
    return v;
}

std::vector<unsigned int> BiallelicData::get_number_of_alleles(unsigned int pattern_index) {
    std::vector<unsigned int> v (2, 0);
    return v;
}

unsigned int BiallelicData::get_pattern_weight(unsigned int pattern_index) {
    return 0;
}

unsigned int BiallelicData::get_number_of_patterns() {
    return 0;
}

unsigned int BiallelicData::get_number_of_taxa() {
    return 0;
}

void BiallelicData::remove_constant_patterns() {
    std::cout << "";
}

void BiallelicData::remove_patterns_with_missing_taxa() {
    std::cout << "";
}
