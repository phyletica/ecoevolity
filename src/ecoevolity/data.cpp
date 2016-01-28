#include "data.hpp"

#define pass ((void)0)

BiallelicData::BiallelicData(const char * path) {
    MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);
    nexusReader.ReadFilepath(path, MultiFormatReader::NEXUS_FORMAT);
    
    int numTaxaBlocks = nexusReader.GetNumTaxaBlocks();
    std::cout << numTaxaBlocks << " TAXA block(s) read.\n";
    for (int i = 0; i < numTaxaBlocks; ++i) {
        NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(i);
        std::string taxaBlockTitle = taxaBlock->GetTitle();
        std::cout << "Taxa block index " << i << " has the Title \"" << taxaBlockTitle << "\"\n";
    
        const unsigned nCharBlocks = nexusReader.GetNumCharactersBlocks(taxaBlock);
        std::cout  <<  nCharBlocks << " CHARACTERS/DATA block(s) refer to this TAXA block\n";
        for (unsigned j = 0; j < nCharBlocks; ++j) {
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
            std::cout << "Char block indix " << j << " data type \"" << (int)dtype << "\" == standard: " << data_is_standard << std::endl;
        }
    }
}


std::vector<unsigned int> BiallelicData::get_num_red_alleles(unsigned int pattern_index) {
    pass;
}

std::vector<unsigned int> BiallelicData::get_num_alleles(unsigned int pattern_index) {
    pass;
}

unsigned int BiallelicData::get_pattern_weight(unsigned int pattern_index) {
    pass;
}

unsigned int BiallelicData::get_number_of_patterns() {
    pass;
}

unsigned int BiallelicData::get_number_of_taxa() {
    pass;
}

void BiallelicData::remove_constant_patterns() {
    pass;
}

void BiallelicData::remove_patterns_with_missing_taxa() {
    pass;
}
