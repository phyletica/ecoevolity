#include "data.hpp"

BiallelicData::BiallelicData(
        const std::string path, 
        const char population_name_delimiter,
        const bool population_name_is_prefix) {
    char pop_name_delimiter = population_name_delimiter;
    if (population_name_delimiter == '_') {
        pop_name_delimiter = ' ';
    }

    MultiFormatReader nexus_reader(-1, NxsReader::WARNINGS_TO_STDERR);
    nexus_reader.ReadFilepath(path.c_str(), MultiFormatReader::NEXUS_FORMAT);
    
    unsigned int num_taxa_blocks = nexus_reader.GetNumTaxaBlocks();
    if (num_taxa_blocks < 1) {
        throw EcoevolityParsingError("No taxa block found", path, 0);
    }
    if (num_taxa_blocks > 1) {
        throw EcoevolityParsingError("More than one taxa block found", path, 0);
    }

    NxsTaxaBlock * taxa_block = nexus_reader.GetTaxaBlock(0);
    std::string taxa_block_title = taxa_block->GetTitle();
    unsigned int num_taxa = taxa_block->GetNTax();
    for (unsigned int taxon_idx = 0; taxon_idx < num_taxa; ++taxon_idx) {
        NxsString seq_label = taxa_block->GetTaxonLabel(taxon_idx);
        std::vector<std::string> seq_label_elements = split(seq_label, pop_name_delimiter);
        ECOEVOLITY_ASSERT(! seq_label_elements.empty());
        std::string pop_label = seq_label_elements.front();
        if (! population_name_is_prefix) {
            pop_label = seq_label_elements.back();
        }
        bool pop_label_found = false;
        for (unsigned int pop_label_idx = 0; pop_label_idx < this->population_labels_.size(); pop_label_idx++) {
            if (this->population_labels_[pop_label_idx] == pop_label) {
                assert(! this->sequence_labels_[pop_label_idx].empty());
                this->sequence_labels_[pop_label_idx].push_back(seq_label);
                pop_label_found = true;
            }
        }
        if (! pop_label_found) {
            this->population_labels_.push_back(pop_label);
            std::vector<std::string> tmp_label_vector = {seq_label};
            this->sequence_labels_.push_back(tmp_label_vector);
        }
    }
        
    ECOEVOLITY_DEBUG(
    std::cerr << "this->populations_labels_:" << std::endl;
    unsigned int pop_idx = 0;
    for (auto p_label: this->population_labels_) {
        std::cerr << p_label << std::endl;
        for (auto s_label: this->sequence_labels_[pop_idx]) {
            std::cerr << "\t" << s_label << std::endl;
        }
        pop_idx += 1;
    }
    )

    const unsigned int num_char_blocks = nexus_reader.GetNumCharactersBlocks(taxa_block);
    if (num_char_blocks < 1) {
        throw EcoevolityParsingError("No character block found", path, 0);
    }
    if (num_char_blocks > 1) {
        throw EcoevolityParsingError("More than one character block found", path, 0);
    }
    const NxsCharactersBlock * char_block = nexus_reader.GetCharactersBlock(taxa_block, 0);
    std::string char_block_title = char_block->GetTitle();
    NxsCharactersBlock::DataTypesEnum dtype;
    dtype = char_block->GetDataType();
    bool data_is_standard;
    data_is_standard = (dtype == NxsCharactersBlock::DataTypesEnum::standard); 

    ECOEVOLITY_DEBUG(
    std::cerr << std::boolalpha; // write booleans as true/false
    std::cerr << "Char block " << char_block_title << " has data type: " << (int)dtype << std::endl;
    std::cerr << "Data type " << (int)dtype << " == standard: " << data_is_standard << std::endl;
    )
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
