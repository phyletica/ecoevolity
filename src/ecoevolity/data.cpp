#include "data.hpp"

BiallelicData::BiallelicData(
        const std::string path, 
        const char population_name_delimiter,
        const bool population_name_is_prefix,
        const bool markers_are_dominant) {
    char pop_name_delimiter = population_name_delimiter;
    if (population_name_delimiter == '_') {
        pop_name_delimiter = ' ';
    }
    this->markers_are_dominant_ = markers_are_dominant;

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
    const unsigned int num_char_blocks = nexus_reader.GetNumCharactersBlocks(taxa_block);
    if (num_char_blocks < 1) {
        throw EcoevolityParsingError("No character block found", path, 0);
    }
    if (num_char_blocks > 1) {
        throw EcoevolityParsingError("More than one character block found", path, 0);
    }

    NxsCharactersBlock * char_block = nexus_reader.GetCharactersBlock(taxa_block, 0);
    std::string char_block_title = char_block->GetTitle();

    unsigned int num_chars = char_block->GetNCharTotal();
    unsigned int num_taxa = char_block->GetNTaxTotal();
    NxsCharactersBlock::DataTypesEnum data_type;
    data_type = char_block->GetDataType();

    ECOEVOLITY_DEBUG(
    std::cerr << "Char block " << char_block_title << " has data type: " <<
            (int)data_type << std::endl;
    )

    for (unsigned int taxon_idx = 0; taxon_idx < num_taxa; ++taxon_idx) {
        NxsString seq_label = char_block->GetTaxonLabel(taxon_idx);
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

    std::vector<const NxsDiscreteDatatypeMapper *> data_type_mappers = char_block->GetAllDatatypeMappers();
    if (data_type_mappers.size() < 1) {
        throw EcoevolityParsingError("No character encoding found", path, 0);
    }
    if (data_type_mappers.size() > 1) {
        throw EcoevolityParsingError("More than one character encoding (i.e., mixed data types) found", path, 0);
    }
    const NxsDiscreteDatatypeMapper * data_type_mapper = data_type_mappers[0];

    if (data_type == NxsCharactersBlock::DataTypesEnum::standard) {
        const NxsDiscreteStateCell highest_state_code = data_type_mapper->GetHighestStateCode();
        if (highest_state_code == 1) {
            this->genotypes_are_diploid_ = false;
        }
        else if (highest_state_code == 2) {
            this->genotypes_are_diploid_ = true;
        }
        else {
            throw EcoevolityParsingError("More than 3 character state codes found", path, 0);
        }

        if (! this->markers_are_dominant_) {
            for (unsigned int site_idx = 0; site_idx < num_chars; ++site_idx) {
                for (unsigned int taxon_idx = 0; taxon_idx < num_taxa; ++taxon_idx) {
                    const NxsDiscreteStateCell state_code = char_block->GetInternalRepresentation(taxon_idx, site_idx);
                    if (state_code >= 0) {
                        const unsigned int num_states = char_block->GetNumStates(taxon_idx, site_idx);
                        if (num_states > 1) {
                            throw EcoevolityInvalidCharacterError(
                                    "Invalid polymorphic character",
                                    path,
                                    char_block->GetTaxonLabel(taxon_idx),
                                    site_idx);
                        }
                        std::cout << state_code << std::endl;
                    }
                }
            }
        }
    }

    nexus_reader.DeleteBlocksFromFactories();
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

unsigned int BiallelicData::get_number_of_populations() {
    return this->population_labels_.size();
}

void BiallelicData::remove_constant_patterns() {
    std::cout << "";
}

void BiallelicData::remove_patterns_with_missing_taxa() {
    std::cout << "";
}
