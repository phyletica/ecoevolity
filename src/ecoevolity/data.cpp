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

#include "data.hpp"

BiallelicData::BiallelicData(
        std::string path, 
        char population_name_delimiter,
        bool population_name_is_prefix,
        bool genotypes_are_diploid,
        bool markers_are_dominant,
        bool validate) {
    this->init(path,
               population_name_delimiter,
               population_name_is_prefix,
               genotypes_are_diploid,
               markers_are_dominant,
               validate);
}

void BiallelicData::init(
        std::string path, 
        char population_name_delimiter,
        bool population_name_is_prefix,
        bool genotypes_are_diploid,
        bool markers_are_dominant,
        bool validate) {
    char pop_name_delimiter = population_name_delimiter;
    if (population_name_delimiter == '_') {
        pop_name_delimiter = ' ';
    }
    this->genotypes_are_diploid_ = genotypes_are_diploid;
    this->markers_are_dominant_ = markers_are_dominant;
    this->path_ = path;

    if ((this->markers_are_dominant_) and (this->genotypes_are_diploid_)) {
        throw EcoevolityBiallelicDataError(
                "Dominant markers must be coded as haploid (i.e., 0/1)",
                this->path_);
    }

    MultiFormatReader nexus_reader(-1, NxsReader::WARNINGS_TO_STDERR);
    nexus_reader.ReadFilepath(this->path_.c_str(), MultiFormatReader::NEXUS_FORMAT);
    
    unsigned int num_taxa_blocks = nexus_reader.GetNumTaxaBlocks();
    if (num_taxa_blocks < 1) {
        throw EcoevolityParsingError("No taxa block found", this->path_, 0);
    }
    if (num_taxa_blocks > 1) {
        throw EcoevolityParsingError("More than one taxa block found", this->path_, 0);
    }

    NxsTaxaBlock * taxa_block = nexus_reader.GetTaxaBlock(0);
    std::string taxa_block_title = taxa_block->GetTitle();
    const unsigned int num_char_blocks = nexus_reader.GetNumCharactersBlocks(taxa_block);
    if (num_char_blocks < 1) {
        throw EcoevolityParsingError("No character block found", this->path_, 0);
    }
    if (num_char_blocks > 1) {
        throw EcoevolityParsingError("More than one character block found", this->path_, 0);
    }

    NxsCharactersBlock * char_block = nexus_reader.GetCharactersBlock(taxa_block, 0);
    std::string char_block_title = char_block->GetTitle();

    unsigned int num_chars = char_block->GetNCharTotal();
    unsigned int num_taxa = char_block->GetNTaxTotal();
    NxsCharactersBlock::DataTypesEnum data_type;
    data_type = char_block->GetDataType();

    // ECOEVOLITY_DEBUG(
    // std::cerr << "Char block " << char_block_title << " has data type: " <<
    //         (int)data_type << std::endl;
    // )

    for (unsigned int taxon_idx = 0; taxon_idx < num_taxa; ++taxon_idx) {
        NxsString seq_label = char_block->GetTaxonLabel(taxon_idx);
        std::vector<std::string> seq_label_elements = string_util::split(
                seq_label,
                pop_name_delimiter);
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
            this->pop_label_to_index_map_[pop_label] = this->population_labels_.size() - 1;
            std::vector<std::string> tmp_label_vector = {seq_label};
            this->sequence_labels_.push_back(tmp_label_vector);
        }
        this->seq_label_to_pop_label_map_[seq_label] = pop_label;
    }

    // ECOEVOLITY_DEBUG(
    // std::cerr << "this->populations_labels_:" << std::endl;
    // unsigned int pop_idx = 0;
    // for (auto p_label: this->population_labels_) {
    //     std::cerr << p_label << std::endl;
    //     for (auto s_label: this->sequence_labels_[pop_idx]) {
    //         std::cerr << "\t" << s_label << std::endl;
    //     }
    //     pop_idx += 1;
    // }
    // )

    std::vector<const NxsDiscreteDatatypeMapper *> data_type_mappers = char_block->GetAllDatatypeMappers();
    if (data_type_mappers.size() < 1) {
        throw EcoevolityParsingError("No character encoding found", this->path_, 0);
    }
    if (data_type_mappers.size() > 1) {
        throw EcoevolityParsingError("More than one character encoding (i.e., mixed data types) found", this->path_, 0);
    }
    const NxsDiscreteDatatypeMapper * data_type_mapper = data_type_mappers[0];

    unsigned int ploidy_multiplier = 1;
    if (this->genotypes_are_diploid_) {
        ploidy_multiplier = 2;
    }
    unsigned int found_pattern_idx = 0;
    bool pattern_was_found = false;
    if (data_type == NxsCharactersBlock::DataTypesEnum::standard) {
        const NxsDiscreteStateCell highest_state_code = data_type_mapper->GetHighestStateCode();
        if (highest_state_code == 1) {
            if (this->genotypes_are_diploid_) {
                throw EcoevolityBiallelicDataError(
                        "Cannot limit diploid data to 0/1 characters",
                        this->path_);
            }
        }
        else if (highest_state_code == 2) {
            if (! this->genotypes_are_diploid_) {
                throw EcoevolityBiallelicDataError(
                        "Haploid data cannot have state codes greater than 1",
                        this->path_);
            }
        }
        else {
            throw EcoevolityParsingError("More than 3 character state codes found", this->path_, 0);
        }

        for (unsigned int site_idx = 0; site_idx < num_chars; ++site_idx) {
            std::vector<unsigned int> allele_cts (this->get_number_of_populations(), 0);
            std::vector<unsigned int> red_allele_cts (this->get_number_of_populations(), 0);
            for (unsigned int taxon_idx = 0; taxon_idx < num_taxa; ++taxon_idx) {
                const NxsDiscreteStateCell state_code = char_block->GetInternalRepresentation(taxon_idx, site_idx);
                if (state_code >= 0) {
                    const unsigned int num_states = char_block->GetNumStates(taxon_idx, site_idx);
                    if (num_states > 1) {
                        throw EcoevolityInvalidCharacterError(
                                "Invalid polymorphic character",
                                this->path_,
                                char_block->GetTaxonLabel(taxon_idx),
                                site_idx);
                    }
                    if ((state_code > 1) && (! this->genotypes_are_diploid_)) {
                        throw EcoevolityInvalidCharacterError(
                                "Invalid diploid character (2) for haploid data",
                                this->path_,
                                char_block->GetTaxonLabel(taxon_idx),
                                site_idx);
                    }
                    const unsigned int& population_idx = this->get_population_index_from_seq_label(char_block->GetTaxonLabel(taxon_idx));
                    red_allele_cts[population_idx] += state_code;
                    allele_cts[population_idx] += 1 * ploidy_multiplier;
                }
            }
            this->get_pattern_index(pattern_was_found, found_pattern_idx,
                    red_allele_cts,
                    allele_cts);
            if (pattern_was_found) {
                this->pattern_weights_[found_pattern_idx] += 1;
            }
            else {
                this->red_allele_counts_.push_back(red_allele_cts);
                this->allele_counts_.push_back(allele_cts);
                this->pattern_weights_.push_back(1);
            }
        }
    }
    else if ((data_type == NxsCharactersBlock::DataTypesEnum::nucleotide) ||
             (data_type == NxsCharactersBlock::DataTypesEnum::dna) ||
             (data_type == NxsCharactersBlock::DataTypesEnum::rna)) {
        if (this->markers_are_dominant_) {
            throw EcoevolityBiallelicDataError(
                    "Dominant data must be coded as 0/1 (not nucleotides)",
                    this->path_);
        }
        for (unsigned int site_idx = 0; site_idx < num_chars; ++site_idx) {
            std::vector<unsigned int> allele_cts (this->get_number_of_populations(), 0);
            std::vector<unsigned int> red_allele_cts (this->get_number_of_populations(), 0);
            NxsDiscreteStateCell red_code = -1;
            NxsDiscreteStateCell green_code = -1;
            for (unsigned int taxon_idx = 0; taxon_idx < num_taxa; ++taxon_idx) {
                const NxsDiscreteStateCell state_code = char_block->GetInternalRepresentation(taxon_idx, site_idx);
                if (state_code >= 0) {
                    const unsigned int num_states = char_block->GetNumStates(taxon_idx, site_idx);
                    if (num_states > 3) {
                        throw EcoevolityInvalidCharacterError(
                                "Invalid polymorphic character with 3 or more states",
                                this->path_,
                                char_block->GetTaxonLabel(taxon_idx),
                                site_idx);
                    }
                    if ((num_states > 1) && (! this->genotypes_are_diploid_)) {
                        throw EcoevolityInvalidCharacterError(
                                "Polymorphic characters are not allowed for haploid data",
                                this->path_,
                                char_block->GetTaxonLabel(taxon_idx),
                                site_idx);
                    }
                    ECOEVOLITY_ASSERT((num_states > 0) && (num_states < 3));
                    std::vector<NxsDiscreteStateCell> states;
                    states.push_back(char_block->GetInternalRepresentation(taxon_idx, site_idx, 0));
                    if (num_states > 1) {
                        states.push_back(char_block->GetInternalRepresentation(taxon_idx, site_idx, 1));
                    }
                    const unsigned int& population_idx = this->get_population_index_from_seq_label(char_block->GetTaxonLabel(taxon_idx));
                    unsigned int pm = ploidy_multiplier;
                    if (num_states > 1) {
                        pm = 1;
                    }
                    for (auto state_iter: states) {
                        if (green_code < 0) {
                            green_code = state_iter;
                            allele_cts[population_idx] += 1 * pm;
                            continue;
                        }
                        else if (green_code == state_iter) {
                            allele_cts[population_idx] += 1 * pm;
                            continue;
                        }
                        else if (red_code < 0) {
                            red_code = state_iter;
                            red_allele_cts[population_idx] += 1 * pm;
                            allele_cts[population_idx] += 1 * pm;
                            continue;
                        }
                        else if (red_code == state_iter) {
                            red_allele_cts[population_idx] += 1 * pm;
                            allele_cts[population_idx] += 1 * pm;
                            continue;
                        }
                        else {
                            throw EcoevolityInvalidCharacterError(
                                    "A third character state was found",
                                    this->path_,
                                    char_block->GetTaxonLabel(taxon_idx),
                                    site_idx);
                        }
                    }
                }
            }
            this->get_pattern_index(pattern_was_found, found_pattern_idx,
                    red_allele_cts,
                    allele_cts);
            if (pattern_was_found) {
                this->pattern_weights_[found_pattern_idx] += 1;
            }
            else {
                this->red_allele_counts_.push_back(red_allele_cts);
                this->allele_counts_.push_back(allele_cts);
                this->pattern_weights_.push_back(1);
            }
        }
    }
    else {
        throw EcoevolityBiallelicDataError("Data type not supported", this->path_);
    }
    nexus_reader.DeleteBlocksFromFactories();
    this->update_max_allele_counts();
    this->update_pattern_booleans();
    if (validate) {
        this->validate();
    }
}

void BiallelicData::update_max_allele_counts() {
    this->max_allele_counts_.assign(this->get_number_of_populations(), 0);
    for (unsigned int pattern_idx = 0; pattern_idx < this->get_number_of_patterns(); ++pattern_idx) {
        const std::vector<unsigned int>& allele_cts = this->get_allele_counts(pattern_idx);
        for (unsigned int pop_idx = 0; pop_idx < allele_cts.size(); ++pop_idx) {
            if (allele_cts.at(pop_idx) > this->max_allele_counts_.at(pop_idx)) {
                this->max_allele_counts_.at(pop_idx) = allele_cts.at(pop_idx);
            }
        }
    }
}

const std::vector<unsigned int>& BiallelicData::get_red_allele_counts(unsigned int pattern_index) const {
    return this->red_allele_counts_.at(pattern_index);
}

const std::vector<unsigned int>& BiallelicData::get_allele_counts(unsigned int pattern_index) const {
    return this->allele_counts_.at(pattern_index);
}
unsigned int BiallelicData::get_red_allele_count(
        unsigned int pattern_index,
        unsigned int population_index) const {
    return this->red_allele_counts_.at(pattern_index).at(population_index);
}
unsigned int BiallelicData::get_allele_count(
        unsigned int pattern_index,
        unsigned int population_index) const {
    return this->allele_counts_.at(pattern_index).at(population_index);
}

unsigned int BiallelicData::get_pattern_weight(unsigned int pattern_index) const {
    return this->pattern_weights_.at(pattern_index);
}

unsigned int BiallelicData::get_max_allele_count(unsigned int population_index) const {
    return this->max_allele_counts_.at(population_index);
}

const std::vector<unsigned int>& BiallelicData::get_max_allele_counts() const {
    return this->max_allele_counts_;
}

unsigned int BiallelicData::get_population_index(std::string population_label) const {
    return this->pop_label_to_index_map_.at(population_label);
}

unsigned int BiallelicData::get_population_index_from_seq_label(std::string seq_label) const {
    const std::string pop_label = this->seq_label_to_pop_label_map_.at(seq_label);
    return this->get_population_index(pop_label);
}

const std::string& BiallelicData::get_population_label(unsigned int population_index) const {
    return this->population_labels_.at(population_index);
}

const std::vector<std::string>& BiallelicData::get_sequence_labels(unsigned int population_index) const {
    return this->sequence_labels_.at(population_index);
}

const std::string& BiallelicData::get_path() const {
    return this->path_;
}

unsigned int BiallelicData::get_number_of_patterns() const {
    return this->pattern_weights_.size();
}

unsigned int BiallelicData::get_number_of_sites() const {
    unsigned int nsites = 0;
    for (auto weight_iter: this->pattern_weights_) {
        nsites += weight_iter;
    }
    return nsites;
}

unsigned int BiallelicData::get_number_of_populations() const {
    return this->population_labels_.size();
}

bool BiallelicData::markers_are_dominant() const {
    return this->markers_are_dominant_;
}

bool BiallelicData::genotypes_are_diploid() const {
    return this->genotypes_are_diploid_;
}

bool BiallelicData::has_constant_patterns() const {
    return this->has_constant_patterns_;
}

bool BiallelicData::has_missing_population_patterns() const {
    return this->has_missing_population_patterns_;
}

bool BiallelicData::has_mirrored_patterns() const {
    return this->has_mirrored_patterns_;
}

bool BiallelicData::patterns_are_folded() const {
    return this->patterns_are_folded_;
}

void BiallelicData::update_has_constant_patterns() {
    this->has_constant_patterns_ = false;
    std::vector<unsigned int> no_red_pattern (this->get_number_of_populations(), 0);
    for (unsigned int pattern_idx = 0; pattern_idx < this->get_number_of_patterns(); ++pattern_idx) {
        if ((this->get_red_allele_counts(pattern_idx) == no_red_pattern) ||
            (this->get_red_allele_counts(pattern_idx) == this->get_allele_counts(pattern_idx))) {
            this->has_constant_patterns_ = true;
            return;
        }
    }
    return;
}

void BiallelicData::update_has_missing_population_patterns() {
    this->has_missing_population_patterns_ = false;
    for (unsigned int pattern_idx = 0; pattern_idx < this->get_number_of_patterns(); ++pattern_idx) {
        for (auto count_iter: this->get_allele_counts(pattern_idx)) {
            if (count_iter == 0) {
                this->has_missing_population_patterns_ = true;
                return;
            }
        }
    }
    return;
}

void BiallelicData::update_has_mirrored_patterns() {
    this->has_mirrored_patterns_ = false;
    unsigned int mirrored_idx = 0;
    bool was_found = false;
    for (unsigned int pattern_idx = 0; pattern_idx < this->get_number_of_patterns(); ++pattern_idx) {
        const std::vector< std::vector<unsigned int> > mirrored_pattern = this->get_mirrored_pattern(pattern_idx);
        this->get_pattern_index(was_found, mirrored_idx,
                mirrored_pattern.at(0),
                mirrored_pattern.at(1));
        if ((was_found) && (mirrored_idx != pattern_idx)) {
            ECOEVOLITY_ASSERT(mirrored_idx > pattern_idx);
            this->has_mirrored_patterns_ = true;
            return;
        }
    }
    return;
}

void BiallelicData::update_patterns_are_folded() {
    this->patterns_are_folded_ = true;
    this->update_has_mirrored_patterns();
    if (this->has_mirrored_patterns()) {
        this->patterns_are_folded_ = false;
        return;
    }
    for (unsigned int pattern_idx = 0; pattern_idx < this->get_number_of_patterns(); ++pattern_idx) {
        const std::vector< std::vector<unsigned int> > mirrored_pattern = this->get_mirrored_pattern(pattern_idx);
        const std::vector<unsigned int>& green_cts = mirrored_pattern.at(0);
        const std::vector<unsigned int>& red_cts = this->get_red_allele_counts(pattern_idx);
        unsigned int red_total = 0;
        unsigned int green_total = 0;
        ECOEVOLITY_ASSERT(green_cts.size() == red_cts.size());
        for (unsigned int pop_idx = 0; pop_idx < red_cts.size(); ++pop_idx) {
            red_total += red_cts.at(pop_idx);
            green_total += green_cts.at(pop_idx);
        }
        if (green_total < red_total) {
            this->patterns_are_folded_ = false;
            return;
        }
    }
    return;
}

void BiallelicData::update_pattern_booleans() {
    this->update_has_constant_patterns();
    this->update_has_missing_population_patterns();
    this->update_has_mirrored_patterns();
    this->update_patterns_are_folded();
}

void BiallelicData::get_pattern_index(
        bool& was_found,
        unsigned int& pattern_index,
        std::vector<unsigned int> red_allele_counts,
        std::vector<unsigned int> allele_counts) const {
    ECOEVOLITY_ASSERT(this->allele_counts_.size() == this->red_allele_counts_.size());
    was_found = false;
    pattern_index = 0;
    for (unsigned int pattern_idx = 0; pattern_idx < this->allele_counts_.size(); ++pattern_idx) {
        if ((this->red_allele_counts_[pattern_idx] == red_allele_counts) &&
            (this->allele_counts_[pattern_idx] == allele_counts)) {
            pattern_index = pattern_idx;
            was_found = true;
            return;
        }
    }
}

// TODO: This method is causing the following error from compiler when using
// flag '-Wstrict-overflow=5':
//   error: assuming signed overflow does not occur when changing X +- C1 cmp C2 to X cmp C2 -+ C1
// This seems to be triggered by the 'erase' method calls and I have no idea
// why. Hacky fix for now is to use -Wstrict-overflow=2.
void BiallelicData::remove_pattern(unsigned int pattern_index) {
    ECOEVOLITY_ASSERT(pattern_index < this->pattern_weights_.size());
    this->pattern_weights_.erase(this->pattern_weights_.begin() + pattern_index);
    this->allele_counts_.erase(this->allele_counts_.begin() + pattern_index);
    this->red_allele_counts_.erase(this->red_allele_counts_.begin() + pattern_index);
    if (this->pattern_weights_.size() < 1) {
        throw EcoevolityBiallelicDataError(
                "Ran out of data while removing patterns",
                this->path_);
    }
    return;
}

void BiallelicData::remove_first_constant_pattern(
        bool& was_removed,
        unsigned int& removed_index) {
    std::vector<unsigned int> no_red_pattern (this->get_number_of_populations(), 0);
    removed_index = 0;
    was_removed = false;
    for (unsigned int pattern_idx = 0; pattern_idx < this->get_number_of_patterns(); ++pattern_idx) {
        if (this->get_red_allele_counts(pattern_idx) == no_red_pattern) {
            this->number_of_constant_green_sites_removed_ += this->get_pattern_weight(pattern_idx);
            this->remove_pattern(pattern_idx);
            removed_index = pattern_idx;
            was_removed = true;
            return;
        }
        if (this->get_red_allele_counts(pattern_idx) == this->get_allele_counts(pattern_idx)) {
            this->number_of_constant_red_sites_removed_ += this->get_pattern_weight(pattern_idx);
            this->remove_pattern(pattern_idx);
            removed_index = pattern_idx;
            was_removed = true;
            return;
        }
    }
}

void BiallelicData::remove_first_missing_population_pattern(
        bool& was_removed,
        unsigned int& removed_index) {
    removed_index = 0;
    was_removed = false;
    for (unsigned int pattern_idx = 0; pattern_idx < this->get_number_of_patterns(); ++pattern_idx) {
        for (unsigned int pop_idx = 0; pop_idx < this->get_number_of_populations(); ++pop_idx) {
            if (this->get_allele_count(pattern_idx, pop_idx) == 0) {
                this->number_of_missing_sites_removed_ += this->get_pattern_weight(pattern_idx);
                this->remove_pattern(pattern_idx);
                removed_index = pattern_idx;
                was_removed = true;
                return;
            }
        }
    }
}

const std::vector< std::vector<unsigned int> > BiallelicData::get_mirrored_pattern(unsigned int pattern_index) const {
    std::vector< std::vector<unsigned int> > mirrored_counts;
    const std::vector<unsigned int>& allele_cts = this->get_allele_counts(pattern_index);
    const std::vector<unsigned int>& red_cts = this->get_red_allele_counts(pattern_index);
    std::vector<unsigned int> green_cts(red_cts.size());
    for (unsigned int pop_idx = 0; pop_idx < red_cts.size(); ++pop_idx) {
        green_cts.at(pop_idx) = allele_cts.at(pop_idx) - red_cts.at(pop_idx);
    }
    mirrored_counts.push_back(green_cts);
    mirrored_counts.push_back(allele_cts);
    return mirrored_counts;
}

void BiallelicData::fold_first_mirrored_pattern(
        bool& was_folded,
        unsigned int& folded_index) {
    was_folded = false;
    folded_index = 0;
    unsigned int mirrored_idx = 0;
    bool was_found = false;
    for (unsigned int pattern_idx = 0; pattern_idx < this->get_number_of_patterns(); ++pattern_idx) {
        const std::vector< std::vector<unsigned int> > mirrored_pattern = this->get_mirrored_pattern(pattern_idx);
        this->get_pattern_index(was_found, mirrored_idx, mirrored_pattern.at(0),
                mirrored_pattern.at(1));
        if (! was_found) {
            continue;
        }
        if (mirrored_idx == pattern_idx) {
            continue;
        }
        ECOEVOLITY_ASSERT(mirrored_idx > pattern_idx);
        this->pattern_weights_.at(pattern_idx) += this->pattern_weights_.at(mirrored_idx);
        this->remove_pattern(mirrored_idx);
        folded_index = mirrored_idx;
        was_folded = true;
        return;
    }
}

unsigned int BiallelicData::remove_constant_patterns(const bool validate) {
    unsigned int number_removed = 0;
    unsigned int return_idx = 0;
    bool was_found = false;
    while (true) {
        this->remove_first_constant_pattern(was_found, return_idx);
        if (! was_found) {
            break;
        }
        number_removed += 1;
    }
    this->update_pattern_booleans();
    this->update_max_allele_counts();
    if (validate) {
        this->validate();
    }
    return number_removed;
}

unsigned int BiallelicData::remove_missing_population_patterns(const bool validate) {
    unsigned int number_removed = 0;
    unsigned int return_idx = 0;
    bool was_found = false;
    while (true) {
        this->remove_first_missing_population_pattern(was_found, return_idx);
        if (! was_found) {
            break;
        }
        number_removed += 1;
    }
    this->update_pattern_booleans();
    this->update_max_allele_counts();
    if (validate) {
        this->validate();
    }
    return number_removed;
}

unsigned int BiallelicData::fold_patterns(const bool validate) {
    if (this->markers_are_dominant()) {
        throw EcoevolityBiallelicDataError(
                "Site patterns cannot be folded for dominant markers",
                this->path_);
    }
    unsigned int number_removed = 0;
    unsigned int return_idx = 0;
    bool was_found = false;
    while (true) {
        this->fold_first_mirrored_pattern(was_found, return_idx);
        if (! was_found) {
            break;
        }
        number_removed += 1;
    }
    for (unsigned int pattern_idx = 0; pattern_idx < this->get_number_of_patterns(); ++pattern_idx) {
        const std::vector< std::vector<unsigned int> > mirrored_pattern = this->get_mirrored_pattern(pattern_idx);
        const std::vector<unsigned int>& green_cts = mirrored_pattern.at(0);
        const std::vector<unsigned int>& red_cts = this->get_red_allele_counts(pattern_idx);
        unsigned int red_total = 0;
        unsigned int green_total = 0;
        ECOEVOLITY_ASSERT(green_cts.size() == red_cts.size());
        for (unsigned int pop_idx = 0; pop_idx < red_cts.size(); ++pop_idx) {
            red_total += red_cts.at(pop_idx);
            green_total += green_cts.at(pop_idx);
        }
        if (green_total < red_total) {
            this->red_allele_counts_.at(pattern_idx) = green_cts;
        }
    }
    this->update_pattern_booleans();
    if (validate) {
        this->validate();
    }
    return number_removed;
}

unsigned int BiallelicData::get_number_of_constant_sites_removed() const {
    return (this->number_of_constant_green_sites_removed_ +
            this->number_of_constant_red_sites_removed_);
}

unsigned int BiallelicData::get_number_of_constant_green_sites_removed() const {
    return this->number_of_constant_green_sites_removed_;
}

unsigned int BiallelicData::get_number_of_constant_red_sites_removed() const {
    return this->number_of_constant_red_sites_removed_;
}

unsigned int BiallelicData::get_number_of_missing_sites_removed() const {
    return this->number_of_missing_sites_removed_;
}

double BiallelicData::get_proportion_of_red_alleles() const {
    unsigned int red_count = 0;
    unsigned int total_count = 0;
    for (unsigned int pattern_idx = 0;
            pattern_idx < this->get_number_of_patterns();
            ++pattern_idx) {
        for (unsigned int pop_idx = 0;
                pop_idx < this->get_number_of_populations();
                ++pop_idx) {
            red_count += this->get_red_allele_count(pattern_idx, pop_idx) * this->get_pattern_weight(pattern_idx);
            total_count += this->get_allele_count(pattern_idx, pop_idx) * this->get_pattern_weight(pattern_idx);
        }
    }
    return (double) red_count / (double) total_count;
}

void BiallelicData::get_empirical_u_v_rates(double& u, double& v) const {
    double p_red = this->get_proportion_of_red_alleles();
    u = 1.0 / (2.0 * p_red);
    v = 1.0 / (2.0 * (1.0 - p_red));
}

void BiallelicData::validate() const {
    if (this->allele_counts_.size() != this->pattern_weights_.size()) {
        throw EcoevolityBiallelicDataError(
                "Different number of allele counts and weights",
                this->path_);
    }
    if (this->allele_counts_.size() != this->red_allele_counts_.size()) {
        throw EcoevolityBiallelicDataError(
                "Different number of allele and red allele counts",
                this->path_);
    }
    if (this->allele_counts_.size() < 1) {
        throw EcoevolityBiallelicDataError(
                "No data found", this->path_);
    }
    std::vector<unsigned int> no_red_pattern (this->get_number_of_populations(), 0);
    unsigned int number_of_pops = this->population_labels_.size();
    bool has_constant = false;
    bool has_missing = false;
    for (unsigned int pattern_idx = 0; pattern_idx < this->get_number_of_patterns(); ++pattern_idx) {
        if (this->allele_counts_.at(pattern_idx).size() != number_of_pops) {
            throw EcoevolityBiallelicDataError(
                    "Different number of populations and allele counts per pattern",
                    this->path_);
        }
        if (this->red_allele_counts_.at(pattern_idx).size() != number_of_pops) {
            throw EcoevolityBiallelicDataError(
                    "Different number of populations and red allele counts per pattern",
                    this->path_);
        }
        if ((this->red_allele_counts_.at(pattern_idx) == no_red_pattern) ||
            (this->red_allele_counts_.at(pattern_idx) == this->allele_counts_.at(pattern_idx))) {
            has_constant = true;
        }
        for (unsigned int pop_idx = 0; pop_idx < number_of_pops; ++pop_idx) {
            if (this->red_allele_counts_.at(pattern_idx).at(pop_idx) > this->allele_counts_.at(pattern_idx).at(pop_idx)) {
                throw EcoevolityBiallelicDataError(
                    "Some patterns have more red alleles than total allele counts",
                    this->path_);
            }
            if (this->allele_counts_.at(pattern_idx).at(pop_idx) == 0) {
                has_missing = true;
            }
        }
    }
    if (has_constant != this->has_constant_patterns()) {
        throw EcoevolityBiallelicDataError(
                "constant pattern boolean is incorrect",
                this->path_);
    }
    if (has_missing != this->has_missing_population_patterns()) {
        throw EcoevolityBiallelicDataError(
                "missing data pattern boolean is incorrect",
                this->path_);
    }
    return;
}

void BiallelicData::write_summary(
        std::ostream& out,
        unsigned int indent_level) const {
    std::string margin = string_util::get_indent(indent_level);
    std::string indent = string_util::get_indent(1);
    out << std::boolalpha 
        << margin <<"Summary of data from \'" << this->get_path() << "\':\n"; 
    if (this->genotypes_are_diploid()) {
        out << margin << indent << "Genotypes: diploid\n";
    }
    else {
        out << margin << indent << "Genotypes: haploid\n";
    }
    out << margin << indent
            << "Markers are dominant? "
            << this->markers_are_dominant() << "\n"
        << margin << indent
            << "Number of populations: "
            << this->get_number_of_populations() << "\n"
        << margin << indent
            << "Number of sites: "
            << this->get_number_of_sites() << "\n"
        << margin << indent
            << "Number of patterns: "
            << this->get_number_of_patterns() << "\n"
        << margin << indent
            << "Patterns folded? "
            << this->patterns_are_folded() << "\n"
        << margin << indent
            << "Population label (max # of alleles sampled):\n";
    for (unsigned int i = 0; i < this->get_number_of_populations(); ++i) {
        out << margin << indent << indent
            << this->get_population_label(i)
            << " (" << this->get_max_allele_count(i) << ")\n";
    }
}
