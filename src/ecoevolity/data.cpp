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
        bool validate,
        bool store_seq_loci_info) {
    this->init(path,
               population_name_delimiter,
               population_name_is_prefix,
               genotypes_are_diploid,
               markers_are_dominant,
               validate,
               store_seq_loci_info);
}

void BiallelicData::init_from_yaml_stream(
        std::istream& stream,
        const std::string& path,
        bool validate) {
    this->path_ = path;
    this->parse_yaml_data(stream, validate);
}

void BiallelicData::init_from_yaml_file(
        const std::string& path,
        bool validate) {
    std::ifstream in_stream;
    in_stream.open(path);
    if (! in_stream.is_open()) {
        throw EcoevolityYamlDataError(
                "Could not open YAML data file",
                path);
    }
    this->init_from_yaml_stream(in_stream, path, validate);
    in_stream.close();
}

void BiallelicData::parse_yaml_data(std::istream& yaml_stream, bool validate) {
    YAML::Node data;
    try {
        data = YAML::Load(yaml_stream);
    }
    catch (...) {
        std::cerr << "ERROR: Problem with YAML-formatting of data\n";
        throw;
    }
    this->parse_yaml_top_level(data);
    this->update_max_allele_counts();
    this->update_pattern_booleans();
    if (validate) {
        this->validate();
    }
}

void BiallelicData::parse_yaml_top_level(const YAML::Node& top_level_node) {
    if (! top_level_node.IsMap()) {
        throw EcoevolityYamlDataError(
                "Expecting top-level of config to be a map, but found: " +
                YamlCppUtils::get_node_type(top_level_node),
                this->path_);
    }
    if (! top_level_node["population_labels"]) {
        throw EcoevolityYamlDataError("No population labels", this->path_);
    }
    if (! top_level_node["allele_count_patterns"]) {
        throw EcoevolityYamlDataError("No allele count patterns", this->path_);
    }
    if (! top_level_node["pattern_weights"]) {
        throw EcoevolityYamlDataError("No pattern weights", this->path_);
    }
    if (! top_level_node["markers_are_dominant"]) {
        this->markers_are_dominant_ = false;
    }
    else {
        this->parse_yaml_marker_dominance(top_level_node["markers_are_dominant"]);
    }
    this->parse_yaml_population_labels(top_level_node["population_labels"]);
    this->parse_yaml_pattern_weights(top_level_node["pattern_weights"]);
    this->parse_yaml_all_allele_count_patterns(top_level_node["allele_count_patterns"]);
}

void BiallelicData::parse_yaml_marker_dominance(const YAML::Node& node) {
    if (! node.IsScalar()) {
        throw EcoevolityYamlDataError(
                "markers_are_dominance node should be a map, but found: " +
                YamlCppUtils::get_node_type(node));
    }
    this->markers_are_dominant_ = node.as<bool>();
}

void BiallelicData::parse_yaml_population_labels(const YAML::Node& node) {
    if (! node.IsSequence()) {
        throw EcoevolityYamlDataError(
                "Expecting population labels to be a sequence, but found: " +
                YamlCppUtils::get_node_type(node),
                this->path_);
    }
    this->population_labels_.clear();
    this->pop_label_to_index_map_.clear();
    unsigned int pop_idx = 0;
    std::set<std::string> label_set;
    for (YAML::const_iterator label = node.begin();
            label != node.end();
            ++label) {
        std::string l = label->as<std::string>();
        if (label_set.count(l) > 0) {
            throw EcoevolityYamlDataError(
                    "Duplicate population label: " + l,
                    this->path_);
        }
        label_set.insert(l);
        this->population_labels_.push_back(l);
        this->pop_label_to_index_map_[l] = pop_idx;
        ++pop_idx;
    }
}

void BiallelicData::parse_yaml_pattern_weights(const YAML::Node& node) {
    if (! node.IsSequence()) {
        throw EcoevolityYamlDataError(
                "Expecting pattern weights to be a sequence, but found: " +
                YamlCppUtils::get_node_type(node),
                this->path_);
    }
    this->pattern_weights_.clear();
    for (YAML::const_iterator w = node.begin();
            w != node.end();
            ++w) {
        unsigned int weight = w->as<unsigned int>();
        this->pattern_weights_.push_back(weight);
    }
}

void BiallelicData::parse_yaml_all_allele_count_patterns(const YAML::Node& node) {
    if (! node.IsSequence()) {
        throw EcoevolityYamlDataError(
                "Expecting allele count patterns to be a sequence, but found: " +
                YamlCppUtils::get_node_type(node),
                this->path_);
    }
    this->allele_counts_.clear();
    this->red_allele_counts_.clear();
    std::vector<unsigned int> tmp_red_allele_counts;
    std::vector<unsigned int> tmp_allele_counts;
    std::vector< std::vector<unsigned int> > pattern;
    std::set< std::vector< std::vector<unsigned int> > > pattern_set;
    for (unsigned int pattern_idx = 0;
            pattern_idx < node.size();
            ++pattern_idx) {
        YAML::Node p = node[pattern_idx];
        this->parse_yaml_allele_count_pattern(p,
                tmp_red_allele_counts,
                tmp_allele_counts);
        ECOEVOLITY_ASSERT(tmp_red_allele_counts.size() == tmp_allele_counts.size());
        if (tmp_allele_counts.size() != this->population_labels_.size()) {
            std::ostringstream message;
            message << "There were "
                    << this->population_labels_.size()
                    << " population labels, so expecting all allele count patterns to consist of counts for "
                    << this->population_labels_.size()
                    << " populations, but found a pattern with counts for "
                    << tmp_allele_counts.size()
                    << " populations\n";
            throw EcoevolityYamlDataError(message.str());
        }
        pattern = {tmp_red_allele_counts, tmp_allele_counts};
        if (pattern_set.count(pattern) > 0) {
            throw EcoevolityYamlDataError("Found duplicate allele count pattern",
                    this->path_);
        }
        pattern_set.insert(pattern);
        this->allele_counts_.push_back(tmp_allele_counts);
        this->red_allele_counts_.push_back(tmp_red_allele_counts);
    }
}

void BiallelicData::parse_yaml_allele_count_pattern(const YAML::Node& node,
        std::vector<unsigned int> & red_allele_counts,
        std::vector<unsigned int> & allele_counts) {
    if (! node.IsSequence()) {
        throw EcoevolityYamlDataError(
                "Expecting each allele count pattern to be a sequence, but found: " +
                YamlCppUtils::get_node_type(node),
                this->path_);
    }
    red_allele_counts.clear();
    allele_counts.clear();
    for (unsigned int pop_idx = 0;
            pop_idx < node.size();
            ++pop_idx) {
        YAML::Node pop_counts = node[pop_idx];
        this->parse_yaml_allele_count(pop_counts,
                red_allele_counts,
                allele_counts);
    }
}

void BiallelicData::parse_yaml_allele_count(const YAML::Node& node,
        std::vector<unsigned int> & red_allele_counts,
        std::vector<unsigned int> & allele_counts) {
    if (! node.IsSequence()) {
        throw EcoevolityYamlDataError(
                "Expecting each population allele count to be a sequence, but found: " +
                YamlCppUtils::get_node_type(node),
                this->path_);
    }
    std::vector<unsigned int> temp_counts;
    for (YAML::const_iterator count = node.begin();
            count != node.end();
            ++count) {
        unsigned int c = count->as<unsigned int>();
        temp_counts.push_back(c);
    }
    if (temp_counts.size() != 2) {
        throw EcoevolityYamlDataError(
                "All population allele counts should be a sequence of 2 integers, but found one with length: " +
                temp_counts.size(),
                this->path_);
    }
    red_allele_counts.push_back(temp_counts.at(0));
    allele_counts.push_back(temp_counts.at(1));
}

void BiallelicData::init(
        std::string path, 
        char population_name_delimiter,
        bool population_name_is_prefix,
        bool genotypes_are_diploid,
        bool markers_are_dominant,
        bool validate,
        bool store_seq_loci_info) {
    char pop_name_delimiter = population_name_delimiter;
    // if (population_name_delimiter == '_') {
    //     pop_name_delimiter = ' ';
    // }
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
    //

    // Store charset info if available
    if (store_seq_loci_info) {
        unsigned int num_assumptions_blocks = nexus_reader.GetNumAssumptionsBlocks(char_block);
        if (num_assumptions_blocks > 1) {
            throw EcoevolityParsingError("More than one assumptions block found", this->path_, 0);
        }
        // If no sets block, do nothing (no error)
        if (num_assumptions_blocks > 0) {
            NxsAssumptionsBlock * assumptions_block = nexus_reader.GetAssumptionsBlock(char_block, 0);
            unsigned int num_charsets = assumptions_block->GetNumCharSets();
            // If there's a sets block but no charsets, let's throw an error
            if (num_charsets < 1) {
                throw EcoevolityParsingError("No charsets found", this->path_, 0);
            }
            std::vector<unsigned int> site_indices_found (num_chars, 0);
            // std::vector<NxsString> charset_names (num_charsets);
            NxsStringVector charset_names (num_charsets);
            assumptions_block->GetCharSetNames(charset_names);
            for (unsigned int cs_idx = 0; cs_idx < num_charsets; ++cs_idx) {
                const NxsUnsignedSet * charset = assumptions_block->GetCharSet(charset_names.at(cs_idx));
                auto mnmx = std::minmax_element(charset->begin(), charset->end());
                // Vet charset to make sure it's a contiguous set of sites
                for (unsigned int site_idx_to_check = *mnmx.first;
                        site_idx_to_check < *mnmx.second + 1;
                        ++site_idx_to_check) {
                    if (charset->count(site_idx_to_check) < 1) {
                        std::ostringstream message;
                        message << "Did not find site "
                                << site_idx_to_check + 1
                                << " in charset \'"
                                << charset_names.at(cs_idx)
                                << "\'\n";
                        throw EcoevolityParsingError(message.str(), this->path_, 0);
                    }
                    if (charset->count(site_idx_to_check) > 1) {
                        std::ostringstream message;
                        message << "Found site "
                                << site_idx_to_check + 1
                                << " multiple times in charset \'"
                                << charset_names.at(cs_idx)
                                << "\'\n";
                        throw EcoevolityParsingError(message.str(), this->path_, 0);
                    }
                    if (site_indices_found.at(site_idx_to_check) > 0) {
                        std::ostringstream message;
                        message << "Site "
                                << site_idx_to_check + 1
                                << " contained in charset \'"
                                << charset_names.at(cs_idx)
                                << "\' was found previously\n";
                        throw EcoevolityParsingError(message.str(), this->path_, 0);
                    }
                    ++site_indices_found.at(site_idx_to_check);
                }
                this->locus_end_indices_.push_back(*mnmx.second);
            }
            std::sort(this->locus_end_indices_.begin(), this->locus_end_indices_.end());
            if (this->locus_end_indices_.back() != (num_chars - 1)) {
                std::ostringstream message;
                message << "Charset error; last locus position is "
                        << this->locus_end_indices_.back() + 1
                        << " but there are "
                        << num_chars
                        << " sites in the alignment\n";
                throw EcoevolityParsingError(message.str(), this->path_, 0);
            }
            // Vet the sets for missing sites
            unsigned int total_sites = 0;
            std::vector<unsigned int> missing_sites;
            for (unsigned int site_idx = 0; site_idx < num_chars; ++ site_idx) {
                if (site_indices_found.at(site_idx) < 1) {
                    missing_sites.push_back(site_idx);
                }
                total_sites += site_indices_found.at(site_idx);
            }
            if (missing_sites.size() > 0) {
                std::ostringstream message;
                message << "The following sites were missing from the charsets:\n"
                        << missing_sites.at(0);
                for (unsigned int i = 1; i < missing_sites.size(); ++i) {
                    message << ", " << missing_sites.at(i);
                }
                message << "\n";
                throw EcoevolityParsingError(message.str(), this->path_, 0);
            }
            ECOEVOLITY_ASSERT(total_sites == num_chars);
            this->storing_seq_loci_info_ = true;
        }
    }


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
                if (this->storing_seq_loci_info_) {
                    this->contiguous_pattern_indices_.push_back(found_pattern_idx);
                }
            }
            else {
                this->red_allele_counts_.push_back(red_allele_cts);
                this->allele_counts_.push_back(allele_cts);
                this->pattern_weights_.push_back(1);
                if (this->storing_seq_loci_info_) {
                    this->contiguous_pattern_indices_.push_back(this->pattern_weights_.size() - 1);
                }
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
            bool triallelic_site = false;
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
                        // Handle 3rd or 4th alleles
                        else {
                            // Code 3rd or 4th alleles as red (1)
                            red_allele_cts[population_idx] += 1 * pm;
                            allele_cts[population_idx] += 1 * pm;
                            triallelic_site = true;
                            /* throw EcoevolityInvalidCharacterError( */
                            /*         "A third character state was found", */
                            /*         this->path_, */
                            /*         char_block->GetTaxonLabel(taxon_idx), */
                            /*         site_idx); */
                        }
                    }
                }
            }
            if (triallelic_site) {
                ++this->number_of_triallelic_sites_recoded_;
            }
            this->get_pattern_index(pattern_was_found, found_pattern_idx,
                    red_allele_cts,
                    allele_cts);
            if (pattern_was_found) {
                ++this->pattern_weights_[found_pattern_idx];
                if (this->storing_seq_loci_info_) {
                    this->contiguous_pattern_indices_.push_back(found_pattern_idx);
                }
            }
            else {
                this->red_allele_counts_.push_back(red_allele_cts);
                this->allele_counts_.push_back(allele_cts);
                this->pattern_weights_.push_back(1);
                if (this->storing_seq_loci_info_) {
                    this->contiguous_pattern_indices_.push_back(this->pattern_weights_.size() - 1);
                }
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

BiallelicData BiallelicData::get_empty_copy() const {
    BiallelicData copy;
    copy.markers_are_dominant_ = this->markers_are_dominant_;
    copy.genotypes_are_diploid_ = this->genotypes_are_diploid_;
    copy.patterns_are_folded_ = this->patterns_are_folded_;
    copy.population_labels_ = this->population_labels_;
    for (auto const & seq_labels: this->sequence_labels_) {
        std::vector<std::string> copied_labels = seq_labels;
        copy.sequence_labels_.push_back(copied_labels);
    }
    copy.seq_label_to_pop_label_map_ = this->seq_label_to_pop_label_map_;
    copy.pop_label_to_index_map_ = this->pop_label_to_index_map_;
    copy.appendable_ = true;
    copy.storing_seq_loci_info_ = this->storing_seq_loci_info_;
    // Do not want pattern indices or locus info to copy over
    // copy.contiguous_pattern_indices_ = this->contiguous_pattern_indices_;
    // copy.locus_end_indices_ = this->locus_end_indices_;
    return copy;
}

bool BiallelicData::add_site(
        const std::vector<unsigned int>& red_allele_counts,
        const std::vector<unsigned int>& allele_counts,
        bool filtering_constant_patterns,
        bool end_of_locus) {
    if (! this->appendable_) {
        throw EcoevolityBiallelicDataError(
                "Patterns cannot be added to a parsed dataset",
                this->path_);
    }
    for (auto cnt: allele_counts) {
        if (cnt == 0) {
            throw EcoevolityBiallelicDataError(
                    "Cannot add missing data patterns",
                    this->path_);
        }
    }
    if (filtering_constant_patterns && this->pattern_is_constant(red_allele_counts, allele_counts)) {
        if (red_allele_counts == allele_counts) {
            ++this->number_of_constant_red_sites_removed_;
        }
        else {
            ++this->number_of_constant_green_sites_removed_;
        }
        return false;
    }
    bool pattern_exists = false;
    unsigned int pattern_index = 0;
    this->get_pattern_index(
            pattern_exists,
            pattern_index,
            red_allele_counts,
            allele_counts);
    if (pattern_exists) {
        ++this->pattern_weights_.at(pattern_index);
        if (this->storing_seq_loci_info_) {
            this->contiguous_pattern_indices_.push_back(pattern_index);
            if (end_of_locus) {
                this->locus_end_indices_.push_back(this->contiguous_pattern_indices_.size() - 1);
            }
        }
        return true;
    }
    this->red_allele_counts_.push_back(red_allele_counts);
    this->allele_counts_.push_back(allele_counts);
    this->pattern_weights_.push_back(1);
    if (this->storing_seq_loci_info_) {
        this->contiguous_pattern_indices_.push_back(this->pattern_weights_.size() - 1);
        if (end_of_locus) {
            this->locus_end_indices_.push_back(this->contiguous_pattern_indices_.size() - 1);
        }
    }
    return true;
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

unsigned int BiallelicData::get_number_of_variable_sites() const {
    unsigned int nsites = 0;
    for (unsigned int i = 0; i < this->get_number_of_patterns(); ++i) {
        if (! this->pattern_is_constant(
                this->get_red_allele_counts(i),
                this->get_allele_counts(i))) {
            nsites += this->get_pattern_weight(i);
        }
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

bool BiallelicData::has_recoded_triallelic_sites() const {
    return (this->number_of_triallelic_sites_recoded_ > 0);
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

bool BiallelicData::pattern_is_constant(
        const std::vector<unsigned int>& red_allele_counts,
        const std::vector<unsigned int>& allele_counts) {
    ECOEVOLITY_ASSERT(allele_counts.size() == red_allele_counts.size());
    bool all_green = true;
    bool all_red = true;
    for (unsigned int i = 0; i < allele_counts.size(); ++i) {
        if (red_allele_counts.at(i) > 0) {
            all_green = false;
        }
        if (red_allele_counts.at(i) < allele_counts.at(i)) {
            all_red = false;
        }
        if ((! all_green) && (! all_red)) {
            break;
        }
    }
    return (all_red || all_green);
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
        const std::vector<unsigned int>& red_allele_counts,
        const std::vector<unsigned int>& allele_counts) const {
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
    // std::cout << "\n";
    // std::cout << "Deleting pattern index " << pattern_index << "\n";
    if (this->storing_seq_loci_info_) {
        std::vector<unsigned int> site_indices_to_erase;
        for (unsigned int site_idx = 0;
                site_idx < this->contiguous_pattern_indices_.size();
                ++site_idx) {
            if (this->contiguous_pattern_indices_.at(site_idx) == pattern_index) {
                site_indices_to_erase.push_back(site_idx);
            } else if (this->contiguous_pattern_indices_.at(site_idx) > pattern_index) {
                --this->contiguous_pattern_indices_.at(site_idx);
            }
        }
        // std::cout << "Sites to delete:\n";
        // for (unsigned int i = 0; i < site_indices_to_erase.size(); ++i) {
        //     std::cout << site_indices_to_erase.at(i) << " ";
        // }
        // std::cout << "\n";
        // std::cout << "Locus ends:\n";
        // for (unsigned int i = 0; i < this->locus_end_indices_.size(); ++i) {
        //     std::cout << this->locus_end_indices_.at(i) << " ";
        // }
        // std::cout << "\n";
        std::vector<unsigned int> ends_to_erase;
        for (int i = (site_indices_to_erase.size() - 1); i >= 0; --i) {
            for (unsigned int locus_idx = 0; locus_idx < this->locus_end_indices_.size(); ++locus_idx) {
                if (this->locus_end_indices_.at(locus_idx) >= site_indices_to_erase.at(i)) {
                    if ((this->locus_end_indices_.at(locus_idx) < 1) ||
                            ((locus_idx > 0) &&
                            (this->locus_end_indices_.at(locus_idx) <=
                            (this->locus_end_indices_.at(locus_idx - 1) + 1)))) {
                        ends_to_erase.push_back(locus_idx);
                    }
                    // } else {
                    //     --this->locus_end_indices_.at(locus_idx);
                    // }
                    if (this->locus_end_indices_.at(locus_idx) > 0) {
                        --this->locus_end_indices_.at(locus_idx);
                    }
                }
            }
        }
        for (int i = (ends_to_erase.size() - 1); i >= 0; --i) {
            this->locus_end_indices_.erase(this->locus_end_indices_.begin() + ends_to_erase.at(i));
        }
        for (int i = (site_indices_to_erase.size() - 1); i >= 0; --i) {
            this->contiguous_pattern_indices_.erase(this->contiguous_pattern_indices_.begin() + site_indices_to_erase.at(i));
        }
        // ECOEVOLITY_ASSERT(this->contiguous_pattern_indices_.size() == this->get_number_of_sites())
        // std::cout << "AFTER: Locus ends:\n";
        // for (unsigned int i = 0; i < this->locus_end_indices_.size(); ++i) {
        //     std::cout << this->locus_end_indices_.at(i) << " ";
        // }
        // std::cout << "\n";
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
        if (this->storing_seq_loci_info_) {
            for (unsigned int site_idx = 0;
                    site_idx < this->contiguous_pattern_indices_.size();
                    ++site_idx) {
                if (this->contiguous_pattern_indices_.at(site_idx) == mirrored_idx) {
                    this->contiguous_pattern_indices_.at(site_idx) = pattern_idx;
                }
            }
        }
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

unsigned int BiallelicData::get_number_of_triallelic_sites_recoded() const {
    return this->number_of_triallelic_sites_recoded_;
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
    if (this->storing_seq_loci_info_) {
        if (this->contiguous_pattern_indices_.size() != (this->locus_end_indices_.back() + 1)) {
            // Debug output
            // std::cout << "NUMBER of sites: "
            //         << this->get_number_of_sites()
            //         << "\n";
            // std::cout << "NUMBER of contiguous pattern indices: "
            //         << this->contiguous_pattern_indices_.size()
            //         << "\n";
            // std::cout << "NUMBER of locus end indices: "
            //         << this->locus_end_indices_.size()
            //         << "\n";
            // std::cout << "LAST locus end index: "
            //         << this->locus_end_indices_.back()
            //         << "\n";
            throw EcoevolityBiallelicDataError(
                    "The number of contiguous pattern indices does not match the end of the last locus",
                    this->path_);
        }
        if (this->has_constant_patterns() && (this->contiguous_pattern_indices_.size() != this->get_number_of_sites())) {
            throw EcoevolityBiallelicDataError(
                    "The number of contiguous pattern indices does not match the number of sites",
                    this->path_);
        }
    }
    return;
}

std::vector< std::vector<std::string> > BiallelicData::get_alignment() const {
    unsigned int npops = this->get_number_of_populations();
    std::vector< std::vector<std::string> > alignment;
    alignment.reserve(npops);
    for (unsigned int pop_idx = 0; pop_idx < npops; ++pop_idx) {
        unsigned int max_allele_count = this->get_max_allele_count(pop_idx);
        std::vector<std::string> row;
        row.reserve(max_allele_count);
        for (unsigned int allele_idx = 0;
                allele_idx < max_allele_count;
                ++allele_idx) {
            row.push_back("");
        }
        alignment.push_back(row);
    }
    if (this->storing_seq_loci_info_) {
        // Maintain the order of the site patterns
        for (auto const &pattern_idx : this->contiguous_pattern_indices_) {
            for (unsigned int pop_idx = 0; pop_idx < npops; ++pop_idx) {
                int n_red_alleles = this->get_red_allele_count(pattern_idx, pop_idx);
                int n_alleles = this->get_allele_count(pattern_idx, pop_idx);
                int n_zero_alleles = n_alleles - n_red_alleles;
                for (auto & row: alignment.at(pop_idx)) {
                    if (n_zero_alleles > 0) {
                        row += "0";
                        --n_zero_alleles;
                        --n_alleles;
                    }
                    else if (n_red_alleles > 0) {
                        row += "1";
                        --n_red_alleles;
                        --n_alleles;
                    }
                    else {
                        ECOEVOLITY_ASSERT(n_alleles == 0);
                        row += "?";
                    }
                }
            }
        }
        return alignment;
    }
    for (unsigned int pattern_idx = 0;
            pattern_idx < this->get_number_of_patterns();
            ++pattern_idx) {
        for (unsigned int pop_idx = 0; pop_idx < npops; ++pop_idx) {
            int n_red_alleles = this->get_red_allele_count(pattern_idx, pop_idx);
            int n_alleles = this->get_allele_count(pattern_idx, pop_idx);
            int n_zero_alleles = n_alleles - n_red_alleles;
            for (auto & row: alignment.at(pop_idx)) {
                if (n_zero_alleles > 0) {
                    row += std::string(this->get_pattern_weight(pattern_idx), '0');
                    --n_zero_alleles;
                    --n_alleles;
                }
                else if (n_red_alleles > 0) {
                    row += std::string(this->get_pattern_weight(pattern_idx), '1');
                    --n_red_alleles;
                    --n_alleles;
                }
                else {
                    ECOEVOLITY_ASSERT(n_alleles == 0);
                    row += std::string(this->get_pattern_weight(pattern_idx), '?');
                }
            }
        }
    }
    return alignment;
}

void BiallelicData::write_alignment(
        std::ostream& out,
        char population_name_delimiter) const {
    unsigned int count = 0;
    std::vector < std::vector<std::string> > alignment = this->get_alignment();
    unsigned int pop_idx = 0;
    for (auto rows: alignment) {
        unsigned int row_idx = 0;
        std::string pop_label = this->get_population_label(pop_idx);
        for (auto row: rows) {
            out << "\'" << pop_label << population_name_delimiter
                << string_util::pad_int(row_idx, 4) << "\'  "
                << row << std::endl;
            ++row_idx;
            ++count;
        }
        ++pop_idx;
    }
}

void BiallelicData::write_nexus(
        std::ostream& out,
        char population_name_delimiter) const {
    unsigned int ntax = 0;
    for (auto num_alleles: this->get_max_allele_counts()) {
        ntax += num_alleles;
    }
    unsigned int nchar = this->get_number_of_sites();
    out << "#NEXUS\n" 
        << "Begin data;\n"
        << "    Dimensions ntax=" << ntax << " nchar=" << nchar << ";\n"
        << "    Format datatype=standard symbols=\"01\" missing=? gap=-;\n"
        << "    Matrix\n";
    this->write_alignment(out, population_name_delimiter);
    out << "    ;\n"
        << "End;\n";
    if (this->storing_seq_loci_info_) {
        out << "\n";
        this->write_charsets(out);
    }
}

void BiallelicData::write_charsets(
        std::ostream& out) const {
    out << "Begin sets;\n";
    for (unsigned int locus_idx = 0;
            locus_idx < this->locus_end_indices_.size();
            ++locus_idx) {
        out << "    Charset locus"
            << locus_idx + 1
            << "=";
        if (locus_idx == 0) {
            out << "1-";
        }
        else {
            out << (this->locus_end_indices_.at(locus_idx - 1) + 2) << "-";
        }
        out << (this->locus_end_indices_.at(locus_idx) + 1) << ";\n";
    }
    out << "End;\n";
}

void BiallelicData::write_yaml(std::ostream& out) const {
    unsigned int pattern_idx;
    unsigned int pop_idx;
    out << std::boolalpha;
    out << "---\n";
    out << "markers_are_dominant: " << this->markers_are_dominant_ << "\n";
    out << "population_labels:\n";
    for (pop_idx = 0;
            pop_idx < this->get_number_of_populations();
            ++pop_idx) {
        out << "    - " << this->get_population_label(pop_idx) << "\n";
    }
    out << "allele_count_patterns:\n";
    for (pattern_idx = 0;
            pattern_idx < this->get_number_of_patterns();
            ++pattern_idx) {
        out << "    - [";
        for (pop_idx = 0;
                pop_idx < this->get_number_of_populations();
                ++pop_idx) {
            if (pop_idx > 0) {
                out << ", ";
            }
            out << "["
                << this->get_red_allele_count(pattern_idx, pop_idx)
                << ","
                << this->get_allele_count(pattern_idx, pop_idx)
                << "]"
        }
        out << "\n";
    }
    out << "pattern_weights:\n";
    for (pattern_idx = 0;
            pattern_idx < this->get_number_of_patterns();
            ++pattern_idx) {
        out << "    - " << this->get_pattern_weight(pattern_idx) << "\n";
    }
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
            << "Number of variable sites: "
            << this->get_number_of_variable_sites() << "\n"
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
    if (this->has_seq_loci_info()) {
        out << margin << indent
                << "Number of loci: "
                << this->locus_end_indices_.size() << "\n";
    }
}

std::map<std::vector<unsigned int>, unsigned int>
BiallelicData::get_unique_allele_counts() const {
    std::map<std::vector<unsigned int>, unsigned int> unique_allele_counts;
    for (unsigned int pattern_idx = 0; pattern_idx < this->get_number_of_patterns(); ++pattern_idx) {
        const std::vector<unsigned int>& allele_cts = this->get_allele_counts(pattern_idx);
        unsigned int w = this->get_pattern_weight(pattern_idx);
        if (unique_allele_counts.count(allele_cts) < 1) {
            unique_allele_counts[allele_cts] = w;
        }
        else {
            unique_allele_counts[allele_cts] += w;
        }
    }
    return unique_allele_counts;
}
