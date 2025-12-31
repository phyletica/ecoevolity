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

#pragma once

#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <set>
#include <string>
#include <vector>

#include "string_util.hpp"
#include "debug.hpp"
#include "error.hpp"


namespace ecoevolity {

/**
 * Class for storing nucleotide sequence data.
 *
 */
    class NucData {
        // TODO: Add `Partition` member and when this gets initiated, store
        // subset alignments in a vector
        public:
                                                    NucData();
                                                    ~NucData();

            void                                    init_from_phylip_stream(
                                                        std::istream& stream);
            void                                    init_from_phylip_path(
                                                        const std::string& path);

            unsigned                                get_num_seqs() const {return _num_seqs;}
            unsigned                                get_num_sites() const {return _num_sites;}
            unsigned                                get_num_states() const {return 4;}
            std::string                             get_seq(const std::string & label) const;

            std::vector< std::string >              get_seq_labels() const;

            void                                    clear();

            void                                    replace_label_underscores_with_spaces();
            
        private:
            std::unordered_map<std::string, std::string>    _sequences;
            unsigned                                        _num_seqs;
            unsigned                                        _num_sites;
            std::set<char> _states = {
                'A', 'a',
                'C', 'c',
                'G', 'g',
                'T', 't',
                'U', 'u',
                'W', 'w',
                'S', 's',
                'M', 'm',
                'K', 'k',
                'R', 'r',
                'Y', 'y',
                'B', 'b',
                'D', 'd',
                'H', 'h',
                'V', 'v',
                'N', 'n',
                'O', 'o',
                'X', 'x',
                '?', '-',
            };

        public:
            typedef std::shared_ptr< NucData >   SharedPtr;
    };

    inline NucData::NucData() {
        //std::cout << "Creating a NucData object" << std::endl;
        this->clear();
    }

    inline NucData::~NucData() {
        //std::cout << "Destroying a NucData object" << std::endl;
    }

    inline void NucData::clear() {
        this->_sequences.clear();
        this->_num_seqs = 0;
        this->_num_sites = 0;
    }

    inline std::string NucData::get_seq(const std::string & label) const {
        return this->_sequences.at(label);
    }

    inline std::vector<std::string> NucData::get_seq_labels() const {
        std::vector<std::string> labels;
        labels.reserve(this->_sequences.size());
        for (const auto& pair: this->_sequences) {
            labels.push_back(pair.first);
        }
        return labels;
    }

    inline void NucData::init_from_phylip_path(const std::string & path) {
        std::ifstream in_stream;
        in_stream.open(path);
        if (! in_stream.is_open()) {
            std::ostringstream msg;
            msg << "Could not open phylip file: "
                << path << std::endl;
            throw EcoevolityError(msg.str());
        }
        try {
            this->init_from_phylip_stream(in_stream);
        }
        catch (...) {
            std::cerr << "Problem parsing phylip file: " << path << std::endl;
            throw;
        }
        in_stream.close();
    }

    inline void NucData::init_from_phylip_stream(std::istream & stream) {
        this->clear();
        std::string line;
        bool parsed_header = false;
        std::string stripped_line;
        std::vector<std::string> line_elements;
        unsigned num_seqs = 0;
        unsigned num_sites = 0;
        unsigned line_idx = 0;
        while (std::getline(stream, line)) {
            stripped_line = string_util::strip(line);
            if (stripped_line.empty()) {
                continue;
            }

            line_elements = string_util::split_wspace(stripped_line);

            if (line_elements.size() != 2) {
                std::ostringstream msg;
                msg << "Line " << line_idx + 1 << " has " << line_elements.size() << " columns; expecting 2" << std::endl;
                throw EcoevolityError(msg.str());
            }

            if (! parsed_header) {
                // Found first 2-column line; should be phylip header
                int nseqs = std::stoi(line_elements.at(0));
                int nsites = std::stoi(line_elements.at(1));
                if (nseqs < 1) {
                    std::ostringstream msg;
                    msg << "Number of sequences must be positive; found: " << nseqs << std::endl;
                    throw EcoevolityError(msg.str());
                }
                if (nsites < 1) {
                    std::ostringstream msg;
                    msg << "Number of sites must be positive; found: " << nsites << std::endl;
                    throw EcoevolityError(msg.str());
                }
                num_seqs = (unsigned)nseqs;
                num_sites = (unsigned)nsites;

                parsed_header = true;

                continue;
            }
            // Should have a non-empty line with label and sequence

            std::string seq_label = line_elements.at(0);
            std::string seq = line_elements.at(1);

            if (seq.size() != num_sites) {
                std::ostringstream msg;
                msg << "Sequence on line " << line_idx + 1 << " has " << seq.size() << " sites; expecting " << num_sites << " sites. Here's the sequence label:" << std::endl << "  " << seq_label << std::endl;
                throw EcoevolityError(msg.str());
            }

            unsigned site_idx = 0;
            for (char state : seq) {
                bool is_valid = _states.find(state) != _states.end();
                if (! is_valid) {
                    std::ostringstream msg;
                    msg << "Found invalid character " << state << " at site " << site_idx + 1 << " of line " << line_idx + 1 << " (label: " << seq_label << ")" << std::endl; 
                    throw(EcoevolityError(msg.str()));
                }
                ++site_idx;
            }

            this->_sequences[seq_label] = seq;
            ++line_idx;
        }
        if (this->_sequences.size() != num_seqs) {
            std::ostringstream msg;
            msg << "Expecting " << num_seqs << " sequences, but found " << _sequences.size() << std::endl;
            throw EcoevolityError(msg.str());
        }
        this->_num_seqs = num_seqs;
        this->_num_sites = num_sites;
    }

    inline void NucData::replace_label_underscores_with_spaces() {
        if (this->_sequences.empty()) {
            return;
        }

        std::unordered_map<std::string, std::string> new_seqs;
        for (const auto& pair: this->_sequences) {
            std::string label = pair.first;
            std::string new_label = string_util::join(string_util::split(label, '_'), " ");
            new_seqs[new_label] = pair.second;
        }
        this->_sequences = new_seqs;
    }
}
