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

#ifndef ECOEVOLITY_SPREADSHEET_HPP
#define ECOEVOLITY_SPREADSHEET_HPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>

#include "assert.hpp"
#include "error.hpp"
#include "string_util.hpp"
#include "stats_util.hpp"

namespace spreadsheet {

inline void parse_header(
        std::istream& in_stream,
        std::vector<std::string>& header,
        char delimiter = '\t') {
    std::string header_line;
    std::getline(in_stream, header_line);
    header = string_util::split(
            header_line, delimiter);
}

inline void parse(
        std::istream& in_stream,
        std::map<std::string, std::vector<std::string> >& column_data,
        std::vector<std::string>& header,
        unsigned int offset = 0,
        char delimiter = '\t') {
    std::string line;
    unsigned int line_index = 0;
    if (header.size() == 0) {
        parse_header(in_stream, header, delimiter);
    }
    else {
        std::vector<std::string> new_header;
        parse_header(in_stream, new_header, delimiter);
        if (new_header != header) {
            throw EcoevolityParsingError(
                    "Headers does not match",
                    line_index + 1);
        }
    }
    if (header.size() == 0) {
        throw EcoevolityParsingError(
                "Could not parse header",
                line_index + 1);
    }
    unsigned int number_of_columns = header.size();
    std::vector<std::string> elements;
    for (const auto &h: header) {
        column_data[h];
    }
    elements.reserve(header.size());
    while (std::getline(in_stream, line)) {
        if (line_index < offset) {
            ++line_index;
            continue;
        }
        elements.clear();
        string_util::split(line, delimiter, elements);
        if (elements.size() != number_of_columns) {
            std::ostringstream message;
            message << "Incorrect number of columns: Expecting "
                    << number_of_columns << ", but found "
                    << elements.size();
            throw EcoevolityParsingError(
                    message.str(),
                    line_index + 1);
        }
        for (unsigned int i = 0; i < elements.size(); ++i) {
            column_data.at(header.at(i)).push_back(elements.at(i));
        }
        ++line_index;
    }
}

inline void parse(
        const std::string& path,
        std::map<std::string, std::vector<std::string> >& column_data,
        std::vector<std::string>& header,
        unsigned int offset = 0,
        char delimiter = '\t') {
    std::ifstream in_stream;
    in_stream.open(path);
    if (! in_stream.is_open()) {
        throw EcoevolityParsingError(
                "Could not open spreadsheet file",
                path);
    }
    try {
        parse(in_stream, column_data, header, offset, delimiter);
    }
    catch (...) {
        std::cerr << "ERROR: Problem parsing spreadsheet \'"
                  << path << "\'\n";
        throw;
    }
    in_stream.close();
}

inline void parse(
        const std::vector<std::string>& paths,
        std::map<std::string, std::vector<std::string> >& column_data,
        std::vector<std::string>& header,
        unsigned int offset = 0,
        char delimiter = '\t') {
    for (const auto &p: paths) {
        parse(p, column_data, header, offset, delimiter);
    }
}

class Spreadsheet {

    protected:

        std::map<std::string, std::vector<std::string> > data_;
        std::vector<std::string> header_;

    public:

        void update(
                std::istream& in_stream,
                unsigned int offset = 0,
                char delimiter = '\t')
        {
            parse(in_stream, this->data_, this->header_, offset, delimiter);

        }

        void update(
                const std::string& path,
                unsigned int offset = 0,
                char delimiter = '\t')
        {
            parse(path, this->data_, this->header_, offset, delimiter);

        }

        void update(
                const std::vector<std::string>& paths,
                unsigned int offset = 0,
                char delimiter = '\t')
        {
            parse(paths, this->data_, this->header_, offset, delimiter);

        }

        const std::vector<std::string>& get_keys() const {
            return this->header_;
        }
        const std::map<std::string, std::vector<std::string> >& get_data() const {
            return this->data_;
        }

        template <typename T>
        void get(
                const std::string& column_label,
                std::vector<T>& target) const
        {
            T value;
            for (const auto &s: this->data_.at(column_label)) {
                std::stringstream converter(s);
                if (! (converter >> value)) {
                    throw EcoevolitySpreadsheetError("could not convert \'" +
                            converter.str() + "\'");
                }
                target.push_back(value);
            }
        }

        template <typename T>
        std::vector<T> get(
                const std::string& column_label) const
        {
            std::vector<T> r;
            r.reserve(this->data_.at(column_label).size());
            this->get<T>(column_label, r);
            return r;
        }

        template <typename T>
        SampleSummarizer<T> summarize(const std::string& column_label) const {
            SampleSummarizer<T> summarizer;
            T value;
            for (const auto &s: this->data_.at(column_label)) {
                std::stringstream converter(s);
                if (! (converter >> value)) {
                    throw EcoevolitySpreadsheetError("could not convert \'" +
                            converter.str() + "\'");
                }
                summarizer.add_sample(value);
            }
            return summarizer;
        }

        bool has_key(const std::string& k) const {
            return (this->data_.count(k) > 0);
        }

};

} // namespace spreadsheet 

#endif
