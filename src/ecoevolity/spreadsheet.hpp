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

namespace spreadsheet {

inline std::vector<std::string> parse_header(
        std::istream& in_stream,
        char delimiter = '\t') {
    std::string header_line;
    std::getline(in_stream, header_line);
    std::vector<std::string> header = string_util::split(
            header_line, delimiter);
    return header;
}

inline void parse(
        std::istream& in_stream,
        std::map<std::string, std::vector<std::string> >& column_data,
        unsigned int offset = 0,
        char delimiter = '\t') {
    std::string line;
    unsigned int line_index = 0;
    std::vector<std::string> header = parse_header(in_stream);
    if (header.size() == 0) {
        throw EcoevolityParsingError(
                "Could not parse header",
                line_index + 1);
    }
    unsigned int number_of_columns = header.size();
    std::vector<std::string> elements;
    if (column_data.size() > 0) {
        if (column_data.size() != header.size()) {
            throw EcoevolityParsingError(
                    "Columns do not match existing column data map",
                    line_index + 1);
        }
        for (unsigned int i = 0; i < header.size(); ++i) {
            if (column_data.count(header.at(i)) < 1) {
                throw EcoevolityParsingError(
                        "Column header \'" + header.at(i) +
                        "\' not found in existing column data map",
                    line_index + 1);
            }
        }
    }
    else {
        for (unsigned int i = 0; i < header.size(); ++i) {
            column_data[header.at(i)];
        }
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
        parse(in_stream, column_data, offset, delimiter);
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
        unsigned int offset = 0,
        char delimiter = '\t') {
    for (unsigned int i = 0; i < paths.size(); ++i) {
        parse(paths.at(i), column_data, offset, delimiter);
    }
}

class Spreadsheet {

    public:
        std::map<std::string, std::vector<std::string> > data;

        void update(
                std::istream& in_stream,
                unsigned int offset = 0,
                char delimiter = '\t')
        {
            parse(in_stream, this->data, offset, delimiter);

        }

        void update(
                const std::string& path,
                unsigned int offset = 0,
                char delimiter = '\t')
        {
            parse(path, this->data, offset, delimiter);

        }

        void update(
                const std::vector<std::string>& paths,
                unsigned int offset = 0,
                char delimiter = '\t')
        {
            parse(paths, this->data, offset, delimiter);

        }

        template <typename T>
        void get(
                const std::string& column_label,
                std::vector<T>& target)
        {
            T value;
            for (unsigned int i = 0;
                    i < this->data.at(column_label).size();
                    ++i) {
                std::stringstream converter(this->data.at(column_label).at(i));
                if (! (converter >> value)) {
                    throw EcoevolitySpreadsheetError("could not convert \'" +
                            converter.str() + "\'");
                }
                target.push_back(value);
            }
        }

        template <typename T>
        std::vector<T> get(
                const std::string& column_label)
        {
            std::vector<T> r;
            r.reserve(this->data.at(column_label).size());
            this->get<T>(column_label, r);
            return r;
        }

        // make sure all columns have same number of samples
        void validate() const {
        }
};

} // namespace spreadsheet 

#endif
