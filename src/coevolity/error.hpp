#ifndef COEVOLITY_EXCEPTION_HPP
#define COEVOLITY_EXCEPTION_HPP

#include <iostream>
#include <stdexcept>
#include <sstream>


/**
 * Base class for errors.
 */
class CoevolityError: public runtime_error {
    public:
        CoevolityError(const std::string message) : message_(message) { }

        virtual const char * what() const throw() {
            return this->get_error_message().c_str();
        }

        string get_error_message() const {
            std::ostringstream error_stream;
            error_stream << "CoevolityError: " << this->message_;
            return error_stream.str();
        }

    private:
        std::string message_;
};

class CoevolityParsingError: public CoevolityError {
    public:
        CoevolityParsingError(
                const std::string message,
                const std::string file_name,
                size_t line_number)
            : message_(message),
              line_number_(line_number),
              file_name_(file_name)
            { }

        string get_error_message() const {
            std::ostringstream error_stream;
            error_stream << "CoevolityParsingError" << std::endl
                         << "File: "  << this->file_name_   << std:endl
                         << "Line: "  << this->line_number_ << std:endl
                         << "Error: " << this->message_     << std:endl;
            return error_stream.str();
        }

    private:
        size_t line_number_;
        std::string file_name_;

#endif
