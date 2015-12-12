#ifndef COEVOLITY_ERROR_HPP
#define COEVOLITY_ERROR_HPP

#include <iostream>
#include <stdexcept>
#include <sstream>


/**
 * Base class for errors.
 */
class CoevolityError: public std::runtime_error {
    public:
        CoevolityError(const std::string message) :
            std::runtime_error(""),
            message_(message)
            { }

        virtual const char * what() const throw() {
            return this->get_error_message().c_str();
        }

        std::string get_error_message() const {
            std::ostringstream error_stream;
            error_stream << "CoevolityError: " << this->message_;
            return error_stream.str();
        }

    protected:
        std::string message_;
};

class CoevolityParsingError: public CoevolityError {
    public:
        CoevolityParsingError(
                const std::string message,
                const std::string file_name,
                size_t line_number) :
            CoevolityError(message),
            line_number_(line_number),
            file_name_(file_name)
            { }

        std::string get_error_message() const {
            std::ostringstream error_stream;
            error_stream << "CoevolityParsingError" << std::endl
                         << "File: "  << this->file_name_   << std::endl
                         << "Line: "  << this->line_number_ << std::endl
                         << "Error: " << this->message_     << std::endl;
            return error_stream.str();
        }

    private:
        size_t line_number_;
        std::string file_name_;
};

#endif
