#ifndef ECOEVOLITY_ERROR_HPP
#define ECOEVOLITY_ERROR_HPP

#include <iostream>
#include <sstream>


/**
 * Base class for errors.
 */
class EcoevolityBaseError: public std::exception {
    public:
        mutable std::string error_report_;

        virtual ~EcoevolityBaseError() throw() { }

        EcoevolityBaseError(
                const std::string & name,
                const std::string & message);
        EcoevolityBaseError(
                const std::string & name,
                const std::string & message,
                const std::string & file_name);
        EcoevolityBaseError(
                const std::string & name,
                const std::string & message,
                const std::string & file_name,
                unsigned int line_number);
        EcoevolityBaseError(
                const std::string & name,
                const std::string & message,
                const std::string & file_name,
                const std::string & taxon_label,
                unsigned int character_index);

        const char * what() const throw() {
            return this->error_report_.empty() ? "EcoevolityError: no message" : this->error_report_.c_str();
        }
};

class EcoevolityError: public EcoevolityBaseError {
    public:
        EcoevolityError(
                const std::string & message) :
            EcoevolityBaseError("EcoevolityError", message) { }
};

class EcoevolityBiallelicDataError: public EcoevolityBaseError {
    public:
        EcoevolityBiallelicDataError(
                const std::string & message,
                const std::string & file_name) :
            EcoevolityBaseError("EcoevolityBiallelicDataError", message,
                    file_name) { }
};

class EcoevolityParsingError: public EcoevolityBaseError {
    public:
        EcoevolityParsingError(
                const std::string & message,
                const std::string & file_name,
                unsigned int line_number) :
            EcoevolityBaseError("EcoevolityParsingError", message, file_name,
                    line_number) { }
};

class EcoevolityInvalidCharacterError: public EcoevolityBaseError {
    public:
        EcoevolityInvalidCharacterError(
                const std::string & message,
                const std::string & file_name,
                const std::string & taxon_label,
                unsigned int character_index) :
            EcoevolityBaseError("EcoevolityInvalidCharacterError", message,
                    file_name, taxon_label, character_index) { }
};

#endif
