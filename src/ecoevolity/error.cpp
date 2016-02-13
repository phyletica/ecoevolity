#include "error.hpp"

EcoevolityBaseError::EcoevolityBaseError(
        const std::string & name,
        const std::string & message) {
    std::ostringstream oss;
    oss << name << ": " << message << std::endl;
    this->error_report_ = oss.str();
}

EcoevolityBaseError::EcoevolityBaseError(
        const std::string & name,
        const std::string & message,
        const std::string & file_name) {
    std::ostringstream oss;
    oss << name << ":" << std::endl
        << "    File: " << file_name << std::endl
        << "    Error: " << message << std::endl;
    this->error_report_ = oss.str();
}

EcoevolityBaseError::EcoevolityBaseError(
        const std::string & name,
        const std::string & message,
        const std::string & file_name,
        unsigned int line_number) {
    std::ostringstream oss;
    oss << name << ":" << std::endl
        << "    File: " << file_name << std::endl
        << "    Line: " << line_number << std::endl
        << "    Error: " << message << std::endl;
    this->error_report_ = oss.str();
}

EcoevolityBaseError::EcoevolityBaseError(
        const std::string & name,
        const std::string & message,
        const std::string & file_name,
        const std::string & taxon_label,
        unsigned int character_index) {
    std::ostringstream oss;
    oss << name << ":" << std::endl
        << "    File: " << file_name << std::endl
        << "    Taxon: " << taxon_label << std::endl
        << "    Character: " << character_index + 1 << std::endl
        << "    Error: " << message << std::endl;
    this->error_report_ = oss.str();
}
