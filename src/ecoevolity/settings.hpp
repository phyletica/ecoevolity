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

// TODO: Error messages could be a lot more informative.

#ifndef ECOEVOLITY_SETTINGS_HPP
#define ECOEVOLITY_SETTINGS_HPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <unordered_map>

#include "yaml-cpp/yaml.h"

#include "error.hpp"
#include "assert.hpp"
#include "string_util.hpp"
#include "path.hpp"
#include "math_util.hpp"
#include "probability.hpp"
#include "data.hpp"


class YamlCppUtils {

    public:

        static std::string get_node_type(const YAML::Node& node) {
            switch(node.Type()) {
                case YAML::NodeType::Null:
                    return "Null";
                case YAML::NodeType::Scalar:
                    return "Scalar";
                case YAML::NodeType::Sequence:
                    return "Sequence";
                case YAML::NodeType::Map:
                    return "Map";
                case YAML::NodeType::Undefined:
                    return "Undefined";
                default:
                    return "Not recognized";
            }
        }
};

class ContinuousDistributionSettings {

    private:
        std::string name_ = "none";
        std::unordered_map<std::string, double> parameters_;

    public:
        ContinuousDistributionSettings() { };
        ContinuousDistributionSettings(
                const std::string& name,
                const std::unordered_map<std::string, double>& parameters) {
            this->name_ = name;
            // Gamma
            if (name == "gamma_distribution") {
                if (parameters.count("shape") < 1) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "shape parameter missing for gamma_distribution"
                            );
                }
                if (parameters.count("scale") < 1) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "scale parameter missing for gamma_distribution"
                            );
                }
                if (parameters.count("offset") < 1) {
                    if (parameters.size() > 2) {
                        throw EcoevolityContinuousDistributionSettingError(
                                "unrecognized parameters for gamma_distribution (recognized parameters: shape, scale, offset)"
                                );
                    }
                }
                else {
                    if (parameters.size() > 3) {
                        throw EcoevolityContinuousDistributionSettingError(
                                "unrecognized parameters for gamma_distribution (recognized parameters: shape, scale, offset)"
                                );
                    }
                }
                if ((parameters.at("shape") <= 0.0) || (parameters.at("scale") <= 0.0)) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "Shape and scale must be greater than zero for gamma distribution");
                }
                this->parameters_["shape"] = parameters.at("shape");
                this->parameters_["scale"] = parameters.at("scale");
                if (parameters.count("offset") == 1) {
                    this->parameters_["offset"] = parameters.at("offset");
                }
            }
            // Exponential
            else if (name == "exponential_distribution") {
                if (parameters.count("rate") < 1) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "rate parameter missing for exponential_distribution"
                            );
                }
                if (parameters.count("offset") < 1) {
                    if (parameters.size() > 1) {
                        throw EcoevolityContinuousDistributionSettingError(
                                "unrecognized parameters for exponential_distribution (recognized parameters: rate, offset)"
                                );
                    }
                }
                else {
                    if (parameters.size() > 2) {
                        throw EcoevolityContinuousDistributionSettingError(
                                "unrecognized parameters for exponential_distribution (recognized parameters: rate, offset)"
                                );
                    }
                }
                if (parameters.at("rate") <= 0.0) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "rate must be greater than 0 for exponential distribution");
                }
                this->parameters_["rate"] = parameters.at("rate");
                if (parameters.count("offset") == 1) {
                    this->parameters_["offset"] = parameters.at("offset");
                }
            }
            // Uniform
            else if (name == "uniform_distribution") {
                if (parameters.count("min") < 1) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "min parameter missing for uniform_distribution"
                            );
                }
                if (parameters.count("max") < 1) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "max parameter missing for uniform_distribution"
                            );
                }
                if (parameters.size() > 2) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "unrecognized parameters for uniform_distribution (recognized parameters: min, max)"
                            );
                }
                if (parameters.at("max") <= parameters.at("min")) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "The upper limit must be greater than lower limit for UniformDistribution");
                }
                this->parameters_["min"] = parameters.at("min");
                this->parameters_["max"] = parameters.at("max");
            }
            else if (name == "none") {
                this->parameters_.clear();
            }
            else {
                std::string message = "unrecognized distribution: " + name;
                throw EcoevolityContinuousDistributionSettingError(message);
            }
        }
        ContinuousDistributionSettings(const YAML::Node& node) {
            if (! node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "continuous distribution node should be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
            }
            if (node.size() != 1) {
                throw EcoevolityYamlConfigError(
                        "continuous distribution node should only have a single key");
            }
            // Gamma
            if (node["gamma_distribution"]) {
                this->name_ = "gamma_distribution";
                YAML::Node parameters = node["gamma_distribution"];
                if (! parameters["shape"]) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "shape parameter missing for gamma_distribution"
                            );
                }
                if (! parameters["scale"]) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "scale parameter missing for gamma_distribution"
                            );
                }
                if (! parameters["offset"]) {
                    if (parameters.size() > 2) {
                        throw EcoevolityContinuousDistributionSettingError(
                                "unrecognized parameters for gamma_distribution (recognized parameters: shape, scale, offset)"
                                );
                    }
                }
                else {
                    if (parameters.size() > 3) {
                        throw EcoevolityContinuousDistributionSettingError(
                                "unrecognized parameters for gamma_distribution (recognized parameters: shape, scale, offset)"
                                );
                    }
                }
                this->parameters_["shape"] = parameters["shape"].as<double>();
                this->parameters_["scale"] = parameters["scale"].as<double>();
                if ((this->parameters_.at("shape") <= 0.0) || (this->parameters_.at("scale") <= 0.0)) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "Shape and scale must be greater than zero for gamma distribution");
                }
                if (parameters["offset"]) {
                    this->parameters_["offset"] = parameters["offset"].as<double>();
                }
            }
            // Exponential
            else if (node["exponential_distribution"]) {
                this->name_ = "exponential_distribution";
                YAML::Node parameters = node["exponential_distribution"];
                if (! parameters["rate"]) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "rate parameter missing for exponential_distribution"
                            );
                }
                if (! parameters["offset"]) {
                    if (parameters.size() > 1) {
                        throw EcoevolityContinuousDistributionSettingError(
                                "unrecognized parameters for exponential_distribution (recognized parameters: rate, offset)"
                                );
                    }
                }
                else {
                    if (parameters.size() > 2) {
                        throw EcoevolityContinuousDistributionSettingError(
                                "unrecognized parameters for exponential_distribution (recognized parameters: rate, offset)"
                                );
                    }
                }
                this->parameters_["rate"] = parameters["rate"].as<double>();
                if (this->parameters_.at("rate") <= 0.0) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "rate must be greater than 0 for exponential distribution");
                }
                if (parameters["offset"]) {
                    this->parameters_["offset"] = parameters["offset"].as<double>();
                }
            }
            // Uniform
            else if (node["uniform_distribution"]) {
                this->name_ = "uniform_distribution";
                YAML::Node parameters = node["uniform_distribution"];
                if (! parameters["min"]) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "min parameter missing for uniform_distribution"
                            );
                }
                if (! parameters["max"]) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "max parameter missing for uniform_distribution"
                            );
                }
                if (parameters.size() > 2) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "unrecognized parameters for uniform_distribution (recognized parameters: min, max)"
                            );
                }
                this->parameters_["min"] = parameters["min"].as<double>();
                this->parameters_["max"] = parameters["max"].as<double>();
                if (this->parameters_.at("max") <= this->parameters_.at("min")) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "The upper limit must be greater than lower limit for UniformDistribution");
                }
            }
            else {
                std::string message = "unrecognized distribution: " + node.begin()->first.as<std::string>();
                throw EcoevolityContinuousDistributionSettingError(message);
            }
        }
        virtual ~ContinuousDistributionSettings() { }
        ContinuousDistributionSettings& operator=(const ContinuousDistributionSettings& other) {
            this->name_ = other.name_;
            this->parameters_ = other.parameters_;
            return * this;
        }

        const std::string& get_name() const {
            return this->name_;
        }

        void nullify() {
            this->name_ = "none";
            this->parameters_.clear();
        }

        std::shared_ptr<ContinuousProbabilityDistribution> get_instance() const {
            std::shared_ptr<ContinuousProbabilityDistribution> p;
            if (this->name_ == "gamma_distribution") {
                if (this->parameters_.count("offset") > 0) { 
                    p = std::make_shared<OffsetGammaDistribution>(
                            this->parameters_.at("shape"),
                            this->parameters_.at("scale"),
                            this->parameters_.at("offset"));
                }
                else {
                    p = std::make_shared<GammaDistribution>(
                            this->parameters_.at("shape"),
                            this->parameters_.at("scale"));
                }
            }
            else if (this->name_ == "exponential_distribution") {
                if (this->parameters_.count("offset") > 0) { 
                    p = std::make_shared<OffsetExponentialDistribution>(
                            this->parameters_.at("rate"),
                            this->parameters_.at("offset"));
                }
                else {
                    p = std::make_shared<ExponentialDistribution>(
                            this->parameters_.at("rate"));
                }
            }
            else if (this->name_ == "uniform_distribution") {
                p = std::make_shared<UniformDistribution>(
                        this->parameters_.at("min"),
                        this->parameters_.at("max"));
            }
            else if (this->name_ == "none") {
                return p;
            }
            else {
                ECOEVOLITY_ASSERT(0 == 1);
            }
            return p;
        }

        std::string to_string(unsigned int indent_level = 0) const {
            if (this->name_ == "none") {
                return "";
            }
            std::ostringstream ss;
            std::string margin = string_util::get_indent(indent_level);
            std::string indent = string_util::get_indent(1);
            ss << margin << this->name_ << ":\n";
            if (this->name_ == "gamma_distribution") {
                ss << margin << indent << "shape: " << this->parameters_.at("shape") << "\n";
                ss << margin << indent << "scale: " << this->parameters_.at("scale") << "\n";
                if (this->parameters_.count("offset") > 0) {
                    ss << margin << indent << "offset: " << this->parameters_.at("offset") << "\n";
                }
            }
            else if (this->name_ == "exponential_distribution") {
                ss << margin << indent << "rate: " << this->parameters_.at("rate") << "\n";
                if (this->parameters_.count("offset") > 0) {
                    ss << margin << indent << "offset: " << this->parameters_.at("offset") << "\n";
                }
            }
            else if (this->name_ == "uniform_distribution") {
                ss << margin << indent << "min: " << this->parameters_.at("min") << "\n";
                ss << margin << indent << "max: " << this->parameters_.at("max") << "\n";
            }
            else {
                ECOEVOLITY_ASSERT(0 == 1);
            }
            return ss.str();
        }
};

class PositiveRealParameterSettings {

    friend class ComparisonSettings;
    friend class CollectionSettings;

    private:
        double value_ = std::numeric_limits<double>::quiet_NaN();
        bool is_fixed_ = false;
        ContinuousDistributionSettings prior_settings_;

    public:
        PositiveRealParameterSettings() { }
        PositiveRealParameterSettings(
                double value,
                bool fix,
                const std::string& prior_name,
                const std::unordered_map<std::string, double>& prior_parameters) {
            if (value < 0.0) {
                throw EcoevolityPositiveRealParameterSettingError(
                        "positive real parameter cannot be less than 0"
                        );
            }
            this->value_ = value;
            this->is_fixed_ = fix;
            if (fix) {
                this->prior_settings_ = ContinuousDistributionSettings("none",
                        prior_parameters);
            }
            else {
                this->prior_settings_ = ContinuousDistributionSettings(prior_name,
                        prior_parameters);
            }

            if ((this->is_fixed_) && (std::isnan(this->value_))) {
                throw EcoevolityPositiveRealParameterSettingError(
                        "cannot fix parameter without a value"
                        );
            }
        }
        PositiveRealParameterSettings(const YAML::Node& node) {
            if (! node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "parameter node should be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
            }
            if (node.size() < 1) {
                throw EcoevolityYamlConfigError(
                        "empty parameter node");
            }
            for (YAML::const_iterator arg = node.begin();
                    arg != node.end();
                    ++arg) {
                if (arg->first.as<std::string>() == "value") {
                    double v = arg->second.as<double>();
                        if (v < 0.0) {
                            throw EcoevolityPositiveRealParameterSettingError(
                                    "positive real parameter cannot be less than 0"
                                    );
                        }
                    this->value_ = v;
                }
                else if (arg->first.as<std::string>() == "estimate") {
                    bool f = arg->second.as<bool>();
                    this->is_fixed_ = (! f);
                }
                else if (arg->first.as<std::string>() == "prior") {
                    this->prior_settings_ = ContinuousDistributionSettings(
                            arg->second);
                }
                else {
                    std::string message = "unrecognized parameter key: " +
                            arg->first.as<std::string>();
                    throw EcoevolityContinuousDistributionSettingError(message);
                }
            }

            if (this->is_fixed_) {
                std::unordered_map<std::string, double> prior_parameters;
                this->prior_settings_ = ContinuousDistributionSettings("none",
                        prior_parameters);
            }
            if ((this->is_fixed_) && (std::isnan(this->value_))) {
                throw EcoevolityPositiveRealParameterSettingError(
                        "cannot fix parameter without a value"
                        );
            }
        }
        virtual ~PositiveRealParameterSettings() { }
        PositiveRealParameterSettings& operator=(const PositiveRealParameterSettings& other) {
            this->value_ = other.value_;
            this->is_fixed_ = other.is_fixed_;
            this->prior_settings_ = other.prior_settings_;
            return * this;
        }

        double get_value() const {
            return this->value_;
        }
        bool is_fixed() const {
            return this->is_fixed_;
        }
        const ContinuousDistributionSettings& get_prior_settings() const {
            return this->prior_settings_;
        }

        virtual std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            ss << std::boolalpha;
            std::string margin = string_util::get_indent(indent_level);
            if (! std::isnan(this->get_value())) {
                ss << margin << "value: " << this->get_value() << "\n";
            }
            ss << margin << "estimate: " << (! this->is_fixed()) << "\n";
            if (! this->is_fixed()) {
                ss << margin << "prior:\n";
                ss << this->prior_settings_.to_string(indent_level + 1);
            }
            return ss.str();
        }
};

class ComparisonSettings {

    friend class CollectionSettings;

    private:
        std::string path_;
        PositiveRealParameterSettings population_size_settings_;
        PositiveRealParameterSettings u_settings_;
        PositiveRealParameterSettings time_multiplier_settings_;

        char population_name_delimiter_ = '_';
        bool population_name_is_prefix_ = true;
        bool genotypes_are_diploid_ = true;
        bool markers_are_dominant_ = false;
        bool constant_sites_removed_ = true;

        bool use_empirical_mutation_rate_starting_values_ = false;
        bool constrain_population_sizes_ = false;
        bool constrain_mutation_rates_ = true;

        void make_consistent() {
            if (this->constrain_mutation_rates_) {
                this->use_empirical_mutation_rate_starting_values_ = false;
                this->u_settings_.prior_settings_.nullify();
                this->u_settings_.value_ = 1.0;
                this->u_settings_.is_fixed_ = true;
            }
        }

        void init_empirical_mutation_rate_starting_values() {
            if (this->use_empirical_mutation_rate_starting_values_) {
                BiallelicData d = BiallelicData(
                        this->path_,
                        this->population_name_delimiter_,
                        this->population_name_is_prefix_,
                        this->genotypes_are_diploid_,
                        this->markers_are_dominant_,
                        true);
                double u;
                double v;
                d.get_empirical_mutation_rates(u, v);
                this->u_settings_.value_ = u;
            }
        }

        void parse_parameter_settings(const YAML::Node& node) {
            if (! node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "comparison parameters node should be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
            }
            if (node.size() < 1) {
                throw EcoevolityYamlConfigError(
                        "empty comparison parameters node");
            }
            for (YAML::const_iterator parameter = node.begin();
                    parameter != node.end();
                    ++parameter) {
                if (parameter->first.as<std::string>() == "population_size") {
                    this->population_size_settings_ = PositiveRealParameterSettings(parameter->second);
                }
                else if (parameter->first.as<std::string>() == "u_rate") {
                    this->u_settings_ = PositiveRealParameterSettings(parameter->second);
                }
                else if (parameter->first.as<std::string>() == "time_multiplier") {
                    this->time_multiplier_settings_ = PositiveRealParameterSettings(parameter->second);
                }
                else {
                    std::string message = "Unrecognized comparison parameter: " +
                            parameter->first.as<std::string>();
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

        void update_from_config(
                const YAML::Node& comparison_node,
                const std::string& config_path,
                bool global_defaults = false) {
            if (! comparison_node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "comparison node should be a map, but found: " +
                        YamlCppUtils::get_node_type(comparison_node));
            }
            if (comparison_node.size() < 1) {
                throw EcoevolityYamlConfigError(
                        "empty comparison node");
            }
            this->path_ = "";
            for (YAML::const_iterator arg = comparison_node.begin();
                    arg != comparison_node.end();
                    ++arg) {
                if (arg->first.as<std::string>() == "path") {
                    std::string p = string_util::strip(arg->second.as<std::string>());
                    this->path_ = path::join(path::dirname(config_path), p);
                }
                else if (arg->first.as<std::string>() == "genotypes_are_diploid") {
                    this->genotypes_are_diploid_ = arg->second.as<bool>();
                }
                else if (arg->first.as<std::string>() == "markers_are_dominant") {
                    this->markers_are_dominant_ = arg->second.as<bool>();
                }
                else if (arg->first.as<std::string>() == "population_name_delimiter") {
                    this->population_name_delimiter_ = arg->second.as<char>();
                }
                else if (arg->first.as<std::string>() == "population_name_is_prefix") {
                    this->population_name_is_prefix_ = arg->second.as<bool>();
                }
                else if (arg->first.as<std::string>() == "constant_sites_removed") {
                    this->constant_sites_removed_ = arg->second.as<bool>();
                }
                else if (arg->first.as<std::string>() == "use_empirical_mutation_rate_starting_values") {
                    this->use_empirical_mutation_rate_starting_values_ = arg->second.as<bool>();
                }
                else if (arg->first.as<std::string>() == "constrain_population_sizes") {
                    this->constrain_population_sizes_ = arg->second.as<bool>();
                }
                else if (arg->first.as<std::string>() == "constrain_mutation_rates") {
                    this->constrain_mutation_rates_ = arg->second.as<bool>();
                }
                else if (arg->first.as<std::string>() == "parameters") {
                    this->parse_parameter_settings(arg->second);
                }
                else {
                    std::string message = "Unrecognized comparison key: " +
                            arg->first.as<std::string>();
                    throw EcoevolityYamlConfigError(message);
                }
            }
            if ((this->path_ == "") && (! global_defaults)) {
                throw EcoevolityYamlConfigError("Every comparison must include a path");
            }
            this->make_consistent();
            if (! global_defaults) {
                this->init_empirical_mutation_rate_starting_values();
            }
        }

    public:
        ComparisonSettings() { }
        ComparisonSettings(
                const std::string& path,
                const PositiveRealParameterSettings& population_size_settings,
                const PositiveRealParameterSettings& u_settings,
                const PositiveRealParameterSettings& time_multiplier_settings,
                char population_name_delimiter = '_',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool use_empirical_mutation_rate_starting_values = false,
                bool constrain_population_sizes = false,
                bool constrain_mutation_rates = true) {

            this->path_ = path;
            this->population_size_settings_ = population_size_settings;
            this->u_settings_ = u_settings;
            this->time_multiplier_settings_ = time_multiplier_settings;
            this->population_name_delimiter_ = population_name_delimiter;
            this->population_name_is_prefix_ = population_name_is_prefix;
            this->genotypes_are_diploid_ = genotypes_are_diploid;
            this->markers_are_dominant_ = markers_are_dominant;
            this->constant_sites_removed_ = constant_sites_removed;
            this->use_empirical_mutation_rate_starting_values_ = use_empirical_mutation_rate_starting_values;
            this->constrain_population_sizes_ = constrain_population_sizes;
            this->constrain_mutation_rates_ = constrain_mutation_rates;
            this->make_consistent();

            // TODO:
            // Not very efficient to parse data just to get rates, but might be
            // more awkward than it's worth to defer it.
            // Make a copy operator for BiallelicData and parse and store here, then
            // can copy it in get_instance method
            this->init_empirical_mutation_rate_starting_values();
        }
        ComparisonSettings(
                const YAML::Node& comparison_node,
                const std::string& config_path,
                bool global_defaults = false) {
            this->update_from_config(comparison_node, config_path, global_defaults);
        }

        virtual ~ComparisonSettings() { }
        ComparisonSettings& operator=(const ComparisonSettings& other) {
            this->path_                                         = other.path_;
            this->population_size_settings_                     = other.population_size_settings_;
            this->u_settings_                                   = other.u_settings_;
            this->time_multiplier_settings_                     = other.time_multiplier_settings_;
            this->population_name_delimiter_                    = other.population_name_delimiter_;
            this->population_name_is_prefix_                    = other.population_name_is_prefix_;
            this->genotypes_are_diploid_                        = other.genotypes_are_diploid_;
            this->markers_are_dominant_                         = other.markers_are_dominant_;
            this->constant_sites_removed_                       = other.constant_sites_removed_;
            this->use_empirical_mutation_rate_starting_values_  = other.use_empirical_mutation_rate_starting_values_;
            this->constrain_population_sizes_                   = other.constrain_population_sizes_;
            this->constrain_mutation_rates_                     = other.constrain_mutation_rates_;
            return * this;
        }

        const std::string& get_path() const {
            return this->path_;
        }
        char get_population_name_delimiter() const {
            return this->population_name_delimiter_;
        }
        bool population_name_is_prefix() const {
            return this->population_name_is_prefix_;
        }
        bool genotypes_are_diploid() const {
            return this->genotypes_are_diploid_;
        }
        bool markers_are_dominant() const {
            return this->markers_are_dominant_;
        }
        bool constant_sites_removed() const {
            return this->constant_sites_removed_;
        }
        bool constrain_mutation_rates() const {
            return this->constrain_mutation_rates_;
        }
        bool constrain_population_sizes() const {
            return this->constrain_population_sizes_;
        }
        bool use_empirical_mutation_rate_starting_values() const {
            return this->use_empirical_mutation_rate_starting_values_;
        }
        const PositiveRealParameterSettings& get_population_size_settings() const {
            return this->population_size_settings_;
        }
        const PositiveRealParameterSettings& get_u_settings() const {
            return this->u_settings_;
        }
        const PositiveRealParameterSettings& get_time_multiplier_settings() const {
            return this->time_multiplier_settings_;
        }

        std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            ss << std::boolalpha;
            std::string margin = string_util::get_indent(indent_level);
            std::string indent = string_util::get_indent(1);
            ss << margin << "path: " << this->path_ << "\n";
            ss << margin << "genotypes_are_diploid: " << this->genotypes_are_diploid_ << "\n";
            ss << margin << "markers_are_dominant: " << this->markers_are_dominant_ << "\n";
            ss << margin << "population_name_delimiter: '" << this->population_name_delimiter_ << "'\n";
            ss << margin << "population_name_is_prefix: " << this->population_name_is_prefix_ << "\n";
            ss << margin << "constant_sites_removed: " << this->constant_sites_removed_ << "\n";
            ss << margin << "use_empirical_mutation_rate_starting_values: " << this->use_empirical_mutation_rate_starting_values_ << "\n";
            ss << margin << "constrain_population_sizes: " << this->constrain_population_sizes_ << "\n";
            ss << margin << "constrain_mutation_rates: " << this->constrain_mutation_rates_ << "\n";
            ss << margin << "parameters:\n";

            ss << margin << indent << "population_size:\n";
            ss << this->population_size_settings_.to_string(indent_level + 2);

            ss << margin << indent << "u_rate:\n";
            ss << this->u_settings_.to_string(indent_level + 2);

            ss << margin << indent << "time_multiplier:\n";
            ss << this->time_multiplier_settings_.to_string(indent_level + 2);

            return ss.str();
        }
};


class OperatorSettings {
    protected:
        double weight_;

    public:
        OperatorSettings() { }
        OperatorSettings(double weight) {
            this->weight_ = weight;
        }
        virtual ~OperatorSettings() { }
        OperatorSettings& operator=(const OperatorSettings& other) {
            this->weight_ = other.weight_;
            return * this;
        }
        double get_weight() const {
            return this->weight_;
        }
        void set_weight(double weight) {
            this->weight_ = weight;
        }
        virtual void update_from_config(const YAML::Node& parameters) {
            if (! parameters.IsMap()) {
                std::string message = (
                        "Expecting operator parameters to be a map, but found: " +
                        YamlCppUtils::get_node_type(parameters));
                throw EcoevolityYamlConfigError(message);
            }

            for (YAML::const_iterator p = parameters.begin();
                    p != parameters.end();
                    ++p) {
                if (p->first.as<std::string>() == "weight") {
                    this->set_weight(p->second.as<double>());
                }
                else {
                    std::string message = (
                            "Unrecognized key in operator parameters: " +
                            p->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

        virtual std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            std::string margin = string_util::get_indent(indent_level);
            ss << margin << "weight: " << this->weight_ << "\n";
            return ss.str();
        }
};

class ModelOperatorSettings : public OperatorSettings {
    protected:
        unsigned int number_of_auxiliary_categories_;

    public:
        ModelOperatorSettings() { }
        ModelOperatorSettings(
                double weight,
                double number_of_auxiliary_categories)
                : OperatorSettings(weight) {
            this->number_of_auxiliary_categories_ = number_of_auxiliary_categories;
        }
        virtual ~ModelOperatorSettings() { }
        ModelOperatorSettings& operator=(const ModelOperatorSettings& other) {
            this->weight_ = other.weight_;
            this->number_of_auxiliary_categories_ = other.number_of_auxiliary_categories_;
            return * this;
        }
        double get_number_of_auxiliary_categories() const {
            return this->number_of_auxiliary_categories_;
        }
        void set_number_of_auxiliary_categories(double number_of_auxiliary_categories) {
            this->number_of_auxiliary_categories_ = number_of_auxiliary_categories;
        }
        virtual void update_from_config(const YAML::Node& parameters) {
            if (! parameters.IsMap()) {
                std::string message = (
                        "Expecting operator parameters to be a map, but found: " +
                        YamlCppUtils::get_node_type(parameters));
                throw EcoevolityYamlConfigError(message);
            }

            for (YAML::const_iterator p = parameters.begin();
                    p != parameters.end();
                    ++p) {
                if (p->first.as<std::string>() == "weight") {
                    this->set_weight(p->second.as<double>());
                }
                else if (p->first.as<std::string>() == "number_of_auxiliary_categories") {
                    this->set_number_of_auxiliary_categories(p->second.as<unsigned int>());
                }
                else {
                    std::string message = (
                            "Unrecognized key in operator parameters: " +
                            p->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

        virtual std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            std::string margin = string_util::get_indent(indent_level);
            ss << margin << "weight: " << this->weight_ << "\n";
            ss << margin << "number_of_auxiliary_categories: " << this->number_of_auxiliary_categories_ << "\n";
            return ss.str();
        }
};

class ScaleOperatorSettings : public OperatorSettings {
    protected:
        double scale_;

    public:
        ScaleOperatorSettings() { }
        ScaleOperatorSettings(double weight, double scale) : OperatorSettings(weight) {
            this->scale_ = scale;
        }
        virtual ~ScaleOperatorSettings() { }
        ScaleOperatorSettings& operator=(const ScaleOperatorSettings& other) {
            this->weight_ = other.weight_;
            this->scale_ = other.scale_;
            return * this;
        }
        double get_scale() const {
            return this->scale_;
        }
        void set_scale(double scale) {
            this->scale_ = scale;
        }
        virtual void update_from_config(const YAML::Node& parameters) {
            if (! parameters.IsMap()) {
                std::string message = (
                        "Expecting operator parameters to be a map, but found: " +
                        YamlCppUtils::get_node_type(parameters));
                throw EcoevolityYamlConfigError(message);
            }

            for (YAML::const_iterator p = parameters.begin();
                    p != parameters.end();
                    ++p) {
                if (p->first.as<std::string>() == "weight") {
                    this->set_weight(p->second.as<double>());
                }
                else if (p->first.as<std::string>() == "scale") {
                    this->set_scale(p->second.as<double>());
                }
                else {
                    std::string message = (
                            "Unrecognized key in operator parameters: " +
                            p->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

        virtual std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            std::string margin = string_util::get_indent(indent_level);
            ss << margin << "weight: " << this->weight_ << "\n";
            ss << margin << "scale: " << this->scale_ << "\n";
            return ss.str();
        }
};

class WindowOperatorSettings : public OperatorSettings {
    protected:
        double window_;

    public:
        WindowOperatorSettings() { }
        WindowOperatorSettings(double weight, double window) : OperatorSettings(weight) {
            this->window_ = window;
        }
        virtual ~WindowOperatorSettings() { }
        WindowOperatorSettings& operator=(const WindowOperatorSettings& other) {
            this->weight_ = other.weight_;
            this->window_ = other.window_;
            return * this;
        }
        double get_window() const {
            return this->window_;
        }
        void set_window(double window) {
            this->window_ = window;
        }
        virtual void update_from_config(const YAML::Node& parameters) {
            if (! parameters.IsMap()) {
                std::string message = (
                        "Expecting operator parameters to be a map, but found: " +
                        YamlCppUtils::get_node_type(parameters));
                throw EcoevolityYamlConfigError(message);
            }

            for (YAML::const_iterator p = parameters.begin();
                    p != parameters.end();
                    ++p) {
                if (p->first.as<std::string>() == "weight") {
                    this->set_weight(p->second.as<double>());
                }
                else if (p->first.as<std::string>() == "window") {
                    this->set_window(p->second.as<double>());
                }
                else {
                    std::string message = (
                            "Unrecognized key in operator parameters: " +
                            p->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

        virtual std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            std::string margin = string_util::get_indent(indent_level);
            ss << margin << "weight: " << this->weight_ << "\n";
            ss << margin << "window: " << this->window_ << "\n";
            return ss.str();
        }
};

class OperatorScheduleSettings {

    friend class CollectionSettings;

    private:
        bool auto_optimize_ = true;
        unsigned int auto_optimize_delay_ = 10000;
        ModelOperatorSettings model_operator_settings_ = ModelOperatorSettings(
                3.0, 4);
        ScaleOperatorSettings concentration_scaler_settings_ = ScaleOperatorSettings(
                1.0, 0.5);
        ScaleOperatorSettings comparison_height_scaler_settings_ = ScaleOperatorSettings(
                1.0, 0.5);
        ScaleOperatorSettings comparison_height_multiplier_scaler_settings_ = ScaleOperatorSettings(
                1.0, 0.3);
        ScaleOperatorSettings root_population_size_scaler_settings_ = ScaleOperatorSettings(
                1.0, 0.5);
        ScaleOperatorSettings child_population_size_scaler_settings_ = ScaleOperatorSettings(
                1.0, 0.5);
        ScaleOperatorSettings mutation_rate_scaler_settings_ = ScaleOperatorSettings(
                1.0, 0.5);

    public:
        OperatorScheduleSettings() { }
        virtual ~OperatorScheduleSettings() { }
        OperatorScheduleSettings& operator=(const OperatorScheduleSettings& other) {
            this->auto_optimize_ = other.auto_optimize_;
            this->auto_optimize_delay_ = other.auto_optimize_delay_;
            this->model_operator_settings_ = other.model_operator_settings_;
            this->concentration_scaler_settings_ = other.concentration_scaler_settings_;
            this->comparison_height_scaler_settings_ = other.comparison_height_scaler_settings_;
            this->comparison_height_multiplier_scaler_settings_ = other.comparison_height_multiplier_scaler_settings_;
            this->root_population_size_scaler_settings_ = other.root_population_size_scaler_settings_;
            this->child_population_size_scaler_settings_ = other.child_population_size_scaler_settings_;
            this->mutation_rate_scaler_settings_ = other.mutation_rate_scaler_settings_;
            return * this;
        }

        bool auto_optimizing() const {
            return this->auto_optimize_;
        }
        unsigned int get_auto_optimize_delay() const {
            return this->auto_optimize_delay_;
        }
        const OperatorSettings& get_model_operator_settings() const {
            return this->model_operator_settings_;
        }
        const ScaleOperatorSettings& get_concentration_scaler_settings() const {
            return this->concentration_scaler_settings_;
        }
        const ScaleOperatorSettings& get_comparison_height_scaler_settings() const {
            return this->comparison_height_scaler_settings_;
        }
        const ScaleOperatorSettings& get_comparison_height_multiplier_scaler_settings() const {
            return this->comparison_height_multiplier_scaler_settings_;
        }
        const ScaleOperatorSettings& get_root_population_size_scaler_settings() const {
            return this->root_population_size_scaler_settings_;
        }
        const ScaleOperatorSettings& get_child_population_size_scaler_settings() const {
            return this->child_population_size_scaler_settings_;
        }
        const ScaleOperatorSettings& get_mutation_rate_scaler_settings() const {
            return this->mutation_rate_scaler_settings_;
        }

        void update_from_config(const YAML::Node& operator_node) {
            if (! operator_node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "Expecting operator_settings node to be a map, but found: " +
                        YamlCppUtils::get_node_type(operator_node));
            }

            for (YAML::const_iterator setting = operator_node.begin();
                    setting != operator_node.end();
                    ++setting) {
                if (setting->first.as<std::string>() == "auto_optimize") {
                    this->auto_optimize_ = setting->second.as<bool>();
                }
                else if (setting->first.as<std::string>() == "auto_optimize_delay") {
                    this->auto_optimize_delay_ = setting->second.as<unsigned int>();
                }
                else if (setting->first.as<std::string>() == "operators") {
                    this->parse_operators(setting->second);
                }
                else {
                    std::string message = (
                            "Unrecognized key in operator_settings: " +
                            setting->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }
        
        void parse_operators(const YAML::Node& operators) {
            if (! operators.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "Expecting operators node to be a map, but found: " +
                        YamlCppUtils::get_node_type(operators));
            }

            for (YAML::const_iterator op = operators.begin();
                    op != operators.end();
                    ++op) {
                if (op->first.as<std::string>() == "ModelOperator") {
                    try {
                        this->model_operator_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing ModelOperator settings\n";
                        throw;
                    }
                }
                else if (op->first.as<std::string>() == "ConcentrationScaler") {
                    try {
                        this->concentration_scaler_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing ConcentrationScaler settings\n";
                        throw;
                    }
                }
                else if (op->first.as<std::string>() == "ComparisonHeightScaler") {
                    try {
                        this->comparison_height_scaler_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing ComparisonHeightScaler settings\n";
                        throw;
                    }
                }
                else if (op->first.as<std::string>() == "ComparisonHeightMultiplierScaler") {
                    try {
                        this->comparison_height_multiplier_scaler_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing ComparisonHeightMultiplierScaler settings\n";
                        throw;
                    }
                }
                else if (op->first.as<std::string>() == "RootPopulationSizeScaler") {
                    try {
                        this->root_population_size_scaler_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing RootPopulationSizeScaler settings\n";
                        throw;
                    }
                }
                else if (op->first.as<std::string>() == "ChildPopulationSizeScaler") {
                    try {
                        this->child_population_size_scaler_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing ChildPopulationSizeScaler settings\n";
                        throw;
                    }
                }
                else if (op->first.as<std::string>() == "MutationRateScaler") {
                    try {
                        this->mutation_rate_scaler_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing MutationRateScaler settings\n";
                        throw;
                    }
                }
                else {
                    std::string message = (
                            "Unrecognized operator: " +
                            op->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

        virtual std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            ss << std::boolalpha;
            std::string margin = string_util::get_indent(indent_level);
            std::string indent = string_util::get_indent(1);
            ss << margin << "operator_settings:\n";
            ss << margin << indent << "auto_optimize: " << this->auto_optimize_ << "\n";
            ss << margin << indent << "auto_optimize_delay: " << this->auto_optimize_delay_ << "\n";
            ss << margin << indent << "operators:\n";
            ss << margin << indent << indent << "ModelOperator:\n";
            ss << this->model_operator_settings_.to_string(indent_level + 3);
            ss << margin << indent << indent << "ConcentrationScaler:\n";
            ss << this->concentration_scaler_settings_.to_string(indent_level + 3);
            ss << margin << indent << indent << "ComparisonHeightScaler:\n";
            ss << this->comparison_height_scaler_settings_.to_string(indent_level + 3);
            ss << margin << indent << indent << "ComparisonHeightMultiplierScaler:\n";
            ss << this->comparison_height_multiplier_scaler_settings_.to_string(indent_level + 3);
            ss << margin << indent << indent << "RootPopulationSizeScaler:\n";
            ss << this->root_population_size_scaler_settings_.to_string(indent_level + 3);
            ss << margin << indent << indent << "ChildPopulationSizeScaler:\n";
            ss << this->child_population_size_scaler_settings_.to_string(indent_level + 3);
            ss << margin << indent << indent << "MutationRateScaler:\n";
            ss << this->mutation_rate_scaler_settings_.to_string(indent_level + 3);
            return ss.str();
        }
};


class CollectionSettings {

    public:

        CollectionSettings() {
            this->init_default_priors();
        }
        CollectionSettings(
                const ContinuousDistributionSettings& time_prior,
                unsigned int chain_length,
                unsigned int sample_frequency,
                const PositiveRealParameterSettings& concentration_settings,
                bool use_dpp)
                : CollectionSettings() {
            this->time_prior_settings_ = time_prior;
            this->chain_length_ = chain_length;
            this->sample_frequency_ = sample_frequency;
            this->concentration_settings_ = concentration_settings;
            this->use_dpp_ = use_dpp;
        }
        CollectionSettings(
                const std::string & yaml_config_path)
                : CollectionSettings() {
            this->init_from_config_file(yaml_config_path);
        }
        CollectionSettings(
                std::istream& yaml_config_stream,
                const std::string& yaml_config_path)
                : CollectionSettings() {
            this->init_from_config_stream(yaml_config_stream, yaml_config_path);
        }
        virtual ~CollectionSettings() { }
        CollectionSettings& operator=(const CollectionSettings& other) {
            this->path_ = other.path_;
            this->operator_schedule_settings_ = other.operator_schedule_settings_;
            this->global_comparison_settings_ = other.global_comparison_settings_;
            this->time_prior_settings_ = other.time_prior_settings_;
            this->chain_length_ = other.chain_length_;
            this->sample_frequency_ = other.sample_frequency_;
            this->concentration_settings_ = other.concentration_settings_;
            this->use_dpp_ = other.use_dpp_;
            this->comparisons_ = other.comparisons_;
            this->default_time_prior_ = other.default_time_prior_;
            this->default_population_size_prior_ = other.default_population_size_prior_;
            this->default_u_prior_ = other.default_u_prior_;
            this->default_time_multiplier_prior_ = other.default_time_multiplier_prior_;
            this->state_log_path_ = other.state_log_path_;
            this->operator_log_path_ = other.operator_log_path_;
            return * this;
        }

        void add_comparison(ComparisonSettings comparison_settings) {
            this->comparisons_.push_back(comparison_settings);
        }

        const std::string& get_path() const {
            return this->path_;
        }
        std::string get_config_directory() const {
            return path::dirname(this->get_path());
        }
        const std::string& get_state_log_path() const {
            return this->state_log_path_;
        }
        const std::string& get_operator_log_path() const {
            return this->operator_log_path_;
        }
        bool using_dpp() const {
            return this->use_dpp_;
        }
        double get_chain_length() const {
            return this->chain_length_;
        }
        double get_sample_frequency() const {
            return this->sample_frequency_;
        }
        unsigned int get_number_of_comparisons() const {
            return this->comparisons_.size();
        }

        unsigned int get_number_of_comparisons_with_free_time_multiplier() const {
            unsigned int nfree = 0;
            for (const ComparisonSettings& comparison : this->comparisons_) {
                if (! comparison.time_multiplier_settings_.is_fixed()) {
                    ++nfree;
                }
            }
            return nfree;
        }

        unsigned int get_number_of_comparisons_with_free_population_size() const {
            unsigned int nfree = 0;
            for (const ComparisonSettings& comparison : this->comparisons_) {
                if (! comparison.population_size_settings_.is_fixed()) {
                    ++nfree;
                }
            }
            return nfree;
        }

        unsigned int get_number_of_comparisons_with_free_u_rate() const {
            unsigned int nfree = 0;
            for (const ComparisonSettings& comparison : this->comparisons_) {
                if (! comparison.u_settings_.is_fixed()) {
                    ++nfree;
                }
            }
            return nfree;
        }

        const ContinuousDistributionSettings& get_time_prior_settings() const {
            return this->time_prior_settings_;
        }

        const PositiveRealParameterSettings& get_concentration_settings() const {
            return this->concentration_settings_;
        }

        const std::vector<ComparisonSettings>& get_comparison_settings() const { 
            return this->comparisons_;
        }

        const OperatorScheduleSettings& get_operator_schedule_settings() const {
            return this->operator_schedule_settings_;
        }

        std::string to_string() const {
            std::ostringstream ss;
            ss << std::boolalpha;
            std::string indent = string_util::get_indent(1);

            ss << "---\n"
               << "event_model_prior:\n";
            if (this->use_dpp_) {
                ss << indent << "dirichlet_process:\n"
                   << indent << indent << "parameters:\n"
                   << indent << indent <<  indent <<"concentration:\n";
                ss << this->concentration_settings_.to_string(4);
            }
            else {
                ss << indent << "uniform:\n";
            }

            ss << "event_time_prior:\n";
            ss << this->time_prior_settings_.to_string(1);

            ss << "mcmc_settings:\n"
               << indent << "chain_length: " << this->chain_length_ << "\n"
               << indent << "sample_frequency: " << this->sample_frequency_ << "\n";

            ss << "output_settings:\n"
               << indent << "state_log_path: " << this->state_log_path_ << "\n"
               << indent << "operator_log_path: " << this->operator_log_path_ << "\n";

            ss << "comparisons:\n";
            for (auto comp : this->comparisons_) {
                ss << "- comparison:\n";
                ss << comp.to_string(1);
            }

            ss << this->operator_schedule_settings_.to_string();
            return ss.str();
        }


    private:

        std::string path_ = "";
        bool use_dpp_ = true;
        unsigned int chain_length_ = 100000;
        unsigned int sample_frequency_ = 100;
        std::string state_log_path_ = "ecoevolity-state.log";
        std::string operator_log_path_ = "ecoevolity-operator.log";

        OperatorScheduleSettings operator_schedule_settings_;

        ContinuousDistributionSettings time_prior_settings_;

        PositiveRealParameterSettings concentration_settings_;

        ComparisonSettings global_comparison_settings_;

        std::vector<ComparisonSettings> comparisons_;

        ContinuousDistributionSettings default_time_prior_;
        ContinuousDistributionSettings default_population_size_prior_;
        ContinuousDistributionSettings default_u_prior_;
        ContinuousDistributionSettings default_time_multiplier_prior_;

        void init_default_priors() {
            std::unordered_map<std::string, double> default_parameters;
            default_parameters["rate"] = 100.0;
            this->default_time_prior_ = ContinuousDistributionSettings(
                    "exponential_distribution",
                    default_parameters);
            default_parameters.clear();
            default_parameters["rate"] = 1000.0;
            this->default_population_size_prior_ = ContinuousDistributionSettings(
                    "exponential_distribution",
                    default_parameters);
            default_parameters.clear();
            default_parameters["rate"] = 2.0;
            default_parameters["offset"] = 0.5;
            this->default_u_prior_ = ContinuousDistributionSettings(
                    "exponential_distribution",
                    default_parameters);
            default_parameters.clear();
            default_parameters["shape"] = 1000.0;
            default_parameters["scale"] = 0.001;
            this->default_time_multiplier_prior_ = ContinuousDistributionSettings(
                    "gamma_distribution",
                    default_parameters);
        }

        void set_output_paths_to_config_directory() {
            this->state_log_path_ = path::join(
                    path::dirname(this->path_),
                    path::basename(this->state_log_path_));
            this->operator_log_path_ = path::join(
                    path::dirname(this->path_),
                    path::basename(this->operator_log_path_));
        }

        void init_from_config_stream(std::istream& stream, const std::string& path) {
            this->path_ = path;
            this->set_output_paths_to_config_directory();
            this->parse_yaml_config(stream);
        }
        void init_from_config_file(const std::string& path) {
            std::ifstream in_stream;
            in_stream.open(path);
            if (! in_stream.is_open()) {
                throw EcoevolityYamlConfigError(
                        "Could not open YAML config file",
                        path);
            }
            this->init_from_config_stream(in_stream, path);
            in_stream.close();
        }

        void parse_yaml_config(std::istream& config_stream) {
            YAML::Node config;
            try {
                config = YAML::Load(config_stream);
            }
            catch (...) {
                std::cerr << "ERROR: Problem with YAML-formatting of config\n";
                throw;
            }
            this->parse_top_level(config);
        }

        void parse_top_level(const YAML::Node& top_level_node) {
            if (! top_level_node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "Expecting top-level of config to be a map, but found: " +
                        YamlCppUtils::get_node_type(top_level_node));
            }

            ///////////////////////////////////////////////////////////////////
            // Set up defaults to be overridden by config settings:
            // Set default time prior
            this->time_prior_settings_ = this->default_time_prior_;
            // Set default time multiplier
            this->global_comparison_settings_.time_multiplier_settings_.value_ = 1.0;
            this->global_comparison_settings_.time_multiplier_settings_.is_fixed_ = true;
            ///////////////////////////////////////////////////////////////////

            if (! top_level_node["comparisons"]) {
                throw EcoevolityYamlConfigError("No comparisons");
            }
            // Need to parse comparisons first, because other setting rely on
            // the knowing the number of comparisons.
            if (top_level_node["global_comparison_settings"]) {
                this->parse_global_comparison_settings(
                        top_level_node["global_comparison_settings"]);
            }
            this->parse_comparisons(top_level_node["comparisons"]);

            ///////////////////////////////////////////////////////////////////
            // More setting up defaults
            // Update unset priors for comparisons
            this->update_default_comparison_priors();
            // Set default concentration prior to be overriden by config
            double default_shape = 2.0;
            double default_prior_mean_num_events = 1.0;
            double default_scale = 1.0;
            if (this->comparisons_.size() > 1) {
                default_prior_mean_num_events = (double)this->comparisons_.size() * 0.75;
                default_scale = get_dpp_gamma_scale(
                        default_prior_mean_num_events,
                        this->comparisons_.size(),
                        default_shape);
            }
            else {
                this->use_dpp_ = false;
            }
            std::unordered_map<std::string, double> default_parameters;
            default_parameters["shape"] = default_shape;
            default_parameters["scale"] = default_scale;
            this->concentration_settings_.prior_settings_ = ContinuousDistributionSettings(
                    "gamma_distribution",
                    default_parameters);
            ///////////////////////////////////////////////////////////////////

            for (YAML::const_iterator top = top_level_node.begin();
                    top != top_level_node.end();
                    ++top) {
                // parse mcmc settings
                if (top->first.as<std::string>() == "mcmc_settings") {
                    this->parse_mcmc_settings(top->second);
                }
                // parse output settings
                else if (top->first.as<std::string>() == "output_settings") {
                    this->parse_output_settings(top->second);
                }
                // parse operator settings
                else if (top->first.as<std::string>() == "operator_settings") {
                    this->operator_schedule_settings_.update_from_config(top->second);
                }
                // parse model prior
                else if (top->first.as<std::string>() == "event_model_prior") {
                    this->parse_model_prior(top->second);
                }
                //parse time prior
                else if (top->first.as<std::string>() == "event_time_prior") {
                    this->parse_time_prior(top->second);
                }
                // parse comparison defaults
                else if (top->first.as<std::string>() == "global_comparison_settings") {
                    // Already parsed, nothing to do here
                    continue;
                }
                // parse comparisons
                else if (top->first.as<std::string>() == "comparisons") {
                    // Already parsed, nothing to do here
                    continue;
                }
                else {
                    throw EcoevolityYamlConfigError(
                            "Unrecognized top level key: '" +
                            top->first.as<std::string>() + "'");
                }
            }

            // Override DPP settings if only one comparison
            if (this->comparisons_.size() < 2) {
                this->use_dpp_ = false;
            }

            // "Turn off" operators that are not needed
            this->update_operator_schedule_settings();
            
        }

        void update_operator_schedule_settings() {
            if ((! this->use_dpp_) || (this->concentration_settings_.is_fixed())) {
                this->operator_schedule_settings_.concentration_scaler_settings_.set_weight(0.0);
            }
            if (this->comparisons_.size() < 2) {
                this->operator_schedule_settings_.model_operator_settings_.set_weight(0.0);
            }
            if (this->get_number_of_comparisons_with_free_time_multiplier() < 1) {
                this->operator_schedule_settings_.comparison_height_multiplier_scaler_settings_.set_weight(0.0);
            }
            if (this->get_number_of_comparisons_with_free_u_rate() < 1) {
                this->operator_schedule_settings_.mutation_rate_scaler_settings_.set_weight(0.0);
            }
            if (this->get_number_of_comparisons_with_free_population_size() < 1) {
                this->operator_schedule_settings_.root_population_size_scaler_settings_.set_weight(0.0);
                this->operator_schedule_settings_.child_population_size_scaler_settings_.set_weight(0.0);
            }
        }

        void update_default_comparison_priors() {
            for (auto&& comp : this->comparisons_) {
                if ((! comp.time_multiplier_settings_.is_fixed()) &&
                        (comp.time_multiplier_settings_.prior_settings_.get_name() == "none")) {
                    comp.time_multiplier_settings_.prior_settings_ = this->default_time_multiplier_prior_;
                }
                if ((! comp.u_settings_.is_fixed()) &&
                        (comp.u_settings_.prior_settings_.get_name() == "none")) {
                    comp.u_settings_.prior_settings_ = this->default_u_prior_;
                }
                if ((! comp.population_size_settings_.is_fixed()) &&
                        (comp.population_size_settings_.prior_settings_.get_name() == "none")) {
                    comp.population_size_settings_.prior_settings_ = this->default_population_size_prior_;
                }
            }
        }

        void parse_mcmc_settings(const YAML::Node& mcmc_node) {
            if (! mcmc_node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "Expecting mcmc_settings to be a map, but found: " +
                        YamlCppUtils::get_node_type(mcmc_node));
            }
            for (YAML::const_iterator mcmc = mcmc_node.begin(); mcmc != mcmc_node.end(); ++mcmc) {
                if (mcmc->first.as<std::string>() == "chain_length") {
                    this->chain_length_ = mcmc->second.as<unsigned int>();

                }
                else if (mcmc->first.as<std::string>() == "sample_frequency") {
                    this->sample_frequency_ = mcmc->second.as<unsigned int>();
                }
                else {
                    throw EcoevolityYamlConfigError(
                            "Unrecognized mcmc_settings key: " +
                            mcmc->first.as<std::string>());
                }
            }
        }

        void parse_output_settings(const YAML::Node& output_node) {
            if (! output_node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "Expecting output_settings to be a map, but found: " +
                        YamlCppUtils::get_node_type(output_node));
            }
            for (YAML::const_iterator output = output_node.begin(); output != output_node.end(); ++output) {
                if (output->first.as<std::string>() == "state_log_path") {
                    this->state_log_path_ = path::join(
                            path::dirname(this->path_),
                            output->second.as<std::string>());

                }
                else if (output->first.as<std::string>() == "operator_log_path") {
                    this->operator_log_path_ = path::join(
                            path::dirname(this->path_),
                            output->second.as<std::string>());
                }
                else {
                    throw EcoevolityYamlConfigError(
                            "Unrecognized output_settings key: " +
                            output->first.as<std::string>());
                }
            }
        }

        void parse_model_prior(const YAML::Node& model_prior_node) {
            if (! model_prior_node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "Expecting event_model_prior to be a map, but found: " +
                        YamlCppUtils::get_node_type(model_prior_node));
            }

            if (model_prior_node.size() > 1) {
                throw EcoevolityYamlConfigError(
                        "event_model_prior node should only have a single key");
            }
            if (model_prior_node["dirichlet_process"]) {
                this->use_dpp_ = true;
                this->parse_dirichlet_process_prior(model_prior_node["dirichlet_process"]);
            }
            else if (model_prior_node["uniform"]) {
                this->use_dpp_ = false;
            }
            else {
                throw EcoevolityYamlConfigError(
                        "Unrecognized event_model_prior key: " +
                        model_prior_node.begin()->first.as<std::string>());
            }
        }

        void parse_time_prior(const YAML::Node& time_prior_node) {
            this->time_prior_settings_ = ContinuousDistributionSettings(time_prior_node);
        }

        void parse_global_comparison_settings(const YAML::Node& global_node) {
            this->global_comparison_settings_.update_from_config(global_node,
                    this->path_, true);
        }

        void parse_comparisons(const YAML::Node& comparisons_node) {
            if (! comparisons_node.IsSequence()) {
                throw EcoevolityYamlConfigError(
                        "Expecting comparisons to be a sequence, but found: " +
                        YamlCppUtils::get_node_type(comparisons_node));
            }
            if (comparisons_node.size() < 1) {
                throw EcoevolityYamlConfigError("No comparsions found");
            }
            for (unsigned int comp_idx = 0;
                    comp_idx < comparisons_node.size();
                    ++comp_idx) {
                YAML::Node comp = comparisons_node[comp_idx];
                if (! comp.IsMap()) {
                    throw EcoevolityYamlConfigError(
                            "Expecting each comparison to be a map, but found: " +
                            YamlCppUtils::get_node_type(comp));
                }
                if (comp.size() != 1) {
                    throw EcoevolityYamlConfigError(
                            "Each comparison should only have one key");
                }
                if (! comp["comparison"]) {
                    throw EcoevolityYamlConfigError(
                            "Each comparison should have a comparison key");
                }
                ComparisonSettings c = this->global_comparison_settings_;
                c.update_from_config(comp["comparison"], this->path_);
                this->comparisons_.push_back(c);
            }
            ECOEVOLITY_ASSERT(comparisons_node.size() == this->comparisons_.size());
        }

        void parse_dirichlet_process_prior(const YAML::Node& dpp_node) {
            if (! dpp_node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "Expecting event_model_prior to be a map, but found: " +
                        YamlCppUtils::get_node_type(dpp_node));
            }

            if (dpp_node.size() != 1) {
                throw EcoevolityYamlConfigError(
                        "dirichlet_process node must have a single key");
            }
            if (! dpp_node["parameters"]) {
                throw EcoevolityYamlConfigError(
                        "dirichlet_process must have a parameters key");
            }
            if (! dpp_node["parameters"]["concentration"]) {
                throw EcoevolityYamlConfigError(
                        "dirichlet_process must specify concentration parameter settings");
            }
            this->parse_concentration_parameter(dpp_node["parameters"]["concentration"]);
        }

        void parse_concentration_parameter(const YAML::Node& node) {
            if (! node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "Expecting concentration node to be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
            }
            if (node.size() < 1) {
                throw EcoevolityYamlConfigError(
                        "empty concentration parameter node");
            }

            for (YAML::const_iterator arg = node.begin();
                    arg != node.end();
                    ++arg) {
                if (arg->first.as<std::string>() == "value") {
                    double v = arg->second.as<double>();
                        if (v < 0.0) {
                            throw EcoevolityPositiveRealParameterSettingError(
                                    "concentration parameter cannot be less than 0"
                                    );
                        }
                    this->concentration_settings_.value_ = v;
                }
                else if (arg->first.as<std::string>() == "estimate") {
                    bool f = arg->second.as<bool>();
                    this->concentration_settings_.is_fixed_ = (! f);
                }
                else if (arg->first.as<std::string>() == "prior") {
                    this->concentration_settings_.prior_settings_ = this->parse_dpp_gamma_hyper_prior(
                            arg->second);
                }
                else {
                    std::string message = "Unrecognized concentration key: " +
                            arg->first.as<std::string>();
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

        ContinuousDistributionSettings parse_dpp_gamma_hyper_prior(const YAML::Node& node) {
            if (! node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "Expecting concentration prior node to be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
            }
            if (node.size() != 1) {
                throw EcoevolityYamlConfigError(
                        "concentration prior node should only have a single key");
            }
            if (! node["gamma_distribution"]) {
                std::string message = "concentration prior must be a gamma_distribution; found: " +
                        node.begin()->first.as<std::string>();
                throw EcoevolityContinuousDistributionSettingError(message);
            }
            YAML::Node parameters = node["gamma_distribution"];
            if (! parameters["shape"]) {
                throw EcoevolityContinuousDistributionSettingError(
                        "shape parameter missing for gamma prior on concentration"
                        );
            }
            if ((! parameters["scale"]) && (! parameters["prior_mean_number_of_events"])) {
                throw EcoevolityContinuousDistributionSettingError(
                        "gamma prior on concentration requires scale or prior_mean_number_of_events parameter"
                        );
            }
            if (parameters.size() > 2) {
                throw EcoevolityContinuousDistributionSettingError(
                        "unrecognized parameters for gamma prior on concentration (should be shape and scale OR shape and prior_mean_number_of_events)"
                        );
            }
            double shape = parameters["shape"].as<double>();
            double scale = -1.0;
            if (parameters["scale"]) {
                scale = parameters["scale"].as<double>();
            }
            else if (parameters["prior_mean_number_of_events"]) {
                double prior_mean_num_events = parameters["prior_mean_number_of_events"].as<double>();
                if (prior_mean_num_events < 1.0) {
                    throw EcoevolityYamlConfigError("prior_mean_number_of_events must be at least 1.0");
                }
                else if (prior_mean_num_events > (double) this->comparisons_.size()) {
                    throw EcoevolityYamlConfigError(
                            "prior_mean_number_of_events cannot be greater than the number of comparisons");
                }
                scale = get_dpp_gamma_scale(
                        prior_mean_num_events,
                        this->comparisons_.size(),
                        shape);
            }
            else {
                throw EcoevolityContinuousDistributionSettingError(
                        "invalid parameter for gamma prior on concentration");
            }
            if ((shape <= 0.0) || (scale <= 0.0)) {
                throw EcoevolityContinuousDistributionSettingError(
                        "Shape and scale must be greater than zero for gamma distribution");
            }
            std::unordered_map<std::string, double> p;
            p["shape"] = shape;
            p["scale"] = scale;
            return ContinuousDistributionSettings("gamma_distribution", p);
        }
};

#endif
