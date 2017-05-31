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
#include <unordered_set>

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
        std::vector<double> concentration_parameters_; // For Dirichlet dist

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
            // Beta 
            else if (name == "beta_distribution") {
                if (parameters.count("alpha") < 1) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "alpha parameter missing for beta_distribution"
                            );
                }
                if (parameters.count("beta") < 1) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "beta parameter missing for beta_distribution"
                            );
                }
                if (parameters.size() > 2) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "unrecognized parameters for beta_distribution (recognized parameters: alpha, beta)"
                            );
                }
                if (parameters.at("alpha") <= 0.0) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "alpha must be greater than 0 for Distribution");
                }
                if (parameters.at("beta") <= 0.0) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "beta must be greater than 0 for BetaDistribution");
                }
                this->parameters_["alpha"] = parameters.at("alpha");
                this->parameters_["beta"] = parameters.at("beta");
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
            // Beta 
            else if (node["beta_distribution"]) {
                this->name_ = "beta_distribution";
                YAML::Node parameters = node["beta_distribution"];
                if (! parameters["alpha"]) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "alpha parameter missing for beta_distribution"
                            );
                }
                if (! parameters["beta"]) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "beta parameter missing for beta_distribution"
                            );
                }
                if (parameters.size() > 2) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "unrecognized parameters for beta_distribution (recognized parameters: alpha, beta)"
                            );
                }
                this->parameters_["alpha"] = parameters["alpha"].as<double>();
                this->parameters_["beta"] = parameters["beta"].as<double>();
                if (this->parameters_.at("alpha") <= 0.0) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "alpha must be greater than 0 for BetaDistribution");
                }
                if (this->parameters_.at("beta") <= 0.0) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "beta must be greater than 0 for BetaDistribution");
                }
            }
            // Dirichlet 
            else if (node["dirichlet_distribution"]) {
                this->name_ = "dirichlet_distribution";
                YAML::Node parameters = node["dirichlet_distribution"];
                if (! parameters["alpha"]) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "alpha parameter missing for dirichlet_distribution"
                            );
                }
                if (parameters.size() > 1) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "unrecognized parameters for dirichlet_distribution (recognized parameters: alpha)"
                            );
                }
                YAML::Node alphas;
                alphas = parameters["alpha"];
                for (YAML::const_iterator a = alphas.begin(); a != alphas.end(); ++a) {
                    double conc_parameter = a->as<double>();
                    if (conc_parameter <= 0.0) {
                        throw EcoevolityContinuousDistributionSettingError(
                                "alphas must be greater than 0 for DirichletDistribution");
                    }
                    this->concentration_parameters_.push_back(conc_parameter);
                }
                if (this->concentration_parameters_.size() < 2) {
                    throw EcoevolityContinuousDistributionSettingError(
                            "at least two alphas must be specified for dirichlet_distribution"
                            );
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
            this->concentration_parameters_ = other.concentration_parameters_;
            return * this;
        }

        const std::string& get_name() const {
            return this->name_;
        }

        void nullify() {
            this->name_ = "none";
            this->parameters_.clear();
        }

        const std::vector<double> & get_concentration_parameters() const {
            return this->concentration_parameters_;
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
            else if (this->name_ == "beta_distribution") {
                p = std::make_shared<BetaDistribution>(
                        this->parameters_.at("alpha"),
                        this->parameters_.at("beta"));
            }
            else if (this->name_ == "dirichlet_distribution") {
                throw EcoevolityError("Call get_dirichlet_distribution_instance");
            }
            else if (this->name_ == "none") {
                return p;
            }
            else {
                ECOEVOLITY_ASSERT(0 == 1);
            }
            return p;
        }

        std::shared_ptr<DirichletDistribution> get_dirichlet_distribution_instance() const {
            std::shared_ptr<DirichletDistribution> p;
            if ((this->name_ != "dirichlet_distribution") && (this->name_ != "none")) {
                throw EcoevolityError("Call get_instance");
            }
            if (this->name_ == "dirichlet_distribution") {
                p = std::make_shared<DirichletDistribution>(
                        this->concentration_parameters_);
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
            else if (this->name_ == "beta_distribution") {
                ss << margin << indent << "alpha: " << this->parameters_.at("alpha") << "\n";
                ss << margin << indent << "beta: " << this->parameters_.at("beta") << "\n";
            }
            else if (this->name_ == "dirichlet_distribution") {
                ss << margin << indent << "alpha: [" << this->concentration_parameters_.at(0);
                for (unsigned int i = 1; i < this->concentration_parameters_.size(); ++i) {
                    ss << ", " << this->concentration_parameters_.at(i);
                }
                ss << "]\n";
            }
            else {
                ECOEVOLITY_ASSERT(0 == 1);
            }
            return ss.str();
        }
};

class PositiveRealParameterSettings {

    template<typename T> friend class BaseComparisonSettings;
    friend class DirichletComparisonSettings;
    friend class RelativeRootComparisonSettings;
    template<typename T> friend class BaseCollectionSettings;

    private:
        double value_ = std::numeric_limits<double>::quiet_NaN();
        std::vector<double> values_; // For Dirichlet distribution
        bool is_fixed_ = false;
        bool is_vector_ = false;
        bool use_empirical_value_ = false;
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
            bool prior_specified = false;
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
                    if (arg->second.IsSequence()) {
                        std::vector<double> temp_values;
                        double sum = 0.0;
                        YAML::Node value_node = arg->second;
                        for (YAML::const_iterator v = value_node.begin();
                                v != value_node.end();
                                ++v) {
                            double val = v->as<double>();
                            if (val <= 0.0) {
                                throw EcoevolityPositiveRealParameterSettingError(
                                        "dirichlet distribution values must be greater than 0"
                                        );
                            }
                            temp_values.push_back(val);
                            sum += val;
                        }
                        if (temp_values.size() < 2) {
                            throw EcoevolityPositiveRealParameterSettingError(
                                    "parameter vector must have at least 2 values"
                                    );
                        }
                        for (auto x : temp_values) {
                            this->values_.push_back(temp_values.size() * (x / sum));
                        }
                        this->is_vector_ = true;
                    }
                    else if (arg->second.as<std::string>() == "empirical") {
                        this->use_empirical_value_ = true;
                    }
                    else {
                        double v = arg->second.as<double>();
                            if (v < 0.0) {
                                throw EcoevolityPositiveRealParameterSettingError(
                                        "positive real parameter cannot be less than 0"
                                        );
                            }
                        this->value_ = v;
                    }
                }
                else if (arg->first.as<std::string>() == "estimate") {
                    bool f = arg->second.as<bool>();
                    this->is_fixed_ = (! f);
                }
                else if (arg->first.as<std::string>() == "prior") {
                    prior_specified = true;
                    this->prior_settings_ = ContinuousDistributionSettings(
                            arg->second);
                }
                else {
                    std::string message = "unrecognized parameter key: " +
                            arg->first.as<std::string>();
                    throw EcoevolityContinuousDistributionSettingError(message);
                }
            }

            if ((this->is_vector_) && (prior_specified)) {
                if (this->prior_settings_.get_name() != "dirichlet_distribution") {
                    throw EcoevolityPositiveRealParameterSettingError(
                            "multiple values requires dirichlet distribution prior"
                            );
                }
                std::vector<double> parameters = this->prior_settings_.get_concentration_parameters();
                if (parameters.size() != this->values_.size()) {
                    throw EcoevolityPositiveRealParameterSettingError(
                            "number of values must match number of dirichlet distribution parameters"
                            );
                }
            }

            if (this->is_fixed_) {
                std::unordered_map<std::string, double> prior_parameters;
                this->prior_settings_ = ContinuousDistributionSettings("none",
                        prior_parameters);
            }
            if ((this->is_fixed_) &&
                    (std::isnan(this->value_)) &&
                    (! this->is_vector_) &&
                    (! this->use_empirical_value_)) {
                throw EcoevolityPositiveRealParameterSettingError(
                        "cannot fix parameter without a value"
                        );
            }
        }
        virtual ~PositiveRealParameterSettings() { }
        PositiveRealParameterSettings& operator=(const PositiveRealParameterSettings& other) {
            this->value_ = other.value_;
            this->values_ = other.values_;
            this->is_fixed_ = other.is_fixed_;
            this->is_vector_ =  other.is_vector_;
            this->use_empirical_value_ = other.use_empirical_value_;
            this->prior_settings_ = other.prior_settings_;
            return * this;
        }

        double get_value() const {
            return this->value_;
        }
        const std::vector<double> & get_values() const {
            return this->values_;
        }
        bool is_fixed() const {
            return this->is_fixed_;
        }
        const ContinuousDistributionSettings& get_prior_settings() const {
            return this->prior_settings_;
        }

        bool use_empirical_value() const {
            return this->use_empirical_value_;
        }

        virtual std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            ss << std::boolalpha;
            std::string margin = string_util::get_indent(indent_level);
            if (this->is_vector_) {
                if (this->values_.size() > 0) {
                    ss << margin << "value: [" << this->values_.at(0);
                    for (unsigned int i = 1; i < this->values_.size(); ++i) {
                        ss << ", " << this->values_.at(i);
                    }
                    ss << "]\n";
                }
            }
            else if (! std::isnan(this->get_value())) {
                ss << margin << "value: " << this->get_value() << "\n";
            }
            else if (this->use_empirical_value_) {
                ss << margin << "value: empirical\n";
            }
            ss << margin << "estimate: " << (! this->is_fixed()) << "\n";
            if (! this->is_fixed()) {
                ss << margin << "prior:\n";
                ss << this->prior_settings_.to_string(indent_level + 1);
            }
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

    template<typename T1> friend class BaseCollectionSettings;

    private:
        bool auto_optimize_ = true;
        unsigned int auto_optimize_delay_ = 10000;
        bool using_population_size_multipliers_ = false;
        ModelOperatorSettings model_operator_settings_ = ModelOperatorSettings(
                3.0, 4);
        ScaleOperatorSettings concentration_scaler_settings_ = ScaleOperatorSettings(
                1.0, 0.5);
        ScaleOperatorSettings composite_time_size_rate_mixer_settings_ = ScaleOperatorSettings(
                1.0, 0.5);
        ScaleOperatorSettings composite_time_size_rate_scaler_settings_ = ScaleOperatorSettings(
                0.0, 0.5);
        ScaleOperatorSettings event_time_scaler_settings_ = ScaleOperatorSettings(
                1.0, 0.5);

    public:
        OperatorScheduleSettings() { }
        virtual ~OperatorScheduleSettings() { }
        OperatorScheduleSettings& operator=(const OperatorScheduleSettings& other) {
            this->auto_optimize_ = other.auto_optimize_;
            this->auto_optimize_delay_ = other.auto_optimize_delay_;
            this->model_operator_settings_ = other.model_operator_settings_;
            this->concentration_scaler_settings_ = other.concentration_scaler_settings_;
            this->composite_time_size_rate_mixer_settings_ = other.composite_time_size_rate_mixer_settings_;
            this->composite_time_size_rate_scaler_settings_ = other.composite_time_size_rate_scaler_settings_;
            this->event_time_scaler_settings_ = other.event_time_scaler_settings_;
            return * this;
        }

        bool using_population_size_multipliers() const {
            return this->using_population_size_multipliers_;
        }

        void turn_on_population_size_multipliers() {
            this->using_population_size_multipliers_ = true;
        }

        bool auto_optimizing() const {
            return this->auto_optimize_;
        }
        unsigned int get_auto_optimize_delay() const {
            return this->auto_optimize_delay_;
        }
        const ModelOperatorSettings& get_model_operator_settings() const {
            return this->model_operator_settings_;
        }
        const ScaleOperatorSettings& get_concentration_scaler_settings() const {
            return this->concentration_scaler_settings_;
        }
        const ScaleOperatorSettings& get_composite_time_size_rate_mixer_settings() const {
            return this->composite_time_size_rate_mixer_settings_;
        }
        const ScaleOperatorSettings& get_composite_time_size_rate_scaler_settings() const {
            return this->composite_time_size_rate_scaler_settings_;
        }
        const ScaleOperatorSettings& get_event_time_scaler_settings() const {
            return this->event_time_scaler_settings_;
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
                else if (op->first.as<std::string>() == "CompositeTimeSizeRateMixer") {
                    try {
                        this->composite_time_size_rate_mixer_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing CompositeTimeSizeRateMixer settings\n";
                        throw;
                    }
                }
                else if (op->first.as<std::string>() == "CompositeTimeSizeRateScaler") {
                    try {
                        this->composite_time_size_rate_scaler_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing CompositeTimeSizeRateScaler settings\n";
                        throw;
                    }
                }
                else if (op->first.as<std::string>() == "EventTimeScaler") {
                    try {
                        this->event_time_scaler_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing EventTimeScaler settings\n";
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
            ss << margin << indent << indent << "CompositeTimeSizeRateMixer:\n";
            ss << this->composite_time_size_rate_mixer_settings_.to_string(indent_level + 3);
            ss << margin << indent << indent << "CompositeTimeSizeRateScaler:\n";
            ss << this->composite_time_size_rate_scaler_settings_.to_string(indent_level + 3);
            ss << margin << indent << indent << "EventTimeScaler:\n";
            ss << this->event_time_scaler_settings_.to_string(indent_level + 3);
            return ss.str();
        }
};


class TreeSpecificOperatorScheduleSettings {

    template<typename T1> friend class BaseCollectionSettings;
    template<typename T> friend class BaseComparisonSettings;
    friend class RelativeRootComparisonSettings;

    protected:
        ScaleOperatorSettings mutation_rate_scaler_settings_ = ScaleOperatorSettings(
                1.0, 0.3);
        WindowOperatorSettings freq_mover_settings_ = WindowOperatorSettings(
                1.0, 0.1);
        ScaleOperatorSettings root_population_size_scaler_settings_ = ScaleOperatorSettings(
                1.0, 0.5);
        ScaleOperatorSettings leaf_population_size_scaler_settings_ = ScaleOperatorSettings(
                1.0, 0.5);

    public:
        TreeSpecificOperatorScheduleSettings() { }
        virtual ~TreeSpecificOperatorScheduleSettings() { }
        TreeSpecificOperatorScheduleSettings& operator=(const TreeSpecificOperatorScheduleSettings& other) {
            this->mutation_rate_scaler_settings_ = other.mutation_rate_scaler_settings_;
            this->freq_mover_settings_ = other.freq_mover_settings_;
            this->root_population_size_scaler_settings_ = other.root_population_size_scaler_settings_;
            this->leaf_population_size_scaler_settings_ = other.leaf_population_size_scaler_settings_;
            return * this;
        }

        virtual bool using_population_size_multipliers() const {
            return false;
        }

        const ScaleOperatorSettings& get_mutation_rate_scaler_settings() const {
            return this->mutation_rate_scaler_settings_;
        }
        const WindowOperatorSettings& get_freq_mover_settings() const {
            return this->freq_mover_settings_;
        }
        virtual const ScaleOperatorSettings& get_root_population_size_scaler_settings() const {
            return this->root_population_size_scaler_settings_;
        }
        virtual const ScaleOperatorSettings& get_leaf_population_size_scaler_settings() const {
            return this->leaf_population_size_scaler_settings_;
        }
        
        virtual const ScaleOperatorSettings& get_mean_population_size_scaler_settings() const {
            throw EcoevolityComparisonSettingError(
                    "cannot call get_mean_population_size_scaler_settings from TreeSpecificOperatorScheduleSettings");
        }
        virtual const ScaleOperatorSettings& get_relative_population_size_mixer_settings() const {
            throw EcoevolityComparisonSettingError(
                    "cannot call get_relative_population_size_mixer_settings from TreeSpecificOperatorScheduleSettings");
        }
        virtual const WindowOperatorSettings& get_root_relative_population_size_mover_settings() const {
            throw EcoevolityComparisonSettingError(
                    "cannot call get_root_relative_population_size_mover_settings from TreeSpecificOperatorScheduleSettings");
        }
        virtual const WindowOperatorSettings& get_leaf_relative_population_size_mover_settings() const {
            throw EcoevolityComparisonSettingError(
                    "cannot call get_leaf_relative_population_size_mover_settings from TreeSpecificOperatorScheduleSettings");
        }

        virtual void update_from_config(const YAML::Node& operators) {
            if (! operators.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "Expecting operators node to be a map, but found: " +
                        YamlCppUtils::get_node_type(operators));
            }

            for (YAML::const_iterator op = operators.begin();
                    op != operators.end();
                    ++op) {
                if (op->first.as<std::string>() == "MutationRateScaler") {
                    try {
                        this->mutation_rate_scaler_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing MutationRateScaler settings\n";
                        throw;
                    }
                }
                else if (op->first.as<std::string>() == "FreqMover") {
                    try {
                        this->freq_mover_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing FreqMover settings\n";
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
                else if (op->first.as<std::string>() == "LeafPopulationSizeScaler") {
                    try {
                        this->leaf_population_size_scaler_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing LeafPopulationSizeScaler settings\n";
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
            ss << margin << "operators:\n";
            ss << margin << indent << "RootPopulationSizeScaler:\n";
            ss << this->root_population_size_scaler_settings_.to_string(indent_level + 2);
            ss << margin << indent << "LeafPopulationSizeScaler:\n";
            ss << this->leaf_population_size_scaler_settings_.to_string(indent_level + 2);
            ss << margin << indent << "MutationRateScaler:\n";
            ss << this->mutation_rate_scaler_settings_.to_string(indent_level + 2);
            ss << margin << indent << "FreqMover:\n";
            ss << this->freq_mover_settings_.to_string(indent_level + 2);
            return ss.str();
        }
};


class DirichletTreeSpecificOperatorScheduleSettings : public TreeSpecificOperatorScheduleSettings {

    template<typename T1> friend class BaseCollectionSettings;
    template<typename T> friend class BaseComparisonSettings;
    friend class DirichletComparisonSettings;

    protected:
        ScaleOperatorSettings mean_population_size_scaler_settings_ = ScaleOperatorSettings(
                1.0, 0.5);
        ScaleOperatorSettings relative_population_size_mixer_settings_ = ScaleOperatorSettings(
                1.0, 0.01);
        WindowOperatorSettings root_relative_population_size_mover_settings_ = WindowOperatorSettings(
                1.0, 0.1);
        WindowOperatorSettings leaf_relative_population_size_mover_settings_ = WindowOperatorSettings(
                1.0, 0.1);

    public:
        DirichletTreeSpecificOperatorScheduleSettings() : TreeSpecificOperatorScheduleSettings() { }
        virtual ~DirichletTreeSpecificOperatorScheduleSettings() { }
        DirichletTreeSpecificOperatorScheduleSettings& operator=(const DirichletTreeSpecificOperatorScheduleSettings& other) {
            this->mutation_rate_scaler_settings_ = other.mutation_rate_scaler_settings_;
            this->freq_mover_settings_ = other.freq_mover_settings_;
            this->mean_population_size_scaler_settings_ = other.mean_population_size_scaler_settings_;
            this->relative_population_size_mixer_settings_ = other.relative_population_size_mixer_settings_;
            this->root_relative_population_size_mover_settings_ = other.root_relative_population_size_mover_settings_;
            this->leaf_relative_population_size_mover_settings_ = other.leaf_relative_population_size_mover_settings_;
            return * this;
        }

        bool using_population_size_multipliers() const {
            return true;
        }

        const ScaleOperatorSettings& get_root_population_size_scaler_settings() const {
            throw EcoevolityComparisonSettingError(
                    "cannot call get_root_population_size_scaler_settings from DirichletTreeSpecificOperatorScheduleSettings");
        }
        const ScaleOperatorSettings& get_leaf_population_size_scaler_settings() const {
            throw EcoevolityComparisonSettingError(
                    "cannot call get_leaf_population_size_scaler_settings from DirichletTreeSpecificOperatorScheduleSettings");
        }

        const ScaleOperatorSettings& get_mean_population_size_scaler_settings() const {
            return this->mean_population_size_scaler_settings_;
        }
        const ScaleOperatorSettings& get_relative_population_size_mixer_settings() const {
            return this->relative_population_size_mixer_settings_;
        }
        const WindowOperatorSettings& get_root_relative_population_size_mover_settings() const {
            return this->root_relative_population_size_mover_settings_;
        }
        const WindowOperatorSettings& get_leaf_relative_population_size_mover_settings() const {
            return this->leaf_relative_population_size_mover_settings_;
        }

        void update_from_config(const YAML::Node& operators) {
            if (! operators.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "Expecting operators node to be a map, but found: " +
                        YamlCppUtils::get_node_type(operators));
            }

            for (YAML::const_iterator op = operators.begin();
                    op != operators.end();
                    ++op) {
                if (op->first.as<std::string>() == "MutationRateScaler") {
                    try {
                        this->mutation_rate_scaler_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing MutationRateScaler settings\n";
                        throw;
                    }
                }
                else if (op->first.as<std::string>() == "FreqMover") {
                    try {
                        this->freq_mover_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing FreqMover settings\n";
                        throw;
                    }
                }
                else if (op->first.as<std::string>() == "MeanPopulationSizeScaler") {
                    try {
                        this->mean_population_size_scaler_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing MeanPopulationSizeScaler settings\n";
                        throw;
                    }
                }
                else if (op->first.as<std::string>() == "RelativePopulationSizeMixer") {
                    try {
                        this->relative_population_size_mixer_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing RelativePopulationSizeMixer settings\n";
                        throw;
                    }
                }
                else if (op->first.as<std::string>() == "RootRelativePopulationSizeMover") {
                    try {
                        this->root_relative_population_size_mover_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing RootRelativePopulationSizeMover settings\n";
                        throw;
                    }
                }
                else if (op->first.as<std::string>() == "LeafRelativePopulationSizeMover") {
                    try {
                        this->leaf_relative_population_size_mover_settings_.update_from_config(op->second);
                    }
                    catch (...) {
                        std::cerr << "ERROR: "
                                  << "Problem parsing LeafRelativePopulationSizeMover settings\n";
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
            ss << margin << "operators:\n";
            ss << margin << indent << "MeanPopulationSizeScaler:\n";
            ss << this->mean_population_size_scaler_settings_.to_string(indent_level + 2);
            ss << margin << indent << "RelativePopulationSizeMixer:\n";
            ss << this->relative_population_size_mixer_settings_.to_string(indent_level + 2);
            ss << margin << indent << "RootRelativePopulationSizeMover:\n";
            ss << this->root_relative_population_size_mover_settings_.to_string(indent_level + 2);
            ss << margin << indent << "LeafRelativePopulationSizeMover:\n";
            ss << this->leaf_relative_population_size_mover_settings_.to_string(indent_level + 2);
            ss << margin << indent << "MutationRateScaler:\n";
            ss << this->mutation_rate_scaler_settings_.to_string(indent_level + 2);
            ss << margin << indent << "FreqMover:\n";
            ss << this->freq_mover_settings_.to_string(indent_level + 2);
            return ss.str();
        }
};


template<class OperatorSettingsType>
class BaseComparisonSettings {

    template<typename T1> friend class BaseCollectionSettings;

    protected:
        std::string path_;
        PositiveRealParameterSettings population_size_settings_;
        PositiveRealParameterSettings freq_1_settings_;
        PositiveRealParameterSettings mutation_rate_settings_;

        OperatorSettingsType operator_settings_;

        char population_name_delimiter_ = ' ';
        bool population_name_is_prefix_ = true;
        bool genotypes_are_diploid_ = true;
        bool markers_are_dominant_ = false;
        bool constant_sites_removed_ = true;

        bool constrain_population_sizes_ = false;
        double ploidy_ = 2.0;

        virtual void make_consistent() {
            if (this->population_size_settings_.use_empirical_value()) {
                throw EcoevolityPositiveRealParameterSettingError(
                        "empirical value not supported for population_size");
            }
            if (this->mutation_rate_settings_.use_empirical_value()) {
                throw EcoevolityPositiveRealParameterSettingError(
                        "empirical value not supported for mutation_rate");
            }
            if (this->genotypes_are_diploid_ && (this->ploidy_ != 2.0)) {
                throw EcoevolityComparisonSettingError(
                        "Genotypes cannot be diploid when ploidy is not 2",
                        this->path_);
            }
        }
        
        virtual void update_operator_settings() {
            if (this->mutation_rate_settings_.is_fixed()) {
                this->operator_settings_.mutation_rate_scaler_settings_.set_weight(0.0);
            }
            if (this->freq_1_settings_.is_fixed()) {
                this->operator_settings_.freq_mover_settings_.set_weight(0.0);
            }
            if (this->population_size_settings_.is_fixed()) {
                this->operator_settings_.root_population_size_scaler_settings_.set_weight(0.0);
                this->operator_settings_.leaf_population_size_scaler_settings_.set_weight(0.0);
            }
            if (this->constrain_population_sizes_) {
                this->operator_settings_.leaf_population_size_scaler_settings_.set_weight(0.0);
            }
        }

        virtual void update_settings_contingent_upon_data() {
            if (this->freq_1_settings_.use_empirical_value_) {
                BiallelicData d = BiallelicData(
                        this->path_,
                        this->population_name_delimiter_,
                        this->population_name_is_prefix_,
                        this->genotypes_are_diploid_,
                        this->markers_are_dominant_,
                        true);
                this->freq_1_settings_.value_ = d.get_proportion_1();
            }
        }

        virtual void parse_parameter_settings(const YAML::Node& node) {
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
                else if (parameter->first.as<std::string>() == "freq_1") {
                    this->freq_1_settings_ = PositiveRealParameterSettings(parameter->second);
                }
                else if (parameter->first.as<std::string>() == "mutation_rate") {
                    this->mutation_rate_settings_ = PositiveRealParameterSettings(parameter->second);
                }
                else {
                    std::string message = "Unrecognized comparison parameter: " +
                            parameter->first.as<std::string>();
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

        void set_path(const std::string& path) {
            this->path_ = path;
        }

        virtual void update_from_config(
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
                else if (arg->first.as<std::string>() == "ploidy") {
                    this->ploidy_ = arg->second.as<double>();
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
                else if (arg->first.as<std::string>() == "equal_population_sizes") {
                    this->constrain_population_sizes_ = arg->second.as<bool>();
                }
                else if (arg->first.as<std::string>() == "parameters") {
                    this->parse_parameter_settings(arg->second);
                }
                else if (arg->first.as<std::string>() == "operators") {
                    this->operator_settings_.update_from_config(arg->second);
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
                this->update_settings_contingent_upon_data();
                this->update_operator_settings();
            }
        }

    public:
        BaseComparisonSettings() { }
        BaseComparisonSettings(
                const std::string& path,
                const PositiveRealParameterSettings& population_size_settings,
                const PositiveRealParameterSettings& freq_1_settings,
                const PositiveRealParameterSettings& mutation_rate_settings,
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool constrain_population_sizes = false,
                double ploidy = 2.0) {

            this->path_ = path;
            this->population_size_settings_ = population_size_settings;
            this->freq_1_settings_ = freq_1_settings;
            this->mutation_rate_settings_ = mutation_rate_settings;
            this->population_name_delimiter_ = population_name_delimiter;
            this->population_name_is_prefix_ = population_name_is_prefix;
            this->genotypes_are_diploid_ = genotypes_are_diploid;
            this->markers_are_dominant_ = markers_are_dominant;
            this->constant_sites_removed_ = constant_sites_removed;
            this->constrain_population_sizes_ = constrain_population_sizes;
            this->ploidy_ = ploidy;
            this->make_consistent();

            // TODO:
            // Not very efficient to parse data just to get rates, but might be
            // more awkward than it's worth to defer it.
            // Make a copy operator for BiallelicData and parse and store here, then
            // can copy it in get_instance method
            this->update_settings_contingent_upon_data();
            this->update_operator_settings();
        }
        BaseComparisonSettings(
                const YAML::Node& comparison_node,
                const std::string& config_path,
                bool global_defaults = false) {
            this->update_from_config(comparison_node, config_path, global_defaults);
        }

        virtual ~BaseComparisonSettings() { }
        BaseComparisonSettings& operator=(const BaseComparisonSettings& other) {
            this->path_                                = other.path_;
            this->population_size_settings_            = other.population_size_settings_;
            this->freq_1_settings_                     = other.freq_1_settings_;
            this->mutation_rate_settings_              = other.mutation_rate_settings_;
            this->population_name_delimiter_           = other.population_name_delimiter_;
            this->population_name_is_prefix_           = other.population_name_is_prefix_;
            this->genotypes_are_diploid_               = other.genotypes_are_diploid_;
            this->markers_are_dominant_                = other.markers_are_dominant_;
            this->constant_sites_removed_              = other.constant_sites_removed_;
            this->constrain_population_sizes_          = other.constrain_population_sizes_;
            this->ploidy_                              = other.ploidy_;
            this->operator_settings_                   = other.operator_settings_;
            return * this;
        }

        virtual bool using_population_size_multipliers() const {
            return false;
        }
        virtual bool using_relative_root_population_size() const {
            return false;
        }

        virtual bool relative_root_population_size_is_fixed() const {
            throw EcoevolityComparisonSettingError(
                    "ComparisonSettings does not support relative root size");
        }

        virtual bool population_size_multipliers_are_fixed() const {
            throw EcoevolityComparisonSettingError(
                    "ComparisonSettings does not support pop size multipliers");
        }

        double get_ploidy() const {
            return this->ploidy_;
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
        virtual bool constrain_population_sizes() const {
            return this->constrain_population_sizes_;
        }
        bool constrain_state_frequencies() const {
            return ((this->freq_1_settings_.is_fixed()) &&
                    (this->freq_1_settings_.get_value() == 0.5));
        }
        const PositiveRealParameterSettings& get_population_size_settings() const {
            return this->population_size_settings_;
        }
        const PositiveRealParameterSettings& get_freq_1_settings() const {
            return this->freq_1_settings_;
        }
        const PositiveRealParameterSettings& get_mutation_rate_settings() const {
            return this->mutation_rate_settings_;
        }
        const OperatorSettingsType& get_operator_settings() const {
            return this->operator_settings_;
        }

        virtual std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            ss << std::boolalpha;
            std::string margin = string_util::get_indent(indent_level);
            std::string indent = string_util::get_indent(1);
            ss << margin << "path: " << this->path_ << "\n";
            ss << margin << "ploidy: " << this->get_ploidy() << "\n";
            ss << margin << "genotypes_are_diploid: " << this->genotypes_are_diploid_ << "\n";
            ss << margin << "markers_are_dominant: " << this->markers_are_dominant_ << "\n";
            ss << margin << "population_name_delimiter: '" << this->population_name_delimiter_ << "'\n";
            ss << margin << "population_name_is_prefix: " << this->population_name_is_prefix_ << "\n";
            ss << margin << "constant_sites_removed: " << this->constant_sites_removed_ << "\n";
            ss << margin << "equal_population_sizes: " << this->constrain_population_sizes_ << "\n";
            ss << margin << "parameters:\n";

            ss << margin << indent << "population_size:\n";
            ss << this->population_size_settings_.to_string(indent_level + 2);

            ss << margin << indent << "mutation_rate:\n";
            ss << this->mutation_rate_settings_.to_string(indent_level + 2);

            ss << margin << indent << "freq_1:\n";
            ss << this->freq_1_settings_.to_string(indent_level + 2);

            ss << this->operator_settings_.to_string(indent_level);

            return ss.str();
        }
};


class ComparisonSettings : public BaseComparisonSettings<TreeSpecificOperatorScheduleSettings> {

    public:
        ComparisonSettings() : BaseComparisonSettings<TreeSpecificOperatorScheduleSettings>() { }
        ComparisonSettings(
                const std::string& path,
                const PositiveRealParameterSettings& population_size_settings,
                const PositiveRealParameterSettings& freq_1_settings,
                const PositiveRealParameterSettings& mutation_rate_settings,
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool constrain_population_sizes = false,
                double ploidy = 2.0)
            : BaseComparisonSettings<TreeSpecificOperatorScheduleSettings>(
                    path,
                    population_size_settings,
                    freq_1_settings,
                    mutation_rate_settings,
                    population_name_delimiter,
                    population_name_is_prefix,
                    genotypes_are_diploid,
                    markers_are_dominant,
                    constant_sites_removed,
                    constrain_population_sizes,
                    ploidy) { }
        ComparisonSettings(
                const YAML::Node& comparison_node,
                const std::string& config_path,
                bool global_defaults = false)
            : BaseComparisonSettings<TreeSpecificOperatorScheduleSettings>(
                    comparison_node,
                    config_path,
                    global_defaults) { }
};

class RelativeRootComparisonSettings : public BaseComparisonSettings<TreeSpecificOperatorScheduleSettings> {

    template<typename T1> friend class BaseCollectionSettings;

    protected:
        PositiveRealParameterSettings relative_root_population_size_settings_;

        void make_consistent() {
            if (this->population_size_settings_.use_empirical_value()) {
                throw EcoevolityPositiveRealParameterSettingError(
                        "empirical value not supported for population_size");
            }
            if (this->mutation_rate_settings_.use_empirical_value()) {
                throw EcoevolityPositiveRealParameterSettingError(
                        "empirical value not supported for mutation_rate");
            }
            if (this->genotypes_are_diploid_ && (this->ploidy_ != 2.0)) {
                throw EcoevolityComparisonSettingError(
                        "Genotypes cannot be diploid when ploidy is not 2",
                        this->path_);
            }
            if (this->constrain_population_sizes_) {
                this->relative_root_population_size_settings_.value_ = 1.0;
                this->relative_root_population_size_settings_.is_fixed_ = true;
                this->relative_root_population_size_settings_.prior_settings_.nullify();
            }
        }
        
        void update_relative_root_population_size_settings() {
            if ((! this->relative_root_population_size_settings_.is_fixed_) &&
                    (this->relative_root_population_size_settings_.prior_settings_.get_name() == "none")) {
                std::stringstream ss;
                ss << "gamma_distribution:\n";
                ss << "    shape: 50.0\n";
                ss << "    scale: 0.02\n";
                YAML::Node n;
                n = YAML::Load(ss);
                this->relative_root_population_size_settings_.prior_settings_ = ContinuousDistributionSettings(n);
            }
        }

        void update_operator_settings() {
            if (this->mutation_rate_settings_.is_fixed()) {
                this->operator_settings_.mutation_rate_scaler_settings_.set_weight(0.0);
            }
            if (this->freq_1_settings_.is_fixed()) {
                this->operator_settings_.freq_mover_settings_.set_weight(0.0);
            }
            if (this->population_size_settings_.is_fixed()) {
                this->operator_settings_.leaf_population_size_scaler_settings_.set_weight(0.0);
            }
            if (this->constrain_population_sizes_) {
                this->operator_settings_.leaf_population_size_scaler_settings_.set_weight(0.0);
            }
            if (this->relative_root_population_size_settings_.is_fixed()) {
                if (! this->constrain_population_sizes_) {
                    this->operator_settings_.root_population_size_scaler_settings_.set_weight(0.0);
                }
                if (this->constrain_population_sizes_ && this->population_size_settings_.is_fixed()) {
                    this->operator_settings_.root_population_size_scaler_settings_.set_weight(0.0);
                }
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
                if (parameter->first.as<std::string>() == "root_relative_population_size") {
                    this->relative_root_population_size_settings_ = PositiveRealParameterSettings(parameter->second);
                }
                else if (parameter->first.as<std::string>() == "population_size") {
                    this->population_size_settings_ = PositiveRealParameterSettings(parameter->second);
                }
                else if (parameter->first.as<std::string>() == "freq_1") {
                    this->freq_1_settings_ = PositiveRealParameterSettings(parameter->second);
                }
                else if (parameter->first.as<std::string>() == "mutation_rate") {
                    this->mutation_rate_settings_ = PositiveRealParameterSettings(parameter->second);
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
                else if (arg->first.as<std::string>() == "ploidy") {
                    this->ploidy_ = arg->second.as<double>();
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
                else if (arg->first.as<std::string>() == "equal_population_sizes") {
                    this->constrain_population_sizes_ = arg->second.as<bool>();
                }
                else if (arg->first.as<std::string>() == "parameters") {
                    this->parse_parameter_settings(arg->second);
                }
                else if (arg->first.as<std::string>() == "operators") {
                    this->operator_settings_.update_from_config(arg->second);
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
                this->update_settings_contingent_upon_data();
                this->update_operator_settings();
            }
            this->update_relative_root_population_size_settings();
        }

    public:
        RelativeRootComparisonSettings() { }
        RelativeRootComparisonSettings(
                const YAML::Node& comparison_node,
                const std::string& config_path,
                bool global_defaults = false) {
            this->update_from_config(comparison_node, config_path, global_defaults);
        }

        RelativeRootComparisonSettings& operator=(const RelativeRootComparisonSettings& other) {
            this->path_                                    = other.path_;
            this->population_size_settings_                = other.population_size_settings_;
            this->relative_root_population_size_settings_  = other.relative_root_population_size_settings_;
            this->freq_1_settings_                         = other.freq_1_settings_;
            this->mutation_rate_settings_                  = other.mutation_rate_settings_;
            this->population_name_delimiter_               = other.population_name_delimiter_;
            this->population_name_is_prefix_               = other.population_name_is_prefix_;
            this->genotypes_are_diploid_                   = other.genotypes_are_diploid_;
            this->markers_are_dominant_                    = other.markers_are_dominant_;
            this->constant_sites_removed_                  = other.constant_sites_removed_;
            this->ploidy_                                  = other.ploidy_;
            this->operator_settings_                       = other.operator_settings_;
            return * this;
        }

        const PositiveRealParameterSettings& get_relative_root_population_size_settings() const {
            return this->relative_root_population_size_settings_;
        }

        bool using_relative_root_population_size() const {
            return true;
        }

        bool relative_root_population_size_is_fixed() const {
            return this->relative_root_population_size_settings_.is_fixed();
        }

        virtual std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            ss << std::boolalpha;
            std::string margin = string_util::get_indent(indent_level);
            std::string indent = string_util::get_indent(1);
            ss << margin << "path: " << this->path_ << "\n";
            ss << margin << "ploidy: " << this->get_ploidy() << "\n";
            ss << margin << "genotypes_are_diploid: " << this->genotypes_are_diploid_ << "\n";
            ss << margin << "markers_are_dominant: " << this->markers_are_dominant_ << "\n";
            ss << margin << "population_name_delimiter: '" << this->population_name_delimiter_ << "'\n";
            ss << margin << "population_name_is_prefix: " << this->population_name_is_prefix_ << "\n";
            ss << margin << "constant_sites_removed: " << this->constant_sites_removed_ << "\n";
            ss << margin << "equal_population_sizes: " << this->constrain_population_sizes_ << "\n";
            ss << margin << "parameters:\n";

            ss << margin << indent << "population_size:\n";
            ss << this->population_size_settings_.to_string(indent_level + 2);

            ss << margin << indent << "root_relative_population_size:\n";
            ss << this->relative_root_population_size_settings_.to_string(indent_level + 2);

            ss << margin << indent << "mutation_rate:\n";
            ss << this->mutation_rate_settings_.to_string(indent_level + 2);

            ss << margin << indent << "freq_1:\n";
            ss << this->freq_1_settings_.to_string(indent_level + 2);

            ss << this->operator_settings_.to_string(indent_level);

            return ss.str();
        }
};


class DirichletComparisonSettings : public BaseComparisonSettings<DirichletTreeSpecificOperatorScheduleSettings> {

    template<typename T1> friend class BaseCollectionSettings;

    protected:
        std::vector<std::string> population_labels_;
        PositiveRealParameterSettings population_size_multiplier_settings_;
        bool pop_size_multipliers_specified_ = false;

        void make_consistent() {
            if (this->population_size_settings_.use_empirical_value()) {
                throw EcoevolityPositiveRealParameterSettingError(
                        "empirical value not supported for population_size");
            }
            if (this->population_size_multiplier_settings_.use_empirical_value()) {
                throw EcoevolityPositiveRealParameterSettingError(
                        "empirical value not supported for population_size_multipliers");
            }
            if (this->mutation_rate_settings_.use_empirical_value()) {
                throw EcoevolityPositiveRealParameterSettingError(
                        "empirical value not supported for mutation_rate");
            }
            if (this->genotypes_are_diploid_ && (this->ploidy_ != 2.0)) {
                throw EcoevolityComparisonSettingError(
                        "Genotypes cannot be diploid when ploidy is not 2",
                        this->path_);
            }
            std::string multiplier_prior_name = this->population_size_multiplier_settings_.prior_settings_.get_name();
            if ((multiplier_prior_name != "none") &&
                (multiplier_prior_name != "dirichlet_distribution")) {
                std::string message = "Invalid prior distribution for population size multipliers: " + multiplier_prior_name;
                throw EcoevolityComparisonSettingError(message);
            }
        }

        void update_operator_settings() {
            if (this->mutation_rate_settings_.is_fixed()) {
                this->operator_settings_.mutation_rate_scaler_settings_.set_weight(0.0);
            }
            if (this->freq_1_settings_.is_fixed()) {
                this->operator_settings_.freq_mover_settings_.set_weight(0.0);
            }
            if (this->population_size_settings_.is_fixed()) {
                this->operator_settings_.mean_population_size_scaler_settings_.set_weight(0.0);
            }
            if (this->population_size_multiplier_settings_.is_fixed()) {
                this->operator_settings_.relative_population_size_mixer_settings_.set_weight(0.0);
                this->operator_settings_.root_relative_population_size_mover_settings_.set_weight(0.0);
                this->operator_settings_.leaf_relative_population_size_mover_settings_.set_weight(0.0);
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
                else if (arg->first.as<std::string>() == "ploidy") {
                    this->ploidy_ = arg->second.as<double>();
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
                else if (arg->first.as<std::string>() == "parameters") {
                    this->parse_parameter_settings(arg->second);
                }
                else if (arg->first.as<std::string>() == "operators") {
                    this->operator_settings_.update_from_config(arg->second);
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
                this->update_settings_contingent_upon_data();
                this->update_operator_settings();
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
                else if (parameter->first.as<std::string>() == "population_size_multipliers") {
                    this->population_size_multiplier_settings_ = PositiveRealParameterSettings(parameter->second);
                    this->pop_size_multipliers_specified_ = true;
                }
                else if (parameter->first.as<std::string>() == "freq_1") {
                    this->freq_1_settings_ = PositiveRealParameterSettings(parameter->second);
                }
                else if (parameter->first.as<std::string>() == "mutation_rate") {
                    this->mutation_rate_settings_ = PositiveRealParameterSettings(parameter->second);
                }
                else {
                    std::string message = "Unrecognized comparison parameter: " +
                            parameter->first.as<std::string>();
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

        void update_settings_contingent_upon_data() {
            BiallelicData d = BiallelicData(
                    this->path_,
                    this->population_name_delimiter_,
                    this->population_name_is_prefix_,
                    this->genotypes_are_diploid_,
                    this->markers_are_dominant_,
                    true);
            unsigned int npops = d.get_number_of_populations();
            if (npops < 1) {
                throw EcoevolityParsingError(
                        "No populations found in nexus file",
                        this->path_);
            }
            else if (npops > 2) {
                throw EcoevolityParsingError(
                        "More than two populations found in nexus file",
                        this->path_);
            }
            for (unsigned int i = 0; i < npops; ++i) {
                this->population_labels_.push_back(d.get_population_label(i));
            }
            if (this->pop_size_multipliers_specified_) {
                unsigned int num_dirichlet_vals = this->population_size_multiplier_settings_.values_.size();
                if (num_dirichlet_vals > 0) {
                    if (num_dirichlet_vals != (this->population_labels_.size() + 1)) {
                        throw EcoevolityParsingError(
                                "Mismatch between number of multiplier values in config and number of populations in nexus file",
                                this->path_);

                    }
                }
                if (this->population_size_multiplier_settings_.prior_settings_.get_name() == "dirichlet_distribution") {
                    if (this->population_size_multiplier_settings_.prior_settings_.get_concentration_parameters().size() != (this->population_labels_.size() + 1)) {
                        throw EcoevolityParsingError(
                                "Mismatch between number of dirichlet alphas in config and number of populations in nexus file",
                                this->path_);
                    }
                }
            }
            else {
                this->population_size_multiplier_settings_.is_vector_ = true;
                this->population_size_multiplier_settings_.values_.clear();
                this->population_size_multiplier_settings_.is_fixed_ = false;
                this->population_size_multiplier_settings_.prior_settings_.nullify();
            }
            if ((! this->population_size_multiplier_settings_.is_fixed_) &&
                    (this->population_size_multiplier_settings_.prior_settings_.get_name() == "none")) {
                std::stringstream ss;
                ss << "dirichlet_distribution:\n";
                ss << "    alpha: [10.0";
                for (unsigned int i = 1; i < (npops + 1); ++i) {
                    ss << ", 10.0";
                }
                ss << "]\n";
                YAML::Node n;
                n = YAML::Load(ss);
                this->population_size_multiplier_settings_.prior_settings_ = ContinuousDistributionSettings(n);
            }
            if (this->freq_1_settings_.use_empirical_value_) {
                this->freq_1_settings_.value_ = d.get_proportion_1();
            }
        }

    public:
        DirichletComparisonSettings() { }
        DirichletComparisonSettings(
                const YAML::Node& comparison_node,
                const std::string& config_path,
                bool global_defaults = false) {
            this->update_from_config(comparison_node, config_path, global_defaults);
        }

        DirichletComparisonSettings& operator=(const DirichletComparisonSettings& other) {
            this->path_                                    = other.path_;
            this->population_size_settings_                = other.population_size_settings_;
            this->population_size_multiplier_settings_     = other.population_size_multiplier_settings_;
            this->population_labels_                       = other.population_labels_;
            this->freq_1_settings_                         = other.freq_1_settings_;
            this->mutation_rate_settings_                  = other.mutation_rate_settings_;
            this->population_name_delimiter_               = other.population_name_delimiter_;
            this->population_name_is_prefix_               = other.population_name_is_prefix_;
            this->genotypes_are_diploid_                   = other.genotypes_are_diploid_;
            this->markers_are_dominant_                    = other.markers_are_dominant_;
            this->constant_sites_removed_                  = other.constant_sites_removed_;
            this->ploidy_                                  = other.ploidy_;
            this->operator_settings_                       = other.operator_settings_;
            return * this;
        }

        bool using_population_size_multipliers() const {
            return true;
        }

        bool population_size_multipliers_are_fixed() const {
            return this->population_size_multiplier_settings_.is_fixed();
        }

        const PositiveRealParameterSettings& get_population_size_multiplier_settings() const {
            return this->population_size_multiplier_settings_;
        }

        bool constrain_population_sizes() const {
            throw EcoevolityError(
                    "Do not use this method for comparison settings that support population size multipliers");
            // if (this->population_size_multiplier_settings_.values_.size() < 1) {
            //     return false;
            // }
            // std::vector<double> v (this->population_size_multiplier_settings_.values_.size(), 1.0);
            // return ((this->population_size_multiplier_settings_.is_fixed()) &&
            //     (this->population_size_multiplier_settings_.values_ == v));
        }

        std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            ss << std::boolalpha;
            std::string margin = string_util::get_indent(indent_level);
            std::string indent = string_util::get_indent(1);
            ss << margin << "path: " << this->path_ << "\n";
            ss << margin << "ploidy: " << this->get_ploidy() << "\n";
            ss << margin << "genotypes_are_diploid: " << this->genotypes_are_diploid_ << "\n";
            ss << margin << "markers_are_dominant: " << this->markers_are_dominant_ << "\n";
            ss << margin << "population_name_delimiter: '" << this->population_name_delimiter_ << "'\n";
            ss << margin << "population_name_is_prefix: " << this->population_name_is_prefix_ << "\n";
            ss << margin << "constant_sites_removed: " << this->constant_sites_removed_ << "\n";
            ss << margin << "parameters:\n";

            ss << margin << indent << "population_size:\n";
            ss << this->population_size_settings_.to_string(indent_level + 2);

            ss << margin << indent << "population_size_multipliers:\n";
            ss << margin << indent << indent
               << "# Multiplier settings map to [" << this->population_labels_.at(0);
            for (unsigned int i = 1; i < this->population_labels_.size(); ++i) {
                ss << ", " << this->population_labels_.at(i);
            }
            ss << ", root]\n";
            ss << this->population_size_multiplier_settings_.to_string(indent_level + 2);

            ss << margin << indent << "mutation_rate:\n";
            ss << this->mutation_rate_settings_.to_string(indent_level + 2);

            ss << margin << indent << "freq_1:\n";
            ss << this->freq_1_settings_.to_string(indent_level + 2);

            ss << this->operator_settings_.to_string(indent_level);

            return ss.str();
        }
};


template<class ComparisonSettingsType>
class BaseCollectionSettings {

    public:

        BaseCollectionSettings() {
            this->init_default_priors();
        }
        BaseCollectionSettings(bool use_population_size_multipliers) {
            if (use_population_size_multipliers) {
                this->operator_schedule_settings_.turn_on_population_size_multipliers();
            }
            this->init_default_priors();
        }
        BaseCollectionSettings(
                const ContinuousDistributionSettings& time_prior,
                unsigned int chain_length,
                unsigned int sample_frequency,
                const PositiveRealParameterSettings& concentration_settings,
                bool use_dpp,
                bool use_population_size_multipliers = false)
                : BaseCollectionSettings(use_population_size_multipliers) {
            this->time_prior_settings_ = time_prior;
            this->chain_length_ = chain_length;
            this->sample_frequency_ = sample_frequency;
            this->concentration_settings_ = concentration_settings;
            this->use_dpp_ = use_dpp;
        }
        BaseCollectionSettings(
                const std::string & yaml_config_path,
                bool use_population_size_multipliers = false)
                : BaseCollectionSettings(use_population_size_multipliers) {
            this->init_from_config_file(yaml_config_path);
        }
        BaseCollectionSettings(
                std::istream& yaml_config_stream,
                const std::string& yaml_config_path,
                bool use_population_size_multipliers = false)
                : BaseCollectionSettings(use_population_size_multipliers) {
            this->init_from_config_stream(yaml_config_stream, yaml_config_path);
        }
        virtual ~BaseCollectionSettings() { }
        BaseCollectionSettings& operator=(const BaseCollectionSettings& other) {
            this->path_ = other.path_;
            this->operator_schedule_settings_ = other.operator_schedule_settings_;
            this->global_comparison_settings_ = other.global_comparison_settings_;
            this->time_prior_settings_ = other.time_prior_settings_;
            this->chain_length_ = other.chain_length_;
            this->sample_frequency_ = other.sample_frequency_;
            this->concentration_settings_ = other.concentration_settings_;
            this->use_dpp_ = other.use_dpp_;
            this->fixed_event_model_indices_ = other.fixed_event_model_indices_;
            this->comparisons_ = other.comparisons_;
            this->default_time_prior_ = other.default_time_prior_;
            this->default_population_size_prior_ = other.default_population_size_prior_;
            this->default_freq_1_prior_ = other.default_freq_1_prior_;
            this->default_mutation_rate_prior_ = other.default_mutation_rate_prior_;
            this->state_log_path_ = other.state_log_path_;
            this->operator_log_path_ = other.operator_log_path_;
            return * this;
        }

        virtual bool using_population_size_multipliers() const {
            return false;
        }

        void add_comparison(ComparisonSettingsType comparison_settings) {
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

        unsigned int get_number_of_comparisons_with_free_mutation_rate() const {
            unsigned int nfree = 0;
            for (const ComparisonSettingsType& comparison : this->comparisons_) {
                if (! comparison.mutation_rate_settings_.is_fixed()) {
                    ++nfree;
                }
            }
            return nfree;
        }

        unsigned int get_number_of_comparisons_with_free_population_size() const {
            unsigned int nfree = 0;
            for (const ComparisonSettingsType& comparison : this->comparisons_) {
                if (! comparison.population_size_settings_.is_fixed() ||
                        (comparison.using_relative_root_population_size() &&
                        (! comparison.relative_root_population_size_is_fixed()))) {
                    ++nfree;
                }
            }
            return nfree;
        }

        unsigned int get_number_of_comparisons_with_free_population_size_multipliers() const {
            unsigned int nfree = 0;
            for (const ComparisonSettingsType& comparison : this->comparisons_) {
                if ((comparison.using_population_size_multipliers()) && (! comparison.population_size_multipliers_are_fixed())) {
                    ++nfree;
                }
            }
            return nfree;
        }

        unsigned int get_number_of_comparisons_with_free_state_frequencies() const {
            unsigned int nfree = 0;
            for (const ComparisonSettingsType& comparison : this->comparisons_) {
                if (! comparison.freq_1_settings_.is_fixed()) {
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

        const std::vector<ComparisonSettingsType>& get_comparison_settings() const { 
            return this->comparisons_;
        }

        const ComparisonSettingsType& get_comparison_setting(
                unsigned int comparison_index) const { 
            ECOEVOLITY_ASSERT(comparison_index < this->comparisons_.size());
            return this->comparisons_.at(comparison_index);
        }

        bool same_comparison_paths(const BaseCollectionSettings& s) const {
            if (this->get_number_of_comparisons() != s.get_number_of_comparisons()) {
                return false;
            }
            std::unordered_set<std::string> paths;
            std::unordered_set<std::string> other_paths;
            for (unsigned int i = 0; i < this->get_number_of_comparisons(); ++i) {
                paths.insert(this->get_comparison_setting(i).get_path());
                other_paths.insert(s.get_comparison_setting(i).get_path());
            }
            if (paths == other_paths) {
                return true;
            }
            return false;
        }

        void replace_comparison_path(
                const std::string & original_path,
                const std::string & new_path) {
            for (auto & comp : this->comparisons_) {
                if (comp.get_path() == original_path) {
                    comp.set_path(new_path);
                    return;
                }
            }
            throw EcoevolityComparisonSettingError(
                    "Comparison path \'" + original_path +
                    "\' was not found in collection settings");
        }

        char get_population_name_delimiter(const std::string& comparison_path) {
            for (auto comp : this->comparisons_) {
                if (comp.get_path() == comparison_path) {
                    return comp.get_population_name_delimiter();
                }
            }
            throw EcoevolityComparisonSettingError(
                    "Comparison path \'" + comparison_path +
                    "\' was not found in collection settings");
        }

        const OperatorScheduleSettings& get_operator_schedule_settings() const {
            return this->operator_schedule_settings_;
        }

        void write_settings(std::ostream& out) const {
            std::string indent = string_util::get_indent(1);
            out << std::boolalpha;

            out << "---\n"
                << "event_model_prior:\n";
            if (this->event_model_is_fixed()) {
                out << indent << "fixed: [" << this->fixed_event_model_indices_.at(0);
                for (unsigned int i = 1; i < this->fixed_event_model_indices_.size(); ++i) {
                    out << ", " << this->fixed_event_model_indices_.at(i);
                }
                out << "]\n";
            }
            else if (this->use_dpp_) {
                out << indent << "dirichlet_process:\n"
                    << indent << indent << "parameters:\n"
                    << indent << indent <<  indent << "concentration:\n";
                out << this->concentration_settings_.to_string(4);
            }
            else {
                out << indent << "uniform:\n"
                    << indent << indent << "parameters:\n"
                    << indent << indent <<  indent << "split_weight:\n";
                out << this->concentration_settings_.to_string(4);
            }

            out << "event_time_prior:\n";
            out << this->time_prior_settings_.to_string(1);

            out << "mcmc_settings:\n"
                << indent << "chain_length: " << this->chain_length_ << "\n"
                << indent << "sample_frequency: " << this->sample_frequency_ << "\n";

            // out << "output_settings:\n"
            //     << indent << "state_log_path: " << this->state_log_path_ << "\n"
            //     << indent << "operator_log_path: " << this->operator_log_path_ << "\n";

            out << "comparisons:\n";
            for (auto comp : this->comparisons_) {
                out << "- comparison:\n";
                out << comp.to_string(1);
            }

            out << this->operator_schedule_settings_.to_string();
        }

        std::string to_string() const {
            std::ostringstream ss;
            this->write_settings(ss);
            return ss.str();
        }

        void blanket_set_population_name_is_prefix(bool b) {
            for (auto comp : this->comparisons_) {
                comp.population_name_is_prefix_ = b;
            }
        }

        bool event_model_is_fixed() const {
            if (this->fixed_event_model_indices_.size() > 0) {
                return true;
            }
            return false;
        }

        const std::vector<unsigned int>& get_fixed_event_model_indices() const {
            return this->fixed_event_model_indices_;
        }


    protected:

        std::string path_ = "";
        bool use_dpp_ = true;
        std::vector<unsigned int> fixed_event_model_indices_;
        unsigned int chain_length_ = 100000;
        unsigned int sample_frequency_ = 100;
        std::string state_log_path_ = "ecoevolity-state-run-1.log";
        std::string operator_log_path_ = "ecoevolity-operator-run-1.log";

        OperatorScheduleSettings operator_schedule_settings_;

        ContinuousDistributionSettings time_prior_settings_;

        PositiveRealParameterSettings concentration_settings_;

        ComparisonSettingsType global_comparison_settings_;

        std::vector<ComparisonSettingsType> comparisons_;

        ContinuousDistributionSettings default_time_prior_;
        ContinuousDistributionSettings default_population_size_prior_;
        ContinuousDistributionSettings default_freq_1_prior_;
        ContinuousDistributionSettings default_mutation_rate_prior_;

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
            default_parameters["alpha"] = 1.0;
            default_parameters["beta"] = 1.0;
            this->default_freq_1_prior_ = ContinuousDistributionSettings(
                    "beta_distribution",
                    default_parameters);
            default_parameters.clear();
            default_parameters["shape"] = 1000.0;
            default_parameters["scale"] = 0.001;
            this->default_mutation_rate_prior_ = ContinuousDistributionSettings(
                    "gamma_distribution",
                    default_parameters);
        }

        void set_output_paths_to_config_directory() {
            std::pair<std::string, std::string> prefix_ext = path::splitext(this->path_);
            this->state_log_path_ = prefix_ext.first + "-state-run-1.log";
            this->operator_log_path_ = prefix_ext.first + "-operator-run-1.log";
        }

        void init_from_config_stream(std::istream& stream, const std::string& path) {
            this->path_ = path;
            this->set_output_paths_to_config_directory();
            this->parse_yaml_config(stream);
            ECOEVOLITY_ASSERT(this->operator_schedule_settings_.using_population_size_multipliers() == this->comparisons_.at(0).using_population_size_multipliers());
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
            if (this->event_model_is_fixed() &&
                    (this->fixed_event_model_indices_.size() != this->comparisons_.size())) {
                throw EcoevolityYamlConfigError(
                        "The number of indices in the fixed event model is incorrect");
            }
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
            // Set default mutation rate
            this->global_comparison_settings_.mutation_rate_settings_.value_ = 1.0;
            this->global_comparison_settings_.mutation_rate_settings_.is_fixed_ = true;
            // Set default freq 1
            this->global_comparison_settings_.freq_1_settings_.value_ = 0.5;
            this->global_comparison_settings_.freq_1_settings_.is_fixed_ = true;
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
                this->concentration_settings_.value_ = 1.0;
                this->concentration_settings_.is_fixed_ = true;
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
                // else if (top->first.as<std::string>() == "output_settings") {
                //     this->parse_output_settings(top->second);
                // }
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
            if ((this->comparisons_.size() < 2) || (this->event_model_is_fixed())) {
                this->operator_schedule_settings_.model_operator_settings_.set_weight(0.0);
            }
            if (this->operator_schedule_settings_.using_population_size_multipliers()) {
                if ((this->get_number_of_comparisons_with_free_mutation_rate() < 1) && 
                    (this->get_number_of_comparisons_with_free_population_size() < 1)) {
                    this->operator_schedule_settings_.composite_time_size_rate_scaler_settings_.set_weight(0.0);
                }
                if ((this->get_number_of_comparisons_with_free_mutation_rate() < 1) && 
                    (this->get_number_of_comparisons_with_free_population_size() < 1) &&
                    (this->get_number_of_comparisons_with_free_population_size_multipliers() < 1)) {
                    this->operator_schedule_settings_.composite_time_size_rate_mixer_settings_.set_weight(0.0);
                }
            }
            else {
                if ((this->get_number_of_comparisons_with_free_mutation_rate() < 1) && 
                    (this->get_number_of_comparisons_with_free_population_size() < 1)) {
                    this->operator_schedule_settings_.composite_time_size_rate_mixer_settings_.set_weight(0.0);
                    this->operator_schedule_settings_.composite_time_size_rate_scaler_settings_.set_weight(0.0);
                }
            }
        }

        void update_default_comparison_priors() {
            for (auto&& comp : this->comparisons_) {
                if ((! comp.mutation_rate_settings_.is_fixed()) &&
                        (comp.mutation_rate_settings_.prior_settings_.get_name() == "none")) {
                    comp.mutation_rate_settings_.prior_settings_ = this->default_mutation_rate_prior_;
                }
                if ((! comp.freq_1_settings_.is_fixed()) &&
                        (comp.freq_1_settings_.prior_settings_.get_name() == "none")) {
                    comp.freq_1_settings_.prior_settings_ = this->default_freq_1_prior_;
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

        // void parse_output_settings(const YAML::Node& output_node) {
        //     if (! output_node.IsMap()) {
        //         throw EcoevolityYamlConfigError(
        //                 "Expecting output_settings to be a map, but found: " +
        //                 YamlCppUtils::get_node_type(output_node));
        //     }
        //     for (YAML::const_iterator output = output_node.begin(); output != output_node.end(); ++output) {
        //         if (output->first.as<std::string>() == "state_log_path") {
        //             this->state_log_path_ = path::join(
        //                     path::dirname(this->path_),
        //                     output->second.as<std::string>());

        //         }
        //         else if (output->first.as<std::string>() == "operator_log_path") {
        //             this->operator_log_path_ = path::join(
        //                     path::dirname(this->path_),
        //                     output->second.as<std::string>());
        //         }
        //         else {
        //             throw EcoevolityYamlConfigError(
        //                     "Unrecognized output_settings key: " +
        //                     output->first.as<std::string>());
        //         }
        //     }
        // }

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
                this->concentration_settings_.value_ = 1.0;
                this->concentration_settings_.is_fixed_ = true;
                this->parse_uniform_model_prior(model_prior_node["uniform"]);
            }
            else if(model_prior_node["fixed"]) {
                this->use_dpp_ = false;
                this->parse_fixed_event_model(model_prior_node["fixed"]);
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
                ComparisonSettingsType c = this->global_comparison_settings_;
                c.update_from_config(comp["comparison"], this->path_);
                this->comparisons_.push_back(c);
            }
            ECOEVOLITY_ASSERT(comparisons_node.size() == this->comparisons_.size());
        }

        void parse_fixed_event_model(const YAML::Node& model_node) {
            if (! model_node.IsSequence()) {
                throw EcoevolityYamlConfigError(
                        "Expecting fixed event model node to be a sequence, but found: " +
                        YamlCppUtils::get_node_type(model_node));
            }
            for (std::size_t i = 0; i < model_node.size(); ++i) {
                this->fixed_event_model_indices_.push_back(model_node[i].as<unsigned int>());
            }
        }

        void parse_dirichlet_process_prior(const YAML::Node& dpp_node) {
            if (! dpp_node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "Expecting dirichlet_process to be a map, but found: " +
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

        void parse_uniform_model_prior(const YAML::Node& uniform_node) {
            if ((uniform_node.IsNull()) || (uniform_node.size() < 1)) {
                return;
            }
            if (! uniform_node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "Expecting event_model_prior -> uniform to be a map, but found: " +
                        YamlCppUtils::get_node_type(uniform_node));
            }
            if (uniform_node.size() != 1) {
                throw EcoevolityYamlConfigError(
                        "event_model_prior -> uniform node must have a single key");
            }
            if (! uniform_node["parameters"]) {
                return;
            }
            if ((uniform_node["parameters"].size() != 1) || (! uniform_node["parameters"]["split_weight"])) {
                throw EcoevolityYamlConfigError(
                        "split_weight is the only parameter setting for event_model_prior -> uniform");
            }
            this->parse_split_weight_parameter(uniform_node["parameters"]["split_weight"]);
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
                    if (arg->second.as<std::string>() == "empirical") {
                        throw EcoevolityPositiveRealParameterSettingError(
                                "empirical value not supported for concentration_parameter");
                    }
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

        void parse_split_weight_parameter(const YAML::Node& node) {
            if (! node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "Expecting split_weight node to be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
            }
            if (node.size() < 1) {
                throw EcoevolityYamlConfigError(
                        "empty split_weight parameter node");
            }

            for (YAML::const_iterator arg = node.begin();
                    arg != node.end();
                    ++arg) {
                if (arg->first.as<std::string>() == "value") {
                    double v = arg->second.as<double>();
                        if (v < 0.0) {
                            throw EcoevolityPositiveRealParameterSettingError(
                                    "split_weight parameter cannot be less than 0"
                                    );
                        }
                    this->concentration_settings_.value_ = v;
                }
                else if (arg->first.as<std::string>() == "estimate") {
                    bool f = arg->second.as<bool>();
                    this->concentration_settings_.is_fixed_ = (! f);
                }
                else if (arg->first.as<std::string>() == "prior") {
                    this->concentration_settings_.prior_settings_ = ContinuousDistributionSettings(arg->second);
                }
                else {
                    std::string message = "Unrecognized split_weight key: " +
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


class CollectionSettings: public BaseCollectionSettings<ComparisonSettings>{
    private:
        typedef BaseCollectionSettings<ComparisonSettings> BaseClass;

    public:
        CollectionSettings() : BaseClass() { }
        CollectionSettings(
                const ContinuousDistributionSettings& time_prior,
                unsigned int chain_length,
                unsigned int sample_frequency,
                const PositiveRealParameterSettings& concentration_settings,
                bool use_dpp)
                : BaseClass(
                        time_prior,
                        chain_length,
                        sample_frequency,
                        concentration_settings,
                        use_dpp) { }
        CollectionSettings(
                const std::string & yaml_config_path)
                : BaseClass(
                        yaml_config_path) { }
        CollectionSettings(
                std::istream& yaml_config_stream,
                const std::string& yaml_config_path)
                : BaseClass(
                        yaml_config_stream,
                        yaml_config_path) { }
};


class RelativeRootCollectionSettings: public BaseCollectionSettings<RelativeRootComparisonSettings>{
    private:
        typedef BaseCollectionSettings<RelativeRootComparisonSettings> BaseClass;

    public:
        RelativeRootCollectionSettings() : BaseClass() { }
        RelativeRootCollectionSettings(
                const ContinuousDistributionSettings& time_prior,
                unsigned int chain_length,
                unsigned int sample_frequency,
                const PositiveRealParameterSettings& concentration_settings,
                bool use_dpp)
                : BaseClass(
                        time_prior,
                        chain_length,
                        sample_frequency,
                        concentration_settings,
                        use_dpp) { }
        RelativeRootCollectionSettings(
                const std::string & yaml_config_path)
                : BaseClass(
                        yaml_config_path) { }
        RelativeRootCollectionSettings(
                std::istream& yaml_config_stream,
                const std::string& yaml_config_path)
                : BaseClass(
                        yaml_config_stream,
                        yaml_config_path) { }
};


class DirichletCollectionSettings: public BaseCollectionSettings<DirichletComparisonSettings>{
    private:
        typedef BaseCollectionSettings<DirichletComparisonSettings> BaseClass;

    public:
        DirichletCollectionSettings() : BaseClass(true) { }
        DirichletCollectionSettings(
                const ContinuousDistributionSettings& time_prior,
                unsigned int chain_length,
                unsigned int sample_frequency,
                const PositiveRealParameterSettings& concentration_settings,
                bool use_dpp)
                : BaseClass(
                        time_prior,
                        chain_length,
                        sample_frequency,
                        concentration_settings,
                        use_dpp,
                        true) {
        }
        DirichletCollectionSettings(
                const std::string & yaml_config_path)
                : BaseClass(
                        yaml_config_path,
                        true) {
        }
        DirichletCollectionSettings(
                std::istream& yaml_config_stream,
                const std::string& yaml_config_path)
                : BaseClass(
                        yaml_config_stream,
                        yaml_config_path,
                        true) {
        }

        bool using_population_size_multipliers() const {
            return true;
        }
};

#endif
