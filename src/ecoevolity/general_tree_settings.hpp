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

#ifndef GENERAL_TREE_SETTINGS_HPP
#define GENERAL_TREE_SETTINGS_HPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <unordered_map>
#include <unordered_set>

#include "yaml_util.hpp"

#include "error.hpp"
#include "assert.hpp"
#include "string_util.hpp"
#include "path.hpp"
#include "math_util.hpp"
#include "probability.hpp"
#include "data.hpp"
#include "options.hpp"


class GeneralTreeOperatorSettings {
    protected:
        std::string name_ = "none";
        double weight_ = 1.0;
        double tuning_parameter_ = 1.0;
        bool auto_optimize_ = true;
        unsigned int auto_optimize_delay_ = 1000;

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
        const std::string & get_name() const {
            return this->name_;
        }
        void set_name(const std::string & name) {
            this->name_ = name;
        }
        double get_weight() const {
            return this->weight_;
        }
        void set_weight(double weight) {
            this->weight_ = weight;
        }
        double get_tuning_parameter() const {
            return this->tuning_parameter_;
        }
        void set_tuning_parameter(double tuning_parameter) {
            this->tuning_parameter_ = tuning_parameter;
        }
        unsigned int get_auto_optimize_delay() const {
            return this->auto_optimize_delay_;
        }
        void set_auto_optimize_delay(unsigned int auto_optimize_delay) {
            this->auto_optimize_delay_ = auto_optimize_delay;
        }
        void set_auto_optimize(bool b) {
            this->auto_optimize_ = b;
        }
        bool auto_optimizing() const {
            return this->auto_optimize_;
        }
        virtual void update_from_config(const YAML::Node& operator_node) {
            if (! operator_node.IsMap()) {
                std::string message = (
                        "Expecting operator parameters to be a map, but found: " +
                        YamlCppUtils::get_node_type(operator_node));
                throw EcoevolityYamlConfigError(message);
            }
            if (operator_node.size() != 1) {
                throw EcoevolityYamlConfigError(
                        "operator node should only have a single key");
            }
            this->name_ = operator_node.begin()->first.as<std::string>();
            YAML::Node parameters = operator_node[this->name_];

            std::unordered_set<std::string> keys;
            for (YAML::const_iterator p = parameters.begin();
                    p != parameters.end();
                    ++p) {
                if (keys.count(p->first.as<std::string>()) > 0) {
                    std::string message = (
                            "Duplicate key in operator parameters: " +
                            p->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
                keys.insert(p->first.as<std::string>());

                if (p->first.as<std::string>() == "weight") {
                    this->set_weight(p->second.as<double>());
                }
                else if (p->first.as<std::string>() == "tuning_parameter") {
                    this->set_tuning_parameter(p->second.as<double>());
                }
                else if (p->first.as<std::string>() == "auto_optimize") {
                    this->set_tuning_parameter(p->second.as<bool>());
                }
                else if (p->first.as<std::string>() == "auto_optimize_delay") {
                    this->set_tuning_parameter(p->second.as<unsigned int>());
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
            std::string indent = string_util::get_indent(1);
            ss << margin << this->name_ << ":\n";
            ss << margin << indent << "weight: " << this->weight_ << "\n";
            ss << margin << indent << "tuning_parameter_: " << this->tuning_parameter_ << "\n";
            ss << margin << indent << "auto_optimize: " << this->auto_optimize_ << "\n";
            ss << margin << indent << "auto_optimize_delay: " << this->auto_optimize_delay_ << "\n";
            return ss.str();
        }
};


class GeneralTreeSettings {
    public:
        GeneralTreeSettings() { }
        GeneralTreeSettings(
                const std::string & yaml_config_path) {
            this->init_from_config_file(yaml_config_path);
        }
        GeneralTreeSettings(
                std::istream& yaml_config_stream,
                const std::string & yaml_config_path) {
            this->init_from_config_stream(yaml_config_stream, yaml_config_path);
        }
        virtual ~BaseCollectionSettings() { }
        GeneralTreeSettings & operator=(const GeneralTreeSettings & other) {
            this->path_ = other.path_;
            this->chain_length_ = other.chain_length_;
            this->sample_frequency_ = other.sample_frequency_;
            this->tree_log_path_ = other.tree_log_path_;
            this->state_log_path_ = other.state_log_path_;
            this->operator_log_path_ = other.operator_log_path_;
            return * this;
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
        double get_chain_length() const {
            return this->chain_length_;
        }
        double get_sample_frequency() const {
            return this->sample_frequency_;
        }
        void write_settings(std::ostream& out) const {
            std::string indent = string_util::get_indent(1);
            out << std::boolalpha;

            out << "---\n"
        }

        std::string to_string() const {
            std::ostringstream ss;
            this->write_settings(ss);
            return ss.str();
        }

    protected:
        std::string path_ = "";
        unsigned int chain_length_ = 100000;
        unsigned int sample_frequency_ = 100;
        std::string tree_log_path_ = "phycoevolity-tree-run-1.log";
        std::string state_log_path_ = "phycoevolity-state-run-1.log";
        std::string operator_log_path_ = "phycoevolity-operator-run-1.log";

        void set_output_paths_to_config_directory() {
            std::pair<std::string, std::string> prefix_ext = path::splitext(this->path_);
            this->tree_log_path_ = prefix_ext.first + "-tree-run-1.log";
            this->state_log_path_ = prefix_ext.first + "-state-run-1.log";
            this->operator_log_path_ = prefix_ext.first + "-operator-run-1.log";
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

            if (! top_level_node["data"]) {
                throw EcoevolityYamlConfigError("No data");
            }
            this->parse_data_settings(top_level_node["data"]);


            std::unordered_set<std::string> keys;
            for (YAML::const_iterator top = top_level_node.begin();
                    top != top_level_node.end();
                    ++top) {
                if (keys.count(top->first.as<std::string>()) > 0) {
                    std::string message = (
                            "Duplicate top-level key: " +
                            top->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
                keys.insert(top->first.as<std::string>());

                // parse mcmc settings
                if (top->first.as<std::string>() == "mcmc_settings") {
                    this->parse_mcmc_settings(top->second);
                }
                // parse operator settings
                else if (top->first.as<std::string>() == "operators") {
                    this->operators_.update_from_config(top->second);
                }
                else if (top->first.as<std::string>() == "data") {
                    // Already parsed, nothing to do here
                    continue;
                }
                else {
                    throw EcoevolityYamlConfigError(
                            "Unrecognized top level key: '" +
                            top->first.as<std::string>() + "'");
                }
            }
        }

        void parse_mcmc_settings(const YAML::Node& mcmc_node) {
            if (! mcmc_node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "Expecting mcmc_settings to be a map, but found: " +
                        YamlCppUtils::get_node_type(mcmc_node));
            }
            std::unordered_set<std::string> keys;
            for (YAML::const_iterator mcmc = mcmc_node.begin(); mcmc != mcmc_node.end(); ++mcmc) {
                if (keys.count(mcmc->first.as<std::string>()) > 0) {
                    std::string message = (
                            "Duplicate mcmc settings key: " +
                            mcmc->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
                keys.insert(mcmc->first.as<std::string>());

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
};

#endif
