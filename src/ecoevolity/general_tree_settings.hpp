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
#include "settings.hpp"


class TreePriorSettings {
    public:
        TreePriorSettings() { }
		enum TreePriorEnum {
            uniform_root_and_betas     = 1,
        };
        
        virtual TreePriorSettings::TreePriorEnum get_type() const = 0;
};

class UniformTreePriorSettings : public TreePriorSettings {
    public:
        UniformTreePriorSettings() { }
};

class UniformRootAndBetasPriorSettings : public UniformTreePriorSettings {
    public:
        UniformRootAndBetasPriorSettings() { }
        UniformRootAndBetasPriorSettings(const UniformRootAndBetasPriorSettings & other) {
            this->root_height_settings_ = other.root_height_settings_;
            this->alpha_of_node_height_beta_settings_ = other.alpha_of_node_height_beta_settings_;
        }
        UniformRootAndBetasPriorSettings& operator=(const UniformRootAndBetasPriorSettings & other) {
            this->root_height_settings_ = other.root_height_settings_;
            this->alpha_of_node_height_beta_settings_ = other.alpha_of_node_height_beta_settings_;
            return * this;
        }
        UniformRootAndBetasPriorSettings(const YAML::Node& node) {
            this->parse_root_and_betas(node);
        }

        TreePriorSettings::TreePriorEnum get_type() const {
            return TreePriorSettings::TreePriorEnum::uniform_root_and_betas;
        }

        void parse_root_and_betas(const YAML::Node& node) {
            if (! node.IsMap()) {
                std::string message = (
                        "Expecting root_and_betas to be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
                throw EcoevolityYamlConfigError(message);
            }

            std::unordered_set<std::string> keys;
            for (YAML::const_iterator arg = node.begin();
                    arg != node.end();
                    ++arg) {
                if (keys.count(arg->first.as<std::string>()) > 0) {
                    std::string message = (
                            "Duplicate key in root_and_betas: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
                keys.insert(arg->first.as<std::string>());

                if (arg->first.as<std::string>() == "root_height") {
                    this->root_height_settings_ = PositiveRealParameterSettings(arg->second);
                }
                else if (arg->first.as<std::string>() == "alpha_of_node_height_beta_prior") {
                    this->alpha_of_node_height_beta_settings_ = PositiveRealParameterSettings(arg->second);
                }
                else {
                    std::string message = (
                            "Unrecognized root_and_betas setting: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

    protected:
        PositiveRealParameterSettings root_height_settings_;
        PositiveRealParameterSettings alpha_of_node_height_beta_settings_;
};


class TreeModelSettings {
    public:
		enum TreeSpaceEnum {
            generalized     = 1,
            bifurcating     = 2,
        };

        std::shared_ptr<TreePriorSettings> tree_prior = nullptr;

        TreeModelSettings() { }
        TreeModelSettings(const TreeModelSettings & other) {
            this->tree_space_ = other.tree_space_;
            this->tree_prior = other.tree_prior;
        }
        TreeModelSettings& operator=(const TreeModelSettings & other) {
            this->tree_space_ = other.tree_space_;
            this->tree_prior = other.tree_prior;
            return * this;
        }
        TreeModelSettings(const YAML::Node& node) {
            this->update_from_config(node);
        }

        TreeModelSettings::TreeSpaceEnum get_tree_space() const {
            return this->tree_space_;
        }

        void set_tree_space(const TreeModelSettings::TreeSpaceEnum tse) {
            this->tree_space_ = tse;
        }

        bool tree_space_generalized() const {
            return (this->get_tree_space() == TreeModelSettings::TreeSpaceEnum::generalized);
        }

        void update_from_config(const YAML::Node& node) {
            this->clear();
            if (! node.IsMap()) {
                std::string message = (
                        "Expecting tree_model node to be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
                throw EcoevolityYamlConfigError(message);
            }

            std::unordered_set<std::string> keys;
            for (YAML::const_iterator arg = node.begin();
                    arg != node.end();
                    ++arg) {
                if (keys.count(arg->first.as<std::string>()) > 0) {
                    std::string message = (
                            "Duplicate key in tree_model: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
                keys.insert(arg->first.as<std::string>());

                if (arg->first.as<std::string>() == "tree_space") {
                    std::string s = arg->second.as<std::string>();
                    if (s == "generalized") {
                        this->set_tree_space(TreeModelSettings::TreeSpaceEnum::generalized);
                    }
                    else if (s == "bifurcating") {
                        this->set_tree_space(TreeModelSettings::TreeSpaceEnum::bifurcating);
                    }
                    else {
                        std::string message = (
                                "Unrecognized tree_space setting: " + s);
                        throw EcoevolityYamlConfigError(message);
                    }
                }
                else if (arg->first.as<std::string>() == "tree_prior") {
                    this->parse_tree_prior(arg->second);
                }
                else {
                    std::string message = (
                            "Unrecognized tree_model setting: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

    protected:
        TreeModelSettings::TreeSpaceEnum tree_space_ = TreeModelSettings::TreeSpaceEnum::generalized;

        void clear() {
            this->tree_space_ = TreeModelSettings::TreeSpaceEnum::generalized;
            this->tree_prior = nullptr;
        }

        void parse_tree_prior(const YAML::Node& node) {
            if (! node.IsMap()) {
                std::string message = (
                        "Expecting tree_prior node to be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
                throw EcoevolityYamlConfigError(message);
            }

            std::unordered_set<std::string> keys;
            for (YAML::const_iterator arg = node.begin();
                    arg != node.end();
                    ++arg) {
                if (keys.count(arg->first.as<std::string>()) > 0) {
                    std::string message = (
                            "Duplicate key in tree_prior: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
                keys.insert(arg->first.as<std::string>());

                if (arg->first.as<std::string>() == "uniform") {
                    this->parse_uniform_node(arg->second);
                }
                else {
                    std::string message = (
                            "Unrecognized tree_prior setting: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

        void parse_uniform_node(const YAML::Node& node) {
            if (! node.IsMap()) {
                std::string message = (
                        "Expecting uniform tree_prior to be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
                throw EcoevolityYamlConfigError(message);
            }

            std::unordered_set<std::string> keys;
            for (YAML::const_iterator arg = node.begin();
                    arg != node.end();
                    ++arg) {
                if (keys.count(arg->first.as<std::string>()) > 0) {
                    std::string message = (
                            "Duplicate key in uniform tree_prior: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
                keys.insert(arg->first.as<std::string>());

                if (arg->first.as<std::string>() == "node_height_prior") {
                    this->parse_node_height_prior(arg->second);
                }
                else {
                    std::string message = (
                            "Unrecognized uniform tree_prior setting: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

        void parse_node_height_prior(const YAML::Node& node) {
            if (! node.IsMap()) {
                std::string message = (
                        "Expecting node_height_prior to be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
                throw EcoevolityYamlConfigError(message);
            }

            std::unordered_set<std::string> keys;
            for (YAML::const_iterator arg = node.begin();
                    arg != node.end();
                    ++arg) {
                if (keys.count(arg->first.as<std::string>()) > 0) {
                    std::string message = (
                            "Duplicate key in node_height_prior: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
                keys.insert(arg->first.as<std::string>());

                if (arg->first.as<std::string>() == "root_and_betas") {
                    this->tree_prior = std::make_shared<UniformRootAndBetasPriorSettings>(arg->second);
                }
                else {
                    std::string message = (
                            "Unrecognized node_height_prior setting: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }
};


class PopSizeSettings : public PositiveRealParameterSettings {
    protected:
        bool population_sizes_are_constrained_ = true;

    public:
        PopSizeSettings() { }
        PopSizeSettings(const YAML::Node& node) {
            YAML::Node n = node;
            if (n["equal_population_sizes"]) {
                this->population_sizes_are_constrained_ = n["equal_population_sizes"].as<bool>();
                n.remove("equal_population_sizes");
            }
            this->init_from_yaml_node(n);
        }
        PopSizeSettings(const PopSizeSettings & other) : PositiveRealParameterSettings(other) {
            this->population_sizes_are_constrained_ = other.population_sizes_are_constrained_;
        }
        PopSizeSettings& operator=(const PopSizeSettings& other) {
            this->value_ = other.value_;
            this->values_ = other.values_;
            this->is_fixed_ = other.is_fixed_;
            this->is_vector_ =  other.is_vector_;
            this->use_empirical_value_ = other.use_empirical_value_;
            this->prior_settings_ = other.prior_settings_;
            this->population_sizes_are_constrained_ = other.population_sizes_are_constrained_;
            return * this;
        }
        bool population_sizes_are_constrained() const {
            return this->population_sizes_are_constrained_;
        }
        void set_population_sizes_are_constrained(bool b) {
            this->population_sizes_are_constrained_ = b;
        }
        virtual std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            ss << std::boolalpha;
            std::string margin = string_util::get_indent(indent_level);
            if (this->population_sizes_are_constrained_) {
                ss << margin << "equal_population_sizes: true\n";
            }
            else {
                ss << margin << "equal_population_sizes: false\n";
            }
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

class GeneralTreeDataSettings {

        friend class GeneralTreeSettings;

    protected:
        std::string path_;
        char population_name_delimiter_ = ' ';
        bool population_name_is_prefix_ = true;
        bool genotypes_are_diploid_ = true;
        bool markers_are_dominant_ = false;
        bool constant_sites_removed_ = true;
        bool using_yaml_data_ = false;
        double ploidy_ = 2.0;

        void parse_yaml_data_settings(const YAML::Node& node,
                const std::string & config_path) {
            this->using_yaml_data_ = true;
            if (! node.IsMap()) {
                std::string message = (
                        "Expecting yaml_allele_counts node to be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
                throw EcoevolityYamlConfigError(message);
            }

            std::unordered_set<std::string> keys;
            for (YAML::const_iterator arg = node.begin();
                    arg != node.end();
                    ++arg) {
                if (keys.count(arg->first.as<std::string>()) > 0) {
                    std::string message = (
                            "Duplicate key in yaml_allele_counts: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
                keys.insert(arg->first.as<std::string>());

                if (arg->first.as<std::string>() == "path") {
                    std::string p = string_util::strip(arg->second.as<std::string>());
                    this->path_ = path::join(path::dirname(config_path), p);
                }
                else {
                    std::string message = (
                            "Unrecognized key in yaml_allele_counts: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

        void parse_alignment_settings(const YAML::Node& node,
                const std::string & config_path) {
            this->using_yaml_data_ = false;
            if (! node.IsMap()) {
                std::string message = (
                        "Expecting alignment node to be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
                throw EcoevolityYamlConfigError(message);
            }

            std::unordered_set<std::string> keys;
            for (YAML::const_iterator arg = node.begin();
                    arg != node.end();
                    ++arg) {
                if (keys.count(arg->first.as<std::string>()) > 0) {
                    std::string message = (
                            "Duplicate key in alignment parameters: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
                keys.insert(arg->first.as<std::string>());

                if (arg->first.as<std::string>() == "genotypes_are_diploid") {
                    this->genotypes_are_diploid_ = arg->second.as<bool>();
                }
                else if (arg->first.as<std::string>() == "markers_are_dominant") {
                    this->markers_are_dominant_ = arg->second.as<bool>();
                }
                else if (arg->first.as<std::string>() == "population_name_is_prefix") {
                    this->population_name_is_prefix_ = arg->second.as<bool>();
                }
                else if (arg->first.as<std::string>() == "population_name_delimiter") {
                    this->population_name_delimiter_ = arg->second.as<char>();
                }
                else if (arg->first.as<std::string>() == "path") {
                    std::string p = string_util::strip(arg->second.as<std::string>());
                    this->path_ = path::join(path::dirname(config_path), p);
                }
                else {
                    std::string message = (
                            "Unrecognized key in alignment parameters: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

    public:
        GeneralTreeDataSettings() { }
        GeneralTreeDataSettings(const GeneralTreeDataSettings & other) {
            this->path_ = other.path_;
            this->population_name_delimiter_ = other.population_name_delimiter_;
            this->population_name_is_prefix_ = other.population_name_is_prefix_;
            this->genotypes_are_diploid_ = other.genotypes_are_diploid_;
            this->markers_are_dominant_ = other.markers_are_dominant_;
            this->constant_sites_removed_ = other.constant_sites_removed_;
            this->using_yaml_data_ = other.using_yaml_data_;
            this->ploidy_ = other.ploidy_;
        }
        GeneralTreeDataSettings& operator=(const GeneralTreeDataSettings& other) {
            this->path_ = other.path_;
            this->population_name_delimiter_ = other.population_name_delimiter_;
            this->population_name_is_prefix_ = other.population_name_is_prefix_;
            this->genotypes_are_diploid_ = other.genotypes_are_diploid_;
            this->markers_are_dominant_ = other.markers_are_dominant_;
            this->constant_sites_removed_ = other.constant_sites_removed_;
            this->using_yaml_data_ = other.using_yaml_data_;
            this->ploidy_ = other.ploidy_;
            return * this;
        }
        virtual ~GeneralTreeDataSettings() { }
        const std::string& get_path() const {
            return this->path_;
        }
        void set_ploidy(double p) {
            this->ploidy_ = p;
        }
        double get_ploidy() const {
            return this->ploidy_;
        }
        void update_from_config(const YAML::Node& node,
                const std::string & config_path) {
            if (! node.IsMap()) {
                std::string message = (
                        "Expecting data node to be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
                throw EcoevolityYamlConfigError(message);
            }

            if ((node["alignment"]) && (node["yaml_allele_counts"])) {
                std::string message(
                        "Data specified as both 'alignment' and 'yaml_allele_counts'");
                throw EcoevolityYamlConfigError(message);
            }
            std::unordered_set<std::string> keys;
            for (YAML::const_iterator p = node.begin();
                    p != node.end();
                    ++p) {
                if (keys.count(p->first.as<std::string>()) > 0) {
                    std::string message = (
                            "Duplicate key in data parameters: " +
                            p->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
                keys.insert(p->first.as<std::string>());

                if (p->first.as<std::string>() == "ploidy") {
                    this->set_ploidy(p->second.as<double>());
                }
                else if (p->first.as<std::string>() == "constant_sites_removed") {
                    this->constant_sites_removed_ = p->second.as<bool>();
                }
                else if (p->first.as<std::string>() == "alignment") {
                    this->parse_alignment_settings(p->second, config_path);
                }
                else if (p->first.as<std::string>() == "yaml_allele_counts") {
                    this->parse_yaml_data_settings(p->second, config_path);
                }
                else {
                    std::string message = (
                            "Unrecognized key in data parameters: " +
                            p->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

        std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            std::string margin = string_util::get_indent(indent_level);
            std::string indent = string_util::get_indent(1);
            ss << margin << "ploidy: " << this->ploidy_ << "\n";
            ss << margin << "constant_sites_removed: " << this->constant_sites_removed_ << "\n";
            if (this->using_yaml_data_) {
                ss << margin << "yaml_allele_counts:\n";
                ss << margin << indent << "path: " << this->path_ << "\n";
            }
            else {
                ss << margin << "alignment:\n";
                ss << margin << indent << "path: " << this->path_ << "\n";
                ss << margin << indent << "genotypes_are_diploid: " << this->genotypes_are_diploid_ << "\n";
                ss << margin << indent << "markers_are_dominant: " << this->markers_are_dominant_ << "\n";
                ss << margin << indent << "population_name_delimiter: '" << this->population_name_delimiter_ << "'\n";
                ss << margin << indent << "population_name_is_prefix: " << this->population_name_is_prefix_ << "\n";
            }

            return ss.str();
        }
};

class GeneralTreeOperatorSettings {
    protected:
        double weight_ = 1.0;

    public:
        GeneralTreeOperatorSettings() { }
        GeneralTreeOperatorSettings(double weight) {
            this->weight_ = weight;
        }
        GeneralTreeOperatorSettings(const GeneralTreeOperatorSettings& other) {
            this->weight_ = other.weight_;
        }
        virtual ~GeneralTreeOperatorSettings() { }
        GeneralTreeOperatorSettings& operator=(const GeneralTreeOperatorSettings& other) {
            this->weight_ = other.weight_;
            return * this;
        }
        virtual bool is_tunable() const {
            return false;
        }
        double get_weight() const {
            return this->weight_;
        }
        void set_weight(double weight) {
            this->weight_ = weight;
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
            std::unordered_set<std::string> keys;
            for (YAML::const_iterator p = operator_node.begin();
                    p != operator_node.end();
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

class GeneralTreeTunableOperatorSettings : public GeneralTreeOperatorSettings {
    protected:
        double tuning_parameter_ = 1.0;
        bool auto_optimize_ = true;
        unsigned int auto_optimize_delay_ = 1000;

    public:
        GeneralTreeTunableOperatorSettings() : GeneralTreeOperatorSettings() { }
        GeneralTreeTunableOperatorSettings(double weight) : GeneralTreeOperatorSettings(weight) { }
        GeneralTreeTunableOperatorSettings(const GeneralTreeTunableOperatorSettings& other) :
                GeneralTreeOperatorSettings(other) {
            this->tuning_parameter_ = other.tuning_parameter_;
            this->auto_optimize_ = other.auto_optimize_;
            this->auto_optimize_delay_ = other.auto_optimize_delay_;
        }
        virtual ~GeneralTreeTunableOperatorSettings() { }
        GeneralTreeTunableOperatorSettings& operator=(const GeneralTreeTunableOperatorSettings& other) {
            this->weight_ = other.weight_;
            this->tuning_parameter_ = other.tuning_parameter_;
            this->auto_optimize_ = other.auto_optimize_;
            this->auto_optimize_delay_ = other.auto_optimize_delay_;
            return * this;
        }
        bool is_tunable() const {
            return true;
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
        void update_from_config(const YAML::Node& operator_node) {
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

            std::unordered_set<std::string> keys;
            for (YAML::const_iterator p = operator_node.begin();
                    p != operator_node.end();
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

        std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            std::string margin = string_util::get_indent(indent_level);
            ss << margin << "weight: " << this->weight_ << "\n";
            ss << margin << "tuning_parameter: " << this->tuning_parameter_ << "\n";
            ss << margin << "auto_optimize: " << this->auto_optimize_ << "\n";
            ss << margin << "auto_optimize_delay: " << this->auto_optimize_delay_ << "\n";
            return ss.str();
        }
};

class GeneralTreeOperatorSettingsCollection {
    protected:
        std::map<std::string, GeneralTreeOperatorSettings> untunable_operators_ {
            {"SplitLumpNodesRevJumpSampler",        GeneralTreeOperatorSettings(10)},
            {"NeighborHeightNodePermute",           GeneralTreeOperatorSettings(6)},
            {"NeighborHeightNodeSwap",              GeneralTreeOperatorSettings(3)},
        };
        std::map<std::string, GeneralTreeTunableOperatorSettings> tunable_operators_ {
            {"MuRateScaler",                        GeneralTreeTunableOperatorSettings(1)},
            {"NodeHeightPriorAlphaScaler",          GeneralTreeTunableOperatorSettings(1)},
            {"NodeHeightPriorAlphaMover",           GeneralTreeTunableOperatorSettings(0)},
            {"TreeScaler",                          GeneralTreeTunableOperatorSettings(1)},
            {"NodeHeightScaler",                    GeneralTreeTunableOperatorSettings(1)},
            {"NodeHeightMover",                     GeneralTreeTunableOperatorSettings(1)},
            {"NodeHeightSlideBumpScaler",           GeneralTreeTunableOperatorSettings(1)},
            {"NodeHeightSlideBumpPermuteScaler",    GeneralTreeTunableOperatorSettings(0)},
            {"NodeHeightSlideBumpSwapScaler",       GeneralTreeTunableOperatorSettings(1)},
            {"NodeHeightSlideBumpMover",            GeneralTreeTunableOperatorSettings(0)},
            {"NodeHeightSlideBumpPermuteMover",     GeneralTreeTunableOperatorSettings(0)},
            {"NodeHeightSlideBumpSwapMover",        GeneralTreeTunableOperatorSettings(0)},
            {"RootHeightScaler",                    GeneralTreeTunableOperatorSettings(1)},
            {"GlobalNodeHeightDirichletOperator",   GeneralTreeTunableOperatorSettings(1)},
            {"NodeHeightDirichletOperator",         GeneralTreeTunableOperatorSettings(1)},
            {"GlobalPopSizeScaler",                 GeneralTreeTunableOperatorSettings(1)},
            {"PopSizeScaler",                       GeneralTreeTunableOperatorSettings(1)},
            {"GlobalHeightSizeMixer",               GeneralTreeTunableOperatorSettings(1)},
            {"HeightSizeMixer",                     GeneralTreeTunableOperatorSettings(1)},
            {"HeightSizeSlideBumpMixer",            GeneralTreeTunableOperatorSettings(0)},
            {"RootHeightSizeMixer",                 GeneralTreeTunableOperatorSettings(1)},
            {"GlobalHeightSizeRateScaler",          GeneralTreeTunableOperatorSettings(1)},
            {"GlobalHeightSizeScaler",              GeneralTreeTunableOperatorSettings(0)},
            {"GlobalHeightRateScaler",              GeneralTreeTunableOperatorSettings(0)},
            {"StateFreqMover",                      GeneralTreeTunableOperatorSettings(1)},
            {"StateFreqDirichletOperator",          GeneralTreeTunableOperatorSettings(1)},
        };

    public:
        GeneralTreeOperatorSettingsCollection() { }
        GeneralTreeOperatorSettingsCollection& operator=(const GeneralTreeOperatorSettingsCollection& other) {
            this->untunable_operators_ = other.untunable_operators_;
            this->tunable_operators_ = other.tunable_operators_;
            return * this;
        }
        GeneralTreeOperatorSettingsCollection(const GeneralTreeOperatorSettingsCollection& other) {
            this->untunable_operators_ = other.untunable_operators_;
            this->tunable_operators_ = other.tunable_operators_;
        }
        bool is_operator(const std::string & name) const {
            if (this->untunable_operators_.count(name) > 0) {
                return true;
            }
            if (this->tunable_operators_.count(name) > 0) {
                return true;
            }
            return false;
        }
        bool is_tunable_operator(const std::string & name) const {
            if (this->tunable_operators_.count(name) > 0) {
                return true;
            }
            return false;
        }
        bool is_untunable_operator(const std::string & name) const {
            if (this->untunable_operators_.count(name) > 0) {
                return true;
            }
            return false;
        }

        void update_from_config(const YAML::Node& operator_node) {
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

            std::unordered_set<std::string> keys;
            for (YAML::const_iterator p = operator_node.begin();
                    p != operator_node.end();
                    ++p) {
                std::string name = p->first.as<std::string>();
                if (keys.count(name) > 0) {
                    std::string message = (
                            "Duplicate key in operator parameters: " +
                            name);
                    throw EcoevolityYamlConfigError(message);
                }
                keys.insert(name);

                if (this->is_untunable_operator(name)) {
                    this->untunable_operators_[name].update_from_config(p->second);
                }
                else if (this->is_tunable_operator(name)) {
                    this->tunable_operators_[name].update_from_config(p->second);
                }
                else {
                    std::string message = (
                            "Unrecognized operator: " +
                            name);
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

        std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            std::string margin = string_util::get_indent(indent_level);
            for (auto op : this->untunable_operators_) {
                ss << op.first << ":\n";
                ss << op.second.to_string(indent_level + 1);
            }
            for (auto op : this->tunable_operators_) {
                ss << op.first << ":\n";
                ss << op.second.to_string(indent_level + 1);
            }
            return ss.str();
        }
};


class GeneralTreeSettings {
    public:
        // Public data members
        PopSizeSettings population_size_settings;
        PositiveRealParameterSettings freq_1_settings;
        PositiveRealParameterSettings mutation_rate_settings;
        GeneralTreeOperatorSettingsCollection operator_settings;
        GeneralTreeDataSettings data_settings;
        TreeModelSettings tree_model_settings;

        // Constructors
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
        GeneralTreeSettings(const GeneralTreeSettings & other) {
            this->config_path_ = other.config_path_;
            this->chain_length_ = other.chain_length_;
            this->sample_frequency_ = other.sample_frequency_;
            this->tree_log_path_ = other.tree_log_path_;
            this->state_log_path_ = other.state_log_path_;
            this->operator_log_path_ = other.operator_log_path_;
            this->population_size_settings = other.population_size_settings;
            this->freq_1_settings = other.freq_1_settings;
            this->mutation_rate_settings = other.mutation_rate_settings;
            this->operator_settings = other.operator_settings;
            this->data_settings = other.data_settings;
        }
        virtual ~GeneralTreeSettings() { }
        GeneralTreeSettings & operator=(const GeneralTreeSettings & other) {
            this->config_path_ = other.config_path_;
            this->chain_length_ = other.chain_length_;
            this->sample_frequency_ = other.sample_frequency_;
            this->tree_log_path_ = other.tree_log_path_;
            this->state_log_path_ = other.state_log_path_;
            this->operator_log_path_ = other.operator_log_path_;
            this->population_size_settings = other.population_size_settings;
            this->freq_1_settings = other.freq_1_settings;
            this->mutation_rate_settings = other.mutation_rate_settings;
            this->operator_settings = other.operator_settings;
            this->data_settings = other.data_settings;
            return * this;
        }

        const std::string& get_config_path() const {
            return this->config_path_;
        }
        std::string get_config_directory() const {
            return path::dirname(this->get_config_path());
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

            out << "---\n";
        }

        std::string to_string() const {
            std::ostringstream ss;
            this->write_settings(ss);
            return ss.str();
        }

    protected:
        std::string config_path_ = "";
        unsigned int chain_length_ = 100000;
        unsigned int sample_frequency_ = 100;
        std::string tree_log_path_ = "phycoevolity-tree-run-1.log";
        std::string state_log_path_ = "phycoevolity-state-run-1.log";
        std::string operator_log_path_ = "phycoevolity-operator-run-1.log";

        void set_output_paths_to_config_directory() {
            std::pair<std::string, std::string> prefix_ext = path::splitext(this->config_path_);
            this->tree_log_path_ = prefix_ext.first + "-tree-run-1.log";
            this->state_log_path_ = prefix_ext.first + "-state-run-1.log";
            this->operator_log_path_ = prefix_ext.first + "-operator-run-1.log";
        }

        void init_from_config_stream(std::istream& stream, const std::string& path) {
            this->config_path_ = path;
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

                // parse data settings
                if (top->first.as<std::string>() == "data") {
                    this->data_settings.update_from_config(top->second, this->config_path_);
                }
                // parse mcmc settings
                else if (top->first.as<std::string>() == "mcmc_settings") {
                    this->parse_mcmc_settings(top->second);
                }
                // parse tree model 
                else if (top->first.as<std::string>() == "tree_model") {
                    this->tree_model_settings.update_from_config(top->second);
                }
                // parse branch parameters
                else if (top->first.as<std::string>() == "branch_parameters") {
                    this->parse_branch_parameters(top->second);
                }
                // parse mutation parameters
                else if (top->first.as<std::string>() == "mutation_parameters") {
                    this->parse_mutation_parameters(top->second);
                }
                else {
                    throw EcoevolityYamlConfigError(
                            "Unrecognized top level key: '" +
                            top->first.as<std::string>() + "'");
                }
            }
        }

        void parse_branch_parameters(const YAML::Node& node) {
            if (! node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "Expecting branch_parameters node to be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
            }
            std::unordered_set<std::string> keys;
            for (YAML::const_iterator sub_node = node.begin();
                    sub_node != node.end();
                    ++sub_node) {
                if (keys.count(sub_node->first.as<std::string>()) > 0) {
                    std::string message = (
                            "Duplicate key in branch_parameters: " +
                            sub_node->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
                keys.insert(sub_node->first.as<std::string>());

                if (sub_node->first.as<std::string>() == "population_size") {
                    this->population_size_settings = PopSizeSettings(sub_node->second);
                    if (this->population_size_settings.use_empirical_value()) {
                        throw EcoevolityPositiveRealParameterSettingError(
                                "empirical value not supported for population_size");
                    }
                }
                else {
                    throw EcoevolityYamlConfigError(
                            "Unrecognized key in branch_parameters: '" +
                            sub_node->first.as<std::string>() + "'");
                }
            }
        }

        void parse_mutation_parameters(const YAML::Node& node) {
            if (! node.IsMap()) {
                throw EcoevolityYamlConfigError(
                        "Expecting mutation_parameters node to be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
            }
            std::unordered_set<std::string> keys;
            for (YAML::const_iterator sub_node = node.begin();
                    sub_node != node.end();
                    ++sub_node) {
                if (keys.count(sub_node->first.as<std::string>()) > 0) {
                    std::string message = (
                            "Duplicate key in mutation_parameters: " +
                            sub_node->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
                keys.insert(sub_node->first.as<std::string>());

                if (sub_node->first.as<std::string>() == "freq_1") {
                    this->freq_1_settings = PositiveRealParameterSettings(sub_node->second);
                }
                else if (sub_node->first.as<std::string>() == "mutation_rate") {
                    this->mutation_rate_settings = PositiveRealParameterSettings(sub_node->second);
                    if (this->mutation_rate_settings.use_empirical_value()) {
                        throw EcoevolityPositiveRealParameterSettingError(
                                "empirical value not supported for mutation_rate");
                    }
                }
                else {
                    throw EcoevolityYamlConfigError(
                            "Unrecognized key in mutation_parameters: '" +
                            sub_node->first.as<std::string>() + "'");
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
                // parse operator settings
                else if (mcmc->first.as<std::string>() == "operators") {
                    this->operator_settings.update_from_config(mcmc->second);
                }
                else {
                    throw EcoevolityYamlConfigError(
                            "Unrecognized mcmc_settings key: " +
                            mcmc->first.as<std::string>());
                }
            }
        }

        virtual void update_settings_contingent_upon_data() {
            if (this->freq_1_settings.use_empirical_value()) {
                BiallelicData d = BiallelicData(
                        this->data_settings.path_,
                        this->data_settings.population_name_delimiter_,
                        this->data_settings.population_name_is_prefix_,
                        this->data_settings.genotypes_are_diploid_,
                        this->data_settings.markers_are_dominant_,
                        true);
                this->freq_1_settings.value_ = d.get_proportion_1();
            }
        }
};

#endif
