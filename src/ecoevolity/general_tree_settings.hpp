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


class StartingTreeSettings {
    public:

        StartingTreeSettings() { }

        enum Source {
            option = 1,
            string = 2,
            path = 3,
        };

		enum Option {
            comb = 1,
            random = 2,
            upgma = 3,
        };

    protected:

        std::string path_ = "";
        std::string string_ = "";
        StartingTreeSettings::Option option_ = StartingTreeSettings::Option::random;
        StartingTreeSettings::Source source_ = StartingTreeSettings::Source::option;

    public:

        StartingTreeSettings::Source get_tree_source() const {
            return this->source_;
        }
        StartingTreeSettings::Option get_tree_option() const {
            return this->option_;
        }
        const std::string& get_tree_path() const {
            return this->path_;
        }
        const std::string& get_tree_string() const {
            return this->string_;
        }

        void update_from_config(const YAML::Node& node,
                const std::string & config_path) {
            if (! node.IsMap()) {
                std::string message = (
                        "Expecting starting_tree to be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
                throw EcoevolityYamlConfigError(message);
            }
            if (node.size() != 1) {
                throw EcoevolityYamlConfigError(
                        "starting_tree node should only have a single key");
            }
            std::unordered_set<std::string> keys;
            for (YAML::const_iterator arg = node.begin();
                    arg != node.end();
                    ++arg) {
                if (keys.count(arg->first.as<std::string>()) > 0) {
                    std::string message = (
                            "Duplicate key in starting_tree: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
                keys.insert(arg->first.as<std::string>());

                if (arg->first.as<std::string>() == "option") {
                    this->source_ = StartingTreeSettings::Source::option;
                    std::string opt = arg->second.as<std::string>();
                    if (opt == "comb") {
                        this->option_ = StartingTreeSettings::Option::comb;
                    }
                    else if (opt == "random") {
                        this->option_ = StartingTreeSettings::Option::random;
                    }
                    else if (opt == "upgma") {
                        this->option_ = StartingTreeSettings::Option::upgma;
                        throw EcoevolityYamlConfigError(
                                "UPGMA starting tree not implemented");
                    }
                    else {
                        throw EcoevolityYamlConfigError(
                                "Unrecognized starting_tree option: " + opt);
                    }
                }
                else if (arg->first.as<std::string>() == "string") {
                    this->source_ = StartingTreeSettings::Source::string;
                    this->string_ = arg->second.as<std::string>();
                }
                else if (arg->first.as<std::string>() == "path") {
                    this->source_ = StartingTreeSettings::Source::path;
                    std::string p = string_util::strip(arg->second.as<std::string>());
                    this->path_ = path::join(path::dirname(config_path), p);
                    if (! path::isfile(this->path_)) {
                        std::string message = (
                                "starting tree path \'" + this->path_ + "\' is not a file");
                        throw EcoevolityYamlConfigError(message);
                    }
                }
                else {
                    std::string message = (
                            "Unrecognized key in starting_tree: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

        std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            ss << std::boolalpha;
            std::string margin = string_util::get_indent(indent_level);
            if (this->source_ == StartingTreeSettings::Source::option) {
                if (this->option_ == StartingTreeSettings::Option::random) {
                    ss << margin << "option: random\n";
                }
                else if (this->option_ == StartingTreeSettings::Option::comb) {
                    ss << margin << "option: comb\n";
                }
                else if (this->option_ == StartingTreeSettings::Option::upgma) {
                    ss << margin << "option: upgma\n";
                }
                else {
                    throw EcoevolityYamlConfigError(
                            "Unexpected StartingTreeSettings::option_");
                }
            }
            else if (this->source_ == StartingTreeSettings::Source::string) {
                ss << margin << "string: " << this->string_ << "\n";
            }
            else if (this->source_ == StartingTreeSettings::Source::path) {
                ss << margin << "path: " << this->path_ << "\n";
            }
            else {
                throw EcoevolityYamlConfigError(
                        "Unexpected StartingTreeSettings::source_");
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
            ss << std::boolalpha;
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
            ss << std::boolalpha;
            std::string margin = string_util::get_indent(indent_level);
            ss << margin << "weight: " << this->weight_ << "\n";
            ss << margin << "tuning_parameter: " << this->tuning_parameter_ << "\n";
            ss << margin << "auto_optimize: " << this->auto_optimize_ << "\n";
            ss << margin << "auto_optimize_delay: " << this->auto_optimize_delay_ << "\n";
            return ss.str();
        }
};


class GeneralTreeOperatorSettingsCollection {
    public:
        std::map<std::string, GeneralTreeOperatorSettings> untunable_operators {
            {"SplitLumpNodesRevJumpSampler",        GeneralTreeOperatorSettings(10)},
            {"NeighborHeightNodePermute",           GeneralTreeOperatorSettings(6)},
            {"NeighborHeightNodeSwap",              GeneralTreeOperatorSettings(3)},
        };
        std::map<std::string, GeneralTreeTunableOperatorSettings> tunable_operators {
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
        };

        GeneralTreeOperatorSettingsCollection() { }
        GeneralTreeOperatorSettingsCollection& operator=(const GeneralTreeOperatorSettingsCollection& other) {
            this->untunable_operators = other.untunable_operators;
            this->tunable_operators = other.tunable_operators;
            return * this;
        }
        GeneralTreeOperatorSettingsCollection(const GeneralTreeOperatorSettingsCollection& other) {
            this->untunable_operators = other.untunable_operators;
            this->tunable_operators = other.tunable_operators;
        }
        bool is_operator(const std::string & name) const {
            if (this->untunable_operators.count(name) > 0) {
                return true;
            }
            if (this->tunable_operators.count(name) > 0) {
                return true;
            }
            return false;
        }
        bool is_tunable_operator(const std::string & name) const {
            if (this->tunable_operators.count(name) > 0) {
                return true;
            }
            return false;
        }
        bool is_untunable_operator(const std::string & name) const {
            if (this->untunable_operators.count(name) > 0) {
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
                    this->untunable_operators[name].update_from_config(p->second);
                }
                else if (this->is_tunable_operator(name)) {
                    this->tunable_operators[name].update_from_config(p->second);
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
            ss << std::boolalpha;
            std::string margin = string_util::get_indent(indent_level);
            for (auto op : this->untunable_operators) {
                ss << margin << op.first << ":\n";
                ss << op.second.to_string(indent_level + 1);
            }
            for (auto op : this->tunable_operators) {
                ss << margin << op.first << ":\n";
                ss << op.second.to_string(indent_level + 1);
            }
            return ss.str();
        }
};

class BetaTreeOperatorSettingsCollection : public GeneralTreeOperatorSettingsCollection {
    public:
        BetaTreeOperatorSettingsCollection() : GeneralTreeOperatorSettingsCollection() {
            this->tunable_operators["NodeHeightPriorAlphaScaler"]  = GeneralTreeTunableOperatorSettings(1);
            this->tunable_operators["NodeHeightPriorAlphaMover"]   = GeneralTreeTunableOperatorSettings(0);
        }
        BetaTreeOperatorSettingsCollection& operator=(const BetaTreeOperatorSettingsCollection& other) {
            this->untunable_operators = other.untunable_operators;
            this->tunable_operators = other.tunable_operators;
            return * this;
        }
        BetaTreeOperatorSettingsCollection(const BetaTreeOperatorSettingsCollection& other) :
            GeneralTreeOperatorSettingsCollection(other) { }
};

class PopulationTreeOperatorSettingsCollection : public GeneralTreeOperatorSettingsCollection {

    public:
        PopulationTreeOperatorSettingsCollection() : GeneralTreeOperatorSettingsCollection() {
            this->tunable_operators["MuRateScaler"]                = GeneralTreeTunableOperatorSettings(1);
            this->tunable_operators["GlobalPopSizeScaler"]         = GeneralTreeTunableOperatorSettings(1);
            this->tunable_operators["PopSizeScaler"]               = GeneralTreeTunableOperatorSettings(1);
            this->tunable_operators["GlobalHeightSizeMixer"]       = GeneralTreeTunableOperatorSettings(1);
            this->tunable_operators["HeightSizeMixer"]             = GeneralTreeTunableOperatorSettings(1);
            this->tunable_operators["HeightSizeSlideBumpMixer"]    = GeneralTreeTunableOperatorSettings(0);
            this->tunable_operators["RootHeightSizeMixer"]         = GeneralTreeTunableOperatorSettings(1);
            this->tunable_operators["GlobalHeightSizeRateScaler"]  = GeneralTreeTunableOperatorSettings(1);
            this->tunable_operators["GlobalHeightSizeScaler"]      = GeneralTreeTunableOperatorSettings(0);
            this->tunable_operators["GlobalHeightRateScaler"]      = GeneralTreeTunableOperatorSettings(0);
            this->tunable_operators["StateFreqMover"]              = GeneralTreeTunableOperatorSettings(1);
            this->tunable_operators["StateFreqDirichletOperator"]  = GeneralTreeTunableOperatorSettings(1);
        }
        PopulationTreeOperatorSettingsCollection& operator=(const PopulationTreeOperatorSettingsCollection& other) {
            this->untunable_operators = other.untunable_operators;
            this->tunable_operators = other.tunable_operators;
            return * this;
        }
        PopulationTreeOperatorSettingsCollection(const PopulationTreeOperatorSettingsCollection& other) :
            GeneralTreeOperatorSettingsCollection(other) { }

};

class BetaPopulationTreeOperatorSettingsCollection : public PopulationTreeOperatorSettingsCollection {

    public:
        BetaPopulationTreeOperatorSettingsCollection() : PopulationTreeOperatorSettingsCollection() {
            this->tunable_operators["NodeHeightPriorAlphaScaler"]  = GeneralTreeTunableOperatorSettings(1);
            this->tunable_operators["NodeHeightPriorAlphaMover"]   = GeneralTreeTunableOperatorSettings(0);
        }
        BetaPopulationTreeOperatorSettingsCollection& operator=(const BetaPopulationTreeOperatorSettingsCollection& other) {
            this->untunable_operators = other.untunable_operators;
            this->tunable_operators = other.tunable_operators;
            return * this;
        }
        BetaPopulationTreeOperatorSettingsCollection(const BetaPopulationTreeOperatorSettingsCollection& other) :
            PopulationTreeOperatorSettingsCollection(other) { }
};



class TreePriorSettings {
    public:
        TreePriorSettings() { }
        
        virtual EcoevolityOptions::TreePrior get_type() const = 0;
        virtual std::string to_string(unsigned int indent_level = 0) const = 0;
        virtual bool root_height_is_fixed() const {
            return false;
        }
        virtual void update_operator_weights(std::shared_ptr<GeneralTreeOperatorSettingsCollection> op_collection) const = 0;
        virtual void set_defaults() = 0;
};

class UniformTreePriorSettings : public TreePriorSettings {
    public:
        UniformTreePriorSettings() : TreePriorSettings() { }
};

class UniformRootAndBetasPriorSettings : public UniformTreePriorSettings {
    public:
        UniformRootAndBetasPriorSettings() : UniformTreePriorSettings() {
            this->set_defaults();
        }
        UniformRootAndBetasPriorSettings(const UniformRootAndBetasPriorSettings & other) {
            this->root_height_settings_ = other.root_height_settings_;
            this->alpha_of_node_height_beta_settings_ = other.alpha_of_node_height_beta_settings_;
        }
        UniformRootAndBetasPriorSettings& operator=(const UniformRootAndBetasPriorSettings & other) {
            this->root_height_settings_ = other.root_height_settings_;
            this->alpha_of_node_height_beta_settings_ = other.alpha_of_node_height_beta_settings_;
            return * this;
        }
        UniformRootAndBetasPriorSettings(const YAML::Node& node) : UniformRootAndBetasPriorSettings() {
            this->init_from_yaml_node(node);
        }

        EcoevolityOptions::TreePrior get_type() const {
            return EcoevolityOptions::TreePrior::uniform_root_and_betas;
        }

        bool root_height_is_fixed() const {
            return this->root_height_settings_.is_fixed();
        }
        bool alpha_of_node_height_beta_prior_is_fixed() const {
            return this->alpha_of_node_height_beta_settings_.is_fixed();
        }

        std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            ss << std::boolalpha;
            std::string margin = string_util::get_indent(indent_level);
            std::string indent = string_util::get_indent(1);
            ss << margin << EcoevolityOptions::get_tree_prior_string(this->get_type()) << ":\n"
               << margin << indent << "parameters:\n"
               << margin << indent << indent << "root_height:\n"
               << this->root_height_settings_.to_string(indent_level + 3)
               << margin << indent << indent << "alpha_of_node_height_beta_prior:\n"
               << this->alpha_of_node_height_beta_settings_.to_string(indent_level + 3);
            return ss.str();
        }

        virtual void set_defaults() {
            YAML::Node n;
            std::stringstream ss;

            ss << "estimate: true\n"
               << "prior:\n"
               << "    exponential_distribution:\n"
               << "        mean: 1.0\n";
            n = YAML::Load(ss);
            this->root_height_settings_ = PositiveRealParameterSettings(n);

            ss.str("");
            ss.clear();
            ss << "value: 1.0\n"
               << "estimate: false\n";
            n = YAML::Load(ss);
            this->alpha_of_node_height_beta_settings_ = PositiveRealParameterSettings(n);
        }

        void update_operator_weights(std::shared_ptr<GeneralTreeOperatorSettingsCollection> op_collection) const {
            if (this->root_height_is_fixed()) {
                op_collection->tunable_operators.at("TreeScaler").set_weight(0);
                op_collection->tunable_operators.at("RootHeightScaler").set_weight(0);
                op_collection->tunable_operators.at("TreeScaler").set_weight(0);
                op_collection->tunable_operators.at("RootHeightScaler").set_weight(0);

                if (op_collection->tunable_operators.count("GlobalHeightSizeRateScaler") > 0) {
                    op_collection->tunable_operators.at("GlobalHeightSizeRateScaler").set_weight(0);
                }
                if (op_collection->tunable_operators.count("GlobalHeightRateScaler") > 0) {
                    op_collection->tunable_operators.at("GlobalHeightRateScaler").set_weight(0);
                }
                if (op_collection->tunable_operators.count("GlobalHeightSizeScaler") > 0) {
                    op_collection->tunable_operators.at("GlobalHeightSizeScaler").set_weight(0);
                }
            }
            if (this->alpha_of_node_height_beta_prior_is_fixed()) {
                if (op_collection->tunable_operators.count("NodeHeightPriorAlphaMover") > 0) {
                    op_collection->tunable_operators.at("NodeHeightPriorAlphaMover").set_weight(0);
                }
                if (op_collection->tunable_operators.count("NodeHeightPriorAlphaScaler") > 0) {
                    op_collection->tunable_operators.at("NodeHeightPriorAlphaScaler").set_weight(0);
                }
            }
        }

    protected:
        PositiveRealParameterSettings root_height_settings_;
        PositiveRealParameterSettings alpha_of_node_height_beta_settings_;

        void init_from_yaml_node(const YAML::Node& node) {
            if (! node.IsMap()) {
                std::string message = (
                        "Expecting uniform_root_and_betas to be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
                throw EcoevolityYamlConfigError(message);
            }

            std::unordered_set<std::string> keys;
            for (YAML::const_iterator arg = node.begin();
                    arg != node.end();
                    ++arg) {
                if (keys.count(arg->first.as<std::string>()) > 0) {
                    std::string message = (
                            "Duplicate key in uniform_root_and_betas: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
                keys.insert(arg->first.as<std::string>());

                if (arg->first.as<std::string>() == "parameters") {
                    this->parse_parameters(arg->second);
                }
                else {
                    std::string message = (
                            "Unrecognized uniform_root_and_betas setting: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

        void parse_parameters(const YAML::Node& node) {
            if (! node.IsMap()) {
                std::string message = (
                        "Expecting uniform_root_and_betas parameters to be a map, but found: " +
                        YamlCppUtils::get_node_type(node));
                throw EcoevolityYamlConfigError(message);
            }

            std::unordered_set<std::string> keys;
            for (YAML::const_iterator arg = node.begin();
                    arg != node.end();
                    ++arg) {
                if (keys.count(arg->first.as<std::string>()) > 0) {
                    std::string message = (
                            "Duplicate key in uniform_root_and_betas parameters: " +
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
                            "Unrecognized uniform_root_and_betas parameters setting: " +
                            arg->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }
};


class TreeModelSettings {
    public:
        std::shared_ptr<TreePriorSettings> tree_prior = std::make_shared<UniformRootAndBetasPriorSettings>();

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

        EcoevolityOptions::TreeSpace get_tree_space() const {
            return this->tree_space_;
        }
        std::string get_tree_space_string() const {
            return EcoevolityOptions::get_tree_space_string(this->tree_space_);
        }

        void set_tree_space(const EcoevolityOptions::TreeSpace ts) {
            this->tree_space_ = ts;
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
                    std::string tree_space_value = arg->second.as<std::string>();

                    if (! EcoevolityOptions::is_tree_space(tree_space_value)) {
                        std::string message = (
                                "Unrecognized tree_space setting: " + tree_space_value);
                        throw EcoevolityYamlConfigError(message);
                    }
                    this->set_tree_space(EcoevolityOptions::get_tree_space(
                                tree_space_value));
                }
                else if (arg->first.as<std::string>() == "tree_prior") {
                    if (! arg->first.IsScalar()) {
                        std::string message = (
                                "Expecting tree_prior node to be a scalar, but found: " +
                                YamlCppUtils::get_node_type(arg->first));
                        throw EcoevolityYamlConfigError(message);
                    }
                    if (arg->second.size() != 1) {
                        throw EcoevolityYamlConfigError(
                                "tree_prior node should have one key");
                    }
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

        void parse_tree_prior(const YAML::Node& node) {
            std::unordered_set<std::string> keys;
            unsigned int key_count = 0;
            for (YAML::const_iterator arg = node.begin();
                    arg != node.end();
                    ++arg) {
                if (key_count > 0) {
                    throw EcoevolityYamlConfigError(
                            "tree_prior node should have at most one key");

                }
                ++key_count;

                std::string tree_prior_key = arg->first.as<std::string>();

                if (! EcoevolityOptions::is_tree_prior(tree_prior_key)) {
                    std::string message = (
                            "Unrecognized tree_prior setting: " + tree_prior_key);
                    throw EcoevolityYamlConfigError(message);
                }
                if (! node.IsMap()) {
                    std::string message = (
                            "Expecting " + tree_prior_key + " node to be a map, but found: " +
                            YamlCppUtils::get_node_type(node));
                    throw EcoevolityYamlConfigError(message);
                }

                if (EcoevolityOptions::get_tree_prior(tree_prior_key) ==
                        EcoevolityOptions::TreePrior::uniform_root_and_betas) {
                    this->tree_prior = std::make_shared<UniformRootAndBetasPriorSettings>(
                            node[tree_prior_key]);
                }
                else {
                    std::string message = (
                            "Unsupported tree_model setting: " + tree_prior_key);
                    throw EcoevolityYamlConfigError(message);
                }
            }
        }

        void update_operator_weights(std::shared_ptr<GeneralTreeOperatorSettingsCollection> op_collection) const {
            if (this->get_tree_space() == EcoevolityOptions::TreeSpace::bifurcating) {
                // Turn off RJ move that moves through generalized tree space
                op_collection->untunable_operators.at("SplitLumpNodesRevJumpSampler").set_weight(0);
            }
            else if (this->get_tree_space() == EcoevolityOptions::TreeSpace::fixed_tree) {
                // Turn off topology changing moves
                op_collection->untunable_operators.at("SplitLumpNodesRevJumpSampler").set_weight(0);
                op_collection->untunable_operators.at("NeighborHeightNodeSwap").set_weight(0);
                op_collection->untunable_operators.at("NeighborHeightNodePermute").set_weight(0);
                op_collection->untunable_operators.at("NodeHeightSlideBumpSwapMover").set_weight(0);
                op_collection->untunable_operators.at("NodeHeightSlideBumpSwapScaler").set_weight(0);
                op_collection->untunable_operators.at("NodeHeightSlideBumpPermuteMover").set_weight(0);
                op_collection->untunable_operators.at("NodeHeightSlideBumpPermuteScaler").set_weight(0);
            }
            this->tree_prior->update_operator_weights(op_collection);
        }

        virtual std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            ss << std::boolalpha;
            std::string margin = string_util::get_indent(indent_level);
            std::string indent = string_util::get_indent(1);
            ss << margin << "tree_space: " << this->get_tree_space_string() << "\n"
               << margin << "tree_prior:\n"
               << this->tree_prior->to_string(indent_level + 1);
            return ss.str();
        }

    protected:
        EcoevolityOptions::TreeSpace tree_space_ = EcoevolityOptions::TreeSpace::generalized;

        void clear() {
            this->tree_space_ = EcoevolityOptions::TreeSpace::generalized;
            this->tree_prior = std::make_shared<UniformRootAndBetasPriorSettings>();
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

        friend class PopulationTreeSettings;

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
                    if (! path::isfile(this->path_)) {
                        std::string message = (
                                "yaml data path \'" + this->path_ + "\' is not a file");
                        throw EcoevolityYamlConfigError(message);
                    }
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
                    if (! path::isfile(this->path_)) {
                        std::string message = (
                                "alignment path \'" + this->path_ + "\' is not a file");
                        throw EcoevolityYamlConfigError(message);
                    }
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
            ss << std::boolalpha;
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


class PopulationTreeSettings {
    public:
        // Public data members
        PopSizeSettings population_size_settings;
        PositiveRealParameterSettings freq_1_settings;
        PositiveRealParameterSettings mutation_rate_settings;
        GeneralTreeDataSettings data_settings;
        TreeModelSettings tree_model_settings;
        std::shared_ptr<GeneralTreeOperatorSettingsCollection> operator_settings;
        StartingTreeSettings starting_tree_settings;

        // Constructors
        PopulationTreeSettings() {
            this->set_defaults();
            this->update_operator_settings_type();
        }
        PopulationTreeSettings(
                const std::string & yaml_config_path) : PopulationTreeSettings() {
            this->init_from_config_file(yaml_config_path);
        }
        PopulationTreeSettings(
                std::istream& yaml_config_stream,
                const std::string & yaml_config_path) : PopulationTreeSettings() {
            this->init_from_config_stream(yaml_config_stream, yaml_config_path);
        }
        PopulationTreeSettings(const PopulationTreeSettings & other) {
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
            this->tree_model_settings = other.tree_model_settings;
            this->starting_tree_settings = this->starting_tree_settings;
        }
        virtual ~PopulationTreeSettings() { }
        PopulationTreeSettings & operator=(const PopulationTreeSettings & other) {
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
            this->tree_model_settings = other.tree_model_settings;
            this->starting_tree_settings = this->starting_tree_settings;
            return * this;
        }

        void set_defaults() {
            YAML::Node n;
            std::stringstream ss;

            ss << "equal_population_sizes: true\n"
               << "estimate: true\n"
               << "prior:\n"
               << "    exponential_distribution:\n"
               << "        mean: 0.001\n";
            n = YAML::Load(ss);
            this->population_size_settings = PopSizeSettings(n);

            ss.str("");
            ss.clear();
            ss << "value: 1.0\n"
               << "estimate: false\n";
            n = YAML::Load(ss);
            this->mutation_rate_settings = PositiveRealParameterSettings(n);

            ss.str("");
            ss.clear();
            ss << "value: 0.5\n"
               << "estimate: false\n";
            n = YAML::Load(ss);
            this->freq_1_settings = PositiveRealParameterSettings(n);
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

        void update_operator_weights() {
            this->tree_model_settings.update_operator_weights(this->operator_settings);

            if (mutation_rate_settings.is_fixed()) {
                this->operator_settings->tunable_operators.at("MuRateScaler").set_weight(0);
            }
            if (freq_1_settings.is_fixed()) {
                this->operator_settings->tunable_operators.at("StateFreqMover").set_weight(0);
                this->operator_settings->tunable_operators.at("StateFreqDirichletOperator").set_weight(0);
            }

            const bool root_height_is_fixed = tree_model_settings.tree_prior->root_height_is_fixed();
            const bool pop_sizes_fixed = population_size_settings.is_fixed();
            const bool pop_sizes_constrained = population_size_settings.population_sizes_are_constrained();

            if (root_height_is_fixed || pop_sizes_constrained || pop_sizes_fixed) {
                this->operator_settings->tunable_operators.at("RootHeightSizeMixer").set_weight(0);
            }
            if (root_height_is_fixed || pop_sizes_fixed) {
                this->operator_settings->tunable_operators.at("GlobalHeightSizeMixer").set_weight(0);
            }
            if (pop_sizes_fixed || pop_sizes_constrained) {
                this->operator_settings->tunable_operators.at("HeightSizeMixer").set_weight(0);
                this->operator_settings->tunable_operators.at("HeightSizeSlideBumpMixer").set_weight(0);
            }
            if (pop_sizes_fixed) {
                this->operator_settings->tunable_operators.at("GlobalPopSizeScaler").set_weight(0);
                this->operator_settings->tunable_operators.at("PopSizeScaler").set_weight(0);
            }
        }

        void write_settings(std::ostream& out) const {
            std::string indent = string_util::get_indent(1);
            out << std::boolalpha;

            out << "---\n"
                << "data:\n"
                << this->data_settings.to_string(1)
                << "tree_model:\n"
                << this->tree_model_settings.to_string(1)
                << "starting_tree:\n"
                << this->starting_tree_settings.to_string(1)
                << "branch_parameters:\n"
                << indent << "population_size:\n"
                << this->population_size_settings.to_string(2)
                << "mutation_parameters:\n"
                << indent << "mutation_rate:\n"
                << this->mutation_rate_settings.to_string(2)
                << indent << "freq_1:\n"
                << this->freq_1_settings.to_string(2)
                << "mcmc_settings:\n"
                << indent << "chain_length: " << this->get_chain_length() << "\n"
                << indent << "sample_frequency: " << this->get_sample_frequency() << "\n"
                << indent << "operators:\n"
                << this->operator_settings->to_string(2);
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
            // Update operators
            this->update_operator_weights();
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

            // parse data settings
            if (! top_level_node["data"]) {
                throw EcoevolityYamlConfigError("No data");
            }
            this->data_settings.update_from_config(
                    top_level_node["data"],
                    this->config_path_);

            // Parse tree model first, because we need it for determining what
            // type off operator collection settings is needed
            if (top_level_node["tree_model"]) {
                this->tree_model_settings.update_from_config(top_level_node["tree_model"]);
                this->update_operator_settings_type();
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

                if (top->first.as<std::string>() == "data") {
                    // Already parsed above
                }
                else if (top->first.as<std::string>() == "tree_model") {
                    // Already parsed above
                }
                // parse mcmc settings
                else if (top->first.as<std::string>() == "mcmc_settings") {
                    this->parse_mcmc_settings(top->second);
                }
                // parse branch parameters
                else if (top->first.as<std::string>() == "branch_parameters") {
                    this->parse_branch_parameters(top->second);
                }
                // parse mutation parameters
                else if (top->first.as<std::string>() == "mutation_parameters") {
                    this->parse_mutation_parameters(top->second);
                }
                else if (top->first.as<std::string>() == "starting_tree") {
                    this->starting_tree_settings.update_from_config(
                            top->second,
                            this->config_path_);

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

        void update_operator_settings_type() {
            EcoevolityOptions::TreePrior tree_prior_type = this->tree_model_settings.tree_prior->get_type();
            if (tree_prior_type == EcoevolityOptions::TreePrior::uniform_root_and_betas) {
                this->operator_settings = std::make_shared<BetaPopulationTreeOperatorSettingsCollection>();
            }
            else {
                std::string message = (
                        "Do not have OperatorSettingsCollection to match tree prior: " +
                        EcoevolityOptions::get_tree_prior_string(tree_prior_type));
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
                    this->operator_settings->update_from_config(mcmc->second);
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
