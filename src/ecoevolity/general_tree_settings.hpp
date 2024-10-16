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
        StartingTreeSettings() { }
        StartingTreeSettings(const std::string & argument,
                const std::string & config_path) {
            std::string arg = string_util::strip(argument);
            std::string potential_path = path::join(
                    path::dirname(config_path),
                    arg);
            // Check if we have one of the starting tree options
            if (arg == "comb") {
                this->source_ = StartingTreeSettings::Source::option;
                this->option_ = StartingTreeSettings::Option::comb;
            }
            else if (arg == "random") {
                this->source_ = StartingTreeSettings::Source::option;
                this->option_ = StartingTreeSettings::Option::random;
            }
            else if (arg == "upgma") {
                this->source_ = StartingTreeSettings::Source::option;
                this->option_ = StartingTreeSettings::Option::upgma;
                throw EcoevolityYamlConfigError(
                        "UPGMA starting tree not implemented");
            }
            // It's not an option, is it a path?
            else if (path::isfile(potential_path)) {
                this->source_ = StartingTreeSettings::Source::path;
                this->path_ = potential_path;
            }
            // It's not an option or a path, so we will treat it as a string
            else {
                this->source_ = StartingTreeSettings::Source::string;
                this->string_ = arg;
            }
        }
        StartingTreeSettings(const StartingTreeSettings& other) {
            this->source_ = other.source_;
            this->option_ = other.option_;
            this->path_ = other.path_;
            this->string_ = other.string_;
        }
        virtual ~StartingTreeSettings() { }
        StartingTreeSettings& operator=(const StartingTreeSettings& other) {
            this->source_ = other.source_;
            this->option_ = other.option_;
            this->path_ = other.path_;
            this->string_ = other.string_;
            return * this;
        }

        StartingTreeSettings::Source get_tree_source() const {
            return this->source_;
        }
        StartingTreeSettings::Option get_tree_option() const {
            return this->option_;
        }
        const std::string& get_tree_path() const {
            return this->path_;
        }
        void set_tree_path(const std::string & path) {
            this->path_ = path;
        }
        const std::string& get_tree_string() const {
            return this->string_;
        }

        std::string to_string() const {
            if (this->source_ == StartingTreeSettings::Source::option) {
                if (this->option_ == StartingTreeSettings::Option::random) {
                    return "random";
                }
                else if (this->option_ == StartingTreeSettings::Option::comb) {
                    return "comb";
                }
                else if (this->option_ == StartingTreeSettings::Option::upgma) {
                    return "upgma";
                }
                else {
                    throw EcoevolityYamlConfigError(
                            "Unexpected StartingTreeSettings::option_");
                }
            }
            else if (this->source_ == StartingTreeSettings::Source::string) {
                return this->string_;
            }
            else if (this->source_ == StartingTreeSettings::Source::path) {
                return this->path_;
            }
            throw EcoevolityYamlConfigError(
                    "Unexpected StartingTreeSettings::source_");
        }
};


class GeneralTreeOperatorSettings {
    protected:
        std::vector<double> weights_;

    public:
        GeneralTreeOperatorSettings() { }
        GeneralTreeOperatorSettings(const std::vector<double> & weights) {
            this->weights_ = weights;
        }
        GeneralTreeOperatorSettings(double weight) {
            this->weights_ = { weight };
        }
        GeneralTreeOperatorSettings(const GeneralTreeOperatorSettings& other) {
            this->weights_ = other.weights_;
        }
        virtual ~GeneralTreeOperatorSettings() { }
        GeneralTreeOperatorSettings& operator=(const GeneralTreeOperatorSettings& other) {
            this->weights_ = other.weights_;
            return * this;
        }

        void init(const std::vector<double> & weights) {
            this->weights_ = weights;
        }

        virtual void clear() {
            this->weights_.clear();
        }
        virtual void turn_off() {
            this->clear();
            this->weights_.push_back(0.0);
        }

        virtual bool is_tunable() const {
            return false;
        }

        double get_weight(unsigned int index) const {
            return this->weights_.at(index);
        }
        const std::vector<double> & get_weights() const {
            return this->weights_;
        }

        unsigned int get_number_of_operators() const {
            return this->weights_.size();
        }
        bool is_empty() const {
            return this->weights_.empty();
        }

        virtual void update_from_config(const YAML::Node& operator_node,
                const bool auto_opt_by_default = false) {
            if (! operator_node.IsMap()) {
                std::string message = (
                        "Expecting operator parameters to be a map, but found: " +
                        YamlCppUtils::get_node_type(operator_node));
                throw EcoevolityYamlConfigError(message);
            }
            if (operator_node.size() != 1) {
                throw EcoevolityYamlConfigError(
                        "untunable operator node should only have a single key");
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
                    if (p->second.IsScalar()) {
                        double wt = p->second.as<double>();
                        if (wt < 0.0) {
                            throw EcoevolityYamlConfigError("operator weight cannot be negative");
                        }
                        this->init( { wt } );
                    }
                    else if (p->second.IsSequence()) {
                        std::vector<double> weights;
                        for (YAML::const_iterator wt_node = p->second.begin();
                                wt_node != p->second.end();
                                ++wt_node) {
                            double wt = wt_node->as<double>();
                            if (wt < 0.0) {
                                throw EcoevolityYamlConfigError("operator weight cannot be negative");
                            }
                            weights.push_back(wt);
                        }
                        this->init(weights);
                    }
                    else {
                        throw EcoevolityYamlConfigError("Operator weight node is not a scalar or sequence");
                    }
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
            if (! this->weights_.empty()) {
                ss << margin << "weight: ";
                if (this->weights_.size() == 1) {
                    ss << this->weights_.at(0) << "\n";
                }
                else {
                    ss << "[" << this->weights_.at(0);
                    for (unsigned int i = 1; i < this->weights_.size(); ++i) {
                        ss << ", " << this->weights_.at(i);
                    }
                    ss << "]\n";
                }
            }
            return ss.str();
        }
};

class GeneralTreeTunableOperatorSettings : public GeneralTreeOperatorSettings {
    protected:
        std::vector<double> tuning_parameters_;
        std::vector<bool> auto_optimize_bools_;
        std::vector<int> auto_optimize_delays_;

    public:
        GeneralTreeTunableOperatorSettings() : GeneralTreeOperatorSettings() { }
        GeneralTreeTunableOperatorSettings(
                const std::vector<double> & weights
                ) : GeneralTreeOperatorSettings() {
            std::vector<double> tuning_parameters(weights.size(), -1.0);
            std::vector<bool> auto_opt_bools(weights.size(), true);
            std::vector<int> auto_opt_delays(weights.size(), -1);
            this->init(weights, tuning_parameters, auto_opt_bools, auto_opt_delays);
        }
        GeneralTreeTunableOperatorSettings(
                double weight
                ) : GeneralTreeOperatorSettings() {
            std::vector<double> weights = { weight };
            std::vector<double> tuning_parameters(weights.size(), -1.0);
            std::vector<bool> auto_opt_bools(weights.size(), true);
            std::vector<int> auto_opt_delays(weights.size(), -1);
            this->init(weights, tuning_parameters, auto_opt_bools, auto_opt_delays);
        }
        GeneralTreeTunableOperatorSettings(
                const std::vector<double> & weights,
                const std::vector<double> & tuning_parameters
                ) : GeneralTreeOperatorSettings() {
            std::vector<bool> auto_opt_bools(weights.size(), true);
            std::vector<int> auto_opt_delays(weights.size(), -1);
            this->init(weights, tuning_parameters, auto_opt_bools, auto_opt_delays);
        }
        GeneralTreeTunableOperatorSettings(
                double weight,
                double tuning_parameter
                ) : GeneralTreeOperatorSettings() {
            std::vector<double> weights = { weight };
            std::vector<double> tuning_parameters = { tuning_parameter };
            std::vector<bool> auto_opt_bools(weights.size(), true);
            std::vector<int> auto_opt_delays(weights.size(), -1);
            this->init(weights, tuning_parameters, auto_opt_bools, auto_opt_delays);
        }
        GeneralTreeTunableOperatorSettings(
                const std::vector<double> & weights,
                const std::vector<double> & tuning_parameters,
                const std::vector<bool> & auto_optimize_bools
                ) : GeneralTreeOperatorSettings() {
            std::vector<int> auto_opt_delays(weights.size(), -1);
            this->init(weights, tuning_parameters, auto_optimize_bools, auto_opt_delays);
        }
        GeneralTreeTunableOperatorSettings(
                double weight,
                double tuning_parameter,
                bool auto_optimize
                ) : GeneralTreeOperatorSettings() {
            std::vector<double> weights = { weight };
            std::vector<double> tuning_parameters = { tuning_parameter };
            std::vector<bool> auto_opt_bools = { auto_optimize };
            std::vector<int> auto_opt_delays(weights.size(), -1);
            this->init(weights, tuning_parameters, auto_opt_bools, auto_opt_delays);
        }
        GeneralTreeTunableOperatorSettings(const GeneralTreeTunableOperatorSettings& other) :
                GeneralTreeOperatorSettings(other) {
            this->tuning_parameters_ = other.tuning_parameters_;
            this->auto_optimize_bools_ = other.auto_optimize_bools_;
            this->auto_optimize_delays_ = other.auto_optimize_delays_;
        }
        virtual ~GeneralTreeTunableOperatorSettings() { }
        GeneralTreeTunableOperatorSettings& operator=(const GeneralTreeTunableOperatorSettings& other) {
            this->weights_ = other.weights_;
            this->tuning_parameters_ = other.tuning_parameters_;
            this->auto_optimize_bools_ = other.auto_optimize_bools_;
            this->auto_optimize_delays_ = other.auto_optimize_delays_;
            return * this;
        }

        void init(const std::vector<double> & weights,
                const std::vector<double> & tuning_parameters,
                const std::vector<bool> & auto_optimize_bools,
                const std::vector<int> auto_optimize_delays) {
            if (weights.size() != tuning_parameters.size()) {
                throw EcoevolityError("Number of tuning parameters does not equal number of weights");
            }
            if (weights.size() != auto_optimize_bools.size()) {
                throw EcoevolityError("Number of auto opt bools does not equal number of weights");
            }
            if (weights.size() != auto_optimize_delays.size()) {
                throw EcoevolityError("Number of auto opt delays does not equal number of weights");
            }
            this->weights_ = weights;
            this->tuning_parameters_ = tuning_parameters;
            this->auto_optimize_bools_ = auto_optimize_bools;
            this->auto_optimize_delays_ = auto_optimize_delays;
        }

        bool is_tunable() const {
            return true;
        }

        virtual void clear() {
            this->weights_.clear();
            this->tuning_parameters_.clear();
            this->auto_optimize_bools_.clear();
            this->auto_optimize_delays_.clear();
        }
        void turn_off() {
            this->clear();
            this->weights_.push_back(0.0);
            this->tuning_parameters_.push_back(-1.0);
            this->auto_optimize_bools_.push_back(false);
            this->auto_optimize_delays_.push_back(-1);
        }

        double get_tuning_parameter(unsigned int index) const {
            return this->tuning_parameters_.at(index);
        }
        const std::vector<double> & get_tuning_parameters() const {
            return this->tuning_parameters_;
        }

        int get_auto_optimize_delay(unsigned int index) const {
            return this->auto_optimize_delays_.at(index);
        }
        const std::vector<int> & get_auto_optimize_delays() const {
            return this->auto_optimize_delays_;
        }

        bool auto_optimizing(unsigned int index) const {
            return this->auto_optimize_bools_.at(index);
        }
        const std::vector<bool> & get_auto_optimize_bools() const {
            return this->auto_optimize_bools_;
        }

        void update_from_config(const YAML::Node& operator_node,
                const bool auto_opt_by_default = true) {
            if (! operator_node.IsMap()) {
                std::string message = (
                        "Expecting operator parameters to be a map, but found: " +
                        YamlCppUtils::get_node_type(operator_node));
                throw EcoevolityYamlConfigError(message);
            }

            std::vector<double> weights;
            std::vector<double> tuning_parameters;
            std::vector<int> opt_delays;
            std::vector<bool> auto_opt_bools;

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
                    if (p->second.IsScalar()) {
                        double wt = p->second.as<double>();
                        if (wt < 0.0) {
                            throw EcoevolityYamlConfigError("operator weight cannot be negative");
                        }
                        weights.push_back(wt);
                    }
                    else if (p->second.IsSequence()) {
                        for (YAML::const_iterator wt_node = p->second.begin();
                                wt_node != p->second.end();
                                ++wt_node) {
                            double wt = wt_node->as<double>();
                            if (wt < 0.0) {
                                throw EcoevolityYamlConfigError("operator weight cannot be negative");
                            }
                            weights.push_back(wt);
                        }
                    }
                    else {
                        throw EcoevolityYamlConfigError("Operator weight node is not a scalar or sequence");
                    }
                }
                else if (p->first.as<std::string>() == "tuning_parameter") {
                    if (p->second.IsScalar()) {
                        double tuning_param = p->second.as<double>();
                        if (tuning_param <= 0.0) {
                            throw EcoevolityYamlConfigError("tuning_parameter must be positive");
                        }
                        tuning_parameters.push_back(tuning_param);
                    }
                    else if (p->second.IsSequence()) {
                        for (YAML::const_iterator tp_node = p->second.begin();
                                tp_node != p->second.end();
                                ++tp_node) {
                            double tp = tp_node->as<double>();
                            if (tp <= 0.0) {
                                throw EcoevolityYamlConfigError("tuning_parameter must be positive");
                            }
                            tuning_parameters.push_back(tp);
                        }
                    }
                    else {
                        throw EcoevolityYamlConfigError("Operator tuning_parameter node is not a scalar or sequence");
                    }
                }
                else if (p->first.as<std::string>() == "auto_optimize") {
                    if (p->second.IsScalar()) {
                        bool auto_opt_bool = p->second.as<bool>();
                        auto_opt_bools.push_back(auto_opt_bool);
                    }
                    else if (p->second.IsSequence()) {
                        for (YAML::const_iterator auto_opt_node = p->second.begin();
                                auto_opt_node != p->second.end();
                                ++auto_opt_node) {
                            bool aob = auto_opt_node->as<bool>();
                            auto_opt_bools.push_back(aob);
                        }
                    }
                    else {
                        throw EcoevolityYamlConfigError("Operator auto_optimize node is not a scalar or sequence");
                    }
                }
                else if (p->first.as<std::string>() == "auto_optimize_delay") {
                    if (p->second.IsScalar()) {
                        int opt_delay = p->second.as<int>();
                        if (opt_delay <= 0) {
                            throw EcoevolityYamlConfigError("auto_optimize_delay must be positive");
                        }
                        opt_delays.push_back(opt_delay);
                    }
                    else if (p->second.IsSequence()) {
                        for (YAML::const_iterator delay_node = p->second.begin();
                                delay_node != p->second.end();
                                ++delay_node) {
                            int d = delay_node->as<int>();
                            if (d <= 0) {
                                throw EcoevolityYamlConfigError("auto_optimize_delay must be positive");
                            }
                            opt_delays.push_back(d);
                        }
                    }
                    else {
                        throw EcoevolityYamlConfigError("Operator auto_optimize_delay node is not a scalar or sequence");
                    }
                }
                else {
                    std::string message = (
                            "Unrecognized key in operator parameters: " +
                            p->first.as<std::string>());
                    throw EcoevolityYamlConfigError(message);
                }
            }

            std::vector<size_t> n_ops_vector =  {
                    weights.size(),
                    tuning_parameters.size(),
                    opt_delays.size(),
                    auto_opt_bools.size()};
            size_t n_ops = *std::max_element(std::begin(n_ops_vector), std::end(n_ops_vector));
            if (n_ops < 1) {
                throw EcoevolityYamlConfigError("operator listed with no options");
            }
            for (auto n : n_ops_vector) {
                if ((n > 0) && (n != n_ops)) {
                    throw EcoevolityYamlConfigError("operator with unequal number of options");
                }
            }
            if (weights.empty()) {
                weights = std::vector<double>(n_ops, -1.0);
            }
            if (tuning_parameters.empty()) {
                tuning_parameters = std::vector<double>(n_ops, -1.0);
            }
            if (opt_delays.empty()) {
                opt_delays = std::vector<int>(n_ops, -1);
            }
            if (auto_opt_bools.empty()) {
                auto_opt_bools = std::vector<bool>(n_ops, true);
                if (! auto_opt_by_default) {
                    auto_opt_bools = std::vector<bool>(n_ops, false);
                }
            }

            this->init(weights, tuning_parameters, auto_opt_bools, opt_delays);
        }

        std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            ss << std::boolalpha;
            std::string margin = string_util::get_indent(indent_level);
            if (! this->weights_.empty()) {
                ss << margin << "weight: ";
                if (this->weights_.size() == 1) {
                    ss << this->weights_.at(0) << "\n";
                }
                else {
                    ss << "[" << this->weights_.at(0);
                    for (unsigned int i = 1; i < this->weights_.size(); ++i) {
                        ss << ", " << this->weights_.at(i);
                    }
                    ss << "]\n";
                }
            }
            if (! this->tuning_parameters_.empty()) {
                ss << margin << "tuning_parameter: ";
                if (this->tuning_parameters_.size() == 1) {
                    ss << this->tuning_parameters_.at(0) << "\n";
                }
                else {
                    ss << "[" << this->tuning_parameters_.at(0);
                    for (unsigned int i = 1; i < this->tuning_parameters_.size(); ++i) {
                        ss << ", " << this->tuning_parameters_.at(i);
                    }
                    ss << "]\n";
                }
            }
            if (! this->auto_optimize_bools_.empty()) {
                ss << margin << "auto_optimize: ";
                if (this->auto_optimize_bools_.size() == 1) {
                    ss << this->auto_optimize_bools_.at(0) << "\n";
                }
                else {
                    ss << "[" << this->auto_optimize_bools_.at(0);
                    for (unsigned int i = 1; i < this->auto_optimize_bools_.size(); ++i) {
                        ss << ", " << this->auto_optimize_bools_.at(i);
                    }
                    ss << "]\n";
                }
            }
            if (! this->auto_optimize_delays_.empty()) {
                ss << margin << "auto_optimize_delay: ";
                if (this->auto_optimize_delays_.size() == 1) {
                    ss << this->auto_optimize_delays_.at(0) << "\n";
                }
                else {
                    ss << "[" << this->auto_optimize_delays_.at(0);
                    for (unsigned int i = 1; i < this->auto_optimize_delays_.size(); ++i) {
                        ss << ", " << this->auto_optimize_delays_.at(i);
                    }
                    ss << "]\n";
                }
            }
            return ss.str();
        }
};


class GeneralTreeOperatorSettingsCollection {
    public:
        std::map<std::string, GeneralTreeOperatorSettings> untunable_operators {
            // If not set in the config, negative weights will be updated to
            // default weight for the operator class and zero weights will
            // prevent operator from being created.
            {"NeighborHeightNodePermute",           GeneralTreeOperatorSettings()},
            {"NeighborHeightNodeSwap",              GeneralTreeOperatorSettings()},
            {"NeighborHeightNodeSwapAll",           GeneralTreeOperatorSettings()},
        };
        std::map<std::string, GeneralTreeTunableOperatorSettings> tunable_operators {
            // If not set in the config, negative weights will be updated to
            // default weight for the operator class and zero weights will
            // prevent operator from being created.
            {"SplitLumpNodesRevJumpSampler",        GeneralTreeTunableOperatorSettings()},
            {"TreeScaler",                          GeneralTreeTunableOperatorSettings()},
            {"NodeHeightScaler",                    GeneralTreeTunableOperatorSettings()},
            {"NodeHeightMover",                     GeneralTreeTunableOperatorSettings()},
            {"NodeHeightSlideBumpScaler",           GeneralTreeTunableOperatorSettings()},
            {"NodeHeightSlideBumpPermuteScaler",    GeneralTreeTunableOperatorSettings()},
            {"NodeHeightSlideBumpSwapScaler",       GeneralTreeTunableOperatorSettings()},
            {"NodeHeightSlideBumpMover",            GeneralTreeTunableOperatorSettings()},
            {"NodeHeightSlideBumpPermuteMover",     GeneralTreeTunableOperatorSettings()},
            {"NodeHeightSlideBumpSwapMover",        GeneralTreeTunableOperatorSettings()},
            {"RootHeightScaler",                    GeneralTreeTunableOperatorSettings()},
            {"GlobalNodeHeightDirichletOperator",   GeneralTreeTunableOperatorSettings()},
            {"NodeHeightDirichletOperator",         GeneralTreeTunableOperatorSettings()},
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
                        "Expecting operators to be a map, but found: " +
                        YamlCppUtils::get_node_type(operator_node));
                throw EcoevolityYamlConfigError(message);
            }

            std::unordered_set<std::string> keys;
            for (YAML::const_iterator p = operator_node.begin();
                    p != operator_node.end();
                    ++p) {
                std::string name = p->first.as<std::string>();
                if (keys.count(name) > 0) {
                    std::string message = (
                            "Duplicate operator: " +
                            name);
                    throw EcoevolityYamlConfigError(message);
                }
                keys.insert(name);

                if (this->is_untunable_operator(name)) {
                    this->untunable_operators[name].update_from_config(p->second);
                }
                else if (this->is_tunable_operator(name)) {
                    if (name == "SplitLumpNodesRevJumpSampler") {
                        this->tunable_operators[name].update_from_config(p->second, false);
                    }
                    else {
                        this->tunable_operators[name].update_from_config(p->second, true);
                    }
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
            // If not set in the config, negative weights will be updated to
            // default weight for the operator class and zero weights will
            // prevent operator from being created.
            this->tunable_operators["NodeHeightPriorAlphaScaler"]  = GeneralTreeTunableOperatorSettings();
            this->tunable_operators["NodeHeightPriorAlphaMover"]   = GeneralTreeTunableOperatorSettings();
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
            // If not set in the config, negative weights will be updated to
            // default weight for the operator class and zero weights will
            // prevent operator from being created.
            this->tunable_operators["MuRateScaler"]                = GeneralTreeTunableOperatorSettings();
            this->tunable_operators["GlobalPopSizeScaler"]         = GeneralTreeTunableOperatorSettings();
            this->tunable_operators["PopSizeScaler"]               = GeneralTreeTunableOperatorSettings();
            this->tunable_operators["GlobalHeightSizeMixer"]       = GeneralTreeTunableOperatorSettings();
            this->tunable_operators["HeightSizeMixer"]             = GeneralTreeTunableOperatorSettings();
            this->tunable_operators["HeightSizeSlideBumpMixer"]    = GeneralTreeTunableOperatorSettings();
            this->tunable_operators["RootHeightSizeMixer"]         = GeneralTreeTunableOperatorSettings();
            this->tunable_operators["GlobalHeightSizeRateScaler"]  = GeneralTreeTunableOperatorSettings();
            this->tunable_operators["GlobalHeightSizeScaler"]      = GeneralTreeTunableOperatorSettings();
            this->tunable_operators["GlobalHeightRateScaler"]      = GeneralTreeTunableOperatorSettings();
            this->tunable_operators["StateFreqMover"]              = GeneralTreeTunableOperatorSettings();
            this->tunable_operators["StateFreqDirichletOperator"]  = GeneralTreeTunableOperatorSettings();
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
            // If not set in the config, negative weights will be updated to
            // default weight for the operator class and zero weights will
            // prevent operator from being created.
            this->tunable_operators["NodeHeightPriorAlphaScaler"]  = GeneralTreeTunableOperatorSettings();
            this->tunable_operators["NodeHeightPriorAlphaMover"]   = GeneralTreeTunableOperatorSettings();
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

        virtual std::map<std::string, PositiveRealParameterSettings>
        get_parameter_settings() const = 0;
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
            this->alpha_of_node_height_settings_ = other.alpha_of_node_height_settings_;
        }
        UniformRootAndBetasPriorSettings& operator=(const UniformRootAndBetasPriorSettings & other) {
            this->root_height_settings_ = other.root_height_settings_;
            this->alpha_of_node_height_settings_ = other.alpha_of_node_height_settings_;
            return * this;
        }
        UniformRootAndBetasPriorSettings(const YAML::Node& node) : UniformRootAndBetasPriorSettings() {
            this->init_from_yaml_node(node);
        }

        std::map<std::string, PositiveRealParameterSettings>
        get_parameter_settings() const {
            std::map<std::string, PositiveRealParameterSettings> p;
            p["root_height"] = this->root_height_settings_;
            p["alpha_of_node_height_beta"] = this->alpha_of_node_height_settings_;
            return p;
        }

        EcoevolityOptions::TreePrior get_type() const {
            return EcoevolityOptions::TreePrior::uniform_root_and_betas;
        }

        bool root_height_is_fixed() const {
            return this->root_height_settings_.is_fixed();
        }
        bool alpha_of_node_height_beta_prior_is_fixed() const {
            return this->alpha_of_node_height_settings_.is_fixed();
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
               << this->alpha_of_node_height_settings_.to_string(indent_level + 3);
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
            this->alpha_of_node_height_settings_ = PositiveRealParameterSettings(n);
        }

        void update_operator_weights(std::shared_ptr<GeneralTreeOperatorSettingsCollection> op_collection) const {
            if (this->root_height_is_fixed()) {
                op_collection->tunable_operators.at("TreeScaler").turn_off();
                op_collection->tunable_operators.at("RootHeightScaler").turn_off();
                op_collection->tunable_operators.at("TreeScaler").turn_off();
                op_collection->tunable_operators.at("RootHeightScaler").turn_off();

                if (op_collection->tunable_operators.count("GlobalHeightSizeRateScaler") > 0) {
                    op_collection->tunable_operators.at("GlobalHeightSizeRateScaler").turn_off();
                }
                if (op_collection->tunable_operators.count("GlobalHeightRateScaler") > 0) {
                    op_collection->tunable_operators.at("GlobalHeightRateScaler").turn_off();
                }
                if (op_collection->tunable_operators.count("GlobalHeightSizeScaler") > 0) {
                    op_collection->tunable_operators.at("GlobalHeightSizeScaler").turn_off();
                }
            }
            if (this->alpha_of_node_height_beta_prior_is_fixed()) {
                if (op_collection->tunable_operators.count("NodeHeightPriorAlphaMover") > 0) {
                    op_collection->tunable_operators.at("NodeHeightPriorAlphaMover").turn_off();
                }
                if (op_collection->tunable_operators.count("NodeHeightPriorAlphaScaler") > 0) {
                    op_collection->tunable_operators.at("NodeHeightPriorAlphaScaler").turn_off();
                }
            }
        }

    protected:
        PositiveRealParameterSettings root_height_settings_;
        PositiveRealParameterSettings alpha_of_node_height_settings_;

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
                    this->root_height_settings_ = PositiveRealParameterSettings(
                            arg->second,
                            true);
                }
                else if (arg->first.as<std::string>() == "alpha_of_node_height_beta_prior") {
                    this->alpha_of_node_height_settings_ = PositiveRealParameterSettings(arg->second);
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
    protected:
        EcoevolityOptions::TreeSpace tree_space_ = EcoevolityOptions::TreeSpace::generalized;

        void clear() {
            this->tree_space_ = EcoevolityOptions::TreeSpace::generalized;
            this->tree_prior = std::make_shared<UniformRootAndBetasPriorSettings>();
            this->starting_tree_settings = StartingTreeSettings();
        }

        void check_settings() const {
            if (
                (this->tree_space_ == EcoevolityOptions::TreeSpace::fixed_tree)
                &&
                (this->starting_tree_settings.get_tree_source() == StartingTreeSettings::Source::option)
               )
            {
                throw EcoevolityYamlConfigError(
                        "tree_space is fixed, but a starting_tree was not provided");
            }
            if (
                (this->tree_space_ == EcoevolityOptions::TreeSpace::bifurcating)
                &&
                (
                 (this->starting_tree_settings.get_tree_source()
                     == StartingTreeSettings::Source::option)
                 &&
                 (this->starting_tree_settings.get_tree_option()
                     == StartingTreeSettings::Option::comb)
                )
               )
            {
                throw EcoevolityYamlConfigError(
                        "tree_space is bifurcating, but a starting_tree was comb");
            }

            std::map<std::string, PositiveRealParameterSettings> tree_parameters;
            tree_parameters = this->tree_prior->get_parameter_settings();

            if (this->tree_prior->get_type() ==
                    EcoevolityOptions::TreePrior::uniform_root_and_betas)
            {
                PositiveRealParameterSettings root_ht = tree_parameters.at("root_height");

                if (std::isnan(root_ht.get_value())
                        && root_ht.is_fixed())
                {
                    if (this->starting_tree_settings.get_tree_source()
                            == StartingTreeSettings::Source::option)
                    {
                        throw EcoevolityYamlConfigError(
                                "root_height fixed, but no value or starting tree provided");
                    }
                }

                if (! std::isnan(root_ht.get_value()))
                {
                    if (this->starting_tree_settings.get_tree_source()
                            != StartingTreeSettings::Source::option)
                    {
                        throw EcoevolityYamlConfigError(
                                "root_height value AND starting tree provided");
                    }
                }
            }
        }

    public:
        std::shared_ptr<TreePriorSettings> tree_prior = std::make_shared<UniformRootAndBetasPriorSettings>();
        StartingTreeSettings starting_tree_settings;

        TreeModelSettings() { this->check_settings(); }
        TreeModelSettings(const TreeModelSettings & other) {
            this->tree_space_ = other.tree_space_;
            this->tree_prior = other.tree_prior;
            this->starting_tree_settings = other.starting_tree_settings;
            this->check_settings();
        }
        TreeModelSettings& operator=(const TreeModelSettings & other) {
            this->tree_space_ = other.tree_space_;
            this->tree_prior = other.tree_prior;
            this->starting_tree_settings = other.starting_tree_settings;
            this->check_settings();
            return * this;
        }
        TreeModelSettings(const YAML::Node& node,
                const std::string & config_path) {
            this->update_from_config(node, config_path);
            this->check_settings();
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

        void update_from_config(const YAML::Node& node,
                const std::string & config_path) {
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
                else if (arg->first.as<std::string>() == "starting_tree") {
                    this->starting_tree_settings = StartingTreeSettings(
                            arg->second.as<std::string>(),
                            config_path);
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
                op_collection->tunable_operators.at("SplitLumpNodesRevJumpSampler").turn_off();
            }
            else if (this->get_tree_space() == EcoevolityOptions::TreeSpace::fixed_tree) {
                // Turn off topology changing moves
                op_collection->tunable_operators.at("SplitLumpNodesRevJumpSampler").turn_off();
                op_collection->untunable_operators.at("NeighborHeightNodeSwap").turn_off();
                op_collection->untunable_operators.at("NeighborHeightNodeSwapAll").turn_off();
                op_collection->untunable_operators.at("NeighborHeightNodePermute").turn_off();

                op_collection->tunable_operators.at("NodeHeightSlideBumpSwapMover").turn_off();
                op_collection->tunable_operators.at("NodeHeightSlideBumpSwapScaler").turn_off();
                op_collection->tunable_operators.at("NodeHeightSlideBumpPermuteMover").turn_off();
                op_collection->tunable_operators.at("NodeHeightSlideBumpPermuteScaler").turn_off();
            }
            this->tree_prior->update_operator_weights(op_collection);
        }

        virtual std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            ss << std::boolalpha;
            std::string margin = string_util::get_indent(indent_level);
            std::string indent = string_util::get_indent(1);
            ss << margin << "tree_space: " << this->get_tree_space_string() << "\n"
               << margin << "starting_tree: " << this->starting_tree_settings.to_string() << "\n"
               << margin << "tree_prior:\n"
               << this->tree_prior->to_string(indent_level + 1);
            return ss.str();
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
        void set_path(const std::string & path) {
            this->path_ = path;
        }
        void set_ploidy(double p) {
            this->ploidy_ = p;
        }
        double get_ploidy() const {
            return this->ploidy_;
        }
        char get_population_name_delimiter() const {
            return this->population_name_delimiter_;
        }
        bool using_yaml_data() const {
            return this->using_yaml_data_;
        }
        void set_using_yaml_data(bool b) {
            this->using_yaml_data_ = b;
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
        void set_constant_sites_removed(const bool b) {
            this->constant_sites_removed_ = b;
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
        const std::string& get_tree_log_path() const {
            return this->tree_log_path_;
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
        bool constrain_state_frequencies() const {
            return ((this->freq_1_settings.is_fixed()) &&
                    (this->freq_1_settings.get_value() == 0.5));
        }

        void update_operator_weights() {
            this->tree_model_settings.update_operator_weights(this->operator_settings);

            if (mutation_rate_settings.is_fixed()) {
                this->operator_settings->tunable_operators.at("MuRateScaler").turn_off();
            }
            if (freq_1_settings.is_fixed()) {
                this->operator_settings->tunable_operators.at("StateFreqMover").turn_off();
                this->operator_settings->tunable_operators.at("StateFreqDirichletOperator").turn_off();
            }

            const bool root_height_is_fixed = tree_model_settings.tree_prior->root_height_is_fixed();
            const bool pop_sizes_fixed = population_size_settings.is_fixed();
            const bool pop_sizes_constrained = population_size_settings.population_sizes_are_constrained();

            if (root_height_is_fixed || pop_sizes_constrained || pop_sizes_fixed) {
                this->operator_settings->tunable_operators.at("RootHeightSizeMixer").turn_off();
            }
            if (root_height_is_fixed || pop_sizes_fixed) {
                this->operator_settings->tunable_operators.at("GlobalHeightSizeMixer").turn_off();
            }
            if (pop_sizes_fixed || pop_sizes_constrained) {
                this->operator_settings->tunable_operators.at("HeightSizeMixer").turn_off();
                this->operator_settings->tunable_operators.at("HeightSizeSlideBumpMixer").turn_off();
            }
            if (pop_sizes_fixed) {
                this->operator_settings->tunable_operators.at("GlobalPopSizeScaler").turn_off();
                this->operator_settings->tunable_operators.at("PopSizeScaler").turn_off();
            }

            if (pop_sizes_constrained && (! pop_sizes_fixed)) {
                // If pop sizes are constrained, we don't want to use
                // PopSizeScaler class' default for the weight, because it is
                // set assuming each branch has an independent pop size (i.e.,
                // the weight would be way too large).
                // If sizes are constrained, PopSizeScaler and
                // GlobalPopSizeScaler are equivalent, so we'll use the latter.
                this->operator_settings->tunable_operators.at("PopSizeScaler").turn_off();
                std::vector<double> global_size_scaler_wts = this->operator_settings->tunable_operators.at("GlobalPopSizeScaler").get_weights();
                bool global_size_scaler_off = true;
                if (global_size_scaler_wts.empty()) {
                    global_size_scaler_off = false;
                }
                else {
                    for (auto w : global_size_scaler_wts) {
                        if (w != 0.0) {
                            global_size_scaler_off = false;
                            break;
                        }
                    }
                }
                if (global_size_scaler_off) {
                    // Make sure we use default weight for GlobalPopSizeScaler
                    this->operator_settings->tunable_operators.at("GlobalPopSizeScaler").clear();
                }
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
        std::string tree_log_path_ = "phycoeval-trees-run-1.nex";
        std::string state_log_path_ = "phycoeval-state-run-1.log";
        std::string operator_log_path_ = "phycoeval-operator-run-1.log";

        void set_output_paths_to_config_directory() {
            std::pair<std::string, std::string> prefix_ext = path::splitext(this->config_path_);
            this->tree_log_path_ = prefix_ext.first + "-trees-run-1.nex";
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
                this->tree_model_settings.update_from_config(
                        top_level_node["tree_model"],
                        this->config_path_);
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
