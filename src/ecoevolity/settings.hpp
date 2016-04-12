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

#ifndef ECOEVOLITY_SETTINGS_HPP
#define ECOEVOLITY_SETTINGS_HPP

#include <iostream>
#include <sstream>
#include <map>
#include <unordered_map>

#include "error.hpp"
#include "assert.hpp"
#include "util.hpp"
#include "probability.hpp"
#include "parameter.hpp"
#include "tree.hpp"

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
                throw EcoevolityContinuousDistributionSettingError(
                        "asked for an instance of a null prior distribution");
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
            std::string margin = get_indent(indent_level);
            std::string indent = get_indent(1);
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
        double value_;
        bool is_fixed_;
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

        virtual std::shared_ptr<PositiveRealParameter> get_instance() const {
            if (this->prior_settings_.get_name() == "none") {
                return std::make_shared<PositiveRealParameter>(
                        this->value_,
                        this->is_fixed_
                        );
            }
            return std::make_shared<PositiveRealParameter>(
                    this->prior_settings_.get_instance(),
                    this->value_,
                    this->is_fixed_
                    );
        }

        virtual std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            std::string margin = get_indent(indent_level);
            if (! std::isnan(this->get_value())) {
                ss << margin << "value: " << this->get_value() << "\n";
            }
            ss << margin << "estimate: " << (! this->is_fixed()) << "\n";
            if (! this->is_fixed()) {
                ss << margin << "prior: " << (! this->is_fixed()) << "\n";
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
        PositiveRealParameterSettings v_settings_;
        PositiveRealParameterSettings time_multiplier_settings_;

        char population_name_delimiter_;
        bool population_name_is_prefix_;
        bool genotypes_are_diploid_;
        bool markers_are_dominant_;
        bool constant_sites_removed_;

        bool use_empirical_mutation_rate_starting_values_;
        bool constrain_population_sizes_;
        bool constrain_mutation_rates_;

        void make_consistent() {
            if (this->constrain_mutation_rates_) {
                this->use_empirical_mutation_rate_starting_values_ = false;
                this->u_settings_.prior_settings_.nullify();
                this->v_settings_.prior_settings_.nullify();
                this->u_settings_.value_ = 1.0;
                this->u_settings_.is_fixed_ = true;
                this->v_settings_.value_ = 1.0;
                this->v_settings_.is_fixed_ = true;
            }
        }

    public:
        ComparisonSettings() { }
        ComparisonSettings(
                const std::string& path,
                const PositiveRealParameterSettings& population_size_settings,
                const PositiveRealParameterSettings& u_settings,
                const PositiveRealParameterSettings& v_settings,
                const PositiveRealParameterSettings& time_multiplier_settings,
                char population_name_delimiter = '_',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = false,
                bool use_empirical_mutation_rate_starting_values = true,
                bool constrain_population_sizes = false,
                bool constrain_mutation_rates = true) {

            this->path_ = path;
            this->population_size_settings_ = population_size_settings;
            this->u_settings_ = u_settings;
            this->v_settings_ = v_settings;
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
            // more awkward than it's worth to defer it
            // Make a copy operator for BiallelicData and parse and store here, then
            // can copy it in get_instance method
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
                this->v_settings_.value_ = v;
            }
        }

        virtual ~ComparisonSettings() { }
        ComparisonSettings& operator=(const ComparisonSettings& other) {
            this->path_                                         = other.path_;
            this->population_size_settings_                     = other.population_size_settings_;
            this->u_settings_                                   = other.u_settings_;
            this->v_settings_                                   = other.v_settings_;
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

        ComparisonPopulationTree get_instance() const {
            ComparisonPopulationTree t = ComparisonPopulationTree(
                    this->path_,
                    this->population_name_delimiter_,
                    this->population_name_is_prefix_,
                    this->genotypes_are_diploid_,
                    this->markers_are_dominant_,
                    this->constant_sites_removed_,
                    true);
            if (this->constrain_mutation_rates_) {
                t.constrain_mutation_rates();
                t.fold_patterns();
            }

            if (this->constrain_population_sizes_) {
                t.constrain_coalescence_rates();
            }
            t.set_population_size_prior(this->population_size_settings_.prior_settings_.get_instance());
            t.set_coalescence_rate(
                    CoalescenceRateParameter::get_rate_from_population_size(
                            this->population_size_settings_.get_value()));
            if (this->population_size_settings_.is_fixed()) {
                t.fix_coalescence_rates();
            }

            t.set_u_parameter(this->u_settings_.get_instance());
            t.set_v_parameter(this->v_settings_.get_instance());
            t.set_node_height_multiplier_parameter(this->time_multiplier_settings_.get_instance());

            return t;
        }

        std::string to_string(unsigned int indent_level = 0) const {
            std::ostringstream ss;
            std::string margin = get_indent(indent_level);
            std::string indent = get_indent(1);
            ss << margin << "path: " << this->path_ << "\n";
            ss << margin << "genotypes_are_diploid: " << this->genotypes_are_diploid_ << "\n";
            ss << margin << "markers_are_dominant: " << this->markers_are_dominant_ << "\n";
            ss << margin << "population_name_delimiter: " << this->population_name_delimiter_ << "\n";
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

            ss << margin << indent << "v_rate:\n";
            ss << this->v_settings_.to_string(indent_level + 2);

            ss << margin << indent << "time_multiplier:\n";
            ss << this->time_multiplier_settings_.to_string(indent_level + 2);

            return ss.str();
        }
};

class CollectionSettings {

    private:

        bool use_dpp_ = true;
        unsigned int chain_length_;
        unsigned int sample_every_;

        ContinuousDistributionSettings time_prior_settings_;

        PositiveRealParameterSettings concentration_;

        std::vector<ComparisonSettings> comparisons_;

    public:

        void add_comparison();

        // TODO: CollectionSettings cs = CollectionSettings::init_from_config_file(path);
        static CollectionSettings init_from_config_file(const std::string& path);
};

#endif
