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

#include "error.hpp"

class ContinuousDistributionSettings {

    private:
        std::string name_;
        std::map<std::string, double> parameters_;

    public:
        ContinuousDistributionSettings(
                const std::string& name,
                const std::unordered_map<std::string, double>& parameters) {
            this->name_ = name;
            // Gamma
            if (name == "gamma_distribution") {
                if (parameters.count("shape") < 1) {
                    throw ContinuousDistributionSettingError(
                            "shape parameter missing for gamma_distribution"
                            );
                }
                if (parameters.count("scale") < 1) {
                    throw ContinuousDistributionSettingError(
                            "scale parameter missing for gamma_distribution"
                            );
                }
                if (parameters.count("offset") < 1) {
                    if (parameters.size() > 2) {
                        throw ContinuousDistributionSettingError(
                                "unrecognized parameters for gamma_distribution (recognized parameters: shape, scale, offset)"
                                );
                    }
                }
                else {
                    if (parameters.size() > 3) {
                        throw ContinuousDistributionSettingError(
                                "unrecognized parameters for gamma_distribution (recognized parameters: shape, scale, offset)"
                                );
                    }
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
                    throw ContinuousDistributionSettingError(
                            "rate parameter missing for exponential_distribution"
                            );
                }
                if (parameters.count("offset") < 1) {
                    if (parameters.size() > 1) {
                        throw ContinuousDistributionSettingError(
                                "unrecognized parameters for exponential_distribution (recognized parameters: rate, offset)"
                                );
                    }
                }
                else {
                    if (parameters.size() > 2) {
                        throw ContinuousDistributionSettingError(
                                "unrecognized parameters for exponential_distribution (recognized parameters: rate, offset)"
                                );
                    }
                }
                this->parameters_["rate"] = parameters.at("rate");
                if (parameters.count("offset") == 1) {
                    this->parameters_["offset"] = parameters.at("offset");
                }
            }
            // Uniform
            else if (name == "uniform_distribution") {
                if (parameters.count("min") < 1) {
                    throw ContinuousDistributionSettingError(
                            "min parameter missing for uniform_distribution"
                            );
                }
                if (parameters.count("max") < 1) {
                    throw ContinuousDistributionSettingError(
                            "max parameter missing for uniform_distribution"
                            );
                }
                if (parameters.size() > 2) {
                    throw ContinuousDistributionSettingError(
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
                throw ContinuousDistributionSettingError(message);
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
            if (this->name_ == "gamma_distribution") {
                if (this->parameters_.count("offset") > 0) { 
                    return std::make_shared<OffsetGammaDistribution>(
                            this->parameters_.at("shape"),
                            this->parameters_.at("scale"),
                            this->parameters_.at("offset"));
                }
                else {
                    return std::make_shared<GammaDistribution>(
                            this->parameters_.at("shape"),
                            this->parameters_.at("scale"));
                }
            }
            else if (this->name_ == "exponential_distribution") {
                if (this->parameters_.count("offset") > 0) { 
                    return std::make_shared<OffsetExponentialDistribution>(
                            this->parameters_.at("rate"),
                            this->parameters_.at("offset"));
                }
                else {
                    return std::make_shared<ExponentialDistribution>(
                            this->parameters_.at("rate"));
                }
            }
            else if (this->name_ == "uniform_distribution") {
                return std::make_shared<UniformDistribution>(
                        this->parameters_.at("min"),
                        this->parameters_.at("max"));
            }
            else if (this->name_ == "none") {
                return std::make_shared<ContinuousProbabilityDistribution>();
            }
            else {
                ECOEVOLITY_ASSERT(0 == 1);
            }
        }

        std::string to_string(unsigned int indent_level = 0) const {
            if (this->name_ == "none") {
                return "\n";
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
        PositiveRealParameterSettings(
                double value,
                bool fix,
                const std::string& prior_name
                const std::unordered_map<std::string, double>& prior_parameters) {
            if (value < 0.0) {
                throw PositiveRealParameterSettingError(
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

        PositiveRealParameter get_instance() const {
            PositiveRealParameter p = PositiveRealParameter(
                    this->prior_settings_.get_instance();
                    this->value_;
                    this->is_fixed_;
                    );
            return p;
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
        bool constant_sites_included_;

        bool use_empirical_mutation_rate_starting_values_;
        bool constrain_population_sizes_;
        bool constrain_mutation_rates_;

        void make_consistent() {
            if (this->constrain_mutation_rates_) {
                this->use_empirical_mutation_rate_starting_values_ = false;
                std::unordered_map<std::string, double> prior_parameters
                this->u_settings_.prior_settings_.nullify();
                this->v_settings_.prior_settings_.nullify();
            }
            if (this->constrain_population_sizes_) {
                this->population_size_settings_.is_fixed_ = true;
                this->population_size_settings_.prior_settings_.nullify();
            }
        }

    public:
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
                bool constant_sites_included = false,
                bool use_empirical_mutation_rate_starting_values = true;
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
            this->constant_sites_included_ = constant_sites_included;
            this->use_empirical_mutation_rate_starting_values_ = use_empirical_mutation_rate_starting_values;
            this->constrain_population_sizes_ = constrain_population_sizes;
            this->constrain_mutation_rates_ = constrain_mutation_rates;
            this->make_consistent();
            // TODO: need to handle empirical mutation rates here.
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
            this->constant_sites_included_                      = other.constant_sites_included_;
            this->use_empirical_mutation_rate_starting_values_  = other.use_empirical_mutation_rate_starting_values_;
            this->constrain_population_sizes_                   = other.constrain_population_sizes_;
            this->constrain_mutation_rates_                     = other.constrain_mutation_rates_;
            return * this;
        }

        ComparisonPopulationTree get_instance() const {
            // TODO:
        }
};

class CollectionSettings {

    private:

        bool use_dpp_ = true;
        unsigned int chain_length_;
        unsigned int sample_every_;

        ContinuousDistributionSettings time_prior_settings_;

        ParameterSettings concentration_;
};

#endif
