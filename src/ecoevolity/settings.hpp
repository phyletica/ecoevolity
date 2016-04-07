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

class ContinuousDistributionSettings {

    private:
        std::string name_;
        std::map<std::string, double> parameters;

    public:
        const std::string& get_name() const {
            return this->name_;
        }

        ContinuousProbabilityDistribution get_instance() const {
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
        double get_value() const {
            return this->value_;
        }
        bool is_fixed() const {
            return this->is_fixed_;
        }

        PositiveRealParameter get_instance() const {
        }

};

class ComparisonSettings {

    friend class CollectionSettings;

    private:
        std::string path_;

        bool markers_are_dominant_;
        bool genotypes_are_diploid_;
        char population_name_delimiter_;
        bool population_name_is_prefix_;
        bool constant_sites_included_;

        bool constrain_population_sizes_;
        bool constrain_mutation_rates_;

        PositiveRealParameterSettings population_size_;
        PositiveRealParameterSettings u_;
        PositiveRealParameterSettings v_;
        PositiveRealParameterSettings time_multiplier_;
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
