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

#ifndef ECOEVOLITY_SETTINGS_IO_HPP
#define ECOEVOLITY_SETTINGS_IO_HPP


#include "general_tree_settings.hpp"
#include "general_tree_operator_schedule.hpp"

template<class TreeType>
inline void write_settings(
        std::ostream & out,
        const PopulationTreeAnalysisSettings & settings,
        const GeneralTreeOperatorSchedule<TreeType> & operator_schedule
        ) {
    std::string indent = string_util::get_indent(1);
    out << std::boolalpha;

    out << "---\n"
        << "data:\n"
        << settings.data_settings.to_string(1)
        << "tree_model:\n"
        << settings.tree_model_settings.to_string(1)
        << "branch_parameters:\n"
        << indent << "population_size:\n"
        << settings.population_size_settings.to_string(2)
        << "mutation_parameters:\n"
        << indent << "mutation_rate:\n"
        << settings.mutation_rate_settings.to_string(2)
        << indent << "freq_1:\n"
        << settings.freq_1_settings.to_string(2)
        << "mcmc_settings:\n"
        << indent << "chain_length: " << settings.get_chain_length() << "\n"
        << indent << "sample_frequency: " << settings.get_sample_frequency() << "\n"
        << indent << "operators:\n";
    std::set<std::string> op_names;
    op_names = operator_schedule.write_op_settings(out, 2);
    std::string margin = string_util::get_indent(2);
    for (auto op : settings.operator_settings->untunable_operators) {
        if (op_names.count(op.first) < 1) {
            out << margin << op.first << ":\n";
            out << margin << indent << "weight: 0\n";
        }
    }
    for (auto op : settings.operator_settings->tunable_operators) {
        if (op_names.count(op.first) < 1) {
            out << margin << op.first << ":\n";
            out << margin << indent << "weight: 0\n";
        }
    }
}

template<class TreeType>
inline void write_settings(
        std::ostream & out,
        const NucTreeAnalysisSettings & settings,
        const GeneralTreeOperatorSchedule<TreeType> & operator_schedule
        ) {
    std::string indent = string_util::get_indent(1);
    out << std::boolalpha;

    out << "---\n"
        << "data:\n"
        << settings.data_settings.to_string(1)
        << "tree_model:\n"
        << settings.tree_model_settings.to_string(1)
        << "mutation_parameters:\n"
        << indent << "state_frequencies:\n"
        << settings.state_freq_settings.to_string(2)
        << indent << "rate_matrix:\n"
        << settings.rate_matrix_settings.to_string(2)
        << indent << "among_site_rate_variation:\n"
        << indent << indent << "discrete_gamma:\n"
        << indent << indent << indent << "number_of_categories: " << settings.asrv_num_cats << "\n"
        << indent << indent << indent << "parameters:\n"
        << indent << indent << indent << indent << "one_over_shape:\n"
        << settings.asrv_one_over_shape_settings.to_string(5)
        << indent << indent << "proportion_invariable_sites:\n"
        << settings.asrv_prop_invar_settings.to_string(3)
        << "mcmc_settings:\n"
        << indent << "chain_length: " << settings.get_chain_length() << "\n"
        << indent << "sample_frequency: " << settings.get_sample_frequency() << "\n"
        << indent << "operators:\n";
    std::set<std::string> op_names;
    op_names = operator_schedule.write_op_settings(out, 2);
    std::string margin = string_util::get_indent(2);
    for (auto op : settings.operator_settings->untunable_operators) {
        if (op_names.count(op.first) < 1) {
            out << margin << op.first << ":\n";
            out << margin << indent << "weight: 0\n";
        }
    }
    for (auto op : settings.operator_settings->tunable_operators) {
        if (op_names.count(op.first) < 1) {
            out << margin << op.first << ":\n";
            out << margin << indent << "weight: 0\n";
        }
    }
}

#endif
