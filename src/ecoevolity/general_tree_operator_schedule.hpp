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

#ifndef ECOEVOLITY_GENERAL_TREE_OPERATOR_SCHEDULE_HPP
#define ECOEVOLITY_GENERAL_TREE_OPERATOR_SCHEDULE_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <memory>

#include "general_tree_operator.hpp"
#include "rng.hpp"
#include "assert.hpp"

template<class TreeType>
class GeneralTreeOperatorSchedule {
    protected:
        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > operators_;
        double total_weight_ = 0.0;
        std::vector<double> cumulative_probs_;
        unsigned int default_auto_optimize_delay_ = 50;

    public:
        GeneralTreeOperatorSchedule() { }
        GeneralTreeOperatorSchedule(
                std::shared_ptr<GeneralTreeOperatorSettingsCollection> settings,
                unsigned int number_of_leaves) {
            for (auto op_settings : settings->untunable_operators) {
                this->_add_untunable_operators(
                        op_settings.first,
                        op_settings.second,
                        number_of_leaves);
            }
            for (auto op_settings : settings->tunable_operators) {
                this->_add_tunable_operators(
                        op_settings.first,
                        op_settings.second,
                        number_of_leaves);
            }
            this->provide_split_lump_rj_move_with_helper_ops();
        }
        virtual ~GeneralTreeOperatorSchedule() { }

        unsigned int get_number_of_operators() const {
            return this->operators_.size();
        }

        void add_operator(std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > o) {
            this->operators_.push_back(o);
            this->total_weight_ += o->get_weight();
            this->cumulative_probs_.push_back(0.0);
            ECOEVOLITY_ASSERT(this->operators_.size() == this->cumulative_probs_.size());
            this->cumulative_probs_.at(0) = this->operators_.at(0)->get_weight() / this->total_weight_;
            for (unsigned int i = 1; i < this->operators_.size(); ++i) {
                this->cumulative_probs_.at(i) =
                        (this->operators_.at(i)->get_weight() /
                        this->total_weight_) + 
                        this->cumulative_probs_.at(i - 1);
            }
        }

        std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > draw_operator(
                RandomNumberGenerator& rng) const {
            double u = rng.uniform_real();
            for (unsigned int i = 0; i < this->cumulative_probs_.size(); ++i) {
                if (u <= this->cumulative_probs_.at(i)) {
                    return this->operators_.at(i);
                }
            }
            return this->operators_.back();
        }

        std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > get_operator(
                unsigned int operator_index) const {
            return this->operators_.at(operator_index);
        }

        bool append_operators(
                const std::string & op_name,
                std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > & ops
                ) const {
            bool found = false;
            for (std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > op: this->operators_) {
                if (op->get_name() == op_name) {
                    ops.push_back(op);
                    found = true;
                }
            }
            return found;
        }

        bool get_operators(
                const std::string & op_name,
                std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > & ops
                ) const {
            bool found = false;
            for (std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > op_iter: this->operators_) {
                if (op_iter->get_name() == op_name) {
                    ops.push_back(op_iter);
                    found = true;
                }
            }
            return found;
        }

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > get_operators(
                const std::string & op_name) const {
            std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > return_ops;
            this->get_operators(op_name, return_ops);
            return return_ops;
        }

        bool get_operators(
                const BaseGeneralTreeOperatorTemplate::OperatorTypeEnum op_type,
                std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > & ops
                ) const {
            bool found = false;
            for (std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > op: this->operators_) {
                if (op->get_type() == op_type) {
                    ops.push_back(op);
                    found = true;
                }
            }
            return found;
        }

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > get_operators(
                const BaseGeneralTreeOperatorTemplate::OperatorTypeEnum op_type
                ) const {
            std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > ops;
            this->get_operators(op_type, ops);
            return ops;
        }

        bool get_operators(
                const BaseGeneralTreeOperatorTemplate::OperatorScopeEnum op_scope,
                std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > & ops
                ) const {
            bool found = false;
            for (std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > op: this->operators_) {
                if (op->get_scope() == op_scope) {
                    ops.push_back(op);
                    found = true;
                }
            }
            return found;
        }

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > get_operators(
                const BaseGeneralTreeOperatorTemplate::OperatorScopeEnum op_scope
                ) const {
            std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > ops;
            this->get_operators(op_scope, ops);
            return ops;
        }

        bool get_split_lump_rj_operators(
                std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > & rj_ops
                ) const {
            return this->get_operators("SplitLumpNodesRevJumpSampler", rj_ops);
        }

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > get_split_lump_rj_operators() const {
            std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > return_ops;
            this->get_split_lump_rj_operators(return_ops);
            return return_ops;
        }

        bool get_root_height_operators(
                std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > & ops
                ) const {
            return this->get_operators(
                    BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::root_height_operator,
                    ops);
        }

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > get_root_height_operators() const {
            std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > ops;
            this->get_root_height_operators(ops);
            return ops;
        }

        void get_node_height_operators(
                std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > & ops
                ) const {
            this->get_root_height_operators(ops);
            this->get_operators(
                    BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::node_height_operator,
                    ops);
            this->get_operators(
                    BaseGeneralTreeOperatorTemplate::OperatorTypeEnum::global_height_operator,
                    ops);
        }

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > get_node_height_operators() const {
            std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > ops;
            this->get_node_height_operators(ops);
            return ops;
        }

        void get_preferred_node_height_operators(
                std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > & ops) const {
            bool got_global_dir_op = false;
            got_global_dir_op = this->get_operators("GlobalNodeHeightDirichletOperator", ops);
            if (got_global_dir_op) {
                this->get_root_height_operators(ops);
                return;
            }
            this->get_node_height_operators(ops);
        }

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > get_preferred_node_height_operators() const {
            std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > ops;
            this->get_preferred_node_height_operators(ops);
            return ops;
        }

        void provide_split_lump_rj_move_with_helper_ops() {
            std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > rj_ops;
            this->get_split_lump_rj_operators(rj_ops);
            if (rj_ops.empty()) {
                return;
            }
            std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > ops;
            this->get_preferred_node_height_operators(ops);
            for (auto rj_op: rj_ops) {
                rj_op->helper_ops = ops;
            }
        }

        double get_total_weight() const {
            return this->total_weight_;
        }

        void write_operator_rates(std::ostream& out) const {
            std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > op = this->get_operator(0);
            out << op->header_string();
            out << op->to_string();
            for (unsigned int i = 1; i < this->operators_.size(); ++i) {
                out << this->get_operator(i)->to_string();
            }
            out << std::flush;
        }

        std::set<std::string> write_op_settings(
                std::ostream & out,
                const unsigned int indent_level = 0) const {
            std::set<std::string> op_names;
            std::map< std::string, std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > > op_map;
            for (auto op : this->operators_) {
                op_names.insert(op->get_name());
                op_map[op->get_name()].push_back(op);
            }

            out << std::boolalpha;
            std::string margin = string_util::get_indent(indent_level);
            std::string indent = string_util::get_indent(1);
            for (auto name_ops : op_map) {
                out << margin << name_ops.first << ":\n";

                out << margin << indent << "weight: ";
                if (name_ops.second.size() == 1) {
                    out << name_ops.second.at(0)->get_weight() << "\n";
                }
                else {
                    out << "[" << name_ops.second.at(0)->get_weight();
                    for (unsigned int i = 1; i < name_ops.second.size(); ++i) {
                        out << ", " << name_ops.second.at(i)->get_weight();
                    }
                    out << "]\n";
                }

                if (std::isnan(name_ops.second.at(0)->get_coercable_parameter_value())) {
                    continue;
                }

                out << margin << indent << "tuning_parameter: ";
                if (name_ops.second.size() == 1) {
                    out << name_ops.second.at(0)->get_coercable_parameter_value() << "\n";
                }
                else {
                    out << "[" << name_ops.second.at(0)->get_coercable_parameter_value();
                    for (unsigned int i = 1; i < name_ops.second.size(); ++i) {
                        out << ", " << name_ops.second.at(i)->get_coercable_parameter_value();
                    }
                    out << "]\n";
                }

                out << margin << indent << "auto_optimize: ";
                if (name_ops.second.size() == 1) {
                    out << name_ops.second.at(0)->auto_optimizing() << "\n";
                }
                else {
                    out << "[" << name_ops.second.at(0)->auto_optimizing();
                    for (unsigned int i = 1; i < name_ops.second.size(); ++i) {
                        out << ", " << name_ops.second.at(i)->auto_optimizing();
                    }
                    out << "]\n";
                }

                out << margin << indent << "auto_optimize_delay: ";
                if (name_ops.second.size() == 1) {
                    out << name_ops.second.at(0)->get_auto_optimize_delay() << "\n";
                }
                else {
                    out << "[" << name_ops.second.at(0)->get_auto_optimize_delay();
                    for (unsigned int i = 1; i < name_ops.second.size(); ++i) {
                        out << ", " << name_ops.second.at(i)->get_auto_optimize_delay();
                    }
                    out << "]\n";
                }
            }
            return op_names;
        }


    protected:
        void _add_untunable_operators(
                const std::string & op_name,
                GeneralTreeOperatorSettings op_settings,
                unsigned int number_of_leaves) {
            if (op_settings.is_empty()) {
                this->_add_untunable_op(
                        op_name,
                        -1.0,
                        1,
                        number_of_leaves);
                return;
            }
            for (unsigned int i = 0; i < op_settings.get_number_of_operators(); ++i) {
                this->_add_untunable_op(
                        op_name,
                        op_settings.get_weight(i),
                        op_settings.get_number_of_operators(),
                        number_of_leaves);
            }
        }

        void _add_untunable_op(
                const std::string & op_name,
                double weight,
                unsigned int operator_count,
                unsigned int number_of_leaves) {
            if (weight == 0.0) {
                return;
            }
            std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > op;
            if (op_name == "NeighborHeightNodePermute") {
                op = std::make_shared<
                                 NeighborHeightNodePermute<TreeType>
                                      >();
            }
            else if (op_name == "NeighborHeightNodeSwap") {
                op = std::make_shared<
                                 NeighborHeightNodeSwap<TreeType>
                                      >();
            }
            else if (op_name == "NeighborHeightNodeSwapAll") {
                op = std::make_shared<
                                 NeighborHeightNodeSwapAll<TreeType>
                                      >();
            }
            else {
                throw EcoevolityError(
                        "GeneralTreeOperatorSchedule: Unrecognized untunable operator: " + op_name);
            }

            if (weight < 0.0) {
                double wt = op->get_default_weight(number_of_leaves);
                if (wt <= 0.0) {
                    return;
                }
                op->set_weight(wt / (double)operator_count);
            }
            else {
                op->set_weight(weight);
            }

            ECOEVOLITY_ASSERT(op->get_weight() > 0.0);

            this->add_operator(op);
        }

        void _add_tunable_operators(
                const std::string & op_name,
                GeneralTreeTunableOperatorSettings op_settings,
                unsigned int number_of_leaves) {
            if (op_settings.is_empty()) {
                if (op_name ==  "SplitLumpNodesRevJumpSampler") {
                    // Special default case of multiple operators for
                    // SplitLumpNodesRevJumpSampler
                    unsigned int n_rj_ops = 3;
                    this->_add_tunable_op(
                            op_name,
                            -1.0,
                            1.0,
                            false,
                            -1,
                            n_rj_ops,
                            number_of_leaves);
                    this->_add_tunable_op(
                            op_name,
                            -1.0,
                            10.0,
                            false,
                            -1,
                            n_rj_ops,
                            number_of_leaves);
                    this->_add_tunable_op(
                            op_name,
                            -1.0,
                            1.0,
                            true,
                            -1,
                            n_rj_ops,
                            number_of_leaves);
                    return;
                }
                this->_add_tunable_op(
                        op_name,
                        -1.0,
                        -1.0,
                        true,
                        -1,
                        1,
                        number_of_leaves);
                return;
            }
            for (unsigned int i = 0; i < op_settings.get_number_of_operators(); ++i) {
                this->_add_tunable_op(
                        op_name,
                        op_settings.get_weight(i),
                        op_settings.get_tuning_parameter(i),
                        op_settings.auto_optimizing(i),
                        op_settings.get_auto_optimize_delay(i),
                        op_settings.get_number_of_operators(),
                        number_of_leaves);
            }
        }

        void _add_tunable_op(
                const std::string & op_name,
                double weight,
                double tuning_parameter,
                bool auto_optimize,
                int auto_optimize_delay,
                unsigned int operator_count,
                unsigned int number_of_leaves) {
            if (weight == 0.0) {
                return;
            }
            std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > op;
            if (op_name == "SplitLumpNodesRevJumpSampler") {
                op = std::make_shared<
                            SplitLumpNodesRevJumpSampler<TreeType>
                                      >();
            }
            else if (op_name == "TreeScaler") {
                op = std::make_shared<
                            TreeScaler<TreeType>
                                      >();
            }
            else if (op_name == "NodeHeightScaler") {
                op = std::make_shared<
                                 NodeHeightScaler<TreeType>
                                      >();
            }
            else if (op_name == "NodeHeightMover") {
                op = std::make_shared<
                                 NodeHeightMover<TreeType>
                                      >();
            }
            else if (op_name == "NodeHeightSlideBumpScaler") {
                op = std::make_shared<
                                 NodeHeightSlideBumpScaler<TreeType>
                                      >();
            }
            else if (op_name == "NodeHeightSlideBumpPermuteScaler") {
                op = std::make_shared<
                                 NodeHeightSlideBumpPermuteScaler<TreeType>
                                      >();
            }
            else if (op_name == "NodeHeightSlideBumpSwapScaler") {
                op = std::make_shared<
                                 NodeHeightSlideBumpSwapScaler<TreeType>
                                      >();
            }
            else if (op_name == "NodeHeightSlideBumpMover") {
                op = std::make_shared<
                                 NodeHeightSlideBumpMover<TreeType>
                                      >();
            }
            else if (op_name == "NodeHeightSlideBumpPermuteMover") {
                op = std::make_shared<
                                 NodeHeightSlideBumpPermuteMover<TreeType>
                                      >();
            }
            else if (op_name == "NodeHeightSlideBumpSwapMover") {
                op = std::make_shared<
                                 NodeHeightSlideBumpSwapMover<TreeType>
                                      >();
            }
            else if (op_name == "RootHeightScaler") {
                op = std::make_shared<
                                 RootHeightScaler<TreeType>
                                      >();
            }
            else if (op_name == "GlobalNodeHeightDirichletOperator") {
                op = std::make_shared<
                                 GlobalNodeHeightDirichletOperator<TreeType>
                                      >();
            }
            else if (op_name == "NodeHeightDirichletOperator") {
                op = std::make_shared<
                                 NodeHeightDirichletOperator<TreeType>
                                      >();
            }
            else if (op_name == "NodeHeightPriorAlphaScaler") {
                op = std::make_shared<
                                 NodeHeightPriorAlphaScaler<TreeType>
                                      >();
            }
            else if (op_name == "NodeHeightPriorAlphaMover") {
                op = std::make_shared<
                                 NodeHeightPriorAlphaMover<TreeType>
                                      >();
            }
            else if (op_name == "MuRateScaler") {
                op = std::make_shared<
                                 MuRateScaler
                                      >();
            }
            else if (op_name == "GlobalPopSizeScaler") {
                op = std::make_shared<
                                 GlobalPopSizeScaler
                                      >();
            }
            else if (op_name == "PopSizeScaler") {
                op = std::make_shared<
                                 PopSizeScaler
                                      >();
            }
            else if (op_name == "GlobalHeightSizeMixer") {
                op = std::make_shared<
                                 GlobalHeightSizeMixer
                                      >();
            }
            else if (op_name == "HeightSizeMixer") {
                op = std::make_shared<
                                 HeightSizeMixer
                                      >();
            }
            else if (op_name == "HeightSizeSlideBumpMixer") {
                op = std::make_shared<
                                 HeightSizeSlideBumpMixer
                                      >();
            }
            else if (op_name == "RootHeightSizeMixer") {
                op = std::make_shared<
                                 RootHeightSizeMixer
                                      >();
            }
            else if (op_name == "GlobalHeightSizeRateScaler") {
                op = std::make_shared<
                                 GlobalHeightSizeRateScaler
                                      >();
            }
            else if (op_name == "GlobalHeightSizeScaler") {
                op = std::make_shared<
                                 GlobalHeightSizeScaler
                                      >();
            }
            else if (op_name == "GlobalHeightRateScaler") {
                op = std::make_shared<
                                 GlobalHeightRateScaler
                                      >();
            }
            else if (op_name == "StateFreqMover") {
                op = std::make_shared<
                                 StateFreqMover
                                      >();
            }
            else if (op_name == "StateFreqDirichletOperator") {
                op = std::make_shared<
                                 StateFreqDirichletOperator
                                      >();
            }
            else {
                throw EcoevolityError(
                        "GeneralTreeOperatorSchedule: Unrecognized tunable operator: " + op_name);
            }

            if (weight < 0.0) {
                double wt = op->get_default_weight(number_of_leaves);
                if (wt <= 0.0) {
                    return;
                }
                op->set_weight(wt / (double)operator_count);
            }
            else {
                op->set_weight(weight);
            }

            op->turn_on_auto_optimize();
            if (! auto_optimize) {
                op->turn_off_auto_optimize();
            }

            if (tuning_parameter > 0.0) {
                op->set_coercable_parameter_value(tuning_parameter);
            }
            else {
                op->set_coercable_parameter_value(
                        op->get_default_coercable_parameter_value());
            }

            if (auto_optimize_delay < 0) {
                op->set_auto_optimize_delay(this->default_auto_optimize_delay_);
            }
            else {
                op->set_auto_optimize_delay(auto_optimize_delay);
            }

            ECOEVOLITY_ASSERT(op->get_weight() > 0.0);

            this->add_operator(op);
        }
};

#endif
