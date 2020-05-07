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

    public:
        GeneralTreeOperatorSchedule() { }
        GeneralTreeOperatorSchedule(
                std::shared_ptr<GeneralTreeOperatorSettingsCollection> settings,
                unsigned int number_of_leaves) {
            for (auto op_settings : settings->untunable_operators) {
                if (op_settings.second.get_weight() != 0.0) {
                    this->_add_untunable_op(
                            op_settings.first,
                            op_settings.second,
                            number_of_leaves);
                }
            }
            for (auto op_settings : settings->tunable_operators) {
                if (op_settings.second.get_weight() != 0.0) {
                    this->_add_tunable_op(
                            op_settings.first,
                            op_settings.second,
                            number_of_leaves);
                }
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

        void append_operator(
                const std::string & op_name,
                std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > & ops
                ) const {
            for (std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > op: this->operators_) {
                if (op->get_name() == op_name) {
                    ops.push_back(op);
                    return;
                }
            }
        }

        void get_operator(
                const std::string & op_name,
                std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > op) const {
            for (std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > op_iter: this->operators_) {
                if (op_iter->get_name() == op_name) {
                    op = op_iter;
                    return;
                }
            }
        }

        std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > get_operator(
                const std::string & op_name) const {
            std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > return_op;
            this->get_operator(op_name, return_op);
            return return_op;
        }

        void get_operators(
                const BaseGeneralTreeOperatorTemplate::OperatorTypeEnum op_type,
                std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > & ops
                ) const {
            for (std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > op: this->operators_) {
                if (op->get_type() == op_type) {
                    ops.push_back(op);
                }
            }
        }

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > get_operators(
                const BaseGeneralTreeOperatorTemplate::OperatorTypeEnum op_type
                ) const {
            std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > ops;
            this->get_operators(op_type);
            return ops;
        }

        void get_operators(
                const BaseGeneralTreeOperatorTemplate::OperatorScopeEnum op_scope,
                std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > & ops
                ) const {
            for (std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > op: this->operators_) {
                if (op->get_scope() == op_scope) {
                    ops.push_back(op);
                }
            }
        }

        std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > get_operators(
                const BaseGeneralTreeOperatorTemplate::OperatorScopeEnum op_scope
                ) const {
            std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > ops;
            this->get_operators(op_scope);
            return ops;
        }

        void get_split_lump_rj_operator(
                std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > rj_op) const {
            this->get_operator("SplitLumpNodesRevJumpSampler", rj_op);
        }

        std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > get_split_lump_rj_operator() const {
            std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > return_op;
            this->get_split_lump_rj_operator(return_op);
            return return_op;
        }

        void get_root_height_operators(
                std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > & ops
                ) const {
            this->get_operators(
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
                std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > ops) const {
            std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> >
                op = this->get_operator("GlobalNodeHeightDirichletOperator");
            if (op) {
                ops.push_back(op);
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
            std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > rj_op = this->get_split_lump_rj_operator();
            if (! rj_op) {
                return;
            }
            std::vector< std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > > ops;
            this->get_preferred_node_height_operators(ops);
            rj_op->helper_ops = ops;
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
            out << std::boolalpha;
            std::string margin = string_util::get_indent(indent_level);
            std::string indent = string_util::get_indent(1);
            for (auto op : this->operators_) {
                op_names.insert(op->get_name());
                out << margin << op->get_name() << ":\n"
                    << margin << indent
                    << "weight: " << op->get_weight() << "\n";
                if (std::isnan(op->get_coercable_parameter_value())) {
                    continue;
                }
                out << margin << indent
                    << "tuning_parameter: "
                    << op->get_coercable_parameter_value() << "\n"
                    << margin << indent
                    << "auto_optimize: " << op->auto_optimizing() << "\n";
                if (op->auto_optimizing()) {
                    out << margin << indent
                        << "auto_optimize_delay: "
                        << op->get_auto_optimize_delay() << "\n";
                }
            }
            return op_names;
        }


    protected:
        void _add_untunable_op(
                const std::string & op_name,
                GeneralTreeOperatorSettings op_settings,
                unsigned int number_of_leaves) {
            std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > op;
            if (op_name == "SplitLumpNodesRevJumpSampler") {
                op = std::make_shared<
                            SplitLumpNodesRevJumpSampler<TreeType>
                                      >();
            }
            else if (op_name == "NeighborHeightNodePermute") {
                op = std::make_shared<
                                 NeighborHeightNodePermute<TreeType>
                                      >();
            }
            else if (op_name == "NeighborHeightNodeSwap") {
                op = std::make_shared<
                                 NeighborHeightNodeSwap<TreeType>
                                      >();
            }
            else {
                throw EcoevolityError(
                        "GeneralTreeOperatorSchedule: Unrecognized untunable operator: " + op_name);
            }

            if (op_settings.get_weight() < 0.0) {
                op->set_default_weight(number_of_leaves);
            }
            else {
                op->set_weight(op_settings.get_weight());
            }

            this->add_operator(op);
        }
        void _add_tunable_op(
                const std::string & op_name,
                GeneralTreeTunableOperatorSettings op_settings,
                unsigned int number_of_leaves) {
            std::shared_ptr< GeneralTreeOperatorTemplate<TreeType> > op;
            if (op_name == "TreeScaler") {
                op = std::make_shared<
                            TreeScaler<TreeType> >();
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

            if (op_settings.get_weight() < 0.0) {
                op->set_default_weight(number_of_leaves);
            }
            else {
                op->set_weight(op_settings.get_weight());
            }

            op->set_auto_optimize_delay(op_settings.get_auto_optimize_delay());

            // If the tuning parameter is negative, it has not been set,
            // and we will keep the default value of the operator class
            if (op_settings.get_tuning_parameter() > 0.0) {
                op->set_coercable_parameter_value(
                        op_settings.get_tuning_parameter());
            }

            op->turn_on_auto_optimize();
            if (! op_settings.auto_optimizing()) {
                op->turn_off_auto_optimize();
            }

            this->add_operator(op);
        }
};

#endif
