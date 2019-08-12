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
#include "settings.hpp"

template<class NodeType>
class GeneralTreeOperatorSchedule {
    protected:
        std::vector< std::shared_ptr< GeneralTreeOperatorInterface<NodeType> > > operators_;
        double total_weight_ = 0.0;
        std::vector<double> cumulative_probs_;

    public:
        GeneralTreeOperatorSchedule() { }
        // GeneralTreeOperatorSchedule(const GeneralTreeSettings& settings) { }
        virtual ~GeneralTreeOperatorSchedule() { }

        void add_operator(std::shared_ptr< GeneralTreeOperatorInterface<NodeType> > o) {
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

        GeneralTreeOperatorInterface<NodeType>& draw_operator(RandomNumberGenerator& rng) const {
            double u = rng.uniform_real();
            for (unsigned int i = 0; i < this->cumulative_probs_.size(); ++i) {
                if (u <= this->cumulative_probs_.at(i)) {
                    return *this->operators_.at(i);
                }
            }
            return *this->operators_.back();
        }

        GeneralTreeOperatorInterface<NodeType>& get_operator(unsigned int operator_index) const {
            return *this->operators_.at(operator_index);
        }

        std::vector< std::shared_ptr< GeneralTreeOperatorInterface<NodeType> > > get_node_height_operators() const {
            std::vector< std::shared_ptr< GeneralTreeOperatorInterface<NodeType> > > ops;
            for (std::shared_ptr< GeneralTreeOperatorInterface<NodeType> > op: this->operators_) {
                if (op->get_type() == GeneralTreeOperatorInterface::OperatorTypeEnum::node_height_operator) {
                    ops.push_back(op);
                }
            }
            return ops;
        }

        double get_total_weight() const {
            return this->total_weight_;
        }

        void write_operator_rates(std::ostream& out) const {
            OperatorInterface& op = this->get_operator(0);
            out << op.header_string();
            out << op.to_string(*this);
            for (unsigned int i = 1; i < this->operators_.size(); ++i) {
                out << this->get_operator(i).to_string(*this);
            }
            out << std::flush;
        }
};

#endif
