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

#ifndef ECOEVOLITY_OPERATOR_SCHEDULE_HPP
#define ECOEVOLITY_OPERATOR_SCHEDULE_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <memory>

#include "rng.hpp"
#include "assert.hpp"
#include "settings.hpp"

class Operator;

class OperatorSchedule {
    protected:
        std::vector< std::shared_ptr<Operator> > operators_;
        double total_weight_ = 0.0;
        std::vector<double> cumulative_probs_;
        unsigned int auto_optimize_delay_ = 10000;
        unsigned int auto_optimize_delay_count_ = 0;
        bool auto_optimize_ = true;

    public:
        OperatorSchedule() { }
        OperatorSchedule(
                const OperatorScheduleSettings& settings,
                bool use_dpp = true);
        virtual ~OperatorSchedule() { }

        void add_operator(std::shared_ptr<Operator> o);

        std::shared_ptr<Operator>& draw_operator(RandomNumberGenerator& rng);

        double calc_delta(std::shared_ptr<const Operator> op, double log_alpha);

        double get_total_weight() const;
        unsigned int get_auto_optimize_delay_count() const;
        unsigned int get_auto_optimize_delay() const;
        void set_auto_optimize_delay(unsigned int delay);

        void write_operator_rates(std::ostream& out) const;

        bool auto_optimizing() const;
        void turn_on_auto_optimize();
        void turn_off_auto_optimize();

        bool using_dpp() const;
};

#endif
