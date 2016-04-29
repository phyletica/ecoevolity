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

#include "operator_schedule.hpp"
#include "operator.hpp"


OperatorSchedule::OperatorSchedule(
        const OperatorScheduleSettings& settings,
        bool use_dpp = true) {
    this->turn_off_auto_optimize();
    if (settings.auto_optimizing()) {
        this->turn_on_auto_optimize();
    }
    this->set_auto_optimize_delay(settings.get_auto_optimize_delay());

    if (settings.get_model_operator_settings().get_weight() > 0.0) {
        if (use_dpp) {
            this->add_operator(std::make_shared<DirichletProcessGibbsSampler>(
                    settings.get_model_operator_settings().get_weight()));
        }
        else {
            this->add_operator(std::make_shared<ReversibleJumpSampler>(
                    settings.get_model_operator_settings().get_weight()));
        }
    }

    if (settings.get_concentration_scaler_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<ConcentrationScaler>(
                settings.get_concentration_scaler_settings().get_weight(),
                settings.get_concentration_scaler_settings().get_scale(),
                ));
    }

    if (settings.get_comparison_height_scaler_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<ComparisonHeightScaler>(
                settings.get_comparison_height_scaler_settings().get_weight(),
                settings.get_comparison_height_scaler_settings().get_scale(),
                ));
    }

    if (settings.get_comparison_height_multiplier_scaler_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<ComparisonHeightMultiplierScaler>(
                settings.get_comparison_height_multiplier_scaler_settings().get_weight(),
                settings.get_comparison_height_multiplier_scaler_settings().get_scale(),
                ));
    }

    if (settings.get_root_coalescence_rate_scaler_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<RootCoalescenceRateScaler>(
                settings.get_root_coalescence_rate_scaler_settings().get_weight(),
                settings.get_root_coalescence_rate_scaler_settings().get_scale(),
                ));
    }

    if (settings.get_child_coalescence_rate_scaler_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<ChildCoalescenceRateScaler>(
                settings.get_child_coalescence_rate_scaler_settings().get_weight(),
                settings.get_child_coalescence_rate_scaler_settings().get_scale(),
                ));
    }

    if (settings.get_mutation_rate_mover_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<MutationRateMover>(
                settings.get_mutation_rate_mover_settings().get_weight(),
                settings.get_mutation_rate_mover_settings().get_window(),
                ));
    }
}

void OperatorSchedule::add_operator(std::shared_ptr<Operator> o) {
    this->operators_.push_back(o);
    o->set_operator_schedule(this);
    this->total_weight_ += o->get_weight();
    this->cumulative_probs_.push_back(0.0);
    ECOEVOLITY_ASSERT(this->operators_.size() == this->cumulative_probs_.size());
    this->cumulative_probs_.at(0) = this->operators_.at(0).get_weight() / this->total_weight_;
    for (unsigned int i = 1; i < this->operators_.size(); ++i) {
        this->cumulative_probs_.at(i) =
                (this->operators_.at(i).get_weight() /
                this->total_weight_) + 
                this->cumulative_probs_.at(i - 1);
    }
}

Operator& OperatorSchedule::draw_operator(RandomNumberGenerator& rng) {
    double u = this->rng.uniform_real();
    for (unsigned int i = 0; i < this->cumulative_probs_.size(); ++i) {
        if (u <= this->cumulative_probs_.at(i)) {
            return this->operators_.at(i);
        }
    }
    return this->operators_.back();
}

double OperatorSchedule::calc_delta(const Operator& op, double log_alpha) {
    if ((this->get_auto_optimize_delay_count() < this->get_auto_optimize_delay()) ||
            (! this->auto_optimize_)) {
        return 0.0;
    }
    double target = op.get_target_acceptance_probability();
    double count = (op.get_number_rejected_for_correction() +
                    op.get_number_accepted_for_correction() +
                    1.0);
    double delta_p = ((1.0 / count) * (std::exp(std::min(log_alpha, 0)) - target));
    double mx = std::numeric_limits<double>::max();
    if ((delta_p > -mx) && (delta_p < mx)) {
        return delta_p;
    }
    return 0.0;
}

double OperatorSchedule::get_total_weight() const {
    return this->total_weight_;
}
unsigned int OperatorSchedule::get_auto_optimize_delay_count() const {
    return this->auto_optimize_delay_count_;
}
unsigned int OperatorSchedule::get_auto_optimize_delay() const {
    return this->auto_optimize_delay_;
}
void OperatorSchedule::set_auto_optimize_delay(unsigned int delay) {
    this->auto_optimize_delay_ = delay;
}

void OperatorSchedule::write_operator_rates(std::ofstream out) {
    const Operator& op = this->operators_.at(0);
    out << op.header_string();
    out << op.to_string();
    for (unsigned int i = 1; i < this->operators_.size(); ++i) {
        out << this->operators_.at(i).to_string();
    }
}

bool OperatorSchedule::auto_optimizing() const {
    return this->auto_optimize_;
}
void OperatorSchedule::turn_on_auto_optimize() {
    this->auto_optimize_ = true;
}
void OperatorSchedule::turn_off_auto_optimize() {
    this->auto_optimize_ = false;
}