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

#include "operator.hpp"

Operator::Operator(double weight) {
    this->set_weight(weight);
}

double Operator::get_target_acceptance_probability() const {
    return 0.234;
}

void Operator::set_operator_schedule(std::shared_ptr<OperatorSchedule> os) {
    this->operator_schedule_ = os;
}

double Operator::get_coercable_parameter_value() {
    return std::numeric_limits<double>::quiet_NaN();
}

void Operator::set_weight(double weight) {
    ECOEVOLITY_ASSERT(weight > 0.0);
    this->weight_ = weight;
}

void Operator::accept() {
    ++this->number_accepted_;
    if (this->operator_schedule_->get_auto_optimize_delay_count() >= this->operator_schedule_->get_auto_optimize_delay()) {
        ++this->number_accepted_for_correction_;
    }
}

void Operator::reject() {
    ++this->number_rejected_;
    if (this->operator_schedule_->get_auto_optimize_delay_count() >= this->operator_schedule_->get_auto_optimize_delay()) {
        ++this->number_rejected_for_correction_;
    }
}

std::string Operator::get_name() const {
    return "BaseOperator";
}

std::string Operator::header_string() const {
    return "name\tnumber_accepted\tnumber_rejected\tweight\tweight_prob\ttuning_parameter\n";
}

std::string Operator::to_string() const {
    std::ostringstream ss;
    ss << this->get_name() << "\t" 
       << this->get_number_accepted() << "\t"
       << this->get_number_rejected() << "\t"
       << this->get_weight() << "\t";

    if (this->operator_schedule_->get_total_weight() > 0.0) {
        ss << this->get_weight() / this->operator_schedule_->get_total_weight() << "\t";
    }
    else {
        ss << "\t";
    }

    double tuning = this->get_coercable_parameter_value();
    if (std::isnan(tuning)) {
        ss << "\t";
    }
    else {
        ss << tuning << "\t";
    }
    ss << "\n";
    return ss.str();
}

double Operator::calc_delta(double log_alpha) const {
    return this->operator_schedule_->calc_delta(this, log_alpha);
}


ScaleOperator::ScaleOperator(double weight, double scale) : Operator(weight) {
    this->set_scale(scale);
}

void ScaleOperator::set_scale(double scale) {
    ECOEVOLITY_ASSERT(scale > 0.0);
    this->scale_ = scale;
}

void ScaleOperator::update(
        RandomNumberGenerator& rng,
        double& parameter_value,
        double& hastings_ratio) const {
    double multiplier = std::exp(this->scale_ * ((2.0 * rng.uniform_real()) - 1.0));
    parameter_value *= multiplier
    hastings_ratio = std::log(multiplier);
}

void ScaleOperator::optimize(double log_alpha) {
    double delta = this->calc_delta(log_alpha);
    delta += std::log(this->scale_);
    this->set_scale(std::exp(delta));
}

double ScaleOperator::get_coercable_parameter_value() {
    return this->scale_;
}

void ScaleOperator::set_coercable_parameter_value(double value) {
    this->set_scale(value);
}

std::string ScaleOperator::get_name() const {
    return "BaseScaleOperator";
}
