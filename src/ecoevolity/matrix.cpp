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

#include "matrix.hpp"

BiallelicPatternProbabilityMatrix::BiallelicPatternProbabilityMatrix(
        unsigned int allele_count) {
    this->reset(allele_count);
}

BiallelicPatternProbabilityMatrix::BiallelicPatternProbabilityMatrix(
                const BiallelicPatternProbabilityMatrix& matrix) {
    this->copy(matrix);
}

BiallelicPatternProbabilityMatrix::BiallelicPatternProbabilityMatrix(
                unsigned int allele_count,
                unsigned int red_allele_count) {
    ECOEVOLITY_ASSERT(red_allele_count <= allele_count);
    this->reset(allele_count);
    if (allele_count > 0) {
        this->set_pattern_probability(allele_count, red_allele_count, 1.0);
    }
    else {
        this->prob_missing_ = 1.0;
        this->pattern_prob_matrix_.clear();
        this->allele_count_ = allele_count;
    }
}

BiallelicPatternProbabilityMatrix::BiallelicPatternProbabilityMatrix(
                unsigned int allele_count,
                const std::vector<double>& pattern_probs) {
    if (allele_count == 0) {
        ECOEVOLITY_ASSERT(pattern_probs.size() == 0);
        this->prob_missing_ = 1.0;
        this->pattern_prob_matrix_.clear();
        this->allele_count_ = allele_count;
        return;
    }
    if (pattern_probs.size() != ((((allele_count + 1) * (allele_count + 2))/2) - 1)) {
        throw EcoevolityError(
                "Allele count and pattern probability vector size mismatch in BiallelicPatternProbabilityMatrix constructor");
    }
    this->allele_count_ = allele_count;
    this->pattern_prob_matrix_ = pattern_probs;
    this->prob_missing_ = 0.0;
}

void BiallelicPatternProbabilityMatrix::resize(unsigned int allele_count) {
    if (allele_count == 0) {
        this->prob_missing_ = 1.0;
        this->pattern_prob_matrix_.clear();
        this->allele_count_ = allele_count;
        return;
    }
    if (this->get_allele_count() == allele_count) {
        ECOEVOLITY_ASSERT(this->pattern_prob_matrix_.size() ==
                ((((allele_count + 1) * (allele_count + 2))/2) - 1));
        return;
    }
    this->allele_count_ = allele_count;
    this->pattern_prob_matrix_.resize(
            ((((allele_count + 1) * (allele_count + 2))/2) - 1),
            0);
    this->prob_missing_ = 0.0;
}

void BiallelicPatternProbabilityMatrix::reset(unsigned int allele_count) {
    if (allele_count == 0) {
        this->prob_missing_ = 1.0;
        this->pattern_prob_matrix_.clear();
        this->allele_count_ = allele_count;
        return;
    }
    if (this->get_allele_count() == allele_count) {
        ECOEVOLITY_ASSERT(this->pattern_prob_matrix_.size() ==
                ((((allele_count + 1) * (allele_count + 2))/2) - 1));
        std::fill(this->pattern_prob_matrix_.begin(), this->pattern_prob_matrix_.end(), 0);
        this->prob_missing_ = 0.0;
        return;
    }
    this->allele_count_ = allele_count;
    this->pattern_prob_matrix_.assign(
            ((((allele_count + 1) * (allele_count + 2))/2) - 1),
            0);
    this->prob_missing_ = 0.0;
}

void BiallelicPatternProbabilityMatrix::copy(
        const BiallelicPatternProbabilityMatrix& m) {
    this->resize(m.get_allele_count());
    ECOEVOLITY_ASSERT(this->get_allele_count() == m.get_allele_count());
    this->pattern_prob_matrix_ = m.get_pattern_prob_matrix();
    this->prob_missing_ = m.prob_missing_;
}

double BiallelicPatternProbabilityMatrix::get_pattern_probability(
        unsigned int allele_count,
        unsigned int red_allele_count) const {
    ECOEVOLITY_ASSERT(red_allele_count <= allele_count);
    if (allele_count == 0) {
        return this->prob_missing_;
    }
    int i = (((allele_count * (allele_count + 1))/2) - 1 + red_allele_count);
    if ((i < 0) || (i >= (int)this->pattern_prob_matrix_.size())) {
        throw std::out_of_range("Allele pattern is out of range of matrix");
    }
    return this->pattern_prob_matrix_.at(i);
}

void BiallelicPatternProbabilityMatrix::set_pattern_probability(
        unsigned int allele_count,
        unsigned int red_allele_count,
        double probability) {
    ECOEVOLITY_ASSERT(red_allele_count <= allele_count);
    if (allele_count == 0) {
        this->prob_missing_ = probability;
        return;
    }
    int i = (((allele_count * (allele_count + 1))/2) - 1 + red_allele_count);
    if ((i < 0) || (i >= (int)this->pattern_prob_matrix_.size())) {
        throw std::out_of_range("Allele pattern is out of range of matrix");
    }
    this->pattern_prob_matrix_.at(i) = probability;
}

unsigned int BiallelicPatternProbabilityMatrix::get_allele_count() const {
    return this->allele_count_;
}

const std::vector<double>& BiallelicPatternProbabilityMatrix::get_pattern_prob_matrix() const {
    return this->pattern_prob_matrix_;
}

std::string BiallelicPatternProbabilityMatrix::to_string() const {
    std::ostringstream ss;
    for (unsigned int i = 0; i <= this->get_allele_count(); ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            if (j > 0) {
                ss << "\t";
            }
            ss << this->get_pattern_probability(i, j);
        }
        ss << std::endl;
    }
    return ss.str();
}
