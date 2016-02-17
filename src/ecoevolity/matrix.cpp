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
    this->reset(allele_count);
    if (allele_count > 0) {
        this->set_pattern_probability(allele_count, red_allele_count, 1.0);
    }
}

BiallelicPatternProbabilityMatrix::BiallelicPatternProbabilityMatrix(
                unsigned int allele_count,
                const std::vector<double>& pattern_probs) {
    if (pattern_probs.size() != ((((allele_count + 1) * (allele_count + 2))/2) - 1)) {
        throw EcoevolityError(
                "Allele count and pattern probability vector size mismatch in BiallelicPatternProbabilityMatrix constructor");
    }
    this->allele_count_ = allele_count;
    this->pattern_prob_matrix_ = pattern_probs;
}

void BiallelicPatternProbabilityMatrix::resize(unsigned int allele_count) {
    if (this->get_allele_count() == allele_count) {
        ECOEVOLITY_ASSERT(this->pattern_prob_matrix_.size() ==
                ((((allele_count + 1) * (allele_count + 2))/2) - 1));
        return;
    }
    this->allele_count_ = allele_count;
    this->pattern_prob_matrix_.resize(
            ((((allele_count + 1) * (allele_count + 2))/2) - 1),
            0);
}

void BiallelicPatternProbabilityMatrix::reset(unsigned int allele_count) {
    if (this->get_allele_count() == allele_count) {
        ECOEVOLITY_ASSERT(this->pattern_prob_matrix_.size() ==
                ((((allele_count + 1) * (allele_count + 2))/2) - 1));
        std::fill(this->pattern_prob_matrix_.begin(), this->pattern_prob_matrix_.end(), 0);
        return;
    }
    this->allele_count_ = allele_count;
    this->pattern_prob_matrix_.assign(
            ((((allele_count + 1) * (allele_count + 2))/2) - 1),
            0);
}

void BiallelicPatternProbabilityMatrix::copy(
        const BiallelicPatternProbabilityMatrix& m) {
    this->resize(m.get_allele_count());
    ECOEVOLITY_ASSERT(this->get_allele_count() == m.get_allele_count());
    this->pattern_prob_matrix_ = m.get_pattern_prob_matrix();
}

const double& BiallelicPatternProbabilityMatrix::get_pattern_probability(
        unsigned int allele_count,
        unsigned int red_allele_count) const {
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
    int i = (((allele_count * (allele_count + 1))/2) - 1 + red_allele_count);
    if ((i < 0) || (i >= (int)this->pattern_prob_matrix_.size())) {
        throw std::out_of_range("Allele pattern is out of range of matrix");
    }
    this->pattern_prob_matrix_.at(i) = probability;
}

const unsigned int& BiallelicPatternProbabilityMatrix::get_allele_count() const {
    return this->allele_count_;
}

const std::vector<double>& BiallelicPatternProbabilityMatrix::get_pattern_prob_matrix() const {
    return this->pattern_prob_matrix_;
}

std::string BiallelicPatternProbabilityMatrix::to_string() const {
    std::ostringstream ss;
    for (unsigned int i = 1; i <= this->get_allele_count(); ++i) {
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
