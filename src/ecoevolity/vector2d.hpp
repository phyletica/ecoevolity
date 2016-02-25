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

/**
 * @brief   1-based 2x2 matrix
 */
class Vector2d {
    private:
        unsigned int nrows_;
        unsigned int ncols_;
        std::vector<double> values_;

    public:
        Vector2d(unsigned int nrows, unsigned int ncols) {
            this->nrows_ = nrows;
            this->ncols_ = ncols;
            this->values_.resize((nrows * ncols), 0);
        }

        Vector2d(unsigned int nrows,
                unsigned int ncols,
                const std::vector<double>& values) {
            ECOEVOLITY_ASSERT((nrows * ncols) == values.size());
            this->nrows_ = nrows;
            this->ncols_ = ncols;
            this->values_ = values;
        }

        /**
         * @brief   Initialize square identitiy matrix
         */
        Vector2d(unsigned int nrows) {
            this->nrows_ = nrows;
            this->ncols_ = nrows;
            this->values_.resize((nrows * nrows), 0);
            for (unsigned int i = 1; i <= nrows; ++i) {
                this->set(i, i, 1.0);
            }
        }

        Vector2d& operator=(const Vector2d& a) {
            this->nrows_ = a.nrows_;
            this->ncols_ = a.ncols_;
            this->values_ = a.values_;
            return * this;
        }

        unsigned int get_nrows() const {
            return this->nrows_;
        }
        unsigned int get_ncols() const {
            return this->ncols_;
        }

        unsigned int get_index(unsigned int i, unsigned int j) const {
            return ((i - 1) * this->ncols_) + j - 1;
        }

        const double& at(unsigned int i, unsigned int j) const {
            return this->values_.at(this->get_index(i, j));
        }

        void set(unsigned int i, unsigned int j, double v) {
            this->values_.at(this->get_index(i,j)) = v;
        }
        void set(const std::vector<double>& zero_based_vector_of_values,
                 unsigned int nrows) {
            ECOEVOLITY_ASSERT((nrows*nrows) == zero_based_vector_of_values.size());
            this->nrows_ = nrows;
            this->ncols_ = nrows;
            this->values_ = zero_based_vector_of_values;
        }

        void divide(unsigned int i, unsigned int j, double v) {
            this->values_.at(this->get_index(i,j)) /= v;
        }

        void multiply(unsigned int i, unsigned int j, double v) {
            this->values_.at(this->get_index(i,j)) *= v;
        }

        void subtract(unsigned int i, unsigned int j, double v) {
            this->values_.at(this->get_index(i,j)) -= v;
        }

        void add(unsigned int i, unsigned int j, double v) {
            this->values_.at(this->get_index(i,j)) += v;
        }
        void add(double multiplier,
                 const std::vector<double>& addends) {
            ECOEVOLITY_ASSERT(this->values_.size() == addends.size());
            for (unsigned int = 0; i < this->values_size(); ++i) {
                this->values_.at(i) += multiplier * addends.at(i);
            }
        }

        /**
         * @brief   Get transposed copy of current matrix
         */
        Vector2d transpose() const {
            Vector2d transpose(this->nrows_, this->ncols_);
            for (unsigned int row_idx = 1; row_idx <= this->nrows_; ++row_idx) {
                for (unsigned int col_idx = 1; col_idx <= this->ncols_; ++col_idx) {
                    transpose.set(col_idx, row_idx, this->at(row_idx, col_idx));
                }
            }
            return transpose;
        }

        /**
         * @brief   Get matrix as zero-based vector
         *
         * Matrix item (i, j) will be at index ((i-1)*ncol)+j-1 in the vector.
         */
        const std::vector<double>& get_zero_based_vector() const {
            return this->values_;
        }

        /**
         * @brief   Get column as 1-based vector
         */
        std::vector<double> get_column(unsigned int column_index) const {
            std::vector<double> c (this->nrows_ + 1, 0);
            for (unsigned int row_idx = 1; row_idx <= this->nrows_; ++row_idx) {
                c.at(row_idx) = this->at(row_idx, column_index);
            }
            return c;
        }
        void get_column(unsigned int column_index,
                        std::vector<double>& column) const {
            column.resize(this->nrows_ + 1);
            for (unsigned int row_idx = 1; row_idx <= this->nrows_; ++row_idx) {
                column.at(row_idx) = this->at(row_idx, column_index);
            }
        }

        void scale(double multiplier) {
            for (unsigned int i = 0; i < this->values_.size(); ++i) {
                this->values_.at(i) *= multiplier;
            }
        }

        void resize(unsigned int nrows, unsigned int ncols) {
            if ((nrows == this->nrows_) && (ncols == this->ncols_)) {
                return;
            }
            std::vector<double> v ((nrows*ncols), 0);
            if (this->nrows_ > nrows) {
                this->nrows_ = nrows;
            }
            for (unsigned int row_idx = 0; row_idx < this->nrows_; ++row_idx) {
                std::copy(this->values_.begin() + (row_idx * this->ncols_),
                          this->values_.begin() + (row_idx * this->ncols_) + std::min(this->ncols_, ncols),
                          v.begin() + (row_idx * ncols));
            }
            this->values_ = v;
            this->nrows_ = nrows;
            this->ncols_ = ncols;
        }

        std::string to_string() const {
            std::ostringstream ss;
            for (unsigned int row_idx = 1; row_idx <= this->nrows_; ++row_idx) {
                for (unsigned int col_idx = 1; col_idx <= this->ncols_; ++col_idx) {
                    if (col_idx > 1) {
                        ss << "\t";
                    }
                    ss << this->values_.at(this->get_index(row_idx, col_idx));
                }
                ss << std::endl;
            }
            return ss.str();
        }
};
