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

#ifndef ECOEVOLITY_QMATRIX_HPP
#define ECOEVOLITY_QMATRIX_HPP

/**
 * @brief Abstract matrix (or linear operator) type; used for sparse matrices.
 */
class AbstractMatrix {
    public:
        AbstractMatrix() { }
        virtual ~AbstractMatrix() { }

        //virtual Cloneable * clone() const = 0;
        virtual unsigned int get_nrows() {return 0;}
        virtual unsigned int get_ncols() {return 0;}
        virtual double inf_norm() {return 0.0;}
        virtual double trace() {return 0.0};

        virtual void multiply(const std::vector<double>& v,
                              const std::vector<double>& av) = 0;

        virtual void solve(const std::vector<Complex>& vc,
                           const Complex& offset,
                           const std::vector<Complex>& xc) = 0;

        virtual void solve(const std::vector<double>& vc_r,
                           const std::vector<double>& vc_i,
                           const double offset_r,
                           const double offset_i,
                           const std::vector<double>& xc_r,
                           const std::vector<double>& xc_i) = 0;

        /**
         * @brief   Computes infinity norm (pg 56 of Golub and van Loan)
         */
        double infinity_norm(Vector2d v) {
            double norm = 0.0;
            unsigned int n = v.get_nrows();
            unsigned int m = v.get_ncols();

            for (unsigned in i = 1; i <= n; ++i) {
                double sum = 0.0;
                for(unsigned int j = 1; j <= m; ++j) {
                    sum += std::abs(v.at(i,j));
                }
                norm = std::max(norm, sum);
            }
            return norm;
        }

        /**
         * @brief   Zero-based vector norm.
         */
        double vector_norm(const std::vector<double>& v) {
            unsigned int n = v.size();
            double sum = 0.0;
            for (unsigned int i = 0; i < n; ++i) {
                sum += v.at(i) * v.at(i);
            }
            return std::sqrt(sum);
        }

        /**
         * @brief   One-based vector norm.
         */
        double vector_norm_1_based(const std::vector<double>& v) {
            unsigned int n = v.size() - 1;
            double sum = 0.0;
            for (unsigned int i = 1; i <= n; ++i) {
                sum += v.at(i) * v.at(i);
            }
            return std::sqrt(sum);
        }
};

#endif
