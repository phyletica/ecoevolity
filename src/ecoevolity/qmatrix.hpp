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

#include "vector2d.hpp"
#include "complex.hpp"

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

class QMatrix : public AbstractMatrix {
    private:
        double u_;
        double v_;
        double coalescence_rate_;
        unsigned int n_;

    public:
        Qmatrix(unsigned int n, double u, double v, double coalescence_rate) {
            this->n_ = n;
            this->u_ = u;
            this->v_ = v;
            this->coalescence_rate_ = coalescence_rate;
        }

        unsigned int get_nrows() const {
            return ((this-n_ + 1) * (this->n_ + 1) + this->n_ - 1)/2;
        }
        unsigned int get_ncols() const {
            return this->get_nrows();
        }

        double inf_norm() {
            if (this->u_ > this->v_) {
                return std::max(
                        2.0 * this->u_ * (this->n_ - 1) + this->coalescence_rate_ * (this->n_ - 1) * (this->n_ - 1),
                        2.0 * this->u_ * (this->n_) + this->coalescence_rate_ * (this->n_) * (this->n_ - 1)/2.0
                        );
            }
            else {
                return std::max(
                        2.0 * this->v_ * (this->n_ - 1) + this->coalescence_rate_ * (this->n_ - 1) * (this->n_ - 1),
                        2.0 * this->v_ * (this->n_) + this->coalescence_rate_ * (this->n_) * (this->n_ - 1)/2.0
                        );
            }
        }

        double trace() {
            return -this->coalescence_rate_ * (this->n_ - 1) * this->n_ * (this->n_ + 1) * (this->n_ + 2)/8 - this->n_ * (this->n_ + 1) * (this->n_ + 2) * (this->u_ + this->v_)/6;
        }

        void multiply(const std::vector<double>& x,
                      const std::vector<double>& ax) {
            unsigned int index = 1;
            for (unsigned int n = 1; n <= this->n_; ++n) {
                double sum = 0.0;
                sum += (- (this->coalescence_rate*(n*(n-1.0)))/2.0 - this->v_*n)*x.at(index);
                sum += n*this->v_*x.at(index+1);
                
                if (n < this->n_) {
                    sum += (n*(n+1.0)/2.0)*this->coalescence_rate_*x.at(index+n+1);
                }
                ax.at(index) = sum;
                index++;
                
                for(unsigned int r = 1; r < n; ++r) {
                    sum = r*this->u_*x.at(index-1);
                    sum +=  (- (this->coalescence_rate_*(n*(n-1.0)))/2.0 - this->v_*(n-r) - this->u_*r)*x.at(index);
                    if (r < n)
                        sum += (n-r)*this->v_*x.at(index+1);
                    
                    if (n < this->n_) {
                        sum += ((n-r)*(n+1.0)/2.0)*this->coalescence_rate_*x.at(index+n+1);
                        sum += (r*(n+1.0)/2.0)*this->coalescence_rate_*x.at(index+n+2);
                    }
                    ax.at(index) = sum;
                    index++;
                }
                sum = n*this->u_*x.at(index-1);
                sum += (- (this->coalescence_rate_*(n*(n-1.0)))/2.0 - this->u_*n)*x.at(index);
            
                if (n < this->n_) {
                    sum += (n*(n+1.0)/2.0)*this->coalescence_rate_*x.at(index+n+2);
                }
                ax.at(index) = sum;
                index++;
            }
        }

        void multiply_orig(const std::vector<double>& x,
                           const std::vector<double>& ax) {
            unsigned int index = 1;
            for(unsigned int n = 1; n <= this->n_; ++n) {
                for(unsigned int r = 0; r <= n; ++r) {
                    double sum = 0.0;
                    if (r > 0) {
                        sum += r*this->u_*x.at(index-1);
                    }
                    sum += (- (this->coalescence_rate_*(n*(n-1.0)))/2.0 - this->v_*(n-r) - this->u_*r)*x.at(index);
                    if (r < n)
                        sum += (n-r)*this->v_*x.at(index+1);
                    
                    if (n < this->n_) {
                        sum += ((n-r)*(n+1.0)/2.0)*this->coalescence_rate_*x.at(index+n+1);
                        sum += (r*(n+1.0)/2.0)*this->coalescence_rate_*x.at(index+n+2);
                    }
                    ax.at(index) = sum;
                    index++;
                }
            }
        }
};

#endif
