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

#ifndef ECOEVOLITY_MATRIX_HPP
#define ECOEVOLITY_MATRIX_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "debug.hpp"
#include "assert.hpp"
#include "error.hpp"
#include "vector2d.hpp"
#include "complex.hpp"

class BiallelicPatternProbabilityMatrix {
    public:
        // Constructor
        BiallelicPatternProbabilityMatrix() { }
        BiallelicPatternProbabilityMatrix(unsigned int allele_count);
        BiallelicPatternProbabilityMatrix(
                const BiallelicPatternProbabilityMatrix& matrix);
        // Construct a leaf likelihood
        BiallelicPatternProbabilityMatrix(
                unsigned int allele_count,
                unsigned int red_allele_count);
        // Construct a top-of-branch likelihood
        BiallelicPatternProbabilityMatrix(
                unsigned int allele_count,
                const std::vector<double>& pattern_probs);
        // Destructor
        // ~AlleleProbabilityMatrix();
        //
        BiallelicPatternProbabilityMatrix& operator=(const BiallelicPatternProbabilityMatrix& m) {
            this->allele_count_ = m.allele_count_;
            this->pattern_prob_matrix_ = m.pattern_prob_matrix_;
            return * this;
        }

        BiallelicPatternProbabilityMatrix* clone() const {
            return new BiallelicPatternProbabilityMatrix(* this);
        }
        
        //Methods
        const double& get_pattern_probability(
                unsigned int allele_count,
                unsigned int red_allele_count) const;
        void set_pattern_probability(
                unsigned int allele_count,
                unsigned int red_allele_count,
                double probability);
        const unsigned int& get_allele_count() const;
        const std::vector<double>& get_pattern_prob_matrix() const;

        void resize(unsigned int allele_count);
        void reset(unsigned int allele_count);
        void copy(const BiallelicPatternProbabilityMatrix& m);

        std::string to_string() const;

    private:
        unsigned int allele_count_ = 0;
        std::vector<double> pattern_prob_matrix_;
};


/**
 * @brief   Abstract matrix (or linear operator) type; used for sparse
 *          matrices.
 *
 * @note    Translated/modified from AbstractMatrix class of the SnAP
 *          package.
 */
class AbstractMatrix {
    public:
        AbstractMatrix() { }
        virtual ~AbstractMatrix() { }

        //virtual Cloneable * clone() const = 0;
        virtual unsigned int get_nrows() const {return 0;}
        virtual unsigned int get_ncols() const {return 0;}
        virtual double inf_norm() const {return 0.0;}
        virtual double trace() const {return 0.0;}

        virtual void multiply(const std::vector<double>& v,
                              std::vector<double>& av) = 0;

        virtual void solve(const std::vector<Complex>& vc,
                           const Complex& offset,
                           std::vector<Complex>& xc) = 0;

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

            for (unsigned int i = 1; i <= n; ++i) {
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

 /**
 * @note    Translated/modified from QMatrix class of the SnAP
 *          package.
 */
class QMatrix : public AbstractMatrix {
    private:
        double u_;
        double v_;
        double coalescence_rate_;
        unsigned int n_;

    public:
        QMatrix(unsigned int n, double u, double v, double coalescence_rate) {
            this->n_ = n;
            this->u_ = u;
            this->v_ = v;
            this->coalescence_rate_ = coalescence_rate;
        }

        unsigned int get_nrows() const {
            return ((this->n_ + 1) * (this->n_ + 1) + this->n_ - 1)/2;
        }
        unsigned int get_ncols() const {
            return this->get_nrows();
        }

        double inf_norm() const {
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

        double trace() const {
            return -this->coalescence_rate_ * (this->n_ - 1) * this->n_ * (this->n_ + 1) * (this->n_ + 2)/8 - this->n_ * (this->n_ + 1) * (this->n_ + 2) * (this->u_ + this->v_)/6;
        }

        void multiply(const std::vector<double>& x,
                      std::vector<double>& ax) {
            unsigned int index = 1;
            for (unsigned int n = 1; n <= this->n_; ++n) {
                double sum = 0.0;
                sum += (- (this->coalescence_rate_*(n*(n-1.0)))/2.0 - this->v_*n)*x.at(index);
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
                           std::vector<double>& ax) {
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

        void check_mrr(double mrr) const {
            // TODO: Comparing floats for equality (!!)
            if (mrr == 0.0) {
                throw EcoevolityError("QMatrix: Error in matrix solve");
            }
        }

        // TODO: check if passed or member variable gets precedence in java;
        // assuming passed gets precedence in translating this from SnAP.
        std::vector<double> solve_central_block_transposed(
                const std::vector<double>& y,
                const double offset,
                const unsigned int n,
                const double u,
                const double v,
                const double coalescence_rate) {

            ECOEVOLITY_DEBUG(
            std::cerr << "ylocal = [";
            for (unsigned int i = 0; i < y.size(); ++i) {
                std::cerr << y.at(i) << " ";
            }
            std::cerr << "];" << std::endl;
            )

            std::vector<double> x (n + 1, 0.0);

            double K = -(coalescence_rate*(n*(n-1.0)))/2.0 - n*v + offset;
            
            // TODO: More comparing floats
            if ((u == 0.0) && (v == 0.0)) { 
                for (unsigned int r = 0; r <= n; ++r) {
                    x.at(r) = y.at(r)/K;
                }
            } else if (u == 0.0) {
                double Mrr = K;
                this->check_mrr(Mrr);
                x.at(0) = y.at(0) / Mrr;
                for(unsigned int r = 1; r <= n; ++r) {
                    Mrr = K+r*(v-u);
                    this->check_mrr(Mrr);
                    x.at(r) = (y.at(r) - ((n-r+1.0)*v)*x.at(r-1))/Mrr;    
                }
            } 
            else if (v == 0.0) {
                double Mrr = K + n*(v-u);
                this->check_mrr(Mrr);
                x.at(n) = y.at(n) / Mrr;
                for (int r = n-1; r >= 0; --r) {
                    Mrr = (K+r*(v-u));
                    x.at(r) = (y.at(r) - ((r+1.0)*u)*x.at(r+1))/Mrr;
                }
            }
            else {
                std::vector<double> d (n+1, 0.0);
                std::vector<double> e (n+1, 0.0);
                d.at(0) = K;
                e.at(0) = y.at(0);
                for(unsigned int r = 1; r <= n; ++r) {
                    this->check_mrr(d.at(r-1));
                    double m = ((n-r+1.0)*v)/d.at(r-1);
                    d.at(r) = K+r*(v-u) - m*r*u;
                    e.at(r) = y.at(r) - m*e.at(r-1);
                }
                
                ECOEVOLITY_DEBUG(
                    std::cerr << "d = [";
                    for (unsigned int i = 0; i < d.size(); ++i) {
                        std::cerr << d.at(i) << " ";
                    }
                    std::cerr << "];" << std::endl;
                    
                    std::cerr << "e = [";
                    for (unsigned int i = 0; i < e.size(); ++i) {
                        std::cerr << e.at(i) << " ";
                    }
                    std::cerr << "];" << std::endl;
                )
                
                //now solve the upper biadiagonal. diagonal is d, upper is same as M
                x.at(n) = e.at(n)/d.at(n);
                for (int r = n-1; r >= 0 ; --r) {
                    this->check_mrr(d.at(r));
                    x.at(r) = (e.at(r) - (r+1.0)*u*x.at(r+1))/d.at(r);
                }
            }
            
            ECOEVOLITY_DEBUG(
                std::cerr << "xlocal = ";
                for(unsigned int i = 0; i < x.size(); ++i) {
                    std::cerr << x.at(i) << " ";
                }
                std::cerr << "];" << std::endl;
                
                double diff = 0.0;
                for(unsigned int r=0; r <= n; ++r) {
                    double sum = 0.0;
                    if (r > 0) {
                        sum += (n+1.0-r)*v*x.at(r-1);
                    }
                    sum += (- (coalescence_rate*(n*(n-1.0)))/2.0 - v*(n-r) - u*r)*x.at(r);
                    if (r < n) {
                        sum+= (r+1.0)*u*x.at(r+1);
                    }
                    
                    std::cerr << sum << std::endl;
                    
                    diff = std::max(diff, std::abs(sum-y.at(r))); //Check that it solves the system
                }
                if (diff > 1e-10) {
                    std::cerr << "QMatrix.solve_central_block_transposed(): Error in solve" << std::endl;
                }
            )
            return x;
        }

        std::vector<double> find_orthogonal_vector() {
            //Finds non-zero vector x such that x'Q' = 0
            /* Want to find non-zero x such that x' Q' = 0.
             Let x = [x_1,x_2,...,x_n]
             We first solve x_1 ' M_1 = 0 manually.
             Then use the recursion
             
             x_{n-1} R_{n-1} + x_n M_n = 0  for n=2,3,...,N
             
             */
            unsigned int nrows = this->get_nrows();
        
            std::vector<double> x (nrows + 1, 0.0);        
            std::vector<double> xn (this->n_ + 1, 0.0); 
            std::vector<double> yn (this->n_ + 1, 0.0);
            
            //First solve the n=1 case.
            //M_1 = [-v   v; -u u] which has solution [u  v]
            xn.at(0) = this->u_; xn.at(1) = this->v_;
            x.at(1) = this->u_;
            x.at(2) = this->v_;
            unsigned int xptr = 3;
            
            //Now the rest.
            for (unsigned int n = 2; n <= this->n_; ++n) {
                
                /***
                 We use the recursion 
                 x_{n-1}' R_{n-1} + x_n' M_n = 0.
                 or equivalently
                 M_n' x_n = - R_{n-1}' x_{n-1}
                 
                 
                 First compute y_n = - R_{n-1}' x_{n-1}  
                 
                 Second, solve M_n' x_n = y_n
                 
                 ****/
                
                ECOEVOLITY_DEBUG(
                    std::cerr << "xn = [";
                    for (unsigned int r = 0; r <= n-1; ++r) {
                        std::cerr << xn.at(r) << " ";
                    }
                    std::cerr << "];" << std::endl;
                )
                
                yn.at(0) = - ((this->coalescence_rate_*(n-1.0)*n)/2.0)*xn.at(0);
                for (unsigned int r = 1; r < n; ++r) {
                    yn.at(r) = - ( (this->coalescence_rate_*(r-1.0)*n)/2.0 )*xn.at(r-1) - ( (this->coalescence_rate_*(n-1.0-r)*n)/2.0 )*xn.at(r);
                }        
                yn.at(n) = - ( (this->coalescence_rate_*(n-1.0)*n)/2.0 )*xn.at(n-1);
                
                ECOEVOLITY_DEBUG(
                    std::cerr << "yn = [";
                    for(unsigned int r = 0;r <= n; ++r) {
                        std::cerr << yn.at(r) << " ";
                    }
                    std::cerr << "];" << std::endl;
                )
                
                
                xn = this->solve_central_block_transposed(
                        yn,
                        0,
                        n,
                        this->u_,
                        this->v_,
                        this->coalescence_rate_);
                
                for (unsigned int i = 0; i < xn.size(); ++i) {
                    x.at(xptr++) = xn.at(i);
                }
                
                ECOEVOLITY_DEBUG(
                    std::cerr << "xn2 = [";
                    for(unsigned int r = 0; r <= n; ++r) {
                        std::cerr << xn.at(r) << " ";
                    }
                    std::cerr << "];" << std::endl;
                )
            }
            return x;
        }

        void solve(const std::vector<Complex>& y,
                   const Complex& offset,
                   std::vector<Complex>& x) {
            /* Suppose that y = [y1',y2',...,yn']' as above. We solve (Q^t + offset*I) x = y.
             This gives us the equations
             
             M_n x_n + R_n x_{n+1} + offset*x_n  = y_n  for n=1,2,...,N - 1
             and
             M_N x_N + offset x_N = y_N
             
             These can be solved in reverse, starting with (M_N + offset * I)x_N = y_N
             and then
             (M_n + offset I) x_n = y_n - R_n x_{n+1}  n = N-1,N_2,...,1
             
             The solution of the tridiaongal system (M_n + offset I) x_n = z 
             is done by the routine solveCentralBlock. The computation of Rn x_{n+1} is done by MultiplyUpperBlock.
             */
            
            if (x.size() != y.size()) {
                throw EcoevolityError("QMatrix.solve(): x & y differ in length");
            }
            for (unsigned int i = 0; i < x.size(); ++i) {
                x.at(i).re_ = 0;
                x.at(i).im_ = 0;
            }
            
            std::vector<Complex> xn (this->n_ + 1);
            std::vector<Complex> yn (this->n_ + 1);
            // for (unsigned int i = 0; i < (this->n_ + 1); ++i) {
            //     xn.at(i) = Complex();
            //     yn.at(i) = Complex();
            // }
            
            
            unsigned int xptr = x.size() - 1 - this->n_;
            unsigned int yptr = y.size() - 1 - this->n_;
            
            //Solve (M_N + offset * I)x_N = y_N and copy solution into xn.
            for (unsigned int i = 0; i < yn.size(); ++i) {
                yn.at(i).re_ = y.at(yptr + i).re_;
                yn.at(i).im_ = y.at(yptr + i).im_;
            }
            this->solve_central_block(
                    yn,
                    offset,
                    this->n_,
                    this->u_,
                    this->v_,
                    this->coalescence_rate_,
                    xn);

            for (unsigned int i = 0; i < xn.size(); ++i) {
                x.at(xptr+i).re_ = xn.at(i).re_;
                x.at(xptr+i).im_ = xn.at(i).im_;
            }
            
            //Solve for the rest
            for (unsigned int n = this->n_-1; n >= 1; --n) {
                xptr = xptr - (n+1);
                this->multiply_upper_block(
                        xn,
                        n,
                        this->coalescence_rate_,
                        yn);
                yptr = yptr - (n+1);
                for (unsigned int r = 0; r <= n; ++r) {
                    yn.at(r).re_ = y.at(yptr + r).re_ - yn.at(r).re_; 
                    yn.at(r).im_ = y.at(yptr + r).im_ - yn.at(r).im_; 
                }
                
                this->solve_central_block(
                        yn,
                        offset,
                        n,
                        this->u_,
                        this->v_,
                        this->coalescence_rate_,
                        xn);
                
                for (unsigned int i = 0; i <= n; ++i) {
                    x.at(xptr + i).re_ = xn.at(i).re_;
                    x.at(xptr + i).im_ = xn.at(i).im_;
                }
            }
        }

        void multiply_upper_block(
                const std::vector<Complex>& x,
                unsigned int n,
                double coalescence_rate,
                std::vector<Complex>& y) {
            if (y.size() < n+1) {
                throw EcoevolityError("QMatrix.multiply_upper_block(): y vector of Complexes shorter than expected");
            }
            for (unsigned int r = 0;r <= n; ++r) { 
              y.at(r).muladd((coalescence_rate*(double)(n-r)*(n+1.0)/2.0), x.at(r), (coalescence_rate*(double)r*(n+1.0)/2.0),x.at(r+1));
            }
        }

        void check_mrr(const Complex& mrr) const {
            if ((std::abs(mrr.re_) < 1e-20) && (std::abs(mrr.im_) < 1e-20)) {
                throw EcoevolityError("QMatrix.check_mrr(): Error in matrix solve");
            }
            
        }

        void solve_central_block(
                const std::vector<Complex>& y,
                const Complex& offset,
                unsigned int n,
                double u,
                double v,
                double coalescence_rate,
                std::vector<Complex>& x) {
            Complex K = Complex((double)(-coalescence_rate*(double)n*(n-1.0)/2.0 - n*v) + offset.re_, offset.im_);
            Complex tmp = Complex();

            if (x.size() < n+1) {
                throw EcoevolityError("QMatrix.solve_central_block(): x vector of Complexes shorter than expected");
            }
            
            Complex Mrr = Complex();
            
            if (v == 0.0) {
                //Lower bidagonal.
                Mrr = K;
                this->check_mrr(Mrr);
                x.at(0).divide(y.at(0), Mrr);
                for (unsigned int r = 1; r <= n; ++r) {
                    Mrr.re_ = K.re_ + (double)r*(v-u);
                    Mrr.im_ = K.im_;
                    this->check_mrr(Mrr);
                    tmp.re_ = y.at(r).re_ - x.at(r-1).re_*(double)(u*r);
                    tmp.im_ = y.at(r).im_ - x.at(r-1).im_*(double)(u*r);
                    x.at(r).divide(tmp, Mrr);
                }
            } //Upper bidiagonal
            else if (u == 0.0) {
                Mrr.re_ = K.re_ + (double)n*(v-u);
                Mrr.im_ = K.im_;
                this->check_mrr(Mrr);
                x.at(n).divide(y.at(n), Mrr);
                for (int r = n-1; r >= 0; --r) {
                    Mrr.re_ = K.re_ + (double)r*(v-u);
                    Mrr.im_ = K.im_;
                    tmp.re_ = y.at(r).re_ - ((double)(n-r)*v*x.at(r+1).re_);
                    tmp.im_ = y.at(r).im_ - ((double)(n-r)*v*x.at(r+1).im_);
                    x.at(r).divide(tmp, Mrr);
                }
            }
            else {
                std::vector<Complex> d (n+1, Complex());
                // for (unsigned int i = 0; i < n+1; ++i) {
                //     d.at(i) = Complex();
                // }
                std::vector<Complex> e (n+1, Complex());
                // for (unsigned int i = 0; i < n+1; ++i) {
                //     e.at(i) = Complex();
                // }
                d.at(0).re_ = K.re_;
                d.at(0).im_ = K.im_;
                e.at(0).re_ = y.at(0).re_;
                e.at(0).im_ = y.at(0).im_;
                for (unsigned int r = 1;r <= n; ++r) {
                    //zero out lower triangular
                    this->check_mrr(d.at(r-1));
                    tmp.re_ = r*u;
                    tmp.im_ = 0.0;
                    Complex m = Complex();
                    m.divide(tmp, d.at(r-1));
                    d.at(r).re_ = K.re_+r*(v-u) - m.re_*(n-r+1.0)*v;
                    d.at(r).im_ = K.im_+        - m.im_*(n-r+1.0)*v;
                    tmp.mul(m, e.at(r-1));
                    e.at(r).re_ = y.at(r).re_ - tmp.re_;
                    e.at(r).im_ = y.at(r).im_ - tmp.im_;
                }
                
                //now solve the upper biadiagonal. diagonal is d, upper is same as M
                x.at(n).divide(e.at(n), d.at(n));
                for (int r = n-1; r >= 0; --r) {
                    this->check_mrr(d.at(r));
                    tmp.re_ = (e.at(r).re_ - (double)(n-r)*v*x.at(r+1).re_);
                    tmp.im_ = (e.at(r).im_ - (double)(n-r)*v*x.at(r+1).im_);
                    x.at(r).divide(tmp, d.at(r));
                }
            }
        }

        void check_mrr(double mrr_r, double mrr_i) const {
            if ((std::abs(mrr_r) < 1e-20) && (std::abs(mrr_i) < 1e-20)) {
               throw EcoevolityError("QMatrix.check_mrr(): Error in matrix solve");
            }
        }

        void multiply_upper_block(
                const std::vector<double>& x_r,
                const std::vector<double>& x_i,
                unsigned int n,
                double coalescence_rate,
                std::vector<double>& y_r,
                std::vector<double>& y_i) {
            if ((y_r.size() < n+1) || (y_i.size() < n+1)) {
                throw EcoevolityError("QMatrix.multiply_upper_block(): y vector shorter than expected");
            }
            for (unsigned int r = 0; r <= n; ++r) { 
                y_r.at(r) = coalescence_rate*(double)(n-r)*(n+1.0)/2.0 * x_r.at(r) + (coalescence_rate*(double)r*(n+1.0)/2.0)*x_r.at(r+1);
                y_i.at(r) = coalescence_rate*(double)(n-r)*(n+1.0)/2.0 * x_i.at(r) + (coalescence_rate*(double)r*(n+1.0)/2.0)*x_i.at(r+1);
            }
            
        }

        void solve_central_block(
                const std::vector<double>& y_r,
                const std::vector<double>& y_i,
                double offset_r,
                double offset_i,
                unsigned int n,
                double u,
                double v,
                double coalescence_rate,
                std::vector<double>& x_r,
                std::vector<double>& x_i) {
            double K_r = (double)(-coalescence_rate*(double)n*(n-1.0)/2.0 - n*v) + offset_r;
            double K_i = offset_i;
            double tmp_r;
            double tmp_i;
            double Mrr_r;
            double Mrr_i;

            if (v == 0.0) {
                //Lower bidagonal.
                Mrr_r = K_r;
                Mrr_i = K_i;
                this->check_mrr(Mrr_r, Mrr_i);
                double f = Mrr_r * Mrr_r + Mrr_i * Mrr_i; 
                x_r.at(0) = (y_r.at(0) * Mrr_r + y_i.at(0) * Mrr_i) / f;
                x_i.at(0) = (y_i.at(0) * Mrr_r - y_r.at(0) * Mrr_i) / f;

                for (unsigned int r = 1; r <= n; ++r) {
                    Mrr_r = K_r + (double)r*(v-u);
                    Mrr_i = K_i;
                    this->check_mrr(Mrr_r, Mrr_i);
                    tmp_r = y_r.at(r) - x_r.at(r-1)*(double)(u*r);
                    tmp_i = y_i.at(r) - x_i.at(r-1)*(double)(u*r);
                    f = Mrr_r * Mrr_r + Mrr_i * Mrr_i; 
                    x_r.at(r) = (tmp_r * Mrr_r + tmp_i * Mrr_i) / f;
                    x_i.at(r) = (tmp_i * Mrr_r - tmp_r * Mrr_i) / f;
                }
            } //Upper bidiagonal
            else if (u == 0.0) {
                Mrr_r = K_r + (double)n*(v-u);
                Mrr_i = K_i;
                this->check_mrr(Mrr_r, Mrr_i);
                double f = Mrr_r * Mrr_r + Mrr_i * Mrr_i; 
                x_r.at(n) = (y_r.at(n) * Mrr_r + y_i.at(n) * Mrr_i) / f;
                x_i.at(n) = (y_i.at(n) * Mrr_r - y_r.at(n) * Mrr_i) / f;
                for (int r = n-1; r >= 0; --r) {
                    Mrr_r = K_r + (double)r*(v-u);
                    Mrr_i = K_i;
                    tmp_r = y_r.at(r) - ((double)(n-r)*v*x_r.at(r+1));
                    tmp_i = y_i.at(r) - ((double)(n-r)*v*x_i.at(r+1));
                    f = Mrr_r * Mrr_r + Mrr_i * Mrr_i; 
                    x_r.at(r) = (tmp_r * Mrr_r + tmp_i * Mrr_i) / f;
                    x_i.at(r) = (tmp_i * Mrr_r - tmp_r * Mrr_i) / f;
                }
            }
            else {
                std::vector<double> d_r (n+1, 0.0);
                std::vector<double> d_i (n+1, 0.0);
                std::vector<double> e_r (n+1, 0.0);
                std::vector<double> e_i (n+1, 0.0);
                d_r.at(0) = K_r;
                d_i.at(0) = K_i;
                e_r.at(0) = y_r.at(0);
                e_i.at(0) = y_i.at(0);
                for (unsigned int r = 1; r <= n; ++r) {
                    this->check_mrr(d_r.at(r-1), d_i.at(r-1));
                    tmp_r = r*u;
                    tmp_i = 0.0;
                    double f = d_r.at(r-1) * d_r.at(r-1) + d_i.at(r-1) * d_i.at(r-1); 
                    double m_r = (tmp_r * d_r.at(r-1) + tmp_i * d_i.at(r-1)) / f;
                    double m_i = (tmp_i * d_r.at(r-1) - tmp_r * d_i.at(r-1)) / f;
                    d_r.at(r) = K_r + r*(v-u) - m_r*(n-r+1.0)*v;
                    d_i.at(r) = K_i +         - m_i*(n-r+1.0)*v;
                    tmp_r = m_r * e_r.at(r-1) - m_i * e_i.at(r-1);
                    tmp_i = m_i * e_r.at(r-1) + m_r * e_i.at(r-1);
                    e_r.at(r) = y_r.at(r) - tmp_r;
                    e_i.at(r) = y_i.at(r) - tmp_i;
                }

                //now solve the upper biadiagonal. diagonal is d, upper is same as M
                double f = d_r.at(n) * d_r.at(n) + d_i.at(n) * d_i.at(n); 
                x_r.at(n) = (e_r.at(n) * d_r.at(n) + e_i.at(n) * d_i.at(n)) / f;
                x_i.at(n) = (e_i.at(n) * d_r.at(n) - e_r.at(n) * d_i.at(n)) / f;
                for (int r = n-1; r >= 0; --r) {
                    this->check_mrr(d_r.at(r), d_i.at(r));
                    tmp_r = (e_r.at(r) - (double)(n-r)*v*x_r.at(r+1));
                    tmp_i = (e_i.at(r) - (double)(n-r)*v*x_i.at(r+1));
                    f = d_r.at(r) * d_r.at(r) + d_i.at(r) * d_i.at(r); 
                    x_r.at(r) = (tmp_r * d_r.at(r) + tmp_i * d_i.at(r)) / f;
                    x_i.at(r) = (tmp_i * d_r.at(r) - tmp_r * d_i.at(r)) / f;
                }
            }
        }

        void solve(
                const std::vector<double>& y_r,
                const std::vector<double>& y_i,
                double offset_r,
                double offset_i,
                const std::vector<double>& x_r,
                const std::vector<double>& x_i) {
            /* Suppose that y = [y1',y2',...,yn']' as above. We solve (Q^t + offset*I) x = y.
             This gives us the equations
             
             M_n x_n + R_n x_{n+1} + offset*x_n  = y_n  for n=1,2,...,N - 1
             and
             M_N x_N + offset x_N = y_N
             
             These can be solved in reverse, starting with (M_N + offset * I)x_N = y_N
             and then
             (M_n + offset I) x_n = y_n - R_n x_{n+1}  n = N-1,N_2,...,1
             
             The solution of the tridiaongal system (M_n + offset I) x_n = z 
             is done by the routine solveCentralBlock. The computation of Rn x_{n+1} is done by MultiplyUpperBlock.
             */
            
            x_r.assign(x_r.size(), 0.0);
            x_i.assign(x_i.size(), 0.0);
            
            std::vector<double> xn_r (this->n_ + 1, 0.0);
            std::vector<double> xn_i (this->n_ + 1, 0.0);
            std::vector<double> yn_r (this->n_ + 1, 0.0);
            std::vector<double> yn_i (this->n_ + 1, 0.0);
            
            int xptr = x_r.size() - 1 - this->n_;
            int yptr = y_r.size() - 1 - this->n_;
            
            for (unsigned int i = 0; i < yn_r.size(); ++i) {
                yn_r.at(i) = y_r.at(yptr + i);
                yn_i.at(i) = y_i.at(yptr + i);
            }
            this->solve_central_block(
                    yn_r, 
                    yn_i,
                    offset_r,
                    offset_i,
                    this->n_,
                    this->u_,
                    this->v_,
                    this->coalescence_rate_,
                    xn_r,
                    xn_i);
            
            
            for (unsigned int i = 0; i < xn_r.size(); ++i) {
                x_r.at(xptr+i) = xn_r.at(i);
                x_i.at(xptr+i) = xn_i.at(i);
            }
            
            //Solve for the rest
            for (unsigned int n = this->n_-1; n >= 1; --n) {
                xptr = xptr - (n+1);
                this->multiply_upper_block(
                        xn_r,
                        xn_i,
                        n,
                        this->coalescence_rate_,
                        yn_r,
                        yn_i);
                yptr = yptr - (n+1);
                for (unsigned int r = 0; r <= n; ++r) {
                    yn_r.at(r) = y_r.at(yptr + r) - yn_r.at(r); 
                    yn_i.at(r) = y_i.at(yptr + r) - yn_i.at(r); 
                }
                
                this->solve_central_block(
                        yn_r,
                        yn_i,
                        offset_r,
                        offset_i,
                        n,
                        this->u_,
                        this->v_,
                        this->coalescence_rate_,
                        xn_r,
                        xn_i);
                
                for (unsigned int i = 0; i <= n; ++i) {
                    x_r.at(xptr + i) = xn_r.at(i);
                    x_i.at(xptr + i) = xn_i.at(i);
                }
            }
        }
};

#endif
