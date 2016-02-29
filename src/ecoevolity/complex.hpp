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

#ifndef ECOEVOLITY_COMPLEX_HPP
#define ECOEVOLITY_COMPLEX_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "debug.hpp"
#include "assert.hpp"
#include "error.hpp"

/**
 * @note    Translated from COMPLEX class of COMPLEX.java of the SnAP package.
 */
class Complex {
    public:
        double re_ = 0.0;
        double im_ = 0.0;

        Complex() { }
        Complex(double re, im) {
            this->re_ = re;
            this->im_ = im;
        }

        void divide(const Complex& numerator, double divisor) {
            this->re_ = numerator.re_ / divisor;
            this->im_ = numerator.im_ / divisor;
        }

        void muladd(const Complex& m1, const Complex& m2) {
            this->re_ += m2.re_ * m1.re_ - m2.im_ * m1.im_;
            this->im_ += m2.re_ * m1.im_ + m2.im_ * m1.re_;
        }

        void muladd(double f1, const Complex& x1, double f2, const Complex& x2) {
            this->re_ = f1 * x1.re_ + f2 * x2.re_;
            this->im_ = f1 * x1.im_ + f2 * x2.im_;
        }

        void divide(const Complex& numerator, const Complex& divisor) {
            double f = divisor.re_ * divisor.re_ + divisor.im_ * divisor.im_;
            this->re_ = (numerator.re_ * divisor.re_ + numerator.im_ * divisor.im_) / f;
            this->im_ = (numerator.im_ * divisor.re_ - numerator.re_ * divisor.im_) / f;
        }

        void mul(const Complex& m1, const Complex& m2) {
            this->re_ = m2.re_ * m1.re_ - m2.im_ * m1.im_;
            this->im_ = m2.re_ * m1.im_ + m2.im_ * m1.re_;
        }

        std::string to_string() {
            std::ostringstream ss;
            ss << "(" << this->re_ << ", " << this->im_ << ")";
            return ss.str();
        }
};

#endif
