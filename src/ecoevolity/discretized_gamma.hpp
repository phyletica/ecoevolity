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

#ifndef DISCRETIZED_GAMMA_HPP
#define DISCRETIZED_GAMMA_HPP

#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>

#include "math_util.hpp"
#include "parameter.hpp"

class DiscretizedGamma {
    protected:
        PositiveRealParameter shape_(1.0, true);
        PositiveRealParameter scale_(1.0, true);
        unsigned int num_categories_ = 4;
        std::vector<double> quantiles_;
        bool use_bin_median_ = false;

        /*!
         * Update the quantiles of the bins.
         *
         * @note    Modified from revbayes (commit c8bf96ec786)
         *          RevBayesCore::DiscretizeGammaFunction::update
         *          (https://github.com/revbayes/revbayes)
         */
        void update() {
            double a = this->get_shape();
            double b = 1.0 / this->get_scale();
            int nCats = (int)this->get_num_categories();
            
            double factor = a / b * nCats;

            if (this->use_bin_median_) {
                /* the median value for each category is used to represent all of the values
                in that category */
                double interval = 1.0 / (2.0 * nCats);
                for (int i=0; i<nCats; i++) 
                    this->quantiles_.at(i) = chi_square_quantile((i * 2.0 + 1.0) * interval, 2.0 * a) / (2.0 * b);
                double t = 0.0;
                for (int i=0; i<nCats; i++) 
                    t += this->quantiles_.at(i);
                for (int i=0; i<nCats; i++)     
                    this->quantiles_.at(i) *= factor / t;
            }
            else
            {
                /* the mean value for each category is used to represent all of the values
                in that category */
                /* calculate the points in the gamma distribution */
                for (int i=0; i<nCats-1; i++) 
                    this->quantiles_.at(i) = chi_square_quantile((i + 1.0) / nCats, 2.0 * a) / (2.0 * b);
                /* calculate the cumulative values */
                double lnGammaValue = std::lgamma(a + 1.0);
                for (int i=0; i<nCats-1; i++) 
                    this->quantiles_.at(i) = incomplete_gamma((*value)[i] * b, a + 1.0, lnGammaValue);
                this->quantiles_.at(nCats-1) = 1.0;
                /* calculate the relative values and rescale */
                for (int i=nCats-1; i>0; i--){
                    this->quantiles_.at(i) -= this->quantiles_.at(i-1);
                    this->quantiles_.at(i) *= factor;
                }
                this->quantiles_.at(0) *= factor;
            }
        }

    public:
        double get_shape() const {
            return this->shape_.get_value();
        }
        double get_scale() const {
            return this->scale_.get_value();
        }
        double get_mean() const {
            return this->get_shape() * this->get_scale();
        }
        unsigned int get_num_categories() const {
            return this->num_categories_;
        }

        void set_shape(double shape) {
            this->shape_.set_value(shape);
            this->update();
        }
        void set_scale(double shape) {
            this->scale_.set_value(scale);
            this->update();
        }

        void store_shape() {
            this->shape_.store();
        }

        void restore_shape() {
            this->shape_.restore();
            this->update();
        }

        void store_scale() {
            this->scale_.store();
        }

        void restore_scale() {
            this->scale_.restore();
            this->update();
        }

        void store() {
            this->shape_.store();
            this->scale_.store();
        }
        void restore() {
            this->shape_.restore();
            this->scale_.restore();
            this->update();
        }

        void set_shape_preserve_mean(double shape) {
            double mean = this->get_mean();
            new_scale = mean / shape;
            this->shape_.set_value(shape);
            this->scale_.set_value(new_scale);
            this->update();
        }
        void set_scale_preserve_mean(double scale) {
            double mean = this->get_mean();
            new_shape = mean / scale;
            this->shape_.set_value(new_shape);
            this->scale_.set_value(scale);
            this->update();
        }

        const std::vector<double>& get_quantiles() const {
            return this->quantiles_;
        }

        double get_quantile(unsigned int bin_index) const {
            return this->quantiles_.at(bin_index);
        }

};

#endif
