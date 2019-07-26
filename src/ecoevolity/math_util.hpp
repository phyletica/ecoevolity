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

#ifndef ECOEVOLITY_MATH_UTIL_HPP
#define ECOEVOLITY_MATH_UTIL_HPP

#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <numeric>

#include "assert.hpp"
#include "error.hpp"

inline void normalize_log_likelihoods(std::vector<double>& v) {
    double mx = v.at(0);
    for (auto v_iter : v) {
        if (v_iter > mx) {
            mx = v_iter;
        }
    }

    for (unsigned int i = 0; i < v.size(); ++i) {
        v.at(i) -= mx;
    }

    double sum = 0.0;
    for (unsigned int i = 0; i < v.size(); ++i) {
        v.at(i) = std::exp(v.at(i));
        sum += v.at(i);
    }

    double t = 0.0;
    for (unsigned int i = 0; i < v.size(); ++i) {
        v.at(i) /= sum;
        t += v.at(i);
    }
    ECOEVOLITY_ASSERT_APPROX_EQUAL(t, 1.0);
}

/**
 * Calculate the expected number of categories under a Dirichlet process.
 *
 * Calculates and returns the expected number of categories for
 * 'number_of_elements' elements under a Dirichlet process controlled by
 * 'concentration' parameter.
 *
 * @note    Modified from `expNumTables` function of `util.h` from
 *          [`DPPDiv`](http://phylo.bio.ku.edu/content/tracy-heath-dppdiv)
 *          version 1.0b (Copyright Tracy Heath, Mark Holder, and John
 *          Huelsenback; licensed under GPL v3;
 *          <http://phylo.bio.ku.edu/content/tracy-heath-dppdiv>).
 */
inline double get_dpp_expected_number_of_categories(
        double concentration,
        unsigned int number_of_elements) {
    double expected_ncats =  0.0;
    for (unsigned int i = 1; i <= number_of_elements; ++i) {
        expected_ncats += (1.0 / (((double) i) - 1.0 + concentration));
    }
    return expected_ncats * concentration;
}


/**
 * Calculate the log of Pochhammer's symbol:
 *
 * (x)_n = x(x + 1)(x + 2)...(x + n - 1) = Gamma(x + n) / Gamma(x)
 *
 */

inline double ln_pochhammer(double x, unsigned int n) {
    return std::lgamma(x + n) - std::lgamma(x);
}


/**
 * Calculate the expected number of categories under a Pitman-Yor process.
 *
 * Calculates and returns the expected number of categories for
 * 'number_of_elements' elements under a Pitman-Yor process controlled by
 * 'concentration' and 'discount' parameters.
 *
 * From Equation 161 of Pitman (2002):
 *
 * Pitman, Jim. 2002. Combinatorial Stochastic Processes. Technical Report No.
 * 621. Lecture notes for St. Flour Course, July 2002.
 * http://www.stat.berkeley.edu/~pitman/621.pdf
 */

inline double get_pyp_expected_number_of_categories(
        double concentration,
        double discount,
        unsigned int number_of_elements) {
    ECOEVOLITY_ASSERT((discount >= 0.0) && (discount < 1.0));
    ECOEVOLITY_ASSERT(concentration > -discount);
    if (discount == 0.0) {
        return get_dpp_expected_number_of_categories(concentration,
                number_of_elements);
    }
    // Calculate first term on log scale, because Pochhammer factorials can be
    // huge. Then exponentiate after taking the ratio of two huge numbers.
    double ln_numerator = ln_pochhammer(concentration + discount, number_of_elements);
    double ln_denom = (std::log(discount) +
            ln_pochhammer(concentration + 1.0, number_of_elements - 1));
    double ln_term_1 = ln_numerator - ln_denom;
    return std::exp(ln_term_1) - (concentration / discount);
}


/**
 * Calculate the Dirichlet-process concentration parameter.
 *
 * Calculates and returns the Dirichlet-process concentration parameter that
 * has an expected number of categories equal to
 * 'expected_number_of_categories' for 'number_of_elements' elements.
 *
 * @note    Modified from `calculateFromPriorMean` function of `util.h` from
 *          [`DPPDiv`](http://phylo.bio.ku.edu/content/tracy-heath-dppdiv)
 *          version 1.0b (Copyright Tracy Heath, Mark Holder, and John
 *          Huelsenback; licensed under GPL v3;
 *          <http://phylo.bio.ku.edu/content/tracy-heath-dppdiv>).
 */
inline double get_dpp_concentration(
        double expected_number_of_categories,
        unsigned int number_of_elements,
        double increment = 0.1,
        double precision = 0.000001,
        double buffer = 0.001) {
    double c = precision;
    const double n = (double) number_of_elements;
    double e = expected_number_of_categories;
    double current_e = get_dpp_expected_number_of_categories(c,
            number_of_elements);

    ECOEVOLITY_ASSERT(e <= n);
    ECOEVOLITY_ASSERT(e >= 1.0);
    if (e > (n - buffer)) {
        e = n - buffer;
    }
    else if (e < (1.0 + buffer)) {
        e = 1.0 + buffer;
    }
    bool increase = false;
    if (current_e < e) {
        increase = true;
    }
    while (fabs(current_e - e) > precision) {
        if ((current_e < e) && increase) {
            c += increment;
        }
        else if ((current_e > e) && (! increase)) {
            c -= increment;
        }
        else if ((current_e < e) && (! increase)) {
            increment /= 2.0;
            increase = true;
            c += increment;
        }
        else {
            increment /= 2.0;
            increase = false;
            c -= increment;
        }
        current_e = get_dpp_expected_number_of_categories(c,
                number_of_elements);
    }
    if (c <= 0.0) {
        return precision;
    }
    return c;
}

/**
 * Calculate the scale parameter of a gamma hyper prior on the concentration
 * parameter of the Dirichlet process.
 */
inline double get_dpp_gamma_scale(
        double expected_number_of_categories,
        unsigned int number_of_elements,
        double shape,
        double increment = 0.1,
        double precision = 0.000001,
        double buffer = 0.001) {
    double concentration = get_dpp_concentration(
            expected_number_of_categories,
            number_of_elements,
            increment,
            precision,
            buffer);
    return concentration / shape;
}


/**
 * Calculate the Pitman-Yor process concentration parameter.
 *
 * Calculates and returns the Pitman-Yor process concentration parameter that
 * has an expected number of categories equal to
 * 'expected_number_of_categories' for 'number_of_elements' elements and a
 * given value of the 'discount' parameter.
 *
 * @note    Modified from `calculateFromPriorMean` function of `util.h` from
 *          [`DPPDiv`](http://phylo.bio.ku.edu/content/tracy-heath-dppdiv)
 *          version 1.0b (Copyright Tracy Heath, Mark Holder, and John
 *          Huelsenback; licensed under GPL v3;
 *          <http://phylo.bio.ku.edu/content/tracy-heath-dppdiv>).
 */
inline double get_pyp_concentration(
        double expected_number_of_categories,
        unsigned int number_of_elements,
        double discount,
        double increment = 0.1,
        double precision = 0.000001,
        double buffer = 0.001) {
    double c = precision;
    const double n = (double) number_of_elements;
    double e = expected_number_of_categories;
    double current_e = get_pyp_expected_number_of_categories(c, discount,
            number_of_elements);

    ECOEVOLITY_ASSERT(e <= n);
    ECOEVOLITY_ASSERT(e >= 1.0);
    if (e > (n - buffer)) {
        e = n - buffer;
    }
    else if (e < (1.0 + buffer)) {
        e = 1.0 + buffer;
    }
    bool increase = false;
    if (current_e < e) {
        increase = true;
    }
    while (fabs(current_e - e) > precision) {
        if ((current_e < e) && increase) {
            c += increment;
        }
        else if ((current_e > e) && (! increase)) {
            c -= increment;
        }
        else if ((current_e < e) && (! increase)) {
            increment /= 2.0;
            increase = true;
            c += increment;
        }
        else {
            increment /= 2.0;
            increase = false;
            c -= increment;
        }
        current_e = get_pyp_expected_number_of_categories(c, discount,
                number_of_elements);
    }
    if (c <= 0.0) {
        return precision;
    }
    return c;
}

/**
 * Calculate the scale parameter of a gamma hyper prior on the concentration
 * parameter of a Pitman-Yor process.
 */
inline double get_pyp_concentration_gamma_scale(
        double expected_number_of_categories,
        unsigned int number_of_elements,
        double discount,
        double shape,
        double increment = 0.1,
        double precision = 0.000001,
        double buffer = 0.001) {
    double concentration = get_pyp_concentration(
            expected_number_of_categories,
            number_of_elements,
            discount,
            increment,
            precision,
            buffer);
    return concentration / shape;
}


template <typename T>
inline double get_dpp_log_prior_probability(
        const std::vector<T>& partition,
        double concentration) {
    ECOEVOLITY_ASSERT(concentration > 0.0);
    double log_concentration = std::log(concentration);
    double log_p = 0.0;
    T current_element;
    std::map<T, unsigned int> subset_counts;
    for (unsigned int i = 0; i < partition.size(); ++i) {
        current_element = partition.at(i);
        if (subset_counts.count(current_element) < 1) {
            log_p += (log_concentration - std::log(concentration + i));
            subset_counts[current_element] = 1;
            continue;
        }
        log_p += (std::log(subset_counts[current_element]) - std::log(concentration + i));
        ++subset_counts[current_element];
    }
    return log_p;
}

inline double get_dpp_log_prior_probability(
        const std::string& partition,
        double concentration) {
    std::vector<char> partition_vector(partition.begin(), partition.end());
    return get_dpp_log_prior_probability<char>(partition_vector, concentration);
}

template <typename T>
inline double get_pyp_log_prior_probability(
        const std::vector<T>& partition,
        double concentration,
        double discount) {
    ECOEVOLITY_ASSERT((discount >= 0.0) && (discount < 1.0));
    ECOEVOLITY_ASSERT(concentration > -discount);
    double log_p = 0.0;
    T current_element;
    std::map<T, unsigned int> subset_counts;
    unsigned ncategories = 0;
    for (unsigned int i = 0; i < partition.size(); ++i) {
        current_element = partition.at(i);
        if (subset_counts.count(current_element) < 1) {
            log_p += (std::log(concentration + (discount * ncategories)) -
                    std::log(concentration + i));
            subset_counts[current_element] = 1;
            ++ncategories;
            continue;
        }
        log_p += (std::log(subset_counts[current_element] - discount) -
                std::log(concentration + i));
        ++subset_counts[current_element];
    }
    return log_p;
}

inline double get_pyp_log_prior_probability(
        const std::string& partition,
        double concentration,
        double discount) {
    std::vector<char> partition_vector(partition.begin(), partition.end());
    return get_pyp_log_prior_probability<char>(partition_vector,
            concentration, discount);
}

template <typename T>
inline double get_wdp_log_prior_probability(
        const std::vector<T>& partition,
        double concentration,
        double discount) {
    ECOEVOLITY_ASSERT((discount >= 0.0) && (discount < 1.0));
    ECOEVOLITY_ASSERT(concentration > -discount);
    double log_p = 0.0;
    T current_element;
    std::map<T, unsigned int> subset_counts;
    unsigned ncategories = 0;
    for (unsigned int i = 0; i < partition.size(); ++i) {
        current_element = partition.at(i);
        if (subset_counts.count(current_element) < 1) {
            log_p += (std::log(concentration + (discount * i)) -
                    std::log(concentration + i));
            subset_counts[current_element] = 1;
            ++ncategories;
            continue;
        }
        log_p += (std::log(subset_counts[current_element] - (discount * subset_counts[current_element])) -
                std::log(concentration + i));
        ++subset_counts[current_element];
    }
    return log_p;
}

inline double get_wdp_log_prior_probability(
        const std::string& partition,
        double concentration,
        double discount) {
    std::vector<char> partition_vector(partition.begin(), partition.end());
    return get_wdp_log_prior_probability<char>(partition_vector,
            concentration, discount);
}


// Recursion is much slower (~ 1000 times slower!) than dynamic programming
// approach below.
// inline unsigned long long stirling2_recurse(int n, int k) {
//     ECOEVOLITY_ASSERT((n > 0) && (k > 0) && (n >= k));
//     if (n == 0 || k == 1 || k == n) {
//         return 1;
//     }
//     return stirling2_recurse(n - 1, k - 1) + k * stirling2_recurse(n - 1, k);
// }

template <typename T>
inline T stirling2_base(int n, int k) {
    ECOEVOLITY_ASSERT((n >= 0) && (k >= 0) && (n >= k));
    if (k == 0) {
        return 0;
    }
    if ((n == k) || (k == 1)) {
        return 1;
    }

    int maxj = n-k;

    std::vector<T> s_numbers (maxj + 1, 1);

    for (int i = 2; i <= k; ++i) {
        for(int j = 1; j <= maxj; ++j) {
            // Check for multiplication overflow
            if (s_numbers.at(j-1) > (std::numeric_limits<T>::max() / i)) {
                throw EcoevolityNumericLimitError("Overflow during multiplication in stirling2_base");
            }
            // Check for addition overflow
            T addend = i * s_numbers.at(j-1);
            if ((addend > 0) &&
                    (s_numbers.at(j) > (std::numeric_limits<T>::max() - addend))) {
                throw EcoevolityNumericLimitError("Overflow during addition in stirling2_base");
            }
            s_numbers.at(j) += addend;
        }
    }

    return s_numbers.at(maxj);
}

inline unsigned long long stirling2(int n, int k) {
    return stirling2_base<unsigned long long>(n, k);
}

inline long double stirling2_float(int n, int k) {
    return stirling2_base<long double>(n, k);
}

template <typename T>
inline T bell_number_base(int n) {
    T b = 0;
    for (int i = 0; i < n; ++i) {
        // Check for addition overflow
        T addend = stirling2_base<T>(n, (i + 1));
        if ((addend > 0) &&
                (b > (std::numeric_limits<T>::max() - addend))) {
            throw EcoevolityNumericLimitError("Overflow during addition in bell_number_base");
        }
        b += addend;
    }
    return b;
}

inline unsigned long long bell_number(int n) {
    // max unsigned long long (and unsigned long): 18446744073709551615
    // Bell number 25:                             4638590332229999353
    // Bell number 26:                             49631246523618756274
    return bell_number_base<unsigned long long>(n);
}

inline long double bell_float(int n) {
    // max double:      1.79769e+308
    // max long double: 1.18973e+4932
    // sage: SetPartitions(218).cardinality() > 1.79769e+308
    // False
    // sage: SetPartitions(219).cardinality() > 1.79769e+308
    // True
    // sage: SetPartitions(2228).cardinality() > 1.18973e+4932
    // False
    // sage: SetPartitions(2229).cardinality() > 1.18973e+4932
    // True
    return bell_number_base<long double>(n);
}

template <typename T>
inline double get_uniform_model_log_prior_probability(
        const unsigned int number_of_elements,
        const unsigned int number_of_categories,
        const double split_weight) {
    ECOEVOLITY_ASSERT(split_weight > 0.0);
    ECOEVOLITY_ASSERT(number_of_elements > 0);
    ECOEVOLITY_ASSERT(number_of_categories > 0);
    ECOEVOLITY_ASSERT(number_of_categories <= number_of_elements);
    T denom = 0;
    for (unsigned int k = 1; k <= number_of_elements; ++k) {
        denom += stirling2_base<T>(number_of_elements, k) * std::pow(split_weight, (k - 1));
    }
    return ((number_of_categories - 1) * std::log(split_weight)) - std::log(denom);
}

inline double get_uniform_model_log_prior_probability(
        const unsigned int number_of_elements,
        const unsigned int number_of_categories,
        const double split_weight) {
    return get_uniform_model_log_prior_probability<long double>(
            number_of_elements,
            number_of_categories,
            split_weight);
}

/**
 * Get all possible integer partitions of n into m parts.
 *
 * @note    This is an implementation of Algorithm H on Page 2 of:
 *          Donald E. Knuth. Generating all partitions, 2004. Pre-fascicle 3B
 *          of The Art of Computer Programming, A draft of sections 7.2.1.4–5
 *          http://www-cs-faculty.stanford.edu/~knuth/fasc3b.ps.gz
 */
inline std::vector< std::vector<unsigned int> > get_integer_partitions(
        int n,
        int m) {
    ECOEVOLITY_ASSERT(n >= m);
    ECOEVOLITY_ASSERT(m > 1);
    std::vector< std::vector<unsigned int> > partitions;
    std::vector<int> a(m + 1, 1);
    a.at(0) = n - m + 1;
    a.back() = -1;
    while (true) {
        std::vector<unsigned int> p(a.begin(), a.end()-1);
        partitions.push_back(p);
        while (a.at(1) < (a.at(0) - 1)) {
            --a.at(0);
            ++a.at(1);
            std::vector<unsigned int> p(a.begin(), a.end()-1);
            partitions.push_back(p);
        }
        int j = 3;
        int s = a.at(0) + a.at(1) - 1;
        while (a.at(j-1) >= (a.at(0) - 1)) {
            s = s + a.at(j-1);
            ++j;
        }
        if (j > m) {
            break;
        }
        int x = a.at(j-1) + 1;
        a.at(j-1) = x;
        --j;
        while (j > 1) {
            a.at(j-1) = x;
            s = s - x;
            --j;
        }
        a.at(0) = s;
    }
    return partitions;
}

/**
 * Calculate log factorial.
 *
 * @note    Adapted from Daniel's post on stackoverflow.com
 *          (http://stackoverflow.com/a/30630392)
 */
inline double log_factorial(unsigned int x) {
    if (x == 0 || x == 1) return 0;
    if (x == 2) {
        return std::log(2); // can add more for efficiency
    }
    if (x > 100) {
        return x * std::log(x) - x; // Stirling's approximation
    }
    std::vector<unsigned int> lx(x);
    std::iota(lx.begin(), lx.end(), 1);
    std::vector<double> tx(x);
    std::transform(lx.cbegin(), lx.cend(), tx.begin(),
                   [] (unsigned int a) { return std::log(static_cast<double>(a)); });
    return std::accumulate(tx.cbegin(), tx.cend(), double {});
}

/**
 * Calculate log multinomial coefficient
 *
 * @note    Adapted from Daniel's post on stackoverflow.com
 *          (http://stackoverflow.com/a/30630392)
 */
inline double log_multinomial_coefficient(const std::vector<unsigned int>& il) {
    ECOEVOLITY_ASSERT(il.size() > 1);
    std::vector<double> denoms(il.size());
    std::transform(il.begin(), il.end(), denoms.begin(), log_factorial);
    unsigned int numerator = std::accumulate(il.begin(), il.end(), 0);
    double denominator = std::accumulate(denoms.begin(), denoms.end(), 0.0);
    return log_factorial(numerator) - denominator;
}


/*!
 * This function calculates the quantiles of a standard normal distribution.
 *
 * \brief Quantile of a standard normal distribution.
 * \param prob is the probability up to the quantile. 
 * \return Returns quantile value. 
 * \throws Does not throw an error.
 * \see Odeh, R. E. and J. O. Evans. 1974. The percentage points of the normal
 *      distribution. Applied Statistics, 22:96-97.
 * \see Wichura, M. J.  1988. Algorithm AS 241: The percentage points of the
 *      normal distribution. 37:477-484.
 * \see Beasley, JD & S. G. Springer. 1977. Algorithm AS 111: The percentage
 *      points of the normal distribution. 26:118-121.
 *
 * @note    From revbayes (commit c8bf96ec786)
 *          RbStatistics::Normal::quantile
 *          (https://github.com/revbayes/revbayes)
 */
inline double normal_quantile(double p) {
    
	double a0 = -0.322232431088;
	double a1 = -1.0;
	double a2 = -0.342242088547;
	double a3 = -0.0204231210245;
 	double a4 = -0.453642210148e-4;
 	double b0 = 0.0993484626060;
 	double b1 = 0.588581570495;
 	double b2 = 0.531103462366; 
 	double b3 = 0.103537752850; 
 	double b4 = 0.0038560700634;
	double p1 = ( p < 0.5 ? p : 1.0 - p);
	if (p1 < 1e-20)
        return (-9999.0);
	double y = std::sqrt( std::log(1.0/(p1*p1)) );   
	double z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
	return ( p < 0.5 ? -z : z );
}


/*!
 * Calculate the incomplete gamma integral.
 *
 * (1) series expansion     if (alpha>x || x<=1)
 * (2) continued fraction   otherwise
 * RATNEST FORTRAN by
 * Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
 * 19: 285-287 (AS32)
 *
 * @note    From revbayes (commit c8bf96ec786)
 *          RbMath::incompleteGamma
 *          (https://github.com/revbayes/revbayes)
 */
inline double incomplete_gamma(double x, double alpha, double scale) {
    
    double accurate = 1e-8, overflow = 1e30;
    double factor, gin, rn, a, b, an, dif, term;
    double pn0, pn1, pn2, pn3, pn4, pn5;
    
    if (x == 0.0) {
        return 0.0;
    }
    if (x < 0.0 || alpha <= 0.0) 
    {
        std::ostringstream s;
        s << "Cannot compute incomplete gamma function for x = " << x << ", alpha = " << alpha << "and scale = " << scale;
        throw EcoevolityError(s.str());
    }
    
    factor = std::exp(alpha * std::log(x) - x - scale);
    
    if (x > 1 && x >= alpha) {
        // continued fraction
        a = 1 - alpha;
        b = a + x + 1;
        term = 0;
        pn0 = 1;
        pn1 = x;
        pn2 = x + 1;
        pn3 = x * b;
        gin = pn2 / pn3;
        
        do {
            a++;
            b += 2;
            term++;
            an = a * term;
            pn4 = b * pn2 - an * pn0;
            pn5 = b * pn3 - an * pn1;
            
            if (pn5 != 0) {
                rn = pn4 / pn5;
                dif = fabs(gin - rn);
                if (dif <= accurate) {
                    if (dif <= accurate * rn) {
                        break;
                    }
                }
                
                gin = rn;
            }
            pn0 = pn2;
            pn1 = pn3;
            pn2 = pn4;
            pn3 = pn5;
            if (fabs(pn4) >= overflow) {
                pn0 /= overflow;
                pn1 /= overflow;
                pn2 /= overflow;
                pn3 /= overflow;
            }
        } while (true);
        gin = 1 - factor * gin;
    } else {
        // series expansion
        gin = 1;
        term = 1;
        rn = alpha;
        do {
            rn++;
            term *= x / rn;
            gin += term;
        }
        while (term > accurate);
        gin *= factor / alpha;
    }
    return gin;
}


/*!
 * This function calculates the quantile of a chi square distribution with v
 * degrees of freedom.
 *
 * \brief Quantile of a chi square distribution.
 * \param df is the degrees of freedom of the chi square.
 * \param prob is the probability up to the quantile. 
 * \return Returns quantile value (or -1 if in error). 
 * \throws Does not throw an error.
 *
 * @note    From revbayes (commit c8bf96ec786)
 *          RbStatistics::ChiSquare::quantile
 *          (https://github.com/revbayes/revbayes)
 */
inline double chi_square_quantile(double prob, double df)
{
    
	double 		e = 0.5e-6, aa = 0.6931471805, p = prob, g,
				xx, c, ch, a = 0.0, q = 0.0, p1 = 0.0, p2 = 0.0, t = 0.0, 
				x = 0.0, b = 0.0, s1, s2, s3, s4, s5, s6;
	
	if (p < 0.000002 || p > 0.999998 || df <= 0.0)
		return (-1.0);
    
	g = std::lgamma(df/2.0);
	xx = df/2.0;
	c = xx - 1.0;
	if (df >= -1.24 * std::log(p))
		goto l1;
	ch = std::pow((p * xx * exp(g + xx * aa)), 1.0/xx);
	if (ch-e < 0) 
		return (ch);
	goto l4;
	l1:
		if (df > 0.32)
			goto l3;
		ch = 0.4;   
		a = std::log(1.0-p);
	l2:
		q = ch;  
		p1 = 1.0+ch*(4.67+ch);  
		p2 = ch*(6.73+ch*(6.66+ch));
		t = -0.5+(4.67+2.0*ch)/p1 - (6.73+ch*(13.32+3.0*ch))/p2;
		ch -= (1.0-exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
		if (fabs(q/ch-1.0)-0.01 <= 0.0) 
			goto l4;
		else                       
			goto l2;
	l3: 
		x = normal_quantile(p);
		p1 = 0.222222/df;
		ch = df*std::pow((x*std::sqrt(p1)+1.0-p1), 3.0);
		if (ch > 2.2*df+6.0)
			ch = -2.0*(std::log(1.0-p)-c*std::log(0.5*ch)+g);
	l4:
        double last_improv = q - ch;
		q = ch; 
		p1 = 0.5*ch;
		if ((t = incomplete_gamma(p1, xx, g)) < 0.0)
        {
            std::cerr<<"\nError incomplete_gamma\n";
			return (-1.0);
		}
		p2 = p-t;
		t = p2*std::exp(xx*aa+g+p1-c*std::log(ch));   
		b = t/ch;  
		a = 0.5*t-b*c;
		s1 = (210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
		s2 = (420.0+a*(735.0+a*(966.0+a*(1141.0+1278.0*a))))/2520.0;
		s3 = (210.0+a*(462.0+a*(707.0+932.0*a)))/2520.0;
		s4 = (252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a)))/5040.0;
		s5 = (84.0+264.0*a+c*(175.0+606.0*a))/2520.0;
		s6 = (120.0+c*(346.0+127.0*c))/5040.0;
		ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
		if (fabs(q/ch-1.0) > e && fabs(q - ch) - fabs(last_improv) < e /* <- against flip-flop */) 
			goto l4;
    
		return (ch);
}


#endif
