/**********************************************************************\
  @file random_functions.hpp

 **********************************************************************

 Copyright (C) 2014, 2015  Bryan R. Gillis

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

\**********************************************************************/

// body file: random_functions.cpp

#ifndef _BRG_RANDOM_FUNCTIONS_HPP_INCLUDED_
#define _BRG_RANDOM_FUNCTIONS_HPP_INCLUDED_

#include <cmath>
#include <random>

#include "brg/global.h"
#include "brg/math/misc_math.hpp"

namespace brgastro
{

extern std::ranlux48 rng; // Initialised in random_functions.cpp

/** Global function declarations **/
#if (1)

// Generates a random int between min and max, inclusive of min, exclusive of max
template< typename T=int, typename T_in=int, typename T_gen=decltype(rng) >
T irand( T_in && min, T_in && max, T_gen & gen=rng )
{
	return std::uniform_int_distribution<T>(std::forward<T_in>(min),std::forward<T_in>(max)-1)(gen);
} // double drand(double min, double max)

// Generates a random double between min and max
template< typename T=double, typename T_gen=decltype(rng) >
T drand( T_gen & gen=rng )
{
	return std::uniform_real_distribution<T>()(gen);
}
template< typename T=double, typename T_in=double, typename T_gen=decltype(rng) >
T drand( T_in && min, T_in && max, T_gen & gen=rng )
{
	return std::uniform_real_distribution<T>(std::forward<T_in>(min),std::forward<T_in>(max))(gen);
} // double drand(double min, double max)

// Returns a random variable from a Gaussian distribution
template< typename T=double, typename T_gen=decltype(rng) >
T Gaus_rand( T_gen & gen=rng )
{
	return std::normal_distribution<T>()(gen);
} // double Gaus_rand()
template< typename T=double, typename T_mean=double, typename T_stddev=double,
		typename T_gen=decltype(rng) >
T Gaus_rand( T_mean && mean, T_stddev && stddev = 1.0, T_gen & gen=rng )
{
	return std::normal_distribution<T>(std::forward<T_mean>(mean),
			std::forward<T_stddev>(stddev))(gen);

} // double Gaus_rand(double mean, double stddev)

// Returns a random variable from a Gaussian distribution, truncated between min and max
template< typename T=double, typename T_mean=double, typename T_stddev=double,
		typename T_min=double, typename T_max=double, typename T_gen=decltype(rng) >
T trunc_Gaus_rand( T_mean && mean, T_stddev && stddev, T_min && min, T_max && max, T_gen & gen=rng )
{
	assert(max>min);

	// Try values until one works
	bool good_value = false;
	int attempt_counter = 0;

	while( (!good_value) and (attempt_counter < 1000) )
	{
		T test_result = Gaus_rand(mean,stddev,gen);
		if((test_result >= min)and(test_result <= max))
		{
			return test_result;
		}
		else
		{
			++attempt_counter;
		}
	}

	// Failsafe
	return (std::forward<T_min>(min)+std::forward<T_max>(max))/2.;

} // T trunc_Gaus_rand( T_in && mean, T_in && stddev, T_in && min, T_in && max, gen_t & gen=rng )

// Returns a random variable from a Gaussian distribution in log space
// Note that "mean" here is the desired mean, NOT the peak of the function (which differ in log space). If you want to use
// the peak, simply use the standard Gaus_rand version instead.
template< typename T=double, typename T_gen=decltype(rng) >
T log10Gaus_rand( T_gen & gen=rng )
{
	const double & fact = std::exp( -square( std::log( 10. ) ) / 2 );

	return ( fact * std::pow(10., Gaus_rand<T>(gen) ) );
} // double log10Gaus_rand()
template< typename T=double, typename T_mean=double, typename T_stddev=double,
		typename T_gen=decltype(rng) >
T log10Gaus_rand( T_mean && mean, T_stddev && stddev = 1., T_gen & gen=rng )
{
	const double & fact = std::exp( -square( stddev*std::log(10.) ) / 2 );

	return ( fact * std::pow(10., Gaus_rand<T,T_mean,T_stddev>(std::forward<T_mean>(mean),
			std::forward<T_stddev>(stddev),gen) ) );
} // double log10Gaus_rand(double mean, double stddev)

// Returns a random variable from a log10_Gaussian distribution, truncated between min and max
template< typename T=double, typename T_mean=double, typename T_stddev=double,
		typename T_min=double, typename T_max=double, typename T_gen=decltype(rng) >
T trunc_log10Gaus_rand( const T_mean & mean, const T_stddev & stddev, T_min && min, T_max && max,
		T_gen & gen=rng )
{
	assert(max>min);

	// Try values until one works
	bool good_value = false;
	int attempt_counter = 0;

	T p10_min = std::pow(10.,min);
	T p10_max = std::pow(10.,max);

	while( (!good_value) and (attempt_counter < 1000) )
	{
		T test_result = log10Gaus_rand(mean,stddev,gen);
		if((test_result >= p10_min)and(test_result <= p10_max))
		{
			return test_result;
		}
		else
		{
			++attempt_counter;
		}
	}

	// Failsafe
	return std::pow(10.,(std::forward<T_min>(min)+std::forward<T_max>(max))/2.);

} // T trunc_Gaus_rand( T_in && mean, T_in && stddev, T_in && min, T_in && max, gen_t & gen=rng )

// Returns a random variable from a Rayleigh distribution
template< typename T=double, typename T_in=double, typename T_gen=decltype(rng) >
T Rayleigh_rand( T_in && sigma=1., T_gen & gen=rng )
{
	return std::forward<T_in>(sigma)*std::sqrt(-2.*std::log(drand<T>(gen)));
}

// Returns a random variable from a Rayleigh distribution, truncated between min and max
template< typename T=double, typename T_sigma=double, typename T_max=double,
		typename T_gen=decltype(rng) >
T trunc_Rayleigh_rand( const T_sigma & sigma, T_max && max, T_gen & gen=rng )
{
	assert(max>0.);

	// Try values until one works
	bool good_value = false;
	int attempt_counter = 0;

	while( (!good_value) and (attempt_counter < 1000) )
	{
		T test_result = Rayleigh_rand(sigma,gen);
		if(test_result <= max)
		{
			return test_result;
		}
		else
		{
			++attempt_counter;
		}
	}

	// Failsafe
	return std::forward<T_max>(max)/2.;

} // T trunc_Rayleigh_rand( T_sigma && sigma, T_max && max, T_gen & gen=rng )

// Get a Rayleigh random variable, smoothly contracted to always be less than max
template< typename T=double, typename T_sigma=double, typename T_max=double,
		typename T_p=double, typename T_gen=decltype(rng) >
T contracted_Rayleigh_rand( T_sigma && sigma, T_max && max, const T_p & p, T_gen & gen=rng )
{
	// Generate an initial random Rayleigh variable
	T first_result = Rayleigh_rand(std::forward<T_sigma>(sigma),gen);

	// Transform it via Bryan's formula to rein in large values to be less than the max_val
	return (first_result / std::pow(1 + std::pow(first_result / std::forward<T_max>(max), p),
			1.0 / p));
}

// Returns a Poisson random variable.
template< typename T=int, typename T_gen=decltype(rng) >
T Pois_rand( T_gen & gen=rng )
{
	return std::poisson_distribution<T>()(gen);
} // T Pois_rand( T_gen & gen=rng )
template< typename T=int, typename T_in=double, typename T_gen=decltype(rng) >
T Pois_rand( T_in && lambda=1., T_gen & gen=rng )
{
	return std::poisson_distribution<T>(std::forward<T_in>(lambda))(gen);
} // T Pois_rand( T_in && lambda=1., T_gen & gen=rng )

#endif // End global function declarations

}

#endif
