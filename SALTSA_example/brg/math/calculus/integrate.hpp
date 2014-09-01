/**********************************************************************\
 @file integrate.hpp
 ------------------

 Functions to be used for integrating functions in various manners.

 **********************************************************************

 Copyright (C) 2014  Bryan R. Gillis

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

#ifndef _BRG_INTEGRATE_HPP_INCLUDED_
#define _BRG_INTEGRATE_HPP_INCLUDED_

#include <cstdlib>
#include <cmath>

#include "brg/global.h"

#include "brg/math/functor/functor_product.hpp"
#include "brg/math/misc_math.hpp"
#include "brg/math/safe_math.hpp"
#include "brg/utility.hpp"
#ifdef _BRG_USE_UNITS
#include "brg_units.h"
#endif

namespace brgastro
{

// Uses trapezoid-rule integration to estimate the integral of a function. Each output parameter is integrated independantly. For multiple input parameters,
// the function works iteratively, using the "passed" parameters seen at the end of the function. These parameters should not be entered by the user unless
// you're sure you know what you're doing. Due to the iterative nature of this function and the overhead involved, it may be unfeasibly slow for large-
// dimensional integration. In the case of num_in_params > ~4, Monte Carlo integration is typically superior.
// Note: For most smooth functions (of num_in_params < ~4), the integrate_Rhomberg function works better. This function is the superior choice for functions with // discontinuities cusps, corners, etc. It also has the benefit that the time spent is predefined by the input parameters, unlike the Rhomberg method which
// must find convergence, so there is no worry about facing a bizarre function which may take surprisingly long to integrate.
//
// Parameters:
// in_params_step: (first version only) The size of the steps used when integrating. Smaller is more accurate, but slower (order 1/in_params_step time).
// num_samples: (second version only) The number of steps used when integrating. Larger is more accurate, but slower (order num_samples time).
// num_passed_in_params & passed_in_params: Ignore these unless you know what you're doing.

// Scalar-in, scaler-out version.
template< typename f, typename T >
inline T integrate_trapezoid( const f * func, const T & min_in_param, const T & max_in_param,
		const T & in_param_step,
		T & out_param, const bool silent = false )
{
	T in_param( 0 );
	T temp_out_param( 0 );
	T last_out_param( 0 );

	bool first_step = true;
	int num_steps;

	// Calculate number of steps for integration
	num_steps = (int)( ( max_in_param - min_in_param )
			/ safe_d( in_param_step ) ) + 1;

	// Standard trapezoid rule integration routine now

	for ( int i = 0; i < num_steps; i++ )
	{
		in_param = min_in_param + in_param_step * i;

		// If we have output params from last time, shift them to the last_out_param array
		last_out_param = temp_out_param;

		// Call function at this value
		temp_out_param = ( *func )( in_param, silent );

		if(first_step)
			first_step = false;
		else
			out_param += ( last_out_param + temp_out_param )
				* in_param_step / 2.;

	} // for( int i = 0; i < num_steps[0]; i++ )

	return out_param;
}

// Vector-in, vector-out version
template< typename f, typename T >
inline std::vector<T> integrate_trapezoid( const f * func,
		const std::vector< T > & min_in_params,
		const std::vector< T > & max_in_params,
		const std::vector< T > & in_params_step,
		const std::vector< T > & passed_in_params = std::vector< T >( 0 ),
		const bool silent = false )
{
	std::vector< T > in_params( 0 );
	std::vector< T > new_min_in_params( 0 );
	std::vector< T > new_max_in_params( 0 );
	std::vector< T > new_in_params_step( 0 );
	std::vector< T > new_passed_in_params( 0 );
	std::vector< T > temp_out_params( 0 );
	std::vector< T > last_out_params( 0 );
	std::vector< T > out_params( 0 );

	typename std::vector<T>::size_type num_in_params=min_in_params.size();
	typename std::vector<T>::size_type num_passed_in_params=passed_in_params.size();
	typename std::vector<T>::size_type num_out_params=(*func)(in_params,true).size();

	bool array_created = false; // So we can only create the out_params array once, after the first step
	std::vector< int > num_steps;
	int param_starting_index;
	int new_num_in_params = 0, new_num_passed_params = 0, num_tot_params =
			num_in_params + num_passed_in_params;

	// Check that we have a sane number of input parameters
	if ( ( num_in_params < 1 ) || ( num_in_params > MAX_STACK_DEPTH )
			|| ( num_in_params != min_in_params.size() )
			|| ( num_in_params != max_in_params.size() )
			|| ( num_in_params != in_params_step.size() )
			|| ( num_passed_in_params != passed_in_params.size() ) )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Bad number of input params passed to integrate().\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	if ( int errcode = make_array( num_steps, num_in_params ) )
		return errcode + LOWER_LEVEL_ERROR;
	if ( int errcode = make_array( in_params, num_tot_params ) )
		return errcode + LOWER_LEVEL_ERROR;

	// Delete out_params array if it exists
	del_array( out_params );

	// Calculate number of steps for integration
	for ( unsigned int i = 0; i < num_in_params; i++ )
	{
		num_steps[i] = (int)( ( max_in_params[i] - min_in_params[i] )
				/ safe_d( in_params_step[i] ) ) + 1;
	}

	// Were any params passed in from a previous iteration?
	if ( num_passed_in_params > 0 )
	{
		// Fill up in_params with the passed parameters
		for ( unsigned int j = 0; j < num_passed_in_params; j++ )
		{
			in_params[j] = passed_in_params[j];
		} // for( int j = 0; j < num_passed_params; j++ )
	} // if ( num_passed_params > 0 )

	param_starting_index = num_passed_in_params; // Set index for parameter we'll be integrating over

	if ( num_in_params == 1 )     // if (num_in_params < 1)
	{
		// Standard trapezoid rule integration routine now

		array_created = false;
		for ( int i = 0; i < num_steps[0]; i++ )
		{
			in_params[param_starting_index] = min_in_params[0]
					+ in_params_step[0] * i;

			// If we have output params from last time, shift them to the last_out_params array
			if ( temp_out_params.size() > 0 )
			{
				for ( unsigned int j = 0; j < num_out_params; j++ )
					last_out_params[j] = temp_out_params[j];
				del_array( temp_out_params );
			}

			// Call function at this value
			temp_out_params = ( *func )( in_params, silent );

			// Create output param arrays if necessary
			if ( !array_created )
			{
				make_array( out_params, num_out_params );
				make_array( last_out_params, num_out_params );
				array_created = true;
			} // If this is the first time, we don't do anything. Wait till next round to start adding in
			else
			{
				// Update the output parameters with those from the function call usind trapezoidal rule
				for ( unsigned int j = 0; j < num_out_params; j++ )
				{
					out_params[j] += ( last_out_params[j] + temp_out_params[j] )
									* in_params_step[0] / 2.;
				}

			}

		} // for( int i = 0; i < num_steps[0]; i++ )

	}
	else if ( num_in_params > 1 )   // else if(num_in_params == 1)
	{
		// In this case, we're going to have to iterate, calling the integration function for each step to integrate other dimensions

		// Set up new passed parameter array
		new_num_passed_params = num_passed_in_params + 1;
		new_num_in_params = num_in_params - 1;
		make_array( new_passed_in_params, new_num_passed_params );
		for ( unsigned int i = 0; i < num_passed_in_params; i++ )
			new_passed_in_params[i] = passed_in_params[i];

		// Set up new in-parameter arrays, excluding this first parameter
		make_array( new_min_in_params, num_in_params - 1 );
		make_array( new_max_in_params, num_in_params - 1 );
		make_array( new_in_params_step, num_in_params - 1 );
		for ( unsigned int i = 0; i < num_in_params - 1; i++ )
		{
			new_min_in_params[i] = min_in_params[i + 1];
			new_max_in_params[i] = max_in_params[i + 1];
			new_in_params_step[i] = in_params_step[i + 1];
		} // for( int i = 0; i < num_in_params-1; i++)

		array_created = false;
		for ( int i = 0; i < num_steps[param_starting_index]; i++ )
		{
			// Determine input param and add it to passed parameters array
			new_passed_in_params[new_num_passed_params - 1] =
					min_in_params[param_starting_index]
							+ in_params_step[param_starting_index] * i;

			// Call integrate on remaining in_params
			integrate_trapezoid( func, new_min_in_params, new_max_in_params, new_in_params_step,
					temp_out_params, new_passed_in_params );

			// Create output param array if necessary
			if ( !array_created )
			{
				make_array( out_params, num_out_params );
				array_created = true;
			}

			// Update the output parameters with those from the integrate call
			for ( unsigned int j = 0; j < num_out_params; j++ )
			{
				out_params[j] += temp_out_params[j]
						* in_params_step[param_starting_index];
			}
		}

	}
	else     // else if (num_in_params > 1)
	{
		throw std::runtime_error("Invalid path!");
	} // else

	return out_params;
}

// Scalar-in, scalar-out version
template< typename f, typename T >
inline T integrate_trapezoid( const f * func, const T & min_in_params, const T & max_in_params, const int num_samples,
		const bool silent = false )
{
	T in_params_step( ( max_in_params - min_in_params )
				/ safe_d( num_samples - 1 ));
	return integrate_trapezoid( func, min_in_params, max_in_params,
			in_params_step, silent );
}

// Vector-in, vector-out version
template< typename f, typename T >
inline std::vector<T> integrate_trapezoid( const f * func,
		const std::vector< T > & min_in_params,
		const std::vector< T > & max_in_params, const int num_samples,
		const std::vector< T > & passed_in_params = std::vector< T >( 0 ),
		const bool silent = false )
{
	const typename std::vector<T>::size_type num_in_params(min_in_params.size());
	std::vector< T > in_params_step( num_in_params );
	for ( unsigned int i = 0; i < num_in_params; i++ )
	{
		in_params_step[i] = ( max_in_params[i] - min_in_params[i] )
				/ safe_d( num_samples - 1 );
	}
	return integrate_trapezoid( func, min_in_params, max_in_params,
			in_params_step, passed_in_params, silent );
}

// Scalar-in, scalar-out version
template< typename f1, typename f2, typename T >
inline T integrate_weighted_trapezoid( const f1 * func, const f2 * func_weight,
		const T & min_in_params, const T & max_in_params, const T & in_params_step,
		const bool silent = false )
{
	functor_product< f1, f2, T > fprod( func, func_weight );
	unsigned int num_prod_out_params = 0, num_weight_out_params = 0;
	T prod_out_params( 0 ), weight_out_params( 0 );

	prod_out_params = integrate_trapezoid( &fprod, min_in_params,
			max_in_params, in_params_step, silent);
	weight_out_params = integrate_trapezoid( func_weight, min_in_params,
			max_in_params, in_params_step, silent);

	return prod_out_params / safe_d( weight_out_params );
}

// Vector-in, vector-out version
template< typename f1, typename f2, typename T >
inline std::vector<T> integrate_weighted_trapezoid( const f1 * func, const f2 * func_weight,
		const std::vector< T > & min_in_params,
		const std::vector< T > & max_in_params,
		const std::vector< T > & in_params_step,
		const std::vector< T > & passed_in_params = std::vector< T >( 0 ),
		const bool silent = false )
{
	functor_product< f1, f2, T > fprod( func, func_weight );
	unsigned int num_prod_out_params = 0, num_weight_out_params = 0;
	std::vector< T > prod_out_params( 0 ), weight_out_params( 0 );

	prod_out_params = integrate_trapezoid( &fprod, min_in_params,
			max_in_params, in_params_step, passed_in_params, silent );
	weight_out_params = integrate_trapezoid( func_weight, min_in_params,
			max_in_params, in_params_step, passed_in_params, silent );

	std::vector<T> out_params( prod_out_params.size() );
	for ( unsigned int i = 0; i < prod_out_params.size(); i++ )
	{
		out_params[i] = prod_out_params[i] / safe_d( weight_out_params[i] );
	}

	return 0;
}

// Scalar-in, scalar-out version
template< typename f1, typename f2, typename T >
inline T integrate_weighted_trapezoid( const f1 * func, const f2 * func_weight,
		const T & min_in_params,
		const T & max_in_params, const int num_samples,
		const bool silent = false )
{
	functor_product< f1, f2, T > fprod( func, func_weight );
	unsigned int num_prod_out_params = 0, num_weight_out_params = 0;
	T prod_out_params( 0 ), weight_out_params( 0 );

	prod_out_params = integrate_trapezoid( &fprod, min_in_params,
			max_in_params, num_samples, silent );
	weight_out_params = integrate_trapezoid( func_weight, min_in_params,
			max_in_params, num_samples, silent );

	return prod_out_params / safe_d( weight_out_params );
}

// Vector-in, vector-out version
template< typename f1, typename f2, typename T >
inline std::vector<T> integrate_weighted_trapezoid( const f1 * func, const f2 * func_weight,
		const std::vector< T > & min_in_params,
		const std::vector< T > & max_in_params, const int num_samples,
		const std::vector< T > & passed_in_params = std::vector< T >( 0 ),
		const bool silent = false )
{
	functor_product< f1, f2, T > fprod( func, func_weight );
	std::vector< T > prod_out_params( 0 ), weight_out_params( 0 );

	prod_out_params = integrate( &fprod, min_in_params,
			max_in_params, num_samples, passed_in_params, silent );
	weight_out_params = integrate( func_weight, min_in_params,
			max_in_params, num_samples, passed_in_params, silent );

	std::vector<T> out_params( prod_out_params.size() );
	for ( unsigned int i = 0; i < prod_out_params.size(); i++ )
	{
		out_params[i] = prod_out_params[i] / safe_d( weight_out_params[i] );
	}

	return 0;
}

// Monte-carlo integration method. NYI
template< typename f, typename T >
inline std::vector< T > integrate_mc( const f * func,
		const std::vector< T > & min_in_params,
		const std::vector< T > & max_in_params, const int num_samples,
		const bool silent = false )
{
	bool first_sample = true;
	std::vector<T> test_in_params, test_out_params, out_params;

	test_in_params = drand(min_in_params,max_in_params);
	out_params = (*func)(test_in_params,silent);

	for(unsigned int i=0; i<num_samples-1; ++i )
	{
		test_in_params = drand(min_in_params,max_in_params);
		test_out_params = (*func)(test_in_params,silent);
		for(unsigned int j=0; j<test_out_params.size(); ++j)
			out_params[j] += test_out_params[j];
	}

	for(unsigned int j=0; j<out_params.size(); ++j)
		out_params[j] /= num_samples;

	return out_params;
}

// Uses Romberg's rule to integrate a function. Each output parameter is integrated independently. For multiple input parameters,
// the function works iteratively, using the "passed" parameters seen at the end of the function. These parameters should not be entered by the user unless
// you're sure you know what you're doing. Due to the iterative nature of this function and the overhead involved, it may be unfeasibly slow for large-
// dimensional integration. In the case of num_in_params > ~4, Monte Carlo integration is typically superior.
// Note: As this integration rule estimates a polynomial form for a function in order to integrate, it has difficulty with functions that have discontinities,
// cusps, corners, etc. The integrate() function above is typically better in these cases.
//
// Parameters:
// precision: The threshold for determining convergence, by comparing the difference of two successive estimates of the integral to their mean. Smaller means
//            more accurate, but longer time to compute.
// loosen_precision: If num_in_params > 1 and this is set to true, the function will accept a lower precision threshold for the sub-integrals performed. This
//                   will dramatically increase speed of the function, but decrease accuracy (whether the trade-off is worth it should be investigated)
// num_passed_in_params & passed_in_params: Ignore these unless you know what you're doing.

// Scalar-in, scalar-out version. !!! Still needs cleanup after testing
template< typename f, typename T >
inline T integrate_Romberg( const f * func,
		T min_in_param, T max_in_param, double precision = 0.00001,
		bool tighten_precision = false, bool silent = false )
{
	T in_param(0);
	std::vector< std::vector< T > > R( 0 );
	std::vector< T > Rn( 0 );
	T Rnm;
	T a0, b0;
	T fa( 0 ), fb( 0 ), ftot( 0 );
	T d;

	// Rhomberg's rule integration routine

	a0 = min_in_param;
	b0 = max_in_param;

	// Get R[0][0] first
	fa = ( *func )( a0, silent );
	fb = ( *func )( b0, silent );

	Rnm = 0.5 * ( b0 - a0 ) * ( fa + fb );

	Rn.push_back( Rnm );
	R.push_back( Rn );
	Rn.resize( 0 );

	for ( int n = 1; n < ROMBERG_N_MAX; n++ )
	{
		// Get R[n][0]

		ftot = 0;

		for ( int k = 1; k <= ipow( 2, n - 1 ); k++ )
		{
			in_param = a0
					+ ( 2 * k - 1 ) * ( b0 - a0 ) / ipow( 2, n );
			ftot += ( *func )( in_param, silent );
		}

		Rnm = 0.5 * R[n - 1][0] + ( b0 - a0 ) / ipow( 2, n ) * ftot;

		Rn.push_back( Rnm );

		for ( int m = 1; m <= n; m++ )
		{
			Rnm = ( ipow( 4, m ) * Rn[m - 1] - R[n - 1][m - 1] )
					/ ( ipow( 4, m ) - 1 );
			Rn.push_back( Rnm );
		}

		R.push_back( Rn );
		Rn.resize( 0 );

		// Check for convergence
		d = ( 2 * fabs( R[n][n] - R[n - 1][n - 1] )
						/ safe_d( fabs( R[n][n] + R[n - 1][n - 1] ) ) );
		if ( d < precision )
		{
			return R[n][n];
		}

	} // for(int n = 0; n < RHOMBERG_N_MAX; n++)

	if(!silent) std::cerr << "WARNING: Integrate_Romberg did not converge.\n";
	return R.back().back();
}

// Vector-in, vector-out version
template< typename f, typename T >
inline std::vector< T > integrate_Romberg( const f * func,
		const std::vector< T > & min_in_params,
		const std::vector< T > & max_in_params,
		const double precision = 0.00001,
		const bool tighten_precision = false,
		const std::vector< T > & passed_in_params = std::vector< T >( 0 ),
		const bool silent = false )
{
	std::vector< T > in_params( 0 );
	std::vector< T > new_in_params( 0 );
	std::vector< T > new_min_in_params( 0 );
	std::vector< T > new_max_in_params( 0 );
	std::vector< T > new_in_params_step( 0 );
	std::vector< T > new_passed_in_params( 0 );
	std::vector< T > temp_out_params( 0 );
	std::vector< T > last_out_params( 0 );
	std::vector< std::vector< std::vector< T > > > R( 0 );
	std::vector< std::vector< T > > Rn( 0 );
	std::vector< T > Rnm;
	std::vector< T > a( 0 ), b( 0 );
	T a0, b0;
	std::vector< T > fa( 0 ), fb( 0 ), ftot( 0 );
	T d;

	unsigned int param_starting_index;
	unsigned int num_in_params = min_in_params.size(), num_passed_in_params = passed_in_params.size(),
			new_num_in_params = 0, new_num_passed_params = 0, num_tot_params =
			num_in_params + num_passed_in_params;

	// Check that we have a sane number of input parameters
	if ( ( num_in_params < 1 ) || ( num_in_params > MAX_STACK_DEPTH )
			|| ( num_in_params != min_in_params.size() )
			|| ( num_in_params != max_in_params.size() )
			|| ( num_passed_in_params != passed_in_params.size() ) )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Bad number of input params passed to integrate_Rhomberg().\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	if ( int errcode = make_array( in_params, num_tot_params ) )
		return errcode + LOWER_LEVEL_ERROR;

	// Were any params passed in from a previous iteration?
	if ( num_passed_in_params > 0 )
	{
		// Fill up in_params with the passed parameters
		for ( unsigned int j = 0; j < num_passed_in_params; j++ )
		{
			in_params[j] = passed_in_params[j];
		} // for( int j = 0; j < num_passed_params; j++ )
	} // if ( num_passed_params > 0 )

	param_starting_index = num_passed_in_params; // Set index for parameter we'll be integrating over

	if ( num_in_params < 1 ) // To catch errors that might have slipped through
	{
		throw std::runtime_error("Invalid path!");
	}
	else if ( num_in_params == 1 )     // if (num_in_params < 1)
	{
		// Rhomberg's rule integration routine now

		a0 = min_in_params[0];
		b0 = max_in_params[0];

		// Get R[0][0] first

		in_params[param_starting_index] = a0;
		fa = ( *func )( in_params, silent );

		in_params[param_starting_index] = b0;
		fb = ( *func )( in_params, silent );

		unsigned int num_out_params = fa.size();

		Rnm.resize( num_out_params );

		for ( unsigned int i = 0; i < num_out_params; i++ )
			Rnm[i] = 0.5 * ( b0 - a0 ) * ( fa[i] + fb[i] );

		Rn.push_back( Rnm );
		R.push_back( Rn );
		Rn.resize( 0 );

		for ( int n = 1; n < ROMBERG_N_MAX; n++ )
		{
			// Get R[n][0]

			make_array( ftot, num_out_params );
			for ( int k = 1; k <= ipow( 2, n - 1 ); k++ )
			{
				in_params[param_starting_index] = a0
						+ ( 2 * k - 1 ) * ( b0 - a0 ) / ipow( 2, n );
				if ( int errcode = ( *func )( in_params, temp_out_params,
						silent ) )
					return errcode + LOWER_LEVEL_ERROR;
				for ( unsigned int i = 0; i < num_out_params; i++ )
					ftot[i] += temp_out_params[i];
			}

			for ( unsigned int i = 0; i < num_out_params; i++ )
				Rnm[i] = 0.5 * R[n - 1][0][i]
						+ ( b0 - a0 ) / ipow( 2, n ) * ftot[i];

			Rn.push_back( Rnm );

			for ( int m = 1; m <= n; m++ )
			{
				for ( unsigned int i = 0; i < num_out_params; i++ )
					Rnm[i] = ( ipow( 4, m ) * Rn[m - 1][i]
							- R[n - 1][m - 1][i] ) / ( ipow( 4, m ) - 1 );
				Rn.push_back( Rnm );
			}

			R.push_back( Rn );
			Rn.resize( 0 );

			// Check for convergence
			d = 0;
			for ( unsigned int i = 0; i < num_out_params; i++ )
			{
				d = quad_add( d,( 2 * fabs( R[n][n][i] - R[n - 1][n - 1][i] )
							/ safe_d( fabs( R[n][n][i] + R[n - 1][n - 1][i] ) ) ) );
			}
			if ( d < precision )
			{
				return R[n][n];
				break;
			}

		} // for(int n = 0; n < RHOMBERG_N_MAX; n++)

	} // else if(num_in_params == 1)
	else if ( num_in_params > 1 )
	{
		// In this case, we're going to have to iterate, calling the integration function for each step to integrate other dimensions

		// Set up new passed parameter array
		new_num_passed_params = num_passed_in_params + 1;
		new_num_in_params = num_in_params - 1;
		double new_precision;
		if ( tighten_precision )
			new_precision = std::pow( precision,
					(double)num_in_params / new_num_in_params );
		else
			new_precision = precision;
		if ( int errcode = make_array( new_passed_in_params,
				new_num_passed_params ) )
			return errcode + LOWER_LEVEL_ERROR;
		for ( unsigned int i = 0; i < num_passed_in_params; i++ )
			new_passed_in_params[i] = passed_in_params[i];

		// Set up new in-parameter arrays, excluding this first parameter
		if ( int errcode = make_array( new_min_in_params, num_in_params - 1 ) )
			return errcode + LOWER_LEVEL_ERROR;
		if ( int errcode = make_array( new_max_in_params, num_in_params - 1 ) )
			return errcode + LOWER_LEVEL_ERROR;
		for ( unsigned int i = 0; i < num_in_params - 1; i++ )
		{
			new_min_in_params[i] = min_in_params[i + 1];
			new_max_in_params[i] = max_in_params[i + 1];
		} // for( int i = 0; i < num_in_params-1; i++)

		a0 = min_in_params[param_starting_index];
		b0 = max_in_params[param_starting_index];

		// Determine input param and add it to passed parameters array
		new_passed_in_params[new_num_passed_params - 1] = a0;
		// Call integrate on remaining in_params
		fa = brgastro::integrate_Romberg( func,
				new_min_in_params, new_max_in_params,
				new_precision, tighten_precision,
				new_passed_in_params, silent );
		// Determine input param and add it to passed parameters array
		new_passed_in_params[new_num_passed_params - 1] = b0;
		// Call integrate on remaining in_params
		fb = brgastro::integrate_Romberg( func,
				new_min_in_params, new_max_in_params,
				new_precision, tighten_precision,
				new_passed_in_params, silent );

		unsigned int num_out_params = fa.size();

		Rnm.resize( num_out_params );

		for ( unsigned int i = 0; i < num_out_params; i++ )
			Rnm[i] = 0.5 * ( b0 - a0 ) * ( fa[i] + fb[i] );

		Rn.push_back( Rnm );
		R.push_back( Rn );
		Rn.resize( 0 );

		for ( int n = 1; n < ROMBERG_N_MAX; n++ )
		{
			// Get R[n][0]

			make_array( ftot, num_out_params );
			for ( int k = 1; k <= ipow( 2, n - 1 ); k++ )
			{
				new_passed_in_params[new_num_passed_params - 1] = a0
						+ ( 2 * k - 1 ) * ( b0 - a0 ) / ipow( 2., n );
				temp_out_params = brgastro::integrate_Romberg( func,
						new_min_in_params,
						new_max_in_params,
						new_precision, tighten_precision,
						new_passed_in_params, silent );
				for ( unsigned int i = 0; i < num_out_params; i++ )
					ftot[i] += temp_out_params[i];
			}

			for ( unsigned int i = 0; i < num_out_params; i++ )
				Rnm[i] = 0.5 * R[n - 1][0][i]
						+ ( b0 - a0 ) / ipow( 2, n ) * ftot[i];

			Rn.push_back( Rnm );

			for ( int m = 1; m <= n; m++ )
			{
				for ( unsigned int i = 0; i < num_out_params; i++ )
					Rnm[i] = ( ipow( 4, m ) * Rn[m - 1][i]
							- R[n - 1][m - 1][i] ) / ( ipow( 4, m ) - 1 );
				Rn.push_back( Rnm );
			}

			R.push_back( Rn );
			Rn.resize( 0 );

			// Check for convergence
			d = 0;
			for ( unsigned int i = 0; i < num_out_params; i++ )
			{
				d = quad_add( d,
						( 2 * fabs( R[n][n][i] - R[n - 1][n - 1][i] )
								/ fabs( R[n][n][i] + R[n - 1][n - 1][i] ) ) );
			}
			if ( d < precision )
			{
				return R[n][n];
			}
		}

	}
	else     // else if (num_in_params > 1)
	{
		throw std::runtime_error("Invalid path!");
	} // else

	if(!silent) std::cerr << "WARNING: Integrate_Romberg did not converge.\n";
	return R.back().back();
}

// Scalar-in, scalar-out version
template< typename f_in_1, typename f_in_2, typename T >
inline T integrate_product_Romberg( const f_in_1 * func1,
		const f_in_2 * func2, const T & min_in_params, const T & max_in_params,
		const double precision = 0.00001, const bool tighten_precision = false,
		const T & passed_in_params = T( 0 ), const bool silent = false )
{
	functor_product< f_in_1, f_in_2, T > fprod( func1, func2 );

	return integrate_Romberg( &fprod, min_in_params, max_in_params,
			precision, tighten_precision, passed_in_params, silent );
}

// Vector-in, vector-out version
template< typename f_in_1, typename f_in_2, typename T >
inline std::vector< T > integrate_product_Romberg( const f_in_1 * func1,
		const f_in_2 * func2, const std::vector< T > & min_in_params,
		const std::vector< T > & max_in_params, const double precision = 0.00001,
		const bool tighten_precision = false,
		const std::vector< T > & passed_in_params = std::vector< T >( 0 ),
		const bool silent = false )
{
	functor_product< f_in_1, f_in_2, T > fprod( func1, func2 );

	return integrate_Romberg( &fprod, min_in_params, max_in_params,
			precision, tighten_precision, passed_in_params, silent );
}

// Scalar-in, scalar-out version
template< typename f_in_1, typename f_in_2, typename T >
inline T integrate_weighted_Romberg( const f_in_1 * func,
		const f_in_2 * func_weight, const T & min_in_params, const T & max_in_params,
		const double precision = 0.00001, const bool tighten_precision = false,
		const bool silent = false )
{
	functor_product< f_in_1, f_in_2, T > fprod( func, func_weight );
	T prod_out_params( 0 ), weight_out_params( 0 );

	prod_out_params = integrate_Romberg( &fprod, min_in_params, max_in_params,
			precision, tighten_precision, silent);
	weight_out_params = integrate_Romberg( func_weight, min_in_params, max_in_params,
			precision, tighten_precision, silent);

	return prod_out_params / safe_d( weight_out_params );
}

// Vector-in, vector-out version
template< typename f_in_1, typename f_in_2, typename T >
inline std::vector< T > integrate_weighted_Romberg( const f_in_1 * func,
		const f_in_2 * func_weight,	const std::vector< T > & min_in_params,
		const std::vector< T > & max_in_params, const double precision = 0.00001,
		const bool tighten_precision = false,
		const std::vector< T > & passed_in_params = std::vector< T >( 0 ),
		const bool silent = false )
{
	functor_product< f_in_1, f_in_2, T > fprod( func, func_weight );
	std::vector< T > prod_out_params( 0 ), weight_out_params( 0 );

	prod_out_params = integrate_Romberg( &fprod, min_in_params, max_in_params,
			precision, tighten_precision, passed_in_params, silent );
	weight_out_params = integrate_Romberg( func_weight, min_in_params, max_in_params,
			precision, tighten_precision, passed_in_params, silent );

	std::vector< T > out_params( prod_out_params.size() );
	for ( unsigned int i = 0; i < prod_out_params.size(); i++ )
	{
		out_params[i] = prod_out_params[i] / safe_d( weight_out_params[i] );
	}

	return out_params;
}

} // namespace brgastro

#endif // _BRG_INTEGRATE_HPP_INCLUDED_
