/*
 * brg_calculus.hpp
 *
 *  Created on: 8 Apr 2014
 *      Author: brg
 */

#ifndef __BRG_CALCULUS_HPP_INCLUDED__
#define __BRG_CALCULUS_HPP_INCLUDED__

#include "brg_global.h"

#include <cstdlib>
#include <cmath>
#include "brg_units.h"
#include "brg_functions.h"

namespace brgastro
{

// Differentiates an arbitrary function numerically. The results are returned in the 2-D vector Jacobian, where Jacobian[i][j] represents the derivative of
// y_i with respect to x_j, with y_i being the output variables and x_j being the input variables, at the position labeled by in_params.
// In current implementation, the size of the differential used is a fraction of the input parameters (SMALL_FACTOR*in_params, where SMALL_FACTOR is defined
// in the brg_global.h header). If any in parameter is zero, the function uses the others as a guide for the size (and if using units, takes the units from
// the passed zero, so make sure your zeros have units if you do this!). If all in_params are zero, the function uses SMALL_FACTOR as the value. Be careful
// about this if evaluating a derivative at zero where the function changes on scales smaller than this - use a value slightly offset from zero instead.
//
// Parameters:
// order: Order of differentiation (1 = first derivative, 2 = second derivative, etc.). Order > 1 is NYI
// power: Used if, instead of the derivative of f(x), you want the derivative of (f(x))^2, for instance (set power = 2 for that), without setting up a
//        different function class.

// Scalar-in, scalar-out version !!! Needs clean-up after testing
template< typename f, typename T >
inline const int differentiate( const f * func,
		const unsigned int num_in_params, const T & in_params,
		unsigned int & num_out_params, T & out_params, T & Jacobian,
		const int order = 1, const double power = 1,
		const bool silent = false )
{

	BRG_UNITS d_in_params( 0 );
	BRG_UNITS base_out_params( 0 );
	BRG_UNITS test_in_params( 0 );
	BRG_UNITS test_out_params( 0 );
	BRG_UNITS small_factor_with_units = SMALL_FACTOR;

	bool power_flag = false;
	bool zero_in_flag = false;

	int order_to_use = (int)max( order, 1 );

	if ( ( order_to_use > 1 ) )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: brgastro::unit_quick_differentiate with order > 1 is not currently supported.\n";
		return UNSPECIFIED_ERROR;
	}

	if ( power != 1 )
		power_flag = true;
	else
		power_flag = false;

	// Check if any in_params are zero. If so, estimate small factor from other in_params
	for ( unsigned int i = 0; i < num_in_params; i++ )
	{
		if ( in_params == 0 )
		{
			zero_in_flag = true;
		}
		else     // if(in_params==0)
		{
			small_factor_with_units = in_params * SMALL_FACTOR;
			d_in_params = small_factor_with_units;
		} // else
	} // for( unsigned int i = 0; i < num_in_params; i++ )

	if ( zero_in_flag )
	{
		if ( small_factor_with_units == 0 )
		{
			// At least try to get the units right
			for ( unsigned int i = 0; i < num_in_params; i++ )
			{
#ifdef _BRG_USE_UNITS_
				d_in_params.set(SMALL_FACTOR,in_params.get_unit_powers());
#else
				d_in_params = SMALL_FACTOR;
#endif
			} // for( unsigned int i = 0; i < num_in_params; i++ )
		}
		else
		{
			for ( unsigned int i = 0; i < num_in_params; i++ )
			{
				if ( in_params == 0 )
				{
#ifdef _BRG_USE_UNITS_
					d_in_params.set(SMALL_FACTOR_units.get_value(),in_params.get_unit_powers());
#else
					d_in_params = small_factor_with_units;
#endif
				} // if(in_params[i]==0)
			} // for( unsigned int i = 0; i < num_in_params; i++ )
		}
	}

	// Get value of function at input parameters
	if ( int errcode = ( *func )( in_params, base_out_params, silent ) )
		return errcode + LOWER_LEVEL_ERROR;

	// Loop over input and output dimensions to get Jacobian
	for ( unsigned int j = 0; j < num_in_params; j++ )
	{
		// Set up test input parameters
		for ( unsigned int j2 = 0; j2 < num_in_params; j2++ )
		{
			if ( j2 == j )
			{
				test_in_params = in_params + d_in_params;
			} // if( j2==j )
			else
			{
				test_in_params = in_params;
			} // else
		}

		// Run the function to get value at test point
		if ( int errcode = ( *func )( test_in_params, test_out_params,
				silent ) )
			return errcode + LOWER_LEVEL_ERROR;

		// Record this derivative
		for ( unsigned int i = 0; i < num_out_params; i++ )
		{
			Jacobian = ( test_out_params - base_out_params ) / d_in_params;
			if ( power_flag )
				Jacobian *= power * safe_pow( base_out_params, power - 1 );
		} // for( int i = 0; i < num_out_params; i++)
	} // for( unsigned int j = 0; j < num_in_params; j++)

	return 0;
}

// Vector-in, vector-out version
template< typename f, typename T >
inline const int differentiate( const f * func,
		const unsigned int num_in_params, const std::vector< T > & in_params,
		unsigned int & num_out_params, std::vector< T > & out_params,
		std::vector< std::vector< T > > & Jacobian, const int order = 1,
		const double power = 1, const bool silent = false )
{

	std::vector< BRG_UNITS > d_in_params( 0 );
	std::vector< BRG_UNITS > base_out_params( 0 );
	std::vector< BRG_UNITS > test_in_params( 0 );
	std::vector< BRG_UNITS > test_out_params( 0 );
	BRG_UNITS small_factor_with_units = SMALL_FACTOR;

	bool power_flag = false;
	bool zero_in_flag = false;

	int order_to_use = (int)max( order, 1 );

	if ( ( order_to_use > 1 ) )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: brgastro::unit_quick_differentiate with order > 1 is not currently supported.\n";
		return UNSPECIFIED_ERROR;
	}

	if ( power != 1 )
		power_flag = true;
	else
		power_flag = false;

	// Delete std::vectors we'll be overwriting in case they previously existed
	del_array( out_params );
	del_array2d( Jacobian, num_in_params );

	// Set up differentials
	if ( int errcode = make_array( d_in_params, num_in_params ) )
		return errcode + LOWER_LEVEL_ERROR;
	if ( int errcode = make_array( test_in_params, num_in_params ) )
		return errcode + LOWER_LEVEL_ERROR;

	// Check if any in_params are zero. If so, estimate small factor from other in_params
	for ( unsigned int i = 0; i < num_in_params; i++ )
	{
		if ( in_params[i] == 0 )
		{
			zero_in_flag = true;
		}
		else     // if(in_params[i]==0)
		{
			small_factor_with_units = in_params[i] * SMALL_FACTOR;
			d_in_params[i] = small_factor_with_units;
		} // else
	} // for( unsigned int i = 0; i < num_in_params; i++ )

	if ( zero_in_flag )
	{
		if ( small_factor_with_units == 0 )
		{
			// At least try to get the units right
			for ( unsigned int i = 0; i < num_in_params; i++ )
			{
#ifdef _BRG_USE_UNITS_
				d_in_params[i].set(SMALL_FACTOR,in_params[i].get_unit_powers());
#else
				d_in_params[i] = SMALL_FACTOR;
#endif
			} // for( unsigned int i = 0; i < num_in_params; i++ )
		}
		else
		{
			for ( unsigned int i = 0; i < num_in_params; i++ )
			{
				if ( in_params[i] == 0 )
				{
#ifdef _BRG_USE_UNITS_
					d_in_params[i].set(SMALL_FACTOR_units.get_value(),in_params[i].get_unit_powers());
#else
					d_in_params[i] = small_factor_with_units;
#endif
				} // if(in_params[i]==0)
			} // for( unsigned int i = 0; i < num_in_params; i++ )
		}
	}

	// Get value of function at input parameters
	if ( int errcode = ( *func )( in_params, base_out_params, silent ) )
		return errcode + LOWER_LEVEL_ERROR;

	// Set up Jacobian
	if ( int errcode = make_array2d( Jacobian, num_out_params,
			num_in_params ) )
		return errcode + LOWER_LEVEL_ERROR;

	// Loop over input and output dimensions to get Jacobian
	for ( unsigned int j = 0; j < num_in_params; j++ )
	{
		// Set up test input parameters
		for ( unsigned int j2 = 0; j2 < num_in_params; j2++ )
		{
			if ( j2 == j )
			{
				test_in_params[j2] = in_params[j2] + d_in_params[j2];
			} // if( j2==j )
			else
			{
				test_in_params[j2] = in_params[j2];
			} // else
		}

		// Run the function to get value at test point
		if ( int errcode = ( *func )( test_in_params, test_out_params,
				silent ) )
			return errcode + LOWER_LEVEL_ERROR;

		// Record this derivative
		for ( unsigned int i = 0; i < num_out_params; i++ )
		{
			Jacobian[i][j] = ( test_out_params[i] - base_out_params[i] )
					/ d_in_params[j];
			if ( power_flag )
				Jacobian[i][j] *= power
						* safe_pow( base_out_params[i], power - 1 );
		} // for( int i = 0; i < num_out_params; i++)
	} // for( unsigned int j = 0; j < num_in_params; j++)

	del_array( base_out_params );
	del_array( test_out_params );
	del_array( d_in_params );
	return 0;
}

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

// Scalar-in, scaler-out version. !!! Needs clean-up after testing
template< typename f, typename T >
inline const int integrate( const f * func, const unsigned int num_in_params,
		T & min_in_params, T & max_in_params, T & in_params_step,
		unsigned int & num_out_params, T & out_params,
		const unsigned int num_passed_in_params = 0,
		const T & passed_in_params = 0, const bool silent = false )
{
	T in_params( 0 );
	T temp_out_params( 0 );
	T last_out_params( 0 );

	bool array_created = false; // So we can only create the out_params array once, after the first step
	int num_steps;

	// Calculate number of steps for integration
	for ( unsigned int i = 0; i < num_in_params; i++ )
	{
		num_steps = (int)( ( max_in_params - min_in_params )
				/ safe_d( in_params_step ) ) + 1;
	}

	// Were any params passed in from a previous iteration?
	if ( num_passed_in_params > 0 )
	{
		// Fill up in_params with the passed parameters
		for ( unsigned int j = 0; j < num_passed_in_params; j++ )
		{
			in_params = passed_in_params;
		} // for( int j = 0; j < num_passed_params; j++ )
	} // if ( num_passed_params > 0 )

	// Standard trapezoid rule integration routine now

	array_created = false;
	for ( int i = 0; i < num_steps; i++ )
	{
		in_params = min_in_params + in_params_step * i;

		// If we have output params from last time, shift them to the last_out_params array
		last_out_params = temp_out_params;

		// Call function at this value
		if ( int errcode = ( *func )( in_params, temp_out_params, silent ) )
			return errcode + LOWER_LEVEL_ERROR;

		// Create output param arrays if necessary
		if ( !array_created )
		{
			array_created = true;
		} // If this is the first time, we don't do anything. Wait till next round to start adding in
		else
		{

			// Update the output parameters with those from the function call using trapezoidal rule
			for ( unsigned int j = 0; j < num_out_params; j++ )
			{
				out_params += ( last_out_params + temp_out_params )
						* in_params_step / 2.;
			}

		}

	} // for( int i = 0; i < num_steps[0]; i++ )

	return 0;
}

// Vector-in, vector-out version
template< typename f, typename T >
inline const int integrate( const f * func, const unsigned int num_in_params,
		const std::vector< T > & min_in_params,
		const std::vector< T > & max_in_params,
		const std::vector< T > & in_params_step, unsigned int & num_out_params,
		std::vector< T > & out_params,
		const unsigned int num_passed_in_params = 0,
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

	if ( num_in_params < 1 ) // To catch errors that might have slipped through
	{
		return errorNOS( silent );
	}
	else if ( num_in_params == 1 )     // if (num_in_params < 1)
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
			if ( int errcode = ( *func )( in_params, temp_out_params,
					silent ) )
				return errcode + LOWER_LEVEL_ERROR;

			// Create output param arrays if necessary
			if ( !array_created )
			{
				if ( int errcode = make_array( out_params, num_out_params ) )
					return errcode + LOWER_LEVEL_ERROR;
				if ( int errcode = make_array( last_out_params,
						num_out_params ) )
					return errcode + LOWER_LEVEL_ERROR;
				array_created = true;
			} // If this is the first time, we don't do anything. Wait till next round to start adding in
			else
			{

				// Update the output parameters with those from the function call usind trapezoidal rule
				for ( unsigned int j = 0; j < num_out_params; j++ )
				{
					out_params[j] +=
							( last_out_params[j] + temp_out_params[j] )
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
		if ( int errcode = make_array( new_in_params_step,
				num_in_params - 1 ) )
			return errcode + LOWER_LEVEL_ERROR;
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
			if ( int errcode = integrate( func, new_num_in_params,
					new_min_in_params, new_max_in_params, new_in_params_step,
					num_out_params, temp_out_params, new_num_passed_params,
					new_passed_in_params ) )
				return errcode + LOWER_LEVEL_ERROR;

			// Create output param array if necessary
			if ( !array_created )
			{
				if ( int errcode = make_array( out_params, num_out_params ) )
					return errcode + LOWER_LEVEL_ERROR;
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
		return errorNOS( silent );
	} // else

	return 0;
}

// Scalar-in, scalar-out version
template< typename f, typename T >
inline const int integrate( const f * func, const unsigned int num_in_params,
		T & min_in_params, T & max_in_params, const int num_samples,
		unsigned int & num_out_params, T & out_params,
		const unsigned int num_passed_in_params = 0,
		const T & passed_in_params = 0, const bool silent = false )
{

	T in_params_step( 0 );

	int result;
	for ( unsigned int i = 0; i < num_in_params; i++ )
	{
		in_params_step = ( max_in_params - min_in_params )
				/ safe_d( num_samples - 1 );
	}
	result = integrate( func, num_in_params, min_in_params, max_in_params,
			in_params_step, num_out_params, out_params, num_passed_in_params,
			passed_in_params );
	return result;
}

// Vector-in, vector-out version
template< typename f, typename T >
inline const int integrate( const f * func, const unsigned int num_in_params,
		const std::vector< T > & min_in_params,
		const std::vector< T > & max_in_params, const int num_samples,
		unsigned int & num_out_params, std::vector< T > & out_params,
		const unsigned int num_passed_in_params = 0,
		const std::vector< T > & passed_in_params = std::vector< T >( 0 ),
		const bool silent = false )
{

	std::vector< T > in_params_step( 0 );

	int result;
	if ( int errcode = make_array( in_params_step, num_in_params ) )
		return errcode + LOWER_LEVEL_ERROR;
	for ( unsigned int i = 0; i < num_in_params; i++ )
	{
		in_params_step[i] = ( max_in_params[i] - min_in_params[i] )
				/ safe_d( num_samples - 1 );
	}
	result = integrate( func, num_in_params, min_in_params, max_in_params,
			in_params_step, num_out_params, out_params, num_passed_in_params,
			passed_in_params );
	del_array( in_params_step );
	return result;
}

// Scalar-in, scalar-out version
template< typename f1, typename f2, typename T >
inline const int integrate_weighted( const f1 * func, const f2 * func_weight,
		const unsigned int num_in_params, const T & min_in_params,
		const T & max_in_params, const T & in_params_step,
		unsigned int & num_out_params, T & out_params,
		const unsigned int num_passed_in_params = 0,
		const T & passed_in_params = T( 0 ), const bool silent = false )
{
	function_product_function< f1, f2, T > fprod( func, func_weight );
	unsigned int num_prod_out_params = 0, num_weight_out_params = 0;
	T prod_out_params( 0 ), weight_out_params( 0 );

	if ( int errcode = integrate( &fprod, num_in_params, min_in_params,
			max_in_params, in_params_step, num_prod_out_params,
			prod_out_params, num_passed_in_params, passed_in_params ) )
		return errcode + LOWER_LEVEL_ERROR;
	if ( int errcode = integrate( func_weight, num_in_params, min_in_params,
			max_in_params, in_params_step, num_weight_out_params,
			weight_out_params, num_passed_in_params, passed_in_params ) )
		return errcode + LOWER_LEVEL_ERROR;

	num_out_params = num_prod_out_params; // By construction of the function_product_function class, this must be the same as num_weight_out_params
	out_params = prod_out_params / safe_d( weight_out_params );

	return 0;
}

// Vector-in, vector-out version
template< typename f1, typename f2, typename T >
inline const int integrate_weighted( const f1 * func, const f2 * func_weight,
		const unsigned int num_in_params,
		const std::vector< T > & min_in_params,
		const std::vector< T > & max_in_params,
		const std::vector< T > & in_params_step, unsigned int & num_out_params,
		std::vector< T > & out_params,
		const unsigned int num_passed_in_params = 0,
		const std::vector< T > & passed_in_params = std::vector< T >( 0 ),
		const bool silent = false )
{
	function_product_function< f1, f2, T > fprod( func, func_weight );
	unsigned int num_prod_out_params = 0, num_weight_out_params = 0;
	std::vector< T > prod_out_params( 0 ), weight_out_params( 0 );

	if ( int errcode = integrate( &fprod, num_in_params, min_in_params,
			max_in_params, in_params_step, num_prod_out_params,
			prod_out_params, num_passed_in_params, passed_in_params ) )
		return errcode + LOWER_LEVEL_ERROR;
	if ( int errcode = integrate( func_weight, num_in_params, min_in_params,
			max_in_params, in_params_step, num_weight_out_params,
			weight_out_params, num_passed_in_params, passed_in_params ) )
		return errcode + LOWER_LEVEL_ERROR;

	num_out_params = num_prod_out_params; // By construction of the function_product_function class, this must be the same as num_weight_out_params
	out_params.resize( num_out_params );
	for ( unsigned int i = 0; i < num_out_params; i++ )
	{
		out_params[i] = prod_out_params[i] / safe_d( weight_out_params[i] );
	}

	return 0;
}

// Scalar-in, scalar-out version
template< typename f1, typename f2, typename T >
inline const int integrate_weighted( const f1 * func, const f2 * func_weight,
		const unsigned int num_in_params, const T & min_in_params,
		const T & max_in_params, const int num_samples,
		unsigned int & num_out_params, T & out_params,
		const unsigned int num_passed_in_params = 0,
		const T & passed_in_params = T( 0 ), const bool silent = false )
{
	function_product_function< f1, f2, T > fprod( func, func_weight );
	unsigned int num_prod_out_params = 0, num_weight_out_params = 0;
	T prod_out_params( 0 ), weight_out_params( 0 );

	if ( int errcode = integrate( &fprod, num_in_params, min_in_params,
			max_in_params, num_samples, num_prod_out_params, prod_out_params,
			num_passed_in_params, passed_in_params ) )
		return errcode + LOWER_LEVEL_ERROR;
	if ( int errcode = integrate( func_weight, num_in_params, min_in_params,
			max_in_params, num_samples, num_weight_out_params,
			weight_out_params, num_passed_in_params, passed_in_params ) )
		return errcode + LOWER_LEVEL_ERROR;

	num_out_params = num_prod_out_params; // By construction of the function_product_function class, this must be the same as num_weight_out_params
	out_params.resize( num_out_params );
	out_params = prod_out_params / safe_d( weight_out_params );

	return 0;
}

// Vector-in, vector-out version
template< typename f1, typename f2, typename T >
inline const int integrate_weighted( const f1 * func, const f2 * func_weight,
		const unsigned int num_in_params,
		const std::vector< T > & min_in_params,
		const std::vector< T > & max_in_params, const int num_samples,
		unsigned int & num_out_params, std::vector< T > & out_params,
		const unsigned int num_passed_in_params = 0,
		const std::vector< T > & passed_in_params = std::vector< T >( 0 ),
		const bool silent = false )
{
	function_product_function< f1, f2, T > fprod( func, func_weight );
	unsigned int num_prod_out_params = 0, num_weight_out_params = 0;
	std::vector< T > prod_out_params( 0 ), weight_out_params( 0 );

	if ( int errcode = integrate( &fprod, num_in_params, min_in_params,
			max_in_params, num_samples, num_prod_out_params, prod_out_params,
			num_passed_in_params, passed_in_params ) )
		return errcode + LOWER_LEVEL_ERROR;
	if ( int errcode = integrate( func_weight, num_in_params, min_in_params,
			max_in_params, num_samples, num_weight_out_params,
			weight_out_params, num_passed_in_params, passed_in_params ) )
		return errcode + LOWER_LEVEL_ERROR;

	num_out_params = num_prod_out_params; // By construction of the function_product_function class, this must be the same as num_weight_out_params
	out_params.resize( num_out_params );
	for ( unsigned int i = 0; i < num_out_params; i++ )
	{
		out_params[i] = prod_out_params[i] / safe_d( weight_out_params[i] );
	}

	return 0;
}

// Monte-carlo integration method. NYI
template< typename f, typename T >
inline const int integrate_mc( const f * func,
		const unsigned int num_in_params,
		const std::vector< T > & min_in_params,
		const std::vector< T > & max_in_params, const int num_samples,
		unsigned int & num_out_params, std::vector< T > & out_params,
		const bool silent = false )
{
	if ( !silent )
		std::cerr
				<< "ERROR: Placeholder function for integrate_mc being used.\n";
	return UNSPECIFIED_ERROR;
}

// Uses Rhomberg's rule to integrate a function. Each output parameter is integrated independently. For multiple input parameters,
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
inline const int integrate_Rhomberg( const f * func,
		const unsigned int num_in_params, const T & min_in_params,
		const T & max_in_params, unsigned int & num_out_params, T & out_params,
		const double precision = 0.00001, const bool tighten_precision = false,
		const unsigned int num_passed_in_params = 0,
		const T & passed_in_params = T( 0 ), const bool silent = false )
{
	T in_params( 0 );
	T temp_out_params( 0 );
	std::vector< std::vector< T > > R( 0 );
	std::vector< T > Rn( 0 );
	T Rnm;
	T a0, b0;
	T fa( 0 ), fb( 0 ), ftot( 0 );
	T d;

	// Check that we have a sane number of input parameters
	if ( ( num_in_params < 1 ) || ( num_in_params > MAX_STACK_DEPTH ) )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Bad number of input params passed to integrate_Rhomberg().\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	if ( num_in_params < 1 ) // To catch errors that might have slipped through
	{
		return errorNOS( silent );
	}
	else if ( num_in_params == 1 )     // if (num_in_params < 1)
	{
		// Rhomberg's rule integration routine now

		a0 = min_in_params;
		b0 = max_in_params;

		// Get R[0][0] first

		in_params = a0;
		if ( int errcode = ( *func )( in_params, temp_out_params, silent ) )
			return errcode + LOWER_LEVEL_ERROR;
		fa = temp_out_params;

		in_params = b0;
		if ( int errcode = ( *func )( in_params, temp_out_params, silent ) )
			return errcode + LOWER_LEVEL_ERROR;
		fb = temp_out_params;

		Rnm = 0.5 * ( b0 - a0 ) * ( fa + fb );

		Rn.push_back( Rnm );
		R.push_back( Rn );
		Rn.resize( 0 );

		for ( int n = 1; n < RHOMBERG_N_MAX; n++ )
		{
			// Get R[n][0]

			for ( int k = 1; k <= std::pow( 2, n - 1 ); k++ )
			{
				in_params = a0
						+ ( 2 * k - 1 ) * ( b0 - a0 ) / std::pow( 2., n );
				if ( int errcode = ( *func )( in_params, temp_out_params,
						silent ) )
					return errcode + LOWER_LEVEL_ERROR;
				ftot += temp_out_params;
			}

			Rnm = 0.5 * R[n - 1][0] + ( b0 - a0 ) / std::pow( 2, n ) * ftot;

			Rn.push_back( Rnm );

			for ( int m = 1; m <= n; m++ )
			{
				Rnm = ( std::pow( 4, m ) * Rn[m - 1] - R[n - 1][m - 1] )
						/ ( std::pow( 4, m ) - 1 );
				Rn.push_back( Rnm );
			}

			R.push_back( Rn );
			Rn.resize( 0 );

			// Check for convergence
			d = 0;
			d = quad_add( d,
					( 2 * fabs( R[n][n] - R[n - 1][n - 1] )
							/ safe_d( fabs( R[n][n] + R[n - 1][n - 1] ) ) ) );
			if ( d < precision )
			{
				out_params = R[n][n];
				break;
			}

		} // for(int n = 0; n < RHOMBERG_N_MAX; n++)

	}
	else
	{
		return errorNOS();
	}

	return 0;
}

// Vector-in, vector-out version
template< typename f, typename T >
inline const int integrate_Rhomberg( const f * func,
		const unsigned int num_in_params,
		const std::vector< T > & min_in_params,
		const std::vector< T > & max_in_params, unsigned int & num_out_params,
		std::vector< T > & out_params, const double precision = 0.00001,
		const bool tighten_precision = false,
		const unsigned int num_passed_in_params = 0,
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

	int param_starting_index;
	int new_num_in_params = 0, new_num_passed_params = 0, num_tot_params =
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

	// Delete out_params array if it exists
	del_array( out_params );

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
		return errorNOS( silent );
	}
	else if ( num_in_params == 1 )     // if (num_in_params < 1)
	{
		// Rhomberg's rule integration routine now

		a0 = min_in_params[0];
		b0 = max_in_params[0];

		// Get R[0][0] first

		in_params[param_starting_index] = a0;
		if ( int errcode = ( *func )( in_params, temp_out_params, silent ) )
			return errcode + LOWER_LEVEL_ERROR;
		fa = temp_out_params;

		in_params[param_starting_index] = b0;
		if ( int errcode = ( *func )( in_params, temp_out_params, silent ) )
			return errcode + LOWER_LEVEL_ERROR;
		fb = temp_out_params;

		Rnm.resize( num_out_params );

		for ( unsigned int i = 0; i < num_out_params; i++ )
			Rnm[i] = 0.5 * ( b0 - a0 ) * ( fa[i] + fb[i] );

		Rn.push_back( Rnm );
		R.push_back( Rn );
		Rn.resize( 0 );

		for ( int n = 1; n < RHOMBERG_N_MAX; n++ )
		{
			// Get R[n][0]

			make_array( ftot, num_out_params );
			for ( int k = 1; k <= std::pow( 2, n - 1 ); k++ )
			{
				in_params[param_starting_index] = a0
						+ ( 2 * k - 1 ) * ( b0 - a0 ) / std::pow( 2., n );
				if ( int errcode = ( *func )( in_params, temp_out_params,
						silent ) )
					return errcode + LOWER_LEVEL_ERROR;
				for ( unsigned int i = 0; i < num_out_params; i++ )
					ftot[i] += temp_out_params[i];
			}

			for ( unsigned int i = 0; i < num_out_params; i++ )
				Rnm[i] = 0.5 * R[n - 1][0][i]
						+ ( b0 - a0 ) / std::pow( 2, n ) * ftot[i];

			Rn.push_back( Rnm );

			for ( int m = 1; m <= n; m++ )
			{
				for ( unsigned int i = 0; i < num_out_params; i++ )
					Rnm[i] = ( std::pow( 4, m ) * Rn[m - 1][i]
							- R[n - 1][m - 1][i] ) / ( std::pow( 4, m ) - 1 );
				Rn.push_back( Rnm );
			}

			R.push_back( Rn );
			Rn.resize( 0 );

			// Check for convergence
			d = 0;
			for ( unsigned int i = 0; i < num_out_params; i++ )
			{
				d =
						quad_add( d,
								( 2 * fabs( R[n][n][i] - R[n - 1][n - 1][i] )
										/ safe_d(
												fabs(
														R[n][n][i]
																+ R[n - 1][n
																		- 1][i] ) ) ) );
			}
			if ( d < precision )
			{
				out_params = R[n][n];
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
		if ( int errcode = brgastro::integrate_Rhomberg( func,
				new_num_in_params, new_min_in_params, new_max_in_params,
				num_out_params, fa, new_precision, tighten_precision,
				new_num_passed_params, new_passed_in_params ) )
			return errcode + LOWER_LEVEL_ERROR;
		// Determine input param and add it to passed parameters array
		new_passed_in_params[new_num_passed_params - 1] = b0;
		// Call integrate on remaining in_params
		if ( int errcode = brgastro::integrate_Rhomberg( func,
				new_num_in_params, new_min_in_params, new_max_in_params,
				num_out_params, fb, new_precision, tighten_precision,
				new_num_passed_params, new_passed_in_params ) )
			return errcode + LOWER_LEVEL_ERROR;

		Rnm.resize( num_out_params );

		for ( unsigned int i = 0; i < num_out_params; i++ )
			Rnm[i] = 0.5 * ( b0 - a0 ) * ( fa[i] + fb[i] );

		Rn.push_back( Rnm );
		R.push_back( Rn );
		Rn.resize( 0 );

		for ( int n = 1; n < RHOMBERG_N_MAX; n++ )
		{
			// Get R[n][0]

			make_array( ftot, num_out_params );
			for ( int k = 1; k <= std::pow( 2, n - 1 ); k++ )
			{
				new_passed_in_params[new_num_passed_params - 1] = a0
						+ ( 2 * k - 1 ) * ( b0 - a0 ) / std::pow( 2., n );
				if ( int errcode = brgastro::integrate_Rhomberg( func,
						new_num_in_params, new_min_in_params,
						new_max_in_params, num_out_params, temp_out_params,
						new_precision, tighten_precision,
						new_num_passed_params, new_passed_in_params ) )
					return errcode + LOWER_LEVEL_ERROR;
				for ( unsigned int i = 0; i < num_out_params; i++ )
					ftot[i] += temp_out_params[i];
			}

			for ( unsigned int i = 0; i < num_out_params; i++ )
				Rnm[i] = 0.5 * R[n - 1][0][i]
						+ ( b0 - a0 ) / std::pow( 2, n ) * ftot[i];

			Rn.push_back( Rnm );

			for ( int m = 1; m <= n; m++ )
			{
				for ( unsigned int i = 0; i < num_out_params; i++ )
					Rnm[i] = ( std::pow( 4, m ) * Rn[m - 1][i]
							- R[n - 1][m - 1][i] ) / ( std::pow( 4, m ) - 1 );
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
				out_params = R[n][n];
				break;
			}
		}

	}
	else     // else if (num_in_params > 1)
	{
		return errorNOS( silent );
	} // else

	return 0;
}

// Scalar-in, scalar-out version
template< typename f_in_1, typename f_in_2, typename T >
inline const int integrate_weighted_Rhomberg( const f_in_1 * func,
		const f_in_2 * func_weight, const unsigned int num_in_params,
		const T & min_in_params, const T & max_in_params,
		unsigned int & num_out_params, T & out_params, const double precision =
				0.00001, const bool tighten_precision = false,
		const unsigned int num_passed_in_params = 0,
		const T & passed_in_params = T( 0 ), const bool silent = false )
{
	function_product_function< f_in_1, f_in_2, T > fprod( func, func_weight );
	unsigned int num_prod_out_params = 0, num_weight_out_params = 0;
	T prod_out_params( 0 ), weight_out_params( 0 );

	if ( int errcode = integrate_Rhomberg( &fprod, num_in_params,
			min_in_params, max_in_params, num_prod_out_params, prod_out_params,
			precision, tighten_precision, num_passed_in_params,
			passed_in_params ) )
		return errcode + LOWER_LEVEL_ERROR;
	if ( int errcode = integrate_Rhomberg( func_weight, num_in_params,
			min_in_params, max_in_params, num_weight_out_params,
			weight_out_params, precision, tighten_precision,
			num_passed_in_params, passed_in_params ) )
		return errcode + LOWER_LEVEL_ERROR;

	num_out_params = num_prod_out_params; // By construction of the function_product_function class, this must be the same as num_weight_out_params
	out_params = prod_out_params / safe_d( weight_out_params );

	return 0;
}

// Vector-in, vector-out version
template< typename f_in_1, typename f_in_2, typename T >
inline const int integrate_weighted_Rhomberg( const f_in_1 * func,
		const f_in_2 * func_weight, const unsigned int num_in_params,
		const std::vector< T > & min_in_params,
		const std::vector< T > & max_in_params, unsigned int & num_out_params,
		std::vector< T > & out_params, const double precision = 0.00001,
		const bool tighten_precision = false,
		const unsigned int num_passed_in_params = 0,
		const std::vector< T > & passed_in_params = std::vector< T >( 0 ),
		const bool silent = false )
{
	function_product_function< f_in_1, f_in_2, T > fprod( func, func_weight );
	unsigned int num_prod_out_params = 0, num_weight_out_params = 0;
	std::vector< T > prod_out_params( 0 ), weight_out_params( 0 );

	if ( int errcode = integrate_Rhomberg( &fprod, num_in_params,
			min_in_params, max_in_params, num_prod_out_params, prod_out_params,
			precision, tighten_precision, num_passed_in_params,
			passed_in_params ) )
		return errcode + LOWER_LEVEL_ERROR;
	if ( int errcode = integrate_Rhomberg( func_weight, num_in_params,
			min_in_params, max_in_params, num_weight_out_params,
			weight_out_params, precision, tighten_precision,
			num_passed_in_params, passed_in_params ) )
		return errcode + LOWER_LEVEL_ERROR;

	num_out_params = num_prod_out_params; // By construction of the function_product_function class, this must be the same as num_weight_out_params
	out_params.resize( num_out_params );
	for ( unsigned int i = 0; i < num_out_params; i++ )
	{
		out_params[i] = prod_out_params[i] / safe_d( weight_out_params[i] );
	}

	return 0;
}

// Leapfrog method for solving a DE. Note that this implementation assumes that the positions and velocities passed to it are already spaced
// out by half a timestep, with velocity at t+t_step/2 (though it does allow phase classes to be passed to it). This method takes a single step,
// using the passed acceleration function. The passed function for this implementation must take in one parameter (the magnitude of distance from
// a centre point) and return one parameter (the magnitude of the acceleration toward this centre point).
template< typename f >
inline const int leapfrog_step( const BRG_DISTANCE &x, const BRG_DISTANCE &y,
		const BRG_DISTANCE &z, const BRG_VELOCITY &vx, const BRG_VELOCITY &vy,
		const BRG_VELOCITY &vz,
		BRG_DISTANCE & new_x, BRG_DISTANCE & new_y, BRG_DISTANCE & new_z,
		BRG_VELOCITY & new_vx, BRG_VELOCITY & new_vy, BRG_VELOCITY & new_vz,
		const BRG_TIME &t_step, const f *accel_func,
		const bool silent = false )
{
	BRG_DISTANCE d;
	BRG_UNITS a;

	d = 0;
	a = 0;

	// Adjust position
	new_x = x + vx * t_step;
	new_y = y + vy * t_step;
	new_z = z + vz * t_step;

	// Calculate acceleration at this new position
	d = dist3d( new_x, new_y, new_z );
	(*accel_func)( d, a, silent );

	// Adjust velocities
	new_vx = vx + a * new_x / d * t_step;
	new_vy = vy + a * new_y / d * t_step;
	new_vz = vz + a * new_z / d * t_step;

	return 0;
}

template< typename f >
inline const int leapfrog_step( BRG_DISTANCE & x, BRG_DISTANCE & y,
		BRG_DISTANCE & z,
		BRG_VELOCITY & vx, BRG_VELOCITY & vy, BRG_VELOCITY & vz,
		const BRG_TIME & t_step, const f *accel_func,
		const bool silent = false )
{
	BRG_DISTANCE new_x, new_y, new_z;
	BRG_VELOCITY new_vx, new_vy, new_vz;

	int result;
	result = leapfrog_step( x, y, z, vx, vy, vz, new_x, new_y, new_z, new_vx,
			new_vy, new_vz, t_step, accel_func, silent );
	x = new_x;
	y = new_y;
	z = new_z;
	vx = new_vx;
	vy = new_vy;
	vz = new_vz;
	return result;
}

template< typename f >
inline const int leapfrog_step( const phase &p, phase & new_p,
		const BRG_TIME &t_step, const f *accel_func,
		const bool silent = false )
{
	return leapfrog_step( p.x, p.y, p.z, p.vx, p.vy, p.vz, new_p.x, new_p.y,
			new_p.z, new_p.vx, new_p.vy, new_p.vz, t_step, accel_func, silent );
}

template< typename f >
inline const int leapfrog_step( phase & p, const BRG_TIME & t_step,
		const f *accel_func, const bool silent = false )
{
	int result;
	phase new_p;
	result = leapfrog_step( p, new_p, t_step, accel_func, silent );
	p = new_p;
	return result;
}

} // namespace brgastro

#endif // __BRG_CALCULUS_HPP_INCLUDED__
