/*
 * brg_calculus.hpp
 *
 *  Created on: 8 Apr 2014
 *      Author: brg
 */

#ifndef __SALTSA_CALCULUS_HPP_INCLUDED__
#define __SALTSA_CALCULUS_HPP_INCLUDED__

#include "SALTSA_global.h"

#include <cstdlib>
#include <cmath>
#include "SALTSA_unitconvs.h"
#include "SALTSA_misc_functions.h"

namespace SALTSA
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

// Scalar-in, scalar-out version !! TO-DO some cleanup here
/**
 *
 * @param func
 * @param num_in_params
 * @param in_params
 * @param num_out_params
 * @param out_params
 * @param Jacobian
 * @param order
 * @param power
 * @param silent
 * @return
 */
template< typename f, typename T >
inline const int differentiate( const f * func,
		const unsigned int num_in_params, const T & in_params,
		unsigned int & num_out_params, T & out_params, T & Jacobian,
		const int order = 1, const double power = 1,
		const bool silent = false )
{

	double d_in_params( 0 );
	double base_out_params( 0 );
	double test_in_params( 0 );
	double test_out_params( 0 );
	double small_factor_with_units = SMALL_FACTOR;

	bool power_flag = false;
	bool zero_in_flag = false;

	int order_to_use = (int)max( order, 1 );

	if ( ( order_to_use > 1 ) )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: SALTSA::unit_quick_differentiate with order > 1 is not currently supported.\n";
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
				d_in_params = SMALL_FACTOR;

			} // for( unsigned int i = 0; i < num_in_params; i++ )
		}
		else
		{
			for ( unsigned int i = 0; i < num_in_params; i++ )
			{
				if ( in_params == 0 )
				{
					d_in_params = small_factor_with_units;
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
/**
 *
 * @param func
 * @param num_in_params
 * @param in_params
 * @param num_out_params
 * @param out_params
 * @param Jacobian
 * @param order
 * @param power
 * @param silent
 * @return
 */
template< typename f, typename T >
inline const int differentiate( const f * func,
		const unsigned int num_in_params, const std::vector< T > & in_params,
		unsigned int & num_out_params, std::vector< T > & out_params,
		std::vector< std::vector< T > > & Jacobian, const int order = 1,
		const double power = 1, const bool silent = false )
{

	std::vector< double > d_in_params( 0 );
	std::vector< double > base_out_params( 0 );
	std::vector< double > test_in_params( 0 );
	std::vector< double > test_out_params( 0 );
	double small_factor_with_units = SMALL_FACTOR;

	bool power_flag = false;
	bool zero_in_flag = false;

	int order_to_use = (int)max( order, 1 );

	if ( ( order_to_use > 1 ) )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: SALTSA::unit_quick_differentiate with order > 1 is not currently supported.\n";
		return UNSPECIFIED_ERROR;
	}

	if ( power != 1 )
		power_flag = true;
	else
		power_flag = false;

	// Delete std::vectors we'll be overwriting in case they previously existed
	out_params.clear();
	Jacobian.clear();

	// Set up differentials
	d_in_params.resize( num_in_params );
	test_in_params.resize( num_in_params );

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
				d_in_params[i] = SMALL_FACTOR;

			} // for( unsigned int i = 0; i < num_in_params; i++ )
		}
		else
		{
			for ( unsigned int i = 0; i < num_in_params; i++ )
			{
				if ( in_params[i] == 0 )
				{
					d_in_params[i] = small_factor_with_units;
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

	base_out_params.clear();
	test_out_params.clear();
	d_in_params.clear();
	return 0;
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
/**
 *
 * @param func
 * @param num_in_params
 * @param min_in_params
 * @param max_in_params
 * @param num_out_params
 * @param out_params
 * @param precision
 * @param tighten_precision
 * @param num_passed_in_params
 * @param passed_in_params
 * @param silent
 * @return
 */
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
/**
 *
 * @param func
 * @param num_in_params
 * @param min_in_params
 * @param max_in_params
 * @param num_out_params
 * @param out_params
 * @param precision
 * @param tighten_precision
 * @param num_passed_in_params
 * @param passed_in_params
 * @param silent
 * @return
 */
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
		if ( int errcode = SALTSA::integrate_Rhomberg( func,
				new_num_in_params, new_min_in_params, new_max_in_params,
				num_out_params, fa, new_precision, tighten_precision,
				new_num_passed_params, new_passed_in_params ) )
			return errcode + LOWER_LEVEL_ERROR;
		// Determine input param and add it to passed parameters array
		new_passed_in_params[new_num_passed_params - 1] = b0;
		// Call integrate on remaining in_params
		if ( int errcode = SALTSA::integrate_Rhomberg( func,
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
				if ( int errcode = SALTSA::integrate_Rhomberg( func,
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
/**
 *
 * @param func
 * @param func_weight
 * @param num_in_params
 * @param min_in_params
 * @param max_in_params
 * @param num_out_params
 * @param out_params
 * @param precision
 * @param tighten_precision
 * @param num_passed_in_params
 * @param passed_in_params
 * @param silent
 * @return
 */
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
/**
 *
 * @param func
 * @param func_weight
 * @param num_in_params
 * @param min_in_params
 * @param max_in_params
 * @param num_out_params
 * @param out_params
 * @param precision
 * @param tighten_precision
 * @param num_passed_in_params
 * @param passed_in_params
 * @param silent
 * @return
 */
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
/**
 *
 * @param x
 * @param y
 * @param z
 * @param vx
 * @param vy
 * @param vz
 * @param new_x
 * @param new_y
 * @param new_z
 * @param new_vx
 * @param new_vy
 * @param new_vz
 * @param t_step
 * @param accel_func
 * @param silent
 * @return
 */
template< typename f >
inline const int leapfrog_step( const double &x, const double &y,
		const double &z, const double &vx, const double &vy,
		const double &vz,
		double & new_x, double & new_y, double & new_z,
		double & new_vx, double & new_vy, double & new_vz,
		const double &t_step, const f *accel_func,
		const bool silent = false )
{
	double d;
	double a;

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

/**
 *
 * @param x
 * @param y
 * @param z
 * @param vx
 * @param vy
 * @param vz
 * @param t_step
 * @param accel_func
 * @param silent
 * @return
 */
template< typename f >
inline const int leapfrog_step( double & x, double & y,
		double & z,
		double & vx, double & vy, double & vz,
		const double & t_step, const f *accel_func,
		const bool silent = false )
{
	double new_x, new_y, new_z;
	double new_vx, new_vy, new_vz;

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

/**
 *
 * @param p
 * @param new_p
 * @param t_step
 * @param accel_func
 * @param silent
 * @return
 */
template< typename f >
inline const int leapfrog_step( const phase &p, phase & new_p,
		const double &t_step, const f *accel_func,
		const bool silent = false )
{
	return leapfrog_step( p.x, p.y, p.z, p.vx, p.vy, p.vz, new_p.x, new_p.y,
			new_p.z, new_p.vx, new_p.vy, new_p.vz, t_step, accel_func, silent );
}

/**
 *
 * @param p
 * @param t_step
 * @param accel_func
 * @param silent
 * @return
 */
template< typename f >
inline const int leapfrog_step( phase & p, const double & t_step,
		const f *accel_func, const bool silent = false )
{
	int result;
	phase new_p;
	result = leapfrog_step( p, new_p, t_step, accel_func, silent );
	p = new_p;
	return result;
}

} // namespace SALTSA

#endif // __SALTSA_CALCULUS_HPP_INCLUDED__
