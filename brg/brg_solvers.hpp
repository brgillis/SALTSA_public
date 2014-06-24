/*
 * brg_solvers.hpp
 *
 *  Created on:     8 April 2014
 *  Last modified:  8 April 2014
 *      Author: brg
 */

#ifndef __BRG_SOLVERS_HPP_INCLUDED__
#define __BRG_SOLVERS_HPP_INCLUDED__

#include "brg_global.h"

#include <cstdlib>
#include <cmath>
#include "brg_units.h"
#include "brg_functions.h"

namespace brgastro
{

// Various methods to solve a function (a child of function_class with an actual function defined must be used):

// solve_iterate: Only valid for functions of the form x = f(x). The defined function should be the f(x), and the solution
// value should be the x where f(x)=x. Note that this solver cannot solve all functions of this form, particularly when the
// f(x) is undefined for certain input values which can occur as output.
//
// Parameters:
// slowdown: Makes the solver slower to converge, but more likely to converge, by using the new test value as the mean of the
//           most recent output and the slowdown last inputs.
// precision: How close successive values must be to each other (in a fraction of the mean) before they're said to have converged.
//            Smaller values make the output more precise, but make the solver take longer.
// max_counter: Maximum number of loops the solver can take before it gives up.
//
// IMPORTANT: This solver returns a value of 0 on failure, and otherwise the solution.

template< typename f, typename T >
inline const T solve_iterate( const f * func, const T &init_param = 0,
		const int slowdown = 1, const double precision = 0.0001,
		const int max_counter = 10000, const bool silent = false )
{
	bool converged_flag;
	int counter = 0;

	BRG_UNITS new_value = init_param;
	BRG_UNITS result = 0;
	BRG_UNITS mean_value = 0;
	std::vector< BRG_UNITS > past_values( 0 );
	int slowdown_to_use;

	if ( slowdown < 0 )
		slowdown_to_use = 0;
	else
		slowdown_to_use = slowdown;
	past_values.resize( slowdown_to_use + 1, 0 );

	while ( counter < max_counter )
	{
		counter++;

		// Shift back values
		for ( int i = slowdown_to_use; i > 0; i-- )
		{
			past_values[i] = past_values[i - 1];
		}
		past_values[0] = new_value;

		// Find new_value, based on standard iteration
		if ( int errcode = ( *func )( new_value, new_value, silent ) )
		{
			throw errcode + LOWER_LEVEL_ERROR;
		}

		// First check if we have an equality. Likely will only happen if zero is a solution
		if ( new_value == past_values[0] )
		{
			counter = 0;
			break;
		}

		if ( slowdown_to_use > 0 )
		{
			mean_value = 0;
			mean_value += new_value / min( slowdown_to_use + 1, counter + 1 );
			for ( int i = 0; i < min( slowdown_to_use, counter ); i++ )
				mean_value += past_values[i]
						/ min( slowdown_to_use + 1, counter + 1 );
			new_value = mean_value;
		}

		converged_flag = true;
		// Check for convergence
		for ( int i = 0; i < slowdown_to_use + 1; i++ )
		{
			if ( std::fabs(
					2 * ( new_value - past_values[i] )
							/ safe_d( new_value + past_values[i] ) )
					> precision )
			{
				converged_flag = false;
				break;
			}
		}

		if ( converged_flag )
		{
			counter = 0;
			break;
		}
	}

	if ( counter != 0 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: solve_brgastro::unit_iterate did not converge.\n";
		return 0;
	}
	result = new_value;
	return result;

}

// solve_sd: Steepest-descent solver. Finds the minimum value of a function. May fail if it reaches input values which give undefined
//           output. Only valid for functions which return only one output value. Returns the best set of input parameters in the
//           result_in_params vector. Returns a value of 0 on success and 1 on failure.
//
// Parameters:
// precision: How close successive values must be to each other (in a fraction of the mean) before they're said to have converged.
//            Smaller values make the output more precise, but make the solver take longer.
// lambda: Relative size of the step size. This implementation scales lambda by the derivative at the starting point over the value
//         at the starting point, so less guesswork is needed in choosing a proper value. Be careful about setting lambda to 1 though,
//         as this will often lead the the function stepping immediately to an input of zero, which will be undefined for many functions.
// cusp_override_power: Normally, a steepest-descent algorithm can only handle functions where the derivative is zero at the minimum. To
//                      handle functions where the derivative is undefined, try changing this value. If the function is expected to have
//                      a corner at the minimum (for instance, you took the absolute value of the distance from a target value), set this
//                      to one. The program will then solve the square of the input function instead. If the function is expected to reach
//                      a cusp at the minimum, set this value to some higher value (the specific value needed will depend on the order of
//                      the cusp).
// max_steps: Maximum number of steps the solver can take before it gives up.

// Scalar-in, scalar-out version !!! needs cleaning
template< typename f, typename T >
const int solve_sd( const f * func, const unsigned int num_in_params,
		const T & init_in_params, T & result_in_params,
		const double precision = 0.00001, const double lambda = 0.1,
		const double cusp_override_power = 0, const int max_steps = 10000,
		const bool silent = false )
{

	BRG_UNITS Jacobian( 0 );
	BRG_UNITS current_in_params( 0 );
	BRG_UNITS last_in_params( 0 );
	BRG_UNITS out_params( 0 );
	BRG_UNITS in_params_to_use = init_in_params;
	BRG_UNITS lambda_norm( 0 );
	BRG_UNITS mag;

	const int lambda_shortening_intervals = 10;
	const double lambda_shortening_factor = 10;
	unsigned int num_out_params = 0;
	int num_steps_taken = 0;
	bool converged_flag = false;

	// Sanity checks on variables
	if ( num_in_params < 1 )
	{
		if ( !silent )
			std::cerr << "ERROR: Invalid num_in_params in solve_sd.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	const double lambda_to_use = max( SMALL_FACTOR,
			min( 1 / SMALL_FACTOR, lambda ) );
	const int max_steps_to_use = (int)max( 1, max_steps );

	// Get normalized step size
	if ( int errcode = differentiate( func, num_in_params, in_params_to_use,
			num_out_params, out_params, Jacobian, 1, cusp_override_power ) )
		return errcode + LOWER_LEVEL_ERROR;
	mag = 0;
	for ( unsigned int i = 0; i < num_in_params; i++ )
	{
		mag = quad_add( mag, Jacobian );
	}
	lambda_norm = lambda_to_use / mag;
	mag = 0;
	for ( unsigned int i = 0; i < num_in_params; i++ )
	{
		mag = quad_add( mag, in_params_to_use );
	}
	if ( mag == 0 )
		mag = 1;
	lambda_norm *= mag;

	// Initialize solver
	num_steps_taken = 0;
	current_in_params = in_params_to_use;

	// Solver

	while ( num_steps_taken < max_steps_to_use )
	{
		num_steps_taken++;
		last_in_params = current_in_params;

		// Get derivative at current location
		if ( int errcode = differentiate( func, num_in_params,
				current_in_params, num_out_params, out_params, Jacobian, 1,
				cusp_override_power ) )
			return errcode + LOWER_LEVEL_ERROR;

		// Make a step
		for ( unsigned int i = 0; i < num_in_params; i++ )
		{
			current_in_params -= lambda_norm * Jacobian;
		}

		// Check if we've converged
		converged_flag = true;

		for ( unsigned int i = 0; i < num_in_params; i++ )
		{
			if ( isnan( current_in_params ) )
			{
				result_in_params = current_in_params;
				if ( !silent )
					std::cerr
							<< "ERROR: Somehow got NaN for in_params in solve_sd.\n";
				return UNSPECIFIED_ERROR;
			}
			if ( std::fabs(
					2 * ( current_in_params - last_in_params )
							/ safe_d( current_in_params + last_in_params ) )
					> precision )
			{
				converged_flag = false;
				break;
			}
		}

		if ( converged_flag )
			break;

		// In case we get to high steps, start decreasing lambda_norm
		if ( divisible( num_steps_taken,
				max_steps_to_use / lambda_shortening_intervals ) )
		{
			lambda_norm /= lambda_shortening_factor;
		}
	}

	// Check if we converged
	if ( !converged_flag )
	{
		if ( !silent )
			std::cerr << "WARNING: Did not converge in solve_sd.\n";
		return UNSPECIFIED_ERROR;
	}
	result_in_params = current_in_params;
	return 0;

}

// Vector-in, vector-out version
template< typename f, typename T >
const int solve_sd( const f * func, const unsigned int num_in_params,
		const std::vector< T > & init_in_params,
		std::vector< T > & result_in_params, const double precision = 0.00001,
		const double lambda = 0.1, const double cusp_override_power = 0,
		const int max_steps = 10000, const bool silent = false )
{

	std::vector< std::vector< BRG_UNITS > > Jacobian( 0 );
	std::vector< BRG_UNITS > current_in_params( 0 );
	std::vector< BRG_UNITS > last_in_params( 0 );
	std::vector< BRG_UNITS > out_params( 0 );
	std::vector< BRG_UNITS > in_params_to_use = init_in_params;
	BRG_UNITS lambda_norm( 0 );
	BRG_UNITS mag;

	const int lambda_shortening_intervals = 10;
	const double lambda_shortening_factor = 10;
	unsigned int num_out_params = 0;
	int num_steps_taken = 0;
	bool converged_flag = false;

	// Sanity checks on variables
	if ( num_in_params < 1 )
	{
		if ( !silent )
			std::cerr << "ERROR: Invalid num_in_params in solve_sd.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	const double lambda_to_use = max( SMALL_FACTOR,
			min( 1 / SMALL_FACTOR, lambda ) );
	const int max_steps_to_use = (int)max( 1, max_steps );

	// Check in_params_to_use. If it's empty, use zero
	if ( in_params_to_use.size() < (unsigned)num_in_params )
	{
		if ( !in_params_to_use.empty() )
		{
			if ( !silent )
				std::cerr
						<< "WARNING: Size mismatch of num_in_params and in_params_to_use.size() in brgastro::unit_solve_sd.\n";
		}
		in_params_to_use.resize( num_in_params, 0 );
		for ( unsigned int i = 0; i < num_in_params; i++ )
		{
			in_params_to_use[i] = 0;
		}
	} // if( in_params_to_use.size() < num_in_params )

	// Get normalized step size
	if ( int errcode = differentiate( func, num_in_params, in_params_to_use,
			num_out_params, out_params, Jacobian, 1, cusp_override_power ) )
		return errcode + LOWER_LEVEL_ERROR;
	mag = 0;
	for ( unsigned int i = 0; i < num_in_params; i++ )
	{
		mag = quad_add( mag, Jacobian[0][i] );
	}
	lambda_norm = lambda_to_use / mag;
	mag = 0;
	for ( unsigned int i = 0; i < num_in_params; i++ )
	{
		mag = quad_add( mag, in_params_to_use[i] );
	}
	if ( mag == 0 )
		mag = 1;
	lambda_norm *= mag;

	// Initialize solver
	num_steps_taken = 0;
	current_in_params = in_params_to_use;

	// Solver

	while ( num_steps_taken < max_steps_to_use )
	{
		num_steps_taken++;
		last_in_params = current_in_params;

		// Get derivative at current location
		if ( int errcode = differentiate( func, num_in_params,
				current_in_params, num_out_params, out_params, Jacobian, 1,
				cusp_override_power ) )
			return errcode + LOWER_LEVEL_ERROR;

		// Make a step
		for ( unsigned int i = 0; i < num_in_params; i++ )
		{
			current_in_params[i] -= lambda_norm * Jacobian[0][i];
		}

		// Check if we've converged
		converged_flag = true;

		for ( unsigned int i = 0; i < num_in_params; i++ )
		{
			if ( isnan( current_in_params[i] ) )
			{
				result_in_params = current_in_params;
				if ( !silent )
					std::cerr
							<< "ERROR: Somehow got NaN for in_params in solve_sd.\n";
				return UNSPECIFIED_ERROR;
			}
			if ( std::fabs(
					2 * ( current_in_params[i] - last_in_params[i] )
							/ safe_d(
									current_in_params[i]
											+ last_in_params[i] ) )
					> precision )
			{
				converged_flag = false;
				break;
			}
		}

		if ( converged_flag )
			break;

		// In case we get to high steps, start decreasing lambda_norm
		if ( divisible( num_steps_taken,
				max_steps_to_use / lambda_shortening_intervals ) )
		{
			lambda_norm /= lambda_shortening_factor;
		}
	}

	// Check if we converged
	if ( !converged_flag )
	{
		if ( !silent )
			std::cerr << "WARNING: Did not converge in solve_sd.\n";
		return UNSPECIFIED_ERROR;
	}
	result_in_params = current_in_params;
	return 0;

}

// solve_grid: Custom solver which searches a grid of input parameters for the ones which give output closest to the target. The best point is then
//             used as the starting point for a narrowing search on progressively finer grids of 3^n points (n = number of input parameters).
//             This solver has the advantage that it can better handle regions of input space which give undefined results, and it is less likely to
//             get stuck in local minima. It is, however, slower than other methods, especially for high numbers of input parameters.
//             This function is set up to handle multiple output parameters. If it is more important that one be close to the target than another,
//             use the out_params_weight vector to assign it a higher weight.
//             Returns 0 on success, 1 on failure (usually only happens if the entire input space gives undefined output or invalid values are
//             passed to the function).
//
// Parameters:
// min/max_in_params: Limits to the search space for a solution
// search_in_params_step: (first version only) Step size used in the initial grid search for the optimum position
// num_search_steps: (second version only) Number of steps to use in the initial grid search (usually around 10 is good).
// target_out_params: Desired output parameters.
// precision: How small must the step size be relative to the steps in the initial grid search before the search is ended.
// search_precision: In the event the initial grid search fails to find any values that return defined results, the algorithm will try one more
//                   time after shrinking the grid step size by this factor.
// out_params_weight: If there are multiple output parameters, this represents how relatively important it is that each be close to the target
//                    value.

// Scalar-in, scalar-out version !!! needs cleaning
template< typename f, typename T >
const int solve_grid( const f * func, const unsigned int num_in_params,
		const T & init_min_in_params, const T & init_max_in_params,
		const T & init_init_in_params_step, const T & target_out_params,
		T & result_in_params, const double init_init_precision = 0.00001,
		const int search_precision = 0.1,
		const double & init_out_params_weight = 0, const bool silent = false )
{

	BRG_UNITS d = 0, d_best = DBL_MAX;
	int i_resid = 0, i_temp = 0, i_best = -1;
	BRG_UNITS init_in_params_step( 0 );
	BRG_UNITS in_params_step( 0 );
	BRG_UNITS test_in_params( 0 );
	BRG_UNITS best_in_params( 0 );
	BRG_UNITS test_out_params( 0 );
	BRG_UNITS min_in_params = init_min_in_params, max_in_params =
			init_max_in_params;
	double out_params_weight = init_out_params_weight;

	const int default_step_number = 10;
	const double grid_shortening_factor = 0.5;

	double init_precision = init_init_precision, precision = init_precision;
	double step_dist = 1;
	int num_test_points = 0;
	int divisor = 1;
	double total_weight = 0;
	int num_in_params_steps( 0 );

	// Check for sanity and set up parameters
	try
	{
		num_test_points = 1;
		total_weight = 0;
		for ( unsigned int i = 0; i < num_in_params; i++ )
		{
			if ( max_in_params < min_in_params )
				std::swap( max_in_params, min_in_params );
			if ( init_init_in_params_step < 0 )
			{
				init_in_params_step = min( -init_init_in_params_step,
						max_in_params - min_in_params );
			}
			else if ( init_init_in_params_step == 0 )
			{
				init_in_params_step = ( max_in_params - min_in_params )
						/ default_step_number;
			}
			else
			{
				init_in_params_step = min( init_init_in_params_step,
						max_in_params - min_in_params );
			}
			num_in_params_steps = (int)floor(
					( max_in_params - min_in_params ) / init_in_params_step )
					+ 1;
			num_test_points *= num_in_params_steps;
			out_params_weight = std::fabs( out_params_weight );
			total_weight += out_params_weight;
		}
		if ( ( init_init_precision > 0 ) && ( init_init_precision <= 1 ) )
			init_precision = init_init_precision;
		if ( total_weight <= 0 )
		{
			out_params_weight = 1;
			total_weight = num_in_params;
		}
	}
	catch ( std::exception & )
	{
		return errorNOS( silent );
	}

	// First step is to search solution space for the best starting point
	i_best = -1;
	d_best = DBL_MAX;
	in_params_step = init_in_params_step;
	for ( int i = 0; i < num_test_points; i++ )
	{
		// Figure out test_in_params for this point
		i_resid = i;
		divisor = num_test_points;
		for ( unsigned int j = 0; j < num_in_params; j++ )
		{
			divisor /= num_in_params_steps;
			i_temp = (int)( i_resid / divisor );
			i_resid -= i_temp * divisor;
			test_in_params = min_in_params + in_params_step * i_temp;
		}
		if ( int errcode = ( *func )( test_in_params, test_out_params,
				silent ) )
			return errcode + LOWER_LEVEL_ERROR;
		d = std::fabs( test_out_params - target_out_params );
		if ( d < d_best )
		{
			d_best = d;
			i_best = i;
		}

	} // for(int i = 0; i < num_test_points; i++ )

	if ( i_best == -1 ) // We didn't find any suitable starting point
	{
		// Try a finer search to see if that can find it
		int search_factor = (int)( 1. / search_precision );
		in_params_step = init_in_params_step;
		for ( unsigned int i = 0; i < num_in_params; i++ )
			in_params_step /= search_factor;
		i_best = -1;
		d_best = DBL_MAX;

		// Recalculate number of test points
		num_test_points = 1;
		for ( unsigned int i = 0; i < num_in_params; i++ )
			num_test_points *= num_in_params_steps * search_factor;

		for ( int i = 0; i < num_test_points; i++ )
		{
			// Figure out test_in_params for this point
			i_resid = i;
			divisor = num_test_points;
			for ( unsigned int j = 0; j < num_in_params; j++ )
			{
				divisor /= num_in_params_steps;
				i_temp = (int)( i_resid / divisor );
				i_resid -= i_temp * divisor;
				test_in_params = min_in_params + in_params_step * i_temp;
			}
			if ( int errcode = ( *func )( test_in_params, test_out_params,
					silent ) )
				return errcode + LOWER_LEVEL_ERROR;
			d = std::fabs( test_out_params - target_out_params );
			if ( d < d_best )
			{
				d_best = d;
				i_best = i;
			}

		} // for(int i = 0; i < num_test_points; i++ )

	} // if(i_best == -1)

	// Check again to see if we've found a good start point
	if ( i_best == -1 )
	{
		// Nope. The function might just be undefined everywhere. Return an warning.
		if ( !silent )
			std::cerr
					<< "WARNING: Could not solve function with solve_grid - no defined points found.\n";
		return UNSPECIFIED_ERROR;
	} // if(i_best == -1)

	// Get best_in_params
	best_in_params = 0;
	i_resid = i_best;
	divisor = num_test_points;
	for ( unsigned int j = 0; j < num_in_params; j++ )
	{
		divisor /= num_in_params_steps;
		i_temp = (int)( i_resid / divisor );
		i_resid -= i_temp * divisor;
		best_in_params = min_in_params + in_params_step * i_temp;
	}

	// Narrowing search
	step_dist = 1;
	num_test_points = round_int( std::pow( 3, num_in_params ) );
	in_params_step = init_in_params_step;
	precision = init_precision;

	while ( step_dist > precision )
	{
		for ( unsigned int i = 0; i < num_in_params; i++ )
		{
			in_params_step *= grid_shortening_factor;
		}

		i_best = -1;
		d_best = DBL_MAX;
		for ( int i = 0; i < num_test_points; i++ )
		{
			// Figure out test_in_params for this point
			i_resid = i;
			divisor = num_test_points;
			for ( unsigned int j = 0; j < num_in_params; j++ )
			{
				divisor /= 3;
				i_temp = (int)( i_resid / divisor );
				i_resid -= i_temp * divisor;
				i_temp -= 1;
				test_in_params = best_in_params + in_params_step * i_temp;
			}
			if ( int errcode = ( *func )( test_in_params, test_out_params,
					silent ) )
				return errcode + LOWER_LEVEL_ERROR;
			d = std::fabs( test_out_params - target_out_params );
			if ( d < d_best )
			{
				d_best = d;
				i_best = i;
			}

		} // for(int i = 0; i < num_test_points; i++ )

		// Figure out test_in_params for this point
		i_resid = i_best;
		divisor = num_test_points;
		for ( unsigned int j = 0; j < num_in_params; j++ )
		{
			divisor /= 3;
			i_temp = (int)( i_resid / divisor );
			i_resid -= i_temp * divisor;
			i_temp -= 1;
			test_in_params = best_in_params + in_params_step * i_temp;
		}
		best_in_params = test_in_params;

		step_dist *= grid_shortening_factor;

	} // while( step_dist > precision )

	result_in_params = best_in_params;
	return 0;
}

// Vector-in, vector-out version
template< typename f, typename T >
const int solve_grid( const f * func, const unsigned int num_in_params,
		const std::vector< T > & init_min_in_params,
		const std::vector< T > & init_max_in_params,
		const std::vector< T > & init_init_in_params_step,
		const std::vector< T > & target_out_params,
		std::vector< T > & result_in_params, const double init_init_precision =
				0.00001, const int search_precision = 0.1,
		const std::vector< double > & init_out_params_weight = std::vector<
				double >( 0 ), const bool silent = false )
{
	unsigned int num_out_params = 0;

	BRG_UNITS d = 0, d_best = DBL_MAX;
	int i_resid = 0, i_temp = 0, i_best = -1;
	std::vector< BRG_UNITS > init_in_params_step( num_in_params, 0 );
	std::vector< BRG_UNITS > in_params_step( num_in_params, 0 );
	std::vector< BRG_UNITS > test_in_params( num_in_params, 0 );
	std::vector< BRG_UNITS > best_in_params( num_in_params, 0 );
	std::vector< BRG_UNITS > test_out_params( num_out_params, 0 );
	std::vector< BRG_UNITS > min_in_params = init_min_in_params,
			max_in_params = init_max_in_params;
	std::vector< double > out_params_weight = init_out_params_weight;

	const int default_step_number = 10;
	const double grid_shortening_factor = 0.5;

	double init_precision = init_init_precision, precision = init_precision;
	double step_dist = 1;
	int num_test_points = 0;
	int divisor = 1;
	double total_weight = 0;
	std::vector< int > num_in_params_steps( num_in_params, 0 );

	// Check for sanity and set up parameters
	try
	{
		num_test_points = 1;
		if ( out_params_weight.size() == 0 )
		{
			out_params_weight.resize( target_out_params.size(), 1 );
		}
		else if ( out_params_weight.size() != target_out_params.size() )
		{
			out_params_weight.resize( target_out_params.size(), 0 );
		}
		total_weight = 0;
		for ( unsigned int i = 0; i < num_in_params; i++ )
		{
			if ( max_in_params.at( i ) < min_in_params.at( i ) )
				std::swap( max_in_params[i], min_in_params[i] );
			if ( init_init_in_params_step.at( i ) < 0 )
			{
				init_in_params_step.at( i ) = min(
						-init_init_in_params_step[i],
						max_in_params.at( i ) - min_in_params.at( i ) );
			}
			else if ( init_init_in_params_step.at( i ) == 0 )
			{
				init_in_params_step.at( i ) = ( max_in_params.at( i )
						- min_in_params.at( i ) ) / default_step_number;
			}
			else
			{
				init_in_params_step.at( i ) = min( init_init_in_params_step[i],
						max_in_params.at( i ) - min_in_params.at( i ) );
			}
			num_in_params_steps.at( i ) = (int)floor(
					( max_in_params.at( i ) - min_in_params.at( i ) )
							/ init_in_params_step.at( i ) ) + 1;
			num_test_points *= num_in_params_steps.at( i );
		}
		for(unsigned int i = 0; i < out_params_weight.size(); i++)
		{
			out_params_weight.at( i ) = fabs( out_params_weight.at( i ) );
			total_weight += out_params_weight.at( i );
		}
		if ( ( init_init_precision > 0 ) && ( init_init_precision <= 1 ) )
			init_precision = init_init_precision;
		if ( total_weight <= 0 )
		{
			out_params_weight.resize( num_in_params, 1 );
			total_weight = num_in_params;
		}
	}
	catch ( std::exception & )
	{

		return errorNOS( silent );
	}

	// First step is to search solution space for the best starting point
	i_best = -1;
	d_best = DBL_MAX;
	in_params_step = init_in_params_step;
	for ( int i = 0; i < num_test_points; i++ )
	{
		// Figure out test_in_params for this point
		i_resid = i;
		divisor = num_test_points;
		for ( unsigned int j = 0; j < num_in_params; j++ )
		{
			divisor /= num_in_params_steps.at( j );
			i_temp = (int)( i_resid / divisor );
			i_resid -= i_temp * divisor;
			test_in_params.at( j ) = min_in_params.at( j )
					+ in_params_step.at( j ) * i_temp;
		}
		if ( int errcode = ( *func )( test_in_params, test_out_params,
				silent ) )
			return errcode + LOWER_LEVEL_ERROR;
		if ( test_out_params.size() != target_out_params.size() )
			throw std::runtime_error("ERROR: target_out_params passed to solve_grid has incorrect size.");
		if ( test_out_params.size() != out_params_weight.size() )
			throw std::runtime_error("ERROR: out_params_weight passed to solve_grid has incorrect size.");
		d = weighted_dist( test_out_params, target_out_params,
				out_params_weight );
		if ( d < d_best )
		{
			d_best = d;
			i_best = i;
		}

	} // for(int i = 0; i < num_test_points; i++ )

	if ( i_best == -1 ) // We didn't find any suitable starting point
	{
		// Try a finer search to see if that can find it
		int search_factor = (int)( 1. / search_precision );
		in_params_step = init_in_params_step;
		for ( unsigned int i = 0; i < in_params_step.size(); i++ )
			in_params_step[i] /= search_factor;
		i_best = -1;
		d_best = DBL_MAX;

		// Recalculate number of test points
		num_test_points = 1;
		for ( unsigned int i = 0; i < num_in_params; i++ )
			num_test_points *= num_in_params_steps.at( i ) * search_factor;

		for ( int i = 0; i < num_test_points; i++ )
		{
			// Figure out test_in_params for this point
			i_resid = i;
			divisor = num_test_points;
			for ( unsigned int j = 0; j < num_in_params; j++ )
			{
				divisor /= num_in_params_steps.at( j );
				i_temp = (int)( i_resid / divisor );
				i_resid -= i_temp * divisor;
				test_in_params.at( j ) = min_in_params.at( j )
						+ in_params_step.at( j ) * i_temp;
			}
			if ( int errcode = ( *func )( test_in_params, test_out_params,
					silent ) )
				return errcode + LOWER_LEVEL_ERROR;
			if ( test_out_params.size() != target_out_params.size() )
				throw std::runtime_error("ERROR: target_out_params passed to solve_grid has incorrect size.");
			if ( test_out_params.size() != out_params_weight.size() )
				throw std::runtime_error("ERROR: out_params_weight passed to solve_grid has incorrect size.");
			d = weighted_dist( test_out_params, target_out_params,
					out_params_weight );
			if ( d < d_best )
			{
				d_best = d;
				i_best = i;
			}

		} // for(int i = 0; i < num_test_points; i++ )

	} // if(i_best == -1)

	// Check again to see if we've found a good start point
	if ( i_best == -1 )
	{
		// Nope. The function might just be undefined everywhere. Return an warning.
		if ( !silent )
			std::cerr
					<< "WARNING: Could not solve function with solve_grid - no defined points found.\n";
		return UNSPECIFIED_ERROR;
	} // if(i_best == -1)

	// Get best_in_params
	best_in_params.resize( num_in_params, 0 );
	i_resid = i_best;
	divisor = num_test_points;
	for ( unsigned int j = 0; j < num_in_params; j++ )
	{
		divisor /= num_in_params_steps.at( j );
		i_temp = (int)( i_resid / divisor );
		i_resid -= i_temp * divisor;
		best_in_params.at( j ) = min_in_params.at( j )
				+ in_params_step.at( j ) * i_temp;
	}

	// Narrowing search
	step_dist = 1;
	num_test_points = round_int( std::pow( 3, num_in_params ) );
	in_params_step = init_in_params_step;
	precision = init_precision;

	while ( step_dist > precision )
	{
		for ( unsigned int i = 0; i < num_in_params; i++ )
		{
			in_params_step[i] *= grid_shortening_factor;
		}

		i_best = -1;
		d_best = DBL_MAX;
		for ( int i = 0; i < num_test_points; i++ )
		{
			// Figure out test_in_params for this point
			i_resid = i;
			divisor = num_test_points;
			for ( unsigned int j = 0; j < num_in_params; j++ )
			{
				divisor /= 3;
				i_temp = (int)( i_resid / divisor );
				i_resid -= i_temp * divisor;
				i_temp -= 1;
				test_in_params.at( j ) = best_in_params.at( j )
						+ in_params_step.at( j ) * i_temp;
			}
			if ( int errcode = ( *func )( test_in_params, test_out_params,
					silent ) )
				return errcode + LOWER_LEVEL_ERROR;
			if ( test_out_params.size() != target_out_params.size() )
				throw std::runtime_error("ERROR: target_out_params passed to solve_grid has incorrect size.");
			if ( test_out_params.size() != out_params_weight.size() )
				throw std::runtime_error("ERROR: out_params_weight passed to solve_grid has incorrect size.");
			d = weighted_dist( test_out_params, target_out_params,
					out_params_weight );
			if ( d < d_best )
			{
				d_best = d;
				i_best = i;
			}

		} // for(int i = 0; i < num_test_points; i++ )

		// Figure out test_in_params for this point
		i_resid = i_best;
		divisor = num_test_points;
		for ( unsigned int j = 0; j < num_in_params; j++ )
		{
			divisor /= 3;
			i_temp = (int)( i_resid / divisor );
			i_resid -= i_temp * divisor;
			i_temp -= 1;
			test_in_params.at( j ) = best_in_params.at( j )
					+ in_params_step.at( j ) * i_temp;
		}
		best_in_params = test_in_params;

		step_dist *= grid_shortening_factor;

	} // while( step_dist > precision )

	result_in_params = best_in_params;
	return 0;
}

// Scalar-in, scalar-out version
template< typename f, typename T >
const int solve_grid( const f * func, const unsigned int num_in_params,
		const T & init_min_in_params, const T & init_max_in_params,
		const int num_search_steps, const T & target_out_params,
		T & result_in_params, const double init_init_precision = 0.00001,
		const int search_precision = 0.1, const double & out_params_weight = 0,
		const bool silent = false )
{
	T in_params_step( 0 );

	const int steps = (int)max( num_search_steps, 1 );

	in_params_step = ( init_max_in_params - init_min_in_params ) / steps;
	return brgastro::solve_grid( func, num_in_params, init_min_in_params,
			init_max_in_params, in_params_step, target_out_params,
			result_in_params, init_init_precision, search_precision,
			out_params_weight );
} // const int solve_grid(...)

// Vector-in, vector-out version
template< typename f, typename T >
const int solve_grid( const f * func, const unsigned int num_in_params,
		const std::vector< T > & init_min_in_params,
		const std::vector< T > & init_max_in_params,
		const int num_search_steps, const std::vector< T > & target_out_params,
		std::vector< T > & result_in_params, const double init_init_precision =
				0.00001, const int search_precision = 0.1,
		const std::vector< double > & out_params_weight =
				std::vector< double >( 0 ), const bool silent = false )
{
	std::vector< T > in_params_step( num_in_params, 0 );

	const int steps = (int)max( num_search_steps, 1 );

	try
	{
		for ( unsigned int i = 0; i < num_in_params; i++ )
			in_params_step.at( i ) = ( init_max_in_params.at( i )
					- init_min_in_params.at( i ) ) / steps;
	}
	catch ( std::exception & )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Incorrect size for an array passed to brgastro::unit_solve_grid.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	return brgastro::solve_grid( func, num_in_params, init_min_in_params,
			init_max_in_params, in_params_step, target_out_params,
			result_in_params, init_init_precision, search_precision,
			out_params_weight );
} // const int solve_grid(...)

/** Attempts to find the minimum output value for the passed function using a Metropolis-Hastings
 *  MCMC algorithm with simulated annealing.
 *
 * @param func
 * @param init_in_param
 * @param init_min_in_param
 * @param init_max_in_param
 * @param init_in_param_step_sigma
 * @param result_in_param
 * @param result_out_param
 * @param max_steps
 * @param annealing_period
 * @param annealing_factor
 * @param silent
 * @return
 */
template< typename f, typename T >
const int solve_MCMC( const f * func, const T init_in_param, const T init_min_in_param,
		const T init_max_in_param, const T init_in_param_step_sigma,
		T & result_in_param, T & result_out_param, const int max_steps=1000000,
		const int annealing_period=100000, const double annealing_factor=4,
		const bool silent = false)
{
	int step_num = 0;
	bool bounds_check = true;

	// Check how bounds were passed in
	if(init_min_in_param==init_max_in_param) bounds_check = false;

	T min_in_param, max_in_param;
	if(init_min_in_param>init_max_in_param)
	{
		// Swap them
		max_in_param = init_min_in_param;
		min_in_param = init_max_in_param;
	}
	else
	{
		min_in_param = init_min_in_param;
		max_in_param = init_max_in_param;
	}
	T in_param_step_sigma;

	// Check step size
	if(init_in_param_step_sigma <= 0)
	{
		if(bounds_check)
		{
			in_param_step_sigma = (max_in_param-min_in_param)/10.;
		}
		else
		{
			in_param_step_sigma = 1; // Default behavior
		}
	}
	else
	{
		in_param_step_sigma = init_in_param_step_sigma;
	}

	// Initialise
	double annealing = 1;
	bool last_cycle = false;
	T current_in_param = init_in_param, test_in_param = init_in_param, best_in_param = init_in_param;
	T out_param, best_out_param;
	T mean_in_param = 0;
	int last_cycle_count = 0;

	// Get value at initial point
	if(func(test_in_param, out_param, silent))
		throw std::runtime_error("Cannot execute solve_MCMC at initial point.");
	best_out_param = out_param;
	double last_likelihood = std::exp(-annealing*out_param/2);

	for(unsigned int step = 0; step < max_steps; step++)
	{
		// Get a new value
		test_in_param = current_in_param + in_param_step_sigma*brgastro::Gaus_rand()/annealing;

		// Check if it's in bounds
		if(bounds_check)
			test_in_param = brgastro::bound(min_in_param,test_in_param,max_in_param);

		// Find the result for this value
		bool good_result = true;
		try
		{
			if(func(test_in_param, out_param, silent))
				good_result = false;
		}
		catch(std::exception &e)
		{
			good_result = false;
		}
		double new_likelihood = std::exp(-annealing*out_param/2);

		// If it's usable, check if we should step to it
		if(good_result)
		{
			bool step_to_it = false;
			if(new_likelihood>=last_likelihood)
			{
				step_to_it = true;
			}
			else
			{
				if(drand48() < new_likelihood/last_likelihood)
					step_to_it = true;
			}
			if(step_to_it)
			{
				last_likelihood = new_likelihood;
				current_in_param = test_in_param;

				// Check if we have a new best point
				if(out_param < best_out_param)
				{
					best_out_param = out_param;
					best_in_param = current_in_param;
				}
			}
		}

		// If we're on the last cycle, add to the mean
		if(last_cycle)
		{
			mean_in_param += current_in_param;
			last_cycle_count += 1;
		}

		if(brgastro::divisible(annealing_period,step))
		{
			annealing *= annealing_factor;

			// Recalculate likelihood
			last_likelihood = std::exp(-annealing*out_param/2);

			// Check if we're going into the last cycle
			if(max_steps-step<=annealing_period)
				last_cycle = true;
		}
	} // for(unsigned int step = 0; step < max_steps; step++)

	// Calculate mean now
	mean_in_param /= last_cycle_count;

	// Check if mean actually gives a better best
	try
	{
		if(!(func(mean_in_param,out_param,silent)))
		{
			if(out_param < best_out_param)
			{
				best_in_param = mean_in_param;
			}
		}
	}
	catch(std::exception &e)
	{
		// Just leave it, no need to do anything
	}

	result_in_param = best_in_param;

	return 0;
}

} // namespace brgastro

#endif // __BRG_SOLVERS_HPP_INCLUDED__
