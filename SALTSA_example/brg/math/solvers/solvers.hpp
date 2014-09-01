/**********************************************************************\
  @file solvers.hpp

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

#ifndef _BRG_SOLVERS_HPP_INCLUDED_
#define _BRG_SOLVERS_HPP_INCLUDED_

#include "brg/global.h"

#include <cstdlib>
#include <cmath>
#include "brg/utility.hpp"
#include "brg/math/random/random_functions.h"
#include "brg/physics/units/unit_obj.h"
#include "brg/vector/elementwise_functions.hpp"
#include "brg/vector/summary_functions.hpp"

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

	T new_value = init_param;
	T mean_value = 0;
	std::vector< T > past_values( 0 );
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
		new_value = ( *func )( new_value, silent );

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
					<< "WARNING: solve_iterate did not converge.\n";
	}
	return new_value;
}

// solve_sd: Steepest-descent solver. Finds the minimum value of a function. May fail if it reaches input values which give undefined
//           output. Only valid for functions which return only one output value. Returns the best set of input parameters in the
//           result_in_params vector.
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
template< typename f, typename T >
T solve_sd( const f * func, const T & init_in_params,
		const double precision = 0.00001, const double lambda = 0.1,
		const double cusp_override_power = 0, const int max_steps = 10000,
		const bool silent = false )
{

	T Jacobian( 0 );
	T current_in_params( 0 );
	T last_in_params( 0 );
	T out_params( 0 );
	T in_params_to_use = init_in_params;
	T lambda_norm( 0 );
	T mag;

	const int lambda_shortening_intervals = 10;
	const double lambda_shortening_factor = 10;
	unsigned int num_out_params = 0;
	int num_steps_taken = 0;
	bool converged_flag = false;

	// Sanity checks on variables
	const double lambda_to_use = max( SMALL_FACTOR,
			min( 1 / SMALL_FACTOR, lambda ) );
	const int max_steps_to_use = (int)max( 1, max_steps );

	// Get normalized step size
	lambda_norm = lambda_to_use / safe_d(differentiate( func, in_params_to_use, 1, cusp_override_power ));

	if ( in_params_to_use != 0 )
		lambda_norm *= in_params_to_use;

	// Initialize solver
	num_steps_taken = 0;
	current_in_params = in_params_to_use;

	// Solver

	while ( num_steps_taken < max_steps_to_use )
	{
		num_steps_taken++;
		last_in_params = current_in_params;

		// Get derivative at current location
		Jacobian = differentiate( func,	current_in_params, 1, cusp_override_power );

		// Make a step
		current_in_params -= lambda_norm * Jacobian;

		// Check if we've converged
		converged_flag = true;

		if ( isbad( current_in_params ) )
		{
			throw std::runtime_error("Somehow got bad value for in_params in solve_sd.");
		}
		if ( std::fabs(
				2 * ( current_in_params - last_in_params )
						/ safe_d( current_in_params + last_in_params ) )
				> precision )
		{
			converged_flag = false;
			break;
		}

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
	}
	return current_in_params;

}

// Vector-in, vector-out version
template< typename f, typename T >
std::vector< T >  solve_sd( const f * func,	const std::vector< T > & init_in_params,
		const double precision = 0.00001,
		const double lambda = 0.1, const double cusp_override_power = 0,
		const int max_steps = 10000, const bool silent = false )
{
	std::vector< std::vector< T > > Jacobian( 0 );
	std::vector< T > current_in_params( 0 );
	std::vector< T > last_in_params( 0 );
	std::vector< T > out_params( 0 );
	std::vector< T > in_params_to_use = init_in_params;
	T lambda_norm( 0 );
	T mag;
	typedef typename std::vector<T>::size_type vsize_t;
	vsize_t num_in_params = init_in_params.size();

	const int lambda_shortening_intervals = 10;
	const double lambda_shortening_factor = 10;
	vsize_t num_out_params = 0;
	int num_steps_taken = 0;
	bool converged_flag = false;

	// Sanity checks on variables
	assert(num_in_params >= 1);
	const double lambda_to_use = max( SMALL_FACTOR,
			min( 1 / SMALL_FACTOR, lambda ) );
	const int max_steps_to_use = (int)max( 1, max_steps );

	// Get normalized step size
	Jacobian = differentiate( func, in_params_to_use, 1, cusp_override_power );
	mag = 0;
	for ( vsize_t i = 0; i < num_in_params; i++ )
	{
		mag = quad_add( mag, Jacobian[0][i] );
	}
	lambda_norm = lambda_to_use / mag;
	mag = 0;
	for ( vsize_t i = 0; i < num_in_params; i++ )
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
		Jacobian = differentiate( func, current_in_params, 1, cusp_override_power );

		// Make a step
		for ( vsize_t i = 0; i < num_in_params; i++ )
		{
			current_in_params[i] -= lambda_norm * Jacobian[0][i];
		}

		// Check if we've converged
		converged_flag = true;

		for ( vsize_t i = 0; i < num_in_params; i++ )
		{
			if ( isbad( current_in_params[i] ) )
			{
				throw std::runtime_error("ERROR: Somehow got NaN for in_params in solve_sd.");
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
	}
	return current_in_params;

}

// solve_grid: Custom solver which searches a grid of input parameters for the ones which give output closest to the target. The best point is then
//             used as the starting point for a narrowing search on progressively finer grids of 3^n points (n = number of input parameters).
//             This solver has the advantage that it can better handle regions of input space which give undefined results, and it is less likely to
//             get stuck in local minima. It is, however, slower than other methods, especially for high numbers of input parameters.
//             This function is set up to handle multiple output parameters. If it is more important that one be close to the target than another,
//             use the out_params_weight vector to assign it a higher weight.
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

// Scalar-in, scalar-out version
template< typename f, typename T >
T solve_grid( const f * func, const T & init_min_in_params, const T & init_max_in_params,
		const T & init_init_in_params_step, const T & target_out_params,
		const double init_init_precision = 0.00001, const int search_precision = 0.1,
		const bool silent = false )
{

	T d = 0, d_best = DBL_MAX;
	unsigned int i_best = -1;
	T init_in_params_step( 0 );
	T in_params_step( 0 );
	T test_in_params( 0 );
	T best_in_params( 0 );
	T test_out_params( 0 );
	T min_in_params = init_min_in_params, max_in_params =
			init_max_in_params;

	const int default_step_number = 10;
	const double grid_shortening_factor = 0.5;

	double init_precision = init_init_precision, precision = init_precision;
	double step_dist = 1;
	unsigned int num_test_points = 0;
	unsigned int num_in_params_steps( 0 );

	// Check for sanity and set up parameters
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
	num_test_points = num_in_params_steps;

	// First step is to search solution space for the best starting point
	i_best = 0;
	d_best = DBL_MAX;
	in_params_step = init_in_params_step;
	bool starting_point_found = false;
	for ( unsigned int i = 0; i < num_test_points; i++ )
	{
		test_in_params = min_in_params + in_params_step * i;

		try
		{
			test_out_params = ( *func )( test_in_params, silent );
			d = std::fabs( test_out_params - target_out_params );
			if ( d < d_best )
			{
				d_best = d;
				i_best = i;
				starting_point_found = true;
			}
		}
		catch (const std::exception &e)
		{
		}

	} // for(int i = 0; i < num_test_points; i++ )

	if ( !starting_point_found )
	{
		// Try a finer search to see if that can find it
		int search_factor = (int)( 1. / search_precision );
		in_params_step *= search_precision;

		// Recalculate number of test points
		num_test_points = num_in_params_steps * search_factor;

		for ( unsigned int i = 0; i < num_test_points; i++ )
		{
			// Figure out test_in_params for this point and get value there
			test_in_params = min_in_params + in_params_step * i;
			try
			{
				test_out_params = ( *func )( test_in_params, silent );
				d = std::fabs( test_out_params - target_out_params );
				if ( d < d_best )
				{
					d_best = d;
					i_best = i;
				}
			}
			catch (const std::exception &e)
			{
			}

		} // for(int i = 0; i < num_test_points; i++ )

	} // if(!starting_point_found)

	// Check again to see if we've found a good start point
	if ( !starting_point_found )
	{
		// Nope. The function might just be undefined everywhere. Throw an exception
		throw std::runtime_error("Solve_grid could not find any defined range points.");
	} // if(!starting_point_found)

	// Get best_in_params
	best_in_params = min_in_params + in_params_step * i_best;

	// Narrowing search
	step_dist = 1;
	num_test_points = 3;
	in_params_step = init_in_params_step;
	precision = init_precision;

	while ( step_dist > precision )
	{
		in_params_step *= grid_shortening_factor;

		i_best = 0;
		d_best = DBL_MAX;
		for ( unsigned int i = 0; i < num_test_points; i++ )
		{
			test_in_params = best_in_params + in_params_step * i - in_params_step;
			try
			{
				test_out_params = ( *func )( test_in_params, silent );
				d = std::fabs( test_out_params - target_out_params );
				if ( d < d_best )
				{
					d_best = d;
					i_best = i;
				}
			}
			catch (const std::exception &e)
			{
			}

		} // for(int i = 0; i < num_test_points; i++ )

		// Figure out best_in_params for next step
		best_in_params += in_params_step * i_best - in_params_step;

		step_dist *= grid_shortening_factor;

	} // while( step_dist > precision )

	return best_in_params;
}

// Vector-in, vector-out version
template< typename f, typename T >
std::vector< T > solve_grid( const f * func,
		const std::vector< T > & init_min_in_params,
		const std::vector< T > & init_max_in_params,
		const std::vector< T > & init_init_in_params_step,
		const std::vector< T > & target_out_params,
		const double init_init_precision = 0.00001, const int search_precision = 0.1,
		const std::vector< double > & init_out_params_weight = std::vector<
				double >( 0 ), const bool silent = false )
{
	typedef typename std::vector<T>::size_type vsize_t;
	vsize_t num_in_params = init_min_in_params.size();
	vsize_t num_out_params = 0;

	T d = 0, d_best = DBL_MAX;
	int i_resid = 0, i_temp = 0, i_best = 0;
	std::vector< T > init_in_params_step( num_in_params, 0 );
	std::vector< T > in_params_step( num_in_params, 0 );
	std::vector< T > test_in_params( num_in_params, 0 );
	std::vector< T > best_in_params( num_in_params, 0 );
	std::vector< T > test_out_params( num_out_params, 0 );
	std::vector< T > min_in_params = init_min_in_params,
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

	// First step is to search solution space for the best starting point
	i_best = 0;
	d_best = DBL_MAX;
	in_params_step = init_in_params_step;
	bool starting_point_found = false;
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
		try
		{
			test_out_params = ( *func )( test_in_params, silent );
			assert(test_out_params.size() == target_out_params.size());
			assert(test_out_params.size() == out_params_weight.size());
			d = weighted_dist( test_out_params, target_out_params,
					out_params_weight );
			if ( d < d_best )
			{
				starting_point_found = true;
				d_best = d;
				i_best = i;
			}
		}
		catch (const std::exception &e)
		{
		}

	} // for(int i = 0; i < num_test_points; i++ )

	if ( !starting_point_found ) // We didn't find any suitable starting point
	{
		// Try a finer search to see if that can find it
		int search_factor = (int)( 1. / search_precision );
		in_params_step = init_in_params_step;
		for ( vsize_t i = 0; i < in_params_step.size(); i++ )
			in_params_step[i] /= search_factor;
		i_best = 0;
		d_best = DBL_MAX;

		// Recalculate number of test points
		num_test_points = 1;
		for ( vsize_t i = 0; i < num_in_params; i++ )
			num_test_points *= num_in_params_steps.at( i ) * search_factor;

		for ( int i = 0; i < num_test_points; i++ )
		{
			// Figure out test_in_params for this point
			i_resid = i;
			divisor = num_test_points;
			for ( vsize_t j = 0; j < num_in_params; j++ )
			{
				divisor /= num_in_params_steps.at( j );
				i_temp = (int)( i_resid / divisor );
				i_resid -= i_temp * divisor;
				test_in_params.at( j ) = min_in_params.at( j )
						+ in_params_step.at( j ) * i_temp;
			}
			try
			{
				test_out_params = ( *func )( test_in_params, silent );

				assert(test_out_params.size() == target_out_params.size());
				assert(test_out_params.size() == out_params_weight.size());

				d = weighted_dist( test_out_params, target_out_params,
						out_params_weight );
				if ( d < d_best )
				{
					starting_point_found = true;
					d_best = d;
					i_best = i;
				}
			}
			catch (const std::exception &e)
			{
			}

		} // for(int i = 0; i < num_test_points; i++ )

	} // if(!starting_point_found)

	// Check again to see if we've found a good start point
	if ( !starting_point_found )
	{
		// Nope. The function might just be undefined everywhere. Throw an exception
		throw std::runtime_error("Solve_grid could not find any defined range points.");
	} // if(!starting_point_found)

	// Get best_in_params
	best_in_params.resize( num_in_params, 0 );
	i_resid = i_best;
	divisor = num_test_points;
	for ( vsize_t j = 0; j < num_in_params; j++ )
	{
		divisor /= num_in_params_steps.at( j );
		i_temp = (unsigned int)( i_resid / divisor );
		i_resid -= i_temp * divisor;
		best_in_params.at( j ) = min_in_params.at( j )
				+ in_params_step.at( j ) * i_temp;
	}

	// Narrowing search
	step_dist = 1;
	num_test_points = ipow( 3, num_in_params );
	in_params_step = init_in_params_step;
	precision = init_precision;

	while ( step_dist > precision )
	{
		for ( vsize_t i = 0; i < num_in_params; i++ )
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
			for ( vsize_t j = 0; j < num_in_params; j++ )
			{
				divisor /= 3;
				i_temp = (int)( i_resid / divisor );
				i_resid -= i_temp * divisor;
				i_temp -= 1;
				test_in_params.at( j ) = best_in_params.at( j )
						+ in_params_step.at( j ) * i_temp;
			}
			try
			{
				test_out_params = ( *func )( test_in_params, silent );

				assert(test_out_params.size() == target_out_params.size());
				assert(test_out_params.size() == out_params_weight.size());

				d = weighted_dist( test_out_params, target_out_params,
						out_params_weight );
				if ( d < d_best )
				{
					d_best = d;
					i_best = i;
				}
			}
			catch (const std::exception &e)
			{
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

	return best_in_params;
}

// Scalar-in, scalar-out version
template< typename f, typename T >
T solve_grid( const f * func, const T & init_min_in_params, const T & init_max_in_params,
		const unsigned int num_search_steps, const T & target_out_params,
		const double init_init_precision = 0.00001,
		const int search_precision = 0.1, const bool silent = false )
{
	const unsigned int steps = (unsigned int)max( num_search_steps, 1 );

	T in_params_step(( init_max_in_params - init_min_in_params ) / steps);
	return brgastro::solve_grid( func, init_min_in_params,
			init_max_in_params, in_params_step, target_out_params,
			init_init_precision, search_precision,
			silent );
} // const int solve_grid(...)

// Vector-in, vector-out version
template< typename f, typename T >
std::vector<T> solve_grid( const f * func,
		const std::vector< T > & init_min_in_params,
		const std::vector< T > & init_max_in_params,
		const unsigned int num_search_steps, const std::vector< T > & target_out_params,
		const double init_init_precision = 0.00001, const int search_precision = 0.1,
		const std::vector< double > & out_params_weight =
				std::vector< double >( 0 ), const bool silent = false )
{
	std::vector< T > in_params_step( init_min_in_params.size(), 0 );

	const unsigned int steps = (unsigned int)max( num_search_steps, 1 );

	assert(init_min_in_params.size()==init_max_in_params.size());

	for ( size_t i = 0; i < init_min_in_params.size(); i++ )
		in_params_step[i] = ( init_max_in_params.at[i]
				- init_min_in_params.at[i] ) / steps;

	return brgastro::solve_grid( func, init_min_in_params,
			init_max_in_params, in_params_step, target_out_params,
			init_init_precision, search_precision,
			out_params_weight, silent );
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
T solve_MCMC( const f * func, const T init_in_param, const T init_min_in_param,
		const T init_max_in_param, const T init_in_param_step_sigma,
		const int max_steps=1000000, const int annealing_period=100000,
		const double annealing_factor=4, const bool silent = false)
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
			in_param_step_sigma = 1; // Default behaviour
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
	try
	{
		out_param = func(test_in_param, silent);
	}
	catch(const std::exception &e)
	{
		throw std::runtime_error("Cannot execute solve_MCMC at initial point.");
	}
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
			out_param = func(test_in_param, silent);
		}
		catch(const std::exception &e)
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
				if(drand() < new_likelihood/last_likelihood)
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
			last_likelihood = std::pow(last_likelihood,annealing_factor);

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
		out_param = func(mean_in_param,silent);
	}
	catch(const std::exception &e)
	{
		// Just leave it, no need to do anything
	}

	return best_in_param;
}

/** Attempts to find the minimum output value for the passed function using a Metropolis-Hastings
 *  MCMC algorithm with simulated annealing. Vector version
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
std::vector<T> solve_MCMC( const f * func, const std::vector<T> & init_in_params,
		const std::vector<T> & init_min_in_params,
		const std::vector<T> & init_max_in_params,
		const std::vector<T> & init_in_param_step_sigmas,
		const int max_steps=1000000,
		const int annealing_period=100000, const double annealing_factor=4,
		const bool silent = false)
{
	bool bounds_check = true;

	// Check how bounds were passed in
	if((init_min_in_params.size()==0) and (init_max_in_params.size()==0))
		bounds_check = false;

	std::vector<T> min_in_params = init_min_in_params;
	std::vector<T> max_in_params = init_max_in_params;

	if(min_in_params.size() < max_in_params.size())
	{
		min_in_params.resize(max_in_params.size(),-DBL_MAX);
	}
	if(max_in_params.size() < min_in_params.size())
	{
		max_in_params.resize(min_in_params.size(),DBL_MAX);
	}

	// Check if any of the params likely have max and min mixed up
	if(bounds_check)
	{
		for(unsigned int i = 0; i < max_in_params.size(); i++)
		{
			if(min_in_params.at(i)>max_in_params.at(i))
			{
				std::swap(min_in_params.at(i),max_in_params.at(i));
			} // if(min_in_params.at(i)>max_in_params.at(i))
		} // for(unsigned int i = 0; i < max_in_params.size(); i++)
	} // if(bounds_check)

	// Check step sizes

	std::vector<T> in_param_step_sigmas = init_in_param_step_sigmas;

	if(in_param_step_sigmas.size() < min_in_params.size())
		in_param_step_sigmas.resize(min_in_params.size(),0);

	for(unsigned int i = 0; i < in_param_step_sigmas.size(); i++)
	{
		if(in_param_step_sigmas.at(i) <= 0)
		{
			if(bounds_check)
			{
				in_param_step_sigmas.at(i) =
						(max_in_params.at(i)-min_in_params.at(i))/10.;
			}
			else
			{
				in_param_step_sigmas.at(i) = 1; // Default behaviour
			}
		}
	}

	// Initialise
	double annealing = 1;
	bool last_cycle = false;
	std::vector<T> current_in_params = init_in_params, test_in_params = init_in_params,
			best_in_params = init_in_params;
	std::vector<T> out_params, best_out_params;
	std::vector<T> mean_in_params(min_in_params.size(),0);
	int last_cycle_count = 0;

	// Get value at initial point
	try
	{
		out_params = ( *func )(test_in_params, silent);
	}
	catch(const std::exception &e)
	{
		throw std::runtime_error("Cannot execute solve_MCMC at initial point.");
	}
	best_out_params = out_params;
	double last_log_likelihood = brgastro::sum(
			brgastro::multiply(-annealing/2,out_params));

	const double (*Gaus_rand)(double,double) = brgastro::Gaus_rand;

	for(int step = 0; step < max_steps; step++)
	{
		// Get a new value
		test_in_params = brgastro::rand_vector(Gaus_rand,
				                               current_in_params,
				                               brgastro::divide(in_param_step_sigmas,annealing));

		// Check if it's in bounds
		if(bounds_check)
		{
			test_in_params = bound(min_in_params,test_in_params,max_in_params);
		}

		// Find the result for this value
		bool good_result = true;
		try
		{
			out_params = ( *func )(test_in_params, silent);
		}
		catch(const std::exception &e)
		{
			good_result = false;
		}
		double new_log_likelihood = brgastro::sum(
				brgastro::multiply(-annealing/2,out_params));

		// If it's usable, check if we should step to it
		bool step_to_it = false;
		if(good_result)
		{
			if(new_log_likelihood>=last_log_likelihood)
			{
				step_to_it = true;
			}
			else
			{
				if(drand() < std::exp(new_log_likelihood - last_log_likelihood))
					step_to_it = true;
			}
			if(step_to_it)
			{
				last_log_likelihood = new_log_likelihood;
				current_in_params = test_in_params;

				// Check if we have a new best point
				if(brgastro::sum(out_params) < brgastro::sum(best_out_params))
				{
					best_out_params = out_params;
					best_in_params = current_in_params;
				}
			}
		}

		// If we're on the last cycle, add to the mean
		if(last_cycle)
		{
			mean_in_params = brgastro::add(mean_in_params,current_in_params);
			last_cycle_count += 1;
		}

		if((step!=0)&&(brgastro::divisible(step,annealing_period)))
		{
			annealing *= annealing_factor;

			// Recalculate likelihood
			last_log_likelihood *= annealing_factor;

			// Check if we're going into the last cycle
			if(max_steps-step<=annealing_period)
				last_cycle = true;
		}
	} // for(unsigned int step = 0; step < max_steps; step++)

	// Calculate mean now
	mean_in_params = brgastro::divide(mean_in_params,brgastro::safe_d(last_cycle_count));

	// Check if mean actually gives a better best
	try
	{
		out_params = ( *func )(mean_in_params,silent);
		if(brgastro::sum(out_params) < brgastro::sum(best_out_params))
		{
			best_in_params = mean_in_params;
			best_out_params = out_params;
		}
	}
	catch(const std::exception &e)
	{
		// Just leave it, no need to do anything
	}

	return best_in_params;
}

} // namespace brgastro

#endif // __BRG_SOLVERS_HPP_INCLUDED__
