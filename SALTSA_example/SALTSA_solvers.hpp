/**********************************************************************\
  @file SALTSA_solvers.hpp

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

#ifndef __SALTSA_SOLVERS_HPP_INCLUDED__
#define __SALTSA_SOLVERS_HPP_INCLUDED__

#include <vector>

#include "SALTSA_global.h"

#include "SALTSA_misc_functions.hpp"
#include "SALTSA_vector.hpp"

namespace SALTSA
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

	double new_value = init_param;
	double result = 0;
	double mean_value = 0;
	std::vector< double > past_values( 0 );
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
					<< "WARNING: solve_SALTSA::unit_iterate did not converge.\n";
		return 0;
	}
	result = new_value;
	return result;

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
		const double  init_out_params_weight = 0, const bool silent = false )
{

	double d = 0, d_best = DBL_MAX;
	int i_resid = 0, i_temp = 0, i_best = -1;
	double init_in_params_step( 0 );
	double in_params_step( 0 );
	double test_in_params( 0 );
	double best_in_params( 0 );
	double test_out_params( 0 );
	double min_in_params = init_min_in_params, max_in_params =
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
	num_test_points = ipow( 3, num_in_params );
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

	double d = 0, d_best = DBL_MAX;
	int i_resid = 0, i_temp = 0, i_best = -1;
	std::vector< double > init_in_params_step( num_in_params, 0 );
	std::vector< double > in_params_step( num_in_params, 0 );
	std::vector< double > test_in_params( num_in_params, 0 );
	std::vector< double > best_in_params( num_in_params, 0 );
	std::vector< double > test_out_params( num_out_params, 0 );
	std::vector< double > min_in_params = init_min_in_params,
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
	num_test_points = ipow( 3, num_in_params );
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
		const int search_precision = 0.1, const double  out_params_weight = 0,
		const bool silent = false )
{
	T in_params_step( 0 );

	const int steps = (int)max( num_search_steps, 1 );

	in_params_step = ( init_max_in_params - init_min_in_params ) / steps;
	return SALTSA::solve_grid( func, num_in_params, init_min_in_params,
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
					<< "ERROR: Incorrect size for an array passed to SALTSA::unit_solve_grid.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	return SALTSA::solve_grid( func, num_in_params, init_min_in_params,
			init_max_in_params, in_params_step, target_out_params,
			result_in_params, init_init_precision, search_precision,
			out_params_weight );
} // const int solve_grid(...)

} // namespace SALTSA

#endif // __SALTSA_SOLVERS_HPP_INCLUDED__
