/**********************************************************************\
 @file differentiate.hpp
 ------------------

 Functions to differentiate a function pointer.

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

#ifndef _BRG_DIFFERENTIATE_HPP_INCLUDED_
#define _BRG_DIFFERENTIATE_HPP_INCLUDED_

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "brg/global.h"

#include "brg/math/misc_math.hpp"
#include "brg/math/safe_math.hpp"
#include "brg/physics/units/unit_obj.h"
#include "brg/utility.hpp"

namespace brgastro {

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

// Scalar-in, scalar-out version
template< typename f, typename T >
inline T differentiate( const f * func, const T & in_param,
		const int order = 1, const double power = 1,
		const bool silent = false )
{

	BRG_UNITS d_in_param( 0 );
	BRG_UNITS base_out_param( 0 );
	BRG_UNITS test_in_param( 0 );
	BRG_UNITS test_out_param( 0 );
	BRG_UNITS small_factor_with_units = SMALL_FACTOR;

	bool power_flag = false;
	bool zero_in_flag = false;

	int order_to_use = max( order, 1 );

	if ( ( order_to_use > 1 ) )
	{
		throw std::logic_error("ERROR: brgastro::differentiate with order > 1 is not currently supported.\n");
	}

	if ( power != 1 )
		power_flag = true;
	else
		power_flag = false;

	// Check if any in_params are zero. If so, estimate small factor from other in_params
	if ( in_param == 0 )
	{
		zero_in_flag = true;
	}
	else     // if(in_params==0)
	{
		small_factor_with_units = in_param * SMALL_FACTOR;
		d_in_param = small_factor_with_units;
	} // else

	if ( zero_in_flag )
	{
		if ( small_factor_with_units == 0 )
		{
#ifdef _BRG_USE_UNITS_
			d_in_param.set(SMALL_FACTOR,in_param.get_unit_powers());
#else
			d_in_param = SMALL_FACTOR;
#endif
		}
		else
		{
			if ( in_param == 0 )
			{
#ifdef _BRG_USE_UNITS_
				d_in_param.set(small_factor_with_units.get_value(),in_param.get_unit_powers());
#else
				d_in_param = small_factor_with_units;
#endif
			} // if(in_params[i]==0)
		}
	}

	// Get value of function at input parameters
	base_out_param = ( *func )( in_param, silent );

	bool bad_function_result = false;
	unsigned int counter = 0;

	T Jacobian=0;

	do {
		counter++;
		bad_function_result = false;

		test_in_param = in_param + d_in_param;

		// Run the function to get value at test point
		try
		{
			test_out_param = ( *func )( test_in_param, silent );
		}
		catch(const std::runtime_error &e)
		{
			bad_function_result = true;
			d_in_param /= 10; // Try again with smaller step
			continue;
		}

		// Record this derivative
		Jacobian = ( test_out_param - base_out_param ) / d_in_param;
		if ( power_flag )
			Jacobian *= power * safe_pow( base_out_param, power - 1 );
		if(isbad(Jacobian))
		{
			bad_function_result = true;
			d_in_param /= 10; // Try again with smaller step
			continue;
		}
	} while ((bad_function_result) && (counter<3));

	if(counter>=3)
		throw std::runtime_error("Cannot differentiate function due to lack of valid nearby points found.");

	return Jacobian;
}

// Vector-in, vector-out version
template< typename f, typename T >
inline std::vector< std::vector< T > > differentiate( const f * func, const std::vector< T > & in_params,
		const int order = 1, const double power = 1, const bool silent = false )
{
	const typename std::vector<T>::size_type num_in_params = in_params.size();
	std::vector< std::vector< T > > Jacobian;

	std::vector< T > d_in_params( 0 );
	std::vector< T > base_out_params( 0 );
	std::vector< T > test_in_params( 0 );
	std::vector< T > test_out_params( 0 );
	T small_factor_with_units = SMALL_FACTOR;

	bool power_flag = false;
	bool zero_in_flag = false;

	int order_to_use = (int)max( order, 1 );

	if ( ( order_to_use > 1 ) )
	{
		throw std::logic_error("brgastro::differentiate with order > 1 is not currently supported.\n");
	}

	if ( power != 1 )
		power_flag = true;
	else
		power_flag = false;

	// Delete std::vectors we'll be overwriting in case they previously existed
	Jacobian.clear();

	// Set up differentials
	make_array( d_in_params, num_in_params );
	make_array( test_in_params, num_in_params );

	// Check if any in_params are zero. If so, estimate small factor from other in_params
	for ( size_t i = 0; i < num_in_params; i++ )
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
	} // for( size_t i = 0; i < num_in_params; i++ )

	if ( zero_in_flag )
	{
		if ( small_factor_with_units == 0 )
		{
			// At least try to get the units right
			for ( size_t i = 0; i < num_in_params; i++ )
			{
#ifdef _BRG_USE_UNITS_
				d_in_params[i].set(SMALL_FACTOR,in_params[i].get_unit_powers());
#else
				d_in_params[i] = SMALL_FACTOR;
#endif
			} // for( size_t i = 0; i < num_in_params; i++ )
		}
		else
		{
			for ( size_t i = 0; i < num_in_params; i++ )
			{
				if ( in_params[i] == 0 )
				{
#ifdef _BRG_USE_UNITS_
					d_in_params[i].set(small_factor_with_units.get_value(),in_params[i].get_unit_powers());
#else
					d_in_params[i] = small_factor_with_units;
#endif
				} // if(in_params[i]==0)
			} // for( size_t i = 0; i < num_in_params; i++ )
		}
	}

	// Get value of function at input parameters
	base_out_params = ( *func )( in_params, silent );
	typename std::vector<T>::size_type num_out_params=base_out_params.size();

	// Set up Jacobian
	make_array2d( Jacobian, num_out_params, num_in_params );

	// Loop over input and output dimensions to get Jacobian

	bool bad_function_result = false;
	unsigned int counter = 0;
	do {
		counter++;
		bad_function_result = false;
		for ( size_t j = 0; j < num_in_params; j++ )
		{
			// Set up test input parameters
			for ( size_t j2 = 0; j2 < num_in_params; j2++ )
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
			try
			{
				test_out_params = ( *func )( test_in_params, silent );
			}
			catch(const std::exception &e)
			{
				bad_function_result = true;
				for(size_t j=0; j< in_params.size(); j++)
					d_in_params[j] /= 10; // Try again with smaller step size
				continue;
			}

			// Record this derivative
			for ( size_t i = 0; i < num_out_params; i++ )
			{
				Jacobian[i][j] = ( test_out_params[i] - base_out_params[i] )
						/ d_in_params[j];
				if ( power_flag )
					Jacobian[i][j] *= power
							* safe_pow( base_out_params[i], power - 1 );
				if(isbad(Jacobian[i][j]))
				{
					bad_function_result = true;
					for(size_t j=0; j< in_params.size(); j++)
						d_in_params[j] /= 10; // Try again with smaller step size
					continue;
				}
			} // for( int i = 0; i < num_out_params; i++)
		} // for( size_t j = 0; j < num_in_params; j++)
	} while (bad_function_result && (counter<3));

	if(counter>=3)
		throw std::runtime_error("Cannot differentiate function due to lack of valid nearby points found.");

	return Jacobian;
}

}



#endif // _BRG_DIFFERENTIATE_HPP_INCLUDED_
