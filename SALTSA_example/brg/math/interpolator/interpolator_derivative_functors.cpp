/**********************************************************************\
  @file interpolator_derivative_functors.cpp

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

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <utility>

#include "brg/math/calculus/differentiate.hpp"
#include "brg/math/calculus/integrate.hpp"
#include "brg/math/random/distributions.hpp"

#include "interpolator_derivative_functors.h"

// brgastro::interpolator_functor class method implementations
#if (1)

// Constructors
brgastro::interpolator_functor::interpolator_functor()
{
	_interpolator_ptr_ = NULL;
	_interpolator_ptr_set_up_ = false;
}
brgastro::interpolator_functor::interpolator_functor(
		const brgastro::interpolator *init_interpolator_ptr )
{
	set_interpolator_ptr( init_interpolator_ptr );
}

// Set functions
void brgastro::interpolator_functor::set_interpolator_ptr(
		const brgastro::interpolator *new_spline_ptr )
{
	_interpolator_ptr_ = new_spline_ptr;
	_interpolator_ptr_set_up_ = true;
}

// Function method
BRG_UNITS brgastro::interpolator_functor::operator()( CONST_BRG_UNITS_REF  in_param,
 const bool silent ) const
{
	if ( !_interpolator_ptr_set_up_ )
	{
		throw std::logic_error("Interpolator pointer must be defined in spline_functor.\n");
	}

	return ( *_interpolator_ptr_ )( in_param );
}
#endif

// brgastro::interpolator_derivative_functor class method implementations
#if (1)

// Constructors
brgastro::interpolator_derivative_functor::interpolator_derivative_functor()
{
	_interpolator_functor_set_up_ = false;
}
brgastro::interpolator_derivative_functor::interpolator_derivative_functor(
		brgastro::interpolator *init_spline_ptr )
{
	set_interpolator_ptr( init_spline_ptr );
}

// Set functions
void brgastro::interpolator_derivative_functor::set_interpolator_ptr(
		const brgastro::interpolator *new_interpolator_ptr )
{
	_interpolator_functor_.set_interpolator_ptr( new_interpolator_ptr );
	_interpolator_functor_set_up_ = true;
}

// Function method
BRG_UNITS brgastro::interpolator_derivative_functor::operator()(
		CONST_BRG_UNITS_REF  in_param, const bool silent ) const
{
	if ( !_interpolator_functor_set_up_ )
	{
		throw std::logic_error("Spline function must be set up in spline_derivative_functor.");
	}
	BRG_UNITS Jacobian(0);

	Jacobian = differentiate( &_interpolator_functor_, in_param );

	return Jacobian;
}

#endif

// brgastro::interpolator_derivative_weight_functor class method implementations
#if (1)

// Constructors
brgastro::interpolator_derivative_weight_functor::interpolator_derivative_weight_functor()
{
	_sample_scale_ = 0;
	_sample_max_width_ = 0;
	_t_max_ = -( 0.99 * DBL_MIN );
	_t_min_ = DBL_MAX;
	_centre_point_ = 0;
}

// Set functions
void brgastro::interpolator_derivative_weight_functor::set_sample_scale(
		double new_sample_scale )
{
	_sample_scale_ = new_sample_scale;
}

void brgastro::interpolator_derivative_weight_functor::set_sample_max_width(
		double new_sample_max_width )
{
	_sample_max_width_ = new_sample_max_width;
}

void brgastro::interpolator_derivative_weight_functor::set_center_point(
		double new_center_point )
{
	_centre_point_ = new_center_point;
}

void brgastro::interpolator_derivative_weight_functor::set_t_min(
		double new_t_min )
{
	_t_min_ = new_t_min;
}

void brgastro::interpolator_derivative_weight_functor::set_t_max(
		double new_t_max )
{
	_t_max_ = new_t_max;
}

// Function method
BRG_UNITS brgastro::interpolator_derivative_weight_functor::operator()(
		CONST_BRG_UNITS_REF in_param, const bool silent ) const
{

	// Check bounds
	if ( std::fabs( in_param - _centre_point_ )
			> _sample_max_width_ * std::fabs( _t_max_ - _t_min_ ) )
	{
		return 0;
	}
	else
	{
		return Gaus_pdf( in_param, _centre_point_,
				_sample_scale_ * std::fabs( _t_max_ - _t_min_ ) );
	}
}
#endif


