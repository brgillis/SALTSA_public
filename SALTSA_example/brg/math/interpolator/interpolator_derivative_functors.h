/**********************************************************************\
  @file interpolator_derivative_functors.h

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

// body file: brg/math/interpolator/interpolator_derivative_functors.cpp

#ifndef _INTERPOLATOR_DERIVATIVE_FUNCTORS_H_INCLUDED_
#define _INTERPOLATOR_DERIVATIVE_FUNCTORS_H_INCLUDED_

#include "brg/global.h"

#include "brg/math/interpolator/interpolator.h"
#include "brg/physics/units/unit_obj.h"

namespace brgastro {

class interpolator_functor
{
	/************************************************************
	 interpolator_functor
	 ---------------

	 This class is used to provide a functor for getting
	 points along a spline (which is used here to differentiate the
	 spline).

	 Use of this class is handled by the spline_derivative class.
	 No need for the end-user to worry too much about how this
	 works.

	 \************************************************************/

private:

	const brgastro::interpolator *_interpolator_ptr_;
	bool _interpolator_ptr_set_up_;

public:

	// Constructors
	interpolator_functor();
	interpolator_functor(const brgastro::interpolator *init_interpolator_ptr );

	// Destructor
	virtual ~interpolator_functor()
	{
	}

	// Operator=
	interpolator_functor& operator=(interpolator_functor other);

	// Set functions
	void set_interpolator_ptr( const brgastro::interpolator *new_interpolator_ptr );

	// Function method
	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param, const bool silent = false ) const;

};
// class interpolator_functor

class interpolator_derivative_functor
{
	/************************************************************
	 interpolator_derivative_functor
	 --------------------------

	 This class is used to provide a functor * for getting
	 the derivative of an interpolator at a given point.

	 Use of this class is handled by the interpolator_derivative class.
	 No need for the end-user to worry too much about how this
	 works.

	 \************************************************************/

private:

	interpolator_functor _interpolator_functor_;
	bool _interpolator_functor_set_up_;

public:

	// Constructors
	interpolator_derivative_functor();
	interpolator_derivative_functor( brgastro::interpolator *init_interpolator_ptr );

	// Destructor
	virtual ~interpolator_derivative_functor()
	{
	}

	// Operator=
	interpolator_derivative_functor& operator=(interpolator_derivative_functor other);

	// Set functions
	void set_interpolator_ptr( const brgastro::interpolator *new_interpolator_ptr );

	// Function method
	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param, const bool silent = false ) const;

};
// class interpolator_derivative_functor

class interpolator_derivative_weight_functor
{
	/************************************************************
	 interpolator_derivative_weight_functor
	 ---------------------------------

	 Child of functor

	 This class is used to provide a functor * for getting
	 the weight of various points in the smoothing kernel for
	 calculating the derivative of an interpolator.

	 Use of this class is handled by the interpolator_derivative class.
	 No need for the end-user to worry too much about how this
	 works.

	 \************************************************************/

private:

	double _sample_scale_, _sample_max_width_;
	double _t_min_, _t_max_, _centre_point_;

public:

	// Constructors
	interpolator_derivative_weight_functor();

	// Destructor
	virtual ~interpolator_derivative_weight_functor()
	{
	}

	// Operator=
	interpolator_derivative_weight_functor & operator=(interpolator_derivative_weight_functor other);

	// Set functions
	void set_sample_scale( double new_sample_scale );

	void set_sample_max_width( double new_sample_max_width );

	void set_center_point( double new_center_point );

	void set_t_min( double new_t_min );

	void set_t_max( double new_t_max );

	// Function method
	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param, const bool silent = false ) const;
};
// class interpolator_derivative_weight_functor

} // end namespace brgastro

#endif /* _INTERPOLATOR_DERIVATIVE_FUNCTORS_H_INCLUDED_ */
