/**********************************************************************\
  @file gabdt.h

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

// body file: brg/physics/astro/SALTSA/gabdt.cpp

#ifndef _GABDT_H_INCLUDED_
#define _GABDT_H_INCLUDED_

#include <vector>

#include "brg/global.h"

#include "brg/physics/density_profile/density_profile.h"
#include "brg/physics/units/unit_obj.h"

namespace brgastro {

class gabdt
{
	/************************************************************
	 gabdt
	 -----

	 This class represents a useful physical construct, of how
	 much particles in halos have been disrupted by tidal shocking.

	 See the g_ab object from equation 10 of Taylor and Babul
	 (2001). This represents that multiplied by the timestep.

	 \************************************************************/

private:
	mutable bool _is_cached_;
	const density_profile *_host_ptr_;

	BRG_DISTANCE _x_, _y_, _z_, _r_;
	BRG_TIME _dt_;
	mutable std::vector< std::vector< long double > > _dv_;

public:

	// Constructors
	gabdt();
	gabdt( const density_profile *init_host, CONST_BRG_DISTANCE_REF init_x,
			CONST_BRG_DISTANCE_REF init_y, CONST_BRG_DISTANCE_REF init_z,
			CONST_BRG_TIME_REF init_dt );

	// Destructor
	virtual ~gabdt()
	{
	}

	// Full clear function
	void clear();

	// Set functions
	void set( const brgastro::density_profile *new_host_ptr,
			CONST_BRG_DISTANCE_REF new_x, CONST_BRG_DISTANCE_REF new_y,
			CONST_BRG_DISTANCE_REF new_z, CONST_BRG_TIME_REF new_dt );
	void set_pos( CONST_BRG_DISTANCE_REF new_x, CONST_BRG_DISTANCE_REF new_y,
			CONST_BRG_DISTANCE_REF new_z );
	void set_dt( CONST_BRG_TIME_REF dt );
	void set_host_ptr( const density_profile *new_host_ptr );
	void override_zero();

	// Calculation function
	void calc_dv( const bool silent = false ) const;

	// Get functions
	const density_profile * host() const;
	BRG_DISTANCE x() const;
	BRG_DISTANCE y() const;
	BRG_DISTANCE z() const;
	BRG_DISTANCE r() const;
	std::vector< std::vector< long double > > dv() const; // involves calculation if necessary
	long double dv( const int x_i, const int y_i ) const; // involves calculation if necessary

	// Operator overloading
	BRG_UNITS operator*( const gabdt & other_gabdt ) const; // Dot-product(ish) operator

	gabdt & operator+=( const gabdt & other_gabdt ); // Addition
	gabdt operator+( const gabdt & other_gabdt ) const;

	gabdt & operator*=( const double scale_fraction ); // Multiplication by a double
	gabdt operator*( const double scale_fraction ) const;

};

class gabdt_functor
{
	/************************************************************
	 gabdt_functor
	 --------------

	 Child of functor

	 This class provides a functor * for getting the 3-D
	 acceleration within a halo.

	 \************************************************************/
public:

	// Constructors
	gabdt_functor();

	// Destructor
	virtual ~gabdt_functor()
	{
	}

	// Host accessor
	const density_profile *host_ptr;

	// Function method
	std::vector< BRG_UNITS > operator()( const std::vector< BRG_UNITS > & in_params,
			const bool silent = false ) const;

};

} // end namespace brgastro

#endif /* _GABDT_H_INCLUDED_ */
