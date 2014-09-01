/**********************************************************************\
 @file DE.hpp
 ------------------

 Functions to be used with differential equations (DEs).

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

// body file: DE.cpp

#ifndef _BRG_DE_HPP_INCLUDED_
#define _BRG_DE_HPP_INCLUDED_

#include "brg/global.h"

#include "brg/math/misc_math.hpp"
#include "brg/physics/phase.hpp"
#include "brg/physics/units/unit_obj.h"

namespace brgastro {

// Leapfrog method for solving a DE. Note that this implementation assumes that the positions and velocities passed to it are already spaced
// out by half a timestep, with velocity at t+t_step/2 (though it does allow phase classes to be passed to it). This method takes a single step,
// using the passed acceleration function. The passed function for this implementation must take in one parameter (the magnitude of distance from
// a centre point) and return one parameter (the magnitude of the acceleration toward this centre point).
template< typename f >
inline const int leapfrog_step( CONST_BRG_DISTANCE_REF x, CONST_BRG_DISTANCE_REF y,
		CONST_BRG_DISTANCE_REF z, CONST_BRG_VELOCITY_REF vx, CONST_BRG_VELOCITY_REF vy,
		CONST_BRG_VELOCITY_REF vz,
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
	a = (*accel_func)( d, silent );

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
		CONST_BRG_TIME_REF t_step, const f *accel_func,
		const bool silent = false )
{
	return leapfrog_step( p.x, p.y, p.z, p.vx, p.vy, p.vz, new_p.x, new_p.y,
			new_p.z, new_p.vx, new_p.vy, new_p.vz, t_step, accel_func, silent );
}

template< typename f >
inline const int leapfrog_step( phase & p, CONST_BRG_TIME_REF t_step,
		const f *accel_func, const bool silent = false )
{
	int result;
	phase new_p(p);
	result = leapfrog_step( p, new_p, t_step, accel_func, silent );
	p = new_p;
	return result;
}

}

#endif // _BRG_DE_HPP_INCLUDED_
