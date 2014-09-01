/**********************************************************************\
  @file phase.hpp

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

#ifndef _BRG_PHASE_HPP_INCLUDED_
#define _BRG_PHASE_HPP_INCLUDED_

#include "brg/global.h"

#ifdef _BRG_USE_UNITS_
#include "brg_units.h"
#endif

namespace brgastro
{

/** Class definitions **/
#if (1)
struct phase
/**********************************************
 phase
 -----

 A structure representing the full phase of an
 object (position and velocity) and also time,
 to help limit the number of variables that need
 to be passed to certain functions.

 All member variables are public, and may be
 accessed directly. Units should be in the
 default set: m for position, m/s for velocity,
 s for time.

 **********************************************/
{
	BRG_DISTANCE x, y, z;

	BRG_VELOCITY vx, vy, vz;

	BRG_TIME t;
	phase( CONST_BRG_DISTANCE_REF init_x = 0, CONST_BRG_DISTANCE_REF init_y = 0,
			CONST_BRG_DISTANCE_REF init_z = 0, CONST_BRG_VELOCITY_REF init_vx = 0,
			CONST_BRG_VELOCITY_REF init_vy = 0, CONST_BRG_VELOCITY_REF init_vz = 0,
	BRG_TIME init_t = 0 )
	{
		x = init_x;
		y = init_y;
		z = init_z;
		vx = init_vx;
		vy = init_vy;
		vz = init_vz;
		t = init_t;
	}
	const int set_phase( CONST_BRG_DISTANCE_REF init_x = 0, CONST_BRG_DISTANCE_REF init_y = 0,
			CONST_BRG_DISTANCE_REF init_z = 0, CONST_BRG_VELOCITY_REF init_vx = 0,
			CONST_BRG_VELOCITY_REF init_vy = 0, CONST_BRG_VELOCITY_REF init_vz = 0,
	BRG_TIME init_t = 0 )
	{
		x = init_x;
		y = init_y;
		z = init_z;
		vx = init_vx;
		vy = init_vy;
		vz = init_vz;
		t = init_t;

		return 0;
	}
};

#endif // end class declarations

}

#endif // __BRG_PHASE_HPP_INCLUDED__
