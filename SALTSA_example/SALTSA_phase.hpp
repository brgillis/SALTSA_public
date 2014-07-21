/**       @file SALTSA_phase.hpp
 *
 *     Project: SALTSA_example
 *        Path: /SALTSA_example/SALTSA_phase.hpp
 *
 *  Created on: 18 Jul 2014
 *      Author: brg
 */

#ifndef SALTSA_PHASE_HPP_
#define SALTSA_PHASE_HPP_
namespace SALTSA {

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
	double x, y, z;
	double vx, vy, vz;
	double t;

	phase( double init_x = 0, double init_y = 0,
	double init_z = 0, double init_vx = 0,
	double init_vy = 0, double init_vz = 0,
	double init_t = 0 )
	{
		x = init_x;
		y = init_y;
		z = init_z;
		vx = init_vx;
		vy = init_vy;
		vz = init_vz;
		t = init_t;
	}
	const int set_phase( double init_x = 0, double init_y = 0,
	double init_z = 0, double init_vx = 0,
	double init_vy = 0, double init_vz = 0,
	double init_t = 0 )
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

} // end namespace SALTSA



#endif /* SALTSA_PHASE_HPP_ */
