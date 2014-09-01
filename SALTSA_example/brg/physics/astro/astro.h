/**********************************************************************\
  @file astro.h

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

// body file: astro.cpp

#ifndef __BRG_ASTRO_H_INCLUDED__
#define __BRG_ASTRO_H_INCLUDED__

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>

#include "brg/global.h"

#include "brg/physics/units/unit_conversions.hpp"
#include "brg/physics/units/unit_obj.h"

/** Constant Definitions **/
#if (1)
/**********************************************************************
 This section stores the definitions of various physical and
 astrophysical constants, stored as variables. Key points:

 -H_0 is set to 70 km/s/Mpc. If this is not changed, all results should
 be assumed to be in h_70 units.
 -Cosmological values are based on best-fit WMAP + priors values
 from Hinshaw et al. 2012.
 -If the _BRG_USE_UNITS_ parameter is set (see the brg_global.h file), values
 that have units will be declared as unit_objs.

 \**********************************************************************/

// Physical Constants
#if (1)
#ifndef _BRG_PI_DEFINED_
#define _BRG_PI_DEFINED_
const double pi = 3.14159265358979323846;
#endif


namespace brgastro
{

#ifndef __PHYS_VARS_DEFINED__
#define __PHYS_VARS_DEFINED__

#ifdef _BRG_USE_UNITS_
const brgastro::unit_obj Gc(6.67384e-11,3,-2,-1,0,0); // In m^3 s^-2 kg^-1
#else
const float Gc = 6.67384e-11; // In m^3 s^-2 kg^-1
#endif
const BRG_VELOCITY c = brgastro::unitconv::ctomps;

#endif

#endif // end Physical Constants

// Cosmological Parameters
#if (1)
#ifndef __COSMO_VARS_DEFINED__
#define __COSMO_VARS_DEFINED__

#ifdef _BRG_USE_UNITS_
const brgastro::unit_obj H_0(70*unitconv::kmtom/unitconv::stos/unitconv::Mpctom,0,-1,0,0,0,0); // So all results will implicitly be in h_70 units
#else
const double H_0 = 70 * brgastro::unitconv::kmtom / brgastro::unitconv::stos / brgastro::unitconv::Mpctom; // So all results will implicitly be in h_70 units
#endif

const float Omega_m = 0.288; // WMAP9 + priors
const float Omega_r = 0.000086; // WMAP9 + priors
const float Omega_k = 0; // Assuming flat space
const float Omega_l = 1 - Omega_k - Omega_m - Omega_r;
const float Omega_b = 0.0472; // WMAP9 + priors
const float sigma_8 = 0.830; // WMAP9 + priors
const float n_s = 0.971;     // WMAP9 + priors

#endif

#endif // end Cosmological Parameters

// Default parameters
#if (1)
#ifndef __BRG_DEFAULT_ASTRO_VARS_DEFINED__
#define __BRG_DEFAULT_ASTRO_VARS_DEFINED__

const float default_c = 6; // To help prevent crashes. Warning will be output
const float default_tau_factor = 2;

#endif
#endif

#endif // End Constant Definitions

/** Function Declarations **/
#if (1)
/**********************************************************************
 This section defines various functions which I find useful for
 astrophysical calculations. All are declared in the namespace brgastro.

 \**********************************************************************/

BRG_UNITS H( const double z );

// Functions to get grid integers or grid boundaries from integers
int get_ra_grid( CONST_BRG_ANGLE_REF ra );
int get_dec_grid( CONST_BRG_ANGLE_REF dec );
int get_z_grid( const double z );

BRG_ANGLE get_ra_grid_lower( const int ra_grid );
BRG_ANGLE get_dec_grid_lower( const int dec_grid );
double get_z_grid_lower( const int z_grid );

BRG_ANGLE get_ra_grid_upper( const int ra_grid );
BRG_ANGLE get_dec_grid_upper( const int dec_grid );
double get_z_grid_upper( const int z_grid );

BRG_ANGLE get_ra_grid_mid( const int ra_grid );
BRG_ANGLE get_dec_grid_mid( const int dec_grid );
double get_z_grid_mid( const int z_grid );

// Functions to get transverse distance (in m) from angle (in rad) or vice-versa
BRG_DISTANCE dfa( CONST_BRG_ANGLE_REF da, const double z );
BRG_DISTANCE dfa( CONST_BRG_ANGLE_REF a1, CONST_BRG_ANGLE_REF a2,
		const double z );
BRG_DISTANCE dfa( CONST_BRG_ANGLE_REF a1x, CONST_BRG_ANGLE_REF a1y,
		CONST_BRG_ANGLE_REF a2x, CONST_BRG_ANGLE_REF a2y, const double z );

BRG_ANGLE afd( CONST_BRG_DISTANCE_REF dd, const double z );
BRG_ANGLE afd( CONST_BRG_DISTANCE_REF d1, CONST_BRG_DISTANCE_REF d2,
		const double z );
BRG_ANGLE afd( CONST_BRG_DISTANCE_REF d1x, CONST_BRG_DISTANCE_REF d1y,
		CONST_BRG_DISTANCE_REF d2x, CONST_BRG_DISTANCE_REF d2y, const double z );

// Functions to work between redshift, scale factor, and time (in s, with zero = present day)
double zfa( const double a );
double afz( const double z );

BRG_TIME tfz( const double z );
BRG_TIME tfa( const double z );
double zft( CONST_BRG_TIME_REF t );
double aft( CONST_BRG_TIME_REF t );

// Functions to integrate out distances
double integrate_add( const double z1, const double z2 );
double integrate_cmd( const double z1, const double z2 );
double integrate_Ld( const double z1, const double z2 );
double integrate_ltd( const double z1, const double z2 );
double integrate_add( const double z );
double integrate_cmd( const double z );
double integrate_Ld( const double z );
double integrate_ltd( const double z );
double integrate_distance( const double z1, const double z2,
		const int mode, const int resolution = 10000 );

// Lensing functions
BRG_DISTANCE ad_distance( double z1, double z2 = 0 );
BRG_UNITS sigma_crit( const double z_lens, const double z_source );

// Like dist2d, but using corrections for spherical geometry
template< typename Tr1, typename Td1, typename Tr2, typename Td2 >
inline const Tr1 skydist2d( const Tr1 ra1, const Td1 dec1,
		const Tr2 ra2, const Td2 dec2 )
{
	return quad_add( ( ra2 - ra1 ) * std::cos( ( dec2 + dec1 ) / 2 ), dec2 - dec1 );
}

#endif // end function declarations

/** Class Definitions **/
#if (1)
class redshift_obj
{
	/**********************************
	 redshift_obj class
	 ------------------

	 A base class for anything with a redshift.

	 Derived classes:

	 sky_obj
	 galaxy
	 group
	 density_profile
	 tNFW_profile
	 point_mass_profile

	 **********************************/
private:
	double _z_, _z_err_;
	mutable BRG_UNITS _H_cache_;
	mutable bool _H_cached_;
	mutable int _z_grid_;
	mutable bool _z_grid_cached_;
	mutable int _local_z_grid_change_num_;
public:

	// Constructor
	redshift_obj( const double init_z = 0, const double init_z_err = 0 )
	{
		_z_ = init_z;
		_z_err_ = init_z_err;
		_H_cache_ = 0;
		_H_cached_ = false;
		_z_grid_ = 0;
		_z_grid_cached_ = false;
		_local_z_grid_change_num_ = -1;
	}

	// Copy constructor
	// (implicit is fine for us)

	// Virtual destructor
	virtual ~redshift_obj()
	{
	}

	// Set functions
#if (1)
	virtual void set_z( const double new_z ) // Sets z
	{
		_z_ = new_z;
		_H_cached_ = false;
		_z_grid_cached_ = false;
	}
	virtual void set_z_err( const double new_z_err ) // Sets z error
	{
		_z_err_ = new_z_err;
	}
#endif

	// Get functions
#if (1)
	double z() const
	{
		return _z_;
	}
	double z_err() const
	{
		return _z_err_;
	}
	int z_grid() const;
#endif

	// Astro calculations

#if (1)
	// H at the redshift of the object, given in units of m/s^2
	BRG_UNITS H() const;

	// Critical density at this redshift
	BRG_UNITS rho_crit() const;

#endif

	// Clone function - Not needed in current implementation
	// virtual redshift_obj *redshift_obj_clone()=0;
};

#endif // end Class Definitions

} // end namespace brgastro

#endif
