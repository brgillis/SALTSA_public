/**********************************************************************\
  @file astro.cpp

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
#include <string>

#include "brg/global.h"

#include "brg/math/cache/cache.hpp"

#include "brg/physics/astro/astro_caches.h"
#include "brg/physics/astro/astro.h"
#include "brg/physics/units/unit_obj.h"

// brgastro::redshift_obj class methods
#if (1)

BRG_UNITS brgastro::redshift_obj::H() const
{
	// If not cached, calculate and cache it
	if(!_H_cached_)
	{
		if(_z_==0)
		{
			_H_cache_ = H_0;
		}
		else
		{
			// Friedmann equation, assuming omega = -1
			_H_cache_ = brgastro::H(_z_);
		}
		_H_cached_ = true;
	}
	return _H_cache_;
}

BRG_UNITS brgastro::redshift_obj::rho_crit() const
{
	return 3*square(H())/(8*pi*Gc);
}

#endif // brgastro::redshift_obj class methods

/** Global Function Definitions **/
#if (1)

BRG_UNITS brgastro::H( const double test_z )
{
	// Friedmann equation, assuming omega = -1
	if(test_z==0) return H_0;
	double zp1 = 1.+test_z;
	return H_0
			* std::sqrt( Omega_r * quart( zp1 )
							+ Omega_m * cube( zp1 )
							+ Omega_k * square( zp1 ) + Omega_l );
}

// dfa and afd functions
#if (1)

BRG_DISTANCE brgastro::dfa( CONST_BRG_ANGLE_REF da, const double z )
{
	return da * dfa_cache().get( z );
}
BRG_DISTANCE brgastro::dfa( CONST_BRG_ANGLE_REF a1, CONST_BRG_ANGLE_REF a2,
		const double z )
{
	return brgastro::dfa( a2 - a1, z );
}
BRG_DISTANCE brgastro::dfa( CONST_BRG_ANGLE_REF a1x, CONST_BRG_ANGLE_REF a1y,
		CONST_BRG_ANGLE_REF a2x, CONST_BRG_ANGLE_REF a2y, const double z )
{
	return brgastro::dfa( skydist2d( a1x, a1y, a2x, a2y ), z );
}
BRG_ANGLE brgastro::afd( CONST_BRG_DISTANCE_REF dd, const double z )
{
	return dd / safe_d( dfa_cache().get( z ) );
}
BRG_ANGLE brgastro::afd( CONST_BRG_DISTANCE_REF d1, CONST_BRG_DISTANCE_REF d2,
		const double z )
{
	return brgastro::afd( fabs( d2 - d1 ), z );
}
BRG_ANGLE brgastro::afd( CONST_BRG_DISTANCE_REF d1x,
		CONST_BRG_DISTANCE_REF d1y, CONST_BRG_DISTANCE_REF d2x,
		CONST_BRG_DISTANCE_REF d2y, const double z )
{
	return brgastro::afd( dist2d( d1x, d1y, d2x, d2y ), z );
}

double brgastro::zfa( const double a )
{
	return 1. / safe_d( a ) - 1.;
}
double brgastro::afz( const double z )
{
	return 1. / safe_d( 1 + z );
}

BRG_TIME brgastro::tfz( const double z )
{
	return brgastro::tfa( afz( z ) );
}
BRG_TIME brgastro::tfa( const double a )
{
	return tfa_cache().get( a );
}
double brgastro::zft( CONST_BRG_TIME_REF t )
{
	return brgastro::zfa( brgastro::aft( t ) );
}
double brgastro::aft( CONST_BRG_TIME_REF t )
{
	brgastro::tfa_cache cache;
	return cache.inverse_get( t );
}

#endif

// Functions to integrate out distances
#if (1)
double brgastro::integrate_add( const double z1, const double z2 )
{
	return brgastro::integrate_distance( z1, z2, 0, 100000 );
}
double brgastro::integrate_cmd( const double z1, const double z2 )
{
	return brgastro::integrate_distance( z1, z2, 1, 10000000 );
}
double brgastro::integrate_Ld( const double z1, const double z2 )
{
	return brgastro::integrate_distance( z1, z2, 2, 10000000 );
}
double brgastro::integrate_ltd( const double z1, const double z2 )
{
	return brgastro::integrate_distance( z1, z2, 3, 10000000 );
}
double brgastro::integrate_add( const double z )
{
	return brgastro::integrate_distance( 0, z, 0, 100000 );
}
double brgastro::integrate_cmd( const double z )
{
	return brgastro::integrate_distance( 0, z, 1, 10000000 );
}
double brgastro::integrate_Ld( const double z )
{
	return brgastro::integrate_distance( 0, z, 2, 10000000 );
}
double brgastro::integrate_ltd( const double z )
{
	return brgastro::integrate_distance( 0, z, 3, 10000000 );
}
double brgastro::integrate_distance( const double z1_init,
		const double z2_init, const int mode, const int n )
{
	// Function that will help calculate cosmic distances, thanks to Richard Powell - http://www.atlasoftheuniverse.com/
	// NOTE: Will return a negative value if z2 < z1. This is by design.

	double OK0;
	double OM0 = Omega_m, OR0 = Omega_r, OL0 = Omega_l;
	double OM, OR, OL, OK;
	double z1 = z1_init, z2 = z2_init;
	double HD; //Hubble distance in billions of lightyears
	double z, a, a1, a2, adot, h1;
	double DC = DBL_MAX, DCC = 0, DT = DBL_MAX, DTT = 0, DM;
	//double age, size;
	int i;
	short int sign = 1;

	if(z1==z2) return 0;

	OK0 = 1 - OM0 - OL0 - OR0;

	if ( z2 < z1 )
	{
		std::swap(z1,z2);
		sign = -1;
	}

	if ( z1 == 0 )
	{
		z = z2;
		h1 = H_0;
		OM = OM0;
		OR = OR0;
		OL = OL0;
		OK = OK0;
	}
	else
	{
		a1 = 1 / ( 1 + z1 );
		a2 = 1 / ( 1 + z2 );
		z = ( a1 / a2 ) - 1;
		h1 = H_0
				* std::sqrt( OR0 * inv_quart( a1 ) + OM0 * inv_cube( a1 )
								+ OK0 * inv_square( a1 ) + OL0 );
		OM = OM0 * square( H_0 / h1 ) * inv_cube( a1 );
		OR = OR0 * square( H_0 / h1 ) * inv_quart( a1 );
		OL = OL0 * square( H_0 / h1 );
		OK = 1 - OM - OR - OL;
	}

	HD = c / h1;

	for ( i = n; i >= 1; i-- )        // This loop is the numerical integration
	{
		a = ( i - 0.5 ) / n;              // Steadily decrease the scale factor
		// Comoving formula (See section 4 of Hogg, but I've added a radiation term too):
		adot = a * std::sqrt( OM * inv_cube(a) + OK * inv_square(a)
								+ OL + OR * inv_quart(a) ); // Note that "a" is equivalent to 1/(1+z)
		 // Collect DC and DT until the correct scale is reached
		DCC = DCC + 1 / ( a * adot ) / n; // Running total of the comoving distance
		DTT = DTT + 1 / adot / n; // Running total of the light travel time (see section 10 of Hogg)
		 // Collect DC and DT until the correct scale is reached
		DC = DCC;                 // Comoving distance DC
		DT = DTT;                 // Light travel time DT
		if ( a <= 1 / ( 1 + z ) ) // Collect DC and DT until the correct scale is reached
		{
			break;
		}
	}

	// Transverse comoving distance DM from section 5 of Hogg:
	if ( OK > 0.0001 )
		DM = ( 1 / std::sqrt( OK ) ) * sinh( std::sqrt( OK ) * DC );
	else if ( OK < -0.0001 )
		DM = ( 1 / std::sqrt( fabs( OK ) ) ) * sin( std::sqrt( fabs( OK ) ) * DC );
	else
		DM = DC;

	//age = HD*DTT;                 // Age of the universe (billions of years)
	//size = HD*DCC;                // Comoving radius of the observable universe

	switch ( mode )
	{
	case 0:
		return sign * HD * DM / ( 1 + z ); // Angular diameter distance (section 6 of Hogg)
		break;
	case 1:
		return sign * HD * DC;             // Comoving distance
		break;
	case 2:
		return sign * HD * DM * ( 1 + z );  // Luminosity distance (section 7 of Hogg)
		break;
	case 3:
		return sign * HD * DT;             // Light travel distance
		break;
	default:
		return sign * HD * DT;             // Light travel distance
	}
}
#endif

#endif // end Global function definitions
