/**********************************************************************\
  @file density_profile.cpp

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

#include <iostream>

#include "brg/global.h"

#include "brg/math/calculus/integrate.hpp"
#include "brg/math/solvers/solvers.hpp"
#include "brg/physics/astro/density_profile/density_profile_functors.h"
#include "brg/physics/units/unit_obj.h"

#include "density_profile.h"

// brgastro::density_profile class methods
#if (1)
BRG_UNITS brgastro::density_profile::Daccel( CONST_BRG_DISTANCE_REF r,
		const bool silent ) const
{
	BRG_DISTANCE dr;
	BRG_UNITS a1, a2;
	// It's simplest, and a ton faster to just manually differentiate here.
	dr = max( r * SMALL_FACTOR, SMALL_FACTOR );
	a1 = accel( r, silent );
	a2 = accel( r + dr, silent );
	return ( a2 - a1 ) / safe_d( dr );
}

BRG_DISTANCE brgastro::density_profile::rhmtot( const bool silent ) const
{
	// If cached, return the cached value
	if ( hmtot_cached )
		return _rhmtot_cache_;

	// Not cached, so we'll have to calculate it
	double target_mass = hmtot();
	solve_rhm_functor func( this, target_mass );
	solve_rhm_functor *func_ptr = &func;

	BRG_UNITS max_r( default_tau_factor * rvir() );
	BRG_UNITS rhm_test( 0 );

	// First check for zero mass/radius/density
	if ( ( mvir() <= 0 ) || ( rvir() <= 0 ) || ( dens( rvir() / 2 ) < 0 ) )
	{
		hmtot_cached = true;
		return _rhmtot_cache_ = 0;
	}

	try
	{
		rhm_test = solve_grid( func_ptr, static_cast<BRG_UNITS>(0), max_r, 10, static_cast<BRG_UNITS>(0.) );
	}
	catch(const std::exception &e)
	{
		if ( !silent )
			std::cerr << "WARNING: Could not solve half-mass radius. Assuming it's zero.\n";

		_rhmtot_cache_ = 0;
		hmtot_cached = true;
		return _rhmtot_cache_;
	}
	_rhmtot_cache_ = std::fabs( rhm_test );
	hmtot_cached = true;
	return _rhmtot_cache_;
}

BRG_DISTANCE brgastro::density_profile::rhmvir( const bool silent ) const
{
	// If cached, return the cached value
	if ( hmvir_cached )
		return _rhmvir_cache_;

	// Not cached, so we'll have to calculate it
	double target_mass = hmvir();
	solve_rhm_functor func( this, target_mass );
	solve_rhm_functor *func_ptr = &func;

	BRG_UNITS max_r( default_tau_factor * rvir() );
	BRG_UNITS rhm_test( 0 );

	// First check for zero mass/radius/density
	if ( ( mvir() <= 0 ) || ( rvir() <= 0 ) || ( dens( rvir() / 2 ) < 0 ) )
	{
		hmvir_cached = true;
		return _rhmvir_cache_ = 0;
	}

	try
	{
		rhm_test = solve_grid( func_ptr, static_cast<BRG_UNITS>(0.), max_r, 10, static_cast<BRG_UNITS>(0.));
	}
	catch(const std::exception &e)
	{
		if ( !silent )
			std::cerr << "WARNING: Could not solve half-mass radius.\n";
		_rhmvir_cache_ = 0;
		return _rhmvir_cache_;
	}
	_rhmvir_cache_ = max(0,std::fabs( rhm_test ));
	hmvir_cached = true;
	return _rhmvir_cache_;
}

BRG_MASS brgastro::density_profile::enc_mass( CONST_BRG_DISTANCE_REF r,
		const bool silent ) const
{
	if ( r == 0 )
		return 0;
	BRG_DISTANCE r_to_use = std::fabs( r );
	brgastro::spherical_density_functor func( this );
	BRG_UNITS min_in_params( 0 ), max_in_params( r_to_use ), out_params( 0 );
	out_params = brgastro::integrate_Romberg( &func,min_in_params,
			max_in_params, 0.00001, false, silent );
	return out_params;
}

#endif

BRG_TIME brgastro::period( const brgastro::density_profile *host,
		CONST_BRG_DISTANCE_REF r, CONST_BRG_VELOCITY_REF vr, CONST_BRG_VELOCITY_REF vt )
{
	BRG_UNITS mu = host->enc_mass( r ) * Gc;
	BRG_VELOCITY v = quad_add( vr, vt );
	BRG_DISTANCE a = -mu / 2 / safe_d( v * v / 2 - mu / safe_d( r ) );
	BRG_TIME result = (
			a > 0 ? 2 * pi * std::sqrt( cube(a) / mu ) : BRG_TIME( 0 ) );
	return result;
}
