/**********************************************************************\
  @file density_profile.h

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

// body file: brg/physics/astro/density_profile/density_profile.cpp

#ifndef _BRG_DENSITY_PROFILE_H_
#define _BRG_DENSITY_PROFILE_H_

#include <stdexcept>
#include <vector>

#include "brg/global.h"

#include "brg/math/misc_math.hpp"
#include "brg/math/safe_math.hpp"
#include "brg/physics/astro/astro.h"
#include "brg/physics/units/unit_obj.h"

namespace brgastro {

#ifndef _BRG_VIRIAL_DENSITY_FACTOR_DEFFED_
#define _BRG_VIRIAL_DENSITY_FACTOR_DEFFED_
const double virial_density_factor = 200;
#endif

class density_profile: public virtual redshift_obj
{
	/**********************************
	 density_profile class
	 -------------

	 An abstract class for anything which has
	 a mass density profile. Meant to be
	 overriden by child classes.

	 The methods for this class are sorted
	 into four groups:
	 -Pure virtual functions (must be implemented
	 by any derived classes)
	 -Virtual functions which must be overwritten if they're going
	 to be used (not pure virtual since they may not be needed)
	 -Virtual functions which should be overwritten if possible and
	 if they're going to be used (those that require integration,
	 but may have an analytic solution for certain profiles)
	 -Virtual functions that don't need to be overwritten in most
	 cases and non-virtual functions non-virtual functions

	 Derived classes:
	 tNFW_profile
	 point_mass profile

	 **********************************/

private:
#if (1)
	char _hm_type_;

	mutable BRG_DISTANCE _rhmvir_cache_, _rhmtot_cache_;

	BRG_VELOCITY _v_from_r( BRG_DISTANCE r ) const
	{
		BRG_UNITS a;

		a = accel( fabs( r ) );
		if ( a >= 0 )
			return 0;
		return std::sqrt( -a * fabs( r ) );
	}
#endif

protected:
	mutable bool hmvir_cached, hmtot_cached;

public:
#if (1)

#if (1) // Constructor
	density_profile()
	{
		_hm_type_ = 1;
		_rhmvir_cache_ = 0;
		hmvir_cached = false;
		_rhmtot_cache_ = 0;
		hmtot_cached = false;
	}
#endif

#if (1) // Destructor
	virtual ~density_profile()
	{
	}
#endif

#if (1) // Pure virtual functions (must be implemented by any derived classes)

	virtual BRG_MASS mvir() const =0; // Virial mass (exact definition can be chosen per profile)
	virtual BRG_UNITS dens( CONST_BRG_DISTANCE_REF r ) const =0; // Local density at radius r

	virtual density_profile *density_profile_clone() const =0; // Creates a clone of this
#endif

#if (1) // Virtual functions which must be overwritten if they're going to be used

#if (1) // Set functions - will return 1 if profile doesn't support this method of setting
	// All take values in default unit set (m, s, kg, K, rad, C)
	virtual void set_mvir( CONST_BRG_MASS_REF new_mvir, bool silent = false )
	{
		throw std::logic_error("density_profile::set_mvir(...) must be overloaded to be used.\n");
	}
	virtual void set_vvir( CONST_BRG_VELOCITY_REF new_vvir, bool silent =
			false )
	{
		throw std::logic_error("density_profile::set_vvir(...) must be overloaded to be used.\n");
	}
	virtual void set_rvir( CONST_BRG_DISTANCE_REF new_rvir, bool silent =
			false )
	{
		throw std::logic_error("density_profile::set_rvir(...) must be overloaded to be used.\n");
	}
	virtual void set_rs( CONST_BRG_DISTANCE_REF new_rs, bool silent = false ) // Scale radius
	{
		throw std::logic_error("density_profile::set_rs(...) must be overloaded to be used.\n");
	}
	virtual void set_rt( CONST_BRG_DISTANCE_REF new_rt, bool silent = false ) // Tidal/truncation radius
	{
		throw std::logic_error("density_profile::set_rt(...) must be overloaded to be used.\n");
	}
	virtual void set_parameters( const std::vector< BRG_UNITS > &parameters, bool silent = false )
	{
		throw std::logic_error("density_profile::set_parameters(...) must be overloaded to be used.\n");
	}
	virtual void set_tau( const double new_tau, bool silent = false ) // Truncation parameter
	{
		throw std::logic_error("density_profile::set_tau(...) must be overloaded to be used.\n");
	}
	virtual void set_c( const double new_c, bool silent = false ) // Concentration
	{
		throw std::logic_error("density_profile::set_c(...) must be overloaded to be used.\n");
	}
	void override_rhmvir( CONST_BRG_DISTANCE_REF new_rhmvir, bool silent =
			false )
	{
		_rhmvir_cache_ = new_rhmvir;
		hmvir_cached = true;
	}
	void override_rhmtot( CONST_BRG_DISTANCE_REF new_rhmtot, bool silent =
			false )
	{
		_rhmtot_cache_ = new_rhmtot;
		hmtot_cached = true;
	}

#endif

#if (1) // Basic get functions
	// All return values in default unit set (m, s, kg, K, rad, C)

	virtual BRG_MASS mtot() const // Total mass
	{
		return 0;
	}

	virtual BRG_DISTANCE rt( const bool silent = false ) const // Tidal/truncation radius
	{
		throw std::logic_error("density_profile::rt() must be overloaded to be used.\n");
	}

#endif // end basic get functions

	virtual unsigned int num_parameters() const
	{
		throw std::logic_error("density_profile::num_parameters() must be overloaded to be used.\n");
	}

	virtual std::vector< BRG_UNITS > get_parameters( const bool silent = false ) const // Returns a set of parameters for this halo (all the variables needed to define it - must be defined for each child)
	{
		throw std::logic_error("density_profile::get_parameters() must be overloaded to be used.\n");
	}

	virtual std::vector< std::string > get_parameter_names( const bool silent = false ) const // Returns a set of names of this halo's parameters
	{
		throw std::logic_error("density_profile::get_parameter_names() must be overloaded to be used.\n");
	}

#endif // end Virtual functions which must be overwritten if they're going to be used

#if (1) // Virtual functions which should be overwritten if at all possible and if they'll be used
	virtual BRG_MASS enc_mass( CONST_BRG_DISTANCE_REF r, const bool silent =
			true ) const; // Mass enclosed with sphere of radius r
#endif

#if (1) // Virtual functions which shouldn't be overwritten in most cases and non-virtual functions

	virtual BRG_DISTANCE rvir() const // Virial radius (exact definition can be chosen per profile if preferred)
	{

		return safe_pow(
				2 * mvir() * Gc / ( square( H() ) * virial_density_factor ),
				1. / 3. );
	}
	BRG_MASS hmvir() const // Half virial mass
	{
		return mvir() / 2;
	}
	BRG_MASS hmtot() const // Half total mass
	{
		return mtot() / 2;
	}
	virtual BRG_MASS hm() const // Half mass (which depends on hm_type value)
	{
		return ( _hm_type_ == 0 ? hmvir() : hmtot() );
	}
	virtual BRG_UNITS enc_dens( CONST_BRG_DISTANCE_REF r,
			const bool silent = false ) const // Mean density enclosed with sphere of radius r
	{
		BRG_DISTANCE r_to_use = max( std::fabs( r ), SMALL_FACTOR );
		return enc_mass( r_to_use, silent )
				/ ( 4. / 3. * pi * cube( std::fabs( r_to_use ) ) );
	}
	virtual BRG_DISTANCE rhmvir( const bool silent = false ) const; // Half-virial-mass radius
	virtual BRG_DISTANCE rhmtot( const bool silent = false ) const; // Half-total-mass radius
	virtual BRG_DISTANCE rhm( const bool silent = false ) const // Half-mass radius
	{
		return ( _hm_type_ == 0 ? rhmvir( silent ) : rhmtot( silent ) );
	}
	virtual BRG_VELOCITY vvir() const // Orbital velocity at rvir
	{
		return _v_from_r( rvir() );
	}
	virtual BRG_VELOCITY vhmvir( const bool silent = false ) const // Orbital velocity at rhmvir
	{
		return _v_from_r( rhmvir( silent ) );
	}
	virtual BRG_VELOCITY vhmtot( const bool silent = false ) const // Orbital velocity at rhmtot
	{
		return _v_from_r( rhmtot( silent ) );
	}
	virtual BRG_VELOCITY vhm( const bool silent = false ) const // Orbital velocity at rhm
	{
		return ( _hm_type_ == 0 ? vhmvir( silent ) : vhmtot( silent ) );
	}
	virtual BRG_VELOCITY vt( const bool silent = false ) const // Orbital velocity at rt
	{
		return _v_from_r( rt( silent ) );
	}
	BRG_TIME otvir() const // Orbital period at rvir
	{
		return 2 * pi * rvir() / vvir();
	}
	BRG_TIME othmvir( const bool silent = false ) const // Orbital period at rhmvir
	{
		return 2 * pi * rhmvir( silent ) / safe_d(vhmvir( silent ));
	}
	BRG_TIME othmtot( const bool silent = false ) const // Orbital period at rhmtot
	{
		return 2 * pi * rhmtot( silent ) / safe_d(vhmtot( silent ));
	}
	BRG_TIME othm( const bool silent = false ) const // Orbital period at rhm
	{
		return 2 * pi * rhm( silent ) / safe_d(vhm( silent ));
	}
	BRG_TIME ott( const bool silent = false ) const // Orbital period at rt
	{
		return 2 * pi * rt( silent ) / safe_d(vt( silent ));
	}

	void set_hm_type( int new_hm_type ) // Whether default half-mass is half virial, half total, or something else
	{
		_hm_type_ = new_hm_type;
	}

#if (1) // advanced get functions

	BRG_UNITS accel( CONST_BRG_DISTANCE_REF r,
			const bool silent = false ) const // Gravitational acceleration at radius r
	{
		if(r==0)
			return 0;
		else
			return -Gc * enc_mass( r, silent ) / square( r );
	}
	virtual BRG_UNITS Daccel( CONST_BRG_DISTANCE_REF r, const bool silent =
			false ) const; // Derivative of acceleration at radius r

#endif // end advanced get functions

#endif // end Virtual functions which shouldn't be overwritten in most cases and non-virtual functions

#if (1) // Other operations

	virtual void truncate_to_fraction( const double fraction,
			bool silent = false ) // Adjusts parameters of this class to decrease the mass to fraction of its previous mass. Must be defined for each child
	{
		throw std::logic_error("density_profile::truncate_to_fraction() must be overloaded to be used.\n");
	}

#endif

#endif
}; // end class density_profile

// Function to estimate orbital period from current position and velocity in a density profile
// Note that this is merely an estimate from analogy to calculations in a Keplerian potential
BRG_TIME period( const density_profile *host, CONST_BRG_DISTANCE_REF r,
		CONST_BRG_VELOCITY_REF vr, CONST_BRG_VELOCITY_REF vt = 0 );

} // end namespace brgastro

#endif /* _BRG_DENSITY_PROFILE_H_ */
