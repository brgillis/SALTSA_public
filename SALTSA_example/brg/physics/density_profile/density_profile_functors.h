/**********************************************************************\
  @file density_profile_functors.h

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

// body file: brg/physics/astro/density_profile/density_profile_functors.cpp

#ifndef _BRG_DENSITY_PROFILE_FUNCTORS_H_
#define _BRG_DENSITY_PROFILE_FUNCTORS_H_

#include "brg/global.h"

#include "brg/physics/density_profile/density_profile.h"
#include "brg/physics/units/unit_obj.h"

namespace brgastro {

class accel_functor
{
	/**********************************
	 accel_functor class
	 -----------------------------

	 Function class for acceleration within a density profile

	 Parent class: function_class (from brg_functors)

	 **********************************/
private:
	const density_profile *_host_ptr_;
public:

	void set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param, const bool silent = false ) const;

	accel_functor();
	accel_functor( const density_profile *init_host_ptr );
	virtual ~accel_functor()
	{
	}
};
// class accel_functor

class solve_rhm_functor
{
	/**********************************
	 solve_rhm_functor class
	 -----------------------------

	 Function class for solving the half-mass
	 radius of a halo.

	 Parent class: function_class (from brg_functors)

	 **********************************/
private:
	const density_profile *_host_ptr_;BRG_MASS _target_mass_;

public:

	void set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}
	void set_target_mass( const BRG_MASS &new_target_mass );
	const BRG_MASS & target_mass()
	{
		return _target_mass_;
	}

	BRG_UNITS operator ()( CONST_BRG_UNITS_REF in_param, const bool silent = false ) const;

	solve_rhm_functor();
	solve_rhm_functor( const density_profile *init_host,
			const BRG_MASS &init_target_mass );

};
// end class unitless_solve_rhm_functor

class spherical_density_functor
{
	/**********************************
	 spherical_density_functor class
	 -----------------------------

	 Function class integrating density in a sphere

	 Parent class: function_class (from brg_functors)

	 **********************************/
private:
	const density_profile *_host_ptr_;

public:

	void set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF in_param, const bool silent = false ) const;

	spherical_density_functor();
	spherical_density_functor( const density_profile *init_host );
	virtual ~spherical_density_functor()
	{
	}
};

} // namespace brgastro

#endif /* _BRG_DENSITY_PROFILE_FUNCTORS_H_ */
