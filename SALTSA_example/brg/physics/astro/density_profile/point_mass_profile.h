/**********************************************************************\
  @file point_mass_profile.h

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

// body file: brg/physics/astro/density_profile/point_mass_profile.cpp

#ifndef _BRG_POINT_MASS_PROFILE_H_
#define _BRG_POINT_MASS_PROFILE_H_

#include <vector>

#include "brg/global.h"

#include "brg/physics/astro/density_profile/density_profile.h"
#include "brg/physics/units/unit_obj.h"

namespace brgastro {

class point_mass_profile: public density_profile
{
	/**********************************
	 point_mass_profile class
	 ------------------------

	 A virtual class for a point mass

	 Defined by two parameters:
	 mass and z

	 Parent class:
	 density profile

	 Derived classes:
	 (none)

	 **********************************/
private:
#if (1) // private member variables specific to this class

	BRG_MASS _mass_;

#endif // end private member variables

public:
#if (1) // public variables and functions

#if (1) // Constructors
	point_mass_profile();

	point_mass_profile( const BRG_MASS init_mass, const double init_z );

#endif // End constructors

	// Destructor
	~point_mass_profile();

#if (1) // Set functions
	virtual void set_mvir( CONST_BRG_MASS_REF new_halo_mass, bool silent =
			false );
	virtual void set_parameters( const std::vector< BRG_UNITS > &new_parameters,
			bool silent = false );
#endif // End set functions

#if (1) // Basic get functions
	BRG_MASS mass() const;

	BRG_MASS mvir() const;
	BRG_MASS mtot() const;

	BRG_DISTANCE rvir() const;
	BRG_DISTANCE rt(const bool silent = false) const;
	BRG_DISTANCE rs() const;

	BRG_VELOCITY vvir() const;

#endif // end basic get functions

#if (1) // advanced get functions
	BRG_UNITS dens( CONST_BRG_DISTANCE_REF r ) const;
	BRG_UNITS enc_dens( CONST_BRG_DISTANCE_REF r,
			const bool silent = false ) const;
	BRG_MASS enc_mass( CONST_BRG_DISTANCE_REF r, const bool silent =
				true ) const; // Mass enclosed with sphere of radius r
	size_t num_parameters() const
	{
		return 2; // Mass and redshift
	}
	std::vector< BRG_UNITS > get_parameters( const bool silent = false ) const;

	std::vector< std::string > get_parameter_names( const bool silent = false ) const;
#endif // end advanced get functions

#if (1) // Other operations

	void truncate_to_fraction( const double fraction,
			bool silent = false );
	virtual density_profile *density_profile_clone() const
	{
		return new point_mass_profile( *this );
	}
	virtual point_mass_profile *point_mass_profile_clone() const
	{
		return new point_mass_profile( *this );
	}

#endif

#endif // end public variables and functions
};
// class point_mass_profile

} // end namespace brgastro

#endif /* _BRG_POINT_MASS_PROFILE_H_ */
