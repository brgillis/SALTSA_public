/**********************************************************************\
  @file point_mass_profile.cpp

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

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "brg/global.h"

#include "brg/physics/units/unit_obj.h"

#include "point_mass_profile.h"

// brgastro::point_mass_profile class methods
#if (1)

#if (1) // Constructors
brgastro::point_mass_profile::point_mass_profile()
{
	_mass_ = 0;
}

brgastro::point_mass_profile::point_mass_profile( const BRG_MASS init_mass,
		const double init_z )
{
	set_mvir( init_mass );
	set_z( init_z );
}

#endif // End constructors

// Destructor
brgastro::point_mass_profile::~point_mass_profile()
{
}

#if (1) // Set functions

void brgastro::point_mass_profile::set_mvir(
		CONST_BRG_MASS_REF new_halo_mass, bool silent )
{
	_mass_ = new_halo_mass;
}
void brgastro::point_mass_profile::set_parameters( const std::vector< BRG_UNITS > &parameters, bool silent )
{
	assert(parameters.size()==2);
	set_mvir( parameters.at( 0 ) );
	set_z( parameters.at( 1 ) );
}

#endif // end set functions

BRG_MASS brgastro::point_mass_profile::mvir() const
{
	return _mass_;
}
BRG_MASS brgastro::point_mass_profile::mass() const
{
	return _mass_;
}

BRG_MASS brgastro::point_mass_profile::mtot() const
{
	return _mass_;
}

BRG_VELOCITY brgastro::point_mass_profile::vvir() const
{
	return std::pow( 10 * Gc * H() * mvir(), 1. / 3. );
}
BRG_DISTANCE brgastro::point_mass_profile::rvir() const
{
	return vvir() / H() / 10;
}
BRG_DISTANCE brgastro::point_mass_profile::rs() const
{
	return 0;
}
BRG_DISTANCE brgastro::point_mass_profile::rt(
		const bool silent) const
{
	return 0;
}

BRG_UNITS brgastro::point_mass_profile::dens(
		CONST_BRG_DISTANCE_REF r ) const
{
#ifdef _BRG_USE_UNITS_
	BRG_UNITS result(0,-3,0,1,0,0,0);
#else
	double result = 0;
#endif

	result = ( r == 0 ? DBL_MAX : 0 );

	return result;
}
BRG_UNITS brgastro::point_mass_profile::enc_dens(
		CONST_BRG_DISTANCE_REF r,
		const bool silent ) const
{
	return enc_mass( r ) / ( 4. / 3. * pi * cube( r ) );
}
BRG_MASS brgastro::point_mass_profile::enc_mass(
		CONST_BRG_DISTANCE_REF r,
		const bool silent ) const
{
	return _mass_;
}

std::vector< BRG_UNITS > brgastro::point_mass_profile::get_parameters( const bool silent ) const
{
	std::vector< BRG_UNITS > parameters(num_parameters());

	parameters[0] = _mass_;
	parameters[1] = z();
	return parameters;
}

std::vector< std::string > brgastro::point_mass_profile::get_parameter_names(
		const bool silent ) const
{
	std::vector< std::string > parameter_names(num_parameters());

	parameter_names.at( 0 ) = "mass";
	parameter_names.at( 1 ) = "z";
	return parameter_names;
}

void brgastro::point_mass_profile::truncate_to_fraction( const double f,
		const bool silent )
{
	if ( ( f <= 0 ) || ( isbad( f ) ) )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Bad, negative, or zero f passed to truncate_to_fraction.\n"
					<< "Will truncate to zero instead.\n";
		_mass_ = 0;
	}
	else
	{
		_mass_ *= f;
	}

}

#endif // end point_mass profile functions
