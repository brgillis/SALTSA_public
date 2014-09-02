/**********************************************************************\
  @file solve_rt_functors.cpp

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

#include "brg/global.h"

#include "brg/physics/density_profile/density_profile.h"
#include "brg/physics/units/unit_obj.h"
#include "brg/utility.hpp"

#include "solve_rt_functors.h"

// brgastro::solve_rt_grid_functor class method implementations
#if (1)

BRG_UNITS brgastro::solve_rt_grid_functor::operator()(
		CONST_BRG_UNITS_REF  in_param, const bool silent ) const
{
	BRG_DISTANCE r;
	BRG_MASS delta_M;
	r = std::fabs( in_param );

	delta_M = sum_delta_rho * 4. / 3. * pi * cube( r );

	if ( r == 0 )
	{
		return DBL_MAX;
	}
	else
	{
		return std::fabs( Gc * ( satellite_ptr->enc_mass( r ) - delta_M )
						/ safe_d(( omega * omega + Daccel ) * cube( r ) ) - 1 );
	}
	return 0;
}

brgastro::solve_rt_grid_functor::solve_rt_grid_functor()
{
	omega = 0;
	Daccel = 0;
	sum_delta_rho = 0;
	satellite_ptr = NULL;
	return;
}
brgastro::solve_rt_grid_functor::solve_rt_grid_functor(
		const BRG_UNITS init_omega, const density_profile *init_satellite,
		const BRG_UNITS init_Daccel, const long double init_sum_delta_rho )
{
	omega = init_omega;
	satellite_ptr = init_satellite;
	Daccel = init_Daccel;
	sum_delta_rho = init_sum_delta_rho;
	return;
}

#endif // end brgastro::solve_rt_grid_functor class function definitions

// brgastro::solve_rt_grid_functor class method implementations
#if(1)

BRG_UNITS brgastro::solve_rt_it_functor::operator()(
		CONST_BRG_UNITS_REF in_param, const bool silent ) const
{
	BRG_DISTANCE r = 0;
	BRG_MASS delta_M = 0;
	BRG_UNITS r3 = 0;

	r = std::fabs( in_param );

	delta_M = sum_delta_rho * 4. / 3. * pi * cube( r );

	r3 = Gc * ( satellite_ptr->enc_mass( r ) - delta_M )
			/ safe_d( omega * omega + Daccel );
	if ( r3 <= 0 )
	{
		return r * 0.9;
	}
	else
	{
		return std::pow( r3, 1. / 3. );
	}
	return 0;
}

brgastro::solve_rt_it_functor::solve_rt_it_functor()
{
	omega = 0;
	satellite_ptr = NULL;
	Daccel = 0;
	sum_delta_rho = 0;
	return;
}
brgastro::solve_rt_it_functor::solve_rt_it_functor(
		const BRG_UNITS init_omega, const density_profile *init_satellite,
		const BRG_UNITS init_Daccel, const long double init_sum_delta_rho )
{
	omega = init_omega;
	satellite_ptr = init_satellite;
	Daccel = init_Daccel;
	sum_delta_rho = init_sum_delta_rho;
	return;
}

#endif // end brgastro::solve_rt_it_functor class function definitions
