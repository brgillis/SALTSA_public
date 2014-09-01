/**********************************************************************\
  @file tNFW_profile_functors.hpp

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

#ifndef _BRG_TNFW_PROFILE_FUNCTORS_HPP_
#define _BRG_TNFW_PROFILE_FUNCTORS_HPP_

#include <cstdlib>

#include "brg/global.h"

#include "brg/physics/astro/density_profile/density_profile.h"
#include "brg/physics/astro/density_profile/tNFW_profile.h"

namespace brgastro
{

class tNFW_solve_rvir_iterative_functor
{
private:
	const tNFW_profile *_halo_;

public:

	tNFW_solve_rvir_iterative_functor(const tNFW_profile *halo)
	: _halo_(halo)
	{
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF in_param,	const bool silent = false ) const
	{
		if(in_param<=0)
		{
			return _halo_->rvir0();
		}
		else
		{
			return in_param * _halo_->enc_dens(in_param)/(virial_density_factor*_halo_->rho_crit());
		}
	}

};

class tNFW_solve_rvir_minimize_functor
{
private:
	const tNFW_profile *_halo_;

public:

	tNFW_solve_rvir_minimize_functor(const tNFW_profile *halo)
	: _halo_(halo)
	{
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param, const bool silent = false ) const
	{
		return std::fabs(virial_density_factor*_halo_->rho_crit()-_halo_->enc_dens(in_param));
	}

};

} // end namespace brgastro

#endif /* _BRG_TNFW_PROFILE_FUNCTORS_HPP_ */
