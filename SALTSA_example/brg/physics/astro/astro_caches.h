/**********************************************************************\
  @file astro_caches.h

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

// body file: brg/physics/astro/astro_caches.cpp

#ifndef _BRG_ASTRO_CACHES_H_INCLUDED_
#define _BRG_ASTRO_CACHES_H_INCLUDED_

#include <string>

#include "brg/global.h"

#include "brg/math/cache/cache.hpp"
#include "brg/physics/units/unit_obj.h"

namespace brgastro
{

class dfa_cache : public brg_cache<dfa_cache>
{
	// "Distance from angle" cache
private:

	DECLARE_BRG_CACHE_STATIC_VARS();

	friend class brg_cache<dfa_cache>;

protected:

	const std::string _name_base() const throw()
	{
		return "dfa";
	}

#ifdef _BRG_USE_UNITS_

	// Tells what units the result should have. Only the units matter in the return, not the value
	const brgastro::unit_obj _units() const throw()
	{
		return brgastro::unit_obj(0,0,1,0,0,0,0);
	}
	const brgastro::unit_obj _inverse_units() const throw()
	{
		return brgastro::unit_obj(0);
	}

#endif // _BRG_USE_UNITS_

	// Long-form calculation function.
	const double _calculate( const double in_param ) const;

public:

};
// class dfa_cache

class tfa_cache : public brg_cache<tfa_cache>
{
	// "Time from a (scale factor)" cache
private:

	DECLARE_BRG_CACHE_STATIC_VARS();

	friend class brg_cache<tfa_cache>;

protected:

	const std::string _name_base() const throw()
	{
		return "tfa";
	}

#ifdef _BRG_USE_UNITS_

	// Tells what units the result should have. Only the units matter in the return, not the value
	const brgastro::unit_obj _units() const throw()
	{
		return brgastro::unit_obj(0,0,1,0,0,0,0);
	}
	const brgastro::unit_obj _inverse_units() const throw()
	{
		return brgastro::unit_obj(0);
	}

#endif // _BRG_USE_UNITS_

	// Long-form calculation function.
	const double _calculate( const double in_params ) const;

public:

};
// class tfa_cache

} // end namespace brgastro

#endif // __BRG_ASTRO_CACHES_H_INCLUDED__

