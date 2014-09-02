/**********************************************************************\
  @file astro_caches.cpp

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
#include <exception>
#include <string>

#include "brg/global.h"

#include "brg/math/cache/cache.hpp"

#include "brg/physics/astro.h"
#include "brg/physics/units/unit_obj.h"

#include "astro_caches.h"

/** Static Class Initialisation **/
#if (1)

// Initialisation for brgastro::dfa_cache
DEFINE_BRG_CACHE_STATIC_VARS( dfa_cache, 0, 5, 0.01 );

// Initialisation for brgastro::tfa_cache
DEFINE_BRG_CACHE_STATIC_VARS( tfa_cache, 0.001, 1.02, 0.001 );

#endif // end Static Class Initialisation

/** Class Method Definitions **/
#if (1)

// brgastro::dfa_cache class methods
#if (1)
const double brgastro::dfa_cache::_calculate( const double in_params ) const
{
	return brgastro::integrate_add( 0, in_params );
}

#endif // end brgastro::dfa_cache functions

const double brgastro::tfa_cache::_calculate( const double in_params ) const
{
	return -brgastro::integrate_ltd( 0, brgastro::zfa( in_params ) ) / c;
}

#endif // end class methods
