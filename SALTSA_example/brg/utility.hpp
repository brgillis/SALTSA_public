/**********************************************************************\
  @file utility.hpp

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

#ifndef _BRG_UTILITY_HPP_INCLUDED_
#define _BRG_UTILITY_HPP_INCLUDED_

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "global.h"

namespace brgastro
{

// Generic functions
#if (1)

// Set_zero function - a way for other template functions to "clear" or initialize a value in various ways
inline void set_zero( int & obj )
{
	obj = 0;
}
inline void set_zero( long int & obj )
{
	obj = 0;
}
inline void set_zero( short int & obj )
{
	obj = 0;
}
inline void set_zero( unsigned int & obj )
{
	obj = 0;
}
inline void set_zero( unsigned long int & obj )
{
	obj = 0;
}
inline void set_zero( unsigned short int & obj )
{
	obj = 0;
}
inline void set_zero( double & obj )
{
	obj = 0;
}
inline void set_zero( long double & obj )
{
	obj = 0;
}
inline void set_zero( float & obj )
{
	obj = 0;
}
#ifdef _BRG_USE_UNITS_
inline void set_zero( unit_obj obj)
{
	obj = unit_obj(0);
}
#endif
inline void set_zero( std::string obj )
{
	obj = "";
}
template< class T >
inline void set_zero( std::vector< T > vec )
{
	vec.clear();
}
template< class T >
inline void set_zero( T *obj )
{
	obj = NULL;
}
template< class obj_type >
inline void set_zero( obj_type obj )
{
	obj = 0;
}

// Various "make" functions, to allocate dynamic memory.
// After allocating memory, these functions initialize the new variables using the
// set_zero function (see above).

template< class obj_type >
inline void make_obj( BRG_UNIQUE_PTR<obj_type> & obj_pointer, const bool silent=false )
{
	obj_pointer = BRG_UNIQUE_PTR<obj_type>(new obj_type);
	set_zero(*obj_pointer);
}

#ifdef _BRG_USE_CPP_11_STD_
template <class array_type>
inline void make_array1d( std::unique_ptr<array_type []> & array_pointer, const unsigned int num_elem,
		const bool silent = false )
{
	array_pointer = std::unique_ptr<array_type []>(new array_type [num_elem]);
	for(int i=0; i<num_elem; i++) set_zero(array_pointer[i]);
}
template <class array_type>
inline void make_array( std::unique_ptr<array_type []> & array_pointer, const unsigned int num_elem,
		const bool silent = false )
{
	make_array1d( array_pointer, num_elem );
}
#endif
template< class array_type >
inline void make_array1d( std::vector< array_type > & array_pointer,
		const unsigned int num_elem, const bool silent = false )
{
	array_pointer.resize( num_elem );
	for ( unsigned int i = 0; i < num_elem; i++ )
		set_zero( array_pointer[i] );
}
template< class array_type >
inline void make_array( std::vector< array_type > & array_pointer,
		const int num_elem, const bool silent = false )
{
	return make_array1d( array_pointer, num_elem );
}

#ifdef _BRG_USE_CPP_11_STD_
template <class array_type>
inline void make_array2d( std::unique_ptr<std::unique_ptr<array_type []> []> & array_pointer, const unsigned int num_elem1, const unsigned int num_elem2, const bool silent=false )
{
	array_pointer = std::unique_ptr<std::unique_ptr<array_type []> []>(new std::unique_ptr<array_type []> [num_elem1]);
	for( unsigned int i = 0; i < num_elem1; i++)
	{
		make_array(array_pointer[i], num_elem2);
	}
}
#endif
template< class array_type >
inline void make_array2d(
		std::vector< std::vector< array_type > > & array_pointer,
		const unsigned int num_elem1, const unsigned int num_elem2, const bool silent = false )
{
	array_pointer.resize( num_elem1 );
	for ( unsigned int i = 0; i < num_elem1; i++ )
	{
		make_array( array_pointer[i], num_elem2 );
	}
}

#ifdef _BRG_USE_CPP_11_STD_
template <class array_type>
inline void make_array3d( std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> & array_pointer, const unsigned int num_elem1,
		const unsigned int num_elem2, const unsigned int num_elem3, const bool silent=false )
{
	array_pointer = std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []>(new
			std::unique_ptr<std::unique_ptr<array_type []> []> [num_elem1]);

	for( unsigned int i = 0; i < num_elem1; i++)
	{
		make_array2d(array_pointer[i], num_elem2, num_elem3);
	}
}
#endif
template< class array_type >
inline void make_array3d(
		std::vector< std::vector< std::vector< array_type > > > & array_pointer,
		const unsigned int num_elem1, const unsigned int num_elem2, const unsigned int num_elem3,
		const bool silent = false )
{
	array_pointer.resize( num_elem1 );
	for ( unsigned int i = 0; i < num_elem1; i++ )
	{
		make_array2d( array_pointer[i], num_elem2, num_elem3 );
	}
}

#ifdef _BRG_USE_CPP_11_STD_
template <class array_type>
inline void make_array4d( std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> []> & array_pointer,
		const unsigned int num_elem1, const unsigned int num_elem2, const unsigned int num_elem3, const unsigned int num_elem4, const bool silent=false )
{
	array_pointer = std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> []>(new
			std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> [num_elem1]);

	for( unsigned int i = 0; i < num_elem1; i++)
	{
		make_array3d(array_pointer[i], num_elem2, num_elem3, num_elem4);
	}
}
#endif
template< class array_type >
inline void make_array4d(
		std::vector< std::vector< std::vector< std::vector< array_type > > > > & array_pointer,
		const unsigned int num_elem1, const unsigned int num_elem2, const unsigned int num_elem3,
		const unsigned int num_elem4, const bool silent = false )
{
	array_pointer.resize( num_elem1 );
	for ( unsigned int i = 0; i < num_elem1; i++ )
	{
		make_array3d( array_pointer[i], num_elem2, num_elem3, num_elem4 );
	}
}

#ifdef _BRG_USE_CPP_11_STD_
template <class array_type>
inline void make_array5d( std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> []> []>
		& array_pointer, const unsigned int num_elem1, const unsigned int num_elem2, const unsigned int num_elem3,
		const unsigned int num_elem4, const unsigned int num_elem5, const bool silent=false )
{
	array_pointer = std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> []> []>(new
			std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> []> [num_elem1]);

	for( int i = 0; i < num_elem1; i++)
	{
		make_array4d(array_pointer[i], num_elem2, num_elem3, num_elem4, num_elem5);
	}
}
#endif
template< class array_type >
inline void make_array5d(
		std::vector<
				std::vector<
						std::vector< std::vector< std::vector< array_type > > > > > & array_pointer,
		const unsigned int num_elem1, const unsigned int num_elem2, const unsigned int num_elem3,
		const unsigned int num_elem4, const unsigned int num_elem5, const bool silent = false )
{
	array_pointer.resize( num_elem1 );
	for ( unsigned int i = 0; i < num_elem1; i++ )
	{
		make_array4d( array_pointer[i], num_elem2,
				num_elem3, num_elem4, num_elem5 );
	}
}
#endif // Ending functions


} // end namespace brgastro

#endif // _BRG_UTILITY_HPP_INCLUDED_
