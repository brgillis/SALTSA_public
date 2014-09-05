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
inline void set_zero( std::string obj )
{
	obj = "";
}
template< typename T >
inline void set_zero( std::vector< T > vec )
{
	vec.clear();
}
template< typename T >
inline void set_zero( T *obj )
{
	obj = NULL;
}
template< typename obj_type >
inline void set_zero( obj_type obj )
{
	obj = obj_type();
}

// Various "make" functions, to allocate dynamic memory.
// After allocating memory, these functions initialize the new variables using the
// set_zero function (see above).

template< typename obj_type >
inline void make_obj( BRG_UNIQUE_PTR<obj_type> & obj_pointer )
{
	obj_pointer = BRG_UNIQUE_PTR<obj_type>(new obj_type);
	set_zero(*obj_pointer);
}

#ifdef _BRG_USE_CPP_11_STD_
template <typename array_type>
inline void make_array1d( std::unique_ptr<array_type []> & array_pointer, const size_t num_elem )
{
	array_pointer = std::unique_ptr<array_type []>(new array_type [num_elem]);
	for(int i=0; i<num_elem; i++) set_zero(array_pointer[i]);
}
template <typename array_type>
inline void make_array( std::unique_ptr<array_type []> & array_pointer, const size_t num_elem )
{
	make_array1d( array_pointer, num_elem );
}
#endif
template< typename array_type >
inline void make_array1d( array_type & array_pointer,
		const size_t num_elem )
{
	array_pointer.resize( num_elem );
	for ( typename array_type::iterator it=array_pointer.begin(); it != array_pointer.end(); ++it)
		set_zero( *it );
}
template< typename array_type >
inline void make_array( array_type & array_pointer,
		const int num_elem )
{
	return make_array1d( array_pointer, num_elem );
}
template< typename array_type, typename other_array_type >
inline void make_array1d( array_type & array_pointer,
		const other_array_type & other_array )
{
	array_pointer.resize( other_array.size() );
	for ( typename array_type::iterator it=array_pointer.begin(); it != array_pointer.end(); ++it)
		set_zero( *it );
}
template< typename array_type, typename other_array_type >
inline void make_array( array_type & array_pointer,
		const other_array_type & other_array )
{
	return make_array1d( array_pointer, other_array );
}

#ifdef _BRG_USE_CPP_11_STD_
template <typename array_type>
inline void make_array2d( std::unique_ptr<std::unique_ptr<array_type []> []> & array_pointer, const size_t num_elem1, const size_t num_elem2, const bool silent=false )
{
	array_pointer = std::unique_ptr<std::unique_ptr<array_type []> []>(new std::unique_ptr<array_type []> [num_elem1]);
	for( size_t i = 0; i < num_elem1; i++)
	{
		make_array1d(array_pointer[i], num_elem2);
	}
}
#endif
template< typename array_type >
inline void make_array2d(
		std::vector< std::vector< array_type > > & array_pointer,
		const size_t num_elem1, const size_t num_elem2, const bool silent = false )
{
	array_pointer.resize( num_elem1 );
	for ( size_t i = 0; i < num_elem1; i++ )
	{
		make_array1d( array_pointer[i], num_elem2 );
	}
}
template< typename array_type, typename other_array_type >
inline void make_array2d( array_type & array_pointer,
		const other_array_type & other_array )
{
	array_pointer.resize( other_array.size() );
	typename other_array_type::const_iterator it_o = other_array.begin();
	for ( typename array_type::iterator it=array_pointer.begin(); it != array_pointer.end(); ++it)
		make_array1d(*it,*it_o);
		++it_o;
}

#ifdef _BRG_USE_CPP_11_STD_
template <typename array_type>
inline void make_array3d( std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> & array_pointer, const size_t num_elem1,
		const size_t num_elem2, const size_t num_elem3, const bool silent=false )
{
	array_pointer = std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []>(new
			std::unique_ptr<std::unique_ptr<array_type []> []> [num_elem1]);

	for( size_t i = 0; i < num_elem1; i++)
	{
		make_array2d(array_pointer[i], num_elem2, num_elem3);
	}
}
#endif
template< typename array_type >
inline void make_array3d(
		std::vector< std::vector< std::vector< array_type > > > & array_pointer,
		const size_t num_elem1, const size_t num_elem2, const size_t num_elem3,
		const bool silent = false )
{
	array_pointer.resize( num_elem1 );
	for ( size_t i = 0; i < num_elem1; i++ )
	{
		make_array2d( array_pointer[i], num_elem2, num_elem3 );
	}
}
template< typename array_type, typename other_array_type >
inline void make_array3d( array_type & array_pointer,
		const other_array_type & other_array )
{
	array_pointer.resize( other_array.size() );
	typename other_array_type::const_iterator it_o = other_array.begin();
	for ( typename array_type::iterator it=array_pointer.begin(); it != array_pointer.end(); ++it)
		make_array2d(*it,*it_o);
		++it_o;
}

#ifdef _BRG_USE_CPP_11_STD_
template <typename array_type>
inline void make_array4d( std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> []> & array_pointer,
		const size_t num_elem1, const size_t num_elem2, const size_t num_elem3, const size_t num_elem4, const bool silent=false )
{
	array_pointer = std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> []>(new
			std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> [num_elem1]);

	for( size_t i = 0; i < num_elem1; i++)
	{
		make_array3d(array_pointer[i], num_elem2, num_elem3, num_elem4);
	}
}
#endif
template< typename array_type >
inline void make_array4d(
		std::vector< std::vector< std::vector< std::vector< array_type > > > > & array_pointer,
		const size_t num_elem1, const size_t num_elem2, const size_t num_elem3,
		const size_t num_elem4, const bool silent = false )
{
	array_pointer.resize( num_elem1 );
	for ( size_t i = 0; i < num_elem1; i++ )
	{
		make_array3d( array_pointer[i], num_elem2, num_elem3, num_elem4 );
	}
}
template< typename array_type, typename other_array_type >
inline void make_array4d( array_type & array_pointer,
		const other_array_type & other_array )
{
	array_pointer.resize( other_array.size() );
	typename other_array_type::const_iterator it_o = other_array.begin();
	for ( typename array_type::iterator it=array_pointer.begin(); it != array_pointer.end(); ++it)
		make_array3d(*it,*it_o);
		++it_o;
}

#ifdef _BRG_USE_CPP_11_STD_
template <typename array_type>
inline void make_array5d( std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> []> []>
		& array_pointer, const size_t num_elem1, const size_t num_elem2, const size_t num_elem3,
		const size_t num_elem4, const size_t num_elem5, const bool silent=false )
{
	array_pointer = std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> []> []>(new
			std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> []> [num_elem1]);

	for( int i = 0; i < num_elem1; i++)
	{
		make_array4d(array_pointer[i], num_elem2, num_elem3, num_elem4, num_elem5);
	}
}
#endif
template< typename array_type >
inline void make_array5d(
		std::vector<
				std::vector<
						std::vector< std::vector< std::vector< array_type > > > > > & array_pointer,
		const size_t num_elem1, const size_t num_elem2, const size_t num_elem3,
		const size_t num_elem4, const size_t num_elem5, const bool silent = false )
{
	array_pointer.resize( num_elem1 );
	for ( size_t i = 0; i < num_elem1; i++ )
	{
		make_array4d( array_pointer[i], num_elem2,
				num_elem3, num_elem4, num_elem5 );
	}
}
template< typename array_type, typename other_array_type >
inline void make_array5d( array_type & array_pointer,
		const other_array_type & other_array )
{
	array_pointer.resize( other_array.size() );
	typename other_array_type::const_iterator it_o = other_array.begin();
	for ( typename array_type::iterator it=array_pointer.begin(); it != array_pointer.end(); ++it)
		make_array4d(*it,*it_o);
		++it_o;
}
#endif // Ending functions


} // end namespace brgastro

#endif // _BRG_UTILITY_HPP_INCLUDED_
