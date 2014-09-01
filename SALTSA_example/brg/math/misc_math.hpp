/**********************************************************************\
  @file misc_math.hpp

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

// body file: random_functions.cpp

#ifndef _BRG_MISC_MATH_HPP_INCLUDED_
#define _BRG_MISC_MATH_HPP_INCLUDED_

#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <vector>

#include "brg/global.h"

namespace brgastro {

// Returns true if val is Not a Number - Personal implementation, to make sure it works for all types
template< typename T >
inline const bool isnan( T val )
{
	return ( val != val );
}

// Returns true if val is infinity - ''
template< typename T >
inline const bool isinf( T val )
{
	using std::fabs;
#ifdef _BRG_USE_UNITS_
	if(std::numeric_limits<T>::max()==0) return fabs( (double)val ) > std::numeric_limits<double>::max();
#endif
	return fabs( val ) > std::numeric_limits<T>::max();
}

// Returns true if val is NaN or Inf
template< typename T >
inline const bool isbad( T val )
{
	using std::isnan;
	using std::isinf; // Use ADL here
	return ( isnan( val ) || isinf( val ) );
}

// Returns true if val is neither NaN nor Inf
template< typename T >
inline const bool isgood( T val )
{
	return !isbad( val );
}

// Min/max: Functions to return the lower/higher of two values.
// The return type is a bit complicated if the variables are of
// different types, so the function goes for the most general
// type of number allowed (unit_obj if units are being used,
// otherwise double).
template< class T1, class T2 >
inline const T1 min( const T1 a, const T2 b )
{
	return ( a < (T1)b ? a : (T1)b );
}
template< class T1, class T2 >
inline const T1 max( const T1 a, const T2 b )
{
	return ( a < (T1)b ? (T1)b : a );
}
template<  class T1, class T2, class T3 >
inline const T2 bound( const T1 lower_bound, const T2 a, const T3 upper_bound)
{
	return min( max( a, lower_bound ) , upper_bound);
}

// The below two variants return by reference, in case you
// want to do something like min_ref( a, b) = 0 to set the
// lower of two values to zero. Note that the values must
// be the same type here.
template< class T >
inline T & min_ref( T &a, T &b )
{
	return ( a > b ? b : a );
}
template< class T >
inline T & max_ref( T &a, T &b )
{
	return ( a < b ? b : a );
}

// Returns true if a is evenly divisible by b
template< typename Ta, typename Tb >
inline const bool divisible( const Ta a, const Tb b )
{
	if(b==0) throw std::runtime_error("It's undefined whether or not value is divisible by zero.");
	return ( a % b == 0 );
}

// Rounds to nearest integer, favoring even
template< typename T >
const int round_int( const T value, const double epsilon=DBL_EPSILON )
{

	if ( value < 0.0 )
		return -round_int( -value, epsilon );

	double ipart;
	std::modf( value, &ipart );

	// If 'value' is exactly halfway between two integers
	if ( fabs( value - ( ipart + 0.5 ) ) < epsilon )
	{
		// If 'ipart' is even then return 'ipart'
		if ( std::fmod( ipart, 2.0 ) < epsilon )
			return (int)ipart;

		// Else return the nearest even integer
		return (int)std::ceil( ipart + 0.5 );
	}

	// Otherwise use the usual round to closest
	// (Either symmetric half-up or half-down will do0
	return (int)std::floor( value + 0.5 );
}

// Inline square
template< typename T >
const T square(T v1)
{
	return v1*v1;
}

// Inline cube
template< typename T >
const T cube(T v1)
{
	return v1*v1*v1;
}

// Inline quart
template< typename T >
const T quart(T v1)
{
	return square(v1)*square(v1);
}

// Inline inverse
template< typename T >
const double inverse(T v1)
{
	return 1./v1;
}
inline const long double inverse(long double v1)
{
	return 1./v1;
}

// Inline square
template< typename T >
const double inv_square(T v1)
{
	return inverse(square(v1));
}
inline const long double inv_square(long double v1)
{
	return inverse(square(v1));
}

// Inline cube
template< typename T >
const double inv_cube(T v1)
{
	return inverse(cube(v1));
}
inline const long double inv_cube(long double v1)
{
	return inverse(cube(v1));
}

// Inline quart
template< typename T >
const double inv_quart(T v1)
{
	return inverse(quart(v1));
}
inline const long double inv_quart(long double v1)
{
	return inverse(quart(v1));
}

// Integer power - use only when you know it'll be an integer, but not the specific value,
// and when it isn't likely to be too high

#if (1)
// This version is optimized so that it won't check for the p==0 and p<0 cases, and generally
// shouldn't be called directly
template< typename T >
const T _ipow( T v, int p )
{
	if(p==1) return v;
	T tmp = _ipow(v,p/2);
	if(p%2==0) return tmp*tmp;
	return v*tmp*tmp;
}
#endif

template< typename T >
const T ipow( T v, int p )
{
	if(p<0) return 1/_ipow(v,-p);
	if(p==0) return 1;
	if(p==1) return v;
	T tmp = _ipow(v,p/2);
	if(p%2==0) return tmp*tmp;
	return v*tmp*tmp;
}

// Returns 1 if a is positive, -1 if it is negative, and 0 if it is 0, NaN if it is NaN
template< class T >
inline const int sign( const T a )
{
	if ( ( a == 0 ) || ( isnan( a ) ) )
		return a;
	if ( a < 0 )
		return -1;
	else
		return 1;
}

// Add two or three values in quadrature.
template < typename T1, typename T2 >
inline const T1 quad_add( const T1 v1, const T2 v2)
{
	using std::sqrt;
	return sqrt( v1 * v1 + v2 * v2);
}
template < typename T1, typename T2, typename T3 >
inline const T1 quad_add( const T1 v1, const T2 v2,
		const T3 v3 )
{
	using std::sqrt;
	return sqrt( v1 * v1 + v2 * v2 + v3 * v3 );
}

// Subtract one value from another in quadrature
template < typename T1, typename T2 >
inline const T1 quad_sub( const T1 v1, const T2 v2 )
{
	using std::sqrt;
	return sqrt( v1 * v1 - v2 * v2 );
}

// Function to calculate the distance between two points in 2-dimensions
template < typename Tx1, typename Ty1, typename Tx2, typename Ty2 >
inline const Tx1 dist2d( const Tx1 x1, const Ty1 y1, const Tx2 x2,
		const Ty2 y2 )
{
	return quad_add( x2 - x1, y2 - y1 );
}
// Function to calculate the distance between a point and (0,0) in 2-dimensions
template < typename Tx1, typename Ty1 >
inline const Tx1 dist2d( const Tx1 x1, const Ty1 y1 )
{
	return quad_add( x1, y1 );
}

// Use the law of cosines to calculate hypotenuse length (lc=Law of Cosines)
template < typename Tx1, typename Ty1, typename Ta1 >
inline const Tx1 lc_add( const Tx1 x1, const Ty1 y1, const Ta1 a1 )
{
	using std::sqrt;
	return sqrt( x1 * x1 + y1 * y1 - 2 * x1 * y1 * cos( a1 ) );
}

// 3-D distance between two points
template < typename Tx1, typename Ty1, typename Tz1, typename Tx2, typename Ty2, typename Tz2 >
inline const Tx1 dist3d( const Tx1 x1, const Ty1 y1, const Tz1 z1,
		const Tx2 x2, const Ty2 y2, const Tz2 z2 )
{
	return quad_add( x2 - x1, y2 - y1, z2 - z1 );
}

// 3-D distance from (0,0,0)
template < typename Tx, typename Ty, typename Tz >
inline const Tx dist3d( const Tx x1, const Ty y1, const Tz z1 )
{
	return quad_add( x1, y1, z1 );
}

// Distance between two vectors, where the dimensions are weighted by vector c
template< typename Ta, typename Tb >
inline const Ta weighted_dist( std::vector< Ta > a,
		std::vector< Tb > b )
{
	using std::sqrt;
	Ta result = 0;
	assert(a.size()==b.size());
	for ( unsigned int i = 0; i < a.size(); i++ )
	{
		result += square( b[i] - a[i] );
	}
	result = sqrt( result );
	return result;
}
template< typename Ta, typename Tb, typename Tc >
inline const Ta weighted_dist( std::vector< Ta > a,
		std::vector< Tb > b, std::vector< Tc > c )
{
	using std::sqrt;
	Ta result = 0;
	assert(a.size()==b.size());
	assert(a.size()==c.size());
	for ( size_t i = 0; i < a.size(); i++ )
	{
		result += square ( ( b[i] - a[i] ) * c[i] );
	}
	result = sqrt( result );
	return result;
}

// Dot-product of two vectors in 3-D space
template< typename Tx1, typename Ty1, typename Tz1, typename Tx2, typename Ty2, typename Tz2 >
inline const Tx1 dot_product( const Tx1 x1, const Ty1 y1,
		const Tz1 z1, const Tx2 x2, const Ty2 y2, const Tz2 z2 )
{
	return x1 * x2 + y1 * y2 + z1 * z2;
}
template< typename T1, typename T2 >
inline const double dot_product( const std::vector< T1 > & a,
		const std::vector< T2 > & b )
{
	assert(a.size()==b.size());
	double result = 0;
	for ( size_t i = 0; i < a.size(); i++ )
	{
		result += a[i] * b[i];
	}
	return result;
}

} // namespace brgastro

#endif /* _BRG_MISC_MATH_HPP_INCLUDED_ */
