/**********************************************************************
 @file SALTSA_misc_functions.hpp
 -------------------------------

 This is a self-contained header file, automatically included along with
 brg_functions.h. This file contains various template and inline
 functions. Note that some functions declared as "inline" cannot actually
 be inlined. They are all declared as such so they'll have local scope,
 and this header can be included in multiple source files without linker
 errors.

 All functions in this file are declared in the namespace SALTSA.

 \**********************************************************************/

#ifndef __SALTSA_MISC_FUNCTIONS_HPP_INCLUDED__
#define __SALTSA_MISC_FUNCTIONS_HPP_INCLUDED__

#include <iostream>
#include <vector>
#include <cstdlib>
#include <limits>

#include "SALTSA_global.h"

namespace SALTSA
{

// Generic functions
#if (1)
/**
 *
 * @param silent
 * @return
 */
inline const int errorNOS( const bool silent = false )
{
	if ( !silent )
		std::cerr
				<< "ERROR-not-otherwise-specified: I'm surprised the code actually found an error here.\n"
				<< "So surprised I didn't bother writing a specific message for it. Sorry; you'll have to\n"
				<< "start debugging and see where it hits this message. Put a breakpoint on the function\n"
				<< "errorNOS() in the brg_functions.hpp file and step out from there.\n";
	return UNSPECIFIED_ERROR;
}

/**
 *
 * @param silent
 * @return
 */
inline const int memory_error( const bool silent = false )
{
	if ( !silent )
		std::cerr << "ERROR: Could not assign sufficient dynamic memory.\n";
	return MEMORY_ERROR;
}

// Returns true if val is Not a Number - Personal implementation, to make sure it's included
template< typename T >
inline const bool isnan( T val )
{
	return ( val != val );
}

// Returns true if val is infinity - ''
template< typename T >
inline const bool isinf( T val )
{
	return std::fabs( val ) > std::numeric_limits<T>::max();
}

// Returns true if val is NaN or Inf
template< typename T >
inline const bool isbad( T val )
{
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
	return ( a < b ? a : (T1)b );
}
template< class T1, class T2 >
inline const T1 max( const T1 a, const T2 b )
{
	return ( a < b ? (T1)b : a );
}
template<  class T1, class T2, class T3 >
inline const T2 bound( const T1 lower_bound, const T2 a, const T3 upper_bound)
{
	return min( max( a, lower_bound ) , upper_bound);
}

// Returns true if a is evenly divisible by b
template< typename Ta, typename Tb >
inline const bool divisible( const Ta a, const Tb b )
{
	if(b==0) return false;
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
		return (int)ceil( ipart + 0.5 );
	}

	// Otherwise use the usual round to closest
	// (Either symmetric half-up or half-down will do0
	return (int)floor( value + 0.5 );
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
template< typename T >
const T ipow( T v, unsigned int p )
{
	if(p>=2) return v*ipow(v,p-1);
	if(p==1) return v;
	if(p<0) return 1/ipow(v,-p);
	return 1;
}

// "Safe" functions - perform the operation specified, but will
// take necessary actions to ensure it won't crash the program
// if the argument is invalid (eg. taking square-root of a negative
// number).
/**
 *
 * @param a
 * @return
 */
template< class T >
const T safe_sqrt( const T a )
{

#ifdef _BRG_WARN_FOR_SAFE_FUNCTIONS_TRIGGERED_
	if(a < 0)
	{
		std::cerr << "WARNING: safe_sqrt() prevented crash from negative input.\n";
	}
#endif
	return sqrt( std::fabs( a ) );
}
/**
 *
 * @param a
 * @return
 */
inline const double safe_sqrt( const int a ) // Special case for integers due to -INT_MIN > INT_MAX
{
	int res;

#ifdef _BRG_WARN_FOR_SAFE_FUNCTIONS_TRIGGERED_
	if(a < 0)
	{
		std::cerr << "WARNING: safe_sqrt() prevented crash from negative input.\n";
	}
#endif

	if ( a == std::numeric_limits<int>::min() )
	{
		res = std::numeric_limits<int>::max();
	}
	else
	{
		res = a < 0 ? -a : a;
	}
	return sqrt( res );
}
template< class Ta, class Tx >
inline const Ta safe_pow( const Ta a, const Tx x )
{
	Ta res;
	double ipart;

	std::modf( a, &ipart );

	if ( ( a < 0 ) && ( ipart != a ) )
	{
#ifdef _BRG_WARN_FOR_SAFE_FUNCTIONS_TRIGGERED_
		std::cerr << "WARNING: safe_pow() prevented error from negative input.\n";
#endif
		res = -a;
	}
	else
	{
		res = a;
	}

	return std::pow( res, x );
}

// Safe_d is used a bit differently: Put it around the denominator to make
// sure you don't hit a divide-by-zero error.
/**
 *
 * @param a
 * @return
 */
template< class T >
const T safe_d( const T a )
{

#ifdef _BRG_WARN_FOR_SAFE_FUNCTIONS_TRIGGERED_
	if( (a == 0) || isbad(a) )
	{
		std::cerr << "WARNING: safe_d() prevented error from zero input or bad input.\n";
	}
#endif
	T result = a;
	T min_d = a; // So it'll have the right units if we're using units here.
	min_d = MIN_DIVISOR;

	if ( isnan( a ) )
		return min_d;
	if ( isinf( a ) )
		return inverse( min_d );

	if (std::fabs(a)>min_d)
		return a;

	if(min_d == 0) return 1; // In case of integers

	return min_d;

}

// Gaussian PDF
template< typename Tx >
inline const double Gaus_pdf( const Tx x )
{
	return Gaus_pdf(x,0.,1.);
}
template< typename Tx, typename Tmean >
inline const double Gaus_pdf( const Tx x, const Tmean mean )
{
	return Gaus_pdf(x,mean,1.);
}
template< typename Tx, typename Tmean, typename Tstddev >
inline const double Gaus_pdf( const Tx x, const Tmean mean = 0,
		const Tstddev std_dev = 1 )
{
	return exp( -square( x - mean ) / ( 2 * std_dev * std_dev ) )
			/ ( std_dev * sqrt( 2 * pi ) );
}

// Add two or three values in quadrature.
template < typename T1, typename T2 >
inline const T1 quad_add( const T1 v1, const T2 v2)
{
	return sqrt( v1 * v1 + v2 * v2);
}
template < typename T1, typename T2, typename T3 >
inline const T1 quad_add( const T1 v1, const T2 v2,
		const T3 v3 )
{
	return sqrt( v1 * v1 + v2 * v2 + v3 * v3 );
}

// Subtract one value from another in quadrature
template < typename T1, typename T2 >
inline const T1 quad_sub( const T1 v1, const T2 v2 )
{
	return safe_sqrt( v1 * v1 - v2 * v2 );
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
// An exception is thrown if the vectors aren't of equal sizes
template< typename Ta, typename Tb >
inline const Ta weighted_dist( std::vector< Ta > a,
		std::vector< Tb > b )
{
	Ta result = 0;
	if ( a.size() != b.size() )
	{
		throw std::runtime_error("ERROR: Vectors in weighted_dist have unequal sizes.\n");
	}
	for ( unsigned int i = 0; i < a.size(); i++ )
	{
		result += square( b[i] - a[i] );
	}
	result = safe_sqrt( result );
	return result;
}
template< typename Ta, typename Tb, typename Tc >
inline const Ta weighted_dist( std::vector< Ta > a,
		std::vector< Tb > b, std::vector< Tc > c )
{
	Ta result = 0;
	if ( ( a.size() != b.size() ) || ( a.size() != c.size() ) )
	{
		throw std::runtime_error("ERROR: Vectors in weighted_dist have unequal sizes.\n");
	}
	for ( size_t i = 0; i < a.size(); i++ )
	{
		result += square ( ( b[i] - a[i] ) * c[i] );
	}
	result = safe_sqrt( result );
	return result;
}

// Dot-product of two vectors in 3-D space or two vectors
// For two vectors, an exception is thrown if they aren't the same size
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
	if ( ( a.size() != b.size() ) )
	{
		throw std::runtime_error("ERROR: Vectors in dot_product have unequal sizes.\n");
	}
	double result = 0;
	for ( size_t i = 0; i < a.size(); i++ )
	{
		result += a[i] * b[i];
	}
	return result;
}

// Generates a random double between min and max
/**
 *
 * @param min
 * @param max
 * @return
 */
inline const double drand( double min, double max )
{

	return min + (max-min)*drand48();

} // double drand(double min, double max)

// Returns 1 if a is positive, -1 if it is negative, and 0 if it is 0, NaN if it is NaN
/**
 *
 * @param a
 * @return
 */
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

// Set_zero function - a way for other template functions to "clear" or initialize a value in various ways
// Types of variables for which the method is defined will return 0, otherwise it will do nothing and
// return 1;
inline const int set_zero( int & obj )
{
	obj = 0;
	return 0;
}
inline const int set_zero( long int & obj )
{
	obj = 0;
	return 0;
}
inline const int set_zero( short int & obj )
{
	obj = 0;
	return 0;
}
inline const int set_zero( unsigned int & obj )
{
	obj = 0;
	return 0;
}
inline const int set_zero( unsigned long int & obj )
{
	obj = 0;
	return 0;
}
inline const int set_zero( unsigned short int & obj )
{
	obj = 0;
	return 0;
}
inline const int set_zero( double & obj )
{
	obj = 0;
	return 0;
}
inline const int set_zero( long double & obj )
{
	obj = 0;
	return 0;
}
inline const int set_zero( float & obj )
{
	obj = 0;
	return 0;
}
#ifdef _BRG_USE_UNITS_
inline const int set_zero( unit_obj obj)
{
	obj = 0;
	return 0;
}
#endif
inline const int set_zero( std::string obj )
{
	obj = "";
	return 0;
}
/**
 *
 * @param vec
 * @return
 */
template< class T >
inline const int set_zero( std::vector< T > vec )
{
	return vec.clear();
}
template< class T >
inline const int set_zero( T *obj )
{
	obj = NULL;
	return 0;
}
template< class obj_type >
inline const int set_zero( obj_type obj )
{
	return INVALID_ARGUMENTS_ERROR;
}

// Various "make" functions, to allocate dynamic memory.
// These functions return 0 on success, 1 on failure (not enough memory available).
// After allocating memory, these functions initialize the new variables using the
// set_zero function (see above).
// Remember that if you use normal pointers, you have to delete the assigned variables
// before they go out of scope. See the del_obj and del_array functions below for that

template< class obj_type >
const int make_obj( obj_type * & obj_pointer, const bool silent = false )
{
	obj_pointer = NULL;
	obj_pointer = new ( std::nothrow ) obj_type;
	if ( obj_pointer == 0 )
		return memory_error( silent );
	set_zero( *obj_pointer );
	return 0;
}

template< class array_type >
const int make_array1d( array_type * & array_pointer, const int num_elem,
		const bool silent = false )
{
	array_pointer = new ( std::nothrow ) array_type[num_elem];
	if ( array_pointer == 0 )
		return memory_error( silent );
	for ( int i = 0; i < num_elem; i++ )
		set_zero( array_pointer[i] );
	return 0;
}
template< class array_type >
const int make_array( array_type * & array_pointer, const int num_elem,
		const bool silent = false )
{
	return make_array1d( array_pointer, num_elem );
}

template< class array_type >
const int make_array1d( std::vector< array_type > & array_pointer,
		const int num_elem, const bool silent = false )
{
	array_pointer.clear();
	array_pointer.resize( num_elem );
	if ( array_pointer.empty() )
		return memory_error( silent );
	for ( int i = 0; i < num_elem; i++ )
		set_zero( array_pointer[i] );
	return 0;
}
template< class array_type >
const int make_array( std::vector< array_type > & array_pointer,
		const int num_elem, const bool silent = false )
{
	return make_array1d( array_pointer, num_elem );
}

template< class array_type >
const int make_array2d( array_type ** & array_pointer, const int num_elem1,
		const int num_elem2, const bool silent = false )
{
	array_pointer = new ( std::nothrow ) array_type *[num_elem1];
	if ( array_pointer == 0 )
		return memory_error( silent );
	for ( int i = 0; i < num_elem1; i++ )
	{
		if ( int errcode = make_array( array_pointer[i], num_elem2 ) )
			return errcode + LOWER_LEVEL_ERROR;
	}
	return 0;
}

template< class array_type >
const int make_array2d(
		std::vector< std::vector< array_type > > & array_pointer,
		const int num_elem1, const int num_elem2, const bool silent = false )
{
	array_pointer.clear();
	array_pointer.resize( num_elem1 );
	if ( array_pointer.empty() )
		return memory_error( silent );
	for ( int i = 0; i < num_elem1; i++ )
	{
		if ( int errcode = make_array( array_pointer[i], num_elem2 ) )
			return errcode + LOWER_LEVEL_ERROR;
	}
	return 0;
}

template< class array_type >
const int make_array3d( array_type *** & array_pointer, const int num_elem1,
		const int num_elem2, const int num_elem3, const bool silent = false )
{
	array_pointer = new ( std::nothrow ) array_type **[num_elem1];
	if ( array_pointer == 0 )
		return memory_error( silent );

	for ( int i = 0; i < num_elem1; i++ )
	{
		if ( int errcode = make_array2d( array_pointer[i], num_elem2,
				num_elem3 ) )
			return errcode + LOWER_LEVEL_ERROR;
	}
	return 0;
}

template< class array_type >
const int make_array3d(
		std::vector< std::vector< std::vector< array_type > > > & array_pointer,
		const int num_elem1, const int num_elem2, const int num_elem3,
		const bool silent = false )
{
	array_pointer.clear();
	array_pointer.resize( num_elem1 );
	if ( array_pointer.empty() )
		return memory_error( silent );
	for ( int i = 0; i < num_elem1; i++ )
	{
		if ( int errcode = make_array2d( array_pointer[i], num_elem2,
				num_elem3 ) )
			return errcode + LOWER_LEVEL_ERROR;
	}
	return 0;
}

template< class array_type >
const int make_array4d( array_type **** & array_pointer, const int num_elem1,
		const int num_elem2, const int num_elem3, const int num_elem4,
		const bool silent = false )
{
	array_pointer = new ( std::nothrow ) array_type ***[num_elem1];
	if ( array_pointer == 0 )
		return memory_error( silent );

	for ( int i = 0; i < num_elem1; i++ )
	{
		if ( int errcode = make_array3d( array_pointer[i], num_elem2,
				num_elem3, num_elem4 ) )
			return errcode + LOWER_LEVEL_ERROR;
	}
	return 0;
}

template< class array_type >
const int make_array4d(
		std::vector< std::vector< std::vector< std::vector< array_type > > > > & array_pointer,
		const int num_elem1, const int num_elem2, const int num_elem3,
		const int num_elem4, const bool silent = false )
{
	array_pointer.clear();
	array_pointer.resize( num_elem1 );
	if ( array_pointer.empty() )
		return memory_error( silent );
	for ( int i = 0; i < num_elem1; i++ )
	{
		if ( int errcode = make_array3d( array_pointer[i], num_elem2,
				num_elem3, num_elem4 ) )
			return errcode + LOWER_LEVEL_ERROR;
	}
	return 0;
}

template< class array_type >
const int make_array5d( array_type ***** & array_pointer, const int num_elem1,
		const int num_elem2, const int num_elem3, const int num_elem4,
		const int num_elem5, const bool silent = false )
{
	array_pointer = new ( std::nothrow ) array_type ****[num_elem1];
	if ( array_pointer == 0 )
		return memory_error( silent );

	for ( int i = 0; i < num_elem1; i++ )
	{
		if ( int errcode = make_array4d( array_pointer[i], num_elem2,
				num_elem3, num_elem4, num_elem5 ) )
			return errcode + LOWER_LEVEL_ERROR;
	}
	return 0;
}

template< class array_type >
const int make_array5d(
		std::vector<
				std::vector<
						std::vector< std::vector< std::vector< array_type > > > > > & array_pointer,
		const int num_elem1, const int num_elem2, const int num_elem3,
		const int num_elem4, const int num_elem5, const bool silent = false )
{
	array_pointer.clear();
	array_pointer.resize( num_elem1 );
	if ( array_pointer.empty() )
		return memory_error( silent );
	for ( int i = 0; i < num_elem1; i++ )
	{
		if ( int errcode = make_array4d( array_pointer[i], num_elem2,
				num_elem3, num_elem4, num_elem5 ) )
			return errcode + LOWER_LEVEL_ERROR;
	}
	return 0;
}

// Delete functions to simply delete dynamically assigned objects, arrays, and multiple-dimension arrays

template< class obj_type >
inline const int del_obj( obj_type * & obj_pointer, const bool silent = false ) throw()
{
	if ( obj_pointer == NULL ) return UNSPECIFIED_ERROR;
	delete obj_pointer;
	obj_pointer = NULL;
	return 0;
}

#endif // Ending functions


} // end namespace SALTSA

#endif // __SALTSA_MISC_FUNCTIONS_HPP_INCLUDED__
