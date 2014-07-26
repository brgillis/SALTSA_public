/**********************************************************************
 brg_misc_functions.hpp
 -----------------

 This is a self-contained header file, automatically included along with
 brg_functions.h. This file contains various template and inline
 functions. Note that some functions declared as "inline" cannot actually
 be inlined. They are all declared as such so they'll have local scope,
 and this header can be included in multiple source files without linker
 errors.

 All functions in this file are declared in the namespace brgastro.

 \**********************************************************************/

#ifndef __BRG_MISC_FUNCTIONS_HPP_INCLUDED__
#define __BRG_MISC_FUNCTIONS_HPP_INCLUDED__

#include <iostream>
#include <vector>
#include <cstdlib>
#include <limits>

#include <boost/math/special_functions/erf.hpp>

#include "brg_global.h"

namespace brgastro
{

// Generic functions
#if (1)
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
	if(b==0) return false;
	return ( a % b == 0 );
}

// "Safe" functions - perform the operation specified, but will
// take necessary actions to ensure it won't crash the program
// if the argument is invalid (eg. taking square-root of a negative
// number).
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
		return pow( min_d, -1 );

	if ( a >= 0 )
		result = max( a, min_d );
	else
		result = a;

	if ( result == 0 )
		result = 1; // In case we're dealing with integers here.

	return result;
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
	return exp( -pow( x - mean, 2 ) / ( 2 * std_dev * std_dev ) )
			/ ( std_dev * sqrt( 2 * pi ) );
}

// Function to integrate a Gaussian from min to max
template< typename Tlo, typename Thi, typename Tmean, typename Tstddev >
inline const double Gaus_int( const Tlo min, const Thi max)
{
	return Gaus_int(min,max,0.,1.);
}
template< typename Tlo, typename Thi, typename Tmean>
inline const double Gaus_int( const Tlo min, const Thi max,
		const Tmean mean)
{
	return Gaus_int(min,max,mean,1.);
}
template< typename Tlo, typename Thi, typename Tmean, typename Tstddev >
inline const double Gaus_int( const Tlo min, const Thi max,
		const Tmean mean, const Tstddev std_dev )
{
	double klo = ( min - mean ) / std_dev;
	double khi = ( max - mean ) / std_dev;

	return fabs( boost::math::erf( khi ) - boost::math::erf( klo ) ) / 2;
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

// Use the law of cosines to calculate hypotenuse length (lc=Law of Cosines)
template < typename Tx1, typename Ty1, typename Ta1 >
inline const Tx1 lc_add( const Tx1 x1, const Ty1 y1, const Ta1 a1 )
{
	return safe_sqrt( x1 * x1 + y1 * y1 - 2 * x1 * y1 * cos( a1 ) );
}

// Like dist2d, but using corrections for spherical geometry
template< typename Tr1, typename Td1, typename Tr2, typename Td2 >
inline const Tr1 skydist2d( const Tr1 ra1, const Td1 dec1,
		const Tr2 ra2, const Td2 dec2 )
{
	return quad_add( ( ra2 - ra1 ) * std::cos( ( dec2 + dec1 ) / 2 ), dec2 - dec1 );
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

template< typename T >
const int round_int( const T value, const double epsilon=DBL_EPSILON )
{

	if ( value < 0.0 )
		return -round_int( -value, epsilon );

	double ipart;
	std::modf( value, &ipart );

	// If 'value' is exctly halfway between two integers
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
		result += std::pow( ( b[i] - a[i] ) , 2 );
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
		result += std::pow( ( b[i] - a[i] ) * c[i], 2 );
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
inline const double drand( double min, double max )
{

	return min + (max-min)*drand48();

} // double drand(double min, double max)

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

template< class obj_type >
inline const int make_obj( BRG_UNIQUE_PTR<obj_type> & obj_pointer, const bool silent=false )
{
	obj_pointer = NULL;
	obj_pointer = new (std::nothrow) obj_type;
	if( obj_pointer == 0 ) return memory_error(silent);
	set_zero(*obj_pointer);
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

#ifdef _BRG_USE_CPP_11_STD_
template <class array_type>
const int make_array1d( std::unique_ptr<array_type []> & array_pointer, const int num_elem )
{
	array_pointer = std::unique_ptr<array_type []>(new (std::nothrow) array_type [num_elem]);
	if( array_pointer == 0 ) return memory_error(silent);
	for(int i=0; i<num_elem; i++) set_zero(array_pointer[i]);
	return 0;
}
template <class array_type>
const int make_array( std::unique_ptr<array_type []> & array_pointer, const int num_elem )
{
	return make_array1d( array_pointer, num_elem );
}
#endif

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

#ifdef _BRG_USE_CPP_11_STD_
template <class array_type>
const int make_array2d( std::unique_ptr<std::unique_ptr<array_type []> []> & array_pointer, const int num_elem1, const int num_elem2, const bool silent=false )
{
	array_pointer = std::unique_ptr<std::unique_ptr<array_type []> []>(new (std::nothrow) std::unique_ptr<array_type []> [num_elem1]);
	if( array_pointer == 0 ) return memory_error(silent);
	for( int i = 0; i < num_elem1; i++)
	{
		if (int errcode = make_array(array_pointer[i], num_elem2)) return errcode+LOWER_LEVEL_ERROR;
	}
	return 0;
}
#endif

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

#ifdef _BRG_USE_CPP_11_STD_
template <class array_type>
const int make_array3d( std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> & array_pointer, const int num_elem1,
		const int num_elem2, const int num_elem3, const bool silent=false )
{
	array_pointer = std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []>(new (std::nothrow)
			std::unique_ptr<std::unique_ptr<array_type []> []> [num_elem1]);
	if( array_pointer == 0 ) return memory_error(silent);

	for( int i = 0; i < num_elem1; i++)
	{
		if (int errcode = make_array2d(array_pointer[i], num_elem2, num_elem3)) return errcode+LOWER_LEVEL_ERROR;
	}
	return 0;
}
#endif

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
#ifdef _BRG_USE_CPP_11_STD_
template <class array_type>
const int make_array4d( std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> []> & array_pointer,
		const int num_elem1, const int num_elem2, const int num_elem3, const int num_elem4, const bool silent=false )
{
	array_pointer = std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> []>(new (std::nothrow)
			std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> [num_elem1]);
	if( array_pointer == 0 ) return memory_error(silent);

	for( int i = 0; i < num_elem1; i++)
	{
		if (int errcode = make_array3d(array_pointer[i], num_elem2, num_elem3, num_elem4)) return errcode+LOWER_LEVEL_ERROR;
	}
	return 0;
}
#endif
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
#ifdef _BRG_USE_CPP_11_STD_
template <class array_type>
const int make_array5d( std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> []> []>
		& array_pointer, const int num_elem1, const int num_elem2, const int num_elem3,
		const int num_elem4, const int num_elem5, const bool silent=false )
{
	array_pointer = std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> []> []>(new (std::nothrow)
			std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> []> [num_elem1]);
	if( array_pointer == 0 ) return memory_error(silent);

	for( int i = 0; i < num_elem1; i++)
	{
		if (int errcode = make_array4d(array_pointer[i], num_elem2, num_elem3, num_elem4, num_elem5)) return errcode+LOWER_LEVEL_ERROR;
	}
	return 0;
}
#endif
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

template< class array_type >
inline const int del_array1d( array_type * & array_pointer,
		const int num_elem = 0, const bool silent = false )
{
	if ( array_pointer == 0 )
		return INVALID_ARGUMENTS_ERROR;
	delete[] array_pointer;
	array_pointer = NULL;
	return 0;
}

template< class array_type >
inline const int del_array2d( array_type ** & array_pointer,
		const int num_elem1, const int num_elem2 = 0,
		const bool silent = false )
{
	if ( array_pointer == 0 )
		return INVALID_ARGUMENTS_ERROR;
	for ( int i = 0; i < num_elem1; i++ )
	{
		del_array( array_pointer[i] );
	}
	delete[] array_pointer;
	array_pointer = NULL;

	return 0;
}

template< class array_type >
inline const int del_array3d( array_type *** & array_pointer,
		const int num_elem1, const int num_elem2, const int num_elem3 = 0,
		const bool silent = false )
{
	if ( array_pointer == 0 )
		return INVALID_ARGUMENTS_ERROR;
	for ( int i = 0; i < num_elem1; i++ )
	{
		if ( int errcode = del_array2d( array_pointer[i], num_elem2 ) )
			return errcode + LOWER_LEVEL_ERROR;
	}
	delete[] array_pointer;
	array_pointer = NULL;

	return 0;
}

template< class array_type >
inline const int del_array4d( array_type **** & array_pointer,
		const int num_elem1, const int num_elem2, const int num_elem3,
		const int num_elem4 = 0, const bool silent = false )
{
	if ( array_pointer == 0 )
		return INVALID_ARGUMENTS_ERROR;
	for ( int i = 0; i < num_elem1; i++ )
	{
		if ( int errcode = del_array3d( array_pointer[i], num_elem2,
				num_elem3 ) )
			return errcode + LOWER_LEVEL_ERROR;
	}
	delete[] array_pointer;
	array_pointer = NULL;

	return 0;
}

template< class array_type >
const int del_array5d( array_type ***** & array_pointer, const int num_elem1,
		const int num_elem2, const int num_elem3, const int num_elem4,
		const int num_elem5 = 0, const bool silent = false )
{
	if ( array_pointer == 0 )
		return INVALID_ARGUMENTS_ERROR;
	for ( int i = 0; i < num_elem1; i++ )
	{
		if ( int errcode = del_array4d( array_pointer[i], num_elem2, num_elem3,
				num_elem4 ) )
			return errcode + LOWER_LEVEL_ERROR;
	}
	delete[] array_pointer;
	array_pointer = NULL;

	return 0;
}
#endif // Ending functions


} // end namespace brgastro

#endif // __BRG_MISC_FUNCTIONS_HPP_INCLUDED__
