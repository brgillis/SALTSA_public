/**********************************************************************
 brg_functions.hpp
 -----------------

 This is a self-contained header file, automatically included along with
 brg_functions.h. This file contains various template and inline
 functions. Note that some functions declared as "inline" cannot actually
 be inlined. They are all declared as such so they'll have local scope,
 and this header can be included in multiple source files without linker
 errors.

 All functions in this file are declared in the namespace brgastro.

 \**********************************************************************/

#ifndef __BRG_FUNCTIONS_HPP_INCLUDED__
#define __BRG_FUNCTIONS_HPP_INCLUDED__

#include "brg_global.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <new>
#include <fstream>
#include <memory>
#include <boost/math/special_functions/erf.hpp>
#ifndef __USE_CPP_11_STD__
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#endif // #ifndef __USE_CPP_11_STD__

#ifdef _BRG_USE_UNITS_
#include "brg_units.h"
#endif

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
inline const bool isnan( double val )
{
	return ( val != val );
}

// Returns true if val is infinity - ''
inline const bool isinf( double val )
{
	return fabs( val ) > DBL_MAX;
}

// Returns true if val is NaN or Inf
inline const bool isbad( double val )
{
	return ( isnan( val ) || isinf( val ) );
}

// Returns true if val is neither NaN nor Inf
inline const bool isgood( double val )
{
	return !isbad( val );
}

// Min/max: Functions to return the lower/higher of two values.
// The return type is a bit complicated if the variables are of
// different types, so the function goes for the most general
// type of number allowed (unit_obj if units are being used,
// otherwise double).
template< class T >
inline const T min( const T &a, const T &b )
{
	return ( a < b ? a : b );
}
template< class T >
inline const T max( const T &a, const T &b )
{
	return ( a < b ? b : a );
}
template< class T >
inline const T bound( const T &lower_bound, const T &a, const T &upper_bound)
{
	return min( max( lower_bound, a ) , upper_bound);
}
#ifdef _BRG_USE_UNITS_
template<class T1, class T2>
inline const brgastro::unit_obj min(const T1 a, const T2 b)
{
	return( (a<b) ? brgastro::unit_obj(a) : brgastro::unit_obj(b) );
}
template<class T1, class T2>
inline const brgastro::unit_obj max(const T1 a, const T2 b)
{
	return( (a<b) ? brgastro::unit_obj(b) : brgastro::unit_obj(a) );
}
template< class T1, class T2, class T3 >
inline const brgastro::unit_obj bound( const T1 &lower_bound, const T2 &a, const T3 &upper_bound)
{
	return min( max( lower_bound, a ) , upper_bound);
}
#else
template< class T1, class T2 >
inline const double min( const T1 a, const T2 b )
{
	return ( a < b ? (double)a : (double)b );
}
template< class T1, class T2 >
inline const double max( const T1 a, const T2 b )
{
	return ( a < b ? (double)b : (double)a );
}
#endif

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
inline const bool divisible( const int a, const int b )
{
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

	if ( a == INT_MIN )
	{
		res = INT_MAX;
	}
	else
	{
		res = a < 0 ? -a : a;
	}
	return sqrt( res );
}
inline const double safe_pow( const int a, const double x )
{
	int res;
	double ipart;

	std::modf( a, &ipart );

	if ( ( a < 0 ) && ( ipart != a ) )
	{
#ifdef _BRG_WARN_FOR_SAFE_FUNCTIONS_TRIGGERED_
		std::cerr << "WARNING: safe_pow() prevented error from negative input.\n";
#endif

		if ( INT_MIN == a )
		{
			res = INT_MAX;
		}
		else
		{
			res = -a;
		}
	}
	else
	{
		res = a;
	}

	return std::pow( res, x );
}
template< class T >
inline const T safe_pow( const T a, const double x )
{
	T res;
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
inline const double Gaus_pdf( const double x, const double mean = 0,
		const double std_dev = 1 )
{
	return exp( -pow( x - mean, 2 ) / ( 2 * std_dev * std_dev ) )
			/ ( std_dev * sqrt( 2 * pi ) );
}

// Function to integrate a Gaussian from min to max
inline const double Gaus_int( const double min, const double max,
		const double mean = 0, const double std_dev = 1 )
{
	double klo = ( min - mean ) / std_dev;
	double khi = ( max - mean ) / std_dev;

	return fabs( boost::math::erf( khi ) - boost::math::erf( klo ) ) / 2;
}

// Add two or three values in quadrature.
inline const double quad_add( const double v1, const double v2,
		const double v3 = 0 )
{
	return sqrt( v1 * v1 + v2 * v2 + v3 * v3 );
}

// Subtract one value from another in quadrature
inline const double quad_sub( const double v1, const double v2 )
{
	return safe_sqrt( v1 * v1 - v2 * v2 );
}

// unit_obj versions of the above
#ifdef _BRG_USE_UNITS_
inline const unit_obj quad_add( const unit_obj v1, const unit_obj v2, const unit_obj v3=0 )
{
	return sqrt( v1*v1 + v2*v2 + v3*v3 );
}
inline const unit_obj quad_sub( const unit_obj v1, const unit_obj v2 )
{
	return safe_sqrt( v1*v1 - v2*v2 );
}
#endif

// Function to calculate the distance between two points in 2-dimensions
#ifdef _BRG_USE_UNITS_
inline const unit_obj dist2d( const unit_obj x1, const unit_obj y1, const unit_obj x2, const unit_obj y2 )
#else
inline const double dist2d( const double x1, const double y1, const double x2,
		const double y2 )
#endif
{
	return quad_add( x2 - x1, y2 - y1 );
}
// Function to calculate the distance between a point and (0,0) in 2-dimensions
#ifdef _BRG_USE_UNITS_
inline const unit_obj dist2d( const unit_obj x1, const unit_obj y1 )
#else
inline const double dist2d( const double x1, const double y1 )
#endif
{
	return quad_add( x1, y1 );
}

// Use the law of cosines to calculate hypotenuse length (lc=Law of Cosines)
#ifdef _BRG_USE_UNITS_
inline const unit_obj lc_add( const unit_obj x1, const unit_obj y1, const unit_obj a1 )
#else
inline const double lc_add( const double x1, const double y1, const double a1 )
#endif
{
	return safe_sqrt( x1 * x1 + y1 * y1 - 2 * x1 * y1 * cos( a1 ) );
}

// Like dist2d, but using corrections for spherical geometry
#ifdef _BRG_USE_UNITS_
inline const unit_angle skydist2d( const unit_angle ra1, const unit_angle dec1, const unit_angle ra2, const unit_angle dec2 )
#else
inline const double skydist2d( const double ra1, const double dec1,
		const double ra2, const double dec2 )
#endif
{
	return quad_add( ( ra2 - ra1 ) * cos( ( dec2 + dec1 ) / 2 ), dec2 - dec1 );
}

// 3-D distance between two points
#ifdef _BRG_USE_UNITS_
inline const unit_obj dist3d( const unit_obj x1, const unit_obj y1, const unit_obj z1, const unit_obj x2, const unit_obj y2, const unit_obj z2 )
#else
inline const double dist3d( const double x1, const double y1, const double z1,
		const double x2, const double y2, const double z2 )
#endif
{
	return quad_add( x2 - x1, y2 - y1, z2 - z1 );
}

// 3-D distance from (0,0,0)
#ifdef _BRG_USE_UNITS_
inline const unit_obj dist3d( const unit_obj x1, const unit_obj y1, const unit_obj z1 )
#else
inline const double dist3d( const double x1, const double y1, const double z1 )
#endif
{
	return quad_add( x1, y1, z1 );
}

// Distance between two vectors, where the dimensions are weighted by vector c
// An exception is thrown if the vectors aren't of equal sizes
inline const double weighted_dist( std::vector< double > a,
		std::vector< double > b, std::vector< double > c =
				std::vector< double >( 0 ) )
{
	double result = 0;
	if ( c.size() == 0 )
		c.resize( a.size(), 1 );
	if ( ( a.size() != b.size() ) || ( a.size() != c.size() ) )
	{
		throw std::runtime_error("ERROR: Vectors in weighted_dist have unequal sizes.\n");
		return -1;
	}
	for ( unsigned int i = 0; i < a.size(); i++ )
	{
		result += std::pow( ( b.at( i ) - a.at( i ) ) * c.at( i ), 2 );
	}
	result = safe_sqrt( result );
	return result;
}
#ifdef _BRG_USE_UNITS_
inline const unit_obj weighted_dist(std::vector<unit_obj> a, std::vector<unit_obj> b, std::vector<double> c = std::vector<double>(0))
{
	unit_obj result = 0;
	if(c.size() == 0) c.resize(a.size());
	if((a.size() != b.size())||(a.size() != c.size()))
	{
		throw std::runtime_error("ERROR: Vectors in weighted_dist have unequal sizes.\n");
		return -1;
	}
	for(unsigned int i = 0; i < a.size(); i++)
	{
		result += pow((b.at(i)-a.at(i))*c.at(i),2);
	}
	result = safe_sqrt(result);
	return result;
}
#endif

// Dot-product of two vectors in 3-D space or two vectors
// For two vectors, an exception is thrown if they aren't the same size
#ifdef _BRG_USE_UNITS_
inline const unit_obj dot_product( const unit_obj x1, const unit_obj y1, const unit_obj z1, const unit_obj x2, const unit_obj y2, const unit_obj z2 )
{
	return x1*x2+y1*y2+z1*z2;
}
#else
inline const double dot_product( const double x1, const double y1,
		const double z1, const double x2, const double y2, const double z2 )
{
	return x1 * x2 + y1 * y2 + z1 * z2;
}
#endif
inline const double dot_product( std::vector< double > a,
		std::vector< double > b )
{
	if ( ( a.size() != b.size() ) )
	{
		throw std::runtime_error("ERROR: Vectors in dot_product have unequal sizes.\n");
		return 0;
	}
	double result = 0;
	for ( unsigned int i = 0; i < a.size(); i++ )
	{
		result += a.at( i ) * b.at( i );
	}
	return result;
}
#ifdef _BRG_USE_UNITS_
inline const unit_obj dot_product(std::vector<unit_obj> a, std::vector<unit_obj> b)
{
	if((a.size() != b.size()))
	{
		throw std::runtime_error("ERROR: Vectors in dot_product have unequal sizes.\n");
		return 0;
	}
	unit_obj result = 0;
	for(unsigned int i = 0; i < a.size(); i++)
	{
		result += a.at(i)*b.at(i);
	}
	return result;
}
#endif

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
inline const int set_zero( int obj )
{
	obj = 0;
	return 0;
}
inline const int set_zero( double obj )
{
	obj = 0;
	return 0;
}
inline const int set_zero( float obj )
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
// These aren't actually necessary for anything but normal pointers, but the other forms are included
// so pointers can be exchanged with vectors and unique/auto_ptrs with minimal issue.

template< class obj_type >
inline const int del_obj( obj_type * & obj_pointer, const bool silent = false )
{
	if ( obj_pointer == 0 )
		return UNSPECIFIED_ERROR;
	delete obj_pointer;
	obj_pointer = NULL;
	return 0;
}
template< class obj_type >
inline const int del_obj( BRG_UNIQUE_PTR<obj_type> & obj_pointer, const bool silent=false )
{
	obj_pointer.reset();
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
inline const int del_array( array_type * & array_pointer, const int num_elem =
		0, const bool silent = false )
{
	return del_array1d( array_pointer );
}
#ifdef _BRG_USE_CPP_11_STD_
template <class array_type>
inline const int del_array1d( std::unique_ptr<array_type []> & array_pointer, const int num_elem=0, const bool silent=false )
{
	if( array_pointer == 0 ) return INVALID_ARGUMENTS_ERROR;
	array_pointer = NULL;
	return 0;
}
template <class array_type>
inline const int del_array( std::unique_ptr<array_type []> & array_pointer, const int num_elem=0, const bool silent=false )
{
	return del_array1d(array_pointer);
}
#endif
template< class array_type >
inline const int del_array1d( std::vector< array_type > & array_pointer,
		const int num_elem = 0, const bool silent = false )
{
	array_pointer.clear();
	return 0;
}
template< class array_type >
inline const int del_array( std::vector< array_type > & array_pointer,
		const int num_elem = 0, const bool silent = false )
{
	return del_array1d( array_pointer );
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
#ifdef _BRG_USE_CPP_11_STD_
template <class array_type>
inline const int del_array2d( std::unique_ptr<std::unique_ptr<array_type []> []> & array_pointer, const int num_elem1, const int num_elem2=0, const bool silent=false )
{
	if( array_pointer == 0 ) return INVALID_ARGUMENTS_ERROR;
	for( int i = 0; i < num_elem1; i++)
	{
		if(int errcode = del_array(array_pointer[i])) return errcode+LOWER_LEVEL_ERROR;
	}
	array_pointer = NULL;

	return 0;
}
#endif
template< class array_type >
inline const int del_array2d(
		std::vector< std::vector< array_type > > & array_pointer,
		const int num_elem1, const int num_elem2 = 0,
		const bool silent = false )
{
	array_pointer.clear();

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
#ifdef _BRG_USE_CPP_11_STD_
template <class array_type>
inline const int del_array3d( std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> & array_pointer,
		const int num_elem1, const int num_elem2, const int num_elem3=0, const bool silent=false )
{
	if( array_pointer == 0 ) return INVALID_ARGUMENTS_ERROR;
	for( int i = 0; i < num_elem1; i++)
	{
		if(int errcode = del_array2d(array_pointer[i], num_elem2)) return errcode+LOWER_LEVEL_ERROR;
	}
	array_pointer = NULL;

	return 0;
}
#endif
template< class array_type >
inline const int del_array3d(
		std::vector< std::vector< std::vector< array_type > > > & array_pointer,
		const int num_elem1, const int num_elem2, const int num_elem3 = 0,
		const bool silent = false )
{
	array_pointer.clear();

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
#ifdef _BRG_USE_CPP_11_STD_
template <class array_type>
const int del_array4d( std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []> []> []> []>
		& array_pointer, const int num_elem1, const int num_elem2, const int num_elem3,
		const int num_elem4=0, const bool silent=false )
{
	if( array_pointer == 0 ) return INVALID_ARGUMENTS_ERROR;
	for( int i = 0; i < num_elem1; i++)
	{
		if(int errcode = del_array3d(array_pointer[i], num_elem2, num_elem3)) return errcode+LOWER_LEVEL_ERROR;
	}
	array_pointer = NULL;

	return 0;
}
#endif
template< class array_type >
const int del_array4d(
		std::vector< std::vector< std::vector< std::vector< array_type > > > > & array_pointer,
		const int num_elem1, const int num_elem2, const int num_elem3,
		const int num_elem4 = 0, const bool silent = false )
{
	array_pointer.clear();

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
#ifdef _BRG_USE_CPP_11_STD_
template <class array_type>
const int del_array5d( std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<std::unique_ptr<array_type []>
		[]> []> []> []> & array_pointer, const int num_elem1, const int num_elem2, const int num_elem3,
		const int num_elem4, const int num_elem5=0, const bool silent=false )
{
	if( array_pointer == 0 ) return INVALID_ARGUMENTS_ERROR;
	for( int i = 0; i < num_elem1; i++)
	{
		if(int errcode = del_array4d(array_pointer[i], num_elem2, num_elem3, num_elem4)) return errcode+LOWER_LEVEL_ERROR;
	}
	array_pointer = NULL;

	return 0;
}
#endif
template< class array_type >
const int del_array5d(
		std::vector<
				std::vector<
						std::vector< std::vector< std::vector< array_type > > > > > & array_pointer,
		const int num_elem1, const int num_elem2, const int num_elem3,
		const int num_elem4, const int num_elem5 = 0,
		const bool silent = false )
{
	array_pointer.clear();

	return 0;
}
#endif // Ending functions


} // end namespace brgastro

#endif
