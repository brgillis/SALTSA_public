#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <new>
#include <fstream>
#include <string>
#include "SALTSA_misc_functions.h"
#include "SALTSA_calculus.hpp"

using namespace std;

/** Class method definitions **/
#if (1)
// SALTSA::phase function implementations
#if (1)
/**
 *
 * @param init_x
 * @param init_y
 * @param init_z
 * @param init_vx
 * @param init_vy
 * @param init_vz
 * @param init_t
 */
SALTSA::phase::phase( double init_x, double init_y,
double init_z,
double init_vx, double init_vy, double init_vz,
double init_t )
{
	x = init_x;
	y = init_y;
	z = init_z;
	vx = init_vx;
	vy = init_vy;
	vz = init_vz;
	t = init_t;
}

/**
 *
 * @param init_x
 * @param init_y
 * @param init_z
 * @param init_vx
 * @param init_vy
 * @param init_vz
 * @param init_t
 * @return
 */
const int SALTSA::phase::set_phase( double init_x, double init_y,
double init_z,
double init_vx, double init_vy, double init_vz,
double init_t )
{
	x = init_x;
	y = init_y;
	z = init_z;
	vx = init_vx;
	vy = init_vy;
	vz = init_vz;
	t = init_t;

	return 0;
}

#endif // end SALTSA::phase function implementations

#endif // end class function definitions

/** Global function implementations **/
#if (1)
/**
 *
 * @param value
 * @param epsilon
 * @return
 */
const int SALTSA::round_int( const double value, const double epsilon )
{

	if ( value < 0.0 )
		return -SALTSA::round_int( -value, epsilon );

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



/**
 *
 * @param mean
 * @param stddev
 * @return
 */
const double SALTSA::Gaus_rand( double mean, double stddev )
{
	double x1, x2, w;

	if ( stddev <= 0 )
		return mean;

	do
	{
		x1 = 2.0 * ( ( (double)rand() ) / RAND_MAX ) - 1.0;
		x2 = 2.0 * ( ( (double)rand() ) / RAND_MAX ) - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = sqrt( ( -2.0 * log( w ) ) / w );
	return ( mean + x1 * w * stddev );

} // double Gaus_rand(double mean, double stddev)

/**
 *
 * @param mean
 * @param stddev
 * @return
 */
const double SALTSA::log10Gaus_rand( double mean, double stddev )
{
	double x1, x2, w, fact;

	if ( stddev <= 0 )
		return mean;

	stddev *= log( 10 ); // Converts dex to natural log
	fact = exp( -std::pow( stddev, 2 ) / 2 );

	do
	{
		x1 = 2.0 * ( ( (double)rand() ) / RAND_MAX ) - 1.0;
		x2 = 2.0 * ( ( (double)rand() ) / RAND_MAX ) - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = sqrt( ( -2.0 * log( w ) ) / w );
	return ( mean * fact * exp( x1 * w * stddev ) );
} // double lnGaus_rand(double mean, double stddev)

/**
 *
 * @param lambda
 * @return
 */
const int SALTSA::Pois_rand( double lambda )
{
	double L, p;
	int k;
	L = exp( -lambda );
	k = 0;
	p = 1;
	do
	{
		k++;
		p *= drand( 0, 1 );
	} while ( p > L );

	return ( k - 1 );
} // int Pois_rand(double lambda)

#endif // end global function implementations
