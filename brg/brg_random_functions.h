/**********************************************************************\
brg_random_functions.h
 -----------

 If this header is used, the source file brg_functions.cpp must be included
 and compiled with the project. This file automatically includes
 brg_functions.hpp, which contains all template and inline functions.
 More complex functions are declared in this file and implemented in
 brg_functions.cpp.

 This file includes various classes and functions for general-purpose
 use. The file is split into two primary sections:

 -Class definitions
 -Global function declarations

 These sections are explained in further detail in their respective
 documentation blocks.

 Everything in this file is declared in the namespace brgastro.

 \**********************************************************************/

#ifndef __BRG_RANDOM_FUNCTIONS_H_INCLUDED__
#define __BRG_RANDOM_FUNCTIONS_H_INCLUDED__

#include "brg_global.h"

namespace brgastro
{

/** Global function declarations **/
#if (1)

// Returns a random variable from a Gaussian distribution
const double Gaus_rand( const double mean = 0, const double stddev = 1 );

// Returns a random variable from a Gaussian distribution in log space
// Note that "mean" here is the desired mean, NOT the peak of the function (which differ in log space). If you want to use
// the peak, simply use the standard Gaus_rand version instead.
const double log10Gaus_rand( const double mean = 0, const double stddev = 1 );

// Returns a Poisson random variable.
const int Pois_rand( const double lambda = 1 );

#endif // End global function declarations

}

#endif
