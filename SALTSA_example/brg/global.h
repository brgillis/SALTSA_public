/**********************************************************************\
  @file global.h

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

#ifndef _BRG_GLOBAL_H_INCLUDED_
#define _BRG_GLOBAL_H_INCLUDED_

// Global compiler directives
// Alter these by switching between #define and #undef

//#define _BRG_WARN_FOR_SAFE_FUNCTIONS_TRIGGERED_ // Warns if a function like "safe_d" prevents an error
// This may be expected or not an issue in some cases though,
// so just undef this for release builds if you're satisfied
// there's no actual problem.

#if(__cplusplus==201103L)
#define _BRG_USE_CPP_11_STD_
#else
#undef _BRG_USE_CPP_11_STD_
#endif

#ifndef NDEBUG
//#define _BRG_USE_UNITS_ // Will use "number-with-units" class for applicable values in code
// This slows things down a bit, but can be useful in debugging. Comment/uncomment or define
// at command line to decide whether or not to use this
#endif // #ifndef NDEBUG

#define _BRG_WARN_FOR_UNIT_MISMATCH_
// Warns in the following scenarios:
// -Adding or subtracting values with incompatible units
// -Setting a variable to something with different units
// Does not warn when:
// -Adding or subtracting to or from zero
// -Adding or subtracting a unitless value to or from an angle
// -Any value in the procedure is not a type with units (eg. it's an int or double)
// -The variable being set is initially unitless

// Set up global parameters

#ifndef MAX_STACK_DEPTH
#define MAX_STACK_DEPTH 100
#endif

#ifndef _BRG_PI_DEFINED_
#define _BRG_PI_DEFINED_
// Defining pi to keep it short, but as a variable so it won't act unusually
const double pi = 3.14159265358979323846;
#endif

#ifndef MIN_DIVISOR
#define MIN_DIVISOR 1e-99
#endif

#ifndef SMALL_FACTOR
#define SMALL_FACTOR 1e-9
#endif

#ifndef INT_MAX
#define INT_MAX std::numeric_limits<int>::max()
#endif

#ifndef INT_MIN
#define INT_MIN std::numeric_limits<int>::min()
#endif

#ifndef DBL_MAX
#define DBL_MAX std::numeric_limits<double>::max()
#endif

#ifndef DBL_MIN
#define DBL_MIN std::numeric_limits<double>::min()
#endif

#ifndef FLT_EPSILON
#define FLT_EPSILON std::numeric_limits<float>::epsilon()
#endif

#ifndef DBL_EPSILON
#define DBL_EPSILON std::numeric_limits<double>::epsilon()
#endif

#ifndef FLT_ROUNDING_EPSILON
#define FLT_ROUNDING_EPSILON 10*FLT_EPSILON
#endif

#ifndef ROUNDING_EPSILON
#define ROUNDING_EPSILON 10*DBL_EPSILON
#endif

#ifndef DBL_MAX_PRECISION
#define DBL_MAX_PRECISION 14
#endif

#ifndef ROMBERG_N_MAX
#define ROMBERG_N_MAX 20
#endif

#ifdef _BRG_USE_UNITS_
#define BRG_UNITS brgastro::unit_obj
#define BRG_DISTANCE brgastro::unit_distance
#define BRG_TIME brgastro::unit_time
#define BRG_MASS brgastro::unit_mass
#define BRG_ANGLE brgastro::unit_angle
#define BRG_CHARGE brgastro::unit_charge
#define BRG_VELOCITY brgastro::unit_velocity

#define CONST_BRG_UNITS_REF const brgastro::unit_obj &
#define CONST_BRG_DISTANCE_REF const brgastro::unit_distance &
#define CONST_BRG_TIME_REF const brgastro::unit_time &
#define CONST_BRG_MASS_REF const brgastro::unit_mass &
#define CONST_BRG_ANGLE_REF const brgastro::unit_angle &
#define CONST_BRG_CHARGE_REF const brgastro::unit_charge &
#define CONST_BRG_VELOCITY_REF const brgastro::unit_velocity &
#else
#define BRG_UNITS double
#define BRG_DISTANCE double
#define BRG_TIME double
#define BRG_MASS double
#define BRG_ANGLE double
#define BRG_CHARGE double
#define BRG_VELOCITY double

#define CONST_BRG_UNITS_REF const double
#define CONST_BRG_DISTANCE_REF const double
#define CONST_BRG_TIME_REF const double
#define CONST_BRG_MASS_REF const double
#define CONST_BRG_ANGLE_REF const double
#define CONST_BRG_CHARGE_REF const double
#define CONST_BRG_VELOCITY_REF const double
#endif // #ifdef _BRG_USE_UNITS_

#ifdef _BRG_USE_CPP_11_STD_
#define BRG_UNIQUE_PTR std::unique_ptr
#define BRG_SHARED_PTR std::shared_ptr
#else
#define BRG_UNIQUE_PTR std::auto_ptr
#define BRG_SHARED_PTR boost::shared_ptr
#endif // #ifdef _BRG_USE_CPP_11_STD_

#ifndef NULL
#ifdef _BRG_USE_CPP_11_STD_
#define NULL nullptr
#else
#define NULL (void *)0
#endif // _BRG_USE_CPP_11_STD_
#endif // #ifndef NULL

// Error code values
#ifndef _BRG_ERR_CODES_DEFINED_
#define _BRG_ERR_CODES_DEFINED_
#define LOWER_LEVEL_ERROR       10
#define UNSPECIFIED_ERROR        1
#define INVALID_ARGUMENTS_ERROR  2
#define OUT_OF_BOUNDS_ERROR      3
#define NOT_SET_UP_ERROR         4
#define MUST_OVERRIDE_ERROR      5
#define MEMORY_ERROR             6
#define FILE_ACCESS_ERROR        7
#define INFINITE_LOOP_ERROR      8
#endif // #ifndef _BRG_ERR_CODES_DEFINED_

#endif
