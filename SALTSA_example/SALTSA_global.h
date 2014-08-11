#ifndef __SALTSA_GLOBAL_H_INCLUDED__
#define __SALTSA_GLOBAL_H_INCLUDED__

// Global compiler directives
// Alter these by switching between #define and #undef

#undef _BRG_WARN_FOR_SAFE_FUNCTIONS_TRIGGERED_ // Warns if a function like "safe_d" prevents an error
// This may be expected or not an issue in some cases though,
// so just undef this for release builds if you're satisfied
// there's no actual problem.

#if(__cplusplus==201103L)
#define _BRG_USE_CPP_11_STD_
#else
#undef _BRG_USE_CPP_11_STD_
#endif

// Set up global parameters

#ifndef MAX_STACK_DEPTH
#define MAX_STACK_DEPTH 100
#endif

#ifndef __PI_DEFINED__
#define __PI_DEFINED__
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
#define INT_MAX 2147483647 // Most likely value
#endif

#ifndef INT_MIN
#define INT_MIN -2147483647 // Most likely value, or possibly one below this
#endif

#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623158e+308 // Most likely value
#endif

#ifndef DBL_MIN
#define DBL_MIN 2.2250738585072014e-308 // Most likely value
#endif

#ifndef FLT_EPSILON
#define FLT_EPSILON 1e-07
#endif

#ifndef DBL_EPSILON
#define DBL_EPSILON 1e-14
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

#ifndef RHOMBERG_N_MAX
#define RHOMBERG_N_MAX 100
#endif

#ifdef _BRG_USE_CPP_11_STD_
#define BRG_UNIQUE_PTR std::unique_ptr
#define BRG_SHARED_PTR std::shared_ptr
#else
#define BRG_UNIQUE_PTR boost::shared_ptr // Next best thing
#define BRG_SHARED_PTR boost::shared_ptr
#endif // #ifdef _BRG_USE_CPP_11_STD_

#ifndef NULL
#define NULL (void *)0
#endif // #ifndef NULL

// Error code values
#ifndef __SALTSA_ERR_CODES_DEFINED__
#define __SALTSA_ERR_CODES_DEFINED__
#define LOWER_LEVEL_ERROR       10
#define UNSPECIFIED_ERROR        1
#define INVALID_ARGUMENTS_ERROR  2
#define OUT_OF_BOUNDS_ERROR      3
#define NOT_SET_UP_ERROR         4
#define MUST_OVERRIDE_ERROR      5
#define MEMORY_ERROR             6
#define FILE_ACCESS_ERROR        7
#define INFINITE_LOOP_ERROR      8
#endif // #ifndef __SALTSA_ERR_CODES_DEFINED__

#endif // __SALTSA_GLOBAL_H_INCLUDED__
