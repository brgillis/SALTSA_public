/**********************************************************************\
brg_functions.h
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

#ifndef __BRG_FUNCTIONS_H_INCLUDED__
#define __BRG_FUNCTIONS_H_INCLUDED__

#include "brg_global.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <memory>

#include "brg_units.h"
#include "brg_functions.hpp"

namespace brgastro
{

/** Class definitions **/
#if (1)
template< typename T >
class functor
/**********************************************
 functor
 -------

 An abstract class representing an arbitrary function,
 which takes num_in_params input parameters
 and returns num_out_params output parameters, with
 an optional bool to silence error messages.

 This class allows polymorphism to be used by
 other functions to, for instance, integrate a
 defined function.

 Must be overridden by child functions for
 specific implementations.

 **********************************************/
{
public:

	// Virtual destructor
	virtual ~functor()
	{
	}

	virtual const int operator()( const T & in_params, T & out_params,
			const bool silent = false ) const =0;
};

template< typename f1, typename f2, typename T >
class function_product_function: public functor< T >
{
	/*****************************************************
	 function_product_function (non-vector implementation)
	 -----------------------------------------------------

	 An example of a function class which returns the
	 product of two other functions (see the
	 integrate_weighted* functions in this file for
	 how this is used).

	 ****************************************************/
private:

	const functor< T > *_f1_ptr_, *_f2_ptr_;
	bool _f1_set_up_, _f2_set_up_;

public:
	// Constructors

	function_product_function()
	{
		_f1_ptr_ = _f2_ptr_ = 0;
		_f1_set_up_ = _f2_set_up_ = false;
	}

	function_product_function( const f1 *new_f1, const f2 *new_f2 )
	{
		_f1_ptr_ = new_f1;
		_f2_ptr_ = new_f2;
		_f1_set_up_ = _f2_set_up_ = true;
	}

	// Virtual destructor
	virtual ~function_product_function()
	{
	}

	// Set methods

	const int set_f1_ptr( const f1*new_f1_ptr )
	{
		_f1_ptr_ = new_f1_ptr;
		_f1_set_up_ = true;
		return 0;
	}

	const int set_f2_ptr( const f2 *new_f2_ptr )
	{
		_f2_ptr_ = new_f2_ptr;
		_f2_set_up_ = true;
		return 0;
	}

	const int set_f1_f2_ptrs( const f1 *new_f1_ptr, const f2 *new_f2_ptr )
	{
		if ( set_f1_ptr( new_f1_ptr ) )
			return 1;
		return set_f2_ptr( new_f2_ptr );
	}

	// Function method

	const int operator()( const T & in_param, T & out_param,
			const bool silent = false ) const
	{
		if ( ( !_f1_set_up_ ) || ( !_f2_set_up_ ) )
		{
			return NOT_SET_UP_ERROR;
		}

		T f1_out_param, f2_out_param;

		if ( ( *_f1_ptr_ )( in_param, f1_out_param, silent ) )
			return 1;
		if ( ( *_f2_ptr_ )( in_param, f2_out_param, silent ) )
			return 1;

		out_param = f1_out_param * f2_out_param;

		return 0;
	}
};

template< typename f1, typename f2, typename T >
class function_product_function< f1, f2, std::vector< T > > : public functor<
		std::vector< T > >
{
	/**********************************************
	 function_product_function
	 -------------------------

	 An example of a function class which returns the
	 product of two other functions (see the
	 integrate_weighted* functions in this file for
	 how this is used).

	 **********************************************/
private:

	const functor< std::vector< T > > *_f1_ptr_, *_f2_ptr_;
	bool _f1_set_up_, _f2_set_up_;

public:
	// Constructors

	function_product_function()
	{
		_f1_ptr_ = _f2_ptr_ = 0;
		_f1_set_up_ = _f2_set_up_ = false;
	}

	function_product_function( const f1 *new_f1, const f2 *new_f2 )
	{
		_f1_ptr_ = new_f1;
		_f2_ptr_ = new_f2;
		_f1_set_up_ = _f2_set_up_ = true;
	}

	// Virtual destructor
	virtual ~function_product_function()
	{
	}

	// Set methods

	const int set_f1_ptr( const f1 *new_f1_ptr )
	{
		_f1_ptr_ = new_f1_ptr;
		_f1_set_up_ = true;
		return 0;
	}

	const int set_f2_ptr( const f2 *new_f2_ptr )
	{
		_f2_ptr_ = new_f2_ptr;
		_f2_set_up_ = true;
		return 0;
	}

	const int set_f1_f2_ptrs( const f1 *new_f1_ptr, const f2 *new_f2_ptr )
	{
		if ( set_f1_ptr( new_f1_ptr ) )
			return 1;
		return set_f2_ptr( new_f2_ptr );
	}

	// Function method

	const int operator()( const std::vector< T > & in_params,
			std::vector< T > & out_params, const bool silent = false ) const
	{
		if ( ( !_f1_set_up_ ) || ( !_f2_set_up_ ) )
		{
			return NOT_SET_UP_ERROR;
		}

		std::vector< T > f1_out_params( 0 ), f2_out_params( 0 );

		if ( _f1_ptr_( in_params, f1_out_params, silent ) )
			return 1;
		if ( _f2_ptr_( in_params, f2_out_params, silent ) )
			return 1;

		if ( f1_out_params.size() != f2_out_params.size() )
		{
			if ( !silent )
				std::cerr
						<< "ERROR: Functions assigned to function_product_function have\n"
						<< "different numbers of output parameters.\n";
			return UNSPECIFIED_ERROR;
		}

		unsigned int num_out_params = f1_out_params.size();
		out_params.resize( num_out_params );

		for ( unsigned int i = 0; i < num_out_params; i++ )
		{
			out_params.at( i ) = f1_out_params.at( i ) * f2_out_params.at( i );
		}

		return 0;
	}
};

struct phase
/**********************************************
 phase
 -----

 A structure representing the full phase of an
 object (position and velocity) and also time,
 to help limit the number of variables that need
 to be passed to certain functions.

 All member variables are public, and may be
 accessed directly. Units should be in the
 default set: m for position, m/s for velocity,
 s for time.

 **********************************************/
{
	BRG_DISTANCE x, y, z;BRG_VELOCITY vx, vy, vz;BRG_TIME t;
	phase( BRG_DISTANCE init_x = 0, BRG_DISTANCE init_y = 0,
	BRG_DISTANCE init_z = 0, BRG_VELOCITY init_vx = 0,
	BRG_VELOCITY init_vy = 0, BRG_VELOCITY init_vz = 0,
	BRG_TIME init_t = 0 );
	const int set_phase( BRG_DISTANCE init_x = 0, BRG_DISTANCE init_y = 0,
	BRG_DISTANCE init_z = 0, BRG_VELOCITY init_vx = 0,
	BRG_VELOCITY init_vy = 0, BRG_VELOCITY init_vz = 0,
	BRG_TIME init_t = 0 );
};

#endif // end class declarations

/** Global function declarations **/
#if (1)
// Rounds to nearest integer, preferring even values if in the middle of two integers (within epsilon distance of the middle)
const int round_int( const double value, const double epsilon =
		ROUNDING_EPSILON );

// Prints a formatted table in the passed stream. header is a vector of strings representing the labels for each column,
// and data is a 2-d vector of the data to be printed, in the format data[c][r], where c is the column index and r is the row index.
// Prints an error message and returns 1 if an error arises (for instance, vectors are too short).
const int print_table( std::ostream & out_stream, const int num_columns,
		const int num_rows, const std::vector< std::string > & header,
		const std::vector< std::vector< std::string > > & data,
		const bool skip_header = false, const bool silent = false );
const int print_table( std::ostream & out_stream, const int num_columns,
		const int num_rows,
		const std::vector< std::vector< std::string > > & data,
		const bool silent = false );

// Load table, either loading in the entire table, or only loading in certain columns into pointed-to
// variables, found by matching header entries to the strings passed
const int load_table( const std::string & table_file_name, std::vector<std::vector<double> > & table_data,
		const bool silent=false);
const int load_table( const std::string & table_file_name, std::vector<std::vector<double> > & table_data,
		std::vector<std::string> & header, const bool silent=false);
const int load_table_columns( const std::string & table_file_name,
		std::vector< std::pair< const std::string, std::vector<double>* > > & header_links,
		const bool case_sensitive=false, const bool silent=false);

// Functions to open a file and check that it's been opened successfully. An error message will be printed if it can't be opened,
// and the functions will return 1. Otherwise they'll return 0.
const int open_file( std::ofstream & stream, const std::string name,
		const bool silent = false );
const int open_file( std::ifstream & stream, const std::string name,
		const bool silent = false );
const int open_file( std::fstream & stream, const std::string name,
		const bool silent = false );

// Functions to get rid of comments lines (those starting with #) from open fstreams and ifstreams. The "one_line" versions will
// trim the next line if it's commented out, and the "all_at_top" versions will trim all commented lines from the current
// position until they find a line which isn't a comment.
// The "one_line" functions will return 1 if the file is already at the end, and 0 otherwise.
// The "all_at_top" functions will return 1 if they reach the end of the file before the run out of comments, and 0 otherwise.
const int trim_comments_one_line( std::ifstream & stream, const bool silent =
		false );
const int trim_comments_one_line( std::fstream & stream, const bool silent =
		false );
const int trim_comments_all_at_top( std::ifstream & stream, const bool silent =
		false );
const int trim_comments_all_at_top( std::fstream & stream, const bool silent =
		false );

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
