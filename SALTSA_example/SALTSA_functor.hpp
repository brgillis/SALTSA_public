/**********************************************************************\
SALTSA_misc_functions.h
 -----------

 If this header is used, the source file SALTSA_misc_functions.cpp must be included
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

 Everything in this file is declared in the namespace SALTSA.

 \**********************************************************************/

#ifndef __SALTSA_MISC_FUNCTIONS_H_INCLUDED__
#define __SALTSA_MISC_FUNCTIONS_H_INCLUDED__

#include <iostream>
#include <vector>

#include "SALTSA_global.h"

namespace SALTSA
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

#endif // end class declarations

} // end namespace SALTSA

#endif // __SALTSA_MISC_FUNCTIONS_H_INCLUDED__
