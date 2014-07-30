/*
 * brg_cache.h
 *
 *  Created on: 25 Apr 2014
 *      Author: brg
 */

#ifndef __BRG_CACHE_HPP_INCLUDED__
#define __BRG_CACHE_HPP_INCLUDED__

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "brg_global.h"
#include "brg_units.h"
#include "brg_functions.h"

// Macro definitions

// SPCP: "Static Polymorphic Const Pointer"
// The constness isn't actually enforced, but this is for the reader's understanding
#define SPCP(name) static_cast<const name*>(this)

// SPP: "Static Polymorphic Pointer"
#define SPP(name) static_cast<name*>(this)

#define DECLARE_BRG_CACHE_STATIC_VARS()		\
	static double _mins_, _maxes_, _steps_; \
	static unsigned int _resolutions_;      \
	static std::vector< double > _results_; \
											\
	static std::string _file_name_;         \
	static std::string _header_string_;     \
											\
	static bool _loaded_, _initialised_;    \
	static short int _is_monotonic_;		\
											\
	static unsigned int _sig_digits_;
#define DEFINE_BRG_CACHE_STATIC_VARS(class_name,init_min,init_max,init_step) \
	double brgastro::class_name::_mins_ = init_min;							\
	double brgastro::class_name::_maxes_ = init_max; 						\
	double brgastro::class_name::_steps_ = init_step; 						\
	bool brgastro::class_name::_loaded_ = false;							\
	bool brgastro::class_name::_initialised_ = false;						\
	short int brgastro::class_name::_is_monotonic_ = 0;						\
	unsigned int brgastro::class_name::_sig_digits_ = 8;					\
	unsigned int brgastro::class_name::_resolutions_ = 0;					\
	std::string brgastro::class_name::_file_name_ = "";						\
	std::string brgastro::class_name::_header_string_ = "";					\
	std::vector<double> brgastro::class_name::_results_;

namespace brgastro
{

template<typename name>
class brg_cache
{
private:

	// Private variables
#if (1)

	DECLARE_BRG_CACHE_STATIC_VARS();

#endif // Private variables

	// Private methods
#if (1)
	const int _init() const throw()
	{
		if(SPCP(name)->_initialised_) return 0;

		SPCP(name)->_resolutions_ = (unsigned int) max( ( ( SPCP(name)->_maxes_ - SPCP(name)->_mins_ ) / safe_d(SPCP(name)->_steps_)) + 1, 1);
		SPCP(name)->_file_name_ = SPCP(name)->_name_base() + "_cache.dat";
		SPCP(name)->_header_string_ = "# " + SPCP(name)->_name_base() + "_cache v1.0";

		SPCP(name)->_initialised_ = true;

		return 0;
	}
	const int _load( const bool silent = false ) const
	{
		std::ifstream in_file;
		std::string file_data;
		bool need_to_calc = false;
		unsigned int i;
		int loop_counter = 0;

		if ( SPCP(name)->_loaded_ )
			return 0;

		do
		{
			if ( loop_counter >= 2 )
			{
				if ( !silent )
					std::cerr
							<< "ERROR: infinite loop detected trying to load " + _file_name_ + " in brgastro::brg_cache.\n";
				return INFINITE_LOOP_ERROR;
			}
			else
			{
				loop_counter++;
			}
			need_to_calc = false;

			if ( open_file( in_file, SPCP(name)->_file_name_, true ) )
			{
				need_to_calc = true;
				if ( SPCP(name)->_calc( silent ) )
					return UNSPECIFIED_ERROR;
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}

			// Check that it has the right header
			getline( in_file, file_data );
			if ( file_data.compare( SPCP(name)->_header_string_ ) )
			{
				need_to_calc = true;
				if ( SPCP(name)->_calc( silent ) )
					return UNSPECIFIED_ERROR;
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}

			// Trim out any other commented lines
			if ( trim_comments_all_at_top( in_file ) )
			{
				need_to_calc = true;
				if ( SPCP(name)->_calc( silent ) )
					return UNSPECIFIED_ERROR;
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}

			// Load range parameters;
			if ( !( in_file >> SPCP(name)->_mins_ >> SPCP(name)->_maxes_ >> SPCP(name)->_steps_ ) )
			{
				need_to_calc = true;
				if ( SPCP(name)->_calc( silent ) )
					return UNSPECIFIED_ERROR;
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}

			// Set up data
			SPCP(name)->_resolutions_ = (int)( ( SPCP(name)->_maxes_ - SPCP(name)->_mins_ ) / SPCP(name)->_steps_ ) + 1;
			make_array( SPCP(name)->_results_, SPCP(name)->_resolutions_ );

			// Read in data

			double temp_data;
			double last_data=0;

			i = 0;
			SPCP(name)->_is_monotonic_ = 0;
			while ( ( !in_file.eof() ) && ( i < SPCP(name)->_resolutions_ ) )
			{
				in_file >> temp_data;
				SPCP(name)->_results_.at(i) = temp_data;
				if(i==1)
				{
					// First monotonic check, so we don't compare to its past values
					if(temp_data > last_data)
					{
						SPCP(name)->_is_monotonic_ = 1;
					}
					else if(temp_data < last_data)
					{
						SPCP(name)->_is_monotonic_ = -1;
					}
					else
					{
						SPCP(name)->_is_monotonic_ = 0;
					}
				}
				else if(i>1)
				{
					// Check for monotonic increase/decrease
					if(temp_data > last_data)
					{
						if(SPCP(name)->_is_monotonic_ != 1)
							SPCP(name)->_is_monotonic_ = 0;
					}
					else if(temp_data < last_data)
					{
						if(SPCP(name)->_is_monotonic_ != -1)
							SPCP(name)->_is_monotonic_ = 0;
					}
					else
					{
						SPCP(name)->_is_monotonic_ = 0;
					}
				}
				last_data = temp_data;
				i++;
			}

			// Check that it was all read properly
			if ( i < SPCP(name)->_resolutions_ )
			{
				need_to_calc = true;
				if ( SPCP(name)->_calc( silent ) )
					return UNSPECIFIED_ERROR;
				SPCP(name)->_unload();
				continue;
			}

		} while ( need_to_calc );

		// Finish up
		in_file.close();
		in_file.clear();
		SPCP(name)->_loaded_ = true;
		return 0;
	}
	const int _unload() const throw()
	{
		SPCP(name)->_loaded_ = false;
		SPCP(name)->_results_.clear();
		return 0;
	}
	const int _calc( const bool silent = false ) const
	{

		// Test that range is sane
		if ( ( SPCP(name)->_maxes_ <= SPCP(name)->_mins_ ) || ( SPCP(name)->_steps_ <= 0 ) )
		{
			if ( !silent )
				std::cerr
						<< "ERROR: Bad range passed to brg_cache::_calc() for " + SPCP(name)->_name_base() + "\n";
			return INVALID_ARGUMENTS_ERROR;
		}

		// Set up data
		SPCP(name)->_resolutions_ = (int)( ( SPCP(name)->_maxes_ - SPCP(name)->_mins_ ) / SPCP(name)->_steps_ ) + 1;
		if ( make_array( SPCP(name)->_results_, SPCP(name)->_resolutions_ ) )
			return 1;

		// Calculate data
		double x = SPCP(name)->_mins_;
		bool bad_result = false;
#pragma omp parallel for
		for ( unsigned int i = 0; i < SPCP(name)->_resolutions_; i++ )
		{
			double result = 0;
			if(SPCP(name)->_calculate(x, result))
			{
				bad_result = true;
			}
			SPCP(name)->_results_.at(i) = result;
			x += SPCP(name)->_steps_;
		}

		if(bad_result) return UNSPECIFIED_ERROR;
		SPCP(name)->_loaded_ = true;

		return 0;
	}
	const int _output( const bool silent = false ) const
	{

		std::ofstream out_file;
		std::string file_data;

		if ( !SPCP(name)->_loaded_ )
		{
			if ( SPCP(name)->_calc( silent ) )
				return 1;
		}

		if ( open_file( out_file, SPCP(name)->_file_name_, true ) )
			return 1;

		// Output header
		out_file << SPCP(name)->_header_string_ << "\n#\n";

		// Set number of significant digits
		out_file.precision(_sig_digits_);

		// Output range
		out_file << SPCP(name)->_mins_ << "\t" << SPCP(name)->_maxes_ << "\t" << SPCP(name)->_steps_ << "\n";

		// Output data
		for ( unsigned int i = 0; i < SPCP(name)->_resolutions_; i++ )
		{
			if ( !( out_file << SPCP(name)->_results_.at(i) << "\n" ) )
				return errorNOS( silent );
		}

		out_file.close();
		out_file.clear();

		return 0;
	}
#endif // Private methods

protected:

	// Protected methods
	// These are made protected instead of private so base classes can overload them
#if (1)
	// Name base - derived classes must overload this to tell how their cache file will be named
	// virtual const std::string _name_base() const throw() =0;

#ifdef _BRG_USE_UNITS_

	// Tells what units the result should have. Only the units matter in the return, not the value
	virtual const brgastro::unit_obj _units() const throw()
	{
		return brgastro::unit_obj(0);
	}
	virtual const brgastro::unit_obj _inverse_units() const throw()
	{
		return brgastro::unit_obj(0);
	}

#endif // _BRG_USE_UNITS_

	// Long-form calculation function. Must be overloaded by child classes
	// virtual const int _calculate( const double in_params, double out_params ) const =0;

#endif // Protected methods

public:

	// Public methods
#if (1)

	const int set_file_name( const std::string new_name )
	{
		if(!SPCP(name)->_initialised_) SPP(name)->_init();
		SPP(name)->_file_name_ = new_name;
		if ( SPCP(name)->_loaded_ )
		{
			return SPCP(name)->_unload();
		}
		return 0;
	} // const int set_file_name()

	const int set_range( const double new_mins, const double new_maxes,
			const double new_steps, const bool silent = false )
	{
		if(!SPCP(name)->_initialised_) SPP(name)->_init();

		// First we try to load, so we can see if there are any changes from
		// the existing cache
		if ( !SPCP(name)->_loaded_ )
			SPCP(name)->_load( true );

		// Go through variables, check if any are actually changed. If so, recalculate cache
		if ( ( SPCP(name)->_mins_ != new_mins ) || ( SPCP(name)->_maxes_ != new_maxes )
				|| ( SPCP(name)->_steps_ != new_steps ) )
		{
			SPP(name)->_mins_ = new_mins;
			SPP(name)->_maxes_ = new_maxes;
			SPP(name)->_steps_ = new_steps;

			if ( SPCP(name)->_unload() )
				return errorNOS( silent );
			if ( SPCP(name)->_calc( silent ) )
				return UNSPECIFIED_ERROR;
		}
		return 0;
	} // const int set_range()

	const int set_precision( const unsigned int new_precision,
			const bool silent = false )
	{
		if(!SPCP(name)->_initialised_) SPP(name)->_init();

		if ( new_precision > 0 )
		{
			SPP(name)->_sig_digits_ = min( new_precision, DBL_MAX_PRECISION );
			return 0;
		}
		else
		{
			if ( !silent )
				std::cerr << "ERROR: Precision for dfa_cache must be > 0.\n";
			return INVALID_ARGUMENTS_ERROR;
		}
	} // const int set_precision()

	const BRG_UNITS get( const double x, const bool silent = false ) const
	{

		double xlo, xhi;
		unsigned int x_i; // Lower nearby array point
#ifdef _BRG_USE_UNITS_
		BRG_UNITS result = SPCP(name)->_units(); // Ensure the result has the proper units
		result = 0;
#else
		double result = 0;
#endif

		if(!SPCP(name)->_initialised_) SPCP(name)->_init();

		// Load if necessary
		if ( !SPCP(name)->_loaded_ )
		{
			// Critical section here, since we can't load multiple times simultaneously
			#pragma omp critical(load_brg_cache)
			{
				if ( SPCP(name)->_load( silent ) )
				{
					result = -1;
				}
			}
			if ( result == -1 )
			{
				std::string err = "ERROR: Could neither load " + SPCP(name)->_file_name_ + " nor calculate in brg_cache::get()\n";
				if ( !silent )
					std::cerr << err;
				throw std::runtime_error(err);
			}
		}

		x_i = (unsigned int)( ( x - SPCP(name)->_mins_ ) / SPCP(name)->_steps_ );
		x_i = max( x_i, (unsigned)0 );
		x_i = min( x_i, SPCP(name)->_resolutions_ - 2 );

		xlo = SPCP(name)->_mins_ + SPCP(name)->_steps_ * x_i;
		xhi = SPCP(name)->_mins_ + SPCP(name)->_steps_ * ( x_i + 1 );

		result = ( ( x - xlo ) * SPCP(name)->_results_.at(x_i + 1) + ( xhi - x ) * SPCP(name)->_results_.at(x_i) )
				/ SPCP(name)->_steps_;

		return result;

	} // get()

	const BRG_UNITS inverse_get( const double y, const bool silent = false ) const
	{
		if(!SPCP(name)->_initialised_) SPCP(name)->_init();

		// Check if it's possible to do an inverse get
		if((SPCP(name)->_is_monotonic_!=1)&&((SPCP(name)->_is_monotonic_!=-1)))
		{
			// Not a monotonic function. Inverse get isn't possible
			std::string err = "ERROR: Attempt to use inverse_get in cache for " + SPCP(name)->_file_name_ + " for function which isn't monotonic.\n";
			if ( !silent )
				std::cerr << err;
			throw std::runtime_error(err);
		}


		double xlo, xhi, ylo, yhi;
#ifdef _BRG_USE_UNITS_
		BRG_UNITS result = SPCP(name)->_inverse_units(); // Ensure the result has the proper units
		result = 0;
#else
		double result = 0;
#endif

		if ( !SPCP(name)->_loaded_ )
		{
			#pragma omp critical(load_brg_cache)
			if ( SPCP(name)->_load( silent ) )
			{
				result = -1;
			}
		}
		if ( result == -1 )
		{
			std::string err = "ERROR: Could neither load " + SPCP(name)->_file_name_ + " nor calculate in brg_cache::inverse_get()\n";
			if ( !silent )
				std::cerr << err;
			throw std::runtime_error(err);
		}

		if(SPCP(name)->_is_monotonic_==1)
		{

			for ( unsigned int x_i = 0; x_i < SPCP(name)->_resolutions_ - 1; x_i++ )
			{
				// Loop through till we find the proper y or reach the end
				yhi = SPCP(name)->_results_.at(x_i);
				if ( ( yhi > y ) || (x_i >= SPCP(name)->_resolutions_ - 2) )
				{
					ylo = SPCP(name)->_results_.at(x_i + 1);

					xlo = SPCP(name)->_mins_ + SPCP(name)->_steps_ * x_i;
					xhi = SPCP(name)->_mins_ + SPCP(name)->_steps_ * ( x_i + 1 );
					result = xlo + ( xhi - xlo ) * ( y - ylo ) / safe_d( yhi - ylo );
					break;
				}
			}
		} // if(_is_monotonic_==1)
		else
		{

			for ( unsigned int x_i = 0; x_i < SPCP(name)->_resolutions_ - 1; x_i++ )
			{
				// Loop through till we find the proper y or reach the end
				ylo = SPCP(name)->_results_.at(x_i);
				if ( ( ylo < y ) || (x_i >= SPCP(name)->_resolutions_ - 2) )
				{
					yhi = SPCP(name)->_results_.at(x_i + 1);

					xlo = SPCP(name)->_mins_ + SPCP(name)->_steps_ * x_i;
					xhi = SPCP(name)->_mins_ + SPCP(name)->_steps_ * ( x_i + 1 );
					result = xlo + ( xhi - xlo ) * ( y - ylo ) / safe_d( yhi - ylo );
					break;
				}
			}

		} // _is_monotonic == -1

		return result;

	}

	// Recalculate function. Call if you want to overwrite a cache when something's changed in the code
	// (for instance, the _calculate() function has been altered)
	const int recalc( const bool silent = false ) const
	{
		SPCP(name)->_unload();
		if(SPCP(name)->_calc(silent)) return LOWER_LEVEL_ERROR;
		if(SPCP(name)->_output(silent)) return LOWER_LEVEL_ERROR;
		SPCP(name)->_unload();
		if(SPCP(name)->_load(silent)) return LOWER_LEVEL_ERROR;
		return 0;
	}

	// Constructor
	brg_cache() throw()
	{
		if(!SPCP(name)->_initialised_) SPP(name)->_init();
	}

	// Deconstructor
	virtual ~brg_cache()
	{
	}

#endif // Public methods

}; // class brg_cache

} // namespace brgastro

// Initialisation of static members
#if (1)

// These ones should be specialised for default values
template<typename name>
double brgastro::brg_cache<name>::_mins_ = 0;

template<typename name>
double brgastro::brg_cache<name>::_maxes_ = 0;

template<typename name>
double brgastro::brg_cache<name>::_steps_ = 1; // To avoid divide-by-zero issues

// These ones don't need to be altered for each child, but they should be copied
template<typename name>
bool brgastro::brg_cache<name>::_loaded_ = false;

template<typename name>
bool brgastro::brg_cache<name>::_initialised_ = false;

template<typename name>
short int brgastro::brg_cache<name>::_is_monotonic_ = 0;

template<typename name>
unsigned int brgastro::brg_cache<name>::_sig_digits_ = 12;

template<typename name>
unsigned int brgastro::brg_cache<name>::_resolutions_ = 0;

template<typename name>
std::string brgastro::brg_cache<name>::_file_name_ = "";

template<typename name>
std::string brgastro::brg_cache<name>::_header_string_ = "";

template<typename name>
std::vector<double> brgastro::brg_cache<name>::_results_;

#endif

// Undef macros
#undef SPP
#undef SPCP

#endif // __BRG_CACHE_HPP_INCLUDED__
