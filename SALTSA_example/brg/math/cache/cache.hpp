/**********************************************************************\
  @file cache.hpp

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

#ifndef _BRG_CACHE_HPP_INCLUDED_
#define _BRG_CACHE_HPP_INCLUDED_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "brg/global.h"

#include "brg/file_functions.h"
#include "brg/math/misc_math.hpp"
#include "brg/math/safe_math.hpp"
#include "brg/physics/units/unit_obj.h"


// Macro definitions

// SPCP: "Static Polymorphic Const Pointer"
// The constness isn't actually enforced, but this is for the reader's understanding
#define SPCP(name) static_cast<const name*>(this)

// SPP: "Static Polymorphic Pointer"
#define SPP(name) static_cast<name*>(this)

#define DECLARE_BRG_CACHE_STATIC_VARS()		\
	static double _min_1_, _max_1_, _step_1_; \
	static unsigned int _resolution_1_;      \
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
	double brgastro::class_name::_min_1_ = init_min;						\
	double brgastro::class_name::_max_1_ = init_max; 						\
	double brgastro::class_name::_step_1_ = init_step; 						\
	bool brgastro::class_name::_loaded_ = false;							\
	bool brgastro::class_name::_initialised_ = false;						\
	short int brgastro::class_name::_is_monotonic_ = 0;						\
	unsigned int brgastro::class_name::_sig_digits_ = 8;					\
	unsigned int brgastro::class_name::_resolution_1_ = 0;					\
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
	void _init() const throw()
	{
		// We check for initialisation twice due to the critical section here.
		// It's expensive to enter, and we don't want to do anything inside it more than once,
		// so we check whether we need to both once outside it and once inside it.
		if(SPCP(name)->_initialised_) return;

		#pragma omp critical(init_brg_cache_nd)
		if(!SPCP(name)->_initialised_)
		{
			SPCP(name)->_resolution_1_ = (unsigned int) max( ( ( SPCP(name)->_max_1_ - SPCP(name)->_min_1_ ) / safe_d(SPCP(name)->_step_1_)) + 1, 2);
			SPCP(name)->_file_name_ = SPCP(name)->_name_base() + "_cache.dat";
			SPCP(name)->_header_string_ = "# " + SPCP(name)->_name_base() + "_cache v1.0";

			SPCP(name)->_initialised_ = true;
		}
	}
	void _load( const bool silent = false ) const
	{
		std::ifstream in_file;
		std::string file_data;
		bool need_to_calc = false;
		unsigned int i;
		int loop_counter = 0;

		if ( SPCP(name)->_loaded_ ) return;

		do
		{
			if ( loop_counter >= 2 )
			{
				throw std::runtime_error("Infinite loop detected trying to load " + SPCP(name)->_file_name_ + " in brgastro::brg_cache.\n");
			}
			else
			{
				loop_counter++;
			}
			need_to_calc = false;

			try
			{
				open_file( in_file, SPCP(name)->_file_name_, true );
			}
			catch(const std::exception &e)
			{
				need_to_calc = true;
				SPCP(name)->_calc( silent );
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}

			// Check that it has the right header
			getline( in_file, file_data );
			if ( file_data.compare( SPCP(name)->_header_string_ ) )
			{
				need_to_calc = true;
				SPCP(name)->_calc( silent );
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}

			// Trim out any other commented lines
			trim_comments_all_at_top( in_file );

			// Load range parameters;
			if ( !( in_file >> SPCP(name)->_min_1_ >> SPCP(name)->_max_1_ >> SPCP(name)->_step_1_ ) )
			{
				need_to_calc = true;
				SPCP(name)->_calc( silent );
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}

			// Set up data
			SPCP(name)->_resolution_1_ = (unsigned int) max( ( ( SPCP(name)->_max_1_ - SPCP(name)->_min_1_ ) / safe_d(SPCP(name)->_step_1_)) + 1, 2);
			SPCP(name)->_results_.resize(SPCP(name)->_resolution_1_ );

			// Read in data

			double temp_data;
			double last_data=0;

			i = 0;
			SPCP(name)->_is_monotonic_ = 0;
			while ( ( !in_file.eof() ) && ( i < SPCP(name)->_resolution_1_ ) )
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
			if ( i < SPCP(name)->_resolution_1_ )
			{
				need_to_calc = true;
				SPCP(name)->_calc( silent );
				SPCP(name)->_unload();
				continue;
			}

		} while ( need_to_calc );

		// Finish up
		in_file.close();
		in_file.clear();
		SPCP(name)->_loaded_ = true;
	}
	void _unload() const throw()
	{
		SPCP(name)->_loaded_ = false;
		SPCP(name)->_results_.clear();
	}
	void _calc( const bool silent = false ) const
	{

		// Test that range is sane
		if ( ( SPCP(name)->_max_1_ <= SPCP(name)->_min_1_ ) || ( SPCP(name)->_step_1_ <= 0 ) )
		{
			throw std::runtime_error("ERROR: Bad range passed to brg_cache::_calc() for " + SPCP(name)->_name_base() + "\n");
		}

		// Set up data
		SPCP(name)->_resolution_1_ = (unsigned int) max( ( ( SPCP(name)->_max_1_ - SPCP(name)->_min_1_ ) / safe_d(SPCP(name)->_step_1_)) + 1, 2);
		SPCP(name)->_results_.resize(SPCP(name)->_resolution_1_ );

		// Calculate data
		double x = SPCP(name)->_min_1_;
		bool bad_result = false;

		#pragma omp parallel for
		for ( unsigned int i = 0; i < SPCP(name)->_resolution_1_; i++ )
		{
			double result = 0;
			try
			{
				result = SPCP(name)->_calculate(x);
			}
			catch(const std::exception &e)
			{
				bad_result = true;
			}
			SPCP(name)->_results_[i] = result;
			x += SPCP(name)->_step_1_;
		}

		if(bad_result) throw std::runtime_error("One or more calculations failed in generating cache " + SPCP(name)->_name_base());
		SPCP(name)->_loaded_ = true;
	}
	void _output( const bool silent = false ) const
	{

		std::ofstream out_file;
		std::string file_data;

		if ( !SPCP(name)->_loaded_ )
		{
			SPCP(name)->_calc( silent );
		}

		open_file( out_file, SPCP(name)->_file_name_, true );

		// Output header
		out_file << SPCP(name)->_header_string_ << "\n#\n";

		// Set number of significant digits
		out_file.precision(SPCP(name)->_sig_digits_);

		// Output range
		out_file << SPCP(name)->_min_1_ << "\t" << SPCP(name)->_max_1_ << "\t" << SPCP(name)->_step_1_ << "\n";

		// Output data
		for ( unsigned int i = 0; i < SPCP(name)->_resolution_1_; i++ )
		{
			out_file << SPCP(name)->_results_.at(i) << "\n";
		}

		out_file.close();
		out_file.clear();
	}
#endif // Private methods

protected:

	// Protected methods
	// These are made protected instead of private so base classes can overload them
#if (1)

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

	// Long calculation function, which is used to generate the cache
	virtual const double _calculate(const double x) const
	{
		return 0;
	}

	// This function should be overloaded to call each cache of the same dimensionality as
	// this cache depends upon in calculation. This is necessary in order to avoid critical
	// sections of the same name being called recursively.
	virtual void _load_cache_dependencies() const
	{
	}

#endif // Protected methods

public:

	// Public methods
#if (1)

	void set_file_name( const std::string & new_name )
	{
		SPP(name)->_file_name_ = new_name;
		if ( SPCP(name)->_loaded_ )
		{
			SPCP(name)->_unload();
		}
	} // set_file_name()

	void set_range( const double new_mins, const double new_maxes,
			const double new_steps, const bool silent = false )
	{
		// First we try to load, so we can see if there are any changes from
		// the existing cache
		if ( !SPCP(name)->_loaded_ )
			SPCP(name)->_load( true );

		// Go through variables, check if any are actually changed. If so, recalculate cache
		if ( ( SPCP(name)->_min_1_ != new_mins ) || ( SPCP(name)->_max_1_ != new_maxes )
				|| ( SPCP(name)->_step_1_ != new_steps ) )
		{
			SPP(name)->_min_1_ = new_mins;
			SPP(name)->_max_1_ = new_maxes;
			SPP(name)->_step_1_ = new_steps;

			SPCP(name)->_unload();
			SPCP(name)->_calc( silent );
		}
	} // const int set_range()

	void set_precision( const unsigned int new_precision,
			const bool silent = false )
	{
		if ( new_precision > 0 )
		{
			SPP(name)->_sig_digits_ = min( new_precision, DBL_MAX_PRECISION );
		}
		else
		{
			throw std::runtime_error("Precision for dfa_cache must be > 0.\n");
		}
	} // const int set_precision()

	void print( std::ostream & out, const bool silent = false ) const
	{
		// Load if necessary
		if ( !SPCP(name)->_loaded_ )
		{
			// Do a test get to make sure it's loaded (and take advantage of the critical section there,
			// so we don't get collisions from loading within two different critical sections at once)
			SPCP(name)->get(SPCP(name)->_min_1_,SPCP(name)->_min_2_,SPCP(name)->_min_3_,
					SPCP(name)->_min_4_,true);
		}

		// Fill up header
		std::vector< std::string > header(3);
		header[0] = "#";
		header[1] = "x_1";
		header[2] = "y";

		// Fill up data
		std::vector< std::vector<std::string> > data(3);
		std::stringstream ss;
		for(unsigned int i_1=0; i_1<SPCP(name)->_resolution_1_; ++i_1)
		{
			data[0].push_back("");
			ss.str("");
			ss << SPCP(name)->_min_1_ + i_1*SPCP(name)->_step_1_;
			data[1].push_back(ss.str());
			ss.str("");
			ss << SPCP(name)->_results_[i_1];
			data[2].push_back(ss.str());
		}

		print_table(out,data,header,silent);
	}

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

		// Load if necessary
		if ( !SPCP(name)->_loaded_ )
		{
			// Load any caches we depend upon before the critical section
			_load_cache_dependencies();

			// Critical section here, since we can't load multiple times simultaneously
			#pragma omp critical(load_brg_cache)
			{
				try
				{
					SPCP(name)->_load( silent );
				}
				catch(const std::exception &e)
				{
					result = -1;
				}
			}
			if ( result == -1 )
			{
				throw std::runtime_error("ERROR: Could neither load " + SPCP(name)->_file_name_ + " nor calculate in brg_cache::get()\n");
			}
		}

		x_i = (unsigned int)bound(0,
				( ( x - SPCP(name)->_min_1_ ) / SPCP(name)->_step_1_ ),
				SPCP(name)->_resolution_1_ - 2 );

		xlo = SPCP(name)->_min_1_ + SPCP(name)->_step_1_ * x_i;
		xhi = SPCP(name)->_min_1_ + SPCP(name)->_step_1_ * ( x_i + 1 );

		result = ( ( x - xlo ) * SPCP(name)->_results_.at(x_i + 1) + ( xhi - x ) * SPCP(name)->_results_.at(x_i) )
				/ SPCP(name)->_step_1_;

		return result;

	} // get()

	const BRG_UNITS inverse_get( const double y, const bool silent = false ) const
	{
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
			// Do a test get to make sure it's loaded (and take advantage of the critical section there,
			// so we don't get collisions from loading within two different critical sections at once)
			SPCP(name)->get(SPCP(name)->_min_1_,true);
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

			for ( unsigned int x_i = 0; x_i < SPCP(name)->_resolution_1_ - 1; x_i++ )
			{
				// Loop through till we find the proper y or reach the end
				yhi = SPCP(name)->_results_.at(x_i);
				if ( ( yhi > y ) || (x_i >= SPCP(name)->_resolution_1_ - 2) )
				{
					ylo = SPCP(name)->_results_.at(x_i + 1);

					xlo = SPCP(name)->_min_1_ + SPCP(name)->_step_1_ * x_i;
					xhi = SPCP(name)->_min_1_ + SPCP(name)->_step_1_ * ( x_i + 1 );
					result = xlo + ( xhi - xlo ) * ( y - ylo ) / safe_d( yhi - ylo );
					break;
				}
			}
		} // if(_is_monotonic_==1)
		else
		{

			for ( unsigned int x_i = 0; x_i < SPCP(name)->_resolution_1_ - 1; x_i++ )
			{
				// Loop through till we find the proper y or reach the end
				ylo = SPCP(name)->_results_.at(x_i);
				if ( ( ylo < y ) || (x_i >= SPCP(name)->_resolution_1_ - 2) )
				{
					yhi = SPCP(name)->_results_.at(x_i + 1);

					xlo = SPCP(name)->_min_1_ + SPCP(name)->_step_1_ * x_i;
					xhi = SPCP(name)->_min_1_ + SPCP(name)->_step_1_ * ( x_i + 1 );
					result = xlo + ( xhi - xlo ) * ( y - ylo ) / safe_d( yhi - ylo );
					break;
				}
			}

		} // _is_monotonic == -1

		return result;

	}

	// Recalculate function. Call if you want to overwrite a cache when something's changed in the code
	// (for instance, the _calculate() function has been altered)
	void recalc( const bool silent = false ) const
	{
		SPCP(name)->_unload();
		SPCP(name)->_calc(silent);
		SPCP(name)->_output(silent);
		SPCP(name)->_unload();
		SPCP(name)->_load(silent);
	}

	// Constructor
	brg_cache()
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

// Undef macros
#undef SPP
#undef SPCP

#endif // __BRG_CACHE_HPP_INCLUDED__
