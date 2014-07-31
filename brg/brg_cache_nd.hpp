/*
 * brg_cache.h
 *
 *  Created on: 25 Apr 2014
 *      Author: brg
 */

// !!! UNTESTED

#ifndef __BRG_CACHE_ND_HPP_INCLUDED__
#define __BRG_CACHE_ND_HPP_INCLUDED__

#define BRG_CACHE_ND_NAME_SIZE 8

#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <exception>
#include "brg_vector.hpp"
#include "brg_vector_functions.hpp"
#include "brg_global.h"
#include "brg_units.h"

// Macro definitions

// SPCP: "Static Polymorphic Const Pointer"
// The constness isn't actually enforced, but this is for the reader's understanding
#define SPCP(name) static_cast<const name*>(this)

// SPP: "Static Polymorphic Pointer"
#define SPP(name) static_cast<name*>(this)

#define DECLARE_BRG_CACHE_ND_STATIC_VARS()		       \
	static brgastro::vector<double> _mins_, _maxes_, _steps_;      \
	static brgastro::vector<unsigned int> _resolutions_;           \
	static brgastro::vector<double> _results_;                     \
											                       \
	static std::string _file_name_;                                \
	static unsigned int _version_number_;                          \
										                           \
	static bool _loaded_, _initialised_;                           \
											                       \
	static unsigned short int _num_dim_;

// Be careful when using this not to use the default constructor for init_steps, which would result in
// divide-by-zero errors
#define DEFINE_BRG_CACHE_ND_STATIC_VARS(class_name,init_mins,init_maxes,init_steps,init_num_dim) \
	brgastro::vector<double> brgastro::class_name::_mins_ = init_mins;	                         \
	brgastro::vector<double> brgastro::class_name::_maxes_ = init_maxes;                         \
	brgastro::vector<double> brgastro::class_name::_steps_ = init_steps;                         \
	bool brgastro::class_name::_loaded_ = false;							                     \
	bool brgastro::class_name::_initialised_ = false;					                         \
	brgastro::vector<unsigned int> brgastro::class_name::_resolutions_ =                         \
		max( (((brgastro::class_name::_maxes_-brgastro::class_name::_mins_) /                    \
				safe_d(brgastro::class_name::_steps_))+1), 1);                                   \
	std::string brgastro::class_name::_file_name_ = "";					     	                 \
	unsigned int brgastro::class_name::_version_number_ = 0;		     		                 \
	brgastro::vector<double> brgastro::class_name::_results_;                                    \
	unsigned short int brgastro::class_name::_num_dim_ = init_num_dim;

namespace brgastro
{

template<typename name>
class brg_cache_nd
{
private:

	// Private variables
#if (1)

	DECLARE_BRG_CACHE_ND_STATIC_VARS();

#endif // Private variables

	// Private methods
#if (1)
	const int _init() const throw()
	{
		if(SPCP(name)->_initialised_) return 0;

		SPCP(name)->_resolutions_ = max( (((SPCP(name)->_maxes_-SPCP(name)->_mins_) / safe_d(SPCP(name)->_steps_))+1), 1);
		SPCP(name)->_file_name_ = SPCP(name)->_name_base() + "_cache.bin";
		SPCP(name)->_version_number_ = 0; // This should be changed when there are changes to this code

		SPCP(name)->_initialised_ = true;

		return 0;
	}
	const int _load( const bool silent = false ) const
	{
		std::ifstream in_file;
		std::string file_data;
		bool need_to_calc = false;
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

			if ( open_bin_file( in_file, SPCP(name)->_file_name_, true ) )
			{
				need_to_calc = true;
				if ( SPCP(name)->_calc( silent ) )
					return UNSPECIFIED_ERROR;
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}

			// Check that it has the right name and version

			char file_name[BRG_CACHE_ND_NAME_SIZE];
			unsigned int file_version = std::numeric_limits<unsigned short int>::max();

			in_file.read(file_name,BRG_CACHE_ND_NAME_SIZE);
			in_file.read((char *)&file_version,sizeof(file_version));

			if( (!in_file) || (((std::string)file_name) != SPCP(name)->_name_base()) ||
					(file_version != SPCP(name)->_version_number_) )
			{
				need_to_calc = true;
				if ( SPCP(name)->_calc( silent ) )
					return UNSPECIFIED_ERROR;
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}

			// Load range parameters;
			SPCP(name)->_mins_.resize(SPCP(name)->_num_dim_,0);
			SPCP(name)->_maxes_.resize(SPCP(name)->_num_dim_,0);
			SPCP(name)->_steps_.resize(SPCP(name)->_num_dim_,0);
			for(unsigned int i = 0; i < SPCP(name)->_num_dim_; i++)
			{
				in_file.read((char *)&(SPCP(name)->_mins_[i]),sizeof(SPCP(name)->_mins_[i]));
				in_file.read((char *)&(SPCP(name)->_maxes_[i]),sizeof(SPCP(name)->_maxes_[i]));
				in_file.read((char *)&(SPCP(name)->_steps_[i]),sizeof(SPCP(name)->_steps_[i]));
			}
			if ( !(in_file) )
			{
				need_to_calc = true;
				if ( SPCP(name)->_calc( silent ) )
					return UNSPECIFIED_ERROR;
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}

			// Set up data
			SPCP(name)->_resolutions_ = max( (((SPCP(name)->_maxes_-SPCP(name)->_mins_) / safe_d(SPCP(name)->_steps_))+1.), 1.);
			SPCP(name)->_results_.reshape(SPCP(name)->_resolutions_.v(),0);

			// Read in data

			// Initialize
			unsigned int i = 0;
			const std::streamsize size = sizeof(SPCP(name)->_results_[0]); // Store the size for speed
			brgastro::vector<unsigned int> position(SPCP(name)->_num_dim_,0);

			while ( ( !in_file.eof() ) && ( i < product(SPCP(name)->_resolutions_) ) && (in_file) )
			{
				in_file.read((char *)&(SPCP(name)->_results_(position)),size);
				i++;

				for(unsigned int d=0; d<SPCP(name)->_num_dim_; d++)
				{
					position[d]++;
					if(position[d] != SPCP(name)->_resolutions_)
						break;
					position[d] = 0;
					// If we get here, we'll go on to increase the next index by 1
				}
			}

			// Check that it was all read properly
			if ( (i != product(SPCP(name)->_resolutions_)) || (!in_file) )
			{
				need_to_calc = true;
				if ( SPCP(name)->_calc( silent ) )
					return UNSPECIFIED_ERROR;
				SPCP(name)->_output();
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
		for(unsigned int i = 0; i < SPCP(name)->_num_dim_; i++)
		{
			try {
				if ( ( SPCP(name)->_maxes_.at(i) <= SPCP(name)->_mins_.at(i) ) || ( SPCP(name)->_steps_.at(i) <= 0 ) )
				{
					if ( !silent )
						std::cerr
								<< "ERROR: Bad range passed to brg_cache::_calc() for " + static_cast<const name*>(this)->_name_base() + "\n";
					return INVALID_ARGUMENTS_ERROR;
				}
			} catch (std::exception &e) {
				if ( !silent )
					std::cerr
							<< "ERROR: Bad range passed to brg_cache::_calc() for " + static_cast<const name*>(this)->_name_base() + "\n";
				return INVALID_ARGUMENTS_ERROR;
			}
		}

		// Set up data
		SPCP(name)->_resolutions_.resize(SPCP(name)->_num_dim_);
		SPCP(name)->_resolutions_ = max( (((SPCP(name)->_maxes_-SPCP(name)->_mins_) / safe_d(SPCP(name)->_steps_))+1.), 1.);
		SPCP(name)->_results_.reshape(SPCP(name)->_resolutions_.v() );

		brgastro::vector<unsigned int> position(SPCP(name)->_num_dim_,0);
		brgastro::vector<double> x(SPCP(name)->_num_dim_,0);
		for ( unsigned int i = 0; i < SPCP(name)->_results_.size(); i++ )
		{
			x = SPCP(name)->_mins_ + SPCP(name)->_steps_*position;
			double result = 0;
			if(SPCP(name)->_calculate(x, result))
			{
				return UNSPECIFIED_ERROR;
			}
			SPCP(name)->_results_(position) = result;

			for(unsigned int d=0; d<SPCP(name)->_num_dim_; d++)
			{
				position[d]++;
				if(position[d] != SPCP(name)->_resolutions_)
					break;
				position[d] = 0;
				// If we get here, we'll go on to increase the next index by 1
			}
		}

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
				return LOWER_LEVEL_ERROR;
		}

		if ( open_bin_file( out_file, SPCP(name)->_file_name_, true ) )
			return LOWER_LEVEL_ERROR;

		// Output name and version

		std::string file_name = SPCP(name)->_name_base();
		unsigned int file_version = SPCP(name)->_version_number_;

		out_file.write(file_name.c_str(),BRG_CACHE_ND_NAME_SIZE);
		out_file.write((char *)&file_version,sizeof(file_version));

		// Output range parameters
		for(unsigned int i = 0; i < SPCP(name)->_num_dim_; i++)
		{
			out_file.write((char *)&(SPCP(name)->_mins_[i]),sizeof(SPCP(name)->_mins_[i]));
			out_file.write((char *)&(SPCP(name)->_maxes_[i]),sizeof(SPCP(name)->_maxes_[i]));
			out_file.write((char *)&(SPCP(name)->_steps_[i]),sizeof(SPCP(name)->_steps_[i]));
		}

		// Output data

		// Initialize
		unsigned int i = 0;
		const std::streamsize size = sizeof(SPCP(name)->_results_[0]);
		brgastro::vector<unsigned int> position(SPCP(name)->_num_dim_,0);

		while ( i < product(SPCP(name)->_resolutions_) )
		{
			out_file.write((char *)&(SPCP(name)->_results_(position)),size);
			i++;

			for(unsigned int d=0; d<SPCP(name)->_num_dim_; d++)
			{
				position[d]++;
				if(position[d] != SPCP(name)->_resolutions_)
					break;
				position[d] = 0;
				// If we get here, we'll go on to increase the next index by 1
			}
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

	const std::string _name_base() const throw()
	{
		char name_base[BRG_CACHE_ND_NAME_SIZE] = "";
		return name_base;
	}

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

	const int set_range( const brgastro::vector<double> & new_mins, const brgastro::vector<double> & new_maxes,
			const brgastro::vector<double> & new_steps, const bool silent = false )
	{
		if(!SPCP(name)->_initialised_) SPP(name)->_init();

		// First we try to load, so we can see if there are any changes from
		// the existing cache
		if ( !SPCP(name)->_loaded_ )
			SPCP(name)->_load( true );

		// Check sizes of passed vectors
		if( (new_mins.size() != SPCP(name)->_num_dim_) || (new_maxes.size() != SPCP(name)->_num_dim_) ||
				(new_steps.size() != SPCP(name)->_num_dim_) )
		{
			throw std::runtime_error("ERROR: Incorrect sizes of vectors passed to set_range.\n");
		}

		// Go through variables, check if any are actually changed. If so, recalculate cache
		for(unsigned int i = 0; i < SPCP(name)->_num_dim_; i++)
		{
			if ( ( SPCP(name)->_mins_.at(i) != new_mins.at(i) ) || ( SPCP(name)->_maxes_.at(i) != new_maxes.at(i) )
					|| ( SPCP(name)->_steps_.at(i) != new_steps.at(i) ) )
			{
				SPP(name)->_mins_ = new_mins;
				SPP(name)->_maxes_ = new_maxes;
				SPP(name)->_steps_ = new_steps;

				if ( SPCP(name)->_unload() )
					return errorNOS( silent );
				if ( SPCP(name)->_calc( silent ) )
					return LOWER_LEVEL_ERROR;
				break;
			}
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

	const BRG_UNITS get( const brgastro::vector<double> & x, const bool silent = false ) const
	{

		brgastro::vector<double> xlo, xhi;
		brgastro::vector<unsigned int> x_i; // Lower nearby array points
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

		x_i = ( x - SPCP(name)->_mins_ ) / SPCP(name)->_steps_ ;
		x_i = brgastro::max( x_i, 0 );
		x_i = brgastro::min( x_i, SPCP(name)->_resolutions_ - 2 );

		xlo = SPCP(name)->_mins_ + SPCP(name)->_steps_ * x_i;
		xhi = SPCP(name)->_mins_ + SPCP(name)->_steps_ * ( x_i + 1 );

		const unsigned int num_surrounding_points = std::pow(2,SPCP(name)->_num_dim_);
		std::vector< std::vector<bool> > use_high(SPCP(name)->_num_dim_);
		for(unsigned int i=0; i < SPCP(name)->_num_dim_; i++)
		{
			int divisor = std::pow(2,i);
			use_high[i].resize(num_surrounding_points);
			for(unsigned int j=0; j<num_surrounding_points; j++)
			{
				use_high[i][j] = divisible(j/divisor,2);
			}
		}

		result = 0;
		double total_weight = 0;
		brgastro::vector<unsigned int> position(SPCP(name)->_num_dim_,0);

		for(unsigned int j=0; j < num_surrounding_points; j++)
		{
			double weight = 1;
			for(unsigned int i=0; i < SPCP(name)->_num_dim_; i++)
			{
				if(use_high[i][j])
				{
					position[i] = x_i[i]+1;
					weight *= x[i]-xlo[i];
				}
				else
				{
					position[i] = x_i[i];
					weight *= xhi[i]-x[i];
				}
			}
			result += SPCP(name)->_results_.at(position) * weight;
			total_weight += weight;
		}

		result /= safe_d(total_weight);

		return result;

	} // get()

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
	brg_cache_nd() throw()
	{
		if(!SPCP(name)->_initialised_) SPP(name)->_init();
	}

	// Deconstructor
	virtual ~brg_cache_nd()
	{
	}

#endif // Public methods

}; // class brg_cache

} // namespace brgastro

// Initialisation of static members
#if (1)

// These ones must be specialised for default values. The zeros in the initialisers and templates
// should be replaced with the number of dimensions.
template<typename name>
brgastro::vector<double> brgastro::brg_cache_nd<name>::_mins_ = {0};

template<typename name>
brgastro::vector<double> brgastro::brg_cache_nd<name>::_maxes_ = {0};

template<typename name>
brgastro::vector<double> brgastro::brg_cache_nd<name>::_steps_ = {2}; // Never let any of the steps be zero

// These ones don't need to be altered for each child, but they must be copied

template<typename name>
brgastro::vector<double> brgastro::brg_cache_nd<name>::_results_;

template<typename name>
bool brgastro::brg_cache_nd<name>::_loaded_ = false;

template<typename name>
bool brgastro::brg_cache_nd<name>::_initialised_ = false;

template<typename name>
brgastro::vector<unsigned int> brgastro::brg_cache_nd<name>::_resolutions_ = {0};

template<typename name>
std::string brgastro::brg_cache_nd<name>::_file_name_ = "";

template<typename name>
unsigned int brgastro::brg_cache_nd<name>::_version_number_ = 0;

#endif

// Undef macros
#undef SPP
#undef SPCP

#endif // __BRG_CACHE_ND_HPP_INCLUDED__
