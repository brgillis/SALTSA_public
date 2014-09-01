/**********************************************************************\
  @file gabdt.cpp

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

#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "brg/global.h"

#include "brg/math/calculus/differentiate.hpp"
#include "brg/physics/astro/density_profile/density_profile.h"
#include "brg/physics/units/unit_obj.h"
#include "brg/utility.hpp"

#include "gabdt.h"

// brgastro::gabdt class method implementations
#if (1)

// Swap functions
void brgastro::gabdt::swap(gabdt &other)
{
	using std::swap;
	swap(_is_cached_, other._is_cached_);
	swap(_host_ptr_, other._host_ptr_);

	swap(_x_, other._x_);
	swap(_y_, other._y_);
	swap(_z_, other._z_);
	swap(_r_, other._r_);
	swap(_dt_, other._dt_);

	swap(_dv_, other._dv_);
}
namespace std
{
	template <>
	void swap(brgastro::gabdt &same,
			brgastro::gabdt &other)
	{
		same.swap(other);
	}
}

brgastro::gabdt::gabdt( void )
{
	clear();
}

brgastro::gabdt::gabdt( const gabdt & other )
{
	_is_cached_ = other._is_cached_;
	_host_ptr_ = other._host_ptr_;

	_x_ = other._x_;
	_y_ = other._y_;
	_z_ = other._z_;
	_r_ = other._r_;
	_dt_ = other._dt_;

	_dv_ = other._dv_;
}

brgastro::gabdt::gabdt( const density_profile *new_host,
		CONST_BRG_DISTANCE_REF new_x, CONST_BRG_DISTANCE_REF new_y,
		CONST_BRG_DISTANCE_REF new_z, CONST_BRG_TIME_REF new_dt )
{
	set( new_host, new_x, new_y, new_z, new_dt );
}

brgastro::gabdt::~gabdt( void )
{
}

const int brgastro::gabdt::clear()
{
	_is_cached_ = false;
	_x_ = 0;
	_y_ = 0;
	_z_ = 0;
	_r_ = 0;
	_dt_ = 0;
	_dv_.clear();
	return 0;
}

const int brgastro::gabdt::set( const density_profile *new_host,
		CONST_BRG_DISTANCE_REF new_x, CONST_BRG_DISTANCE_REF new_y,
		CONST_BRG_DISTANCE_REF new_z, CONST_BRG_TIME_REF new_dt )
{
	_is_cached_ = false;
	_host_ptr_ = new_host;
	_dt_ = new_dt;
	_dv_.clear();
	return 0;
}
const int brgastro::gabdt::set_host_ptr( const density_profile *new_host )
{
	_is_cached_ = false;
	_host_ptr_ = new_host;
	return 0;
}
const int brgastro::gabdt::set_pos( CONST_BRG_DISTANCE_REF new_x,
		CONST_BRG_DISTANCE_REF new_y, CONST_BRG_DISTANCE_REF new_z )
{
	_is_cached_ = false;
	_x_ = new_x;
	_y_ = new_y;
	_z_ = new_z;
	_r_ = dist3d( _x_, _y_, _z_ );
	_dv_.clear();
	return 0;
}
const int brgastro::gabdt::set_dt( CONST_BRG_TIME_REF new_dt )
{
	if ( _is_cached_ )
	{
		// Multiply dv by ratio of new_dt to old
		for ( int i = 0; i < 3; i++ )
			for ( int j = 0; j < 3; j++ )
			{
				_dv_[i][j] *= new_dt / _dt_;
			}
	}
	_dt_ = new_dt;
	return 0;
}
const int brgastro::gabdt::calc_dv( const bool silent ) const
{
	gabdt_functor gabdtf;
	gabdtf.host_ptr = _host_ptr_;
	unsigned int num_in_params = 3, num_out_params = 3;
	std::vector< BRG_UNITS > in_params( num_in_params, 0 ), out_params(
			num_out_params, 0 );
	std::vector< std::vector< BRG_UNITS > > Jacobian;
	make_array2d( _dv_, num_out_params, num_in_params );

	in_params[0] = _x_;
	in_params[1] = _y_;
	in_params[2] = _z_;

	try
	{
		Jacobian = differentiate( &gabdtf, in_params);
	}
	catch(const std::exception &e)
	{
		if ( !silent )
		{
			std::cerr << "ERROR: Cannot differentiate in gabdt::calc_dv():\n";
			std::cerr << e.what() << std::endl;
			std::cerr.flush();
		}
		for ( unsigned int i = 0; i < num_out_params; i++ )
			for ( unsigned int j = 0; j < num_in_params; j++ )
				_dv_[i][j] = 0; // To be safe;
		_is_cached_ = true;
		return 1;
	}

	// Multiply by dt
	bool bad_value = false;
	for ( unsigned int i = 0; i < num_out_params; i++ )
		for ( unsigned int j = 0; j < num_in_params; j++ )
		{
			if(isbad(Jacobian[i][j]))
			{
				_dv_[i][j] = 0;
				bad_value = true;
			}
			else;
			{
				_dv_[i][j] = Jacobian[i][j]*_dt_;
			}
		}

	if(bad_value && (!silent))
	{
		std::cerr << "WARNING: Bad value in Jacobian in calculation of gabdt. Treating as zero to be safe.\n";
		std::cerr.flush();
	}

	_is_cached_ = true;

	return 0;
}
const int brgastro::gabdt::override_zero()
{
	make_array2d( _dv_, 3, 3 );
	for ( int x_i = 0; x_i < 3; x_i++ )
	{
		for ( int y_i = 0; y_i < 3; y_i++ )
		{
			_dv_[x_i][y_i] = 0;
		}
	}
	_is_cached_ = true;
	return 0;
}
const brgastro::density_profile * brgastro::gabdt::host() const
{
	return _host_ptr_;
}
const BRG_DISTANCE brgastro::gabdt::x() const
{
	return _x_;
}
const BRG_DISTANCE brgastro::gabdt::y() const
{
	return _y_;
}
const BRG_DISTANCE brgastro::gabdt::z() const
{
	return _z_;
}
const BRG_DISTANCE brgastro::gabdt::r() const
{
	return _r_;
}
const std::vector< std::vector< long double > > brgastro::gabdt::dv() const
{
	if ( !_is_cached_ )
	{
		if ( calc_dv() )
			_dv_.clear();
	}
	return _dv_;
}
const long double brgastro::gabdt::dv( const int x_i, const int y_i ) const
{
	return dv().at(x_i).at(y_i);
}

brgastro::gabdt & brgastro::gabdt::operator=( gabdt other )
{
	swap(other);
	return *this;
}

const BRG_UNITS brgastro::gabdt::operator*( const gabdt & other_gabdt ) const // "Dot-product" operator
{
	if ( !_is_cached_ )
		if ( calc_dv() )
			_dv_.clear();
	double result = 0;
	for ( int x_i = 0; x_i < 3; x_i++ )
	{
		for ( int y_i = 0; y_i < 3; y_i++ )
		{
			result += dv().at(x_i).at(y_i)*other_gabdt.dv().at(x_i).at(y_i);
			if(isbad(result)) std::cerr << "WARNING: Bad result in gabdt dot-product when multiplying "
					<< dv().at(x_i).at(y_i) << " and " << other_gabdt.dv().at(x_i).at(y_i) << std::endl;
		}
	}
	return result;
}
brgastro::gabdt & brgastro::gabdt::operator+=( const gabdt & other_gabdt )
{
	if ( !_is_cached_ )
		if ( calc_dv() )
		{
			_dv_.clear();
			return *this;
		}

	for ( int x_i = 0; x_i < 3; x_i++ )
	{
		for ( int y_i = 0; y_i < 3; y_i++ )
		{
			long double tmp = _dv_.at(x_i).at(y_i);
			_dv_.at(x_i).at(y_i) += other_gabdt.dv().at(x_i).at(y_i);
			if(isbad(_dv_[x_i][y_i])) std::cerr << "WARNING: Bad result in gabdt dot-product when multiplying "
					<< tmp << " and " << other_gabdt.dv().at(x_i).at(y_i) << std::endl;
		}
	}
	return *this;
}
brgastro::gabdt brgastro::gabdt::operator+( const gabdt & other_gabdt ) const
{
	if ( !_is_cached_ )
		if ( calc_dv() )
		{
			_dv_.clear();
			return *this;
		}
	gabdt result = gabdt( *this );
	result += other_gabdt;
	return result;
}

brgastro::gabdt & brgastro::gabdt::operator*=( const double scale_fraction )
{
	if(isbad(scale_fraction))
	{
		std::cerr << "WARNING: Bad scale fraction passed to gabdt*=: " << scale_fraction << std::endl;
		return *this;
	}
	if ( !_is_cached_ )
		if ( calc_dv() )
		{
			_dv_.clear();
			return *this;
		}

	for ( int x_i = 0; x_i < 3; x_i++ )
	{
		for ( int y_i = 0; y_i < 3; y_i++ )
		{
			_dv_[x_i][y_i] *= scale_fraction;
		}
	}
	return *this;

}

brgastro::gabdt brgastro::gabdt::operator*( const double scale_fraction ) const
{
	if ( !_is_cached_ )
		if ( calc_dv() )
		{
			_dv_.clear();
			return *this;
		}
	gabdt result = gabdt( *this );
	result *= scale_fraction;
	return result;
}

#endif // end brgastro::gabdt class function definitions

// brgastro::gabdt_functor class method implementations
#if (1)

// Swap functions
void brgastro::gabdt_functor::swap(gabdt_functor &other)
{
	using std::swap;
	swap(host_ptr,other.host_ptr);
}
namespace std
{
	template <>
	void swap(brgastro::gabdt_functor &same,
			brgastro::gabdt_functor &other)
	{
		same.swap(other);
	}
}

// Constructors
brgastro::gabdt_functor::gabdt_functor()
{
	host_ptr = NULL;
}
brgastro::gabdt_functor::gabdt_functor(const gabdt_functor &other)
{
	host_ptr = other.host_ptr;
}

// Operator=
brgastro::gabdt_functor & brgastro::gabdt_functor::operator=(gabdt_functor other)
{
	swap(other);
	return *this;
}

// Operator()
std::vector< BRG_UNITS > brgastro::gabdt_functor::operator()(
		const std::vector< BRG_UNITS > & in_params, const bool silent ) const
{
	std::vector< BRG_UNITS > out_params(3,0);

	if ( in_params.size() != 3 )
	{
		throw std::runtime_error("ERROR: in_params.size() and num_in_params must == 3 in gabdt_functor.\n");
	}
	double R = dist3d( in_params[0], in_params[1], in_params[2] );

	unsigned int num_out_params = 3;
	make_array( out_params, num_out_params );
	if ( R == 0 )
	{
		// Special handling for this case, returning a zero result
		for ( unsigned int i = 0; i < num_out_params; i++ )
		{
			out_params[i] = 0;
		}

	}
	else
	{
		double accel_mag = host_ptr->accel( R );

		for ( unsigned int i = 0; i < num_out_params; i++ )
		{
			out_params[i] = in_params[i] / R * accel_mag;
		}
	}
	return out_params;
}

#endif // end brgastro::gabdt_functor class function definitions
