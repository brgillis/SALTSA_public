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
#include "brg/physics/density_profile/density_profile.h"
#include "brg/physics/units/unit_obj.h"
#include "brg/utility.hpp"

#include "gabdt.h"

// brgastro::gabdt class method implementations
#if (1)

brgastro::gabdt::gabdt( void )
{
	clear();
}

brgastro::gabdt::gabdt( const density_profile *new_host,
		CONST_BRG_DISTANCE_REF new_x, CONST_BRG_DISTANCE_REF new_y,
		CONST_BRG_DISTANCE_REF new_z, CONST_BRG_TIME_REF new_dt )
{
	set( new_host, new_x, new_y, new_z, new_dt );
}

void brgastro::gabdt::set( const density_profile *new_host,
		CONST_BRG_DISTANCE_REF new_x, CONST_BRG_DISTANCE_REF new_y,
		CONST_BRG_DISTANCE_REF new_z, CONST_BRG_TIME_REF new_dt )
{
	_is_cached_ = false;
	_host_ptr_ = new_host;
	_dt_ = new_dt;
	_dv_.clear();
}
void brgastro::gabdt::clear()
{
	_is_cached_ = false;
	_host_ptr_ =  NULL;

	_x_ = _y_ = _z_ = _r_ = 0;
	_dt_ = 0;
	_dv_.clear();
}
void brgastro::gabdt::set_host_ptr( const density_profile *new_host )
{
	_is_cached_ = false;
	_host_ptr_ = new_host;
}
void brgastro::gabdt::set_pos( CONST_BRG_DISTANCE_REF new_x,
		CONST_BRG_DISTANCE_REF new_y, CONST_BRG_DISTANCE_REF new_z )
{
	_is_cached_ = false;
	_x_ = new_x;
	_y_ = new_y;
	_z_ = new_z;
	_r_ = dist3d( _x_, _y_, _z_ );
	_dv_.clear();
}
void brgastro::gabdt::set_dt( CONST_BRG_TIME_REF new_dt )
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
}
void brgastro::gabdt::calc_dv( const bool silent ) const
{
	gabdt_functor gabdtf;
	gabdtf.host_ptr = _host_ptr_;
	size_t num_in_params = 3, num_out_params = 3;
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
		}
		for ( size_t i = 0; i < num_out_params; i++ )
			for ( size_t j = 0; j < num_in_params; j++ )
				_dv_[i][j] = 0; // To be safe;
		_is_cached_ = true;
		throw;
	}

	// Multiply by dt
	bool bad_value = false;
	for ( size_t i = 0; i < num_out_params; i++ )
		for ( size_t j = 0; j < num_in_params; j++ )
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
}
void brgastro::gabdt::override_zero()
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
}
const brgastro::density_profile * brgastro::gabdt::host() const
{
	return _host_ptr_;
}
BRG_DISTANCE brgastro::gabdt::x() const
{
	return _x_;
}
BRG_DISTANCE brgastro::gabdt::y() const
{
	return _y_;
}
BRG_DISTANCE brgastro::gabdt::z() const
{
	return _z_;
}
BRG_DISTANCE brgastro::gabdt::r() const
{
	return _r_;
}
std::vector< std::vector< long double > > brgastro::gabdt::dv() const
{
	if ( !_is_cached_ )
	{
		calc_dv();
	}
	return _dv_;
}
long double brgastro::gabdt::dv( const int x_i, const int y_i ) const
{
	return dv().at(x_i).at(y_i);
}

BRG_UNITS brgastro::gabdt::operator*( const gabdt & other_gabdt ) const // "Dot-product" operator
{
	if ( !_is_cached_ )
		calc_dv();
	double result = 0;
	for ( size_t x_i = 0; x_i < 3; x_i++ )
	{
		for ( size_t y_i = 0; y_i < 3; y_i++ )
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
		calc_dv();

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
		calc_dv();
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
		calc_dv();

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
		calc_dv();
	gabdt result = gabdt( *this );
	result *= scale_fraction;
	return result;
}

#endif // end brgastro::gabdt class function definitions

// brgastro::gabdt_functor class method implementations
#if (1)

// Constructors
brgastro::gabdt_functor::gabdt_functor()
{
	host_ptr = NULL;
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

	size_t num_out_params = 3;
	make_array( out_params, num_out_params );
	if ( R == 0 )
	{
		// Special handling for this case, returning a zero result
		for ( size_t i = 0; i < num_out_params; i++ )
		{
			out_params[i] = 0;
		}

	}
	else
	{
		double accel_mag = host_ptr->accel( R );

		for ( size_t i = 0; i < num_out_params; i++ )
		{
			out_params[i] = in_params[i] / R * accel_mag;
		}
	}
	return out_params;
}

#endif // end brgastro::gabdt_functor class function definitions
