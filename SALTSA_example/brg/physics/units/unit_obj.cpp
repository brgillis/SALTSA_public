/**********************************************************************\
  @file unit_obj.cpp

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

#include "brg/global.h"

#ifdef _BRG_USE_UNITS_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <sstream>
#include <memory>

#include "brg/physics/units/unit_conversions.hpp"

#include "unit_obj.h"

using namespace unitconv;
using std::cout;
using std::cerr;
using std::string;
using std::stringstream;

//-------------------------
// Full functions (brgastro::unit_obj)
//-------------------------

brgastro::unit_obj::unit_obj(const double init__value_, const float d_units_power, const float t_units_power, const float m_units_power,
		const float T_units_power, const float a_units_power, const float c_units_power)
{
	_fix_unit_powers_ = false;
	_unit_powers_.resize(NUM_UNIT_TYPES,0);

	set(init__value_,d_units_power,t_units_power,m_units_power,
			T_units_power,a_units_power,c_units_power);
}
brgastro::unit_obj::unit_obj(const double init__value_,
		const double d_units, const float d_units_power, const double t_units, const float t_units_power,
		const double m_units, const float m_units_power, const double T_units, const float T_units_power,
		const double a_units, const float a_units_power, const double c_units, const float c_units_power)
{
	_fix_unit_powers_ = false;
	_unit_powers_.resize(NUM_UNIT_TYPES,0);

	set(init__value_,d_units,d_units_power,t_units,t_units_power,m_units,m_units_power,
			T_units,T_units_power,a_units,a_units_power,c_units,c_units_power);
}
brgastro::unit_obj::~unit_obj()
{
	_unit_powers_.clear(); // Not strictly necessary with a vector, but good practice in case it's redefined
}

const int brgastro::unit_obj::round_powers()
{
	if(_fix_unit_powers_) return 0;
	for(int i=0; i < NUM_UNIT_TYPES; i++)
	{
		_unit_powers_[i] = round_int(_unit_powers_[i]*UNIT_POWER_ROUND_PRECISION_FACTOR)/(float)UNIT_POWER_ROUND_PRECISION_FACTOR;
	}

	return 0;
}

// Reset variable and set units (distance, time, mass, temperature, angle, all taken as doubles)
const int brgastro::unit_obj::reset(const float dp, const float tp, const float mp, const float Tp, const float ap, const float cp)
{
	_value_ = 0;
	if(_fix_unit_powers_) return 0;
	_unit_powers_[DIST_UNIT_INDEX] = dp;
	_unit_powers_[TIME_UNIT_INDEX] = tp;
	_unit_powers_[MASS_UNIT_INDEX] = mp;
	_unit_powers_[TEMP_UNIT_INDEX] = Tp;
	_unit_powers_[ANGL_UNIT_INDEX] = ap;
	_unit_powers_[CHRG_UNIT_INDEX] = cp;

	return 0;
}
// Set units, but don't change _value_ue in default units - use with care!
const int brgastro::unit_obj::set_unit_powers(const float dp, const float tp, const float mp, const float Tp, const float ap, const float cp)
{
	if(!_fix_unit_powers_)
	{
		_unit_powers_[DIST_UNIT_INDEX] = dp;
		_unit_powers_[TIME_UNIT_INDEX] = tp;
		_unit_powers_[MASS_UNIT_INDEX] = mp;
		_unit_powers_[TEMP_UNIT_INDEX] = Tp;
		_unit_powers_[ANGL_UNIT_INDEX] = ap;
		_unit_powers_[CHRG_UNIT_INDEX] = cp;
		return 0;
	}
	else
	{
#ifdef _BRG_WARN_FOR_UNIT_MISMATCH_
		std::cerr << "WARNING: Attempt to change unit powers of brgastro::unit_obj with fixed unit powers.\n";
#endif
		return INVALID_ARGUMENTS_ERROR;
	}
}
// Set units, but don't change _value_ue in default units - use with care!
const int brgastro::unit_obj::set_unit_powers(std::vector<float> new__unit_powers_)
{
	if(_fix_unit_powers_==false)
	{
		if(new__unit_powers_.size()==NUM_UNIT_TYPES)
		{
			_unit_powers_[DIST_UNIT_INDEX] = new__unit_powers_[DIST_UNIT_INDEX];
			_unit_powers_[TIME_UNIT_INDEX] = new__unit_powers_[TIME_UNIT_INDEX];
			_unit_powers_[MASS_UNIT_INDEX] = new__unit_powers_[MASS_UNIT_INDEX];
			_unit_powers_[TEMP_UNIT_INDEX] = new__unit_powers_[TEMP_UNIT_INDEX];
			_unit_powers_[ANGL_UNIT_INDEX] = new__unit_powers_[ANGL_UNIT_INDEX];
			_unit_powers_[CHRG_UNIT_INDEX] = new__unit_powers_[CHRG_UNIT_INDEX];
		}
		else
		{
			return INVALID_ARGUMENTS_ERROR;
		}
	}
	else
	{
#ifdef _BRG_WARN_FOR_UNIT_MISMATCH_
		std::cerr << "WARNING: Attempt to change unit powers of brgastro::unit_obj with fixed unit powers.\n";
#endif
		return INVALID_ARGUMENTS_ERROR;
	}
	return 0;
}
// Set _value_ue in default or specified units
const int brgastro::unit_obj::set_value(const double init__value_, const double conv_factor)
{
	_value_ = init__value_*conv_factor;
	return 0;
}
// Set _value_ue in specified units, but don't change units
// IMPORTANT: Always use the unitconv::pctom (units used here to default) form for units)
// EXAMPLE 1: To load 4 kpc, use set_value_ue(4, kpctom, 0, 0, 0, 0, 0)
// EXAMPLE 2: To load 1 m/s^2, use set_value_ue(1, mtom, stos, 0, 0, 0, 0) (the function will automatically figure out powers from the unit)
const int brgastro::unit_obj::set_value(const double new__value_, const double d_units, const double t_units, const double m_units, const double T_units,
		const double a_units, const double c_units)
{
	// Check if any units were given as 0. If so, change them to 1 so we can multiply
	double d_units_use = ( d_units == 0 ? 1 : d_units);
	double t_units_use = ( t_units == 0 ? 1 : t_units);
	double m_units_use = ( m_units == 0 ? 1 : m_units);
	double T_units_use = ( T_units == 0 ? 1 : T_units);
	double a_units_use = ( a_units == 0 ? 1 : a_units);
	double c_units_use = ( c_units == 0 ? 1 : c_units);
	_value_ = new__value_ * std::pow(d_units_use,_unit_powers_[DIST_UNIT_INDEX]) * std::pow(t_units_use,_unit_powers_[TIME_UNIT_INDEX]) *
	std::pow(m_units_use,_unit_powers_[MASS_UNIT_INDEX]) * std::pow(T_units_use,_unit_powers_[TEMP_UNIT_INDEX]) * std::pow(a_units_use,_unit_powers_[ANGL_UNIT_INDEX]) *
	std::pow(c_units_use,_unit_powers_[CHRG_UNIT_INDEX]);
	return 0;
}

// Set _value_ue and changes units - used by overloaded operators mostly
const int brgastro::unit_obj::set(const double new__value_, const std::vector<float> new__unit_powers_)
{
	if(new__unit_powers_.size() != NUM_UNIT_TYPES) return INVALID_ARGUMENTS_ERROR;
	_value_ = new__value_; // We allow the _value_ue to be changed even if there's a problem with changing unit powers, so it will behave as expected
	if(_fix_unit_powers_==false)
	{
		_unit_powers_[DIST_UNIT_INDEX] = new__unit_powers_[DIST_UNIT_INDEX];
		_unit_powers_[TIME_UNIT_INDEX] = new__unit_powers_[TIME_UNIT_INDEX];
		_unit_powers_[MASS_UNIT_INDEX] = new__unit_powers_[MASS_UNIT_INDEX];
		_unit_powers_[TEMP_UNIT_INDEX] = new__unit_powers_[TEMP_UNIT_INDEX];
		_unit_powers_[ANGL_UNIT_INDEX] = new__unit_powers_[ANGL_UNIT_INDEX];
		_unit_powers_[CHRG_UNIT_INDEX] = new__unit_powers_[CHRG_UNIT_INDEX];
	}
	else
	{
		// Check if any unit powers are actually being changed here
		if((_unit_powers_[DIST_UNIT_INDEX] != new__unit_powers_[DIST_UNIT_INDEX]) ||
				(_unit_powers_[TIME_UNIT_INDEX] != new__unit_powers_[TIME_UNIT_INDEX]) ||
				(_unit_powers_[MASS_UNIT_INDEX] != new__unit_powers_[MASS_UNIT_INDEX]) ||
				(_unit_powers_[TEMP_UNIT_INDEX] != new__unit_powers_[TEMP_UNIT_INDEX]) ||
				(_unit_powers_[ANGL_UNIT_INDEX] != new__unit_powers_[ANGL_UNIT_INDEX]) ||
				(_unit_powers_[CHRG_UNIT_INDEX] != new__unit_powers_[CHRG_UNIT_INDEX]))
		{
			std::cerr << "ERROR: Attempt to change unit powers of brgastro::unit_obj with fixed unit powers.\n";
			return INVALID_ARGUMENTS_ERROR;
		}
	}
	return 0;
}
const int brgastro::unit_obj::set(const double new__value_, const float dp, const float tp, const float mp, const float Tp, const float ap, const float cp)
{
	return set(new__value_, 1, dp, 1, tp, 1, mp, 1, Tp, 1, ap, 1, cp);
}

// Set _value_ue in specified units, and change units of variable
// IMPORTANT: Always use the unitconv::pctom (units used here to default) form for units)
// EXAMPLE 1: To load 4 kpc, use set_value_ue(4, kpctom, 1, 0, 0, 0, 0, 0, 0, 0, 0)
// EXAMPLE 2: To load 1 m/s^2, use set_value_ue(1, mtom, 1, stos, -2, 0, 0, 0, 0, 0, 0)
const int brgastro::unit_obj::set(const double new__value_, const double d_units, const float d_units_power, const double t_units,
		const float t_units_power, const double m_units, const float m_units_power, const double T_units,
		const float T_units_power, const double a_units, const float a_units_power, const double c_units,
		const float c_units_power)
{
	int return__value_ue = 0;

	if(_fix_unit_powers_==false)
	{
		_unit_powers_[DIST_UNIT_INDEX] = d_units_power;
		_unit_powers_[TIME_UNIT_INDEX] = t_units_power;
		_unit_powers_[MASS_UNIT_INDEX] = m_units_power;
		_unit_powers_[TEMP_UNIT_INDEX] = T_units_power;
		_unit_powers_[ANGL_UNIT_INDEX] = a_units_power;
		_unit_powers_[CHRG_UNIT_INDEX] = c_units_power;
	}
#ifdef _BRG_WARN_FOR_UNIT_MISMATCH_
	else
	{

		// Check if any unit powers are actually being changed here
		if((_unit_powers_[DIST_UNIT_INDEX] != d_units_power) ||
				(_unit_powers_[TIME_UNIT_INDEX] != t_units_power) ||
				(_unit_powers_[MASS_UNIT_INDEX] != m_units_power) ||
				(_unit_powers_[TEMP_UNIT_INDEX] != T_units_power) ||
				(_unit_powers_[ANGL_UNIT_INDEX] != a_units_power) ||
				(_unit_powers_[CHRG_UNIT_INDEX] != c_units_power))
		{
			std::cerr << "ERROR: Attempt to change unit powers of brgastro::unit_obj with fixed unit powers.\n";
			return__value_ue = 1;
		}
	}
#endif

	// Check if any units were given as 0. If so, change them to 1 so we can multiply safely
	double d_units_use = ( d_units == 0 ? 1 : d_units);
	double t_units_use = ( t_units == 0 ? 1 : t_units);
	double m_units_use = ( m_units == 0 ? 1 : m_units);
	double T_units_use = ( T_units == 0 ? 1 : T_units);
	double a_units_use = ( a_units == 0 ? 1 : a_units);
	double c_units_use = ( c_units == 0 ? 1 : c_units);
	_value_ = new__value_ * std::pow(d_units_use,_unit_powers_[DIST_UNIT_INDEX]) * std::pow(t_units_use,_unit_powers_[TIME_UNIT_INDEX]) *
	std::pow(m_units_use,_unit_powers_[MASS_UNIT_INDEX]) * std::pow(T_units_use,_unit_powers_[TEMP_UNIT_INDEX]) * std::pow(a_units_use,_unit_powers_[ANGL_UNIT_INDEX]) *
	std::pow(c_units_use,_unit_powers_[CHRG_UNIT_INDEX]);
	return return__value_ue;
}

// Get _value_ue in default units
const double brgastro::unit_obj::get_value() const
{
	return _value_;
}

const double brgastro::unit_obj::get_value(const double d_units, const double t_units, const double m_units, const double T_units, const double a_units,
		const double c_units) const
{
	// Get _value_ue in specified units
	// IMPORTANT: Always use the unitconv::pctom (units used here to default) form for units)
	// EXAMPLE 1: To load 4 kpc, use set_value_ue(4, kpctom, 0, 0, 0, 0)
	// EXAMPLE 2: To load 1 m/s^2, use set_value_ue(1, mtom, stos, 0, 0, 0) (the function will automatically figure out powers from the unit)

	// Check if any units were given as 0. If so, change them to 1 so we can multiply
	double d_units_use = ( d_units == 0 ? 1 : d_units);
	double t_units_use = ( t_units == 0 ? 1 : t_units);
	double m_units_use = ( m_units == 0 ? 1 : m_units);
	double T_units_use = ( T_units == 0 ? 1 : T_units);
	double a_units_use = ( a_units == 0 ? 1 : a_units);
	double c_units_use = ( c_units == 0 ? 1 : c_units);
	return(_value_ * std::pow(d_units_use,-_unit_powers_[DIST_UNIT_INDEX]) * std::pow(t_units_use,-_unit_powers_[TIME_UNIT_INDEX]) *
			std::pow(m_units_use,-_unit_powers_[MASS_UNIT_INDEX]) * std::pow(T_units_use,-_unit_powers_[TEMP_UNIT_INDEX]) * std::pow(a_units_use,-_unit_powers_[ANGL_UNIT_INDEX]) *
			std::pow(c_units_use,-_unit_powers_[CHRG_UNIT_INDEX]));
	// Note the negative sign for powers here, to flip the sign on the unit conversion for output
}

const double brgastro::unit_obj::get_value(const double d_units, const float d_units_power, const double t_units, const float t_units_power,
		const double m_units, const float m_units_power, const double T_units, const float T_units_power,
		const double a_units, const float a_units_power, const double c_units, const float c_units_power) const
{
	// Get _value_ue in specified units, overriding units of variable
	// IMPORTANT: Always use the unitconv::mtopc (units used here to default) form for units)
	// EXAMPLE 1: To load 4 kpc, use set_value_ue(4, kpctom, 1, 0, 0, 0, 0, 0, 0, 0, 0)
	// EXAMPLE 2: To load 1 m/s^2, use set_value_ue(1, mtom, 1, stos, -2, 0, 0, 0, 0, 0, 0)

	// Check if any units were given as 0. If so, change them to 1 so we can multiply
	double d_units_use = ( d_units == 0 ? 1 : d_units);
	double t_units_use = ( t_units == 0 ? 1 : t_units);
	double m_units_use = ( m_units == 0 ? 1 : m_units);
	double T_units_use = ( T_units == 0 ? 1 : T_units);
	double a_units_use = ( a_units == 0 ? 1 : a_units);
	double c_units_use = ( c_units == 0 ? 1 : c_units);
	return(_value_ * std::pow(d_units_use,-d_units_power) * std::pow(t_units_use,-t_units_power) *
			std::pow(m_units_use,-m_units_power) * std::pow(T_units_use,-T_units_power) * std::pow(a_units_use,-a_units_power) *
			std::pow(c_units_use,-c_units_power));
	// Note the negative sign for powers here, to flip the sign on the unit conversion for output
}

const std::vector<float> brgastro::unit_obj::get_unit_powers() const
{
	// Returns a 6-component std::vector<float> in the format [d_units_power,t_units_power,m_units_power,T_units_power,a_units_power,c_units_power]
	std::vector<float> _unit_powers_copy(_unit_powers_);
	return _unit_powers_copy;
}

const std::string brgastro::unit_obj::get_string() const
{
	// Returns a std::string detailing the _value_ue and units
	// Good for debugging purposes or output
	stringstream ss;
	bool firstunitflag = true;
	ss.str("");
	ss << _value_;
	if(_unit_powers_[DIST_UNIT_INDEX] != 0)
	{
		if(firstunitflag) firstunitflag = false;
		else ss << " ";
		if(_unit_powers_[DIST_UNIT_INDEX]==1)
		{
			ss << "m";
		}
		else
		{
			ss << "m^" << _unit_powers_[DIST_UNIT_INDEX];
		}
	}
	if(_unit_powers_[TIME_UNIT_INDEX] != 0)
	{
		if(firstunitflag) firstunitflag = false;
		else ss << " ";
		if(_unit_powers_[TIME_UNIT_INDEX]==1)
		{
			ss << "s";
		}
		else
		{
			ss << "s^" << _unit_powers_[TIME_UNIT_INDEX];
		}
	}
	if(_unit_powers_[MASS_UNIT_INDEX] != 0)
	{
		if(firstunitflag) firstunitflag = false;
		else ss << " ";
		if(_unit_powers_[MASS_UNIT_INDEX]==1)
		{
			ss << "kg";
		}
		else
		{
			ss << "kg^" << _unit_powers_[MASS_UNIT_INDEX];
		}
	}
	if(_unit_powers_[TEMP_UNIT_INDEX] != 0)
	{
		if(firstunitflag) firstunitflag = false;
		else ss << " ";
		if(_unit_powers_[TEMP_UNIT_INDEX]==1)
		{
			ss << "K";
		}
		else
		{
			ss << "K^" << _unit_powers_[TEMP_UNIT_INDEX];
		}
	}
	if(_unit_powers_[ANGL_UNIT_INDEX] != 0)
	{
		if(firstunitflag) firstunitflag = false;
		else ss << " ";
		if(_unit_powers_[ANGL_UNIT_INDEX]==1)
		{
			ss << "rad";
		}
		else
		{
			ss << "rad^" << _unit_powers_[ANGL_UNIT_INDEX];
		}
	}
	if(_unit_powers_[CHRG_UNIT_INDEX] != 0)
	{
		if(firstunitflag) firstunitflag = false;
		else ss << " ";
		if(_unit_powers_[CHRG_UNIT_INDEX]==1)
		{
			ss << "C";
		}
		else
		{
			ss << "C^" << _unit_powers_[CHRG_UNIT_INDEX];
		}
	}
	return ss.str();
}

// distance Units
const double brgastro::unit_obj::m() const
{
	return _value_;
}

const double brgastro::unit_obj::mm() const
{
	return _value_*std::pow(unitconv::mtomm,_unit_powers_[DIST_UNIT_INDEX]);
}

const double brgastro::unit_obj::um() const
{
	return _value_*std::pow(unitconv::mtoum,_unit_powers_[DIST_UNIT_INDEX]);
}

const double brgastro::unit_obj::nm() const
{
	return _value_*std::pow(unitconv::mtonm,_unit_powers_[DIST_UNIT_INDEX]);
}

const double brgastro::unit_obj::cm() const
{
	return _value_*std::pow(unitconv::mtocm,_unit_powers_[DIST_UNIT_INDEX]);
}

const double brgastro::unit_obj::angstrom() const
{
	return _value_*std::pow(unitconv::mtoangstrom,_unit_powers_[DIST_UNIT_INDEX]);
}

const double brgastro::unit_obj::km() const
{
	return _value_*std::pow(unitconv::mtokm,_unit_powers_[DIST_UNIT_INDEX]);
}

const double brgastro::unit_obj::ltyr() const
{
	return _value_*std::pow(unitconv::mtoltyr,_unit_powers_[DIST_UNIT_INDEX]);
}

const double brgastro::unit_obj::pc() const
{
	return _value_*std::pow(unitconv::mtopc,_unit_powers_[DIST_UNIT_INDEX]);
}

const double brgastro::unit_obj::kpc() const
{
	return _value_*std::pow(unitconv::mtokpc,_unit_powers_[DIST_UNIT_INDEX]);
}

const double brgastro::unit_obj::Mpc() const
{
	return _value_*std::pow(unitconv::mtoMpc,_unit_powers_[DIST_UNIT_INDEX]);
}

const double brgastro::unit_obj::mi() const
{
	return _value_*std::pow(unitconv::mtomi,_unit_powers_[DIST_UNIT_INDEX]);
}

const double brgastro::unit_obj::Mmi() const
{
	return _value_*std::pow(unitconv::mtoMmi,_unit_powers_[DIST_UNIT_INDEX]);
}

const double brgastro::unit_obj::ft() const
{
	return _value_*std::pow(unitconv::mtoft,_unit_powers_[DIST_UNIT_INDEX]);
}

const double brgastro::unit_obj::yd() const
{
	return _value_*std::pow(unitconv::mtoyd,_unit_powers_[DIST_UNIT_INDEX]);
}

// Time units

const double brgastro::unit_obj::s() const
{
	return _value_;
}

const double brgastro::unit_obj::ms() const
{
	return _value_*std::pow(unitconv::stoms,_unit_powers_[TIME_UNIT_INDEX]);
}

const double brgastro::unit_obj::cs() const
{
	return _value_*std::pow(unitconv::stocs,_unit_powers_[TIME_UNIT_INDEX]);
}

const double brgastro::unit_obj::ns() const
{
	return _value_*std::pow(unitconv::stons,_unit_powers_[TIME_UNIT_INDEX]);
}

const double brgastro::unit_obj::us() const
{
	return _value_*std::pow(unitconv::stous,_unit_powers_[TIME_UNIT_INDEX]);
}

const double brgastro::unit_obj::min() const
{
	return _value_*std::pow(unitconv::stomin,_unit_powers_[TIME_UNIT_INDEX]);
}

const double brgastro::unit_obj::hr() const
{
	return _value_*std::pow(unitconv::stohr,_unit_powers_[TIME_UNIT_INDEX]);
}

const double brgastro::unit_obj::day() const
{
	return _value_*std::pow(unitconv::stoday,_unit_powers_[TIME_UNIT_INDEX]);
}

const double brgastro::unit_obj::week() const
{
	return _value_*std::pow(unitconv::stoweek,_unit_powers_[TIME_UNIT_INDEX]);
}

const double brgastro::unit_obj::month() const
{
	return _value_*std::pow(unitconv::stomonth,_unit_powers_[TIME_UNIT_INDEX]);
}

const double brgastro::unit_obj::yr() const
{
	return _value_*std::pow(unitconv::stoyr,_unit_powers_[TIME_UNIT_INDEX]);
}

const double brgastro::unit_obj::kyr() const
{
	return _value_*std::pow(unitconv::stokyr,_unit_powers_[TIME_UNIT_INDEX]);
}

const double brgastro::unit_obj::Myr() const
{
	return _value_*std::pow(unitconv::stoMyr,_unit_powers_[TIME_UNIT_INDEX]);
}

const double brgastro::unit_obj::Gyr() const
{
	return _value_*std::pow(unitconv::stoGyr,_unit_powers_[TIME_UNIT_INDEX]);
}

const double brgastro::unit_obj::kg() const
{
	return _value_;
}

const double brgastro::unit_obj::gm() const
{
	return _value_*std::pow(unitconv::kgtogm,_unit_powers_[MASS_UNIT_INDEX]);
}

const double brgastro::unit_obj::Mearth() const
{
	return _value_*std::pow(unitconv::kgtoMearth,_unit_powers_[MASS_UNIT_INDEX]);
}

const double brgastro::unit_obj::Msun() const
{
	return _value_*std::pow(unitconv::kgtoMsun,_unit_powers_[MASS_UNIT_INDEX]);
}

const double brgastro::unit_obj::ttMsun() const
{
	return _value_*std::pow(unitconv::kgtottMsun,_unit_powers_[MASS_UNIT_INDEX]);
}

// Temperature

const double brgastro::unit_obj::K() const
{
	return _value_;
}

const double brgastro::unit_obj::degC() const
{
	// Check if this is a pure unit_temperature. If so, do proper conversion. Otherwise just scale
	for(int i=0; i<NUM_UNIT_TYPES; i++)
	{
		if((_unit_powers_[i]!=0)&&(i!=TEMP_UNIT_INDEX))
		{
			return _value_*std::pow(unitconv::KtodegC,_unit_powers_[TEMP_UNIT_INDEX]);
		}
	}
	if(_unit_powers_[TEMP_UNIT_INDEX] != 1)
	{
		return _value_*std::pow(unitconv::KtodegC,_unit_powers_[TEMP_UNIT_INDEX]);
	}
	return _value_*KtodegC+ABS_ZERO_C;
}

const double brgastro::unit_obj::degR() const
{
	return _value_*std::pow(unitconv::KtodegR,_unit_powers_[TEMP_UNIT_INDEX]);
}

const double brgastro::unit_obj::degF() const
{
	// Check if this is a pure unit_temperature. If so, do proper conversion. Otherwise just scale
	for(int i=0; i<NUM_UNIT_TYPES; i++)
	{
		if((_unit_powers_[i]!=0)&&(i!=TEMP_UNIT_INDEX))
		{
			return _value_*std::pow(unitconv::KtodegF,_unit_powers_[TEMP_UNIT_INDEX]);
		}
	}
	if(_unit_powers_[TEMP_UNIT_INDEX] != 1)
	{
		return _value_*std::pow(unitconv::KtodegF,_unit_powers_[TEMP_UNIT_INDEX]);
	}
	return _value_*KtodegF+ABS_ZERO_F;
}

// Velocity

const double brgastro::unit_obj::mps() const
{
	return _value_;
}

const double brgastro::unit_obj::kmps() const
{
	return _value_*mpstokmps;
}

const double brgastro::unit_obj::c() const
{
	return _value_*mpstoc;
}

const double brgastro::unit_obj::miphr() const
{
	return _value_*mpstomiphr;
}

const double brgastro::unit_obj::kpcpGyr() const
{
	return _value_*mtokpc/stoGyr;
}

// Acceleration

const double brgastro::unit_obj::kmpspGyr() const
{
	return _value_*unitconv::mtokm/unitconv::stos/unitconv::stoGyr;
}
const double brgastro::unit_obj::kpcpGyr2() const
{
	return _value_*unitconv::mtokpc/unitconv::stoGyr/unitconv::stoGyr;
}

// Charge

const double brgastro::unit_obj::C() const
{
	return _value_;
}

const double brgastro::unit_obj::esu() const
{
	return _value_*std::pow(unitconv::Ctoesu,_unit_powers_[CHRG_UNIT_INDEX]);
}

// Operator overloading

// Copy constructor
brgastro::unit_obj::unit_obj(const brgastro::unit_obj& old_unit_obj, bool maintain_unit_fix)
{
	make_array(_unit_powers_, NUM_UNIT_TYPES);
	for(int i = 0; i < NUM_UNIT_TYPES; i++) _unit_powers_[i]=old_unit_obj._unit_powers_[i];
	_value_ = old_unit_obj._value_;
	if(maintain_unit_fix)
	{
		_fix_unit_powers_ = old_unit_obj._fix_unit_powers_;
	}
	else
	{
		_fix_unit_powers_ = false;
	}
}

// Assignment
brgastro::unit_obj & brgastro::unit_obj::operator=(const brgastro::unit_obj &old_unit_obj)
{
	if (this != &old_unit_obj)
	{
#ifdef _BRG_WARN_FOR_UNIT_MISMATCH_
		bool no_problem = true;
		// First check if the variable is unitless. If so, allow changing it
		for(int i = 0; i < NUM_UNIT_TYPES; i++)
		{
			if(_unit_powers_[i]!=0) no_problem = false;
		}

		if(!no_problem)
		{
			no_problem = true;
			for(int i = 0; i < NUM_UNIT_TYPES; i++)
			{
				if((_unit_powers_[i]!=old_unit_obj._unit_powers_[i])&&(i!=ANGL_UNIT_INDEX)) // We don't check unit_angle here, as that might make sense
				{
					std::cerr << "WARNING: Unit mismatch in assignment.\n";
					no_problem = false;
					break;
				}
			}
			// Check angle now, make sure we're only allowing proper assignments with it
			if(no_problem && (_unit_powers_[ANGL_UNIT_INDEX] != old_unit_obj._unit_powers_[ANGL_UNIT_INDEX]))
			{
				if(((_unit_powers_[ANGL_UNIT_INDEX] == 0 ) && (old_unit_obj._unit_powers_[ANGL_UNIT_INDEX] == 1)) ||
						((_unit_powers_[ANGL_UNIT_INDEX] == 1 ) && (old_unit_obj._unit_powers_[ANGL_UNIT_INDEX] == 0)))
				{
					_unit_powers_[ANGL_UNIT_INDEX] = 1;
				}
				else
				{
					std::cerr << "WARNING: Unit mismatch in assignment.\n";
					no_problem = false;
				}
			}
		}
#endif // #ifdef _BRG_WARN_FOR_UNIT_MISMATCH_

		for(int i = 0; i < NUM_UNIT_TYPES; i++) _unit_powers_[i]=old_unit_obj._unit_powers_[i];
		_value_ = old_unit_obj._value_;
	}
	return *this;
}
brgastro::unit_obj & brgastro::unit_obj::operator=(const int &new__value_ue)
{
	_value_ = new__value_ue;
	return *this;
}
brgastro::unit_obj & brgastro::unit_obj::operator=(const double &new__value_ue)
{
	_value_ = new__value_ue;
	return *this;
}
brgastro::unit_obj & brgastro::unit_obj::operator=(const float &new__value_ue)
{
	_value_ = new__value_ue;
	return *this;
}

// Addition
brgastro::unit_obj & brgastro::unit_obj::operator+=(const brgastro::unit_obj &other_unit_obj)
{
#ifndef _BRG_WARN_FOR_UNIT_MISMATCH_

	_value_ += other_unit_obj._value_;
	return *this;

#else

	// If current _value_ue is zero, simply overwrite with other _value_ue
	if(_value_==0)
	{
		return (*this = other_unit_obj);
	}

	// If other _value_ is zero, return
	if(other_unit_obj._value_ == 0)
	{
		return *this;
	}

	bool no_problem = true;
	for(int i = 0; i < NUM_UNIT_TYPES; i++)
	{
		if((_unit_powers_[i]!=other_unit_obj._unit_powers_[i])&&(i!=ANGL_UNIT_INDEX)) // We don't check unit_angle here, as that might make sense
		{
			std::cerr << "WARNING: Unit mismatch in addition.\n";
			no_problem = false;
			break;
		}
	}
	// Check angle now, make sure we're only allowing proper assignments with it
	if(no_problem && (_unit_powers_[ANGL_UNIT_INDEX] != other_unit_obj._unit_powers_[ANGL_UNIT_INDEX]))
	{
		if(((_unit_powers_[ANGL_UNIT_INDEX] == 0 ) && (other_unit_obj._unit_powers_[ANGL_UNIT_INDEX] == 1)) ||
				(((_unit_powers_[ANGL_UNIT_INDEX] == 1 ) && (other_unit_obj._unit_powers_[ANGL_UNIT_INDEX] == 0))))
		{
			_unit_powers_[ANGL_UNIT_INDEX] = 1;
		}
		else
		{
			std::cerr << "WARNING: Unit mismatch in addition.\n";
			no_problem = false;
		}
	}

	_value_ += other_unit_obj._value_;
	return *this;

#endif // #ifndef _BRG_WARN_FOR_UNIT_MISMATCH_ -> #else
}
const brgastro::unit_obj brgastro::unit_obj::operator+(const brgastro::unit_obj &other_unit_obj) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result += other_unit_obj;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator+=(const double &new__value_ue)
{
	_value_ += new__value_ue;
	return *this;
}
const brgastro::unit_obj brgastro::unit_obj::operator+(const double &new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result += new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator+=(const float &new__value_ue)
{
	_value_ += new__value_ue;
	return *this;
}
const brgastro::unit_obj brgastro::unit_obj::operator+(const float &new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result += new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator+=(const int &new__value_ue)
{
	_value_ += new__value_ue;
	return *this;
}
const brgastro::unit_obj brgastro::unit_obj::operator+(const int &new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result += new__value_ue;
	return result;
}

// Subtraction
brgastro::unit_obj & brgastro::unit_obj::operator-=(const brgastro::unit_obj &other_unit_obj)
{
#ifndef _BRG_WARN_FOR_UNIT_MISMATCH_

	_value_ -= other_unit_obj._value_;
	return *this;

#else

	// If current _value_ue is zero, simply overwrite with other _value_ue
	if(_value_==0)
	{
		return (*this = -other_unit_obj);
	}

	// If other _value_ is zero, return
	if(other_unit_obj._value_ == 0)
	{
		return *this;
	}

	bool no_problem = true;
	for(int i = 0; i < NUM_UNIT_TYPES; i++)
	{
		if((_unit_powers_[i]!=other_unit_obj._unit_powers_[i])&&(i!=ANGL_UNIT_INDEX)) // We don't check unit_angle here, as that might make sense
		{
			std::cerr << "WARNING: Unit mismatch in subtraction.\n";
			no_problem = false;
			break;
		}
	}
	// Check angle now, make sure we're only allowing proper assignments with it
	if(no_problem && (_unit_powers_[ANGL_UNIT_INDEX] != other_unit_obj._unit_powers_[ANGL_UNIT_INDEX]))
	{
		if(((_unit_powers_[ANGL_UNIT_INDEX] == 0 ) && (other_unit_obj._unit_powers_[ANGL_UNIT_INDEX] == 1)) ||
				(((_unit_powers_[ANGL_UNIT_INDEX] == 1 ) && (other_unit_obj._unit_powers_[ANGL_UNIT_INDEX] == 0))))
		{
			_unit_powers_[ANGL_UNIT_INDEX] = 1;
		}
		else
		{
			std::cerr << "WARNING: Unit mismatch in subtraction.\n";
			no_problem = false;
		}
	}

	_value_ -= other_unit_obj._value_;
	return *this;

#endif // #ifndef _BRG_WARN_FOR_UNIT_MISMATCH_ -> #else
}
const brgastro::unit_obj brgastro::unit_obj::operator-(const brgastro::unit_obj &other_unit_obj) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result -= other_unit_obj;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator-=(const double &new__value_ue)
{
	_value_ -= new__value_ue;
	return *this;
}
const brgastro::unit_obj brgastro::unit_obj::operator-(const double &new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result -= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator-=(const float &new__value_ue)
{
	_value_ -= new__value_ue;
	return *this;
}
const brgastro::unit_obj brgastro::unit_obj::operator-(const float &new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result -= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator-=(const int &new__value_ue)
{
	_value_ -= new__value_ue;
	return *this;
}
const brgastro::unit_obj brgastro::unit_obj::operator-(const int &new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result -= new__value_ue;
	return result;
}

// Multiplication
brgastro::unit_obj & brgastro::unit_obj::operator*=(const brgastro::unit_obj &other_unit_obj)
{
	for(int i = 0; i < NUM_UNIT_TYPES; i++)
	{
		_unit_powers_[i] += other_unit_obj._unit_powers_[i];
	} // for(int i = 0; i < NUM_UNIT_TYPES; i++)
	_value_ *= other_unit_obj._value_;
	return *this;
}
const brgastro::unit_obj brgastro::unit_obj::operator*(const brgastro::unit_obj &other_unit_obj) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result *= other_unit_obj;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator*=(const double &new__value_ue)
{
	_value_ *= new__value_ue;
	return *this;
}
const brgastro::unit_obj brgastro::unit_obj::operator*(const double &new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result *= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator*=(const float &new__value_ue)
{
	_value_ *= new__value_ue;
	return *this;
}
const brgastro::unit_obj brgastro::unit_obj::operator*(const float &new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result *= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator*=(const int &new__value_ue)
{
	_value_ *= new__value_ue;
	return *this;
}
const brgastro::unit_obj brgastro::unit_obj::operator*(const int &new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result *= new__value_ue;
	return result;
}

// Division
brgastro::unit_obj & brgastro::unit_obj::operator/=(const brgastro::unit_obj &other_unit_obj)
{
	for(int i = 0; i < NUM_UNIT_TYPES; i++)
	{
		_unit_powers_[i] -= other_unit_obj._unit_powers_[i];
	} // for(int i = 0; i < NUM_UNIT_TYPES; i++)
	_value_ /= other_unit_obj._value_;
	return *this;
}
const brgastro::unit_obj brgastro::unit_obj::operator/(const brgastro::unit_obj &other_unit_obj) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result /= other_unit_obj;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator/=(const double &new__value_ue)
{
	_value_ /= new__value_ue;
	return *this;
}
const brgastro::unit_obj brgastro::unit_obj::operator/(const double &new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result /= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator/=(const float &new__value_ue)
{
	_value_ /= new__value_ue;
	return *this;
}
const brgastro::unit_obj brgastro::unit_obj::operator/(const float &new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result /= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator/=(const int &new__value_ue)
{
	_value_ /= new__value_ue;
	return *this;
}
const brgastro::unit_obj brgastro::unit_obj::operator/(const int &new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result /= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator++()   //Prefix
{
	_value_ += 1;
	return *this;
}
const brgastro::unit_obj brgastro::unit_obj::operator++(int)   //Postfix
{
	brgastro::unit_obj obj_copy = brgastro::unit_obj(*this);
	_value_ += 1;
	return obj_copy;
}
brgastro::unit_obj & brgastro::unit_obj::operator--()   //Prefix
{
	_value_ -= 1;
	return *this;
}
const brgastro::unit_obj brgastro::unit_obj::operator--(int)   //Postfix
{
	brgastro::unit_obj obj_copy = brgastro::unit_obj(*this);
	_value_ -= 1;
	return obj_copy;
}
const brgastro::unit_obj brgastro::unit_obj::operator-() const   //Prefix
{
	brgastro::unit_obj result = brgastro::unit_obj(*this);
	return result *= -1;
}
const bool brgastro::unit_obj::operator<(const brgastro::unit_obj & other_unit_obj) const
{
	return (_value_ < other_unit_obj._value_);
}
const bool brgastro::unit_obj::operator<(const int &comp__value_) const
{
	return (_value_ < comp__value_);
}
const bool brgastro::unit_obj::operator<(const double &comp__value_) const
{
	return (_value_ < comp__value_);
}
const bool brgastro::unit_obj::operator<(const float &comp__value_) const
{
	return (_value_ < comp__value_);
}
const bool brgastro::unit_obj::operator>(const brgastro::unit_obj & other_unit_obj) const
{
	return (_value_ > other_unit_obj._value_);
}
const bool brgastro::unit_obj::operator>(const int &comp__value_) const
{
	return (_value_ > comp__value_);
}
const bool brgastro::unit_obj::operator>(const double &comp__value_) const
{
	return (_value_ > comp__value_);
}
const bool brgastro::unit_obj::operator>(const float &comp__value_) const
{
	return (_value_ > comp__value_);
}
const bool brgastro::unit_obj::operator==(const brgastro::unit_obj & other_unit_obj) const
{
	return (_value_ == other_unit_obj._value_);
}
const bool brgastro::unit_obj::operator==(const int &comp__value_) const
{
	return (_value_ == comp__value_);
}
const bool brgastro::unit_obj::operator==(const double &comp__value_) const
{
	return (_value_ == comp__value_);
}
const bool brgastro::unit_obj::operator==(const float &comp__value_) const
{
	return (_value_ == comp__value_);
}
const bool brgastro::unit_obj::operator<=(const brgastro::unit_obj & other_unit_obj) const
{
	return (_value_ <= other_unit_obj._value_);
}
const bool brgastro::unit_obj::operator<=(const int &comp__value_) const
{
	return (_value_ <= comp__value_);
}
const bool brgastro::unit_obj::operator<=(const double &comp__value_) const
{
	return (_value_ <= comp__value_);
}
const bool brgastro::unit_obj::operator<=(const float &comp__value_) const
{
	return (_value_ <= comp__value_);
}
const bool brgastro::unit_obj::operator>=(const brgastro::unit_obj & other_unit_obj) const
{
	return (_value_ >= other_unit_obj._value_);
}
const bool brgastro::unit_obj::operator>=(const int &comp__value_) const
{
	return (_value_ >= comp__value_);
}
const bool brgastro::unit_obj::operator>=(const double &comp__value_) const
{
	return (_value_ >= comp__value_);
}
const bool brgastro::unit_obj::operator>=(const float &comp__value_) const
{
	return (_value_ >= comp__value_);
}
const bool brgastro::unit_obj::operator!=(const brgastro::unit_obj & other_unit_obj) const
{
	return (_value_ != other_unit_obj._value_);
}
const bool brgastro::unit_obj::operator!=(const int &comp__value_) const
{
	return (_value_ != comp__value_);
}
const bool brgastro::unit_obj::operator!=(const double &comp__value_) const
{
	return (_value_ != comp__value_);
}
const bool brgastro::unit_obj::operator!=(const float &comp__value_) const
{
	return (_value_ != comp__value_);
}

brgastro::unit_obj::operator double() const
{
	return _value_;
}

// Overloaded operators relating to unit_objs

const brgastro::unit_obj operator+(const int lhs, const brgastro::unit_obj & rhs)
{
	return rhs+lhs;
}
const brgastro::unit_obj operator+(const double lhs, const brgastro::unit_obj & rhs)
{
	return rhs+lhs;
}
const brgastro::unit_obj operator+(const float lhs, const brgastro::unit_obj & rhs)
{
	return rhs+lhs;
}
const brgastro::unit_obj operator-(const int lhs, const brgastro::unit_obj & rhs)
{
	return -(rhs-lhs);
}
const brgastro::unit_obj operator-(const double lhs, const brgastro::unit_obj & rhs)
{
	return -(rhs-lhs);
}
const brgastro::unit_obj operator-(const float lhs, const brgastro::unit_obj & rhs)
{
	return -(rhs-lhs);
}
const brgastro::unit_obj operator*(const int lhs, const brgastro::unit_obj & rhs)
{
	return rhs*lhs;
}
const brgastro::unit_obj operator*(const double lhs, const brgastro::unit_obj & rhs)
{
	return rhs*lhs;
}
const brgastro::unit_obj operator*(const float lhs, const brgastro::unit_obj & rhs)
{
	return rhs*lhs;
}
const brgastro::unit_obj operator/(const int lhs, const brgastro::unit_obj & rhs)
{
	return brgastro::unit_obj(lhs)/rhs;
}
const brgastro::unit_obj operator/(const double lhs, const brgastro::unit_obj & rhs)
{
	return brgastro::unit_obj(lhs)/rhs;
}
const brgastro::unit_obj operator/(const float lhs, const brgastro::unit_obj & rhs)
{
	return brgastro::unit_obj(lhs)/rhs;
}

const bool operator<(const double lhs, const brgastro::unit_obj & rhs)
{
	return rhs > lhs;
}
const bool operator<(const float lhs, const brgastro::unit_obj & rhs)
{
	return rhs > lhs;
}
const bool operator<(const int lhs, const brgastro::unit_obj & rhs)
{
	return rhs > lhs;
}
const bool operator>(const double lhs, const brgastro::unit_obj & rhs)
{
	return rhs < lhs;
}
const bool operator>(const float lhs, const brgastro::unit_obj & rhs)
{
	return rhs < lhs;
}
const bool operator>(const int lhs, const brgastro::unit_obj & rhs)
{
	return rhs < lhs;
}
const bool operator<=(const double lhs, const brgastro::unit_obj & rhs)
{
	return rhs >= lhs;
}
const bool operator<=(const float lhs, const brgastro::unit_obj & rhs)
{
	return rhs >= lhs;
}
const bool operator<=(const int lhs, const brgastro::unit_obj & rhs)
{
	return rhs >= lhs;
}
const bool operator>=(const double lhs, const brgastro::unit_obj & rhs)
{
	return rhs <= lhs;
}
const bool operator>=(const float lhs, const brgastro::unit_obj & rhs)
{
	return rhs <= lhs;
}
const bool operator>=(const int lhs, const brgastro::unit_obj & rhs)
{
	return rhs <= lhs;
}
const bool operator==(const double lhs, const brgastro::unit_obj & rhs)
{
	return rhs == lhs;
}
const bool operator==(const float lhs, const brgastro::unit_obj & rhs)
{
	return rhs == lhs;
}
const bool operator==(const int lhs, const brgastro::unit_obj & rhs)
{
	return rhs == lhs;
}
const bool operator!=(const double lhs, const brgastro::unit_obj & rhs)
{
	return rhs != lhs;
}
const bool operator!=(const float lhs, const brgastro::unit_obj & rhs)
{
	return rhs != lhs;
}
const bool operator!=(const int lhs, const brgastro::unit_obj & rhs)
{
	return rhs != lhs;
}

std::ostream & operator<< (std::ostream &out, brgastro::unit_obj &obj)
{
	out << obj.get_value();
	return out;
}

const brgastro::unit_obj std::pow(const brgastro::unit_obj & lhs, const double rhs)
{
	brgastro::unit_obj result = brgastro::unit_obj(lhs);
	result.set_value(std::pow(lhs.get_value(),rhs));
	for(int i = 0; i < NUM_UNIT_TYPES; i++) result._unit_powers_[i]*=(float)rhs;
	result.round_powers();
	return result;
}
const brgastro::unit_obj std::pow(const brgastro::unit_obj & lhs, const float rhs)
{
	brgastro::unit_obj result = brgastro::unit_obj(lhs);
	result.set_value(std::pow(lhs.get_value(),rhs));
	for(int i = 0; i < NUM_UNIT_TYPES; i++) result._unit_powers_[i]*=rhs;
	result.round_powers();
	return result;
}
const brgastro::unit_obj std::pow(const brgastro::unit_obj & lhs, const int rhs)
{
	brgastro::unit_obj result = brgastro::unit_obj(lhs);
	result.set_value(std::pow(lhs.get_value(),rhs));
	for(int i = 0; i < NUM_UNIT_TYPES; i++) result._unit_powers_[i]*=rhs;
	result.round_powers();
	return result;
}
const brgastro::unit_obj std::sqrt(const brgastro::unit_obj & obj)
{
	brgastro::unit_obj result = brgastro::unit_obj(obj);
	result.set_value(std::sqrt(obj.get_value()));
	for(int i = 0; i < NUM_UNIT_TYPES; i++) result._unit_powers_[i]/=2.;
	result.round_powers();
	return result;
}
const brgastro::unit_obj std::fabs(const brgastro::unit_obj & obj)
{
	brgastro::unit_obj result = brgastro::unit_obj(obj);
	if(result<0) result *= -1;
	return result;
}

//unit_distance full functions
// Default units: meters (m)
brgastro::unit_distance::unit_distance()
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[DIST_UNIT_INDEX] = 1;
	_value_ = 0;
	_fix_unit_powers_ = true;
}

brgastro::unit_distance::unit_distance(double init__value_, double conv_factor)
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[DIST_UNIT_INDEX] = 1;
	_value_ = init__value_*conv_factor;
	_fix_unit_powers_ = true;
}

brgastro::unit_distance::unit_distance(const brgastro::unit_obj & other_unit_obj)
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[DIST_UNIT_INDEX] = 1;
	_value_ = other_unit_obj.get_value();
	_fix_unit_powers_ = true;
}

//unit_time full functions
brgastro::unit_time::unit_time()
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[TIME_UNIT_INDEX] = 1;
	_value_ = 0;
	_fix_unit_powers_ = true;
}

brgastro::unit_time::unit_time(double init__value_, double conv_factor)
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[TIME_UNIT_INDEX] = 1;
	_value_ = init__value_*conv_factor;
	_fix_unit_powers_ = true;
}

brgastro::unit_time::unit_time(const brgastro::unit_obj & other_unit_obj)
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[TIME_UNIT_INDEX] = 1;
	_value_ = other_unit_obj.get_value();
	_fix_unit_powers_ = true;
}

//Velocity full functions
// Default units: meters per second (m/s)
brgastro::unit_velocity::unit_velocity()
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[DIST_UNIT_INDEX] = 1;
	_unit_powers_[TIME_UNIT_INDEX] = -1;
	_value_ = 0;
	_fix_unit_powers_ = true;
}

brgastro::unit_velocity::unit_velocity(double init__value_, double conv_factor)
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[DIST_UNIT_INDEX] = 1;
	_unit_powers_[TIME_UNIT_INDEX] = -1;
	_value_ = init__value_*conv_factor;
	_fix_unit_powers_ = true;
}

brgastro::unit_velocity::unit_velocity(const brgastro::unit_obj & other_unit_obj)
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[DIST_UNIT_INDEX] = 1;
	_unit_powers_[TIME_UNIT_INDEX] = -1;
	_value_ = other_unit_obj.get_value();
	_fix_unit_powers_ = true;
}

//unit_mass full functions
// Default units: kilograms (kg)

brgastro::unit_mass::unit_mass()
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[MASS_UNIT_INDEX] = 1;
	_value_ = 0;
	_fix_unit_powers_ = true;
}

brgastro::unit_mass::unit_mass(double init__value_, double conv_factor)
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[MASS_UNIT_INDEX] = 1;
	_value_ = init__value_*conv_factor;
	_fix_unit_powers_ = true;
}

brgastro::unit_mass::unit_mass(const brgastro::unit_obj & other_unit_obj)
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[MASS_UNIT_INDEX] = 1;
	_value_ = other_unit_obj.get_value();
	_fix_unit_powers_ = true;
}

//unit_temperature full functions
// Default units: Kelvins

brgastro::unit_temperature::unit_temperature()
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[TEMP_UNIT_INDEX] = 1;
	_value_ = 0;
	_fix_unit_powers_ = true;
}

brgastro::unit_temperature::unit_temperature(double init__value_, double conv_factor)
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[TEMP_UNIT_INDEX] = 1;
	_value_ = init__value_*conv_factor;
	_fix_unit_powers_ = true;
}

brgastro::unit_temperature::unit_temperature(const brgastro::unit_obj & other_unit_obj)
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[TEMP_UNIT_INDEX] = 1;
	_value_ = other_unit_obj.get_value();
	_fix_unit_powers_ = true;
}

const double brgastro::unit_temperature::K() const
{
	return _value_;
}

const double brgastro::unit_temperature::C() const
{
	return _value_+ABS_ZERO_C;
}

const double brgastro::unit_temperature::R() const
{
	return _value_*KtodegR;
}

const double brgastro::unit_temperature::F() const
{
	return _value_*KtodegR+ABS_ZERO_F;
}

//unit_angle full functions
brgastro::unit_angle::unit_angle()
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[ANGL_UNIT_INDEX] = 1;
	_value_ = 0;
	_fix_unit_powers_ = true;
}
brgastro::unit_angle::unit_angle(double init__value_, double conv_factor)
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[ANGL_UNIT_INDEX] = 1;
	_value_ = init__value_*conv_factor;
	_fix_unit_powers_ = true;
}
brgastro::unit_angle::unit_angle(const brgastro::unit_obj & other_unit_obj)
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[ANGL_UNIT_INDEX] = 1;
	_value_ = other_unit_obj.get_value();
	_fix_unit_powers_ = true;
}

//unit_charge full functions
brgastro::unit_charge::unit_charge()
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[ANGL_UNIT_INDEX] = 1;
	_value_ = 0;
	_fix_unit_powers_ = true;
}
brgastro::unit_charge::unit_charge(double init__value_, double conv_factor)
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[ANGL_UNIT_INDEX] = 1;
	_value_ = init__value_*conv_factor;
	_fix_unit_powers_ = true;
}
brgastro::unit_charge::unit_charge(const brgastro::unit_obj & other_unit_obj)
{
	_unit_powers_.resize(NUM_UNIT_TYPES,0);
	_unit_powers_[ANGL_UNIT_INDEX] = 1;
	_value_ = other_unit_obj.get_value();
	_fix_unit_powers_ = true;
}

#endif // #ifdef _BRG_USE_UNITS_
