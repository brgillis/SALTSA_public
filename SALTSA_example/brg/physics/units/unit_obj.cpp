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
#include <stdexcept>
#include <memory>

#include "brg/math/misc_math.hpp"
#include "brg/physics/units/unit_conversions.hpp"

#include "unit_obj.h"

using namespace brgastro::unitconv;
using std::cout;
using std::cerr;
using std::string;
using std::stringstream;

//-------------------------
// Full functions (brgastro::unit_obj)
//-------------------------

brgastro::unit_obj::unit_obj(double init__value_, float d_units_power, float t_units_power, float m_units_power,
		float T_units_power, float a_units_power, float c_units_power)
{
	_fix_unit_powers_ = false;
	_unit_powers_.resize(NUM_UNIT_TYPES,0);

	set(init__value_,d_units_power,t_units_power,m_units_power,
			T_units_power,a_units_power,c_units_power);
}
brgastro::unit_obj::unit_obj(double init__value_,
		double d_units, float d_units_power, double t_units, float t_units_power,
		double m_units, float m_units_power, double T_units, float T_units_power,
		double a_units, float a_units_power, double c_units, float c_units_power)
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

void brgastro::unit_obj::round_powers()
{
	if(_fix_unit_powers_) return;
	for(size_t i=0; i < NUM_UNIT_TYPES; i++)
	{
		_unit_powers_[i] = round_int(_unit_powers_[i]*UNIT_POWER_ROUND_PRECISION_FACTOR)/(float)UNIT_POWER_ROUND_PRECISION_FACTOR;
	}
}

// Reset variable and set units (distance, time, mass, temperature, angle, all taken as doubles)
void brgastro::unit_obj::reset(float dp, float tp, float mp, float Tp, float ap, float cp)
{
	_value_ = 0;
	if(_fix_unit_powers_) return;
	_unit_powers_[DIST_UNIT_INDEX] = dp;
	_unit_powers_[TIME_UNIT_INDEX] = tp;
	_unit_powers_[MASS_UNIT_INDEX] = mp;
	_unit_powers_[TEMP_UNIT_INDEX] = Tp;
	_unit_powers_[ANGL_UNIT_INDEX] = ap;
	_unit_powers_[CHRG_UNIT_INDEX] = cp;
}
// Set units, but don't change _value_ue in default units - use with care!
void brgastro::unit_obj::set_unit_powers(float dp, float tp, float mp, float Tp, float ap, float cp)
{
	if(!_fix_unit_powers_)
	{
		_unit_powers_[DIST_UNIT_INDEX] = dp;
		_unit_powers_[TIME_UNIT_INDEX] = tp;
		_unit_powers_[MASS_UNIT_INDEX] = mp;
		_unit_powers_[TEMP_UNIT_INDEX] = Tp;
		_unit_powers_[ANGL_UNIT_INDEX] = ap;
		_unit_powers_[CHRG_UNIT_INDEX] = cp;
	}
	else
	{
#ifdef _BRG_WARN_FOR_UNIT_MISMATCH_
		std::cerr << "WARNING: Attempt to change unit powers of brgastro::unit_obj with fixed unit powers.\n";
#endif
	}
}
// Set units, but don't change _value_ue in default units - use with care!
void brgastro::unit_obj::set_unit_powers(std::vector<float> new__unit_powers_)
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
			return;
		}
	}
	else
	{
#ifdef _BRG_WARN_FOR_UNIT_MISMATCH_
		std::cerr << "WARNING: Attempt to change unit powers of brgastro::unit_obj with fixed unit powers.\n";
#endif
	}
}
// Set _value_ue in default or specified units
void brgastro::unit_obj::set_value(double init__value_, double conv_factor)
{
	_value_ = init__value_*conv_factor;
}
// Set _value_ue in specified units, but don't change units
// IMPORTANT: Always use the unitconv::pctom (units used here to default) form for units)
// EXAMPLE 1: To load 4 kpc, use set_value_ue(4, kpctom, 0, 0, 0, 0, 0)
// EXAMPLE 2: To load 1 m/s^2, use set_value_ue(1, mtom, stos, 0, 0, 0, 0) (the function will automatically figure out powers from the unit)
void brgastro::unit_obj::set_value(double new__value_, double d_units, double t_units, double m_units, double T_units,
		double a_units, double c_units)
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
}

// Set _value_ue and changes units - used by overloaded operators mostly
void brgastro::unit_obj::set(double new__value_, const std::vector<float> & new__unit_powers_)
{
	if(new__unit_powers_.size() != NUM_UNIT_TYPES)
		throw std::logic_error("Invalide unit powers vector size.");
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
			std::cerr << "WARNING: Attempt to change unit powers of brgastro::unit_obj with fixed unit powers.\n";
		}
	}
	return;
}
void brgastro::unit_obj::set(double new__value_, float dp, float tp, float mp, float Tp, float ap, float cp)
{
	return set(new__value_, 1, dp, 1, tp, 1, mp, 1, Tp, 1, ap, 1, cp);
}

// Set _value_ue in specified units, and change units of variable
// IMPORTANT: Always use the unitconv::pctom (units used here to default) form for units)
// EXAMPLE 1: To load 4 kpc, use set_value_ue(4, kpctom, 1, 0, 0, 0, 0, 0, 0, 0, 0)
// EXAMPLE 2: To load 1 m/s^2, use set_value_ue(1, mtom, 1, stos, -2, 0, 0, 0, 0, 0, 0)
void brgastro::unit_obj::set(double new__value_, double d_units, float d_units_power, double t_units,
		float t_units_power, double m_units, float m_units_power, double T_units,
		float T_units_power, double a_units, float a_units_power, double c_units,
		float c_units_power)
{
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
}

// Get _value_ue in default units
double brgastro::unit_obj::get_value() const
{
	return _value_;
}

double brgastro::unit_obj::get_value(double d_units, double t_units, double m_units, double T_units, double a_units,
		double c_units) const
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

double brgastro::unit_obj::get_value(double d_units, float d_units_power, double t_units, float t_units_power,
		double m_units, float m_units_power, double T_units, float T_units_power,
		double a_units, float a_units_power, double c_units, float c_units_power) const
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

std::vector<float> brgastro::unit_obj::get_unit_powers() const
{
	// Returns a 6-component std::vector<float> in the format [d_units_power,t_units_power,m_units_power,T_units_power,a_units_power,c_units_power]
	return _unit_powers_;
}

std::string brgastro::unit_obj::get_string() const
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
double brgastro::unit_obj::m() const
{
	return _value_;
}

double brgastro::unit_obj::mm() const
{
	return _value_*std::pow(unitconv::mtomm,_unit_powers_[DIST_UNIT_INDEX]);
}

double brgastro::unit_obj::um() const
{
	return _value_*std::pow(unitconv::mtoum,_unit_powers_[DIST_UNIT_INDEX]);
}

double brgastro::unit_obj::nm() const
{
	return _value_*std::pow(unitconv::mtonm,_unit_powers_[DIST_UNIT_INDEX]);
}

double brgastro::unit_obj::cm() const
{
	return _value_*std::pow(unitconv::mtocm,_unit_powers_[DIST_UNIT_INDEX]);
}

double brgastro::unit_obj::angstrom() const
{
	return _value_*std::pow(unitconv::mtoangstrom,_unit_powers_[DIST_UNIT_INDEX]);
}

double brgastro::unit_obj::km() const
{
	return _value_*std::pow(unitconv::mtokm,_unit_powers_[DIST_UNIT_INDEX]);
}

double brgastro::unit_obj::ltyr() const
{
	return _value_*std::pow(unitconv::mtoltyr,_unit_powers_[DIST_UNIT_INDEX]);
}

double brgastro::unit_obj::pc() const
{
	return _value_*std::pow(unitconv::mtopc,_unit_powers_[DIST_UNIT_INDEX]);
}

double brgastro::unit_obj::kpc() const
{
	return _value_*std::pow(unitconv::mtokpc,_unit_powers_[DIST_UNIT_INDEX]);
}

double brgastro::unit_obj::Mpc() const
{
	return _value_*std::pow(unitconv::mtoMpc,_unit_powers_[DIST_UNIT_INDEX]);
}

double brgastro::unit_obj::mi() const
{
	return _value_*std::pow(unitconv::mtomi,_unit_powers_[DIST_UNIT_INDEX]);
}

double brgastro::unit_obj::Mmi() const
{
	return _value_*std::pow(unitconv::mtoMmi,_unit_powers_[DIST_UNIT_INDEX]);
}

double brgastro::unit_obj::ft() const
{
	return _value_*std::pow(unitconv::mtoft,_unit_powers_[DIST_UNIT_INDEX]);
}

double brgastro::unit_obj::yd() const
{
	return _value_*std::pow(unitconv::mtoyd,_unit_powers_[DIST_UNIT_INDEX]);
}

// Time units

double brgastro::unit_obj::s() const
{
	return _value_;
}

double brgastro::unit_obj::ms() const
{
	return _value_*std::pow(unitconv::stoms,_unit_powers_[TIME_UNIT_INDEX]);
}

double brgastro::unit_obj::cs() const
{
	return _value_*std::pow(unitconv::stocs,_unit_powers_[TIME_UNIT_INDEX]);
}

double brgastro::unit_obj::ns() const
{
	return _value_*std::pow(unitconv::stons,_unit_powers_[TIME_UNIT_INDEX]);
}

double brgastro::unit_obj::us() const
{
	return _value_*std::pow(unitconv::stous,_unit_powers_[TIME_UNIT_INDEX]);
}

double brgastro::unit_obj::min() const
{
	return _value_*std::pow(unitconv::stomin,_unit_powers_[TIME_UNIT_INDEX]);
}

double brgastro::unit_obj::hr() const
{
	return _value_*std::pow(unitconv::stohr,_unit_powers_[TIME_UNIT_INDEX]);
}

double brgastro::unit_obj::day() const
{
	return _value_*std::pow(unitconv::stoday,_unit_powers_[TIME_UNIT_INDEX]);
}

double brgastro::unit_obj::week() const
{
	return _value_*std::pow(unitconv::stoweek,_unit_powers_[TIME_UNIT_INDEX]);
}

double brgastro::unit_obj::month() const
{
	return _value_*std::pow(unitconv::stomonth,_unit_powers_[TIME_UNIT_INDEX]);
}

double brgastro::unit_obj::yr() const
{
	return _value_*std::pow(unitconv::stoyr,_unit_powers_[TIME_UNIT_INDEX]);
}

double brgastro::unit_obj::kyr() const
{
	return _value_*std::pow(unitconv::stokyr,_unit_powers_[TIME_UNIT_INDEX]);
}

double brgastro::unit_obj::Myr() const
{
	return _value_*std::pow(unitconv::stoMyr,_unit_powers_[TIME_UNIT_INDEX]);
}

double brgastro::unit_obj::Gyr() const
{
	return _value_*std::pow(unitconv::stoGyr,_unit_powers_[TIME_UNIT_INDEX]);
}

double brgastro::unit_obj::kg() const
{
	return _value_;
}

double brgastro::unit_obj::gm() const
{
	return _value_*std::pow(unitconv::kgtogm,_unit_powers_[MASS_UNIT_INDEX]);
}

double brgastro::unit_obj::Mearth() const
{
	return _value_*std::pow(unitconv::kgtoMearth,_unit_powers_[MASS_UNIT_INDEX]);
}

double brgastro::unit_obj::Msun() const
{
	return _value_*std::pow(unitconv::kgtoMsun,_unit_powers_[MASS_UNIT_INDEX]);
}

double brgastro::unit_obj::ttMsun() const
{
	return _value_*std::pow(unitconv::kgtottMsun,_unit_powers_[MASS_UNIT_INDEX]);
}

// Temperature

double brgastro::unit_obj::K() const
{
	return _value_;
}

double brgastro::unit_obj::degC() const
{
	// Check if this is a pure unit_temperature. If so, do proper conversion. Otherwise just scale
	for(size_t i=0; i<NUM_UNIT_TYPES; i++)
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

double brgastro::unit_obj::degR() const
{
	return _value_*std::pow(unitconv::KtodegR,_unit_powers_[TEMP_UNIT_INDEX]);
}

double brgastro::unit_obj::degF() const
{
	// Check if this is a pure unit_temperature. If so, do proper conversion. Otherwise just scale
	for(size_t i=0; i<NUM_UNIT_TYPES; i++)
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

double brgastro::unit_obj::mps() const
{
	return _value_;
}

double brgastro::unit_obj::kmps() const
{
	return _value_*mpstokmps;
}

double brgastro::unit_obj::c() const
{
	return _value_*mpstoc;
}

double brgastro::unit_obj::miphr() const
{
	return _value_*mpstomiphr;
}

double brgastro::unit_obj::kpcpGyr() const
{
	return _value_*mtokpc/stoGyr;
}

// Acceleration

double brgastro::unit_obj::kmpspGyr() const
{
	return _value_*unitconv::mtokm/unitconv::stos/unitconv::stoGyr;
}
double brgastro::unit_obj::kpcpGyr2() const
{
	return _value_*unitconv::mtokpc/unitconv::stoGyr/unitconv::stoGyr;
}

// Charge

double brgastro::unit_obj::C() const
{
	return _value_;
}

double brgastro::unit_obj::esu() const
{
	return _value_*std::pow(unitconv::Ctoesu,_unit_powers_[CHRG_UNIT_INDEX]);
}

// Operator overloading

// Copy constructor
brgastro::unit_obj::unit_obj(const brgastro::unit_obj& old_unit_obj, bool maintain_unit_fix)
{
	_unit_powers_ = old_unit_obj._unit_powers_;
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
		for(size_t i = 0; i < NUM_UNIT_TYPES; i++)
		{
			if(_unit_powers_[i]!=0) no_problem = false;
		}

		if(!no_problem)
		{
			no_problem = true;
			for(size_t i = 0; i < NUM_UNIT_TYPES; i++)
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

		for(size_t i = 0; i < NUM_UNIT_TYPES; i++) _unit_powers_[i]=old_unit_obj._unit_powers_[i];
		_value_ = old_unit_obj._value_;
	}
	return *this;
}
brgastro::unit_obj & brgastro::unit_obj::operator=(int new__value_ue)
{
	_value_ = new__value_ue;
	return *this;
}
brgastro::unit_obj & brgastro::unit_obj::operator=(long int new__value_ue)
{
	_value_ = new__value_ue;
	return *this;
}
brgastro::unit_obj & brgastro::unit_obj::operator=(unsigned int new__value_ue)
{
	_value_ = new__value_ue;
	return *this;
}
brgastro::unit_obj & brgastro::unit_obj::operator=(short unsigned int new__value_ue)
{
	_value_ = new__value_ue;
	return *this;
}
brgastro::unit_obj & brgastro::unit_obj::operator=(long unsigned int new__value_ue)
{
	_value_ = new__value_ue;
	return *this;
}
brgastro::unit_obj & brgastro::unit_obj::operator=(double new__value_ue)
{
	_value_ = new__value_ue;
	return *this;
}
brgastro::unit_obj & brgastro::unit_obj::operator=(long double new__value_ue)
{
	_value_ = new__value_ue;
	return *this;
}
brgastro::unit_obj & brgastro::unit_obj::operator=(float new__value_ue)
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
	for(size_t i = 0; i < NUM_UNIT_TYPES; i++)
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
brgastro::unit_obj brgastro::unit_obj::operator+(const brgastro::unit_obj &other_unit_obj) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result += other_unit_obj;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator+=(double new__value_ue)
{
	_value_ += new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator+(double new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result += new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator+=(long double new__value_ue)
{
	_value_ += new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator+(long double new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result += new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator+=(float new__value_ue)
{
	_value_ += new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator+(float new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result += new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator+=(int new__value_ue)
{
	_value_ += new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator+(int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result += new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator+=(long int new__value_ue)
{
	_value_ += new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator+(long int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result += new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator+=(unsigned int new__value_ue)
{
	_value_ += new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator+(unsigned int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result += new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator+=(short unsigned int new__value_ue)
{
	_value_ += new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator+(short unsigned int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result += new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator+=(long unsigned int new__value_ue)
{
	_value_ += new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator+(long unsigned int new__value_ue) const
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
	for(size_t i = 0; i < NUM_UNIT_TYPES; i++)
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
brgastro::unit_obj brgastro::unit_obj::operator-(const brgastro::unit_obj &other_unit_obj) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result -= other_unit_obj;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator-=(double new__value_ue)
{
	_value_ -= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator-(double new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result -= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator-=(long double new__value_ue)
{
	_value_ -= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator-(long double new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result -= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator-=(float new__value_ue)
{
	_value_ -= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator-(float new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result -= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator-=(int new__value_ue)
{
	_value_ -= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator-(int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result -= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator-=(long int new__value_ue)
{
	_value_ -= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator-(long int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result -= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator-=(unsigned int new__value_ue)
{
	_value_ -= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator-(unsigned int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result -= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator-=(short unsigned int new__value_ue)
{
	_value_ -= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator-(short unsigned int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result -= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator-=(long unsigned int new__value_ue)
{
	_value_ -= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator-(long unsigned int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result -= new__value_ue;
	return result;
}

// Multiplication
brgastro::unit_obj & brgastro::unit_obj::operator*=(const brgastro::unit_obj &other_unit_obj)
{
	for(size_t i = 0; i < NUM_UNIT_TYPES; i++)
	{
		_unit_powers_[i] += other_unit_obj._unit_powers_[i];
	} // for(size_t i = 0; i < NUM_UNIT_TYPES; i++)
	_value_ *= other_unit_obj._value_;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator*(const brgastro::unit_obj &other_unit_obj) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result *= other_unit_obj;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator*=(double new__value_ue)
{
	_value_ *= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator*(double new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result *= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator*=(long double new__value_ue)
{
	_value_ *= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator*(long double new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result *= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator*=(float new__value_ue)
{
	_value_ *= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator*(float new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result *= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator*=(int new__value_ue)
{
	_value_ *= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator*(int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result *= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator*=(long int new__value_ue)
{
	_value_ *= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator*(long int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result *= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator*=(unsigned int new__value_ue)
{
	_value_ *= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator*(unsigned int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result *= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator*=(short unsigned int new__value_ue)
{
	_value_ *= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator*(short unsigned int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result *= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator*=(long unsigned int new__value_ue)
{
	_value_ *= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator*(long unsigned int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result *= new__value_ue;
	return result;
}

// Division
brgastro::unit_obj & brgastro::unit_obj::operator/=(const brgastro::unit_obj &other_unit_obj)
{
	for(size_t i = 0; i < NUM_UNIT_TYPES; i++)
	{
		_unit_powers_[i] -= other_unit_obj._unit_powers_[i];
	} // for(size_t i = 0; i < NUM_UNIT_TYPES; i++)
	_value_ /= other_unit_obj._value_;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator/(const brgastro::unit_obj &other_unit_obj) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result /= other_unit_obj;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator/=(double new__value_ue)
{
	_value_ /= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator/(double new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result /= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator/=(long double new__value_ue)
{
	_value_ /= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator/(long double new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result /= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator/=(float new__value_ue)
{
	_value_ /= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator/(float new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result /= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator/=(int new__value_ue)
{
	_value_ /= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator/(int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result /= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator/=(long int new__value_ue)
{
	_value_ /= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator/(long int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result /= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator/=(unsigned int new__value_ue)
{
	_value_ /= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator/(unsigned int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result /= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator/=(short unsigned int new__value_ue)
{
	_value_ /= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator/(short unsigned int new__value_ue) const
{
	brgastro::unit_obj result = brgastro::unit_obj(*this,true);
	result /= new__value_ue;
	return result;
}
brgastro::unit_obj & brgastro::unit_obj::operator/=(long unsigned int new__value_ue)
{
	_value_ /= new__value_ue;
	return *this;
}
brgastro::unit_obj brgastro::unit_obj::operator/(long unsigned int new__value_ue) const
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
brgastro::unit_obj brgastro::unit_obj::operator++(int)   //Postfix
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
brgastro::unit_obj brgastro::unit_obj::operator--(int)   //Postfix
{
	brgastro::unit_obj obj_copy = brgastro::unit_obj(*this);
	_value_ -= 1;
	return obj_copy;
}
brgastro::unit_obj brgastro::unit_obj::operator-() const   //Prefix
{
	brgastro::unit_obj result = brgastro::unit_obj(*this);
	return result *= -1;
}
bool brgastro::unit_obj::operator<(const brgastro::unit_obj & other_unit_obj) const
{
	return (_value_ < other_unit_obj._value_);
}
bool brgastro::unit_obj::operator<(int comp__value_) const
{
	return (_value_ < comp__value_);
}
bool brgastro::unit_obj::operator<(long int comp__value_) const
{
	return (_value_ < comp__value_);
}
bool brgastro::unit_obj::operator<(unsigned int comp__value_) const
{
	return (_value_ < comp__value_);
}
bool brgastro::unit_obj::operator<(short unsigned int comp__value_) const
{
	return (_value_ < comp__value_);
}
bool brgastro::unit_obj::operator<(long unsigned int comp__value_) const
{
	return (_value_ < comp__value_);
}
bool brgastro::unit_obj::operator<(double comp__value_) const
{
	return (_value_ < comp__value_);
}
bool brgastro::unit_obj::operator<(long double comp__value_) const
{
	return (_value_ < comp__value_);
}
bool brgastro::unit_obj::operator<(float comp__value_) const
{
	return (_value_ < comp__value_);
}
bool brgastro::unit_obj::operator>(const brgastro::unit_obj & other_unit_obj) const
{
	return (_value_ > other_unit_obj._value_);
}
bool brgastro::unit_obj::operator>(int comp__value_) const
{
	return (_value_ > comp__value_);
}
bool brgastro::unit_obj::operator>(long int comp__value_) const
{
	return (_value_ > comp__value_);
}
bool brgastro::unit_obj::operator>(unsigned int comp__value_) const
{
	return (_value_ > comp__value_);
}
bool brgastro::unit_obj::operator>(short unsigned int comp__value_) const
{
	return (_value_ > comp__value_);
}
bool brgastro::unit_obj::operator>(long unsigned int comp__value_) const
{
	return (_value_ > comp__value_);
}
bool brgastro::unit_obj::operator>(double comp__value_) const
{
	return (_value_ > comp__value_);
}
bool brgastro::unit_obj::operator>(long double comp__value_) const
{
	return (_value_ > comp__value_);
}
bool brgastro::unit_obj::operator>(float comp__value_) const
{
	return (_value_ > comp__value_);
}
bool brgastro::unit_obj::operator==(const brgastro::unit_obj & other_unit_obj) const
{
	return (_value_ == other_unit_obj._value_);
}
bool brgastro::unit_obj::operator==(int comp__value_) const
{
	return (_value_ == comp__value_);
}
bool brgastro::unit_obj::operator==(long int comp__value_) const
{
	return (_value_ == comp__value_);
}
bool brgastro::unit_obj::operator==(unsigned int comp__value_) const
{
	return (_value_ == comp__value_);
}
bool brgastro::unit_obj::operator==(short unsigned int comp__value_) const
{
	return (_value_ == comp__value_);
}
bool brgastro::unit_obj::operator==(long unsigned int comp__value_) const
{
	return (_value_ == comp__value_);
}
bool brgastro::unit_obj::operator==(double comp__value_) const
{
	return (_value_ == comp__value_);
}
bool brgastro::unit_obj::operator==(long double comp__value_) const
{
	return (_value_ == comp__value_);
}
bool brgastro::unit_obj::operator==(float comp__value_) const
{
	return (_value_ == comp__value_);
}
bool brgastro::unit_obj::operator<=(const brgastro::unit_obj & other_unit_obj) const
{
	return (_value_ <= other_unit_obj._value_);
}
bool brgastro::unit_obj::operator<=(int comp__value_) const
{
	return (_value_ <= comp__value_);
}
bool brgastro::unit_obj::operator<=(long int comp__value_) const
{
	return (_value_ <= comp__value_);
}
bool brgastro::unit_obj::operator<=(unsigned int comp__value_) const
{
	return (_value_ <= comp__value_);
}
bool brgastro::unit_obj::operator<=(short unsigned int comp__value_) const
{
	return (_value_ <= comp__value_);
}
bool brgastro::unit_obj::operator<=(long unsigned int comp__value_) const
{
	return (_value_ <= comp__value_);
}
bool brgastro::unit_obj::operator<=(double comp__value_) const
{
	return (_value_ <= comp__value_);
}
bool brgastro::unit_obj::operator<=(long double comp__value_) const
{
	return (_value_ <= comp__value_);
}
bool brgastro::unit_obj::operator<=(float comp__value_) const
{
	return (_value_ <= comp__value_);
}
bool brgastro::unit_obj::operator>=(const brgastro::unit_obj & other_unit_obj) const
{
	return (_value_ >= other_unit_obj._value_);
}
bool brgastro::unit_obj::operator>=(int comp__value_) const
{
	return (_value_ >= comp__value_);
}
bool brgastro::unit_obj::operator>=(long int comp__value_) const
{
	return (_value_ >= comp__value_);
}
bool brgastro::unit_obj::operator>=(unsigned int comp__value_) const
{
	return (_value_ >= comp__value_);
}
bool brgastro::unit_obj::operator>=(short unsigned int comp__value_) const
{
	return (_value_ >= comp__value_);
}
bool brgastro::unit_obj::operator>=(long unsigned int comp__value_) const
{
	return (_value_ >= comp__value_);
}
bool brgastro::unit_obj::operator>=(double comp__value_) const
{
	return (_value_ >= comp__value_);
}
bool brgastro::unit_obj::operator>=(long double comp__value_) const
{
	return (_value_ >= comp__value_);
}
bool brgastro::unit_obj::operator>=(float comp__value_) const
{
	return (_value_ >= comp__value_);
}
bool brgastro::unit_obj::operator!=(const brgastro::unit_obj & other_unit_obj) const
{
	return (_value_ != other_unit_obj._value_);
}
bool brgastro::unit_obj::operator!=(int comp__value_) const
{
	return (_value_ != comp__value_);
}
bool brgastro::unit_obj::operator!=(long int comp__value_) const
{
	return (_value_ != comp__value_);
}
bool brgastro::unit_obj::operator!=(unsigned int comp__value_) const
{
	return (_value_ != comp__value_);
}
bool brgastro::unit_obj::operator!=(short unsigned int comp__value_) const
{
	return (_value_ != comp__value_);
}
bool brgastro::unit_obj::operator!=(long unsigned int comp__value_) const
{
	return (_value_ != comp__value_);
}
bool brgastro::unit_obj::operator!=(double comp__value_) const
{
	return (_value_ != comp__value_);
}
bool brgastro::unit_obj::operator!=(long double comp__value_) const
{
	return (_value_ != comp__value_);
}
bool brgastro::unit_obj::operator!=(float comp__value_) const
{
	return (_value_ != comp__value_);
}

brgastro::unit_obj::operator double() const
{
	return _value_;
}

// Overloaded operators relating to unit_objs

brgastro::unit_obj operator+(int lhs, const brgastro::unit_obj & rhs)
{
	return rhs+lhs;
}
brgastro::unit_obj operator+(long int lhs, const brgastro::unit_obj & rhs)
{
	return rhs+lhs;
}
brgastro::unit_obj operator+(double lhs, const brgastro::unit_obj & rhs)
{
	return rhs+lhs;
}
brgastro::unit_obj operator+(long double lhs, const brgastro::unit_obj & rhs)
{
	return rhs+lhs;
}
brgastro::unit_obj operator+(float lhs, const brgastro::unit_obj & rhs)
{
	return rhs+lhs;
}
brgastro::unit_obj operator-(int lhs, const brgastro::unit_obj & rhs)
{
	return -(rhs-lhs);
}
brgastro::unit_obj operator-(long int lhs, const brgastro::unit_obj & rhs)
{
	return -(rhs-lhs);
}
brgastro::unit_obj operator-(double lhs, const brgastro::unit_obj & rhs)
{
	return -(rhs-lhs);
}
brgastro::unit_obj operator-(long double lhs, const brgastro::unit_obj & rhs)
{
	return -(rhs-lhs);
}
brgastro::unit_obj operator-(float lhs, const brgastro::unit_obj & rhs)
{
	return -(rhs-lhs);
}
brgastro::unit_obj operator*(int lhs, const brgastro::unit_obj & rhs)
{
	return rhs*lhs;
}
brgastro::unit_obj operator*(long int lhs, const brgastro::unit_obj & rhs)
{
	return rhs*lhs;
}
brgastro::unit_obj operator*(double lhs, const brgastro::unit_obj & rhs)
{
	return rhs*lhs;
}
brgastro::unit_obj operator*(long double lhs, const brgastro::unit_obj & rhs)
{
	return rhs*lhs;
}
brgastro::unit_obj operator*(float lhs, const brgastro::unit_obj & rhs)
{
	return rhs*lhs;
}
brgastro::unit_obj operator/(int lhs, const brgastro::unit_obj & rhs)
{
	return brgastro::unit_obj(lhs)/rhs;
}
brgastro::unit_obj operator/(long int lhs, const brgastro::unit_obj & rhs)
{
	return brgastro::unit_obj(lhs)/rhs;
}
brgastro::unit_obj operator/(double lhs, const brgastro::unit_obj & rhs)
{
	return brgastro::unit_obj(lhs)/rhs;
}
brgastro::unit_obj operator/(long double lhs, const brgastro::unit_obj & rhs)
{
	return brgastro::unit_obj(lhs)/rhs;
}
brgastro::unit_obj operator/(float lhs, const brgastro::unit_obj & rhs)
{
	return brgastro::unit_obj(lhs)/rhs;
}

bool operator<(double lhs, const brgastro::unit_obj & rhs)
{
	return rhs > lhs;
}
bool operator<(long double lhs, const brgastro::unit_obj & rhs)
{
	return rhs > lhs;
}
bool operator<(float lhs, const brgastro::unit_obj & rhs)
{
	return rhs > lhs;
}
bool operator<(int lhs, const brgastro::unit_obj & rhs)
{
	return rhs > lhs;
}
bool operator<(unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs > lhs;
}
bool operator<(short unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs > lhs;
}
bool operator<(long unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs > lhs;
}
bool operator>(double lhs, const brgastro::unit_obj & rhs)
{
	return rhs < lhs;
}
bool operator>(long double lhs, const brgastro::unit_obj & rhs)
{
	return rhs < lhs;
}
bool operator>(float lhs, const brgastro::unit_obj & rhs)
{
	return rhs < lhs;
}
bool operator>(int lhs, const brgastro::unit_obj & rhs)
{
	return rhs < lhs;
}
bool operator>(unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs < lhs;
}
bool operator>(short unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs < lhs;
}
bool operator>(long unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs < lhs;
}
bool operator<=(double lhs, const brgastro::unit_obj & rhs)
{
	return rhs >= lhs;
}
bool operator<=(long double lhs, const brgastro::unit_obj & rhs)
{
	return rhs >= lhs;
}
bool operator<=(float lhs, const brgastro::unit_obj & rhs)
{
	return rhs >= lhs;
}
bool operator<=(int lhs, const brgastro::unit_obj & rhs)
{
	return rhs >= lhs;
}
bool operator<=(unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs >= lhs;
}
bool operator<=(short unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs >= lhs;
}
bool operator<=(long unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs >= lhs;
}
bool operator>=(double lhs, const brgastro::unit_obj & rhs)
{
	return rhs <= lhs;
}
bool operator>=(long double lhs, const brgastro::unit_obj & rhs)
{
	return rhs <= lhs;
}
bool operator>=(float lhs, const brgastro::unit_obj & rhs)
{
	return rhs <= lhs;
}
bool operator>=(int lhs, const brgastro::unit_obj & rhs)
{
	return rhs <= lhs;
}
bool operator>=(unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs <= lhs;
}
bool operator>=(short unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs <= lhs;
}
bool operator>=(long unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs <= lhs;
}
bool operator==(double lhs, const brgastro::unit_obj & rhs)
{
	return rhs == lhs;
}
bool operator==(long double lhs, const brgastro::unit_obj & rhs)
{
	return rhs == lhs;
}
bool operator==(float lhs, const brgastro::unit_obj & rhs)
{
	return rhs == lhs;
}
bool operator==(int lhs, const brgastro::unit_obj & rhs)
{
	return rhs == lhs;
}
bool operator==(unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs == lhs;
}
bool operator==(short unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs == lhs;
}
bool operator==(long unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs == lhs;
}
bool operator!=(double lhs, const brgastro::unit_obj & rhs)
{
	return rhs != lhs;
}
bool operator!=(long double lhs, const brgastro::unit_obj & rhs)
{
	return rhs != lhs;
}
bool operator!=(float lhs, const brgastro::unit_obj & rhs)
{
	return rhs != lhs;
}
bool operator!=(int lhs, const brgastro::unit_obj & rhs)
{
	return rhs != lhs;
}
bool operator!=(unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs != lhs;
}
bool operator!=(short unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs != lhs;
}
bool operator!=(long unsigned int lhs, const brgastro::unit_obj & rhs)
{
	return rhs != lhs;
}

std::ostream & operator<< (std::ostream &out, brgastro::unit_obj &obj)
{
	out << obj.get_value();
	return out;
}

brgastro::unit_obj std::pow(const brgastro::unit_obj & lhs, double rhs)
{
	brgastro::unit_obj result = brgastro::unit_obj(lhs);
	result.set_value(std::pow(lhs.get_value(),rhs));
	for(size_t i = 0; i < NUM_UNIT_TYPES; i++) result._unit_powers_[i]*=(float)rhs;
	result.round_powers();
	return result;
}
brgastro::unit_obj std::pow(const brgastro::unit_obj & lhs, long double rhs)
{
	brgastro::unit_obj result = brgastro::unit_obj(lhs);
	result.set_value(std::pow(lhs.get_value(),rhs));
	for(size_t i = 0; i < NUM_UNIT_TYPES; i++) result._unit_powers_[i]*=(float)rhs;
	result.round_powers();
	return result;
}
brgastro::unit_obj std::pow(const brgastro::unit_obj & lhs, float rhs)
{
	brgastro::unit_obj result = brgastro::unit_obj(lhs);
	result.set_value(std::pow(lhs.get_value(),rhs));
	for(size_t i = 0; i < NUM_UNIT_TYPES; i++) result._unit_powers_[i]*=rhs;
	result.round_powers();
	return result;
}
brgastro::unit_obj std::pow(const brgastro::unit_obj & lhs, int rhs)
{
	brgastro::unit_obj result = brgastro::unit_obj(lhs);
	result.set_value(std::pow(lhs.get_value(),rhs));
	for(size_t i = 0; i < NUM_UNIT_TYPES; i++) result._unit_powers_[i]*=rhs;
	result.round_powers();
	return result;
}
brgastro::unit_obj std::pow(const brgastro::unit_obj & lhs, unsigned int rhs)
{
	brgastro::unit_obj result = brgastro::unit_obj(lhs);
	result.set_value(std::pow(lhs.get_value(),rhs));
	for(size_t i = 0; i < NUM_UNIT_TYPES; i++) result._unit_powers_[i]*=rhs;
	result.round_powers();
	return result;
}
brgastro::unit_obj std::pow(const brgastro::unit_obj & lhs, short unsigned int rhs)
{
	brgastro::unit_obj result = brgastro::unit_obj(lhs);
	result.set_value(std::pow(lhs.get_value(),rhs));
	for(size_t i = 0; i < NUM_UNIT_TYPES; i++) result._unit_powers_[i]*=rhs;
	result.round_powers();
	return result;
}
brgastro::unit_obj std::pow(const brgastro::unit_obj & lhs, long unsigned int rhs)
{
	brgastro::unit_obj result = brgastro::unit_obj(lhs);
	result.set_value(std::pow(lhs.get_value(),rhs));
	for(size_t i = 0; i < NUM_UNIT_TYPES; i++) result._unit_powers_[i]*=rhs;
	result.round_powers();
	return result;
}
brgastro::unit_obj std::sqrt(const brgastro::unit_obj & obj)
{
	brgastro::unit_obj result = brgastro::unit_obj(obj);
	result.set_value(std::sqrt(obj.get_value()));
	for(size_t i = 0; i < NUM_UNIT_TYPES; i++) result._unit_powers_[i]/=2.;
	result.round_powers();
	return result;
}
brgastro::unit_obj std::fabs(const brgastro::unit_obj & obj)
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

double brgastro::unit_temperature::K() const
{
	return _value_;
}

double brgastro::unit_temperature::C() const
{
	return _value_+ABS_ZERO_C;
}

double brgastro::unit_temperature::R() const
{
	return _value_*KtodegR;
}

double brgastro::unit_temperature::F() const
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
