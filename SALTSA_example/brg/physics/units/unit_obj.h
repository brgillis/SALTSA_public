/**********************************************************************\
  @file unit_obj.h

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

// body file: brg/physics/units/unit_obj.cpp

/**********************************************************************\
 unit_obj.h
 -----------

 brgastro::unit_obj declarations: A class for a double with units. This
 section is only included if the header directive _BRG_USE_UNITS_ is defined.
 This can be changed in the file global.h. If this is the case, the
 file units.cpp must be included and compiled with the project. This
 declares the following classes:

 brgastro::unit_obj - Double with generic units

 brgastro::unit_distance    - Child of unit_obj, with units fixed to
 distance
 brgastro::unit_time        - " fixed to time
 brgastro::unit_mass        - " fixed to mass
 brgastro::unit_temperature - " fixed to temp
 brgastro::unit_angle       - " fixed to angle
 brgastro::unit_charge      - " fixed to charge
 brgastro::unit_velocity    - " fixed to distance/time


 The class internally stores the value of the variable in default units
 and the powers of each unit type in a size-6 vector. All operations are
 performed in the default unit set, except for "set" and "get" functions
 which specify alternate units to use.

 I've found that using this class is most useful for debugging, as unit
 mismatches raise an obvious flag when something is wrong. As use of this
 class does slow the program down (order of ~50% longer for some programs
 I've run), it might be recommended to disable its use once debugging is
 complete, by toggling the _BRG_USE_UNITS_ flag in the brg_global.h file.

 As programmed, the class allows only basic units of distance, time,
 mass, temperature, angle, and charge. This can be expanded if necessary
 or convenient, but be aware of possible backwards-compatibility issues
 from changing the function formats.


 The class defaults to kms units for operations with unitless values. For
 instance, the expression:
 unit_distance d = 2;
 will assign d a value of 2 metres. The full default unit set is:

 Distance:    Metres    (m)
 Time:        Seconds   (s)
 Mass:        Kilograms (kg)
 Temperature: Kelvins   (K)
 Angle:       Radians   (rad)
 Charge:      Coulombs  (C)


 Functions for unit_objs:


 Key things to remember for these functions:
 -When in doubt, use default units (see list above), and expect
 unit_objs to behave like doubles.
 -For any unit conversions, use the "specified-to-default" format
 (unitconv::kmtom, not unitconv::mtokm). Do not raise the unit
 conversion factor to any power when you're giving it as a
 function argument, except for when using the
 set_value(new_val, conv_factor) function.
 -Values for angle and charge units/unit powers can always be left
 off if not needed; they'll take default values for all functions.


 brgastro::unit_obj x;

 Initializes x with a value of zero and unitless.


 brgastro::unit_obj x(const double init_val=0,
 const double d_units_power=0, const double t_units_power=0,
 const double m_units_power=0, const double T_units_power=0,
 const double a_units_power=0, const double c_units_power=0);

 Initializes x with a value of init_val and specified unit powers. For
 instance, brgastro::unit_obj x(4,2) will initializes x with as 4 m^2.


 brgastro::unit_obj x(const double init_val,
 const double d_units, const float d_units_power,
 const double t_units, const float t_units_power,
 const double m_units, const float m_units_power,
 const double T_units, const float T_units_power,
 const double a_units=0, const float a_units_power=0,
 const double c_units=0, const float c_units_power=0);

 Initializes x with a value of init_val in the specified units, with the
 specified unit powers. For instance,
 brgastro::unit_obj x( 4, unitconv::kmtom, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0 )
 will initialize x with a value of 4 km^2, which is internally stored as
 4e6 m^2. Always use the "kmtom"-like form of the unit conversion, with
 the default unit on the right, for this and all in all other class
 functions which use specified units.


 brgastro::unit_obj x(const unit_obj& other_unit_obj,
 const bool maintain_unit_fix=false);

 Initializes x as a copy of other_unit_obj. If maintain_unit_fix is set
 to true, and other_unit_obj is a value with fixed units (eg. it was
 declared as a unit_distance), x will have fixed units as well. Generally
 this is not necessary, and this value can be ignored for any explicit
 declarations.


 Set functions:


 const int brgastro::unit_obj::reset(const float d_units_power=0,
 const float t_units_power=0, const float m_units_power=0,
 const float T_units_power=0, const float a_units_power=0,
 const float c_units_power=0);

 Resets the value and optionally alters the unit powers, otherwise
 setting them to zero. Always returns a value of 0, but
 can be defined otherwise if necessary.


 const int brgastro::unit_obj::set_unit_powers(
 const float d_units_power=0, const float t_units_power=0,
 const float m_units_power=0, const float T_units_power=0,
 const float a_units_power=0, const float c_units_power=0);
 const int brgastro::unit_obj::set_unit_powers(
 std::vector<float> new_unitpowers);

 Set the unit powers of the variable, but DO NOT affect the stored value.
 Not recommended to be used, but available if necessary. Either function
 will return 1 if the variable's unit powers are fixed, and the latter
 will return 1 if the vector does not have length equal to the number
 of unit types.


 const int brgastro::unit_obj::set_value(const double new_val,
 const double conv_factor=1);

 Set the value of the variable in default units, or in the units
 specified by the conv_factor parameter. For instance,
 x.set_value(4, unitconv::kmtom/squarew(unitconv::hrtos)) will assign
 x a value of 4 km/hr^2, but it WILL NOT change the unit powers of
 x correspondingly. Use this only if you are sure x already has
 the correct unit powers. Always returns a value of 0.


 const int brgastro::unit_obj::set_value(const double new_val,
 const double d_units, const double t_units, const double m_units,
 const double T_units, const double a_units=1,
 const double c_units=1);

 As above, but units are set individually. DO NOT raise units to
 the appropriate power; the function determines it automatically
 from the stored unit powers. For instance,
 x.set_value(4, unitconv::kmtom, unitconv::hrtos, 1, 1, 1, 1)
 will assign x a value of 4 km/hr if x has unit powers of distance/time,
 but it will assign x a value of 4 km/hr^2 if x has unit powers of
 distance/time^2.


 const int brgastro::unit_obj::set(const double new_val,
 const std::vector<float> new_unitpowers);
 const int brgastro::unit_obj::set(const double new_val,
 const float d_units_power, const float t_units_power,
 const float m_units_power, const float T_units_power,
 const float a_units_power=0, const float c_units_power=0);

 Sets the value in default units and changes the unit powers of the
 variable.


 const int brgastro::unit_obj::set(const double new_val,
 const double d_units, const float d_units_power,
 const double t_units, const float t_units_power,
 const double m_units, const float m_units_power,
 const double T_units, const float T_units_power,
 const double a_units=1, const float a_units_power=0,
 const double c_units=1, const float c_units_power=0);

 Sets the value in the specified units and changes the unit powers of the
 variable. DO NOT raise units to the appropriate power; the function
 determines it automatically from the specified unit powers.


 Get functions:


 const double brgastro::unit_obj::get_value() const;

 Returns the value in default units.


 const double brgastro::unit_obj::get_value(const double d_units,
 const double t_units, const double m_units, const double T_units,
 const double a_units=1, const double C_units=1) const;

 Returns the value in specified units. DO NOT raise units to
 the appropriate power; the function determines it automatically
 from the stored unit powers.


 const double brgastro::unit_obj::get_value(
 const double d_units, const float d_units_power,
 const double t_units, const float t_units_power,
 const double m_units, const float m_units_power,
 const double T_units, const float T_units_power,
 const double a_units=1, const float a_units_power=0,
 const double c_units=1, const float c_units_power=0) const;

 Returns the value in specified units, using specified unit powers
 instead of the stored powers. Typically should not be needed.


 const std::vector<float> brgastro::unit_obj::get_unit_powers() const;

 Returns a vector of the unit powers of the variable.


 const std::string brgastro::unit_obj::get_string() const;

 Returns a string writing out the variable's value in default units,
 along with the powers of those units.


 Specific get functions:
 These functions automatically adjust for the unit powers. For instance,
 if x has a value of 100 m^2, x.km() will return a value of 0.0001 (km^2).

 Get distance in units of...
 const double brgastro::unit_obj::m() const;
 const double brgastro::unit_obj::mm() const;
 const double brgastro::unit_obj::um() const;
 const double brgastro::unit_obj::nm() const;
 const double brgastro::unit_obj::cm() const;
 const double brgastro::unit_obj::angstrom() const;
 const double brgastro::unit_obj::km() const;
 const double brgastro::unit_obj::ltyr() const;
 const double brgastro::unit_obj::pc() const;
 const double brgastro::unit_obj::kpc() const;
 const double brgastro::unit_obj::Mpc() const;
 const double brgastro::unit_obj::mi() const;
 const double brgastro::unit_obj::Mmi() const;
 const double brgastro::unit_obj::ft() const;
 const double brgastro::unit_obj::yd() const;

 Get time in units of...
 const double brgastro::unit_obj::s() const;
 const double brgastro::unit_obj::ms() const;
 const double brgastro::unit_obj::cs() const;
 const double brgastro::unit_obj::ns() const;
 const double brgastro::unit_obj::us() const;
 const double brgastro::unit_obj::min() const;
 const double brgastro::unit_obj::hr() const;
 const double brgastro::unit_obj::day() const;
 const double brgastro::unit_obj::week() const;
 const double brgastro::unit_obj::month() const;
 const double brgastro::unit_obj::yr() const;
 const double brgastro::unit_obj::kyr() const;
 const double brgastro::unit_obj::Myr() const;
 const double brgastro::unit_obj::Gyr() const;

 Get mass in units of...
 const double brgastro::unit_obj::kg() const;
 const double brgastro::unit_obj::gm() const;
 const double brgastro::unit_obj::Mearth() const;
 const double brgastro::unit_obj::Msun() const;
 const double brgastro::unit_obj::ttMsun() const; // 10^10 Msun

 Get temperature in units of...
 Note: These automatically adjust for the different
 zero-values of the different temperature unit scales
 if the value is a pure temperature. If scaling only is
 preferred, simply use only the K() and degR() functions.
 const double brgastro::unit_obj::K() const ;
 const double brgastro::unit_obj::degC() const ;
 const double brgastro::unit_obj::degF() const;
 const double brgastro::unit_obj::degR() const;

 Get angle in units of...
 const double brgastro::unit_obj::rad() const;
 const double brgastro::unit_obj::deg() const;
 const double brgastro::unit_obj::amin() const;
 const double brgastro::unit_obj::asec() const;

 Get charge in units of...
 const double brgastro::unit_obj::C() const;
 const double brgastro::unit_obj::esu() const;

 Get other values in units of...
 Note: These functions do not automatically adjust for the variable's
 unit powers, since they involve more than one unit type.

 const double brgastro::unit_obj::mps() const; // For velocity
 const double brgastro::unit_obj::kmps() const; // ''
 const double brgastro::unit_obj::c() const; // ''
 const double brgastro::unit_obj::miphr() const; // ''
 const double brgastro::unit_obj::kpcpGyr() const; // ''

 const double brgastro::unit_obj::kmpspGyr() const; // For acceleration
 const double brgastro::unit_obj::kpcpGyr2() const; // ''


 Other functions:

 const int round_powers();

 Rounds the unit powers. Precision is determined by the
 power_round_precision value. By default, it rounds so that no unit
 power will have a denominator greater than 12. Use this if round-off
 error is a concern after raising a variable to a fractional power.

 \**********************************************************************/

#ifndef _BRG_UNIT_OBJ_H_INCLUDED_
#define _BRG_UNIT_OBJ_H_INCLUDED_

#ifdef _BRG_USE_UNITS_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "brg/global.h"

#include "brg/physics/units/unit_conversions.hpp"

#define DIST_UNIT_INDEX 0
#define TIME_UNIT_INDEX 1
#define MASS_UNIT_INDEX 2
#define TEMP_UNIT_INDEX 3
#define ANGL_UNIT_INDEX 4
#define CHRG_UNIT_INDEX 5
#define NUM_UNIT_TYPES 6
#define ABS_ZERO_C -273.15
#define ABS_ZERO_F -459.67
// Power round precision - Set to zero to enforce no rounding. Higher is finer rounding
#define UNIT_POWER_ROUND_PRECISION_FACTOR 27720 // Lowest number that's evenly divisible by every integer 1 through 12

namespace brgastro
{
	class unit_obj;
}

namespace std
{

// Some specialisations to add to the std:: namespace
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const int rhs);
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const short int rhs);
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const long int rhs);
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const unsigned int rhs);
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const short unsigned int rhs);
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const long unsigned int rhs);
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const double rhs);
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const long double rhs);
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const float rhs);
brgastro::unit_obj sqrt(const brgastro::unit_obj &obj);
brgastro::unit_obj fabs(const brgastro::unit_obj &obj);

}

namespace brgastro
{

	class unit_obj
	{
	protected:

		// Default units: m, s, kg, K, rad, C
		std::vector<float> _unit_powers_;// Distance, Time, Mass, Temperature, Angle
		double _value_;// Value in default units

		bool _fix_unit_powers_;// Used by derived classes to keep unit powers from being altered

	public:

		//------------------------------
		// Function prototypes (unit_obj)
		//------------------------------

		// Constructors
		unit_obj(double init_val=0, float d_units_power=0, float t_units_power=0, float m_units_power=0,
				float T_units_power=0, float a_units_power=0, float c_units_power=0);
		unit_obj(double init_val,
				double d_units, float d_units_power, double t_units, float t_units_power,
				double m_units, float m_units_power, double T_units, float T_units_power,
				 double a_units=1, float a_units_power=0, double c_units=1, float c_units_power=0);

		// Copy constructor
		unit_obj(const unit_obj& other_unit_obj, bool maintain_unit_fix=false);

		// Virtual Destructor (in case a derived class might need this functionality)
		virtual ~unit_obj();

		// Reset function
		void reset(float d_units_power=0, float t_units_power=0, float m_units_power=0,
				float T_units_power=0, float a_units_power=0, float c_units_power=0);

		// Set functions
		void set_unit_powers(float d_units_power=0, float t_units_power=0, float m_units_power=0,
				float T_units_power=0, float a_units_power=0, float c_units_power=0);// Use with care
		void set_unit_powers( std::vector<float> new_unitpowers);// Use with care
		void set_value(double new_val, double conv_factor=1);
		void set_value(double new_val, double d_units, double t_units, double m_units, double T_units,
				double a_units=1, double c_units=1);
		void set(double new_val, const std::vector<float> & new_unitpowers);
		void set(double new_val, float d_units_power, float t_units_power, float m_units_power,
				float T_units_power, float a_units_power=0, float c_units_power=0);
		void set(double new_val, double d_units, float d_units_power, double t_units, float t_units_power,
				double m_units, float m_units_power, double T_units, float T_units_power, double a_units=1,
				float a_units_power=0, double c_units=1, float c_units_power=0);

		// Get functions
		double get_value() const;
		double get_value(double d_units, double t_units, double m_units, double T_units, double a_units=1,
				double C_units=1) const;
		double get_value(double d_units, float d_units_power, double t_units, float t_units_power, double m_units,
				float m_units_power, double T_units, float T_units_power, double a_units=1,
				float a_units_power=0, double c_units=1, float c_units_power=0) const;// Use with care
		std::vector<float> get_unit_powers() const;
		std::string get_string() const;

		// Misc functions
		void round_powers();

		// Get distance in units of...
		double m() const;
		double mm() const;
		double um() const;
		double nm() const;
		double cm() const;
		double angstrom() const;
		double km() const;
		double ltyr() const;
		double pc() const;
		double kpc() const;
		double Mpc() const;
		double mi() const;
		double Mmi() const;
		double ft() const;
		double yd() const;

		// Get time in units of...
		double s() const;
		double ms() const;
		double cs() const;
		double ns() const;
		double us() const;
		double min() const;
		double hr() const;
		double day() const;
		double week() const;
		double month() const;
		double yr() const;
		double kyr() const;
		double Myr() const;
		double Gyr() const;

		// Get mass in units of...
		double kg() const;
		double gm() const;
		double Mearth() const;
		double Msun() const;
		double ttMsun() const;

		// Get temperature in units of...
		double K() const;
		double degC() const;
		double degF() const;
		double degR() const;

		// Get angle in units of...
		double rad() const;
		double deg() const;
		double amin() const;
		double asec() const;

		// Get charge in units of...
		double C() const;
		double esu() const;

		// Get other values in units of...

		double mps() const;// For velocity
		double kmps() const;// ''
		double c() const;// ''
		double miphr() const;// ''
		double kpcpGyr() const;// ''

		double kmpspGyr() const;// For acceleration
		double kpcpGyr2() const;// ''

		// Operator overloading
		unit_obj & operator=(const unit_obj & );
		unit_obj & operator=(int);
		unit_obj & operator=(short int);
		unit_obj & operator=(long int);
		unit_obj & operator=(unsigned int);
		unit_obj & operator=(short unsigned int);
		unit_obj & operator=(long unsigned int);
		unit_obj & operator=(double);
		unit_obj & operator=(long double);
		unit_obj & operator=(float);
		unit_obj operator+(const unit_obj & ) const;
		unit_obj operator+(int) const;
		unit_obj operator+(short int) const;
		unit_obj operator+(long int) const;
		unit_obj operator+(unsigned int) const;
		unit_obj operator+(short unsigned int) const;
		unit_obj operator+(long unsigned int) const;
		unit_obj operator+(double) const;
		unit_obj operator+(long double) const;
		unit_obj operator+(float) const;
		unit_obj operator-(const unit_obj & ) const;
		unit_obj operator-(int) const;
		unit_obj operator-(short int) const;
		unit_obj operator-(long int) const;
		unit_obj operator-(unsigned int) const;
		unit_obj operator-(short unsigned int) const;
		unit_obj operator-(long unsigned int) const;
		unit_obj operator-(double) const;
		unit_obj operator-(long double) const;
		unit_obj operator-(float) const;
		unit_obj operator*(const unit_obj & ) const;
		unit_obj operator*(int) const;
		unit_obj operator*(short int) const;
		unit_obj operator*(long int) const;
		unit_obj operator*(unsigned int) const;
		unit_obj operator*(short unsigned int) const;
		unit_obj operator*(long unsigned int) const;
		unit_obj operator*(double) const;
		unit_obj operator*(long double) const;
		unit_obj operator*(float) const;
		unit_obj operator/(const unit_obj & ) const;
		unit_obj operator/(int) const;
		unit_obj operator/(short int) const;
		unit_obj operator/(long int) const;
		unit_obj operator/(unsigned int) const;
		unit_obj operator/(short unsigned int) const;
		unit_obj operator/(long unsigned int) const;
		unit_obj operator/(double) const;
		unit_obj operator/(long double) const;
		unit_obj operator/(float) const;
		unit_obj & operator+=(const unit_obj & );
		unit_obj & operator+=(int);
		unit_obj & operator+=(short int);
		unit_obj & operator+=(long int);
		unit_obj & operator+=(unsigned int);
		unit_obj & operator+=(short unsigned int);
		unit_obj & operator+=(long unsigned int);
		unit_obj & operator+=(double);
		unit_obj & operator+=(long double);
		unit_obj & operator+=(float);
		unit_obj & operator-=(const unit_obj & );
		unit_obj & operator-=(int);
		unit_obj & operator-=(short int);
		unit_obj & operator-=(long int);
		unit_obj & operator-=(unsigned int);
		unit_obj & operator-=(short unsigned int);
		unit_obj & operator-=(long unsigned int);
		unit_obj & operator-=(double);
		unit_obj & operator-=(long double);
		unit_obj & operator-=(float);
		unit_obj & operator*=(const unit_obj & );
		unit_obj & operator*=(int);
		unit_obj & operator*=(short int);
		unit_obj & operator*=(long int);
		unit_obj & operator*=(unsigned int);
		unit_obj & operator*=(short unsigned int);
		unit_obj & operator*=(long unsigned int);
		unit_obj & operator*=(double);
		unit_obj & operator*=(long double);
		unit_obj & operator*=(float);
		unit_obj & operator/=(const unit_obj & );
		unit_obj & operator/=(int);
		unit_obj & operator/=(short int);
		unit_obj & operator/=(long int);
		unit_obj & operator/=(unsigned int);
		unit_obj & operator/=(short unsigned int);
		unit_obj & operator/=(long unsigned int);
		unit_obj & operator/=(double);
		unit_obj & operator/=(long double);
		unit_obj & operator/=(float);
		unit_obj & operator++();
		unit_obj operator++(int);
		unit_obj & operator--();
		unit_obj operator--(int);
		unit_obj operator-() const;
		bool operator<(const unit_obj & ) const;
		bool operator<(int) const;
		bool operator<(short int) const;
		bool operator<(long int) const;
		bool operator<(unsigned int) const;
		bool operator<(short unsigned int) const;
		bool operator<(long unsigned int) const;
		bool operator<(double) const;
		bool operator<(long double) const;
		bool operator<(float) const;
		bool operator>(const unit_obj & ) const;
		bool operator>(int) const;
		bool operator>(short int) const;
		bool operator>(long int) const;
		bool operator>(unsigned int) const;
		bool operator>(short unsigned int) const;
		bool operator>(long unsigned int) const;
		bool operator>(double) const;
		bool operator>(long double) const;
		bool operator>(float) const;
		bool operator==(const unit_obj & ) const;
		bool operator==(int) const;
		bool operator==(short int) const;
		bool operator==(long int) const;
		bool operator==(unsigned int) const;
		bool operator==(short unsigned int) const;
		bool operator==(long unsigned int) const;
		bool operator==(double) const;
		bool operator==(long double) const;
		bool operator==(float) const;
		bool operator<=(const unit_obj & ) const;
		bool operator<=(int) const;
		bool operator<=(short int) const;
		bool operator<=(long int) const;
		bool operator<=(unsigned int) const;
		bool operator<=(short unsigned int) const;
		bool operator<=(long unsigned int) const;
		bool operator<=(double) const;
		bool operator<=(long double) const;
		bool operator<=(float) const;
		bool operator>=(const unit_obj & ) const;
		bool operator>=(int) const;
		bool operator>=(short int) const;
		bool operator>=(long int) const;
		bool operator>=(unsigned int) const;
		bool operator>=(short unsigned int) const;
		bool operator>=(long unsigned int) const;
		bool operator>=(double) const;
		bool operator>=(long double) const;
		bool operator>=(float ) const;
		bool operator!=(const unit_obj & ) const;
		bool operator!=(int) const;
		bool operator!=(short int) const;
		bool operator!=(long int) const;
		bool operator!=(unsigned int) const;
		bool operator!=(short unsigned int) const;
		bool operator!=(long unsigned int) const;
		bool operator!=(double) const;
		bool operator!=(long double) const;
		bool operator!=(float) const;
		operator double() const;

		// Friends

		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, int rhs);
		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, short int rhs);
		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, long int rhs);
		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, unsigned int rhs);
		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, short unsigned int rhs);
		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, long unsigned int rhs);
		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, double rhs);
		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, long double rhs);
		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, float rhs);
		friend brgastro::unit_obj std::sqrt(const brgastro::unit_obj &obj);
		friend brgastro::unit_obj std::fabs(const brgastro::unit_obj &obj);

	}; // unit_obj class

// Distance class
	class unit_distance: public unit_obj
	{

	public:
		unit_distance();
		unit_distance(double init_val, double conv_factor=1);
		unit_distance(const unit_obj &other_unit_obj);

	}; // class distance

	class unit_time: public unit_obj
	{

	public:
		unit_time();
		unit_time(double init_val, double conv_factor=1);
		unit_time(const unit_obj &other_unit_obj);

	}; // class time

	class unit_velocity: public unit_obj
	{
	public:
		unit_velocity();
		unit_velocity(double init_val, double conv_factor=1);
		unit_velocity(const unit_obj &other_unit_obj);

	}; // class unit_velocity

	class unit_mass: public unit_obj
	{
	public:
		unit_mass();
		unit_mass(double init_val, double conv_factor=1);
		unit_mass(const unit_obj &other_unit_obj);

	}; // class mass

	class unit_temperature: public unit_obj
	{
	public:
		unit_temperature();
		unit_temperature(double init_val, double conv_factor=1);
		unit_temperature(const unit_obj &other_unit_obj);

		double K() const;
		double C() const;
		double R() const;
		double F() const;

	}; // class unit_temperature

	class unit_angle: public unit_obj
	{

	public:
		unit_angle();
		unit_angle(double init_val, double conv_factor=1);
		unit_angle(const unit_obj &other_unit_obj);

	}; // class unit_angle

	class unit_charge: public unit_obj
	{

	public:
		unit_charge();
		unit_charge(double init_val, double conv_factor=1);
		unit_charge(const unit_obj &other_unit_obj);

	}; // class unit_charge

	// brgastro function overloads for unit_objs
#if (1)

	inline void set_zero( brgastro::unit_obj & v1)
	{
		v1 = unit_obj(0);
	}
	inline void set_zero( brgastro::unit_distance & v1)
	{
		v1 = unit_obj(0);
	}
	inline void set_zero( brgastro::unit_time & v1)
	{
		v1 = unit_obj(0);
	}
	inline void set_zero( brgastro::unit_velocity & v1)
	{
		v1 = unit_obj(0);
	}
	inline void set_zero( brgastro::unit_mass & v1)
	{
		v1 = unit_obj(0);
	}
	inline void set_zero( brgastro::unit_temperature & v1)
	{
		v1 = unit_obj(0);
	}
	inline void set_zero( brgastro::unit_angle & v1)
	{
		v1 = unit_obj(0);
	}
	inline void set_zero( brgastro::unit_charge & v1)
	{
		v1 = unit_obj(0);
	}

	inline bool isinf( unit_obj val )
	{
		return std::fabs( val.get_value() ) > std::numeric_limits<double>::max();
	}

	inline brgastro::unit_obj square( const brgastro::unit_obj & v1)
	{
		return v1*v1;
	}
	inline brgastro::unit_obj square( const brgastro::unit_distance & v1)
	{
		return v1*v1;
	}
	inline brgastro::unit_obj square( const brgastro::unit_time & v1)
	{
		return v1*v1;
	}
	inline brgastro::unit_obj square( const brgastro::unit_velocity & v1)
	{
		return v1*v1;
	}
	inline brgastro::unit_obj square( const brgastro::unit_mass & v1)
	{
		return v1*v1;
	}
	inline brgastro::unit_obj square( const brgastro::unit_temperature & v1)
	{
		return v1*v1;
	}
	inline brgastro::unit_obj square( const brgastro::unit_angle & v1)
	{
		return v1*v1;
	}
	inline brgastro::unit_obj square( const brgastro::unit_charge & v1)
	{
		return v1*v1;
	}
	inline brgastro::unit_obj cube( const brgastro::unit_obj & v1)
	{
		return v1*square(v1);
	}
	inline brgastro::unit_obj cube( const brgastro::unit_distance & v1)
	{
		return v1*square(v1);
	}
	inline brgastro::unit_obj cube( const brgastro::unit_time & v1)
	{
		return v1*square(v1);
	}
	inline brgastro::unit_obj cube( const brgastro::unit_velocity & v1)
	{
		return v1*square(v1);
	}
	inline brgastro::unit_obj cube( const brgastro::unit_mass & v1)
	{
		return v1*square(v1);
	}
	inline brgastro::unit_obj cube( const brgastro::unit_temperature & v1)
	{
		return v1*square(v1);
	}
	inline brgastro::unit_obj cube( const brgastro::unit_angle & v1)
	{
		return v1*square(v1);
	}
	inline brgastro::unit_obj cube( const brgastro::unit_charge & v1)
	{
		return v1*square(v1);
	}
	inline brgastro::unit_obj quart( const brgastro::unit_obj & v1)
	{
		return square(v1)*square(v1);
	}
	inline brgastro::unit_obj quart( const brgastro::unit_distance & v1)
	{
		return square(v1)*square(v1);
	}
	inline brgastro::unit_obj quart( const brgastro::unit_time & v1)
	{
		return square(v1)*square(v1);
	}
	inline brgastro::unit_obj quart( const brgastro::unit_velocity & v1)
	{
		return square(v1)*square(v1);
	}
	inline brgastro::unit_obj quart( const brgastro::unit_mass & v1)
	{
		return square(v1)*square(v1);
	}
	inline brgastro::unit_obj quart( const brgastro::unit_temperature & v1)
	{
		return square(v1)*square(v1);
	}
	inline brgastro::unit_obj quart( const brgastro::unit_angle & v1)
	{
		return square(v1)*square(v1);
	}
	inline brgastro::unit_obj quart( const brgastro::unit_charge & v1)
	{
		return square(v1)*square(v1);
	}
	inline brgastro::unit_obj inverse( const brgastro::unit_obj & v1)
	{
		return 1/v1;
	}
	inline brgastro::unit_obj inverse( const brgastro::unit_distance & v1)
	{
		return 1/v1;
	}
	inline brgastro::unit_obj inverse( const brgastro::unit_time & v1)
	{
		return 1/v1;
	}
	inline brgastro::unit_obj inverse( const brgastro::unit_velocity & v1)
	{
		return 1/v1;
	}
	inline brgastro::unit_obj inverse( const brgastro::unit_mass & v1)
	{
		return 1/v1;
	}
	inline brgastro::unit_obj inverse( const brgastro::unit_temperature & v1)
	{
		return 1/v1;
	}
	inline brgastro::unit_obj inverse( const brgastro::unit_angle & v1)
	{
		return 1/v1;
	}
	inline brgastro::unit_obj inverse( const brgastro::unit_charge & v1)
	{
		return 1/v1;
	}
	inline brgastro::unit_obj inv_square( const brgastro::unit_obj & v1)
	{
		return inverse(square(v1));
	}
	inline brgastro::unit_obj inv_square( const brgastro::unit_distance & v1)
	{
		return inverse(square(v1));
	}
	inline brgastro::unit_obj inv_square( const brgastro::unit_time & v1)
	{
		return inverse(square(v1));
	}
	inline brgastro::unit_obj inv_square( const brgastro::unit_velocity & v1)
	{
		return inverse(square(v1));
	}
	inline brgastro::unit_obj inv_square( const brgastro::unit_mass & v1)
	{
		return inverse(square(v1));
	}
	inline brgastro::unit_obj inv_square( const brgastro::unit_temperature & v1)
	{
		return inverse(square(v1));
	}
	inline brgastro::unit_obj inv_square( const brgastro::unit_angle & v1)
	{
		return inverse(square(v1));
	}
	inline brgastro::unit_obj inv_square( const brgastro::unit_charge & v1)
	{
		return inverse(square(v1));
	}
	inline brgastro::unit_obj inv_cube( const brgastro::unit_obj & v1)
	{
		return inverse(cube(v1));
	}
	inline brgastro::unit_obj inv_cube( const brgastro::unit_distance & v1)
	{
		return inverse(cube(v1));
	}
	inline brgastro::unit_obj inv_cube( const brgastro::unit_time & v1)
	{
		return inverse(cube(v1));
	}
	inline brgastro::unit_obj inv_cube( const brgastro::unit_velocity & v1)
	{
		return inverse(cube(v1));
	}
	inline brgastro::unit_obj inv_cube( const brgastro::unit_mass & v1)
	{
		return inverse(cube(v1));
	}
	inline brgastro::unit_obj inv_cube( const brgastro::unit_temperature & v1)
	{
		return inverse(cube(v1));
	}
	inline brgastro::unit_obj inv_cube( const brgastro::unit_angle & v1)
	{
		return inverse(cube(v1));
	}
	inline brgastro::unit_obj inv_cube( const brgastro::unit_charge & v1)
	{
		return inverse(cube(v1));
	}
	inline brgastro::unit_obj inv_quart( const brgastro::unit_obj & v1)
	{
		return inverse(quart(v1));
	}
	inline brgastro::unit_obj inv_quart( const brgastro::unit_distance & v1)
	{
		return inverse(quart(v1));
	}
	inline brgastro::unit_obj inv_quart( const brgastro::unit_time & v1)
	{
		return inverse(quart(v1));
	}
	inline brgastro::unit_obj inv_quart( const brgastro::unit_velocity & v1)
	{
		return inverse(quart(v1));
	}
	inline brgastro::unit_obj inv_quart( const brgastro::unit_mass & v1)
	{
		return inverse(quart(v1));
	}
	inline brgastro::unit_obj inv_quart( const brgastro::unit_temperature & v1)
	{
		return inverse(quart(v1));
	}
	inline brgastro::unit_obj inv_quart( const brgastro::unit_angle & v1)
	{
		return inverse(quart(v1));
	}
	inline brgastro::unit_obj inv_quart( const brgastro::unit_charge & v1)
	{
		return inverse(quart(v1));
	}
	inline brgastro::unit_obj ipow( const brgastro::unit_obj & v1, int p)
	{
		return ipow((brgastro::unit_obj)v1,p);
	}
	inline brgastro::unit_obj ipow( const brgastro::unit_distance & v1, int p)
	{
		return ipow((brgastro::unit_obj)v1,p);
	}
	inline brgastro::unit_obj ipow( const brgastro::unit_time & v1, int p)
	{
		return ipow((brgastro::unit_obj)v1,p);
	}
	inline brgastro::unit_obj ipow( const brgastro::unit_velocity & v1, int p)
	{
		return ipow((brgastro::unit_obj)v1,p);
	}
	inline brgastro::unit_obj ipow( const brgastro::unit_mass & v1, int p)
	{
		return ipow((brgastro::unit_obj)v1,p);
	}
	inline brgastro::unit_obj ipow( const brgastro::unit_temperature & v1, int p)
	{
		return ipow((brgastro::unit_obj)v1,p);
	}
	inline brgastro::unit_obj ipow( const brgastro::unit_angle & v1, int p)
	{
		return ipow((brgastro::unit_obj)v1,p);
	}
	inline brgastro::unit_obj ipow( const brgastro::unit_charge & v1, int p)
	{
		return ipow((brgastro::unit_obj)v1,p);
	}
#endif

} // namespace brgastro

// Overloaded operators relating to unit_objs

brgastro::unit_obj operator+(int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator+(short int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator+(long int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator+(unsigned int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator+(short unsigned int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator+(long unsigned int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator+(double lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator+(long double lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator+(float lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(short int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(long int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(unsigned int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(short unsigned int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(long unsigned int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(double lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(long double lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(float lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(short int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(long int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(unsigned int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(short unsigned int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(long unsigned int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(double lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(long double lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(float lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(short int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(long int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(unsigned int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(short unsigned int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(long unsigned int lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(double lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(long double lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(float lhs, const brgastro::unit_obj &rhs);

bool operator<(int lhs, const brgastro::unit_obj &rhs);
bool operator<(short int lhs, const brgastro::unit_obj &rhs);
bool operator<(long int lhs, const brgastro::unit_obj &rhs);
bool operator<(unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator<(short unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator<(long unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator<(double lhs, const brgastro::unit_obj &rhs);
bool operator<(long double lhs, const brgastro::unit_obj &rhs);
bool operator<(float lhs, const brgastro::unit_obj &rhs);
bool operator>(int lhs, const brgastro::unit_obj &rhs);
bool operator>(short int lhs, const brgastro::unit_obj &rhs);
bool operator>(long int lhs, const brgastro::unit_obj &rhs);
bool operator>(unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator>(short unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator>(long unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator>(double lhs, const brgastro::unit_obj &rhs);
bool operator>(long double lhs, const brgastro::unit_obj &rhs);
bool operator>(float lhs, const brgastro::unit_obj &rhs);
bool operator<=(int lhs, const brgastro::unit_obj &rhs);
bool operator<=(short int lhs, const brgastro::unit_obj &rhs);
bool operator<=(long int lhs, const brgastro::unit_obj &rhs);
bool operator<=(unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator<=(short unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator<=(long unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator<=(double lhs, const brgastro::unit_obj &rhs);
bool operator<=(long double lhs, const brgastro::unit_obj &rhs);
bool operator<=(float lhs, const brgastro::unit_obj &rhs);
bool operator>=(int lhs, const brgastro::unit_obj &rhs);
bool operator>=(short int lhs, const brgastro::unit_obj &rhs);
bool operator>=(long int lhs, const brgastro::unit_obj &rhs);
bool operator>=(unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator>=(short unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator>=(long unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator>=(double lhs, const brgastro::unit_obj &rhs);
bool operator>=(long double lhs, const brgastro::unit_obj &rhs);
bool operator>=(float lhs, const brgastro::unit_obj &rhs);
bool operator==(int lhs, const brgastro::unit_obj &rhs);
bool operator==(short int lhs, const brgastro::unit_obj &rhs);
bool operator==(long int lhs, const brgastro::unit_obj &rhs);
bool operator==(unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator==(short unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator==(long unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator==(double lhs, const brgastro::unit_obj &rhs);
bool operator==(long double lhs, const brgastro::unit_obj &rhs);
bool operator==(float lhs, const brgastro::unit_obj &rhs);
bool operator!=(int lhs, const brgastro::unit_obj &rhs);
bool operator!=(short int lhs, const brgastro::unit_obj &rhs);
bool operator!=(long int lhs, const brgastro::unit_obj &rhs);
bool operator!=(unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator!=(short unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator!=(long unsigned int lhs, const brgastro::unit_obj &rhs);
bool operator!=(double lhs, const brgastro::unit_obj &rhs);
bool operator!=(long double lhs, const brgastro::unit_obj &rhs);
bool operator!=(float lhs, const brgastro::unit_obj &rhs);

std::ostream & operator<<(std::ostream &out, brgastro::unit_obj &obj);

#endif

#endif // _BRG_UNIT_OBJ_H_INCLUDED_
