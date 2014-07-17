/**********************************************************************\
brg_units.h
 -----------

 Header file including two components:


 unitconv namespace: Namespace containing various unit conversions.


 brgastro::unit_obj declarations: A class for a double with units. This
 section is only included if the header directive _BRG_USE_UNITS_ is defined.
 This can be changed in the file brg_global.h. If this is the case, the
 file brg_units.cpp must be included and compiled with the project. This
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
 x.set_value(4, unitconv::kmtom/pow(unitconv::hrtos,2)) will assign
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

#ifndef __SALTSA_UNITCONVS_H_INCLUDED__
#define __SALTSA_UNITCONVS_H_INCLUDED__

#include "SALTSA_global.h"

namespace unitconv
{
// All unit conversions are exact unless noted

// Distance
// Default unit: meter (m)
const double mtom = 1;
const double mtomm = 1e3;
const double mmtom = 1 / mtomm;
const double mtocm = 1e2;
const double cmtom = 1 / mtocm;
const double mtoum = 1e6;
const double umtom = 1 / mtoum;
const double mtonm = 1e9;
const double nmtom = 1 / mtonm;
const double mtoangstrom = 1e10;
const double angstromtom = 1 / mtoangstrom;
const double mtokm = 1e-3;
const double kmtom = 1 / mtokm;
const double ltyrtom = 9460730472580800;
const double mtoltyr = 1 / ltyrtom;
const double AUtom = 149597870700;
const double mtoAU = 1 / AUtom;
const double pctom = AUtom * 648000 / pi;
const double mtopc = 1 / pctom;
const double kpctom = 1000 * pctom;
const double mtokpc = 1 / kpctom;
const double Mpctom = 1000000 * pctom;
const double mtoMpc = 1 / Mpctom;
const double mitom = 1609.344;
const double mtomi = 1 / mitom;
const double Mmitom = 1e6 * mitom;
const double mtoMmi = 1 / Mmitom;
const double fttom = 0.3048;
const double mtoft = 1 / fttom;
const double intom = .0254;
const double mtoin = 1 / intom;
const double ydtom = 0.9144;
const double mtoyd = 1 / ydtom;

// Time
// Default unit: second (s)
const double stos = 1;
const double stocs = 1e2;
const double cstos = 1 / stocs;
const double stoms = 1e3;
const double mstos = 1 / stoms;
const double stous = 1e6;
const double ustos = 1 / stous;
const double stons = 1e9;
const double nstos = 1 / stons;
const double mintos = 60;
const double stomin = 1 / mintos;
const double hrtos = mintos * 60;
const double stohr = 1 / hrtos;
const double daytos = hrtos * 24;
const double stoday = 1 / daytos;
const double weektos = daytos * 7;
const double stoweek = 1 / weektos;
const double yrtos = daytos * 365.24219; // Approximate tropical year
const double stoyr = 1 / yrtos; // Approximate
const double monthtos = yrtos / 12; // Mean month length for tropical year
const double stomonth = 1 / monthtos; // Approximate
const double kyrtos = yrtos * 1e3; // Approximate
const double stokyr = 1 / kyrtos; // Approximate
const double Myrtos = yrtos * 1e6; // Approximate
const double stoMyr = 1 / Myrtos; // Approximate
const double Gyrtos = yrtos * 1e9; // Approximate
const double stoGyr = 1 / Gyrtos; // Approximate

// Velocity
// Default units: meters per second (mps)
const double mpstomps = 1;
const double mpstokmps = 1e-3;
const double kmpstomps = 1 / mpstokmps;
const double ctomps = 299792458;
const double mpstoc = 1 / ctomps;
const double mpstomiphr = mtomi / stohr;
const double miphr = 1 / mpstomiphr;

// Mass
// Default unit: kilogram (kg)
const double kgtokg = 1;
const double kgtogm = 1e3;
const double gmtokg = 1 / kgtogm;
const double Mearthtokg = 5.9736e24; // Approximate
const double kgtoMearth = 1 / Mearthtokg; // Approximate
const double Msuntokg = 1.9891e30; // Approximate
const double kgtoMsun = 1 / Msuntokg; // Approximate
const double kgtottMsun = kgtoMsun * 1e-10; // Approximate
const double ttMsuntokg = 1 / kgtottMsun; // Approximate

// Temperature
// Default unit: Kelvin (K)
const double KtoK = 1;
const double KtodegF = 1.8;
const double degCtoK = KtoK;
const double KtodegC = 1 / degCtoK;
const double degFtoK = 1 / KtodegF;
const double degCtodegF = KtodegF;
const double degFtodegC = degFtoK;
const double KtodegR = KtodegF;
const double degRtoK = degFtoK;
const double degCtodegR = KtodegF;
const double degRtodegC = degFtoK;

// Angle
// Default unit: radian (rad)
const double radtorad = 1;
const double degtorad = pi / 180;
const double radtodeg = 1 / degtorad;
const double amintodeg = 60;
const double degtoamin = 1 / amintodeg;
const double amintoasec = 60;
const double asectoamin = 1 / amintoasec;
const double asectodeg = asectoamin * amintodeg;
const double degtoasec = 1 / asectodeg;
const double amintorad = amintodeg * degtorad;
const double radtoamin = 1 / amintorad;
const double asectorad = asectodeg * degtorad;
const double radtoasec = 1 / asectorad;

// Charge
// Default unit: Coulomb (C)
const double CtoC = 1;
const double Ctoesu = 6.241509324e18; // Approximate
const double esutoC = 1 / Ctoesu;

}

#endif // __SALTSA_UNITCONVS_H_INCLUDED__
