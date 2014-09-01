/**********************************************************************\
  @file unit_conversions.hpp

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

#ifndef _BRG_UNIT_CONVERSIONS_HPP_INCLUDED_
#define _BRG_UNIT_CONVERSIONS_HPP_INCLUDED_

namespace brgastro{ namespace unitconv {
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

} } // namespace brgastro::unitconv



#endif /* _BRG_UNIT_CONVERSIONS_HPP_INCLUDED_ */
