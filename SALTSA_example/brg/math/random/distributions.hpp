/**********************************************************************\
  @file distributions.hpp

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

#ifndef _BRG_DISTRIBUTIONS_HPP_INCLUDED_
#define _BRG_DISTRIBUTIONS_HPP_INCLUDED_

#include <cstdlib>

#include "brg/global.h"

#include "brg/math/misc_math.hpp"

namespace brgastro {

// Error function - Thanks to John D. Cook for this code
template< typename Tx >
inline Tx erf(Tx x)
{
	using std::fabs;
	using std::exp;

    // constants
    const double a1 =  0.254829592;
    const double a2 = -0.284496736;
    const double a3 =  1.421413741;
    const double a4 = -1.453152027;
    const double a5 =  1.061405429;
    const double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    Tx y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-square(x));

    return y*=sign;
}

// Gaussian PDF
template< typename Tx >
inline const double Gaus_pdf( const Tx x )
{
	return Gaus_pdf(x,0.,1.);
}
template< typename Tx, typename Tmean >
inline const double Gaus_pdf( const Tx x, const Tmean mean )
{
	return Gaus_pdf(x,mean,1.);
}
template< typename Tx, typename Tmean, typename Tstddev >
inline const double Gaus_pdf( const Tx x, const Tmean mean,
		const Tstddev std_dev )
{
	return std::exp( -square( x - mean ) / ( 2 * std_dev * std_dev ) )
			/ ( std_dev * std::sqrt( 2 * pi ) );
}

// Function to integrate a Gaussian from min to max
template< typename Tlo, typename Thi, typename Tmean, typename Tstddev >
inline const double Gaus_int( const Tlo min, const Thi max)
{
	return Gaus_int(min,max,0.,1.);
}
template< typename Tlo, typename Thi, typename Tmean>
inline const double Gaus_int( const Tlo min, const Thi max,
		const Tmean mean)
{
	return Gaus_int(min,max,mean,1.);
}
template< typename Tlo, typename Thi, typename Tmean, typename Tstddev >
inline const double Gaus_int( const Tlo min, const Thi max,
		const Tmean mean, const Tstddev std_dev )
{
	double klo = ( min - mean ) / std_dev;
	double khi = ( max - mean ) / std_dev;

	return std::fabs( erf( khi ) - ( klo ) ) / 2;
}

} // namespace brgastro

#endif /* _BRG_DISTRIBUTIONS_HPP_INCLUDED_ */
