/**********************************************************************\
  @file interpolator.h

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

// body file: interpolator.cpp

#ifndef _BRG_INTERPOLATOR_H_INCLUDED_
#define _BRG_INTERPOLATOR_H_INCLUDED_

#include <cstdlib>
#include <vector>
#include <stdexcept>

#include "brg/global.h"

#include "brg/math/interpolator/tk_spline.h"

namespace brgastro {

bool p1first_lt_p2first(std::pair<double,double> pair1, std::pair<double,double> pair2);
bool p1first_lt_v2(std::pair<double,double> pair1, double v2);

class interpolator {
public:

	enum allowed_interpolation_type {
		LOWER,
		UPPER,
		LINEAR,
		SPLINE
	}; // end enum allowed_interpolation_type

private:
#if (1)

	std::vector< std::pair<double,double> > _data_;

	mutable std::vector< std::pair<double,double> > _sorted_data_;
	mutable bool _sorted_data_cached_;

	mutable tk::spline _spline_;
	mutable bool _spline_cached_;

	static allowed_interpolation_type _default_interpolation_type_;
	allowed_interpolation_type _interpolation_type_;

	void _set_spline_points() const;

#endif // end private members and methods
public:
#if (1)

	// Swap functions
	void swap(interpolator &other);
	friend void swap(interpolator &same, interpolator &other) {same.swap(other);}

	// Constructors
	interpolator();
	interpolator(const interpolator &other);

	// Destructor
	virtual ~interpolator()
	{
	}

	// Operator=
	interpolator & operator=(interpolator other);

	// Clearing functions
	void clear();
	void clear_points();

	// Accessors to current and default interpolation types
	static const allowed_interpolation_type default_interpolation_type()
	{return _default_interpolation_type_;}
	const allowed_interpolation_type interpolation_type() const
	{return _interpolation_type_;}

	// Set functions for the current and default interpolation types
	static void set_default_interpolation_type(const allowed_interpolation_type new_default_type);
	void set_default_interpolation_type(const allowed_interpolation_type new_default_type,
			const bool override_current);
	void set_interpolation_type(const allowed_interpolation_type new_type);

	// This version doesn't check for duplicate x values, but if one does exist, an exception will
	// eventually be thrown
	void add_point(const double x, const double y);

	// This version checks if there's a point with a duplicate x value. If so, it throws an
	// exception
	void try_add_point(const double x, const double y);

	unsigned int size() const
	{
		return sorted_data().size();
	}

	std::vector< std::pair<double,double> > & sorted_data() const;

	const double operator()(const double x) const;

#endif // public

}; // end class interpolator

} // end namespace brgastro



#endif //__BRG_INTERPOLATOR_H_INCLUDED__
