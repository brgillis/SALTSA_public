/**       @file brg_interpolator.h
 *
 *     Project: brg_example
 *        Path: /brg_example/brg_interpolator.h
 *
 *  Created on: 16 Jul 2014
 *      Author: brg
 */

#ifndef __BRG_INTERPOLATOR_H_INCLUDED__
#define __BRG_INTERPOLATOR_H_INCLUDED__

#include <cstdlib>
#include <vector>
#include <stdexcept>
#include "tk_spline.h"

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
