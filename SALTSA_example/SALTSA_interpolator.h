/**       @file SALTSA_interpolator.h
 *
 *     Project: SALTSA_example
 *        Path: /SALTSA_example/SALTSA_interpolator.h
 *
 *  Created on: 16 Jul 2014
 *      Author: brg
 */

#ifndef __SALTSA_INTERPOLATOR_H_INCLUDED__
#define __SALTSA_INTERPOLATOR_H_INCLUDED__

#include "tk_spline.h"

namespace SALTSA {

bool p1first_lt_p2first(std::pair<double,double> pair1, std::pair<double,double> pair2);
bool p1first_lt_v2(std::pair<double,double> pair1, double v2);

	class interpolator {
	private:
#if (1)

		std::vector< std::pair<double,double> > _data_;

		mutable std::vector< std::pair<double,double> > _sorted_data_;
		mutable bool _sorted_data_cached_;

		mutable tk::spline _spline_;
		mutable bool _spline_cached_;

		void _set_spline_points() const;

#endif // end private members and methods
	public:
		enum allowed_interpolation_type {
			LOWER,
			UPPER,
			LINEAR,
			SPLINE
		}; // end enum allowed_interpolation_type

		static allowed_interpolation_type default_interpolation_type;
		allowed_interpolation_type interpolation_type;

		interpolator();

		void clear();
		void clear_points();

		void add_point(const double x, const double y);

		unsigned int size() const
		{
			return sorted_data().size();
		}

		std::vector< std::pair<double,double> > & sorted_data() const;

		const double operator()(const double x) const;

	}; // end class interpolator

} // end namespace SALTSA



#endif //__SALTSA_INTERPOLATOR_H_INCLUDED__
