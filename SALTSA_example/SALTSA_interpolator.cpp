/**       @file SALTSA_interpolator.cpp
 *
 *     Project: SALTSA_example
 *        Path: /SALTSA_example/SALTSA_interpolator.cpp
 *
 *  Created on: 16 Jul 2014
 *      Author: brg
 */

#include "SALTSA_interpolator.h"

// Global function implementations
bool SALTSA::p1first_lt_p2first(std::pair<double,double> pair1, std::pair<double,double> pair2)
{
	if(pair1.first == pair2.first)
		throw std::runtime_error("ERROR: Two points passed to interpolator have same domain value.\n");
	return (pair1.first < pair2.first);
}

// Global function implementations
bool SALTSA::p1first_lt_v2(std::pair<double,double> pair1, double v2)
{
	return (pair1.first < v2);
}

// Implement interpolator static variables
SALTSA::interpolator::allowed_interpolation_type SALTSA::interpolator::_default_interpolation_type_ = SPLINE;

// Implement interpolator methods

// Set functions for the current and default interpolation types
void SALTSA::interpolator::set_default_interpolation_type(
		const allowed_interpolation_type new_default_type)
{
	_default_interpolation_type_ = new_default_type;
}
void SALTSA::interpolator::set_default_interpolation_type(
		const allowed_interpolation_type new_default_type,
		const bool override_current)
{
	_default_interpolation_type_ = new_default_type;
	set_interpolation_type(_default_interpolation_type_);
}
void SALTSA::interpolator::set_interpolation_type(
		const allowed_interpolation_type new_type)
{
	_interpolation_type_ = new_type;
}

void SALTSA::interpolator::_set_spline_points() const
{
	std::vector<double> x_points(0), y_points(0);
	sorted_data(); // Ensure it's cached
	for(unsigned int i = 0; i < _sorted_data_.size(); i++)
	{
		x_points.push_back(_sorted_data_[i].first);
		y_points.push_back(_sorted_data_[i].second);
	}

	_spline_.set_points(x_points,y_points);
	_spline_cached_ = true;
}

SALTSA::interpolator::interpolator()
{
	_interpolation_type_ = _default_interpolation_type_;
	_spline_cached_ = false;
	_sorted_data_cached_ = false;
}

void SALTSA::interpolator::clear()
{
	clear_points();
	_interpolation_type_ = _default_interpolation_type_;
}

void SALTSA::interpolator::clear_points()
{
	_data_.clear();
	_sorted_data_.clear();
	_spline_cached_ = false;
	_sorted_data_cached_ = false;
}

void SALTSA::interpolator::add_point(const double x, const double y)
{
	_data_.push_back(std::make_pair(x,y));
	_spline_cached_ = false;
	_sorted_data_cached_ = false;
}

void SALTSA::interpolator::try_add_point(const double x, const double y)
{
	for (unsigned int i=0;i<_data_.size();i++)
	{
		if(x==_data_[i].first)
		{
			throw std::runtime_error("ERROR: Attempt to add point to interpolator with duplicate domain value.");
		}
	}
	_data_.push_back(std::make_pair(x,y));
	_spline_cached_ = false;
	_sorted_data_cached_ = false;
}

std::vector< std::pair<double,double> > & SALTSA::interpolator::sorted_data() const
{
	if(_sorted_data_cached_)
		return _sorted_data_;

	_sorted_data_ = _data_;
	std::sort(_sorted_data_.begin(),_sorted_data_.end(),SALTSA::p1first_lt_p2first);

	_sorted_data_cached_ = true;
	return _sorted_data_;
}

const double SALTSA::interpolator::operator()(const double x) const
{
	if(_interpolation_type_==SPLINE)
	{
		if(_data_.size() < 2)
			throw std::runtime_error("ERROR: Interpolator called before at least 2 points were loaded.\n");

		if(!_spline_cached_)
			_set_spline_points();

		return _spline_(x);
	}
	else if(_interpolation_type_==LINEAR)
	{
		if(_data_.size() < 2)
			throw std::runtime_error("ERROR: Interpolator called before at least 2 points were loaded.\n");

		double xlo, xhi, ylo, yhi;
		if(x<=sorted_data().front().first)
		{
			xlo = _sorted_data_[0].first;
			ylo = _sorted_data_[0].second;
			xhi = _sorted_data_[1].first;
			yhi = _sorted_data_[1].second;
		}
		else if(x>=_sorted_data_.back().first)
		{
			xlo = _sorted_data_[_sorted_data_.size()-2].first;
			ylo = _sorted_data_[_sorted_data_.size()-2].second;
			xhi = _sorted_data_.back().first;
			yhi = _sorted_data_.back().second;
		}
		else
		{
			std::vector< std::pair<double,double> >::iterator it=_sorted_data_.end();

			it = std::lower_bound(_sorted_data_.begin(),_sorted_data_.end(),
					x,p1first_lt_v2);

			if(it==_sorted_data_.end())
				throw std::runtime_error("ERROR: Could not find x value in interpolator. Check it's a valid number.");

			xlo = (it-1)->first;
			ylo = (it-1)->second;
			xhi = it->first;
			yhi = it->second;
		}

		return ylo + (yhi-ylo)/(xhi-xlo);
	}
	else if(_interpolation_type_==LOWER)
	{
		if(_data_.size() < 1)
			throw std::runtime_error("ERROR: Interpolator called before at least 1 point was loaded.\n");

		if(x<=sorted_data().front().first)
		{
			return _sorted_data_.front().second;
		}
		else if(x>=_sorted_data_.back().first)
		{
			return _sorted_data_.back().second;
		}
		else
		{
			std::vector< std::pair<double,double> >::iterator it=_sorted_data_.end();

			it = std::lower_bound(_sorted_data_.begin(),_sorted_data_.end(),
					x,p1first_lt_v2);

			if(it==_sorted_data_.end())
				throw std::runtime_error("ERROR: Could not find x value in interpolator. Check it's a valid number.");

			return (it-1)->second;
		}
	}
	else if(_interpolation_type_==UPPER)
	{
		if(_data_.size() < 1)
			throw std::runtime_error("ERROR: Interpolator called before at least 1 point was loaded.\n");

		if(x<=sorted_data().front().first)
		{
			return _sorted_data_.front().second;
		}
		else if(x>=_sorted_data_.back().first)
		{
			return _sorted_data_.back().second;
		}
		else
		{
			std::vector< std::pair<double,double> >::iterator it=_sorted_data_.end();

			it = std::lower_bound(_sorted_data_.begin(),_sorted_data_.end(),
					x,p1first_lt_v2);

			if(it==_sorted_data_.end())
				throw std::runtime_error("ERROR: Could not find x value in interpolator. Check it's a valid number.");

			return it->second;
		}
	}
	// Should never get here
	throw std::runtime_error("ERROR: Bad path reached in interpolator, which may be the result of code changes\nto the allow interpolator types.\n");
	return -1; // Just to keep editor from giving a warning

}
