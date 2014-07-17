/**********************************************************************
 brg_orbit.cpp
 -------------

 This file contains implementations of the classes and functions declared
 in the brg_orbit.h file. For the casual user, it is not necessary to
 understand everything going here, and so I have not put as much effort
 into documenting it.

 \**********************************************************************/

#include <cstdlib>
#include <cmath>
#include <vector>
#include "SALTSA.h"
#include "SALTSA_orbit.h"
#include "SALTSA_solvers.hpp"
#include "SALTSA_calculus.hpp"

using namespace std;

/** Static member initialisations **/
#if (1)

// spline_derivative static member initialisations
#if (1)
double SALTSA::spline_derivative::_default_sample_scale_ = 0.02; // As a fraction of t_max-t_min
double SALTSA::spline_derivative::_default_sample_max_width_ = 0.05; // As a fraction of t_max-t_min
double SALTSA::spline_derivative::_default_sample_precision_ = 0.01;
#endif

// stripping_orbit static member initialisations
#if (1)
// Default integration parameters
#if(1)
// Default number of steps for which stripping is calculated
int SALTSA::stripping_orbit::_default_spline_resolution_ = 100;

// Variable step length tweaking: Time step length is proportional to (v_0/v)^(step_length_power)
// This gives smaller steps when the satellite is moving faster.
// If you want to turn off adaptive step size, set step_length_power to 0
// Alternatively, set step_length_power to 1 for even steps in position
double SALTSA::stripping_orbit::_default_v_0_ = 400 * unitconv::kmpstomps; // 400 km/s
double SALTSA::stripping_orbit::_default_r_0_ = 400 * unitconv::kpctom; // 400 kpc
double SALTSA::stripping_orbit::_default_step_length_power_ = 1.5;
double SALTSA::stripping_orbit::_default_step_factor_max_ = 10; // Maximum allowed value of (v_0/v)^(step_length_power)
double SALTSA::stripping_orbit::_default_step_factor_min_ = 0.01; // Minimum allowed value of (v_0/v)^(step_length_power)
#endif

// Default tuning values

#if(1)
// Tuning parameters, for how strong stripping and shocking are and when shocking is active
double SALTSA::stripping_orbit::_default_tidal_stripping_amplification_ = 0.675; // Tuned
double SALTSA::stripping_orbit::_default_tidal_stripping_deceleration_ = 0.1; // Tuned
double SALTSA::stripping_orbit::_default_tidal_shocking_amplification_ = 3.0; // Tuned
double SALTSA::stripping_orbit::_default_tidal_shocking_persistance_ = 1.0; // How long shocking is active for
double SALTSA::stripping_orbit::_default_tidal_shocking_power_ = -1.5; // Affects interplay of stripping and satellite halo profile
#endif

#endif

#endif

/** Class method implementations **/
#if (1)
// SALTSA::spline_function class method implementations
#if (1)
// Constructors
SALTSA::spline_function::spline_function()
{
	_spline_ptr_ = 0;
	_spline_ptr_set_up_ = false;
}
SALTSA::spline_function::spline_function(
		SALTSA::interpolator *init_spline_ptr )
{
	set_spline_ptr( init_spline_ptr );
}

// Set functions
const int SALTSA::spline_function::set_spline_ptr(
		SALTSA::interpolator *new_spline_ptr )
{
	_spline_ptr_ = new_spline_ptr;
	_spline_ptr_set_up_ = true;
	return 0;
}

// Function method
const int SALTSA::spline_function::operator()( const double & in_param,
double & out_param, const bool silent ) const
{
	if ( !_spline_ptr_set_up_ )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Spline pointer must be defined in spline_function.\n";
		return NOT_SET_UP_ERROR;
	}

	out_param = ( *_spline_ptr_ )( in_param );
	return 0;
}
#endif

// SALTSA::spline_derivative_function class method implementations
#if (1)

// Constructors
SALTSA::spline_derivative_function::spline_derivative_function()
{
	_spline_function_set_up_ = false;
}
SALTSA::spline_derivative_function::spline_derivative_function(
		SALTSA::interpolator *init_spline_ptr )
{
	set_spline_ptr( init_spline_ptr );
}

// Set functions
const int SALTSA::spline_derivative_function::set_spline_ptr(
		SALTSA::interpolator *new_spline_ptr )
{
	_spline_function_.set_spline_ptr( new_spline_ptr );
	_spline_function_set_up_ = true;
	return 0;
}

// Function method
const int SALTSA::spline_derivative_function::operator()(
		const double & in_param,
		double & out_param, const bool silent ) const
{
	if ( !_spline_function_set_up_ )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Spline function must be set up in spline_derivative_function.\n";
		return NOT_SET_UP_ERROR;
	}
	double temp_out_params;
	double Jacobian;
	unsigned int temp_num_out_params;

	if ( differentiate( &_spline_function_, 1, in_param, temp_num_out_params,
			temp_out_params, Jacobian ) )
		return 1;

	out_param = Jacobian;
	return 0;
}

#endif

// SALTSA::spline_derivative_weight_function class method implementations
#if (1)
// Constructors
SALTSA::spline_derivative_weight_function::spline_derivative_weight_function()
{
	_sample_scale_ = 0;
	_sample_max_width_ = 0;
	_t_max_ = -( 0.99 * DBL_MIN );
	_t_min_ = DBL_MAX;
	_centre_point_ = 0;
}

// Set functions
const int SALTSA::spline_derivative_weight_function::set_sample_scale(
		double new_sample_scale )
{
	_sample_scale_ = new_sample_scale;
	return 0;
}

const int SALTSA::spline_derivative_weight_function::set_sample_max_width(
		double new_sample_max_width )
{
	_sample_max_width_ = new_sample_max_width;
	return 0;
}

const int SALTSA::spline_derivative_weight_function::set_center_point(
		double new_center_point )
{
	_centre_point_ = new_center_point;
	return 0;
}

const int SALTSA::spline_derivative_weight_function::set_t_min(
		double new_t_min )
{
	_t_min_ = new_t_min;
	return 0;
}

const int SALTSA::spline_derivative_weight_function::set_t_max(
		double new_t_max )
{
	_t_max_ = new_t_max;
	return 0;
}

// Function method
const int SALTSA::spline_derivative_weight_function::operator()(
		const double & in_param,
		double & out_param, const bool silent ) const
{

	double result;

	// Check bounds
	if ( std::fabs( in_param - _centre_point_ )
			> _sample_max_width_ * std::fabs( _t_max_ - _t_min_ ) )
	{
		result = 0;
	}
	else
	{
		result = Gaus_pdf( in_param, _centre_point_,
				_sample_scale_ * std::fabs( _t_max_ - _t_min_ ) );
	}

	out_param = result;
	return 0;
}
#endif

// SALTSA::spline_derivative class method implementations
#if (1)
// Constructors
SALTSA::spline_derivative::spline_derivative()
{
	clear();
}
SALTSA::spline_derivative::spline_derivative(
		SALTSA::interpolator *init_spline_ptr )
{
	clear();
	set_spline_ptr( init_spline_ptr );
}

// Set functions
const int SALTSA::spline_derivative::set_spline_ptr(
		SALTSA::interpolator *new_spline_ptr )
{
	if ( _spline_func_.set_spline_ptr( new_spline_ptr ) )
		return 1;
	_spline_ptr_ = new_spline_ptr;
	_spline_ptr_set_up_ = true;
	return 0;
}

const int SALTSA::spline_derivative::clear_spline_ptr()
{
	_spline_func_.set_spline_ptr( 0 );
	_spline_ptr_ = 0;
	_spline_ptr_set_up_ = false;
	return 0;
}

const int SALTSA::spline_derivative::set_default_sample_scale(
		double new_default_sample_scale )
{
	_default_sample_scale_ = new_default_sample_scale;
	return 0;
}

const int SALTSA::spline_derivative::set_default_sample_max_width(
		double new_default_sample_max_width )
{
	_default_sample_max_width_ = new_default_sample_max_width;
	return 0;
}

const int SALTSA::spline_derivative::set_sample_scale(
		double new_sample_scale )
{
	_sample_scale_ = new_sample_scale;
	_calculated_ = false;
	return 0;
}

const int SALTSA::spline_derivative::set_sample_max_width(
		double new_sample_max_width )
{
	_sample_max_width_ = new_sample_max_width;
	_calculated_ = false;
	return 0;
}

const int SALTSA::spline_derivative::reset_sample_scale() // Sets it to default
{
	_sample_scale_ = _default_sample_scale_;
	_calculated_ = false;
	return 0;
}

const int SALTSA::spline_derivative::reset_sample_max_width() // Sets it to default
{
	_sample_max_width_ = _default_sample_max_width_;
	_calculated_ = false;
	return 0;
}

const int SALTSA::spline_derivative::clear_known_points()
{
	_known_spline_.clear();
	_calculated_ = false;
	return 0;
}

const int SALTSA::spline_derivative::clear_unknown_points()
{
	_unknown_t_list_.clear();
	_calculated_ = false;
	return 0;
}

const int SALTSA::spline_derivative::clear_points()
{
	_known_spline_.clear();
	_unknown_t_list_.clear();
	_calculated_ = false;
	return 0;
}

const int SALTSA::spline_derivative::clear()
{
	_spline_ptr_ = 0;
	_spline_ptr_set_up_ = false;
	_known_spline_.clear();
	_unknown_t_list_.clear();
	_t_max_ = ( -DBL_MAX );
	_t_min_ = DBL_MAX;
	_calculated_ = false;
	_sample_scale_ = _default_sample_scale_;
	_sample_max_width_ = _default_sample_max_width_;
	_sample_precision_ = _default_sample_precision_;
	return 0;
}

const int SALTSA::spline_derivative::add_point( const double t,
		const double x )
{
	_known_spline_.add_point( t, x );
	_calculated_ = false;
	return 0;
}

const int SALTSA::spline_derivative::add_unknown_point( const double t )
{
	_unknown_t_list_.push_back( t );
	_calculated_ = false;
	return 0;
}

// Get functions
const double SALTSA::spline_derivative::operator()( double xval ) const
{
	if ( !_spline_ptr_set_up_ )
	{
		if ( _known_spline_.size() >= 2 ) // We can use the known spline for everything
		{
			return _known_spline_( xval );
		} // if(known_spline.size() >= 2)
		else // We don't know enough to get any points
		{
			throw std::runtime_error("ERROR: Spline_derivative called without spline assigned to it.\n");
			return DBL_MIN;
		} //  if(known_spline.size() >= 2) ... else
	} // if(!spline_set_up)

	if ( _calculated_ )
	{
		return _estimated_spline_( xval );
	} // if(calculated)
	else // We'll have to calculate
	{
		// Get t_min and t_max
		_t_min_ = DBL_MAX;
		_t_max_ = ( -DBL_MAX );

		for ( unsigned int i = 0; i < _known_spline_.sorted_data().size(); i++ )
		{
			if ( _known_spline_.sorted_data().at( i ).first < _t_min_ )
				_t_min_ = _known_spline_.sorted_data().at( i ).first;
			if ( _known_spline_.sorted_data().at( i ).first > _t_max_ )
				_t_max_ = _known_spline_.sorted_data().at( i ).first;
		}

		for ( unsigned int i = 0; i < _unknown_t_list_.size(); i++ )
		{
			if ( _unknown_t_list_[i] < _t_min_ )
				_t_min_ = _unknown_t_list_[i];
			if ( _unknown_t_list_[i] > _t_max_ )
				_t_max_ = _unknown_t_list_[i];
		}

		// Set up the estimated spline, starting by making a copy of the known spline
		_estimated_spline_ = _known_spline_;
		unsigned int num_points_to_calculate = _unknown_t_list_.size();
		spline_derivative_function spline_derivative_function_val(
				_spline_ptr_ );
		spline_derivative_weight_function spline_derivative_weight_function_val;

		spline_derivative_weight_function_val.set_sample_scale(
				_sample_scale_ );
		spline_derivative_weight_function_val.set_sample_max_width(
				_sample_max_width_ );
		spline_derivative_weight_function_val.set_t_min( _t_min_ );
		spline_derivative_weight_function_val.set_t_max( _t_max_ );

		double delta_t = fabs( _t_max_ - _t_min_ ) * _sample_max_width_;

		for ( unsigned int i = 0; i < num_points_to_calculate; i++ ) // For each point we need to calculate
		{
			double t = _unknown_t_list_[i];
			spline_derivative_weight_function_val.set_center_point( t );
			unsigned int num_in_params = 1, num_out_params = 1;
			double min_in_params( t - delta_t );
			double max_in_params( t + delta_t );
			double out_params( 0 );
			double Jacobian( 0 );

			if ( delta_t <= 0 )
			{
				if ( differentiate( &spline_derivative_function_val,
						num_in_params, min_in_params, num_out_params,
						out_params, Jacobian ) )
					return 1;

				_estimated_spline_.add_point( t, Jacobian );
			}
			else
			{

				if ( integrate_weighted_Rhomberg(
						&spline_derivative_function_val,
						&spline_derivative_weight_function_val, num_in_params,
						min_in_params, max_in_params, num_out_params,
						out_params, _sample_precision_ ) )
					return 1;

				_estimated_spline_.add_point( t, out_params );

			}

		} // for(unsigned int i = 0; i < num_points_to_calculate; i++ )
	} // if(calculated) ... else

	_calculated_ = true;

	return _estimated_spline_( xval );
} // const double operator()(double xval)

#endif

// SALTSA::stripping_orbit class method implementations
#if (1)

SALTSA::stripping_orbit::stripping_orbit()
{
	_init();
}

const int SALTSA::stripping_orbit::_init()
{
	return clear();
}

const std::vector< SALTSA::stripping_orbit_segment >::iterator SALTSA::stripping_orbit::_final_good_segment() const
{
	if(!_calculated_)
	{
		_final_good_segment_ = _orbit_segments_.end();
		calc(); // It might not work, but running it will at least give us the most sensible result
	}
	return _final_good_segment_;
}

SALTSA::stripping_orbit::stripping_orbit(
		const stripping_orbit &other_stripping_orbit )
{
	// Integration parameters
#if(1)
	_spline_resolution_ = other_stripping_orbit._spline_resolution_;
	_v_0_ = other_stripping_orbit._v_0_;
	_r_0_ = other_stripping_orbit._r_0_;
	_step_length_power_ = other_stripping_orbit._step_length_power_;
	_step_factor_max_ = other_stripping_orbit._step_factor_max_;
	_step_factor_min_ = other_stripping_orbit._step_factor_min_;
#endif

	// Tuning values
#if(1)
	_tidal_stripping_amplification_ = other_stripping_orbit._tidal_stripping_amplification_;
	_tidal_stripping_deceleration_ = other_stripping_orbit._tidal_stripping_deceleration_;
	_tidal_shocking_amplification_ = other_stripping_orbit._tidal_shocking_amplification_;
	_tidal_shocking_persistance_ = other_stripping_orbit._tidal_shocking_persistance_;
	_tidal_shocking_power_ = other_stripping_orbit._tidal_shocking_power_;
#endif

	_num_patameters_ = other_stripping_orbit._num_patameters_;
	_satellite_parameter_unitconvs_ =
			other_stripping_orbit._satellite_parameter_unitconvs_;
	_satellite_output_parameters_ =
			other_stripping_orbit._satellite_output_parameters_;
	_host_parameter_unitconvs_ =
			other_stripping_orbit._host_parameter_unitconvs_;
	_host_output_parameters_ = other_stripping_orbit._host_output_parameters_;
	_record_full_data_ = other_stripping_orbit._record_full_data_;
	_num_segments_ = other_stripping_orbit._num_segments_;

	_t_min_natural_value_ = other_stripping_orbit._t_min_natural_value_;
	_t_max_natural_value_ = other_stripping_orbit._t_max_natural_value_;
	_t_min_override_value_ = other_stripping_orbit._t_min_override_value_;
	_t_max_override_value_ = other_stripping_orbit._t_max_override_value_;
	_override_t_min_ = other_stripping_orbit._override_t_min_;
	_override_t_max_ = other_stripping_orbit._override_t_max_;

	_final_fmret_list_ = other_stripping_orbit._final_fmret_list_;

	_x_spline_points_ = other_stripping_orbit._x_spline_points_;
	_y_spline_points_ = other_stripping_orbit._y_spline_points_;
	_z_spline_points_ = other_stripping_orbit._z_spline_points_;
	_test_mass_spline_points_ =
			other_stripping_orbit._test_mass_spline_points_;
	_vx_spline_points_ = other_stripping_orbit._vx_spline_points_;
	_vy_spline_points_ = other_stripping_orbit._vy_spline_points_;
	_vz_spline_points_ = other_stripping_orbit._vz_spline_points_;
	_vx_spline_unknown_points_ =
			other_stripping_orbit._vx_spline_unknown_points_;
	_vy_spline_unknown_points_ =
			other_stripping_orbit._vy_spline_unknown_points_;
	_vz_spline_unknown_points_ =
			other_stripping_orbit._vz_spline_unknown_points_;
	_host_parameter_spline_points_ =
			other_stripping_orbit._host_parameter_spline_points_;
	_discontinuity_times_ = other_stripping_orbit._discontinuity_times_;
	_cleaned_discontinuity_times_ =
			other_stripping_orbit._cleaned_discontinuity_times_;
	_num_discontinuities_ = other_stripping_orbit._num_discontinuities_;
	_num_cleaned_discontinuities_ =
			other_stripping_orbit._num_cleaned_discontinuities_;

	_init_host_ptr_ = other_stripping_orbit._init_host_ptr_;
	_init_satellite_ptr_ = other_stripping_orbit._init_satellite_ptr_;
	_host_is_evolving_ = other_stripping_orbit._host_is_evolving_;
	_calculated_ = other_stripping_orbit._calculated_;
	_bad_result_ = other_stripping_orbit._bad_result_;
	_host_loaded_ = other_stripping_orbit._host_loaded_;
	_satellite_loaded_ = other_stripping_orbit._satellite_loaded_;
	_using_private_init_host_ =
			other_stripping_orbit._using_private_init_host_;
	_using_private_init_satellite_ =
			other_stripping_orbit._using_private_init_satellite_;
	_private_tNFW_init_host_ = other_stripping_orbit._private_tNFW_init_host_;
	_private_tNFW_init_satellite_ =
			other_stripping_orbit._private_tNFW_init_satellite_;

	_orbit_segments_ = other_stripping_orbit._orbit_segments_;

	if( (other_stripping_orbit._final_good_segment_ == other_stripping_orbit._orbit_segments_.end()) || ( _orbit_segments_.empty() ) )
	{
		_final_good_segment_ = _orbit_segments_.end();
	}
	else
	{
		// We'll have to find the segment
		std::vector< SALTSA::stripping_orbit_segment >::iterator test_iterator_old = other_stripping_orbit._orbit_segments_.begin();
		std::vector< SALTSA::stripping_orbit_segment >::iterator test_iterator_new = _orbit_segments_.begin();
		bool found = false;
		while(test_iterator_old != other_stripping_orbit._orbit_segments_.end())
		{
			if(test_iterator_old == other_stripping_orbit._final_good_segment_)
			{
				found = true;
				_final_good_segment_ = test_iterator_new;
				break;
			}
			else
			{
				test_iterator_old++;
				test_iterator_new++;
			}
		}
		if(!found) _final_good_segment_ = _orbit_segments_.end();
	}

	if ( _using_private_init_host_ )
	{
		_init_host_ptr_ = &_private_tNFW_init_host_;
	}
	if ( _using_private_init_satellite_ )
	{
		_init_satellite_ptr_ = &_private_tNFW_init_satellite_;
	}
}

SALTSA::stripping_orbit & SALTSA::stripping_orbit::operator=(
		const stripping_orbit &other_stripping_orbit )
{
	if ( this != &other_stripping_orbit )
	{
		clear();

		// Integration parameters
#if(1)
		_spline_resolution_ = other_stripping_orbit._spline_resolution_;
		_v_0_ = other_stripping_orbit._v_0_;
		_r_0_ = other_stripping_orbit._r_0_;
		_step_length_power_ = other_stripping_orbit._step_length_power_;
		_step_factor_max_ = other_stripping_orbit._step_factor_max_;
		_step_factor_min_ = other_stripping_orbit._step_factor_min_;
#endif

		// Tuning values
#if(1)
		_tidal_stripping_amplification_ = other_stripping_orbit._tidal_stripping_amplification_;
		_tidal_stripping_deceleration_ = other_stripping_orbit._tidal_stripping_deceleration_;
		_tidal_shocking_amplification_ = other_stripping_orbit._tidal_shocking_amplification_;
		_tidal_shocking_persistance_ = other_stripping_orbit._tidal_shocking_persistance_;
		_tidal_shocking_power_ = other_stripping_orbit._tidal_shocking_power_;
#endif

		_num_patameters_ = other_stripping_orbit._num_patameters_;
		_satellite_parameter_unitconvs_ =
				other_stripping_orbit._satellite_parameter_unitconvs_;
		_satellite_output_parameters_ =
				other_stripping_orbit._satellite_output_parameters_;
		_host_parameter_unitconvs_ =
				other_stripping_orbit._host_parameter_unitconvs_;
		_host_output_parameters_ =
				other_stripping_orbit._host_output_parameters_;
		_record_full_data_ = other_stripping_orbit._record_full_data_;
		_num_segments_ = other_stripping_orbit._num_segments_;

		_t_min_natural_value_ = other_stripping_orbit._t_min_natural_value_;
		_t_max_natural_value_ = other_stripping_orbit._t_max_natural_value_;
		_t_min_override_value_ = other_stripping_orbit._t_min_override_value_;
		_t_max_override_value_ = other_stripping_orbit._t_max_override_value_;
		_override_t_min_ = other_stripping_orbit._override_t_min_;
		_override_t_max_ = other_stripping_orbit._override_t_max_;

		_final_fmret_list_ = other_stripping_orbit._final_fmret_list_;

		_x_spline_points_ = other_stripping_orbit._x_spline_points_;
		_y_spline_points_ = other_stripping_orbit._y_spline_points_;
		_z_spline_points_ = other_stripping_orbit._z_spline_points_;
		_test_mass_spline_points_ =
				other_stripping_orbit._test_mass_spline_points_;
		_vx_spline_points_ = other_stripping_orbit._vx_spline_points_;
		_vy_spline_points_ = other_stripping_orbit._vy_spline_points_;
		_vz_spline_points_ = other_stripping_orbit._vz_spline_points_;
		_vx_spline_unknown_points_ =
				other_stripping_orbit._vx_spline_unknown_points_;
		_vy_spline_unknown_points_ =
				other_stripping_orbit._vy_spline_unknown_points_;
		_vz_spline_unknown_points_ =
				other_stripping_orbit._vz_spline_unknown_points_;
		_host_parameter_spline_points_ =
				other_stripping_orbit._host_parameter_spline_points_;
		_discontinuity_times_ = other_stripping_orbit._discontinuity_times_;
		_cleaned_discontinuity_times_ =
				other_stripping_orbit._cleaned_discontinuity_times_;
		_num_discontinuities_ = other_stripping_orbit._num_discontinuities_;
		_num_cleaned_discontinuities_ =
				other_stripping_orbit._num_cleaned_discontinuities_;

		_init_host_ptr_ = other_stripping_orbit._init_host_ptr_;
		_init_satellite_ptr_ = other_stripping_orbit._init_satellite_ptr_;
		_host_is_evolving_ = other_stripping_orbit._host_is_evolving_;
		_calculated_ = other_stripping_orbit._calculated_;
		_bad_result_ = other_stripping_orbit._bad_result_;
		_host_loaded_ = other_stripping_orbit._host_loaded_;
		_satellite_loaded_ = other_stripping_orbit._satellite_loaded_;
		_using_private_init_host_ =
				other_stripping_orbit._using_private_init_host_;
		_using_private_init_satellite_ =
				other_stripping_orbit._using_private_init_satellite_;
		_private_tNFW_init_host_ =
				other_stripping_orbit._private_tNFW_init_host_;
		_private_tNFW_init_satellite_ =
				other_stripping_orbit._private_tNFW_init_satellite_;

		_orbit_segments_ = other_stripping_orbit._orbit_segments_;

		if( (other_stripping_orbit._final_good_segment_ == other_stripping_orbit._orbit_segments_.end()) || ( _orbit_segments_.empty() ) )
		{
			_final_good_segment_ = _orbit_segments_.end();
		}
		else
		{
			// We'll have to find the segment
			std::vector< SALTSA::stripping_orbit_segment >::iterator test_iterator_old = other_stripping_orbit._orbit_segments_.begin();
			std::vector< SALTSA::stripping_orbit_segment >::iterator test_iterator_new = _orbit_segments_.begin();
			bool found = false;
			while(test_iterator_old != other_stripping_orbit._orbit_segments_.end())
			{
				if(test_iterator_old == other_stripping_orbit._final_good_segment_)
				{
					found = true;
					_final_good_segment_ = test_iterator_new;
					break;
				}
				else
				{
					test_iterator_old++;
					test_iterator_new++;
				}
			}
			if(!found) _final_good_segment_ = _orbit_segments_.end();
		}

		if ( _using_private_init_host_ )
		{
			_init_host_ptr_ = &_private_tNFW_init_host_;
		}
		if ( _using_private_init_satellite_ )
		{
			_init_satellite_ptr_ = &_private_tNFW_init_satellite_;
		}
	}
	return *this;
}

SALTSA::stripping_orbit::stripping_orbit( density_profile *init_init_host,
		density_profile *init_init_satellite, const int init_resolution )
{
	clear();
	_init_host_ptr_ = init_init_host;
	_init_satellite_ptr_ = init_init_satellite;
	_spline_resolution_ = init_resolution;
	_host_loaded_ = true;
	_satellite_loaded_ = true;

}

SALTSA::stripping_orbit::~stripping_orbit()
{
}

SALTSA::stripping_orbit *SALTSA::stripping_orbit::stripping_orbit_clone()
{
	return new stripping_orbit( *this );
}

const int SALTSA::stripping_orbit::clear()
{
	clear_calcs();

	// Integration parameters
#if(1)
	_spline_resolution_ = _default_spline_resolution_;
	_v_0_ = _default_v_0_;
	_r_0_ = _default_r_0_;
	_step_length_power_ = _default_step_length_power_;
	_step_factor_max_ = _default_step_factor_max_;
	_step_factor_min_ = _default_step_factor_min_;
#endif

	// Tuning values
#if(1)
	_tidal_stripping_amplification_ = _default_tidal_stripping_amplification_;
	_tidal_stripping_deceleration_ = _default_tidal_stripping_deceleration_;
	_tidal_shocking_amplification_ = _default_tidal_shocking_amplification_;
	_tidal_shocking_persistance_ = _default_tidal_shocking_persistance_;
	_tidal_shocking_power_ = _default_tidal_shocking_power_;
#endif

	_satellite_parameter_unitconvs_.clear();
	_satellite_output_parameters_.clear();
	_host_parameter_unitconvs_.clear();
	_host_output_parameters_.clear();

	_x_spline_points_.clear();
	_y_spline_points_.clear();
	_z_spline_points_.clear();
	_vx_spline_points_.clear();
	_vy_spline_points_.clear();
	_vz_spline_points_.clear();
	_vx_spline_unknown_points_.clear();
	_vy_spline_unknown_points_.clear();
	_vz_spline_unknown_points_.clear();
	_test_mass_spline_points_.clear();
	_host_parameter_spline_points_.clear();

	_t_min_natural_value_ = DBL_MAX;
	_t_max_natural_value_ = ( -DBL_MAX );
	_t_min_override_value_ = DBL_MAX;
	_t_max_override_value_ = ( -DBL_MAX );
	_override_t_min_ = false;
	_override_t_max_ = false;
	_num_patameters_ = 0;
	_record_full_data_ = false;
	_num_discontinuities_ = 0;
	_num_cleaned_discontinuities_ = 0;

	_init_host_ptr_ = 0;
	_init_satellite_ptr_ = 0;

	_host_loaded_ = false;
	_satellite_loaded_ = false;
	_using_private_init_host_ = false;
	_using_private_init_satellite_ = false;
	_host_is_evolving_ = false;

	return 0;
}

const int SALTSA::stripping_orbit::clear_calcs() const
{
	_orbit_segments_.clear();
	_num_segments_ = 0;
	_final_good_segment_ = _orbit_segments_.end();
	_final_fmret_list_.clear();

	_calculated_ = false;
	_bad_result_ = false;
	return 0;
}

const int SALTSA::stripping_orbit::add_point( const double &x,
		const double &y, const double &z, const double &t,
		const double new_mass )
{
	_calculated_ = false;
	try
	{
		_x_spline_points_.push_back( std::pair< double, double >( t, x ) );
		_y_spline_points_.push_back( std::pair< double, double >( t, y ) );
		_z_spline_points_.push_back( std::pair< double, double >( t, z ) );
		_vx_spline_unknown_points_.push_back( t );
		_vy_spline_unknown_points_.push_back( t );
		_vz_spline_unknown_points_.push_back( t );
		_test_mass_spline_points_.push_back(
				std::pair< double, double >( t, new_mass ) );
		if ( t < _t_min_natural_value_ )
			_t_min_natural_value_ = t;
		if ( t > _t_max_natural_value_ )
			_t_max_natural_value_ = t;
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit::add_point().\n"
				<< "Exception: " << e.what();
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	catch (...)
	{
		std::cerr << "ERROR: Exception in stripping_orbit::add_point().\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

const int SALTSA::stripping_orbit::add_point( const double &x,
		const double &y, const double &z, const double &vx,
		const double &vy, const double &vz, const double &t,
		const double new_test_mass )
{
	_calculated_ = false;
	try
	{
		_x_spline_points_.push_back( std::pair< double, double >( t, x ) );
		_y_spline_points_.push_back( std::pair< double, double >( t, y ) );
		_z_spline_points_.push_back( std::pair< double, double >( t, z ) );
		_vx_spline_points_.push_back( std::pair< double, double >( t, vx ) );
		_vy_spline_points_.push_back( std::pair< double, double >( t, vy ) );
		_vz_spline_points_.push_back( std::pair< double, double >( t, vz ) );
		_test_mass_spline_points_.push_back(
				std::pair< double, double >( t, new_test_mass ) );
		if ( t < _t_min_natural_value_ )
			_t_min_natural_value_ = t;
		if ( t > _t_max_natural_value_ )
			_t_max_natural_value_ = t;
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit::add_point().\n"
				<< "Exception: " << e.what();
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	catch (...)
	{
		std::cerr << "ERROR: Exception in stripping_orbit::add_point().\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

const int SALTSA::stripping_orbit::add_host_parameter_point(
		const unsigned int num_parameters,
		const std::vector< double > &parameters, const double &t,
		const bool silent )
{
	// Check num_parameters matches vector size
	if ( num_parameters != parameters.size() )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: num_parameters must == parameters.size() in stripping_orbit::add_host_parameter_point.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	if ( num_parameters <= 0 )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: num_parameters must be > 0 in stripping_orbit::add_host_parameter_point.\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	_host_parameter_spline_points_.push_back(
			std::pair< double, std::vector< double > >( t, parameters ) );

	if ( _host_parameter_spline_points_.size() >= 2 )
		_host_is_evolving_ = true;

	_calculated_ = false;

	return 0;
}

const int SALTSA::stripping_orbit::add_discontinuity_time(
		const double &t )
{
	try
	{
		_discontinuity_times_.push_back( t );
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit::add_discontinuity_time().\n"
				<< "Exception: " << e.what();
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	catch (...)
	{
		std::cerr << "ERROR: Exception in stripping_orbit::add_discontinuity_time().\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

const int SALTSA::stripping_orbit::clear_points()
{
	_x_spline_points_.clear();
	_y_spline_points_.clear();
	_z_spline_points_.clear();
	_vx_spline_points_.clear();
	_vy_spline_points_.clear();
	_vz_spline_points_.clear();
	_test_mass_spline_points_.clear();
	_calculated_ = false;
	return 0;
}
const int SALTSA::stripping_orbit::clear_discontinuity_times()
{
	_discontinuity_times_.clear();
	_calculated_ = false;
	return 0;
}
const int SALTSA::stripping_orbit::clear_host_parameter_points()
{
	_host_parameter_spline_points_.clear();
	_calculated_ = false;
	return 0;
}

const int SALTSA::stripping_orbit::set_init_host(
		const density_profile *new_init_host )
{
	_init_host_ptr_ = new_init_host;
	_host_loaded_ = true;
	_using_private_init_host_ = false;
	_calculated_ = false;
	return 0;
}

const int SALTSA::stripping_orbit::set_init_satellite(
		const density_profile *new_init_satellite )
{
	_init_satellite_ptr_ = new_init_satellite;
	_satellite_loaded_ = true;
	_using_private_init_satellite_ = false;
	_calculated_ = false;
	return 0;
}

const int SALTSA::stripping_orbit::set_tNFW_init_satellite(
		const double &new_init_mvir0, const double z,
		const double new_init_c, const double new_init_tau )
{
	_using_private_init_satellite_ = true;
	_private_tNFW_init_satellite_ = tNFW_profile( new_init_mvir0, z,
			new_init_c, new_init_tau );
	_init_satellite_ptr_ = &_private_tNFW_init_satellite_;

	if ( _using_private_init_host_ )
	{
		if ( _private_tNFW_init_host_.z() != z )
			_private_tNFW_init_host_.set_z( z );
	}

	_satellite_loaded_ = true;
	_calculated_ = false;
	return 0;
}

const int SALTSA::stripping_orbit::set_tNFW_host(
		const double &new_init_mvir0, const double z,
		const double new_init_c, const double new_init_tau )
{
	_using_private_init_host_ = true;
	_private_tNFW_init_host_ = tNFW_profile( new_init_mvir0, z, new_init_c,
			new_init_tau );
	_init_host_ptr_ = &_private_tNFW_init_host_;

	if ( _using_private_init_satellite_ )
	{
		if ( _private_tNFW_init_satellite_.z() != z )
		{
			_private_tNFW_init_satellite_.set_z( z );
		}
	}

	_host_loaded_ = true;
	_calculated_ = false;
	return 0;
}

const int SALTSA::stripping_orbit::clear_init_satellite()
{
	_init_satellite_ptr_ = 0;
	_satellite_loaded_ = false;
	_using_private_init_satellite_ = false;
	_calculated_ = false;
	return 0;
}
const int SALTSA::stripping_orbit::clear_init_host()
{
	_init_host_ptr_ = 0;
	_host_loaded_ = false;
	_using_private_init_host_ = false;
	_calculated_ = false;
	return 0;
}

// Setting default integration parameters
#if(1)
const int SALTSA::stripping_orbit::set_default_resolution( const int new_default_resolution,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_resolution == _default_spline_resolution_ )
		return 0;

	if ( new_default_resolution < 2 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set default resolution to value below minimum of 2.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	_default_spline_resolution_ = new_default_resolution;

	reset_resolution();
	return 0;
}
const int SALTSA::stripping_orbit::set_default_v_0( const double new_default_v_0,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_v_0 == _default_v_0_ )
		return 0;

	if ( new_default_v_0 <= 0 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set default v_0 to value at or below minimum of 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	_default_v_0_ = new_default_v_0;

	reset_v_0();
	return 0;
}
const int SALTSA::stripping_orbit::set_default_r_0( const double new_default_r_0,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_r_0 == _default_r_0_ )
		return 0;

	if ( new_default_r_0 <= 0 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set default r_0 to value at or below minimum of 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	_default_r_0_ = new_default_r_0;

	reset_r_0();
	return 0;
}
const int SALTSA::stripping_orbit::set_default_step_length_power(
		const double new_default_step_length_power,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_step_length_power == _default_step_length_power_ )
		return 0;
	_default_step_length_power_ = new_default_step_length_power;

	reset_step_length_power();
	return 0;
}
const int SALTSA::stripping_orbit::set_default_step_factor_max(
		const double new_default_step_factor_max,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_step_factor_max == _default_step_factor_max_ )
		return 0;

	if ( new_default_step_factor_max < 1 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set default step_factor_max to value below minimum of 1.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	_default_step_factor_max_ = new_default_step_factor_max;

	reset_step_factor_max();
	return 0;
}
const int SALTSA::stripping_orbit::set_default_step_factor_min(
		const double new_default_step_factor_min,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_step_factor_min == _default_step_factor_min_ )
		return 0;

	if ( new_default_step_factor_min > 1 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set default step_factor_min to value above minimum of 1.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	_default_step_factor_min_ = new_default_step_factor_min;

	reset_step_factor_min();
	return 0;
}
#endif

// Setting default tuning parameters
#if(1)

const int SALTSA::stripping_orbit::set_default_tidal_stripping_amplification(
		const double new_default_tidal_stripping_amplification,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_tidal_stripping_amplification == _default_tidal_stripping_amplification_ )
		return 0;

	if ( new_default_tidal_stripping_amplification < 0 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set default tidal_stripping_amplification to value below minimum of 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	_default_tidal_stripping_amplification_ = new_default_tidal_stripping_amplification;

	reset_tidal_stripping_amplification();
	return 0;
}
const int SALTSA::stripping_orbit::set_default_tidal_stripping_deceleration(
		const double new_default_tidal_stripping_deceleration,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_tidal_stripping_deceleration == _default_tidal_stripping_deceleration_ )
		return 0;
	_default_tidal_stripping_deceleration_ = new_default_tidal_stripping_deceleration;

	reset_tidal_stripping_deceleration();
	return 0;
}
const int SALTSA::stripping_orbit::set_default_tidal_shocking_amplification(
		const double new_default_tidal_shocking_amplification,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_tidal_shocking_amplification == _default_tidal_shocking_amplification_ )
		return 0;

	if ( new_default_tidal_shocking_amplification < 0 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set default tidal_shocking_amplification to value below minimum of 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	_default_tidal_shocking_amplification_ = new_default_tidal_shocking_amplification;

	reset_tidal_shocking_amplification();

	return 0;
}
const int SALTSA::stripping_orbit::set_default_tidal_shocking_persistance(
		const double new_default_tidal_shocking_persistance,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_tidal_shocking_persistance == _default_tidal_shocking_persistance_ )
		return 0;

	if ( new_default_tidal_shocking_persistance <= 0 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set default tidal_stripping_persistance to value at or below minimum of 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	_default_tidal_shocking_persistance_ = new_default_tidal_shocking_persistance;

	reset_tidal_shocking_persistance();

	return 0;
}
const int SALTSA::stripping_orbit::set_default_tidal_shocking_power(
		const double new_default_tidal_shocking_power,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_tidal_shocking_power == _default_tidal_shocking_power_ )
		return 0;
	_default_tidal_shocking_power_ = new_default_tidal_shocking_power;

	reset_tidal_shocking_power();

	return 0;
}
#endif

// Setting integration parameters
#if(1)
const int SALTSA::stripping_orbit::set_resolution( const int new_resolution,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_resolution == _spline_resolution_ )
		return 0;

	if ( new_resolution < 2 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set resolution to value below minimum of 2.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	clear_calcs();
	_spline_resolution_ = new_resolution;
	return 0;
}
const int SALTSA::stripping_orbit::set_v_0( const double new_v_0,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_v_0 == _v_0_ )
		return 0;

	if ( new_v_0 <= 0 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set v_0 to value at or below minimum of 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	clear_calcs();
	_v_0_ = new_v_0;
	return 0;
}
const int SALTSA::stripping_orbit::set_r_0( const double new_r_0,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_r_0 == _r_0_ )
		return 0;

	if ( new_r_0 <= 0 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set r_0 to value at or below minimum of 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	clear_calcs();
	_r_0_ = new_r_0;
	return 0;
}
const int SALTSA::stripping_orbit::set_step_length_power( const double new_step_length_power,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_step_length_power == _step_length_power_ )
		return 0;

	clear_calcs();
	_step_length_power_ = new_step_length_power;
	return 0;
}
const int SALTSA::stripping_orbit::set_step_factor_max( const double new_step_factor_max,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_step_factor_max == _step_factor_max_ )
		return 0;

	if ( new_step_factor_max < 1 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set step_factor_max to value below minimum of 1.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	clear_calcs();
	_step_factor_max_ = new_step_factor_max;
	return 0;
}
const int SALTSA::stripping_orbit::set_step_factor_min( const double new_step_factor_min,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_step_factor_min == _step_factor_min_ )
		return 0;

	if ( new_step_factor_min > 1 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set step_factor_min to value above minimum of 1.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	clear_calcs();
	_step_factor_min_ = new_step_factor_min;
	return 0;
}
#endif

// Setting tuning parameters
#if(1)

const int SALTSA::stripping_orbit::set_tidal_stripping_amplification(
		const double new_tidal_stripping_amplification,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_stripping_amplification == _tidal_stripping_amplification_ )
		return 0;

	if ( new_tidal_stripping_amplification < 0 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set tidal_stripping_amplification to value below minimum of 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	clear_calcs();
	_tidal_stripping_amplification_ = new_tidal_stripping_amplification;
	return 0;
}
const int SALTSA::stripping_orbit::set_tidal_stripping_deceleration(
		const double new_tidal_stripping_deceleration,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_stripping_deceleration == _tidal_stripping_deceleration_ )
		return 0;

	clear_calcs();
	_tidal_stripping_deceleration_ = new_tidal_stripping_deceleration;
	return 0;
}
const int SALTSA::stripping_orbit::set_tidal_shocking_amplification(
		const double new_tidal_shocking_amplification,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_shocking_amplification == _tidal_shocking_amplification_ )
		return 0;

	if ( new_tidal_shocking_amplification < 0 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set tidal_shocking_amplification to value below minimum of 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	clear_calcs();
	_tidal_shocking_amplification_ = new_tidal_shocking_amplification;
	return 0;
}
const int SALTSA::stripping_orbit::set_tidal_shocking_persistance(
		const double new_tidal_shocking_persistance,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_shocking_persistance == _tidal_shocking_persistance_ )
		return 0;

	if ( new_tidal_shocking_persistance <= 0 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set tidal_shocking_persistance to value at or below minimum of 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	clear_calcs();
	_tidal_shocking_persistance_ = new_tidal_shocking_persistance;
	return 0;
}
const int SALTSA::stripping_orbit::set_tidal_shocking_power(
		const double new_tidal_shocking_power,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_shocking_power == _tidal_shocking_power_ )
		return 0;

	clear_calcs();
	_tidal_shocking_power_ = new_tidal_shocking_power;
	return 0;
}
#endif

// Setting integration parameters
#if(1)
const int SALTSA::stripping_orbit::reset_resolution()
{
	// Check if anything is actually changing here
	if ( _spline_resolution_ == _default_spline_resolution_ )
		return 0;

	clear_calcs();
	_spline_resolution_ = _default_spline_resolution_;
	return 0;
}
const int SALTSA::stripping_orbit::reset_v_0()
{
	// Check if anything is actually changing here
	if ( _default_v_0_ == _v_0_ )
		return 0;

	clear_calcs();
	_v_0_ = _default_v_0_;
	return 0;
}
const int SALTSA::stripping_orbit::reset_r_0()
{
	// Check if anything is actually changing here
	if ( _default_r_0_ == _r_0_ )
		return 0;

	clear_calcs();
	_r_0_ = _default_r_0_;
	return 0;
}
const int SALTSA::stripping_orbit::reset_step_length_power()
{
	// Check if anything is actually changing here
	if ( _default_step_length_power_ == _step_length_power_ )
		return 0;

	clear_calcs();
	_step_length_power_ = _default_step_length_power_;
	return 0;
}
const int SALTSA::stripping_orbit::reset_step_factor_max()
{
	// Check if anything is actually changing here
	if ( _default_step_factor_max_ == _step_factor_max_ )
		return 0;

	clear_calcs();
	_step_factor_max_ = _default_step_factor_max_;
	return 0;
}
const int SALTSA::stripping_orbit::reset_step_factor_min()
{
	// Check if anything is actually changing here
	if ( _default_step_factor_min_ == _step_factor_min_ )
		return 0;

	clear_calcs();
	_step_factor_min_ = _default_step_factor_min_;
	return 0;
}
#endif

// Resetting tuning parameters
#if(1)

const int SALTSA::stripping_orbit::reset_tidal_stripping_amplification()
{
	// Check if anything is actually changing here
	if ( _default_tidal_stripping_amplification_ == _tidal_stripping_amplification_ )
		return 0;

	clear_calcs();
	_tidal_stripping_amplification_ = _default_tidal_stripping_amplification_;
	return 0;
}
const int SALTSA::stripping_orbit::reset_tidal_stripping_deceleration()
{
	// Check if anything is actually changing here
	if ( _default_tidal_stripping_deceleration_ == _tidal_stripping_deceleration_ )
		return 0;

	clear_calcs();
	_tidal_stripping_deceleration_ = _default_tidal_stripping_deceleration_;
	return 0;
}
const int SALTSA::stripping_orbit::reset_tidal_shocking_amplification()
{
	// Check if anything is actually changing here
	if ( _default_tidal_shocking_amplification_ == _tidal_shocking_amplification_ )
		return 0;

	clear_calcs();
	_tidal_shocking_amplification_ = _default_tidal_shocking_amplification_;
	return 0;
}
const int SALTSA::stripping_orbit::reset_tidal_shocking_persistance()
{
	// Check if anything is actually changing here
	if ( _default_tidal_shocking_persistance_ == _tidal_shocking_persistance_ )
		return 0;

	clear_calcs();
	_tidal_shocking_persistance_ = _default_tidal_shocking_persistance_;
	return 0;
}
const int SALTSA::stripping_orbit::reset_tidal_shocking_power()
{
	// Check if anything is actually changing here
	if ( _default_tidal_shocking_power_ == _tidal_shocking_power_ )
		return 0;

	clear_calcs();
	_tidal_shocking_power_ = _default_tidal_shocking_power_;
	return 0;
}
#endif

const int SALTSA::stripping_orbit::set_t_min( const double &new_t_min )
{
	_t_min_override_value_ = new_t_min;
	_override_t_min_ = true;
	return 0;
}

const int SALTSA::stripping_orbit::set_t_max( const double &new_t_max )
{
	_t_max_override_value_ = new_t_max;
	_override_t_max_ = true;
	return 0;
}

const int SALTSA::stripping_orbit::reset_t_min()
{
	_t_min_natural_value_ = DBL_MAX;
	for ( unsigned int i = 0; i < _x_spline_points_.size(); i++ )
	{
		if ( _x_spline_points_.at( i ).first < _t_min_natural_value_ )
			_t_min_natural_value_ = _x_spline_points_.at( i ).first;
	}
	_override_t_min_ = false;
	return 0;
}
const int SALTSA::stripping_orbit::reset_t_max()
{
	_t_max_natural_value_ = ( -DBL_MAX );
	for ( unsigned int i = 0; i < _x_spline_points_.size(); i++ )
	{
		if ( _x_spline_points_.at( i ).first > _t_max_natural_value_ )
			_t_max_natural_value_ = _x_spline_points_.at( i ).first;
	}
	_override_t_max_ = false;
	return 0;
}

// Functions for determining how calc() will be called
const int SALTSA::stripping_orbit::set_record_full_data(
		const bool new_record_full_data ) const
{
	// Check if anything is actually changing here
	if ( new_record_full_data == _record_full_data_ )
		return 0;

	clear_calcs();
	_record_full_data_ = new_record_full_data;
	return 0;
}

// Function to calculate stripping
const int SALTSA::stripping_orbit::calc( const bool silent ) const
{
	std::vector< int > segment_resolutions( 0 );
	std::vector< bool > segments_to_skip( 0 );
	int num_segments_to_skip = 0;
	double t_min_to_use = (
			_override_t_min_ ? _t_min_override_value_ : _t_min_natural_value_ );
	double t_max_to_use = (
			_override_t_max_ ? _t_max_override_value_ : _t_max_natural_value_ );
	clear_calcs();

	// First, determine number of segments we'll be using

	// Start by generating the cleaned_discontinuity_times list
	_cleaned_discontinuity_times_.clear();
	for ( unsigned int i = 0; i < _discontinuity_times_.size(); i++ )
	{
		// Check if this discontinuity is actually between t_min and t_max
		if ( ( _discontinuity_times_.at( i ) > t_min_to_use )
				&& ( _discontinuity_times_.at( i ) < t_max_to_use ) )
		{
			// It is, so add it to the cleaned list
			_cleaned_discontinuity_times_.push_back(
					_discontinuity_times_.at( i ) );
		}
	} // for(int i = 0; i < discontinuity_times.size(); i++)

	// Sort list of discontinuities
	std::sort( _cleaned_discontinuity_times_.begin(),
			_cleaned_discontinuity_times_.end() );

	_num_segments_ = _cleaned_discontinuity_times_.size() + 1;
	_orbit_segments_.resize( _num_segments_ );
	segment_resolutions.resize( _num_segments_, 0 );
	segments_to_skip.resize( _num_segments_, false );

	// Add each point to its proper segment

	for ( unsigned int i = 0; i < _x_spline_points_.size(); i++ )
	{
		double t = _x_spline_points_.at( i ).first;
		double x = _x_spline_points_.at( i ).second;
		int segment_counter = 0;
		for ( int j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_x_point( x, t );
	}

	for ( unsigned int i = 0; i < _y_spline_points_.size(); i++ )
	{
		double t = _y_spline_points_.at( i ).first;
		double y = _y_spline_points_.at( i ).second;
		int segment_counter = 0;
		for ( int j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_y_point( y, t );
	}

	for ( unsigned int i = 0; i < _z_spline_points_.size(); i++ )
	{
		double t = _z_spline_points_.at( i ).first;
		double z = _z_spline_points_.at( i ).second;
		int segment_counter = 0;
		for ( int j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_z_point( z, t );
	}

	for ( unsigned int i = 0; i < _vx_spline_points_.size(); i++ )
	{
		double t = _vx_spline_points_.at( i ).first;
		double vx = _vx_spline_points_.at( i ).second;
		int segment_counter = 0;
		for ( int j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_vx_point( vx, t );
	}

	for ( unsigned int i = 0; i < _vy_spline_points_.size(); i++ )
	{
		double t = _vy_spline_points_.at( i ).first;
		double vy = _vy_spline_points_.at( i ).second;
		int segment_counter = 0;
		for ( int j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_vy_point( vy, t );
	}

	for ( unsigned int i = 0; i < _vz_spline_points_.size(); i++ )
	{
		double t = _vz_spline_points_.at( i ).first;
		double vz = _vz_spline_points_.at( i ).second;
		int segment_counter = 0;
		for ( int j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_vz_point( vz, t );
	}

	for ( unsigned int i = 0; i < _vx_spline_unknown_points_.size(); i++ )
	{
		double t = _vx_spline_unknown_points_.at( i );
		int segment_counter = 0;
		for ( int j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_unknown_vx_point( t );
	}

	for ( unsigned int i = 0; i < _vy_spline_unknown_points_.size(); i++ )
	{
		double t = _vy_spline_unknown_points_.at( i );
		int segment_counter = 0;
		for ( int j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_unknown_vy_point( t );
	}

	for ( unsigned int i = 0; i < _vz_spline_unknown_points_.size(); i++ )
	{
		double t = _vz_spline_unknown_points_.at( i );
		int segment_counter = 0;
		for ( int j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_unknown_vz_point( t );
	}

	for ( unsigned int i = 0; i < _test_mass_spline_points_.size(); i++ )
	{
		double t = _test_mass_spline_points_.at( i ).first;
		double test_mass = _test_mass_spline_points_.at( i ).second;
		int segment_counter = 0;
		for ( int j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_test_mass_point( test_mass,
				t );
	}

	for ( unsigned int i = 0; i < _host_parameter_spline_points_.size(); i++ )
	{
		double t = _host_parameter_spline_points_.at( i ).first;
		std::vector< double > host_parameters =
				_host_parameter_spline_points_.at( i ).second;
		int segment_counter = 0;
		for ( int j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_host_parameter_point(
				host_parameters.size(), host_parameters, t );
	}

	// Adjust t_min and t_max for each segment to leave no gaps, get segment resolution,
	// and check for segments with too few spline points
	if ( t_max_to_use <= t_min_to_use )
	{
		if ( !silent )
			std::cerr << "ERROR: t_max is <= t_min for stripping_orbit.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	for ( int i = 0; i < _num_segments_; i++ )
	{
		// Adjust t_min and t_max
		double new_t_min;
		double new_t_max;
		if ( i == 0 )
		{
			new_t_min = t_min_to_use;
		}
		else
		{
			new_t_min = _cleaned_discontinuity_times_.at( i - 1 );
		}
		if ( i == _num_segments_ - 1 )
		{
			new_t_max = t_max_to_use;
		}
		else
		{
			new_t_max = _cleaned_discontinuity_times_.at( i );
		}
		_orbit_segments_.at( i ).set_t_min( new_t_min );
		_orbit_segments_.at( i ).set_t_max( new_t_max );

		// Adjust resolution for each segment
		int segment_resolution = (int)( _spline_resolution_
				* std::fabs( new_t_max - new_t_min )
				/ std::fabs( t_max_to_use - t_min_to_use ) ) + 1;
		segment_resolutions.at( i ) = segment_resolution;

		// Check if the segment has too few spline points or is too short
		if ( ( _orbit_segments_.at( i ).length() < 2 ) || ( segment_resolutions.at( i ) < 2 ) )
		{
			segments_to_skip.at( i ) = true;
			num_segments_to_skip++;
		}
	} // for(int i = 0; i < num_segments; i++ )

	// Check to make sure we have at least one segment we can calculate stripping for
	if ( num_segments_to_skip == _num_segments_ )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Cannot calculate stripping for any segments of orbit.\n";
		return UNSPECIFIED_ERROR;
	}

	// Now loop through and calculate stripping for each segment in turn

	density_profile *temp_satellite = NULL;
	density_profile *temp_host = NULL;

	int last_good_segment = -1;

	try
	{
		for ( int i = 0; i < _num_segments_; i++ )
		{
			double fmret = 1;
			if ( !( segments_to_skip.at( i ) ) )
			{

				if ( last_good_segment == -1 ) // Special handling for first good segment
				{
					temp_satellite = _init_satellite_ptr_->density_profile_clone();
					temp_host = _init_host_ptr_->density_profile_clone();
				}
				else
				{
					delete temp_satellite;
					delete temp_host;
					_orbit_segments_.at( last_good_segment ).clone_final_satellite(
							temp_satellite );
					_orbit_segments_.at( last_good_segment ).clone_final_host(
							temp_host );
					try
					{
						_orbit_segments_.at( i ).set_init_sum_deltarho(
								_orbit_segments_.at( last_good_segment ).final_sum_deltarho() );
						_orbit_segments_.at( i ).set_init_sum_gabdt(
								_orbit_segments_.at( last_good_segment ).final_sum_gabdt() );
					}
					catch ( ... )
					{
						if ( !silent )
							std::cerr
									<< "ERROR: Could not connect orbit segments properly.\n";
						std::cerr.flush();
						throw std::runtime_error("ERROR: Could not calculate stripping for orbit segment.\n");
					}
				}
				last_good_segment = i;

				// Pass parameters to the segment
				if ( _orbit_segments_.at( i ).set_init_satellite(
						temp_satellite ) )
					throw std::runtime_error("ERROR: Could not calculate stripping for orbit segment.\n");
				if ( _orbit_segments_.at( i ).set_init_host( temp_host ) )
					throw std::runtime_error("ERROR: Could not calculate stripping for orbit segment.\n");

				if ( _orbit_segments_.at( i ).set_resolution(
						segment_resolutions.at( i ) ) )
					throw std::runtime_error("ERROR: Could not calculate stripping for orbit segment.\n");
				if ( _orbit_segments_.at( i ).set_v_0(
						_v_0_ ) )
					throw std::runtime_error("ERROR: Could not calculate stripping for orbit segment.\n");
				if ( _orbit_segments_.at( i ).set_r_0(
						_r_0_ ) )
					throw std::runtime_error("ERROR: Could not calculate stripping for orbit segment.\n");
				if ( _orbit_segments_.at( i ).set_step_length_power(
						_step_length_power_ ) )
					throw std::runtime_error("ERROR: Could not calculate stripping for orbit segment.\n");
				if ( _orbit_segments_.at( i ).set_step_factor_min(
						_step_factor_min_ ) )
					throw std::runtime_error("ERROR: Could not calculate stripping for orbit segment.\n");
				if ( _orbit_segments_.at( i ).set_step_factor_max(
						_step_factor_max_ ) )
					throw std::runtime_error("ERROR: Could not calculate stripping for orbit segment.\n");
				if ( _orbit_segments_.at( i ).set_tidal_stripping_amplification(
						_tidal_stripping_amplification_ ) )
					throw std::runtime_error("ERROR: Could not calculate stripping for orbit segment.\n");
				if ( _orbit_segments_.at( i ).set_tidal_stripping_deceleration(
						_tidal_stripping_deceleration_ ) )
					throw std::runtime_error("ERROR: Could not calculate stripping for orbit segment.\n");
				if ( _orbit_segments_.at( i ).set_tidal_shocking_amplification(
						_tidal_shocking_amplification_ ) )
					throw std::runtime_error("ERROR: Could not calculate stripping for orbit segment.\n");
				if ( _orbit_segments_.at( i ).set_tidal_shocking_persistance(
						_tidal_shocking_persistance_ ) )
					throw std::runtime_error("ERROR: Could not calculate stripping for orbit segment.\n");
				if ( _orbit_segments_.at( i ).set_tidal_shocking_power(
						_tidal_shocking_power_ ) )
					throw std::runtime_error("ERROR: Could not calculate stripping for orbit segment.\n");
				if ( _orbit_segments_.at( i ).set_record_full_data(
						_record_full_data_ ) )
					throw std::runtime_error("ERROR: Could not calculate stripping for orbit segment.\n");
				if ( _orbit_segments_.at( i ).calc() )
				{
					throw std::runtime_error("ERROR: Could not calculate stripping for orbit segment.\n");
				}
				if (_orbit_segments_.at( i ).bad_result() )
				{
					if ( !silent )
						std::cerr
								<< "ERROR: Could not calculate stripping for orbit segment.\n";
					std::cerr.flush();
					throw std::runtime_error("ERROR: Could not calculate stripping for orbit segment.\n");
				}

				_orbit_segments_.at( i ).get_final_fmret( fmret );

				// Record this as the final good segment, in case we don't find any more afterward
				_final_good_segment_ = _orbit_segments_.begin() + i;
			}
			else
			{
				fmret = 1;
			}

			if ( _final_fmret_list_.size() > 0 )
				_final_fmret_list_.push_back( fmret * _final_fmret_list_.back() );
			else
				_final_fmret_list_.push_back( fmret );
		}

		delete temp_satellite;
		delete temp_host;

		_calculated_ = true;
	}
	catch (const std::exception &e)
	{
		clear_calcs();
		del_obj(temp_satellite);
		del_obj(temp_host);

		_calculated_ = true;
		_bad_result_ = true;

		return UNSPECIFIED_ERROR;
	}

	return 0;

}

const int SALTSA::stripping_orbit::set_satellite_output_parameters(
		const unsigned int num_out_parameters,
		const std::vector< bool > & new_output_parameters )
{
	_satellite_output_parameters_ = new_output_parameters;
	return 0;
}
const int SALTSA::stripping_orbit::set_satellite_parameter_unitconvs(
		const unsigned int num_out_parameters,
		const std::vector< double > & new_parameter_unitconvs )
{
	_satellite_parameter_unitconvs_ = new_parameter_unitconvs;
	return 0;
}

const int SALTSA::stripping_orbit::set_host_output_parameters(
		const unsigned int num_out_parameters,
		const std::vector< bool > & new_output_parameters )
{
	_host_output_parameters_ = new_output_parameters;
	return 0;
}
const int SALTSA::stripping_orbit::set_host_parameter_unitconvs(
		const unsigned int num_out_parameters,
		const std::vector< double > & new_parameter_unitconvs )
{
	_host_parameter_unitconvs_ = new_parameter_unitconvs;
	return 0;
}

const int SALTSA::stripping_orbit::clear_satellite_output_parameters()
{
	_satellite_output_parameters_.clear();
	return 0;
}
const int SALTSA::stripping_orbit::clear_satellite_parameter_unitconvs()
{
	_satellite_parameter_unitconvs_.clear();
	return 0;
}

const int SALTSA::stripping_orbit::clear_host_output_parameters()
{
	_host_output_parameters_.clear();
	return 0;
}
const int SALTSA::stripping_orbit::clear_host_parameter_unitconvs()
{
	_host_parameter_unitconvs_.clear();
	return 0;
}

const int SALTSA::stripping_orbit::print_full_data( std::ostream *out ) const
{
	if ( ( !_calculated_ ) || ( !_record_full_data_ ) )
	{
		if ( int errcode = set_record_full_data( true ) )
			return errcode + LOWER_LEVEL_ERROR;
		if ( calc() )
			return UNSPECIFIED_ERROR;
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}

	// Loop through segments, printing each of them in turn
	bool first_good_segment = true;
	for ( int i = 0; i < _num_segments_; i++ )
	{
		if ( _orbit_segments_.at( i ).length() >= 2 )
		{
			_orbit_segments_.at( i ).set_satellite_output_parameters(
					_satellite_output_parameters_.size(),
					_satellite_output_parameters_ );
			_orbit_segments_.at( i ).set_satellite_parameter_unitconvs(
					_satellite_parameter_unitconvs_.size(),
					_satellite_parameter_unitconvs_ );
			_orbit_segments_.at( i ).set_host_output_parameters(
					_host_output_parameters_.size(),
					_host_output_parameters_ );
			_orbit_segments_.at( i ).set_host_parameter_unitconvs(
					_host_parameter_unitconvs_.size(),
					_host_parameter_unitconvs_ );
			if ( first_good_segment ) // Special handling for first segment
			{
				if ( _orbit_segments_.at( i ).print_full_data( out, true, 1 ) )
					return 1; // Print header, mret_multiplier = 1
				first_good_segment = false;
			}
			else
			{
				if ( _orbit_segments_.at( i ).print_full_data( out, false,
						_final_fmret_list_.at( i - 1 ) ) )
					return 1; // Don't print header, mret_multiplier = fmret after last segment
			}
		}
	}
	return 0;
}

const int SALTSA::stripping_orbit::get_final_mret( double & mret ) const
{
	if ( !_calculated_ )
	{
		if ( calc() )
			return UNSPECIFIED_ERROR;
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}
	mret = _init_satellite_ptr_->mvir()*_final_fmret_list_.back();
	return 0;
}
const int SALTSA::stripping_orbit::get_final_fmret( double & fmret ) const
{
	if ( !_calculated_ )
	{
		if ( calc() )
			return UNSPECIFIED_ERROR;
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}
	fmret = _final_fmret_list_.back();
	return 0;
}

const int SALTSA::stripping_orbit::get_final_sum_deltarho(
double & final_sum_deltarho ) const
{
	if ( !_calculated_ )
	{
		if ( calc() )
			return UNSPECIFIED_ERROR;
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}
	_final_good_segment()->get_final_sum_deltarho( final_sum_deltarho );
	return 0;
}

const int SALTSA::stripping_orbit::get_final_sum_gabdt(
		gabdt & final_sum_gabdt ) const
{
	if ( !_calculated_ )
	{
		if ( calc() )
			return UNSPECIFIED_ERROR;
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}
	_final_good_segment()->get_final_sum_gabdt( final_sum_gabdt );
	return 0;
}

const int SALTSA::stripping_orbit::get_last_infall_time(
		double & t ) const
{
	if ( !_calculated_ )
	{
		if ( calc() )
			return UNSPECIFIED_ERROR;
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}
	t = _final_good_segment()->t_min_natural_value();
	return 0;
}

const int SALTSA::stripping_orbit::clone_final_satellite(
		SALTSA::density_profile * & final_satellite_clone ) const
{
	if ( !_calculated_ )
	{
		if ( calc() )
		{
			// If this function is being called, something being assigned is expected.

			if(_satellite_loaded_)
			{
				// Next most logical is to use the initial satellite.
				final_satellite_clone = _init_satellite_ptr_->density_profile_clone();
				return UNSPECIFIED_ERROR;
			}
			else
			{
				// Failsafe, we'll assign a new tNFW profile, just so it can be deleted
				// if that's attempted.
				final_satellite_clone = new SALTSA::tNFW_profile;
				return NOT_SET_UP_ERROR;
			}
		}
	}
	if( (_bad_result_) || (_final_good_segment()==_orbit_segments_.end()) )
	{
		final_satellite_clone = _init_satellite_ptr_->density_profile_clone(); // Sanest option in this case
		return UNSPECIFIED_ERROR;
	}
	else
	{
		final_satellite_clone = _final_good_segment()->final_satellite()->density_profile_clone();
		return 0;
	}

	return 0;
}

const int SALTSA::stripping_orbit::clone_final_host(
		SALTSA::density_profile * & final_host_clone ) const
{
	if ( !_calculated_ )
	{
		if ( calc() )
		{
			// If this function is being called, something being assigned is expected.

			if(_host_loaded_)
			{
				// Next most logical is to use the initial host.
				final_host_clone = _init_host_ptr_->density_profile_clone();
				return UNSPECIFIED_ERROR;
			}
			else
			{
				// Failsafe, we'll assign a new tNFW profile, just so it can be deleted
				// if that's attempted.
				final_host_clone = new SALTSA::tNFW_profile;
				return NOT_SET_UP_ERROR;
			}
		}
	}
	if( (_bad_result_) || (_final_good_segment()==_orbit_segments_.end()) )
	{
		final_host_clone = _init_host_ptr_->density_profile_clone(); // Sanest option in this case
		return UNSPECIFIED_ERROR;
	}
	else
	{
		final_host_clone = _final_good_segment()->final_host()->density_profile_clone();
		return 0;
	}
}

// Get final data (throws exception on failure)
const double SALTSA::stripping_orbit::final_mret() const
{
	double result = -1;

	if ( get_final_mret( result ) )
	{
		throw std::runtime_error("ERROR: Could not calculate stripping in stripping_orbit::final_mret.\n");
	}
	return result;
}

const double SALTSA::stripping_orbit::final_sum_deltarho() const
{
	double result = -1;

	if ( get_final_sum_deltarho( result ) )
	{
		throw std::runtime_error("ERROR: Could not calculate stripping in stripping_orbit::final_sum_deltarho.\n");
	}
	return result;
}

const double SALTSA::stripping_orbit::final_fmret() const
{
	double result;

	if ( get_final_fmret( result ) )
	{
		throw std::runtime_error("ERROR: Could not calculate stripping in stripping_orbit::final_fmret.\n");
	}
	return result;
}

const SALTSA::gabdt SALTSA::stripping_orbit::final_sum_gabdt() const
{
	SALTSA::gabdt result;

	if ( get_final_sum_gabdt( result ) )
	{
		throw std::runtime_error("ERROR: Could not calculate stripping in stripping_orbit::final_sum_gabdt.\n");
	}
	return result;
}

const double SALTSA::stripping_orbit::last_infall_time() const
{
	double result;

	if ( get_last_infall_time( result ) )
	{
		throw std::runtime_error("ERROR: Could not calculate stripping in stripping_orbit::last_infall_time.\n");
	}
	return result;
}

const SALTSA::density_profile * SALTSA::stripping_orbit::final_satellite() const
{
	if ( !_calculated_ )
	{
		if(_satellite_loaded_)
		{
			if ( calc() )
			{
				return _init_satellite_ptr_;
			}
		}
		else
		{
			throw std::runtime_error("ERROR: Attempt to call stripping_orbit::final_satellite() without init_satellite assigned.\n");
			return 0;
		}
	}
	if( (_bad_result_) || (_final_good_segment()==_orbit_segments_.end()) )
		return _init_satellite_ptr_; // Sanest option in this case
	else
		return _final_good_segment()->final_satellite();
}

const SALTSA::density_profile * SALTSA::stripping_orbit::final_host() const
{
	if ( !_calculated_ )
	{
		if(_host_loaded_)
		{
			if ( calc() )
			{
				return _init_host_ptr_;
			}
		}
		else
		{
			throw std::runtime_error("ERROR: Attempt to call stripping_orbit::final_host() without init_host assigned.\n");
			return 0;
		}
	}
	if( (_bad_result_) || (_final_good_segment()==_orbit_segments_.end()) )
		return _init_host_ptr_; // Sanest option in this case
	else
		return _final_good_segment()->final_host();
}

#endif // end SALTSA::stripping_orbit_segment class function definitions

// SALTSA::stripping_orbit_segment class method implementations
#if (1)

SALTSA::stripping_orbit_segment::stripping_orbit_segment()
{
	_init();
}

const int SALTSA::stripping_orbit_segment::_init()
{
	_current_satellite_in_use_ = false;
	_current_host_in_use_ = false;
	_host_loaded_ = false;
	_satellite_loaded_ = false;
	return clear();
}

SALTSA::stripping_orbit_segment::stripping_orbit_segment(
		const stripping_orbit_segment &other_orbit_spline )
{

	// Integration parameters
#if(1)
	_spline_resolution_ = other_orbit_spline._spline_resolution_;
	_v_0_ = other_orbit_spline._v_0_;
	_r_0_ = other_orbit_spline._r_0_;
	_step_length_power_ = other_orbit_spline._step_length_power_;
	_step_factor_max_ = other_orbit_spline._step_factor_max_;
	_step_factor_min_ = other_orbit_spline._step_factor_min_;
#endif

	// Tuning values
#if(1)
	_tidal_stripping_amplification_ = other_orbit_spline._tidal_stripping_amplification_;
	_tidal_stripping_deceleration_ = other_orbit_spline._tidal_stripping_deceleration_;
	_tidal_shocking_amplification_ = other_orbit_spline._tidal_shocking_amplification_;
	_tidal_shocking_persistance_ = other_orbit_spline._tidal_shocking_persistance_;
	_tidal_shocking_power_ = other_orbit_spline._tidal_shocking_power_;
#endif

	_rt_list_ = other_orbit_spline._rt_list_;
	_rt_ratio_list_ = other_orbit_spline._rt_ratio_list_;
	_delta_rho_list_ = other_orbit_spline._delta_rho_list_;
	_sum_delta_rho_list_ = other_orbit_spline._sum_delta_rho_list_;
	_x_data_ = other_orbit_spline._x_data_;
	_y_data_ = other_orbit_spline._y_data_;
	_z_data_ = other_orbit_spline._z_data_;
	_vx_data_ = other_orbit_spline._vx_data_;
	_vy_data_ = other_orbit_spline._vy_data_;
	_vz_data_ = other_orbit_spline._vz_data_;
	_t_min_natural_value_ = other_orbit_spline._t_min_natural_value_;
	_t_max_natural_value_ = other_orbit_spline._t_max_natural_value_;
	_satellite_parameter_data_ = other_orbit_spline._satellite_parameter_data_;
	_host_parameter_data_ = other_orbit_spline._host_parameter_data_;
	_num_parameters_ = other_orbit_spline._num_parameters_;
	mret_list = other_orbit_spline.mret_list;
	_satellite_parameter_unitconvs_ =
			other_orbit_spline._satellite_parameter_unitconvs_;
	_satellite_output_parameters_ =
			other_orbit_spline._satellite_output_parameters_;
	_host_parameter_unitconvs_ = other_orbit_spline._host_parameter_unitconvs_;
	_host_output_parameters_ = other_orbit_spline._host_output_parameters_;
	_gabdt_list_ = other_orbit_spline._gabdt_list_;
	_sum_gabdt_list_ = other_orbit_spline._sum_gabdt_list_;
	_phase_list_ = other_orbit_spline._phase_list_;
	_phase_output_list_ = other_orbit_spline._phase_output_list_;
	_record_full_data_ = other_orbit_spline._record_full_data_;
	_x_spline_ = other_orbit_spline._x_spline_;
	_y_spline_ = other_orbit_spline._y_spline_;
	_z_spline_ = other_orbit_spline._z_spline_;
	_vx_spline_ = other_orbit_spline._vx_spline_;
	_vy_spline_ = other_orbit_spline._vy_spline_;
	_vz_spline_ = other_orbit_spline._vz_spline_;
	_test_mass_spline_ = other_orbit_spline._test_mass_spline_;
	_host_parameter_splines_ = other_orbit_spline._host_parameter_splines_;
	_init_host_ptr_ = other_orbit_spline._init_host_ptr_;
	_init_satellite_ptr_ = other_orbit_spline._init_satellite_ptr_;
	_init_sum_delta_rho_ = other_orbit_spline._init_sum_delta_rho_;
	_init_sum_gabdt_ = other_orbit_spline._init_sum_gabdt_;
	_calculated_ = other_orbit_spline._calculated_;
	_bad_result_ = other_orbit_spline._bad_result_;
	_current_satellite_in_use_ = other_orbit_spline._current_satellite_in_use_;
	_current_host_in_use_ = other_orbit_spline._current_host_in_use_;
	_host_loaded_ = other_orbit_spline._host_loaded_;
	_satellite_loaded_ = other_orbit_spline._satellite_loaded_;
	_using_private_init_host_ = other_orbit_spline._using_private_init_host_;
	_using_private_init_satellite_ =
			other_orbit_spline._using_private_init_satellite_;
	_private_tNFW_init_host_ = other_orbit_spline._private_tNFW_init_host_;
	_private_tNFW_init_satellite_ =
			other_orbit_spline._private_tNFW_init_satellite_;
	_evolving_host_ = other_orbit_spline._evolving_host_;
	_override_t_max_ = other_orbit_spline._override_t_max_;
	_override_t_min_ = other_orbit_spline._override_t_min_;
	_t_max_override_val_ = other_orbit_spline._t_max_override_val_;
	_t_min_override_val_ = other_orbit_spline._t_min_override_val_;

	if ( _current_satellite_in_use_ )
	{
		_current_satellite_ptr_ =
				other_orbit_spline._current_satellite_ptr_->density_profile_clone();
	}
	else
	{
		_current_satellite_ptr_ = 0;
	}
	if ( _current_host_in_use_ )
	{
		_current_host_ptr_ =
				other_orbit_spline._current_host_ptr_->density_profile_clone();
	}
	else
	{
		_current_host_ptr_ = 0;
	}
	if ( _using_private_init_host_ )
	{
		_init_host_ptr_ = &_private_tNFW_init_host_;
	}
	if ( _using_private_init_satellite_ )
	{
		_init_satellite_ptr_ = &_private_tNFW_init_satellite_;
	}
}

SALTSA::stripping_orbit_segment & SALTSA::stripping_orbit_segment::operator=(
		const stripping_orbit_segment &other_orbit_spline )
{
	if ( this != &other_orbit_spline )
	{
		clear();

		// Integration parameters
#if(1)
		_spline_resolution_ = other_orbit_spline._spline_resolution_;
		_v_0_ = other_orbit_spline._v_0_;
		_r_0_ = other_orbit_spline._r_0_;
		_step_length_power_ = other_orbit_spline._step_length_power_;
		_step_factor_max_ = other_orbit_spline._step_factor_max_;
		_step_factor_min_ = other_orbit_spline._step_factor_min_;
#endif

	// Tuning parameters
#if(1)
		_tidal_stripping_amplification_ = other_orbit_spline._tidal_stripping_amplification_;
		_tidal_stripping_deceleration_ = other_orbit_spline._tidal_stripping_deceleration_;
		_tidal_shocking_amplification_ = other_orbit_spline._tidal_shocking_amplification_;
		_tidal_shocking_persistance_ = other_orbit_spline._tidal_shocking_persistance_;
		_tidal_shocking_power_ = other_orbit_spline._tidal_shocking_power_;
#endif

		_rt_list_ = other_orbit_spline._rt_list_;
		_rt_ratio_list_ = other_orbit_spline._rt_ratio_list_;
		_delta_rho_list_ = other_orbit_spline._delta_rho_list_;
		_sum_delta_rho_list_ = other_orbit_spline._sum_delta_rho_list_;
		_x_data_ = other_orbit_spline._x_data_;
		_y_data_ = other_orbit_spline._y_data_;
		_z_data_ = other_orbit_spline._z_data_;
		_vx_data_ = other_orbit_spline._vx_data_;
		_vy_data_ = other_orbit_spline._vy_data_;
		_vz_data_ = other_orbit_spline._vz_data_;
		_t_min_natural_value_ = other_orbit_spline._t_min_natural_value_;
		_t_max_natural_value_ = other_orbit_spline._t_max_natural_value_;
		_satellite_parameter_data_ =
				other_orbit_spline._satellite_parameter_data_;
		_host_parameter_data_ = other_orbit_spline._host_parameter_data_;
		_num_parameters_ = other_orbit_spline._num_parameters_;
		mret_list = other_orbit_spline.mret_list;
		_satellite_parameter_unitconvs_ =
				other_orbit_spline._satellite_parameter_unitconvs_;
		_satellite_output_parameters_ =
				other_orbit_spline._satellite_output_parameters_;
		_host_parameter_unitconvs_ =
				other_orbit_spline._host_parameter_unitconvs_;
		_host_output_parameters_ = other_orbit_spline._host_output_parameters_;
		_gabdt_list_ = other_orbit_spline._gabdt_list_;
		_sum_gabdt_list_ = other_orbit_spline._sum_gabdt_list_;
		_phase_list_ = other_orbit_spline._phase_list_;
		_phase_output_list_ = other_orbit_spline._phase_output_list_;
		_record_full_data_ = other_orbit_spline._record_full_data_;
		_x_spline_ = other_orbit_spline._x_spline_;
		_y_spline_ = other_orbit_spline._y_spline_;
		_z_spline_ = other_orbit_spline._z_spline_;
		_vx_spline_ = other_orbit_spline._vx_spline_;
		_vy_spline_ = other_orbit_spline._vy_spline_;
		_vz_spline_ = other_orbit_spline._vz_spline_;
		_test_mass_spline_ = other_orbit_spline._test_mass_spline_;
		_host_parameter_splines_ = other_orbit_spline._host_parameter_splines_;
		_init_host_ptr_ = other_orbit_spline._init_host_ptr_;
		_init_satellite_ptr_ = other_orbit_spline._init_satellite_ptr_;
		_init_sum_delta_rho_ = other_orbit_spline._init_sum_delta_rho_;
		_init_sum_gabdt_ = other_orbit_spline._init_sum_gabdt_;
		_calculated_ = other_orbit_spline._calculated_;
		_bad_result_ = other_orbit_spline._bad_result_;
		_current_satellite_in_use_ =
				other_orbit_spline._current_satellite_in_use_;
		_current_host_in_use_ = other_orbit_spline._current_host_in_use_;
		_host_loaded_ = other_orbit_spline._host_loaded_;
		_satellite_loaded_ = other_orbit_spline._satellite_loaded_;
		_using_private_init_host_ =
				other_orbit_spline._using_private_init_host_;
		_using_private_init_satellite_ =
				other_orbit_spline._using_private_init_satellite_;
		_private_tNFW_init_host_ = other_orbit_spline._private_tNFW_init_host_;
		_private_tNFW_init_satellite_ =
				other_orbit_spline._private_tNFW_init_satellite_;
		_evolving_host_ = other_orbit_spline._evolving_host_;
		_override_t_max_ = other_orbit_spline._override_t_max_;
		_override_t_min_ = other_orbit_spline._override_t_min_;
		_t_max_override_val_ = other_orbit_spline._t_max_override_val_;
		_t_min_override_val_ = other_orbit_spline._t_min_override_val_;

		if ( _current_satellite_in_use_ )
		{
			_current_satellite_ptr_ =
					other_orbit_spline._current_satellite_ptr_->density_profile_clone();
		}
		else
		{
			_current_satellite_ptr_ = 0;
		}
		if ( _current_host_in_use_ )
		{
			_current_host_ptr_ =
					other_orbit_spline._current_host_ptr_->density_profile_clone();
		}
		else
		{
			_current_host_ptr_ = 0;
		}
		if ( _using_private_init_host_ )
		{
			_init_host_ptr_ = &_private_tNFW_init_host_;
		}
		if ( _using_private_init_satellite_ )
		{
			_init_satellite_ptr_ = &_private_tNFW_init_satellite_;
		}
	}
	return *this;
}

SALTSA::stripping_orbit_segment::stripping_orbit_segment(
		const density_profile *init_init_host,
		const density_profile *init_init_satellite, const int init_resolution )
{
	_init();
	_init_host_ptr_ = init_init_host;
	_init_satellite_ptr_ = init_init_satellite;
	_spline_resolution_ = init_resolution;
	_current_satellite_ptr_ = _init_satellite_ptr_->density_profile_clone();
	_current_satellite_in_use_ = true;
	_current_host_ptr_ = _init_host_ptr_->density_profile_clone();
	_current_host_in_use_ = true;
	_host_loaded_ = true;
	_satellite_loaded_ = true;

}

SALTSA::stripping_orbit_segment::~stripping_orbit_segment()
{
	if ( _current_satellite_in_use_ )
	{
		delete _current_satellite_ptr_;
	}
	if ( _current_host_in_use_ )
	{
		delete _current_host_ptr_;
	}
}

SALTSA::stripping_orbit_segment *SALTSA::stripping_orbit_segment::stripping_orbit_spline_clone() const
{
	return new stripping_orbit_segment( *this );
}

const int SALTSA::stripping_orbit_segment::clear()
{
	clear_calcs();

	// Integration parameters
#if(1)
	_spline_resolution_ = SALTSA::stripping_orbit::default_spline_resolution();
	_v_0_ = SALTSA::stripping_orbit::default_spline_resolution();
	_r_0_ = SALTSA::stripping_orbit::default_spline_resolution();
	_step_length_power_ = SALTSA::stripping_orbit::default_spline_resolution();
	_step_factor_max_ = SALTSA::stripping_orbit::default_spline_resolution();
	_step_factor_min_ = SALTSA::stripping_orbit::default_spline_resolution();
#endif

	// Tuning values
#if(1)
	_tidal_stripping_amplification_ = SALTSA::stripping_orbit::default_tidal_stripping_amplification();
	_tidal_stripping_deceleration_ = SALTSA::stripping_orbit::default_tidal_stripping_deceleration();
	_tidal_shocking_amplification_ = SALTSA::stripping_orbit::default_tidal_shocking_amplification();
	_tidal_shocking_persistance_ = SALTSA::stripping_orbit::default_tidal_shocking_persistance();
	_tidal_shocking_power_ = SALTSA::stripping_orbit::default_tidal_shocking_power();
#endif


	_phase_list_.clear();

	_x_data_.clear();
	_y_data_.clear();
	_z_data_.clear();
	_vx_data_.clear();
	_vy_data_.clear();
	_vz_data_.clear();

	_satellite_parameter_data_.clear();
	_host_parameter_data_.clear();
	_satellite_parameter_unitconvs_.clear();
	_satellite_output_parameters_.clear();
	_host_parameter_unitconvs_.clear();
	_host_output_parameters_.clear();

	_x_spline_.clear();
	_y_spline_.clear();
	_z_spline_.clear();
	_vx_spline_.clear();
	_vy_spline_.clear();
	_vz_spline_.clear();
	_test_mass_spline_.clear();
	_host_parameter_splines_.clear();

	_t_min_natural_value_ = DBL_MAX;
	_t_max_natural_value_ = ( -DBL_MAX );
	_num_parameters_ = 0;
	_record_full_data_ = false;

	_init_sum_delta_rho_ = 0;
	_init_sum_gabdt_.override_zero();

	_init_host_ptr_ = 0;
	_init_satellite_ptr_ = 0;

	if ( _current_satellite_in_use_ )
	{
		delete _current_satellite_ptr_;
	}
	if ( _current_host_in_use_ )
	{
		delete _current_host_ptr_;
	}

	_current_satellite_ptr_ = 0;
	_current_satellite_in_use_ = false;
	_current_host_ptr_ = 0;
	_current_host_in_use_ = false;
	_host_loaded_ = false;
	_satellite_loaded_ = false;
	_using_private_init_host_ = false;
	_using_private_init_satellite_ = false;
	_evolving_host_ = false;

	return 0;
}

const int SALTSA::stripping_orbit_segment::clear_calcs() const
{
	_rt_list_.clear();
	_rt_ratio_list_.clear();
	_phase_output_list_.clear();
	mret_list.clear();
	_delta_rho_list_.clear();
	_sum_delta_rho_list_.clear();
	_gabdt_list_.clear();
	_sum_gabdt_list_.clear();

	_calculated_ = false;
	_bad_result_ = false;
	return 0;
}

// Setting integration parameters
#if(1)
const int SALTSA::stripping_orbit_segment::set_resolution( const int new_resolution,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_resolution == _spline_resolution_ )
		return 0;

	if ( new_resolution < 2 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set resolution to value below minimum of 2.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	clear_calcs();
	_spline_resolution_ = new_resolution;
	return 0;
}
const int SALTSA::stripping_orbit_segment::set_v_0( const double new_v_0,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_v_0 == _v_0_ )
		return 0;

	if ( new_v_0 <= 0 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set v_0 to value at or below minimum of 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	clear_calcs();
	_v_0_ = new_v_0;
	return 0;
}
const int SALTSA::stripping_orbit_segment::set_r_0( const double new_r_0,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_r_0 == _r_0_ )
		return 0;

	if ( new_r_0 <= 0 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set r_0 to value at or below minimum of 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	clear_calcs();
	_r_0_ = new_r_0;
	return 0;
}
const int SALTSA::stripping_orbit_segment::set_step_length_power( const double new_step_length_power,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_step_length_power == _step_length_power_ )
		return 0;

	clear_calcs();
	_step_length_power_ = new_step_length_power;
	return 0;
}
const int SALTSA::stripping_orbit_segment::set_step_factor_max( const double new_step_factor_max,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_step_factor_max == _step_factor_max_ )
		return 0;

	if ( new_step_factor_max < 1 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set step_factor_max to value below minimum of 1.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	clear_calcs();
	_step_factor_max_ = new_step_factor_max;
	return 0;
}
const int SALTSA::stripping_orbit_segment::set_step_factor_min( const double new_step_factor_min,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_step_factor_min == _step_factor_min_ )
		return 0;

	if ( new_step_factor_min > 1 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set step_factor_min to value above minimum of 1.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	clear_calcs();
	_step_factor_min_ = new_step_factor_min;
	return 0;
}
#endif

// Setting tuning parameters
#if(1)

const int SALTSA::stripping_orbit_segment::set_tidal_stripping_amplification(
		const double new_tidal_stripping_amplification,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_stripping_amplification == _tidal_stripping_amplification_ )
		return 0;

	if ( new_tidal_stripping_amplification < 0 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set tidal_stripping_amplification to value below minimum of 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	clear_calcs();
	_tidal_stripping_amplification_ = new_tidal_stripping_amplification;
	return 0;
}
const int SALTSA::stripping_orbit_segment::set_tidal_stripping_deceleration(
		const double new_tidal_stripping_deceleration,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_stripping_deceleration == _tidal_stripping_deceleration_ )
		return 0;

	clear_calcs();
	_tidal_stripping_deceleration_ = new_tidal_stripping_deceleration;
	return 0;
}
const int SALTSA::stripping_orbit_segment::set_tidal_shocking_amplification(
		const double new_tidal_shocking_amplification,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_shocking_amplification == _tidal_shocking_amplification_ )
		return 0;

	if ( new_tidal_shocking_amplification < 0 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set tidal_shocking_amplification to value below minimum of 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	clear_calcs();
	_tidal_shocking_amplification_ = new_tidal_shocking_amplification;
	return 0;
}
const int SALTSA::stripping_orbit_segment::set_tidal_shocking_persistance(
		const double new_tidal_shocking_persistance,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_shocking_persistance == _tidal_shocking_persistance_ )
		return 0;

	if ( new_tidal_shocking_persistance <= 0 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to set tidal_shocking_persistance to value at or below minimum of 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	clear_calcs();
	_tidal_shocking_persistance_ = new_tidal_shocking_persistance;
	return 0;
}
const int SALTSA::stripping_orbit_segment::set_tidal_shocking_power(
		const double new_tidal_shocking_power,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_shocking_power == _tidal_shocking_power_ )
		return 0;

	clear_calcs();
	_tidal_shocking_power_ = new_tidal_shocking_power;
	return 0;
}
#endif


const int SALTSA::stripping_orbit_segment::add_point( const double &x,
		const double &y, const double &z, const double &t,
		const double new_mass )
{
	_calculated_ = false;
	try
	{
		_x_spline_.add_point( t, x );
		_y_spline_.add_point( t, y );
		_z_spline_.add_point( t, z );
		_vx_spline_.add_unknown_point( t );
		_vy_spline_.add_unknown_point( t );
		_vz_spline_.add_unknown_point( t );
		_test_mass_spline_.add_point( t, new_mass );
		if ( t < _t_min_natural_value_ )
			_t_min_natural_value_ = t;
		if ( t > _t_max_natural_value_ )
			_t_max_natural_value_ = t;
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_point().\n"
				<< "Exception: " << e.what();
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	catch (...)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_point().\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

const int SALTSA::stripping_orbit_segment::add_point( const double &x,
		const double &y, const double &z, const double &vx,
		const double &vy, const double &vz, const double &t,
		const double new_test_mass )
{
	_calculated_ = false;
	try
	{
		_x_spline_.add_point( t, x );
		_y_spline_.add_point( t, y );
		_z_spline_.add_point( t, z );
		_vx_spline_.add_point( t, vx );
		_vy_spline_.add_point( t, vy );
		_vz_spline_.add_point( t, vz );
		_test_mass_spline_.add_point( t, new_test_mass );
		if ( t < _t_min_natural_value_ )
			_t_min_natural_value_ = t;
		if ( t > _t_max_natural_value_ )
			_t_max_natural_value_ = t;
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_point().\n"
				<< "Exception: " << e.what();
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	catch (...)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_point().\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

const int SALTSA::stripping_orbit_segment::add_x_point(
		const double &x, const double &t )
{
	_calculated_ = false;
	try
	{
		_x_spline_.add_point( t, x );
		if ( t < _t_min_natural_value_ )
			_t_min_natural_value_ = t;
		if ( t > _t_max_natural_value_ )
			_t_max_natural_value_ = t;
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_x_point().\n"
				<< "Exception: " << e.what();
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	catch (...)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_x_point().\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

const int SALTSA::stripping_orbit_segment::add_y_point(
		const double &y, const double &t )
{
	_calculated_ = false;
	try
	{
		_y_spline_.add_point( t, y );
		if ( t < _t_min_natural_value_ )
			_t_min_natural_value_ = t;
		if ( t > _t_max_natural_value_ )
			_t_max_natural_value_ = t;
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_y_point().\n"
				<< "Exception: " << e.what();
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	catch (...)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_y_point().\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

const int SALTSA::stripping_orbit_segment::add_z_point(
		const double &z, const double &t )
{
	_calculated_ = false;
	try
	{
		_z_spline_.add_point( t, z );
		if ( t < _t_min_natural_value_ )
			_t_min_natural_value_ = t;
		if ( t > _t_max_natural_value_ )
			_t_max_natural_value_ = t;
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_z_point().\n"
				<< "Exception: " << e.what();
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	catch (...)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_z_point().\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

const int SALTSA::stripping_orbit_segment::add_vx_point(
		const double &vx, const double &t )
{
	_calculated_ = false;
	try
	{
		_vx_spline_.add_point( t, vx );
		if ( t < _t_min_natural_value_ )
			_t_min_natural_value_ = t;
		if ( t > _t_max_natural_value_ )
			_t_max_natural_value_ = t;
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_vx_point().\n"
				<< "Exception: " << e.what();
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	catch (...)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_vx_point().\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

const int SALTSA::stripping_orbit_segment::add_vy_point(
		const double &vy, const double &t )
{
	_calculated_ = false;
	try
	{
		_vy_spline_.add_point( t, vy );
		if ( t < _t_min_natural_value_ )
			_t_min_natural_value_ = t;
		if ( t > _t_max_natural_value_ )
			_t_max_natural_value_ = t;
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_vy_point().\n"
				<< "Exception: " << e.what();
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	catch (...)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_vy_point().\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

const int SALTSA::stripping_orbit_segment::add_vz_point(
		const double &vz, const double &t )
{
	_calculated_ = false;
	try
	{
		_vz_spline_.add_point( t, vz );
		if ( t < _t_min_natural_value_ )
			_t_min_natural_value_ = t;
		if ( t > _t_max_natural_value_ )
			_t_max_natural_value_ = t;
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_vz_point().\n"
				<< "Exception: " << e.what();
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	catch (...)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_vz_point().\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

const int SALTSA::stripping_orbit_segment::add_unknown_vx_point(
		const double &t )
{
	_calculated_ = false;
	try
	{
		_vx_spline_.add_unknown_point( t );
		if ( t < _t_min_natural_value_ )
			_t_min_natural_value_ = t;
		if ( t > _t_max_natural_value_ )
			_t_max_natural_value_ = t;
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_unknown_vx_point().\n"
				<< "Exception: " << e.what();
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	catch (...)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_unknown_vx_point().\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

const int SALTSA::stripping_orbit_segment::add_unknown_vy_point(
		const double &t )
{
	_calculated_ = false;
	try
	{
		_vy_spline_.add_unknown_point( t );
		if ( t < _t_min_natural_value_ )
			_t_min_natural_value_ = t;
		if ( t > _t_max_natural_value_ )
			_t_max_natural_value_ = t;
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_unknown_vy_point().\n"
				<< "Exception: " << e.what();
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	catch (...)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_unknown_vy_point().\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

const int SALTSA::stripping_orbit_segment::add_unknown_vz_point(
		const double &t )
{
	_calculated_ = false;
	try
	{
		_vz_spline_.add_unknown_point( t );
		if ( t < _t_min_natural_value_ )
			_t_min_natural_value_ = t;
		if ( t > _t_max_natural_value_ )
			_t_max_natural_value_ = t;
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_unknown_vz_point().\n"
				<< "Exception: " << e.what();
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	catch (...)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_unknown_vz_point().\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

const int SALTSA::stripping_orbit_segment::add_test_mass_point(
		const double test_mass, const double &t )
{
	_calculated_ = false;
	try
	{
		_test_mass_spline_.add_point( t, test_mass );
		if ( t < _t_min_natural_value_ )
			_t_min_natural_value_ = t;
		if ( t > _t_max_natural_value_ )
			_t_max_natural_value_ = t;
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_test_mass_point().\n"
				<< "Exception: " << e.what();
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	catch (...)
	{
		std::cerr << "ERROR: Exception in stripping_orbit_segment::add_test_mass_point().\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

const int SALTSA::stripping_orbit_segment::add_host_parameter_point(
		const unsigned int num_parameters,
		const std::vector< double > & parameters, const double &t,
		const bool silent )
{
	// Check num_parameters matches vector size
	if ( num_parameters != parameters.size() )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: num_parameters must == parameters.size() in stripping_orbit_segment::add_host_parameter_point.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	if ( num_parameters <= 0 )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: num_parameters must be > 0 in stripping_orbit_segment::add_host_parameter_point.\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	if ( _host_parameter_splines_.size() == 0 )
		_host_parameter_splines_.resize( num_parameters );

	if ( _host_parameter_splines_.size() != num_parameters )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: All parameter lists passed to stripping_orbit_segment must have same size.\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	for ( unsigned int i = 0; i < num_parameters; i++ )
	{
		_host_parameter_splines_[i].add_point( t, parameters[i] );
	}

	if ( _host_parameter_splines_[0].size() >= 2 )
		_evolving_host_ = true;
	if ( t < _t_min_natural_value_ )
		_t_min_natural_value_ = t;
	if ( t > _t_max_natural_value_ )
		_t_max_natural_value_ = t;

	return 0;
}

const int SALTSA::stripping_orbit_segment::_reserve( const int n,
		const bool silent ) const
{
	if ( n < 1 )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Attempt to reserve length of 0 or less in stripping_orbit_segment::reserve.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	_phase_list_.reserve( n );
	mret_list.reserve( n );
	_phase_output_list_.reserve( n );
	_rt_list_.reserve( n );
	_rt_ratio_list_.reserve( n );
	_delta_rho_list_.reserve( n );
	_sum_delta_rho_list_.reserve( n );
	_gabdt_list_.reserve( n );
	_sum_gabdt_list_.reserve( n );
	return 0;
}

const int SALTSA::stripping_orbit_segment::set_init_host(
		const density_profile *new_init_host )
{
	_init_host_ptr_ = new_init_host;
	if ( _current_host_in_use_ )
	{
		delete _current_host_ptr_;
	}
	_current_host_ptr_ = _init_host_ptr_->density_profile_clone();
	_current_host_in_use_ = true;
	_host_loaded_ = true;
	_using_private_init_host_ = false;
	_calculated_ = false;
	return 0;
}

const int SALTSA::stripping_orbit_segment::set_init_satellite(
		const density_profile *new_init_satellite )
{
	_init_satellite_ptr_ = new_init_satellite;
	if ( _current_satellite_in_use_ )
	{
		delete _current_satellite_ptr_;
	}
	_current_satellite_ptr_ = _init_satellite_ptr_->density_profile_clone();
	_current_satellite_in_use_ = true;
	_satellite_loaded_ = true;
	_using_private_init_satellite_ = false;
	_calculated_ = false;
	return 0;
}

const int SALTSA::stripping_orbit_segment::set_init_sum_deltarho(
		const double &new_init_sum_deltarho )
{
	_init_sum_delta_rho_ = new_init_sum_deltarho;
	return 0;
}

const int SALTSA::stripping_orbit_segment::set_init_sum_gabdt(
		const gabdt &new_init_sum_gabdt )
{
	_init_sum_gabdt_ = new_init_sum_gabdt;
	return 0;
}

const int SALTSA::stripping_orbit_segment::set_tNFW_init_satellite(
		const double &new_init_mvir0, const double z,
		const double new_init_c, const double new_init_tau )
{
	_using_private_init_satellite_ = true;
	_private_tNFW_init_satellite_ = tNFW_profile( new_init_mvir0, z,
			new_init_c, new_init_tau );
	_init_satellite_ptr_ = &_private_tNFW_init_satellite_;

	if ( _current_satellite_in_use_ )
	{
		delete _current_satellite_ptr_;
	}
	_current_satellite_ptr_ = _init_satellite_ptr_->density_profile_clone();
	_current_satellite_in_use_ = true;

	if ( _using_private_init_host_ )
	{
		if ( _private_tNFW_init_host_.z() != z )
			_private_tNFW_init_host_.set_z( z );
	}

	_satellite_loaded_ = true;
	_calculated_ = false;
	return 0;
}

const int SALTSA::stripping_orbit_segment::set_tNFW_host(
		const double &new_init_mvir0, const double z,
		const double new_init_c, const double new_init_tau )
{
	_using_private_init_host_ = true;
	_private_tNFW_init_host_ = tNFW_profile( new_init_mvir0, z, new_init_c,
			new_init_tau );
	_init_host_ptr_ = &_private_tNFW_init_host_;

	if ( _using_private_init_satellite_ )
	{
		if ( _private_tNFW_init_satellite_.z() != z )
		{
			_private_tNFW_init_satellite_.set_z( z );

			if ( _current_satellite_in_use_ )
			{
				delete _current_satellite_ptr_;
			}
			_current_satellite_ptr_ =
					_init_satellite_ptr_->density_profile_clone();
			_current_satellite_in_use_ = true;
		}
	}

	_host_loaded_ = true;
	_calculated_ = false;
	return 0;
}

const int SALTSA::stripping_orbit_segment::set_t_min(
		const double &new_t_min )
{
	_t_min_override_val_ = new_t_min;
	_override_t_min_ = true;
	return 0;
}

const int SALTSA::stripping_orbit_segment::set_t_max(
		const double &new_t_max )
{
	_t_max_override_val_ = new_t_max;
	_override_t_max_ = true;
	return 0;
}

const int SALTSA::stripping_orbit_segment::reset_t_min()
{
	_t_min_natural_value_ = DBL_MAX;
	for ( unsigned int i = 0; i < _x_spline_.size(); i++ )
	{
		if ( _x_spline_.sorted_data().at( i ).first < _t_min_natural_value_ )
			_t_min_natural_value_ = _x_spline_.sorted_data().at( i ).first;
	}
	_override_t_min_ = false;
	return 0;
}

const int SALTSA::stripping_orbit_segment::reset_t_max()
{
	_t_max_natural_value_ = ( -DBL_MAX );
	for ( unsigned int i = 0; i < _x_spline_.size(); i++ )
	{
		if ( _x_spline_.sorted_data().at( i ).first > _t_max_natural_value_ )
			_t_max_natural_value_ = _x_spline_.sorted_data().at( i ).first;
	}
	_override_t_min_ = false;
	return 0;
}

const int SALTSA::stripping_orbit_segment::clear_points()
{
	_x_spline_.clear();
	_y_spline_.clear();
	_z_spline_.clear();
	_vx_spline_.clear();
	_vy_spline_.clear();
	_vz_spline_.clear();
	_test_mass_spline_.clear();
	_calculated_ = false;
	return 0;
}
const int SALTSA::stripping_orbit_segment::clear_host_parameter_points()
{
	_host_parameter_splines_.clear();
	_calculated_ = false;
	return 0;
}
const int SALTSA::stripping_orbit_segment::clear_init_sum_deltarho()
{
	_init_sum_delta_rho_ = 0;
	_calculated_ = false;
	return 0;
}
const int SALTSA::stripping_orbit_segment::clear_init_sum_gabdt()
{
	_init_sum_gabdt_.clear();
	_init_sum_gabdt_.override_zero();
	_calculated_ = false;
	return 0;
}
const int SALTSA::stripping_orbit_segment::clear_init_satellite()
{
	_init_satellite_ptr_ = 0;
	if ( _current_satellite_in_use_ )
	{
		delete _current_satellite_ptr_;
	}
	_current_satellite_ptr_ = 0;
	_current_satellite_in_use_ = false;
	_satellite_loaded_ = false;
	_using_private_init_satellite_ = false;
	_calculated_ = false;
	return 0;
}
const int SALTSA::stripping_orbit_segment::clear_init_host()
{
	_init_host_ptr_ = 0;
	if ( _current_host_in_use_ )
	{
		delete _current_host_ptr_;
	}
	_current_host_ptr_ = 0;
	_current_host_in_use_ = false;
	_host_loaded_ = false;
	_using_private_init_host_ = false;
	_calculated_ = false;
	return 0;
}

const unsigned int SALTSA::stripping_orbit_segment::length() const
{
	unsigned int result = INT_MAX;
	if ( _x_spline_.size() < result )
		result = _x_spline_.size();
	if ( _y_spline_.size() < result )
		result = _y_spline_.size();
	if ( _z_spline_.size() < result )
		result = _z_spline_.size();
	return result;
}

// Functions for determining how calc() will be called
const int SALTSA::stripping_orbit_segment::set_record_full_data(
		const bool new_record_full_data ) const
{
	// Check if anything is actually changing here
	if ( new_record_full_data == _record_full_data_ )
		return 0;

	clear_calcs();
	_record_full_data_ = new_record_full_data;
	return 0;
}

// Function to calculate stripping
const int SALTSA::stripping_orbit_segment::calc( const bool silent ) const
{
	double r;
	double v, vt, vr;
	double t, t_step, t_min_to_use, t_max_to_use;
	double current_deltarho;
	std::vector< double > satellite_parameters;
	std::vector< double > host_parameters;
	int num_satellite_parameters = 0;
	int num_host_parameters = 0;
	int counter = 0;
	double step_length_factor = 1;
	double mret = 1;

	if ( ( !_host_loaded_ ) || ( !_satellite_loaded_ ) )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host or satellite not loaded in stripping_orbit_segment::calc.\n";
		return NOT_SET_UP_ERROR;
	}

	if ( length() < 2 )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Too few data points to calculate stripping in stripping_orbit_segment::calc().\n";
		return UNSPECIFIED_ERROR;
	}

	// Initialisation

	if ( _override_t_min_ )
	{
		t_min_to_use = _t_min_override_val_;
	}
	else
	{
		t_min_to_use = _t_min_natural_value_;
	}

	if ( _override_t_max_ )
	{
		t_max_to_use = _t_max_override_val_;
	}
	else
	{
		t_max_to_use = _t_max_natural_value_;
	}

	t_step = ( t_max_to_use - t_min_to_use ) / _spline_resolution_;
	if ( t_step <= 0 )
	{
		if ( !silent )
			std::cerr << "ERROR t_max <= t_min for stripping_orbit_segment!\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	gabdt current_gabdt = _init_sum_gabdt_;
	t = 0;
	r = 0;
	v = vt = vr = 0;
	current_deltarho = 0;

	clear_calcs();
	_reserve( 2 * _spline_resolution_ + 1 );

	// Setup

	// Velocity splines setup
	_vx_spline_.set_spline_ptr( &_x_spline_ );
	_vy_spline_.set_spline_ptr( &_y_spline_ );
	_vz_spline_.set_spline_ptr( &_z_spline_ );

	mret_list.push_back( 1 );
	_delta_rho_list_.push_back( current_deltarho );
	_sum_delta_rho_list_.push_back( _init_sum_delta_rho_ );

	_sum_gabdt_list_.push_back( current_gabdt );
	current_gabdt.set_host_ptr( _init_host_ptr_ );
	current_gabdt.override_zero();
	_gabdt_list_.push_back( current_gabdt );
	current_gabdt.override_zero();

	// Calculate radius and velocity at init position
	r = SALTSA::dist3d( _x_spline_( t_min_to_use ),
			_y_spline_( t_min_to_use ), _z_spline_( t_min_to_use ) );
	v = SALTSA::dist3d( _vx_spline_( t_min_to_use ),
			_vy_spline_( t_min_to_use ), _vz_spline_( t_min_to_use ) );
	vr = SALTSA::dot_product( _vx_spline_( t_min_to_use ),
			_vy_spline_( t_min_to_use ), _vz_spline_( t_min_to_use ),
			_x_spline_( t_min_to_use ), _y_spline_( t_min_to_use ),
			_z_spline_( t_min_to_use ) ) / r;
	if ( std::fabs( vr ) >= std::fabs( v ) )
		vt = 0;
	else
		vt = SALTSA::quad_sub( v, vr );

	// Make sure current_satellite and current_host are set up correctly
	if ( !_current_satellite_in_use_ )
	{
		_current_satellite_ptr_ =
				_init_satellite_ptr_->density_profile_clone();
		_current_satellite_in_use_ = true;
	}
	if ( !_current_host_in_use_ )
	{
		_current_host_ptr_ = _init_host_ptr_->density_profile_clone();
		_current_host_in_use_ = true;
	}

	step_length_factor = _step_length_factor(v, r);

	_rt_list_.push_back(
			get_rt( _current_host_ptr_, _current_satellite_ptr_, r,
					vr, vt, t_step, current_deltarho ) );
	_rt_ratio_list_.push_back( _rt_list_.at( 0 ) / _rvir() );

	// Get satellite parameters at init position
	if ( _record_full_data_ )
	{
		if ( _current_satellite_ptr_->get_parameters( num_satellite_parameters,
				satellite_parameters ) )
		{
			if ( !silent )
				std::cerr
						<< "WARNING: Cannot get satellite parameters. No extra information can be output.\n";
		}
		else
		{
			_satellite_parameter_data_.push_back( satellite_parameters );
		}
		if ( _current_host_ptr_->get_parameters( num_host_parameters,
				host_parameters ) )
		{
			if ( !silent )
				std::cerr
						<< "WARNING: Cannot get host parameters. No extra information can be output.\n";
		}
		else
		{
			_host_parameter_data_.push_back( host_parameters );
		}
	}

	// Loop over time
	for ( t = t_min_to_use; t < t_max_to_use;
			t += t_step * step_length_factor )
	{
		if ( _record_full_data_ )
			_phase_output_list_.push_back(
					phase( _x_spline_( t ), _y_spline_( t ), _z_spline_( t ),
							_vx_spline_( t ), _vy_spline_( t ),
							_vz_spline_( t ), t ) );
		counter += 1;
		// Calculate radius and velocity at this position
		try
		{
			r = SALTSA::dist3d( _x_spline_( t ), _y_spline_( t ),
					_z_spline_( t ) );
			v = SALTSA::quad_add( _vx_spline_( t ), _vy_spline_( t ),
					_vz_spline_( t ) );
			vr = SALTSA::dot_product( _vx_spline_( t ), _vy_spline_( t ),
					_vz_spline_( t ), _x_spline_( t ), _y_spline_( t ),
					_z_spline_( t ) ) / r;
			if ( std::fabs( vr ) >= std::fabs( v ) )
				vt = 0;
			else
				vt = SALTSA::quad_sub( v, vr );
			step_length_factor = _step_length_factor(v, r);

			// Update host for this timestep if necessary
			if ( _evolving_host_ )
			{
				num_host_parameters = _host_parameter_splines_.size();
				host_parameters.resize( num_host_parameters );
				for ( int i = 0; i < num_host_parameters; i++ )
				{
					host_parameters.at( i ) = _host_parameter_splines_.at( i )(
							t );
				}
				if ( _current_host_ptr_->set_parameters( num_host_parameters,
						host_parameters ) )
				{
					if ( !silent )
						std::cerr
								<< "ERROR: Cannot update host in stripping_orbit.\n";
					return 1;
				}
			}

			// Calculate effects of tidal stripping and shocks
			// Effects of stripping
			mret = tidal_strip_retained( _current_host_ptr_,
					_current_satellite_ptr_, r, vr, vt,
					t_step * step_length_factor,
					_sum_delta_rho_list_.at( counter - 1 ) );
			mret_list.push_back( mret_list.at( counter - 1 ) * mret );
			_current_satellite_ptr_->truncate_to_fraction( mret );

			//mass_list.push_back( mass_list.at(counter-1) * SALTSA::tidal_strip_retained( host, init_mass, tau_list.at(counter-1), r, vr, vt, t_step*step_length_factor, sum_deltarho_list.at(counter-1) ) );
			//tau_list.push_back( SALTSA::taufm( mass_list.at(counter)/safe_d(mass_list.at(0)), conc ) );

			// Effects of shocking

			double t_shock = r / safe_d(v);
			double t_recover = _current_satellite_ptr_->rhmvir()
					/ safe_d(_current_satellite_ptr_->vhmvir());
			double x;
			double gabdt_scaling_factor;

			if(isbad(t_recover) or (t_recover<=0))
			{
				x = 0;
				gabdt_scaling_factor = 0;
				_bad_result_ = true;
			}
			else
			{
				x = max(0,t_shock / safe_d(t_recover) / ( 2 * pi ));
				gabdt_scaling_factor = max(0,1
									- ( t_step * step_length_factor
											/ safe_d( t_recover * _tidal_shocking_persistance_ ) ) );
			}

			current_gabdt.set_host_ptr( _current_host_ptr_ );
			current_gabdt.set_pos( _x_spline_( t ), _y_spline_( t ),
					_z_spline_( t ) );
			current_gabdt.set_dt( t_step * step_length_factor );
			current_gabdt.calc_dv();
			_gabdt_list_.push_back( current_gabdt );
			_sum_gabdt_list_.push_back(
					_sum_gabdt_list_.back() * gabdt_scaling_factor
							+ current_gabdt );

			current_deltarho = _delta_rho( counter, x,
                                    t_step * step_length_factor );

			_delta_rho_list_.push_back( current_deltarho );

			_sum_delta_rho_list_.push_back(
							_sum_delta_rho_list_.back() + current_deltarho);

			_rt_list_.push_back(
					get_rt( _current_host_ptr_,
							_current_satellite_ptr_, r, vr, vt,
							t_step * step_length_factor,
							_sum_delta_rho_list_.at( counter - 1 ) ) );
			_rt_ratio_list_.push_back( _rt_list_.at( counter ) / _rvir( 0 ) );

			if ( !_current_satellite_ptr_->get_parameters(
					num_satellite_parameters, satellite_parameters ) )
			{
				_satellite_parameter_data_.push_back( satellite_parameters );
			}

			if ( !_current_host_ptr_->get_parameters( num_host_parameters,
					host_parameters ) )
			{
				_host_parameter_data_.push_back( host_parameters );
			}

		}
		catch ( std::exception & e )
		{
			_calculated_ = true; // We know there's an issue, so we won't calculate again unless something changes
			_bad_result_ = true;
			return UNSPECIFIED_ERROR;
		}
	}
	if ( _record_full_data_ )
		_phase_output_list_.push_back(
				phase( _x_spline_( t ), _y_spline_( t ), _z_spline_( t ),
						_vx_spline_( t ), _vy_spline_( t ), _vz_spline_( t ),
						t ) );
	_calculated_ = true;
	return 0;
}

const int SALTSA::stripping_orbit_segment::set_satellite_output_parameters(
		const unsigned int num_out_parameters,
		const std::vector< bool > & new_satellite_output_parameters )
{
	_satellite_output_parameters_ = new_satellite_output_parameters;
	return 0;
}
const int SALTSA::stripping_orbit_segment::set_satellite_parameter_unitconvs(
		const unsigned int num_out_parameters,
		const std::vector< double > &new_satellite_parameter_unitconvs )
{
	_satellite_parameter_unitconvs_ = new_satellite_parameter_unitconvs;
	return 0;
}

const int SALTSA::stripping_orbit_segment::set_host_output_parameters(
		const unsigned int num_out_parameters,
		const std::vector< bool > & new_host_output_parameters )
{
	_host_output_parameters_ = new_host_output_parameters;
	return 0;
}
const int SALTSA::stripping_orbit_segment::set_host_parameter_unitconvs(
		const unsigned int num_out_parameters,
		const std::vector< double > & new_host_parameter_unitconvs )
{
	_host_parameter_unitconvs_ = new_host_parameter_unitconvs;
	return 0;
}

const int SALTSA::stripping_orbit_segment::clear_satellite_output_parameters()
{
	_satellite_output_parameters_.clear();
	return 0;
}
const int SALTSA::stripping_orbit_segment::clear_satellite_parameter_unitconvs()
{
	_satellite_parameter_unitconvs_.clear();
	return 0;
}

const int SALTSA::stripping_orbit_segment::clear_host_output_parameters()
{
	_host_output_parameters_.clear();
	return 0;
}
const int SALTSA::stripping_orbit_segment::clear_host_parameter_unitconvs()
{
	_host_parameter_unitconvs_.clear();
	return 0;
}

const int SALTSA::stripping_orbit_segment::print_full_data(
		std::ostream *out, const bool include_header,
		const double mret_multiplier, const bool silent ) const
{
	const int num_columns_base = 16;
	int num_extra_satellite_columns = 0;
	int num_extra_host_columns = 0;
	int num_rows = 0;
	std::vector< std::string > header;
	std::vector< std::vector< std::string > > data;
	stringstream ss;

	// Calculate data if necessary
	if ( ( !_calculated_ ) || ( !_record_full_data_ ) )
	{
		if ( int errcode = set_record_full_data( true ) )
			return errcode + LOWER_LEVEL_ERROR;
		if ( calc() )
			return 1;
	}

	num_rows = _phase_output_list_.size();
	if ( num_rows < 2 )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: No stripping data to output in stripping_orbit_segment::print_full_data.\n";
		return UNSPECIFIED_ERROR;
	}

	// Check that vectors for density_profile params are sane.
	// If so, set up extra columns. If not, no extra columns.
	num_extra_satellite_columns = 0;
	if ( !( ( _satellite_output_parameters_.size() == 0 )
			&& ( _satellite_parameter_unitconvs_.size() == 0 ) ) )
	{
		// If output_parameter_unitconvs was never initialized, set it to all 1
		if ( _satellite_parameter_unitconvs_.size() == 0 )
			_satellite_parameter_unitconvs_.resize(
					_satellite_output_parameters_.size(), 1 );

		// Check both vectors are the same size
		if ( _satellite_output_parameters_.size()
				!= _satellite_parameter_unitconvs_.size() )
		{
			if ( !silent )
				std::cerr
						<< "WARNING: satellite_output_parameters and satellite_output_parameter_unitconvs have different sizes.\n"
						<< "No extra parameters will be output.\n";
		}
		else
		{
			// Count up extra columns to print out
			try
			{
				for ( unsigned int i = 0;
						i < _satellite_output_parameters_.size(); i++ )
				{
					if ( _satellite_output_parameters_.at( i ) )
						num_extra_satellite_columns++;
				}
			}
			catch ( std::exception & )
			{
				if ( !silent )
					std::cerr
							<< "WARNING: Error reading satellite_output_parameters.\n"
							<< "No extra parameters will be output.\n";
				num_extra_satellite_columns = 0;
			}
		}
	}

	// As before, for host now
	num_extra_host_columns = 0;
	if ( !( ( _host_output_parameters_.size() == 0 )
			&& ( _host_parameter_unitconvs_.size() == 0 ) ) )
	{
		// If output_parameter_unitconvs was never initialized, set it to all 1
		if ( _host_parameter_unitconvs_.size() == 0 )
			_host_parameter_unitconvs_.resize( _host_output_parameters_.size(),
					1 );

		// Check both vectors are the same size
		if ( _host_output_parameters_.size()
				!= _host_parameter_unitconvs_.size() )
		{
			if ( !silent )
				std::cerr
						<< "WARNING: host_output_parameters and host_output_parameter_unitconvs have different sizes.\n"
						<< "No extra parameters will be output.\n";
		}
		else
		{
			// Count up extra columns to print out
			try
			{
				for ( unsigned int i = 0; i < _host_output_parameters_.size();
						i++ )
				{
					if ( _host_output_parameters_.at( i ) )
						num_extra_host_columns++;
				}
			}
			catch ( std::exception & )
			{
				if ( !silent )
					std::cerr
							<< "WARNING: Error reading host_output_parameters.\n"
							<< "No extra parameters will be output.\n";
				num_extra_host_columns = 0;
			}
		}
	}

	int num_columns = num_columns_base + num_extra_satellite_columns
			+ num_extra_host_columns;

	if ( make_array1d( header, num_columns ) )
		return 1;
	if ( make_array2d( data, num_columns, num_rows ) )
		return 1;

	header[0] = "#";
	header[1] = "t";
	header[2] = "x";
	header[3] = "y";
	header[4] = "z";
	header[5] = "vx";
	header[6] = "vy";
	header[7] = "vz";
	header[8] = "d";
	header[9] = "v";
	header[10] = "m_ret";
	header[11] = "m_frac_lost";
	header[12] = "sim_m_ret";
	header[13] = "sim_m_frac_lost";
	header[14] = "rt";
	header[15] = "rt/rvir";

	if ( num_extra_satellite_columns > 0 )
	{
		int extra_column_counter = 0;
		int num_parameter_names = 0;
		std::vector< std::string > parameter_names( 0 );

		if ( _current_satellite_ptr_->get_parameter_names( num_parameter_names,
				parameter_names ) )
			return 1;

		for ( unsigned int i = 0; i < _satellite_output_parameters_.size();
				i++ )
		{
			if ( _satellite_output_parameters_.at( i ) )
			{
				try
				{
					ss.str( "" );
					ss << "Satellite_" << parameter_names.at( i );
					header[num_columns_base + extra_column_counter] = ss.str();
				}
				catch ( std::exception & )
				{
					ss.str( "" );
					ss << "Unknown_parameter_" << ( extra_column_counter + 1 );
					header[num_columns_base + extra_column_counter] = ss.str();
				}
				extra_column_counter++;
			}
		}
	}

	if ( num_extra_host_columns > 0 )
	{
		int extra_column_counter = 0;
		int num_parameter_names = 0;
		std::vector< std::string > parameter_names( 0 );

		if ( _current_host_ptr_->get_parameter_names( num_parameter_names,
				parameter_names ) )
			return 1;

		for ( unsigned int i = 0; i < _host_output_parameters_.size(); i++ )
		{
			if ( _host_output_parameters_.at( i ) )
			{
				try
				{
					ss.str( "" );
					ss << "Host_" << parameter_names.at( i );
					header[num_columns_base + num_extra_satellite_columns
							+ extra_column_counter] = ss.str();
				}
				catch ( std::exception & )
				{
					ss.str( "" );
					ss << "Unknown_parameter_"
							<< ( num_extra_satellite_columns
									+ extra_column_counter + 1 );
					header[num_columns_base + num_extra_satellite_columns
							+ extra_column_counter] = ss.str();
				}
				extra_column_counter++;
			}
		}
	}

	double t_min_to_use, t_max_to_use;

	if ( _override_t_min_ )
	{
		t_min_to_use = _t_min_override_val_;
	}
	else
	{
		t_min_to_use = _t_min_natural_value_;
	}

	if ( _override_t_max_ )
	{
		t_max_to_use = _t_max_override_val_;
	}
	else
	{
		t_max_to_use = _t_max_natural_value_;
	}

	if ( t_max_to_use <= t_min_to_use )
	{
		if ( !silent )
			std::cerr << "ERROR t_max <= t_min for stripping_orbit_segment!\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	for ( int i = 0; i < num_rows; i++ )
	{
		data[0][i] = "";
		ss.str( "" );
		ss << _phase_output_list_.at( i ).t * unitconv::stoGyr;
		data[1][i] = ss.str();
		ss.str( "" );
		ss << _phase_output_list_.at( i ).x * unitconv::mtokpc;
		data[2][i] = ss.str();
		ss.str( "" );
		ss << _phase_output_list_.at( i ).y * unitconv::mtokpc;
		data[3][i] = ss.str();
		ss.str( "" );
		ss << _phase_output_list_.at( i ).z * unitconv::mtokpc;
		data[4][i] = ss.str();
		ss.str( "" );
		ss
				<< _phase_output_list_.at( i ).vx * unitconv::mtokpc
						/ unitconv::stoGyr;
		data[5][i] = ss.str();
		ss.str( "" );
		ss
				<< _phase_output_list_.at( i ).vy * unitconv::mtokpc
						/ unitconv::stoGyr;
		data[6][i] = ss.str();
		ss.str( "" );
		ss
				<< _phase_output_list_.at( i ).vz * unitconv::mtokpc
						/ unitconv::stoGyr;
		data[7][i] = ss.str();
		ss.str( "" );
		ss
				<< SALTSA::dist3d( _phase_output_list_.at( i ).x,
						_phase_output_list_.at( i ).y,
						_phase_output_list_.at( i ).z ) * unitconv::mtokpc;
		data[8][i] = ss.str();
		ss.str( "" );
		ss
				<< quad_add( _vx_spline_( _phase_output_list_.at( i ).t ),
						_vy_spline_( _phase_output_list_.at( i ).t ),
						_vz_spline_( _phase_output_list_.at( i ).t ) )
						* unitconv::mtokpc / unitconv::stoGyr;
		data[9][i] = ss.str();
		ss.str( "" );
		ss << mret_list.at( i ) * mret_multiplier;
		data[10][i] = ss.str();
		ss.str( "" );
		if ( i > 0 )
			ss
					<< ( 1 - mret_list.at( i ) / mret_list.at( i - 1 ) )
							/ ( _phase_output_list_.at( i ).t
									- _phase_output_list_.at( i - 1 ).t );
		else
			ss << 0;
		data[11][i] = ss.str();
		ss.str( "" );
		ss << _test_mass_spline_( _phase_output_list_.at( i ).t );
		data[12][i] = ss.str();
		ss.str( "" );
		if ( i > 0 )
			ss
					<< ( 1
							- _test_mass_spline_(
									_phase_output_list_.at( i ).t )
									/ safe_d(
											_test_mass_spline_(
													_phase_output_list_.at(
															i - 1 ).t ) ) )
							/ ( _phase_output_list_.at( i ).t
									- _phase_output_list_.at( i - 1 ).t );
		else
			ss << 0;
		data[13][i] = ss.str();
		ss.str( "" );
		ss << _rt_list_.at( i ) * unitconv::mtokpc;
		data[14][i] = ss.str();
		ss.str( "" );
		ss << _rt_ratio_list_.at( i );
		data[15][i] = ss.str();

		if ( num_extra_satellite_columns > 0 )
		{
			int extra_column_counter = 0;
			for ( unsigned int j = 0; j < _satellite_output_parameters_.size();
					j++ )
			{
				if ( _satellite_output_parameters_.at( j ) )
				{
					try
					{
						ss.str( "" );
						ss
								<< ( _satellite_parameter_data_.at( i ).at( j )
										/ safe_d(
												_satellite_parameter_unitconvs_.at(
														j ) ) );
						data[num_columns_base + extra_column_counter][i] =
								ss.str();
					}
					catch ( std::exception & )
					{
						ss.str( "" );
						ss << DBL_MIN;
						data[num_columns_base + extra_column_counter][i] =
								ss.str();
					}
					extra_column_counter++;
				}
			}
		}

		if ( num_extra_host_columns > 0 )
		{
			int extra_column_counter = 0;
			for ( unsigned int j = 0; j < _host_output_parameters_.size();
					j++ )
			{
				if ( _host_output_parameters_.at( j ) )
				{
					try
					{
						ss.str( "" );
						ss
								<< ( _host_parameter_data_.at( i ).at( j )
										/ safe_d(
												_host_parameter_unitconvs_.at(
														j ) ) );
						data[num_columns_base + num_extra_satellite_columns
								+ extra_column_counter][i] = ss.str();
					}
					catch ( std::exception & )
					{
						ss.str( "" );
						ss << DBL_MIN;
						data[num_columns_base + num_extra_satellite_columns
								+ extra_column_counter][i] = ss.str();
					}
					extra_column_counter++;
				}
			}
		}
	}

	if ( include_header )
		return print_table( *out, num_columns, num_rows, header, data, false );
	else
		return print_table( *out, num_columns, num_rows, header, data, true ); // Format for header, but don't print it again
}

const double SALTSA::stripping_orbit_segment::_delta_rho(
		const int index, const double x, const double &t_step,
		const bool silent ) const
{
	if ( index < 1 )
		return 0;
	double result = ( 1. / 6. / Gc * std::pow( 1 + x, _tidal_shocking_power_ )
			* _tidal_shocking_amplification_
			* ( 2 * ( _gabdt_list_[index] * _sum_gabdt_list_[index - 1] )
					+ _gabdt_list_[index] * _gabdt_list_[index] ) );
	if ( isbad( result ) )
	{
		if(!silent)
		{
			std::cerr << "ERROR: Could not calculate delta_rho in stripping_orbit_segment::_delta_rho.\n";
			std::cerr.flush();
		}
		throw std::runtime_error("ERROR: Could not calculate delta_rho in stripping_orbit_segment::_delta_rho.\n");
	}
	return max( result, 0 );
}

const double SALTSA::stripping_orbit_segment::_step_length_factor( const double & v, const double & r ) const
{
	/* Method to determine the appropriate factor for the step length.
	 * It uses the host's virial velocity and radius so that parts of
	 * the orbit with either high velocity or low distance from the host
	 * will take small steps (and the inverse as well).
	 */
	double v_0;
	double r_0;

	// Determine proper scaling from the virial velocities, using the
	// most appropriate value depending on what's been loaded.
	if(_current_host_in_use_)
	{
		v_0 = _current_host_ptr_->vvir();
		r_0 = _current_host_ptr_->rvir();
	}
	else
	{
		if(_host_loaded_)
		{
			v_0 = _init_host_ptr_->vvir();
			r_0 = _init_host_ptr_->rvir();
		}
		else
		{
			v_0 = _v_0_;
			r_0 = _r_0_;
		}
	}

	// Figure out what the step factor would be from each of v and r
	// alone.

	// Low-v -> large step factor (and the inverse), bounded by max and min
	double step_length_factor_v = bound( _step_factor_min_,
			std::pow( std::fabs(v_0 / safe_d( v )), _step_length_power_ ), _step_factor_max_);

	// Low-r -> small step factor (and the inverse), bounded by max and min
	double step_length_factor_r = bound( _step_factor_min_,
			std::pow( std::fabs(r / safe_d( r_0 )), _step_length_power_ ), _step_factor_max_);
	// Use the smaller of the two calculated factors
	return min(step_length_factor_v,step_length_factor_r);
}

const double SALTSA::stripping_orbit_segment::_rvir(
		const int index ) const
{
	try
	{
		return _current_satellite_ptr_->rvir();
	}
	catch ( std::exception & )
	{
		return -1;
	}
}

const int SALTSA::stripping_orbit_segment::get_final_mret(
double & mret, const bool silent ) const
{
	if ( !_calculated_ )
	{
		if ( calc() )
			return UNSPECIFIED_ERROR;
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}
	mret = _init_satellite_ptr_->mvir()*mret_list.back();
	return 0;
}
const int SALTSA::stripping_orbit_segment::get_final_fmret( double & fmret,
		const bool silent ) const
{
	if ( !_calculated_ )
	{
		if ( calc() )
			return UNSPECIFIED_ERROR;
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}
	fmret = mret_list.back() / safe_d(mret_list.front());
	return 0;
}

const int SALTSA::stripping_orbit_segment::get_final_sum_deltarho(
double & final_sum_deltarho, const bool silent ) const
{
	if ( !_calculated_ )
	{
		if ( calc() )
			return UNSPECIFIED_ERROR;
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}
	if(_bad_result_) return UNSPECIFIED_ERROR;
	final_sum_deltarho = _sum_delta_rho_list_.back();
	return 0;
}

const int SALTSA::stripping_orbit_segment::get_final_sum_gabdt(
		gabdt & final_sum_gabdt, const bool silent ) const
{
	if ( !_calculated_ )
	{
		if ( calc() )
			return UNSPECIFIED_ERROR;
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}
	final_sum_gabdt = _sum_gabdt_list_.back();
	return 0;
}

const int SALTSA::stripping_orbit_segment::clone_final_satellite(
		density_profile * & final_satellite_clone, const bool silent ) const
{
	if ( !_calculated_ )
	{
		if ( calc() )
		{
			// We can't calculate. Next most logical option is to use the initial
			// host
			if(_host_loaded_)
			{
				final_satellite_clone = _init_satellite_ptr_->density_profile_clone();
				return UNSPECIFIED_ERROR;
			}
			else
			{
				// If this function is being called, something being assigned is expected,
				// so we'll assign a new tNFW profile as a default
				final_satellite_clone = new SALTSA::tNFW_profile;
				return NOT_SET_UP_ERROR;
			}
		}
	}
	if( _bad_result_ )
	{
		final_satellite_clone = _init_satellite_ptr_->density_profile_clone(); // Sanest option in this case
		return UNSPECIFIED_ERROR;
	}
	else
	{
		final_satellite_clone = _current_satellite_ptr_->density_profile_clone(); // Successful!
		return 0;
	}
}

const int SALTSA::stripping_orbit_segment::clone_final_host(
		density_profile * & final_host_clone, const bool silent ) const
{
	if ( !_calculated_ )
	{
		if ( calc() )
		{
			// We can't calculate. Next most logical option is to use the initial
			// host
			if(_host_loaded_)
			{
				final_host_clone = _init_host_ptr_->density_profile_clone();
				return UNSPECIFIED_ERROR;
			}
			else
			{
				// If this function is being called, something being assigned is expected,
				// so we'll assign a new tNFW profile as a default
				final_host_clone = new SALTSA::tNFW_profile;
				return NOT_SET_UP_ERROR;
			}
		}
	}
	if( _bad_result_ )
	{
		final_host_clone = _init_host_ptr_->density_profile_clone(); // Sanest option in this case
		return UNSPECIFIED_ERROR;
	}
	else
	{
		final_host_clone = _current_host_ptr_->density_profile_clone(); // Successful!
		return 0;
	}
}

#if (1) // Get final data (throws exception on error)

const double SALTSA::stripping_orbit_segment::final_mret() const
{
	double result = -1;

	if ( get_final_mret( result ) )
	{
		std::cerr << "ERROR: Could not calculate in stripping_orbit_segment::final_mret.\n";
		std::cerr.flush();
		throw std::runtime_error("ERROR: Could not calculate in stripping_orbit_segment::final_mret.\n");
	}
	return result;
}

const double SALTSA::stripping_orbit_segment::final_sum_deltarho() const
{
	double result;

	if ( get_final_sum_deltarho( result ) )
	{
		std::cerr << "ERROR: Could not calculate in stripping_orbit_segment::final_sum_deltarho.\n";
		std::cerr.flush();
		throw std::runtime_error("ERROR: Could not calculate in stripping_orbit_segment::final_sum_deltarho.\n");
		result = -1;
	}
	return result;
}

const double SALTSA::stripping_orbit_segment::final_fmret() const
{
	double result = -1;

	if ( get_final_fmret( result ) )
	{
		std::cerr << "ERROR: Could not calculate in stripping_orbit_segment::final_fmret.\n";
		std::cerr.flush();
		throw std::runtime_error("ERROR: Could not calculate in stripping_orbit_segment::final_fmret.\n");
	}
	return result;
}

const SALTSA::gabdt SALTSA::stripping_orbit_segment::final_sum_gabdt() const
{
	gabdt result;

	if ( get_final_sum_gabdt( result ) )
	{
		std::cerr << "ERROR: Could not calculate in stripping_orbit_segment::final_sum_gabdt.\n";
		std::cerr.flush();
		throw std::runtime_error("ERROR: Could not calculate in stripping_orbit_segment::final_sum_gabdt.\n");
	}
	return result;
}

const SALTSA::density_profile * SALTSA::stripping_orbit_segment::final_satellite() const
{
	if ( !_calculated_ )
	{
		if(_satellite_loaded_)
		{
			if ( calc() )
			{
				return _init_satellite_ptr_;
			}
		}
		else
		{
			throw std::runtime_error("ERROR: Attempt to call stripping_orbit::final_satellite() without init_satellite assigned.\n");
			return 0;
		}
	}
	if( _bad_result_ )
		return _init_satellite_ptr_; // Sanest option in this case
	else
		return _current_satellite_ptr_;
}

const SALTSA::density_profile * SALTSA::stripping_orbit_segment::final_host() const
{
	if ( !_calculated_ )
	{
		if(_host_loaded_)
		{
			if ( calc() )
			{
				return _init_host_ptr_;
			}
		}
		else
		{
			throw std::runtime_error("ERROR: Attempt to call stripping_orbit::final_host() without init_host assigned.\n");
			return 0;
		}
	}
	if( _bad_result_ )
		return _init_host_ptr_; // Sanest option in this case
	else
		return _current_host_ptr_;
}

#endif

#endif // end SALTSA::stripping_orbit_segment class function definitions

// SALTSA::solve_rt_sd_function class method implementations
#if (1)

const int SALTSA::solve_rt_grid_function::operator()(
		const double & in_param,
		double & out_param, const bool silent ) const
{
	double r;
	double delta_M;
	r = std::fabs( in_param );

	delta_M = sum_rho * 4. / 3. * pi * std::pow( r, 3 );

	if ( r == 0 )
	{
		out_param = DBL_MAX;
	}
	else
	{
		out_param = std::fabs(
				Gc * ( satellite_ptr->enc_mass( r ) - delta_M )
						/ safe_d(
								( omega * omega + Daccel ) * std::pow( r, 3 ) )
						- 1 );
	}
	return 0;
}

SALTSA::solve_rt_grid_function::solve_rt_grid_function()
{
	omega = 0;
	Daccel = 0;
	sum_rho = 0;
	satellite_ptr = 0;
	return;
}
SALTSA::solve_rt_grid_function::solve_rt_grid_function(
		const double init_omega, const density_profile *init_satellite,
		const double init_Daccel, const double init_sum_rho )
{
	omega = init_omega;
	satellite_ptr = init_satellite;
	Daccel = init_Daccel;
	sum_rho = init_sum_rho;
	return;
}

const int SALTSA::solve_rt_it_function::operator()(
		const double & in_param,
		double & out_param, const bool silent ) const
{
	double r = 0;
	double delta_M = 0;
	double r3 = 0;

	r = std::fabs( in_param );

	delta_M = sum_rho * 4. / 3. * pi * std::pow( r, 3 );

	r3 = Gc * ( satellite_ptr->enc_mass( r ) - delta_M )
			/ safe_d( omega * omega + Daccel );
	if ( r3 <= 0 )
	{
		out_param = r * 0.9;
	}
	else
	{
		out_param = std::pow( r3, 1. / 3. );
	}
	return 0;
}

SALTSA::solve_rt_it_function::solve_rt_it_function()
{
	omega = 0;
	satellite_ptr = 0;
	Daccel = 0;
	sum_rho = 0;
	return;
}
SALTSA::solve_rt_it_function::solve_rt_it_function(
		const double init_omega, const density_profile *init_satellite,
		const double init_Daccel, const double init_sum_rho )
{
	omega = init_omega;
	satellite_ptr = init_satellite;
	Daccel = init_Daccel;
	sum_rho = init_sum_rho;
	return;
}

#endif // end SALTSA::solve_rt_sd_function class function definitions

// SALTSA::gabdt class method implementations
#if (1)

SALTSA::gabdt::gabdt( void )
{
	clear();
}

SALTSA::gabdt::gabdt( const gabdt & other_gabdt )
{
	_is_cached_ = other_gabdt._is_cached_;
	_host_ptr_ = other_gabdt._host_ptr_;

	_x_ = other_gabdt._x_;
	_y_ = other_gabdt._y_;
	_z_ = other_gabdt._z_;
	_r_ = other_gabdt._r_;
	_dt_ = other_gabdt._dt_;

	_dv_ = other_gabdt._dv_;
}

SALTSA::gabdt::gabdt( const density_profile *new_host,
		const double &new_x, const double &new_y,
		const double &new_z, const double &new_dt )
{
	set( new_host, new_x, new_y, new_z, new_dt );
}

SALTSA::gabdt::~gabdt( void )
{
}

const int SALTSA::gabdt::clear()
{
	_is_cached_ = false;
	_x_ = 0;
	_y_ = 0;
	_z_ = 0;
	_r_ = 0;
	_dt_ = 0;
	_dv_.clear();
	return 0;
}

const int SALTSA::gabdt::set( const density_profile *new_host,
		const double &new_x, const double &new_y,
		const double &new_z, const double &new_dt )
{
	_is_cached_ = false;
	_host_ptr_ = new_host;
	_dt_ = new_dt;
	_dv_.clear();
	return 0;
}
const int SALTSA::gabdt::set_host_ptr( const density_profile *new_host )
{
	_is_cached_ = false;
	_host_ptr_ = new_host;
	return 0;
}
const int SALTSA::gabdt::set_pos( const double &new_x,
		const double &new_y, const double &new_z )
{
	_is_cached_ = false;
	_x_ = new_x;
	_y_ = new_y;
	_z_ = new_z;
	_r_ = dist3d( _x_, _y_, _z_ );
	_dv_.clear();
	return 0;
}
const int SALTSA::gabdt::set_dt( const double &new_dt )
{
	if ( _is_cached_ )
	{
		// Multiply dv by ratio of new_dt to old
		for ( int i = 0; i < 3; i++ )
			for ( int j = 0; j < 3; j++ )
			{
				_dv_[i][j] *= new_dt / _dt_;
			}
	}
	_dt_ = new_dt;
	return 0;
}
const int SALTSA::gabdt::calc_dv( const bool silent ) const
{
	gabdt_function gabdtf;
	gabdtf.host_ptr = _host_ptr_;
	unsigned int num_in_params = 3, num_out_params = 3;
	std::vector< double > in_params( num_in_params, 0 ), out_params(
			num_out_params, 0 );
	make_array2d( _dv_, num_out_params, num_in_params );

	in_params[0] = _x_;
	in_params[1] = _y_;
	in_params[2] = _z_;

	if ( differentiate( &gabdtf, num_in_params, in_params, num_out_params,
			out_params, _dv_ ) )
	{
		if ( !silent )
			std::cerr << "ERROR: Cannot differentiate in gabdt::calc_dv().\n";
		return 1;
	}

	// Multiply by dt
	for ( unsigned int i = 0; i < num_out_params; i++ )
		for ( unsigned int j = 0; j < num_in_params; j++ )
		{
			_dv_[i][j] *= _dt_;
		}

	_is_cached_ = true;

	return 0;
}
const int SALTSA::gabdt::override_zero()
{
	make_array2d( _dv_, 3, 3 );
	for ( int x_i = 0; x_i < 3; x_i++ )
	{
		for ( int y_i = 0; y_i < 3; y_i++ )
		{
			_dv_[x_i][y_i] = 0;
		}
	}
	_is_cached_ = true;
	return 0;
}
const SALTSA::density_profile * SALTSA::gabdt::host() const
{
	return _host_ptr_;
}
const double SALTSA::gabdt::x() const
{
	return _x_;
}
const double SALTSA::gabdt::y() const
{
	return _y_;
}
const double SALTSA::gabdt::z() const
{
	return _z_;
}
const double SALTSA::gabdt::r() const
{
	return _r_;
}
const std::vector< std::vector< double > > SALTSA::gabdt::dv() const
{
	if ( !_is_cached_ )
	{
		if ( calc_dv() )
			_dv_.clear();
	}
	return _dv_;
}
const double SALTSA::gabdt::dv( const int x_i, const int y_i ) const
{
	return dv()[x_i][y_i];
}

SALTSA::gabdt & SALTSA::gabdt::operator=( const gabdt &other_gabdt )
{
	if ( &other_gabdt != this )
	{

		_is_cached_ = other_gabdt._is_cached_;
		_host_ptr_ = other_gabdt._host_ptr_;

		_x_ = other_gabdt._x_;
		_y_ = other_gabdt._y_;
		_z_ = other_gabdt._z_;
		_r_ = other_gabdt._r_;
		_dt_ = other_gabdt._dt_;

		_dv_ = other_gabdt._dv_;

	}
	return *this;
}

const double SALTSA::gabdt::operator*( const gabdt & other_gabdt ) const // "Dot-product" operator
{
	if ( !_is_cached_ )
		if ( calc_dv() )
			_dv_.clear();
	double result = 0;
	for ( int x_i = 0; x_i < 3; x_i++ )
	{
		for ( int y_i = 0; y_i < 3; y_i++ )
		{
			result += dv().at(x_i).at(y_i)*other_gabdt.dv().at(x_i).at(y_i);
		}
	}
	return result;
}
SALTSA::gabdt & SALTSA::gabdt::operator+=( const gabdt & other_gabdt )
{
	if ( !_is_cached_ )
		if ( calc_dv() )
		{
			_dv_.clear();
			return *this;
		}

	for ( int x_i = 0; x_i < 3; x_i++ )
	{
		for ( int y_i = 0; y_i < 3; y_i++ )
		{
			_dv_[x_i][y_i] += other_gabdt.dv()[x_i][y_i];
		}
	}
	return *this;
}
SALTSA::gabdt SALTSA::gabdt::operator+( const gabdt & other_gabdt ) const
{
	if ( !_is_cached_ )
		if ( calc_dv() )
		{
			_dv_.clear();
			return *this;
		}
	gabdt result = gabdt( *this );
	result += other_gabdt;
	return result;
}

SALTSA::gabdt & SALTSA::gabdt::operator*=( const double scale_fraction )
{
	if ( !_is_cached_ )
		if ( calc_dv() )
		{
			_dv_.clear();
			return *this;
		}

	for ( int x_i = 0; x_i < 3; x_i++ )
	{
		for ( int y_i = 0; y_i < 3; y_i++ )
		{
			_dv_[x_i][y_i] *= scale_fraction;
		}
	}
	return *this;

}

SALTSA::gabdt SALTSA::gabdt::operator*( const double scale_fraction ) const
{
	if ( !_is_cached_ )
		if ( calc_dv() )
		{
			_dv_.clear();
			return *this;
		}
	gabdt result = gabdt( *this );
	result *= scale_fraction;
	return result;
}

#endif // end SALTSA::gabdt class function definitions

// SALTSA::gabdt_function class method implementations
#if (1)

SALTSA::gabdt_function::gabdt_function()
{
	host_ptr = 0;
}

const int SALTSA::gabdt_function::operator()(
		const std::vector< double > & in_params,
		std::vector< double > & out_params, const bool silent ) const
{
	if ( in_params.size() != 3 )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: in_params.size() and num_in_params must == 3 in gabdt_function.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	double R = dist3d( in_params[0], in_params[1], in_params[2] );

	unsigned int num_out_params = 3;
	if ( make_array( out_params, num_out_params ) )
		return 1;
	if ( R == 0 )
	{
		// Special handling for this case, returning a zero result
		for ( unsigned int i = 0; i < num_out_params; i++ )
		{
			out_params[i] = 0;
		}

	}
	else
	{
		double accel_mag = host_ptr->accel( R );

		for ( unsigned int i = 0; i < num_out_params; i++ )
		{
			out_params[i] = in_params[i] / R * accel_mag;
		}
	}
	return 0;
}

#endif // end SALTSA::gabdt_function class function definitions

#endif // end class function definitions

// SALTSA function definitions
#if (1)
const double SALTSA::stripping_orbit_segment::tidal_strip_retained( const density_profile *host_group,
		const density_profile *satellite, const double &r,
		const double &vr, const double &vt,
		const double &time_step, const double &sum_rho ) const
{
	double new_rt;
	double inst_orbital_period, hm_period, stripping_period;
	double mass_frac_retained, mass_frac_lost_total;
	inst_orbital_period = 2 * pi * r / safe_d(vt);
	hm_period = satellite->othm();

	if(isbad(hm_period) or (hm_period<=0))
		stripping_period = inst_orbital_period;
	else
		stripping_period = inst_orbital_period *
			std::pow(hm_period/safe_d(inst_orbital_period), _tidal_stripping_deceleration_);

	new_rt = get_rt( host_group, satellite, r, vr, vt, time_step, sum_rho );

	if ( !( new_rt > 0 ) )
		mass_frac_lost_total = 0;
	else
	{
		mass_frac_lost_total = max(
				1 - satellite->enc_mass( new_rt ) / safe_d(satellite->mtot()), 0 );
	}
	mass_frac_retained = max(
			min(
					1
							- mass_frac_lost_total * time_step / stripping_period
									* _tidal_stripping_amplification_, 1 ), 0 );

	return mass_frac_retained;
}

const double SALTSA::stripping_orbit_segment::get_rt( const density_profile *host_group,
		const density_profile *satellite, const double &r,
		const double &vr, const double &vt,
		const double &time_step, const double &sum_rho,
		const bool silent ) const
{
	double omega;
	double new_rt, old_rt;
	double max_rt = 0;

	omega = vt / safe_d( r );

	// Check for null case
	if ( satellite->mtot() <= 0 )
	{
		return 0;
	}

	solve_rt_grid_function rt_grid_solver( omega, satellite,
			host_group->Daccel( r ), sum_rho );
	solve_rt_it_function rt_it_solver( omega, satellite,
			host_group->Daccel( r ), sum_rho );

	if ( satellite->mtot() <= 0 )
		return 0;

	old_rt = satellite->rt();

	// First, we try solving iteratively
	new_rt = solve_iterate( &rt_it_solver, old_rt, 1, 0.0001, 1000, true );
	if ( ( new_rt == 0 ) || ( isbad( new_rt ) ) )
	{
		// Iteratively didn't work, so we go to the grid option

		max_rt = 2 * default_tau_factor * satellite->rvir();
		old_rt = new_rt;
		if ( solve_grid( &rt_grid_solver, 1, 0., max_rt, 100, 0., new_rt ) )
		{
			if ( !silent )
				std::cerr << "WARNING: Could not solve rt.\n";
			new_rt = 0; // Most likely value in the case where we can't solve it
		}
	}

	return ( max( new_rt, 0 ) );
}

#endif // end SALTSA function definitions
